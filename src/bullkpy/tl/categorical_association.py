from __future__ import annotations

from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp

from scipy.stats import kruskal, f_oneway, chi2_contingency

import anndata as ad

from ..logging import info, warn


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    m = np.isfinite(p)
    if m.sum() == 0:
        return out
    pv = p[m]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[m] = q[inv]
    return out


def _as_str_cat(s: pd.Series) -> pd.Categorical:
    return pd.Categorical(s.astype(str))


def _epsilon2_from_kruskal(H: float, n: int, k: int) -> float:
    # common effect size for Kruskal–Wallis
    if n <= 1 or k <= 1:
        return np.nan
    return max(0.0, (H - k + 1) / (n - k))


def _eta2_from_anova(F: float, df_between: int, df_within: int) -> float:
    # approximate eta^2 from F and dof (one-way)
    if df_between <= 0 or df_within <= 0:
        return np.nan
    return (F * df_between) / (F * df_between + df_within)


def _get_X_for_genes(
    adata: ad.AnnData,
    genes: Sequence[str],
    *,
    layer: str | None = "log1p_cpm",
) -> tuple[np.ndarray, list[str]]:
    genes = [str(g) for g in genes]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found in adata.var_names (first 10): {missing[:10]}")
    idx = [adata.var_names.get_loc(g) for g in genes]

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    M = X[:, idx]
    if sp.issparse(M):
        M = M.toarray()
    M = np.asarray(M, dtype=float)
    return M, genes


def gene_categorical_association(
    adata: ad.AnnData,
    *,
    groupby: str,
    genes: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    method: Literal["kruskal", "anova"] = "kruskal",
    effect_size: Literal["epsilon2", "eta2"] | None = "epsilon2",
    min_group_size: int = 2,
    adjust: Literal["fdr_bh", "none"] = "fdr_bh",
) -> pd.DataFrame:
    """
    A) Association between gene expression (numeric) and a categorical obs column.

    Returns tidy df with columns:
      groupby, gene, statistic, pval, qval, effect, n_groups, n, means (optional columns)
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    g = _as_str_cat(adata.obs[groupby])
    cats = list(g.categories)
    if len(cats) < 2:
        raise ValueError(f"groupby='{groupby}' has <2 groups.")

    if genes is None:
        genes = list(adata.var_names)

    M, genes = _get_X_for_genes(adata, genes, layer=layer)
    info(f"gene_categorical_association: {len(genes)} genes vs {groupby} ({len(cats)} groups), method={method}")

    rows = []
    for j, gene in enumerate(genes):
        y = M[:, j]
        # split by group
        groups = []
        means = {}
        n_total = 0
        for c in cats:
            vals = y[(g == c)]
            vals = vals[np.isfinite(vals)]
            if vals.size >= min_group_size:
                groups.append(vals)
                means[str(c)] = float(np.mean(vals)) if vals.size else np.nan
                n_total += vals.size

        if len(groups) < 2:
            continue

        if method == "kruskal":
            stat, p = kruskal(*groups)
            eff = _epsilon2_from_kruskal(float(stat), int(n_total), int(len(groups))) if effect_size == "epsilon2" else np.nan
        elif method == "anova":
            stat, p = f_oneway(*groups)
            # df between = k-1, within = n-k
            dfb = len(groups) - 1
            dfw = n_total - len(groups)
            eff = _eta2_from_anova(float(stat), int(dfb), int(dfw)) if effect_size == "eta2" else np.nan
        else:
            raise ValueError(f"Unknown method='{method}'")

        rows.append(
            {
                "groupby": str(groupby),
                "gene": str(gene),
                "statistic": float(stat),
                "pval": float(p),
                "effect": float(eff) if np.isfinite(eff) else np.nan,
                "n_groups": int(len(groups)),
                "n": int(n_total),
                # optional: keep means for reference (can be used later)
                "group_means": means,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["qval"] = _bh_fdr(out["pval"].to_numpy()) if adjust == "fdr_bh" else np.nan
    out = out.sort_values(["qval", "pval"], ascending=True, na_position="last").reset_index(drop=True)
    return out


def obs_categorical_association(
    adata: ad.AnnData,
    *,
    groupby: str,
    obs_keys: Sequence[str] | None = None,
    method: Literal["kruskal", "anova"] = "kruskal",
    effect_size: Literal["epsilon2", "eta2"] | None = "epsilon2",
    min_group_size: int = 2,
    adjust: Literal["fdr_bh", "none"] = "fdr_bh",
) -> pd.DataFrame:
    """
    B) Association between numeric obs columns and a categorical obs column.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    g = _as_str_cat(adata.obs[groupby])
    cats = list(g.categories)
    if len(cats) < 2:
        raise ValueError(f"groupby='{groupby}' has <2 groups.")

    if obs_keys is None:
        # all numeric obs
        obs_keys = list(adata.obs.select_dtypes(include=[np.number]).columns)

    obs_keys = [str(k) for k in obs_keys if k in adata.obs.columns]
    if len(obs_keys) == 0:
        raise ValueError("No obs_keys found / selected.")

    info(f"obs_categorical_association: {len(obs_keys)} numeric obs vs {groupby} ({len(cats)} groups), method={method}")

    rows = []
    for k in obs_keys:
        y = pd.to_numeric(adata.obs[k], errors="coerce").to_numpy(dtype=float)

        groups = []
        means = {}
        n_total = 0
        for c in cats:
            vals = y[(g == c)]
            vals = vals[np.isfinite(vals)]
            if vals.size >= min_group_size:
                groups.append(vals)
                means[str(c)] = float(np.mean(vals)) if vals.size else np.nan
                n_total += vals.size

        if len(groups) < 2:
            continue

        if method == "kruskal":
            stat, p = kruskal(*groups)
            eff = _epsilon2_from_kruskal(float(stat), int(n_total), int(len(groups))) if effect_size == "epsilon2" else np.nan
        elif method == "anova":
            stat, p = f_oneway(*groups)
            dfb = len(groups) - 1
            dfw = n_total - len(groups)
            eff = _eta2_from_anova(float(stat), int(dfb), int(dfw)) if effect_size == "eta2" else np.nan
        else:
            raise ValueError(f"Unknown method='{method}'")

        rows.append(
            {
                "groupby": str(groupby),
                "obs": str(k),
                "statistic": float(stat),
                "pval": float(p),
                "effect": float(eff) if np.isfinite(eff) else np.nan,
                "n_groups": int(len(groups)),
                "n": int(n_total),
                "group_means": means,
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["qval"] = _bh_fdr(out["pval"].to_numpy()) if adjust == "fdr_bh" else np.nan
    out = out.sort_values(["qval", "pval"], ascending=True, na_position="last").reset_index(drop=True)
    return out


def cat_cat_association(
    adata: ad.AnnData,
    *,
    key1: str,
    key2: str,
    adjust: Literal["fdr_bh", "none"] = "none",
) -> pd.DataFrame:
    """
    C) Association between two categorical obs columns.

    Returns:
      key1, key2, chi2, dof, pval, qval, cramers_v, n
    """
    if key1 not in adata.obs.columns:
        raise KeyError(f"'{key1}' not in adata.obs")
    if key2 not in adata.obs.columns:
        raise KeyError(f"'{key2}' not in adata.obs")

    s1 = adata.obs[key1].astype(str)
    s2 = adata.obs[key2].astype(str)
    tab = pd.crosstab(s1, s2)

    chi2, p, dof, _ = chi2_contingency(tab.to_numpy())
    n = float(tab.to_numpy().sum())
    r, k = tab.shape

    # Cramér’s V
    if n <= 0 or min(r - 1, k - 1) <= 0:
        v = np.nan
    else:
        v = np.sqrt((chi2 / n) / min(r - 1, k - 1))

    out = pd.DataFrame(
        [
            {
                "key1": str(key1),
                "key2": str(key2),
                "chi2": float(chi2),
                "dof": int(dof),
                "pval": float(p),
                "cramers_v": float(v) if np.isfinite(v) else np.nan,
                "n": int(n),
            }
        ]
    )
    out["qval"] = _bh_fdr(out["pval"].to_numpy()) if adjust == "fdr_bh" else np.nan
    return out