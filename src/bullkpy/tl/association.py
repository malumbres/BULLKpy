from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp

from scipy.stats import (
    chi2_contingency,
    kruskal,
    f_oneway,
    mannwhitneyu,
    ttest_ind,
)

try:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
except Exception:  # pragma: no cover
    adjusted_rand_score = None
    normalized_mutual_info_score = None

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
    out_idx = np.empty_like(order)
    out_idx[order] = np.arange(n)
    out[m] = q[out_idx]
    return out


def _as_array_layer(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _get_gene_vector(adata: ad.AnnData, gene: str, layer: str | None) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names.")
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    j = int(adata.var_names.get_loc(gene))
    if sp.issparse(X):
        v = X[:, j].toarray().ravel()
    else:
        v = np.asarray(X[:, j], dtype=float).ravel()
    return v


def _cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 2 or y.size < 2:
        return np.nan
    sx = np.nanstd(x, ddof=1)
    sy = np.nanstd(y, ddof=1)
    s = np.sqrt(((x.size - 1) * sx**2 + (y.size - 1) * sy**2) / max(x.size + y.size - 2, 1))
    if s == 0:
        return np.nan
    return float((np.nanmean(x) - np.nanmean(y)) / s)


def _eta2_from_anova(groups: list[np.ndarray]) -> float:
    # eta^2 = SS_between / SS_total
    allv = np.concatenate(groups)
    grand = np.mean(allv)
    ss_total = np.sum((allv - grand) ** 2)
    ss_between = 0.0
    for g in groups:
        if g.size == 0:
            continue
        ss_between += g.size * (np.mean(g) - grand) ** 2
    if ss_total == 0:
        return np.nan
    return float(ss_between / ss_total)


def _epsilon2_from_kruskal(H: float, n: int, k: int) -> float:
    # Common effect size for Kruskal–Wallis
    # epsilon^2 = (H - k + 1) / (n - k)
    denom = (n - k)
    if denom <= 0:
        return np.nan
    return float((H - k + 1) / denom)


def _cramers_v(table: np.ndarray, chi2: float) -> float:
    n = table.sum()
    if n == 0:
        return np.nan
    r, c = table.shape
    denom = n * max(min(r - 1, c - 1), 1)
    return float(np.sqrt(chi2 / denom))


def _is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s)


def _as_categorical(s: pd.Series) -> pd.Categorical:
    if pd.api.types.is_categorical_dtype(s):
        return s.astype("category")
    return pd.Categorical(s.astype(str))


# -----------------------------
# A) gene ↔ categorical
# -----------------------------
def gene_categorical_association(
    adata: ad.AnnData,
    *,
    genes: list[str] | None = None,
    groupby: str,
    layer: str | None = "log1p_cpm",
    method: Literal["auto", "anova", "kruskal", "ttest", "wilcoxon"] = "auto",
    effect_size: Literal["eta2", "epsilon2", "log2fc", "cohen_d", "none"] = "epsilon2",
    adjust: Literal["fdr_bh", "none"] = "fdr_bh",
    min_group_size: int = 2,
) -> pd.DataFrame:
    """
    Association of gene expression with a categorical obs column.
    Returns a tidy dataframe (one row per gene).
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    cats = _as_categorical(adata.obs[groupby])
    groups = list(pd.Categorical(cats).categories)

    if genes is None:
        genes = list(map(str, adata.var_names))
    else:
        genes = [str(g) for g in genes]

    info(f"gene_categorical_association: {len(genes)} genes vs {groupby} ({len(groups)} groups)")

    rows = []
    for g in genes:
        y = _get_gene_vector(adata, g, layer=layer)
        df = pd.DataFrame({"y": y, "grp": cats})
        df = df.dropna(subset=["y", "grp"])
        if df.empty:
            continue

        # group vectors
        gv = []
        means = {}
        for lv in groups:
            v = df.loc[df["grp"].astype(str) == str(lv), "y"].to_numpy(dtype=float)
            if v.size >= 1:
                means[f"mean_{lv}"] = float(np.mean(v))
            else:
                means[f"mean_{lv}"] = np.nan
            gv.append(v)

        # filter by size
        ok = [v for v in gv if v.size >= min_group_size]
        k_eff = len(ok)
        if k_eff < 2:
            stat = np.nan
            pval = np.nan
        else:
            m = method
            if m == "auto":
                m = "kruskal" if k_eff > 2 else "wilcoxon"

            if m == "anova":
                stat, pval = f_oneway(*ok)
            elif m == "kruskal":
                stat, pval = kruskal(*ok)
            elif m == "ttest":
                if k_eff != 2:
                    stat, pval = np.nan, np.nan
                else:
                    stat, pval = ttest_ind(ok[0], ok[1], equal_var=False, nan_policy="omit")
            elif m == "wilcoxon":
                if k_eff != 2:
                    stat, pval = np.nan, np.nan
                else:
                    # mann-whitney is the usual "wilcoxon rank-sum" for independent groups
                    stat, pval = mannwhitneyu(ok[0], ok[1], alternative="two-sided")
            else:
                raise ValueError(f"Unknown method='{method}'")

        # effect
        eff = np.nan
        if effect_size != "none" and k_eff >= 2 and np.isfinite(stat):
            if effect_size == "eta2" and (method in ("anova", "auto")):
                eff = _eta2_from_anova(ok)
            elif effect_size == "epsilon2" and (method in ("kruskal", "auto")):
                n = int(sum(v.size for v in ok))
                eff = _epsilon2_from_kruskal(float(stat), n=n, k=k_eff)
            elif effect_size == "cohen_d" and k_eff == 2:
                eff = _cohen_d(ok[0], ok[1])
            elif effect_size == "log2fc" and k_eff == 2:
                # log2(mean(group1)+eps) - log2(mean(group2)+eps)
                eps = 1e-9
                eff = float(np.log2(np.mean(ok[0]) + eps) - np.log2(np.mean(ok[1]) + eps))

        row = {
            "gene": g,
            "groupby": groupby,
            "n_groups": int(len(groups)),
            "statistic": float(stat) if np.isfinite(stat) else np.nan,
            "pval": float(pval) if np.isfinite(pval) else np.nan,
            "effect": eff,
            "method": (method if method != "auto" else ("kruskal" if len(groups) > 2 else "wilcoxon")),
            "effect_size": effect_size,
        }
        row.update(means)
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    if adjust == "fdr_bh":
        out["qval"] = _bh_fdr(out["pval"].to_numpy())
    else:
        out["qval"] = np.nan

    out = out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)
    return out


# -----------------------------
# B) numeric obs ↔ categorical
# -----------------------------
def obs_categorical_association(
    adata: ad.AnnData,
    *,
    obs_keys: list[str] | None = None,
    groupby: str,
    method: Literal["auto", "anova", "kruskal", "ttest", "wilcoxon"] = "auto",
    effect_size: Literal["eta2", "epsilon2", "cohen_d", "none"] = "epsilon2",
    adjust: Literal["fdr_bh", "none"] = "fdr_bh",
    min_group_size: int = 2,
) -> pd.DataFrame:
    """
    Association of numeric obs columns with a categorical obs column.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    cats = _as_categorical(adata.obs[groupby])
    groups = list(pd.Categorical(cats).categories)

    if obs_keys is None:
        obs_keys = [c for c in adata.obs.columns if _is_numeric_series(adata.obs[c])]
    else:
        obs_keys = [str(k) for k in obs_keys]

    info(f"obs_categorical_association: {len(obs_keys)} numeric obs vs {groupby}")

    rows = []
    for k in obs_keys:
        if k not in adata.obs.columns:
            continue
        s = adata.obs[k]
        if not _is_numeric_series(s):
            continue

        df = pd.DataFrame({"y": pd.to_numeric(s, errors="coerce"), "grp": cats})
        df = df.dropna(subset=["y", "grp"])
        if df.empty:
            continue

        gv = []
        means = {}
        for lv in groups:
            v = df.loc[df["grp"].astype(str) == str(lv), "y"].to_numpy(dtype=float)
            means[f"mean_{lv}"] = float(np.mean(v)) if v.size else np.nan
            gv.append(v)

        ok = [v for v in gv if v.size >= min_group_size]
        k_eff = len(ok)
        if k_eff < 2:
            stat = np.nan
            pval = np.nan
        else:
            m = method
            if m == "auto":
                m = "kruskal" if k_eff > 2 else "wilcoxon"

            if m == "anova":
                stat, pval = f_oneway(*ok)
            elif m == "kruskal":
                stat, pval = kruskal(*ok)
            elif m == "ttest":
                if k_eff != 2:
                    stat, pval = np.nan, np.nan
                else:
                    stat, pval = ttest_ind(ok[0], ok[1], equal_var=False, nan_policy="omit")
            elif m == "wilcoxon":
                if k_eff != 2:
                    stat, pval = np.nan, np.nan
                else:
                    stat, pval = mannwhitneyu(ok[0], ok[1], alternative="two-sided")
            else:
                raise ValueError(f"Unknown method='{method}'")

        eff = np.nan
        if effect_size != "none" and k_eff >= 2 and np.isfinite(stat):
            if effect_size == "eta2" and (method in ("anova", "auto")):
                eff = _eta2_from_anova(ok)
            elif effect_size == "epsilon2" and (method in ("kruskal", "auto")):
                n = int(sum(v.size for v in ok))
                eff = _epsilon2_from_kruskal(float(stat), n=n, k=k_eff)
            elif effect_size == "cohen_d" and k_eff == 2:
                eff = _cohen_d(ok[0], ok[1])

        row = {
            "obs": k,
            "groupby": groupby,
            "n_groups": int(len(groups)),
            "statistic": float(stat) if np.isfinite(stat) else np.nan,
            "pval": float(pval) if np.isfinite(pval) else np.nan,
            "qval": np.nan,
            "effect": eff,
            "method": (method if method != "auto" else ("kruskal" if len(groups) > 2 else "wilcoxon")),
            "effect_size": effect_size,
        }
        row.update(means)
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    if adjust == "fdr_bh":
        out["qval"] = _bh_fdr(out["pval"].to_numpy())
    out = out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)
    return out


# -----------------------------
# C) categorical ↔ categorical
# -----------------------------
def categorical_association(
    adata: ad.AnnData,
    *,
    key1: str,
    key2: str,
    metrics: Iterable[str] = ("chi2", "cramers_v", "ari", "nmi"),
    dropna: bool = True,
) -> dict:
    """
    Categorical association between two obs columns.
    Returns dict with contingency table + requested metrics.
    """
    if key1 not in adata.obs.columns:
        raise KeyError(f"'{key1}' not in adata.obs")
    if key2 not in adata.obs.columns:
        raise KeyError(f"'{key2}' not in adata.obs")

    s1 = adata.obs[key1]
    s2 = adata.obs[key2]

    df = pd.DataFrame({key1: s1, key2: s2})
    if dropna:
        df = df.dropna(subset=[key1, key2])

    a = df[key1].astype(str)
    b = df[key2].astype(str)

    tab = pd.crosstab(a, b)
    table = tab.to_numpy()

    out: dict = {"table": tab}

    metrics = [m.lower() for m in metrics]

    if "chi2" in metrics or "cramers_v" in metrics:
        chi2, p, dof, _ = chi2_contingency(tab, correction=False)
        if "chi2" in metrics:
            out["chi2"] = {"statistic": float(chi2), "pval": float(p), "dof": int(dof)}
        if "cramers_v" in metrics:
            out["cramers_v"] = _cramers_v(table, float(chi2))

    # ARI / NMI (clustering similarity)
    if ("ari" in metrics or "nmi" in metrics) and (adjusted_rand_score is None or normalized_mutual_info_score is None):
        warn("sklearn not available: cannot compute ARI/NMI.")
    else:
        if "ari" in metrics:
            out["ari"] = float(adjusted_rand_score(a, b))
        if "nmi" in metrics:
            out["nmi"] = float(normalized_mutual_info_score(a, b))

    return out


# -----------------------------
# Dispatcher: association(x, y)
# -----------------------------
def association(
    adata: ad.AnnData,
    *,
    x: str,
    y: str,
    layer: str | None = "log1p_cpm",
    method: str = "auto",
) -> object:
    """
    Unified association dispatcher.
      - x gene, y categorical -> gene_categorical_association (single-gene)
      - x numeric obs, y categorical -> obs_categorical_association (single-obs)
      - x categorical, y categorical -> categorical_association
    """
    x_is_gene = x in adata.var_names
    y_is_gene = y in adata.var_names
    x_in_obs = x in adata.obs.columns
    y_in_obs = y in adata.obs.columns

    if x_is_gene and y_in_obs:
        # gene vs categorical (or numeric, but here we focus categorical)
        if _is_numeric_series(adata.obs[y]):
            raise ValueError("gene↔numeric handled elsewhere (use your numeric correlation utilities).")
        return gene_categorical_association(adata, genes=[x], groupby=y, layer=layer, method=method)

    if x_in_obs and y_in_obs:
        if _is_numeric_series(adata.obs[x]) and (not _is_numeric_series(adata.obs[y])):
            return obs_categorical_association(adata, obs_keys=[x], groupby=y, method=method)
        if (not _is_numeric_series(adata.obs[x])) and (not _is_numeric_series(adata.obs[y])):
            return categorical_association(adata, key1=x, key2=y)
        raise ValueError("This dispatcher covers numeric↔categorical and categorical↔categorical only.")

    if y_is_gene and x_in_obs:
        # allow swapped order
        if _is_numeric_series(adata.obs[x]):
            raise ValueError("gene↔numeric handled elsewhere.")
        return gene_categorical_association(adata, genes=[y], groupby=x, layer=layer, method=method)

    raise KeyError("Could not resolve x/y as gene or obs columns.")