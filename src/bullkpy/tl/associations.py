from __future__ import annotations

from typing import Literal, Sequence, Mapping, Any, Optional, Tuple
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from anndata import AnnData
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

try:
    from scipy import stats
except Exception:  # pragma: no cover
    stats = None

from ..logging import info, warn



# -----------------------------
# Utilities
# -----------------------------
def _require_scipy() -> None:
    if stats is None:
        raise ImportError("associations requires scipy. Please install scipy.")


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR correction. Keeps NaNs."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)

    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out

    p0 = p[ok]
    n = p0.size
    order = np.argsort(p0)
    ranked = p0[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    out_ok = np.empty_like(p0)
    out_ok[order] = q
    out[ok] = out_ok
    return out


def _get_matrix(adata: ad.AnnData, layer: str | None = None) -> Any:
    """Return X/layer as (n_obs x n_vars), sparse or dense."""
    if layer is not None and layer in adata.layers:
        return adata.layers[layer]
    return adata.X


def _get_gene_vector(adata: ad.AnnData, gene: str, layer: str | None = None) -> np.ndarray:
    """Return expression vector for one gene (n_obs,)."""
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names.")
    X = _get_matrix(adata, layer=layer)
    j = int(adata.var_names.get_loc(gene))
    if sp.issparse(X):
        return np.asarray(X[:, j].toarray()).ravel().astype(float)
    return np.asarray(X[:, j], dtype=float).ravel()


def _as_categorical(s: pd.Series) -> pd.Categorical:
    if pd.api.types.is_categorical_dtype(s):
        return s.astype("category")
    return pd.Categorical(s.astype(str))


def _is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s)


def _eta2_from_anova(groups: list[np.ndarray]) -> float:
    """eta^2 = SS_between / SS_total"""
    allv = np.concatenate([g for g in groups if g.size > 0], axis=0)
    if allv.size == 0:
        return np.nan
    grand = float(np.mean(allv))
    ss_total = float(np.sum((allv - grand) ** 2))
    if ss_total == 0:
        return np.nan

    ss_between = 0.0
    for g in groups:
        if g.size == 0:
            continue
        ss_between += g.size * float((np.mean(g) - grand) ** 2)
    return float(ss_between / ss_total)


def _epsilon2_from_kruskal(H: float, n: int, k: int) -> float:
    """Common effect size for Kruskal–Wallis: epsilon^2 = (H - k + 1) / (n - k)."""
    denom = (n - k)
    if denom <= 0:
        return np.nan
    return float((H - k + 1.0) / denom)


# -----------------------------
# A) Gene ↔ categorical obs (global test across all groups)
# -----------------------------
def gene_categorical_association(
    adata: ad.AnnData,
    *,
    groupby: str,
    genes: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    test: Literal["auto", "kruskal", "anova"] = "auto",
    effect_size: Literal["auto", "epsilon2", "eta2", "none"] = "auto",
    min_group_size: int = 2,
    adjust: Literal["bh", "none"] = "bh",
) -> pd.DataFrame:
    """
    Scan genes for association with a categorical obs column (multi-class).

    This answers: “Which genes vary across the categories in adata.obs[groupby]?”

    Parameters
    ----------
    adata
        AnnData with samples in `.obs` and genes in `.var`.
    groupby
        Categorical column in `adata.obs` defining groups (e.g. subtype, cluster).
    genes
        Optional subset of genes. If None, uses all `adata.var_names`.
    layer
        Expression layer to use. If None or missing, uses `adata.X`.
        Typical: "log1p_cpm".
    test
        - "kruskal": Kruskal–Wallis (non-parametric; robust default)
        - "anova": one-way ANOVA (parametric)
        - "auto": uses "kruskal"
    effect_size
        - "epsilon2": for Kruskal–Wallis
        - "eta2": for ANOVA
        - "auto": epsilon2 if kruskal/auto, else eta2
        - "none": skip effect
    min_group_size
        Minimum samples per group to include that group in the test.
    adjust
        Multiple testing correction:
        - "bh": Benjamini–Hochberg FDR
        - "none": no correction

    Returns
    -------
    DataFrame with one row per gene, columns:
        ['gene','groupby','n_groups','test','statistic','pval','qval','effect']
    plus mean per group columns: mean_<group>.
    """
    _require_scipy()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    cats = _as_categorical(adata.obs[groupby])
    levels = list(pd.Categorical(cats).categories)

    genes_use = list(map(str, adata.var_names)) if genes is None else [str(g) for g in genes]
    if len(genes_use) == 0:
        raise ValueError("No genes provided.")

    test_use = "kruskal" if test in ("auto", "kruskal") else "anova"
    eff_use = effect_size
    if eff_use == "auto":
        eff_use = "epsilon2" if test_use == "kruskal" else "eta2"

    info(f"gene_categorical_association: {len(genes_use)} genes vs '{groupby}' ({len(levels)} groups), test={test_use}")

    rows: list[dict[str, Any]] = []
    for g in genes_use:
        y = _get_gene_vector(adata, g, layer=layer)
        df = pd.DataFrame({"y": y, "grp": cats})
        df = df.dropna(subset=["y", "grp"])

        # group vectors (keep group order stable)
        gv: list[np.ndarray] = []
        means: dict[str, float] = {}
        for lv in levels:
            v = df.loc[df["grp"].astype(str) == str(lv), "y"].to_numpy(dtype=float)
            means[f"mean_{lv}"] = float(np.mean(v)) if v.size else np.nan
            if v.size >= int(min_group_size):
                gv.append(v)

        if len(gv) < 2:
            stat_val, p_val, eff = np.nan, np.nan, np.nan
        else:
            if test_use == "kruskal":
                stat_val, p_val = stats.kruskal(*gv)
                eff = (
                    _epsilon2_from_kruskal(float(stat_val), n=int(sum(v.size for v in gv)), k=int(len(gv)))
                    if eff_use == "epsilon2"
                    else np.nan
                )
            else:
                stat_val, p_val = stats.f_oneway(*gv)
                eff = _eta2_from_anova(gv) if eff_use == "eta2" else np.nan

        row = {
            "gene": str(g),
            "groupby": str(groupby),
            "n_groups": int(len(levels)),
            "test": test_use,
            "statistic": float(stat_val) if np.isfinite(stat_val) else np.nan,
            "pval": float(p_val) if np.isfinite(p_val) else np.nan,
            "effect": float(eff) if np.isfinite(eff) else np.nan,
        }
        row.update(means)
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["qval"] = _bh_fdr(out["pval"].to_numpy()) if adjust == "bh" else np.nan
    out = out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)
    return out


# -----------------------------
# B) Fast “rank_genes_groups”-like for one group vs reference/rest
# -----------------------------
def rank_genes_groups_fast(
    adata: ad.AnnData,
    *,
    groupby: str,
    group: str,
    reference: str | None = None,  # if None -> rest
    layer: str | None = "log1p_cpm",
    genes: Sequence[str] | None = None,
    method: Literal["mwu", "ttest"] = "mwu",
    effect: Literal["delta_mean", "delta_median"] = "delta_mean",
    adjust: Literal["bh", "none"] = "bh",
    min_n: int = 2,
) -> pd.DataFrame:
    """
    Fast DE-like scan for one group vs reference/rest.

    This answers: “Which genes are most associated with THIS category label?”

    - method="mwu": Mann–Whitney U (two-sided) — robust default
    - method="ttest": Welch t-test (unequal variance)

    Returns a DataFrame with:
        gene, effect, pval, qval, mean_group, mean_ref, median_group, median_ref
    """
    _require_scipy()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    grp = adata.obs[groupby].astype(str)
    m1 = grp.eq(str(group)).to_numpy()

    if reference is None:
        m2 = ~m1
        ref_name = "rest"
    else:
        ref_name = str(reference)
        m2 = grp.eq(ref_name).to_numpy()

    if int(m1.sum()) < int(min_n) or int(m2.sum()) < int(min_n):
        raise ValueError(f"Not enough samples in '{group}' and '{ref_name}' for rank_genes_groups_fast.")

    genes_use = list(map(str, adata.var_names)) if genes is None else [str(g) for g in genes]
    if len(genes_use) == 0:
        raise ValueError("No genes provided.")

    X = _get_matrix(adata, layer=layer)
    idx = adata.var_names.get_indexer(genes_use)
    if np.any(idx < 0):
        miss = [genes_use[i] for i in np.where(idx < 0)[0][:10]]
        raise KeyError(f"Genes not in adata.var_names (first 10): {miss}")

    M = X[:, idx]
    if sp.issparse(M):
        M = M.toarray()
    else:
        M = np.asarray(M, dtype=float)

    pvals = np.full(len(genes_use), np.nan, dtype=float)
    effs = np.full(len(genes_use), np.nan, dtype=float)
    mean1 = np.full(len(genes_use), np.nan, dtype=float)
    mean2 = np.full(len(genes_use), np.nan, dtype=float)
    med1 = np.full(len(genes_use), np.nan, dtype=float)
    med2 = np.full(len(genes_use), np.nan, dtype=float)

    info(f"rank_genes_groups_fast: '{group}' vs '{ref_name}' ({len(genes_use)} genes), method={method}, effect={effect}")

    for k in range(len(genes_use)):
        x1 = M[m1, k]
        x2 = M[m2, k]

        mean1[k] = float(np.nanmean(x1))
        mean2[k] = float(np.nanmean(x2))
        med1[k] = float(np.nanmedian(x1))
        med2[k] = float(np.nanmedian(x2))

        # test
        try:
            if method == "ttest":
                p = stats.ttest_ind(x1, x2, equal_var=False, nan_policy="omit").pvalue
            else:
                p = stats.mannwhitneyu(x1, x2, alternative="two-sided").pvalue
        except Exception:
            p = np.nan
        pvals[k] = float(p) if np.isfinite(p) else np.nan

        # effect
        if effect == "delta_median":
            effs[k] = med1[k] - med2[k]
        else:
            effs[k] = mean1[k] - mean2[k]

    qvals = _bh_fdr(pvals) if adjust == "bh" else np.full_like(pvals, np.nan, dtype=float)

    df = pd.DataFrame(
        {
            "gene": genes_use,
            "groupby": str(groupby),
            "group": str(group),
            "reference": str(ref_name),
            "effect": effs,
            "pval": pvals,
            "qval": qvals,
            "mean_group": mean1,
            "mean_ref": mean2,
            "median_group": med1,
            "median_ref": med2,
        }
    )
    df = df.sort_values("qval", na_position="last").reset_index(drop=True)
    return df


# -----------------------------
# C) Categorical ↔ categorical (chi-square)
# -----------------------------
def categorical_association(
    adata: ad.AnnData,
    *,
    key1: str,
    key2: str,
) -> dict[str, Any]:
    """
    Association between two categorical obs columns using chi-square contingency test.

    Returns a dict with:
        - chi2, dof, pval
        - table (contingency DataFrame)
        - expected (numpy array)
    """
    _require_scipy()

    if key1 not in adata.obs.columns:
        raise KeyError(f"'{key1}' not found in adata.obs")
    if key2 not in adata.obs.columns:
        raise KeyError(f"'{key2}' not found in adata.obs")

    s1 = adata.obs[key1].astype(str)
    s2 = adata.obs[key2].astype(str)
    table = pd.crosstab(s1, s2)

    if table.size == 0:
        raise ValueError("Empty contingency table.")

    chi2, pval, dof, expected = stats.chi2_contingency(table.to_numpy())
    return {
        "key1": key1,
        "key2": key2,
        "chi2": float(chi2),
        "dof": int(dof),
        "pval": float(pval),
        "table": table,
        "expected": expected,
    }


# -----------------------------
# D) Convenience: store results in adata.uns
# -----------------------------
def store_gene_categorical_association(
    adata: ad.AnnData,
    *,
    groupby: str,
    key: str = "associations",
    subkey: str | None = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Run `gene_categorical_association` and store results in `adata.uns[key][subkey]`.

    If subkey is None, uses f"gene~{groupby}".
    """
    subkey = subkey or f"gene~{groupby}"
    res = gene_categorical_association(adata, groupby=groupby, **kwargs)

    adata.uns.setdefault(key, {})
    adata.uns[key][subkey] = res
    info(f"Stored gene_categorical_association results in adata.uns['{key}']['{subkey}']")

    return res


def store_rank_genes_groups_fast(
    adata: ad.AnnData,
    *,
    groupby: str,
    group: str,
    reference: str | None = None,
    key: str = "rank_genes_groups_fast",
    subkey: str | None = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Run `rank_genes_groups_fast` and store results in `adata.uns[key][subkey]`.

    If subkey is None, uses f"{groupby}:{group}_vs_{reference or 'rest'}".
    """
    ref_name = reference if reference is not None else "rest"
    subkey = subkey or f"{groupby}:{group}_vs_{ref_name}"

    res = rank_genes_groups_fast(adata, groupby=groupby, group=group, reference=reference, **kwargs)

    adata.uns.setdefault(key, {})
    adata.uns[key][subkey] = res
    info(f"Stored rank_genes_groups_fast results in adata.uns['{key}']['{subkey}']")

    return res


# -----------------------------
# E) Very small dispatcher (optional)
# -----------------------------
def association(
    adata: ad.AnnData,
    *,
    x: str,
    y: str,
    layer: str | None = "log1p_cpm",
) -> Any:
    """
    Minimal association dispatcher (categorical focus).

    - gene vs categorical obs -> rank_genes_groups_fast-like is not applicable (needs a target group),
      so we run the global scan (gene_categorical_association) for that categorical obs.
    - categorical vs categorical -> categorical_association

    For numeric correlations, keep using correlations.py utilities.
    """
    x_is_gene = x in adata.var_names
    y_is_gene = y in adata.var_names
    x_in_obs = x in adata.obs.columns
    y_in_obs = y in adata.obs.columns

    if x_is_gene and y_in_obs:
        if _is_numeric_series(adata.obs[y]):
            raise ValueError("Use correlations utilities for gene↔numeric obs.")
        return gene_categorical_association(adata, groupby=y, genes=[x], layer=layer)

    if y_is_gene and x_in_obs:
        if _is_numeric_series(adata.obs[x]):
            raise ValueError("Use correlations utilities for gene↔numeric obs.")
        return gene_categorical_association(adata, groupby=x, genes=[y], layer=layer)

    if x_in_obs and y_in_obs:
        sx_num = _is_numeric_series(adata.obs[x])
        sy_num = _is_numeric_series(adata.obs[y])
        if (not sx_num) and (not sy_num):
            return categorical_association(adata, key1=x, key2=y)
        raise ValueError("This dispatcher only covers categorical↔categorical and gene↔categorical.")

    raise KeyError("Could not resolve x/y as gene or obs columns.")


def posthoc_per_gene(
    adata: ad.AnnData,
    *,
    gene: str,
    groupby: str,
    layer: str | None = "log1p_cpm",
    method: Literal["mwu", "ttest"] = "mwu",
    adjust: Literal["bh", "none"] = "bh",
    min_n: int = 2,
) -> pd.DataFrame:
    """
    Pairwise post-hoc comparisons for ONE gene across all categories.

    Use AFTER selecting a gene of interest.

    Returns
    -------
    DataFrame with columns:
        ['gene','group1','group2','effect','pval','qval']
    """
    _require_scipy()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")
    if gene not in adata.var_names:
        raise KeyError(f"gene='{gene}' not found in adata.var_names")

    grp = adata.obs[groupby].astype(str)
    cats = list(pd.Categorical(grp).categories)

    y = _get_gene_vector(adata, gene, layer=layer)

    rows = []
    for i, g1 in enumerate(cats):
        for g2 in cats[i + 1 :]:
            x1 = y[grp == g1]
            x2 = y[grp == g2]

            if x1.size < min_n or x2.size < min_n:
                continue

            try:
                if method == "ttest":
                    p = stats.ttest_ind(x1, x2, equal_var=False, nan_policy="omit").pvalue
                else:
                    p = stats.mannwhitneyu(x1, x2, alternative="two-sided").pvalue
            except Exception:
                p = np.nan

            eff = float(np.nanmean(x1) - np.nanmean(x2))

            rows.append(
                {
                    "gene": gene,
                    "group1": g1,
                    "group2": g2,
                    "effect": eff,
                    "pval": float(p) if np.isfinite(p) else np.nan,
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    df["qval"] = _bh_fdr(df["pval"].to_numpy()) if adjust == "bh" else np.nan
    return df.sort_values("qval", na_position="last").reset_index(drop=True)



def gene_metadata_association_scan(
    adata: AnnData,
    *,
    metadata_key: str,
    layer: str = "log1p_cpm",
    genes: Optional[Sequence[str]] = None,
    groupby: Optional[str] = None,   # optional: run separately per tumor type
    method: str = "spearman",        # "spearman" (v1)
    adjust: str = "fdr_bh",
    out_key: str = "gene_meta_assoc",
) -> pd.DataFrame:
    """
    Scan gene–metadata associations (v1: Spearman correlation for continuous metadata).

    Returns a tidy DataFrame with columns:
      group (if groupby), gene, n, rho, pval, qval

    Stores in `adata.uns[out_key]`.
    """
    if metadata_key not in adata.obs:
        raise KeyError(f"{metadata_key!r} not found in adata.obs")
    if layer not in adata.layers:
        raise KeyError(f"layer {layer!r} not found in adata.layers")

    y_all = adata.obs[metadata_key].to_numpy(dtype=float)

    X = adata.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)

    var_names = np.asarray(adata.var_names)
    if genes is not None:
        genes = [g for g in genes if g in adata.var_names]
        idx = np.isin(var_names, genes)
        X = X[:, idx]
        var_names = var_names[idx]

    def run_block(y, Xblk, genes_blk, group_label=None):
        rows = []
        for j, g in enumerate(genes_blk):
            x = Xblk[:, j].astype(float)
            ok = np.isfinite(x) & np.isfinite(y)
            n = int(ok.sum())
            if n < 10:
                rows.append({"group": group_label, "gene": g, "n": n, "rho": np.nan, "pval": np.nan})
                continue
            rho, p = spearmanr(x[ok], y[ok])
            rows.append({"group": group_label, "gene": g, "n": n, "rho": float(rho), "pval": float(p)})
        df = pd.DataFrame(rows)
        # adjust p-values within this block
        pvals = df["pval"].to_numpy(dtype=float)
        q = np.full_like(pvals, np.nan, dtype=float)
        okp = np.isfinite(pvals)
        if okp.any():
            q[okp] = multipletests(pvals[okp], method=adjust)[1]
        df["qval"] = q
        return df

    if groupby is None:
        out = run_block(y_all, X, var_names, group_label=None)
        out = out.drop(columns=["group"])
    else:
        if groupby not in adata.obs:
            raise KeyError(f"{groupby!r} not found in adata.obs")
        out_parts = []
        for grp, idx in adata.obs.groupby(groupby).groups.items():
            idx = np.asarray(list(idx))
            out_parts.append(run_block(y_all[idx], X[idx, :], var_names, group_label=str(grp)))
        out = pd.concat(out_parts, ignore_index=True)

    adata.uns[out_key] = out
    return out





def signature_axis(
    adata: AnnData,
    *,
    pos_key: str,
    neg_key: str,
    standardize: str = "global",  # "global" | "groupby"
    groupby: Optional[str] = None,
    out_key: str = "signature_axis",
) -> None:
    """
    Create a two-sided axis from two signature score columns:
        axis = z(pos) - z(neg)

    Stores in `adata.obs[out_key]`.
    """
    obs = adata.obs
    if pos_key not in obs or neg_key not in obs:
        raise KeyError("pos_key and/or neg_key not found in adata.obs")

    pos = obs[pos_key].to_numpy(dtype=float)
    neg = obs[neg_key].to_numpy(dtype=float)

    axis = np.full_like(pos, np.nan, dtype=float)

    if standardize == "global":
        def z(x):
            mu = np.nanmean(x)
            sd = np.nanstd(x)
            return (x - mu) / (sd if sd != 0 else np.nan)

        axis = z(pos) - z(neg)

    elif standardize == "groupby":
        if groupby is None or groupby not in obs:
            raise KeyError("groupby must be provided and exist in adata.obs when standardize='groupby'")
        axis[:] = np.nan
        for grp, idx in obs.groupby(groupby).groups.items():
            idx = np.asarray(list(idx))
            p = pos[idx]; n = neg[idx]
            def z_local(x):
                mu = np.nanmean(x)
                sd = np.nanstd(x)
                return (x - mu) / (sd if sd != 0 else np.nan)
            axis[idx] = z_local(p) - z_local(n)
    else:
        raise ValueError("standardize must be 'global' or 'groupby'")

    adata.obs[out_key] = axis


def quantile_groups(
    adata: AnnData,
    *,
    key: str,
    q: Tuple[float, float] = (0.25, 0.75),
    labels: Tuple[str, str, str] = ("low", "mid", "high"),
    out_key: Optional[str] = None,
) -> None:
    """
    Assign samples into low/mid/high based on quantiles of a continuous `.obs[key]`.

    Stores categorical groups in `.obs[out_key or f'{key}_qgrp']`.
    """
    if key not in adata.obs:
        raise KeyError(f"{key!r} not found in adata.obs")

    x = adata.obs[key].to_numpy(dtype=float)
    lo, hi = np.nanquantile(x, q)

    grp = np.array([labels[1]] * len(x), dtype=object)
    grp[np.isfinite(x) & (x <= lo)] = labels[0]
    grp[np.isfinite(x) & (x >= hi)] = labels[2]
    grp[~np.isfinite(x)] = pd.NA

    target = out_key or f"{key}_qgrp"
    adata.obs[target] = pd.Categorical(grp, categories=list(labels), ordered=True)