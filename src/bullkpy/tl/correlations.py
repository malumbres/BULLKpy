from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Iterable, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp

from scipy.stats import t as t_dist
from scipy.stats import spearmanr

import anndata as ad

try:
    import seaborn as sns  # optional
except Exception:  # pragma: no cover
    sns = None

import matplotlib.pyplot as plt

from ..logging import info, warn


Method = Literal["pearson", "spearman"]
BatchMode = Literal["none", "residualize", "within_batch", "meta"]


# ---------------------------
# Core helpers
# ---------------------------

def _as_list(x: str | Sequence[str] | None) -> list[str]:
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


def _get_layer_matrix(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _get_gene_vector(adata: ad.AnnData, gene: str, *, layer: str | None) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names")
    X = _get_layer_matrix(adata, layer)
    j = adata.var_names.get_loc(gene)
    return np.asarray(X[:, j], dtype=float).ravel()


def _numeric_obs_df(adata: ad.AnnData, *, obs_keys: list[str] | None = None) -> pd.DataFrame:
    if obs_keys is None:
        # auto-pick numeric columns
        df = adata.obs.copy()
        num_cols = []
        for c in df.columns:
            s = df[c]
            # pandas numeric?
            if pd.api.types.is_numeric_dtype(s):
                num_cols.append(c)
        return df[num_cols].copy()
    else:
        missing = [k for k in obs_keys if k not in adata.obs.columns]
        if missing:
            raise KeyError(f"Obs keys not found: {missing}")
        df = adata.obs[obs_keys].copy()
        # ensure numeric
        for c in df.columns:
            if not pd.api.types.is_numeric_dtype(df[c]):
                raise TypeError(f"Obs column '{c}' is not numeric (dtype={df[c].dtype})")
        return df


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR correction."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)

    m = np.isfinite(p)
    if not np.any(m):
        return q

    pv = p[m]
    n = pv.size
    order = np.argsort(pv)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)

    qv = pv * (n / ranks)
    # monotone
    qv_sorted = qv[order]
    qv_sorted = np.minimum.accumulate(qv_sorted[::-1])[::-1]
    qv[order] = np.clip(qv_sorted, 0.0, 1.0)

    q[m] = qv
    return q


def _corr_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return np.nan, np.nan, n

    xa, ya = x[m], y[m]
    # centered
    xa = xa - xa.mean()
    ya = ya - ya.mean()

    denom = np.sqrt(np.sum(xa * xa) * np.sum(ya * ya))
    if denom == 0:
        return np.nan, np.nan, n

    r = float(np.sum(xa * ya) / denom)

    # p-value from t-stat
    df = n - 2
    if not np.isfinite(r) or df <= 0 or abs(r) >= 1:
        p = 0.0 if abs(r) == 1 else np.nan
    else:
        t = r * np.sqrt(df / (1.0 - r * r))
        p = float(2.0 * t_dist.sf(np.abs(t), df))
    return r, p, n


def _corr_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return np.nan, np.nan, n
    r, p = spearmanr(x[m], y[m], nan_policy="omit")
    return float(r), float(p), n


def _corr(x: np.ndarray, y: np.ndarray, *, method: Method) -> tuple[float, float, int]:
    if method == "pearson":
        return _corr_pearson(x, y)
    elif method == "spearman":
        return _corr_spearman(x, y)
    else:
        raise ValueError(f"Unknown method: {method}")


def _design_matrix_from_obs(
    adata: ad.AnnData,
    covariates: list[str],
    *,
    add_intercept: bool = True,
) -> np.ndarray:
    """
    Build a design matrix from obs columns (numeric + categorical one-hot).
    """
    if len(covariates) == 0:
        return np.ones((adata.n_obs, 1), dtype=float) if add_intercept else np.empty((adata.n_obs, 0), dtype=float)

    parts = []
    for c in covariates:
        if c not in adata.obs.columns:
            raise KeyError(f"Covariate '{c}' not in adata.obs")
        s = adata.obs[c]

        if pd.api.types.is_numeric_dtype(s):
            v = pd.to_numeric(s, errors="coerce").to_numpy(dtype=float).reshape(-1, 1)
            parts.append(v)
        else:
            # categorical / object → one-hot
            cat = pd.Categorical(s.astype(str))
            dummies = pd.get_dummies(cat, drop_first=True)  # avoid collinearity
            parts.append(dummies.to_numpy(dtype=float))

    X = np.concatenate(parts, axis=1) if len(parts) else np.empty((adata.n_obs, 0), dtype=float)
    if add_intercept:
        X = np.concatenate([np.ones((adata.n_obs, 1), dtype=float), X], axis=1)
    return X


def _residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """
    Residualize y with respect to design matrix X via least squares.
    Handles NaNs by fitting on finite rows and returning NaNs elsewhere.
    """
    y = np.asarray(y, dtype=float).ravel()
    if X.size == 0:
        return y.copy()

    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    out = np.full_like(y, np.nan, dtype=float)
    if m.sum() < X.shape[1] + 2:
        return out

    Xm = X[m, :]
    ym = y[m]
    beta, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
    out[m] = ym - Xm @ beta
    return out


def _batch_groups(adata: ad.AnnData, batch_key: str) -> pd.Series:
    if batch_key not in adata.obs.columns:
        raise KeyError(f"batch_key='{batch_key}' not in adata.obs")
    return adata.obs[batch_key].astype("category")


def _within_batch_zscore(v: np.ndarray, batch: pd.Series) -> np.ndarray:
    v = np.asarray(v, dtype=float).ravel()
    out = np.full_like(v, np.nan, dtype=float)
    for b in batch.cat.categories:
        m = (batch == b).to_numpy()
        x = v[m]
        mm = np.isfinite(x)
        if mm.sum() < 2:
            continue
        mu = np.nanmean(x)
        sd = np.nanstd(x, ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        out[m] = (x - mu) / sd
    return out


def _fisher_meta_cor(r_list: list[float], n_list: list[int]) -> tuple[float, float]:
    """
    Combine correlations via Fisher z (weighted by n-3).
    Returns (r_meta, z_se) where z_se is SE in Fisher-z space.
    """
    rs = np.asarray(r_list, dtype=float)
    ns = np.asarray(n_list, dtype=float)

    m = np.isfinite(rs) & np.isfinite(ns) & (ns >= 4) & (np.abs(rs) < 1)
    if not np.any(m):
        return np.nan, np.nan

    rs = rs[m]
    ns = ns[m]
    w = ns - 3.0
    z = np.arctanh(rs)
    zbar = np.sum(w * z) / np.sum(w)
    se = np.sqrt(1.0 / np.sum(w))
    r_meta = float(np.tanh(zbar))
    return r_meta, float(se)


def _corr_batch_aware(
    x: np.ndarray,
    y: np.ndarray,
    *,
    adata: ad.AnnData,
    method: Method,
    batch_key: str | None,
    batch_mode: BatchMode,
    covariates: list[str] | None = None,
) -> tuple[float, float, int]:
    """
    Compute correlation under batch-aware mode.

    - none: standard correlation
    - residualize: regress out batch_key (+covariates) from x and y, then correlate residuals
    - within_batch: z-score x and y within each batch, then correlate pooled
    - meta: correlate within each batch, then Fisher-z combine, approximate p-value from z/se
    """
    covariates = covariates or []

    if batch_mode == "none" or batch_key is None:
        return _corr(x, y, method=method)

    batch = _batch_groups(adata, batch_key)

    if batch_mode == "residualize":
        X = _design_matrix_from_obs(adata, covariates=[batch_key] + covariates, add_intercept=True)
        xr = _residualize(x, X)
        yr = _residualize(y, X)
        return _corr(xr, yr, method=method)

    if batch_mode == "within_batch":
        xz = _within_batch_zscore(x, batch)
        yz = _within_batch_zscore(y, batch)
        return _corr(xz, yz, method=method)

    if batch_mode == "meta":
        r_list: list[float] = []
        p_list: list[float] = []
        n_list: list[int] = []
        for b in batch.cat.categories:
            m = (batch == b).to_numpy()
            rb, pb, nb = _corr(np.asarray(x)[m], np.asarray(y)[m], method=method)
            r_list.append(rb)
            p_list.append(pb)
            n_list.append(nb)

        r_meta, se = _fisher_meta_cor(r_list, n_list)
        if not np.isfinite(r_meta) or not np.isfinite(se) or se <= 0:
            return np.nan, np.nan, int(np.nansum(n_list))

        # z-stat in Fisher space
        z = np.arctanh(r_meta) / se
        p = float(2.0 * (1.0 - 0.5 * (1.0 + np.math.erf(np.abs(z) / np.sqrt(2.0)))))  # 2*sf(|z|) approx
        return r_meta, p, int(np.nansum(n_list))

    raise ValueError(f"Unknown batch_mode='{batch_mode}'")


# ---------------------------
# A) Strongest gene–gene pairs
# ---------------------------

def top_gene_gene_correlations(
    adata: ad.AnnData,
    *,
    genes: list[str] | None = None,
    layer: str | None = "log1p_cpm",
    method: Method = "pearson",
    top_n: int = 200,
    min_abs_r: float | None = None,
    use_abs: bool = True,
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
) -> pd.DataFrame:
    """
    Find strongest gene–gene correlations among a gene set (or all genes if provided).
    Returns long table with (gene1, gene2, r, pval, qval, n).
    """
    X = _get_layer_matrix(adata, layer)
    if genes is None:
        raise ValueError("For safety, provide genes=... for top_gene_gene_correlations (all genes is too big).")
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found: {missing[:10]}")

    idx = [adata.var_names.get_loc(g) for g in genes]
    M = np.asarray(X[:, idx], dtype=float)  # n_obs x n_genes
    names = [adata.var_names[i] for i in idx]

    rows = []
    pvals = []

    # O(G^2) loop; OK for moderate gene panels
    G = len(names)
    for i in range(G):
        xi = M[:, i]
        for j in range(i + 1, G):
            yj = M[:, j]
            r, p, n = _corr_batch_aware(
                xi, yj, adata=adata, method=method,
                batch_key=batch_key, batch_mode=batch_mode,
                covariates=covariates,
            )
            if min_abs_r is not None and np.isfinite(r) and abs(r) < float(min_abs_r):
                continue
            rows.append((names[i], names[j], r, p, n))
            pvals.append(p)

    df = pd.DataFrame(rows, columns=["gene1", "gene2", "r", "pval", "n"])
    df["qval"] = _bh_fdr(df["pval"].to_numpy())

    score = np.abs(df["r"].to_numpy()) if use_abs else df["r"].to_numpy()
    df = df[np.isfinite(score)].copy()
    if df.empty:
        return df

    k = min(int(top_n), df.shape[0])
    top_idx = np.argpartition(score[np.isfinite(score)], -k)[-k:]
    # safer: just sort
    df = df.sort_values("r", key=lambda s: np.abs(s), ascending=False).head(k).reset_index(drop=True)

    df["method"] = method
    df["batch_key"] = batch_key
    df["batch_mode"] = batch_mode
    return df


# ---------------------------
# D) One gene vs genes
# ---------------------------

def gene_gene_correlations(
    adata: ad.AnnData,
    *,
    gene: str,
    genes: list[str] | None = None,
    layer: str | None = "log1p_cpm",
    method: Method = "pearson",
    top_n: int = 50,
    min_abs_r: float | None = None,
    use_abs: bool = True,
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
) -> pd.DataFrame:
    """
    Correlate one gene with all other genes or a subset.
    Returns table: (gene, r, pval, qval, n).
    """
    x0 = _get_gene_vector(adata, gene, layer=layer)
    X = _get_layer_matrix(adata, layer)

    if genes is None:
        idx = [i for i in range(adata.n_vars) if adata.var_names[i] != gene]
        names = [adata.var_names[i] for i in idx]
    else:
        missing = [g for g in genes if g not in adata.var_names]
        if missing:
            raise KeyError(f"Genes not found: {missing[:10]}")
        idx = [adata.var_names.get_loc(g) for g in genes if g != gene]
        names = [adata.var_names[i] for i in idx]

    rows = []
    for i, nm in zip(idx, names):
        xi = X[:, i]
        r, p, n = _corr_batch_aware(
            x0, xi, adata=adata, method=method,
            batch_key=batch_key, batch_mode=batch_mode,
            covariates=covariates,
        )
        rows.append((nm, r, p, n))

    df = pd.DataFrame(rows, columns=["gene", "r", "pval", "n"])
    df["qval"] = _bh_fdr(df["pval"].to_numpy())

    if min_abs_r is not None:
        df = df[np.abs(df["r"]) >= float(min_abs_r)]

    score = np.abs(df["r"].to_numpy()) if use_abs else df["r"].to_numpy()
    df = df[np.isfinite(score)].copy()
    if df.empty:
        return df

    df = df.sort_values("r", key=lambda s: np.abs(s), ascending=False).head(int(top_n)).reset_index(drop=True)
    df.insert(0, "query_gene", gene)
    df["method"] = method
    df["batch_key"] = batch_key
    df["batch_mode"] = batch_mode
    return df


# ---------------------------
# B) Gene vs numeric obs
# ---------------------------

def top_gene_obs_correlations(
    adata: ad.AnnData,
    *,
    gene: str | list[str],
    obs: str | list[str] | None = None,
    obs_keys: list[str] | None = None,
    layer: str | None = "log1p_cpm",
    method: Method = "pearson",
    top_n: int = 50,
    min_abs_r: float | None = None,
    use_abs: bool = True,
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
) -> pd.DataFrame:
    """
    Correlate one or more genes with numeric obs columns.

    - obs: restrict to these obs columns
    - obs_keys: alternative explicit numeric obs selection
    """
    genes = _as_list(gene)
    for g in genes:
        if g not in adata.var_names:
            raise KeyError(f"Gene '{g}' not in adata.var_names")

    obs_df = _numeric_obs_df(adata, obs_keys=obs_keys)
    if obs is not None:
        obs_list = _as_list(obs)
        missing = [k for k in obs_list if k not in obs_df.columns]
        if missing:
            raise KeyError(f"Obs keys not found or not numeric: {missing}")
        obs_df = obs_df[obs_list].copy()

    rows = []
    for g in genes:
        xg = _get_gene_vector(adata, g, layer=layer)
        for col in obs_df.columns:
            y = obs_df[col].to_numpy(dtype=float)
            r, p, n = _corr_batch_aware(
                xg, y, adata=adata, method=method,
                batch_key=batch_key, batch_mode=batch_mode,
                covariates=covariates,
            )
            if min_abs_r is not None and np.isfinite(r) and abs(r) < float(min_abs_r):
                continue
            rows.append((g, col, r, p, n))

    df = pd.DataFrame(rows, columns=["gene", "obs", "r", "pval", "n"])
    df["qval"] = _bh_fdr(df["pval"].to_numpy())

    score = np.abs(df["r"].to_numpy()) if use_abs else df["r"].to_numpy()
    df = df[np.isfinite(score)].copy()
    if df.empty:
        return df

    df = df.sort_values("r", key=lambda s: np.abs(s), ascending=False).head(int(top_n)).reset_index(drop=True)
    df["method"] = method
    df["batch_key"] = batch_key
    df["batch_mode"] = batch_mode
    return df


# ---------------------------
# C) Numeric obs × numeric obs
# ---------------------------

def top_obs_obs_correlations(
    adata: ad.AnnData,
    *,
    focus: str | list[str],
    against: list[str] | None = None,
    obs_keys: list[str] | None = None,
    method: Method = "pearson",
    top_n: int = 100,
    min_abs_r: float | None = None,
    use_abs: bool = True,
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
) -> pd.DataFrame:
    """
    Strongest correlations between numeric obs columns.
    If against is None: compare focus vs all other numeric obs.
    Returns long table of pairs.
    """
    obs_df = _numeric_obs_df(adata, obs_keys=obs_keys)
    focus_l = _as_list(focus)

    for k in focus_l:
        if k not in obs_df.columns:
            raise KeyError(f"focus '{k}' not numeric / not found in obs")

    if against is None:
        against_l = [c for c in obs_df.columns if c not in focus_l]
    else:
        against_l = _as_list(against)
        missing = [k for k in against_l if k not in obs_df.columns]
        if missing:
            raise KeyError(f"against keys not numeric / not found: {missing}")

    rows = []
    for f in focus_l:
        x = obs_df[f].to_numpy(dtype=float)
        for a in against_l:
            y = obs_df[a].to_numpy(dtype=float)
            r, p, n = _corr_batch_aware(
                x, y, adata=adata, method=method,
                batch_key=batch_key, batch_mode=batch_mode,
                covariates=covariates,
            )
            if min_abs_r is not None and np.isfinite(r) and abs(r) < float(min_abs_r):
                continue
            rows.append((f, a, r, p, n))

    df = pd.DataFrame(rows, columns=["focus", "against", "r", "pval", "n"])
    df["qval"] = _bh_fdr(df["pval"].to_numpy())

    score = np.abs(df["r"].to_numpy()) if use_abs else df["r"].to_numpy()
    df = df[np.isfinite(score)].copy()
    if df.empty:
        return df

    df = df.sort_values("r", key=lambda s: np.abs(s), ascending=False).head(int(top_n)).reset_index(drop=True)
    df["method"] = method
    df["batch_key"] = batch_key
    df["batch_mode"] = batch_mode
    return df


def obs_obs_corr_matrix(
    adata: ad.AnnData,
    *,
    focus: str | list[str],
    against: list[str] | None = None,
    obs_keys: list[str] | None = None,
    method: Method = "pearson",
    drop_na: bool = True,
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
    return_p: bool = False,
    return_q: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, pd.DataFrame] | tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Correlation matrix between numeric obs columns (focus × against).
    Returns wide matrix (rows=focus, cols=against).
    Optionally returns p-values and q-values matrices.
    """
    obs_df = _numeric_obs_df(adata, obs_keys=obs_keys)
    focus_l = _as_list(focus)

    if against is None:
        against_l = [c for c in obs_df.columns if c not in focus_l]
    else:
        against_l = _as_list(against)

    for k in focus_l + against_l:
        if k not in obs_df.columns:
            raise KeyError(f"'{k}' not found among numeric obs")

    # compute pairwise
    R = pd.DataFrame(index=focus_l, columns=against_l, dtype=float)
    P = pd.DataFrame(index=focus_l, columns=against_l, dtype=float)

    for f in focus_l:
        x = obs_df[f].to_numpy(dtype=float)
        for a in against_l:
            y = obs_df[a].to_numpy(dtype=float)
            r, p, _ = _corr_batch_aware(
                x, y, adata=adata, method=method,
                batch_key=batch_key, batch_mode=batch_mode,
                covariates=covariates,
            )
            R.loc[f, a] = r
            P.loc[f, a] = p

    if not return_p and not return_q:
        return R

    # q-values across all tested entries
    qvals = _bh_fdr(P.to_numpy().ravel())
    Q = pd.DataFrame(qvals.reshape(P.shape), index=P.index, columns=P.columns)

    if return_p and return_q:
        return R, P, Q
    if return_p:
        return R, P
    return R, Q


# ---------------------------
# Partial correlation convenience
# ---------------------------

def partial_corr(
    adata: ad.AnnData,
    *,
    x: np.ndarray,
    y: np.ndarray,
    covariates: list[str],
    method: Method = "pearson",
) -> tuple[float, float, int]:
    """
    Partial correlation: residualize x and y on covariates (obs columns) then correlate.
    """
    X = _design_matrix_from_obs(adata, covariates=covariates, add_intercept=True)
    xr = _residualize(x, X)
    yr = _residualize(y, X)
    return _corr(xr, yr, method=method)


# ---------------------------
# Plotting wrappers
# ---------------------------

def plot_corr_scatter(
    adata: ad.AnnData,
    *,
    x: str,
    y: str,
    hue: str | None = None,
    kind: Literal["obs_obs", "gene_obs"] = "obs_obs",
    layer: str | None = "log1p_cpm",
    method: Method = "pearson",
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
    figsize: tuple[float, float] = (5.5, 4.5),
    s: float = 18.0,
    alpha: float = 0.8,
    cmap: str = "viridis",
    palette: str = "tab20",
    title: str | None = None,
    show: bool = True,
):
    """
    Scatter plot + correlation in title.
    kind="obs_obs": x,y are obs numeric columns
    kind="gene_obs": x is gene, y is obs numeric (or vice versa)
    """
    if kind == "obs_obs":
        if x not in adata.obs.columns or y not in adata.obs.columns:
            raise KeyError("x and y must be obs columns for kind='obs_obs'")
        xv = pd.to_numeric(adata.obs[x], errors="coerce").to_numpy(dtype=float)
        yv = pd.to_numeric(adata.obs[y], errors="coerce").to_numpy(dtype=float)
    elif kind == "gene_obs":
        # allow x=gene, y=obs or vice versa
        if x in adata.var_names and y in adata.obs.columns:
            xv = _get_gene_vector(adata, x, layer=layer)
            yv = pd.to_numeric(adata.obs[y], errors="coerce").to_numpy(dtype=float)
        elif y in adata.var_names and x in adata.obs.columns:
            xv = pd.to_numeric(adata.obs[x], errors="coerce").to_numpy(dtype=float)
            yv = _get_gene_vector(adata, y, layer=layer)
        else:
            raise KeyError("For kind='gene_obs', one of (x,y) must be a gene and the other a numeric obs.")
    else:
        raise ValueError(f"Unknown kind={kind}")

    r, p, n = _corr_batch_aware(
        xv, yv, adata=adata, method=method,
        batch_key=batch_key, batch_mode=batch_mode,
        covariates=covariates,
    )

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    if hue is None:
        ax.scatter(xv, yv, s=s, alpha=alpha, edgecolors="none")
    else:
        if hue not in adata.obs.columns:
            raise KeyError(f"hue='{hue}' not in adata.obs")
        h = adata.obs[hue]
        if pd.api.types.is_numeric_dtype(h):
            vals = pd.to_numeric(h, errors="coerce").to_numpy(dtype=float)
            sc = ax.scatter(xv, yv, c=vals, s=s, alpha=alpha, cmap=cmap, edgecolors="none")
            cb = fig.colorbar(sc, ax=ax, pad=0.01)
            cb.set_label(hue)
        else:
            cat = pd.Categorical(h.astype(str))
            cats = list(cat.categories)
            if sns is not None:
                pal = sns.color_palette(palette, n_colors=len(cats))
                colors = {c: pal[i] for i, c in enumerate(cats)}
            else:
                cm = plt.get_cmap(palette)
                colors = {c: cm(i % cm.N) for i, c in enumerate(cats)}
            for c in cats:
                m = (cat == c)
                ax.scatter(xv[m], yv[m], s=s, alpha=alpha, edgecolors="none", label=str(c), color=colors[c])
            ax.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ttl = title or f"{method} r={r:.3f}, p={p:.2e}, n={n}"
    ax.set_title(ttl)
    if show:
        plt.show()
    return fig, ax


def plot_corr_heatmap(
    mat: pd.DataFrame,
    *,
    figsize: tuple[float, float] = (7.5, 4.5),
    cmap: str = "vlag",
    center: float = 0.0,
    vmin: float | None = None,
    vmax: float | None = None,
    annot: bool = False,
    fmt: str = ".2f",
    title: str | None = None,
    show: bool = True,
):
    """
    Heatmap for a wide correlation matrix (e.g., from obs_obs_corr_matrix).
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    if sns is not None:
        sns.heatmap(mat, ax=ax, cmap=cmap, center=center, vmin=vmin, vmax=vmax, annot=annot, fmt=fmt)
    else:
        im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto")
        fig.colorbar(im, ax=ax)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns, rotation=90)
        ax.set_yticklabels(mat.index)

    if title is not None:
        ax.set_title(title)

    if show:
        plt.show()
    return fig, ax