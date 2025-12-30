from __future__ import annotations
import warnings
from ._nb_utils import _as_dense, deseq2_size_factors, estimate_dispersion_mom, shrink_dispersion_to_trend
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import anndata as ad

from ..logging import info, warn

from pathlib import Path
from typing import Sequence
import matplotlib.pyplot as plt

def de(
    adata: ad.AnnData,
    *,
    groupby: str,
    group: str,
    reference: str,
    method: Literal["welch_ttest", "nb_glm"] = "welch_ttest",
    layer_counts: str = "counts",
    layer_expr: str = "log1p_cpm",
    key_added: str = "de",
    alpha: float = 0.05,
    shrink_dispersion: bool = True,
    prior_df: float = 10.0,
    ) -> None:

    """
    Two-group differential expression.

    Stores results in:
      adata.uns[key_added][f"{groupby}:_{group}_vs_{reference}"] = DataFrame-like dict

    Methods
    -------
    welch_ttest:
        Fast, exploratory; runs Welch t-test on `layer_expr` (default log1p_cpm).
    nb_glm:
        Proper bulk model; NB GLM on raw counts (layer_counts) with log(libsize) offset.
        Requires statsmodels.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    labels = adata.obs[groupby].astype(str)
    g_mask = labels == str(group)
    r_mask = labels == str(reference)

    if g_mask.sum() < 2 or r_mask.sum() < 2:
        raise ValueError(
            f"Need >=2 samples per group. Got {g_mask.sum()} in {group}, {r_mask.sum()} in {reference}"
        )

    contrast_key = f"{groupby}_{group}_vs_{reference}"
    adata.uns.setdefault(key_added, {})

    if method == "welch_ttest":
        res = _de_welch_ttest(adata, g_mask, r_mask, layer=layer_expr)
        adata.uns[key_added][contrast_key] = {
            "method": method,
            "groupby": groupby,
            "group": str(group),
            "reference": str(reference),
            "layer": layer_expr,
            "results": res,
        }
        info(f"Stored DE results (welch_ttest) in adata.uns['{key_added}']['{contrast_key}'].")

    elif method == "nb_glm":
        res = _de_nb_glm(adata, g_mask, r_mask, layer_counts=layer_counts)
        adata.uns[key_added][contrast_key] = {
            "method": method,
            "groupby": groupby,
            "group": str(group),
            "reference": str(reference),
            "layer_counts": layer_counts,
            "results": res,
        }
        info(f"Stored DE results (nb_glm) in adata.uns['{key_added}']['{contrast_key}'].")

    else:
        raise ValueError("method must be 'welch_ttest' or 'nb_glm'")


def _get_matrix(adata: ad.AnnData, layer: str) -> np.ndarray:
    if layer is None:
        X = adata.X
    else:
        if layer not in adata.layers:
            raise KeyError(f"layer='{layer}' not found in adata.layers")
        X = adata.layers[layer]
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def _de_welch_ttest(
    adata: ad.AnnData,
    g_mask: np.ndarray,
    r_mask: np.ndarray,
    *,
    layer: str = "log1p_cpm",
) -> pd.DataFrame:
    info(f"DE (welch_ttest) on layer='{layer}'")
    X = _get_matrix(adata, layer)

    Xg = X[g_mask, :]
    Xr = X[r_mask, :]

    # logFC as difference in means on log scale (interpretable as log fold-change approx)
    mean_g = Xg.mean(axis=0)
    mean_r = Xr.mean(axis=0)
    logfc = mean_g - mean_r

    # Welch t-test per gene
    t, p = stats.ttest_ind(Xg, Xr, axis=0, equal_var=False, nan_policy="omit")
    p = np.nan_to_num(p, nan=1.0)

    q = _bh_fdr(p)

    res = pd.DataFrame(
        {
            "gene": adata.var_names.to_numpy(),
            "log2FC": logfc,
            "t": t,
            "pval": p,
            "qval": q,
            "mean_group": mean_g,
            "mean_ref": mean_r,
        }
    ).sort_values("qval", ascending=True)

    return res


def _de_nb_glm(
    adata: ad.AnnData,
    g_mask: np.ndarray,
    r_mask: np.ndarray,
    *,
    layer_counts: str = "counts",
    shrink_dispersion: bool = True,
    prior_df: float = 10.0,
) -> pd.DataFrame:
    """
    DESeq2-like NB Wald test:
      - size factors (median-of-ratios)
      - gene-wise dispersion (MoM) + optional shrinkage to trend
      - NB GLM with offset(log(size_factor))
    """
    try:
        import statsmodels.api as sm
    except Exception as e:
        raise ImportError(
            "nb_glm requires statsmodels. Install with:\n"
            "  pip install statsmodels\n"
            f"Original error: {e}"
        )

    info(f"DE (nb_glm, deseq2-like) on raw counts layer='{layer_counts}'")
    Yfull = _get_matrix(adata, layer_counts)  # samples x genes

    # subset samples
    mask = g_mask | r_mask
    Y = Yfull[mask, :]
    grp = np.where(g_mask[mask], 1.0, 0.0)  # 1=group, 0=reference

    # Size factors
    sf = deseq2_size_factors(Y)
    offset = np.log(sf)

    # Normalized counts for dispersion estimation
    Yn = Y / sf[:, None]

    mu, alpha_hat = estimate_dispersion_mom(Yn)
    if shrink_dispersion:
        alpha_use = shrink_dispersion_to_trend(mu, alpha_hat, prior_df=prior_df)
    else:
        alpha_use = alpha_hat

    # Design
    X = np.column_stack([np.ones_like(grp), grp])

    n_genes = Y.shape[1]
    coef = np.zeros(n_genes, dtype=float)
    se = np.zeros(n_genes, dtype=float)
    pval = np.ones(n_genes, dtype=float)

    # avoid warning spam globally
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Negative binomial dispersion parameter alpha not set.*",
        )

        for j in range(n_genes):
            y = Y[:, j]
            if y.sum() == 0:
                continue
            a = float(alpha_use[j])

            try:
                model = sm.GLM(
                    y,
                    X,
                    family=sm.families.NegativeBinomial(alpha=a),
                    offset=offset,
                )
                fit = model.fit(maxiter=100, disp=0)

                # group coefficient (log fold-change on natural log scale)
                coef[j] = fit.params[1]
                se[j] = fit.bse[1]

                z = coef[j] / se[j] if se[j] > 0 else 0.0
                pval[j] = 2 * (1 - stats.norm.cdf(abs(z)))
            except Exception:
                continue

    qval = _bh_fdr(pval)

    res = pd.DataFrame(
        {
            "gene": adata.var_names.to_numpy(),
            "log2FC": coef / np.log(2),
            "se": se,
            "pval": pval,
            "qval": qval,
            "dispersion": alpha_use,
            "mean_norm": mu,
        }
    ).sort_values("qval", ascending=True)

    return res

