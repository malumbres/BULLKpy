from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def highly_variable_genes(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    n_top_genes: int = 2000,
    n_bins: int = 20,
    min_mean: float = 0.0,
    max_mean: float = np.inf,
    min_disp: float = -np.inf,
    key_added: str = "highly_variable",
) -> None:
    """
    

    Select highly variable genes (bulk-friendly version of HVGs).

    This computes mean and variance across samples on the chosen layer
    and identifies genes with high dispersion relative to genes of similar mean
    (mean-binned z-score of log-dispersion, Scanpy/Seurat spirit).

    Stores results in `adata.var`:
      - mean, variance, dispersion, dispersion_norm
      - boolean flag `adata.var[key_added]` (default: 'highly_variable')
    """
    X = adata.layers[layer] if layer is not None else adata.X
    info(f"Computing HVGs from layer={layer!r} with n_top_genes={n_top_genes}")

    if sp.issparse(X):
        X = X.tocsr()
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
        var = mean_sq - mean**2
        var[var < 0] = 0.0  # numerical safety
    else:
        X = np.asarray(X, dtype=float)
        mean = X.mean(axis=0)
        var = X.var(axis=0, ddof=0)

    # Dispersion (variance / mean) in log space (stable)
    with np.errstate(divide="ignore", invalid="ignore"):
        disp = var / mean
    disp[~np.isfinite(disp)] = 0.0

    # Filter by mean/disp thresholds (optional)
    valid = (
        (mean >= min_mean)
        & (mean <= max_mean)
        & (disp >= min_disp)
    )

    if valid.sum() == 0:
        raise ValueError("No genes pass mean/disp filters; relax thresholds.")

    # Bin genes by mean (on log scale to spread bins better)
    mean_for_bins = np.log1p(mean)
    bins = pd.qcut(mean_for_bins[valid], q=min(n_bins, int(valid.sum())), duplicates="drop")

    disp_log = np.log1p(disp)  # log(1+disp)
    disp_norm = np.full_like(disp_log, fill_value=np.nan, dtype=float)

    # Z-score within mean bins
    valid_idx = np.where(valid)[0]
    for b in bins.categories:
        mask = np.zeros_like(valid, dtype=bool)
        # bins is only for valid genes; map back to global indices

        in_bin = np.asarray(bins == b)
        mask[valid_idx[in_bin]] = True

        vals = disp_log[mask]
        if vals.size < 2:
            continue
        mu = np.mean(vals)
        sd = np.std(vals, ddof=0)
        if sd == 0:
            sd = 1.0
        disp_norm[mask] = (vals - mu) / sd

    # Rank genes by normalized dispersion (higher = more variable)
    # Fill NaNs with -inf so they never get selected
    score = np.nan_to_num(disp_norm, nan=-np.inf)
    n_top = int(min(n_top_genes, score.size))
    top_idx = np.argsort(score)[::-1][:n_top]

    hvg = np.zeros(score.size, dtype=bool)
    hvg[top_idx] = True

    # Store in adata.var
    adata.var["means"] = mean
    adata.var["variances"] = var
    adata.var["dispersions"] = disp
    adata.var["dispersions_norm"] = disp_norm
    adata.var[key_added] = hvg

    info(f"Marked {hvg.sum()} highly variable genes in adata.var['{key_added}']")