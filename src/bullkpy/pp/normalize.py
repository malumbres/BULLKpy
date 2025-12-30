from __future__ import annotations

from typing import Literal

import numpy as np
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def set_raw_counts(
    adata: ad.AnnData,
    *,
    layer: str = "counts",
    overwrite: bool = False,
) -> None:
    """
    Store current adata.X into adata.layers[layer] as raw counts.

    Call this once right after reading counts (before any normalization).
    """
    if layer in adata.layers and not overwrite:
        warn(f"adata.layers['{layer}'] already exists; not overwriting.")
        return
    adata.layers[layer] = adata.X.copy()
    info(f"Stored raw counts in adata.layers['{layer}'].")


def normalize_cpm(
    adata: ad.AnnData,
    *,
    layer: str | None = "counts",
    target_sum: float = 1e6,
    out_layer: str = "cpm",
    inplace_X: bool = False,
    eps: float = 1e-12,
) -> None:
    """
    CPM normalize counts per sample.

    Parameters
    ----------
    layer
        Input layer. If None, uses adata.X.
    target_sum
        Scale factor (1e6 = CPM).
    out_layer
        Where to write normalized values (adata.layers[out_layer]).
    inplace_X
        If True, also write normalized values into adata.X.
    eps
        Small constant to avoid division by zero.
    """
    X = adata.layers[layer] if layer is not None else adata.X

    if sp.issparse(X):
        X = X.tocsr()
        libsize = np.asarray(X.sum(axis=1)).ravel()
        libsize = np.maximum(libsize, eps)
        scale = target_sum / libsize
        X_norm = X.multiply(scale[:, None])
    else:
        libsize = X.sum(axis=1)
        libsize = np.maximum(libsize, eps)
        scale = (target_sum / libsize).astype(float)
        X_norm = (X.astype(float).T * scale).T

    adata.layers[out_layer] = X_norm
    info(f"Wrote CPM-normalized data to adata.layers['{out_layer}'].")

    if inplace_X:
        adata.X = X_norm
        info("Updated adata.X with CPM-normalized data.")

    # Store library size for reference (bulk QC)
    adata.obs["libsize"] = libsize


def log1p(
    adata: ad.AnnData,
    *,
    layer: str = "cpm",
    out_layer: str = "log1p_cpm",
    inplace_X: bool = False,
) -> None:
    """
    log1p-transform a layer (default: CPM) and store as a new layer.

    Works with dense or sparse matrices.
    """
    X = adata.layers[layer] if layer is not None else adata.X

    if sp.issparse(X):
        X_log = X.copy()
        X_log.data = np.log1p(X_log.data)
    else:
        X_log = np.log1p(X.astype(float))

    adata.layers[out_layer] = X_log
    info(f"Wrote log1p-transformed data to adata.layers['{out_layer}'].")

    if inplace_X:
        adata.X = X_log
        info("Updated adata.X with log1p-transformed data.")