from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from .logging import info, warn


def qc_metrics(
    adata: ad.AnnData,
    *,
    mt_prefix: str | Iterable[str] = ("MT-", "mt-"),
    ribo_prefix: str | Iterable[str] = ("RPS", "RPL"),
) -> None:
    """
    Compute basic QC metrics for bulk RNA-seq data.

    Adds metrics to:
    - adata.obs (per-sample)
    - adata.var (per-gene)

    Parameters
    ----------
    adata
        AnnData object with raw counts in `.X`
    mt_prefix
        Prefix(es) identifying mitochondrial genes
    ribo_prefix
        Prefix(es) identifying ribosomal genes
    """
    info("Computing bulk RNA-seq QC metrics")

    X = adata.X

    if sp.issparse(X):
        X = X.tocsr()

    # --------------------
    # Sample-level metrics
    # --------------------
    adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    adata.obs["n_genes_detected"] = np.asarray((X > 0).sum(axis=1)).ravel()

    # --------------------
    # Gene-level metrics
    # --------------------
    adata.var["total_counts"] = np.asarray(X.sum(axis=0)).ravel()
    adata.var["n_samples_detected"] = np.asarray((X > 0).sum(axis=0)).ravel()

    # --------------------
    # Mitochondrial genes
    # --------------------
    if isinstance(mt_prefix, str):
        mt_prefix = (mt_prefix,)

    mt_mask = adata.var_names.str.startswith(tuple(mt_prefix))

    if mt_mask.any():
        mt_counts = X[:, mt_mask].sum(axis=1)
        adata.obs["pct_counts_mt"] = (
            np.asarray(mt_counts).ravel()
            / adata.obs["total_counts"]
            * 100
        )
    else:
        warn("No mitochondrial genes detected using mt_prefix")
        adata.obs["pct_counts_mt"] = 0.0

    # --------------------
    # Ribosomal genes
    # --------------------
    if isinstance(ribo_prefix, str):
        ribo_prefix = (ribo_prefix,)

    ribo_mask = adata.var_names.str.startswith(tuple(ribo_prefix))

    if ribo_mask.any():
        ribo_counts = X[:, ribo_mask].sum(axis=1)
        adata.obs["pct_counts_ribo"] = (
            np.asarray(ribo_counts).ravel()
            / adata.obs["total_counts"]
            * 100
        )
    else:
        warn("No ribosomal genes detected using ribo_prefix")
        adata.obs["pct_counts_ribo"] = 0.0

    info("QC metrics added to adata.obs and adata.var")