from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn

def _is_integerish(X, max_check: int = 50000, tol: float = 1e-8) -> bool:
    """Heuristic: does X look like raw counts (non-negative integers)?"""
    if sp.issparse(X):
        data = X.data
        if data.size == 0:
            return True  # all zeros
        if data.size > max_check:
            data = data[:max_check]
        if np.min(data) < 0:
            return False
        frac = np.abs(data - np.round(data))
        return bool(np.max(frac) < tol)
    else:
        A = np.asarray(X)
        if A.size == 0:
            return True
        flat = A.ravel()
        if flat.size > max_check:
            flat = flat[:max_check]
        if np.min(flat) < 0:
            return False
        frac = np.abs(flat - np.round(flat))
        return bool(np.max(frac) < tol)


def qc_metrics(
    adata: ad.AnnData,
    *,
    layer: str | None = "counts",
    mt_prefix: str = "MT-",
    mt_var_key: str | None = None,
    ribo_prefixes: tuple[str, ...] = ("RPS", "RPL"),
    ribo_var_key: str | None = None,
    detection_threshold: float = 0.0,
    compute_total_counts: bool = True,
    compute_n_genes: bool = True,
    compute_pct_mt: bool = True,
    compute_pct_ribo: bool = True,
) -> None:
    """
    Compute bulk QC metrics in ``adata.obs``.

    Adds, where applicable: ``total_counts``, ``n_genes_detected``,
    ``pct_counts_mt`` and ``pct_counts_ribo``.

    Smart behavior:
      - If `layer` looks like raw counts: compute all metrics.
      - If `layer` looks non-integer (log/normalized): compute only
        ``n_genes_detected`` using `detection_threshold`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    layer
        Layer to read counts from; falls back to ``adata.X`` when absent.
    mt_prefix
        Gene-name prefix identifying mitochondrial genes.
    mt_var_key
        Optional boolean column in ``adata.var`` flagging mitochondrial genes.
        Takes precedence over `mt_prefix`.
    ribo_prefixes
        Gene-name prefixes identifying ribosomal protein genes.
    ribo_var_key
        Optional boolean column in ``adata.var`` flagging ribosomal genes.
        Takes precedence over `ribo_prefixes`.
    detection_threshold
        Expression above which a gene counts as detected (non-count layers).
    compute_total_counts, compute_n_genes, compute_pct_mt, compute_pct_ribo
        Toggles for the individual metrics.
    """
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    is_counts = _is_integerish(X)

    if not is_counts:
        # Logged/normalized: only detection QC makes sense
        if compute_total_counts or compute_pct_mt:
            warn(
                f"qc_metrics: layer='{layer}' does not look like raw integer counts. "
                "Skipping total_counts and pct_counts_mt; computing only n_genes_detected "
                f"using detection_threshold>{detection_threshold}."
            )
        compute_total_counts = False
        compute_pct_mt = False
        compute_pct_ribo = False

    # --- total_counts ---
    if compute_total_counts:
        if sp.issparse(X):
            adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
        else:
            adata.obs["total_counts"] = np.sum(np.asarray(X), axis=1)

    # --- n_genes_detected ---
    if compute_n_genes:
        thr = float(detection_threshold) if not is_counts else 0.0
        if sp.issparse(X):
            adata.obs["n_genes_detected"] = np.asarray((X > thr).sum(axis=1)).ravel()
        else:
            adata.obs["n_genes_detected"] = (np.asarray(X) > thr).sum(axis=1)

    gene_names = pd.Series(adata.var_names.astype(str)).str.upper()

    if sp.issparse(X):
        tot = np.asarray(X.sum(axis=1)).ravel()
    else:
        tot = np.sum(np.asarray(X), axis=1)
    tot_safe = np.where(tot == 0, np.nan, tot)

    def _pct_for(mask: np.ndarray) -> np.ndarray:
        Xs = X[:, mask]
        sub = np.asarray(Xs.sum(axis=1)).ravel() if sp.issparse(X) else np.sum(np.asarray(Xs), axis=1)
        return 100.0 * (sub / tot_safe)

    added = []
    if compute_total_counts:
        added.append("total_counts")
    if compute_n_genes:
        added.append("n_genes_detected")

    # --- pct_counts_mt ---
    if compute_pct_mt:
        if mt_var_key is not None and mt_var_key in adata.var.columns:
            mt_mask = adata.var[mt_var_key].astype(bool).to_numpy()
        else:
            mt_mask = gene_names.str.startswith(mt_prefix.upper()).to_numpy()

        if mt_mask.sum() == 0:
            warn(
                "qc_metrics: no mitochondrial genes detected with given "
                "mt_prefix/mt_var_key; pct_counts_mt not computed."
            )
        else:
            adata.obs["pct_counts_mt"] = _pct_for(mt_mask)
            added.append("pct_counts_mt")

    # --- pct_counts_ribo ---
    if compute_pct_ribo:
        if ribo_var_key is not None and ribo_var_key in adata.var.columns:
            ribo_mask = adata.var[ribo_var_key].astype(bool).to_numpy()
        else:
            prefixes = tuple(p.upper() for p in ribo_prefixes)
            ribo_mask = gene_names.str.startswith(prefixes).to_numpy()

        if ribo_mask.sum() == 0:
            warn(
                "qc_metrics: no ribosomal genes detected with given "
                "ribo_prefixes/ribo_var_key; pct_counts_ribo not computed."
            )
        else:
            adata.obs["pct_counts_ribo"] = _pct_for(ribo_mask)
            added.append("pct_counts_ribo")

    info(f"QC metrics added to adata.obs: {', '.join(added) if added else '(none)'}.")