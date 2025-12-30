from __future__ import annotations

from typing import Literal

import numpy as np
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def _get_X(adata: ad.AnnData, layer: str | None) -> np.ndarray | sp.spmatrix:
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"layer='{layer}' not found in adata.layers")
    return adata.layers[layer]


def _sum_axis0(X):
    return np.asarray(X.sum(axis=0)).ravel() if sp.issparse(X) else np.asarray(X).sum(axis=0)


def _sum_axis1(X):
    return np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else np.asarray(X).sum(axis=1)


def _count_nonzero_axis0(X, threshold: float):
    if sp.issparse(X):
        return np.asarray((X > threshold).sum(axis=0)).ravel()
    return (np.asarray(X) > threshold).sum(axis=0)


def _count_nonzero_axis1(X, threshold: float):
    if sp.issparse(X):
        return np.asarray((X > threshold).sum(axis=1)).ravel()
    return (np.asarray(X) > threshold).sum(axis=1)


def _infer_gene_mask(
    adata: ad.AnnData,
    kind: Literal["mt", "ribo"],
    *,
    prefer_var_key: str | None = None,
) -> np.ndarray | None:
    """
    Return boolean mask over genes for mt/ribo, or None if cannot infer.
    """
    if prefer_var_key is not None and prefer_var_key in adata.var.columns:
        m = adata.var[prefer_var_key]
        if m.dtype == bool:
            return m.to_numpy()
        # allow 0/1
        try:
            return m.to_numpy(dtype=bool)
        except Exception:
            pass

    names = adata.var_names.astype(str)

    if kind == "mt":
        # common conventions
        # human: "MT-" prefix; mouse: "mt-"
        return names.str.upper().str.startswith("MT-").to_numpy()

    if kind == "ribo":
        # common ribosomal protein prefixes
        up = names.str.upper()
        return (up.str.startswith("RPS") | up.str.startswith("RPL")).to_numpy()

    return None


def filter_genes(
    adata: ad.AnnData,
    *,
    layer: str | None = "counts",
    min_samples: int = 3,
    min_expr: float = 0.0,
    inplace: bool = True,
) -> ad.AnnData:
    """
    Filter genes by detection across samples.

    A gene is "detected" in a sample if expression > min_expr in `layer`.
    Keeps genes detected in at least `min_samples` samples.

    Works for counts (min_expr=0) or log layers (use e.g. min_expr=0.1).
    """
    X = _get_X(adata, layer)
    det = _count_nonzero_axis0(X, threshold=float(min_expr))
    keep = det >= int(min_samples)

    info(
        f"Filtering genes: layer={layer!r}, min_expr>{min_expr}, min_samples={min_samples} "
        f"-> keep {int(keep.sum())}/{adata.n_vars}"
    )

    # store metrics
    adata.var["n_samples_detected"] = det

    adata.uns.setdefault("pp", {})
    adata.uns["pp"]["filter_genes"] = {
        "layer": layer,
        "min_expr": float(min_expr),
        "min_samples": int(min_samples),
        "n_vars_before": int(adata.n_vars),
        "n_vars_after": int(keep.sum()),
    }

    if inplace:
        adata._inplace_subset_var(keep)
        return adata
    return adata[:, keep].copy()


def filter_samples(
    adata: ad.AnnData,
    *,
    # thresholds
    min_genes: int | None = None,
    min_counts: float | None = None,
    max_pct_mt: float | None = None,
    max_pct_ribo: float | None = None,
    # where to compute from if obs columns missing
    layer: str | None = "counts",
    expr_threshold_for_genes: float = 0.0,
    mt_var_key: str | None = None,
    ribo_var_key: str | None = None,
    inplace: bool = True,
) -> ad.AnnData:
    """
    Filter samples using QC metrics in `adata.obs` if present; otherwise compute from matrix.

    No hard requirement that any of these columns exist:
      - total_counts
      - n_genes_detected
      - pct_counts_mt
      - pct_counts_ribo

    If a requested threshold is provided but cannot be computed, it will be skipped with a warning.
    """
    n_before = adata.n_obs
    keep = np.ones(adata.n_obs, dtype=bool)

    # ---------- total_counts ----------
    if min_counts is not None:
        if "total_counts" in adata.obs.columns:
            tc = adata.obs["total_counts"].to_numpy(dtype=float)
        else:
            try:
                X = _get_X(adata, layer)
                tc = _sum_axis1(X).astype(float)
                adata.obs["total_counts"] = tc
                info(f"Computed adata.obs['total_counts'] from layer={layer!r}")
            except Exception as e:
                warn(f"min_counts set but could not compute total_counts (layer={layer!r}): {e}. Skipping min_counts.")
                tc = None

        if tc is not None:
            keep &= tc >= float(min_counts)

    # ---------- n_genes_detected ----------
    if min_genes is not None:
        if "n_genes_detected" in adata.obs.columns:
            ng = adata.obs["n_genes_detected"].to_numpy(dtype=float)
        else:
            try:
                X = _get_X(adata, layer)
                ng = _count_nonzero_axis1(X, threshold=float(expr_threshold_for_genes)).astype(float)
                adata.obs["n_genes_detected"] = ng
                info(
                    f"Computed adata.obs['n_genes_detected'] from layer={layer!r} "
                    f"(expr > {expr_threshold_for_genes})"
                )
            except Exception as e:
                warn(f"min_genes set but could not compute n_genes_detected (layer={layer!r}): {e}. Skipping min_genes.")
                ng = None

        if ng is not None:
            keep &= ng >= float(min_genes)

    # ---------- pct_counts_mt ----------
    if max_pct_mt is not None:
        if "pct_counts_mt" in adata.obs.columns:
            pm = adata.obs["pct_counts_mt"].to_numpy(dtype=float)
        else:
            try:
                X = _get_X(adata, layer)
                tc = adata.obs["total_counts"].to_numpy(dtype=float) if "total_counts" in adata.obs.columns else _sum_axis1(X)
                mt_mask = _infer_gene_mask(adata, "mt", prefer_var_key=mt_var_key)
                if mt_mask is None or mt_mask.sum() == 0:
                    raise ValueError("Could not infer mitochondrial genes (mt_mask empty).")
                mt_counts = _sum_axis1(X[:, mt_mask])
                pm = (mt_counts / np.maximum(tc, 1e-12)) * 100.0
                adata.obs["pct_counts_mt"] = pm
                info(f"Computed adata.obs['pct_counts_mt'] from layer={layer!r}")
            except Exception as e:
                warn(f"max_pct_mt set but could not compute pct_counts_mt (layer={layer!r}): {e}. Skipping max_pct_mt.")
                pm = None

        if pm is not None:
            keep &= pm <= float(max_pct_mt)

    # ---------- pct_counts_ribo ----------
    if max_pct_ribo is not None:
        if "pct_counts_ribo" in adata.obs.columns:
            pr = adata.obs["pct_counts_ribo"].to_numpy(dtype=float)
        else:
            try:
                X = _get_X(adata, layer)
                tc = adata.obs["total_counts"].to_numpy(dtype=float) if "total_counts" in adata.obs.columns else _sum_axis1(X)
                ribo_mask = _infer_gene_mask(adata, "ribo", prefer_var_key=ribo_var_key)
                if ribo_mask is None or ribo_mask.sum() == 0:
                    raise ValueError("Could not infer ribosomal genes (ribo_mask empty).")
                ribo_counts = _sum_axis1(X[:, ribo_mask])
                pr = (ribo_counts / np.maximum(tc, 1e-12)) * 100.0
                adata.obs["pct_counts_ribo"] = pr
                info(f"Computed adata.obs['pct_counts_ribo'] from layer={layer!r}")
            except Exception as e:
                warn(f"max_pct_ribo set but could not compute pct_counts_ribo (layer={layer!r}): {e}. Skipping max_pct_ribo.")
                pr = None

        if pr is not None:
            keep &= pr <= float(max_pct_ribo)

    n_after = int(keep.sum())
    info(f"Filtering samples -> keep {n_after}/{n_before}")

    adata.uns.setdefault("pp", {})
    adata.uns["pp"]["filter_samples"] = {
        "min_genes": None if min_genes is None else int(min_genes),
        "min_counts": None if min_counts is None else float(min_counts),
        "max_pct_mt": None if max_pct_mt is None else float(max_pct_mt),
        "max_pct_ribo": None if max_pct_ribo is None else float(max_pct_ribo),
        "layer": layer,
        "expr_threshold_for_genes": float(expr_threshold_for_genes),
        "n_obs_before": int(n_before),
        "n_obs_after": int(n_after),
    }

    if inplace:
        adata._inplace_subset_obs(keep)
        return adata
    return adata[keep, :].copy()