from __future__ import annotations

from typing import Iterable, Optional, Sequence
import re
import numpy as np
import pandas as pd
from anndata import AnnData


def _mad(x: np.ndarray, axis: int = 0) -> np.ndarray:
    med = np.nanmedian(x, axis=axis, keepdims=True)
    return np.nanmedian(np.abs(x - med), axis=axis)


def sample_dispersion(
    adata: AnnData,
    *,
    layer: str = "log1p_cpm",
    genes: Optional[Sequence[str]] = None,
    groupby: Optional[str] = "Project_ID",
    standardize_within_group: bool = True,
    method: str = "mad",  # "mad" | "iqr" | "std"
    out_key: str = "dispersion_sample",
    min_genes: int = 5000,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Compute per-sample transcriptomic dispersion (robust 'heterogeneity' proxy).

    If `standardize_within_group=True`, dispersion is computed on z-scored genes
    within each `groupby` group (recommended for cross-tumor comparability).

    Writes `adata.obs[out_key]`.
    """
    ad = adata.copy() if copy else adata

    if layer not in ad.layers:
        raise KeyError(f"layer {layer!r} not found in adata.layers")

    X = ad.layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)

    var_names = np.asarray(ad.var_names)
    if genes is not None:
        genes = [g for g in genes if g in ad.var_names]
        if len(genes) < min(50, min_genes):
            raise ValueError(f"Too few genes found in var_names ({len(genes)}).")
        idx = np.isin(var_names, genes)
        X = X[:, idx]
    else:
        if X.shape[1] < min_genes:
            raise ValueError(f"Too few genes ({X.shape[1]}) for robust dispersion.")

    obs = ad.obs
    if groupby is not None and groupby not in obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    disp = np.full(X.shape[0], np.nan, dtype=float)

    if groupby is None or not standardize_within_group:
        Z = X
        if method == "mad":
            disp = _mad(Z, axis=1)
        elif method == "iqr":
            disp = np.nanpercentile(Z, 75, axis=1) - np.nanpercentile(Z, 25, axis=1)
        elif method == "std":
            disp = np.nanstd(Z, axis=1)
        else:
            raise ValueError("method must be one of {'mad','iqr','std'}")
    else:
        for grp, idx_rows in obs.groupby(groupby).groups.items():
            idx_rows = np.asarray(list(idx_rows))
            Z = X[idx_rows, :]

            # z-score per gene within group
            mu = np.nanmean(Z, axis=0, keepdims=True)
            sd = np.nanstd(Z, axis=0, keepdims=True)
            sd[sd == 0] = np.nan
            Z = (Z - mu) / sd

            if method == "mad":
                disp[idx_rows] = _mad(Z, axis=1)
            elif method == "iqr":
                disp[idx_rows] = np.nanpercentile(Z, 75, axis=1) - np.nanpercentile(Z, 25, axis=1)
            elif method == "std":
                disp[idx_rows] = np.nanstd(Z, axis=1)
            else:
                raise ValueError("method must be one of {'mad','iqr','std'}")

    ad.obs[out_key] = disp
    return ad if copy else None


def signature_dispersion_by_group(
    adata: AnnData,
    *,
    score_keys: Sequence[str],
    groupby: str = "Project_ID",
    method: str = "mad",
    out_key: str = "signature_dispersion",
) -> pd.DataFrame:
    """
    Compute dispersion of signature scores across samples within each group.

    Returns a tidy DataFrame: [group, signature, dispersion, n].
    Also stores in `adata.uns[out_key]`.
    """
    obs = adata.obs
    if groupby not in obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    for k in score_keys:
        if k not in obs:
            raise KeyError(f"score key {k!r} not found in adata.obs")

    rows = []
    for grp, df in obs.groupby(groupby):
        for sig in score_keys:
            x = df[sig].to_numpy(dtype=float)
            if method == "mad":
                d = float(np.nanmedian(np.abs(x - np.nanmedian(x))))
            elif method == "iqr":
                d = float(np.nanpercentile(x, 75) - np.nanpercentile(x, 25))
            elif method == "std":
                d = float(np.nanstd(x))
            else:
                raise ValueError("method must be one of {'mad','iqr','std'}")
            rows.append({"group": grp, "signature": sig, "dispersion": d, "n": int(np.isfinite(x).sum())})

    out = pd.DataFrame(rows)
    adata.uns[out_key] = out
    return out

