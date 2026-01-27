from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple
import numpy as np
import pandas as pd
from anndata import AnnData


def standardize_tcga_obs(
    adata: AnnData,
    *,
    keys: Optional[Sequence[str]] = None,
    na_values: Tuple[str, ...] = ("NA", "N/A", "Unknown", "unknown", ""),
    strip: bool = True,
    lower: bool = False,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Standardize TCGA-like metadata in `adata.obs` (missing values, whitespace, casing).

    Parameters
    ----------
    adata
        AnnData with sample metadata in `.obs`.
    keys
        Which columns to standardize. If None, applies to all object/string columns.
    na_values
        Values to convert to pd.NA.
    strip
        Strip leading/trailing whitespace.
    lower
        Convert strings to lowercase (often false for TCGA labels).
    copy
        If True, returns a copy of AnnData.

    Returns
    -------
    AnnData | None
        Returns AnnData if copy=True else modifies in place and returns None.
    """
    ad = adata.copy() if copy else adata

    obs = ad.obs.copy()
    if keys is None:
        keys = [c for c in obs.columns if obs[c].dtype == object or pd.api.types.is_string_dtype(obs[c])]

    for k in keys:
        s = obs[k]
        if not (s.dtype == object or pd.api.types.is_string_dtype(s)):
            continue
        s = s.astype("string")

        if strip:
            s = s.str.strip()

        # normalize missingness
        s = s.replace(list(na_values), pd.NA)

        if lower:
            s = s.str.lower()

        obs[k] = s

    ad.obs = obs
    return ad if copy else None


def obs_map_categories(
    adata: AnnData,
    *,
    key: str,
    mapping: dict,
    out_key: Optional[str] = None,
    inplace: bool = True,
) -> Optional[pd.Series]:
    """
    Map categories in `adata.obs[key]` to harmonized labels.

    Parameters
    ----------
    adata
        AnnData.
    key
        Column in `.obs` to map.
    mapping
        Dict mapping old values -> new values.
    out_key
        If provided, writes to `.obs[out_key]` instead of overwriting `key`.
    inplace
        If True, store in `.obs`; otherwise return mapped Series.

    Returns
    -------
    pd.Series | None
    """
    if key not in adata.obs:
        raise KeyError(f"{key!r} not found in adata.obs")

    s = adata.obs[key].astype("string").map(lambda x: mapping.get(x, x))
    if not inplace:
        return s

    target = out_key or key
    adata.obs[target] = pd.Categorical(s)
    return None