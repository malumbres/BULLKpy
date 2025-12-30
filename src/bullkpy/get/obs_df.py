from __future__ import annotations

import pandas as pd
import anndata as ad

from .vector import vector


def obs_df(
    adata: ad.AnnData,
    keys: list[str],
    *,
    layer: str | None = None,
) -> pd.DataFrame:
    """
    Scanpy-like helper: return a DataFrame with requested keys.

    Each key can be:
    - an obs column
    - a gene name (from adata.var_names), taken from `layer` if provided

    Returns
    -------
    pd.DataFrame indexed by adata.obs_names.
    """
    data = {}
    for k in keys:
        data[k] = vector(adata, k, layer=layer)
    return pd.DataFrame(data, index=adata.obs_names)