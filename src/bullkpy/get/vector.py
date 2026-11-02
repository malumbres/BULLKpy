from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import anndata as ad
from .._genes import gene_loc as _gene_loc


def vector(
    adata: ad.AnnData,
    key: str,
    *,
    layer: str | None = None,
    allow_obs: bool = True,
    allow_gene: bool = True,
) -> np.ndarray:
    """
    Return a 1D numpy vector from AnnData.

    - If key is in adata.obs and allow_obs=True -> returns obs column
    - If key is in adata.var_names and allow_gene=True -> returns gene expression
      from `layer` (if provided and present) otherwise adata.X.

    Notes
    -----
    - For categorical obs, this returns the raw pandas array; cast outside as needed.
    - For gene expression, result is float array shape (n_obs,).
    """
    if allow_obs and key in adata.obs.columns:
        v = adata.obs[key].to_numpy()
        return np.asarray(v)

    if allow_gene and key in adata.var_names:
        if layer is not None and layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found. Available: {list(adata.layers)}")
        X = adata.layers[layer] if layer is not None else adata.X
        j = int(_gene_loc(adata, key))
        v = X[:, j]
        if sp.issparse(v):
            v = v.toarray()
        return np.asarray(v, dtype=float).ravel()

    raise KeyError(f"'{key}' not found in adata.obs or adata.var_names")