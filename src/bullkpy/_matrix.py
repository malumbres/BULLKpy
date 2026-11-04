"""Expression-matrix and argument helpers shared across modules."""
from __future__ import annotations

from typing import Any, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = ["get_matrix", "as_list"]


def get_matrix(
    adata,
    layer: str | None = None,
    *,
    dense: bool = True,
    use: Literal["samples", "genes"] = "samples",
) -> Any:
    """Expression matrix of ``adata.X`` or a layer, samples x genes.

    A layer that is requested but missing raises ``KeyError``. Falling back to
    ``adata.X`` instead would silently analyse a different quantity, such as raw
    counts where log-CPM was meant.

    Parameters
    ----------
    dense
        Return a float ``ndarray``. With ``False`` the matrix is returned as stored,
        sparse or dense.
    use
        ``"genes"`` returns the transpose (genes x samples).
    """
    if layer is None:
        X = adata.X
    elif layer in adata.layers:
        X = adata.layers[layer]
    else:
        raise KeyError(f"layer {layer!r} not found in adata.layers; available: {list(adata.layers)}")
    if dense:
        X = X.toarray() if sp.issparse(X) else X
        X = np.asarray(X, dtype=float)
    return X.T if use == "genes" else X


def as_list(x, *, empty_if_none: bool = False) -> list | None:
    """Wrap a scalar in a list; list-like values become lists; ``None`` is kept.

    With ``empty_if_none=True``, ``None`` becomes ``[]``.
    """
    if x is None:
        return [] if empty_if_none else None
    if isinstance(x, (list, tuple, set, np.ndarray, pd.Index)):
        return list(x)
    return [x]
