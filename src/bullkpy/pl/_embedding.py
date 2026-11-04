"""Scatter of a 2-D embedding (PCA, UMAP, ...) coloured by an obs column or a gene."""
from __future__ import annotations

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .._compat import is_categorical_like
from .._matrix import as_list, get_gene_vector
from ..logging import warn
from ._colors import categorical_palette

__all__ = ["plot_embedding_one"]


def _scatter(ax, x, y, point_size, alpha, **kw):
    return ax.scatter(x, y, s=point_size, alpha=alpha, edgecolors="none", **kw)


def plot_embedding_one(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    adata: ad.AnnData,
    color: str | None,
    layer: str | None,
    point_size: float,
    alpha: float,
    palette: str,
    cmap: str,
    highlight: str | list[str] | None,
    grey_color: str,
    title: str | None,
    xlabel: str | None = None,
    ylabel: str | None = None,
) -> None:
    """Draw one panel.

    ``color`` may be None, a gene (continuous, from ``layer``) or an obs column
    (categorical or numeric). With ``highlight``, only the listed classes of a
    categorical column are coloured and the rest are grey; for a numeric column,
    zero values are greyed out.
    """
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    if color is None:
        _scatter(ax, x, y, point_size, alpha)
        return

    if color in adata.var_names:
        vals = get_gene_vector(adata, color, layer)
        sc = _scatter(ax, x, y, point_size, alpha, c=vals, cmap=cmap)
        plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05).set_label(str(color))
        return

    if color not in adata.obs.columns:
        warn(f"color='{color}' not found in adata.obs or adata.var_names; plotting without color.")
        _scatter(ax, x, y, point_size, alpha)
        return

    s = adata.obs[color]
    hl = as_list(highlight)
    hl = None if hl is None else [str(v) for v in hl]

    if pd.api.types.is_numeric_dtype(s.dtype) and not is_categorical_like(s):
        vals = s.to_numpy(dtype=float)
        if hl is not None:
            _scatter(ax, x, y, point_size, alpha, color=grey_color)
            keep = np.isfinite(vals) & (vals != 0)
            sc = _scatter(ax, x[keep], y[keep], point_size, alpha, c=vals[keep], cmap=cmap)
        else:
            sc = _scatter(ax, x, y, point_size, alpha, c=vals, cmap=cmap)
        plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05)
        return

    labels = s.astype(str)
    names = pd.Categorical(labels).categories.tolist()
    shown = hl if hl is not None else names
    colors = categorical_palette(shown, palette=palette)
    if hl is not None:
        _scatter(ax, x, y, point_size, alpha, color=grey_color)
    for name in shown:
        mask = labels.values == name
        _scatter(ax, x[mask], y[mask], point_size, alpha, color=colors.get(str(name), "k"), label=str(name))
    legend_title = f"{color} (highlight)" if hl is not None else str(color)
    ax.legend(title=legend_title, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
