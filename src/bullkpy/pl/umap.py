from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import anndata as ad

from ._style import set_style, _savefig
from .._matrix import as_list as _as_list
from ._embedding import plot_embedding_one as _plot_embedding_one

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None


def umap(
    adata: ad.AnnData,
    *,
    basis: str = "X_umap",
    color: str | list[str] | None = None,
    layer: str | None = "log1p_cpm",
    point_size: float = 25.0,
    alpha: float = 0.85,
    figsize: tuple[float, float] = (6.5, 5.5),
    title: str | None = None,
    palette: str = "Set1",
    cmap: str = "viridis",
    highlight: str | list[str] | None = None,
    grey_color: str = "#D3D3D3",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scanpy-like UMAP plotting.

    - color can be:
        - None
        - obs column (categorical or numeric)
        - gene name (continuous, from `layer`)
        - list of any of the above -> multiple panels in one row
    - highlight (categoricals): show only selected classes in color, all others grey.
    """
    set_style()

    if basis not in adata.obsm:
        raise KeyError(f"adata.obsm['{basis}'] not found. Run bk.tl.umap(adata) first.")

    X = np.asarray(adata.obsm[basis], dtype=float)
    x = X[:, 0]
    y = X[:, 1]

    colors = _as_list(color) if color is not None else [None]
    n = len(colors)

    fig_w = figsize[0] * n if n > 1 else figsize[0]
    fig_h = figsize[1]
    fig, axes = plt.subplots(1, n, figsize=(fig_w, fig_h), constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, c in zip(axes, colors):
        _plot_embedding_one(
            ax,
            x,
            y,
            adata=adata,
            color=c,
            layer=layer,
            point_size=float(point_size),
            alpha=float(alpha),
            palette=palette,
            cmap=cmap,
            highlight=highlight,
            grey_color=grey_color,
            title=(title if (title is not None and n == 1) else ("UMAP" if c is None else str(c))),
            xlabel="UMAP1",
            ylabel="UMAP2",
        )

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes