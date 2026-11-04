from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import anndata as ad

from ..pl._style import set_style, _savefig
from .._matrix import as_list as _as_list
from ._embedding import plot_embedding_one as _plot_embedding_one

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None




def _pc_label(adata: ad.AnnData, pc_index_0: int, *, key: str = "pca") -> str:
    pc = pc_index_0 + 1
    try:
        vr = np.asarray(adata.uns[key]["variance_ratio"], dtype=float)
        if pc_index_0 < len(vr) and np.isfinite(vr[pc_index_0]):
            return f"PC{pc} ({vr[pc_index_0]*100:.1f}%)"
    except Exception:
        pass
    return f"PC{pc}"


def pca_scatter(
    adata: ad.AnnData,
    *,
    basis: str = "X_pca",
    components: tuple[int, int] = (1, 2),
    color: str | list[str] | None = None,
    layer: str | None = "log1p_cpm",
    point_size: float = 20.0,
    alpha: float = 0.85,
    figsize: tuple[float, float] = (6.5, 5.0),
    title: str | None = None,
    palette: str = "Set1",
    cmap: str = "viridis",
    highlight: str | list[str] | None = None,
    grey_color: str = "#D3D3D3",
    key: str = "pca",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scanpy-like PCA scatter plot.

    - color can be:
        - None
        - obs column (categorical or numeric)
        - gene name (continuous)
        - list of any of the above -> multiple panels in one row
    - highlight (categoricals): show only selected classes in color, all others grey.
    """
    set_style()

    if basis not in adata.obsm:
        raise KeyError(f"adata.obsm['{basis}'] not found. Run bk.tl.pca() first.")

    X = np.asarray(adata.obsm[basis], dtype=float)

    pcx = int(components[0]) - 1
    pcy = int(components[1]) - 1
    if pcx < 0 or pcy < 0 or pcx >= X.shape[1] or pcy >= X.shape[1]:
        raise ValueError(f"components={components} out of range for {basis} with shape {X.shape}")

    x = X[:, pcx]
    y = X[:, pcy]

    colors = _as_list(color) if color is not None else [None]
    n = len(colors)

    # scale width with number of panels (scanpy-like behavior)
    fig_w = figsize[0] * n if n > 1 else figsize[0]
    fig_h = figsize[1]
    fig, axes = plt.subplots(1, n, figsize=(fig_w, fig_h), constrained_layout=True)
    if n == 1:
        axes = [axes]

    xlabel = _pc_label(adata, pcx, key=key)
    ylabel = _pc_label(adata, pcy, key=key)

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
            title=(title if (title is not None and n == 1) else ("PCA" if c is None else str(c))),
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes

def pca_variance_ratio(
    adata: ad.AnnData,
    *,
    key: str = "pca",
    n_comps: int | None = None,
    cumulative: bool = True,
    figsize: tuple[float, float] = (6.5, 4.5),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scree plot of PCA variance ratio (+ optional cumulative curve).
    Requires adata.uns[key]['variance_ratio'].
    """
    set_style()

    if key not in adata.uns or "variance_ratio" not in adata.uns[key]:
        raise KeyError(f"Missing adata.uns['{key}']['variance_ratio']. Run bk.tl.pca() first.")

    vr = np.asarray(adata.uns[key]["variance_ratio"], dtype=float)
    if n_comps is not None:
        vr = vr[: int(n_comps)]

    xs = np.arange(1, len(vr) + 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xs, vr, marker="o")
    ax.set_xlabel("PC")
    ax.set_ylabel("Explained variance ratio")
    ax.set_title("PCA variance ratio")

    if cumulative:
        ax2 = ax.twinx()
        ax2.plot(xs, np.cumsum(vr), marker="o")
        ax2.set_ylabel("Cumulative variance")

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax