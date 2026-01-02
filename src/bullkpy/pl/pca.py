from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import to_hex
import anndata as ad

from ..logging import info, warn
from .._settings import settings
from ..pl._style import set_style, _savefig

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None




def _as_list(x):
    if x is None:
        return None
    if isinstance(x, (list, tuple, set)):
        return list(x)
    return [x]


def _pc_label(adata: ad.AnnData, pc_index_0: int, *, key: str = "pca") -> str:
    pc = pc_index_0 + 1
    try:
        vr = np.asarray(adata.uns[key]["variance_ratio"], dtype=float)
        if pc_index_0 < len(vr) and np.isfinite(vr[pc_index_0]):
            return f"PC{pc} ({vr[pc_index_0]*100:.1f}%)"
    except Exception:
        pass
    return f"PC{pc}"


def _categorical_palette(categories: Sequence[str], palette: str = "Set1") -> dict[str, str]:
    cats = [str(c) for c in categories]
    n = len(cats)

    # seaborn palettes are usually nicer and can go > 20 colors (e.g., husl)
    if sns is not None:
        try:
            cols = sns.color_palette(palette, n_colors=n)
            return {cats[i]: to_hex(cols[i]) for i in range(n)}
        except Exception:
            pass

    # fallback to matplotlib colormap
    cmap = mpl.cm.get_cmap(palette) if palette in plt.colormaps() else mpl.cm.get_cmap("tab20")
    return {cats[i]: to_hex(cmap(i / max(n - 1, 1))) for i in range(n)}


def _get_gene_vector(adata: ad.AnnData, gene: str, *, layer: str | None) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names")

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    j = adata.var_names.get_loc(gene)
    v = X[:, j]

    if sp.issparse(v):
        v = v.toarray().ravel()
    else:
        v = np.asarray(v).ravel()
    return v.astype(float)


def _is_categorical_series(s: pd.Series) -> bool:
    return pd.api.types.is_categorical_dtype(s.dtype) or (s.dtype == object)


def _is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s.dtype)


def _plot_embedding_one(
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
):
    # no color
    if color is None:
        ax.scatter(x, y, s=point_size, alpha=alpha, edgecolors="none")
        if title:
            ax.set_title(title)
        return

    # gene expression -> continuous
    if color in adata.var_names:
        vals = _get_gene_vector(adata, color, layer=layer)
        sc = ax.scatter(x, y, c=vals, cmap=cmap, s=point_size, alpha=alpha, edgecolors="none")
        cb = plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05)
        cb.set_label(str(color))
        if title:
            ax.set_title(title)
        return

    # obs column
    if color not in adata.obs.columns:
        warn(f"color='{color}' not found in adata.obs or adata.var_names; plotting without color.")
        ax.scatter(x, y, s=point_size, alpha=alpha, edgecolors="none")
        if title:
            ax.set_title(title)
        return

    s = adata.obs[color]
    hl = _as_list(highlight)
    if hl is not None:
        hl = [str(v) for v in hl]

    # categorical
    if _is_categorical_series(s):
        cats = s.astype(str)
        names = pd.Categorical(cats).categories.tolist()

        # highlight mode: plot grey background, then highlighted classes
        if hl is not None:
            ax.scatter(x, y, s=point_size, alpha=alpha, color=grey_color, edgecolors="none")
            cmap_map = _categorical_palette(hl, palette=palette)
            for name in hl:
                mask = (cats.values == name)
                ax.scatter(
                    x[mask], y[mask],
                    s=point_size,
                    alpha=alpha,
                    color=cmap_map.get(name, "k"),
                    edgecolors="none",
                    label=name,
                )
            ax.legend(
                title=f"{color} (highlight)",
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                frameon=False,
            )
            if title:
                ax.set_title(title)
            return

        # normal categorical mode: each category colored
        cmap_map = _categorical_palette(names, palette=palette)
        for name in names:
            mask = (cats.values == name)
            ax.scatter(
                x[mask], y[mask],
                s=point_size,
                alpha=alpha,
                edgecolors="none",
                color=cmap_map[str(name)],
                label=str(name),
            )
        ax.legend(title=str(color), bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        if title:
            ax.set_title(title)
        return

    # numeric (continuous)
    if _is_numeric_series(s):
        vals = s.to_numpy(dtype=float)

        # highlight for numeric is optional; simplest useful behavior:
        # if highlight is provided -> only show nonzero (or finite) values in color, rest grey
        if hl is not None:
            ax.scatter(x, y, s=point_size, alpha=alpha, color=grey_color, edgecolors="none")
            mask = np.isfinite(vals) & (vals != 0)
            sc = ax.scatter(x[mask], y[mask], c=vals[mask], cmap=cmap, s=point_size, alpha=alpha, edgecolors="none")
            plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05)
            if title:
                ax.set_title(title)
            return

        sc = ax.scatter(x, y, c=vals, cmap=cmap, s=point_size, alpha=alpha, edgecolors="none")
        plt.colorbar(sc, ax=ax, pad=0.01, fraction=0.05)
        if title:
            ax.set_title(title)
        return

    # fallback: treat as categorical strings
    cats = s.astype(str)
    names = pd.Categorical(cats).categories.tolist()

    if hl is not None:
        ax.scatter(x, y, s=point_size, alpha=alpha, color=grey_color, edgecolors="none")
        cmap_map = _categorical_palette(hl, palette=palette)
        for name in hl:
            mask = (cats.values == name)
            ax.scatter(x[mask], y[mask], s=point_size, alpha=alpha, color=cmap_map.get(name, "k"), edgecolors="none", label=name)
        ax.legend(title=f"{color} (highlight)", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
        if title:
            ax.set_title(title)
        return

    cmap_map = _categorical_palette(names, palette=palette)
    for name in names:
        mask = (cats.values == name)
        ax.scatter(x[mask], y[mask], s=point_size, alpha=alpha, edgecolors="none", color=cmap_map[str(name)], label=str(name))
    ax.legend(title=str(color), bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    if title:
        ax.set_title(title)


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