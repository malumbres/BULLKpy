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

from ..logging import warn
from ._style import set_style, _savefig

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


def _categorical_palette(categories: Sequence[str], palette: str = "Set1") -> dict[str, str]:
    cats = [str(c) for c in categories]
    n = len(cats)

    if sns is not None:
        try:
            cols = sns.color_palette(palette, n_colors=n)
            return {cats[i]: to_hex(cols[i]) for i in range(n)}
        except Exception:
            pass

    cmap = mpl.cm.get_cmap(palette) if palette in plt.colormaps() else mpl.cm.get_cmap("tab20")
    return {cats[i]: to_hex(cmap(i / max(n - 1, 1))) for i in range(n)}


def _is_categorical_series(s: pd.Series) -> bool:
    return pd.api.types.is_categorical_dtype(s.dtype) or (s.dtype == object)


def _is_numeric_series(s: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(s.dtype)


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
    xlabel: str,
    ylabel: str,
):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

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

        # optional numeric highlight behavior (simple, useful default):
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