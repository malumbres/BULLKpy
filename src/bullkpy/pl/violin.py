from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad

from ._colors import displayed_categories, positional_palette
from ._style import set_style, _savefig
from .._genes import gene_loc as _gene_loc
from .._matrix import as_list as _as_list_impl


def _as_list(x):
    return _as_list_impl(x, empty_if_none=True)


def _get_matrix_for_layer(adata: ad.AnnData, layer: str | None):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"layer='{layer}' not found in adata.layers. Available: {list(adata.layers.keys())}")
    return adata.layers[layer]


def violin(
    adata: ad.AnnData,
    *,
    keys: list[str],
    groupby: str,
    layer: str | None = "log1p_cpm",
    figsize: tuple[float, float] = (8, 4),
    panel_size: tuple[float, float] | None = None,
    show_points: bool = True,
    point_size: float = 2.0,
    point_alpha: float = 0.35,
    palette: str | None = None,
    order: list[str] | None = None,
    rotate_xticks: float = 45,
    inner: str = "quartile",
    cut: float = 0.0,
    save: str | Path | None = None,
    show: bool = True,
):

    """
    Violin plots of sample-level variables and/or gene expression across groups.

    Notes
    -----
    - Each entry in ``keys`` is interpreted as an ``adata.obs`` column if present,
      otherwise as a gene in ``adata.var_names``.
    - Gene expression is taken from ``layer`` (or ``adata.X`` if ``layer=None``).
    """

    set_style()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    keys = [str(k) for k in keys]
    if len(keys) == 0:
        raise ValueError("keys must be a non-empty list")

    # Determine which keys are obs vs genes
    obs_keys: list[str] = []
    gene_keys: list[str] = []
    missing: list[str] = []
    for k in keys:
        if k in adata.obs.columns:
            obs_keys.append(k)
        elif k in adata.var_names:
            gene_keys.append(k)
        else:
            missing.append(k)

    if missing:
        raise KeyError(
            "Some keys were not found in adata.obs or adata.var_names: "
            f"{missing}. (obs keys available: {len(adata.obs.columns)}, genes: {adata.n_vars})"
        )

    # Build dataframe
    df = adata.obs[[groupby]].copy()

    # enforce category order
    if order is not None:
        cats = [str(x) for x in order]
        df[groupby] = pd.Categorical(df[groupby].astype(str), categories=cats, ordered=True)
    else:
        df[groupby] = df[groupby].astype("category")

    # add obs columns
    for k in obs_keys:
        df[k] = adata.obs[k].values

    # add gene expression columns
    if gene_keys:
        X = _get_matrix_for_layer(adata, layer)
        gidx = [_gene_loc(adata, g) for g in gene_keys]
        if sp.issparse(X):
            M = X[:, gidx].toarray()
        else:
            M = np.asarray(X[:, gidx], dtype=float)
        for j, g in enumerate(gene_keys):
            df[g] = M[:, j]

    # default palette choice
    # - if palette is None, let seaborn decide
    # - if many categories, husl avoids repeating as quickly as tab10/tab20
    if palette is None:
        n_cats = df[groupby].nunique(dropna=False)
        palette = "husl" if n_cats > 20 else "Set2"

    # panel sizing: either use figsize or derive from panel_size
    n = len(keys)
    if panel_size is not None:
        w = float(panel_size[0]) * n
        h = float(panel_size[1])
        fig_size = (w, h)
    else:
        fig_size = (float(figsize[0]), float(figsize[1]))

    fig, axes = plt.subplots(
        1, n,
        figsize=fig_size,
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes.ravel()

    # seaborn maps palettes to hue levels correctly only for string x with explicit order/hue_order
    shown = [str(c) for c in displayed_categories(df[groupby], order)]
    df_plot = df.assign(**{groupby: df[groupby].astype(str)})

    # plotting
    for ax, k in zip(axes, keys):
        sns.violinplot(
            data=df_plot,
            x=groupby,
            y=k,
            hue=groupby,
            hue_order=shown,
            dodge=False,
            legend=False,
            ax=ax,
            inner=inner,
            cut=cut,
            palette=positional_palette(palette, shown),
            order=shown,
        )
        if show_points:
            sns.stripplot(
                data=df_plot,
                x=groupby,
                y=k,
                ax=ax,
                color="k",
                size=float(point_size),
                alpha=float(point_alpha),
                order=shown,
            )
        ax.set_title(k)
        ax.tick_params(axis="x", rotation=float(rotate_xticks))

    # if keys < axes (shouldn't happen), hide extras
    for ax in axes[len(keys):]:
        ax.axis("off")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, axes