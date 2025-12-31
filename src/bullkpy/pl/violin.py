from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad

from ._style import set_style, _savefig


def _get_matrix_for_layer(adata: ad.AnnData, layer: str | None):
    if layer is None:
        return adata.X
    if layer not in adata.layers:
        raise KeyError(f"layer='{layer}' not found in adata.layers. Available: {list(adata.layers.keys())}")
    return adata.layers[layer]


def _as_list(x) -> list:
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return list(x)
    return [x]


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

    This function behaves Scanpy-like: each entry in `keys` is interpreted as
    (i) an `adata.obs` column if present, otherwise (ii) a gene in `adata.var_names`.

    Parameters
    ----------
    adata
        AnnData object with samples in `.obs` and genes in `.var_names`.
    keys
        List of variables to plot. Each key can be either:
        - name of a column in `adata.obs` (QC/clinical/signature scores), or
        - a gene name found in `adata.var_names` (expression will be taken from `layer`).
    groupby
        Categorical column in `adata.obs` used to define groups on the x-axis.
    layer
        Layer to use for gene expression keys (default: `"log1p_cpm"`).
        If `None`, uses `adata.X`.
    figsize
        Base figure size `(width, height)` in inches.
        If multiple keys are provided, the final width is scaled by the number of panels.
    panel_size
        Alternative to `figsize`: size per panel `(width, height)` in inches.
        If provided, overrides `figsize`.
    show_points
        Whether to overlay individual samples as points (strip plot).
    point_size
        Point size for the overlaid samples.
    point_alpha
        Transparency for the overlaid points.
    palette
        Categorical palette name (matplotlib/seaborn). If `None`, uses global defaults.
    order
        Explicit order of categories for `groupby`. If `None`, uses category order in `adata.obs[groupby]`.
    rotate_xticks
        Rotation angle (degrees) for x-axis tick labels.
    inner
        Passed to `seaborn.violinplot(inner=...)` (e.g. `"quartile"`, `"box"`, `None`).
    cut
        Passed to `seaborn.violinplot(cut=...)`.
    save
        If provided, path to save the figure.
    show
        If True, displays the plot (matplotlib `plt.show()`).

    Returns
    -------
    fig
        Matplotlib Figure.
    axes
        Array/list of Axes for each panel.

    Notes
    -----
    - If `groupby` is numeric with many unique values, consider converting it to a categorical.
    - For gene keys, missing genes are ignored/raised depending on implementation; see error message.

    See Also
    --------
    bullkpy.pl.gene_association : association tests for genes vs categories
    bullkpy.tl.score_genes      : compute signature scores for plotting in violin

    Examples
    --------
    QC variables:

    >>> bk.pl.violin(adata, keys=["total_counts", "pct_counts_mt"], groupby="Project_ID")

    Gene expression:

    >>> bk.pl.violin(adata, keys=["DLL3", "SOX10"], groupby="Subtype_PAM50", layer="log1p_cpm")

    Control category order and tick rotation:

    >>> bk.pl.violin(
    ...     adata,
    ...     keys=["CDC20"],
    ...     groupby="Project_ID",
    ...     order=["LUAD", "LUSC", "BRCA"],
    ...     rotate_xticks=90,
    ... )
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
        gidx = [adata.var_names.get_loc(g) for g in gene_keys]
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

    # plotting
    for ax, k in zip(axes, keys):
        sns.violinplot(
            data=df,
            x=groupby,
            y=k,
            ax=ax,
            inner=inner,
            cut=cut,
            palette=palette,
            order=order,
        )
        if show_points:
            sns.stripplot(
                data=df,
                x=groupby,
                y=k,
                ax=ax,
                color="k",
                size=float(point_size),
                alpha=float(point_alpha),
                order=order,
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