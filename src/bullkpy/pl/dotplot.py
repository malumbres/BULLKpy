from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
except Exception:  # pragma: no cover
    linkage = None
    dendrogram = None

import anndata as ad

from ._style import set_style, _savefig


def _scipy_leafpos_to_index(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return (v - 5.0) / 10.0


def _plot_row_dendrogram_aligned(
    ax: plt.Axes,
    Z: np.ndarray,
    n_leaves: int,
    *,
    row_spacing: float = 1.0,
    color: str = "0.4",
    lw: float = 1.2,
    invert_y: bool = True,
    mirror_x: bool = False,
) -> None:
    """
    Draw a row dendrogram aligned to row centers:
      rows at y = [0, row_spacing, 2*row_spacing, ...]
    """
    dd = dendrogram(Z, orientation="right", no_plot=True)

    max_x = 0.0
    if mirror_x:
        for x in dd["dcoord"]:
            max_x = max(max_x, float(np.max(x)))

    for x, y in zip(dd["dcoord"], dd["icoord"]):
        yy = _scipy_leafpos_to_index(np.asarray(y))
        xx = np.asarray(x, dtype=float)
        if mirror_x:
            xx = max_x - xx
        ax.plot(xx, yy, color=color, lw=lw)

    
    ax.set_ylim(-0.5, n_leaves - 0.5)
    if invert_y:
        ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])
    for spn in ax.spines.values():
        spn.set_visible(False)


def _plot_col_dendrogram_aligned(
    ax: plt.Axes,
    Z: np.ndarray,
    n_leaves: int,
    *,
    color: str = "0.4",
    lw: float = 1.2,
) -> None:
    dd = dendrogram(Z, orientation="top", no_plot=True)
    for x, y in zip(dd["icoord"], dd["dcoord"]):
        xx = _scipy_leafpos_to_index(np.asarray(x))
        ax.plot(xx, y, color=color, lw=lw)

    ax.set_xlim(-0.5, n_leaves - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spn in ax.spines.values():
        spn.set_visible(False)


def dotplot(
    adata: ad.AnnData,
    *,
    var_names: Sequence[str] | None = None,
    var_groups: dict[str, Sequence[str]] | None = None,
    groupby: str | Sequence[str] = "leiden",
    layer: str | None = "log1p_cpm",
    fraction_layer: str | None = "counts",
    expr_threshold: float = 0.0,
    standard_scale: str | None = None,
    swap_axes: bool = False,
    row_spacing: float = 0.7,
    dendrogram_top: bool = False,
    dendrogram_rows: bool = False,
    row_dendrogram_position: Literal["right", "left", "outer_left"] = "right",
    cluster_rows: bool | None = None,
    cluster_cols: bool | None = None,
    cmap: str = "Reds",
    vmin: float | None = None,
    vmax: float | None = None,
    dot_min: float = 0.0,
    dot_max: float = 1.0,
    gamma: float = 0.5,
    smallest_dot: float = 12.0,
    largest_dot: float = 260.0,
    figsize: tuple[float, float] | None = None,
    invert_yaxis: bool = True,
    title: str | None = None,
    size_title: str = "Fraction of samples\nin group (%)",
    colorbar_title: str = "Mean expression\nin group",
    save: str | Path | None = None,
    show: bool = True,
):
    set_style()

    if (dendrogram_top or dendrogram_rows) and (linkage is None or dendrogram is None):
        raise ImportError("Dendrograms require scipy (scipy.cluster.hierarchy).")

    if cluster_rows is None:
        cluster_rows = dendrogram_rows
    if cluster_cols is None:
        cluster_cols = dendrogram_top

    # ---- groupby ----
    if isinstance(groupby, (list, tuple)):
        for g in groupby:
            if g not in adata.obs.columns:
                raise KeyError(f"groupby='{g}' not found in adata.obs")
        grp_df = adata.obs[list(groupby)].copy()
        grp_key = grp_df.astype(str).agg(" | ".join, axis=1)
        groups = pd.Categorical(grp_key)
    else:
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not found in adata.obs")
        groups = adata.obs[groupby].astype("category")

    cats = list(pd.Categorical(groups).categories)

    # ---- genes ----
    if var_groups is not None:
        ordered_genes: list[str] = []
        for _, genes in var_groups.items():
            ordered_genes.extend([str(g) for g in genes])
        var_names = ordered_genes
    elif var_names is None:
        raise ValueError("Provide either var_names=... or var_groups={...}.")

    var_names = [str(v) for v in var_names]
    missing = [g for g in var_names if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found in adata.var_names (first 10): {missing[:10]}")

    # ---- matrices ----
    X_mean = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    X_frac = (
        adata.layers[fraction_layer]
        if (fraction_layer is not None and fraction_layer in adata.layers)
        else X_mean
    )

    gidx = [adata.var_names.get_loc(g) for g in var_names]
    M_mean = X_mean[:, gidx].toarray() if sp.issparse(X_mean) else np.asarray(X_mean[:, gidx], dtype=float)
    M_frac = X_frac[:, gidx].toarray() if sp.issparse(X_frac) else np.asarray(X_frac[:, gidx], dtype=float)

    # ---- aggregate (groups x genes) ----
    groups_cat = pd.Categorical(groups, categories=cats, ordered=True)
    mean_expr = np.zeros((len(cats), len(var_names)), dtype=float)
    frac_expr = np.zeros((len(cats), len(var_names)), dtype=float)

    for i, c in enumerate(cats):
        mask = (groups_cat == c)
        mean_expr[i, :] = M_mean[mask, :].mean(axis=0)
        frac_expr[i, :] = (M_frac[mask, :] > expr_threshold).mean(axis=0)

    disp = mean_expr.copy()
    if standard_scale == "var":
        mu = disp.mean(axis=0, keepdims=True)
        sd = disp.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        disp = (disp - mu) / sd
    elif standard_scale == "group":
        mu = disp.mean(axis=1, keepdims=True)
        sd = disp.std(axis=1, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        disp = (disp - mu) / sd

    # ---- sizes ----
    f = np.clip(frac_expr, dot_min, dot_max)
    u = np.clip((f - dot_min) / (dot_max - dot_min), 0, 1) ** float(gamma)
    sizes = smallest_dot + (largest_dot - smallest_dot) * u

    # ---- plotted matrix ----
    if swap_axes:
        plot_vals = disp.T
        plot_sizes = sizes.T
        row_labels = list(var_names)  # genes
        col_labels = list(cats)       # groups
    else:
        plot_vals = disp
        plot_sizes = sizes
        row_labels = list(cats)
        col_labels = list(var_names)

    plot_vals = np.asarray(plot_vals, dtype=float)
    plot_sizes = np.asarray(plot_sizes, dtype=float)

    # ---- clustering (orders) ----
    row_order = np.arange(plot_vals.shape[0])
    col_order = np.arange(plot_vals.shape[1])
    Z_row = None
    Z_col = None

    if cluster_rows and plot_vals.shape[0] > 2:
        Z_row = linkage(plot_vals, method="average", metric="euclidean")
        row_order = np.array(dendrogram(Z_row, no_plot=True)["leaves"], dtype=int)

    if cluster_cols and plot_vals.shape[1] > 2:
        Z_col = linkage(plot_vals.T, method="average", metric="euclidean")
        col_order = np.array(dendrogram(Z_col, no_plot=True)["leaves"], dtype=int)

    plot_vals = plot_vals[row_order, :][:, col_order]
    plot_sizes = plot_sizes[row_order, :][:, col_order]
    row_labels = [row_labels[i] for i in row_order]
    col_labels = [col_labels[i] for i in col_order]

    # ---- autosize ----
    if figsize is None:
        n_x = len(col_labels)
        n_y = len(row_labels)
        w = max(5.2, 0.55 * n_x + 3.4)
        h = max(3.8, 0.55 * n_y + 2.0)
        figsize = (w, h)

    # ---- color scaling ----
    vmin_eff = vmin if vmin is not None else float(np.nanmin(plot_vals))
    vmax_eff = vmax if vmax is not None else float(np.nanmax(plot_vals))
    norm = mpl.colors.Normalize(vmin=vmin_eff, vmax=vmax_eff)
    cmap_obj = mpl.cm.get_cmap(cmap)

    has_row_dendro = bool(dendrogram_rows and (Z_row is not None) and (len(row_labels) > 2))
    has_top_dendro = bool(dendrogram_top and (Z_col is not None) and (len(col_labels) > 2))

    # ---- layout ----
    outer_left = 0.20 if (has_row_dendro and row_dendrogram_position == "outer_left") else 0.001
    inner_left = 0.18 if (has_row_dendro and row_dendrogram_position == "left") else 0.001
    inner_right = 0.18 if (has_row_dendro and row_dendrogram_position == "right") else 0.001
    legends_w = 0.58

    fig = plt.figure(figsize=figsize, constrained_layout=False)
    gs = GridSpec(
        nrows=2,
        ncols=5,
        figure=fig,
        height_ratios=[0.22, 1.0] if has_top_dendro else [0.001, 1.0],
        width_ratios=[outer_left, inner_left, 1.0, inner_right, legends_w],
        hspace=0.05,
        wspace=0.14,
    )

    ax_top = fig.add_subplot(gs[0, 2]) if has_top_dendro else None
    ax = fig.add_subplot(gs[1, 2])

    ax_outer_left = fig.add_subplot(gs[1, 0]) if (has_row_dendro and row_dendrogram_position == "outer_left") else None
    ax_inner_left = fig.add_subplot(gs[1, 1]) if (has_row_dendro and row_dendrogram_position == "left") else None
    ax_inner_right = fig.add_subplot(gs[1, 3]) if (has_row_dendro and row_dendrogram_position == "right") else None

    gs_leg = gs[:, 4].subgridspec(nrows=2, ncols=1, height_ratios=[0.58, 0.42], hspace=0.25)
    ax_leg = fig.add_subplot(gs_leg[0, 0])
    ax_cbar = fig.add_subplot(gs_leg[1, 0])

    # ---- main dots ----
    xs = np.arange(len(col_labels))
    ys = np.arange(len(row_labels), dtype=float)

    for yi in range(plot_vals.shape[0]):
        ax.scatter(
            xs,
            np.full(xs.shape, ys[yi], dtype=float),
            s=plot_sizes[yi, :],
            c=cmap_obj(norm(plot_vals[yi, :])),
            edgecolors="0.2",
            linewidths=0.35,
        )

    ax.set_xlim(-0.5, len(col_labels) - 0.5)
    ax.set_ylim(-0.5, len(row_labels) - 0.5)
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.set_xticklabels([str(x) for x in col_labels], rotation=90)
    ax.set_yticklabels([str(y) for y in row_labels])

    ax.tick_params(axis="y", pad=40 if swap_axes else 18)       # MM:adapt pad=XX to gene names

    if has_top_dendro:
        ax.tick_params(axis="x", labeltop=False, labelbottom=True)
    else:
        ax.tick_params(axis="x", labeltop=True, labelbottom=False)

    if invert_yaxis:
        ax.invert_yaxis()

    ax.set_xlabel("")
    ax.set_ylabel("")
    if title is not None:
        ax.set_title(title)


    # ---- dendrograms ----
    if ax_top is not None:
        _plot_col_dendrogram_aligned(ax_top, Z_col, n_leaves=len(col_labels))

    if ax_outer_left is not None:
        _plot_row_dendrogram_aligned(
            ax_outer_left, Z_row, n_leaves=len(row_labels),
            row_spacing=row_spacing, invert_y=invert_yaxis, mirror_x=True
        )
    if ax_inner_left is not None:
        _plot_row_dendrogram_aligned(
            ax_inner_left, Z_row, n_leaves=len(row_labels),
            row_spacing=row_spacing, invert_y=invert_yaxis, mirror_x=True
        )
    if ax_inner_right is not None:
        _plot_row_dendrogram_aligned(
            ax_inner_right, Z_row, n_leaves=len(row_labels),
            row_spacing=row_spacing, invert_y=invert_yaxis, mirror_x=False
        )

    # ---- legends ----
    ax_leg.axis("off")
    ax_leg.text(0.0, 1.00, size_title, ha="left", va="top", transform=ax_leg.transAxes)

    ref = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
    ref_u = np.clip((ref - dot_min) / (dot_max - dot_min), 0, 1) ** float(gamma)
    ref_s = smallest_dot + (largest_dot - smallest_dot) * ref_u

    x0, y0, dx = 0.12, 0.55, 0.16
    for j, rs in enumerate(ref_s):
        ax_leg.scatter(
            [x0 + j * dx],
            [y0],
            s=rs,
            color="0.55",
            edgecolors="0.2",
            linewidths=0.3,
            transform=ax_leg.transAxes,
        )
    ax_leg.text(0.12, 0.25, "20  40  60  80  100", ha="left", va="center", transform=ax_leg.transAxes)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    cb = fig.colorbar(sm, cax=ax_cbar, orientation="horizontal")
    cb.set_label(colorbar_title)


    # ---- margins ----
    if swap_axes:
        # genes on y-axis -> need more space on the left
        left = 0.42 if has_row_dendro else 0.30
        bottom = 0.18
    else:
        left = 0.32 if has_row_dendro else 0.20
        bottom = 0.14

    fig.subplots_adjust(left=left, right=0.98, top=0.92, bottom=bottom)

    # --- shrink the dot panel height to match row_spacing (so no big empty area) ---
    def _shrink_axis_height(a: plt.Axes, factor: float) -> None:
        if a is None:
            return
        pos = a.get_position()
        new_h = pos.height * float(factor)
        # anchor at the top so it shrinks downward (Scanpy-like)
        a.set_position([pos.x0, pos.y0 + (pos.height - new_h), pos.width, new_h])

    def _shrink_axis_box(a: plt.Axes, height_factor: float = 0.55) -> None:
        if a is None:
            return
        pos = a.get_position()
        new_h = pos.height * float(height_factor)
        a.set_position([pos.x0, pos.y0 + (pos.height - new_h), pos.width, new_h])

    # Make colorbar shorter (Scanpy-like)
    _shrink_axis_box(ax_cbar, height_factor=0.2)      # check height factor

    if row_spacing != 1.0:
        _shrink_axis_height(ax, row_spacing)

        # shrink whichever row dendrogram axis is being used, so it stays aligned
        if ax_outer_left is not None:
            _shrink_axis_height(ax_outer_left, row_spacing)
        if ax_inner_left is not None:
            _shrink_axis_height(ax_inner_left, row_spacing)
        if ax_inner_right is not None:
            _shrink_axis_height(ax_inner_right, row_spacing)


    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax