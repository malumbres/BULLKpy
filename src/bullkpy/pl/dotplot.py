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


# -----------------------------------------------------------------------------
# Dendrogram helpers (aligned to 0.5-centered matrix coordinates)
# -----------------------------------------------------------------------------
def _scipy_leafpos_to_index(v: np.ndarray) -> np.ndarray:
    """
    SciPy dendrogram leaf positions come in a 5,15,25,... spacing.
    Convert to 0,1,2,... indices.
    """
    v = np.asarray(v, dtype=float)
    return (v - 5.0) / 10.0


def _plot_row_dendrogram_aligned(
    ax: plt.Axes,
    Z: np.ndarray,
    n_leaves: int,
    *,
    color: str = "0.4",
    lw: float = 1.2,
    invert_y: bool = True,
    mirror_x: bool = False,
) -> None:
    """
    Draw a row dendrogram aligned to row centers at y = 0.5, 1.5, ..., n-0.5.
    """
    dd = dendrogram(Z, orientation="right", no_plot=True)

    max_x = 0.0
    if mirror_x:
        for x in dd["dcoord"]:
            max_x = max(max_x, float(np.max(x)))

    for x, y in zip(dd["dcoord"], dd["icoord"]):
        yy = _scipy_leafpos_to_index(np.asarray(y)) + 0.5  # align to centers
        xx = np.asarray(x, dtype=float)
        if mirror_x:
            xx = max_x - xx
        ax.plot(xx, yy, color=color, lw=lw)

    ax.set_ylim(0.0, float(n_leaves))
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
    """
    Draw a column dendrogram aligned to column centers at x = 0.5, 1.5, ..., n-0.5.
    """
    dd = dendrogram(Z, orientation="top", no_plot=True)
    for x, y in zip(dd["icoord"], dd["dcoord"]):
        xx = _scipy_leafpos_to_index(np.asarray(x)) + 0.5  # align to centers
        ax.plot(xx, y, color=color, lw=lw)

    ax.set_xlim(0.0, float(n_leaves))
    ax.set_xticks([])
    ax.set_yticks([])
    for spn in ax.spines.values():
        spn.set_visible(False)


# -----------------------------------------------------------------------------
# Data extraction
# -----------------------------------------------------------------------------
def _get_feature_matrix(
    adata: ad.AnnData,
    *,
    features: Sequence[str],
    layer: str | None,
    use_obs: bool,
) -> np.ndarray:
    """
    Returns (n_obs x n_features) numeric matrix for either:
      - var/features from X/layers (use_obs=False)
      - obs columns (use_obs=True)
    """
    feats = [str(f) for f in features]

    if use_obs:
        missing = [f for f in feats if f not in adata.obs.columns]
        if missing:
            raise KeyError(f"obs features not found in adata.obs (first 10): {missing[:10]}")
        M = adata.obs[feats].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        return M

    missing = [g for g in feats if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found in adata.var_names (first 10): {missing[:10]}")

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    gidx = [adata.var_names.get_loc(g) for g in feats]
    M = X[:, gidx].toarray() if sp.issparse(X) else np.asarray(X[:, gidx], dtype=float)
    return M


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def dotplot(
    adata: ad.AnnData,
    *,
    # ---- what to plot ----
    var_names: Sequence[str] | None = None,
    var_groups: dict[str, Sequence[str]] | None = None,
    # allow plotting obs (e.g., metaprograms stored in adata.obs["MP_*"])
    obs_names: Sequence[str] | None = None,
    obs_groups: dict[str, Sequence[str]] | None = None,
    use_obs: bool = False,
    # ---- grouping ----
    groupby: str | Sequence[str] = "leiden",
    # ---- expression sources ----
    layer: str | None = "log1p_cpm",
    fraction_layer: str | None = "counts",
    # Scanpy-style naming:
    expression_cutoff: float = 0.0,
    mean_only_expressed: bool = False,
    # Back-compat alias (will override expression_cutoff if provided as non-default):
    expr_threshold: float | None = None,
    # Scanpy-style standard_scale: "var" or "group" performs MIN-MAX scaling [0,1]
    standard_scale: Literal["var", "group", "zscore_var", "zscore_group"] | None = None,
    # ---- layout ----
    swap_axes: bool = False,
    row_spacing: float = 1.0,
    dendrogram_top: bool = False,
    dendrogram_rows: bool = False,
    row_dendrogram_position: Literal["right", "left", "outer_left"] = "right",
    cluster_rows: bool | None = None,
    cluster_cols: bool | None = None,
    # ---- coloring ----
    cmap: str = "Reds",
    vmin: float | None = None,
    vmax: float | None = None,
    # ---- size encoding ----
    dot_min: float | None = None,
    dot_max: float | None = None,
    # Scanpy-style name:
    size_exponent: float = 1.5,
    # Back-compat alias: if gamma is provided, it overrides size_exponent
    gamma: float | None = None,
    smallest_dot: float = 0.0,
    largest_dot: float = 200.0,
    # AUTO scaling to figsize/axes box
    scale_dots_to_fig: bool = True,
    dot_scale: float = 1.0,
    # Padding (Scanpy-like; units are "tick distance" where 1.0 is one cell)
    x_padding: float = 0.8,
    y_padding: float = 1.0,
    # ---- figure ----
    figsize: tuple[float, float] | None = None,
    invert_yaxis: bool = True,
    title: str | None = None,
    size_title: str = "Fraction of samples\nin group (%)",
    colorbar_title: str = "Mean expression\nin group",
    # NEW: allow size to represent an obs column directly (e.g. purity)
    size_obs_key: str | None = None,
    size_clip: tuple[float, float] | None = None,
    # ---- IO ----
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Dotplot like Scanpy (0.5-centered coordinates, dot scaling, minmax standard_scale),
    with BULLKpy additions:
      - can plot var genes OR numeric adata.obs columns (use_obs=True)
      - dot sizes can auto-scale to figsize/axes (scale_dots_to_fig=True)
      - optional dendrograms (top + rows)
      - optional size encoding from an obs column (size_obs_key)
    """
    set_style()

    if (dendrogram_top or dendrogram_rows) and (linkage is None or dendrogram is None):
        raise ImportError("Dendrograms require scipy (scipy.cluster.hierarchy).")

    if cluster_rows is None:
        cluster_rows = dendrogram_rows
    if cluster_cols is None:
        cluster_cols = dendrogram_top

    # cutoff alias handling
    if expr_threshold is not None:
        expression_cutoff = float(expr_threshold)

    # exponent alias handling
    if gamma is not None:
        size_exponent = float(gamma)

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

    # ---- features (var or obs) ----
    if use_obs:
        if obs_groups is not None:
            ordered: list[str] = []
            for _, feats in obs_groups.items():
                ordered.extend([str(f) for f in feats])
            obs_names = ordered
        elif obs_names is None:
            raise ValueError("With use_obs=True, provide obs_names=... or obs_groups={...}.")
        feature_names = [str(v) for v in obs_names]
    else:
        if var_groups is not None:
            ordered = []
            for _, genes in var_groups.items():
                ordered.extend([str(g) for g in genes])
            var_names = ordered
        elif var_names is None:
            raise ValueError("Provide either var_names=... or var_groups={...}.")
        feature_names = [str(v) for v in var_names]

    # ---- matrices ----
    # mean-expression matrix
    M_mean = _get_feature_matrix(adata, features=feature_names, layer=layer, use_obs=use_obs)

    # size matrix: either fraction > cutoff OR an obs key
    if size_obs_key is not None:
        if size_obs_key not in adata.obs.columns:
            raise KeyError(f"size_obs_key='{size_obs_key}' not found in adata.obs")
        sraw = pd.to_numeric(adata.obs[size_obs_key], errors="coerce").to_numpy(float)
        if size_clip is not None:
            lo, hi = map(float, size_clip)
            sraw = np.clip(sraw, lo, hi)
        # broadcast per-feature (same size for all features within each sample)
        M_size = np.repeat(sraw[:, None], repeats=len(feature_names), axis=1)
        size_mode = "mean"
    else:
        if use_obs:
            # for obs features, define "fraction expressed" as proportion above cutoff
            M_frac = M_mean.copy()
        else:
            X_frac = (
                adata.layers[fraction_layer]
                if (fraction_layer is not None and fraction_layer in adata.layers)
                else (adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X)
            )
            gidx = [adata.var_names.get_loc(g) for g in feature_names]
            M_frac = X_frac[:, gidx].toarray() if sp.issparse(X_frac) else np.asarray(X_frac[:, gidx], dtype=float)
        M_size = M_frac
        size_mode = "fraction"

    # ---- aggregate (groups x features) ----
    groups_cat = pd.Categorical(groups, categories=cats, ordered=True)
    mean_expr = np.zeros((len(cats), len(feature_names)), dtype=float)
    frac_expr = np.zeros((len(cats), len(feature_names)), dtype=float)

    for i, c in enumerate(cats):
        mask = (groups_cat == c)
        Xg = M_mean[mask, :]

        # mean expression: optionally over expressed only
        if mean_only_expressed and size_obs_key is None:
            # "expressed" defined by expression_cutoff on the SAME value used for mean
            # (if you want fraction_layer to define expression, keep mean_only_expressed=False)
            expressed = Xg > float(expression_cutoff)
            with np.errstate(invalid="ignore", divide="ignore"):
                num = np.nansum(np.where(expressed, Xg, 0.0), axis=0)
                den = np.nansum(expressed.astype(float), axis=0)
                den = np.where(den == 0, np.nan, den)
                mean_expr[i, :] = num / den
            mean_expr[i, :] = np.nan_to_num(mean_expr[i, :], nan=0.0)
        else:
            mean_expr[i, :] = np.nanmean(Xg, axis=0)

        # size encoding
        if size_mode == "mean":
            frac_expr[i, :] = np.nanmean(M_size[mask, :], axis=0)
        else:
            frac_expr[i, :] = np.nanmean(M_size[mask, :] > float(expression_cutoff), axis=0)

    # ---- standard_scale (Scanpy-like minmax) ----
    disp = mean_expr.copy()

    if standard_scale == "group":
        # per-row (group) minmax to [0,1]
        mn = np.nanmin(disp, axis=1, keepdims=True)
        mx = np.nanmax(disp, axis=1, keepdims=True)
        denom = np.where((mx - mn) == 0, 1.0, (mx - mn))
        disp = (disp - mn) / denom
        disp = np.nan_to_num(disp, nan=0.0)

    elif standard_scale == "var":
        # per-col (feature) minmax to [0,1]
        mn = np.nanmin(disp, axis=0, keepdims=True)
        mx = np.nanmax(disp, axis=0, keepdims=True)
        denom = np.where((mx - mn) == 0, 1.0, (mx - mn))
        disp = (disp - mn) / denom
        disp = np.nan_to_num(disp, nan=0.0)

    # Optional legacy z-score modes (if you still want them)
    elif standard_scale == "zscore_var":
        mu = np.nanmean(disp, axis=0, keepdims=True)
        sd = np.nanstd(disp, axis=0, ddof=0, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        disp = (disp - mu) / sd

    elif standard_scale == "zscore_group":
        mu = np.nanmean(disp, axis=1, keepdims=True)
        sd = np.nanstd(disp, axis=1, ddof=0, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        disp = (disp - mu) / sd

    elif standard_scale is None:
        pass
    else:
        raise ValueError("standard_scale must be one of: None, 'var', 'group', 'zscore_var', 'zscore_group'")

    # ---- plotted matrix (swap axes) ----
    if swap_axes:
        plot_vals = disp.T
        plot_frac = frac_expr.T
        row_labels = list(feature_names)  # features
        col_labels = list(cats)           # groups
    else:
        plot_vals = disp
        plot_frac = frac_expr
        row_labels = list(cats)
        col_labels = list(feature_names)

    plot_vals = np.asarray(plot_vals, dtype=float)
    plot_frac = np.asarray(plot_frac, dtype=float)

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
    plot_frac = plot_frac[row_order, :][:, col_order]
    row_labels = [row_labels[i] for i in row_order]
    col_labels = [col_labels[i] for i in col_order]

    # ---- autosize figure ----
    if figsize is None:
        n_x = len(col_labels)
        n_y = len(row_labels)
        w = max(5.8, 0.50 * n_x + 3.6)
        h = max(4.2, 0.42 * n_y + 2.6)
        figsize = (w, h)

    # ---- color scaling ----
    vmin_eff = vmin if vmin is not None else float(np.nanmin(plot_vals))
    vmax_eff = vmax if vmax is not None else float(np.nanmax(plot_vals))
    norm = mpl.colors.Normalize(vmin=vmin_eff, vmax=vmax_eff)
    cmap_obj = mpl.colormaps.get_cmap(cmap) if hasattr(mpl, "colormaps") else mpl.cm.get_cmap(cmap)

    has_row_dendro = bool(dendrogram_rows and (Z_row is not None) and (len(row_labels) > 2))
    has_top_dendro = bool(dendrogram_top and (Z_col is not None) and (len(col_labels) > 2))

    # ---- dot_min/dot_max (Scanpy style) ----
    frac = np.asarray(plot_frac, float)
    frac_flat = frac[np.isfinite(frac)]
    if frac_flat.size == 0:
        frac_flat = np.array([0.0], dtype=float)

    if dot_max is None:
        # ceil to nearest 0.1 like Scanpy for nicer legends
        dot_max_eff = float(np.ceil(np.nanmax(frac_flat) * 10.0) / 10.0)
    else:
        dot_max_eff = float(dot_max)

    if dot_min is None:
        dot_min_eff = 0.0
    else:
        dot_min_eff = float(dot_min)

    if not (0.0 <= dot_min_eff <= 1.0 and 0.0 <= dot_max_eff <= 1.0 and dot_min_eff <= dot_max_eff):
        raise ValueError("dot_min/dot_max must satisfy 0<=dot_min<=dot_max<=1")

    # clip + rescale to [0,1]
    f = np.clip(frac, dot_min_eff, dot_max_eff)
    if (dot_max_eff - dot_min_eff) > 0:
        u = (f - dot_min_eff) / (dot_max_eff - dot_min_eff)
    else:
        u = np.zeros_like(f, dtype=float)

    u = np.clip(u, 0.0, 1.0) ** float(size_exponent)

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

    # ---- dot sizes (AUTO scale to figure/axes) ----
    if scale_dots_to_fig:
        # compute available point space for the dot grid from the axes box
        fig.canvas.draw_idle()
        pos = ax.get_position()
        ax_w_pts = figsize[0] * 72.0 * pos.width
        ax_h_pts = figsize[1] * 72.0 * pos.height
        nx = max(1, len(col_labels))
        ny = max(1, len(row_labels))
        cell_w = ax_w_pts / nx
        cell_h = ax_h_pts / ny
        max_diam = 0.85 * min(cell_w, cell_h)
        min_diam = 0.22 * min(cell_w, cell_h)

        largest_dot_eff = (max_diam ** 2) * float(dot_scale)
        smallest_dot_eff = (min_diam ** 2) * float(dot_scale)
        sizes = smallest_dot_eff + (largest_dot_eff - smallest_dot_eff) * u
    else:
        smallest_dot_eff = float(smallest_dot)
        largest_dot_eff = float(largest_dot)
        sizes = smallest_dot_eff + (largest_dot_eff - smallest_dot_eff) * u

    # ---- main dots (Scanpy-style coordinates) ----
    n_rows, n_cols = plot_vals.shape
    yy, xx = np.indices((n_rows, n_cols))
    x = xx.ravel().astype(float) + 0.5
    y = yy.ravel().astype(float) + 0.5

    s = sizes.ravel().astype(float)
    c = cmap_obj(norm(plot_vals.ravel().astype(float)))

    ax.scatter(
        x,
        y,
        s=s,
        c=c,
        edgecolors="0.2",
        linewidths=0.35,
    )

    # ticks at centers
    x_ticks = np.arange(n_cols, dtype=float) + 0.5
    y_ticks = np.arange(n_rows, dtype=float) + 0.5
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)

    ax.set_xticklabels([str(v) for v in col_labels], rotation=90, ha="center")
    ax.set_yticklabels([str(v) for v in row_labels])

    # axis limits (Scanpy-style padding when color_on="dot")
    xpad = float(x_padding) - 0.5
    ypad = float(y_padding) - 0.5
    ax.set_xlim(-xpad, float(n_cols) + xpad)
    ax.set_ylim(float(n_rows) + ypad, -ypad)  # inverted like Scanpy

    if not invert_yaxis:
        # user wants non-inverted y: flip back
        ax.set_ylim(-ypad, float(n_rows) + ypad)

    # reduce y tick label padding so it doesn't crash into dendrogram area
    ax.tick_params(axis="y", pad=10 if has_row_dendro else 6)

    # x tick placement: keep bottom labels by default (more typical for bulk figures)
    ax.tick_params(axis="x", labeltop=False, labelbottom=True)

    ax.set_xlabel("")
    ax.set_ylabel("")
    if title is not None:
        ax.set_title(title)

    # ---- dendrograms ----
    if ax_top is not None:
        _plot_col_dendrogram_aligned(ax_top, Z_col, n_leaves=len(col_labels))

    # mirror_x=True flips the row dendrogram horizontally (so it “points” to the matrix)
    if ax_outer_left is not None:
        _plot_row_dendrogram_aligned(ax_outer_left, Z_row, n_leaves=len(row_labels), invert_y=invert_yaxis, mirror_x=True)
    if ax_inner_left is not None:
        _plot_row_dendrogram_aligned(ax_inner_left, Z_row, n_leaves=len(row_labels), invert_y=invert_yaxis, mirror_x=True)
    if ax_inner_right is not None:
        _plot_row_dendrogram_aligned(ax_inner_right, Z_row, n_leaves=len(row_labels), invert_y=invert_yaxis, mirror_x=False)

    # ---- legends ----
    ax_leg.axis("off")
    ax_leg.text(0.0, 1.00, size_title, ha="left", va="top", transform=ax_leg.transAxes)

    # legend reference points (fractions)
    ref = np.array([0.2, 0.4, 0.6, 0.8, 1.0], dtype=float)
    ref = np.clip(ref, dot_min_eff, dot_max_eff)
    if (dot_max_eff - dot_min_eff) > 0:
        ref_u = ((ref - dot_min_eff) / (dot_max_eff - dot_min_eff)) ** float(size_exponent)
    else:
        ref_u = np.zeros_like(ref)

    ref_s = smallest_dot_eff + (largest_dot_eff - smallest_dot_eff) * ref_u

    x0, y0, dx = 0.12, 0.55, 0.16
    for j, rs in enumerate(ref_s):
        ax_leg.scatter(
            [x0 + j * dx],
            [y0],
            s=float(rs),
            color="0.55",
            edgecolors="0.2",
            linewidths=0.3,
            transform=ax_leg.transAxes,
        )

    if size_obs_key is None:
        ax_leg.text(0.12, 0.25, "20  40  60  80  100", ha="left", va="center", transform=ax_leg.transAxes)
    else:
        ax_leg.text(0.12, 0.25, "low          high", ha="left", va="center", transform=ax_leg.transAxes)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    cb = fig.colorbar(sm, cax=ax_cbar, orientation="horizontal")
    cb.set_label(colorbar_title)

    # ---- margins (make room for long y-labels) ----
    if swap_axes:
        left = 0.42 if has_row_dendro else 0.30
        bottom = 0.18
    else:
        left = 0.32 if has_row_dendro else 0.20
        bottom = 0.14

    fig.subplots_adjust(left=left, right=0.98, top=0.92, bottom=bottom)

    # --- shrink legend colorbar height ---
    def _shrink_axis_box(a: plt.Axes, height_factor: float = 0.2) -> None:
        if a is None:
            return
        pos = a.get_position()
        new_h = pos.height * float(height_factor)
        a.set_position([pos.x0, pos.y0 + (pos.height - new_h), pos.width, new_h])

    _shrink_axis_box(ax_cbar, height_factor=0.2)

    # ---- optional row spacing compression (kept for compatibility) ----
    def _shrink_axis_height(a: plt.Axes, factor: float) -> None:
        if a is None:
            return
        pos = a.get_position()
        new_h = pos.height * float(factor)
        a.set_position([pos.x0, pos.y0 + (pos.height - new_h), pos.width, new_h])

    if float(row_spacing) != 1.0:
        _shrink_axis_height(ax, float(row_spacing))
        _shrink_axis_height(ax_outer_left, float(row_spacing))
        _shrink_axis_height(ax_inner_left, float(row_spacing))
        _shrink_axis_height(ax_inner_right, float(row_spacing))

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax