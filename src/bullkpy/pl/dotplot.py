from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
except Exception:  # pragma: no cover
    linkage = None
    dendrogram = None

import anndata as ad

from ..logging import warn
from ._style import set_style, _savefig
from .._genes import gene_loc as _gene_loc


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
    gidx = [_gene_loc(adata, g) for g in feats]
    M = X[:, gidx].toarray() if sp.issparse(X) else np.asarray(X[:, gidx], dtype=float)
    return M


def _font_pt(size) -> float:
    """Font size in points; rcParams may hold names such as "medium"."""
    from matplotlib.font_manager import FontProperties

    return float(FontProperties(size=size).get_size_in_points())


def _text_extent_in(texts, fontsize: float) -> float:
    """Widest of ``texts`` in inches when drawn at ``fontsize`` (no pyplot state touched)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    texts = [str(t) for t in texts if str(t)]
    if not texts:
        return 0.0
    fig = Figure(figsize=(1, 1))
    renderer = FigureCanvasAgg(fig).get_renderer()
    widths = [fig.text(0, 0, t, fontsize=fontsize).get_window_extent(renderer).width / fig.dpi for t in texts]
    return float(max(widths))


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
    cell_size: float = 0.32,
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

    Layout is computed in inches from the number of rows and columns and the label lengths:
    ``cell_size`` is the width and height of one dot cell (default 0.32 in); the panel,
    dendrograms and legends are sized from it and the dots scale with the cell. Pass
    ``figsize`` to fix the total size instead; the panel then takes the space left after the
    labels and legends.
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
            gidx = [_gene_loc(adata, g) for g in feature_names]
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

    # ---- color scaling ----
    vmin_eff = vmin if vmin is not None else float(np.nanmin(plot_vals))
    vmax_eff = vmax if vmax is not None else float(np.nanmax(plot_vals))
    norm = mpl.colors.Normalize(vmin=vmin_eff, vmax=vmax_eff)
    cmap_obj = mpl.colormaps.get_cmap(cmap)

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

    # ---- layout in physical units (inches) ----
    n_rows, n_cols = plot_vals.shape
    xpad = float(x_padding) - 0.5
    ypad = float(y_padding) - 0.5
    tick_fs = _font_pt(mpl.rcParams["ytick.labelsize"])
    label_fs = _font_pt(mpl.rcParams["axes.labelsize"])
    line_in = 1.3 * label_fs / 72.0  # height of one line of text

    pad = 0.08
    dendro_w = 0.55
    dendro_h = 0.50
    ylab_w = _text_extent_in(row_labels, tick_fs) + 0.12
    xlab_h = _text_extent_in(col_labels, tick_fs) + 0.12
    size_title_lines = (str(size_title).count("\n") + 1) if size_title else 0
    cbar_title_lines = (str(colorbar_title).count("\n") + 1) if colorbar_title else 0

    left_dendro = has_row_dendro and row_dendrogram_position == "left"
    outer_dendro = has_row_dendro and row_dendrogram_position == "outer_left"
    right_dendro = has_row_dendro and row_dendrogram_position == "right"

    left_margin = ylab_w + (dendro_w + pad if outer_dendro else 0.0) + (dendro_w + pad if left_dendro else 0.0)
    right_of_panel = (dendro_w + pad) if right_dendro else 0.0
    bottom_margin = xlab_h
    title_h = (line_in + pad) if title else 0.0
    top_extra = ((dendro_h + pad) if has_top_dendro else 0.0) + title_h
    top_pad = 0.10
    gap_legend = 0.30
    row_h = float(row_spacing)

    def _dot_diameters(cell_w_in: float, cell_h_in: float) -> tuple[float, float]:
        if not scale_dots_to_fig:
            return float(np.sqrt(smallest_dot)), float(np.sqrt(largest_dot))
        base = min(cell_w_in, cell_h_in) * 72.0
        return 0.22 * base * float(np.sqrt(dot_scale)), 0.85 * base * float(np.sqrt(dot_scale))

    def _legend_geometry(max_diam_pt: float) -> dict:
        d_in = max_diam_pt / 72.0
        spacing = max(1.35 * d_in, 0.30)
        dots_w = 4.0 * spacing + d_in
        text_w = max(
            _text_extent_in(str(size_title or "").split("\n"), label_fs),
            _text_extent_in(str(colorbar_title or "").split("\n"), label_fs),
        )
        w = max(dots_w + 0.1, text_w, 1.30)
        size_h = size_title_lines * line_in + (0.08 if size_title_lines else 0.0) + d_in + 0.06 + line_in
        cbar_h = cbar_title_lines * line_in + (0.08 if cbar_title_lines else 0.0) + 0.13 + 1.1 * line_in
        return dict(d_in=d_in, spacing=spacing, w=w, size_h=size_h, cbar_h=cbar_h, stack_h=size_h + 0.25 + cbar_h)

    if figsize is None:
        cell_in = float(cell_size)
        cell_w_in, cell_h_in = cell_in, cell_in * row_h
        panel_w = (n_cols + 2 * xpad) * cell_w_in
        panel_h = (n_rows + 2 * ypad) * cell_h_in
        leg = _legend_geometry(_dot_diameters(cell_w_in, cell_h_in)[1])
        fig_w = left_margin + panel_w + right_of_panel + gap_legend + leg["w"] + 0.05
        fig_h = max(bottom_margin + panel_h + top_extra, bottom_margin + leg["stack_h"]) + top_pad
    else:
        fig_w, fig_h = float(figsize[0]), float(figsize[1])
        leg = _legend_geometry(_dot_diameters(float(cell_size), float(cell_size) * row_h)[1])
        # figsize is a request: grow it when labels and legends cannot fit around a readable panel
        min_cell = 0.20
        need_w = left_margin + right_of_panel + gap_legend + leg["w"] + 0.05 + (n_cols + 2 * xpad) * min_cell
        need_h = (
            max(bottom_margin + top_extra + (n_rows + 2 * ypad) * min_cell * row_h, bottom_margin + leg["stack_h"])
            + top_pad
        )
        if fig_w < need_w - 1e-6 or fig_h < need_h - 1e-6:
            new_w, new_h = max(fig_w, need_w), max(fig_h, need_h)
            warn(
                f"dotplot: figsize=({fig_w:g}, {fig_h:g}) is too small for the labels and legends; "
                f"using ({new_w:.1f}, {new_h:.1f}). Omit figsize to size the plot automatically."
            )
            fig_w, fig_h = new_w, new_h
        for _ in range(2):  # legend width depends on dot size, which depends on the panel
            panel_w = max(0.5, fig_w - left_margin - right_of_panel - gap_legend - leg["w"] - 0.05)
            panel_h = max(0.5, fig_h - bottom_margin - top_extra - top_pad)
            cell_w_in = panel_w / (n_cols + 2 * xpad)
            cell_h_in = panel_h / (n_rows + 2 * ypad)
            leg = _legend_geometry(_dot_diameters(cell_w_in, cell_h_in)[1])

    smallest_d, largest_d = _dot_diameters(cell_w_in, cell_h_in)
    smallest_dot_eff, largest_dot_eff = smallest_d**2, largest_d**2
    sizes = smallest_dot_eff + (largest_dot_eff - smallest_dot_eff) * u

    fig = plt.figure(figsize=(fig_w, fig_h), constrained_layout=False)

    def _rect(x: float, y: float, w: float, h: float) -> list[float]:
        return [x / fig_w, y / fig_h, w / fig_w, h / fig_h]

    x0, y0 = left_margin, bottom_margin
    ax = fig.add_axes(_rect(x0, y0, panel_w, panel_h))
    ax_top = fig.add_axes(_rect(x0, y0 + panel_h + pad, panel_w, dendro_h)) if has_top_dendro else None
    ax_outer_left = fig.add_axes(_rect(0.0, y0, dendro_w, panel_h)) if outer_dendro else None
    ax_inner_left = fig.add_axes(_rect(x0 - pad - dendro_w, y0, dendro_w, panel_h)) if left_dendro else None
    ax_inner_right = fig.add_axes(_rect(x0 + panel_w + pad, y0, dendro_w, panel_h)) if right_dendro else None

    # ---- main dots (Scanpy-style coordinates) ----
    yy, xx = np.indices((n_rows, n_cols))
    x = xx.ravel().astype(float) + 0.5
    y = yy.ravel().astype(float) + 0.5
    s = sizes.ravel().astype(float)
    c = cmap_obj(norm(plot_vals.ravel().astype(float)))
    ax.scatter(x, y, s=s, c=c, edgecolors="0.2", linewidths=0.35)
    ax.grid(False)

    ax.set_xticks(np.arange(n_cols, dtype=float) + 0.5)
    ax.set_yticks(np.arange(n_rows, dtype=float) + 0.5)
    ax.set_xticklabels([str(v) for v in col_labels], rotation=90, ha="center")
    ax.set_yticklabels([str(v) for v in row_labels])

    ax.set_xlim(-xpad, float(n_cols) + xpad)
    if invert_yaxis:
        ax.set_ylim(float(n_rows) + ypad, -ypad)  # first group on top, like Scanpy
    else:
        ax.set_ylim(-ypad, float(n_rows) + ypad)

    # an inner-left dendrogram sits between the labels and the panel
    ax.tick_params(axis="y", pad=(dendro_w + pad) * 72.0 + 4.0 if left_dendro else 4.0)
    ax.tick_params(axis="x", labeltop=False, labelbottom=True, pad=3.0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    if title is not None:
        top_y = y0 + panel_h + ((dendro_h + pad) if has_top_dendro else 0.0) + 0.04
        fig.text((x0 + panel_w / 2) / fig_w, top_y / fig_h, title, ha="center", va="bottom", fontsize=label_fs)

    # ---- dendrograms (limits copied from the panel so leaves sit on the dots) ----
    if ax_top is not None:
        _plot_col_dendrogram_aligned(ax_top, Z_col, n_leaves=len(col_labels))
        ax_top.set_xlim(ax.get_xlim())
    for dax, mirror in ((ax_outer_left, True), (ax_inner_left, True), (ax_inner_right, False)):
        if dax is not None:
            _plot_row_dendrogram_aligned(dax, Z_row, n_leaves=len(row_labels), invert_y=invert_yaxis, mirror_x=mirror)
            dax.set_ylim(ax.get_ylim())

    # ---- legends: one column to the right, size legend above colour bar ----
    x_leg = x0 + panel_w + right_of_panel + gap_legend
    top_y = fig_h - top_pad
    d_in, spacing = leg["d_in"], leg["spacing"]

    if size_title:
        fig.text(x_leg / fig_w, top_y / fig_h, size_title, ha="left", va="top", fontsize=label_fs, linespacing=1.15)
    dots_y = top_y - size_title_lines * line_in - (0.08 if size_title_lines else 0.0) - d_in / 2.0

    ref = np.clip(np.round(np.linspace(0.2, 1.0, 5) * dot_max_eff * 20.0) / 20.0, dot_min_eff, dot_max_eff)
    if (dot_max_eff - dot_min_eff) > 0:
        ref_u = ((ref - dot_min_eff) / (dot_max_eff - dot_min_eff)) ** float(size_exponent)
    else:
        ref_u = np.zeros_like(ref)
    ref_s = smallest_dot_eff + (largest_dot_eff - smallest_dot_eff) * ref_u

    ax_leg = fig.add_axes([0.0, 0.0, 1.0, 1.0], facecolor="none")  # full-figure overlay in inches
    ax_leg.set_xlim(0.0, fig_w)
    ax_leg.set_ylim(0.0, fig_h)
    ax_leg.axis("off")
    for j, (val, rs) in enumerate(zip(ref, ref_s)):
        cx = x_leg + d_in / 2.0 + j * spacing
        ax_leg.scatter([cx], [dots_y], s=float(rs), color="0.55", edgecolors="0.2", linewidths=0.3)
        if size_obs_key is None:
            label = f"{int(round(val * 100))}"
        else:
            label = "low" if j == 0 else ("high" if j == len(ref) - 1 else "")
        ax_leg.text(cx, dots_y - d_in / 2.0 - 0.05, label, ha="center", va="top", fontsize=tick_fs)

    cb_top = top_y - leg["size_h"] - 0.25
    if colorbar_title:
        fig.text(x_leg / fig_w, cb_top / fig_h, colorbar_title, ha="left", va="top", fontsize=label_fs, linespacing=1.15)
    bar_y = cb_top - cbar_title_lines * line_in - (0.08 if cbar_title_lines else 0.0) - 0.13
    ax_cbar = fig.add_axes(_rect(x_leg, bar_y, max(1.1, min(leg["w"], 1.5)), 0.13))
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj), cax=ax_cbar, orientation="horizontal")
    cb.ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
    cb.ax.tick_params(labelsize=tick_fs, length=2.5, pad=2)
    cb.outline.set_linewidth(0.5)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax
