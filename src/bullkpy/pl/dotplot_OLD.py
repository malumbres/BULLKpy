from __future__ import annotations

from pathlib import Path
from typing import Sequence

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


def dotplot(
    adata: ad.AnnData,
    *,
    var_names: Sequence[str],
    groupby: str,
    layer: str | None = "log1p_cpm",
    expr_threshold: float = 0.0,
    standard_scale: str | None = None,  # "var" (per gene across groups) | "group" | None
    dendrogram_rows: bool = False,
    dendrogram_cols: bool = False,
    var_group_positions: Sequence[tuple[int, int, str]] | None = None,
    # e.g. [(0,2,"B"), (2,5,"Myeloid"), ...] indices in var_names
    cmap: str = "Reds",
    vmin: float | None = None,
    vmax: float | None = None,
    dot_min: float = 0.0,
    dot_max: float = 1.0,
    smallest_dot: float = 10.0,
    largest_dot: float = 350.0,
    figsize: tuple[float, float] = (6.8, 4.8),
    size_title: str = "Fraction of samples\nin group (%)",
    colorbar_title: str = "Mean expression\nin group",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scanpy-like dotplot for bulk data.
      - dot color = mean expression per group
      - dot size  = fraction of samples expressing (> expr_threshold) per group

    Notes:
      - groupby must be categorical-ish in adata.obs
      - var_names must be in adata.var_names
      - dendrogram_* requires scipy
    """
    set_style()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    var_names = [str(v) for v in var_names]
    missing = [g for g in var_names if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not found in adata.var_names (first 10): {missing[:10]}")

    # --- extract matrix ---
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    gidx = [adata.var_names.get_loc(g) for g in var_names]
    if sp.issparse(X):
        M = X[:, gidx].toarray()
    else:
        M = np.asarray(X[:, gidx], dtype=float)

    groups = adata.obs[groupby].astype("category")
    cats = list(groups.cat.categories)

    mean_expr = np.zeros((len(cats), len(var_names)), dtype=float)
    frac_expr = np.zeros((len(cats), len(var_names)), dtype=float)

    for i, c in enumerate(cats):
        mask = (groups == c).to_numpy()
        Xi = M[mask, :]
        mean_expr[i, :] = Xi.mean(axis=0)
        frac_expr[i, :] = (Xi > expr_threshold).mean(axis=0)

    # --- scale like scanpy options ---
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

    # --- optional dendrogram ordering ---
    row_order = np.arange(len(cats))
    col_order = np.arange(len(var_names))

    if (dendrogram_rows or dendrogram_cols) and (linkage is None or dendrogram is None):
        raise ImportError("dendrogram_rows/cols requires scipy (scipy.cluster.hierarchy).")

    if dendrogram_rows and len(cats) > 2:
        Zr = linkage(disp, method="average", metric="euclidean")
        row_order = dendrogram(Zr, no_plot=True)["leaves"]

    if dendrogram_cols and len(var_names) > 2:
        Zc = linkage(disp.T, method="average", metric="euclidean")
        col_order = dendrogram(Zc, no_plot=True)["leaves"]

    disp = disp[row_order, :][:, col_order]
    frac_expr = frac_expr[row_order, :][:, col_order]
    cats_ord = [cats[i] for i in row_order]
    genes_ord = [var_names[i] for i in col_order]

    # --- dot size scaling (fraction -> area) ---
    f = np.clip(frac_expr, dot_min, dot_max)
    # map dot_min..dot_max to [smallest_dot..largest_dot]
    if dot_max <= dot_min:
        raise ValueError("dot_max must be > dot_min")
    s = smallest_dot + (largest_dot - smallest_dot) * ((f - dot_min) / (dot_max - dot_min))

    # --- color scaling ---
    if vmin is None:
        vmin = np.nanmin(disp)
    if vmax is None:
        vmax = np.nanmax(disp)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap_obj = mpl.cm.get_cmap(cmap)

    # --- layout (Scanpy-like): main + row dendro + size legend + colorbar ---
    fig = plt.figure(figsize=figsize, constrained_layout=False)
    gs = GridSpec(
        nrows=2,
        ncols=3,
        figure=fig,
        height_ratios=[1.0, 0.08],   # main + colorbar row
        width_ratios=[1.0, 0.18, 0.34],  # main, dendro (optional), legend
        wspace=0.05,
        hspace=0.25,
    )

    ax = fig.add_subplot(gs[0, 0])

    ax_dendro = fig.add_subplot(gs[0, 1]) if dendrogram_rows else None
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_cbar = fig.add_subplot(gs[1, 2])

    # --- draw dots ---
    xs = np.arange(len(genes_ord))
    ys = np.arange(len(cats_ord))
    for i in range(len(cats_ord)):
        ax.scatter(
            xs,
            np.full_like(xs, ys[i]),
            s=s[i, :],
            c=cmap_obj(norm(disp[i, :])),
            edgecolors="none",
        )

    ax.set_xlim(-0.5, len(genes_ord) - 0.5)
    ax.set_ylim(-0.5, len(cats_ord) - 0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(genes_ord, rotation=90)
    ax.set_yticks(ys)
    ax.set_yticklabels([str(c) for c in cats_ord])

    # Scanpy-like “categories on top”
    ax.tick_params(axis="x", labeltop=True, labelbottom=True)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # --- optional row dendrogram on the right ---
    if dendrogram_rows and ax_dendro is not None and len(cats_ord) > 2:
        Zr = linkage(disp, method="average", metric="euclidean")
        dendrogram(Zr, orientation="right", ax=ax_dendro, color_threshold=None, no_labels=True)
        ax_dendro.set_xticks([])
        ax_dendro.set_yticks([])
        for spine in ax_dendro.spines.values():
            spine.set_visible(False)

    # --- size legend (Scanpy-like) ---
    ax_leg.axis("off")
    ax_leg.set_title(size_title, loc="left")

    # show size legend as 5 reference dots
    ref = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
    ref_s = smallest_dot + (largest_dot - smallest_dot) * ((ref - dot_min) / (dot_max - dot_min))
    x0, y0 = 0.15, 0.75
    dx = 0.16
    for j, (r, rs) in enumerate(zip(ref, ref_s)):
        ax_leg.scatter([x0 + j * dx], [y0], s=rs, color="gray", edgecolors="none")
        ax_leg.text(x0 + j * dx, y0 - 0.10, f"{int(round(r*100))}", ha="center", va="top")
    ax_leg.text(x0 - 0.03, y0 - 0.22, "%", ha="right", va="top")

    # --- colorbar ---
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    cb = fig.colorbar(sm, cax=ax_cbar, orientation="horizontal")
    cb.set_label(colorbar_title)

    # --- gene group brackets (Scanpy-ish) ---
    if var_group_positions is not None:
        # var_group_positions are in original var_names indices; convert to new col_order
        # mapping: original index -> new index
        inv = {orig_i: new_i for new_i, orig_i in enumerate(col_order)}
        topy = len(cats_ord) - 0.2  # near top
        for start, end, label in var_group_positions:
            # end is exclusive like python slicing
            start_n = inv.get(int(start), None)
            end_n = inv.get(int(end - 1), None)
            if start_n is None or end_n is None:
                continue
            # bracket coordinates in new order
            x1 = start_n - 0.5
            x2 = end_n + 0.5
            y = len(cats_ord) - 0.5 + 0.25
            ax.plot([x1, x1, x2, x2], [y, y + 0.2, y + 0.2, y], lw=1, color="black", clip_on=False)
            ax.text((x1 + x2) / 2, y + 0.25, label, ha="center", va="bottom", rotation=90)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax