from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from ._style import set_style, _savefig


def gsea_bubbleplot(
    df_gsea: pd.DataFrame,
    *,
    pathways: Mapping[str, Sequence[str]] | Sequence[str],
    comparison_col: str = "comparison",
    term_col: str = "Term",
    nes_col: str = "NES",
    fdr_col: str = "FDR q-val",
    # ordering
    comparison_order: Sequence[str] | None = None,
    drop_empty_comparisons: bool = True,
    # size mapping
    size_from: str = "fdr",  # "fdr" or "pval" (if you want)
    min_q: float = 1e-300,
    size_min: float = 10.0,
    size_max: float = 350.0,
    fdr_floor: float = 1e-50,
    size_clip_quantile: float | None = 0.99,   # None disables
    # color mapping
    cmap: str = "RdBu_r",
    center: float = 0.0,
    vmin: float | None = None,
    vmax: float | None = None,
    # cosmetics
    figsize: tuple[float, float] | None = None,
    row_spacing: float = 1.0,
    col_spacing: float = 1.0,
    row_height: float = 0.32,                  # inches per row (Scanpy-ish)
    col_width: float = 0.32,                   # inches per column
    dot_edgecolor: str = "0.15",
    dot_linewidth: float = 0.35,
    show_grid: bool = False,
    group_label_rotation: float = 90,
    xtick_rotation: float = 90,
    title: str | None = None,
    # output
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Bubble plot matrix for GSEA results.

    Rows: comparisons (contrasts)
    Cols: pathways (terms)
    Color: NES (diverging, centered at `center`)
    Size: -log10(FDR q-val) with floor & optional clipping

    `pathways` can be:
      - dict: {"Immune": [term1, term2], "Metabolism": [term3]}
      - list: [term1, term2, ...]
    """
    set_style()

    if not isinstance(df_gsea, pd.DataFrame):
        raise TypeError("df_gsea must be a pandas DataFrame.")
    for c in (comparison_col, term_col, nes_col, fdr_col):
        if c not in df_gsea.columns:
            raise KeyError(
                f"'{c}' not found in df_gsea columns. Available: {list(df_gsea.columns)[:30]} ..."
            )

    # ---- flatten pathways + keep group spans ----
    if isinstance(pathways, Mapping):
        terms: list[str] = []
        spans: list[tuple[int, int, str]] = []
        start = 0
        for grp, lst in pathways.items():
            lst = [str(x) for x in lst]
            terms.extend(lst)
            end = start + len(lst)
            spans.append((start, end, str(grp)))
            start = end
    else:
        terms = [str(x) for x in pathways]
        spans = []

    if len(terms) == 0:
        raise ValueError("No pathways provided.")

    # ---- subset and make pivot matrices ----
    sub = df_gsea[df_gsea[term_col].astype(str).isin(terms)].copy()
    sub[term_col] = sub[term_col].astype(str)
    sub[comparison_col] = sub[comparison_col].astype(str)

    if comparison_order is None:
        comps = list(pd.Categorical(sub[comparison_col]).categories) if sub.shape[0] else []
        if not comps:
            comps = sorted(df_gsea[comparison_col].astype(str).unique().tolist())
    else:
        comps = [str(x) for x in comparison_order]

    nes_mat = sub.pivot_table(index=comparison_col, columns=term_col, values=nes_col, aggfunc="mean")
    q_mat = sub.pivot_table(index=comparison_col, columns=term_col, values=fdr_col, aggfunc="mean")

    nes_mat = nes_mat.reindex(index=comps, columns=terms)
    q_mat = q_mat.reindex(index=comps, columns=terms)

    if drop_empty_comparisons:
        keep = ~(nes_mat.isna().all(axis=1))
        nes_mat = nes_mat.loc[keep]
        q_mat = q_mat.loc[keep]

    comps_final = nes_mat.index.tolist()
    if len(comps_final) == 0:
        raise ValueError("No comparisons have any of the selected pathways in df_gsea (after filtering).")

    # ---- size: -log10(q) with floor & clipping ----
    q_vals = q_mat.to_numpy(dtype=float)
    q_vals = np.where(np.isfinite(q_vals), q_vals, np.nan)

    # floor: treat 0 or negative as fdr_floor
    q_vals = np.where(q_vals <= 0, fdr_floor, q_vals)
    q_vals = np.clip(q_vals, fdr_floor, 1.0)

    size_signal = -np.log10(q_vals)

    if size_clip_quantile is not None:
        cap = np.nanquantile(size_signal[np.isfinite(size_signal)], float(size_clip_quantile))
        size_signal = np.minimum(size_signal, cap)

    finite = np.isfinite(size_signal)
    if finite.any():
        smin = float(np.nanmin(size_signal[finite]))
        smax = float(np.nanmax(size_signal[finite]))
        if smax == smin:
            smax = smin + 1.0
        u = (size_signal - smin) / (smax - smin)
        sizes = size_min + (size_max - size_min) * np.clip(u, 0, 1)
    else:
        smin, smax = 0.0, 1.0
        sizes = np.full_like(size_signal, size_min, dtype=float)

    # missing NES → do not draw dot
    nes_vals = nes_mat.to_numpy(dtype=float)
    sizes[~np.isfinite(nes_vals)] = 0.0

    # ---- color scaling (CRITICAL: same norm used by scatter + colorbar) ----
    nes_all = nes_vals[np.isfinite(nes_vals)]
    if nes_all.size == 0:
        # fallback
        vmin_eff, vmax_eff = -1.0, 1.0
    else:
        if vmin is None or vmax is None:
            vmax_abs = float(np.nanmax(np.abs(nes_all)))
            vmax_abs = vmax_abs if vmax_abs > 0 else 1.0
            vmin_eff, vmax_eff = center - vmax_abs, center + vmax_abs
        else:
            vmin_eff, vmax_eff = float(vmin), float(vmax)

            # if user gives asymmetric bounds but wants centering, you can still
            # keep TwoSlopeNorm; it will work, but not symmetric.
            # We'll keep as given.
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin_eff, vcenter=float(center), vmax=vmax_eff)
    cmap_obj = mpl.cm.get_cmap(cmap)

    # ---- autosize figure ----
    if figsize is None:
        w = max(6.0, col_width * len(terms) + 3.2)          # room for legends
        h = max(3.5, row_height * len(comps_final) + 2.0)   # room for labels
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=False)
    fig.subplots_adjust(right=0.78)

    # coords
    xs = np.arange(len(terms), dtype=float) * float(col_spacing)
    ys = np.arange(len(comps_final), dtype=float) * float(row_spacing)

    # ---- draw dots (IMPORTANT: c must be numeric NES, not RGBA) ----
    # We draw all points in one scatter so colorbar matches exactly.
    XX, YY = np.meshgrid(xs, ys)
    XX = XX.ravel()
    YY = YY.ravel()
    S = sizes.ravel()
    C = nes_vals.ravel()

    mask = np.isfinite(C) & (S > 0)
    sc = ax.scatter(
        XX[mask],
        YY[mask],
        s=S[mask],
        c=C[mask],                 # numeric NES
        cmap=cmap_obj,
        norm=norm,
        edgecolors=dot_edgecolor,
        linewidths=dot_linewidth,
    )

    # axes ticks
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.set_xticklabels(terms, rotation=xtick_rotation, ha="right")
    ax.set_yticklabels(comps_final)

    # limits
    if len(xs) == 1:
        ax.set_xlim(xs[0] - 0.6 * col_spacing, xs[0] + 0.6 * col_spacing)
    else:
        ax.set_xlim(xs.min() - 0.6 * col_spacing, xs.max() + 0.6 * col_spacing)

    if len(ys) == 1:
        ax.set_ylim(ys[0] - 0.6 * row_spacing, ys[0] + 0.6 * row_spacing)
    else:
        ax.set_ylim(ys.min() - 0.6 * row_spacing, ys.max() + 0.6 * row_spacing)

    ax.invert_yaxis()  # scanpy-like

    if show_grid:
        ax.grid(True, linewidth=0.4, alpha=0.35)
    else:
        ax.grid(False)

    ax.set_xlabel("")
    ax.set_ylabel("")
    if title:
        ax.set_title(title)

    # ---- pathway group brackets (if dict provided) ----
    if spans:
        top_y = ys.min() - 1.1 * row_spacing
        for start, end, label in spans:
            x1 = xs[start] - 0.5 * col_spacing
            x2 = xs[end - 1] + 0.5 * col_spacing
            ax.plot(
                [x1, x1, x2, x2],
                [top_y + 0.2 * row_spacing, top_y, top_y, top_y + 0.2 * row_spacing],
                lw=1.0,
                color="0.2",
                clip_on=False,
            )
            ax.text(
                (x1 + x2) / 2,
                top_y - 0.15 * row_spacing,
                label,
                ha="center",
                va="bottom",
                rotation=group_label_rotation,
                clip_on=False,
            )

    # ---- colorbar (now matches bubbles) ----
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cax = inset_axes(
        ax,
        width="3%",
        height="55%",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.25, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label("NES")

    # ---- size legend (uses the same scaling, including floor/clipping) ----
    ref_q = np.array([0.05, 0.01, 0.001, fdr_floor], dtype=float)
    ref_q = np.clip(ref_q, fdr_floor, 1.0)
    ref_sig = -np.log10(ref_q)

    if size_clip_quantile is not None and finite.any():
        cap = np.nanquantile((-np.log10(q_vals[np.isfinite(q_vals)])), float(size_clip_quantile))
        ref_sig = np.minimum(ref_sig, cap)

    if finite.any():
        uref = (ref_sig - smin) / (smax - smin)
        sref = size_min + (size_max - size_min) * np.clip(uref, 0, 1)
    else:
        sref = np.full_like(ref_sig, size_min, dtype=float)

    labels = [f"q={q:g}" for q in ref_q[:-1]] + [f"q≤{fdr_floor:g}"]
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="0.6",
            markeredgecolor=dot_edgecolor,
            markersize=float(np.sqrt(sr)),   # area->approx marker size
            label=lab,
        )
        for sr, lab in zip(sref, labels)
    ]

    ax.legend(
        handles=handles,
        title="-log10(FDR)",
        bbox_to_anchor=(1.02, 0.22),
        loc="upper left",
        frameon=False,
    )

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax