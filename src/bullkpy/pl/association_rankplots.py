from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from ..logging import warn
from .dotplot import _font_pt, _text_extent_in
from ._style import set_style, _savefig
from .._deprecated import renamed_kwargs


def rankplot_association(
    *,
    res: pd.DataFrame,
    gene_col: str = "gene",
    effect_col: str = "effect",      # e.g. log2FC, eta2, mean_diff
    sort_by: str = "qval",           # "qval" | "pval" | effect_col
    direction: str = "both",         # "up" | "down" | "both"
    n_items: int = 20,
    figsize: tuple[float, float] = (7, 6),
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Ranked barplot for categorical association results.

    Behavior mirrors bk.pl.rankplot():
      • Upregulated / positive = red
      • Downregulated / negative = blue
      • Strongest up at top
      • Most downregulated LAST
    """
    set_style()

    if gene_col not in res.columns:
        raise KeyError(f"'{gene_col}' not found in result table.")
    if effect_col not in res.columns:
        raise KeyError(f"'{effect_col}' not found in result table.")
    if sort_by not in res.columns:
        raise KeyError(f"sort_by='{sort_by}' not found in result table.")

    df = res.copy()
    df = df.dropna(subset=[gene_col, effect_col, sort_by]).copy()
    df[effect_col] = pd.to_numeric(df[effect_col], errors="coerce")
    df = df.dropna(subset=[effect_col]).copy()

    n_items = int(n_items)
    if n_items <= 0:
        raise ValueError("n_items must be > 0")

    # Base ordering
    if sort_by == effect_col:
        df_sorted = df
    else:
        df_sorted = df.sort_values(sort_by, ascending=True)

    if direction == "up":
        sub = df_sorted[df_sorted[effect_col] > 0].copy()
        if sort_by == effect_col:
            sub = sub.sort_values(effect_col, ascending=False).head(n_items)
        else:
            sub = sub.head(n_items).sort_values(effect_col, ascending=False)

        # strongest up at top
        sub = sub.sort_values(effect_col, ascending=False)

    elif direction == "down":
        sub = df_sorted[df_sorted[effect_col] < 0].copy()
        if sort_by == effect_col:
            sub = sub.sort_values(effect_col, ascending=True).head(n_items)
        else:
            sub = sub.head(n_items)

        # most downregulated LAST
        sub = sub.sort_values(effect_col, ascending=False)

    elif direction == "both":
        n_up = n_items // 2
        n_down = n_items - n_up

        up = df_sorted[df_sorted[effect_col] > 0].copy()
        down = df_sorted[df_sorted[effect_col] < 0].copy()

        if sort_by == effect_col:
            up = up.sort_values(effect_col, ascending=False).head(n_up)
            down = down.sort_values(effect_col, ascending=True).head(n_down)
        else:
            up = up.head(n_up)
            down = down.head(n_down)

        # enforce final display order
        up = up.sort_values(effect_col, ascending=False)
        down = down.sort_values(effect_col, ascending=False)

        sub = pd.concat([up, down], axis=0)

    else:
        raise ValueError("direction must be 'up', 'down', or 'both'.")

    labels = sub[gene_col].astype(str).tolist()
    vals = sub[effect_col].to_numpy(dtype=float)

    colors = np.where(vals >= 0, "#D62728", "#1F77B4")  # red / blue

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    ax.barh(labels, vals, color=colors)
    ax.axvline(0, linewidth=1)

    ax.set_xlabel(effect_col)
    ax.set_ylabel(gene_col)
    ax.invert_yaxis()

    ax.set_title(title or "Ranked associations")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax

@renamed_kwargs(q_col="qval_col")
def dotplot_association(
    df: pd.DataFrame,
    *,
    feature_col: str,            # "gene" or "obs"
    groupby_col: str = "groupby",
    effect_col: str = "effect",
    qval_col: str = "qval",
    top_n: int = 50,
    figsize: tuple[float, float] | None = None,
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    size_min: float = 10.0,
    size_max: float = 250.0,
    pad: float = 0.7,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Scanpy-like dotplot for association results across multiple contrasts/groupby runs.
      - dot color = effect
      - dot size  = -log10(qval)
    Works when df contains multiple groupby/contrasts (groupby_col), otherwise still works.

    ``pad`` is the empty space around the outermost dots, in units of one dot cell. ``size_max`` is an
    upper bound: dots are also capped to the size of their cell so they never overlap or touch the frame.
    A size legend shows -log10(q) for four reference dots; the default title is omitted because the
    colour bar and legend already say what colour and size mean.
    """
    set_style()

    d = df.copy()
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=[feature_col, groupby_col, effect_col, qval_col])

    # take top_n per groupby_col
    d = d.sort_values([groupby_col, qval_col], ascending=[True, True])
    d = d.groupby(groupby_col, group_keys=False).head(int(top_n))

    piv_eff = d.pivot_table(index=feature_col, columns=groupby_col, values=effect_col, aggfunc="first")
    piv_q = d.pivot_table(index=feature_col, columns=groupby_col, values=qval_col, aggfunc="first")

    eff = piv_eff.to_numpy(dtype=float)
    qq = piv_q.to_numpy(dtype=float)
    siz = -np.log10(np.maximum(qq, 1e-300))

    n_rows, n_cols = piv_eff.shape
    smin, smax = np.nanmin(siz), np.nanmax(siz)

    own_axes = ax is None
    if own_axes:
        # smallest canvas that holds the labels, a readable panel, the colour bar and the size legend
        font = _font_pt(mpl.rcParams["ytick.labelsize"])
        ylab_w = _text_extent_in([str(r) for r in piv_eff.index], font) + 0.15
        xlab_h = _text_extent_in([str(c) for c in piv_eff.columns], font) + 0.15
        panel_w = (n_cols - 1 + 2 * float(pad)) * 0.28
        panel_h = (n_rows - 1 + 2 * float(pad)) * 0.28
        min_w = ylab_w + panel_w + 0.30 + 0.95 + 0.45 + max(_text_extent_in([f"-log10({qval_col})"], font), 0.7)
        min_h = max(xlab_h + panel_h + 0.25, 0.9 + 1.9)  # colour bar above the size legend
        if figsize is None:
            figsize = (max(min_w, 0.40 * n_cols + 2.6 + ylab_w), max(min_h, 0.28 * n_rows + xlab_h + 0.4))
        elif figsize[0] < min_w - 1e-6 or figsize[1] < min_h - 1e-6:
            new = (max(figsize[0], min_w), max(figsize[1], min_h))
            warn(
                f"dotplot_association: figsize=({figsize[0]:g}, {figsize[1]:g}) is too small for the labels, colour "
                f"bar and size legend; using ({new[0]:.1f}, {new[1]:.1f}). Omit figsize to size the plot automatically."
            )
            figsize = new
    elif figsize is None:
        figsize = (4.0, 3.0)

    if own_axes:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    if vmin is None:
        vmin = float(np.nanmin(eff))
    if vmax is None:
        vmax = float(np.nanmax(eff))
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    norm = (
        mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        if (vmin < 0 < vmax)
        else mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    )

    xs = np.arange(n_cols)
    ys = np.arange(n_rows)
    ax.set_xlim(-float(pad), n_cols - 1 + float(pad))
    ax.set_ylim(n_rows - 1 + float(pad), -float(pad))  # first row on top

    # Cap the dots to their cell: draw once so the (constrained) layout is final, then measure.
    fig.canvas.draw()
    box = ax.get_window_extent()
    cell_pt = (
        min(
            box.width / max(n_cols - 1 + 2 * float(pad), 1.0),
            box.height / max(n_rows - 1 + 2 * float(pad), 1.0),
        )
        * 72.0
        / fig.dpi
    )
    size_max_eff = min(float(size_max), (0.9 * cell_pt) ** 2)
    size_min_eff = min(float(size_min), 0.25 * size_max_eff)

    def _area(v):
        v = np.asarray(v, dtype=float)
        if not (np.isfinite(smin) and np.isfinite(smax)) or smax == smin:
            return np.full_like(v, size_max_eff)
        frac = np.clip((v - smin) / (smax - smin), 0.0, 1.0)
        return size_min_eff + (size_max_eff - size_min_eff) * frac

    X, Y = np.meshgrid(xs, ys)
    show_cell = np.isfinite(eff) & np.isfinite(siz)
    sc = ax.scatter(
        X[show_cell],
        Y[show_cell],
        s=_area(siz[show_cell]),
        c=eff[show_cell],
        cmap=mpl.colormaps.get_cmap(cmap),
        norm=norm,
        edgecolors="0.2",
        linewidths=0.3,
    )

    ax.set_xticks(xs)
    ax.set_xticklabels([str(c) for c in piv_eff.columns], rotation=90)
    ax.set_yticks(ys)
    ax.set_yticklabels([str(r) for r in piv_eff.index])
    ax.grid(False)

    cbar_shrink = min(0.6, 1.5 / max(box.height / fig.dpi, 1e-6))  # at most 1.5 in tall
    cbar = fig.colorbar(sc, ax=ax, pad=0.02, shrink=cbar_shrink, aspect=12)
    cbar.set_label(effect_col)

    # size legend: reference dots spanning the observed -log10(q) range
    if np.isfinite(smin) and np.isfinite(smax):
        ref = np.unique(np.round(np.linspace(smin, smax, 4), 1))
    else:
        ref = np.array([])
    handles = [
        ax.scatter([], [], s=float(_area(v)), color="0.6", edgecolors="0.2", linewidths=0.3, label=f"{v:g}")
        for v in ref
    ]
    if handles:
        legend_kw = dict(handles=handles, title=f"-log10({qval_col})", frameon=False, labelspacing=1.1, borderpad=0.4)
        if own_axes:
            fig.legend(loc="outside right lower", **legend_kw)
        else:
            ax.legend(loc="upper left", bbox_to_anchor=(1.25, 0.0), **legend_kw)

    if title:
        ax.set_title(title)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


def heatmap_association(
    df: pd.DataFrame,
    *,
    feature_col: str,
    groupby_col: str = "groupby",
    value_col: str = "effect",
    top_n: int = 60,
    cmap: str = "RdBu_r",
    center: float = 0.0,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Heatmap of association values (effect by default), selecting top_n rows by best qval per column.
    """
    set_style()

    d = df.copy()
    if "qval" in d.columns:
        d = d.sort_values(["qval"], ascending=True)
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=[feature_col, groupby_col, value_col])

    # pick top_n unique features overall (simple + robust)
    feats = d[feature_col].astype(str).drop_duplicates().head(int(top_n)).tolist()
    d = d[d[feature_col].astype(str).isin(feats)]

    piv = d.pivot_table(index=feature_col, columns=groupby_col, values=value_col, aggfunc="first")
    mat = piv.to_numpy(dtype=float)

    if figsize is None:
        figsize = (max(6, 0.45 * piv.shape[1] + 3), max(4, 0.22 * piv.shape[0] + 1.8))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    im = ax.imshow(mat, aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels([str(c) for c in piv.columns], rotation=90)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels([str(r) for r in piv.index])
    ax.invert_yaxis()

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label(value_col)

    ax.set_title(title or f"Association heatmap ({value_col})")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax