from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from ._style import set_style, _savefig


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
    show: bool = True,
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

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
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

def dotplot_association(
    df: pd.DataFrame,
    *,
    feature_col: str,            # "gene" or "obs"
    groupby_col: str = "groupby",
    effect_col: str = "effect",
    q_col: str = "qval",
    top_n: int = 50,
    figsize: tuple[float, float] | None = None,
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    size_min: float = 10.0,
    size_max: float = 250.0,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scanpy-like dotplot for association results across multiple contrasts/groupby runs.
      - dot color = effect
      - dot size  = -log10(qval)
    Works when df contains multiple groupby/contrasts (groupby_col), otherwise still works.
    """
    set_style()

    d = df.copy()
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=[feature_col, groupby_col, effect_col, q_col])

    # take top_n per groupby_col
    d = d.sort_values([groupby_col, q_col], ascending=[True, True])
    d = d.groupby(groupby_col, group_keys=False).head(int(top_n))

    piv_eff = d.pivot_table(index=feature_col, columns=groupby_col, values=effect_col, aggfunc="first")
    piv_q = d.pivot_table(index=feature_col, columns=groupby_col, values=q_col, aggfunc="first")

    eff = piv_eff.to_numpy(dtype=float)
    qq = piv_q.to_numpy(dtype=float)
    siz = -np.log10(np.maximum(qq, 1e-300))

    # scale size
    smin, smax = np.nanmin(siz), np.nanmax(siz)
    if not np.isfinite(smin) or not np.isfinite(smax) or smax == smin:
        u = np.zeros_like(siz)
    else:
        u = (siz - smin) / (smax - smin)
    sizes = size_min + (size_max - size_min) * u

    if figsize is None:
        figsize = (max(6, 0.45 * piv_eff.shape[1] + 2.8), max(4, 0.22 * piv_eff.shape[0] + 1.6))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    if vmin is None:
        vmin = float(np.nanmin(eff))
    if vmax is None:
        vmax = float(np.nanmax(eff))
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax) if (vmin < 0 < vmax) else mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cm = mpl.cm.get_cmap(cmap)

    xs = np.arange(piv_eff.shape[1])
    ys = np.arange(piv_eff.shape[0])

    for i in range(piv_eff.shape[0]):
        ax.scatter(
            xs,
            np.full_like(xs, ys[i]),
            s=sizes[i, :],
            c=cm(norm(eff[i, :])),
            edgecolors="0.2",
            linewidths=0.3,
        )

    ax.set_xticks(xs)
    ax.set_xticklabels([str(c) for c in piv_eff.columns], rotation=90)
    ax.set_yticks(ys)
    ax.set_yticklabels([str(r) for r in piv_eff.index])
    ax.invert_yaxis()

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cm)
    cbar = fig.colorbar(sm, ax=ax, pad=0.01)
    cbar.set_label(effect_col)

    ax.set_title(title or f"Association dotplot: color={effect_col}, size=-log10({q_col})")

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
    show: bool = True,
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

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
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