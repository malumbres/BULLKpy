from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad


from ._style import set_style, _savefig
from .._compat import group_labels as _group_labels


def association_heatmap(
    df: pd.DataFrame,
    *,
    index: str,
    columns: str,
    values: str,
    agg: Literal["mean", "max", "min"] = "mean",
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Generic heatmap for association outputs.

    Parameters
    ----------
    df
        Long-form association table (tidy dataframe).
    index, columns, values
        Column names used to build the pivot (heatmap matrix).
    agg
        Aggregation used when multiple rows map to the same (index, column) cell.
    cmap, vmin, vmax
        Matplotlib colormap and scaling.
    figsize
        Figure size in inches. If None, auto-sized from matrix shape.
    title
        Plot title.
    save
        Path to save the figure (png/pdf/svg).
    show
        Whether to display the figure.

    Returns
    -------
    fig, ax
    """
    set_style()

    piv = pd.pivot_table(df, index=index, columns=columns, values=values, aggfunc=agg)
    mat = piv.to_numpy(dtype=float)

    if figsize is None:
        figsize = (max(6, 0.25 * piv.shape[1] + 3), max(4, 0.25 * piv.shape[0] + 2))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels([str(c) for c in piv.columns], rotation=90)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels([str(r) for r in piv.index])

    ax.set_title(title or f"Heatmap: {values}")

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label(values)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


def boxplot_with_stats(
    adata: ad.AnnData,
    *,
    y: str,
    groupby: str,
    figsize: tuple[float, float] = (7, 3.5),
    kind: Literal["box", "violin"] = "violin",
    show_points: bool = True,
    point_size: float = 2.0,
    point_alpha: float = 0.3,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Box/violin plot of a numeric obs column across categorical groups, with a
    simple global p-value annotation.

    - If 2 groups: Mann–Whitney U
    - If >2 groups: Kruskal–Wallis

    Returns (fig, ax).
    """
    set_style()

    if y not in adata.obs.columns:
        raise KeyError(f"'{y}' not in adata.obs")
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' not in adata.obs")

    s = pd.to_numeric(adata.obs[y], errors="coerce")
    g = _group_labels(adata.obs[groupby])
    df = pd.DataFrame({"y": s, "g": g}).dropna()

    cats = list(pd.Categorical(df["g"]).categories)
    groups = [df.loc[df["g"] == c, "y"].to_numpy(dtype=float) for c in cats]
    k_eff = sum(v.size > 0 for v in groups)

    pval = np.nan
    if k_eff >= 2:
        if len(cats) == 2:
            from scipy.stats import mannwhitneyu
            pval = mannwhitneyu(groups[0], groups[1], alternative="two-sided").pvalue
        else:
            from scipy.stats import kruskal
            pval = kruskal(*[v for v in groups if v.size > 0]).pvalue

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    positions = np.arange(len(cats))
    if kind == "box":
        ax.boxplot(groups, positions=positions, showfliers=False)
    else:
        parts = ax.violinplot(groups, positions=positions, showmeans=False, showextrema=False, showmedians=True)
        for pc in parts.get("bodies", []):
            pc.set_alpha(0.8)

    if show_points:
        rng = np.random.RandomState(0)
        for i, v in enumerate(groups):
            if v.size == 0:
                continue
            xj = rng.normal(loc=i, scale=0.06, size=v.size)
            ax.scatter(xj, v, s=point_size, alpha=point_alpha, edgecolors="none")

    ax.set_xticks(positions)
    ax.set_xticklabels(cats, rotation=90)
    ax.set_ylabel(y)

    ttl = title or f"{y} by {groupby}"
    if np.isfinite(pval):
        ttl += f"  (p={pval:.2e})"
    ax.set_title(ttl)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


