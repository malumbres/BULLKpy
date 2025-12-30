from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata as ad

from ._style import set_style, _savefig


def library_size_vs_genes(
    adata: ad.AnnData,
    *,
    x: str = "total_counts",
    y: str = "n_genes_detected",
    groupby: str | None = None,
    # thresholds
    min_counts: float | None = None,
    max_counts: float | None = None,
    min_genes: float | None = None,
    max_genes: float | None = None,
    # display
    logx: bool = True,
    logy: bool = True,
    s: float = 18.0,
    alpha: float = 0.85,
    linewidth: float = 0.25,
    edgecolor: str = "0.15",
    # highlight outliers
    show_outliers: bool = True,
    outlier_color: str = "crimson",
    outlier_marker: str = "x",
    outlier_size: float = 28.0,
    # labels
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    # figure
    figsize: tuple[float, float] = (5.5, 4.5),
    legend: bool = True,
    legend_loc: str = "best",
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    QC scatter: library size vs detected genes with optional thresholds.

    Requires columns in adata.obs:
      - x (default: total_counts)
      - y (default: n_genes_detected)

    Typical usage:
      bk.pl.library_size_vs_genes(adata, min_counts=1e6, min_genes=12000, groupby="Subtype")
    """
    set_style()

    if x not in adata.obs.columns:
        raise KeyError(f"adata.obs['{x}'] not found")
    if y not in adata.obs.columns:
        raise KeyError(f"adata.obs['{y}'] not found")

    df = adata.obs[[x, y]].copy()
    df[x] = pd.to_numeric(df[x], errors="coerce")
    df[y] = pd.to_numeric(df[y], errors="coerce")
    df = df.dropna(subset=[x, y])

    # compute pass/fail mask
    ok = np.ones(len(df), dtype=bool)

    if min_counts is not None:
        ok &= df[x].to_numpy() >= float(min_counts)
    if max_counts is not None:
        ok &= df[x].to_numpy() <= float(max_counts)
    if min_genes is not None:
        ok &= df[y].to_numpy() >= float(min_genes)
    if max_genes is not None:
        ok &= df[y].to_numpy() <= float(max_genes)

    # handle grouping/palette via seaborn set_palette in set_style()
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"adata.obs['{groupby}'] not found")
        df[groupby] = adata.obs.loc[df.index, groupby].astype(str)

    fig, ax = plt.subplots(figsize=figsize)

    # scatter (inliers)
    if groupby is None:
        ax.scatter(
            df.loc[ok, x],
            df.loc[ok, y],
            s=s,
            alpha=alpha,
            linewidths=linewidth,
            edgecolors=edgecolor,
        )
    else:
        # plot per category (keeps legend manageable and stable)
        cats = pd.Categorical(df[groupby]).categories
        for c in cats:
            m = ok & (df[groupby] == c).to_numpy()
            if m.sum() == 0:
                continue
            ax.scatter(
                df.loc[m, x],
                df.loc[m, y],
                s=s,
                alpha=alpha,
                linewidths=linewidth,
                edgecolors=edgecolor,
                label=str(c),
            )

    # outliers on top
    if show_outliers and (~ok).any():
        ax.scatter(
            df.loc[~ok, x],
            df.loc[~ok, y],
            s=outlier_size,
            alpha=0.95,
            marker=outlier_marker,
            c=outlier_color,
            linewidths=0.8,
            label="QC fail" if groupby is not None else None,
        )

    # log scales
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    # threshold lines
    def _vline(val):
        ax.axvline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    def _hline(val):
        ax.axhline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    if min_counts is not None:
        _vline(min_counts)
    if max_counts is not None:
        _vline(max_counts)
    if min_genes is not None:
        _hline(min_genes)
    if max_genes is not None:
        _hline(max_genes)

    # labels
    ax.set_xlabel(xlabel if xlabel is not None else x)
    ax.set_ylabel(ylabel if ylabel is not None else y)

    if title is None:
        n_fail = int((~ok).sum())
        title = f"{y} vs {x} (QC fail: {n_fail})" if (min_counts or max_counts or min_genes or max_genes) else f"{y} vs {x}"
    ax.set_title(title)

    if legend and groupby is not None:
        ax.legend(loc=legend_loc, frameon=False, ncol=1)

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax