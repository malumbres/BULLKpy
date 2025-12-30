from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import anndata as ad

from ._style import set_style, _savefig
from ..tl.association import categorical_association


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
    Example: genes x categories using effect or -log10(qval).
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
    y: str,                       # numeric obs key
    groupby: str,                 # categorical obs key
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
    Simple scanpy-like box/violin plot annotated with global p-value (Kruskal for >2 groups, MWU for 2).
    """
    set_style()

    if y not in adata.obs.columns:
        raise KeyError(f"'{y}' not in adata.obs")
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' not in adata.obs")

    s = pd.to_numeric(adata.obs[y], errors="coerce")
    g = adata.obs[groupby].astype(str)
    df = pd.DataFrame({"y": s, "g": g}).dropna()

    cats = list(pd.Categorical(df["g"]).categories)
    groups = [df.loc[df["g"] == c, "y"].to_numpy(dtype=float) for c in cats]
    k_eff = sum(v.size > 0 for v in groups)

    # compute p-value
    pval = np.nan
    if k_eff >= 2:
        if len(cats) == 2:
            from scipy.stats import mannwhitneyu
            pval = mannwhitneyu(groups[0], groups[1], alternative="two-sided").pvalue
        else:
            from scipy.stats import kruskal
            pval = kruskal(*[v for v in groups if v.size > 0]).pvalue

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # lightweight plotting without seaborn dependency
    # (if you prefer seaborn, we can switch)
    positions = np.arange(len(cats))
    if kind == "box":
        ax.boxplot(groups, positions=positions, showfliers=False)
    else:
        parts = ax.violinplot(groups, positions=positions, showmeans=False, showextrema=False, showmedians=True)
        for pc in parts["bodies"]:
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


def categorical_confusion(
    adata: ad.AnnData,
    *,
    key1: str,
    key2: str,
    normalize: Literal["none", "row", "col", "all"] = "row",
    cmap: str = "Blues",
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Confusion-style heatmap for two categorical obs columns + Cramér’s V / ARI / NMI.
    """
    set_style()

    res = categorical_association(adata, key1=key1, key2=key2, metrics=("chi2", "cramers_v", "ari", "nmi"))
    tab: pd.DataFrame = res["table"].copy()

    mat = tab.to_numpy(dtype=float)

    if normalize == "row":
        mat = mat / np.maximum(mat.sum(axis=1, keepdims=True), 1.0)
    elif normalize == "col":
        mat = mat / np.maximum(mat.sum(axis=0, keepdims=True), 1.0)
    elif normalize == "all":
        mat = mat / np.maximum(mat.sum(), 1.0)

    if figsize is None:
        figsize = (max(6, 0.25 * tab.shape[1] + 3), max(4, 0.25 * tab.shape[0] + 2))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1 if normalize != "none" else None)

    ax.set_xticks(np.arange(tab.shape[1]))
    ax.set_xticklabels(tab.columns.tolist(), rotation=90)
    ax.set_yticks(np.arange(tab.shape[0]))
    ax.set_yticklabels(tab.index.tolist())

    crv = res.get("cramers_v", np.nan)
    ari = res.get("ari", np.nan)
    nmi = res.get("nmi", np.nan)

    ttl = title or f"{key1} vs {key2}"
    ttl += f" | Cramér’s V={crv:.3f}"
    if np.isfinite(ari):
        ttl += f" | ARI={ari:.3f}"
    if np.isfinite(nmi):
        ttl += f" | NMI={nmi:.3f}"
    ax.set_title(ttl)

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label("fraction" if normalize != "none" else "count")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax