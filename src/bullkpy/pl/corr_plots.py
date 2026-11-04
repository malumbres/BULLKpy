from __future__ import annotations

from typing import Literal

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..tl.correlations import BatchMode, Method, _corr_batch_aware
from .._matrix import get_gene_vector as _get_gene_vector
import seaborn as sns


def plot_corr_scatter(
    adata: ad.AnnData,
    *,
    x: str,
    y: str,
    hue: str | None = None,
    kind: Literal["obs_obs", "gene_obs"] = "obs_obs",
    layer: str | None = "log1p_cpm",
    method: Method = "pearson",
    batch_key: str | None = None,
    batch_mode: BatchMode = "none",
    covariates: list[str] | None = None,
    figsize: tuple[float, float] = (5.5, 4.5),
    s: float = 18.0,
    alpha: float = 0.8,
    cmap: str = "viridis",
    palette: str = "tab20",
    title: str | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Scatter plot + correlation in title.
    kind="obs_obs": x,y are obs numeric columns
    kind="gene_obs": x is gene, y is obs numeric (or vice versa)
    """
    if kind == "obs_obs":
        if x not in adata.obs.columns or y not in adata.obs.columns:
            raise KeyError("x and y must be obs columns for kind='obs_obs'")
        xv = pd.to_numeric(adata.obs[x], errors="coerce").to_numpy(dtype=float)
        yv = pd.to_numeric(adata.obs[y], errors="coerce").to_numpy(dtype=float)
    elif kind == "gene_obs":
        # allow x=gene, y=obs or vice versa
        if x in adata.var_names and y in adata.obs.columns:
            xv = _get_gene_vector(adata, x, layer=layer)
            yv = pd.to_numeric(adata.obs[y], errors="coerce").to_numpy(dtype=float)
        elif y in adata.var_names and x in adata.obs.columns:
            xv = pd.to_numeric(adata.obs[x], errors="coerce").to_numpy(dtype=float)
            yv = _get_gene_vector(adata, y, layer=layer)
        else:
            raise KeyError("For kind='gene_obs', one of (x,y) must be a gene and the other a numeric obs.")
    else:
        raise ValueError(f"Unknown kind={kind}")

    r, p, n = _corr_batch_aware(
        xv, yv, adata=adata, method=method,
        batch_key=batch_key, batch_mode=batch_mode,
        covariates=covariates,
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    if hue is None:
        ax.scatter(xv, yv, s=s, alpha=alpha, edgecolors="none")
    else:
        if hue not in adata.obs.columns:
            raise KeyError(f"hue='{hue}' not in adata.obs")
        h = adata.obs[hue]
        if pd.api.types.is_numeric_dtype(h):
            vals = pd.to_numeric(h, errors="coerce").to_numpy(dtype=float)
            sc = ax.scatter(xv, yv, c=vals, s=s, alpha=alpha, cmap=cmap, edgecolors="none")
            cb = fig.colorbar(sc, ax=ax, pad=0.01)
            cb.set_label(hue)
        else:
            cat = pd.Categorical(h.astype(str))
            cats = list(cat.categories)
            if sns is not None:
                pal = sns.color_palette(palette, n_colors=len(cats))
                colors = {c: pal[i] for i, c in enumerate(cats)}
            else:
                cm = plt.get_cmap(palette)
                colors = {c: cm(i % cm.N) for i, c in enumerate(cats)}
            for c in cats:
                m = (cat == c)
                ax.scatter(xv[m], yv[m], s=s, alpha=alpha, edgecolors="none", label=str(c), color=colors[c])
            ax.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ttl = title or f"{method} r={r:.3f}, p={p:.2e}, n={n}"
    ax.set_title(ttl)
    if show:
        plt.show()
    return fig, ax
def plot_corr_heatmap(
    mat: pd.DataFrame,
    *,
    figsize: tuple[float, float] = (7.5, 4.5),
    cmap: str = "vlag",
    center: float = 0.0,
    vmin: float | None = None,
    vmax: float | None = None,
    annot: bool = False,
    fmt: str = ".2f",
    title: str | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Heatmap for a wide correlation matrix (e.g., from obs_obs_corr_matrix).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure
    if sns is not None:
        sns.heatmap(mat, ax=ax, cmap=cmap, center=center, vmin=vmin, vmax=vmax, annot=annot, fmt=fmt)
    else:
        im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto")
        fig.colorbar(im, ax=ax)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns, rotation=90)
        ax.set_yticklabels(mat.index)

    if title is not None:
        ax.set_title(title)

    if show:
        plt.show()
    return fig, ax
