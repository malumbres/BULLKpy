from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

from ._style import set_style, _savefig, _apply_clustergrid_style

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None

import anndata as ad

from ._style import set_style, _savefig


def _get_matrix(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def corr_heatmap(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    method: Literal["pearson", "spearman"] = "pearson",
    use: Literal["samples", "genes"] = "samples",
    groupby: str | None = None,
    groups: Sequence[str] | None = None,
    col_colors: str | Sequence[str] | None = None,
    cmap: str = "vlag",
    center: float = 0.0,
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] | None = None,
    show_labels: bool = False,
    dendrogram: bool = True,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Correlation heatmap for sample QC (or gene-gene if use="genes").

    Parameters
    ----------
    use
        "samples" -> correlation between samples (obs x obs) [default, QC]
        "genes"   -> correlation between genes (var x var)
    groupby / groups
        If provided, subset and order samples by adata.obs[groupby].
    col_colors
        obs key(s) used to annotate columns (and rows, since it's symmetric) when plotting samples.
        Values are mapped to colors automatically.
    dendrogram
        If True, uses seaborn clustermap. If False, uses heatmap without clustering.
    """
    set_style()
    if sns is None:
        raise ImportError("corr_heatmap requires seaborn. Please install seaborn.")

    X = _get_matrix(adata, layer)

    # --- choose axis for correlation ---
    if use == "samples":
        data = X  # samples x genes
        names = adata.obs_names.astype(str).tolist()
        axis_name = "samples"
    else:
        data = X.T  # genes x samples
        names = adata.var_names.astype(str).tolist()
        axis_name = "genes"

    # --- optional subsetting/ordering by groupby (samples only) ---
    order = np.arange(data.shape[0])

    if use == "samples" and groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not found in adata.obs")

        g = adata.obs[groupby].astype(str)

        mask = np.ones(adata.n_obs, dtype=bool)
        if groups is not None:
            groups = [str(x) for x in groups]
            mask = g.isin(groups).to_numpy()

        data = data[mask, :]
        g = g[mask]
        names = adata.obs_names[mask].astype(str).tolist()

        # group order: use provided groups, else alphabetical
        if groups is None:
            cat_order = sorted(pd.unique(g))
        else:
            cat_order = [x for x in groups if x in set(g)]

        order = np.argsort(pd.Categorical(g, categories=cat_order, ordered=True))
        data = data[order, :]
        names = [names[i] for i in order]
        g = g.iloc[order]
    else:
        g = None

    # --- correlation ---
    if method == "spearman":
        df = pd.DataFrame(data, index=names)
        corr = df.T.corr(method="spearman")
    else:
        # pearson fast path
        corr = np.corrcoef(data)
        corr = pd.DataFrame(corr, index=names, columns=names)

    # --- col_colors mapping for samples ---
    col_colors_df = None
    if use == "samples" and col_colors is not None:
        if isinstance(col_colors, str):
            col_colors = [col_colors]

        ann = {}
        # rebuild obs aligned to corr index
        obs_sub = adata.obs.loc[corr.index]

        for key in col_colors:
            if key not in obs_sub.columns:
                raise KeyError(f"col_colors obs key '{key}' not found in adata.obs")

            vals = obs_sub[key].astype(str)
            cats = pd.Categorical(vals).categories
            pal = sns.color_palette("tab20", n_colors=len(cats))
            lut = {cat: pal[i] for i, cat in enumerate(cats)}
            ann[key] = vals.map(lut)

        col_colors_df = pd.DataFrame(ann, index=corr.index)

    # --- autosize ---
    if figsize is None:
        n = corr.shape[0]
        # clamp: big matrices become huge otherwise
        w = min(max(6.0, 0.18 * n + 2.0), 18.0)
        h = min(max(6.0, 0.18 * n + 2.0), 18.0)
        figsize = (w, h)

    # --- plot ---
    if dendrogram:
        cg = sns.clustermap(
            corr,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            row_cluster=True,
            col_cluster=True,
            xticklabels=show_labels,
            yticklabels=show_labels,
            figsize=figsize,
            col_colors=col_colors_df,
            cbar_kws={"label": f"{method} correlation ({axis_name})"},
        )

        _apply_clustergrid_style(cg)
        cg.fig.tight_layout()   
     
        fig = cg.fig

    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            corr,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            square=True,
            xticklabels=show_labels,
            yticklabels=show_labels,
            cbar_kws={"label": f"{method} correlation ({axis_name})"},
            ax=ax,
        )
        ax.set_title(f"{method} correlation ({axis_name})")

    _apply_clustergrid_style(cg)
    cg.fig.tight_layout()

    if save is not None:
        _savefig(cg.fig, save)
    if show:
        plt.show()

    return cg