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



def gene_panel_correlation_heatmap(
    adata,
    *,
    genes: Sequence[str],
    layer: str | None = "log1p_cpm",
    method: Literal["pearson", "spearman"] = "pearson",
    use_abs: bool = False,
    figsize: tuple[float, float] | None = None,
    cmap: str = "vlag",
    vmin: float | None = -1.0,
    vmax: float | None = 1.0,
    annotate: bool = False,
    annot_fontsize: float = 7.0,
    tick_fontsize: float = 8.0,
    title: str | None = None,
    cluster: bool = False,          # uses seaborn.clustermap if True
    dendrogram_ratio: float = 0.12, # only for cluster=True
    cbar_pos: tuple[float, float, float, float] = (0.92, 0.15, 0.02, 0.7),
    save: str | None = None,
    show: bool = True,
):
    """
    Pairwise correlation heatmap for a gene panel.

    Returns
    -------
    dfC : pd.DataFrame (gene x gene correlation)
    fig, ax_or_grid : matplotlib figure + axis (or seaborn ClusterGrid if cluster=True)
    """
    # --- expression matrix ---
    genes = [str(g) for g in genes if str(g) in adata.var_names]
    if len(genes) < 2:
        raise ValueError("Need at least 2 genes present in adata.var_names.")

    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    gidx = adata.var_names.get_indexer(genes)
    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # --- correlation ---
    if method == "spearman":
        Mr = np.apply_along_axis(lambda x: pd.Series(x).rank(method="average").to_numpy(), 0, M)
        C = np.corrcoef(Mr, rowvar=False)
    elif method == "pearson":
        C = np.corrcoef(M, rowvar=False)
    else:
        raise ValueError("method must be 'pearson' or 'spearman'.")

    C = np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)
    if use_abs:
        C = np.abs(C)
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = 1.0

    dfC = pd.DataFrame(C, index=genes, columns=genes)

    # --- size ---
    if figsize is None:
        s = max(4.5, 0.28 * len(genes) + 2.0)
        figsize = (s, s)


    # --- clustered heatmap (optional) ---
    if cluster:
        try:
            import seaborn as sns
            import matplotlib.pyplot as plt
        except Exception as e:
            raise ImportError(f"cluster=True requires seaborn/matplotlib. ({e})")

        g = sns.clustermap(
            dfC,
            cmap=cmap,
            row_cluster=True,
            col_cluster=True,
            xticklabels=True,
            yticklabels=True,
            figsize=figsize,
            dendrogram_ratio=float(dendrogram_ratio),
            cbar_pos=cbar_pos,  # initial hint
        )

        # --- title ---
        if title is None:
            title = f"Gene–gene correlation ({method}{', abs' if use_abs else ''})"
        g.ax_heatmap.set_title(title, pad=12)

        # --- ticks ---
        for lab in g.ax_heatmap.get_xticklabels():
            lab.set_rotation(90)
            lab.set_fontsize(float(tick_fontsize))
        for lab in g.ax_heatmap.get_yticklabels():
            lab.set_fontsize(float(tick_fontsize))

        # --- spacing (leave room for labels + the cbar) ---
        # (do this BEFORE final cbar placement)
        g.fig.subplots_adjust(bottom=0.22)

        # IMPORTANT in seaborn 0.11.x:
        # Force layout to settle, THEN force the colorbar position at the end.
        g.fig.canvas.draw()
        g.cax.set_position(cbar_pos)
        g.fig.canvas.draw()

        # if you want df in plotted order:
        row_order = [dfC.index[i] for i in g.dendrogram_row.reordered_ind]
        col_order = [dfC.columns[i] for i in g.dendrogram_col.reordered_ind]
        dfC_plot_order = dfC.loc[row_order, col_order]

        if save is not None:
            # NOTE: bbox_inches="tight" can override manual axis positions.
            # Use bbox_inches=None to preserve exact cbar_pos.
            g.fig.savefig(save, dpi=300)  # <-- no bbox_inches="tight"
        if show:
            plt.show()

        return dfC_plot_order, g.fig, g

    # --- plain matplotlib heatmap ---
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(dfC.to_numpy(), aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)

    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=float(tick_fontsize))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=float(tick_fontsize))

    if title is None:
        title = f"Gene–gene correlation ({method}{', abs' if use_abs else ''})"
    ax.set_title(title, pad=12)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f"{'| ' if use_abs else ''}{method} corr")

    if annotate and len(genes) <= 35:
        arr = dfC.to_numpy()
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                ax.text(j, i, f"{arr[i, j]:.2f}", ha="center", va="center", fontsize=float(annot_fontsize))

    fig.tight_layout()

    if save is not None:
        fig.savefig(save, bbox_inches="tight", dpi=300)
    if show:
        plt.show()

    return dfC, fig, ax