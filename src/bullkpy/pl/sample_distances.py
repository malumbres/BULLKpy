from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp

import anndata as ad
from matplotlib.colors import to_hex
from ._colors import categorical_colors_array, get_categorical_colors


from ._style import set_style, _savefig

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None

try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import linkage
except Exception:  # pragma: no cover
    pdist = None
    squareform = None
    linkage = None


def _get_matrix(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    use: Literal["samples", "genes"] = "samples",
) -> np.ndarray:
    """
    Returns dense matrix for distance computation.
    - use="samples": rows=samples, cols=genes  (default; typical)
    - use="genes":   rows=genes, cols=samples  (gene-gene distances)
    """
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    if use == "genes":
        X = X.T
    return X


def _metadata_colors(adata, *, columns, palette="tab20"):
    if columns is None or len(columns) == 0:
        return None, {}

    out = []
    legend_maps = {}
    for c in columns:
        # stable mapping stored in bk.settings
        cmap = get_categorical_colors(adata, key=c, where="obs", palette=palette)
        legend_maps[c] = cmap
        out.append(categorical_colors_array(adata, key=c, where="obs", palette=palette))

    col_colors_df = pd.concat(out, axis=1)
    col_colors_df.columns = list(columns)
    col_colors_df.index = adata.obs_names
    return col_colors_df, legend_maps

def sample_distances(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    metric: str = "euclidean",
    method: str = "average",
    use: Literal["samples", "genes"] = "samples",
    col_colors: Sequence[str] | None = None,
    palette: str = "tab20",
    z_score: bool = False,
    figsize: tuple[float, float] | None = None,
    show_labels: bool = False,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Sample (or gene) distance clustergram.

    - Computes pairwise distances (pdist) on X (samples x genes).
    - Uses seaborn.clustermap with hierarchical clustering.
    - Optionally annotate samples with metadata columns via col_colors.

    Notes:
      - For sample QC, use metric="correlation" (distance = 1-corr) often works well.
      - z_score=True will z-score genes across samples before distance computation.
    """
    set_style()
    if sns is None:
        raise ImportError("sample_distances requires seaborn. Please install seaborn.")
    if pdist is None or squareform is None or linkage is None:
        raise ImportError("sample_distances requires scipy (distance + hierarchy).")

    X = _get_matrix(adata, layer=layer, use=use)

    # Optional z-score across features (genes) for sample distance
    if z_score:
        mu = X.mean(axis=0, keepdims=True)
        sd = X.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        X = (X - mu) / sd

    # Distance matrix
    d = pdist(X, metric=metric)
    D = squareform(d)

    labels = list(adata.obs_names) if use == "samples" else list(adata.var_names)
    dfD = pd.DataFrame(D, index=labels, columns=labels)

    # Clustering on distances (linkage expects condensed distances)
    Z = linkage(d, method=method)

    # Metadata colors only meaningful for samples
    col_colors_df = None
    legend_maps = {}
    if use == "samples" and col_colors is not None and len(col_colors) > 0:
        col_colors_df, legend_maps = _metadata_colors(adata, columns=col_colors, palette=palette)
        # reorder to match dfD labels
        col_colors_df = col_colors_df.loc[dfD.index]

    # autosize
    if figsize is None:
        n = len(labels)
        w = max(6.0, min(16.0, 0.18 * n + 4.0))
        h = w
        figsize = (w, h)

    cg = sns.clustermap(
        dfD,
        row_linkage=Z,
        col_linkage=Z,
        cmap="viridis",
        figsize=figsize,
        xticklabels=show_labels,
        yticklabels=show_labels,
        col_colors=col_colors_df,
        cbar_kws={"label": f"{metric} distance"},
    )

    cg.ax_heatmap.set_title("Sample distances" if use == "samples" else "Gene distances", pad=10)

    # Add metadata legends on the right
    if legend_maps:
        ax = cg.ax_heatmap
        # Place small legends outside
        x0 = 1.02
        y0 = 1.0
        dy = 0.06
        for j, (col, cmap) in enumerate(legend_maps.items()):
            y = y0 - j * (dy * (len(cmap) + 1))
            ax.text(x0, y, col, transform=ax.transAxes, ha="left", va="top", fontsize=plt.rcParams["font.size"])
            y -= dy
            for k, v in cmap.items():
                ax.scatter([x0], [y], transform=ax.transAxes, s=40, c=v, clip_on=False)
                ax.text(x0 + 0.03, y, str(k), transform=ax.transAxes, ha="left", va="center")
                y -= dy

    if save is not None:
        _savefig(cg.fig, save)
    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return cg


def sample_correlation_clustergram(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    method: Literal["pearson", "spearman"] = "pearson",
    linkage_method: str = "average",
    col_colors: Sequence[str] | None = None,
    palette: str = "tab20",
    figsize: tuple[float, float] | None = None,
    show_labels: bool = False,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Sample correlation clustergram (more interpretable than raw distances for bulk QC).
    Displays correlation, clusters by (1 - correlation).

    Heatmap values: correlation in [-1, 1]
    Clustering: on distance = 1 - correlation
    """
    set_style()
    if sns is None:
        raise ImportError("sample_correlation_clustergram requires seaborn.")
    if pdist is None or squareform is None or linkage is None:
        raise ImportError("sample_correlation_clustergram requires scipy.")

    X = _get_matrix(adata, layer=layer, use="samples")
    if method == "spearman":
        Xr = pd.DataFrame(X, index=adata.obs_names).T.rank(axis=0).to_numpy().T
        X = Xr

    # correlation matrix
    C = np.corrcoef(X)
    labels = list(adata.obs_names)
    dfC = pd.DataFrame(C, index=labels, columns=labels)

    # distance for clustering
    d = squareform(1.0 - C, checks=False)
    Z = linkage(d, method=linkage_method)

    col_colors_df = None
    legend_maps = {}
    if col_colors is not None and len(col_colors) > 0:
        col_colors_df, legend_maps = _metadata_colors(adata, columns=col_colors, palette=palette)
        col_colors_df = col_colors_df.loc[dfC.index]

    if figsize is None:
        n = len(labels)
        w = max(6.0, min(16.0, 0.18 * n + 4.0))
        figsize = (w, w)

    cg = sns.clustermap(
        dfC,
        row_linkage=Z,
        col_linkage=Z,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        center=0,
        figsize=figsize,
        xticklabels=show_labels,
        yticklabels=show_labels,
        col_colors=col_colors_df,
        cbar_kws={"label": f"{method} correlation"},
    )

    cg.ax_heatmap.set_title(f"Sample correlation ({method})", pad=10)

    if legend_maps:
        ax = cg.ax_heatmap
        x0 = 1.02
        y0 = 1.0
        dy = 0.06
        import matplotlib.pyplot as plt
        for j, (col, cmap) in enumerate(legend_maps.items()):
            y = y0 - j * (dy * (len(cmap) + 1))
            ax.text(x0, y, col, transform=ax.transAxes, ha="left", va="top", fontsize=plt.rcParams["font.size"])
            y -= dy
            for k, v in cmap.items():
                ax.scatter([x0], [y], transform=ax.transAxes, s=40, c=v, clip_on=False)
                ax.text(x0 + 0.03, y, str(k), transform=ax.transAxes, ha="left", va="center")
                y -= dy

    if save is not None:
        _savefig(cg.fig, save)
    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return cg