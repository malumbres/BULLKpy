from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad

from ._style import set_style, _savefig


def corr_heatmap(
    adata: ad.AnnData,
    *,
    layer: str = "log1p_cpm",
    method: str = "pearson",     # "pearson" or "spearman"
    groupby: str | None = None,  # annotate by group
    figsize: tuple[float, float] = (7, 6),
    cmap: str = "vlag",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Sample-sample correlation heatmap for QC.
    """
    set_style()
    X = adata.layers[layer] if (layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X, dtype=float)

    # correlation across genes -> samples x samples
    if method == "pearson":
        C = np.corrcoef(X, rowvar=True)
    elif method == "spearman":
        # rank transform per gene dimension
        from scipy.stats import rankdata
        R = np.apply_along_axis(rankdata, 0, X)  # rank each gene across samples
        C = np.corrcoef(R, rowvar=True)
    else:
        raise ValueError("method must be 'pearson' or 'spearman'")

    C = np.asarray(C, dtype=float)
    mat = pd.DataFrame(C, index=adata.obs_names, columns=adata.obs_names)

    # annotation
    col_colors = None
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not in adata.obs")
        grp = adata.obs[groupby].astype("category")
        palette = sns.color_palette("Set2", n_colors=len(grp.cat.categories))
        lut = {cat: palette[i] for i, cat in enumerate(grp.cat.categories)}
        col_colors = grp.map(lut)

    g = sns.clustermap(
        mat,
        cmap=cmap,
        center=0,
        row_colors=col_colors,
        col_colors=col_colors,
        figsize=figsize,
        xticklabels=False,
        yticklabels=False,
    )
    g.fig.suptitle(f"Sample correlation ({method})", y=1.02)

    if save is not None:
        _savefig(g.fig, save)
    if show:
        plt.show()
    return g