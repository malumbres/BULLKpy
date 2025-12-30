from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import matplotlib.pyplot as plt
import anndata as ad

from ..logging import warn
from ._style import set_style, _savefig


def qc_pairplot(
    adata: ad.AnnData,
    *,
    keys: Sequence[str] = ("total_counts", "n_genes_detected", "pct_counts_mt"),
    color: str | None = "pct_counts_mt",
    log1p: Sequence[str] = ("total_counts",),
    point_size: float = 14.0,
    alpha: float = 0.7,
    figsize: tuple[float, float] = (8, 8),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scatter-matrix (pairplot) of QC metrics in `adata.obs`.

    Diagonal: histograms
    Off-diagonal: scatter plots
    """
    set_style()

    for k in keys:
        if k not in adata.obs.columns:
            raise KeyError(f"Missing '{k}' in adata.obs. Run bk.pp.qc_metrics(adata) first.")

    C = None
    if color is not None:
        if color in adata.obs.columns:
            C = adata.obs[color].to_numpy()
        else:
            warn(f"color='{color}' not found in adata.obs; coloring disabled.")

    # prepare transformed columns
    data = {}
    for k in keys:
        v = adata.obs[k].to_numpy(dtype=float)
        if k in log1p:
            v = np.log1p(v)
        data[k] = v

    n = len(keys)
    fig, axes = plt.subplots(n, n, figsize=figsize, constrained_layout=True)

    for i, yi in enumerate(keys):
        for j, xj in enumerate(keys):
            ax = axes[i, j]
            if i == j:
                ax.hist(data[xj], bins=30)
                ax.set_ylabel("")
            else:
                sc = ax.scatter(
                    data[xj], data[yi],
                    c=C, s=point_size, alpha=alpha, edgecolors="none"
                )
            # labels only on left and bottom
            if i < n - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(_label(xj, log1p))
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(_label(yi, log1p))

    if C is not None:
        cbar = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.01)
        cbar.set_label(color)

    fig.suptitle("QC pairplot", y=1.02)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes


def _label(k: str, log1p: Sequence[str]) -> str:
    return f"log1p({k})" if k in log1p else k