from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import matplotlib.pyplot as plt
import anndata as ad

from ._style import set_style, _savefig


def qc_by_group(
    adata: ad.AnnData,
    *,
    groupby: str,
    keys: Sequence[str] = ("total_counts", "n_genes_detected", "pct_counts_mt", "pct_counts_ribo"),
    kind: Literal["violin", "box"] = "violin",
    log1p: Sequence[str] = ("total_counts",),
    figsize: tuple[float, float] = (11, 4),
    rotate_xticks: int = 45,
    save: str | Path | None = None,
    show: bool = True,
    show_n: bool = True,
):
    """
    Plot QC metrics grouped by a metadata column in `adata.obs` (e.g. batch, cohort).
    """
    set_style()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    for k in keys:
        if k not in adata.obs.columns:
            raise KeyError(f"Missing '{k}' in adata.obs. Run bk.pp.qc_metrics(adata) first.")

    groups = adata.obs[groupby].astype("category")
    cat = groups.cat.categories.tolist()

    counts = groups.value_counts().reindex(cat)

    fig, axes = plt.subplots(1, len(keys), figsize=figsize, constrained_layout=True)
    if len(keys) == 1:
        axes = [axes]

    for ax, k in zip(axes, keys):
        vals = adata.obs[k].to_numpy(dtype=float)
        if k in log1p:
            vals = np.log1p(vals)

        data = [vals[groups == g] for g in cat]

        if kind == "violin":
            parts = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
            # keep default coloring; just clean edges
            for pc in parts.get("bodies", []):
                pc.set_alpha(0.8)
        else:
            ax.boxplot(data, showfliers=False)

        ax.set_title(f"{_label(k, log1p)}")

        ax.set_xticks(range(1, len(cat) + 1))

        if show_n:
            labels = [f"{c} (n={counts[c]})" for c in cat]
        else:
            labels = cat

        ax.set_xticklabels(labels, rotation=rotate_xticks, ha="right")


        ax.set_ylabel("value")

    fig.suptitle(f"QC by {groupby}", y=1.05)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes


def _label(k: str, log1p: Sequence[str]) -> str:
    return f"log1p({k})" if k in log1p else k