from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata as ad

from ._style import set_style, _savefig
from ..get.rank_genes_groups_df_all import rank_genes_groups_df_all


def rank_genes_groups(
    adata: ad.AnnData,
    *,
    key: str = "rank_genes_groups",
    groups: Sequence[str] | None = None,
    n_genes: int = 10,
    sort_by: Literal["scores", "logfoldchanges", "pvals_adj", "pvals"] = "scores",
    show: bool = True,
    save: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Scanpy-like quick view: show top genes per group as a compact table plot.
    """
    set_style()

    df = rank_genes_groups_df_all(adata, key=key, groups=groups, sort_by=sort_by)
    if df.shape[0] == 0:
        raise ValueError("No rank_genes_groups results found to plot.")

    # keep top per group
    df = df.groupby("group", group_keys=False).head(int(n_genes)).copy()

    # format table text
    df["txt"] = (
        df["gene"].astype(str)
        + "  "
        + df["log2FC"].map(lambda x: f"{x:+.2f}")
        + "  "
        + df["qval"].map(lambda x: f"q={x:.2g}")
    )

    groups_order = list(pd.unique(df["group"]))

    # build a rectangular table (rows = rank, cols = groups)
    max_r = int(n_genes)
    mat = []
    for r in range(max_r):
        row = []
        for g in groups_order:
            sub = df[df["group"] == g].iloc[r : r + 1]
            row.append(sub["txt"].values[0] if len(sub) else "")
        mat.append(row)

    col_labels = groups_order
    row_labels = [f"{i+1}" for i in range(max_r)]

    if figsize is None:
        figsize = (max(6.5, 1.6 * len(col_labels) + 2.5), max(3.0, 0.35 * max_r + 1.8))

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    tbl = ax.table(
        cellText=mat,
        colLabels=col_labels,
        rowLabels=row_labels,
        loc="center",
        cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(max(8, plt.rcParams.get("font.size", 12) - 2))
    tbl.scale(1.0, 1.15)

    if title is None:
        title = f"Top {n_genes} ranked genes per group"
    ax.set_title(title, pad=12)

    plt.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax