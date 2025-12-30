from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata as ad

from ._style import set_style, _savefig


def pca_loadings_bar(
    adata: ad.AnnData,
    *,
    pc: int = 1,
    n_top: int = 15,
    loadings_key: str = "PCs",
    use_abs: bool = False,
    show_negative: bool = True,
    gene_symbol_key: str | None = None,  # e.g. "gene_symbol" in adata.var
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Plot top PCA loadings for a single PC.

    - If use_abs=True: shows top |loading| (all positive bars).
    - Else: shows top positive and (optionally) top negative loadings.
    """
    set_style()

    if loadings_key not in adata.varm:
        raise KeyError(f"adata.varm['{loadings_key}'] not found. Run bk.tl.pca first.")

    PCs = np.asarray(adata.varm[loadings_key], dtype=float)  # (n_vars x n_comps)
    n_vars, n_comps = PCs.shape
    pc0 = int(pc) - 1
    if pc0 < 0 or pc0 >= n_comps:
        raise ValueError(f"pc must be in 1..{n_comps}, got {pc}")

    load = PCs[:, pc0]
    ok = np.isfinite(load)
    load = load[ok]

    if gene_symbol_key is not None and gene_symbol_key in adata.var.columns:
        genes_all = adata.var[gene_symbol_key].astype(str).to_numpy()
    else:
        genes_all = adata.var_names.astype(str).to_numpy()
    genes = genes_all[ok]

    df = pd.DataFrame({"gene": genes, "loading": load})

    if use_abs:
        top = (
            df.assign(abs_loading=df["loading"].abs())
            .sort_values("abs_loading", ascending=False)
            .head(int(n_top))
            .copy()
        )
        top = top.sort_values("loading")  # nice ordering in plot
        top["group"] = "abs"
        plot_df = top

    else:
        pos = df[df["loading"] > 0].sort_values("loading", ascending=False).head(int(n_top)).copy()
        neg = df[df["loading"] < 0].sort_values("loading", ascending=True).head(int(n_top)).copy()

        pos["group"] = "pos"
        neg["group"] = "neg"

        plot_df = pos
        if show_negative:
            plot_df = pd.concat([neg, pos], axis=0)

        plot_df = plot_df.sort_values("loading")

    if figsize is None:
        h = max(3.2, 0.22 * plot_df.shape[0] + 1.2)
        figsize = (6.8, h)

    fig, ax = plt.subplots(figsize=figsize)

    # color convention: negatives one tone, positives another; abs as neutral
    colors = []
    for _, r in plot_df.iterrows():
        if r["group"] == "neg":
            colors.append("0.55")
        elif r["group"] == "pos":
            colors.append("0.15")
        else:
            colors.append("0.25")

    ax.barh(plot_df["gene"], plot_df["loading"], color=colors)
    ax.axvline(0, lw=1, color="0.6")

    ax.set_xlabel("Loading")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=max(1, plt.rcParams.get("font.size", 12) - 1))

    if title is None:
        title = f"PCA loadings: PC{pc}"
    ax.set_title(title)

    plt.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax