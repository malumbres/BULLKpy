from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata as ad

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None

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


def pca_loadings_heatmap(
    adata: ad.AnnData,
    *,
    pcs: Sequence[int] = (1, 2, 3),
    n_top: int = 15,
    loadings_key: str = "PCs",
    use_abs: bool = False,
    show_negative: bool = True,
    gene_symbol_key: str | None = None,
    z_score: bool = False,  # z-score per gene across PCs (like Scanpy option)
    cluster_genes: bool = True,
    cluster_pcs: bool = False,
    cmap: str = "vlag",
    center: float = 0.0,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap of PCA loadings for union of top genes across selected PCs.

    - Selects top positive and (optional) top negative genes for each PC.
    - Builds a matrix [genes x PCs] of loadings.
    - Optional clustering and z-scoring.
    """
    set_style()
    if sns is None:
        raise ImportError("pca_loadings_heatmap requires seaborn. Please install seaborn.")

    if loadings_key not in adata.varm:
        raise KeyError(f"adata.varm['{loadings_key}'] not found. Run bk.tl.pca first.")

    PCs = np.asarray(adata.varm[loadings_key], dtype=float)  # (n_vars x n_comps)
    n_vars, n_comps = PCs.shape

    pcs = [int(p) for p in pcs]
    for p in pcs:
        if p < 1 or p > n_comps:
            raise ValueError(f"PC {p} not available; choose within 1..{n_comps}")

    if gene_symbol_key is not None and gene_symbol_key in adata.var.columns:
        genes_all = adata.var[gene_symbol_key].astype(str).to_numpy()
    else:
        genes_all = adata.var_names.astype(str).to_numpy()

    # Build union gene set
    selected = []
    for p in pcs:
        v = PCs[:, p - 1]
        ok = np.isfinite(v)
        df = pd.DataFrame({"gene": genes_all[ok], "loading": v[ok]})

        if use_abs:
            df = df.assign(abs_loading=df["loading"].abs()).sort_values("abs_loading", ascending=False)
            selected.extend(df["gene"].head(int(n_top)).tolist())
        else:
            pos = df[df["loading"] > 0].sort_values("loading", ascending=False).head(int(n_top))
            selected.extend(pos["gene"].tolist())
            if show_negative:
                neg = df[df["loading"] < 0].sort_values("loading", ascending=True).head(int(n_top))
                selected.extend(neg["gene"].tolist())

    # unique while preserving order
    seen = set()
    genes = [g for g in selected if not (g in seen or seen.add(g))]

    if len(genes) == 0:
        raise ValueError("No genes selected for heatmap (all NaN or filtered).")

    # Map gene symbols back to indices (if using gene_symbol_key, may have duplicates)
    # We choose first occurrence for each symbol.
    symbol_to_idx = {}
    for i, g in enumerate(genes_all):
        if g not in symbol_to_idx:
            symbol_to_idx[g] = i
    idx = [symbol_to_idx[g] for g in genes if g in symbol_to_idx]

    mat = PCs[idx, :][:, [p - 1 for p in pcs]]  # genes x pcs
    dfm = pd.DataFrame(mat, index=[genes_all[i] for i in idx], columns=[f"PC{p}" for p in pcs])

    if z_score:
        mu = dfm.mean(axis=1)
        sd = dfm.std(axis=1, ddof=0).replace(0, 1.0)
        dfm = dfm.sub(mu, axis=0).div(sd, axis=0)

    if figsize is None:
        w = max(4.8, 0.7 * len(pcs) + 3.0)
        h = max(4.0, 0.22 * dfm.shape[0] + 2.0)
        figsize = (w, h)

    # clustermap gives Scanpy-like dendrograms; if not clustering, use heatmap
    if cluster_genes or cluster_pcs:
        cg = sns.clustermap(
            dfm,
            cmap=cmap,
            center=center,
            row_cluster=bool(cluster_genes),
            col_cluster=bool(cluster_pcs),
            yticklabels=True,
            xticklabels=True,
            figsize=figsize,
            cbar_kws={"label": "loading" + (" (z)" if z_score else "")},
        )
        if title is None:
            title = "PCA loadings heatmap"
        cg.ax_heatmap.set_title(title, pad=10)

        fig = cg.fig
    else:
        fig, ax = plt.subplots(figsize=figsize)
        sns.heatmap(
            dfm,
            cmap=cmap,
            center=center,
            yticklabels=True,
            xticklabels=True,
            cbar_kws={"label": "loading" + (" (z)" if z_score else "")},
            ax=ax,
        )
        if title is None:
            title = "PCA loadings heatmap"
        ax.set_title(title)
        plt.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig