from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import anndata as ad

from ._style import set_style, _savefig
from ..logging import warn


def _get_expr_matrix(
    adata: ad.AnnData,
    genes: Sequence[str],
    *,
    layer: str | None = "log1p_cpm",
) -> pd.DataFrame:
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not in adata.var_names (first 10): {missing[:10]}")

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    gidx = [adata.var_names.get_loc(g) for g in genes]
    M = X[:, gidx].toarray() if sp.issparse(X) else np.asarray(X[:, gidx], dtype=float)
    return pd.DataFrame(M, index=adata.obs_names.astype(str), columns=[str(g) for g in genes])


def _binary_from_obs(adata: ad.AnnData, cols: Sequence[str]) -> pd.DataFrame:
    miss = [c for c in cols if c not in adata.obs.columns]
    if miss:
        raise KeyError(f"Mutation columns not found in adata.obs (first 10): {miss[:10]}")

    df = adata.obs[list(cols)].copy()
    out = pd.DataFrame(index=df.index.astype(str))
    for c in cols:
        s = df[c]
        if str(s.dtype) == "bool":
            v = s.astype(int)
        else:
            x = pd.to_numeric(s, errors="coerce")
            if x.isna().mean() > 0.2:
                v = (
                    s.astype(str)
                    .str.lower()
                    .map({"1": 1, "0": 0, "true": 1, "false": 0})
                    .fillna(0)
                    .astype(int)
                )
            else:
                v = x.fillna(0).astype(int)
        out[c] = v.clip(lower=0, upper=1).astype(int)
    return out


def _sort_genes_by_freq(mut_samp_x_gene: pd.DataFrame) -> list[str]:
    freq = mut_samp_x_gene.mean(axis=0).sort_values(ascending=False)
    return freq.index.astype(str).tolist()


def _sort_samples_mutation_first(mut_samp_x_gene: pd.DataFrame, gene_order: Sequence[str]) -> pd.Index:
    """
    cBioPortal-like:
      Sort samples so that mutants in gene_order[0] are left,
      then within those, mutants in gene_order[1] are left, etc.
    Equivalent to lex sort on the mutation bit-vectors (mutants first).
    """
    M = mut_samp_x_gene.loc[:, list(gene_order)].to_numpy(dtype=int)  # (n_samples x n_genes)
    burden = M.sum(axis=1)

    # lexsort: last key is primary, so build keys reversed.
    # We want: gene0 mut first, then gene1 mut first, ... and higher burden first.
    keys = [-(burden.astype(int))]
    for j in range(M.shape[1] - 1, -1, -1):
        keys.append(-M[:, j])
    order = np.lexsort(keys)  # stable, deterministic
    return mut_samp_x_gene.index[order]


def _apply_group_contiguity(
    mut: pd.DataFrame,
    groups: pd.Series,
    *,
    group_order: Sequence[str] | None = None,
    within_group_sort: Literal["mut_first", "burden", "none"] = "mut_first",
) -> tuple[pd.DataFrame, pd.Series]:
    g = groups.loc[mut.index].astype(str)

    if group_order is not None:
        cats = [str(x) for x in group_order]
        gcat = pd.Categorical(g, categories=cats, ordered=True)
    else:
        gcat = pd.Categorical(g)

    blocks = []
    g_out = []
    for level in list(gcat.categories):
        idx = mut.index[gcat == level]
        if len(idx) == 0:
            continue
        sub = mut.loc[idx]
        if within_group_sort == "mut_first":
            gene_order = list(sub.columns)
            idx2 = _sort_samples_mutation_first(sub, gene_order)
            sub = sub.loc[idx2]
        elif within_group_sort == "burden":
            burden = sub.sum(axis=1).sort_values(ascending=False)
            sub = sub.loc[burden.index]
        elif within_group_sort == "none":
            pass
        else:
            raise ValueError("within_group_sort must be 'mut_first', 'burden', or 'none'.")
        blocks.append(sub)
        g_out.append(pd.Series([str(level)] * sub.shape[0], index=sub.index))

    mut2 = pd.concat(blocks, axis=0) if blocks else mut
    g2 = pd.concat(g_out) if g_out else g
    return mut2, g2


def oncoprint(
    adata: ad.AnnData,
    *,
    mut_cols: Sequence[str],
    expr_genes: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    # ordering
    sort_genes: bool = True,
    sort_samples: Literal["mut_first", "burden", "none"] = "mut_first",
    drop_all_wt: bool = True,
    max_samples: int | None = None,
    # group blocks
    groupby: str | None = None,
    group_order: Sequence[str] | None = None,
    group_blocks: bool = True,
    within_group_sort: Literal["mut_first", "burden", "none"] = "mut_first",
    # styling
    show_sample_labels: bool = False,
    mut_color: str = "#222222",
    wt_color: str = "#FFFFFF",
    grid_color: str = "0.85",
    expr_cmap: str = "viridis",
    expr_vmin: float | None = None,
    expr_vmax: float | None = None,
    expr_zscore: bool = False,
    cell_size: float | None = None,     # inches per sample (auto if None)
    row_height: float = 0.35,
    expr_row_height: float = 0.25,
    top_annotation_height: float = 0.35,
    title: str | None = None,
    # output safety
    save_dpi: int | None = None,        # auto-reduced if needed
    max_pixels: int = 60000,            # keep < 2^16 per dimension
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Binary oncoprint (mut vs wt) from adata.obs 0/1 columns.
    Implements mutation-first sample ordering and optional group blocks.
    """
    set_style()

    mut = _binary_from_obs(adata, mut_cols)  # (samples x genes)
    mut.index = mut.index.astype(str)

    # drop all-wt samples (C)
    if drop_all_wt:
        keep = mut.sum(axis=1) > 0
        mut = mut.loc[keep]
        if mut.shape[0] == 0:
            raise ValueError("After drop_all_wt=True, no samples have any mutations in mut_cols.")

    # gene order (B)
    if sort_genes:
        gene_order = _sort_genes_by_freq(mut)
        mut = mut[gene_order]
    else:
        gene_order = list(mut.columns)

    # sample cap
    if max_samples is not None and mut.shape[0] > int(max_samples):
        mut = mut.iloc[: int(max_samples)].copy()
        warn(f"Truncated to first max_samples={max_samples} samples for plotting.")

    # group info
    groups = None
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not found in adata.obs")
        groups = adata.obs.loc[mut.index, groupby].astype(str)

    # sample ordering (B) + group blocks (A)
    if groups is not None and group_blocks:
        # keep groups contiguous, sort inside each block
        mut, groups = _apply_group_contiguity(
            mut,
            groups,
            group_order=group_order,
            within_group_sort=within_group_sort,
        )
    else:
        if sort_samples == "none":
            pass
        elif sort_samples == "burden":
            burden = mut.sum(axis=1).sort_values(ascending=False)
            mut = mut.loc[burden.index]
            if groups is not None:
                groups = groups.loc[mut.index]
        elif sort_samples == "mut_first":
            idx = _sort_samples_mutation_first(mut, gene_order=gene_order)
            mut = mut.loc[idx]
            if groups is not None:
                groups = groups.loc[mut.index]
        else:
            raise ValueError("sort_samples must be 'mut_first', 'burden', or 'none'.")

    # expression tracks (optional)
    expr_df = None
    if expr_genes is not None and len(expr_genes) > 0:
        expr_df = _get_expr_matrix(adata, expr_genes, layer=layer).loc[mut.index]
        if expr_zscore:
            arr = expr_df.to_numpy(dtype=float)
            mu = np.nanmean(arr, axis=0, keepdims=True)
            sd = np.nanstd(arr, axis=0, keepdims=True)
            sd[sd == 0] = 1.0
            expr_df = pd.DataFrame((arr - mu) / sd, index=expr_df.index, columns=expr_df.columns)

    n_samples = mut.shape[0]
    n_genes = mut.shape[1]
    n_expr = 0 if expr_df is None else expr_df.shape[1]

    # auto cell size to avoid gigantic images (D)
    if cell_size is None:
        # target a manageable width in inches even for big n
        # (still readable, but prevents exploding pixel sizes)
        cell_size = float(np.clip(18.0 / max(n_samples, 1), 0.02, 0.22))

    # figure size
    width = max(6.0, cell_size * n_samples + 2.8)
    height = (
        (top_annotation_height if groups is not None else 0.0)
        + row_height * n_genes
        + (0.35 if n_expr > 0 else 0.0)
        + expr_row_height * n_expr
        + 1.2
    )

    fig = plt.figure(figsize=(width, height), constrained_layout=False)

    heights = []
    if groups is not None:
        heights.append(top_annotation_height)
    heights.append(row_height * n_genes)
    if n_expr > 0:
        heights.append(expr_row_height * n_expr + 0.35)

    gs = GridSpec(
        nrows=len(heights),
        ncols=2,
        figure=fig,
        width_ratios=[1.0, 0.10],
        height_ratios=heights,
        hspace=0.15,
        wspace=0.08,
    )

    row_i = 0

    # --- top group strip + separators (A) ---
    ax_top = None
    group_boundaries = []
    if groups is not None:
        ax_top = fig.add_subplot(gs[row_i, 0])
        row_i += 1

        cats = list(pd.Categorical(groups).categories)
        cmap = mpl.cm.get_cmap("tab20" if len(cats) <= 20 else "hsv")
        cols = [mpl.colors.to_hex(cmap(i / max(1, len(cats) - 1))) for i in range(len(cats))]
        m = {str(c): cols[i] for i, c in enumerate(cats)}

        g_list = groups.astype(str).tolist()
        for i, g in enumerate(g_list):
            ax_top.add_patch(mpl.patches.Rectangle((i, 0), 1, 1, facecolor=m.get(g, "#cccccc"), edgecolor="none"))

        # group boundaries for separators
        for i in range(1, len(g_list)):
            if g_list[i] != g_list[i - 1]:
                group_boundaries.append(i)

        ax_top.set_xlim(0, n_samples)
        ax_top.set_ylim(0, 1)
        ax_top.set_xticks([])
        ax_top.set_yticks([])
        for spn in ax_top.spines.values():
            spn.set_visible(False)

        # legend
        handles = [mpl.patches.Patch(color=m[str(c)], label=str(c)) for c in cats]
        ax_top.legend(
            handles=handles,
            title=groupby,
            bbox_to_anchor=(1.22, 0.0),  # moved right so it doesn't overlap freq bar
            loc="lower left",
            frameon=False,
            borderaxespad=0.0,
    )

    # --- mutation panel ---
    ax = fig.add_subplot(gs[row_i, 0])
    row_i += 1

    ax.set_facecolor(wt_color)

    A = mut.to_numpy(dtype=int).T  # genes x samples
    for gi in range(n_genes):
        y = n_genes - 1 - gi  # top gene at top
        muts = np.where(A[gi, :] == 1)[0]
        if muts.size:
            # rasterize rectangles for speed on large cohorts
            for x in muts:
                ax.add_patch(
                    mpl.patches.Rectangle((x, y), 1, 1, facecolor=mut_color, edgecolor="none", rasterized=True)
                )

    # grid
    ax.set_xlim(0, n_samples)
    ax.set_ylim(0, n_genes)
    ax.set_xticks(np.arange(0, n_samples + 1, 1), minor=True)
    ax.set_yticks(np.arange(0, n_genes + 1, 1), minor=True)
    ax.grid(which="minor", color=grid_color, linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    # separators between group blocks
    for b in group_boundaries:
        ax.axvline(b, color="0.5", lw=1.2, alpha=0.8)

    # labels
    ax.set_yticks(np.arange(n_genes) + 0.5)
    ax.set_yticklabels(list(reversed(mut.columns.astype(str).tolist())))
    ax.set_xticks(np.arange(n_samples) + 0.5)
    ax.set_xticklabels(mut.index.astype(str).tolist(), rotation=90, fontsize=7 if show_sample_labels else 7)
    if not show_sample_labels:
        ax.set_xticklabels([])
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)

    ax.set_xlabel("Samples")
    ax.set_ylabel("Mutations")
    ax.set_title(title or "Oncoprint", pad=8)

    # side bar: frequency
    ax_bar = fig.add_subplot(gs[(0 if ax_top is None else 1), 1])
    freq = mut.mean(axis=0).to_numpy(dtype=float)[::-1]
    ax_bar.barh(np.arange(n_genes) + 0.5, freq, height=0.85)
    ax_bar.set_ylim(0, n_genes)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Freq")
    ax_bar.xaxis.set_label_position("top")
    ax_bar.xaxis.tick_top()
    for spn in ["right", "bottom", "left"]:
        ax_bar.spines[spn].set_visible(False)

    # expression panel
    if expr_df is not None and n_expr > 0:
        ax_expr = fig.add_subplot(gs[row_i, 0])

        X = expr_df.to_numpy(dtype=float).T  # genes x samples
        vmin = float(np.nanmin(X)) if expr_vmin is None else float(expr_vmin)
        vmax = float(np.nanmax(X)) if expr_vmax is None else float(expr_vmax)

        im = ax_expr.imshow(
            X,
            aspect="auto",
            interpolation="nearest",
            cmap=expr_cmap,
            vmin=vmin,
            vmax=vmax,
        )
        for b in group_boundaries:
            ax_expr.axvline(b - 0.5, color="0.5", lw=1.2, alpha=0.8)

        ax_expr.set_yticks(np.arange(n_expr))
        ax_expr.set_yticklabels(expr_df.columns.astype(str).tolist())
        ax_expr.set_xticks([])
        ax_expr.tick_params(axis="y", length=0)

        # compact colorbar
        cax = fig.add_axes([0.92, 0.12, 0.015, 0.16])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(("Z-score" if expr_zscore else "Expression"), rotation=90)

    fig.subplots_adjust(left=0.26, right=0.86, top=0.92, bottom=0.10)

    # ---- SAFE SAVE (D): avoid matplotlib 2^16 px limit ----
    if save is not None:
        # estimate current dpi
        current_dpi = int(mpl.rcParams.get("savefig.dpi", mpl.rcParams.get("figure.dpi", 150)))
        dpi_use = current_dpi if save_dpi is None else int(save_dpi)

        # enforce pixel limits
        px_w = int(width * dpi_use)
        px_h = int(height * dpi_use)
        if px_w > int(max_pixels) or px_h > int(max_pixels):
            dpi_w = int(max(20, max_pixels // max(width, 1e-6)))
            dpi_h = int(max(20, max_pixels // max(height, 1e-6)))
            dpi_use = int(max(20, min(dpi_use, dpi_w, dpi_h)))
            warn(f"Reducing save DPI to {dpi_use} to avoid huge image ({px_w}x{px_h}px).")

        # save directly with dpi override
        save = Path(save)
        save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save, dpi=dpi_use, bbox_inches="tight")

    if show:
        plt.show()

    return fig, ax