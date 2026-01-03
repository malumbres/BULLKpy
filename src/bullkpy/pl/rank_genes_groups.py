from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import anndata as ad

from ._style import set_style, _savefig
from ..get.rank_genes_groups_df_all import rank_genes_groups_df_all
from .dotplot import dotplot

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

def rank_genes_groups_dotplot(
    adata: ad.AnnData,
    *,
    groupby: str,
    groups: Sequence[str] | None = None,
    key: str = "rank_genes_groups",
    n_genes: int = 5,
    sort_by: Literal["scores", "logfoldchanges", "pvals_adj", "pvals"] = "scores",
    unique: bool = True,
    use_abs: bool = False,

    # --- NEW 1) what to color by ---
    values_to_plot: Literal["expression", "logfoldchanges"] = "expression",

    # expression parameters (passed to dotplot)
    layer: str | None = "log1p_cpm",
    fraction_layer: str | None = "counts",
    expr_threshold: float = 0.0,

    # --- NEW 2) fraction filters for selecting genes ---
    min_in_group_fraction: float | None = None,  # e.g. 0.2
    max_in_group_fraction: float | None = None,  # e.g. 0.95

    # --- NEW 3) auto standard_scale for expression mode ---
    standard_scale: str | None = "auto",  # "auto" | "var" | "group" | None

    # dotplot layout options
    swap_axes: bool = True,
    dendrogram_top: bool = True,
    dendrogram_rows: bool = False,
    row_dendrogram_position: Literal["right", "left", "outer_left"] = "right",
    row_spacing: float = 0.75,
    cmap: str = "Reds",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scanpy-like dotplot of ranked genes per group.

    Enhancements:
      1) values_to_plot:
         - "expression"      -> dot color = mean expression per group (your dotplot default)
         - "logfoldchanges"  -> dot color = log2FC from rank_genes_groups
      2) min/max fraction filters (computed from fraction_layer + expr_threshold)
      3) standard_scale="auto" -> uses "var" for expression mode (Scanpy-like)
    """
    set_style()

    # ---- gather DE/ranking table ----
    df = rank_genes_groups_df_all(adata, key=key, groups=groups, sort_by=sort_by).copy()
    if df.shape[0] == 0:
        raise ValueError("No rank_genes_groups results found.")

    if use_abs:
        df["__rankval"] = df["scores"].abs()
        df = df.sort_values(["group", "__rankval"], ascending=[True, False])
    else:
        df = df.sort_values(["group", "scores"], ascending=[True, False])

    # ---- optionally compute fraction per (group, gene) and filter ----
    if (min_in_group_fraction is not None) or (max_in_group_fraction is not None):
        # compute fractions for ONLY genes that might be used (top-ish chunk)
        # (keeps it efficient)
        # take a little buffer so filtering doesn't leave empty sets
        buffer_n = max(int(n_genes) * 5, int(n_genes) + 10)
        df_buf = df.groupby("group", group_keys=False).head(buffer_n)
        genes_needed = pd.unique(df_buf["gene"].astype(str))

        # pull fraction matrix from AnnData
        Xf = (
            adata.layers[fraction_layer]
            if (fraction_layer is not None and fraction_layer in adata.layers)
            else (adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X)
        )
        if hasattr(Xf, "toarray"):
            Xf = Xf[:, adata.var_names.get_indexer(genes_needed)].toarray()
        else:
            Xf = np.asarray(Xf[:, adata.var_names.get_indexer(genes_needed)], dtype=float)

        # group labels
        g = adata.obs[groupby].astype(str)
        cats = pd.unique(df["group"].astype(str)) if groups is None else [str(x) for x in groups]
        cats = [c for c in cats if c in set(g)]

        frac_map = {}  # (group, gene)-> frac
        for i, c in enumerate(cats):
            mask = (g == c).to_numpy()
            if mask.sum() == 0:
                continue
            Xi = Xf[mask, :]
            frac = (Xi > float(expr_threshold)).mean(axis=0)
            for j, gene in enumerate(genes_needed):
                frac_map[(c, str(gene))] = float(frac[j])

        def _passes(row) -> bool:
            f = frac_map.get((str(row["group"]), str(row["gene"])), np.nan)
            if not np.isfinite(f):
                return False
            if (min_in_group_fraction is not None) and (f < float(min_in_group_fraction)):
                return False
            if (max_in_group_fraction is not None) and (f > float(max_in_group_fraction)):
                return False
            return True

        df = df[df.apply(_passes, axis=1)].copy()
        if df.shape[0] == 0:
            raise ValueError("After fraction filtering, no genes remain. Relax min/max fraction.")

    # ---- select top genes per group ----
    df_top = df.groupby("group", group_keys=False).head(int(n_genes)).copy()

    # ---- var_groups dict: group -> [genes...] ----
    var_groups: dict[str, list[str]] = {}
    used = set()
    for gname, sub in df_top.groupby("group", sort=False):
        genes = []
        for gene in sub["gene"].astype(str).tolist():
            if gene not in adata.var_names:
                continue
            if (not unique) or (gene not in used):
                genes.append(gene)
                used.add(gene)
        if genes:
            var_groups[str(gname)] = genes

    if not var_groups:
        raise ValueError("No genes selected (all filtered or not found).")

    # ---- standard_scale auto behavior ----
    if standard_scale == "auto":
        standard_scale_eff = "var" if values_to_plot == "expression" else None
    else:
        standard_scale_eff = standard_scale

    # ---- 1) expression mode: call dotplot normally ----
    if values_to_plot == "expression":
        return dotplot(
            adata,
            var_groups=var_groups,
            groupby=groupby,
            layer=layer,
            fraction_layer=fraction_layer,
            expr_threshold=expr_threshold,
            standard_scale=standard_scale_eff,
            swap_axes=swap_axes,
            dendrogram_top=dendrogram_top,
            dendrogram_rows=dendrogram_rows,
            row_dendrogram_position=row_dendrogram_position,
            row_spacing=row_spacing,
            cmap=cmap,
            save=save,
            show=show,
            title=None,  # avoids overlap with dendrograms
        )

    # ---- 1) logfoldchanges mode: build a temporary matrix where "expression"=logFC ----
    # dotplot aggregates mean per group. We want mean per group to equal the logFC,
    # so we construct a per-sample matrix where each sample in group g has constant value logFC(g, gene).
    # Fraction still computed from fraction_layer on original adata.
    # This keeps dot size (fraction) meaningful and dot color = logFC.
    # Note: requires that group names in df match adata.obs[groupby].
    all_genes = []
    for gname, glist in var_groups.items():
        all_genes.extend(glist)
    all_genes = list(dict.fromkeys(all_genes))  # unique, preserve order
    gidx = adata.var_names.get_indexer(all_genes)

    # Build logFC lookup: (group,gene)->log2FC
    lfc_lookup = {}
    for _, r in df_top.iterrows():
        lfc_lookup[(str(r["group"]), str(r["gene"]))] = float(r["log2FC"])

    # Create matrix (n_obs x n_genes_selected)
    tmp = np.zeros((adata.n_obs, len(all_genes)), dtype=float)
    gobs = adata.obs[groupby].astype(str).to_numpy()
    for i in range(adata.n_obs):
        gi = gobs[i]
        for j, gene in enumerate(all_genes):
            tmp[i, j] = lfc_lookup.get((gi, gene), 0.0)

    ad_tmp = ad.AnnData(
        X=tmp,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(all_genes, name=adata.var_names.name)),
    )
    # carry fraction layer across for dot sizes
    if fraction_layer is not None and fraction_layer in adata.layers:
        Xf = adata.layers[fraction_layer][:, gidx]
        ad_tmp.layers[fraction_layer] = Xf

    return dotplot(
        ad_tmp,
        var_groups=var_groups,
        groupby=groupby,
        layer=None,  # use ad_tmp.X (logFC)
        fraction_layer=fraction_layer,
        expr_threshold=expr_threshold,
        standard_scale=None,  # don't z-score logFC unless user really wants it
        swap_axes=swap_axes,
        dendrogram_top=dendrogram_top,
        dendrogram_rows=dendrogram_rows,
        row_dendrogram_position=row_dendrogram_position,
        row_spacing=row_spacing,
        cmap=cmap,
        save=save,
        show=show,
        title=None,
        colorbar_title="log2FC",
    )