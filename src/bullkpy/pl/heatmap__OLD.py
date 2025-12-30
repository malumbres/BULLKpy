from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

try:
    import seaborn as sns
except Exception as e:  # pragma: no cover
    sns = None

import anndata as ad

from ._style import set_style, _savefig


def _get_matrix(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def heatmap_de(
    adata: ad.AnnData,
    *,
    de_key: str = "de",
    contrast: str = "LumB_vs_LumA",
    results_key: str = "results",
    groupby: str | None = None,
    groups: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    top_n: int = 50,
    mode: Literal["up", "down", "both", "abs"] = "both",
    sort_by: Literal["pval_adj", "pval", "log2fc", "score"] = "pval_adj",
    ascending: bool | None = None,
    z_score: Literal["row", "none"] = "row",
    clip_z: float | None = 3.0,
    center: float = 0.0,
    cmap: str = "vlag",
    figsize: tuple[float, float] | None = None,
    show_sample_labels: bool = False,
    col_colors: str | Sequence[str] | None = None,  # obs key(s) to annotate columns
    dendrogram_rows: bool = True,
    dendrogram_cols: bool = True,
    linewidths: float = 0.0,
    linecolor: str = "white",
    cbar_label: str = "z-scored expression" ,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap of DE genes from adata.uns[de_key][contrast][results_key].

    Expects `results` as a DataFrame-like with at least:
      - gene (or index as gene names)
      - log2FC (or log2fc)
      - pval_adj (or padj) [optional but recommended]

    Parameters
    ----------
    groupby / groups
        If provided, subsets samples (obs) to those groups and orders columns by group.
    layer
        Expression layer used for heatmap (typically log1p_cpm).
    mode
        - "up": top genes with positive log2FC
        - "down": top genes with negative log2FC
        - "both": top_n/2 up + top_n/2 down
        - "abs": top_n by absolute log2FC
    z_score
        "row" z-scores each gene across samples.
    col_colors
        obs column name(s) to annotate samples in seaborn clustermap.
    """
    set_style()

    if sns is None:
        raise ImportError("heatmap_de requires seaborn. Please install seaborn.")

    # --- fetch DE results ---
    try:
        res = adata.uns[de_key][contrast][results_key]
    except KeyError as e:
        raise KeyError(
            f"Could not find DE results at adata.uns['{de_key}']['{contrast}']['{results_key}']"
        ) from e

    if isinstance(res, dict):
        res = pd.DataFrame(res)
    elif not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)

    # --- infer gene column ---
    gene_col = None
    for c in ["gene", "genes", "symbol", "feature", "name"]:
        if c in res.columns:
            gene_col = c
            break

    if gene_col is None:
        # fall back to index
        if res.index.dtype == object:
            genes = res.index.astype(str).tolist()
        else:
            raise ValueError("DE results has no obvious gene column and non-string index.")
    else:
        genes = res[gene_col].astype(str).tolist()

    res = res.copy()
    if gene_col is None:
        res["gene"] = genes
        gene_col = "gene"

    # --- infer log2fc and p columns ---
    lfc_col = "log2FC" if "log2FC" in res.columns else ("log2fc" if "log2fc" in res.columns else None)
    if lfc_col is None:
        raise ValueError("DE results must contain 'log2FC' (or 'log2fc').")

    p_adj_col = None
    for c in ["pval_adj", "padj", "adj_pval", "qval", "fdr"]:
        if c in res.columns:
            p_adj_col = c
            break
    p_col = "pval" if "pval" in res.columns else ("p_value" if "p_value" in res.columns else None)

    if sort_by == "pval_adj" and p_adj_col is None:
        # degrade gracefully
        sort_by = "pval" if p_col is not None else "log2fc"

    if sort_by == "log2fc" and lfc_col is not None:
        sort_col = lfc_col
    else:
        sort_col = sort_by if sort_by in res.columns else (p_adj_col or p_col or lfc_col)

    if ascending is None:
        # p-values ascending, log2fc descending
        if sort_col in [p_adj_col, p_col, "pval_adj", "padj", "pval"]:
            ascending = True
        else:
            ascending = False

    # --- select genes ---
    res = res.dropna(subset=[lfc_col])
    if mode == "up":
        sub = res[res[lfc_col] > 0].sort_values(sort_col, ascending=ascending)
        pick = sub.head(top_n)
    elif mode == "down":
        sub = res[res[lfc_col] < 0].sort_values(sort_col, ascending=ascending)
        pick = sub.head(top_n)
    elif mode == "abs":
        sub = res.copy()
        sub["_abs_lfc"] = np.abs(sub[lfc_col].values)
        pick = sub.sort_values("_abs_lfc", ascending=False).head(top_n)
    else:  # both
        n1 = top_n // 2
        n2 = top_n - n1
        up = res[res[lfc_col] > 0].sort_values(sort_col, ascending=ascending).head(n1)
        down = res[res[lfc_col] < 0].sort_values(sort_col, ascending=ascending).head(n2)
        pick = pd.concat([up, down], axis=0)

    gene_list = pick[gene_col].astype(str).tolist()

    # keep only present genes
    gene_list = [g for g in gene_list if g in adata.var_names]
    if len(gene_list) == 0:
        raise ValueError("None of the selected DE genes are present in adata.var_names.")

    # --- subset samples by groupby/groups (optional) ---
    obs_mask = np.ones(adata.n_obs, dtype=bool)
    col_order = np.arange(adata.n_obs)

    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not found in adata.obs")

        g = adata.obs[groupby].astype(str)
        if groups is not None:
            groups = [str(x) for x in groups]
            obs_mask = g.isin(groups).to_numpy()
        g_sub = g[obs_mask]
        # order by group then by sample name
        col_order = np.argsort(pd.Categorical(g_sub, categories=sorted(g_sub.unique()), ordered=True))
    else:
        g_sub = None

    # --- build expression matrix (genes x samples) ---
    X = _get_matrix(adata, layer)
    gidx = [adata.var_names.get_loc(g) for g in gene_list]
    Xg = X[:, gidx]  # samples x genes
    Xg = Xg[obs_mask, :]
    Xh = Xg.T  # genes x samples

    # --- z-score per gene ---
    if z_score == "row":
        mu = Xh.mean(axis=1, keepdims=True)
        sd = Xh.std(axis=1, keepdims=True, ddof=0)
        sd[sd == 0] = 1.0
        Xh = (Xh - mu) / sd
        if clip_z is not None:
            Xh = np.clip(Xh, -float(clip_z), float(clip_z))

    # --- reorder columns if grouping requested ---
    if groupby is not None:
        Xh = Xh[:, col_order]
        g_sub = g_sub.iloc[col_order]

    # --- dataframe for seaborn ---
    colnames = adata.obs_names[obs_mask].astype(str).tolist()
    if groupby is not None:
        colnames = [colnames[i] for i in col_order]

    df = pd.DataFrame(Xh, index=gene_list, columns=colnames)

    # --- column colors (annotations) ---
    col_colors_df = None
    if col_colors is not None:
        if isinstance(col_colors, str):
            col_colors = [col_colors]
        ann = {}
        for key in col_colors:
            if key not in adata.obs.columns:
                raise KeyError(f"col_colors obs key '{key}' not found")
            v = adata.obs.loc[df.columns, key].astype(str)
            ann[key] = v
        col_colors_df = pd.DataFrame(ann)

    # --- autosize ---
    if figsize is None:
        w = max(6.0, 0.15 * df.shape[1] + 3.0)
        h = max(5.0, 0.18 * df.shape[0] + 2.5)
        figsize = (w, h)

    # --- clustermap ---
    cg = sns.clustermap(
        df,
        cmap=cmap,
        center=center if z_score == "row" else None,
        row_cluster=bool(dendrogram_rows),
        col_cluster=bool(dendrogram_cols),
        col_colors=col_colors_df,
        linewidths=linewidths,
        linecolor=linecolor,
        xticklabels=show_sample_labels,
        yticklabels=True,
        figsize=figsize,
        cbar_kws={"label": cbar_label},
    )

    # title
    cg.ax_heatmap.set_title(f"DE heatmap: {contrast} (top_n={top_n}, mode={mode})", pad=12)

    # If groupby present, add a colorbar-like label row above
    if groupby is not None and g_sub is not None:
        # put a small label on the top
        cg.ax_col_dendrogram.text(
            0.0, 1.02, f"Columns ordered by {groupby}",
            transform=cg.ax_col_dendrogram.transAxes,
            ha="left", va="bottom",
        )

    if save is not None:
        _savefig(cg.fig, save)
    if show:
        plt.show()

    return cg