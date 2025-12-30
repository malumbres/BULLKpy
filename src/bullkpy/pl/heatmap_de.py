from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import to_hex

from ._style import set_style, _savefig, _apply_clustergrid_style

try:
    import seaborn as sns
except Exception:  # pragma: no cover
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
    contrast: str,
    de_key: str = "de",
    results_key: str = "results",
    layer: str | None = "log1p_cpm",
    groupby: str | None = None,
    groups: Sequence[str] | None = None,
    top_n: int = 50,
    mode: Literal["up", "down", "both", "abs"] = "both",
    sort_by: Literal["qval", "pval", "log2FC", "t"] = "qval",
    z_score: Literal["row", "none"] = "row",
    clip_z: float | None = 3.0,
    cmap: str = "vlag",
    center: float = 0.0,
    col_colors: str | Sequence[str] | None = None,  # obs key(s) for column annotation
    dendrogram_rows: bool = True,
    dendrogram_cols: bool = True,
    show_sample_labels: bool = False,
    figsize: tuple[float, float] | None = None,
    cbar_label: str = "z-scored expression",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap of top DE genes using results stored at:
        adata.uns[de_key][contrast][results_key]

    Expected DE columns:
        ['gene', 'log2FC', 't', 'pval', 'qval', 'mean_group', 'mean_ref']

    - Heatmap values come from `layer` (default: log1p_cpm)
    - Gene selection uses DE table (default sort: qval)
    - Optionally subset/order samples by `groupby` and `groups`
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

    if not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)

    required = {"gene", "log2FC"}
    missing = required - set(res.columns)
    if missing:
        raise ValueError(f"DE results missing required columns: {sorted(missing)}")

    # --- gene selection ---
    res = res.copy()
    res["gene"] = res["gene"].astype(str)

    # pick sorting column
    if sort_by not in res.columns:
        # fallback
        sort_by = "qval" if "qval" in res.columns else ("pval" if "pval" in res.columns else "log2FC")

    # p/q ascending, effect size descending
    if sort_by in {"qval", "pval"}:
        res = res.sort_values(sort_by, ascending=True)
    else:
        res = res.sort_values(sort_by, ascending=False)

    if mode == "up":
        pick = res[res["log2FC"] > 0].head(top_n)
    elif mode == "down":
        pick = res[res["log2FC"] < 0].head(top_n)
    elif mode == "abs":
        pick = res.assign(_abs=np.abs(res["log2FC"])).sort_values("_abs", ascending=False).head(top_n)
    else:  # both
        n1 = top_n // 2
        n2 = top_n - n1
        up = res[res["log2FC"] > 0].head(n1)
        down = res[res["log2FC"] < 0].head(n2)
        pick = pd.concat([up, down], axis=0)

    genes = [g for g in pick["gene"].tolist() if g in adata.var_names]
    if len(genes) == 0:
        raise ValueError("None of the selected DE genes are present in adata.var_names.")

    # --- subset/order samples ---
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
        # keep group order as provided if groups given; else alphabetical
        if groups is None:
            cat_order = sorted(g_sub.unique())
        else:
            cat_order = [x for x in groups if x in set(g_sub)]
        col_order = np.argsort(pd.Categorical(g_sub, categories=cat_order, ordered=True))
    else:
        g_sub = None

    # --- expression matrix: genes x samples ---
    X = _get_matrix(adata, layer)
    gidx = [adata.var_names.get_loc(g) for g in genes]
    Xg = X[:, gidx]          # samples x genes
    Xg = Xg[obs_mask, :]     # subset samples
    Xh = Xg.T                # genes x samples

    # --- z-score per gene ---
    if z_score == "row":
        mu = Xh.mean(axis=1, keepdims=True)
        sd = Xh.std(axis=1, keepdims=True, ddof=0)
        sd[sd == 0] = 1.0
        Xh = (Xh - mu) / sd
        if clip_z is not None:
            Xh = np.clip(Xh, -float(clip_z), float(clip_z))

    # order columns
    if groupby is not None:
        Xh = Xh[:, col_order]
        g_sub = g_sub.iloc[col_order]

    # dataframe
    colnames = adata.obs_names[obs_mask].astype(str).tolist()
    if groupby is not None:
        colnames = [colnames[i] for i in col_order]
    df = pd.DataFrame(Xh, index=genes, columns=colnames)

    # --- column annotations ---
    # --- column annotations (map categories -> colors) ---
    col_colors_df = None
    if col_colors is not None:
        if isinstance(col_colors, str):
            col_colors = [col_colors]

        ann_colors = {}
        for key in col_colors:
            if key not in adata.obs.columns:
                raise KeyError(f"col_colors obs key '{key}' not found")

            vals = adata.obs.loc[df.columns, key]

            # categorical/stringify for stable mapping
            vals_str = vals.astype(str)
            cats = pd.Categorical(vals_str).categories

            # build a palette (tab20 is good default)
            palette = sns.color_palette("tab20", n_colors=len(cats))
            lut = {cat: to_hex(palette[i]) for i, cat in enumerate(cats)}

            ann_colors[key] = vals_str.map(lut)

        col_colors_df = pd.DataFrame(ann_colors, index=df.columns)

    # --- autosize ---
    if figsize is None:
        w = max(6.0, 0.15 * df.shape[1] + 3.0)
        h = max(4.8, 0.18 * df.shape[0] + 2.2)
        figsize = (w, h)

    # --- plot ---
    cg = sns.clustermap(
        df,
        cmap=cmap,
        center=center if z_score == "row" else None,
        row_cluster=bool(dendrogram_rows),
        col_cluster=bool(dendrogram_cols),
        col_colors=col_colors_df,
        xticklabels=show_sample_labels,
        yticklabels=True,
        figsize=figsize,
        cbar_kws={"label": cbar_label},
    )

    _apply_clustergrid_style(cg)
    cg.fig.tight_layout()

    # title
    cg.ax_heatmap.set_title(f"{contrast} — top {top_n} ({mode})", pad=10)

    if save is not None:
        _savefig(cg.fig, save)
    if show:
        plt.show()

    return cg