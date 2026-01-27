from __future__ import annotations

from typing import Optional, Sequence, Literal, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData

# Optional scipy for clustering/dendrogram ordering
try:
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False



def axis_vs_dispersion(
    adata: AnnData,
    *,
    axis_key: str,
    dispersion_key: str,
    groupby: Optional[str] = None,
    max_groups: int = 12,
    title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Scatter plot of signature axis vs per-sample dispersion.
    If `groupby` is provided, colors by group (limited to max_groups).
    """
    if axis_key not in adata.obs or dispersion_key not in adata.obs:
        raise KeyError("axis_key and/or dispersion_key not found in adata.obs")

    x = adata.obs[axis_key].to_numpy(dtype=float)
    y = adata.obs[dispersion_key].to_numpy(dtype=float)

    fig, ax = plt.subplots()
    if groupby is None:
        ax.scatter(x, y, s=10, alpha=0.7)
    else:
        if groupby not in adata.obs:
            raise KeyError(f"{groupby!r} not found in adata.obs")
        cats = adata.obs[groupby].astype("category")
        levels = list(cats.cat.categories)[:max_groups]
        for lvl in levels:
            idx = (cats == lvl).to_numpy()
            ax.scatter(x[idx], y[idx], s=10, alpha=0.7, label=str(lvl))
        ax.legend(frameon=False, fontsize=7, ncol=2)

    ax.set_xlabel(axis_key)
    ax.set_ylabel(dispersion_key)
    ax.set_title(title or f"{axis_key} vs {dispersion_key}")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def dispersion_summary(
    adata: AnnData,
    *,
    dispersion_key: str,
    groupby: str,
    title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Boxplot of per-sample dispersion across groups (e.g., tumor types).
    """
    if dispersion_key not in adata.obs:
        raise KeyError(f"{dispersion_key!r} not found in adata.obs")
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    df = adata.obs[[groupby, dispersion_key]].dropna()
    groups = df[groupby].astype("category")
    order = df.groupby(groupby)[dispersion_key].median().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(8, 3))
    data = [df.loc[df[groupby] == g, dispersion_key].to_numpy() for g in order]
    ax.boxplot(data, labels=[str(g) for g in order], showfliers=False)
    ax.set_ylabel(dispersion_key)
    ax.set_title(title or f"{dispersion_key} by {groupby}")
    ax.tick_params(axis="x", labelrotation=90)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def _get_mp_matrix(
    adata,
    *,
    obsm_key: str = "X_mp",
    mp_names: Sequence[str] | None = None,
    use_obs_prefix: str | None = None,
) -> tuple[np.ndarray, list[str]]:
    """
    Returns (M, names) where M is (n_samples, n_mps).
    Priority:
      1) adata.obsm[obsm_key] + mp_names (or adata.uns['metaprograms']['mp_names'])
      2) adata.obs columns starting with use_obs_prefix
    """
    if use_obs_prefix is not None:
        cols = [c for c in adata.obs.columns if str(c).startswith(str(use_obs_prefix))]
        if len(cols) == 0:
            raise KeyError(f"No adata.obs columns start with prefix='{use_obs_prefix}'.")
        df = adata.obs[cols].apply(pd.to_numeric, errors="coerce")
        M = df.to_numpy(float)
        names = list(df.columns.astype(str))
        return M, names

    if obsm_key not in getattr(adata, "obsm", {}):
        raise KeyError(
            f"Neither adata.obsm['{obsm_key}'] found nor use_obs_prefix provided."
        )
    M = np.asarray(adata.obsm[obsm_key], float)

    if mp_names is None:
        mp_names = (
            adata.uns.get("metaprograms", {}).get("mp_names", None)
            if hasattr(adata, "uns")
            else None
        )
    if mp_names is None:
        mp_names = [f"MP{i+1}" for i in range(M.shape[1])]
    mp_names = [str(x) for x in list(mp_names)]

    if len(mp_names) != M.shape[1]:
        raise ValueError(
            f"mp_names length ({len(mp_names)}) does not match M.shape[1] ({M.shape[1]})."
        )
    return M, mp_names


def _cluster_order_from_corr(C: np.ndarray, method: str = "average") -> np.ndarray:
    """
    C is (p,p) correlation matrix.
    Returns an ordering of indices. Uses hierarchical clustering if scipy available.
    """
    p = C.shape[0]
    if p <= 2:
        return np.arange(p)

    # distance = 1 - corr (clip to [0,2])
    D = 1.0 - np.clip(C, -1.0, 1.0)
    D = np.clip(D, 0.0, 2.0)

    if _HAVE_SCIPY:
        # squareform expects condensed distance
        d = squareform(D, checks=False)
        Z = linkage(d, method=method)
        return leaves_list(Z)

    # fallback: order by first principal component of corr
    # (very rough, but deterministic)
    vals, vecs = np.linalg.eigh(C)
    pc1 = vecs[:, np.argmax(vals)]
    return np.argsort(pc1)


def metaprogram_corr(
    adata,
    *,
    obsm_key: str = "X_mp",
    mp_names: Sequence[str] | None = None,
    use_obs_prefix: str | None = None,
    method: Literal["pearson", "spearman"] = "pearson",
    cluster: bool = True,
    cluster_method: str = "average",
    vmin: float = -1.0,
    vmax: float = 1.0,
    figsize: tuple[float, float] = (7.5, 6.5),
    title: str | None = None,
    ax=None,
    show_colorbar: bool = True,
    rotate_xticks: int = 90,
):
    """
    Plot MP-MP correlation matrix heatmap (optionally clustered).

    Returns
    -------
    df_corr : pd.DataFrame
        Correlation matrix (possibly re-ordered).
    """
    M, names = _get_mp_matrix(
        adata, obsm_key=obsm_key, mp_names=mp_names, use_obs_prefix=use_obs_prefix
    )

    df = pd.DataFrame(M, columns=names)
    if str(method).lower() == "spearman":
        C = df.corr(method="spearman").to_numpy(float)
    else:
        C = df.corr(method="pearson").to_numpy(float)

    order = np.arange(len(names))
    if bool(cluster):
        order = _cluster_order_from_corr(C, method=str(cluster_method))

    names_ord = [names[i] for i in order]
    C_ord = C[np.ix_(order, order)]
    df_corr = pd.DataFrame(C_ord, index=names_ord, columns=names_ord)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(C_ord, vmin=float(vmin), vmax=float(vmax), aspect="auto")
    ax.set_xticks(np.arange(len(names_ord)))
    ax.set_yticks(np.arange(len(names_ord)))
    ax.set_xticklabels(names_ord, rotation=int(rotate_xticks), ha="right")
    ax.set_yticklabels(names_ord)

    ax.set_xlabel("Metaprograms")
    ax.set_ylabel("Metaprograms")
    if title is None:
        title = f"Metaprogram correlation ({method})"
    ax.set_title(title)

    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Correlation")

    plt.tight_layout()
    return df_corr


def metaprogram_heatmap(
    adata,
    *,
    obsm_key="X_mp",
    mp_names=None,
    use_obs_prefix=False,
    groupby="Project_ID",
    agg="mean",                      # "mean" or "median"
    scale="zscore_cols",             # None, "zscore_cols", "zscore_rows"
    order_groups="size",             # "size", "alpha", None, "cluster"
    order_mps="cluster",             # "alpha", None, "cluster"
    cluster_method="average",        # linkage method
    cluster_metric="correlation",    # distance metric
    row_dendrogram=True,
    col_dendrogram=True,
    dendrogram_ratio=(0.18, 0.18),   # (row_dendro_width, col_dendro_height) in figure fraction-ish
    cmap="viridis",
    vmin=None,
    vmax=None,
    row_dendro_width=0.28,
    col_dendro_height=0.18,
    dendro_gap=0.03,
    ytick_pad=6,
    figsize=(12, 6),
    title=None,
    ax=None,                         # if provided, dendrograms are disabled unless you pass ax=None
    show_colorbar=True,
    rotate_xticks=45,
    min_n_per_group=1,
):
    """
    Heatmap of metaprogram (MP) scores aggregated by group, with optional row/col dendrograms.

    Returns:
      df_heat: aggregated (and reordered/scaled) DataFrame [groups x MPs]
      df_group_info: DataFrame with group sizes and ordering
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # ---- optional scipy for dendrograms ----
    if row_dendrogram or col_dendrogram or order_groups == "cluster" or order_mps == "cluster":
        try:
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import pdist
        except Exception as e:
            raise ImportError(f"metaprogram_heatmap: scipy is required for clustering/dendrograms. ({e})")

    # ---------- pull MP matrix ----------
    if obsm_key not in adata.obsm_keys():
        raise KeyError(f"obsm_key='{obsm_key}' not found in adata.obsm")

    X = np.asarray(adata.obsm[obsm_key])
    if X.ndim != 2:
        raise ValueError(f"adata.obsm['{obsm_key}'] must be 2D")

    n, p = X.shape

    # infer mp_names
    if mp_names is None:
        if use_obs_prefix:
            # attempt to infer from obs columns like "MP_*" (not recommended unless you store them)
            mp_names = [f"MP{i+1}" for i in range(p)]
        else:
            mp_names = [f"MP{i+1}" for i in range(p)]
    else:
        mp_names = [str(m) for m in mp_names]
        if len(mp_names) != p:
            raise ValueError(f"mp_names length {len(mp_names)} != X_mp columns {p}")

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    # ---------- build df2 ----------
    df2 = pd.DataFrame(X, columns=mp_names, index=adata.obs_names)
    df2[groupby] = adata.obs[groupby].astype("object").to_numpy()

    # group sizes (robust)
    gcounts = (
        df2[groupby]
        .value_counts()
        .rename_axis("group")
        .reset_index(name="n")
    )
    keep_groups = set(gcounts.loc[gcounts["n"] >= int(min_n_per_group), "group"].tolist())
    df2 = df2[df2[groupby].isin(keep_groups)].copy()
    gcounts = gcounts[gcounts["group"].isin(keep_groups)].reset_index(drop=True)

    if df2.shape[0] == 0:
        raise ValueError("No samples left after min_n_per_group filtering.")

    # ---------- aggregate ----------
    if str(agg).lower() == "median":
        df_heat = df2.groupby(groupby)[mp_names].median()
    else:
        df_heat = df2.groupby(groupby)[mp_names].mean()

    # update group info
    df_group_info = gcounts.copy()
    df_group_info = df_group_info.rename(columns={"group": groupby})

    # ---------- scaling ----------
    H = df_heat.to_numpy(float)

    if scale is None or str(scale).lower() == "none":
        Hs = H.copy()
    elif str(scale).lower() == "zscore_cols":
        mu = np.nanmean(H, axis=0, keepdims=True)
        sd = np.nanstd(H, axis=0, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        Hs = (H - mu) / sd
    elif str(scale).lower() == "zscore_rows":
        mu = np.nanmean(H, axis=1, keepdims=True)
        sd = np.nanstd(H, axis=1, keepdims=True)
        sd = np.where(sd == 0, 1.0, sd)
        Hs = (H - mu) / sd
    else:
        raise ValueError("scale must be None/'none', 'zscore_cols', or 'zscore_rows'")

    df_heat_scaled = pd.DataFrame(Hs, index=df_heat.index.astype(str), columns=df_heat.columns.astype(str))

    # ---------- ordering helpers ----------
    def _cluster_order(mat_2d, axis="rows"):
        # mat_2d is dataframe
        A = mat_2d.to_numpy(float)
        if axis == "cols":
            A = A.T
        # pdist expects rows are observations
        D = pdist(A, metric=str(cluster_metric))
        Z = linkage(D, method=str(cluster_method))
        leaves = dendrogram(Z, no_plot=True)["leaves"]
        return leaves, Z

    # initial order
    row_order = list(range(df_heat_scaled.shape[0]))
    col_order = list(range(df_heat_scaled.shape[1]))
    Zr = Zc = None

    # order columns
    if order_mps == "alpha":
        col_order = list(np.argsort(df_heat_scaled.columns.to_numpy()))
    elif order_mps == "cluster":
        col_order, Zc = _cluster_order(df_heat_scaled, axis="cols")

    # order rows
    if order_groups == "alpha":
        row_order = list(np.argsort(df_heat_scaled.index.to_numpy()))
    elif order_groups == "size":
        # use df_group_info
        size_map = dict(zip(df_group_info[groupby].astype(str), df_group_info["n"].astype(int)))
        sizes = np.array([size_map.get(str(g), 0) for g in df_heat_scaled.index], dtype=int)
        row_order = list(np.argsort(-sizes))  # descending
    elif order_groups == "cluster":
        row_order, Zr = _cluster_order(df_heat_scaled, axis="rows")

    # apply ordering
    df_heat_scaled = df_heat_scaled.iloc[row_order, :].iloc[:, col_order]

    # reorder group info to match
    df_group_info = df_group_info.copy()
    df_group_info[groupby] = df_group_info[groupby].astype(str)
    df_group_info = df_group_info.set_index(groupby).reindex(df_heat_scaled.index).reset_index()

    # ---------- plotting ----------
    # If ax is provided, we draw only the heatmap (no dendrograms)
    if ax is not None:
        row_dendrogram = False
        col_dendrogram = False

    if ax is None:
        fig = plt.figure(figsize=figsize)

        # user-tunable layout knobs (add these as function params if you want)
        row_dendro_width = row_dendro_width
        col_dendro_height = col_dendro_height
        dendro_gap = dendro_gap         # spacer between row dendrogram and heatmap
        cb_w = 0.03 if show_colorbar else 0.001

        top_h = col_dendro_height if col_dendrogram else 0.001
        left_w = row_dendro_width if row_dendrogram else 0.001

        import matplotlib.gridspec as gridspec

        # NOTE: 4 columns now: [row dendro | GAP | heatmap | colorbar]
        gs = gridspec.GridSpec(
            2, 4,
            height_ratios=[top_h, 1.0],
            width_ratios=[left_w, dendro_gap, 1.0, cb_w],
            hspace=0.02,
            wspace=0.02,   # small global spacing; gap column does the main separation
        )

        ax_col = fig.add_subplot(gs[0, 2]) if col_dendrogram else None
        ax_row = fig.add_subplot(gs[1, 0]) if row_dendrogram else None
        ax_gap = fig.add_subplot(gs[1, 1])  # spacer column
        ax_hm  = fig.add_subplot(gs[1, 2])
        ax_cb  = fig.add_subplot(gs[1, 3]) if show_colorbar else None

        ax_gap.axis("off")
        if ax_col is not None:
            ax_col.axis("off")
        if ax_row is not None:
            ax_row.axis("off")
    else:
        ax_hm = ax
        ax_cb = None
        ax_col = ax_row = None


    # ---------- dendrograms ----------
    if col_dendrogram:
        if Zc is None:
            _, Zc = _cluster_order(df_heat_scaled, axis="cols")
        dd = dendrogram(Zc, ax=ax_col, no_labels=True, color_threshold=None)
        ax_col.set_xticks([]); ax_col.set_yticks([])
        # clip lines to axis
        for line in ax_col.get_lines():
            line.set_clip_on(True)

    if row_dendrogram:
        if Zr is None:
            _, Zr = _cluster_order(df_heat_scaled, axis="rows")

        dendrogram(
            Zr,
            ax=ax_row,
            orientation="right",
            no_labels=True,
            color_threshold=None,
        )

        ax_row.invert_yaxis()
        ax_row.set_xlim(ax_row.get_xlim()[::-1])
        ax_row.set_xticks([])
        ax_row.set_yticks([])

        for line in ax_row.get_lines():
            line.set_clip_on(True)

    # Heatmap
    im = ax_hm.imshow(
        df_heat_scaled.to_numpy(float),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )

    ax_hm.set_yticks(np.arange(df_heat_scaled.shape[0]))
    ax_hm.set_yticklabels(df_heat_scaled.index.tolist(), fontsize=9)

    # push y labels left so they don't collide visually with dendrogram area
    ax_hm.tick_params(axis="y", pad=ytick_pad)
    ax_hm.set_xticks(np.arange(df_heat_scaled.shape[1]))
    ax_hm.set_xticklabels(df_heat_scaled.columns.tolist(), fontsize=9, rotation=rotate_xticks, ha="right")

    ax_hm.set_xlabel("Metaprograms")
    ax_hm.set_ylabel("")

    if title is not None:
        ax_hm.set_title(title)

    if show_colorbar and ax_cb is not None:
        plt.colorbar(im, cax=ax_cb)

    # return the scaled/reordered heat matrix and group info
    return df_heat_scaled, df_group_info



# Backward-compatible names for your requested API
def metaprogram_heatmap_plot(*args, **kwargs):
    return metaprogram_heatmap(*args, **kwargs)

def metaprogram_corr_plot(*args, **kwargs):
    return metaprogram_corr(*args, **kwargs)


##########################
### NEW METAPROGRAMS
##########################


def metaprogram_dispersion_heatmap(
    adata: AnnData,
    *,
    uns_key: str = "metaprogram_dispersion",
    sort_groups_by: str = "median",   # "median" | "name"
    sort_programs_by: str = "median", # "median" | "name"
    title: str = "Metaprogram dispersion by tumor type",
    cmap="viridis",
    figsize=(8, 6),
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Heatmap of per-group dispersion (rows: metaprograms, cols: groups).

    Requires adata.uns[uns_key] with columns: group, metaprogram, dispersion.
    """
    if uns_key not in adata.uns:
        raise KeyError(f"{uns_key!r} not found in adata.uns")
    df = adata.uns[uns_key]
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    req = {"group", "metaprogram", "dispersion"}
    if not req.issubset(df.columns):
        raise ValueError(f"{uns_key} must contain columns {req}")

    mat = df.pivot(index="metaprogram", columns="group", values="dispersion")

    # sorting
    if sort_groups_by == "median":
        mat = mat.loc[:, mat.median(axis=0, skipna=True).sort_values(ascending=False).index]
    elif sort_groups_by == "name":
        mat = mat.loc[:, sorted(mat.columns)]
    if sort_programs_by == "median":
        mat = mat.loc[mat.median(axis=1, skipna=True).sort_values(ascending=False).index, :]
    elif sort_programs_by == "name":
        mat = mat.loc[sorted(mat.index), :]

    arr = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(arr, aspect="auto", cmap=cmap)
    ax.set_title(title)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index.tolist(), fontsize=7)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns.tolist(), fontsize=7, rotation=90)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Dispersion (MAD)", fontsize=8)

    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def metaprogram_metrics_summary(
    adata: AnnData,
    *,
    groupby: str = "Project_ID",
    heterogeneity_key: str = "mp_heterogeneity_entropy",
    dispersion_key: str = "mp_dispersion_mad_zwithin_Project_ID",
    title: str = "Metaprogram heterogeneity and dispersion by tumor type",
    figsize=(12, 3),
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Two-panel summary (boxplots):
      - heterogeneity by group
      - dispersion by group
    """
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")
    if heterogeneity_key not in adata.obs:
        raise KeyError(f"{heterogeneity_key!r} not found in adata.obs")
    if dispersion_key not in adata.obs:
        raise KeyError(f"{dispersion_key!r} not found in adata.obs")

    df = adata.obs[[groupby, heterogeneity_key, dispersion_key]].dropna()
    order = df.groupby(groupby)[dispersion_key].median().sort_values(ascending=False).index.tolist()

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=False)
    fig.suptitle(title, y=1.02, fontsize=10)

    # heterogeneity
    data_h = [df.loc[df[groupby] == g, heterogeneity_key].to_numpy() for g in order]
    axes[0].boxplot(data_h, labels=order, showfliers=False)
    axes[0].set_title("Per-sample heterogeneity (entropy)")
    axes[0].set_ylabel(heterogeneity_key)
    axes[0].tick_params(axis="x", labelrotation=90)

    # dispersion
    data_d = [df.loc[df[groupby] == g, dispersion_key].to_numpy() for g in order]
    axes[1].boxplot(data_d, labels=order, showfliers=False)
    axes[1].set_title("Per-sample dispersion (MAD of z-scored MPs)")
    axes[1].set_ylabel(dispersion_key)
    axes[1].tick_params(axis="x", labelrotation=90)

    plt.tight_layout()

    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return axes


def metaprogram_ne_scatter(
    adata: AnnData,
    *,
    ne_key: str,
    x_key: str = "mp_heterogeneity_entropy",
    groupby: Optional[str] = "Project_ID",
    max_groups: int = 10,
    title: Optional[str] = None,
    figsize=(5, 4),
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Scatter: NE score vs metaprogram heterogeneity (or dispersion).
    Colors by group (subset of groups for readability).
    """
    if ne_key not in adata.obs:
        raise KeyError(f"{ne_key!r} not found in adata.obs")
    if x_key not in adata.obs:
        raise KeyError(f"{x_key!r} not found in adata.obs")

    x = adata.obs[x_key].to_numpy(dtype=float)
    y = adata.obs[ne_key].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=figsize)
    if groupby is None:
        ax.scatter(x, y, s=10, alpha=0.7)
    else:
        if groupby not in adata.obs:
            raise KeyError(f"{groupby!r} not found in adata.obs")
        cats = adata.obs[groupby].astype("category")
        levels = list(cats.cat.categories)[:max_groups]
        for lvl in levels:
            idx = (cats == lvl).to_numpy()
            ax.scatter(x[idx], y[idx], s=10, alpha=0.7, label=str(lvl))
        ax.legend(frameon=False, fontsize=7, ncol=2)

    ax.set_xlabel(x_key)
    ax.set_ylabel(ne_key)
    ax.set_title(title or f"{ne_key} vs {x_key}")

    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def metaprogram_dominance_ridgeplot_like(
    adata,
    *,
    groupby: str = "Project_ID",
    topk_key: str = "mp_topk",
    rank: int = 1,
    weight_key: str = "weight",
    figsize=(10, 4),
    ylim=None,
    y_text=None,
    label_rotation=90,
    fontsize=8,
    show=True,
    save=None,
):
    """
    Ridge/violin-style plot of dominance of the top metaprogram (rank-1 weight)
    across groups (e.g. tumor types).

    Adds:
      - dominant metaprogram label on top of each violin
      - configurable y-limits and text height
      - figsize control
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    if topk_key not in adata.uns:
        raise KeyError(f"{topk_key!r} not found in adata.uns")

    df = adata.uns[topk_key].copy()
    df = df[df["rank"] == rank].copy()

    # map group
    df[groupby] = adata.obs.loc[df["sample"], groupby].astype(str).values

    # order groups by median dominance
    order = (
        df.groupby(groupby)[weight_key]
          .median()
          .sort_values(ascending=False)
          .index.tolist()
    )

    data = [df.loc[df[groupby] == g, weight_key].values for g in order]

    fig, ax = plt.subplots(figsize=figsize)

    parts = ax.violinplot(
        data,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    # style violins
    for pc in parts["bodies"]:
        pc.set_facecolor("#d95f5f")
        pc.set_alpha(0.35)
        pc.set_edgecolor("none")

    # medians
    parts["cmedians"].set_color("#b22222")
    parts["cmedians"].set_linewidth(2)

    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90)
    ax.set_ylabel("Top metaprogram weight (softmax)")
    ax.set_title("Dominance of the top metaprogram across tumor types (rank-1 weight)")

    # y-limits
    ymax = max(np.nanmax(d) for d in data if len(d) > 0)
    if ylim is None:
        ax.set_ylim(0, ymax * 1.25)
    else:
        ax.set_ylim(*ylim)

    # text height
    if y_text is None:
        y_text = ax.get_ylim()[1] * 0.97

    # ---- add dominant metaprogram labels
    for i, g in enumerate(order, start=1):
        sub = df[df[groupby] == g]
        if sub.empty:
            continue
        top_mp = sub["metaprogram"].value_counts().idxmax()
        ax.text(
            i,
            y_text,
            top_mp,
            ha="center",
            va="top",
            rotation=label_rotation,
            fontsize=fontsize,
        )

    plt.tight_layout()
    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    return ax