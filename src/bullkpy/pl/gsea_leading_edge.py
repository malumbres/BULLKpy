from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd

try:
    import seaborn as sns  # type: ignore
except Exception:
    sns = None

import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    import scipy.sparse as sp  # type: ignore
except Exception:
    sp = None

try:
    from scipy.cluster.hierarchy import linkage, fcluster  # type: ignore
    from scipy.spatial.distance import squareform  # type: ignore
except Exception:
    linkage = None
    fcluster = None
    squareform = None

try:
    from sklearn.manifold import MDS
except Exception:
    MDS = None

# ---------------------------------------------------------------------
# Your existing style helpers (adapt as needed)
# ---------------------------------------------------------------------
def set_style():
    """Minimal style hook (use your existing bullkpy set_style if available)."""
    try:
        plt.rcParams["figure.dpi"] = 120
    except Exception:
        pass


def _savefig(fig, save: str | Path):
    save = Path(save)
    save.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save, bbox_inches="tight")


def _rotate_gene_labels(ax: plt.Axes, fontsize: float = 7.0):
    for lab in ax.get_xticklabels():
        lab.set_rotation(90)
        lab.set_ha("center")
        lab.set_va("top")
        lab.set_fontsize(float(fontsize))


# ---------------------------------------------------------------------
# GSEA result handling (supports pre_res object OR DataFrame)
# ---------------------------------------------------------------------
def _get_res2d(pre_res) -> pd.DataFrame:
    """
    Accepts either:
      - gseapy prerank result object with .res2d
      - pandas DataFrame already (res2d-like)
    """
    if pre_res is None:
        raise ValueError("pre_res cannot be None")
    if isinstance(pre_res, pd.DataFrame):
        return pre_res
    if hasattr(pre_res, "res2d"):
        df = pre_res.res2d
        if not isinstance(df, pd.DataFrame):
            raise TypeError("pre_res.res2d must be a pandas DataFrame")
        return df
    raise TypeError("pre_res must be a gseapy.prerank result object with .res2d OR a pandas DataFrame.")


def _find_col(df: pd.DataFrame, candidates: Sequence[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of the candidate columns found: {list(candidates)}")


def _split_leading_edge(x) -> list[str]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return []
    # common separators: '/', ';', ',', whitespace
    for sep in [";", ",", "/", "|"]:
        s = s.replace(sep, " ")
    genes = [g for g in s.split() if g]
    return genes


def _normalize_term_selector(
    df: pd.DataFrame,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Filter by either:
      - terms (explicit names), OR
      - term_idx (slice | int | Sequence[int]) using iloc positional indices.
    """
    if terms is not None and len(terms) > 0:
        term_col = _find_col(df, ["Term", "term", "pathway", "Pathway", "NAME", "name"])
        terms_set = set(map(str, terms))
        out = df[df[term_col].astype(str).isin(terms_set)].copy()
        if out.shape[0] == 0:
            raise ValueError("No rows matched the provided `terms`.")
        return out

    if term_idx is None:
        return df

    if isinstance(term_idx, int):
        out = df.iloc[[int(term_idx)]].copy()
    elif isinstance(term_idx, slice):
        out = df.iloc[term_idx].copy()
    elif isinstance(term_idx, (list, tuple, np.ndarray, pd.Index)):
        idx = [int(i) for i in list(term_idx)]
        n = df.shape[0]
        bad = [i for i in idx if i < 0 or i >= n]
        if bad:
            raise IndexError(f"term_idx contains out-of-range indices (0..{n-1}): {bad[:10]}")
        out = df.iloc[idx].copy()
    else:
        raise TypeError("term_idx must be None | int | slice | Sequence[int]")

    if out.shape[0] == 0:
        raise ValueError("term_idx selection produced an empty result table.")
    return out


def _leading_edge_sets(
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
) -> tuple[list[str], dict[str, set[str]]]:
    """
    Returns:
      - term_names (ordered)
      - dict term -> set(leading-edge genes)
    """
    df = _get_res2d(pre_res)
    df = _normalize_term_selector(df, term_idx=term_idx, terms=terms)

    term_col = _find_col(df, ["Term", "term", "pathway", "Pathway", "NAME", "name"])
    le_col = _find_col(df, ["Lead_genes", "lead_genes", "ledge_genes", "ledge", "Lead_genes "])

    term_names: list[str] = []
    le_sets: dict[str, set[str]] = {}

    for _, row in df.iterrows():
        t = str(row[term_col])
        genes = _split_leading_edge(row[le_col])
        term_names.append(t)
        le_sets[t] = set(map(str, genes))

    term_names = [t for t in term_names if len(le_sets.get(t, set())) > 0]
    le_sets = {t: le_sets[t] for t in term_names}

    if len(term_names) == 0:
        raise ValueError("No leading-edge genes found for the selected terms/indices.")
    return term_names, le_sets


# ---------------------------------------------------------------------
# Core computation helpers
# ---------------------------------------------------------------------
def _jaccard_matrix(term_names: list[str], le_sets: dict[str, set[str]]) -> pd.DataFrame:
    n = len(term_names)
    J = np.zeros((n, n), dtype=float)
    for i, ti in enumerate(term_names):
        Ai = le_sets[ti]
        for j, tj in enumerate(term_names):
            Aj = le_sets[tj]
            inter = len(Ai & Aj)
            union = len(Ai | Aj)
            J[i, j] = (inter / union) if union > 0 else 0.0
    return pd.DataFrame(J, index=term_names, columns=term_names)


def _cohesion(dfJ: pd.DataFrame, labels: np.ndarray) -> pd.DataFrame:
    """
    Simple cluster cohesion metrics from similarity matrix.
    - within_mean: mean similarity within cluster (off-diagonal)
    - between_mean: mean similarity to outside cluster
    """
    terms = dfJ.index.to_list()
    out = []
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        if idx.size <= 1:
            within = np.nan
        else:
            sub = dfJ.to_numpy()[np.ix_(idx, idx)]
            # exclude diagonal
            within = float((sub.sum() - np.trace(sub)) / (idx.size * (idx.size - 1)))
        # between
        other = np.where(labels != c)[0]
        if other.size == 0:
            between = np.nan
        else:
            subb = dfJ.to_numpy()[np.ix_(idx, other)]
            between = float(np.nanmean(subb)) if subb.size else np.nan
        out.append({"cluster": int(c), "n_terms": int(idx.size), "within_mean": within, "between_mean": between})
    return pd.DataFrame(out).sort_values(["n_terms", "within_mean"], ascending=[False, False]).reset_index(drop=True)


# ---------------------------------------------------------------------
# NEW 1) Pathway clustering (nodules) + metrics
# ---------------------------------------------------------------------
def leading_edge_pathway_clusters(
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
    method: str = "average",
    threshold: float | None = None,
    n_clusters: int | None = None,
    min_shared_genes: int = 0,
) -> dict:
    """
    Cluster pathways into 'nodules' using leading-edge Jaccard similarity.

    You can specify:
      - threshold: cut dendrogram by distance (distance = 1 - Jaccard), OR
      - n_clusters: request a fixed number of clusters.

    Returns
    -------
    dict with:
      - term_names
      - le_sets
      - dfJ (Jaccard similarity)
      - clusters (pd.Series: pathway -> cluster_id)
      - metrics (cluster cohesion table)
    """
    if linkage is None or fcluster is None or squareform is None:
        raise ImportError("leading_edge_pathway_clusters requires scipy (cluster.hierarchy + distance).")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)
    dfJ = _jaccard_matrix(term_names, le_sets)

    if int(min_shared_genes) > 0:
        # de-noise: zero similarities if intersection too small
        for i, ti in enumerate(term_names):
            for j, tj in enumerate(term_names):
                if i == j:
                    continue
                if len(le_sets[ti] & le_sets[tj]) < int(min_shared_genes):
                    dfJ.iloc[i, j] = 0.0

    D = 1.0 - dfJ.to_numpy()
    # ensure valid condensed form
    d = squareform(D, checks=False)
    Z = linkage(d, method=method)

    if (threshold is None) == (n_clusters is None):
        raise ValueError("Provide exactly one of: threshold OR n_clusters.")

    if threshold is not None:
        # distance threshold (smaller threshold -> more clusters)
        labels = fcluster(Z, t=float(threshold), criterion="distance")
    else:
        labels = fcluster(Z, t=int(n_clusters), criterion="maxclust")

    clusters = pd.Series(labels, index=term_names, name="cluster").astype(int)
    metrics = _cohesion(dfJ, labels)

    return {
        "term_names": term_names,
        "le_sets": le_sets,
        "dfJ": dfJ,
        "clusters": clusters,
        "metrics": metrics,
        "linkage": Z,
    }


# ---------------------------------------------------------------------
# NEW 2) High-level leading-edge genes per pathway nodule
# ---------------------------------------------------------------------
def leading_edge_cluster_driver_genes(
    le_sets: dict[str, set[str]],
    clusters: pd.Series,
    *,
    top_n: int = 25,
    min_gene_freq_in_cluster: int = 2,
) -> pd.DataFrame:
    """
    Compute 'driver genes' for each pathway cluster (nodule).

    Scoring:
      - freq_in_cluster: in how many pathways (within the cluster) gene appears
      - coverage_in_cluster: freq / n_terms_in_cluster
      - freq_global: how many pathways overall contain the gene
      - tfidf_like: coverage_in_cluster * log( (1 + n_terms_total) / (1 + freq_global) )

    Returns
    -------
    DataFrame with columns:
      cluster, gene, freq_in_cluster, coverage_in_cluster, freq_global, tfidf_like
    """
    clusters = clusters.astype(int)
    term_to_cluster = clusters.to_dict()

    # global freq
    global_freq: dict[str, int] = {}
    for t, gs in le_sets.items():
        for g in gs:
            global_freq[g] = global_freq.get(g, 0) + 1
    n_terms_total = len(le_sets)

    # cluster sizes
    cl_terms: dict[int, list[str]] = {}
    for t, c in term_to_cluster.items():
        cl_terms.setdefault(int(c), []).append(t)

    rows = []
    for c, terms_in_c in cl_terms.items():
        # cluster freq
        freq_c: dict[str, int] = {}
        for t in terms_in_c:
            for g in le_sets.get(t, set()):
                freq_c[g] = freq_c.get(g, 0) + 1

        n_terms_c = max(1, len(terms_in_c))
        for g, f in freq_c.items():
            if int(f) < int(min_gene_freq_in_cluster):
                continue
            cov = float(f) / float(n_terms_c)
            fg = int(global_freq.get(g, 0))
            tfidf = cov * float(np.log((1.0 + n_terms_total) / (1.0 + fg)))
            rows.append(
                {
                    "cluster": int(c),
                    "gene": str(g),
                    "freq_in_cluster": int(f),
                    "coverage_in_cluster": cov,
                    "freq_global": int(fg),
                    "tfidf_like": tfidf,
                    "n_terms_in_cluster": int(n_terms_c),
                }
            )

    df = pd.DataFrame(rows)
    if df.shape[0] == 0:
        return df

    # rank within cluster
    df = df.sort_values(["cluster", "tfidf_like", "freq_in_cluster", "gene"], ascending=[True, False, False, True])
    df = df.groupby("cluster", group_keys=False).head(int(top_n)).reset_index(drop=True)
    return df




# ---------------------------------------------------------------------
# Plot 1: Leading-edge expression heatmap (adata required)
# ---------------------------------------------------------------------
def gsea_leading_edge_heatmap(
    adata,
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    use: str = "samples",  # "samples" | "group_mean"
    groupby: str | None = None,
    min_gene_freq: int = 1,
    max_genes: int | None = 200,
    z_score: str | None = "row",  # "row" | None
    clip_z: float | None = 3.0,
    row_cluster: bool = True,
    col_cluster: bool = True,
    cmap: str = "vlag",
    figsize: tuple[float, float] | None = None,
    show_labels: bool = False,
    gene_label_fontsize: float = 7.0,
    label_fontsize: float = 9.0,
    dendrogram_ratio: tuple[float, float] = (0.06, 0.12),
    title: str | None = None,
    show_title: bool = True,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap of expression for leading-edge genes from selected pathways.
    Supports pre_res object OR DataFrame.

    Enhancements:
      - dendrogram_ratio to shrink left dendrogram
      - label_fontsize controls axis labels
      - show_title toggles title
    """
    set_style()
    if sns is None:
        raise ImportError("gsea_leading_edge_heatmap requires seaborn. Please install seaborn.")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)

    # gene frequency across pathways
    freq: dict[str, int] = {}
    for s in le_sets.values():
        for g in s:
            freq[g] = freq.get(g, 0) + 1

    genes = [g for g, f in freq.items() if f >= int(min_gene_freq)]
    if len(genes) == 0:
        raise ValueError(f"No genes pass min_gene_freq={min_gene_freq}.")

    genes = sorted(genes, key=lambda g: (-freq[g], g))
    if max_genes is not None:
        genes = genes[: int(max_genes)]

    # keep only genes present
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) == 0:
        raise ValueError("None of the leading-edge genes are present in adata.var_names.")

    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    gidx = [adata.var_names.get_loc(g) for g in genes]
    M = X[:, gidx].toarray() if (sp is not None and sp.issparse(X)) else np.asarray(X[:, gidx], dtype=float)

    if use not in {"samples", "group_mean"}:
        raise ValueError("use must be 'samples' or 'group_mean'")

    if use == "group_mean":
        if groupby is None:
            raise ValueError("groupby must be provided when use='group_mean'")
        if groupby not in adata.obs.columns:
            raise KeyError(f"groupby='{groupby}' not found in adata.obs")

        s = adata.obs[groupby].astype("category")
        cats = list(s.cat.categories)

        out = np.zeros((len(cats), len(genes)), dtype=float)
        for i, c in enumerate(cats):
            mask = (s == c).to_numpy()
            out[i, :] = np.nan if mask.sum() == 0 else M[mask, :].mean(axis=0)

        df = pd.DataFrame(out, index=[str(c) for c in cats], columns=genes)
    else:
        df = pd.DataFrame(M, index=adata.obs_names.astype(str), columns=genes)

    if z_score == "row":
        arr = df.to_numpy(dtype=float)
        mu = np.nanmean(arr, axis=0, keepdims=True)
        sd = np.nanstd(arr, axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        arr = (arr - mu) / sd
        if clip_z is not None:
            arr = np.clip(arr, -float(clip_z), float(clip_z))
        df = pd.DataFrame(arr, index=df.index, columns=df.columns)
    elif z_score is None:
        pass
    else:
        raise ValueError("z_score must be 'row' or None")

    if figsize is None:
        w = max(7.0, 0.10 * df.shape[1] + 4.5)
        h = max(5.0, 0.12 * df.shape[0] + 3.0)
        figsize = (w, h)

    g = sns.clustermap(
        df,
        cmap=cmap,
        row_cluster=bool(row_cluster),
        col_cluster=bool(col_cluster),
        xticklabels=bool(show_labels),
        yticklabels=True,
        figsize=figsize,
        dendrogram_ratio=dendrogram_ratio,
        cbar_kws={"label": "Z-score" if z_score == "row" else "Expression"},
    )

    ax = g.ax_heatmap
    if show_title:
        ax.set_title(
            title or f"Leading-edge expression heatmap (terms={len(term_names)}, genes={df.shape[1]})",
            pad=10,
        )
    ax.set_xlabel("Leading-edge genes", fontsize=float(label_fontsize))
    ax.set_ylabel("Samples" if use == "samples" else str(groupby), fontsize=float(label_fontsize))

    for lab in ax.get_yticklabels():
        lab.set_fontsize(float(label_fontsize))

    if show_labels:
        _rotate_gene_labels(ax, fontsize=float(gene_label_fontsize))
        g.fig.subplots_adjust(bottom=0.30)

    if save is not None:
        _savefig(g.fig, save)
    if show:
        plt.show()

    return df, g


# ---------------------------------------------------------------------
# Plot 2: Jaccard similarity heatmap
# ---------------------------------------------------------------------
def leading_edge_jaccard_heatmap(
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
    min_shared_genes: int = 0,
    row_cluster: bool = True,
    col_cluster: bool = True,
    cmap: str = "viridis",
    vmin: float = 0.0,
    vmax: float = 1.0,
    figsize: tuple[float, float] | None = None,
    show_labels: bool = True,
    label_fontsize: float = 9.0,
    dendrogram_ratio: tuple[float, float] = (0.06, 0.12),
    title: str | None = None,
    show_title: bool = True,
    save: str | Path | None = None,
    show: bool = True,
):
    set_style()
    if sns is None:
        raise ImportError("leading_edge_jaccard_heatmap requires seaborn. Please install seaborn.")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)
    dfJ = _jaccard_matrix(term_names, le_sets)

    if int(min_shared_genes) > 0:
        for i, ti in enumerate(term_names):
            for j, tj in enumerate(term_names):
                if i == j:
                    continue
                if len(le_sets[ti] & le_sets[tj]) < int(min_shared_genes):
                    dfJ.iloc[i, j] = 0.0

    n = len(term_names)
    if figsize is None:
        s = max(6.0, 0.35 * n + 3.0)
        figsize = (s, s)

    g = sns.clustermap(
        dfJ,
        cmap=cmap,
        vmin=float(vmin),
        vmax=float(vmax),
        row_cluster=bool(row_cluster),
        col_cluster=bool(col_cluster),
        xticklabels=bool(show_labels),
        yticklabels=bool(show_labels),
        figsize=figsize,
        dendrogram_ratio=dendrogram_ratio,
    )

    ax = g.ax_heatmap
    if show_title:
        ax.set_title(title or "Leading-edge Jaccard similarity (pathway × pathway)", pad=10)

    if show_labels:
        for lab in ax.get_xticklabels():
            lab.set_rotation(90)
            lab.set_ha("center")
            lab.set_va("top")
            lab.set_fontsize(float(label_fontsize))
        for lab in ax.get_yticklabels():
            lab.set_fontsize(float(label_fontsize))
        g.fig.subplots_adjust(bottom=0.30)

    if save is not None:
        _savefig(g.fig, save)
    if show:
        plt.show()

    return dfJ, g


# ---------------------------------------------------------------------
# Plot 3: Pathway × gene overlap matrix (binary)
# ---------------------------------------------------------------------
def leading_edge_overlap_matrix(
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
    min_gene_freq: int = 2,
    sort_genes_by: str = "freq",  # "freq" | "alpha"
    row_cluster: bool = True,
    col_cluster: bool = False,
    cmap: str = "Greys",
    figsize: tuple[float, float] | None = None,
    show_gene_labels: bool = True,
    gene_label_fontsize: float = 8.0,
    show_term_labels: bool = True,
    label_fontsize: float = 9.0,
    dendrogram_ratio: tuple[float, float] = (0.06, 0.12),
    title: str | None = None,
    show_title: bool = True,
    grid_every: int = 5,
    grid_color: str = "0.90",
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Pathway × gene binary matrix for leading-edge membership.

    Enhancements:
      - grid lines every N rows/cols for readability
      - label_fontsize controls pathway labels too
      - dendrogram_ratio shrinks left dendrogram
      - when returning df, returns it *in the clustered order* (rows/cols)
    """
    set_style()
    if sns is None:
        raise ImportError("leading_edge_overlap_matrix requires seaborn. Please install seaborn.")
    if int(min_gene_freq) < 1:
        raise ValueError("min_gene_freq must be >= 1")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)

    all_genes = sorted(set().union(*le_sets.values()))
    if len(all_genes) == 0:
        raise ValueError("No genes in leading-edge sets.")

    mat = np.zeros((len(term_names), len(all_genes)), dtype=int)
    for i, t in enumerate(term_names):
        s = le_sets[t]
        for j, g in enumerate(all_genes):
            mat[i, j] = 1 if g in s else 0

    df = pd.DataFrame(mat, index=term_names, columns=all_genes)

    gene_freq = df.sum(axis=0).astype(int)
    keep = gene_freq[gene_freq >= int(min_gene_freq)].index.tolist()
    df = df.loc[:, keep]
    if df.shape[1] == 0:
        raise ValueError(f"No genes pass min_gene_freq={min_gene_freq}.")

    if sort_genes_by == "freq":
        df = df.loc[:, df.sum(axis=0).sort_values(ascending=False).index]
    elif sort_genes_by == "alpha":
        df = df.loc[:, sorted(df.columns)]
    else:
        raise ValueError("sort_genes_by must be 'freq' or 'alpha'")

    if figsize is None:
        w = max(7.0, 0.18 * df.shape[1] + 4.5)
        h = max(4.5, 0.22 * df.shape[0] + 2.5)
        figsize = (w, h)

    # ensure 0 is white even if user uses a different cmap
    cm = mpl.cm.get_cmap(cmap)
    if hasattr(cm, "with_extremes"):
        cm = cm.with_extremes(under="white")
    # but easiest: vmin=0 and df is 0/1
    g = sns.clustermap(
        df,
        cmap=cm,
        vmin=0,
        vmax=1,
        row_cluster=bool(row_cluster),
        col_cluster=bool(col_cluster),
        linewidths=0.0,
        xticklabels=bool(show_gene_labels),
        yticklabels=bool(show_term_labels),
        figsize=figsize,
        dendrogram_ratio=dendrogram_ratio,
        cbar_pos=None,
    )

    ax = g.ax_heatmap
    ax.set_xlabel(f"Leading-edge genes (kept={df.shape[1]}, freq≥{min_gene_freq})", fontsize=float(label_fontsize))
    ax.set_ylabel("Pathways", fontsize=float(label_fontsize))
    if show_title:
        ax.set_title(title or "Leading-edge overlap (pathway × gene)", pad=10)

    # grid every N
    if int(grid_every) > 0:
        nrows, ncols = df.shape
        for r in range(0, nrows + 1, int(grid_every)):
            ax.axhline(r, color=grid_color, lw=0.8, zorder=10)
        for c in range(0, ncols + 1, int(grid_every)):
            ax.axvline(c, color=grid_color, lw=0.8, zorder=10)

    # apply fonts
    for lab in ax.get_yticklabels():
        lab.set_fontsize(float(label_fontsize))
    if show_gene_labels:
        _rotate_gene_labels(ax, fontsize=float(gene_label_fontsize))
        g.fig.subplots_adjust(bottom=0.30)

    # IMPORTANT: return df in clustered order (rows + cols)
    row_order = g.dendrogram_row.reordered_ind if g.dendrogram_row is not None else list(range(df.shape[0]))
    col_order = g.dendrogram_col.reordered_ind if g.dendrogram_col is not None else list(range(df.shape[1]))
    df_ord = df.iloc[row_order, :].iloc[:, col_order].copy()

    if save is not None:
        _savefig(g.fig, save)
    if show:
        plt.show()

    return df_ord, g

# ---------------------------------------------------------------------
# Plot 4: Pathway cluster bubbles
# ---------------------------------------------------------------------
def leading_edge_cluster_bubbles(
    res: dict,
    *,
    min_overlap_genes: int = 1,
    layout: Literal["mds", "spring", "pcoa"] = "pcoa",
    figsize: tuple[float, float] = (8.0, 7.0),
    bubble_alpha: float = 0.35,
    edge_alpha: float = 0.40,
    edge_width_scale: float = 6.0,
    radius_scale: float = 0.10,
    radius_mode: Literal["sqrt", "log"] = "log",
    dist_scale: float = 4.0,
    normalize_xy: bool = True,          # NEW: can disable if you want raw embedding scale
    show_labels: bool = True,
    label_fontsize: float = 10.0,
    label_style: Literal["id", "size"] = "id",
    show_edges: bool = True,
    seed: int = 0,
    title: str | None = "Pathway clusters bubble map",
    repel: bool = True,
    repel_iters: int = 200,
    repel_strength: float = 0.8,
    repel_padding: float = 0.02,
    return_metrics: bool = False,        # NEW: optionally return embedding metrics as 8th object
):
    """
    Bubble map for pathway clusters.

    - Bubble size ~ #unique leading-edge genes in the cluster (union of member pathways).
    - Bubble distance ~ attempts to match D = 1 - Jaccard between cluster unions.
    - Optional edges connect clusters sharing >= min_overlap_genes genes (width ~ shared genes).
    - dist_scale spreads the embedding; radius_mode/radius_scale control bubble sizes.
    - repel performs a simple overlap-reduction pass (no extra deps).

    Returns
    -------
    (fig, ax, df_xy, dfJc, dfOc, cluster_gene_sets, cluster_terms)
    If return_metrics=True, adds an extra 8th return dict with:
      {"stress": ..., "dist_corr": ...}
    """
    # ---- extract ----
    if "clusters" not in res or "le_sets" not in res:
        raise ValueError("res must contain 'clusters' and 'le_sets' from leading_edge_pathway_clusters().")

    clusters: pd.Series = res["clusters"]
    le_sets: dict = res["le_sets"]

    # ---- build cluster unions: cluster_id -> set(genes) ----
    cluster_ids = sorted(pd.unique(clusters.astype(int)))
    cluster_gene_sets: dict[int, set[str]] = {}
    cluster_terms: dict[int, list[str]] = {}

    for term, cid in clusters.items():
        cid = int(cid)
        term = str(term)
        cluster_terms.setdefault(cid, []).append(term)
        cluster_gene_sets.setdefault(cid, set()).update(set(map(str, le_sets.get(term, set()))))

    cluster_ids = [cid for cid in cluster_ids if len(cluster_gene_sets.get(cid, set())) > 0]
    if len(cluster_ids) < 2:
        raise ValueError("Need at least 2 non-empty clusters for bubble plot.")

    genesets = [cluster_gene_sets[cid] for cid in cluster_ids]
    n = len(cluster_ids)

    # ---- overlap + cluster-level Jaccard ----
    O = np.zeros((n, n), dtype=int)
    J = np.zeros((n, n), dtype=float)

    for i in range(n):
        Ai = genesets[i]
        for j in range(n):
            Aj = genesets[j]
            inter = len(Ai & Aj)
            union = len(Ai | Aj)
            O[i, j] = inter
            J[i, j] = (inter / union) if union > 0 else 0.0

    D = 1.0 - J
    np.fill_diagonal(D, 0.0)

    # ---- helpers for embedding quality ----
    stress = np.nan
    dist_corr = np.nan

    def _pairwise_euclid(XY: np.ndarray) -> np.ndarray:
        # returns square distance matrix
        XY = np.asarray(XY, dtype=float)
        n0 = XY.shape[0]
        DD = np.zeros((n0, n0), dtype=float)
        for i in range(n0):
            for j in range(i + 1, n0):
                d = float(np.hypot(XY[i, 0] - XY[j, 0], XY[i, 1] - XY[j, 1]))
                DD[i, j] = d
                DD[j, i] = d
        return DD

    def _dist_corr(Dtrue: np.ndarray, Dplot: np.ndarray) -> float:
        m = np.triu_indices_from(Dtrue, k=1)
        a = Dtrue[m].astype(float)
        b = Dplot[m].astype(float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 3:
            return float("nan")
        a = a[ok]
        b = b[ok]
        if np.std(a) == 0 or np.std(b) == 0:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    # ---- layout ----
    dv = float(np.nanstd(D[np.isfinite(D)])) if np.isfinite(D).any() else 0.0
    if dv < 1e-8:
        # all distances ~ equal -> place on circle
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        XY = np.c_[np.cos(ang), np.sin(ang)]
    else:
        if layout == "pcoa":
            # Classical MDS / PCoA
            D2 = D**2
            H = np.eye(n) - np.ones((n, n)) / n
            B = -0.5 * H @ D2 @ H
            evals, evecs = np.linalg.eigh(B)
            order = np.argsort(evals)[::-1]
            evals = evals[order]
            evecs = evecs[:, order]
            evals = np.maximum(evals[:2], 0.0)
            XY = evecs[:, :2] * np.sqrt(evals)

        elif layout == "mds":
            if MDS is None:
                raise ImportError("layout='mds' requires scikit-learn. Install with: pip install scikit-learn")
            mds = MDS(
                n_components=2,
                dissimilarity="precomputed",
                metric=True,
                random_state=int(seed),
                n_init=20,
                max_iter=3000,
                eps=1e-12,
            )
            XY = mds.fit_transform(D)
            try:
                stress = float(mds.stress_)
            except Exception:
                stress = np.nan

        elif layout == "spring":
            try:
                import networkx as nx
            except Exception as e:
                raise ImportError(f"layout='spring' requires networkx. Install with: pip install networkx. ({e})")

            G = nx.Graph()
            for cid in cluster_ids:
                G.add_node(cid)
            for i in range(n):
                for j in range(i + 1, n):
                    w = float(J[i, j])
                    if w > 0:
                        G.add_edge(cluster_ids[i], cluster_ids[j], weight=w)

            pos = nx.spring_layout(G, seed=int(seed), weight="weight")
            XY = np.array([[pos[c][0], pos[c][1]] for c in cluster_ids], dtype=float)

        else:
            raise ValueError("layout must be 'pcoa' or 'mds' or 'spring'")

    # ---- center/scale coordinates ----
    x = XY[:, 0].astype(float)
    y = XY[:, 1].astype(float)

    # center always
    x -= np.nanmean(x)
    y -= np.nanmean(y)

    if normalize_xy:
        s = float(np.nanstd(np.r_[x, y]))
        if not np.isfinite(s) or s == 0:
            s = 1.0
        x /= s
        y /= s

    # expand distances
    x *= float(dist_scale)
    y *= float(dist_scale)

    # ---- radii ----
    n_genes = np.array([len(cluster_gene_sets[cid]) for cid in cluster_ids], dtype=float)
    if radius_mode == "sqrt":
        r = np.sqrt(np.maximum(n_genes, 1.0))
    elif radius_mode == "log":
        r = np.log1p(np.maximum(n_genes, 1.0))
    else:
        raise ValueError("radius_mode must be 'sqrt' or 'log'")
    r = float(radius_scale) * r

    # ---- repel (push circles apart) ----
    if repel and n >= 2:
        rng = np.random.default_rng(int(seed))
        x = x + rng.normal(0, 1e-4, size=n)
        y = y + rng.normal(0, 1e-4, size=n)

        for _ in range(int(repel_iters)):
            moved = False
            for i in range(n):
                for j in range(i + 1, n):
                    dx = x[j] - x[i]
                    dy = y[j] - y[i]
                    dist = float(np.hypot(dx, dy))
                    target = float(r[i] + r[j] + repel_padding)

                    if dist < 1e-12:
                        ang = float(rng.uniform(0, 2 * np.pi))
                        dx, dy = np.cos(ang), np.sin(ang)
                        dist = 1e-6

                    if dist < target:
                        push = (target - dist) / dist
                        step = float(repel_strength) * push * 0.5
                        x[i] -= dx * step
                        y[i] -= dy * step
                        x[j] += dx * step
                        y[j] += dy * step
                        moved = True
            if not moved:
                break

    # ---- compute fidelity metric ----
    Dplot = _pairwise_euclid(np.c_[x, y])
    dist_corr = _dist_corr(D, Dplot)

    # ---- plot ----
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal", adjustable="datalim")

    if show_edges:
        maxO = max(1, int(O.max()))
        for i in range(n):
            for j in range(i + 1, n):
                shared = int(O[i, j])
                if shared < int(min_overlap_genes):
                    continue
                lw = float(edge_width_scale) * (shared / maxO)
                ax.plot(
                    [x[i], x[j]], [y[i], y[j]],
                    linewidth=lw, alpha=float(edge_alpha),
                    color="0.40", zorder=1
                )

    for i, cid in enumerate(cluster_ids):
        circ = plt.Circle(
            (x[i], y[i]), float(r[i]),
            alpha=float(bubble_alpha),
            ec="0.25", lw=1.0, zorder=2
        )
        ax.add_patch(circ)

    if show_labels:
        for i, cid in enumerate(cluster_ids):
            if label_style == "size":
                lab = f"C{cid}\n{int(n_genes[i])} genes"
            else:
                lab = f"C{cid}"
            ax.text(
                x[i], y[i], lab,
                ha="center", va="center",
                fontsize=float(label_fontsize),
                zorder=3
            )

    pad = float(np.max(r)) * 1.6 if np.isfinite(r).any() else 0.5
    ax.set_xlim(float(np.min(x) - pad), float(np.max(x) + pad))
    ax.set_ylim(float(np.min(y) - pad), float(np.max(y) + pad))
    ax.set_xticks([])
    ax.set_yticks([])

    if title:
        # attach fidelity info lightly (optional)
        extra = ""
        if np.isfinite(dist_corr):
            extra += f"  (dist corr={dist_corr:.2f})"
        ax.set_title(str(title) + extra, pad=12)

    # ---- tables ----
    df_xy = pd.DataFrame(
        {
            "cluster": cluster_ids,
            "x": x,
            "y": y,
            "n_le_genes": n_genes,
            "radius": r,
            "n_terms": [len(cluster_terms[cid]) for cid in cluster_ids],
        }
    ).set_index("cluster")

    dfJc = pd.DataFrame(J, index=cluster_ids, columns=cluster_ids)
    dfOc = pd.DataFrame(O, index=cluster_ids, columns=cluster_ids)

    if return_metrics:
        metrics = {"stress": float(stress) if np.isfinite(stress) else np.nan,
                   "dist_corr": float(dist_corr) if np.isfinite(dist_corr) else np.nan}
        return fig, ax, df_xy, dfJc, dfOc, cluster_gene_sets, cluster_terms, metrics

    return fig, ax, df_xy, dfJc, dfOc, cluster_gene_sets, cluster_terms

# ---------------------------------------------------------------------
# Export for Cytoscape
# ---------------------------------------------------------------------
def export_leading_edge_clusters_cytoscape(
    res: dict,
    *,
    min_overlap_genes: int = 1,
    min_jaccard: float | None = None,
    include_genes: bool = True,
) -> dict[str, pd.DataFrame]:
    """
    Export pathway clusters and their overlaps for Cytoscape.

    Returns a dict of DataFrames:
      - nodes_clusters: cluster node table
      - edges_clusters: cluster-cluster edge table (overlap + jaccard)
      - (optional) nodes_genes: gene node table
      - (optional) edges_cluster_gene: bipartite edges cluster->gene

    Notes
    -----
    - Nodes are labeled "C{cluster_id}" for clusters, and gene symbols for genes.
    """
    if "clusters" not in res or "le_sets" not in res:
        raise ValueError("res must contain 'clusters' and 'le_sets'.")

    clusters: pd.Series = res["clusters"].astype(int)
    le_sets: dict = res["le_sets"]

    # build unions
    cluster_ids = sorted(pd.unique(clusters))
    cluster_gene_sets: dict[int, set[str]] = {int(cid): set() for cid in cluster_ids}
    cluster_terms: dict[int, list[str]] = {int(cid): [] for cid in cluster_ids}

    for term, cid in clusters.items():
        cid = int(cid)
        cluster_terms[cid].append(str(term))
        cluster_gene_sets[cid].update(set(map(str, le_sets.get(term, set()))))

    # drop empties
    cluster_ids = [cid for cid in cluster_ids if len(cluster_gene_sets[cid]) > 0]
    n = len(cluster_ids)
    if n == 0:
        raise ValueError("No non-empty clusters found.")

    # nodes: clusters
    nodes_clusters = pd.DataFrame(
        {
            "id": [f"C{cid}" for cid in cluster_ids],
            "cluster": cluster_ids,
            "n_terms": [len(cluster_terms[cid]) for cid in cluster_ids],
            "n_le_genes": [len(cluster_gene_sets[cid]) for cid in cluster_ids],
            "terms": [";".join(cluster_terms[cid]) for cid in cluster_ids],
        }
    )

    # edges: cluster-cluster
    rows = []
    for i in range(n):
        ci = cluster_ids[i]
        Ai = cluster_gene_sets[ci]
        for j in range(i + 1, n):
            cj = cluster_ids[j]
            Aj = cluster_gene_sets[cj]
            inter = len(Ai & Aj)
            union = len(Ai | Aj)
            jac = (inter / union) if union > 0 else 0.0

            if inter < int(min_overlap_genes):
                continue
            if (min_jaccard is not None) and (jac < float(min_jaccard)):
                continue

            rows.append(
                {
                    "source": f"C{ci}",
                    "target": f"C{cj}",
                    "shared_genes": int(inter),
                    "jaccard": float(jac),
                    "union_genes": int(union),
                }
            )

    edges_clusters = pd.DataFrame(rows)

    out = {"nodes_clusters": nodes_clusters, "edges_clusters": edges_clusters}

    # optional bipartite network cluster-gene
    if include_genes:
        all_genes = sorted(set().union(*[cluster_gene_sets[c] for c in cluster_ids]))
        nodes_genes = pd.DataFrame({"id": all_genes, "type": "gene"})
        nodes_clusters2 = nodes_clusters.copy()
        nodes_clusters2["type"] = "cluster"

        # edges cluster->gene
        e2 = []
        for cid in cluster_ids:
            for g in cluster_gene_sets[cid]:
                e2.append({"source": f"C{cid}", "target": str(g), "type": "cluster_gene"})
        edges_cluster_gene = pd.DataFrame(e2)

        out["nodes_all"] = pd.concat([nodes_clusters2[["id", "type"]], nodes_genes], axis=0).reset_index(drop=True)
        out["edges_cluster_gene"] = edges_cluster_gene

    return out