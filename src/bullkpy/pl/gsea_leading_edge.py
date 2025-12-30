from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import warn
from ._style import set_style, _savefig

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None


# ----------------------------
# Helpers
# ----------------------------
def _get_res2d(pre_res) -> pd.DataFrame:
    """
    GSEApy returns an object with .res2d (DataFrame).
    We keep this as a helper to be resilient to different versions.
    """
    if pre_res is None or not hasattr(pre_res, "res2d"):
        raise ValueError("pre_res must be a gseapy.prerank result object with attribute .res2d")
    df = pre_res.res2d
    if not isinstance(df, pd.DataFrame):
        raise TypeError("pre_res.res2d must be a pandas DataFrame")
    return df


def _find_col(df: pd.DataFrame, candidates: Sequence[str]) -> str:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        c = cols_lower.get(cand.lower())
        if c is not None:
            return c
    raise KeyError(f"None of the expected columns found: {candidates}. Available: {list(df.columns)}")


def _normalize_term_selector(df: pd.DataFrame, term_idx=None, terms: Sequence[str] | None = None) -> pd.DataFrame:
    if terms is not None:
        terms_set = set(map(str, terms))
        term_col = _find_col(df, ["Term", "term", "pathway", "Pathway", "NAME", "name"])
        out = df[df[term_col].astype(str).isin(terms_set)].copy()
        if out.empty:
            raise ValueError("No terms matched `terms=`.")
        return out

    if term_idx is None:
        return df.copy()

    if isinstance(term_idx, (slice, list, tuple, np.ndarray)):
        return df.iloc[term_idx].copy()

    if isinstance(term_idx, int):
        return df.iloc[[term_idx]].copy()

    raise TypeError("term_idx must be None, int, slice, or list-like. Or use terms=[...].")


def _split_leading_edge(s: str) -> list[str]:
    """
    Split GSEApy Lead_genes string into list of genes.
    Most common delimiter is ';', but we also support ',' and whitespace.
    """
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    s = str(s).strip()
    if not s:
        return []
    parts = re.split(r"[;,\s]+", s)
    parts = [p.strip() for p in parts if p.strip()]
    return parts


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
    le_col = _find_col(df, ["Lead_genes", "lead_genes", "ledge_genes", "ledge"])

    term_names: list[str] = []
    le_sets: dict[str, set[str]] = {}

    for _, row in df.iterrows():
        t = str(row[term_col])
        genes = _split_leading_edge(row[le_col])
        term_names.append(t)
        le_sets[t] = set(map(str, genes))

    # drop empties (sometimes GSEApy can have missing lead genes)
    term_names = [t for t in term_names if len(le_sets.get(t, set())) > 0]
    le_sets = {t: le_sets[t] for t in term_names}

    if len(term_names) == 0:
        raise ValueError("No leading-edge genes found for the selected terms/indices.")

    return term_names, le_sets


def _rotate_gene_labels(ax, fontsize: float = 8.0) -> None:
    # robust label formatting (works with seaborn clustermap axes)
    for lab in ax.get_xticklabels():
        lab.set_rotation(90)
        lab.set_ha("center")
        lab.set_va("top")
        lab.set_fontsize(fontsize)


# ----------------------------
# 1) Overlap matrix
# ----------------------------
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
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Pathway × gene binary matrix for leading-edge membership.

    This is the closest to the "Leading Edge Viewer" idea:
    you can see if a small set of genes drives many enriched pathways.
    """
    set_style()
    if sns is None:
        raise ImportError("leading_edge_overlap_matrix requires seaborn. Please install seaborn.")

    if int(min_gene_freq) < 1:
        raise ValueError("min_gene_freq must be >= 1")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)

    # union genes
    all_genes = sorted(set().union(*le_sets.values()))
    if len(all_genes) == 0:
        raise ValueError("No genes in leading-edge sets.")

    mat = np.zeros((len(term_names), len(all_genes)), dtype=int)
    for i, t in enumerate(term_names):
        s = le_sets[t]
        for j, g in enumerate(all_genes):
            mat[i, j] = 1 if g in s else 0

    df = pd.DataFrame(mat, index=term_names, columns=all_genes)

    # filter genes by frequency
    gene_freq = df.sum(axis=0).astype(int)
    keep = gene_freq[gene_freq >= int(min_gene_freq)].index.tolist()
    df = df.loc[:, keep]

    if df.shape[1] == 0:
        raise ValueError(
            f"No genes pass min_gene_freq={min_gene_freq}. "
            f"Try lowering it (e.g. 1) or selecting fewer terms."
        )

    # sort genes
    if sort_genes_by == "freq":
        df = df.loc[:, df.sum(axis=0).sort_values(ascending=False).index]
    elif sort_genes_by == "alpha":
        df = df.loc[:, sorted(df.columns)]
    else:
        raise ValueError("sort_genes_by must be 'freq' or 'alpha'")

    # auto figsize
    if figsize is None:
        w = max(7.0, 0.18 * df.shape[1] + 4.5)
        h = max(4.5, 0.22 * df.shape[0] + 2.5)
        figsize = (w, h)

    g = sns.clustermap(
        df,
        cmap=cmap,
        row_cluster=bool(row_cluster),
        col_cluster=bool(col_cluster),
        linewidths=0.0,
        xticklabels=bool(show_gene_labels),
        yticklabels=bool(show_term_labels),
        figsize=figsize,
        cbar_pos=None,
    )

    ax = g.ax_heatmap
    ax.set_xlabel(f"Leading-edge genes (kept={df.shape[1]}, freq≥{min_gene_freq})")
    ax.set_ylabel("Pathways")
    ax.set_title("Leading-edge overlap (pathway × gene)", pad=10)

    if show_gene_labels:
        _rotate_gene_labels(ax, fontsize=float(gene_label_fontsize))
        # Ensure enough room for rotated labels
        g.fig.subplots_adjust(bottom=0.30)

    if save is not None:
        _savefig(g.fig, save)

    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return df, g


# ----------------------------
# 2) Jaccard heatmap (term × term)
# ----------------------------
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
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Jaccard similarity between leading-edge gene sets for each pathway.

    Jaccard(A,B) = |A ∩ B| / |A ∪ B|

    Useful to see clusters of pathways driven by similar leading-edge genes.
    """
    set_style()
    if sns is None:
        raise ImportError("leading_edge_jaccard_heatmap requires seaborn. Please install seaborn.")

    term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)

    n = len(term_names)
    J = np.zeros((n, n), dtype=float)

    for i, ti in enumerate(term_names):
        Ai = le_sets[ti]
        for j, tj in enumerate(term_names):
            Aj = le_sets[tj]
            inter = len(Ai & Aj)
            union = len(Ai | Aj)
            J[i, j] = (inter / union) if union > 0 else 0.0

    dfJ = pd.DataFrame(J, index=term_names, columns=term_names)

    if int(min_shared_genes) > 0:
        # zero out similarities if intersection too small (visual de-noising)
        for i, ti in enumerate(term_names):
            for j, tj in enumerate(term_names):
                if i == j:
                    continue
                if len(le_sets[ti] & le_sets[tj]) < int(min_shared_genes):
                    dfJ.iloc[i, j] = 0.0

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
    )

    ax = g.ax_heatmap
    ax.set_title("Leading-edge Jaccard similarity (pathway × pathway)", pad=10)

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
        import matplotlib.pyplot as plt

        plt.show()

    return dfJ, g


# ----------------------------
# 3) Expression heatmap of leading-edge genes
# ----------------------------
def gsea_leading_edge_heatmap(
    adata: ad.AnnData,
    pre_res,
    *,
    term_idx=None,
    terms: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    use: str = "samples",  # "samples" | "group_mean"
    groupby: str | None = None,  # required if use="group_mean"
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
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap of expression for leading-edge genes from selected pathways.

    - Builds a union of leading-edge genes across selected terms
    - Optionally keeps only genes that appear in >= min_gene_freq pathways
    - Plots expression across:
        * samples (use="samples"), OR
        * group means (use="group_mean", requires groupby)

    This helps you see whether the shared leading-edge genes form coherent
    expression patterns across samples/subtypes.
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
        raise ValueError(
            f"No genes pass min_gene_freq={min_gene_freq}. "
            f"Try lowering it or selecting different terms."
        )

    # prioritize by frequency (and cap)
    genes = sorted(genes, key=lambda g: (-freq[g], g))
    if max_genes is not None:
        genes = genes[: int(max_genes)]

    # keep only genes present in adata
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) == 0:
        raise ValueError("None of the leading-edge genes are present in adata.var_names.")

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    gidx = [adata.var_names.get_loc(g) for g in genes]
    M = X[:, gidx].toarray() if sp.issparse(X) else np.asarray(X[:, gidx], dtype=float)

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
            if mask.sum() == 0:
                out[i, :] = np.nan
            else:
                out[i, :] = M[mask, :].mean(axis=0)

        df = pd.DataFrame(out, index=[str(c) for c in cats], columns=genes)
    else:
        df = pd.DataFrame(M, index=adata.obs_names.astype(str), columns=genes)

    # z-score per gene (row = gene in clustermap is axis=0, but here genes are columns)
    # We'll z-score columns so each gene has mean 0/var 1 across samples/groups.
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
        cbar_kws={"label": "Z-score" if z_score == "row" else "Expression"},
    )

    ax = g.ax_heatmap
    ax.set_title(
        f"Leading-edge expression heatmap (terms={len(term_names)}, genes={df.shape[1]})",
        pad=10,
    )
    ax.set_xlabel("Leading-edge genes")
    ax.set_ylabel("Samples" if use == "samples" else str(groupby))

    if show_labels:
        _rotate_gene_labels(ax, fontsize=float(gene_label_fontsize))
        g.fig.subplots_adjust(bottom=0.30)

    if save is not None:
        _savefig(g.fig, save)

    if show:
        import matplotlib.pyplot as plt

        plt.show()

    return df, g