from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


try:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
except Exception:  # pragma: no cover
    linkage = fcluster = squareform = None


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
