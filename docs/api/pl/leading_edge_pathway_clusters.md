# Leading-edge pathway clusters

```{eval-rst}
.. autofunction:: bullkpy.pl.leading_edge_pathway_clusters

```

Groups GSEA pathways into clusters ("nodules") based on the overlap of their
leading-edge gene sets, using hierarchical clustering on pairwise Jaccard similarity.

This is a **data transformation function** (not a plot): it returns a dict that
feeds downstream visualizations such as `pl.leading_edge_cluster_bubbles` or
`pl.leading_edge_cluster_driver_genes`.

## What it does

1. Extracts the leading-edge gene sets for the selected pathways.
2. Builds a pairwise Jaccard similarity matrix across pathways.
3. Optionally zeroes out pairs that share fewer than `min_shared_genes` genes (de-noising).
4. Performs hierarchical clustering (scipy `linkage`) on the distance matrix (`1 − Jaccard`).
5. Cuts the dendrogram either by a distance `threshold` or by a fixed `n_clusters`.
6. Returns cluster assignments and per-cluster cohesion metrics.

Pathways in the same cluster share a large fraction of their leading-edge genes,
suggesting they are driven by the same biological mechanism.

## Expected input

`pre_res` is the result object returned by `bk.tl.gsea_preranked()` (a `gseapy.GSEA`
or `gseapy.Prerank` object). It must contain leading-edge gene information.

Select pathways with **either** `term_idx` (integer positions) **or** `terms` (pathway
names) — not both.

## Parameters

**pre_res**
GSEA result object (from `bk.tl.gsea_preranked`). Must contain leading-edge gene sets.

**term_idx**
List of integer indices selecting which pathways to include.
Mutually exclusive with `terms`.

**terms**
List of pathway name strings to include.
Mutually exclusive with `term_idx`.

**method**
Linkage method passed to `scipy.cluster.hierarchy.linkage`.
Common options: `"average"` (default), `"complete"`, `"ward"`.

**threshold**
Distance threshold for cutting the dendrogram (`distance = 1 − Jaccard`).
Smaller values → more clusters. Mutually exclusive with `n_clusters`.

**n_clusters**
Fixed number of clusters to request.
Mutually exclusive with `threshold`.

**min_shared_genes**
Minimum number of genes two pathways must share for their Jaccard similarity to be
retained; pairs below this are zeroed out before clustering. Default: `0` (no filtering).

## Returns

```python
{
    "term_names": list[str],        # ordered pathway names
    "le_sets":    dict[str, set],   # pathway -> leading-edge gene set
    "dfJ":        pd.DataFrame,     # Jaccard similarity matrix (pathways × pathways)
    "clusters":   pd.Series,        # pathway -> cluster_id (int)
    "metrics":    pd.DataFrame,     # per-cluster cohesion table
    "linkage":    np.ndarray,       # scipy linkage matrix Z
}
```

## Examples

1) Cluster by distance threshold
```python
pre_res = bk.tl.gsea_preranked(adata, gene_sets="MSigDB_Hallmark_2020")

result = bk.pl.leading_edge_pathway_clusters(
    pre_res,
    threshold=0.6,   # distance = 1 - Jaccard; lower -> tighter clusters
)

print(result["clusters"].value_counts())
```

2) Request a fixed number of clusters
```python
result = bk.pl.leading_edge_pathway_clusters(
    pre_res,
    n_clusters=5,
)
```

3) Restrict to a subset of pathways and apply de-noising
```python
result = bk.pl.leading_edge_pathway_clusters(
    pre_res,
    terms=["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNF_SIGNALING_VIA_NFKB",
           "HALLMARK_IL6_JAK_STAT3_SIGNALING"],
    threshold=0.5,
    min_shared_genes=3,
)
```

4) Pass result to downstream visualizations
```python
result = bk.pl.leading_edge_pathway_clusters(pre_res, n_clusters=6)

bk.pl.leading_edge_cluster_bubbles(pre_res, result)
```

## Notes
- Exactly one of `threshold` or `n_clusters` must be provided; passing both (or neither) raises a `ValueError`.
- Requires `scipy` (`cluster.hierarchy` and `spatial.distance`).
- The returned `linkage` matrix can be used directly with `scipy.cluster.hierarchy.dendrogram` for custom visualizations.
- Jaccard similarity of a pathway with itself is 1 (diagonal), but clustering uses only off-diagonal distances.

## See also
- `tl.gsea_preranked`
- `pl.leading_edge_cluster_driver_genes`
- `pl.leading_edge_cluster_bubbles`
- `pl.leading_edge_jaccard_heatmap`
- `pl.gsea_leading_edge_heatmap`
