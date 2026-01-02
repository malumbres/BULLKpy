# Cluster metrics

```{eval-rst}
.. autofunction:: bullkpy.tl.cluster_metrics

```

Compute **agreement metrics** between a clustering assignment and a known
ground-truth label, plus a **silhouette score** to summarize cluster separation.

This is useful to evaluate clustering solutions (e.g. Leiden) against curated
annotations, and to compare alternative resolutions / representations.

## What it computes

Returns a dictionary with:
	•	n_used: number of samples used after filtering (if dropna=True)
	•	ari: Adjusted Rand Index between true_key and cluster_key
	•	nmi: Normalized Mutual Information between true_key and cluster_key
	•	cramers_v: Cramér’s V association between true_key and cluster_key
	•	silhouette: Silhouette score measuring how well clusters separate in the
chosen feature space (rep, X, or layer). 

Agreement metrics quantify how similar the clustering labels are to the known
labels; silhouette measures intrinsic cluster compactness/separation.

## Parameters

#### Required

**adata**  
AnnData object containing both ground-truth and clustering labels in adata.obs.

**true_key**   
Column name in adata.obs with the reference/known labels (e.g. "cell_type").

#### Optional: clustering labels

**cluster_key** (default: "leiden")    
Column name in adata.obs with cluster assignments to evaluate.

#### Optional: silhouette configuration

Silhouette can be computed on different data representations:  

**silhouette_on** (default: "rep")   
Which data to use for silhouette:
- "rep": use adata.obsm[use_rep] (recommended; e.g. PCA space)
- "X": use adata.X
- "layer": use adata.layers[layer]

**use_rep** (default: "X_pca")    
Key in adata.obsm used when silhouette_on="rep".

**layer** (default: None)    
Layer name used when silhouette_on="layer".

**n_pcs** (default: None)    
If provided and silhouette_on="rep", only the first n_pcs columns of
adata.obsm[use_rep] are used.

**metric** (default: "euclidean")    
Distance metric passed to the silhouette score (scikit-learn). Common choices:
"euclidean", "cosine".

#### Optional: missing values

**dropna** (default: True)   
If True, samples with missing true_key or cluster_key are removed before
computing metrics. If False, missing values will raise an error.

## Returns

**dict[str, float]**
Dictionary with keys: n_used, ari, nmi, cramers_v, silhouette.

Notes:  
- If the number of clusters is <2 or >= n_used, silhouette is not defined and
the function returns silhouette = NaN (and emits a warning).
- Requires scikit-learn for ARI/NMI and silhouette.

## Examples

1) Evaluate Leiden vs known labels (PCA silhouette)
```python
m = bk.tl.cluster_metrics(
    adata,
    true_key="cell_type",
    cluster_key="leiden",
    silhouette_on="rep",
    use_rep="X_pca",
    n_pcs=20,
)
m
```

2) Silhouette on expression layer (e.g. log1p CPM)
```python
m = bk.tl.cluster_metrics(
    adata,
    true_key="cell_type",
    cluster_key="leiden",
    silhouette_on="layer",
    layer="log1p_cpm",
)
M
```

3) Compare multiple clusterings
```python
for key in ["leiden_0.5", "leiden_1.0", "leiden_2.0"]:
    m = bk.tl.cluster_metrics(adata, true_key="cell_type", cluster_key=key)
    print(key, m["ari"], m["silhouette"])
```

## Interpretation tips
- ARI / NMI: higher is better (1.0 = perfect agreement).
- Cramér’s V: association strength between categorical variables (0..1).
- Silhouette: higher means better separation (typically between -1 and 1).
Use silhouette mainly to compare representations or resolution choices on
the same dataset.

## See also
	•	tl.leiden_resolution_scan – scan Leiden resolutions and score vs ground truth
	•	pl.ari_resolution_heatmap – visualize resolution scan metrics