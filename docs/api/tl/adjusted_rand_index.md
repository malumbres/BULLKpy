# Adjusted rand index

```{eval-rst}
.. autofunction:: bullkpy.tl.adjusted_rand_index

```

Compute the **Adjusted Rand Index (ARI)** between two categorical annotations
stored in `adata.obs`.

This is a lightweight convenience wrapper around
`sklearn.metrics.adjusted_rand_score` that handles missing values in AnnData.


## What it does
- Compares two categorical labelings (true_key vs pred_key)
- Automatically filters out samples with missing values in either column
- Returns a single float ARI score

ARI measures similarity between two partitions, corrected for chance:
- 1.0 → perfect agreement
- 0.0 → random agreement
- < 0.0 → worse than random

## Parameters

**adata**   
AnnData object containing both label vectors in adata.obs.

**true_key**  
Column in adata.obs with the reference / ground-truth labels
(e.g. "cell_type").

**pred_key**   
Column in adata.obs with predicted labels
(e.g. "leiden", "kmeans", "predicted_type").

## Returns

float. 
Adjusted Rand Index computed on samples where both labels are present.

**Raises**:

ValueError.  
If no samples have non-missing values in both true_key and pred_key.

## Examples

Basic usage
```python
ari = bk.tl.adjusted_rand_index(
    adata,
    true_key="cell_type",
    pred_key="leiden",
)
print(f"ARI = {ari:.3f}")
```

Compare multiple clusterings
```python
for key in ["leiden_0.5", "leiden_1.0", "leiden_2.0"]:
    ari = bk.tl.adjusted_rand_index(
        adata,
        true_key="cell_type",
        pred_key=key,
    )
    print(key, ari)
```

## Notes
- This function only computes ARI.  
For a richer evaluation (NMI, Cramér’s V, silhouette score), see:
	•	tl.cluster_metrics
	•	Requires scikit-learn.

## See also
	•	tl.cluster_metrics – ARI, NMI, Cramér’s V, silhouette in one call
	•	pl.ari_resolution_heatmap – visualize ARI across Leiden resolutions