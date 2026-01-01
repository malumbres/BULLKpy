# Leiden resolution scan

```{eval-rst}
.. autofunction:: bullkpy.tl.leiden_resolution_scan

```

Scan Leiden clustering resolutions and score against a ground-truth annotation.

This utility evaluates **multiple Leiden resolutions** and quantifies how well
the resulting clusters match a known categorical label
(e.g. cell type, subtype, batch) using **external validation metrics**.

## What it does

For each Leiden resolution, the function:  

1. Computes (or reuses) a kNN graph.  

2. Runs Leiden clustering.  

3. Compares cluster labels to a reference annotation.  

4. Computes agreement metrics:   
	•	ARI (Adjusted Rand Index)
	•	NMI (Normalized Mutual Information)
	•	Cramér’s V (always available). 

5. Returns a tidy table summarizing clustering performance.  

Optionally stores results in adata.uns.

## When to use

Use leiden_resolution_scan when you want to:
	•	Choose an optimal Leiden resolution
	•	Benchmark clustering against known labels
	•	Detect over- or under-clustering
	•	Quantify clustering robustness across resolutions

Typical use cases:
	•	Validate clustering against known subtypes
	•	Tune resolution before downstream analysis
	•	Compare embeddings or representations

## Parameters

**adata**   
AnnData object with an embedding or representation

**true_key**   
Column in adata.obs containing ground-truth labels  
(e.g. "Subtype", "Batch", "CellType").  

**resolutions**    
Sequence of Leiden resolutions to scan
	•	Default: np.linspace(0.2, 2.0, 10)

**base_key**    
Base name for clustering columns in adata.obs.    
Each resolution creates a column like "{base_key}_{resolution}".  

**store_key**
Key under adata.uns to store results (if inplace=True)

**use_rep**   
Representation used for neighbor graph construction    
(e.g. "X_pca")

**n_pcs**   
Number of principal components to use (if applicable)

**n_neighbors**   
Number of neighbors for kNN graph

**metric**   
Distance metric for neighbor search

**recompute_neighbors**   
Force recomputation of the neighbors graph

**inplace**   
If True, store results in adata.uns[store_key]

**random_state**.   
Random seed for Leiden clustering

## Metrics explained

**Adjusted Rand Index (ARI)** 
	•	Measures similarity between two clusterings
	•	Adjusted for chance
	•	Range: –1 to 1
	•	Higher = better agreement

**Normalized Mutual Information (NMI)**
	•	Information-theoretic similarity
	•	Range: 0 to 1
	•	Higher = better agreement
	•	Requires scikit-learn

**Cramér’s V**
	•	Association measure for categorical variables
	•	Range: 0 to 1
	•	Always computed (no external dependencies)

## Output

Returns a tidy DataFrame with one row per resolution:

| Column | Description |
| ---------- | -------------------- |
| resolution | Leiden resolution |
| key | obs column storing cluster labels |
| n_clusters | Number of clusters |
| ARI | Adjusted Rand Index |
| NMI | Normalized Mutual Information |
| cramers_v | Cramér’s V |
| n_used | Samples used in scoring |
| n_missing | Samples excluded (NA labels) |

Sorted by increasing resolution.

## Examples

**Basic resolution scan**   

```python
df = bk.tl.leiden_resolution_scan(
    adata,
    true_key="Subtype",
)
df
```

**Custom resolution grid**    

```python
df = bk.tl.leiden_resolution_scan(
    adata,
    true_key="CellType",
    resolutions=[0.2, 0.4, 0.6, 0.8, 1.0],
)
```

**Plot ARI vs resolution**   

```python
import matplotlib.pyplot as plt

plt.plot(df["resolution"], df["ARI"], marker="o")
plt.xlabel("Leiden resolution")
plt.ylabel("ARI")
plt.show()
```

## Notes
- The neighbors graph is computed once unless recompute_neighbors=True
- Missing labels in true_key are automatically excluded from scoring
- If scikit-learn is unavailable:
	•	ARI and NMI are returned as NaN
	•	Cramér’s V is still computed
- Clustering columns (leiden_<resolution>) remain in adata.obs

## See also
	•	tl.neighbors
	•	tl.cluster
	•	pl.rankplot
	•	pl.violin
	•	scanpy.tl.leiden