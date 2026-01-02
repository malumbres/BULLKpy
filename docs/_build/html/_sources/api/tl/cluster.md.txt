# Clustering

```{eval-rst}
.. autofunction:: bullkpy.tl.cluster

```

Cluster samples stored in an `AnnData` object using either a **graph-based Leiden algorithm** (recommended) or a **K-means fallback** in PCA space.

This function mirrors Scanpy’s clustering workflow while remaining **explicit, bulk-friendly, and dependency-aware**.

## Overview

cluster assigns each sample (adata.obs) to a discrete cluster and stores the result as a categorical column in adata.obs.  

Two methods are supported:  
- Leiden (default, recommended)
	- Graph-based clustering
	- Operates on adata.obsp["connectivities"]
	- Produces resolution-controlled communities
- K-means (fallback)
	- Runs directly in PCA space
	- Requires no extra dependencies
	- Useful when graph-based clustering is unavailable

## Requirements

**For method="leiden"**    

You must have:  
1. A neighbors graph:  
```python
bk.tl.neighbors(adata)
```

This creates adata.obsp["connectivities"].  

2. Optional but recommended:  
- PCA via bk.tl.pca(adata). 

3. Python packages:  
- igraph
- leidenalg.  

If any of these are missing, an informative error is raised.  

**For method="kmeans"**   

You must have:
- A representation in adata.obsm[use_rep] (default: "X_pca"). 

No additional dependencies are required.

## Parameters

#### Common parameters

**method**   
Clustering algorithm to use:  
- "leiden" (default)
- "kmeans".  

**key_added**   
Name of the column written to adata.obs.  
The column is stored as a categorical variable.  
Default: "clusters".  

**random_state**   
Random seed for reproducibility.  
Default: 0.  

#### Leiden-specific parameters

**resolution**.  
Controls cluster granularity:  
- Higher → more clusters. 
- Lower → fewer clusters. 
Typical range: 0.2–2.0. 
Default: 1.0.  

#### K-means–specific parameters

**n_clusters**
Number of clusters (k).  
Default: 8.  

**use_rep**   
Representation in adata.obsm used for clustering.  
Default: "X_pca".  

**n_pcs**  
Number of dimensions to use from use_rep.  
Default: 20.  

## Method details

**Leiden clustering**  
- Converts the connectivities graph into an undirected weighted igraph object
- Runs Leiden using RBConfigurationVertexPartition
- Cluster labels are integers, stored as strings for categorical safety
- Results are deterministic given random_state. 

This is the preferred method for most analyses.  

**K-means clustering**   
- Runs standard K-means in PCA space
- Uses only the first n_pcs dimensions
- Faster and dependency-free, but:
- Ignores graph structure
- Less robust for complex manifolds

## Stored results

**In adata.obs**.   
- adata.obs[key_added]
	- Categorical cluster labels
	- Labels are strings ("0", "1", …). 

**In adata.uns["clusters"]**   

Each clustering run is recorded with its parameters:  

Leiden example
```python
adata.uns["clusters"]["clusters"] = {
    "method": "leiden",
    "resolution": 1.0,
    "random_state": 0,
}
```

K-means example
```python
adata.uns["clusters"]["clusters"] = {
    "method": "kmeans",
    "n_clusters": 8,
    "random_state": 0,
    "use_rep": "X_pca",
    "n_pcs": 20,
}
```

## Typical workflow

```python
# 1. Dimensionality reduction
bk.tl.pca(adata)

# 2. Build neighbor graph
bk.tl.neighbors(adata)

# 3. Graph-based clustering
bk.tl.cluster(
    adata,
    method="leiden",
    resolution=1.0,
    key_added="leiden",
)
```

Fallback without graph:

```python
bk.tl.cluster(
    adata,
    method="kmeans",
    n_clusters=10,
    key_added="kmeans",
)
```

## Notes and caveats
- Leiden clustering requires a neighbors graph.
- Re-running cluster with the same key_added will overwrite the column.
- Cluster labels are stored as categorical strings, not integers.
- For resolution tuning, see:
	- bk.tl.leiden_resolution_scan
	- bk.pl.ari_resolution_heatmap

## See also
	•	bk.tl.neighbors – build kNN graph
	•	bk.tl.pca – PCA representation
	•	bk.tl.leiden_resolution_scan – resolution benchmarking
	•	bk.tl.cluster_metrics – clustering quality metrics
