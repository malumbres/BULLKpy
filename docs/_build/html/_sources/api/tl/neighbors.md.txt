# Neighbors
```{eval-rst}
.. autofunction:: bullkpy.tl.neighbors

```

Compute a **k-nearest-neighbor (kNN) graph between samples** using a low-dimensional representation (typically PCA), closely mirroring Scanpy’s neighbors workflow but implemented in a **bulk-friendly, explicit** way.  

This function is the backbone for downstream steps such as **Leiden clustering**, **graph-based UMAP**, and **clustering quality scans**.

## Overview

neighbors builds two sample–sample graphs:  
- Distances graph.   
Sparse matrix of pairwise distances to the k nearest neighbors.

- Connectivities graph.  
A symmetrized, locally scaled Gaussian kernel graph, suitable for:  
	- Leiden / Louvain clustering
	- Graph-based UMAP
	- Graph-based metrics. 

Both graphs are stored in adata.obsp.  

## Requirements

Before calling neighbors, you must have:
- A low-dimensional representation in adata.obsm[use_rep]
	- Typically created by bk.tl.pca(adata)
- No missing samples (all rows in use_rep must be valid). 

If use_rep is missing, a KeyError is raised.

## Parameters

#### Core parameters

**n_neighbors**   
Number of nearest neighbors (k).  
Will be clipped to n_obs - 1 automatically.  
Default: 15.  

**use_rep**   
Key in adata.obsm containing the representation used for neighbor search.  
Default: "X_pca".  

**n_pcs**
Number of dimensions from use_rep to use.  
If None, uses all available dimensions.  
Default: 20.  

**metric**    
Distance metric used to find neighbors:  
- "euclidean" (default)
- "cosine" (implemented via L2 normalization + Euclidean distance). 

**key_added**  
Key under which parameters are stored in adata.uns.  
Default: "neighbors".  

## Method details

**Neighbor search**   
- Uses scipy.spatial.cKDTree for fast kNN queries
- Queries k + 1 neighbors and removes self-neighbors
- Distances are symmetrized using the minimum distance between pairs

**Local scaling → connectivities**   

Distances are converted to connectivities using a locally scaled Gaussian kernel:
- For each sample i, a local scale
(\sigma_i =) distance to its k-th nearest neighbor
- Connectivity between i and j:  
[  
w_{ij} = \exp\left(-\frac{d_{ij}^2}{2 \sigma_i \sigma_j}\right)   
]   

This produces a smooth, robust graph suitable for clustering and visualization.

## Stored results

After running neighbors, the following fields are populated:  

**In adata.obsp**   
- adata.obsp["distances"]
	- CSR sparse matrix
	- Shape: (n_obs, n_obs)
	- Symmetric nearest-neighbor distances
- adata.obsp["connectivities"]. 
	- CSR sparse matrix
	- Gaussian-kernel–weighted neighbor graph. 

**In adata.uns[key_added]**   
```python
adata.uns["neighbors"] = {
    "params": {
        "n_neighbors": 15,
        "n_pcs": 20,
        "use_rep": "X_pca",
        "metric": "euclidean",
    }
}
```

## Typical workflow

```python
# 1. Dimensionality reduction
bk.tl.pca(adata)

# 2. Build neighbors graph
bk.tl.neighbors(
    adata,
    n_neighbors=15,
    n_pcs=20,
    metric="euclidean",
)

# 3. Downstream analyses
bk.tl.cluster(adata, method="leiden")
bk.tl.umap_graph(adata)
```

## Notes and caveats
- This function operates on samples (obs), not genes.
- For cosine distance, vectors are L2-normalized first; results are equivalent to cosine similarity ranking.
- The connectivities graph is symmetric and weighted, unlike a raw directed kNN graph.
- Re-running neighbors will overwrite existing distances and connectivities.

## See also
	•	bk.tl.pca – compute PCA representation
	•	bk.tl.cluster – Leiden/Louvain clustering
	•	bk.tl.umap_graph – UMAP embedding from the neighbor graph
	•	bk.tl.leiden_resolution_scan – resolution benchmarking

