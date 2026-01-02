# UMAP embedding

```{eval-rst}
.. autofunction:: bullkpy.tl.umap

```

Compute a **UMAP embedding** from an existing low-dimensional representation
(default: PCA), following Scanpy’s standard workflow.

This function is intentionally lightweight and bulk-friendly: it **does not**
recompute neighbors internally, but instead embeds samples directly from a
representation such as `X_pca`.

## When to use

Use umap when you want to:  
- visualize samples in 2D or 3D after PCA
- mirror Scanpy-style workflows (pca → umap)
- control UMAP parameters explicitly (neighbors, distance, metric, etc.). 

This function assumes you have already computed PCA (or another representation)
and stored it in adata.obsm.

## Parameters

#### Input / representation

**adata**   
AnnData object with samples in rows.

**use_rep**   
Key in adata.obsm containing the representation to embed.  
Default: "X_pca".  
If the key is missing, a KeyError is raised.  

**n_pcs**   
Number of dimensions from use_rep to use.  
If None, all dimensions are used.  
Default: 20.  

#### UMAP parameters

**n_neighbors**   
Size of the local neighborhood (controls local vs global structure).  
Typical values: 10–50.  
Default: 15.  

**min_dist**    
Minimum distance between embedded points.  
Smaller values → tighter clusters.  
Default: 0.5.  

**spread**   
Controls the overall scale of the embedding.  
Default: 1.0.  

**n_components**   
Output dimensionality of the embedding.  
- 2 → 2D UMAP
- 3 → 3D UMAP. 
Default: 2.  

**metric**   
Distance metric used in the high-dimensional space.  
Default: "euclidean".  

**random_state**   
Random seed for reproducibility.    
Default: 0.  

**init**   
Initialization method for UMAP.  
Common values: "spectral", "random".  
Default: "spectral".  

## Output

After running, the following fields are populated:  

**Embedding**.   
- adata.obsm["X_umap"].  
NumPy array of shape (n_obs, n_components) containing the UMAP coordinates.

**Metadata**  
- adata.uns["umap"].  
Dictionary storing UMAP parameters for provenance and reproducibility:

```python
{
  "params": {
    "use_rep": "X_pca",
    "n_neighbors": 15,
    "n_pcs": 20,
    "min_dist": 0.5,
    "spread": 1.0,
    "n_components": 2,
    "metric": "euclidean",
    "random_state": 0,
    "init": "spectral",
  }
}


## Dependencies
- Requires umap-learn:

```python
pip install umap-learn

If umap-learn is not installed, an informative ImportError is raised.

## Examples

Standard PCA → UMAP workflow
```python
bk.tl.pca(adata, layer="log1p_cpm", n_comps=30)
bk.tl.umap(adata)
```

Use more neighbors and tighter clusters
```python
bk.tl.umap(
    adata,
    n_neighbors=30,
    min_dist=0.1,
)
```

3D UMAP from the first 50 PCs
```python
bk.tl.umap(
    adata,
    n_components=3,
    n_pcs=50,
)
```

## Notes / caveats   
- This function embeds directly from the chosen representation; it does not
rebuild or reuse a neighbors graph.
- For consistency across analyses, you should usually use the same n_pcs
here as in PCA-based downstream steps (clustering, neighbors, etc.).
- Results are stochastic unless random_state is fixed.

