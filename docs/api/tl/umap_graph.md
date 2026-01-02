# UMAP graph

```{eval-rst}
.. autofunction:: bullkpy.tl.umap_graph

```

Compute a **UMAP embedding strictly from a precomputed neighbor graph** (e.g. Leiden/UMAP workflow where the graph is fixed and you want the embedding to reflect that exact graph).

Unlike `bk.tl.umap`, which embeds directly from a representation (e.g. PCA), **`umap_graph` uses `adata.obsp[graph_key]` as the source of neighborhood structure**. The chosen representation in `adata.obsm[use_rep]` is used **only for initialization**.


## When to use

Use umap_graph when you want:
- a UMAP embedding that matches a specific neighbors graph already computed (e.g. after tuning n_neighbors, metric, batch correction, etc.)
- reproducible embeddings from a stored graph (the structure is fixed in obsp)
- Scanpy-like consistency: neighbors → clustering → umap where UMAP is “based on the graph”

## Requirements

This function requires:
- adata.obsp[graph_key] exists.  
Typically created by bk.tl.neighbors(adata) as "connectivities".
- adata.obsm[use_rep] exists.  
Typically created by bk.tl.pca(adata) as "X_pca" (used only for init).
- umap-learn installed:
```python
pip install umap-learn
```

If any requirement is missing, a clear KeyError/ImportError is raised.

## Parameters

#### Graph and initialization

**graph_key**   
Key in adata.obsp containing the neighbor graph.  
Default: "connectivities".  

**use_rep**   
Key in adata.obsm used only to initialize the embedding.  
Default: "X_pca".  

**n_pcs**   
Number of dimensions from use_rep to use for initialization.  
Default: 20.  

#### UMAP layout parameters

**min_dist**   
Minimum spacing between points in the embedding.  
Smaller → tighter clusters.  
Default: 0.5.  

**spread**   
Overall scale of the embedding.  
Default: 1.0.  

**n_components**   
Output dimensionality (2D or 3D).  
Default: 2.   

**init**.   
Initialization method:    
- "spectral" (recommended, graph/structure-aware)
- "random". 
Default: "spectral".

**random_state**   
Random seed for reproducibility.  
Default: 0.  

#### Optimization controls

**negative_sample_rate**   
UMAP optimization parameter; larger can improve separation at cost of speed.  
Default: 5.  

**n_epochs**    
Number of training epochs.  
If None, umap-learn chooses a default; this implementation falls back to a safe default in some versions.  
Default: None.  

## Output

After running, the following are created:
- adata.obsm["X_umap_graph"]. 
Array of shape (n_obs, n_components) with graph-based UMAP coordinates.   
- adata.uns["umap_graph"]. 
Stores parameters for provenance, e.g.:  

```python
{
  "params": {
    "mode": "graph",
    "graph_key": "connectivities",
    "use_rep_init": "X_pca",
    "n_pcs_init": 20,
    "min_dist": 0.5,
    "spread": 1.0,
    "n_components": 2,
    "random_state": 0,
    "init": "spectral",
    "negative_sample_rate": 5,
    "n_epochs": None,
  }
}
```

## Notes / caveats
- The graph is treated as fixed: changing use_rep or n_pcs changes only the initialization, not which neighbors exist.
- This implementation uses some umap-learn internals (_fit_embed_data, graph_) for compatibility across versions. If you upgrade umap-learn and see errors, pinning a known working version may help.
- Ensure adata.obsp[graph_key] is a proper sparse connectivities matrix (CSR recommended). The function will coerce to CSR and remove explicit zeros.

## Examples

Standard neighbors → UMAP-from-graph
```python
bk.tl.pca(adata)
bk.tl.neighbors(adata, n_neighbors=15, n_pcs=20)
bk.tl.umap_graph(adata)
```

Tighter embedding, more optimization
```python
bk.tl.umap_graph(
    adata,
    min_dist=0.1,
    negative_sample_rate=10,
    n_epochs=800,
)
```

3D graph UMAP
```python
bk.tl.umap_graph(adata, n_components=3)
```