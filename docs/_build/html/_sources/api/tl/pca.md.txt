# Principal Component Analysis

```{eval-rst}
.. autofunction:: bullkpy.tl.pca

```

Principal Component Analysis (PCA) for **bulk-friendly** AnnData objects.

This function performs a classical PCA using SVD on an expression matrix
(`adata.X` or a selected layer), stores sample scores in `adata.obsm`,
gene loadings in `adata.varm`, and PCA metadata in `adata.uns`.

It is intentionally simple, transparent, and Scanpy-compatible in terms of
where results are stored, while remaining suitable for bulk RNA-seq data.

## What it does

1. Selects an expression matrix:  
- adata.layers[layer] if layer is not None
- otherwise adata.X.  

2. Optionally restricts genes to those marked as
adata.var["highly_variable"].  

3. Optionally centers and/or scales genes.  

4. Computes PCA using SVD:  
[   
X_c = U S V^\top  
]  

5. Stores:  
- Scores (samples × PCs) in adata.obsm["X_pca"]
- Loadings (genes × PCs) in adata.varm["PCs"]
- Explained variance and parameters in adata.uns[key_added]

## Parameters

**adata**   
AnnData object containing expression data.

**layer**   
Which matrix to use for PCA.
- Default: "log1p_cpm"
- If None, uses adata.X.   

**n_comps**  
Number of principal components to compute.  
The effective number is capped at:  
```python
min(n_comps, n_obs - 1, n_genes_used)
```
**center**   
If True (default), center each gene by subtracting its mean.

**scale**   
If True, scale each gene to unit variance after centering
(z-scoring per gene).

**use_highly_variable**  
If True, restrict PCA to genes where
adata.var["highly_variable"] == True.   

Notes:
- If the column does not exist, all genes are used with a warning.
- If it exists but has zero True values, an error is raised.

**key_added**   
Key under adata.uns where PCA metadata is stored.  
Default: "pca".  

## Output

**Sample scores**   

Stored in:
```python
adata.obsm["X_pca"]
```
Shape: (n_samples, n_comps)

These are the PCA coordinates of samples, suitable for:
- QC plots
- clustering
- downstream embeddings (UMAP, t-SNE, etc.)

**Gene loadings**   

Stored in:
```python
adata.varm["PCs"]
```
Shape: (n_genes, n_comps)

Notes:
- If use_highly_variable=True, genes not used in PCA receive NaN loadings.
- This keeps alignment with the full gene set, as expected by AnnData.

**PCA metadata**   

Stored in:
```python
adata.uns[key_added]
```

Contains:
- params
	- layer
	- n_comps
	- center
	- scale
	- use_highly_variable
	- n_vars_used.  
- variance: Eigenvalues (explained variance) for each PC.
- variance_ratio: Fraction of total variance explained by each PC.
- used_genes_mask: Boolean mask of genes actually used in PCA (length = n_vars).
- mean: Per-gene mean used for centering (only for used genes).
- std: Per-gene standard deviation used for scaling (or None if scale=False).

## Notes and best practices

**Recommended input**    
Use normalized, approximately Gaussian data such as:
- log1p_cpm
- log1p_tpm
- other log-transformed normalized expression

**Scaling**   
- scale=False (default) is often appropriate for bulk RNA-seq.
- scale=True may be useful when genes have very different variances.

**Highly variable genes** 
- Optional but recommended for large gene sets.
- Works well after running pp.highly_variable_genes.

**Sparse matrices**  
- Sparse inputs are densified internally.
- For very large datasets, ensure memory availability.

## Examples

Basic PCA on log1p CPM
```python
bk.tl.pca(adata)
```

PCA using only highly variable genes
```python
bk.pp.highly_variable_genes(adata)
bk.tl.pca(adata, use_highly_variable=True)
```

PCA with scaling
```python
bk.tl.pca(
    adata,
    layer="log1p_cpm",
    n_comps=30,
    center=True,
    scale=True,
)
```

Access results
```python
adata.obsm["X_pca"]              # sample scores
adata.varm["PCs"]                # gene loadings
adata.uns["pca"]["variance_ratio"]
```

## See also
	•	pp.highly_variable_genes – select informative genes
	•	tl.neighbors – build neighbor graph from PCA
	•	pl.corr_heatmap – sample correlation QC using PCA or expression
	•	scanpy.tl.pca – Scanpy’s reference implementation
