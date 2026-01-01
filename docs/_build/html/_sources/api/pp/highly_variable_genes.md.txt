# Highly variable genes

```{eval-rst}
.. autofunction:: bullkpy.pp.highly_variable_genes



```

Identify highly variable genes (HVGs) in bulk RNA-seq data.


## Description

highly_variable_genes selects genes whose expression varies strongly across samples,
relative to genes of similar mean expression.  

This function follows the Scanpy / Seurat spirit, but is adapted for bulk RNA-seq:
variation is computed across samples, not cells, and works directly on normalized
expression layers.  

The method:
	1.	Computes per-gene mean and variance across samples
	2.	Calculates dispersion (variance / mean)
	3.	Log-transforms dispersion for numerical stability
	4.	Bins genes by mean expression
	5.	Computes a z-score of dispersion within each mean bin
	6.	Selects the top-ranked genes by normalized dispersion

Results are stored in adata.var.  

## Parameters

**adata**  
AnnData object with samples in `.obs` and genes in `.var`.

**layer**.   
  - `str | None`
  - Layer used to compute variability (default: `"log1p_cpm"`).  
    If `None`, uses `adata.X`.  

**n_top_genes**  
  - `int`
  - Number of highly variable genes to select (default: `2000`).

**n_bins**   
  - `int`
  - Number of mean-expression bins used for dispersion normalization.

*min_mean**   
  - `float`
  - Minimum mean expression for a gene to be considered.

**max_mean**   
  - `float`
  - Maximum mean expression for a gene to be considered.

**min_disp**   
  - `float`
  - Minimum dispersion threshold.

**key_added**   
  - `str`
  - Column name in `adata.var` where the HVG boolean flag is stored.


## Returns   

Results are stored in `adata.var`.

The following columns are added to adata.var:

Column --> Description.   
means --> Mean expression across samples.   
variances --> Expression variance across samples.   
dispersions --> Variance / mean.  
dispersions_norm --> Mean-binned z-score of log-dispersion.   
<key_added> --> Boolean flag indicating highly variable genes.   

## Example usage

```python
#Basic HVG selection
bk.pp.highly_variable_genes(adata)

# Using a different layer
bk.pp.highly_variable_genes(
    adata,
    layer="log2_tpm",
    n_top_genes=3000,
)

# Applying expression filters
bk.pp.highly_variable_genes(
    adata,
    min_mean=0.1,
    max_mean=5.0,
    min_disp=0.5,
)

# Subset AnnData to HVGs
adata = adata[:, adata.var["highly_variable"]].copy()
```

## Notes  
	•	Designed for bulk RNA-seq, not single-cell data
	•	Works with dense or sparse matrices
	•	Uses population variance (ddof=0)
	•	Mean binning is performed on log1p(mean) to improve stability
	•	Genes failing mean/dispersion thresholds are excluded from ranking

## Related functions   
	•	filter_genes – Remove lowly detected genes
	•	filter_samples – Sample-level QC filtering
	•	pca – Dimensionality reduction using HVGs

See also
	•	Scanpy: scanpy.pp.highly_variable_genes
	•	Seurat: FindVariableFeatures