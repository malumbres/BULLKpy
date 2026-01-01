# Filter genes

```{eval-rst}
.. autofunction:: bullkpy.pp.filter_genes

```

Filter genes by detection across samples.

`bk.pp.filter_genes` removes genes that are detected in too few samples.
A gene is considered *detected* in a sample if its expression value exceeds
a user-defined threshold (`min_expr`) in the selected expression layer.

This function is designed to work seamlessly with:
- raw count matrices, and
- normalized or log-transformed layers (e.g. `log1p_cpm`).

---

## What it does  

For each gene, the function counts in how many samples its expression is above
`min_expr`. Genes detected in fewer than `min_samples` samples are removed.   

Detection rule:  

A gene is detected in a sample if `expression > min_expr` in the chosen layer.

---

## When to use it  

Typical use cases include:  

- Removing genes expressed in only a handful of samples
- Reducing noise before PCA / clustering
- Speeding up downstream differential expression
- Making bulk RNA-seq analyses more robust

This is especially important for large cohorts (e.g. TCGA) where many genes
are barely expressed.  

---

## Basic usage

### Filter genes detected in at least 3 samples (counts)

```python
bk.pp.filter_genes(
    adata,
    min_samples=3,
    layer="counts",
)
```
This is the most common use for raw count data.  

### Filter genes using a log-transformed layer   

```python
bk.pp.filter_genes(
    adata,
    min_samples=10,
    min_expr=0.1,
    layer="log1p_cpm",
)
```

This is recommended when your data are already normalized or log-transformed.    


## Parameters    

**layer**   

Expression matrix to use for detection:  
	•	"counts" (default): raw counts
	•	any key in adata.layers
	•	None: uses adata.X

**min_expr**   

Minimum expression value required to count a gene as detected in a sample.  

Typical values:
	•	counts: min_expr = 0
	•	log-normalized data: min_expr ≈ 0.1


**min_samples**    

Minimum number of samples in which a gene must be detected to be kept.  

**In-place vs copy**   

```python
# By default, filtering is done in-place:
bk.pp.filter_genes(adata, inplace=True)

# To keep the original object unchanged:
adata_filt = bk.pp.filter_genes(adata, inplace=False)

```

## Stored metadata

The function records filtering details for reproducibility:  

```python
# In adata.var

adata.var["n_samples_detected"]
# Number of samples in which each gene was detected.

# In adata.uns
adata.uns["pp"]["filter_genes"]
```
Contains:  
	•	layer used  
	•	detection threshold  
	•	minimum samples  
	•	number of genes before and after filtering  


## Interaction with other preprocessing steps

A common preprocessing order is:  

```python
bk.pp.filter_genes(adata)
bk.pp.filter_samples(adata)
bk.pp.qc_metrics(adata)

```

Filtering genes early reduces noise and improves QC and dimensionality reduction.  
