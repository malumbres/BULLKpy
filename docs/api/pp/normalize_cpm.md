# Normalize cpm

```{eval-rst}
.. autofunction:: bullkpy.pp.normalize_cpm

```

Counts-per-million (CPM) normalization for bulk RNA-seq data.

`normalize_cpm` rescales gene expression values so that each sample has the same
total library size (`target_sum`, default 1e6). This is a standard normalization
step for bulk RNA-seq prior to log transformation, visualization, or batch
correction.

## What it does   
	•	Computes per-sample library sizes
	•	Scales expression values so that each sample sums to target_sum
	•	Writes normalized values to a new layer by default
	•	Optionally replaces adata.X with the normalized matrix
	•	Stores library sizes in adata.obs["libsize"]

## Parameters

**adata**  
AnnData object containing raw or unnormalized expression data.

**layer**   
Layer containing the input matrix.  
If None, uses adata.X.  
Default: "counts".  

**target_sum**   
Total counts per sample after normalization.  
Default: 1e6 (CPM).  

**out_layer**   
Name of the layer where normalized values will be stored.  
Default: "cpm".  

**inplace_X**   
If True, normalized values are also written to adata.X.  
Default: False.  

**eps**   
Small constant added to library sizes to avoid division by zero.  
Default: 1e-12.  

## Recommended workflow.  

CPM normalization should be applied after raw counts are stored and
before log transformation or batch correction:

```python
bk.pp.set_raw_counts(adata)
bk.pp.normalize_cpm(adata)
bk.pp.log1p(adata, layer="cpm", out_layer="log1p_cpm")
```

## Examples

#### Basic CPM normalization

```python
bk.pp.normalize_cpm(adata)
```

Normalized values are stored in:

```python
adata.layers["cpm"]
```
#### Use a custom input layer

```python
bk.pp.normalize_cpm(
    adata,
    layer="counts",
    out_layer="cpm"
)
```

#### Overwrite adata.X with CPM values
```python
bk.pp.normalize_cpm(
    adata,
    inplace_X=True
)
```

## Notes
	•	CPM normalization assumes that most genes are not differentially expressed
and is appropriate for exploratory analysis and visualization.
	•	For differential expression testing, raw counts (or proper model-based
normalization) should still be used.
	•	This function is bulk-RNA-seq oriented; it does not perform cell-level scaling.

## See also
	•	pp.set_raw_counts
	•	pp.log1p
	•	pp.batch_correct_combat
	•	pp.qc_metrics