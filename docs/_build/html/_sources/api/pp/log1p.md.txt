# log1p

```{eval-rst}
.. autofunction:: bullkpy.pp.log1p

```


Log-transform expression values using `log(1 + x)`.

`log1p` applies a natural log–transformation to an expression matrix, which
reduces the influence of highly expressed genes and makes expression values
more suitable for visualization, clustering, and batch correction.   

This function is typically applied **after CPM normalization**.  

## What it does  
	•	Applies log(1 + x) to expression values
	•	Works with both dense and sparse matrices
	•	Writes results to a new layer by default
	•	Optionally overwrites adata.X

## Parameters   

**adata**  
AnnData object containing normalized expression data.

**layer**   
Input layer to transform.  
If None, uses adata.X.  
Default: "cpm".  

**out_layer**   
Name of the layer where log-transformed values are stored.  
Default: "log1p_cpm".  

**inplace_X**   
If True, also writes the transformed matrix to adata.X. 
Default: False. 

## Recommended workflow

A standard bulk RNA-seq preprocessing pipeline:

```python
bk.pp.set_raw_counts(adata)
bk.pp.normalize_cpm(adata)
bk.pp.log1p(adata)
```

After this, log-transformed values are available in:

```python
adata.layers["log1p_cpm"]
```
## Examples

#### Default usage (CPM → log1p)

```python
bk.pp.log1p(adata)
```

#### Transform a custom layer 

```python
bk.pp.log1p(
    adata,
    layer="cpm",
    out_layer="log_expr"
)
```

#### Overwrite adata.X

```python
bk.pp.log1p(
    adata,
    inplace_X=True
)
```

## Notes
	•	log1p is preferred over log(x) because it is well-defined for zero counts.
	•	This transformation assumes input data are non-negative.
	•	For batch correction and visualization (e.g. PCA, violin plots), log1p
expression is recommended.

## See also
	•	pp.normalize_cpm
	•	pp.set_raw_counts
	•	pp.highly_variable_genes
	•	pp.batch_correct_combat

