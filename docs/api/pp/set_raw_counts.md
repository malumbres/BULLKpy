# Set raw counts

```{eval-rst}
.. autofunction:: bullkpy.pp.set_raw_counts

```

Store the current expression matrix as raw counts in an AnnData layer.

`set_raw_counts` is a convenience function to **preserve the original count matrix**
before any normalization or transformation steps are applied. This mirrors the
recommended Scanpy workflow of keeping raw data accessible for QC, filtering,
and reference.

## What it does
	•	Copies the current adata.X matrix into adata.layers[layer]
	•	Intended to be called immediately after reading count data
	•	Prevents accidental overwriting of existing raw layers by default


## Parameters  

**adata**  
AnnData object containing raw counts in adata.X.

**layer**   
Name of the layer where raw counts will be stored  
(default: "counts").  

**overwrite**   
If False (default), the function will not overwrite an existing layer
and will emit a warning instead.  
If True, any existing layer with the same name will be replaced.  

## Recommended usage

Call this function once right after loading the count matrix and before
any of the following steps:  
	•	normalization (CPM / TPM)
	•	log transformation
	•	batch correction
	•	filtering

This ensures that downstream steps can always reference unmodified counts.  

## Examples    

Basic usage after reading counts

```python
import bullkpy as bk

adata = bk.io.read_counts("counts.tsv")
bk.pp.set_raw_counts(adata)
```

Raw counts are now available as:

```python
adata.layers["counts"]
```
Overwrite an existing raw layer (not recommended unless intentional)

```python
bk.pp.set_raw_counts(adata, overwrite=True)
```
## Notes
	•	This function does not check whether the data in adata.X are truly raw.
It assumes the user calls it at the appropriate time.
	•	Downstream functions such as filter_samples and filter_genes can
reference the stored raw counts via layer="counts".

## See also
	•	io.read_counts
	•	pp.qc_metrics
	•	pp.filter_samples
	•	pp.filter_genes