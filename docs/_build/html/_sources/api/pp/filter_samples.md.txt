# Filter samples

```{eval-rst}
.. autofunction:: bullkpy.pp.filter_samples

```

Filter samples based on QC thresholds (with graceful fallback computation).

`bk.pp.filter_samples` removes low-quality samples using QC metrics stored in
`adata.obs` **if available**, and **computes missing metrics from the expression matrix**
when needed (Scanpy-like behavior, but adapted for bulk RNA-seq).

A key design choice: **there are no hard requirements** that any QC columns exist.
If a threshold is requested but the corresponding metric cannot be computed, the
function **skips that threshold with a warning** instead of failing.

---

## What it does  

Given one or more thresholds, the function builds a boolean mask over samples and keeps
only samples that pass all requested filters.  

Supported thresholds:  

- `min_counts` → keeps samples with `total_counts >= min_counts`  
- `min_genes` → keeps samples with `n_genes_detected >= min_genes`  
- `max_pct_mt` → keeps samples with `pct_counts_mt <= max_pct_mt`   
- `max_pct_ribo` → keeps samples with `pct_counts_ribo <= max_pct_ribo`    

If these columns are missing in `adata.obs`, the function attempts to compute them from
`adata.layers[layer]` (default `"counts"`) or `adata.X` if `layer=None`.  

Computed metrics are **stored back into `adata.obs`** so they can be reused by other functions.

---

## Inputs and requirements

### Expression matrix source

By default, metrics are computed from:

- `adata.layers["counts"]` (recommended for raw counts), otherwise
- `adata.X` if `layer=None`

If your data are already normalized or log-transformed, the filters still run,
but interpretation changes (see below).

### Gene masks for MT and ribosomal genes. 

To compute mitochondrial/ribosomal percentages, the function needs to identify
mitochondrial and ribosomal genes. It does this using:  

- provided masks in `adata.var` (recommended), e.g. `mt_var_key="mt"`  
- otherwise, heuristics based on gene names (implementation-dependent)   

If gene masks cannot be inferred, `%mt` / `%ribo` filters are skipped with a warning.

---

## Basic usage

### Filter only by detected genes

```python
bk.pp.filter_samples(
    adata,
    min_genes=5000,
    layer="counts",
)
```

### Filter by multiple QC criteria
```python
bk.pp.filter_samples(
    adata,
    min_counts=1e6,
    min_genes=8000,
    max_pct_mt=20,
    max_pct_ribo=50,
    layer="counts",
)
```

### Using precomputed QC metrics

If you already ran:
```python
bk.pp.qc_metrics(adata)
```

then filter_samples will reuse the columns in adata.obs instead of recomputing.  

This is faster and guarantees consistency with your QC plots.  


### Working with non-count data

If your matrix is already log-transformed (e.g. "log1p_cpm"), you can still filter,
but some thresholds become less meaningful:
	•	min_genes remains meaningful (counts non-zero genes above threshold)
	•	min_counts depends on the scale of the data (sum of transformed values)
	•	%mt and %ribo become relative fractions of transformed values

For normalized/log data, filtering is best done using relative thresholds
or distribution inspection (e.g. percentile cutoffs).

### Custom definition of “detected genes”

The number of detected genes is computed as:

number of genes with expression > expr_threshold_for_genes
```python
bk.pp.filter_samples(
    adata,
    min_genes=6000,
    expr_threshold_for_genes=0.1,
    layer="log1p_cpm",
)
```

### In-place vs copy


```python
# By default, filtering is done in-place (faster):
bk.pp.filter_samples(adata, min_genes=5000, inplace=True)

# To return a new AnnData without modifying the original:
adata_filt = bk.pp.filter_samples(adata, min_genes=5000, inplace=False)

```

### Stored metadata

This function records filter parameters and before/after counts in:
```python
adata.uns["pp"]["filter_samples"]
```
This is useful for reproducibility and reporting.

