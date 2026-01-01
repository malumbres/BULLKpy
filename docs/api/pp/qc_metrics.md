# QC metrics

```{eval-rst}
.. autofunction:: bullkpy.pp.qc_metrics

```

Compute basic quality-control (QC) metrics for bulk RNA-seq samples.

`bk.pp.qc_metrics` calculates commonly used per-sample QC statistics
and stores them in `adata.obs`, enabling downstream filtering and QC plots.

## Purpose

This function computes **sample-level quality metrics** such as:

- total library size
- number of detected genes
- fraction of mitochondrial reads
- fraction of ribosomal reads

These metrics are essential for:
- identifying low-quality samples
- detecting outliers
- visual QC inspection before normalization and PCA

## Metrics computed

Depending on gene annotations available in `adata.var`, the following
columns may be added to `adata.obs`:

| Column | Description |
|------|------------|
| `total_counts` | Sum of counts across all genes per sample |
| `n_genes_detected` | Number of genes with non-zero counts per sample |
| `pct_counts_mt` | Percentage of counts from mitochondrial genes |
| `pct_counts_ribo` | Percentage of counts from ribosomal genes |

## Basic usage

```python
import bullkpy as bk

bk.pp.qc_metrics(adata)
```

After running this function, QC metrics are available in adata.obs:

```python
adata.obs[["total_counts", "n_genes_detected"]].head()
```

## Basic usage

### Gene annotation requirements

#### Mitochondrial genes

To compute mitochondrial fractions, gene names must contain a
mitochondrial prefix (by default "MT-"):

```python
adata.var_names[:5]
# ['MT-ND1', 'MT-CO1', ...]
```

#### Ribosomal genes

Ribosomal fractions are computed from gene names starting with
"RPS" or "RPL".

If no mitochondrial or ribosomal genes are detected, the corresponding
metrics are silently skipped.

## Handling transformed data

If your expression matrix is already log-transformed or normalized,
QC metrics can still be computed, but interpretation changes:
	•	total_counts reflects summed expression values (not raw counts)
	•	n_genes_detected remains meaningful
	•	percentages reflect relative expression contributions

For TCGA-style data that is already normalized, QC metrics are still
useful for relative sample comparison, not absolute thresholds.

## Example

```python
bk.pp.qc_metrics(adata)

bk.pl.qc_metrics(
    adata,
    color="n_genes_detected",
)
```

This produces:
	•	library size vs detected genes scatter
	•	QC histograms
	•	optional color-coded overlays


## Downstream usage

QC metrics are commonly used for:

### Sample filtering

```python
adata = adata[
    (adata.obs["n_genes_detected"] > 5000) &
    (adata.obs["pct_counts_mt"] < 20)
]
```

### PCA diagnostics
```python
bk.pl.pca_scatter(
    adata,
    color="pct_counts_mt",
)
```

## Output
	•	QC metrics are added in-place to adata.obs
	•	The function returns the modified AnnData object
	•	No changes are made to .X, .layers, or .var


## Notes and best practices
	•	Run qc_metrics before normalization and PCA
	•	For large cohorts (e.g. TCGA), avoid hard thresholds; inspect distributions
	•	QC metrics are compatible with log-transformed data but should be interpreted carefully
