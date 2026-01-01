# Top Gene-Obs correlations

```{eval-rst}
.. autofunction:: bullkpy.tl.top_gene_obs_correlations

```

Correlate one or more genes with **numeric** sample-level variables in `adata.obs`.

This function scans numeric columns in `adata.obs` (e.g. QC metrics, clinical variables,
scores, signatures) and returns the **strongest gene↔obs correlations**.

Typical uses:

- correlate genes with QC metrics (library size, % mito, etc.)
- associate genes with continuous phenotypes (age, tumor purity, signature score)
- quickly identify top gene–trait relationships


## Parameters

**adata**   
AnnData object with expression and sample metadata.

**gene**   
A gene name (string) or list of gene names. All must exist in adata.var_names.

**obs**   
Optional obs key (string) or list of obs keys to restrict the analysis to these columns.
Only numeric columns are allowed.

**obs_keys**   
Alternative explicit list of numeric obs columns to consider.   
If both obs and obs_keys are given, obs acts as a final filter.   

**layer**   
Expression layer to use (default: "log1p_cpm"). If None, uses adata.X.

**method**   
Correlation method:
	•	"pearson" (linear correlation)
	•	"spearman" (rank correlation)

**top_n**   
Maximum number of gene–obs pairs returned after ranking (default: 50).

**min_abs_r**    
Optional minimum absolute correlation threshold. Pairs with |r| < min_abs_r
are skipped before ranking.

**use_abs**   
If True (default), rank by |r| (strongest correlations regardless of sign).   
If False, rank by signed r.  

**batch_key, batch_mode, covariates**   
Optional batch-aware correlation controls (same behavior as other correlation utilities):
	•	batch_key: obs column defining batches
	•	batch_mode: "none", "within", or "residual"
	•	covariates: numeric obs columns to regress out before correlation

## Returned value

A DataFrame with one row per tested gene–obs pair, ranked by correlation strength.

Columns:

| column | description |
| --------------- | -------------------- |
| gene | Gene name |
| obs | Numeric obs column |
| r | Correlation coefficient |
| pval | Raw p-value |
| qval | FDR-adjusted p-value (BH) |
| n | Number of samples used |
| method | Correlation method |
| batch_key | Batch column used (or None) |
| batch_mode | Batch handling strategy |

## Examples

Correlate one gene vs all numeric obs

```python
df = bk.tl.top_gene_obs_correlations(
    adata,
    gene="TP53",
)
df.head()
```

Restrict to specific obs variables

```python
df = bk.tl.top_gene_obs_correlations(
    adata,
    gene="MKI67",
    obs=["purity", "signature_IFNG", "pct_counts_mt"],
    method="spearman",
)
```

Multiple genes

```python
df = bk.tl.top_gene_obs_correlations(
    adata,
    gene=["MYC", "CDKN1A", "MKI67"],
    top_n=30,
)
```

With batch adjustment

```python
df = bk.tl.top_gene_obs_correlations(
    adata,
    gene="EPCAM",
    batch_key="Batch",
    batch_mode="residual",
)
```

## Interpretation notes
	•	Correlation does not imply causation; treat results as exploratory.
	•	Batch effects and confounders can inflate correlations — use batch_mode="residual"
and/or covariates=[...] when appropriate.
	•	For visualization, follow up with a scatter plot of gene expression vs the selected obs.

## See also
	•	tl.gene_gene_correlations
	•	tl.top_gene_gene_correlations
	•	tl.association
	•	pl.scatter


