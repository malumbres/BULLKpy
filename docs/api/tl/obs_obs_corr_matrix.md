# Obs-Obs correlation matrix

```{eval-rst}
.. autofunction:: bullkpy.tl.obs_obs_corr_matrix

```

Compute a **correlation matrix between numeric sample-level variables** stored in `adata.obs`.

This function returns a **wide matrix** (rows = `focus`, columns = `against`) with correlation
coefficients, and can optionally also return the corresponding **p-value** and **q-value (FDR)**
matrices.

## What it does
- Selects numeric columns from adata.obs (optionally restricted by obs_keys)
- Computes correlations for all pairs focus × against
- Uses batch-aware correlation if batch_key / batch_mode / covariates are provided
- Optionally returns:
	•	raw p-values (return_p=True)
	•	BH-FDR q-values across all tested entries (return_q=True)

## Parameters

**adata**   
AnnData object with sample-level metadata in .obs.

**focus**   
One obs column name or a list of obs column names. These define the rows of the matrix.

**against**   
Optional list of obs column names defining the columns of the matrix.  
If None (default), uses all other numeric obs columns not in focus.  

**obs_keys**    
Optional list of obs column names used to define the numeric universe.  
If provided, only these columns are considered.  

**method**   
Correlation method:
	•	"pearson" (default)
	•	"spearman"

**drop_na**   
Whether to drop missing values before correlation (pairwise).  
(If your internal _corr_batch_aware always handles NaNs, this may have no effect.)

**batch_key, batch_mode, covariates**   
Batch-aware correlation settings:
	•	batch_key: obs column defining batch membership
	•	batch_mode: "none", "within", or "residual"
	•	covariates: numeric obs columns to regress out

**return_p**   
If True, also return a p-value matrix.

**return_q**  
If True, also return a q-value (BH-FDR) matrix computed across all tested pairs.

## Returns
	•	If return_p=False and return_q=False (default):  
R — DataFrame of correlation coefficients with shape (len(focus), len(against)).
	•	If return_p=True:  
(R, P) — correlation matrix and p-value matrix.
	•	If return_q=True:  
(R, Q) — correlation matrix and q-value matrix.
	•	If both are True:  
(R, P, Q)

## Examples

1) One focus variable vs all other numeric obs

```python
R = bk.tl.obs_obs_corr_matrix(
    adata,
    focus="pct_counts_mt",
)
R
```

2) Focus vs a restricted set of variables

```python
R, P, Q = bk.tl.obs_obs_corr_matrix(
    adata,
    focus=["libsize", "n_genes_detected"],
    against=["pct_counts_mt", "pct_counts_ribo", "purity"],
    method="spearman",
    return_p=True,
    return_q=True,
)
```

3) Batch-adjusted correlation matrix

```python
R = bk.tl.obs_obs_corr_matrix(
    adata,
    focus=["libsize", "pct_counts_mt"],
    against=["signature_IFNG", "signature_Tcell"],
    batch_key="Batch",
    batch_mode="residual",
)
```

## Interpretation notes
	•	Correlations can reflect technical confounding (library size, batch).
Use batch_mode="residual" or include relevant covariates when appropriate.
	•	If you request q-values, they are adjusted across all tested pairs in the matrix,
so the threshold is interpretable as a global FDR across the full grid.

## See also
	•	tl.top_obs_obs_correlations (returns a ranked long table instead of a matrix)
	•	tl.top_gene_obs_correlations
	•	pl.heatmap / your correlation heatmap plotting utilities