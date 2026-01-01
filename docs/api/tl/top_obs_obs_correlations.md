# Top Obs-Obs correlations

```{eval-rst}
.. autofunction:: bullkpy.tl.top_obs_obs_correlations

```

Identify the **strongest correlations between numeric sample-level variables**
stored in `adata.obs`.

This function is useful to explore relationships between QC metrics, clinical
variables, signatures, and continuous phenotypes, either:

- one or more *focus* variables vs all other numeric obs, or  
- *focus* variables vs a restricted set of *against* variables.

## Parameters

**adata**   
AnnData object containing sample metadata in .obs.

**focus**   
One obs key (string) or a list of obs keys.  
These variables are correlated against other numeric obs columns.  

**against**   
Optional list of obs keys to correlate against.  
If None (default), all other numeric obs columns are used.  

**obs_keys**   
Optional list of obs columns used to define the numeric universe.  
If provided, only these columns are considered numeric candidates.  

**method**   
Correlation method:
	•	"pearson" (default)
	•	"spearman"

**top_n**   
Maximum number of obs–obs pairs returned after ranking (default: 100).

**min_abs_r**   
Optional minimum absolute correlation threshold.  
Pairs with |r| < min_abs_r are discarded before ranking.  

**use_abs**   
If True (default), ranking is based on |r|.    
If False, ranking is based on signed r.   

**batch_key, batch_mode, covariates**   
Optional batch-aware correlation settings:
	•	batch_key: obs column defining batches
	•	batch_mode: "none", "within", or "residual"
	•	covariates: numeric obs columns to regress out

## Returned value

A long-format DataFrame with one row per correlated obs–obs pair.

Columns:

| column | description |
| ---------- | -------------------- |
| focus | Focus obs variable |
| against | Other obs variable |
| r | Correlation coefficient |
| pval | Raw p-value |
| qval | FDR-adjusted p-value (BH) |
| n | Number of samples used |
| method | Correlation method |
| batch_key | Batch column used (or None) |
| batch_mode | Batch handling strategy |

## Examples

Correlate one QC metric vs all others

df = bk.tl.top_obs_obs_correlations(
    adata,
    focus="pct_counts_mt",
)
df.head()

Multiple focus variables

```python
df = bk.tl.top_obs_obs_correlations(
    adata,
    focus=["libsize", "n_genes_detected"],
    top_n=50,
)
```

Restrict comparisons to selected variables

```python
df = bk.tl.top_obs_obs_correlations(
    adata,
    focus="purity",
    against=["signature_IFNG", "signature_Tcell", "pct_counts_mt"],
    method="spearman",
)
```

With batch adjustment

```python
df = bk.tl.top_obs_obs_correlations(
    adata,
    focus="libsize",
    batch_key="Batch",
    batch_mode="residual",
)
```

## Interpretation notes
	•	Correlations may reflect confounding technical effects (batch, library size).
Use batch-aware options when appropriate.
	•	Strong correlations can guide feature selection or downstream modeling,
but should be validated visually (e.g. scatter plots).
	•	For categorical obs associations, use tl.categorical_association instead.

## See also
	•	tl.top_gene_obs_correlations
	•	tl.association
	•	pl.scatter
	•	pl.correlation_heatmap