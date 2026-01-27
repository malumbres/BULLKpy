# Partial correlations

```{eval-rst}
.. autofunction:: bullkpy.tl.partial_corr

```

Compute a **partial correlation** between two variables by first regressing out
one or more covariates taken from `adata.obs`.

This is useful when you want to quantify the association between `x` and `y`
*after removing the linear effects* of known confounders (e.g. batch, library size,
purity scores).

## What it does
1. Builds a design matrix from adata.obs[covariates] (with an intercept). 
2. Residualizes x and y with respect to these covariates. 
3. Computes a correlation between the residuals. 

Conceptually:

[
x^* = x - \hat{x}(\text{covariates}), \quad
y^* = y - \hat{y}(\text{covariates})
]

[
\text{partial corr}(x, y \mid \text{covariates}) = \text{corr}(x^*, y^*)
]

## Parameters

**adata**   
AnnData object containing covariates in .obs.

**x**   
1D numeric array (length = number of observations).

**y**   
1D numeric array (length = number of observations).

**covariates**    
List of obs column names to regress out.  
These should typically be numeric (categorical covariates should be encoded
beforehand).

**method**   
Correlation method:
	•	"pearson" (default)
	•	"spearman"

## Returns

A tuple:

```python
(r, pval, n)
```
- r — partial correlation coefficient
- pval — two-sided p-value
- n — number of observations used in the correlation

## Examples

1) Partial correlation between two genes controlling for library size

```python
x = adata.layers["log1p_cpm"][:, adata.var_names.get_loc("TP53")]
y = adata.layers["log1p_cpm"][:, adata.var_names.get_loc("MDM2")]

r, p, n = bk.tl.partial_corr(
    adata,
    x=x,
    y=y,
    covariates=["libsize"],
)
```

2) Gene–obs partial correlation adjusting for batch and purity

```python
x = bk.tl._get_gene_vector(adata, "PDCD1", layer="log1p_cpm")
y = adata.obs["signature_Tcell"].to_numpy()

r, p, n = bk.tl.partial_corr(
    adata,
    x=x,
    y=y,
    covariates=["Batch", "purity"],
    method="spearman",
)
```

## Notes
- Covariates are treated linearly (ordinary least squares residualization).
- This function does not handle missing values explicitly; NaNs in x, y,
or covariates should be removed or handled upstream.
- For large-scale analyses, prefer higher-level wrappers such as:
	•	top_gene_obs_correlations
	•	top_gene_gene_correlations
which can call partial correlation internally via batch/covariate logic.

## See also
	•	tl.top_gene_obs_correlations
	•	tl.top_obs_obs_correlations
	•	tl.gene_gene_correlations
	•	statsmodels.api.OLS (for explicit modeling)