# Negative binomial GLM differential expression

```{eval-rst}
.. autofunction:: bullkpy.tl.de_glm

```

Run **negative binomial GLM differential expression** (bulk-friendly) using a **formula-defined design matrix** and a **two-level contrast**. This is intended for raw (integer-like) counts and supports adding covariates (e.g., batch) in the model.

The implementation is “DESeq2-like” in spirit:  
- it computes **size factors** and uses a **log(size factor)** offset
- estimates **gene-wise dispersion** from normalized counts
- optionally **shrinks dispersion** toward a trend
- fits a **gene-wise NB-GLM** and performs a **Wald test** for a specified contrast
- reports log2FC, Wald z, p-values, and BH FDR q-values

## When to use

Use de_glm when you want:   
- two-group DE while controlling for additional covariates (e.g., batch, sex, purity)
- a count model (NB) instead of a t-test on log-expression
- a workflow similar to bulk RNA-seq DE tools (DESeq2/edgeR style). 

You should provide raw counts in adata.layers["counts"] (or another layer via layer_counts).

## Parameters

#### Core inputs

**adata**  
AnnData with samples in rows (n_obs) and genes in columns (n_vars).

**formula**   
A Patsy-like formula string defining the design matrix using adata.obs columns.
Examples:  
- "~ Subtype"
- "~ Subtype + Batch"
- "~ Subtype + Batch + Sex". 
The design matrix is constructed via an internal helper (_design_matrix_from_formula).

**contrast**   
Required in the current implementation.  
A tuple (var, level_a, level_b) meaning level_a − level_b for a categorical variable.  
Example:
- ("Subtype", "Basal", "Luminal") → Basal vs Luminal (Basal − Luminal)
Internally a contrast vector c is built from the coefficient names (coef_names) using _contrast_vector.

**layer_counts**   
Layer name containing the raw counts matrix (samples × genes).    
Default: "counts".    

#### Dispersion options

**shrink_dispersion**  
If True, shrink gene-wise dispersion estimates toward a trend (stabilizes estimates when n is small).  
Default: True.

**prior_df**  
Prior degrees of freedom controlling the strength of dispersion shrinkage (higher → stronger shrinkage).  
Default: 10.0.  

#### Output storage

**key_added**   
Where results are stored in adata.uns.  
Default: "de_glm".  

**add_intercept**   
Whether to include an intercept term in the design matrix.  
Default: True.  

## Output

Results are stored as:
```python
adata.uns[key_added][contrast_name] = {
  "method": ...,
  "formula": ...,
  "contrast": ...,
  "coef_names": ...,
  "layer_counts": ...,
  "shrink_dispersion": ...,
  "prior_df": ...,
  "results": <pd.DataFrame>
}
```

Where:
- contrast_name is constructed as:
	- f"{var}:{level_a}_vs_{level_b}"
	- Example: "Subtype:Basal_vs_Luminal"

**Result table columns**.   

The stored DataFrame (results) contains:
- gene : gene name (from adata.var_names)
- log2FC : estimated log2 fold-change for the contrast (level_a − level_b)
- wald_z : Wald z-statistic for the contrast
- pval : two-sided p-value from Wald z
- qval : BH-adjusted p-value (FDR)
- dispersion : final dispersion used per gene (after optional shrink)
- mean_norm : mean of normalized counts per gene (used in dispersion estimation). 

The table is sorted by qval ascending.

## Statistical details (high level)

1. Design matrix. 
Built from adata.obs using the provided formula (plus optional intercept).  
	
2. Counts and offsets. 
- counts: Y = adata.layers[layer_counts]
- size factors: sf = deseq2_size_factors(Y)
- offset: log(sf) included in the GLM. 

3. Dispersion. 
- estimates dispersion from normalized counts (Y / sf)
- optional shrinkage to a dispersion trend. 

4. Model fitting (per gene). 
Fits:  
- NB(mu, alpha) with log link
- log(mu) = Xβ + offset.  

5. Wald test for contrast
Tests cᵀβ = 0 using:  
- z = (cᵀβ) / sqrt(cᵀCov(β)c)
- p-value from standard normal approximation

## Requirements / dependencies
- statsmodels is required:
	- installed via pip install statsmodels
- scipy is used for the normal CDF in Wald p-values.  

If statsmodels is missing, the function raises an ImportError with install instructions.

## Examples

Basic DE with batch covariate
```python
bk.tl.de_glm(
    adata,
    formula="~ Subtype + Batch",
    contrast=("Subtype", "Basal", "Luminal"),
    layer_counts="counts",
)
```

Retrieve results
```python
res = adata.uns["de_glm"]["Subtype:Basal_vs_Luminal"]["results"]
res.head()
```

Turn off dispersion shrinkage
```python
bk.tl.de_glm(
    adata,
    formula="~ Subtype + Batch",
    contrast=("Subtype", "Basal", "Luminal"),
    shrink_dispersion=False,
)
```

## Notes / caveats
- The current implementation requires contrast (no “test all coefficients” mode yet).
- This is a gene-wise GLM loop; it can be slower than vectorized approaches for very large gene counts.
- Counts with total sum 0 for a gene are skipped (remain default p=1 / NaNs for some stats).
- The contrast mechanism assumes a fairly standard dummy encoding for categorical variables (handled by the internal design/contrast helpers). If coefficient naming differs from expectations, _contrast_vector may fail or target the wrong coefficient—double-check coef_names stored in adata.uns.