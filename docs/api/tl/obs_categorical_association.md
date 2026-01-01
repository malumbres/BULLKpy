# Obs-Category association

```{eval-rst}
.. autofunction:: bullkpy.tl.obs_categorical_association

```

Association between numeric sample-level variables and a categorical variable.

This function tests whether **numeric observation columns**
(e.g. QC metrics, clinical variables, scores)
**differ across categories** of another observation.

It is the **obs-level analogue** of
`gene_categorical_association`.

## What it does

For each numeric column in adata.obs, the function:     
1. Splits samples by categories in adata.obs[groupby]   
2. Tests whether values differ across groups using:   
	•	Kruskal–Wallis (default, non-parametric)  
	•	One-way ANOVA (parametric)   
3. Computes an optional effect size    
4. Applies multiple-testing correction    
5. Returns a tidy results table (one row per obs variable)   

## When to use

Use obs_categorical_association when you want to:
	•	Test QC metrics across conditions (e.g. library size vs batch)
	•	Assess clinical variables across subtypes
	•	Screen many numeric obs columns at once
	•	Perform a global multi-group test (not pairwise)

Typical examples:
	•	Are QC metrics different across batches?
	•	Do clinical scores vary across tumor subtypes?
	•	Is tumor purity associated with molecular class?

## Parameters

**adata**   
AnnData object containing observations in .obs

**groupby**  
Categorical column in adata.obs defining groups

**obs_keys**   
Numeric obs columns to test
	•	None (default): all numeric columns in adata.obs

**method**   
Global test:
	•	"kruskal" – Kruskal–Wallis H-test (default, non-parametric)
	•	"anova" – One-way ANOVA

**effect_size**   
Effect size to compute:
	•	"epsilon2" – Kruskal–Wallis effect size (recommended)
	•	"eta2" – ANOVA effect size
	•	None – skip effect size

**min_group_size**   
Minimum number of samples per group required to include that group

**adjust**   
Multiple-testing correction:
	•	"fdr_bh" – Benjamini–Hochberg (default)
	•	"none" – no correction

## Output

Returns a tidy DataFrame with one row per obs variable and columns:

| Column | Description |
| ---------- | -------------------- |
| groupby | Name of grouping variable |
| obs | Obs column name |
| statistic | Test statistic (H or F) |
| pval | Raw p-value |
| qval | BH-adjusted p-value |
| effect | Effect size (ε² or η²) |
| n_groups | Number of groups tested |
| n | Total samples used |
| group_means | Mean value per group (dict)| 

Results are sorted by qval, then pval.

## Effect sizes
- ε² (epsilon-squared) – Kruskal–Wallis
	•	Range: 0–1
	•	Proportion of variance explained by group membership
	•	Robust, non-parametric (recommended)
- η² (eta-squared) – ANOVA
	•	Parametric analogue of ε²

## Examples

**Test all numeric obs columns across groups**.   

```python
res = bk.tl.obs_categorical_association(
    adata,
    groupby="Batch",
)
res.head()
```

Test selected QC metrics only

```python
res = bk.tl.obs_categorical_association(
    adata,
    groupby="Subtype",
    obs_keys=["total_counts", "pct_counts_mt", "libsize"],
)
```

Use ANOVA instead of Kruskal–Wallis

```python
res = bk.tl.obs_categorical_association(
    adata,
    groupby="Project_ID",
    method="anova",
    effect_size="eta2",
)
```

## Notes
- Groups with fewer than min_group_size samples are ignored
- At least two valid groups are required per variable
- This is a global association test:
	•	Use posthoc_per_gene or pairwise_posthoc for pairwise comparisons
- Works naturally with QC metrics computed via pp.qc_metrics

## See also
	•	tl.gene_categorical_association
	•	tl.rank_genes_categorical
	•	tl.posthoc_per_gene
	•	tl.pairwise_posthoc
	•	pl.violin
	•	pl.rankplot
