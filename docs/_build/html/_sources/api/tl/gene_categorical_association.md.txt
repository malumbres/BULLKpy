# Gene-Category association

```{eval-rst}
.. autofunction:: bullkpy.tl.gene_categorical_association

```

Association between gene expression and a categorical variable.

This function tests whether **gene expression differs across multiple
categories** of an observation (e.g. subtype, condition, project),
using **global statistical tests** rather than pairwise contrasts.

It is the natural **multi-group complement** to
`rank_genes_categorical`.

## What it does

For each gene, the function:

1. Splits samples by categories in adata.obs[groupby].  

2. Tests whether expression distributions differ across groups using:  
	•	Kruskal–Wallis (default, non-parametric)
	•	One-way ANOVA (parametric)

3. Computes an optional effect size.   

4. Applies multiple-testing correction.  

5. Returns a tidy results table (one row per gene).   

## When to use

Use gene_categorical_association when:  
	•	groupby has more than two categories
	•	You want a global test (not pairwise)
	•	You want to screen many genes at once
	•	Data may be non-Gaussian (bulk RNA-seq is often skewed)

Typical examples:   
	•	Tumor subtype association
	•	Batch or project effects
	•	Clinical categories (stage, grade, response)

## Parameters. 

**adata**   
AnnData object containing expression data and sample annotations

**groupby**   
Categorical column in adata.obs defining groups

**genes**    
Genes to test
	•	None (default): test all genes in adata.var_names

**layer**     
Expression layer to use (default: "log1p_cpm")

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

Returns a tidy DataFrame with one row per gene and columns:

| Column | Description |
| ---------- | -------------------- |
| groupby | Name of grouping variable |
|gene | Gene name |
|statistic | Test statistic (H or F) |
| pval | Raw p-value |
| qval | BH-adjusted p-value |
| effect | Effect size (ε² or η²) |
| n_groups | Number of groups tested |
| n | Total samples used |
| group_means | Mean expression per group (dict) |

Results are sorted by qval, then pval.

## Effect sizes

- ε² (epsilon-squared) – Kruskal–Wallis
	•	Range: 0–1
	•	Interpretable as proportion of variance explained
	•	Robust and recommended for bulk RNA-seq
- η² (eta-squared) – ANOVA
	•	Parametric analogue of ε²

## Examples

**Test gene–subtype association**   

```python
res = bk.tl.gene_categorical_association(
    adata,
    groupby="Subtype",
)
res.head()
```

**Restrict to selected genes**

```python
res = bk.tl.gene_categorical_association(
    adata,
    groupby="Project_ID",
    genes=["TP53", "RB1", "EGFR"],
)
```
   
**Use ANOVA instead of Kruskal–Wallis**

```python
res = bk.tl.gene_categorical_association(
    adata,
    groupby="Subtype",
    method="anova",
    effect_size="eta2",
)
```

## Notes
	•	Groups with fewer than min_group_size samples are ignored
	•	At least two valid groups are required per gene
	•	This is a global test:
	•	use posthoc_per_gene or pairwise_posthoc for pairwise differences
	•	Effect sizes are approximate but informative for ranking

## See also
	•	tl.rank_genes_categorical
	•	tl.posthoc_per_gene
	•	tl.pairwise_posthoc
	•	pl.violin
	•	pl.rankplot