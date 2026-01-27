# Pairwise posthoc

```{eval-rst}
.. autofunction:: bullkpy.tl.pairwise_posthoc

```

Pairwise post-hoc statistical tests between groups.

This function performs **all pairwise comparisons** between categories in a
grouping column, using either a **non-parametric** or **parametric** test.
It is typically used after a global association test to identify **which
groups differ from each other**.

## What it does
- Computes pairwise group comparisons for a numeric variable
- Supports:
	•	Mann–Whitney U test (default, non-parametric)
	•	Welch’s t-test (parametric)
- Reports:
	•	p-values
	•	BH-adjusted q-values
	•	effect sizes
	•	mean and median differences
- Returns a tidy DataFrame, easy to plot or export

## When to use

Typical use cases include:
- Post-hoc analysis after:
	•	rank_genes_categorical
	•	Kruskal–Wallis / ANOVA–like tests
- Pairwise comparison of:
	•	gene expression across multiple groups
	•	signature scores
	•	QC metrics

## Parameters

**df**   
Input DataFrame containing:
	•	one column with group labels
	•	one column with numeric values

**group_col**  
Column defining groups (default: "grp")

**value_col**   
Numeric column to test (default: "y")

**method**   
Statistical test:
	•	"mwu" – Mann–Whitney U test (two-sided, default)
	•	"ttest" – Welch two-sample t-test

**correction**   
Multiple-testing correction method:
	•	"bh" – Benjamini–Hochberg FDR (default)

**dropna**     
Whether to drop rows with missing group or value

## Output

Returns a DataFrame with one row per pairwise comparison and columns:

| Column | Description |
|---------- | -------------------- |
| group1, group2 | Compared groups |
| n1, n2 | Sample sizes |
| pval | Raw p-value |
| qval | BH-adjusted p-value |
| effect | Effect size (rank-biserial or Cohen’s d) |
| delta_mean | Mean(group1) − Mean(group2) |
| delta_median | Median(group1) − Median(group2) |

Results are sorted by increasing qval.

## Effect sizes

**MWU**  
Rank-biserial correlation
	•	Range: −1 to +1
	•	Sign indicates direction of shift

**t-test**    
Cohen’s d (approximate, pooled SD)

## Examples

**Pairwise tests for a single gene**    

```python
from bullkpy.tl import pairwise_posthoc

df = pd.DataFrame({
    "grp": adata.obs["Project_ID"],
    "y": adata[:, "TP53"].X.ravel(),
})

post = pairwise_posthoc(df, method="mwu")
post
```

**Post-hoc after categorical association**

```python
res = bk.tl.rank_genes_categorical(
    adata,
    groupby="Subtype",
)

posthoc = {}
for g in ["TP53", "RB1"]:
    posthoc[g] = bk.tl.posthoc_per_gene(
        adata,
        genes=[g],
        groupby="Subtype",
    )[g]
```

## Notes
	•	Requires ≥2 samples per group
	•	MWU is recommended for:
	•	small sample sizes
	•	non-Gaussian distributions
	•	t-test assumes approximate normality
	•	q-values are computed across all pairwise tests

## See also
	•	tl.rank_genes_categorical
	•	tl.posthoc_per_gene
	•	pl.violin
	•	pl.rankplot