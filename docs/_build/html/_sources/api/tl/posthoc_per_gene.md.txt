# Posthoc per gene

```{eval-rst}
.. autofunction:: bullkpy.tl.posthoc_per_gene

```

Run pairwise posthoc tests for selected genes across all categories of a grouping variable.

This function performs **gene-wise pairwise comparisons** between all levels of a
categorical variable in `adata.obs`, returning a separate result table per gene.

It is typically used **after a global categorical association test**
(e.g. Kruskal–Wallis) to identify **which specific group pairs differ**.

## What it does

For each gene in genes, the function:  
	1.	Extracts expression values (from layer or adata.X)
	2.	Groups samples by groupby
	3.	Performs all pairwise group comparisons
	4.	Returns results as a dictionary:

```python
{
    "GENE1": DataFrame,
    "GENE2": DataFrame,
    ...
}
```
Each DataFrame contains pairwise statistics between categories.

## Statistical methods

Supported pairwise tests:  

|Method | Test | Notes |
|---------- | ---------- | ----------|
|"mwu" (default) | Mann–Whitney U | Non-parametric, robust|    
|"ttest" | Welch’s t-test | Assumes approximate normality|   

The actual pairwise testing is delegated to:

```python
bk.tl.pairwise_posthoc(df, method=...)
```
## Returned format

Each value in the returned dictionary is a DataFrame with pairwise results.

Typical columns include (depending on pairwise_posthoc implementation):

| Column | Description |
|-------- | ---------- |
| group1 | First group |
| group2 | Second group |
| pval | Raw p-value |
| qval | FDR-corrected p-value |
| effect_size | Pairwise effect size |
| mean_1 | Mean expression in group1 |
| mean_2 | Mean expression in group2 |

## Parameters

#### Gene selection

**genes**
List of gene names to test.  
All genes must be present in adata.var_names.  

### Group definition

**groupby**  
Categorical column in adata.obs defining groups.

#### Expression source

**layer**
Expression layer to use (e.g. "log1p_cpm").  
If None, uses adata.X.  

#### Statistical test

**method**
Pairwise test to apply:
	•	"mwu" (default)
	•	"ttest"

## Examples

#### Pairwise testing after global association

```python
# Global test
res = bk.tl.rank_genes_categorical(
    adata,
    groupby="Subtype",
    group="Basal",
)

# Posthoc comparisons for selected genes
posthoc = bk.tl.posthoc_per_gene(
    adata,
    genes=["TP53", "E2F1", "MYC"],
    groupby="Subtype",
)
```

Inspect pairwise results for one gene

```python
posthoc["TP53"]
```

Use t-test instead of MWU

```python
posthoc = bk.tl.posthoc_per_gene(
    adata,
    genes=["CDKN2A"],
    groupby="Project_ID",
    method="ttest",
)
```

## Typical workflow

**1. Global test**.   
rank_genes_categorical or cat_cat_association

**2. Select genes of interest**.  
Based on q-value and effect size

**3. Posthoc testing**.  
posthoc_per_gene

**4. Visualization**
	•	boxplots / violins
	•	effect size heatmaps
	•	pairwise significance tables

## Notes
	•	This function does not apply multiple-testing correction across genes
(only within each gene’s pairwise comparisons).
	•	For many genes × many groups, computation can be expensive.
	•	Non-parametric tests are recommended for heterogeneous bulk cohorts.

## See also
	•	tl.rank_genes_categorical
	•	tl.pairwise_posthoc
	•	pl.violin
	•	pl.rankplot