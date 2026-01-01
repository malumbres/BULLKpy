# Differential expression

```{eval-rst}
.. autofunction:: bullkpy.tl.de

```

Two-group differential expression analysis for bulk RNA-seq data.

This function performs **gene-wise differential expression** between two groups
defined by a categorical column in `adata.obs`. It is designed for **bulk RNA-seq**
and integrates tightly with the BULLKpy workflow and `AnnData`.

Results are stored in `adata.uns` in a Scanpy-like structure and can be visualized
with downstream plotting functions such as `pl.volcano`, `pl.rankplot`, or used
for GSEA.

## What it does

	•	Compares two groups of samples (group vs reference)
	•	Supports fast exploratory and proper statistical models
	•	Computes per-gene statistics (effect size, p-value, q-value, log2FC)
	•	Stores results in a contrast-aware location in adata.uns

## Methods

**welch_ttest** (default)
	•	Welch two-sample t-test
	•	Operates on normalized expression (default: log1p_cpm)
	•	Fast and suitable for:
	•	QC
	•	exploratory analysis
	•	visualization
	•	Assumes approximately Gaussian data

**Effect size**
	•	Mean difference
	•	log2 fold change

**nb_glm**
	•	Negative binomial GLM (bulk RNA-seq–appropriate)
	•	Operates on raw counts
	•	Includes library-size offset
	•	Requires statsmodels
	•	Slower but statistically rigorous

Recommended for:  
	•	publication-quality DE
	•	datasets with strong mean–variance relationships

## Parameters

**groupby**   
Column in adata.obs defining sample groups

**group / reference**   
Two categories to compare

**method**  
"welch_ttest" or "nb_glm"

**layer_expr**   
Expression layer for Welch test (default: "log1p_cpm")

**layer_counts**   
Raw count layer for NB GLM (default: "counts")

**key_added**   
Key under adata.uns where results are stored

## Output structure

Results are stored in:

```python
adata.uns[key_added][f"{groupby}_{group}_vs_{reference}"]
```

Each entry contains:

```python
{
  "method": "...",
  "groupby": "...",
  "group": "...",
  "reference": "...",
  "layer": "...",           # or layer_counts
  "results": DataFrame
}
```

## Results table columns

Typical columns include:
	•	gene
	•	log2FC
	•	pval
	•	qval
	•	mean_group
	•	mean_ref
	•	t (Welch only)
	•	stat (NB GLM only)

(Exact columns depend on the method)

## Examples

**Welch t-test (exploratory)**    

```python
bk.tl.de(
    adata,
    groupby="Project_ID",
    group="LUAD",
    reference="LUSC",
)
```

**NB GLM (rigorous bulk DE)**    

```python
bk.tl.de(
    adata,
    groupby="RB1_mut",
    group="1.0",
    reference="0.0",
    method="nb_glm",
)
```

## Access results

```python
res = adata.uns["de"]["RB1_mut_1.0_vs_0.0"]["results"]
res.head()
```

## Plot results

```python
bk.pl.volcano(
    res,
    title="RB1 mutant vs WT",
)

bk.pl.rankplot(
    res=res,
    direction="both",
)
```

## Notes
	•	Requires ≥2 samples per group
	•	welch_ttest assumes approximately normal expression
	•	nb_glm is preferred when raw counts and library sizes are available
	•	Multiple contrasts can be stored under the same key_added

## See also
	•	pl.volcano
	•	pl.rankplot
	•	tl.rank_genes_categorical
	•	tl.gsea_preranked