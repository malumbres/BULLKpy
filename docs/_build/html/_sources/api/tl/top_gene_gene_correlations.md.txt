# Top gene-gene correlations

```{eval-rst}
.. autofunction:: bullkpy.tl.top_gene_gene_correlations

```

Identify the strongest gene–gene correlations within a selected gene set.

This function computes pairwise correlations between genes and returns the
**top-ranked gene pairs** according to correlation strength.  
It is designed for **targeted panels, pathways, or signatures**, not genome-wide scans.

## Purpose

top_gene_gene_correlations helps answer questions such as:
	•	Which genes are most tightly co-regulated?
	•	Are genes within a pathway strongly correlated?
	•	Do correlations persist after accounting for batch effects?
	•	Which gene–gene relationships are strongest in my dataset?

Because the computation scales as O(G²), the function requires an explicit gene list
to avoid accidental genome-wide scans.

## Parameters

**adata**   
AnnData object containing expression data.

**genes**    
List of gene names to test.    
Required for safety (all-vs-all is intentionally disallowed).   

**layer**   
Expression layer to use (default: "log1p_cpm").

**method**  
Correlation method:
	•	"pearson"
	•	"spearman"

**top_n**    
Number of strongest gene pairs to return (default: 200).

**min_abs_r**   
Optional minimum absolute correlation threshold.    
Pairs below this value are discarded early.    

**use_abs**   
If True (default), ranking is based on |r|.   
If False, ranking uses signed correlation.   

**batch_key**   
Optional obs column specifying batch labels.

**batch_mode**  
How to handle batch effects:
	•	"none": ignore batches
	•	"within": compute correlations within batches
	•	"residual": regress out batch effects before correlation

**covariates**  
Optional obs columns to regress out before correlation
(e.g. library size, QC metrics).

## What is computed

For each gene pair (gene1, gene2):
	•	correlation coefficient r
	•	p-value
	•	Benjamini–Hochberg FDR (qval)
	•	number of samples used (n)

Batch-aware correlation is applied if batch_key is provided.

## Returned value

A tidy DataFrame with one row per gene pair:

| column | description |
| ---------- | -------------------- |
| gene1 | First gene |
| gene2 | Second gene |
| r | Correlation coefficient |
| pval | Raw p-value |
| qval | FDR-adjusted p-value |
| n | Number of samples used |
| method | Correlation method |
| batch_key | Batch column used (if any) |
| batch_mode | Batch handling strategy |

Rows are sorted by correlation strength (strongest first).

## Examples

**Basic usage**   

```python
bk.tl.top_gene_gene_correlations(
    adata,
    genes=["TP53", "MDM2", "CDKN1A", "BAX"],
)
```

**With batch correction**    

```python
bk.tl.top_gene_gene_correlations(
    adata,
    genes=hallmark_genes,
    batch_key="Batch",
    batch_mode="residual",
)
```

**Enforce minimum correlation strength**   

```python
bk.tl.top_gene_gene_correlations(
    adata,
    genes=genes_of_interest,
    min_abs_r=0.4,
    top_n=50,
)
```

## Interpretation tips
- High |r| suggests co-regulation or shared biology, not causality
- Strong correlations can arise from:
	•	shared pathway activity
	•	technical confounders
	•	cell-type composition
- Use batch_key and covariates to reduce confounding
- Consider visualizing top hits with:
	•	scatter plots
	•	correlation heatmaps
	•	network graphs

## Performance notes
	•	Runtime scales as O(G²)
→ keep genes lists reasonably small (tens to a few hundred).
	•	For genome-wide correlation analysis, use specialized methods instead.

## See also
	•	tl.association
	•	tl.gene_categorical_association
	•	pl.scatter
	•	pl.heatmap
	•	pl.network
