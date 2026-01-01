# Gene-gene correlations

```{eval-rst}
.. autofunction:: bullkpy.tl.gene_gene_correlations

```

Correlate one gene against all other genes (or a specified subset).

This function computes correlations between a **single query gene** and many
other genes, returning the **strongest associations** ranked by correlation
strength. It is useful for:

- identifying co-expressed genes
- exploring pathway membership
- validating gene modules
- hypothesis-driven correlation analysis

## Purpose

gene_gene_correlations answers the question:  

*Which genes are most correlated with this gene of interest?*  

Compared to top_gene_gene_correlations, this function is:
	•	O(G) instead of O(G²)
	•	focused on a single anchor gene
	•	safe to run on genome-wide data

## Parameters

**adata**   
AnnData object containing expression data.

**gene**   
Query gene name. Must be present in adata.var_names.

**genes**   
Optional list of target genes to correlate against.  
If None, all genes except the query gene are used.   

**layer**
Expression layer to use (default: "log1p_cpm").

**method**    
Correlation method:
	•	"pearson"
	•	"spearman"

**top_n**   
Maximum number of correlated genes to return (default: 50).

**min_abs_r**     
Optional minimum absolute correlation threshold.

**use_abs**    
If True (default), ranking is based on |r|.    
If False, ranking uses signed correlation.    

**batch_key**   
Optional obs column specifying batch labels.

**batch_mode**   
Batch handling strategy:
	•	"none": ignore batch structure
	•	"within": compute correlations within batches
	•	"residual": regress out batch effects before correlation

**covariates**   
Optional numeric obs columns to regress out prior to correlation.

## What is computed

For each gene:
	•	correlation coefficient (r)
	•	p-value
	•	FDR-adjusted q-value
	•	number of samples used (n)

If batch handling is enabled, correlations are computed using
batch-aware residualization or within-batch aggregation.

## Returned value

A DataFrame sorted by correlation strength:

| column | description |
| ------------ | -------------------- |
| query_gene | Anchor gene |
| gene | Correlated gene |
| r | Correlation coefficient |
| pval | Raw p-value | 
| qval | FDR-adjusted p-value |
| n | Number of samples used |
| method | Correlation method |
| batch_key | Batch column used |
| batch_mode | Batch handling strategy |

## Examples

Genome-wide correlation for a single gene

```python
bk.tl.gene_gene_correlations(
    adata,
    gene="TP53",
)
```

Restrict to a pathway

```python
bk.tl.gene_gene_correlations(
    adata,
    gene="MYC",
    genes=hallmark_cell_cycle,
    top_n=20,
)
```

With batch correction

```python
bk.tl.gene_gene_correlations(
    adata,
    gene="CDKN1A",
    batch_key="Batch",
    batch_mode="residual",
)
```

Require strong correlations only

```python
bk.tl.gene_gene_correlations(
    adata,
    gene="BRCA1",
    min_abs_r=0.5,
)
```

## Interpretation notes
	•	High correlation suggests co-regulation, not causation
	•	Batch effects can strongly inflate correlations
	•	Always inspect results visually (scatter plots recommended)
	•	Use alongside:
	•	tl.top_gene_gene_correlations
	•	pl.scatter
	•	pl.heatmap

## See also
	•	tl.top_gene_gene_correlations
	•	tl.association
	•	pl.scatter
	•	pl.rankplot
