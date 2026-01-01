# Rank genes categorical

```{eval-rst}
.. autofunction:: bullkpy.tl.rank_genes_categorical

```

Rank genes associated with a categorical sample annotation.

This function performs a **bulk-friendly, Scanpy-like gene ranking** for categorical
variables stored in `adata.obs` (e.g. subtype, response group, mutation status).

It supports **two-group** and **multi-group** comparisons and reports
both statistical significance and effect size.

## What it does
Given a categorical variable in adata.obs, this function:  
1. Splits samples into:  
	•	a target group
	•	a reference group (either another category or “rest”). 

2.	Tests each gene for association with group membership. 

3.	Computes:  
	•	p-value
	•	FDR-corrected q-value
	•	effect size
	•	group and reference means
	•	log2 fold-change. 

4.	Returns a ranked DataFrame, optionally storing results in adata.uns

The behavior closely mirrors scanpy.tl.rank_genes_groups, but is designed for bulk RNA-seq.

## Supported comparisons

#### Two-group comparison

Used when:
	•	group is provided, or
	•	group is None and groupby has exactly 2 categories

Supported methods:

**Method - Test - Effect size**  
"mwu" (default) - Mann–Whitney U - Rank-biserial correlation  
"ttest" - Welch’s t-test - Cohen’s d (approx.)    
"kruskal" - Kruskal–Wallis - η² (rough)    
"anova" - One-way ANOVA - η² (rough)

#### Multi-group comparison

Used when groupby has >2 categories and a global test is requested:     

*Method - Test - Effect size*    
"kruskal" - Kruskal–Wallis - η² (rough)     
"anova" - One-way ANOVA - η² (rough)    

## Returned columns

The returned DataFrame contains:   

Column: Description    
- gene: Gene name.   
- pval: Raw p-value.   
- qval: Benjamini–Hochberg FDR.   
- effect_size: Method-dependent effect size.   
- mean_group: Mean expression in target group.   
- mean_ref: Mean expression in reference group.  
- log2FC: log2(mean_group / mean_ref).  

Results are sorted by qval, then pval.  

## Parameters

#### Group definition

**groupby**   
Column in adata.obs defining categories.    

**group**   
Target category.   
If None and groupby has exactly 2 categories, the first is used.  

**reference**   
Reference group:   
	•	"rest" (default): all other samples
	•	or a specific category name

####Expression source. 

**layer**    
Expression layer to use (e.g. "log1p_cpm").  
If None, uses adata.X.   

**genes**   
Optional list of genes to test.   
Default: all genes in adata.var_names.   

## Statistical method

**method**   
One of: "mwu", "ttest", "kruskal", "anova".

#### Storing results

**store_key**  
If provided, results are stored in:

```python
adata.uns["assoc"][store_key]
```

along with metadata describing the test.

## Examples

#### Binary comparison (default MWU)

```python
res = bk.tl.rank_genes_categorical(
    adata,
    groupby="Subtype",
    group="Basal",
    reference="rest",
    layer="log1p_cpm",
)
```

#### Two-group comparison with explicit reference

```python
res = bk.tl.rank_genes_categorical(
    adata,
    groupby="Response",
    group="Responder",
    reference="Non-responder",
    method="ttest",
)
```

#### Restrict to selected genes

```python
res = bk.tl.rank_genes_categorical(
    adata,
    groupby="RB1_mut",
    group="1",
    genes=["TP53", "CDKN2A", "E2F1"],
)
```

#### Store results in AnnData

```python
bk.tl.rank_genes_categorical(
    adata,
    groupby="Project_ID",
    group="ACC",
    store_key="assoc:Project_ID:ACC_vs_rest",
)
```

Later access:

```python
adata.uns["assoc"]["assoc:Project_ID:ACC_vs_rest"]["results"]
```

## Notes
	•	MWU is recommended for robust, non-parametric testing in bulk RNA-seq.
	•	Effect sizes are always reported and should be used alongside p/q-values.
	•	For very small groups (n < 2), the function raises an error.

## See also
	•	tl.rank_genes_groups (Scanpy, single-cell)
	•	pl.volcano
	•	pl.rankplot