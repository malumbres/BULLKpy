# Score genes from a dictionary of genes

```{eval-rst}
.. autofunction:: bullkpy.tl.score_genes_dict

```

Compute **multiple gene signature scores** in one call and store each score as a column
in `adata.obs`.

This is a convenience wrapper around [`score_genes`](score_genes.md), designed for
cases where you want to score **many pathways / signatures** in a consistent way
(e.g. Hallmark sets, cell-type programs, custom modules).


## What it does

For each entry in gene_sets:  
1. Calls score_genes(...).   
2. Computes a per-sample signature score:   
```python
mean(signature genes) − mean(matched control genes)
```
3. Stores the result in adata.obs under a standardized column name.

Each gene set is processed independently, but all share the same:
	•	expression layer
	•	control-gene pool
	•	binning strategy
	•	random seed
	•	scaling option

## Parameters

#### Required

**adata**   
AnnData object containing expression data.

**gene_sets**  
Mapping of {name → genes} defining multiple signatures.  
Example:
```python
{
    "T_cell": ["CD3D", "CD3E", "CD247"],
    "B_cell": ["MS4A1", "CD79A"],
    "IFN":    ["IFIT1", "ISG15", "MX1"],
}
```

#### Optional

**layer** (default: None)   
Expression layer to use.
	•	None → adata.X
	•	recommended: "log1p_cpm"

**prefix** (default: "")
Prefix added to each output column name.

**suffix**
suffix (default: "_score").  
Suffix added to each output column name.   

**gene_pool** (default: all genes)   
Pool of genes used for control-gene sampling.  
Signature genes are automatically excluded.  

**ctrl_size** (default: 50)   
Number of control genes sampled per signature gene.

**n_bins** (default: 25)
Number of expression bins for matching control genes.

**random_state** (default: 0)  
Random seed for reproducible control selection.

**scale** (default: False)    
If True, z-score genes before computing scores.

##  Returns
None.  

All results are written in-place to adata.obs.


## Output columns

For each gene set name name, the output column is:
```python
f"{prefix}{name}{suffix}"
```

Example:
```python
prefix="sig_"
suffix="_score"

→ sig_T_cell_score
→ sig_B_cell_score
→ sig_IFN_score
```

## Examples

1) Score multiple cell-type signatures
```python
gene_sets = {
    "Tcell": ["CD3D", "CD3E", "CD247"],
    "Bcell": ["MS4A1", "CD79A", "CD79B"],
    "Myeloid": ["LYZ", "CTSD", "LGALS3"],
}

bk.tl.score_genes_dict(
    adata,
    gene_sets,
    layer="log1p_cpm",
)  
```
Resulting columns in adata.obs:
```python
Tcell_score
Bcell_score
Myeloid_score
```

2) Use a prefix for tidy downstream plotting
```python
bk.tl.score_genes_dict(
    adata,
    gene_sets,
    prefix="sig_",
    layer="log1p_cpm",
)
```

3) Z-scored signatures with a custom gene pool
```python
bk.tl.score_genes_dict(
    adata,
    gene_sets,
    gene_pool=expressed_genes,
    scale=True,
    random_state=42,
)
```

## Notes & best practices

**Normalization first**.   
Use normalized, log-scale data (e.g. CPM → log1p).

**Consistent settings**.    
Because all signatures share the same parameters, this function is ideal
for comparative analyses across pathways or cell programs.

**Large gene sets**.   
For very large signatures, consider reducing ctrl_size or n_bins
to avoid overfitting control selection.

**Interpretation**.   
Scores are relative measures—use them for comparisons across samples or groups,
not as absolute expression levels.

## See also
	•	tl.score_genes – single-signature scoring
	•	pl.boxplot_with_stats – visualize scores by group
	•	pl.plot_corr_scatter – correlate scores with phenotypes
	•	pl.corr_heatmap – inspect relationships between multiple scores