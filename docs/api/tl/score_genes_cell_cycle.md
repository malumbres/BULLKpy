# Score cell cycle genes 

```{eval-rst}
.. autofunction:: bullkpy.tl.score_genes_cell_cycle

```

Compute **cell-cycle scores and phase assignment** (G1 / S / G2M) using a
Scanpy-compatible strategy, adapted to **bulk RNA-seq** (or pseudobulk) data.  

This function is a thin, opinionated wrapper around
[`score_genes`](score_genes.md) that reproduces the behavior of
`scanpy.tl.score_genes_cell_cycle`, while remaining robust for non-single-cell
use cases.

## What it does

1. Computes two gene signature scores per sample:  
	•	S phase score. 
	•	G2/M phase score. 
using matched control genes (expression-bin matched).

2. Assigns a discrete cell-cycle phase per sample using Scanpy’s rule:  

| Condition | Phase |
|------------------------ | ----------- |
| S_score > G2M_score and S_score > 0 | S |
| G2M_score > S_score and G2M_score > 0 | G2M |
|otherwise | G1 |

3. Stores results directly in adata.obs.  

## Output

The function writes three columns to adata.obs:

| Column | Type | Description |
|---------------- | ----------- | --------------------- |
| s_score | float | S-phase signature score |
| g2m_score | float | G2/M-phase signature score | 
| phase | categorical | "G1", "S", "G2M" (ordered) |

Default column names:
```python
S_score, G2M_score, phase
```
## Parameters

#### Required

**adata**    
AnnData object containing expression data.

**s_genes**   
Gene symbols defining the S-phase signature.

**g2m_genes**    
Gene symbols defining the G2/M-phase signature.

#### Optional (scoring behavior)

**layer** (default: None)   
Expression layer to use.
	•	None → adata.X.   
	•	recommended: "log1p_cpm".  

**gene_pool** (default: all genes)   
Background gene pool used to sample control genes.  
Signature genes are automatically excluded.  

**ctrl_size** (default: 50)    
Number of control genes sampled per signature gene.

**n_bins** (default: 25)   
Number of expression bins used to match control genes.

**random_state** (default: 0)   
Random seed for reproducible control selection.  
(Internally, G2M uses random_state + 1.). 

**scale** (default: False)    
If True, z-score genes before computing scores.

#### Optional (output names)

**s_score** (default: "S_score")   
Column name for S-phase score.

**g2m_score** (default: "G2M_score")    
Column name for G2/M score.  

**phase** (default: "phase")   
Column name for phase assignment.

## Returns
None.  

All results are written in-place to adata.obs.

## Example

1) Standard cell-cycle scoring
```python
bk.tl.score_genes_cell_cycle(
    adata,
    s_genes=S_GENES,
    g2m_genes=G2M_GENES,
    layer="log1p_cpm",
)
```

Creates:
```python
adata.obs["S_score"]
adata.obs["G2M_score"]
adata.obs["phase"]
```

2) Custom output column names
```python
bk.tl.score_genes_cell_cycle(
    adata,
    s_genes=S_GENES,
    g2m_genes=G2M_GENES,
    s_score="cc_S",
    g2m_score="cc_G2M",
    phase="cell_cycle",
)
```

3) Z-scored signatures with restricted gene pool
```python
bk.tl.score_genes_cell_cycle(
    adata,
    s_genes=S_GENES,
    g2m_genes=G2M_GENES,
    gene_pool=expressed_genes,
    scale=True,
)
```

## Notes & best practices

**Normalization matters**   
Use normalized, log-transformed data (e.g. log1p CPM or TPM).   
Raw counts are not recommended.   

**Bulk vs single-cell**   
While originally designed for single-cell data, this implementation
works well for bulk and pseudobulk, where phases should be interpreted
as dominant cell-cycle programs rather than discrete cell states.

**Mostly G1 warning**.  
If >95% of samples are labeled "G1", a warning is emitted.  
Common causes:  
	•	gene symbols not matching adata.var_names
	•	inappropriate layer (e.g. raw counts)
	•	very low proliferation signal

**Interpretation**.  
Cell-cycle phase here is a relative classification useful for QC,
covariate adjustment, or exploratory analysis—not a definitive cell state.

## See also
	•	tl.score_genes – single gene-set scoring
	•	tl.score_genes_dict – score multiple signatures
	•	pl.boxplot_with_stats – visualize S/G2M scores by group
	•	pl.corr_heatmap – inspect correlations between scores