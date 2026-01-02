# Score genes

```{eval-rst}
.. autofunction:: bullkpy.tl.score_genes

```

Compute a **signature score per sample** using a Scanpy-like approach, adapted for **bulk RNA-seq**.

This function assigns each sample a score defined as the **mean expression of a gene set** minus the **mean expression of matched control genes**, where controls are selected to match the expression distribution of the signature genes.  

The result is stored directly in `adata.obs`.

## What it does

For each sample:

[   
\text{score} = \mathrm{mean}(\text{expression of signature genes})    
- \mathrm{mean}(\text{expression of control genes})   
]    

Control genes are:  
	•	drawn from a gene pool (default: all genes),  
	•	excluding the signature genes,  
	•	matched by expression level using quantile bins (similar to scanpy.tl.score_genes).  

This makes the score robust to global expression biases and library-size effects.

## Parameters

#### Required

**adata** 
AnnData object containing expression data.

**genes**    
Sequence of gene names defining the signature to score.  
All genes must be present in adata.var_names.  

#### Optional

**score_name**   (default: "score")    
Name of the column written to adata.obs.

**layer** (default: None)   
Expression layer to use.
	•	None → use adata.X
	•	recommended: log-normalized layer (e.g. "log1p_cpm")

**gene_pool** (default: all genes)    
Pool of genes from which control genes are sampled.  
Signature genes are automatically removed from this pool.  

**ctrl_size** (default: 50)    
Number of control genes sampled per signature gene (before deduplication).

**n_bins** (default: 25)    
Number of expression bins used to match control genes to signature genes.  
Larger values give finer matching but may fail with small gene pools.  

**random_state** (default: 0)   
Random seed for reproducible control-gene sampling.  
Set to None for non-deterministic behavior.  

**scale**  (default: False)    
If True, gene expression is z-scored per gene before computing the score:  
[  
z = (x - \mu) / \sigma   
]    
This can be useful when genes differ strongly in variance.

## Returns

None.  

The result is written in-place to:
```python
adata.obs[score_name]
```

After running:
```python
adata.obs[score_name]
```

## Examples

1) Basic gene signature score
```python
bk.tl.score_genes(
    adata,
    genes=["CD3D", "CD3E", "CD247"],
    score_name="T_cell_score",
    layer="log1p_cpm",
)
```

2) Custom gene pool and reproducible sampling

```python
bk.tl.score_genes(
    adata,
    genes=my_signature,
    gene_pool=expressed_genes,
    ctrl_size=100,
    n_bins=20,
    random_state=42,
)
```

3) Z-scored signature (variance-normalized)
```python
bk.tl.score_genes(
    adata,
    genes=metabolic_genes,
    layer="log1p_cpm",
    scale=True,
    score_name="metabolism_score",
)
```

## Notes & recommendations
- Use normalized data.  
Raw counts are not recommended. Use log-scale normalized expression
(e.g. CPM → log1p).   
- Bulk-friendly by design.  
Unlike single-cell settings, bulk samples have fewer observations but more stable
gene means—this implementation is optimized for that regime.
- Control gene sampling. 
If you get errors about empty bins, try:
	•	lowering n_bins
	•	lowering ctrl_size
	•	providing a larger gene_pool
- Interpretation
	•	Scores are relative, not absolute.
	•	Best used for comparisons across samples or groups, not as standalone quantities.

## See also
	•	scanpy.tl.score_genes – original single-cell implementation
	•	tl.normalize_cpm, tl.log1p – recommended preprocessing
	•	pl.boxplot_with_stats – visualize scores by group
	•	pl.plot_corr_scatter – correlate scores with phenotypes
