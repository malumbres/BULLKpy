# PCA loadings heatmap

```{eval-rst}
.. autofunction:: bullkpy.pl.pca_loadings_heatmap

```

Heatmap of PCA loadings for the **union of top-loading genes** across one or more principal components.

This plot helps interpret multiple PCs at once by showing which genes load strongly on each selected
PC, optionally separating positive/negative contributors and clustering genes (and/or PCs).


## What it does. 

1. Reads PCA loadings from adata.varm[loadings_key] (shape: n_genes × n_pcs).    
2. For each PC in pcs, selects top genes:    
- If use_abs=True: top n_top genes by |loading|.  
- Else: top n_top positive genes, plus (optionally) top n_top negative genes.   
3. Forms the union of selected genes across PCs (preserving first-seen order).  
4. Builds a matrix: genes × PCs of loadings for those genes.  
5. Optionally z-scores each gene across PCs.   
6. Plots as:  
- seaborn.clustermap if clustering is enabled (Scanpy-like dendrograms), else
- seaborn.heatmap.   

## Parameters

#### Core inputs

**adata**: AnnData.   
Must contain PCA loadings at adata.varm[loadings_key] (typically from bk.tl.pca()).   

**pcs**: Sequence[int], default (1, 2, 3).  
PCs to include (1-based indexing, e.g. (1,2,3) means PC1–PC3).    

**n_top**: int, default 15.   
Number of genes selected per PC (per sign if signed mode).    

**loadings_key**: str, default "PCs".  
Key in adata.varm where loadings are stored.   

#### Gene selection mode. 

**use_abs**: bool, default False.  
- False: signed selection (positive + optional negative).
- True: select genes by absolute loading magnitude only.   

**show_negative**: bool, default True. 
Only used when use_abs=False. If True, includes top negative loadings per PC.  

####Gene labels 

**gene_symbol_key**: str | None, default None.  
If provided and present in adata.var, uses this column for gene labels instead of
adata.var_names.   
Note: If symbols are duplicated, the function uses the first occurrence of each symbol
when mapping back to indices (a pragmatic choice for plotting).   

#### Transformations

**z_score**: bool, default False. 
If True, z-scores each gene across PCs:   
[    
z_{gene,pc} = \frac{loading_{gene,pc} - \mu_{gene}}{\sigma_{gene}}   
]   
Useful to emphasize relative PC preference per gene rather than absolute magnitude.     

####Clustering / plotting.  

**cluster_genes**: bool, default True.  
Cluster rows (genes) with hierarchical clustering (uses seaborn.clustermap).   

**cluster_pcs**: bool, default False.   
Cluster columns (PCs). Often off by default because PC order is meaningful.    

**cmap**: str, default "vlag".  
Diverging colormap suited for signed loadings.    

**center**: float, default 0.0.   
Center value for diverging colormap normalization (0 is typical for loadings).    

**figsize**: tuple[float, float] | None. 
Auto-sized if None based on number of genes and PCs.    

**title**: str | None.  
Default: "PCA loadings heatmap".  

####Output controls

**save**: str | Path | None  
If provided, saves the figure via the project’s _savefig() utility.   

**show**: bool, default True.  
Whether to display the plot with plt.show().    

## Returns

**fig**: matplotlib.figure.Figure    
If clustering is enabled, this is the clustermap figure (cg.fig); otherwise the standard
heatmap figure.   

## Raises
- ImportError:  If seaborn is not installed.
- KeyError: If adata.varm[loadings_key] is missing.
- ValueError. If any requested PC is out of range, or if no genes are selected (e.g. all loadings are NaN).

## Examples

1) Signed loadings for PC1–PC3 (pos + neg)
```python
bk.pl.pca_loadings_heatmap(
    adata,
    pcs=(1, 2, 3),
    n_top=15,
    use_abs=False,
    show_negative=True,
)
```

2) Magnitude-only drivers across PCs
```python
bk.pl.pca_loadings_heatmap(
    adata,
    pcs=(1, 2, 3, 4),
    n_top=25,
    use_abs=True,
)
```

3) Emphasize per-gene PC preference (z-score across PCs)
```python
bk.pl.pca_loadings_heatmap(
    adata,
    pcs=(1, 2, 3),
    n_top=20,
    z_score=True,
)
```

4) No clustering (simple heatmap, preserves order)
```python
bk.pl.pca_loadings_heatmap(
    adata,
    pcs=(1, 2, 3),
    cluster_genes=False,
    cluster_pcs=False,
)
```

5) Use gene symbols (if available in adata.var)
```python
bk.pl.pca_loadings_heatmap(
    adata,
    pcs=(1, 2),
    gene_symbol_key="gene_symbol",
)
```

## Notes & tips
- Signed mode (use_abs=False) is best when you want to interpret opposing programs on a PC
(genes loading positive vs negative).
- Absolute mode (use_abs=True) is best to identify the strongest overall contributors.
- If you see repeated labels with gene_symbol_key, consider deduplicating upstream or switching
to adata.var_names to avoid ambiguity.
- For exporting per-PC gene sets (pos/neg/abs) for enrichment, use bk.tl.pca_loadings() instead.