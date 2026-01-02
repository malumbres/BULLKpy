# GSEA Leading edge heatmap

```{eval-rst}
.. autofunction:: bullkpy.pl.gsea_leading_edge_heatmap

```

Plot a heatmap of **expression values for leading-edge genes** from one or more GSEA pathways.

This plot is designed to answer: *“Are the genes driving enrichment (leading-edge) coherently up/down across my samples or groups?”*

It works by:
1) extracting leading-edge gene sets for selected terms,  
2) taking the **union** of those genes (optionally filtering by frequency),  
3) plotting expression across **samples** or **group means**, optionally **z-scoring per gene**.


## Parameters

#### Required

**adata**: anndata.AnnData    
Contains expression values (adata.X or adata.layers[layer]) and sample annotations (adata.obs).  

**pre_res**  
GSEApy prerank result object (e.g. returned from gseapy.prerank). Used to obtain leading-edge gene sets via _leading_edge_sets(...).  

#### Term selection

**term_idx**: optional   
Index-like selection understood by _leading_edge_sets. Useful to pick “top N” terms by rank.  

**terms**: Sequence[str] | None.   
Explicit list of pathway/term names to include.  

If neither is provided, the helper typically uses all terms (depending on _leading_edge_sets).  

#### Expression source

**layer**: str | None (default "log1p_cpm")   
Which matrix to plot:  
- If layer exists in adata.layers, uses adata.layers[layer]
- Otherwise falls back to adata.X

#### What to plot on rows

- **use**: "samples" | "group_mean" (default "samples")
- **"samples"**: rows are samples (adata.obs_names)
- **"group_mean"**: rows are group means; requires groupby
- **groupby**: str | None. Required when use="group_mean". Categorical adata.obs[groupby] defines the groups.

#### Gene set construction

**min_gene_freq**: int (default 1)  
Keep only genes that appear in at least min_gene_freq leading-edge sets among selected terms.   

**max_genes**: int | None (default 200)   
Cap the number of genes shown (after sorting). If None, no cap.  

**Gene prioritization rule**  
Genes are sorted by:  
	1.	decreasing frequency across selected pathways, then
	2.	alphabetical (tie-breaker)

**Important filtering**   
After selection, genes are filtered to those present in adata.var_names.  
If none remain, an error is raised.  

**Scaling / normalization for visualization**  

**z_score**: "row" | None (default "row")   
Despite the name, this implementation z-scores per gene across rows (i.e., z-score each column):  
- subtract mean and divide by std for each gene across samples/groups
- makes patterns comparable across genes with different dynamic ranges. 

**clip_z**: float | None (default 3.0)    
If provided, clamps z-scores into [-clip_z, +clip_z] to reduce outlier domination.  

#### Clustering / appearance

**row_cluster**: bool (default True)   
Cluster rows in seaborn clustermap (samples/groups).  

**col_cluster**: bool (default True).    
Cluster columns (genes).   

**cmap**: str (default "vlag").   
Colormap for heatmap (diverging recommended for z-scores).   

**figsize**: (w, h) | None.  
Auto-sized if None:   
- w = max(7.0, 0.10 * n_genes + 4.5)
- h = max(5.0, 0.12 * n_rows + 3.0).  

**show_labels**: bool (default False). 
Show gene labels on x-axis (can be slow/cluttered for many genes).   

**gene_label_fontsize**: float (default 7.0).   
Font size for gene labels if enabled.   

#### Output control
**save**: str | Path | None.  
Save figure via _savefig(...).   

**show**: bool (default True).  
Calls plt.show().  

## What it does

1. Dependency check. 
Requires seaborn (sns.clustermap). Raises ImportError if unavailable.  

2. Extract leading-edge gene sets. 

```python
term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)
```

le_sets is expected to be term -> set(genes).  

3. Build gene list
- Counts each gene’s frequency across selected leading-edge sets.
- Keeps genes with freq >= min_gene_freq.
- Sorts by (-freq, gene) and applies max_genes cap.
- Filters to genes present in adata.var_names.   

4. Build expression matrix. 
- Subsets expression to selected genes.
- If use="samples" → DataFrame rows = samples
- If use="group_mean" → DataFrame rows = group categories, values = mean expression within each group. 
	
5. Optional per-gene z-scoring. 
If z_score == "row", standardizes each gene across rows and optionally clips.  

6. Plot clustermap.  
- Uses sns.clustermap(df, ...)
- Colorbar label is "Z-score" if z-scored, else "Expression".

## Returns
- **df**: pd.DataFrame. 
The matrix used for plotting:  
	-- rows: samples or groups
	-- columns: leading-edge genes
	-- values: expression or z-scores. 
- **g**: seaborn ClusterGrid. 
The clustermap object (for further customization/export).  

## Raises
- ImportError if seaborn is not installed.
- ValueError if:
	-- no genes pass min_gene_freq
	-- none of the selected leading-edge genes are present in adata.var_names
	-- invalid use or z_score value
	-- use="group_mean" but groupby is missing. 
- KeyError if groupby not found in adata.obs (when required). 

## Notes / tips
- For interpretability, z-scoring is usually recommended (z_score="row"), especially when genes have different baselines.
- If you want to emphasize absolute expression rather than relative patterns, set z_score=None.
- For many genes, use:
	•	max_genes=50..150
	•	show_labels=False
	•	or increase figsize

## Examples

1) Heatmap across samples (default)
```python
df, g = bk.pl.gsea_leading_edge_heatmap(adata, pre_res, terms=[
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
])
```

2) Heatmap of group means (e.g., subtype averages)
```python
df, g = bk.pl.gsea_leading_edge_heatmap(
    adata, pre_res,
    terms=["HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE"],
    use="group_mean",
    groupby="Subtype",
    min_gene_freq=2,
    max_genes=120,
    z_score="row",
)
```

3) No z-scoring, show gene labels, save figure
```python
df, g = bk.pl.gsea_leading_edge_heatmap(
    adata, pre_res,
    term_idx=range(20),
    z_score=None,
    show_labels=True,
    gene_label_fontsize=6,
    save="leading_edge_expr_heatmap.png",
)
```