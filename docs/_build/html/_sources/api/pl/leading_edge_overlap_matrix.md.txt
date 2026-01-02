# Leading edge overlap matrix

```{eval-rst}
.. autofunction:: bullkpy.pl.leading_edge_overlap_matrix

```

Build and plot a **pathway × gene binary matrix** indicating **leading-edge membership** from a GSEApy prerank result. This is meant to mimic a “Leading Edge Viewer”: you can quickly see whether a small set of genes drives many enriched pathways.

- **Rows** = pathways/terms   
- **Columns** = genes   
- **Cell value** = `1` if the gene is in the pathway’s leading-edge set, else `0`   
- Optionally filters genes by how often they appear across leading edges (`min_gene_freq`), and clusters rows/columns via seaborn `clustermap`.


```{figure} /_static/leading_edge_overlap_matrix.png
:alt: Violin plot example
:width: 1000px
:align: center
```

Example Leading edge overlap matrix

## Parameters

#### Required

**pre_res**  
A GSEApy prerank result object (the pre_res returned by gseapy.prerank(...)).  
This function expects that pre_res contains enough information to recover each term’s leading-edge genes, via the internal helper _leading_edge_sets(pre_res, ...).  

**Term selection**  

**term_idx**: None or index-like   
Optional selection mechanism understood by _leading_edge_sets. Use this when you want to select terms by index/rank rather than by name.  

**terms**: Sequence[str] | None    
Explicit list of term names to include (subset of available terms in pre_res).  

If both are provided, the helper decides precedence (typically: explicit terms wins; otherwise term_idx; otherwise all).  

#### Gene filtering / ordering

**min_gene_freq**: int (default 2)    
Keep only genes that appear in at least this many leading-edge sets across the selected terms.
- 1 keeps everything
- larger values focus on shared drivers

**sort_genes_by**: "freq" or "alpha"   
- "freq" (default): sort columns by decreasing gene frequency across pathways
- "alpha": sort columns alphabetically. 

#### Clustering / plotting

**row_cluster**: bool (default True)  
Cluster pathways (rows).  

**col_cluster**: bool (default False)    
Cluster genes (columns). Often disabled because there may be many genes.  

**cmap**: colormap (default "Greys")   
Binary heatmap color scale.  

**figsize**: (w, h) | None    
If None, size is chosen automatically based on matrix dimensions.   

**show_gene_labels**: bool (default True)   
Show gene names on the x-axis.  

**gene_label_fontsize**: float (default 8.0)   
Font size for gene labels (use smaller values for large matrices).   

**show_term_labels**: bool (default True)   
Show pathway names on the y-axis.   

**save**: str | Path | None. 
If provided, saves the figure (not just the matrix) using _savefig.   

**show**: bool (default True).   
Calls plt.show().   

## What it does. 
1. Validates dependencies
- Requires seaborn (uses sns.clustermap).
- Requires min_gene_freq >= 1.  

2. Extracts leading-edge sets. 
Calls:
```python
term_names, le_sets = _leading_edge_sets(pre_res, term_idx=term_idx, terms=terms)
```
le_sets is expected to be a dict: term -> set(genes).

3. Builds the binary matrix. 
- Creates union of all leading-edge genes across selected terms.
- Constructs a matrix mat[i, j] = 1 if gene j is in term i’s leading-edge set.  

4. Filters genes by frequency. 
- Computes gene_freq = df.sum(axis=0).
- Keeps only genes where gene_freq >= min_gene_freq.  

5. Sorts columns  
Either by frequency (descending) or alphabetically.  

6. Plots with seaborn clustermap. 
- Optional clustering for rows/cols.
- Removes the colorbar (cbar_pos=None) because it’s binary.
- Adds axis labels and title.
- If gene labels are shown, rotates them via _rotate_gene_labels and increases bottom margin.  

7. Optionally saves and shows. 

## Returns
- **df**: pd.DataFrame. 
The final binary pathway × gene matrix after filtering and sorting.  

- **g**: seaborn ClusterGrid.  
The object returned by sns.clustermap (useful for fine-grained figure edits).   

## Raises
- ImportError if seaborn is missing.
- ValueError if:
	-- min_gene_freq < 1
	-- no leading-edge genes are found for the selected terms
	-- filtering removes all genes (df.shape[1] == 0)
	-- invalid sort_genes_by
- Potential errors from _leading_edge_sets if terms cannot be resolved in pre_res.  

## Notes / tips
- Start with min_gene_freq=2 to highlight “shared driver” genes across enriched pathways.
- If you have many terms, consider:
	-- selecting a smaller subset via terms=[...]
	-- disabling gene labels (show_gene_labels=False)
	-- increasing min_gene_freq
- If you want clustering among genes, set col_cluster=True, but this can be slow for large matrices.

## Examples

1) Default: show shared leading-edge genes across all enriched terms
```python
df_le, g = bk.pl.leading_edge_overlap_matrix(pre_res)
```

2) Focus on specific Hallmark pathways
```python
terms = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1",
]
df_le, g = bk.pl.leading_edge_overlap_matrix(
    pre_res,
    terms=terms,
    min_gene_freq=1,
    col_cluster=True,
)
```

3) Emphasize only genes shared across many pathways
```python
df_le, g = bk.pl.leading_edge_overlap_matrix(
    pre_res,
    min_gene_freq=4,
    show_gene_labels=False,
)
```
