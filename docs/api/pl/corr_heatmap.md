# Correlation heatmap

```{eval-rst}
.. autofunction:: bullkpy.pl.corr_heatmap

```

Correlation heatmap for **sample–sample** similarity (QC) or **gene–gene** similarity.

This function computes a correlation matrix from an expression matrix (typically a log-normalized layer)
and visualizes it as a heatmap. When `use="samples"` (default), it is a convenient QC plot to detect:
- batch structure
- outlier samples
- mislabeled groups
- strong technical gradients

When `use="genes"`, it produces a gene–gene correlation heatmap (useful for small gene panels).

## What it does

- Extracts a matrix from adata.layers[layer] (or adata.X if layer=None)
- Computes a correlation matrix using Pearson or Spearman
- Optionally subsets and orders samples using groupby / groups
- Optionally annotates columns/rows with color bars from one or more obs keys (col_colors)
- Plots either:
	•	a clustered heatmap (seaborn clustermap) if dendrogram=True
	•	a standard heatmap if dendrogram=False

Dependency: requires seaborn for plotting.

## Parameters

#### Core

**adata**   
AnnData object containing expression matrix and annotations.

**layer** (default: "log1p_cpm")    
Layer to use for correlations.   
Use log-scale normalized expression for best behavior.   

**method** (default: "pearson")    
Correlation metric:
	•	"pearson" – linear correlation (fast)
	•	"spearman" – rank correlation (more robust to non-linearity / outliers)

**use** (default: "samples")    
Which correlation matrix to compute:
	•	"samples" – sample × sample correlation (QC use-case)
	•	"genes" – gene × gene correlation (for smaller panels)

#### Subsetting / ordering (samples only)

**groupby**   
Categorical obs key used to subset and order samples.

**groups**    
Optional list of groups to keep (and define ordering).  
If provided, only those groups are used and ordered as given.  

#### Annotations (samples only)

**col_colors**    
One or more obs keys used to draw color annotations aligned to samples.  
Values are converted to categorical strings and mapped to colors automatically.  
Examples:  
	•	col_colors="batch"
	•	col_colors=["batch", "sex", "tumor_type"]

#### Plot appearance

**cmap, center, vmin, vmax**   
Control heatmap colormap and scaling.

**figsize**    
Figure size in inches. If None, chosen automatically based on matrix size.

**show_labels**   
Whether to show row/column labels (off by default for readability).

**dendrogram** (default: True)     
If True, performs hierarchical clustering and shows dendrograms.

**save**   
Path to save the plot.   

**show**   
Whether to display the plot.

## Returns

When dendrogram=True (default):
- cg — seaborn ClusterGrid object (contains cg.fig, cg.ax_heatmap, etc.)

When dendrogram=False:
- (fig, ax) — matplotlib Figure and Axes

Tip: if you want a consistent return type, keep dendrogram=True.

## Examples

1) Sample–sample QC correlation heatmap (default)
```python
bk.pl.corr_heatmap(adata)
```

2) Spearman correlation (robust to outliers)
```python
bk.pl.corr_heatmap(adata, method="spearman")
```

3) Subset/order samples by group and annotate batch
```python
bk.pl.corr_heatmap(
    adata,
    groupby="condition",
    groups=["control", "treated"],
    col_colors="batch",
    show_labels=False,
)
```

4) Multiple annotation bars (batch + clinical)
```python
bk.pl.corr_heatmap(
    adata,
    col_colors=["batch", "sex", "tumor_stage"],
)
```

5) Gene–gene correlation heatmap (small gene panel)
```python
bk.pl.corr_heatmap(
    adata,
    use="genes",
    method="pearson",
    show_labels=True,
)
```

For use="genes", consider filtering adata beforehand to a manageable number of genes.

6) Disable clustering (plain heatmap)
```python
bk.pl.corr_heatmap(
    adata,
    dendrogram=False,
    show_labels=False,
)
```

## Notes / recommendations
- Use log-scale expression (e.g. log1p_cpm) for meaningful correlations.
- Spearman can be preferable when expression distributions are heavy-tailed or have outliers.
- For large datasets, correlation matrices can become visually dense; consider:
	•	subsetting to relevant groups
	•	turning labels off (show_labels=False)
	•	using groupby to enforce interpretable ordering
	•	reducing to representative samples (or averaging within groups)

## See also
	•	pl.plot_corr_heatmap (for plotting a precomputed correlation matrix)
	•	tl.leiden_resolution_scan + pl.ari_resolution_heatmap (cluster-quality diagnostics)
	•	pl.categorical_confusion (compare categorical labels)







