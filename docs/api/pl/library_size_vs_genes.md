# Library size versus genes

```{eval-rst}
.. autofunction:: bullkpy.pl.library_size_vs_genes

```

QC scatter plot of **library size** versus **number of detected genes** per sample, with optional threshold lines and outlier highlighting.

This is a bulk/sample QC helper: it reads two numeric columns from `adata.obs` (by default `total_counts` and `n_genes_detected`), optionally colors points by a categorical annotation, and marks samples that fail any user-provided thresholds.

## What it plots
- **x-axis**: adata.obs[x] (default: "total_counts"). 
- **y-axis**: adata.obs[y] (default: "n_genes_detected"). 
- Optionally:
	-- color/legend by adata.obs[groupby]
	-- log-scale for either axis (logx, logy)
	-- dashed threshold lines for min/max bounds
	-- outliers (QC fails) overplotted with a distinct marker

**Required adata.obs columns**  
- x must exist in adata.obs and be coercible to numeric. 
- y must exist in adata.obs and be coercible to numeric. 

If groupby is provided, it must exist in adata.obs (treated as categorical via .astype(str)).

**QC thresholding logic**.  

A sample is considered passing (inlier) if it satisfies all provided bounds:  
- min_counts: x >= min_counts
- max_counts: x <= max_counts
- min_genes:  y >= min_genes
- max_genes:  y <= max_genes

If a threshold is None, that bound is ignored.  

Samples failing any bound are outliers (QC fail) and can be highlighted if show_outliers=True.  

## Parameters

#### Data selection
**x**: str (default "total_counts")  
obs column for library size / total counts.  

**y**: str (default "n_genes_detected").   
obs column for detected genes.   

**groupby**: str | None.  
If set, points are plotted per category with a stable legend.   

#### Thresholds (optional). 

**min_counts, max_counts**: float | None.   
Lower/upper bounds on x.   

**min_genes, max_genes**: float | None.  
Lower/upper bounds on y.   

#### Display / styling. 

**logx, logy**: bool.   
Apply logarithmic scaling to axes.   

**s**: float.   
Marker size for inliers.    

**alpha**: float.  
Marker transparency for inliers.   

**linewidth**: float.  
Edge line width for inliers.   

**edgecolor**: str.  
Edge color for inliers (matplotlib color spec).   

#### Outlier overlay. 

**show_outliers**: bool.  
If True, overlays failed samples on top.   

**outlier_color**: str. 

####Color for outliers.  

** outlier_marker**: str. 
Marker style for outliers (default "x").  

**outlier_size**: float. 
Marker size for outliers.  

#### Labels and layout 

**title**: str | None. 
If None, auto-generates:  
- with any threshold set: "{y} vs {x} (QC fail: {n_fail})"
- otherwise: "{y} vs {x}".   

**xlabel, ylabel**: str | None. 
If None, defaults to x and y.  

**figsize**: (float, float). 
Figure size in inches.  

**legend**: bool. 
Show legend when groupby is provided.  

**legend_loc**: str. 
Matplotlib legend location string (e.g. "best", "upper right").  

#### Output. 

**save**: str | Path | None. 
If provided, saves the figure via _savefig(fig, save).  

**show**: bool. 
If True, calls plt.show().  

## Returns
- **fig**: matplotlib.figure.Figure
- **ax**: matplotlib.axes.Axes

## Raises
- KeyError if x, y, or groupby (when provided) are missing from adata.obs.
- Drops rows where x or y cannot be coerced to numeric (via pd.to_numeric(..., errors="coerce")).

## Examples

1) Basic QC scatter
```python
bk.pl.library_size_vs_genes(adata)
```

2) Add thresholds and color by subtype
```python
bk.pl.library_size_vs_genes(
    adata,
    min_counts=1e6,
    min_genes=12000,
    groupby="Subtype",
)
```

3) Linear axes (no log scaling)
```python
bk.pl.library_size_vs_genes(adata, logx=False, logy=False)
```

4) Save to file
```python
bk.pl.library_size_vs_genes(
    adata,
    min_counts=5e5,
    min_genes=8000,
    save="qc_library_vs_genes.png",
    show=False,
)
```