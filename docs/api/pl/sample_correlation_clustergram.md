# Sample correlation clustergram
```{eval-rst}
.. autofunction:: bullkpy.pl.sample_correlation_clustergram

```

Sample–sample **correlation clustergram** for bulk expression QC and exploratory analysis.  

This plot is often more interpretable than raw distance heatmaps because:  
- the **heatmap shows correlations directly** (range `[-1, 1]`), and
- hierarchical clustering is performed on **distance = `1 − correlation`**.


## What it does. 

1. Extracts a sample × gene matrix
```python 
X = _get_matrix(adata, layer=layer, use="samples")
```
2. Optional Spearman transformation.   
If method="spearman":   
- Each gene is ranked across samples.
- Pearson correlation is then computed on the ranked data.   
This matches the standard definition of Spearman correlation.  

3. Computes the sample–sample correlation matrix
```python 
C = np.corrcoef(X)   # shape: (n_samples, n_samples)
```
Values range from -1 (anti-correlated) to +1 (perfectly correlated).  
  
4. Builds a clustering distance
```python 
distance = 1.0 - correlation
```
- Converted to condensed form with squareform.
- Hierarchical clustering is performed with scipy.cluster.hierarchy.linkage.  

5. Plots a seaborn clustergram.  
- Heatmap values = correlation
- Color map = "vlag" (diverging, centered at 0)
- Dendrograms reflect clustering on 1 − correlation.   

6. Optionally annotates samples with metadata.   
- col_colors adds color bars for adata.obs columns.
- Legends are drawn manually to the right of the plot.

## Parameters

#### Core computation

**adata** (AnnData): Input object.  

**layer** (str | None, default "log1p_cpm"): Expression layer used for correlations. Passed to _get_matrix.  

**method** ("pearson" | "spearman", default "pearson").  
Correlation type:   
- "pearson": linear correlation.
- "spearman": rank-based correlation (robust to outliers, monotonic trends).  

**linkage_method** (str, default "average"):  Linkage method used for hierarchical clustering.  
Common choices: "average", "complete", "single".  

#### Metadata annotations

**col_colors** (Sequence[str] | None).   
List of adata.obs columns to show as color annotations above the heatmap.    
 Only categorical metadata are meaningful here.   

**palette** (str, default "tab20"): Palette used for mapping metadata categories to colors.

#### Display

**figsize** ((w, h) | None): If None, auto-sized based on number of samples:
```python
w = max(6.0, min(16.0, 0.18 * n_samples + 4.0))
figsize = (w, w)
```
show_labels (bool, default False).  Whether to show sample names on axes.  

### Output
**save** (str | Path | None): If provided, saves the figure using _savefig.    
**show (bool, default True). Whether to display the figure with plt.show().   

## Returns
**cg**: seaborn.matrix.ClusterGrid. 
- Main heatmap axis: cg.ax_heatmap
- Figure: cg.fig

## Requirements
- seaborn
- scipy (pdist, squareform, linkage)
- Raises ImportError if dependencies are missing.

## Interpretation guide
- Red (positive) values → samples with similar expression profiles.
- Blue (negative) values → anti-correlated samples.
- Block-diagonal structure → coherent sample groups (often biological subtypes).
- Mixed blocks or striping → batch effects or gradual expression gradients.

## Best practices
For bulk RNA-seq QC, this plot is often preferable to distance heatmaps.   
Use:  
- method="pearson" for general similarity.
- method="spearman" when outliers or non-linear monotonic trends are suspected.  

Combine with metadata annotations:
```python
col_colors=["Subtype", "Batch"]
```
to visually assess confounding effects.

## Examples

1) Pearson correlation QC
```python
bk.pl.sample_correlation_clustergram(
    adata,
    layer="log1p_cpm",
    method="pearson",
    col_colors=["Subtype", "Batch"],
)
```

2) Spearman correlation (robust)
```python
bk.pl.sample_correlation_clustergram(
    adata,
    layer="log1p_cpm",
    method="spearman",
    col_colors=["Patient"],
)
```

3) Compact overview (no labels)
```python
bk.pl.sample_correlation_clustergram(
    adata,
    show_labels=False,
)
```