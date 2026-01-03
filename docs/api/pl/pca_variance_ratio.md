# PCA variance ratio

```{eval-rst}
.. autofunction:: bullkpy.pl.pca_variance_ratio

```

Scree plot for **PCA explained variance**, with an optional cumulative variance curve.
This is the standard diagnostic plot used to decide how many principal components
to retain.


```{figure} /_static/pca_svariance_ratio.png
:alt: PCA variance ratio
:width: 500px
:align: center
```

Example PCA variance ratio

## What it does

Reads per-PC explained variance ratios from adata.uns[key]["variance_ratio"] (as written by bk.tl.pca).  

Plots a scree plot:
- x-axis: principal component index (PC1, PC2, …)
- y-axis: explained variance ratio.  

Optionally overlays a cumulative variance curve on a secondary y-axis.  

## Requirements

PCA must have been run beforehand:  
```python
bk.tl.pca(adata)
```

The AnnData object must contain:
```python
adata.uns[key]["variance_ratio"]
```

## Parameters

#### Data source

**key** (str, default "pca").  
Key in adata.uns where PCA metadata is stored.   
The function expects adata.uns[key]["variance_ratio"].   

**n_comps** (int | None, default None). 
If provided, plot only the first n_comps principal components.  

#### Plot options

**cumulative** (bool, default True).  
If True, adds a cumulative explained variance curve on a secondary y-axis.   

**figsize** (tuple[float, float], default (6.5, 4.5)).  
Figure size in inches.  

#### Output

**save** (str | Path | None, default None). 
If provided, saves the figure to this path via _savefig.  

**show** (bool, default True).  
If True, displays the plot with plt.show().   

## Returns
- **fig** (matplotlib.figure.Figure)

- **ax** (matplotlib.axes.Axes). 

The primary axes (explained variance ratio).  
If cumulative=True, the cumulative curve is drawn on a secondary y-axis
created via ax.twinx().  

## Interpretation
- The scree plot (variance ratio per PC) shows how much variance each PC explains.
- The cumulative curve helps choose a cutoff:  
 e.g. retain PCs until cumulative variance ≥ 0.8 or 0.9. 
- A sharp “elbow” in the scree plot often indicates a natural dimensionality.

## Examples

1) Default scree plot with cumulative variance
```python
bk.pl.pca_variance_ratio(adata)
```

2) Show only the first 20 PCs
```python
bk.pl.pca_variance_ratio(adata, n_comps=20)
```

3) Scree plot without cumulative curve
```python
bk.pl.pca_variance_ratio(adata, cumulative=False)
```

4) Save to file
```python
bk.pl.pca_variance_ratio(adata, save="pca_variance_ratio.pdf")
```

