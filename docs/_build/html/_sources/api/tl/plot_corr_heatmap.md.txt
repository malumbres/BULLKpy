# Plot Correlation heatmap

```{eval-rst}
.. autofunction:: bullkpy.tl.plot_corr_scatter

```

Heatmap visualization for correlation matrices.

This function is designed to visualize **wide correlation matrices**
(e.g. obs × obs, gene × gene) produced by functions such as
`obs_obs_corr_matrix`, `gene_gene_correlations`, or similar utilities.

## What it does
	•	Takes a wide pandas DataFrame (rows × columns)
	•	Displays a heatmap of correlation values
	•	Uses seaborn if available, otherwise falls back to pure matplotlib
	•	Supports optional value annotations and color scaling

This function is intentionally lightweight and does not recompute
correlations—it only visualizes an existing matrix.

## Parameters

**mat**   
A pandas DataFrame where rows and columns correspond to variables
and values are correlation coefficients.

**figsize**   
Size of the figure in inches.

**cmap**   
Colormap used for the heatmap (default: "vlag").

**center**   
Value at which to center the colormap (default: 0.0,
appropriate for correlations).

**vmin, vmax**   
Optional minimum and maximum values for the color scale.  
Useful to enforce symmetric limits (e.g. -1, 1).  

**annot**.  
If True, annotate each cell with its numeric value.

**fmt**   
String formatting for annotations (e.g. ".2f").

**title**   
Optional title for the plot.

**show**      
Whether to immediately display the figure using plt.show().

## Returns

```python
(fig, ax)
```
	•	fig — matplotlib Figure
	•	ax — matplotlib Axes

Returning the objects allows further customization or saving by the caller.

## Examples

1) Heatmap from an obs–obs correlation matrix

```python
R = bk.tl.obs_obs_corr_matrix(
    adata,
    focus=["libsize", "pct_counts_mt", "pct_counts_ribo"],
)

bk.pl.plot_corr_heatmap(
    R,
    title="QC metric correlations",
)
```

2) Fixed correlation scale and annotations

```python
bk.pl.plot_corr_heatmap(
    R,
    vmin=-1,
    vmax=1,
    annot=True,
    fmt=".2f",
    title="Correlation matrix (annotated)",
)
```

3) Large matrices (disable annotations)

```python
bk.pl.plot_corr_heatmap(
    R,
    figsize=(10, 8),
    annot=False,
)
```

## Notes
	•	If seaborn is installed, it is used for nicer defaults and labels.
	•	If seaborn is not available, the function gracefully falls back to
matplotlib.imshow.
	•	For very large matrices, consider:
	•	disabling annot
	•	pre-filtering rows/columns
	•	clustering externally before plotting

## See also
	•	tl.obs_obs_corr_matrix
	•	tl.gene_gene_correlations
	•	pl.plot_corr_scatter