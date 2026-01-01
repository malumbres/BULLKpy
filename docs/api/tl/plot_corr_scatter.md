# Plot Correlation scatter

```{eval-rst}
.. autofunction:: bullkpy.tl.plot_corr_scatter

```

Scatter plot for visualizing the relationship between two variables together
with their correlation statistic.

This function supports **obs–obs** and **gene–obs** correlations, optional
batch/covariate adjustment, and coloring points by a third variable.

## What it does
	1.	Extracts vectors x and y from adata depending on kind
	2.	Computes a correlation (optionally batch- / covariate-aware)
	3.	Draws a scatter plot
	4.	Displays the correlation coefficient, p-value, and sample size in the title

This is intended as a diagnostic / exploratory visualization to accompany
the numerical correlation utilities in bullkpy.tl.

## Parameters

**adata**   
AnnData object containing expression data and/or obs annotations.

**x, y**   
Variables to plot:
	•	For kind="obs_obs": both must be numeric columns in adata.obs
	•	For kind="gene_obs": one must be a gene name (adata.var_names),
the other a numeric obs column

**hue**   
Optional obs column to color points by.
	•	Numeric → continuous color scale
	•	Categorical → discrete colors + legend

**kind**   
Type of relationship:
	•	"obs_obs" (default): obs vs obs
	•	"gene_obs": gene vs obs (order does not matter)

**layer**   
Expression layer to use when a gene is involved (default: "log1p_cpm").

**method**   
Correlation method: "pearson" or "spearman".

**batch_key, batch_mode, covariates**   
Optional batch-aware / partial correlation settings.  
These are forwarded to the internal correlation engine.  

**figsize**   
Figure size in inches.

**s**   
Marker size.

**alpha**  
Marker transparency.

**cmap**  
Colormap for numeric hue.

**palette**   
Color palette for categorical hue.

**title**   
Custom title. If None, a title with correlation statistics is generated.

**show**   
Whether to call plt.show().

## Returns

```python
(fig, ax)
```
	•	fig — matplotlib Figure
	•	ax — matplotlib Axes

This allows further customization or saving by the caller.

## Examples

1) Obs–obs correlation

```python
bk.pl.plot_corr_scatter(
    adata,
    x="libsize",
    y="pct_counts_mt",
    kind="obs_obs",
)
```

2) Gene–obs correlation

```python
bk.pl.plot_corr_scatter(
    adata,
    x="PDCD1",
    y="signature_Tcell",
    kind="gene_obs",
)
```

3) Coloring by a categorical variable

```python
bk.pl.plot_corr_scatter(
    adata,
    x="TP53",
    y="MDM2",
    kind="gene_obs",
    hue="Batch",
)
```

4) Batch-aware correlation visualization

```python
bk.pl.plot_corr_scatter(
    adata,
    x="GZMB",
    y="cytotoxic_score",
    kind="gene_obs",
    batch_key="Batch",
    batch_mode="within",
)
```

## Notes
	•	Correlation statistics are computed before plotting, using the same
batch/covariate logic as the numerical correlation functions.
	•	Missing values are handled internally by the correlation routine.
	•	For large datasets, consider subsampling before plotting to improve
rendering performance.

## See also
	•	tl.top_gene_obs_correlations
	•	tl.gene_gene_correlations
	•	tl.partial_corr
	•	pl.rankplot

