# Dotplot association

```{eval-rst}
.. autofunction:: bullkpy.pl.dotplot_association

```

Scanpy-like dot plot for association results across multiple contrasts or groupings.

This visualization summarizes **association analyses run across multiple
groupby variables or contrasts**, encoding:

- **dot color** → effect size  
- **dot size** → statistical significance (−log10 q-value)

It is especially useful when you have **multiple categorical association
runs** (e.g. different `groupby` variables, contrasts, or experimental factors)
and want to compare them side by side.

## What it does
- Produces a feature × contrast dot plot
- Mimics Scanpy’s dotplot semantics:
	•	color = effect size
	•	size = significance
- Automatically selects the top N most significant features per contrast
- Works with:
	•	gene–categorical associations
	•	obs–categorical associations
	•	mixed results stored in a single tidy table

This plot is ideal for high-level comparative summaries and figure panels.

## Expected input format

A tidy pandas DataFrame with (at minimum):

| column | description |
| ------------ | ------------------ |
| feature_col | gene or obs feature name |
| groupby_col | contrast / groupby identifier |
| effect_col | signed effect size |
| q_col | adjusted p-value (q-value) |

The exact column names are configurable.

## Parameters

**df**   
Association results table (long format).

**feature_col**   
Column identifying features (e.g. "gene" or "obs").

**groupby_col**   
Column identifying contrasts or groupby runs.
	
**effect_col**    
Column used for dot color (signed effect size).

**q_col**   
Column used for dot size (−log10(q)).

**top_n**   
Number of most significant features shown per groupby/contrast.

**figsize**   
Figure size in inches. If None, computed automatically.

**cmap**   
Colormap for effect sizes (default: "RdBu_r").

**vmin, vmax**   
Color scale limits. Defaults inferred from data.

**size_min, size_max**   
Minimum and maximum dot sizes.

**title**    
Optional plot title.

**save**    
Path to save the figure.

**show**   
Whether to display the plot immediately.

## Returns

```python
(fig, ax)
```
- **fig**: matplotlib Figure
- **ax**: matplotlib Axes.  

## Examples

1) Dotplot across multiple categorical associations

```python
bk.pl.dotplot_association(
    df=assoc_df,
    feature_col="gene",
    groupby_col="groupby",
    effect_col="effect",
)
```

2) Limit to top 20 features per contrast
```python
bk.pl.dotplot_association(
    df=assoc_df,
    feature_col="gene",
    top_n=20,
)
```

3) Custom color scale and dot sizes
```python
bk.pl.dotplot_association(
    df=assoc_df,
    feature_col="gene",
    cmap="coolwarm",
    size_min=20,
    size_max=300,
)
```

4) Save without displaying
```python
bk.pl.dotplot_association(
    df=assoc_df,
    feature_col="gene",
    save="association_dotplot.pdf",
    show=False,
)
```

## Interpretation guide
- Large, dark red dots: strong positive association with high confidence
- Large, dark blue dots: strong negative association with high confidence
- Small dots: weak or non-significant associations

Rows are ordered by feature name; columns correspond to contrasts/groupby runs.

## Notes
- Dot sizes are scaled globally across the plot for comparability.
- The plot works even with a single contrast, but is most informative with ≥2.
- Effect sizes should be signed for meaningful color interpretation.

## See also
	•	pl.rankplot_association
	•	pl.rankplot
	•	pl.volcano
	•	pl.ma
	•	tl.gene_categorical_association
	•	tl.obs_categorical_association







