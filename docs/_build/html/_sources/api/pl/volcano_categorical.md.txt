# Volcano for categorial comparisons

```{eval-rst}
.. autofunction:: bullkpy.pl.volcano_categorical

```

Volcano plot for categorical associations.

This function visualizes effect sizes (e.g. log2 fold change or regression
coefficients) against statistical significance for categorical comparisons,
highlighting significant features and optionally labeling top hits.


## What it does  
	•	Creates a volcano plot with:
	•	x-axis: effect size (e.g. log2FC, effect_size)
	•	y-axis: -log10(q-value)
	•	Highlights statistically significant features
	•	Optionally draws effect-size and significance thresholds
	•	Labels the top significant genes
	•	Returns (fig, ax) for further customization

## Expected input

results must be a pandas.DataFrame containing at least:
	•	gene — feature / gene name
	•	q_col — multiple-testing corrected p-values (default: qval)
	•	x — effect size (e.g. log2FC, regression coefficient)


## Parameters. 

**results**  
DataFrame with association results.

**x**   
Column name containing effect sizes.  
Default: "effect_size".  

**q_col**  
Column name containing q-values.  
Default: "qval".  

**label_top**   
Number of top (lowest q-value) genes to label.  
Default: 12.  

**q_thr**   
Significance threshold on q-values.  
Default: 0.05.  

**effect_thr**  
Optional effect size threshold. If provided, vertical dashed lines are drawn
at ±effect_thr.

**figsize**   
Figure size in inches (width, height).  

**title**   
Optional plot title.   

**save**    
File path to save the figure (e.g. "volcano.png").  

**show**   
Whether to display the plot using plt.show().

## Returns
	•	fig — matplotlib.figure.Figure
	•	ax — matplotlib.axes.Axes

## Example  

#### Basic volcano plot

```python
fig, ax = bk.pl.volcano_categorical(
    results=df_results,
    x="log2FC",
    q_col="qval",
)
```
#### With thresholds and labels
```python
fig, ax = bk.pl.volcano_categorical(
    results=df_results,
    x="effect_size",
    q_col="qval",
    q_thr=0.01,
    effect_thr=1.0,
    label_top=20,
    title="Treatment vs Control",
)
```

#### Save without displaying

```python
bk.pl.volcano_categorical(
    results=df_results,
    save="volcano_plot.pdf",
    show=False,
)
```

## Notes. 
	•	Infinite and missing values are automatically removed.
	•	The y-axis uses a clipped minimum q-value (1e-300) to avoid numerical issues.
	•	For continuous covariates, consider regression-based effect sizes rather than
fold changes.

## See also
	•	pl.volcano_continuous
	•	tl.differential_expression
	•	pl.rankplot
