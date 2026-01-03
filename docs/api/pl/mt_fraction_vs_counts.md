# Mitochondrial fraction vs. counts

```{eval-rst}
.. autofunction:: bullkpy.pl.mt_fraction_vs_counts

```

Scatter plot for **QC inspection of mitochondrial fraction vs library size**.

This plot is commonly used to identify samples with **high mitochondrial
content**, which often indicates RNA degradation, low-quality input, or
technical artifacts in bulk RNA-seq (and sample-level single-cell summaries). 

## What it does

Plots a scatter of samples:  
- x-axis: library size (default: total_counts)
- y-axis: mitochondrial fraction (default: pct_counts_mt). 

Optionally:  
- Colors samples by a categorical variable
- Applies QC thresholds for counts and/or mitochondrial fraction
- Highlights QC-failing samples. 

Automatically reports the number of QC failures in the title when
thresholds are used. 

## Requirements

adata.obs must contain:  
- the column specified by x (default: "total_counts")
- the column specified by y (default: "pct_counts_mt").  

If groupby is provided, that column must also exist in adata.obs.   

## Parameters

#### Axes

**x** (str, default "total_counts"): Column in adata.obs used for the x-axis (library size).  
**y** (str, default "pct_counts_mt"): Column in adata.obs used for the y-axis (mitochondrial fraction, typically in %).

#### Grouping

**groupby** (str | None, default None). 
Optional categorical column in adata.obs used to color samples. 
(e.g. batch, subtype).  

#### QC thresholds

Samples are considered QC-pass only if all specified conditions are met.  

**min_counts** (float | None): Minimum allowed value for x.  
**max_counts** (float | None): Maximum allowed value for x.   
**min_mt** (float | None): Minimum allowed value for y (rarely used).  

**max_mt** (float | None):  Maximum allowed mitochondrial fraction (commonly used).  

When any threshold is provided:  
- QC-failing samples are visually marked
- The plot title includes the number of failing samples

#### Scaling
**logx** (bool, default True): Apply log scaling to the x-axis (recommended for library size).  
**logy** (bool, default False): Apply log scaling to the y-axis (usually unnecessary for percentages).  

#### Figure and output

**figsize** (tuple[float, float], default (5.5, 4.5)). 
Figure size in inches.  

**save** (str | Path | None, default None). 
If provided, saves the figure to this path via _savefig.  

**show** (bool, default True)
If True, calls plt.show().

## Returns

- **fig** (matplotlib.figure.Figure). The created figure.  
- **ax** (matplotlib.axes.Axes). The scatter plot axis.  

## Interpretation tips

**Upper region (high pct_counts_mt)**. 
Samples dominated by mitochondrial reads → often low quality.   

**Low counts + high mt%**  
Strong indicator of degraded or failed libraries.  

**Typical bulk RNA-seq cutoffs**.  
- pct_counts_mt < 5–10% (context-dependent)
- Combined with minimum total_counts

This plot is typically used before filtering, to define reasonable
mitochondrial thresholds.  

## Examples

Basic mitochondrial QC
```python
bk.pl.mt_fraction_vs_counts(adata)
```

Apply mitochondrial cutoff
```python
bk.pl.mt_fraction_vs_counts(
    adata,
    max_mt=10.0,
)
```

Combine with library size threshold
```python
bk.pl.mt_fraction_vs_counts(
    adata,
    min_counts=1e6,
    max_mt=8.0,
)
```

Color by subtype
```python
bk.pl.mt_fraction_vs_counts(
    adata,
    groupby="Subtype",
    max_mt=10.0,
)
```

Save to file
```python
bk.pl.mt_fraction_vs_counts(
    adata,
    max_mt=7.5,
    save="mt_fraction_qc.png",
    show=False,
)
```