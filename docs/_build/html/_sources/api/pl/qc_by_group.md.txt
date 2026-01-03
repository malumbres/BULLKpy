# QC by group

```{eval-rst}
.. autofunction:: bullkpy.pl.qc_by_group

```

Grouped **QC metric distribution plots** for bulk RNA-seq (or pseudo-bulk) samples.

This function visualizes how standard QC metrics vary **across levels of a
categorical metadata variable** (e.g. batch, cohort, condition), making it
easy to detect **systematic quality differences** between groups.

```{figure} /_static/qc_by_group.png
:alt: QC by group
:width: 700px
:align: center
```
Example QC by group

## What it does

For each QC metric in keys, this function:  
1. Groups samples by adata.obs[groupby].  
2. Plots one panel per metric.   
3. Shows the distribution per group as either:   
- violins (default), or
- boxplots.  
4. Optionally log-transforms selected metrics.  
5. Annotates group labels with sample counts.  

This is especially useful for:
- Batch effect diagnostics
- Comparing cohorts or experimental conditions
- Detecting QC-driven confounding before downstream analysis

## Requirements
- groupby must be a categorical column in adata.obs
- Each entry in keys must exist in adata.obs. 

Typically, these metrics are created by:  
```python
bk.pp.qc_metrics(adata)
```

## Parameters

#### Grouping

**groupby ** (str, required). 
Categorical column in adata.obs used to group samples
(e.g. "Batch", "Cohort", "Platform").  

#### QC metrics
**keys** (Sequence[str], default.  
("total_counts", "n_genes_detected", "pct_counts_mt", "pct_counts_ribo"))
QC metrics to plot, one panel per metric.   

Each key must exist in adata.obs.  

#### Plot type

**kind** ("violin" | "box", default "violin").  
Type of distribution plot:   
- "violin" : shows full distribution shape + median
- "box" →: compact summary (quartiles, median). 

#### Transformations

**log1p** (Sequence[str], default ("total_counts",)). 
Metrics that should be transformed using log1p before plotting.  

Useful for highly skewed variables such as library size.  

#### Layout and labels

**figsize** (tuple[float, float], default (11, 4)). 
Overall figure size.  

**rotate_xticks** (int, default 45). 
Rotation angle for group labels.  

**show_n** (bool, default True). 
If True, appends sample counts to group labels
(e.g. Batch1 (n=24)).  

#### Output

**save** (str | Path | None): Path to save the figure.
**show** (bool, default True): If True, calls plt.show().

## Returns
- **fig** (matplotlib.figure.Figure). The figure object.  
- **axes** (list[matplotlib.axes.Axes]). One axis per QC metric, in the same order as keys.  

## Interpretation guide

Shifted distributions between groups: potential batch or cohort effects.  

Higher mt% or lower gene counts in a group: degraded RNA or sample preparation issues.  

Broader distributions in one group: increased technical variability.  

These plots are most informative before filtering, to guide threshold
selection or batch-aware QC decisions.  

## Examples

Default QC comparison by batch
```python
bk.pl.qc_by_group(adata, groupby="Batch")
```

Boxplots instead of violins
```python
bk.pl.qc_by_group(
    adata,
    groupby="Cohort",
    kind="box",
)
```

Custom metrics and transformations
```python
bk.pl.qc_by_group(
    adata,
    groupby="Platform",
    keys=("total_counts", "pct_counts_mt"),
    log1p=("total_counts",),
)
```

Save without displaying
```python
bk.pl.qc_by_group(
    adata,
    groupby="Batch",
    save="qc_by_batch.png",
    show=False,
)
```

## Related functions
	•	qc_scatter_panel
	•	qc_metrics
	•	library_size_vs_genes
	•	mt_fraction_vs_counts
