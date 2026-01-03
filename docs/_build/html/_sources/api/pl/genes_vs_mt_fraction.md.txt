# Genes vs. mitochondrial fraction

```{eval-rst}
.. autofunction:: bullkpy.pl.genes_vs_mt_fraction

```

Scatter plot for **QC inspection of gene complexity vs mitochondrial fraction**.

This plot complements *library size–based QC* by focusing on the relationship
between **detected gene count** and **mitochondrial fraction**, which is
particularly useful for identifying low-quality or degraded samples.


## What it does

Plots each sample as a point:  
- x-axis: mitochondrial fraction (default: pct_counts_mt)
- y-axis: number of detected genes (default: n_genes_detected).  

Optionally:  
- Colors samples by a categorical variable
- Applies QC thresholds on mitochondrial fraction and/or gene count
- Highlights QC-failing samples.   

Automatically reports the number of QC failures in the title when
thresholds are provided. 

This view is especially useful when:  
- Samples with high mt fraction show reduced gene complexity
- Library-size effects have already been inspected separately

## Requirements

adata.obs must contain:  
- the column specified by x (default: "pct_counts_mt")
- the column specified by y (default: "n_genes_detected").  

If groupby is provided, it must also exist in adata.obs.   

## Parameters

#### Axes

**x** (str, default "pct_counts_mt").  
Column in adata.obs used for the x-axis (mitochondrial fraction).   

**y** (str, default "n_genes_detected").   
Column in adata.obs used for the y-axis (gene complexity).   

#### Grouping
**groupby** (str | None, default None). 
Optional categorical column in adata.obs used to color samples
(e.g. batch, subtype).  

#### QC thresholds

Samples are considered QC-pass only if all specified conditions are met.  

**min_mt** (float | None): Minimum allowed mitochondrial fraction (rarely used).  
**max_mt** (float | None). Maximum allowed mitochondrial fraction (commonly used).  
**min_genes** (float | None). Minimum allowed number of detected genes.  
**max_genes** (float | None). Maximum allowed number of detected genes (occasionally useful to remove extreme outliers).  

When any threshold is provided:  
- QC-failing samples are visually highlighted
- The plot title includes the number of failing samples

#### Scaling

**logx** (bool, default False). Apply log scaling to the x-axis (usually unnecessary for percentages).  
**logy** (bool, default True). Apply log scaling to the y-axis (recommended for gene counts).  

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

**High mt fraction + low gene count**.   
Strong indicator of degraded RNA or failed library prep.   

**Low mt fraction + high gene count**.   
Typically high-quality samples.    

**Diagonal trend**.  
Expected in many datasets; extreme deviations are worth inspecting.   

This plot is best interpreted together with:  
- library_size_vs_genes
- mt_fraction_vs_counts
- qc_metrics

## Examples

Basic QC plot
```python
bk.pl.genes_vs_mt_fraction(adata)
```

Apply mitochondrial and gene-count thresholds
```python
bk.pl.genes_vs_mt_fraction(
    adata,
    max_mt=10.0,
    min_genes=12000,
)
```

Color by subtype
```python
bk.pl.genes_vs_mt_fraction(
    adata,
    groupby="Subtype",
    max_mt=8.0,
)
```

Save figure without displaying
```python
bk.pl.genes_vs_mt_fraction(
    adata,
    max_mt=7.5,
    save="genes_vs_mt_qc.png",
    show=False,
)
```

## Related functions
	•	library_size_vs_genes — library size vs gene complexity QC
	•	mt_fraction_vs_counts — mitochondrial fraction vs library size QC
	•	qc_metrics — combined QC overview