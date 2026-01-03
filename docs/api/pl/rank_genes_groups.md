# Rank genes groups

```{eval-rst}
.. autofunction:: bullkpy.pl.rank_genes_groups

```

Scanpy-like **quick summary plot** for differential expression results produced by
`sc.tl.rank_genes_groups` (or compatible outputs stored in `adata.uns`).

This function creates a **compact table-style figure** showing the top-ranked genes
per group, suitable for quick inspection, slides, or reports.


## What it does. 

1. Loads rank_genes_groups results. 
Uses rank_genes_groups_df_all(...) to convert Scanpy-style results stored in `adata.uns[key]`.  

Supports standard Scanpy fields:
- gene
- group
- scores
- logfoldchanges (shown as log2FC)
- pvals
- pvals_adj (shown as qval).    

2. Selects top genes per group
- Results are grouped by group
- Keeps the top n_genes rows per group, after sorting.   

3. Formats a readable table.  
Each cell contains: GENE_NAME | +log2FC | q=FDR.  
Example:  
```python 
EPCAM  +2.34  q=1.2e-05
```. 

4. Renders a matplotlib table
- Columns = groups
- Rows = rank (1, 2, 3, …)
- No axes; designed for presentation

## Parameters

#### Core

**adata** (AnnData): Annotated data matrix containing DE results in adata.uns.  

**key** (str, default "rank_genes_groups"): Key in adata.uns where results are stored.  

**groups** (Sequence[str] | None): Subset of groups to include.  
If None: include all groups found in the results.   

**n_genes** (int, default 10): Number of top genes to display per group.   

**sort_by** ("scores" | "logfoldchanges" | "pvals_adj" | "pvals", default "scores"): Metric used to rank genes within each group.  

#### Display

**figsize** ((width, height) | None): Auto-scaled if None:   
```python 
width  ~ number_of_groups
height ~ n_genes
```

 **title** (str | None): Plot title.   
Defaults to: "Top {n_genes} ranked genes per group".   

#### Output

**save** (str | Path | None): If provided, saves the figure via _savefig.   

**show** (bool, default True):  Whether to display the plot.

## Returns
- **fig**: matplotlib.figure.Figure
- **ax**: matplotlib.axes.Axes
(axis is turned off; contains only the table)

## Requirements
- matplotlib
- Results must already exist in adata.uns[key]
- Raises ValueError if no results are found

## Interpretation guide
Each column → one group (cluster / condition). 
Each row → rank within that group.  
Displayed metrics   
- +log2FC: direction and magnitude of change
- q=: Benjamini–Hochberg FDR.   

This view is not intended for statistical inference, but for fast comparison and reporting.

## Best practices

Use this after:  
```python
sc.tl.rank_genes_groups(adata, groupby="leiden")
```
For publication-quality plots, combine with:
- dotplot
- heatmap_de
- gene_plot

## Examples

1) Default view (all groups)
```python
bk.pl.rank_genes_groups(adata)
```

2) Specific groups only
```python
bk.pl.rank_genes_groups(
    adata,
    groups=["Basal", "Luminal"],
    n_genes=15,
)
```

3) Sort by effect size
```python
bk.pl.rank_genes_groups(
    adata,
    sort_by="logfoldchanges",
)
```

4) Save for a report
```python
bk.pl.rank_genes_groups(
    adata,
    save="ranked_genes_table.pdf",
    show=False,
)
```

## When to use this vs other plots

| Goal | Recommended plot |
| --------------- | ------------------- |
| Quick DE summary | rank_genes_groups |
| Expression patterns | heatmap_de |
| Marker comparison | dotplot |
| Single gene validation | gene_plot |




