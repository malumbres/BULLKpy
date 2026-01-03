# Rank genes groups dotplot

```{eval-rst}
.. autofunction:: bullkpy.pl.rank_genes_groups_dotplot

```

Dotplot for visualizing **top ranked genes per group** from
`adata.uns[key]` (typically created by `rank_genes_groups`). This is a convenience
wrapper that:  

1) selects the top `n_genes` per group from the ranking table,    
2) optionally filters genes by **within-group detection fraction**, and    
3) calls `dotplot()` to render a multi-group gene summary.  
  
It supports two color modes:  
- **expression** (mean expression per group; classic Scanpy-style dotplot).  
- **logfoldchanges** (dot color encodes DE log2FC instead of expression).  

## What it does. 

1) Read ranking results. 
Loads a long-form table via:  
- ank_genes_groups_df_all(adata, key=key, groups=groups, sort_by=sort_by).  

Expected columns include at least:  
- group, gene, scores, log2FC (and often pval, qval).   

2) Rank genes per group.  
- Default ranking uses descending scores.
- If use_abs=True, ranks by abs(scores).   

3) Optional fraction-based gene filtering.  

If min_in_group_fraction and/or max_in_group_fraction are set:   
- It computes, for candidate genes, the fraction of samples within each group with expression > expr_threshold.
- Fraction is computed from: fraction_layer if present in adata.layers, else falls back to layer or adata.X.
- Rows failing the fraction constraints are removed before selecting top genes.

This helps remove:  
- genes expressed in too few samples (dropouts / unstable markers)
- genes expressed in almost all samples (uninformative housekeeping-like).  

4) Select top genes per group.  

Keeps the first n_genes per group after sorting (and filtering).  

If unique=True, the same gene will not be reused across multiple groups
(first group “claims” it).  

5) Render a dotplot (two modes). 

A) values_to_plot="expression" (default): Calls dotplot() on the original adata:
- Dot color = mean expression per group (from layer)
- Dot size  = fraction of samples expressing the gene (from fraction_layer + expr_threshold)
- standard_scale="auto" maps to "var" (Scanpy-like) for expression mode.   

B) values_to_plot="logfoldchanges": Constructs a temporary AnnData where:
- ad_tmp.X[i, gene] = log2FC(group_of_sample_i, gene) (constant within each group).  

Then calls dotplot() so that:  
- Dot color behaves like “mean value per group” → equals log2FC
- Dot size still uses detection fraction from the original counts/frac layer.  

This keeps the size meaningful (detection) while coloring by effect size.  

## Parameters

#### Inputs / selection

**groupby** (str, required).   
adata.obs[groupby] defines the groups (clusters / conditions) used for plotting
and for matching log2FC values.   

**groups** (Sequence[str] | None). 
Subset of groups to include. None uses what’s in the ranking table.    

**key** (str, default "rank_genes_groups").  
Location of ranking results in adata.uns[key].    

**n_genes** (int, default 5).   
Number of genes selected per group.   

**sort_by** ("scores" | "logfoldchanges" | "pvals_adj" | "pvals").  
Passed through to rank_genes_groups_df_all(...) for initial ordering.  

**unique** (bool, default True). 
If True, prevents the same gene from appearing in multiple group panels.   

**use_abs** (bool, default False).  
If True, ranks by abs(scores) rather than signed scores.    

#### Color mode

**values_to_plot** ("expression" | "logfoldchanges", default "expression"). 
Controls what drives dot color:  
- "expression": mean expression per group
- "logfoldchanges": log2FC from ranking results

#### Expression / fraction settings (used by dotplot)

**layer** (str | None, default "log1p_cpm").  
Expression layer used for mean expression when values_to_plot="expression".   

**fraction_layer** (str | None, default "counts").  
Layer used to compute detection fraction (dot size).  
Falls back to layer or adata.X if missing.   

**expr_threshold** (float, default 0.0). 
A sample “expresses” the gene if value > expr_threshold in fraction_layer
(or fallback).  

#### Fraction filters (gene selection)

**min_in_group_fraction** (float | None). 
Keep genes expressed in at least this fraction of samples within the group.  

**max_in_group_fraction** (float | None). 
Keep genes expressed in at most this fraction of samples within the group.  

If filtering removes all candidates, a ValueError is raised advising to relax thresholds.  

#### Scaling

**standard_scale** ("auto" | "var" | "group" | None, default "auto"). 
Passed to dotplot() as standard_scale:  
- "auto" → "var" when plotting expression, else None for logFC mode
- "var"  → z-score per gene across groups
- "group"→ z-score per group across genes
- None   → no scaling

#### Layout (forwarded to dotplot)

**swap_axes** (bool, default True). 
Scanpy-like orientation: genes on y-axis, groups on x-axis.  

**dendrogram_top** (bool, default True). 
Cluster columns (typically groups) and show top dendrogram.  

**dendrogram_rows** (bool, default False).   
Cluster rows (typically genes) and show row dendrogram.   

**row_dendrogram_position** ("right" | "left" | "outer_left").  
Placement for row dendrogram if enabled.   

**row_spacing** (float, default 0.75).  
Compresses/expands vertical spacing between rows (useful for long gene lists).   

**cmap** (str, default "Reds").  
Colormap for dot color (expression or log2FC).   

#### Output

**save** (str | Path | None). 
Save the plot via _savefig.  

**show** (bool, default True).  
Whether to display the figure.   

## Returns

Returns whatever dotplot() returns: (fig, ax). 
(Exact ax type depends on your dotplot implementation.). 


## Notes
- Requires matching group labels: logFC mode assumes the group labels in the
ranking table correspond to adata.obs[groupby] categories.
- fraction_layer should be comparable across samples (counts or normalized counts).
If using raw counts, expr_threshold=0 is typical.
- In logfoldchanges mode, genes without a stored log2FC for a given group are
assigned 0.0 for that group (so they appear neutral).

## Examples

1) Classic Scanpy-style: expression colored, fraction sized
```python
bk.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    n_genes=5,
)
```

2) Color by log2FC instead of expression
```python
bk.pl.rank_genes_groups_dotplot(
    adata,
    groupby="Subtype",
    n_genes=8,
    values_to_plot="logfoldchanges",
    cmap="RdBu_r",
)
```

3) Filter to robust markers (expressed in ≥20% of samples in-group)
```python
bk.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    n_genes=6,
    min_in_group_fraction=0.20,
    fraction_layer="counts",
    expr_threshold=0.0,
)
```

4) Remove “too ubiquitous” genes (expressed in >95% within group)
```python
bk.pl.rank_genes_groups_dotplot(
    adata,
    groupby="leiden",
    n_genes=6,
    min_in_group_fraction=0.10,
    max_in_group_fraction=0.95,
)
```


