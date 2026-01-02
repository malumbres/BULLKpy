# Rank genes groups

```{eval-rst}
.. autofunction:: bullkpy.tl.rank_genes_groups

```

Scanpy-like per-group gene ranking for **bulk-friendly** AnnData objects.  

`rank_genes_groups` tests, for each group in `adata.obs[groupby]`, which genes are
most different between that group and a reference (either `"rest"` or a specific
group). Results are stored in `adata.uns[key_added]` in a Scanpy-compatible
“dict-of-arrays per group” structure. 

## What it does

For each selected group g:  
1. Split samples into:  
- group: samples where adata.obs[groupby] == g
- reference:
	- "rest": all other samples, or
	- a specific group name.  

2. Compute per-gene:  
- mean expression in group and reference
- a test statistic (scores)
- p-values and FDR-adjusted p-values
- a fold-change-like quantity (logfoldchanges).   

3. Rank genes by scores (or abs(scores) if use_abs=True) and store the top
n_genes.

## Statistics

**method="t-test" (default)**   
- Performs a Welch t-test per gene (unequal variances).
- score = t statistic
- p-value = two-sided t-test p-value using Welch–Satterthwaite df.  

**method="wilcoxon"**   
- Performs a Mann–Whitney U test per gene.
- score = z-like approximation derived from U (no tie correction)
- p-value = two-sided MWU p-value.  

Notes:  
- wilcoxon is usually slower because it loops gene-by-gene.
- For bulk-sized datasets (hundreds of samples, thousands of genes), this is
typically still OK, but for very large matrices prefer "t-test".

## Parameters

**adata**  
AnnData with expression data in adata.X or adata.layers[layer].

**groupby**  
Categorical column in adata.obs defining groups (e.g. "condition").

**groups**  
Subset of group names to analyze. If None, uses all categories in
adata.obs[groupby].

**reference**  
Reference for each group comparison:
- "rest" (default): all samples not in the target group
- "SomeGroup": compare each group vs that explicit group

**layer**   
Expression matrix to use:
- if not None, uses adata.layers[layer]
- else uses adata.X

**method**  
"t-test" (Welch) or "wilcoxon" (MWU).

**corr_method**  
Multiple testing correction. Currently only:
- "benjamini-hochberg" (BH FDR)

**use_abs**   
If True, ranks by abs(scores) (largest magnitude changes first).  
If False, ranks by scores descending (positive effects first).  

**n_genes**   
Number of top genes stored per group.

**key_added**   
Key in adata.uns where results are stored.

## Output

Results are stored in:  
```python
adata.uns[key_added]
```

with:  
- params: dict of the parameters used
- names[group]: array of top gene names
- scores[group]: array of test statistics (t or z-like)
- logfoldchanges[group]: array of log2 fold-change-like values
- pvals[group]: array of raw p-values
- pvals_adj[group]: array of BH-adjusted p-values
- mean_group[group]: mean expression in the target group (bulk-friendly extra)
- mean_ref[group]: mean expression in the reference (bulk-friendly extra)

This layout is intentionally similar to Scanpy so downstream plotting/helpers can
consume it easily.

### Interpretation of logfoldchanges

logfoldchanges is computed as:

```python
(mean_group - mean_ref) / ln(2)
```

This is exact log2 fold change only if the input layer is on a natural-log
scale. If you use layer="log1p_cpm" (default), it behaves like a log-scale
mean difference expressed in log2 units, which is often still useful for
ranking/visualization.

If you need strict log2FC from counts, compute it from raw means on a linear
scale (or use your DE routine that models counts).

### Raises
- ImportError if scipy is not available (scipy.stats is required).
- KeyError if groupby, requested groups, or reference are not found.
- ValueError if an unsupported method or corr_method is provided.

Groups or references with 0 samples are skipped with a warning.

## Examples

Rank genes for all groups vs rest
```python
bk.tl.rank_genes_groups(
    adata,
    groupby="condition",
    reference="rest",
    layer="log1p_cpm",
    method="t-test",
    n_genes=50,
)
```

Compare each group vs an explicit reference group
```python
bk.tl.rank_genes_groups(
    adata,
    groupby="treatment",
    groups=["drugA", "drugB"],
    reference="vehicle",
    method="wilcoxon",
    n_genes=100,
)
```

Access results
```python
rg = adata.uns["rank_genes_groups"]
rg["names"]["drugA"][:10]
rg["pvals_adj"]["drugA"][:10]
```

## See also
	•	tl.de – two-group DE with Welch t-test or NB GLM (counts-aware)
	•	pl.volcano, pl.rankplot – visualization of DE-like tables
	•	tl.rank_genes_categorical – single contrast categorical association (bulk-friendly)