# GSEA preranked

```{eval-rst}
.. autofunction:: bullkpy.tl.gsea_preranked

```

Run preranked GSEA using **GSEApy** and return a tidy results table.

This function performs a **GSEA preranked analysis** starting from a differential
expression (or association) table and is designed to integrate cleanly with
the BULLKpy pipeline and `AnnData`.

It wraps `gseapy.prerank` but adds several bulk-friendly and workflow-oriented
improvements.

## What it does
	•	Builds a preranked gene list from a results table (e.g. DE or association results)
	•	Runs GSEApy prerank
	•	Writes results to a structured output folder:

```python
outdir/<comparison>/
```
	•	Returns a tidy DataFrame with NES, p-values, FDR, leading-edge genes, etc.
	•	Optionally returns the raw gseapy result object (pre_res) for downstream plots
	•	Optionally stores results in adata.uns (tables only, h5ad-safe)

## Key features

**1. Flexible gene-set input**   

The gene_sets argument accepts:

**Enrichr library names**   
e.g. "MSigDB_Hallmark_2020"

**GMT files**   
e.g. "path/to/genesets.gmt"

**Dictionaries of gene sets**  
{ "PathwayA": [...], "PathwayB": [...] }

**Convenience aliases**  
	•	"hallmark" / "hallmarks"
	•	"c2" / "curated"

**Lists combining any of the above**   

**2. Comparison-aware output**   
Results are written to:

```python
outdir/<comparison>/
```

	•	A "comparison" column is added to the output table (optional)
	•	Makes it easy to concatenate GSEA results across many contrasts

**3. Duplicate gene handling**   

If duplicated gene names are found in the ranking:
	•	A warning is issued
	•	The entry with the largest absolute score is kept (standard GSEA practice)

**4. Safe AnnData integration**   
	•	Only tables are stored in adata.uns
	•	The raw gseapy object is not stored (not h5ad-serializable)
	•	Provenance and parameters are recorded for reproducibility

## Parameters

**res**   
Differential expression or association results table

**gene_sets**  
Gene sets to test (see supported formats above)

**score_col**  
Ranking statistic (recommended: "t")

**comparison**  
Name of the contrast (used for folder names and table column)

**return_pre_res**  
Whether to return the raw GSEApy object for plotting

**add_comparison_column**   
Add a "comparison" column to the results table

## Returns

**pd.DataFrame**    
Tidy GSEA results table

**(pd.DataFrame, pre_res)** (optional)   
If return_pre_res=True, also returns the raw GSEApy result object

## Output table columns

Typical columns include:
	•	comparison
	•	Term
	•	NES
	•	pval
	•	FDR q-val
	•	Lead_genes
	•	ES
	•	Tag %
	•	Gene %

(Exact columns depend on GSEApy version)

## Examples

**Basic preranked GSEA**   

```python
df_gsea = bk.tl.gsea_preranked(
    adata,
    res=de_res,
    gene_sets="MSigDB_Hallmark_2020",
    comparison="LumB_vs_LumA",
)
```

**Use t-statistic ranking and keep raw GSEA object**   

```python
df_gsea, pre_res = bk.tl.gsea_preranked(
    adata,
    res=de_res,
    score_col="t",
    gene_sets=["hallmark", "c2"],
    comparison="RB1_mut_vs_WT",
    return_pre_res=True,
)
```

**Concatenate multiple GSEA runs**   

```python
dfs = []
for comp, res in all_results.items():
    df = bk.tl.gsea_preranked(
        adata,
        res=res,
        gene_sets="hallmark",
        comparison=comp,
    )
    dfs.append(df)

df_all = pd.concat(dfs, axis=0)
```

**Plot later using the raw GSEA object**     

```python
pre_res.plot(
    term="HALLMARK_INTERFERON_GAMMA_RESPONSE",
    ofname="ifng_enrichment.png",
)
```

## Notes
	•	Recommended ranking statistic: "t"
(best behaved for preranked GSEA)
	•	Input data should be approximately Gaussian
(e.g. log1p CPM, voom, variance-stabilized counts)
	•	The returned pre_res object is intentionally not stored in AnnData

## See also
	•	tl.list_enrichr_libraries
	•	pl.gsea_bubbleplot
	•	tl.gsea_leading_edge_heatmap
	•	tl.leading_edge_overlap_matrix
	•	tl.leading_edge_jaccard_heatmap