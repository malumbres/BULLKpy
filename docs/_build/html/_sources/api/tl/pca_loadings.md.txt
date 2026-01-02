# PCA loadings

```{eval-rst}
.. autofunction:: bullkpy.tl.pca_loadings

```

Rank genes by **PCA loadings** for one or more principal components (PCs), and optionally export the results for downstream enrichment analysis (e.g., GSEA/Enrichr).

This is a bulk-friendly utility that assumes PCA has already been computed with `bk.tl.pca`, producing:

- `adata.varm["PCs"]` (or your chosen `loadings_key`) with shape `(n_genes, n_comps)`
- (optionally) `adata.uns["pca"]` (or your chosen `key`) containing PCA metadata

The function returns per-PC gene rankings (positive / negative or absolute), stores them in `adata.uns`, and can export them as TSV/CSV and/or GMT.


## What it does

For each requested PC:
- Extract the loading vector (one value per gene).
- Optionally drop NaNs and/or filter by min_abs_loading.
- Rank genes by:
	- absolute loading (use_abs=True) → one list per PC, or
	- signed loading (use_abs=False) → top positive list, plus optional top negative list
- Store results in adata.uns[store_key].
- Optionally export:
	- a long table (PC, sign, rank, gene, loading)
	- a wide table (one column per PC/sign, rows are genes)
	- a GMT file (gene sets per PC/sign).   

## Parameters

#### Inputs / keys

**adata**   
AnnData with PCA loadings stored in adata.varm[loadings_key].

**key**   
PCA metadata key in adata.uns (not strictly required, mostly used for provenance strings).  
Default: "pca".  

**loadings_key**    
Where PCA loadings are stored in adata.varm.  
Default: "PCs".  

#### PC selection and ranking behavior

**pcs**   
PCs to process using 1-based indexing (e.g. [1, 2, 3]).  
If None, processes all PCs available in the loadings matrix.  

**n_top**   
Number of genes to return per PC.
- If use_abs=True: returns n_top genes per PC.
- If use_abs=False and include_negative=True: returns up to n_top positive + up to n_top negative per PC.

**use_abs**   
If True, rank by abs(loading) and return a single list per PC.  
Output keys look like: PC1_abs, PC2_abs, ….  

**include_negative**    
Only used when use_abs=False.  
If True, returns both positive and negative gene lists per PC:  
- PC1_pos = highest positive loadings
- PC1_neg = most negative loadings

#### Filtering / formatting

**dropna**   
If True, drop genes with NaN loadings before ranking.

**gene_col**   
Column name used for the gene identifier in output tables.  
Default: "gene".  

**min_abs_loading**    
Optional threshold: keep only genes with abs(loading) >= min_abs_loading.  
Useful to avoid exporting near-zero loadings.  

#### Storage and export

**store_key**  
Where results are stored in adata.uns.  
Default: "pca_loadings".  

**export**   
If provided, writes:  
- 1. a long-format table to this path, and
- 2. a wide-format table next to it named:
"<stem>.wide<suffix>"

**export_sep**   
Delimiter for export. Either "\t" or ",".  
Default: "\t".  

**export_gmt**   
If provided, writes a GMT file where each PC/sign is a gene set:
- If use_abs=False: PC1_pos, PC1_neg, …
- If use_abs=True: PC1_abs, PC2_abs, …

## Returns

A dictionary mapping set keys → DataFrames.  

Keys are one of:
- Signed mode (use_abs=False):
	- "PC1_pos", "PC1_neg", "PC2_pos", …
- Absolute mode (use_abs=True):
	- "PC1_abs", "PC2_abs", …. 

Each DataFrame contains:  
- pc (e.g., "PC1")
- sign ("pos", "neg", or "abs")
- rank (1..N)
- <gene_col> (gene name)
- loading (raw loading value)

#### Stored in adata.uns

Results are stored under:
```python
adata.uns[store_key]
```
Structure:  
- params: run parameters (pcs, n_top, use_abs, include_negative, etc.)
- results: dict of per-set tables (same as return value)
- table: concatenated long-format table across all PCs/signs. 

## Exports

**Long table (export)**.  

Columns: pc, sign, rank, gene, loading.   

This is convenient for:
- auditing loadings
- making plots
- joining with annotations

**Wide table (<stem>.wide<suffix>)**.   

Each column is a PC/sign set and rows are ranked gene names:
- PC1_pos, PC1_neg, …

Useful for quick inspection or manual copy/paste.

**GMT (export_gmt)**.   

Each line is a gene set:
```python
set_name<TAB>description<TAB>gene1<TAB>gene2...
```
Description is set to "{key}:{loadings_key}" (e.g. "pca:PCs").

## Examples

Get top positive/negative genes for first 3 PCs
```python
out = bk.tl.pca_loadings(
    adata,
    pcs=[1, 2, 3],
    n_top=50,
    use_abs=False,
    include_negative=True,
)
```

Rank by absolute loadings (single list per PC)
```python
out = bk.tl.pca_loadings(
    adata,
    pcs=[1, 2, 3],
    use_abs=True,
    n_top=100,
)
```

Filter weak loadings and export TSV + GMT
```python
bk.tl.pca_loadings(
    adata,
    pcs=[1, 2, 3],
    n_top=100,
    include_negative=True,
    min_abs_loading=0.02,
    export="results/pca_loadings.tsv",
    export_gmt="results/pca_loadings.gmt",
)
```

Access the stored long table
```python
adata.uns["pca_loadings"]["table"].head()
```

## Notes / tips
- PC indexing is 1-based (pcs=[1,2] means PC1 and PC2).
- If you ran PCA using only HVGs, non-used genes may have NaN loadings in adata.varm["PCs"].
	- Keep dropna=True (default) to avoid empty/garbage ranks.
- For enrichment, signed sets (PC1_pos vs PC1_neg) are often more interpretable than absolute sets.
