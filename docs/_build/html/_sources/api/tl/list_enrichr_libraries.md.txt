# List enrichr libraries

```{eval-rst}
.. autofunction:: bullkpy.tl.list_enrichr_libraries

```

List available Enrichr gene-set libraries.

This function queries **Enrichr** (via `gseapy`) and returns the exact names of all
gene-set libraries currently available for enrichment analysis.

It is mainly intended as a **discovery helper**, so users know which library names
can be passed to Enrichr-based workflows (e.g. GSEApy preranked analysis).

## What it does

	•	Connects to the Enrichr API through gseapy
	•	Fetches the list of registered gene-set libraries
	•	Returns them as a Python list of strings

These names must be used exactly as returned when running Enrichr analyses.

## Returns

**list[str]**  
List of Enrichr library identifiers.

Examples include:
	•	"MSigDB_Hallmark_2020"
	•	"KEGG_2021_Human"
	•	"GO_Biological_Process_2023"
	•	"Reactome_2022"

## Requirements.  
	•	Internet access
	•	gseapy installed
	•	Enrichr service available

If Enrichr is unreachable, a RuntimeError is raised.

## Examples

List all available Enrichr libraries

```python
libs = bk.tl.list_enrichr_libraries()
print(libs[:10])
```

Check if a library exists before running GSEA

```python
if "MSigDB_Hallmark_2020" in bk.tl.list_enrichr_libraries():
    print("Hallmark gene sets available")
```

Use with Enrichr-based GSEA

```python
df_gsea = bk.tl.gsea_preranked(
    adata,
    res=de_res,
    gene_sets="MSigDB_Hallmark_2020",
)
```

## Notes
	•	Library availability may change over time as Enrichr updates its database.
	•	For reproducible analyses, it is recommended to record the library name and date.
	•	This function does not download gene sets; it only lists metadata.

## See also
	•	tl.gsea_preranked
	•	tl.gsea_leading_edge_heatmap
	•	tl.leading_edge_overlap_matrix
	•	tl.leading_edge_jaccard_heatmap
