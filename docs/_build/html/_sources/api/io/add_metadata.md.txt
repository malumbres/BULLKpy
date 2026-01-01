# Add metadata

```{eval-rst}
.. autofunction:: bullkpy.io.add_metadata

```

Add sample-level metadata to an existing `AnnData` object.

`bk.io.add_metadata` merges clinical, experimental, or annotation data
into `adata.obs`, aligning rows by sample identifiers.


## Purpose

This function is used to attach **sample annotations** (clinical variables,
phenotypes, batch information, mutation status, etc.) to an `AnnData` object
created from a count matrix.

It is typically run **immediately after** `read_counts`.

---

## Supported metadata formats

- Tab-separated files (`.tsv`)
- Comma-separated files (`.csv`)
- Excel files (`.xls`, `.xlsx`)

The metadata file must contain **one column identifying samples**, which
will be matched to `adata.obs_names`.

---

## Basic usage

```python
import bullkpy as bk

adata = bk.io.read_counts("counts.tsv")

adata = bk.io.add_metadata(
    adata,
    "metadata.tsv",
    index_col="Sample_ID",
)
```

After running this function, all metadata columns will be available in
adata.obs.


## Basic usage

The column specified by index_col must contain sample IDs identical
to adata.obs_names.

Example metadata table:
Sample_ID, Project_ID, Age, Sex  
TCGA-01, LUAD, 64, F   
TCGA-02, LUSC, 71, M   

```python
adata = bk.io.add_metadata(
    adata,
    "clinical.tsv",
    index_col="Sample_ID",
)
```

### Merge strategies

The how argument controls how samples are retained during the merge.

#### Keep all samples (default)

```python
adata = bk.io.add_metadata(
    adata,
    "metadata.tsv",
    index_col="Sample_ID",
    how="left",
)
```
	•	Keeps all samples in adata
	•	Samples without metadata will contain NaN values


#### Keep only samples with metadata

```python
adata = bk.io.add_metadata(
    adata,
    "metadata.tsv",
    index_col="Sample_ID",
    how="inner",
)
```
	•	Drops samples not present in the metadata file
	•	Useful when metadata completeness is required



## Warnings and diagnostics

The function provides informative warnings:
	•	Duplicated sample IDs in metadata
	•	Missing metadata for some samples
	•	Number of samples successfully matched

Example warning:
```python
WARNING: 12 samples in AnnData are missing metadata (showing up to 5)
```

## Output

The function modifies and returns the same AnnData object:
	•	Metadata columns are appended to adata.obs
	•	Existing .obs columns are preserved
	•	Index order remains consistent with sample order

```python
adata.obs.head()
```

## Notes and best practices
	•	Metadata columns with mixed types may be imported as object
	•	For large metadata tables, consider sanitizing dtypes before saving
to .h5ad (e.g. convert categorical variables to category)
	•	This function does not modify .X, .var, or .layers


