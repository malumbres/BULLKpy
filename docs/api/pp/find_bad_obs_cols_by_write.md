# Find bad Obs columns

```{eval-rst}
.. autofunction:: bullkpy.pp.find_bad_obs_cols_by_write

```

Detect problematic columns in adata.obs that prevent saving an AnnData object to disk.

## Description

find_bad_obs_cols_by_write is a diagnostic utility that identifies columns in
adata.obs which cause adata.write() to fail (typically when writing .h5ad files).  

This problem commonly arises when metadata columns contain:  
	•	mixed Python object types
	•	lists, dictionaries, or sets
	•	byte strings
	•	partially parsed datetimes
	•	other non-HDF5-serializable values

Unlike simple dtype checks, this function actually attempts to write each column
into a minimal AnnData object, closely mimicking the real adata.write() behavior.



## Why this function is needed.  

Standard checks such as:   

```python
adata.obs.dtypes
```

or even select_dtypes() are often insufficient, because:  
	•	a column may look like object or string
	•	only some rows contain problematic values
	•	the first row is often clean, while later rows break serialization

This function explicitly tests multiple informative rows per column to reliably
expose hidden issues.  

## Parameters  

**adata**   
AnnData object whose `.obs` table will be tested.

**n_rows**  
  - `int`, optional
  - Maximum number of rows per column to test (default: 200).  
    Rows are selected to maximize the chance of exposing problematic values.

**include_index_test**   
  - `bool`, optional
  - Whether to also test whether `adata.obs_names` (index) can be written safely.

## Returns    

**bad**:   
  - `list[tuple[str, str]]`
  - List of `(column_name, error_message)` pairs for columns that failed to write.

**index_error**   
  - `str | None`
  - Error message if the obs index itself fails to write, otherwise `None`.

## Typical usage   

```python
bad, index_err = bk.io.find_bad_obs_cols_by_write(adata)
bad
```
 Example output:
```python
[
    ("clinical_notes", "TypeError: Can't implicitly convert non-string objects to strings"),
    ("mutation_list", "Object dtype dtype('O') has no native HDF5 equivalent"),
]
```

## Inspecting problematic columns

Once an offending column is identified, inspect its contents:  
```python
col = bad[0][0]
adata.obs[col].map(type).value_counts()
```

## Recommended fixes

Depending on the column: 

**Issue type -->  Suggested fix**.  
Lists / dicts --> Convert to string or explode.   
Mixed numeric + strings --> Use pd.to_numeric(errors="coerce").   
High-cardinality strings --> Drop column or convert to categorical.  
Free-text fields --> Keep as pandas "string" dtype.  
Metadata table --> Use sanitize_metadata() before merging.  

## Related functions
`sanitize_metadata`   
Automatically cleans metadata tables before adding them to adata.obs.

`adata.write()` (AnnData).  
Final serialization step that this function helps debug.

## Notes
	•	This function does not modify adata
	•	Only .obs is tested (not .var, .uns, .obsm, etc.)
	•	If no bad columns are found but adata.write() still fails, the issue is likely in:
	•	adata.uns
	•	adata.var
	•	adata.obsm / adata.varm
	•	adata.layers

If needed, similar diagnostic functions can be applied to those components.