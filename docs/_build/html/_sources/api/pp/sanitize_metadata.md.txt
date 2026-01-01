# Sanitize metadata

```{eval-rst}
.. autofunction:: bullkpy.pp.sanitize_metadata

```

Clean and standardize sample metadata before adding it to `AnnData`.

## When should I use this?

You should run sanitize_metadata() before assigning a metadata table to
adata.obs, especially when:  
	•	Metadata comes from TCGA / GEO / Excel
	•	Many columns are object dtype
	•	.h5ad writing fails
	•	Correlation or association functions ignore metadata columns

## Typical workflow. 

```python
meta = pd.read_csv("metadata.tsv", sep="\t")
meta = bk.pp.sanitize_metadata(meta, index_col="Sample_ID")
adata.obs = meta
```

##Key behaviors  
	•	Converts numeric-looking strings to numeric
	•	Converts low-cardinality strings to categories
	•	Converts date-like columns to datetime
	•	Converts remaining text to pandas string
	•	Optionally drops ID-like columns


## Parameters   

**df**   
    Metadata table with samples in rows and annotations in columns.

**index_col**   
    Optional column name to set as the DataFrame index (e.g. sample ID).

**numeric_min_frac**   
    Minimum fraction of non-missing values that must successfully parse as
    numeric for a column to be converted to numeric (default: 0.9).  

**category_max_unique**  
    Maximum number of unique values for a string column to be converted to
    a categorical dtype (default: 50).

**category_max_frac**   
    Maximum fraction of unique values (relative to number of rows) for a column
    to be converted to categorical (default: 0.2).

**datetime_min_frac**  
    Minimum fraction of values that must successfully parse as datetimes for
    a column to be converted to datetime64 (default: 0.9).

**drop_high_cardinality_strings**   
    If True, drop string columns that appear ID-like (very high cardinality).  

**high_cardinality_frac**   
    Fraction of unique values above which a column is considered high-cardinality
    and ID-like (default: 0.5).  

**verbose**    
    If True, print a per-column report describing how each column was treated.

## Returns

pandas.DataFrame  
    A sanitized copy of the input DataFrame with improved dtypes.

## Notes  

- Common missing-value tokens (e.g. `"NA"`, `"None"`, `"."`, `"-"`) are normalized
  to `NaN` before type inference.
- Numeric inference uses `pandas.to_numeric(errors="coerce")`.
- Datetime inference uses `pandas.to_datetime(errors="coerce")`.
- Categorical conversion is preferred for low-cardinality strings.
- Remaining text columns are converted to pandas `string` dtype (not Python
  `object`) to improve HDF5 compatibility.

## Examples   

Sanitize a TCGA metadata table before adding it to AnnData:

```python
meta = pd.read_csv("TCGA_metadata.tsv", sep="\t")
meta = sanitize_metadata(meta, index_col="Sample_ID")
adata.obs = meta
```


Drop ID-like columns automatically:  

```python

meta = sanitize_metadata(
           meta,
           drop_high_cardinality_strings=True,
)
```
