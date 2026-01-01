# Categorical association

```{eval-rst}
.. autofunction:: bullkpy.tl.categorical_association

```

Association analysis between two categorical annotations.

This function quantifies the relationship between **two categorical columns**
in `adata.obs` using contingency tables and multiple association metrics.
It is useful for comparing clusterings, annotations, batches, or any pair
of categorical labels.

## What it does

Given two categorical variables (key1, key2), the function:
	1.	Builds a contingency table
	2.	Computes one or more association metrics
	3.	Returns all results in a single dictionary

Supported metrics include:
	•	Chi-squared test
	•	Cramér’s V
	•	Adjusted Rand Index (ARI)
	•	Normalized Mutual Information (NMI)

## When to use

Use categorical_association when you want to:
	•	Compare cluster labels vs known annotations
	•	Quantify agreement between two clustering solutions
	•	Assess batch effects
	•	Explore relationships between categorical metadata fields

Examples:
	•	Leiden clusters vs cell type
	•	Batch vs condition
	•	Manual annotation vs automated labels

## Parameters

**adata**   
AnnData object containing the annotations in .obs

**key1**    
First categorical column in adata.obs

**key2**   
Second categorical column in adata.obs

**metrics**   
Iterable of metrics to compute. Supported values:
	•	"chi2" – Chi-squared test of independence
	•	"cramers_v" – Effect size for categorical association
	•	"ari" – Adjusted Rand Index (requires scikit-learn)
	•	"nmi" – Normalized Mutual Information (requires scikit-learn)

**dropna**   
If True, rows with missing values in either column are dropped before analysis

## Returned value

Returns a dict with the following entries:

**table**    

A pandas DataFrame representing the contingency table:

```python
key2 categories → (columns)
key1 categories ↓ (rows)
```

**Optional metric entries**   

Depending on metrics, the dictionary may also include:

**chi2**   

```python
{
  "statistic": float,
  "pval": float,
  "dof": int
}
```

**cramers_v**.  

```python
float
```

**ari**. 

```python
float
```

**nmi**.   

```python
float
```

## Metrics explained

**Chi-squared test**.   
	•	Tests independence between categories
	•	Sensitive to sample size
	•	Returns statistic and p-value

**Cramér’s V**.  
	•	Effect size for categorical association
	•	Range: 0–1
	•	Interpretable regardless of table size

**Adjusted Rand Index (ARI)**.   
	•	Measures clustering similarity
	•	Adjusted for chance
	•	Range: –1 to 1
	•	Requires scikit-learn

**Normalized Mutual Information (NMI)**.   
	•	Information-theoretic similarity
	•	Range: 0–1
	•	Requires scikit-learn

## Examples

**Compare clustering vs annotation**    

```python
out = bk.tl.categorical_association(
    adata,
    key1="leiden_1.0",
    key2="CellType",
)

out["cramers_v"]
```

**Full metric set**   

```python
out = bk.tl.categorical_association(
    adata,
    key1="Batch",
    key2="Condition",
    metrics=("chi2", "cramers_v", "ari", "nmi"),
)
```

**Inspect contingency table**   

```python
out["table"]
```

## Notes
- ARI and NMI require scikit-learn
- If scikit-learn is unavailable:
	•	ARI/NMI are skipped with a warning
- Chi-squared uses no Yates correction
- All categories are cast to strings before comparison

## See also
	•	tl.leiden_resolution_scan
	•	tl.rank_genes_categorical
	•	tl.gene_categorical_association
	•	tl.obs_categorical_association
	•	scanpy.tl.leiden



