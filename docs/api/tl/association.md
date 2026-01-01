# Association

```{eval-rst}
.. autofunction:: bullkpy.tl.association

```

Unified dispatcher for association analyses between genes and annotations.

This function provides a **single entry point** to test associations between:

- a **gene** and a **categorical annotation**
- a **numeric obs column** and a **categorical annotation**
- two **categorical annotations**

It automatically detects the types of `x` and `y` and routes the analysis
to the appropriate lower-level function.

## What it does

Depending on the nature of x and y, association dispatches to:

| x | y | Dispatched function |
| ---------- | ---------- | -------------------- |
| gene | categorical obs | gene_categorical_association |
| numeric obs | categorical obs | obs_categorical_association |
| categorical obs | categorical obs | categorical_association |

The goal is to let users write:

```python
bk.tl.association(adata, x="TP53", y="Subtype")
```

instead of worrying about which specific association function to call.

## Parameters

**adata**   
AnnData object containing expression data and annotations

**x**   
Name of a gene (adata.var_names) or   
name of an obs column (adata.obs)   

**y**   
Name of a gene (adata.var_names) or    
name of an obs column (adata.obs)   

**layer**    
Expression layer to use when x or y refers to a gene    
(default: "log1p_cpm")    

**method**   
Statistical method to use (passed through to the underlying function).     
Common values:    
	•	"kruskal", "anova" (multi-group)
	•	"mwu", "ttest" (two-group)
	•	"auto" (use defaults of the dispatched function)

## Dispatch logic (conceptual)

```python
gene ↔ categorical obs
    → gene_categorical_association (single gene)

numeric obs ↔ categorical obs
    → obs_categorical_association (single variable)

categorical obs ↔ categorical obs
    → categorical_association
```

Numeric–numeric and gene–numeric associations are intentionally excluded
from this dispatcher and should be handled by correlation utilities instead.

⸻

## Returned value

The return type depends on the dispatched function:


| Case | Return type |
| ---------- | ---------- |
|gene ↔ categorical | pd.DataFrame |
| numeric obs ↔ categorical | pd.DataFrame |
| categorical ↔ categorical | dict |

See the documentation of the underlying function for exact structure.

⸻

## Examples

**Gene vs categorical annotation**   

```python
bk.tl.association(
    adata,
    x="TP53",
    y="Subtype",
)
```

Equivalent to:

```python
bk.tl.gene_categorical_association(
    adata,
    genes=["TP53"],
    groupby="Subtype",
)
```

**Numeric obs vs categorical annotation**  

```python
bk.tl.association(
    adata,
    x="age",
    y="Subtype",
)
```

Equivalent to:

```python
bk.tl.obs_categorical_association(
    adata,
    obs_keys=["age"],
    groupby="Subtype",
)
```

**Categorical vs categorical**    

```python
bk.tl.association(
    adata,
    x="Batch",
    y="Subtype",
)
```

Equivalent to:

```python
bk.tl.categorical_association(
    adata,
    key1="Batch",
    key2="Subtype",
)
```

## Error handling

The function raises informative errors when:
	•	x or y cannot be resolved to gene or obs
	•	unsupported combinations are requested (e.g. gene ↔ numeric)
	•	required columns are missing

## Design notes
	•	This is a dispatcher, not a statistical method itself
	•	It favors explicitness over magic:
unsupported cases are rejected
	•	Intended for interactive exploration and pipelines
	•	Keeps the public API minimal while preserving flexibility

## See also
	•	tl.gene_categorical_association
	•	tl.obs_categorical_association
	•	tl.categorical_association
	•	tl.rank_genes_categorical
	•	pl.violin


