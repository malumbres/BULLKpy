# Batch correction with Combat

```{eval-rst}
.. autofunction:: bullkpy.pp.batch_correct_combat

```

ComBat batch correction (Johnson et al.) for bulk expression matrices.

`batch_correct_combat` removes unwanted batch-associated variation while optionally **preserving biological covariates** (e.g., subtype, Project_ID). It is intended for **approximately Gaussian** expression values, so use **log-normalized** layers (e.g. `log1p_cpm`) rather than raw counts.

## What it does
	•	Takes an expression matrix (default: adata.layers["log1p_cpm"]) and a batch label in adata.obs[batch_key].
	•	Fits a linear model with:
	•	intercept
	•	optional covariates (biological variables to preserve)
	•	batch indicators
	•	Applies empirical Bayes shrinkage to batch effect parameters and returns a corrected matrix.
	•	By default, stores the corrected values as a new layer: adata.layers[key_added] (default: "combat").

## Inputs

### Required
	•	adata: AnnData (samples in .obs, genes in .var_names)
	•	batch_key: column in adata.obs with batch labels (categorical recommended)

### Optional
	•	layer: which matrix to correct ("log1p_cpm" recommended). If None, uses adata.X.
	•	covariates: list of adata.obs columns to include in the design matrix (preserved effects).
	•	Numeric covariates are included as-is
	•	Categorical covariates are one-hot encoded
	•	key_added: layer name for corrected matrix if overwrite=False
	•	overwrite:
	•	False (default): write to adata.layers[key_added]
	•	True: overwrite the selected layer (or .X if layer=None)
	•	inplace:
	•	True (default): write into adata, return None
	•	False: return corrected matrix as a NumPy array (n_samples, n_genes)


## Returns
	•	If inplace=True: returns None and stores corrected matrix in .layers (or overwrites)
	•	If inplace=False: returns a NumPy array with corrected expression values (samples × genes)

## Notes and recommendations
	•	Do NOT run ComBat on raw counts. Use log-normalized expression (e.g. CPM/TPM + log1p).
	•	If your batch variable has < 2 batches, the function will skip correction (returns input matrix if inplace=False).
	•	Use covariates to avoid removing true biological effects correlated with batch.

Typical covariates:
	•	tumor subtype, tissue type, sex, Project_ID (if biological)
	•	avoid adding covariates that are actually batch proxies unless you explicitly want to preserve them

## Examples

#### 1) Correct by sequencing center / batch and store to a new layer

```python
import bullkpy as bk

bk.pp.batch_correct_combat(
    adata,
    batch_key="Center",       # e.g., sequencing center
    layer="log1p_cpm",
    key_added="combat",
)

# downstream: use corrected layer
bk.tl.pca(adata, layer="combat", use_highly_variable=True)
bk.pl.pca_scatter(adata, color="Center")
```

#### 2) Preserve a biological covariate (e.g. subtype) while correcting batch
```python
bk.pp.batch_correct_combat(
    adata,
    batch_key="Center",
    layer="log1p_cpm",
    covariates=["Project_ID"],   # preserve biological signal
    key_added="combat_cov",
)
```

#### 3) Overwrite the original layer
```python
bk.pp.batch_correct_combat(
    adata,
    batch_key="Center",
    layer="log1p_cpm",
    overwrite=True,      # writes back to adata.layers["log1p_cpm"]
)
```

#### 4) Get the corrected matrix without modifying adata
```python
Xc = bk.pp.batch_correct_combat(
    adata,
    batch_key="Center",
    layer="log1p_cpm",
    inplace=False,
)
print(Xc.shape)  # (n_samples, n_genes)
```

## See also
	•	pp.qc_metrics, pp.filter_samples, pp.highly_variable_genes
	•	tl.pca, pl.pca_scatter, pl.umap