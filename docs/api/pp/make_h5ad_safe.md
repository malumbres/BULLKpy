# Make h5ad safe

```{eval-rst}
.. autofunction:: bullkpy.pp.make_h5ad_safe

```

Make an AnnData writable to `.h5ad`, fixing both key names and column dtypes.

Writing fails for two unrelated reasons, and a dataset with real clinical
metadata usually hits both.

## Problem 1 — illegal key names

HDF5 uses `/` as a group separator, so it can never appear in a key. Any
`.obs` or `.var` column whose name contains one makes `adata.write()` fail:

```python
adata.write("results.h5ad", compression="gzip")
```

```text
ValueError: Forward slashes are not allowed in keys in <class 'h5py._hl.group.Group'>
Error raised while writing key 'GP1_Proliferation/DNA_repair' of <class 'h5py._hl.group.Group'> to /obs
```

Score and signature columns are the usual culprits, because names like
`GP1_Proliferation/DNA_repair` or `GP4_MES/ECM` read naturally in an
analysis but are illegal on disk. The same applies to `.uns`, `.obsm`,
`.varm`, `.layers`, `.obsp` and `.varp` keys.

## Problem 2 — un-serialisable dtypes

pandas represents a missing string as `float('nan')`, so a column can end up
holding a mix of Python types. h5py refuses those:

```text
TypeError: Can't implicitly convert non-string objects to strings
Error raised while writing key 'samples.is_ffpe' of <class 'h5py._hl.group.Group'> to /obs
```

anndata copes with more than it may appear: numeric, `bool`, datetime,
`category` and pandas' `str` dtype all write cleanly, and an `object` column
holding only strings is converted to a categorical on write with missing
values preserved. The one case that fails is an **`object` column containing
non-string values** — booleans or numbers alongside the `nan` used for missing
entries. Only those are repaired, so nothing else in your metadata is touched.

Repairs are made by intent rather than by blanket stringification:

| Column | Becomes |
| --- | --- |
| numeric apart from its missing values | numeric, with `NaN` preserved |
| anything else | `category` with string categories, `NaN` preserved |

Missing values stay missing — they never become the literal string `"nan"`.

## The fix

Call it immediately before writing:

```python
bk.pp.make_h5ad_safe(adata)
adata.write("results.h5ad", compression="gzip")
```

Only **names** change — values, dtypes and column order are untouched, and
nothing is dropped.

## What it renames

| Location | Renamed |
| --- | --- |
| `adata.obs` / `adata.var` | column names |
| `adata.uns` | keys, recursively through nested dicts |
| `adata.obsm` / `varm` / `layers` / `obsp` / `varp` | keys |

It also renames anything using anndata's reserved `_index` key, and appends a
numeric suffix if two names would collide after renaming — so
`A/B` alongside an existing `A_B` becomes `A_B_1`, never a silent merge.

## Inspecting the changes

The return value maps each location to its renames:

```python
report = bk.pp.make_h5ad_safe(adata)
report["obs"]
```

```text
{'GP1_Proliferation/DNA_repair': 'GP1_Proliferation_DNA_repair',
 'GP4_MES/ECM': 'GP4_MES_ECM'}
```

Use `copy=True` to get a sanitised copy and leave the original alone:

```python
safe = bk.pp.make_h5ad_safe(adata, copy=True)
```

## Choosing the replacement

`_` is the default. Pass `replacement=""` to delete the character instead, or
any other string:

```python
bk.pp.make_h5ad_safe(adata, replacement="-")
```

## Skipping the dtype pass

Pass `dtypes=False` to rename keys only and leave every column as it is:

```python
bk.pp.make_h5ad_safe(adata, dtypes=False)
```

## If writing still fails

Find the offending column by trial write, which names it directly:

```python
bad, index_err = bk.pp.find_bad_obs_cols_by_write(adata)
```

## See also

- [`find_bad_obs_cols_by_write`](find_bad_obs_cols_by_write.md)
- [`make_obs_h5ad_safe_strict`](make_obs_h5ad_safe_strict.md)
- [`make_var_h5ad_safe_strict`](make_var_h5ad_safe_strict.md)
- [`sanitize_metadata`](sanitize_metadata.md)
