# Obs-Category association

```{eval-rst}
.. autofunction:: bullkpy.tl.obs_categorical_association

```

Association between **numeric sample-level variables** and a **categorical** variable.

This function tests whether numeric observation columns (QC metrics, clinical
variables, signature scores) **differ across the categories** of another
observation. It is the obs-level analogue of
[`gene_categorical_association`](gene_categorical_association.md).

## What it does

For each numeric column in `adata.obs`, the function:

1. Splits samples by the categories in `adata.obs[groupby]`.
2. Tests whether values differ across groups using Kruskal–Wallis (default,
   non-parametric) or one-way ANOVA (parametric).
3. Computes an effect size.
4. Applies multiple-testing correction.
5. Returns a tidy table with one row per obs variable.

## When to use

Use `obs_categorical_association` when you want to:

- test QC metrics across conditions (e.g. library size vs batch),
- assess clinical variables across subtypes,
- screen many numeric obs columns at once,
- perform a global multi-group test (not pairwise).

## Parameters

**adata**
AnnData object containing observations in `.obs`.

**groupby**
Categorical column in `adata.obs` defining the groups.

**obs_keys**
Numeric obs columns to test. If `None` (default), every numeric column in
`adata.obs` is used, excluding `groupby` itself.

**test**
Global test to apply:

- `"kruskal"` – Kruskal–Wallis H-test (default, non-parametric)
- `"anova"` – one-way ANOVA
- `"auto"` – resolves to `"kruskal"`

**method**
Alias for `test`, accepted for convenience. If both are supplied, `method` wins.

**effect_size**

- `"epsilon2"` – Kruskal–Wallis effect size (recommended)
- `"eta2"` – ANOVA effect size
- `"auto"` – picks the one matching `test`
- `"none"` – skip effect size

**min_group_size**
Minimum number of samples a group needs to be included in the test.

**adjust**
Multiple-testing correction: `"bh"` (Benjamini–Hochberg, default) or `"none"`.

## Output

Returns a tidy `DataFrame` with one row per obs variable:

| Column | Description |
| --- | --- |
| `obs` | Obs column name |
| `groupby` | Name of the grouping variable |
| `n_groups` | Number of groups |
| `test` | Test actually used |
| `statistic` | Test statistic (H or F) |
| `pval` | Raw p-value |
| `qval` | BH-adjusted p-value |
| `effect` | Effect size (ε² or η²) |
| `mean_<group>` | Mean value in each group |

Results are sorted by `qval`, then `pval`.

## Effect sizes

- **ε² (epsilon-squared)** – Kruskal–Wallis; range 0–1; proportion of variance
  explained by group membership. Robust and non-parametric (recommended).
- **η² (eta-squared)** – the parametric ANOVA analogue.

## Examples

Test all numeric obs columns across groups:

```python
res = bk.tl.obs_categorical_association(adata, groupby="Batch")
res.head()
```

Test selected QC metrics only:

```python
res = bk.tl.obs_categorical_association(
    adata,
    groupby="Subtype",
    obs_keys=["total_counts", "pct_counts_mt", "pct_counts_ribo"],
)
```

Use ANOVA instead of Kruskal–Wallis:

```python
res = bk.tl.obs_categorical_association(
    adata,
    groupby="Project_ID",
    obs_keys=["purity", "NE25_score"],
    method="anova",
    effect_size="eta2",
)
```

## Notes

- Groups with fewer than `min_group_size` samples are ignored.
- At least two valid groups are required per variable; otherwise the row is
  returned with `NaN` statistics.
- This is a **global** association test. For pairwise comparisons use
  [`pairwise_posthoc`](pairwise_posthoc.md).
- Works naturally with QC metrics computed by `bk.pp.qc_metrics`.

## See also

- [`gene_categorical_association`](gene_categorical_association.md)
- [`categorical_association`](categorical_association.md)
- [`pairwise_posthoc`](pairwise_posthoc.md)
- [`posthoc_per_gene`](posthoc_per_gene.md)
