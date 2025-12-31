# Violin plots

Use `bk.pl.violin` to plot QC/clinical variables from `adata.obs` or gene expression
(from `adata.var_names`) across groups.

```{figure} ../../_static/violin_genes_example.png
:alt: Violin plot example
:width: 500px
:align: center
```

Example violin plot showing gene expression across tumor types.


## Basic usage

```python
bk.pl.violin(
    adata,
    keys=["CDC20", "ASCL1", "CD3D"],
    groupby="Project_ID",
    figsize=(10, 3),
)
```
## Notes

- If a key is found in `adata.obs`, it is treated as metadata.
- If a key is found in `adata.var_names`, expression is extracted from `layer`.
- Supports mixed metadata + genes in a single call.

```{automodule} bullkpy.pl.violin
:members:
:undoc-members: false
:show-inheritance: false
```