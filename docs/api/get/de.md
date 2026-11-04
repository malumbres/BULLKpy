# de

```{eval-rst}
.. autofunction:: bullkpy.get.de
```

Return the results table of one differential-expression contrast, whichever function produced it.

`bk.tl.de` stores results in `adata.uns["de"]` under keys such as `Project_ID_LUAD_vs_LUSC`, while
`bk.tl.de_glm` uses `adata.uns["de_glm"]` and keys such as `Subtype:LumB_vs_LumA`. `bk.get.de` searches
both stores and treats `:` and `_` as the same separator.

## Examples

```python
bk.tl.de(adata, groupby="Project_ID", group="LUAD", reference="LUSC")
res = bk.get.de(adata, "Project_ID_LUAD_vs_LUSC")
res = bk.get.de(adata, "Project_ID:LUAD_vs_LUSC")   # same table

bk.tl.de_glm(adata, formula="~ Subtype", contrast=("Subtype", "LumB", "LumA"))
res = bk.get.de(adata, "Subtype:LumB_vs_LumA")
```

If nothing matches, the `KeyError` lists every contrast that exists.
