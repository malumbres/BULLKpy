# Rank genes groups (fast)

```{eval-rst}
.. autofunction:: bullkpy.tl.rank_genes_groups_fast
```

## Storing the result

Pass `key_added="rank_genes_groups_fast"` (and optionally `subkey=`) to also store `{"results", "params"}` in `adata.uns`. This replaces the former `store_rank_genes_groups_fast()` function, which still works but is deprecated.
