# Correlation scatter

```{eval-rst}
.. autofunction:: bullkpy.pl.plot_corr_scatter

```

Plotting counterpart to the correlation utilities in
[`bk.tl`](../tl/index.md).

This function is implemented in `tl/correlations.py`, next to the statistics it
draws, and is exported from both namespaces. `bk.pl.plot_corr_scatter` is the
preferred spelling — it renders a figure, so it belongs with the other plotting
functions — while `bk.tl.plot_corr_scatter` continues to work for existing code.

See the [`bk.tl.plot_corr_scatter`](../tl/plot_corr_scatter.md) page for the full parameter
reference; the two names refer to the same function object.

## See also

- [`bk.pl.corrplot`](corrplot.md)
- [`bk.pl.corr_heatmap`](corr_heatmap.md)
- [`bk.tl.gene_gene_correlations`](../tl/gene_gene_correlations.md)
