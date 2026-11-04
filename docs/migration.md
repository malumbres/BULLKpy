# Migration guide

This page lists every public name or parameter that changed in 0.3. Nothing was removed: each old spelling still works and raises a `DeprecationWarning`, and will be removed in 1.0. To see the warnings in a notebook run `import warnings; warnings.simplefilter('default')`.

## Renamed parameters

Same behaviour, one spelling for each concept across the package.

| Old | New | Functions |
|---|---|---|
| `de_key=` | `key=` | `pl.heatmap_de` |
| `event_key=` | `event_col=` | `pl.km_2x2_interaction`, `pl.km_univariate`, `tl.cox_interaction`, `tl.cox_univariate`, `tl.run_cox_per_group`, `tl.surv_1d_bins`, `tl.surv_2x2_bins` |
| `fdr_col=` | `qval_col=` | `pl.gsea_bubbleplot` |
| `group_col=` | `groupby=` | `tl.run_cox_per_group`, `tl.summarize_by_group` |
| `out_key=` | `key_added=` | `pl.km_2x2_interaction`, `pl.km_univariate`, `pl.metaprogram_rank1_composition_stackedbar`, `pp.obs_map_categories`, `tl.cox_interaction`, `tl.cox_univariate`, `tl.gene_metadata_association_scan`, `tl.metaprogram_dispersion_by_group`, `tl.metaprogram_ne_enrichment`, `tl.metaprogram_topk_contribution`, `tl.quantile_groups`, `tl.sample_dispersion`, `tl.signature_axis`, `tl.signature_dispersion_by_group` |
| `p_col=` | `pval_col=` | `pl.volcano` |
| `q_col=` | `qval_col=` | `pl.dotplot_association`, `pl.volcano_categorical`, `pl.volcano` |
| `random_state=` | `seed=` | `tl.cluster`, `tl.leiden_resolution_scan`, `tl.score_genes_cell_cycle`, `tl.score_genes_dict`, `tl.score_genes`, `tl.umap` |
| `store_key=` | `key=` | `pl.ari_resolution_heatmap` |
| `store_key=` | `key_added=` | `tl.gsea_preranked`, `tl.leiden_resolution_scan`, `tl.pca_loadings` |
| `time_key=` | `time_col=` | `pl.km_2x2_interaction`, `pl.km_univariate`, `tl.cox_interaction`, `tl.cox_univariate`, `tl.run_cox_per_group`, `tl.surv_1d_bins`, `tl.surv_2x2_bins` |
| `uns_key=` | `key=` | `pl.cox_forest_from_uns`, `pl.metaprogram_dispersion_heatmap`, `pl.rankplot`, `pl.signature_dispersion_heatmap` |
| `uns_key=` | `key_added=` | `tl.score_metaprograms` |

Concept map: survival time/event are `time_col`/`event_col`; a column to group by is `groupby`; the name under which a result is **written** is `key_added`; the name of a stored result to **read** is `key`; the random seed is `seed`; p- and q-value columns are `pval_col`/`qval_col`.

## Functions that moved

| Old | New | Why |
|---|---|---|
| `pl.bootstrap_pr_bands` | `tl.bootstrap_pr_bands` | it computes bands; the plot is `pl.plot_pr_curves_signatures_bootstrap` |
| `pl.export_leading_edge_clusters_cytoscape` | `tl.export_leading_edge_clusters_cytoscape` | writes files |
| `pl.leading_edge_cluster_driver_genes` | `tl.leading_edge_cluster_driver_genes` | computation, not plotting |
| `pl.leading_edge_pathway_clusters` | `tl.leading_edge_pathway_clusters` | clustering, not plotting |
| `pl.run_cox_per_group` | `tl.run_cox_per_group` | it fits models and stores results; nothing is drawn |
| `tl.load_signature_bank` | `io.load_signature_bank` | reads a file |
| `tl.metaprogram_scores_get` | `get.metaprogram_scores` | it is an accessor |
| `tl.mp_dispersion_metrics` | `tl.metaprogram_dispersion_metrics` | one `metaprogram_` prefix for all metaprogram functions |
| `tl.plot_corr_heatmap` | `pl.plot_corr_heatmap` | it draws a figure |
| `tl.plot_corr_scatter` | `pl.plot_corr_scatter` | it draws a figure |
| `tl.read_metaprograms_xlsx` | `io.read_metaprograms_xlsx` | reads a file |
| `tl.save_signature_bank` | `io.save_signature_bank` | writes a file |
| `tl.score_signature_on_PPV` | `tl.score_signature_on_ppv` | lower-case function name |
| `tl.store_gene_categorical_association` | `tl.gene_categorical_association(key_added=...)` | now `bk.tl.gene_categorical_association(..., key_added=...)` |
| `tl.store_rank_genes_groups_fast` | `tl.rank_genes_groups_fast(key_added=...)` | now `bk.tl.rank_genes_groups_fast(..., key_added=...)` |
| `tl.umap_graph` | `tl.umap(method="graph")` | now `bk.tl.umap(..., method="graph")` |

## New

- `bk.get.de(adata, contrast)` returns a DE results table from either `tl.de` or `tl.de_glm`, with `:` and `_` treated alike.
- `ax=` on 43 single-axes plotting functions, for building multi-panel figures. Clustered heatmaps (`cluster=True`, `dendrogram=True`) create their own figure and raise if `ax=` is also given.
- `bk.settings.reset()`; a warning when a palette name is unknown.

## Behaviour changes (not just names)

- **`tl.de` `log2FC`** was a difference of natural logs; it is now log2, as the column name says. Fold-change cut-offs chosen on old output select different genes.
- **Untestable genes** get `NaN` p/q-values (and do not count towards the number of tests) in `tl.de` and `tl.de_glm`; they were set to p = 1.
- **Missing group labels** (`NaN` in the `groupby` column) are excluded from comparisons. With `reference="rest"` they used to be counted as the comparison group.
- **A named layer that does not exist** raises `KeyError`. Many functions used to fall back to `adata.X` silently.
- **Figures** are saved at `bk.settings.plot_dpi` through one code path (a few used 300 dpi or matplotlib's default).
- **Importing `bullkpy` in Jupyter** no longer runs `%matplotlib inline`.
- **matplotlib >= 3.9** is required.

## Considered and deliberately not changed

These were on the list and were left alone because changing them cannot be done safely with an alias, or would alter results. Each is a candidate for 1.0.

| Item | Why it was left as it is |
|---|---|
| One storage layout and key format for DE results (`uns["de"]["Project_ID_LUAD_vs_LUSC"]` vs `uns["de_glm"]["Subtype:LumB_vs_LumA"]`) | Code and notebooks index `adata.uns` directly, so a new key format cannot be aliased and would break them silently. Use `bk.get.de()`, which accepts both. |
| One return type for plotting functions (now a mix of `(fig, ax)`, `fig`, `ax`, `ClusterGrid`, tuples with a table) | `fig, ax = ...` unpacking of a plot call in user code would break with no way to warn. |
| Merging the four PR-curve functions, `stability_barplot`/`stability_freq_bar`, `rankplot`/`rankplot_association`, `association_heatmap`/`heatmap_association`, the three redundancy filters, `cox_forest`/`cox_forest_from_uns` | They take different inputs or implement different algorithms, so merging changes behaviour. |
| `make_obs_h5ad_safe_strict`, `make_var_h5ad_safe_strict`, `find_bad_*` vs `make_h5ad_safe` | The strict variants coerce more aggressively and return a report; `make_h5ad_safe` is the recommended default. |
| `layer_expr` vs `layer_counts` in `tl.de`, `copy` vs `inplace` | Each pair means two different things. |
| `pl.*` calling `set_style()` (changes global matplotlib settings) | By design; a scoped style would touch 47 call sites. |

