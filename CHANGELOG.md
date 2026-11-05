# Changelog

All notable changes to **BULLKpy** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [0.3.0] - 2026-11-05

### Highlights
- Three tutorial notebooks on public TCGA lung data (LUAD + LUSC) replace the single pan-cancer notebook:
  1. build the AnnData object, QC and 2-D embedding; 2. clustering, signatures and statistics;
  3. differential expression, pathway enrichment and survival.
- `bk.datasets.tcga_lung()` downloads the tutorial data (integer counts, gene annotation, clinical
  metadata) from Zenodo ([10.5281/zenodo.23966459](https://doi.org/10.5281/zenodo.23966459)), verifies
  checksums and builds the AnnData object.
- A plain `pip install bullkpy` now installs everything the tutorials need (UMAP, Leiden, networkx);
  the old `umap`, `leiden` and `network` extras remain as no-ops. Requires matplotlib >= 3.9,
  seaborn >= 0.13, pandas >= 2.
- Parameter names are consistent across the package and several functions moved to the namespace that
  matches what they do. Every old name keeps working with a `DeprecationWarning`; see `docs/migration.md`.
- Plot layouts that clipped or overlapped (dotplot, dotplot_association, gene_panel_correlation_heatmap,
  volcano, cox_forest_from_uns) were rebuilt, and 43 plotting functions accept `ax=`.
- Tests and CI run on Linux, macOS and Windows.

### Results may change
- `tl.de()` reported `log2FC` as a difference of natural logs, about 1.44 times too
  small on `pp.log1p` layers. It is now on the log2 scale, matching the column name,
  `tl.de_glm()` and `tl.rank_genes_groups()`.
- Genes that cannot be tested now get `NaN` for `pval` and `qval` in `tl.de()` and
  `tl.de_glm()`, and no longer count towards the number of tests. They were previously
  set to `p = 1`, which made q-values slightly too conservative.

- Samples with a missing group label are now excluded from analyses. With
  `reference="rest"` they used to be counted in the comparison group, and on pandas 2
  they could also form a group called `"nan"` in tests and plots.

### Changed
- `pl.dotplot()` lays out in inches from the number of rows/columns and the label lengths, so
  figures stay compact and uncluttered at any size (a 2 x 3 plot no longer gets a 5.8 x 4.2 in canvas,
  a 25 x 30 one is smaller than before). Labels sit under their legend dots, dendrograms line up with
  the dots, and the colour bar and size legend are compact. New `cell_size=` (inches per dot cell).
  A `figsize=` too small for the labels and legends is enlarged, with a warning, instead of squashed.
  The size legend now follows `dot_max` and no longer shows the faint grid behind the dots.

- `pl.cox_forest_from_uns()` leaves empty space around the points and confidence bars (the first and last
  rows no longer touch the frame, and the x range extends past the CIs), puts the p-values in a column
  outside the frame, sizes the figure from the number of rows, wraps and centres the title above the plot,
  and accepts `ax=`. The two near-identical drawing blocks (tidy table and single model) are now one helper.
- `pl.volcano()` puts the legend to the right of the plot (`legend_outside=True`) instead of over the
  points, adds headroom so labels at the extremes stay inside the frame, and moves gene labels apart so
  none overlaps another or leaves the axes (`avoid_label_overlap=True`); a thin line joins a moved label
  to its point. Both can be switched off.
- `pl.dotplot_association()` leaves empty space around the outermost dots (`pad=`, default 0.7 cells),
  caps dots to the size of their cell so none is clipped or overlaps, adds a size legend for
  -log10(q) (it had none), keeps the colour bar to a sensible length, skips missing cells, and drops
  the long default title. A `figsize` too small for the labels, colour bar and legend is enlarged with
  a warning.
- `pl.gene_panel_correlation_heatmap(cluster=True)` places the colour bar to the right of the row
  labels (measured after drawing, widening the figure if needed) instead of at a fixed position that
  long gene names ran under, and it now passes `vmin`/`vmax` to the clustermap. `cbar_pos=` is `None`
  by default (automatic); a tuple still sets it by hand. The title sits above the dendrogram.

### Fixed
- `tl.gene_gene_correlations(batch_mode="meta")` crashed on NumPy 2 (`np.math` was removed).
- `pl.violin()` and `pl.gene_association()` no longer trigger seaborn's "palette without hue"
  FutureWarning; colours are unchanged. `pl.metaprogram_heatmap()` no longer uses a deprecated AnnData call.
- `io.add_metadata()` returns an unfragmented `.obs`, so adding QC metrics or group columns to a wide
  metadata table no longer prints a long "DataFrame is highly fragmented" warning.
- `pl.gene_plot()` and other functions no longer fail on pandas 3 when the grouping
  column has missing values.
- Replaced `boxplot(labels=)`, removed in matplotlib 3.11, in the dispersion and
  metaprogram summary plots. BULLKpy now requires matplotlib 3.9 or newer.
- `pl.km_2x2_interaction()` failed whenever the Z variable needed quantile binning.
- Survival binning and other functions that write to `.obs` now work on AnnData views
  (`adata[mask]` without `.copy()`); new columns were being lost.
- One Benjamini-Hochberg implementation (`_stats.bh_fdr`) replaces eight copies that
  treated `NaN` differently, including an incorrect fallback in `rank_genes_groups`.
- `store_rank_genes_groups_fast()` and `store_gene_categorical_association()` now store
  `{"results", "params"}`, the layout plotting functions read.
- `tl.metaprogram_heterogeneity()` no longer fails when the grouping column has missing
  values; those samples are skipped.
- `get.vector()` raises on a missing layer instead of silently using `.X`.
- Added `settings.reset()`, and made integer verbosity levels mean the same everywhere.
- `lxml` is now installed, as gseapy's MSigDB utilities need it; `scikit-survival` is
  available through the `survival-ml` extra.

### Changed
- Functions that read an expression layer now raise `KeyError` when the named layer does
  not exist. Several used to fall back to `adata.X` silently, which analysed a different
  quantity (for example raw counts instead of log-CPM) without any message.
- All plotting functions save through one helper, so `bk.settings.plot_dpi` applies
  everywhere (several used the matplotlib default, a few a fixed 300 dpi). Library code
  logs through `bullkpy.logging` instead of `print`.
- Removed dead modules (`tl/cluster_metrics.py`, `plotting/`, the shadowed `_settings.py`)
  and duplicate definitions that were never exported.
- Importing `bullkpy` in Jupyter no longer runs `%matplotlib inline` or changes the inline
  backend's settings. An unknown palette name now logs a warning before falling back.
- Added `.zenodo.json` and the author ORCID, so the Zenodo archive gets the right metadata.
- Seven copies of the matrix getter and five of `_as_list` are now one implementation each
  (`bullkpy._matrix`).

### Deprecated (still work, with a `DeprecationWarning`; see docs/migration.md)
- Parameter names normalised across 38 functions: `time_key`/`event_key` -> `time_col`/`event_col`,
  `out_key`/`store_key`/`uns_key` -> `key_added` (written) or `key` (read), `de_key` -> `key`,
  `random_state` -> `seed`, `group_col` -> `groupby`, `q_col`/`fdr_col` -> `qval_col`, `p_col` -> `pval_col`.
- Functions moved to the namespace that matches what they do: `tl.plot_corr_*` -> `pl`;
  `pl.run_cox_per_group`, `pl.bootstrap_pr_bands`, `pl.leading_edge_pathway_clusters`,
  `pl.leading_edge_cluster_driver_genes`, `pl.export_leading_edge_clusters_cytoscape` -> `tl`;
  `tl.read_metaprograms_xlsx`, `tl.save_signature_bank`, `tl.load_signature_bank` -> `io`;
  `tl.metaprogram_scores_get` -> `get.metaprogram_scores`.
- `tl.umap_graph` -> `tl.umap(method="graph")`; `tl.store_gene_categorical_association` and
  `tl.store_rank_genes_groups_fast` -> `key_added=` on the functions they wrapped.
- `tl.mp_dispersion_metrics` -> `tl.metaprogram_dispersion_metrics`;
  `tl.score_signature_on_PPV` -> `tl.score_signature_on_ppv`.

### Added
- `ax=` on 43 single-axes plotting functions. `bk.get.de()` finds DE results from `tl.de` or `tl.de_glm`.
- `docs/migration.md`, generated from the deprecation registry.

### Removed
- `settings.figdir`, `autoshow` and `autosave`, which no code ever read.

### Added
- `scripts/make_tcga_lung_subset.py` builds the TCGA lung tutorial data with integer counts.

## [0.2.0] - 2026-11-02

Repository-wide audit and cleanup ahead of the first fully public release.

### Fixed
- **`io.read_counts()` silently truncated non-integer matrices.** The default
  `dtype="int64"` was applied unconditionally, so a float expression matrix was
  cast without warning: `10.7698` became `10` and `0.9999` became `0`. This
  matters because several public resources — UCSC Xena's GDC Pan-Cancer table
  among them — distribute **log2(count + 1)** values under a `*_counts.tsv`
  filename. Loading one produced a matrix of small integers that looked
  plausible but had lost the data, and every downstream CPM, log1p and
  fold-change was computed from it. `read_counts` now checks whether the values
  are integral before casting, and if they are not it loads them unchanged and
  explains how to recover counts. Genuine integer matrices are still cast.
- `tl.adjusted_rand_index()` was entirely non-functional: it called
  `adjusted_rand_score` without importing it, raising `NameError` on every call.
- Nine `NameError`-at-runtime defects caused by notebook code pasted into the
  library (`bk.` alias used inside modules, `plt`/`ad` used without import,
  and a QC title referencing a variable before it was assigned).
- `pl.gene_panel_correlation_heatmap(cluster=False)` and
  `pl.sample_distances(col_colors=...)` raised `UnboundLocalError`: a
  function-local `import matplotlib.pyplot as plt` inside a conditional branch
  shadowed the module-level import for the whole function.
- `io.add_metadata()` crashed on `.xls`/`.xlsx` input because `low_memory` was
  forwarded to `pandas.read_excel`, which does not accept it. It is now applied
  only to delimited files, where it is valid.
- pandas 3 compatibility: string columns now use the dedicated `str` dtype
  rather than `object`, so `s.dtype == object` no longer identified them.
  `tl.de_glm()` consequently pushed string covariates into the numeric branch
  and failed with `Unable to parse string`. Fixed here and at five other sites
  via a new internal `_compat` helper.
- `pp.qc_metrics()` ignored `compute_pct_mt=False` and computed the metric anyway.
- `pl.qc_by_group()` raised `KeyError` on any dataset lacking `pct_counts_ribo`;
  missing QC metrics are now skipped with a warning.
- `tl.cluster_metrics()` defaulted to reading `obs["leiden"]` while
  `tl.cluster()` writes `obs["clusters"]`, so the documented sequence always
  failed. The defaults now match.
- `reference="rest"` — the spelling used throughout the docs and tutorial — was
  not recognised by `tl.rank_genes_groups_fast()`; only `reference=None` worked.
  Both are now equivalent, unless a real category is named `rest`.
- Removed a duplicate `add_signature_from_genes()` definition in `tl/signature.py`
  where the second silently shadowed the first.
- `src/bullkpy/io.py` used legacy CR-only line endings, which broke diffs and
  text tooling. Normalised, with `.gitattributes` added to prevent recurrence.
- Corrected three broken image references and six references to non-existent
  functions in the API documentation.
- Removed a duplicate PyPI publishing workflow that fired twice per release.
- Replaced `mpl.cm.get_cmap`, removed in matplotlib 3.9, which made plots raise
  `AttributeError` on unknown palette names and in several other code paths.
- Gene lookups now work when `var_names` holds duplicates or `NaN`, unless a
  requested gene is itself ambiguous, in which case the error names it.
  Previously `get_loc` could return a slice or mask instead of one position.
- `pl.heatmap_de()` lists the available contrasts when it cannot find one, and
  the `tl.de()` docstring now gives the real results key format.
- License metadata now matches the BSD 3-Clause `LICENSE` file; `pyproject.toml`,
  the README and `CITATION.cff` wrongly said MIT.

### Breaking
- Parameters that defaulted to column names from one particular study
  (`Project_ID`, `OS.time`, `OS`, `PFS_6m`, `NR_6m`, `Neuroendocrine_score`,
  `mp_heterogeneity_entropy`) no longer carry those defaults. 47 parameters
  across the Cox, signature, dispersion and metaprogram functions are now either
  required keyword arguments or default to `None` where "no grouping" is
  meaningful. A missing-argument `TypeError` is self-explanatory; a wrong default
  produced a `KeyError` about a column the user had never heard of.
  `tl.tcga_define_groups(project_key="Project_ID")` keeps its default, being
  TCGA-specific by definition.

### Added
- `tl.obs_categorical_association()` — tests numeric `.obs` columns for
  association with a categorical variable (the obs-level analogue of
  `tl.gene_categorical_association`). It was referenced by the tutorial, the
  `tl.association` dispatcher and the docs, but had no implementation.
- A test suite (`tests/`) covering the public API, including API/doc contract
  checks that fail when exports, documentation pages and toctrees drift apart.
- `Tests` CI workflow running pytest on Python 3.10–3.13, plus lint, docs build
  and distribution metadata validation.
- Ribosomal QC: `pp.qc_metrics()` now computes `pct_counts_ribo`, which the
  plotting layer already expected.
- Docstrings for the 16 exported functions that had none, including `tl.pca`;
  their Read the Docs pages previously rendered as empty stubs. A contract test
  now fails if any export loses its docstring.
- `CITATION.cff`, so GitHub offers a "Cite this repository" button pointing at
  the bioRxiv preprint.
- `py.typed` marker, so the existing annotations are visible to type checkers in
  downstream projects.
- Notebook contract tests that bind every `bk.*` call in the tutorial against the
  real signatures, catching drift that stored outputs otherwise hide.
- `pp.make_h5ad_safe()` — renames keys that `.h5ad` writing cannot represent.
  HDF5 uses `/` as a group separator, so a column such as
  `GP1_Proliferation/DNA_repair` made `adata.write()` fail with
  *"Forward slashes are not allowed in keys"*. None of the existing sanitise
  helpers fixed this: `find_bad_obs_cols_by_write()` correctly identified the
  offending columns, but nothing renamed them. The new function covers `.obs`,
  `.var`, `.uns` (recursively), `.obsm`, `.varm`, `.layers`, `.obsp` and
  `.varp`, de-duplicates names that would collide after renaming, and reports
  what it changed.

  It also repairs the column dtypes h5py cannot serialise. pandas stores a
  missing string as `float('nan')`, so an `object` column mixing booleans or
  numbers with missing entries fails with *"Can't implicitly convert
  non-string objects to strings"*. Only genuinely un-writable columns are
  touched — anndata handles numeric, `bool`, datetime, `category`, pandas'
  `str` dtype and all-string `object` columns by itself. A column that is
  numeric apart from its gaps becomes numeric; anything else becomes a string
  categorical. Missing values stay missing rather than becoming `"nan"`.

### Changed
- `tl.association()` now covers every combination of gene and `.obs` inputs
  (gene/gene, gene/categorical, gene/numeric, categorical/categorical,
  numeric/categorical, numeric/numeric) instead of raising for most pairs, and
  reports which name it could not resolve.
- `pl.ma()` accepts the results frame as `res` — positionally or by keyword —
  matching `pl.volcano()` and `pl.rankplot()`. The old `result=` keyword still
  works. It also infers the mean-expression column instead of defaulting to
  `mean_norm`, which `tl.de()` does not produce.
- `pl.plot_corr_scatter()` and `pl.plot_corr_heatmap()` are now exported from
  `pl`, where the rest of the plotting API lives. The `tl` names still resolve to
  the same functions.
- Tutorial notebook paths repointed at this machine's `BioDATA` root, with the
  provenance of the pan-cancer matrix documented in the notebook rather than
  assumed present. The notebook now loads it with `dtype=None` and
  back-transforms `2**x - 1` to recover counts, since Xena ships log2(count+1).
- Removed the empty `bullkpy.tools` submodule, which shipped in every wheel.
- `requires-python` raised from `>=3.9` to `>=3.10`, verified against the
  test suite on 3.10, 3.11, 3.12 and 3.13.
- `nbsphinx` removed from runtime dependencies (it is a docs tool, and the docs
  use `myst-nb`); `openpyxl` added, since Excel metadata reading requires it.
- `__version__` is now derived solely from package metadata instead of being
  hardcoded a second time and overriding it.
- Documentation dependencies consolidated into the `docs` extra so Read the Docs
  and local builds no longer drift.
- Sphinx no longer mocks `anndata`, `lifelines` and `gseapy`; mocking the core
  `AnnData` type produced misleading API signatures.
- Private helpers (`_is_integerish`, `_savefig`, `_get_series`) removed from
  public `__all__`; `io` now declares one.
- README rewritten around a runnable quick start, replacing a schematic call
  list that was fenced as shell, contained Python, and could never have run as
  written (`bk.tl.de(adata)` has always required `groupby`). The remaining call
  list is relabelled as a function index. A test now binds every `bk.*` call in
  the README against the real signatures.
- `docs/install.md` documents the Python floor, the optional extras and the
  upgrade path; the README no longer claims Python 3.9 in one place and 3.10 in
  another.
- Dropped four unused image assets (1.6 MB), including a stale 1 MB logo.

---

## [0.1.0] - 2026-01-27

### Added
- First public, documented API for bulk OMICs analysis.
- Integrated AnnData-based workflow for bulk RNA-seq data.
- Core plotting utilities (e.g. correlation heatmaps, oncoprints).
- Initial tutorials and example notebooks.
- ReadTheDocs documentation with API reference.
- Flexible `x_source` / `y_source` resolution in corrplot()

### Changed
- Improved metadata handling and validation.
- More robust plotting defaults and styling.
- Documentation structure reorganized.
- Renamed `corrplot_obs()` → `corrplot()`
- `corrplot()` now supports obs–obs, gene–gene, and gene–obs correlations

### Fixed
- Multiple `.h5ad` serialization edge cases.
- Metadata coercion issues when writing AnnData objects.

---

## [0.0.1] - 2025-12-23

### Added
- Initial internal release.
- Core data structures and early utilities.