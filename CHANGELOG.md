# Changelog

All notable changes to **BULLKpy** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [0.2.0] - 2026-08-10

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