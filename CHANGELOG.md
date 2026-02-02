# Changelog

All notable changes to **BULLKpy** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

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