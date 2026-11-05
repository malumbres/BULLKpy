# BULLKpy 🧬
<img src="https://raw.githubusercontent.com/malumbres/BULLKpy/main/docs/images/BULLKpy_logo.png" width="300">

**BULLKpy** is a Python pipeline for **bulk OMICs analysis**, based on AnnData objects and inspired by Scanpy/scverse but adapted for bulk transcriptomics. It integrates QC, normalization, clustering, correlation and association utilities, differential expression, gene set enrichment analysis (GSEA), metaprograms, and rich visualization utilities (oncoprints, etc.).  

---

## 📄 Documentation

BULLKpy documentation in Read The Docs:

https://bullkpy.readthedocs.io/en/latest/    

--- 

## 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy
```

Install from Pypi:  
https://pypi.org/project/bullkpy/  

```bash
pip install bullkpy
```

BULLKpy requires **Python 3.10 or newer**. A plain `pip install bullkpy` installs everything the tutorials
need (UMAP, Leiden clustering, GSEA, survival analysis, networks). The one optional extra is
`bullkpy[survival-ml]` for `pl.time_dependent_roc` (scikit-survival).

Check an installation end to end with [`scripts/smoke_test.py`](scripts/smoke_test.py).

--- 

## 🛠️ Development

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"

pytest                                  # run the test suite
ruff check src/                         # lint
python -m sphinx -b html docs docs/_build/html   # build the docs
```

--- 

## 🚀 Tutorials

The tutorials use public TCGA lung cancer RNA-seq (LUAD and LUSC) and follow one analysis from raw counts to
survival. Run them in order:

1. [Data, QC and 2-D embedding](https://bullkpy.readthedocs.io/en/latest/notebooks/BULLKpy_nb1.html): build the
   AnnData object, quality control, normalisation, PCA and UMAP.
2. [Clustering, signatures and statistics](https://bullkpy.readthedocs.io/en/latest/notebooks/BULLKpy_nb2.html):
   Leiden and k-means, gene signatures, correlations and associations.
3. [Differential expression, pathways and survival](https://bullkpy.readthedocs.io/en/latest/notebooks/BULLKpy_nb3.html):
   DE, GSEA and enrichment, metaprograms, Cox models and Kaplan-Meier curves.

The data are downloaded once by `bk.datasets.tcga_lung()` (see *Data* below), or build them yourself in notebook 1.

--- 
## 🧪 Data

`bk.datasets.tcga_lung()` returns the TCGA lung cancer data used by the tutorials (1,108 open-access samples:
567 LUAD and 541 LUSC) as an AnnData object with **integer counts**. The files are downloaded once from Zenodo
into a cache folder and checked against their SHA-256. Data DOI:
[10.5281/zenodo.23966459](https://doi.org/10.5281/zenodo.23966459) (always the latest version); version 1.0, the one
downloaded by this release, is [10.5281/zenodo.23966460](https://doi.org/10.5281/zenodo.23966460).

```python
import bullkpy as bk

adata = bk.datasets.tcga_lung()                   # all samples
adata = bk.datasets.tcga_lung(primary_only=True)  # primary tumours only
```

*Provenance.* UCSC Xena distributes the GDC Pan-Cancer RNA-seq table (`GDC-PANCAN.htseq_counts.tsv`) as
log2(count + 1). Integer counts were recovered as `2**x - 1` with
[`scripts/make_tcga_lung_subset.py`](scripts/make_tcga_lung_subset.py). Every value was an integer to within
floating-point precision (largest deviation 2.4e-7), except for 10 samples whose values were half-integers; those
were rounded half-to-even and are listed in the dataset record. The recovered counts are consistent with the original
GDC HTSeq counts but were not compared with the raw GDC files. Clinical annotation is from the TCGA Pan-Cancer
Clinical Data Resource (Liu *et al.*, Cell 2018).

--- 
## 📦 Project structure

```bash
bullkpy-skeleton/
├── src/                # BULLKpy Python package
│   └── bullkpy/
|       ├── io.py.      # input/output tools
│       ├── pp/         # preprocessing
│       ├── tl/         # tools (DE, clustering, GSEA, associations)
│       ├── pl/         # plotting
│       └── settings.py
│
├── notebooks/          # analysis notebooks (examples, use cases)
├── data/               # large input datasets (NOT tracked by git)
├── docs/		# Read the Docs at `https://bullkpy.readthedocs.io/en/latest/` 
├── results/            # analysis outputs (NOT tracked by git)
│
├── pyproject.toml      # package configuration
├── README.md
├── CHANGELOG.md
├── LICENSE
├── .gitignore
└── .readthedocs.yaml

```
---

## 🧪 Quick start

```python
import bullkpy as bk

# ---- load -------------------------------------------------------------
adata = bk.io.read_counts("counts.tsv", sep="\t", orientation="genes_by_samples")
bk.io.add_metadata(adata, "metadata.tsv", sep="\t", index_col="sample")

# ---- QC and normalisation --------------------------------------------
bk.pp.set_raw_counts(adata)
bk.pp.qc_metrics(adata)
bk.pl.qc_metrics(adata)
bk.pp.normalize_cpm(adata)
bk.pp.log1p(adata)

# ---- dimensionality reduction and clustering -------------------------
bk.pp.highly_variable_genes(adata, n_top_genes=2000)
bk.tl.pca(adata, n_comps=30, use_highly_variable=True)
bk.pl.pca_scatter(adata, color="Subtype")
bk.tl.neighbors(adata, n_neighbors=15)
bk.tl.cluster(adata, method="leiden", resolution=1.0)   # -> adata.obs["clusters"]
bk.tl.cluster_metrics(adata, true_key="Subtype")

# ---- differential expression -----------------------------------------
bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
res = adata.uns["de"]["Subtype_Basal_vs_Luminal"]["results"]
bk.pl.volcano(res)
bk.pl.heatmap_de(adata, contrast="Subtype_Basal_vs_Luminal", groupby="Subtype")

# ---- associations and survival ---------------------------------------
bk.tl.gene_categorical_association(adata, groupby="Subtype", genes=["MKI67", "TP53"])
bk.tl.cox_gene_association(adata, time_col="OS.time", event_col="OS", genes=["MKI67"])

# ---- save -------------------------------------------------------------
bk.pp.make_h5ad_safe(adata)          # rename keys HDF5 cannot represent
adata.write("results.h5ad", compression="gzip")
```

Every call above is checked against the real signatures by the test suite.

---
## 📚 Function index

A map of what is available, grouped by task. See the
[API reference](https://bullkpy.readthedocs.io/en/latest/api/index.html) for the
arguments each one takes.

```text
# Load data
adata = bk.io.read_counts("counts.tsv", sep="\t")
# Load metadata
adata = bk.add_metadata(adata, "metadata.tsv", sep="\t")
# QC
bk.pp.qc_metrics(adata)
bk.pl.qc_metrics(adata)
bk.pp.filter_genes(adata)
bk.pp.filter_samples(adata)
# PCA + UMAP
bk.pp.highly_variable_genes(adata)
bk.tl.pca(adata)
bk.pl.pca_scatter(adata)
bk.pl.pca_variance_ratio(adata)
bk.tl.pca_loadings(adata)
bk.pl.pca_loadings_bar(adata)
bk.pl.pca_loadings_heatmap(adata)
bk.tl.neighbors(adata)
bk.tl.umap(adata)
bk.tl.umap(adata, method="graph")
bk.pl.umap(adata)
# Clustering
bk.tl.leiden_resolution_scan(adata)
bk.pl.ari_resolution_heatmap(adata)
bk.tl.cluster(adata, method="leiden")
bk.tl.cluster(adata, method="means")
bk.tl.cluster_metrics(adata)
# Genes and gene signatures
bk.tl.score_genes(adata, signature)
bk.tl.score_genes_cell_cycle(adata)
# Correlations and associations
bk.pl.corr_heatmap(adata)
bk.tl.gene_gene_correlations(adata)
bk.tl.gene_gene_correlations(adata)
bk.tl.top_gene_obs_correlations(adata)
bk.tl.obs_obs_corr_matrix(adata)
bk.pl.corrplot(adata)
bk.pl.plot_corr_scatter(adata)
bk.tl.gene_categorical_association(adata)
bk.pl.association_heatmap(dfg)
bk.tl.categorical_association(adata)
bk.pl.boxplot_with_stats(adata)
bk.pl.categorical_confusion(adata)
bk.pl.gene_association(adata)
bk.pl.gene_association_volcano(adata)
bk.tl.pairwise_posthoc(y, method="mwu")
bk.pl.dotplot_association(df_all)
bk.pl.heatmap_association(df_all)
bk.tl.rank_genes_groups_fast(adata)
bk.pl.rankplot_association(dfo)
bk.pl.volcano_categorical(res)
bk.tl.posthoc_per_gene(adata)
# Marker genes and Differential expression
bk.tl.de(adata)
bk.tl.de_glm(adata)  # DESeq-like   
bk.pl.volcano(res)
bk.pl.rankplot(res)
bk.pl.ma(res)
# GSEA, genesets and pathway analysis
bk.tl.gsea_preranked(adata)
bk.pl.gsea_bubbleplot(df_gsea)
bk.pl.gsea_leading_edge_heatmap(adata)
bk.pl.leading_edge_jaccard_heatmap(pre_res)
bk.pl.leading_edge_overlap_matrix(pre_res)
bk.tl.list_enrichr_libraries()
# Metaprograms
bk.tl.score_metaprograms(adata)
bk.tl.metaprogram_heterogeneity(adata)
bk.pl.metaprogram_corr(adata)
bk.pl.metaprogram_heatmap(adata)
bk.tl.metaprogram_sample_metrics(adata)
bk.tl.metaprogram_dispersion_by_group(adata)
bk.tl.metaprogram_ne_enrichment(adata)
bk.tl.metaprogram_topk_contribution(adata)
bk.pl.metaprogram_dispersion_heatmap(adata)
bk.pl.metaprogram_metrics_summary(adata)
bk.pl.metaprogram_ne_scatter(adata)
bk.pl.metaprogram_dominance_ridgeplot_like(adata)
bk.pl.metaprogram_rank1_composition_stackedbar(adata)
# Survival analysis
bk.tl.cox_univariate(adata)
bk.pl.cox_forest_from_uns(adata)
bk.tl.run_cox_per_group(adata)
bk.tl.cox_interaction(adata)
bk.pl.km_univariate(adata)
bk.pl.km_2x2_interaction(adata)
# Plots
bk.pl.violin(adata)
bk.pl.dotplot(adata)
bk.pl.heatmap_de(adata)
bk.pl.sample_distances(adata)
bk.pl.sample_correlation_clustergram(adata)
bk.pl.gene_plot(adata)
bk.pl.oncoprint(adata)
# Other Utilities
bk.pp.make_h5ad_safe(adata)
bk.pp.find_bad_obs_cols_by_write(adata)
bk.pp.find_bad_var_cols_by_write(adata)
bk.pp.make_obs_h5ad_safe_strict(adata)
bk.pp.make_var_h5ad_safe_strict(adata)
bk.pp.batch_correct_combat(adata)
```

---
## 📊 Features

	•	Bulk RNA-seq, small and large projects. QC & filtering
	•	PCA, UMAP, Leiden, k-means clustering
	•	Gene scores and signatures
	•	Gene–obs and obs–obs associations and correlations
	•	Differential expression from counts or log data
	•	GSEA preranked pipeline (GSEApy)
	•	Leading-edge GSEA analysis
	•	Oncoprint-style mutation plots
	•	Scanpy-like API (pp, tl, pl)
---
## ⚠️ Notes

	•	data/ and results/ are not versioned
	•	Designed for small or large datasets (TCGA-scale)
	•	Requires Python ≥ 3.10

---
## Upgrading to 0.3.0

0.3.0 renames parameters for consistency and moves some functions to the namespace that matches what they do. **Nothing
was removed**: every old name still works and raises a `DeprecationWarning`. The full table is in the
[migration guide](https://bullkpy.readthedocs.io/en/latest/migration.html). Behaviour changes that can alter results:

- `tl.de()` reported `log2FC` on the natural-log scale (about 1.44 times too small); it is now log2.
- Genes that cannot be tested get `NaN` p/q-values and no longer count towards the number of tests.
- Samples with a missing group label are excluded from comparisons.
- A named layer that does not exist raises an error instead of silently using `adata.X`.
- Requires matplotlib >= 3.9, seaborn >= 0.13 and pandas >= 2.

---
## Upgrading to 0.2.0

0.2.0 fixes several defects that silently produced wrong results, so re-running
existing analyses is worthwhile. Two changes affect existing code:

- **`io.read_counts()` no longer casts non-integer matrices to `int64`.** It
  previously truncated them without warning (`10.77` became `10`), which
  corrupts any file distributed as log2(count+1) — UCSC Xena's tables among
  them. Such files are now loaded unchanged, with a warning.
- **Parameters that defaulted to one study's column names are now required.**
  Calls relying on defaults such as `time_key="OS.time"`, `groupby="Project_ID"`
  or `label_col="PFS_6m"` must now name the column explicitly. A missing
  argument raises `TypeError` naming the parameter.

See [CHANGELOG.md](CHANGELOG.md) for the full list.

---
## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a full list of changes.

---
## 📄 Publication & License

BSD 3-Clause License

Please cite the software (use *Cite this repository* on GitHub, or the Zenodo record of the release you used) and the preprint:
Malumbres M. (2026) BULLKpy: An AnnData-Inspired Unified Framework for Comprehensive Bulk OMICs Analysis. BioRxiv 10.64898/2026.01.26.701768v1. doi: https://doi.org/10.64898/2026.01.26.701768. 

https://www.biorxiv.org/content/10.64898/2026.01.26.701768v1. 


