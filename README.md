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
--- 

## 🚀 Tutorials

Use this initial Jupyter notebook for a practical presentation of BULLKpy functions:  

https://bullkpy.readthedocs.io/en/latest/notebooks/BULLKpy_TCGA_example_notebook.html

More functions and field-specific tutorials to come...

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

## 🧪 Typical workflow

```bash
import bullkpy as bk
import pandas
import seaborn as sns
import anndata as ad

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
bk.tl.pca_variance_ratio(adata)
bk.tl.pca_loadings(adata)
bk.pl.pca_loadings_bar(adata)
bk.pl.pca_loadings_heatmap(adata)
bk.tl.neighbors(adata)
bk.tl.umap(adata)
bk.tl.umap_graph(adata)
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
bk.pl.corrplot_obs(adata)
bk.tl.plot_corr_scatter(adata)
bk.tl.gene_categorical_association(adata)
bk.pl.association_heatmap(dfg)
bk.tl.obs_categorical_association(adata)
bk.pl.boxplot_with_stats(adata)
bk.pl.categorical_confusion(adata)
bk.pl.gene_association(adata)
bk.pl.gene_association_volcano(adata)
bk.tl.pairwise_posthoc(y, method="mwu")
bk.tl.cat_cat_association(adata)
bk.pl.dotplot_association(df_all)
bk.pl.heatmap_association(df_all)
bk.tl.rank_genes_categorical(adata)
bk.pl.rankplot_association(dfo)
bk.pl.volcano_categorical(res)
bk.tl.posthoc_per_gene(adata)

# Marker genes and Differential expression
bk.tl.de(adata)
bk.tl.de_glm(data).  # DESeq-like   
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
bk.pl.run_cox_per_group(adata)
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
bk.pp.find_bad_obs_cols_by_write(adata)
bk.pp.find_bad_var_cols_by_write(adata)
bk.pp.make_obs_h5ad_safe_strict(adata)
bbk.pp.make_var_h5ad_safe_strict(adata
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
	•	Requires Python ≥ 3.9

---
## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a full list of changes.

---
## 📄 Publication & License

MIT License

Please cite:   
Malumbres M. (2026) BULLKpy: An AnnData-Inspired Unified Framework for Comprehensive Bulk OMICs Analysis. BioRxiv 10.64898/2026.01.26.701768v1. doi: https://doi.org/10.64898/2026.01.26.701768. 

https://www.biorxiv.org/content/10.64898/2026.01.26.701768v1. 


