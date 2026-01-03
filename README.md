# BULLKpy 🧬
<img width="617" height="592" alt="Captura de pantalla 2026-01-03 a las 19 13 16" src="https://github.com/user-attachments/assets/d48c768e-ff71-4fce-859c-8ce165bee71d" />

**BULLKpy** is a Python pipeline for **bulk RNA-seq analysis**, inspired by Scanpy but adapted for
bulk transcriptomics. It integrates QC, normalization, clustering, differential expression,
gene set enrichment analysis (GSEA), and rich visualization utilities.

Developed and used for TCGA and large-scale cancer transcriptomics analyses.

---

## 📄 Documentation

BULLKpy documentation in Read The Docs:

https://bullkpy.readthedocs.io/en/latest/    

--- 

## 📦 Project structure

```bash
bullkpy-skeleton/
├── src/                # BULLKpy Python package
│   └── bullkpy/
│       ├── pp/         # preprocessing
│       ├── tl/         # tools (DE, clustering, GSEA, associations)
│       ├── pl/         # plotting
│       ├── io.py
│       └── settings.py
│
├── notebooks/          # analysis notebooks (examples, use cases)
├── data/               # large input datasets (NOT tracked by git)
├── docs/		# Read the Docs at `https://bullkpy.readthedocs.io/en/latest/` 
├── results/            # analysis outputs (NOT tracked by git)
│
├── pyproject.toml      # package configuration
├── README.md
├── LICENSE
└── .gitignore
```
---

## 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy
```

Install in editable mode:
```bash
pip install -e .
```

## 🧪 Typical workflow

```bash
import bullkpy as bk
import pandas
import seaborn as sns

# Load data
adata = bk.read_counts("counts.tsv")

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
bk.tl.cluster(adata, method="leiden")
bk.tl.umap(adata)
bk.tl.umap_graph(adata)
bk.pl.umap(adata)

# Clustering
bk.tl.leiden_resolution_scan(adata)
bk.pl.ari_resolution_heatmap(adata)
bk.tl.cluster(adata)
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
res = bk.tl.de(adata)
bk.tl.de_glm(data)
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

# Plots
bk.pl.violin(adata)
bk.pl.dotplot(adata)
bk.pl.heatmap_de(adata)
bk.pl.sample_distances(adata)
bk.pl.sample_correlation_clustergram(adata)
bk.pl.gene_plot(adata)
bk.pl.oncoprint(adata)
```

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

## ⚠️ Notes

	•	data/ and results/ are not versioned
	•	Designed for small or large datasets (TCGA-scale)
	•	Requires Python ≥ 3.9

## 📄 License

MIT License
