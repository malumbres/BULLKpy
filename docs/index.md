# BULLKpy 🧬

<img src="https://raw.githubusercontent.com/malumbres/BULLKpy/main/docs/images/BULLKpy_logo.png" width="300">

**BULLKpy** is a Python framework for **comprehensive bulk OMICs data analysis**,  
with a strong focus on **biomedical and cancer research**.

It provides a unified, AnnData-inspired workflow to perform:

- Quality control and preprocessing
- Dimensionality reduction and clustering
- Differential expression analysis
- Pathway and gene set enrichment
- Metaprograms and tumor heterogeneity analysis
- Survival analysis and clinical associations
- Publication-ready visualization

[BULLKpy on GitHub](https://github.com/malumbres/BULLKpy)
[BULLKpy on Pypi](https://pypi.org/project/bullkpy/)

BULLKpy is based on AnnData structures and is designed to integrate seamlessly with the **scverse ecosystem**,   and to help **standardize and democratize** bulk OMICs analysis in Python.

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

## 🚀 Getting started

### 📘 Table of contents
```{toctree}
:maxdepth: 2
:caption: Contents

install
api/index
```

## 🚀 Tutorials

Step-by-step tutorials

```{toctree}
:maxdepth: 2
:caption: Tutorial

notebooks/index
```

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

## 📘 Quickstart

```bash
import bullkpy as bk
import pandas as pd                        # for better manipulation of metadata
import seaborn as sns                      # for some plots
import gseapy as gspy                      # Additional GSEA tools

## Data input and loading
adata = bk.io.read_counts("counts.tsv").      # when using counts
bk.io.add_metadata(adata, metadata)

## Preprocessing
bk.pp.qc_metrics(adata)
bk.pl.qc_metrics(
    adata,
    color="n_genes_detected",
    vars_to_plot=("n_genes_detected",),
)
bk.pp.filter_samples(adata)
bk.pp.filter_genes(adata)

## PCA and bidimensional representation
bk.pp.highly_variable_genes(adata, layer="log1p_cpm", n_top_genes=2000)
bk.tl.pca(adata, layer="log1p_cpm", n_comps=20, 
         use_highly_variable=True,)  

bk.pl.pca_scatter(adata, color="Project_ID", palette="tab20", figsize=(7,7), point_size=8)
bk.pl.pca_variance_ratio(adata)
bk.tl.pca_loadings(adata, pcs=[1,2], n_top=10)

### Neighbors and UMAP
bk.tl.neighbors(adata, n_neighbors=15, n_pcs=20, metric="euclidean")
bk.tl.umap(adata, n_neighbors=8, n_pcs=30, min_dist=0.3, random_state=0)
bk.pl.umap(adata, color="Project_ID", 
           figsize=(8,7),
           point_size=10, palette="tab20")

## Clustering and groups
bk.tl.leiden_resolution_scan(adata=
bk.tl.cluster(adata, method="leiden", key_added="leiden_2.2", resolution=2.2)

## Genes and signatures
bk.tl.score_genes(adata, ne_score_25, score_name="NE25_score", layer="log1p_cpm")

## Data exploration
bk.pl.violin(
    adata,
    keys=["MYC_amplified", "n_genes_detected"],
    groupby="Project_ID",
)

bk.pl.corr_heatmap(
    adata_lung,
    layer="log1p_cpm",
    method="pearson",
    groupby="Project_ID",
)

## Markers and Differential Expression
bk.tl.de(adata, groupby="condition", groupA="A", groupB="B")
bk.pl.volcano(adata)

##Pathway and Gene Set Enrichment Analysis
bk.tl.gsea_preranked(
    adata,
    res=res,
    comparison=comparison,
    score_col="t",                                       
    gene_sets=["hallmark"], 
)
```
⸻

## 🔗 Links

BULLKpy is available on GitHub ([https://github.com/malumbres/BULLKpy](https://github.com/malumbres/BULLKpy)).

Issue tracker: ([https://github.com/malumbres/BULLKpy/issues](https://github.com/malumbres/BULLKpy/issues))

---
## 📄 Citation

Please refer to: Malumbres M. BULLKpy: An AnnData-Inspired Unified Framework for Comprehensive Bulk OMICs Analysis. BioRxiv 10.64898/2026.01.26.701768v1. doi: https://doi.org/10.64898/2026.01.26.701768. 

https://www.biorxiv.org/content/10.64898/2026.01.26.701768v1. 

