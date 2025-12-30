# BULLKpy 🧬

**BULLKpy** is a Python pipeline for **bulk RNA-seq analysis**, inspired by Scanpy but adapted for
bulk transcriptomics. It integrates QC, normalization, clustering, differential expression,
gene set enrichment analysis (GSEA), and rich visualization utilities.

Developed and used for TCGA and large-scale cancer transcriptomics analyses.

---

## 📦 Project structure

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
├── results/            # analysis outputs (NOT tracked by git)
│
├── pyproject.toml      # package configuration
├── README.md
├── LICENSE
└── .gitignore

---

## 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy


Install in editable mode:

pip install -e .


🧪 Typical workflow
import bullkpy as bk

# Load data
adata = bk.read_counts("counts.tsv")

# QC
bk.pp.qc_metrics(adata)
bk.pl.qc_metrics(adata)

# PCA + clustering
bk.tl.pca(adata)
bk.tl.cluster(adata, method="leiden")

# Differential expression
res = bk.tl.de(
    adata,
    groupby="Project_ID",
    groupA="LUAD",
    groupB="LUSC",
)

# Volcano plot
bk.pl.volcano(res)

# GSEA
df_gsea, pre_res = bk.tl.gsea_preranked(
    adata,
    res=res,
    gene_sets=["Hallmark_2020"],
)
bk.pl.gsea_bubbleplot(df_gsea)

📊 Features

	•	Bulk RNA-seq QC & filtering
	•	PCA, UMAP, Leiden clustering
	•	Differential expression
	•	GSEA preranked pipeline (GSEApy)
	•	Gene–obs and obs–obs associations
	•	Leading-edge GSEA analysis
	•	Oncoprint-style mutation plots
	•	Scanpy-like API (pp, tl, pl)

⚠️ Notes

	•	data/ and results/ are not versioned
	•	Designed for large datasets (TCGA-scale)
	•	Requires Python ≥ 3.9

📄 License
MIT License