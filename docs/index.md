#  🧬 BULLKpy 🧬

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

BULLKpy is based on **AnnData structures** and is designed to integrate seamlessly with the **scverse ecosystem**,   and to help **standardize and democratize** bulk OMICs analysis in Python.

---

## 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy
```

Install from Pypi:   
([https://pypi.org/project/bullkpy/](https://pypi.org/project/bullkpy/)) 

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
migration
api/index
```
---

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

## 🔗 Links

BULLKpy is available on GitHub ([https://github.com/malumbres/BULLKpy](https://github.com/malumbres/BULLKpy)).

Issue tracker: ([https://github.com/malumbres/BULLKpy/issues](https://github.com/malumbres/BULLKpy/issues))

**Malumbreslab.org**: ([http://malumbreslab.org/](http://malumbreslab.org/))

---
## 📄 Citation

Software (Zenodo, all versions): [10.5281/zenodo.23966752](https://doi.org/10.5281/zenodo.23966752); version 0.3.0: [10.5281/zenodo.23966753](https://doi.org/10.5281/zenodo.23966753).

Preprint:  
Malumbres M. (2026) BULLKpy: An AnnData-Inspired Unified Framework for Comprehensive Bulk OMICs Analysis. BioRxiv 10.64898/2026.01.26.701768v1. doi: https://doi.org/10.64898/2026.01.26.701768. 

**BioRxiv**: ([https://www.biorxiv.org/content/10.64898/2026.01.26.701768v1](https://www.biorxiv.org/content/10.64898/2026.01.26.701768v1)

