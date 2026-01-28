# BULLKpy

<img src="docs/images/BULLKpy_logo.png" width="300">


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

BULLKpy is designed to integrate seamlessly with the **scverse ecosystem**,  
and to help **standardize and democratize** bulk transcriptomics analysis in Python.

---

## 🚀 Getting started

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
```

⸻

📘 Tutorials & notebooks

```{toctree}
:maxdepth: 2
:caption: Tutorials

notebooks/index
```
⸻

📚 API reference

```{toctree}
:maxdepth: 2
:caption: API reference

api/index
```

🔗 Links

- GitHub: https://github.com/malumbres/BULLKpy
- Issue tracker: https://github.com/malumbres/BULLKpy/issues

---

## ✅ Required companion files (quick checklist)

Make sure these exist, otherwise RTD will silently skip sections:

### 1️⃣ `docs/notebooks/index.md` (or `index.rst`)
Example:

```md
# Tutorials

```{toctree}
:maxdepth: 1

260123_TGA_analysis_in_BULLKpy
```


(Notebook filename **without `.ipynb`**)

---

### 2️⃣ `docs/api/index.rst` (autosummary root)
Example:

```rst
API reference
=============

.. toctree::
   :maxdepth: 2

   pl
   tl
   pp
   io
```

