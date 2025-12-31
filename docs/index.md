# BULLKpy

Scanpy-inspired pipeline for bulk RNA-seq analysis built around **AnnData**.

BULLKpy is a Scanpy-inspired Python toolkit for **bulk RNA-seq analysis**,
built on top of **AnnData** and designed for:

- Differential expression analysis
- QC and filtering of bulk RNA-seq
- Clustering and dimensionality reduction
- Gene set enrichment analysis (GSEA)
- Rich publication-quality plots

## Quickstart

```python
import bullkpy as bk

adata = bk.read_counts("counts.tsv")
bk.pp.qc_metrics(adata)
bk.tl.de(adata, groupby="condition", groupA="A", groupB="B")
bk.pl.volcano(adata)
```

```{toctree}
:maxdepth: 2
:caption: Contents

install
api/index
notebooks/index
```