# BULLKpy

```{figure} /_static/BULLKpy_logo.png
:alt: BULLKpy
:width: 400px
:align: left
``` 

Scanpy-inspired pipeline for bulk RNA-seq analysis built around **AnnData**.

BULLKpy is a Scanpy-inspired Python toolkit for **bulk RNA-seq analysis**,
built on top of **AnnData** and designed for:

- QC and filtering of bulk RNA-seq
- Clustering and dimensionality reduction
- Gene set enrichment analysis (GSEA)
- Differential expression analysis
- Oncoprint and other tools for cancer research
- Rich publication-quality plots

## Quickstart

```python
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

```{toctree}
:maxdepth: 2
:caption: Contents

install
api/index
notebooks/index
```