# Tutorials

Three notebooks follow one analysis of **TCGA lung cancer RNA-seq** (LUAD and LUSC) from raw counts to survival. Run
them in order; each one starts from the file the previous one wrote. The data come from
`bk.datasets.tcga_lung()`, which downloads them once (see the *Data* section of the README for how the counts were
obtained).

```{toctree}
:maxdepth: 1

Notebook 1: data, QC and 2-D embedding <BULLKpy_nb1>
Notebook 2: clustering, signatures and statistics <BULLKpy_nb2>
Notebook 3: differential expression, pathways and survival <BULLKpy_nb3>
```
