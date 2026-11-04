# tcga_lung

```{eval-rst}
.. autofunction:: bullkpy.datasets.tcga_lung
```

The TCGA lung cancer data used by the three tutorials: 1,108 open-access samples (567 LUAD and 541 LUSC) with
**integer counts**, gene annotation and the clinical metadata.

## Where the counts come from

UCSC Xena distributes the GDC Pan-Cancer RNA-seq table (`GDC-PANCAN.htseq_counts.tsv`) as
log2(count + 1). Integer counts were recovered as `2**x - 1` by `scripts/make_tcga_lung_subset.py` in the
repository. After the back-transform every value was an integer to within floating-point precision
(largest deviation 2.4e-7), except for 10 samples whose values were all half-integers; those were rounded half
to even, and their barcodes are listed in `TCGA_lung_half_integer_samples.txt` in the data record (you can drop them
without changing any conclusion). The data are therefore consistent with the original GDC HTSeq counts, but were
not compared with the raw GDC files.

## Files

| File | Contents |
| --- | --- |
| `TCGA_lung_counts.tsv.gz` | integer counts, genes x samples, Ensembl IDs |
| `TCGA_lung_genes.tsv` | Ensembl ID, gene name, gene type, chromosome |
| `TCGA_lung_metadata.tsv` | 387 clinical and annotation columns (listed in `scripts/tcga_lung_metadata_columns.txt`) |
| `TCGA_lung_half_integer_samples.txt` | barcodes of the 10 samples whose counts were rounded from half-integers |

The files are downloaded once into `~/.cache/bullkpy` (set `BULLKPY_DATA` to change it) and verified with SHA-256.

## Examples

```python
import bullkpy as bk

adata = bk.datasets.tcga_lung()                     # downloads once, then reads the cache
adata = bk.datasets.tcga_lung(primary_only=True)    # primary tumours only
adata = bk.datasets.tcga_lung(path="my/folder")     # from files you already have
```

`Sample_type` (Primary, Recurrent, Normal) is derived from the sample barcode. Use primary tumours for tumour
comparisons and for survival analysis, since survival is recorded per patient.
