#!/usr/bin/env python
"""End-to-end smoke test of an installed BULLKpy, on synthetic data.

Run after ``pip install bullkpy`` (or an installed wheel) to check that the package and the
dependencies that ship with it work together: load counts, QC, normalise, PCA, UMAP, Leiden,
differential expression, a plot and an ``.h5ad`` round trip.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import bullkpy as bk  # noqa: E402


def main() -> None:
    rng = np.random.default_rng(0)
    n_genes, n_samples = 300, 60
    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:03d}" for i in range(n_samples)]
    base = rng.lognormal(4.0, 1.0, size=(n_genes, 1))
    counts = rng.poisson(base * rng.uniform(0.6, 1.6, size=(1, n_samples))).astype(np.int64)
    counts[:40, : n_samples // 2] = (counts[:40, : n_samples // 2] * 3.5).astype(np.int64)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        counts_path = tmp / "counts.tsv"
        pd.DataFrame(counts, index=genes, columns=samples).to_csv(counts_path, sep="\t")
        meta_path = tmp / "meta.tsv"
        pd.DataFrame(
            {
                "sample": samples,
                "group": ["A"] * (n_samples // 2) + ["B"] * (n_samples - n_samples // 2),
                "time": rng.uniform(1, 100, n_samples),
                "event": rng.integers(0, 2, n_samples),
            }
        ).to_csv(meta_path, sep="\t", index=False)

        adata = bk.io.read_counts(counts_path, orientation="genes_by_samples")
        bk.io.add_metadata(adata, meta_path, index_col="sample")
        bk.pp.set_raw_counts(adata)
        bk.pp.qc_metrics(adata)
        bk.pp.normalize_cpm(adata)
        bk.pp.log1p(adata)
        bk.tl.pca(adata, n_comps=10)
        bk.tl.neighbors(adata, n_neighbors=8)
        bk.tl.umap(adata)
        bk.tl.cluster(adata, method="leiden", resolution=0.8)
        bk.tl.de(adata, groupby="group", group="A", reference="B")
        res = bk.get.de(adata, "group_A_vs_B")
        assert (res["qval"] < 0.05).sum() > 10, "differential expression found nothing"

        fig, _ = bk.pl.volcano(res, show=False)
        fig.savefig(tmp / "volcano.png")
        bk.pp.make_h5ad_safe(adata)
        adata.write(tmp / "out.h5ad")
        assert ad.read_h5ad(tmp / "out.h5ad").shape == adata.shape

    print(f"bullkpy {bk.__version__}: smoke test passed")


if __name__ == "__main__":
    sys.exit(main())
