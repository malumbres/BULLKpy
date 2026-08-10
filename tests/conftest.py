"""Shared fixtures: a small synthetic bulk RNA-seq dataset with clinical metadata."""
from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")  # never open windows during tests

N_SAMPLES = 60
N_GENES = 300


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(0)


@pytest.fixture
def counts_df():
    """Raw integer count matrix, genes x samples (the read_counts default layout)."""
    rng = np.random.default_rng(0)
    genes = [f"GENE{i:04d}" for i in range(N_GENES)]
    # a few recognisable gene families used by qc_metrics
    genes[:5] = ["MT-CO1", "MT-ND1", "MT-ATP6", "RPS3", "RPL7"]
    samples = [f"S{i:03d}" for i in range(N_SAMPLES)]

    base = rng.lognormal(mean=4.0, sigma=1.0, size=(N_GENES, 1))
    lib = rng.uniform(0.6, 1.6, size=(1, N_SAMPLES))
    lam = base * lib
    X = rng.poisson(lam).astype(np.int64)

    # inject real group structure into the first 40 genes so DE/clustering has signal
    half = N_SAMPLES // 2
    X[:40, :half] = (X[:40, :half] * 3.5).astype(np.int64)

    return pd.DataFrame(X, index=genes, columns=samples)


@pytest.fixture
def counts_file(tmp_path, counts_df):
    p = tmp_path / "counts.tsv"
    counts_df.to_csv(p, sep="\t")
    return p


@pytest.fixture
def metadata_file(tmp_path, counts_df):
    rng = np.random.default_rng(1)
    samples = list(counts_df.columns)
    half = len(samples) // 2
    meta = pd.DataFrame(
        {
            "sample_id": samples,
            "Subtype": ["Basal"] * half + ["Luminal"] * (len(samples) - half),
            "Batch": rng.choice(["b1", "b2"], size=len(samples)),
            "Stage": rng.choice(["I", "II", "III"], size=len(samples)),
            "age": rng.integers(35, 85, size=len(samples)),
            "purity": rng.uniform(0.3, 0.95, size=len(samples)),
            "OS_time": rng.uniform(1.0, 100.0, size=len(samples)),
            "OS_event": rng.integers(0, 2, size=len(samples)),
            # binary response label used by the signature/PR helpers
            "resp": rng.choice(["R", "NR"], size=len(samples)),
        }
    )
    p = tmp_path / "meta.tsv"
    meta.to_csv(p, sep="\t", index=False)
    return p


@pytest.fixture
def adata_raw(counts_file, metadata_file):
    """AnnData straight out of io.read_counts + io.add_metadata."""
    import bullkpy as bk

    adata = bk.io.read_counts(counts_file, orientation="genes_by_samples")
    bk.io.add_metadata(adata, metadata_file, index_col="sample_id")
    return adata


@pytest.fixture
def adata(adata_raw):
    """Fully preprocessed AnnData: QC + CPM + log1p, ready for tl/pl functions."""
    import bullkpy as bk

    a = adata_raw
    bk.pp.set_raw_counts(a)
    bk.pp.qc_metrics(a)
    bk.pp.normalize_cpm(a)
    bk.pp.log1p(a)
    return a


@pytest.fixture
def adata_pca(adata):
    """Preprocessed + PCA + neighbors, for embedding/clustering tests."""
    import bullkpy as bk

    bk.tl.pca(adata, n_comps=10)
    bk.tl.neighbors(adata, n_neighbors=8)
    return adata


@pytest.fixture
def genes(adata):
    return list(adata.var_names[:10])
