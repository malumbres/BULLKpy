"""bk.datasets.tcga_lung builds the tutorial AnnData from the three published files."""
from __future__ import annotations

import gzip
import os
import shutil

import numpy as np
import pandas as pd
import pytest

import bullkpy as bk
from bullkpy import datasets

SAMPLES = ["TCGA-AA-0001-01A", "TCGA-AA-0002-01A", "TCGA-AA-0003-11A", "TCGA-AA-0004-02A"]
ENSG = ["ENSG00000000001.1", "ENSG00000000002.2", "ENSG00000000003.3", "ENSG00000000004.4"]


@pytest.fixture
def lung_folder(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    counts = pd.DataFrame(rng.integers(0, 500, (4, 4)), index=ENSG, columns=SAMPLES)
    counts.index.name = "Ensembl_ID"
    with gzip.open(tmp_path / "TCGA_lung_counts.tsv.gz", "wt") as fh:
        counts.to_csv(fh, sep="\t")
    pd.DataFrame(
        {
            "ensembl_id": ENSG,
            "gene_name": ["TP53", "Y_RNA", "Y_RNA", None],  # a duplicate and a gene without a symbol
            "gene_type": "protein_coding",
            "chromosome": "1",
        }
    ).to_csv(tmp_path / "TCGA_lung_genes.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"sample": SAMPLES, "Project_ID": ["LUAD", "LUAD", "LUSC", "LUSC"], "OS": [1, 0, 1, 0]}
    ).to_csv(tmp_path / "TCGA_lung_metadata.tsv", sep="\t", index=False)
    # the real checksums belong to the published files; point the module at these small ones
    monkeypatch.setattr(datasets, "_FILES", {n: datasets._sha256(tmp_path / n) for n in datasets._FILES})
    return tmp_path


def test_builds_an_anndata_with_integer_counts_and_clean_names(lung_folder):
    a = bk.datasets.tcga_lung(path=lung_folder)
    assert a.shape == (4, 4)
    assert np.issubdtype(a.X.dtype, np.integer)
    assert a.var_names.is_unique
    assert list(a.var_names) == ["TP53", "Y_RNA", "Y_RNA-1", "ENSG00000000004.4"]  # fallback to the Ensembl ID
    assert list(a.var["ensembl_id"]) == ENSG
    assert list(a.obs["Project_ID"]) == ["LUAD", "LUAD", "LUSC", "LUSC"]
    assert list(a.obs["Sample_type"]) == ["Primary", "Primary", "Normal", "Recurrent"]


def test_primary_only_drops_normal_and_recurrent(lung_folder):
    a = bk.datasets.tcga_lung(path=lung_folder, primary_only=True)
    assert list(a.obs_names) == SAMPLES[:2]


def test_missing_files_in_path_are_named(tmp_path):
    with pytest.raises(FileNotFoundError, match="TCGA_lung_counts"):
        bk.datasets.tcga_lung(path=tmp_path)


def test_download_is_refused_until_a_record_is_set(tmp_path, monkeypatch):
    monkeypatch.setenv("BULLKPY_DATA", str(tmp_path))
    monkeypatch.setattr(datasets, "ZENODO_RECORD", None)
    with pytest.raises(RuntimeError, match="manually"):
        bk.datasets.tcga_lung()


def test_download_goes_through_the_cache_and_checks_the_checksum(lung_folder, tmp_path_factory, monkeypatch):
    cache = tmp_path_factory.mktemp("cache")
    monkeypatch.setenv("BULLKPY_DATA", str(cache))
    monkeypatch.setattr(datasets, "ZENODO_RECORD", "123456")
    seen = []

    def fake_download(url, dest):
        seen.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(lung_folder / dest.name, dest)

    monkeypatch.setattr(datasets, "_download", fake_download)
    a = bk.datasets.tcga_lung()
    assert a.n_obs == 4 and len(seen) == 3 and all("/records/123456/files/" in u for u in seen)
    seen.clear()
    bk.datasets.tcga_lung()  # second call reads the cache
    assert not seen

    # a corrupted download is discarded
    for f in (cache / "tcga_lung").iterdir():
        f.unlink()

    def corrupt(url, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"not the data")

    monkeypatch.setattr(datasets, "_download", corrupt)
    with pytest.raises(RuntimeError, match="Checksum mismatch"):
        bk.datasets.tcga_lung()
    assert not list((cache / "tcga_lung").glob("TCGA_lung_counts*"))


def test_a_zenodo_record_is_set_so_the_download_works_out_of_the_box():
    assert datasets.ZENODO_RECORD and datasets.ZENODO_RECORD.isdigit()


@pytest.mark.skipif(not os.environ.get("BULLKPY_TEST_NETWORK"), reason="set BULLKPY_TEST_NETWORK=1 to download from Zenodo")
def test_the_published_zenodo_files_download_and_match_their_checksums(tmp_path, monkeypatch):
    monkeypatch.setenv("BULLKPY_DATA", str(tmp_path))
    a = bk.datasets.tcga_lung()
    assert a.shape[0] == 1108 and set(a.obs["Project_ID"]) == {"LUAD", "LUSC"}
