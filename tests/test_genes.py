"""Gene lookups must fail clearly -- and only when needed -- on non-unique var_names."""
from __future__ import annotations

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

import bullkpy as bk
from bullkpy._genes import gene_indexer, gene_loc

matplotlib.use("Agg")


@pytest.fixture
def dup_adata():
    rng = np.random.default_rng(0)
    names = ["TP53", "MKI67", "Y_RNA", "Y_RNA", np.nan, "CDK1", np.nan, "EGFR"]
    a = ad.AnnData(rng.normal(size=(20, len(names))).astype("float32"))
    a.var_names = pd.Index(names, dtype=object)
    return a


def test_indexer_resolves_unique_genes_despite_duplicates(dup_adata):
    idx = gene_indexer(dup_adata, ["TP53", "EGFR", "CDK1", "NOT_A_GENE"])
    assert list(idx) == [0, 7, 5, -1]


def test_indexer_names_the_ambiguous_gene(dup_adata):
    with pytest.raises(ValueError, match=r"Y_RNA") as e:
        gene_indexer(dup_adata, ["TP53", "Y_RNA"])
    assert "var_names_make_unique" in str(e.value)
    assert "missing (NaN)" in str(e.value)


def test_indexer_matches_pandas_on_unique_index():
    a = ad.AnnData(np.ones((3, 4)), var=pd.DataFrame(index=list("abcd")))
    assert list(gene_indexer(a, ["d", "a", "zz"])) == [3, 0, -1]


def test_loc_unique_gene_ok_and_ambiguous_gene_raises(dup_adata):
    assert gene_loc(dup_adata, "CDK1") == 5
    with pytest.raises(ValueError, match="more than once"):
        gene_loc(dup_adata, "Y_RNA")
    with pytest.raises(KeyError):
        gene_loc(dup_adata, "NOT_A_GENE")


def test_public_function_works_for_unique_genes_and_explains_duplicates(dup_adata):
    bk.pl.gene_panel_correlation_heatmap(
        dup_adata, genes=["TP53", "MKI67", "CDK1"], cluster=False, show=False
    )
    with pytest.raises(ValueError, match="Y_RNA"):
        bk.pl.gene_panel_correlation_heatmap(
            dup_adata, genes=["TP53", "Y_RNA", "CDK1"], cluster=False, show=False
        )
