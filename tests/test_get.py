"""get: accessors that pull tidy frames out of AnnData."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import bullkpy as bk

LAYER = "log1p_cpm"


def test_obs_df_selects_columns(adata):
    df = bk.get.obs_df(adata, keys=["age", "purity"])
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["age", "purity"]
    assert len(df) == adata.n_obs


def test_vector_returns_gene_expression(adata, genes):
    v = bk.get.vector(adata, key=genes[0], layer=LAYER)
    assert len(v) == adata.n_obs
    expected = np.asarray(adata.layers[LAYER])[:, adata.var_names.get_loc(genes[0])]
    assert np.allclose(np.asarray(v, dtype=float), expected)


def test_vector_returns_obs_column(adata):
    v = bk.get.vector(adata, key="purity")
    assert len(v) == adata.n_obs


def test_rank_genes_groups_df(adata):
    bk.tl.rank_genes_groups(adata, groupby="Subtype", layer=LAYER)
    df = bk.get.rank_genes_groups_df(adata, group="Basal")
    assert isinstance(df, pd.DataFrame) and len(df) > 0


def test_rank_genes_groups_df_all(adata):
    bk.tl.rank_genes_groups(adata, groupby="Subtype", layer=LAYER)
    df = bk.get.rank_genes_groups_df_all(adata)
    assert isinstance(df, pd.DataFrame) and len(df) > 0


def test_rank_genes_groups_df_requires_results(adata):
    with pytest.raises(KeyError):
        bk.get.rank_genes_groups_df(adata, group="Basal")


def test_get_de_finds_results_in_either_store_with_either_separator(adata):
    import bullkpy as bk

    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
    bk.tl.de_glm(adata, formula="~ Subtype", contrast=("Subtype", "Basal", "Luminal"))
    a = bk.get.de(adata, "Subtype_Basal_vs_Luminal")
    b = bk.get.de(adata, "Subtype:Basal_vs_Luminal", key="de")
    c = bk.get.de(adata, "Subtype:Basal_vs_Luminal", key="de_glm")
    assert a.equals(b) and "wald_z" in c.columns
    import pytest

    with pytest.raises(KeyError, match="Available:.*Subtype"):
        bk.get.de(adata, "nope")
