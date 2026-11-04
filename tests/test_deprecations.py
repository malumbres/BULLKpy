"""Renamed parameters and moved functions keep working, with a DeprecationWarning."""
from __future__ import annotations

import importlib
import inspect
import warnings

import numpy as np
import pytest

import bullkpy as bk
from bullkpy._deprecated import MOVED, RENAMED_KWARGS


def _resolve(qualified: str):
    module, _, name = qualified.rpartition(".")
    return getattr(importlib.import_module(module), name)


@pytest.mark.parametrize("qualified", sorted(RENAMED_KWARGS))
def test_signature_has_new_names_and_not_old(qualified):
    func = _resolve(qualified)
    params = set(inspect.signature(func).parameters)
    for old, new in RENAMED_KWARGS[qualified].items():
        assert new in params, f"{qualified}: new parameter {new!r} missing"
        assert old not in params, f"{qualified}: old parameter {old!r} still in the signature"
        assert func.__deprecated_kwargs__[old] == new


def test_renamed_parameters_cover_the_documented_set():
    # Guards against the decorator silently not being applied.
    assert len(RENAMED_KWARGS) >= 35


def test_old_keyword_works_and_warns(adata):
    new = bk.tl.cox_univariate(
        adata, time_col="OS_time", event_col="OS_event", x_keys=["age"], groupby=None, strata=None
    )
    with pytest.warns(DeprecationWarning, match="time_key.*time_col"):
        old = bk.tl.cox_univariate(
            adata, time_key="OS_time", event_key="OS_event", x_keys=["age"], groupby=None, strata=None
        )
    assert type(old) is type(new)


def test_old_and_new_together_is_an_error(adata):
    with pytest.raises(TypeError, match="both"):
        bk.tl.cox_univariate(
            adata, time_key="OS_time", time_col="OS_time", event_col="OS_event", x_keys=["age"]
        )


def test_random_state_is_now_seed(adata_pca):
    bk.tl.cluster(adata_pca, method="kmeans", n_clusters=3, use_rep="X_pca", n_pcs=5, seed=1)
    a = adata_pca.obs["clusters"].astype(str).to_numpy().copy()
    with pytest.warns(DeprecationWarning, match="random_state"):
        bk.tl.cluster(adata_pca, method="kmeans", n_clusters=3, use_rep="X_pca", n_pcs=5, random_state=1)
    assert (adata_pca.obs["clusters"].astype(str).to_numpy() == a).all()


@pytest.mark.parametrize("old", sorted(MOVED))
def test_moved_functions_resolve_and_warn(old):
    namespace, _, name = old.rpartition(".")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        via_old = getattr(getattr(bk, namespace), name)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert via_old is _resolve("bullkpy." + MOVED[old])
    assert name not in getattr(getattr(bk, namespace), "__all__", [])


def test_umap_graph_is_now_a_method_of_umap(adata_pca):
    pytest.importorskip("umap")
    bk.tl.umap(adata_pca, method="graph", n_epochs=20)
    new = adata_pca.obsm["X_umap_graph"].copy()
    assert new.shape == (adata_pca.n_obs, 2)
    del adata_pca.obsm["X_umap_graph"]
    with pytest.warns(DeprecationWarning, match="umap_graph"):
        legacy = bk.tl.umap_graph
    legacy(adata_pca, n_epochs=20)
    np.testing.assert_allclose(adata_pca.obsm["X_umap_graph"], new)
    with pytest.raises(ValueError, match="method"):
        bk.tl.umap(adata_pca, method="nope")


def test_rank_genes_groups_fast_stores_when_key_added_is_given(adata):
    bk.tl.rank_genes_groups_fast(adata, groupby="Subtype", group="Basal", layer="log1p_cpm")
    assert "rank_genes_groups_fast" not in adata.uns
    res = bk.tl.rank_genes_groups_fast(
        adata, groupby="Subtype", group="Basal", layer="log1p_cpm", key_added="rgg"
    )
    stored = adata.uns["rgg"]["Subtype:Basal_vs_rest"]
    assert stored["results"] is res
    assert stored["params"]["groupby"] == "Subtype" and "reference" not in stored["params"]


def test_legacy_store_functions_match_key_added(adata):
    with pytest.warns(DeprecationWarning):
        legacy = bk.tl.store_gene_categorical_association
    legacy(adata, groupby="Subtype", layer="log1p_cpm", genes=list(adata.var_names[:5]))
    assert "gene~Subtype" in adata.uns["associations"]
    assert set(adata.uns["associations"]["gene~Subtype"]) == {"results", "params"}


def test_moved_and_renamed_functions_have_new_homes():
    for name in ("plot_corr_scatter", "plot_corr_heatmap"):
        assert name in bk.pl.__all__ and name not in bk.tl.__all__
    for name in ("run_cox_per_group", "bootstrap_pr_bands", "leading_edge_pathway_clusters"):
        assert name in bk.tl.__all__ and name not in bk.pl.__all__
    assert {"read_metaprograms_xlsx", "save_signature_bank", "load_signature_bank"} <= set(bk.io.__all__)
    assert "metaprogram_scores" in bk.get.__all__
    assert "metaprogram_dispersion_metrics" in bk.tl.__all__ and "score_signature_on_ppv" in bk.tl.__all__
