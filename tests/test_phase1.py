"""Regression tests for statistical and API correctness fixes."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

import bullkpy as bk
from bullkpy._stats import bh_fdr


def test_bh_fdr_matches_statsmodels():
    p = np.random.default_rng(0).uniform(size=200)
    assert np.allclose(bh_fdr(p), multipletests(p, method="fdr_bh")[1])


def test_bh_fdr_excludes_nan_from_the_number_of_tests():
    q = bh_fdr([0.01, np.nan, 0.04])
    assert np.isnan(q[1])
    assert np.allclose(q[[0, 2]], multipletests([0.01, 0.04], method="fdr_bh")[1])


def _two_group_adata():
    # Group A is exactly 4-fold higher on the count scale for gene g0.
    rng = np.random.default_rng(1)
    x = rng.normal(5, 0.2, size=(40, 4))
    x[:20, 0] += np.log(4.0)
    x[:, 3] = 1.0  # constant: cannot be tested
    a = ad.AnnData(x, var=pd.DataFrame(index=["g0", "g1", "g2", "g3"]))
    a.layers["log1p_cpm"] = x.copy()
    a.obs["grp"] = pd.Categorical(["A"] * 20 + ["B"] * 20)
    return a


def test_de_log2fc_is_on_the_log2_scale():
    a = _two_group_adata()
    bk.tl.de(a, groupby="grp", group="A", reference="B", method="welch_ttest", layer_expr="log1p_cpm")
    res = a.uns["de"]["grp_A_vs_B"]["results"].set_index("gene")
    assert res.loc["g0", "log2FC"] == pytest.approx(2.0, abs=0.15)


def test_de_untestable_gene_is_nan_and_not_counted():
    a = _two_group_adata()
    bk.tl.de(a, groupby="grp", group="A", reference="B", method="welch_ttest", layer_expr="log1p_cpm")
    res = a.uns["de"]["grp_A_vs_B"]["results"].set_index("gene")
    assert np.isnan(res.loc["g3", "pval"]) and np.isnan(res.loc["g3", "qval"])
    tested = res.drop("g3")
    assert np.allclose(tested["qval"], multipletests(tested["pval"], method="fdr_bh")[1])


def test_stored_rank_genes_results_use_the_plotting_layout(adata):
    bk.tl.store_rank_genes_groups_fast(adata, groupby="Subtype", group="Basal", reference="rest", layer="log1p_cpm")
    stored = adata.uns["rank_genes_groups_fast"]["Subtype:Basal_vs_rest"]
    assert isinstance(stored["results"], pd.DataFrame)


def test_settings_reset_restores_defaults():
    bk.settings.plot_dpi = 42
    bk.settings.set_palette_for("x", palette="Set2")
    bk.settings.reset()
    assert bk.settings.plot_dpi == 150
    assert bk.settings.categorical_palette_overrides == {}
    bk.settings.verbosity = 1
    assert bk.settings.verbosity.level == 20


def test_vector_raises_on_missing_layer():
    a = _two_group_adata()
    with pytest.raises(KeyError, match="nope"):
        bk.get.vector(a, "g0", layer="nope")


def test_metaprogram_heterogeneity_skips_missing_groups():
    a = _two_group_adata()
    a.obs["grp"] = a.obs["grp"].astype(object)
    a.obs.iloc[0, a.obs.columns.get_loc("grp")] = np.nan
    a.obsm["X_mp"] = np.random.default_rng(0).normal(size=(40, 3))
    import inspect

    params = inspect.signature(bk.tl.metaprogram_heterogeneity).parameters
    kw = {"obsm_key": "X_mp"} if "obsm_key" in params else {}
    df = bk.tl.metaprogram_heterogeneity(a, groupby="grp", stat="sd", **kw)
    assert set(df[[c for c in df.columns if df[c].astype(str).isin(["A"]).any()][0]]) == {"A", "B"}
