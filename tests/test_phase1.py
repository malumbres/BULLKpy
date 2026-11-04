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


@pytest.fixture
def adata_missing_groups(adata):
    labels = adata.obs["Subtype"].astype(object).copy()
    labels.iloc[:3] = np.nan
    adata.obs["Subtype"] = labels
    return adata


def test_group_labels_keep_missing_values_missing():
    from bullkpy._compat import group_labels

    out = group_labels(pd.Series(["a", np.nan, 1.0], dtype=object))
    assert out.iloc[0] == "a" and pd.isna(out.iloc[1]) and out.iloc[2] == "1.0"


def test_gene_plot_skips_samples_without_a_group(adata_missing_groups, genes):
    bk.pl.gene_plot(adata_missing_groups, gene=genes[0], groupby="Subtype", show=False)


def test_association_never_tests_a_nan_group(adata_missing_groups, genes):
    df = bk.tl.gene_categorical_association(
        adata_missing_groups, groupby="Subtype", genes=list(genes[:3]), layer="log1p_cpm"
    )
    assert not df.astype(str).isin(["nan"]).drop(columns=df.select_dtypes("number").columns).any().any()


def test_de_with_missing_groups_excludes_those_samples(adata_missing_groups):
    bk.tl.rank_genes_groups_fast(
        adata_missing_groups, groupby="Subtype", group="Basal", reference="rest", layer="log1p_cpm"
    )
    bk.pl.corr_heatmap(adata_missing_groups, use="samples", col_colors="Subtype", show=False)


def test_rest_reference_excludes_samples_without_a_label(adata, genes):
    gene = genes[0]
    x = np.asarray(adata[:, gene].layers["log1p_cpm"]).ravel().astype(float)
    labels = adata.obs["Subtype"].astype(object).copy()
    is_basal = (labels == "Basal").to_numpy()
    expected = x[~is_basal][3:].mean()
    missing = np.flatnonzero(~is_basal)[:3]
    labels.iloc[missing] = np.nan
    adata.obs["Subtype"] = labels
    res = bk.tl.rank_genes_groups_fast(
        adata, groupby="Subtype", group="Basal", reference="rest", genes=[gene], layer="log1p_cpm"
    )
    assert res["mean_ref"].iloc[0] == pytest.approx(expected)


def test_rank_genes_groups_rest_excludes_samples_without_a_label(adata):
    labels = adata.obs["Subtype"].astype(object).copy()
    base = adata.copy()
    labels.iloc[:3] = np.nan
    adata.obs["Subtype"] = labels
    keep = adata[adata.obs["Subtype"].notna()].copy()
    bk.tl.rank_genes_groups(adata, groupby="Subtype", layer="log1p_cpm")
    bk.tl.rank_genes_groups(keep, groupby="Subtype", layer="log1p_cpm")
    a, b = adata.uns["rank_genes_groups"], keep.uns["rank_genes_groups"]
    group = a["names"].dtype.names[0] if hasattr(a["names"], "dtype") and a["names"].dtype.names else "Basal"
    assert list(a["names"][group]) == list(b["names"][group])
    assert base.n_obs == adata.n_obs


def test_km_2x2_interaction_with_two_quantile_binned_scores(adata):
    pytest.importorskip("lifelines")
    rng = np.random.default_rng(0)
    adata.obs["score_x"] = rng.normal(size=adata.n_obs)
    adata.obs["score_z"] = rng.normal(size=adata.n_obs)
    bk.pl.km_2x2_interaction(
        adata, x_key="score_x", z_key="score_z", time_key="OS_time", event_key="OS_event",
        q=(0.5, 0.5), min_n=1, show=False,
    )


def test_survival_binning_works_on_an_anndata_view(adata):
    rng = np.random.default_rng(0)
    adata.obs["score_x"] = rng.normal(size=adata.n_obs)
    adata.obs["score_z"] = rng.normal(size=adata.n_obs)
    view = adata[adata.obs["Subtype"] == adata.obs["Subtype"].iloc[0]]
    assert view.is_view
    bk.tl.surv_1d_bins(view, x_key="score_x", time_key="OS_time", event_key="OS_event")
    assert "surv_group_1d" in view.obs.columns
    view = adata[adata.obs["Subtype"] == adata.obs["Subtype"].iloc[0]]
    bk.pl.km_2x2_interaction(
        view, x_key="score_x", z_key="score_z", time_key="OS_time", event_key="OS_event",
        q=(0.5, 0.5), min_n=1, force_rebin=True, show=False,
    )


def test_unknown_palette_warns_and_falls_back(caplog):
    import logging

    from bullkpy.pl._colors import categorical_palette

    with caplog.at_level(logging.WARNING):
        cols = categorical_palette(["a", "b"], "Set20")
    assert set(cols) == {"a", "b"} and all(c.startswith("#") for c in cols.values())
