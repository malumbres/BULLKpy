"""tl: dimensionality reduction, clustering, DE, correlations, associations, survival."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import bullkpy as bk

LAYER = "log1p_cpm"


# ------------------------------------------------------- embeddings / graph
def test_pca_populates_obsm_varm_uns(adata):
    bk.tl.pca(adata, n_comps=10)
    assert adata.obsm["X_pca"].shape == (adata.n_obs, 10)
    assert adata.varm["PCs"].shape[1] == 10
    assert "variance_ratio" in adata.uns["pca"]


def test_pca_variance_ratio_is_decreasing(adata):
    bk.tl.pca(adata, n_comps=10)
    vr = np.asarray(adata.uns["pca"]["variance_ratio"])
    assert np.all(np.diff(vr) <= 1e-12)
    assert 0 < vr.sum() <= 1.0 + 1e-9


def test_pca_loadings_returns_per_pc_tables(adata_pca):
    out = bk.tl.pca_loadings(adata_pca)
    assert isinstance(out, dict) and out
    assert "PC1_pos" in out and "PC1_neg" in out
    assert isinstance(out["PC1_pos"], pd.DataFrame)
    assert {"gene", "loading", "rank"} <= set(out["PC1_pos"].columns)


def test_neighbors_builds_graph(adata_pca):
    assert "connectivities" in adata_pca.obsp
    assert adata_pca.obsp["connectivities"].shape == (adata_pca.n_obs, adata_pca.n_obs)


def test_cluster_assigns_labels(adata_pca):
    bk.tl.cluster(adata_pca, resolution=1.0)
    assert "clusters" in adata_pca.obs.columns
    assert adata_pca.obs["clusters"].notna().all()


def test_umap_embedding_shape(adata_pca):
    pytest.importorskip("umap")
    bk.tl.umap(adata_pca)
    assert adata_pca.obsm["X_umap"].shape == (adata_pca.n_obs, 2)


# --------------------------------------------------------------------- DE
def test_de_stores_results_under_contrast_key(adata):
    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
    key = "Subtype_Basal_vs_Luminal"
    assert key in adata.uns["de"]
    df = adata.uns["de"][key]["results"]
    assert isinstance(df, pd.DataFrame)
    assert {"gene", "log2FC", "pval", "qval"} <= set(df.columns)


def test_de_recovers_planted_signal(adata):
    """The fixture inflates the first 40 genes in the Basal half."""
    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
    df = adata.uns["de"]["Subtype_Basal_vs_Luminal"]["results"]

    planted = set(adata.var_names[:40])
    top = df.reindex(df["log2FC"].abs().sort_values(ascending=False).index).head(40)
    hits = len(planted & set(top["gene"].astype(str)))
    assert hits >= 25, f"only {hits}/40 planted DE genes in the top 40 by |log2FC|"


def test_de_planted_genes_are_significant(adata):
    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
    df = adata.uns["de"]["Subtype_Basal_vs_Luminal"]["results"].set_index("gene")
    planted = [g for g in adata.var_names[:40] if g in df.index]
    assert (df.loc[planted, "qval"] < 0.05).mean() > 0.8


def test_de_glm_handles_string_dtype_covariates(adata):
    """pandas 3 stores string columns as `str`, not `object`; de_glm must still
    treat them as categorical rather than coercing them to numbers."""
    assert not pd.api.types.is_numeric_dtype(adata.obs["Subtype"])
    bk.tl.de_glm(adata, formula="~ Subtype", contrast=("Subtype", "Basal", "Luminal"))
    assert "de_glm" in adata.uns


def test_rank_genes_groups_stores_results(adata):
    bk.tl.rank_genes_groups(adata, groupby="Subtype", layer=LAYER)
    assert "rank_genes_groups" in adata.uns


def test_rank_genes_groups_fast_returns_ranked_frame(adata):
    df = bk.tl.rank_genes_groups_fast(adata, groupby="Subtype", group="Basal", layer=LAYER)
    assert isinstance(df, pd.DataFrame) and len(df) > 0
    assert {"gene", "pval", "qval"} <= set(df.columns)
    assert df["qval"].is_monotonic_increasing or df["pval"].is_monotonic_increasing


# ------------------------------------------------------- cluster comparison
def test_adjusted_rand_index_identical_labels_is_one(adata):
    adata.obs["copy_of_subtype"] = adata.obs["Subtype"].astype(str)
    ari = bk.tl.adjusted_rand_index(adata, true_key="Subtype", pred_key="copy_of_subtype")
    assert ari == pytest.approx(1.0)


def test_adjusted_rand_index_rejects_unknown_key(adata):
    with pytest.raises(KeyError):
        bk.tl.adjusted_rand_index(adata, true_key="Subtype", pred_key="nope")


def test_cluster_metrics_default_key_matches_cluster_output(adata_pca):
    """cluster() writes obs['clusters']; cluster_metrics() must read the same key."""
    bk.tl.cluster(adata_pca, resolution=1.0)
    out = bk.tl.cluster_metrics(adata_pca, true_key="Subtype")
    assert isinstance(out, dict) and out
    assert {"ari", "nmi"} <= set(out)


def test_categorical_confusion(adata):
    out = bk.tl.categorical_confusion(adata, key1="Subtype", key2="Batch")
    assert isinstance(out, dict)


def test_leiden_resolution_scan(adata_pca):
    df = bk.tl.leiden_resolution_scan(
        adata_pca, true_key="Subtype", resolutions=(0.5, 1.0), n_pcs=10
    )
    assert isinstance(df, pd.DataFrame) and len(df) == 2


# ------------------------------------------------------------- gene scoring
def test_score_genes_adds_obs_column(adata, genes):
    bk.tl.score_genes(adata, genes, score_name="my_score", layer=LAYER)
    assert "my_score" in adata.obs.columns
    assert adata.obs["my_score"].notna().any()


def test_score_genes_dict_adds_one_column_per_set(adata, genes):
    bk.tl.score_genes_dict(adata, {"setA": genes[:5], "setB": genes[5:]}, layer=LAYER)
    assert "setA_score" in adata.obs.columns
    assert "setB_score" in adata.obs.columns


def test_signature_score(adata, genes):
    w = pd.DataFrame({"gene": genes, "beta": np.ones(len(genes))})
    bk.tl.signature_score(adata, weights=w, layer=LAYER)
    assert "signature_score" in adata.obs.columns


# ------------------------------------------------------------ correlations
def test_gene_gene_correlations(adata, genes):
    df = bk.tl.gene_gene_correlations(adata, gene=genes[0], genes=genes, layer=LAYER)
    assert isinstance(df, pd.DataFrame) and len(df) > 0


def test_top_gene_obs_correlations(adata, genes):
    df = bk.tl.top_gene_obs_correlations(
        adata, gene=genes[0], obs=["purity", "age"], layer=LAYER
    )
    assert isinstance(df, pd.DataFrame)


def test_obs_obs_corr_matrix(adata):
    df = bk.tl.obs_obs_corr_matrix(adata, focus="purity", against=["age", "OS_time"])
    assert isinstance(df, pd.DataFrame)


def test_partial_corr_returns_triple(adata):
    r, p, n = bk.tl.partial_corr(
        adata,
        x=adata.obs["age"].to_numpy(float),
        y=adata.obs["purity"].to_numpy(float),
        covariates=["OS_time"],
    )
    assert -1.0 <= r <= 1.0
    assert 0.0 <= p <= 1.0
    assert n > 0


# ------------------------------------------------------------ associations
def test_categorical_association(adata):
    out = bk.tl.categorical_association(adata, key1="Subtype", key2="Batch")
    assert isinstance(out, dict) and "table" in out


def test_gene_categorical_association(adata, genes):
    df = bk.tl.gene_categorical_association(
        adata, groupby="Subtype", genes=genes, layer=LAYER
    )
    assert isinstance(df, pd.DataFrame) and len(df) == len(genes)


@pytest.mark.parametrize(
    "x,y",
    [
        ("Subtype", "Batch"),      # categorical x categorical
        ("Subtype", "purity"),     # categorical x numeric obs
        ("purity", "Subtype"),     # numeric obs x categorical (order independent)
        ("age", "purity"),         # numeric x numeric
    ],
)
def test_association_dispatcher_covers_obs_quadrants(adata, x, y):
    out = bk.tl.association(adata, x=x, y=y)
    assert isinstance(out, (pd.DataFrame, dict))


@pytest.mark.parametrize("swap", [False, True])
def test_association_dispatcher_handles_genes(adata, genes, swap):
    x, y = (genes[0], "Subtype") if not swap else ("Subtype", genes[0])
    out = bk.tl.association(adata, x=x, y=y, layer=LAYER)
    assert isinstance(out, pd.DataFrame) and len(out) > 0


def test_association_dispatcher_gene_vs_numeric_obs(adata, genes):
    out = bk.tl.association(adata, x=genes[0], y="purity", layer=LAYER)
    assert isinstance(out, pd.DataFrame)


def test_association_dispatcher_rejects_unknown_names(adata):
    with pytest.raises(KeyError, match="neither a gene"):
        bk.tl.association(adata, x="not_a_thing", y="Subtype")


def test_posthoc_per_gene(adata, genes):
    df = bk.tl.posthoc_per_gene(adata, gene=genes[0], groupby="Stage", layer=LAYER)
    assert isinstance(df, pd.DataFrame)


def test_pairwise_posthoc(adata):
    df = pd.DataFrame(
        {"grp": adata.obs["Stage"].astype(str).to_numpy(), "y": adata.obs["purity"].to_numpy()}
    )
    out = bk.tl.pairwise_posthoc(df)
    assert isinstance(out, pd.DataFrame)


def test_gene_metadata_association_scan(adata, genes):
    df = bk.tl.gene_metadata_association_scan(
        adata, metadata_key="purity", genes=genes, layer=LAYER
    )
    assert isinstance(df, pd.DataFrame)


# --------------------------------------------------------------- survival
def test_cox_gene_association(adata, genes):
    df = bk.tl.cox_gene_association(
        adata, time_col="OS_time", event_col="OS_event", genes=genes, layer=LAYER
    )
    assert isinstance(df, pd.DataFrame) and len(df) > 0
    assert {"gene", "pval"} <= set(df.columns)


def test_cox_multivariate(adata, genes):
    df = bk.tl.cox_multivariate(
        adata, genes=genes, time_col="OS_time", event_col="OS_event", layer=LAYER
    )
    assert isinstance(df, pd.DataFrame)


def test_cox_univariate(adata):
    out = bk.tl.cox_univariate(
        adata,
        time_key="OS_time",
        event_key="OS_event",
        x_keys=["age", "purity"],
        groupby=None,
        strata=None,
    )
    assert out is not None


def test_surv_1d_bins_adds_group_column(adata):
    bk.tl.surv_1d_bins(adata, x_key="purity", time_key="OS_time", event_key="OS_event")
    assert "surv_group_1d" in adata.obs.columns


def test_surv_2x2_bins_adds_group_column(adata):
    bk.tl.surv_2x2_bins(
        adata, x_key="purity", z_key="age", time_key="OS_time", event_key="OS_event"
    )
    assert "surv_group_2x2" in adata.obs.columns


# --------------------------------------------------------------- signature
def test_filter_genes_var(adata):
    out = bk.tl.filter_genes_var(adata, label_col="resp", label="NR", layer=LAYER)
    assert out is not None


def test_rank_genes_univariate_pr(adata, genes):
    out = bk.tl.rank_genes_univariate_pr(
        adata, genes, label_col="resp", label="NR", layer=LAYER
    )
    assert out is not None


# --------------------------------------- obs vs categorical (restored in 0.1.1)
def test_obs_categorical_association_scans_all_numeric_columns(adata):
    df = bk.tl.obs_categorical_association(adata, groupby="Subtype")
    assert isinstance(df, pd.DataFrame) and len(df) > 0
    assert {"obs", "groupby", "test", "statistic", "pval", "qval", "effect"} <= set(df.columns)
    # groupby itself must not be tested, and no non-numeric column may appear
    assert "Subtype" not in set(df["obs"])
    for key in df["obs"]:
        assert pd.api.types.is_numeric_dtype(adata.obs[key])


def test_obs_categorical_association_respects_obs_keys(adata):
    df = bk.tl.obs_categorical_association(
        adata, groupby="Subtype", obs_keys=["purity", "age"]
    )
    assert set(df["obs"]) == {"purity", "age"}


def test_obs_categorical_association_detects_planted_difference(adata):
    """A column built to differ by group must come out significant."""
    is_basal = (adata.obs["Subtype"].astype(str) == "Basal").to_numpy()
    adata.obs["planted"] = np.where(is_basal, 10.0, 0.0) + np.random.default_rng(0).normal(
        0, 0.1, adata.n_obs
    )
    adata.obs["noise"] = np.random.default_rng(1).normal(0, 1, adata.n_obs)

    df = bk.tl.obs_categorical_association(
        adata, groupby="Subtype", obs_keys=["planted", "noise"]
    ).set_index("obs")

    assert df.loc["planted", "qval"] < 0.01
    assert df.loc["planted", "effect"] > df.loc["noise", "effect"]
    # the planted column ranks first
    assert df.index[0] == "planted"


def test_obs_categorical_association_group_means(adata):
    df = bk.tl.obs_categorical_association(
        adata, groupby="Subtype", obs_keys=["purity"]
    )
    for lvl in adata.obs["Subtype"].astype(str).unique():
        assert f"mean_{lvl}" in df.columns


def test_obs_categorical_association_anova_alias(adata):
    df = bk.tl.obs_categorical_association(
        adata, groupby="Subtype", obs_keys=["purity"], method="anova"
    )
    assert df.loc[0, "test"] == "anova"


def test_obs_categorical_association_rejects_unknown_keys(adata):
    with pytest.raises(KeyError):
        bk.tl.obs_categorical_association(adata, groupby="Subtype", obs_keys=["nope"])


def test_obs_categorical_association_rejects_unknown_groupby(adata):
    with pytest.raises(KeyError):
        bk.tl.obs_categorical_association(adata, groupby="nope")


def test_rank_genes_groups_fast_accepts_rest_sentinel(adata):
    """`reference="rest"` is what docs and notebooks use; it must equal reference=None."""
    a = bk.tl.rank_genes_groups_fast(
        adata, groupby="Subtype", group="Basal", reference=None, layer=LAYER
    )
    b = bk.tl.rank_genes_groups_fast(
        adata, groupby="Subtype", group="Basal", reference="rest", layer=LAYER
    )
    assert len(a) == len(b)
    assert np.allclose(
        a.set_index("gene")["pval"].to_numpy(),
        b.set_index("gene").loc[a["gene"], "pval"].to_numpy(),
        equal_nan=True,
    )


def test_rank_genes_groups_fast_prefers_a_real_group_named_rest(adata):
    """If a category is genuinely called 'rest', it wins over the sentinel."""
    lab = adata.obs["Subtype"].astype(str).to_numpy().copy()
    lab[lab == "Luminal"] = "rest"
    adata.obs["with_rest"] = lab
    df = bk.tl.rank_genes_groups_fast(
        adata, groupby="with_rest", group="Basal", reference="rest", layer=LAYER
    )
    assert len(df) > 0


def test_store_rank_genes_groups_fast_writes_uns(adata):
    bk.tl.store_rank_genes_groups_fast(
        adata, groupby="Subtype", group="Basal", reference="rest",
        layer=LAYER, subkey="Subtype:Basal_vs_rest",
    )
    assert "rank_genes_groups_fast" in adata.uns
