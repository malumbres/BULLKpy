"""pl: every plotting entry point must build a figure without raising."""
from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import pytest

import bullkpy as bk

matplotlib.use("Agg")
LAYER = "log1p_cpm"


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def de_table(adata):
    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal")
    return adata.uns["de"]["Subtype_Basal_vs_Luminal"]["results"]


# ------------------------------------------------------------------ QC plots
def test_qc_metrics_panel(adata):
    fig, axes = bk.pl.qc_metrics(adata, show=False)
    assert fig is not None and len(axes) >= 2


def test_library_size_vs_genes(adata):
    bk.pl.library_size_vs_genes(adata, show=False)


def test_mt_fraction_vs_counts(adata):
    bk.pl.mt_fraction_vs_counts(adata, show=False)


def test_genes_vs_mt_fraction(adata):
    bk.pl.genes_vs_mt_fraction(adata, show=False)


def test_mt_fraction_vs_counts_with_thresholds(adata):
    """Passing thresholds used to raise NameError while building the title."""
    fig, ax = bk.pl.mt_fraction_vs_counts(adata, min_counts=1.0, max_mt=99.0, show=False)
    assert "QC fail" in ax.get_title()


def test_genes_vs_mt_fraction_with_thresholds(adata):
    fig, ax = bk.pl.genes_vs_mt_fraction(adata, min_mt=0.0, max_genes=1e9, show=False)
    assert "QC fail" in ax.get_title()


def test_qc_scatter_panel(adata):
    bk.pl.qc_scatter_panel(adata, show=False)


def test_qc_by_group(adata):
    bk.pl.qc_by_group(adata, groupby="Subtype", show=False)


def test_qc_by_group_skips_missing_metrics(adata):
    """Missing QC columns should be skipped, not raise."""
    del adata.obs["pct_counts_ribo"]
    bk.pl.qc_by_group(adata, groupby="Subtype", show=False)


def test_qc_pairplot(adata):
    bk.pl.qc_pairplot(adata, show=False)


# ----------------------------------------------------------------- embeddings
def test_pca_scatter_categorical_color(adata_pca):
    bk.pl.pca_scatter(adata_pca, color="Subtype", show=False)


def test_pca_scatter_numeric_color(adata_pca):
    bk.pl.pca_scatter(adata_pca, color="purity", show=False)


def test_pca_variance_ratio(adata_pca):
    bk.pl.pca_variance_ratio(adata_pca, show=False)


def test_pca_loadings_bar(adata_pca):
    bk.pl.pca_loadings_bar(adata_pca, show=False)


def test_pca_loadings_heatmap(adata_pca):
    bk.pl.pca_loadings_heatmap(adata_pca, show=False)


def test_umap_plot(adata_pca):
    pytest.importorskip("umap")
    bk.tl.umap(adata_pca)
    bk.pl.umap(adata_pca, color="Subtype", show=False)


# ------------------------------------------------------------- expression
def test_violin(adata, genes):
    bk.pl.violin(adata, keys=genes[:3], groupby="Subtype", layer=LAYER, show=False)


def test_dotplot(adata, genes):
    bk.pl.dotplot(adata, var_names=genes[:5], groupby="Subtype", layer=LAYER, show=False)


def test_gene_plot(adata, genes):
    bk.pl.gene_plot(adata, gene=genes[0], groupby="Subtype", layer=LAYER, show=False)


def test_corrplot_gene_vs_gene(adata, genes):
    bk.pl.corrplot(
        adata, x=genes[0], y=genes[1], x_source="gene", y_source="gene",
        layer=LAYER, show=False,
    )


def test_corr_heatmap(adata):
    bk.pl.corr_heatmap(adata, layer=LAYER, use="samples", show=False)


def test_gene_panel_correlation_heatmap_clustered(adata, genes):
    bk.pl.gene_panel_correlation_heatmap(adata, genes=genes, layer=LAYER, cluster=True, show=False)


def test_gene_panel_correlation_heatmap_unclustered(adata, genes):
    """cluster=False previously raised UnboundLocalError on `plt`."""
    bk.pl.gene_panel_correlation_heatmap(adata, genes=genes, layer=LAYER, cluster=False, show=False)


# --------------------------------------------------------------- distances
def test_sample_distances(adata):
    bk.pl.sample_distances(adata, layer=LAYER, show=False)


def test_sample_distances_with_metadata_colors(adata):
    """The legend branch previously raised UnboundLocalError on `plt`."""
    bk.pl.sample_distances(adata, layer=LAYER, col_colors=["Subtype"], show=False)


def test_sample_correlation_clustergram(adata):
    bk.pl.sample_correlation_clustergram(adata, layer=LAYER, show=False)


# ---------------------------------------------------------------- DE plots
def test_volcano(de_table):
    bk.pl.volcano(de_table, show=False)


def test_rankplot_from_results_frame(de_table):
    bk.pl.rankplot(res=de_table, show=False)


def test_rankplot_from_adata(adata, de_table):
    bk.pl.rankplot(adata, contrast="Subtype_Basal_vs_Luminal", show=False)


def test_ma(de_table):
    bk.pl.ma(result=de_table, mean_col="mean_group", show=False)


def test_heatmap_de(adata, de_table):
    bk.pl.heatmap_de(
        adata, contrast="Subtype_Basal_vs_Luminal", groupby="Subtype",
        layer=LAYER, top_n=20, show=False,
    )


# ------------------------------------------------------------- categorical
def test_categorical_confusion_plot(adata):
    bk.pl.categorical_confusion(adata, key1="Subtype", key2="Batch", show=False)


# -------------------------------------------------------------- colour API
def test_get_palette_returns_n_colors():
    pal = bk.pl.get_palette(5)
    assert len(pal) == 5


def test_get_categorical_colors_covers_all_levels(adata):
    cmap = bk.pl.get_categorical_colors(adata, key="Subtype", where="obs")
    assert set(cmap) == set(adata.obs["Subtype"].astype(str).unique())


def test_categorical_colors_array_matches_n_obs(adata):
    arr = bk.pl.categorical_colors_array(adata, key="Subtype", where="obs")
    assert len(arr) == adata.n_obs


def test_set_style_is_idempotent():
    bk.pl.set_style()
    bk.pl.set_style()
