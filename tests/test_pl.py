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


def test_categorical_palette_unknown_name_falls_back():
    # An unrecognised palette name takes the matplotlib fallback branch, which
    # called mpl.cm.get_cmap -- removed in matplotlib 3.9.
    from bullkpy.pl._colors import categorical_palette

    cols = categorical_palette(["a", "b", "c"], "NotAPalette")
    assert set(cols) == {"a", "b", "c"}


def test_pca_scatter_unknown_palette_falls_back(adata_pca):
    bk.pl.pca_scatter(adata_pca, color="Subtype", palette="NotAPalette", show=False)


def test_no_removed_matplotlib_get_cmap():
    # mpl.cm.get_cmap was removed in matplotlib 3.9; use mpl.colormaps.get_cmap.
    import pathlib
    import re

    root = pathlib.Path(bk.__file__).parent
    bad = [
        str(p.relative_to(root))
        for p in root.rglob("*.py")
        if re.search(r"(?:mpl|matplotlib)\.cm\.get_cmap\(", p.read_text())
    ]
    assert not bad, f"removed API mpl.cm.get_cmap still used in: {bad}"


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


def test_heatmap_de_missing_contrast_lists_available(adata):
    adata.uns["de"] = {"grp_A_vs_B": {"results": None}}
    with pytest.raises(KeyError, match=r"Available contrasts.*grp_A_vs_B"):
        bk.pl.heatmap_de(adata, contrast="A_vs_B", groupby="grp", show=False)


def test_no_removed_boxplot_labels_argument():
    # Axes.boxplot(labels=) was removed in matplotlib 3.11; use tick_labels.
    import pathlib
    import re

    root = pathlib.Path(bk.__file__).parent
    bad = [
        str(p.relative_to(root))
        for p in root.rglob("*.py")
        if re.search(r"\.boxplot\([^\n]*\blabels=", p.read_text().replace("tick_labels=", ""))
    ]
    assert not bad, f"removed boxplot(labels=) still used in: {bad}"


def test_every_figure_is_saved_through_the_shared_helper():
    # One save path means one dpi rule (bk.settings.plot_dpi) and one place to fix bugs.
    import pathlib
    import re

    root = pathlib.Path(bk.__file__).parent
    bad = [
        str(p.relative_to(root))
        for p in root.rglob("*.py")
        if p.name != "_style.py" and re.search(r"\.savefig\(", p.read_text())
    ]
    assert not bad, f"call _savefig from pl._style instead of savefig in: {bad}"


def test_library_code_does_not_print():
    import pathlib
    import re

    root = pathlib.Path(bk.__file__).parent
    bad = [
        str(p.relative_to(root))
        for p in root.rglob("*.py")
        if re.search(r"^\s+print\(", p.read_text(), flags=re.M)
    ]
    assert not bad, f"use bullkpy.logging (info/warn) instead of print in: {bad}"


def test_plot_save_writes_a_file(adata_pca, tmp_path):
    out = tmp_path / "sub" / "pca.png"
    bk.pl.pca_scatter(adata_pca, color="Subtype", save=out, show=False)
    assert out.exists() and out.stat().st_size > 0


def _fresh_ax():
    fig, ax = plt.subplots()
    return fig, ax


def test_single_axes_plots_draw_on_a_given_axes(adata, de_table, genes):
    fig, axes = plt.subplots(1, 4, figsize=(12, 3))
    bk.pl.volcano(de_table, ax=axes[0], show=False)
    bk.pl.ma(de_table, ax=axes[1], show=False)
    bk.pl.gene_plot(adata, gene=genes[0], groupby="Subtype", ax=axes[2], show=False)
    bk.pl.library_size_vs_genes(adata, ax=axes[3], show=False)
    for ax in axes:
        assert ax.has_data(), "nothing was drawn on the axes that was passed in"
    assert len(plt.get_fignums()) == 1, "a second figure was created"
    plt.close("all")


def test_ax_argument_is_added_to_the_single_axes_plots():
    import inspect

    for name in ("volcano", "ma", "rankplot", "gene_plot", "pca_variance_ratio", "boxplot_with_stats",
                 "km_univariate", "cox_forest", "gsea_bubbleplot"):
        assert "ax" in inspect.signature(getattr(bk.pl, name)).parameters, name


def test_clustered_heatmaps_refuse_an_external_axes(adata, genes):
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="ax="):
        bk.pl.gene_panel_correlation_heatmap(adata, genes=genes[:4], cluster=True, ax=ax, show=False)
    with pytest.raises(ValueError, match="ax="):
        bk.pl.corr_heatmap(adata, use="samples", dendrogram=True, ax=ax, show=False)
    plt.close("all")


def _violin_colour_by_position(ax):
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import to_hex

    found = {}
    for coll in ax.collections:
        if isinstance(coll, PolyCollection):
            for path in coll.get_paths():
                if len(path.vertices) > 10:
                    found[round(float(path.vertices[:, 0].mean()))] = to_hex(coll.get_facecolor()[0])
    return [found[i] for i in sorted(found)]


def test_violin_colours_follow_the_drawing_order(adata):
    # Palette colours go to categories in the order they are drawn (first drawn = first colour),
    # for default order, an explicit order and a categorical column with its own category order.
    import seaborn as sns
    from matplotlib.colors import to_hex

    adata.obs["Subtype"] = adata.obs["Subtype"].astype("category").cat.reorder_categories(["Luminal", "Basal"])
    first, second = (to_hex(sns.desaturate(c, 0.75)) for c in sns.color_palette("Set2", 2))
    fig, axes = bk.pl.violin(adata, keys=["purity"], groupby="Subtype", show=False)
    assert _violin_colour_by_position(axes[0]) == [first, second]
    fig, axes = bk.pl.violin(adata, keys=["purity"], groupby="Subtype", order=["Basal", "Luminal"], show=False)
    assert _violin_colour_by_position(axes[0]) == [first, second]
    plt.close("all")


def test_seaborn_plots_do_not_use_the_deprecated_palette_without_hue(adata, genes):
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        bk.pl.violin(adata, keys=["purity"], groupby="Subtype", show=False)
        bk.pl.gene_association(adata, gene=[genes[0]], groupby="Subtype", layer="log1p_cpm", show=False)
        bk.pl.gene_association(adata, gene=[genes[0]], groupby="Subtype", layer="log1p_cpm", kind="box", show=False)
    plt.close("all")


# ----------------------------------------------------------------- dotplot layout
def _dot_adata(n_groups, n_genes, n=240, seed=0):
    import anndata as ad
    import numpy as np
    import pandas as pd

    r = np.random.default_rng(seed)
    X = r.gamma(2.0, 1.0, size=(n, n_genes)).astype("float32")
    a = ad.AnnData(X, var=pd.DataFrame(index=[f"GENE{i + 1}" for i in range(n_genes)]))
    a.obs["grp"] = pd.Categorical([f"A fairly long group name {i}" for i in r.integers(0, n_groups, n)])
    a.layers["log1p_cpm"] = np.log1p(X)
    a.layers["counts"] = X
    return a


def _inside_figure(fig, tol=0.03):
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    w, h = fig.get_size_inches()
    return bb.x0 >= -tol and bb.y0 >= -tol and bb.x1 <= w + tol and bb.y1 <= h + tol


@pytest.mark.parametrize(
    "n_groups,n_genes,kwargs",
    [
        (2, 3, {}),
        (3, 8, {"dendrogram_top": True, "dendrogram_rows": True}),
        (10, 6, {"dendrogram_top": True, "dendrogram_rows": True}),
        (12, 25, {}),
        (25, 12, {"swap_axes": True, "dendrogram_top": True, "dendrogram_rows": True}),
        (6, 6, {"dendrogram_rows": True, "row_dendrogram_position": "outer_left"}),
        (6, 6, {"dendrogram_rows": True, "row_dendrogram_position": "left"}),
        (6, 6, {"title": "A title", "dendrogram_top": True}),
    ],
)
def test_dotplot_layout_fits_the_figure_for_any_size(n_groups, n_genes, kwargs):
    a = _dot_adata(n_groups, n_genes)
    fig, ax = bk.pl.dotplot(a, var_names=list(a.var_names), groupby="grp", layer="log1p_cpm", show=False, **kwargs)
    assert _inside_figure(fig), "labels, dendrograms or legends are cut off or outside the figure"
    # legends sit to the right of the panel and dendrograms, never on top of them
    panel = ax.get_position()
    others = [x for x in fig.axes if x is not ax and x.get_visible()]
    for other in others:
        pos = other.get_position()
        overlap_x = min(panel.x1, pos.x1) - max(panel.x0, pos.x0)
        overlap_y = min(panel.y1, pos.y1) - max(panel.y0, pos.y0)
        is_overlay = pos.width > 0.99 and pos.height > 0.99  # the transparent full-figure legend canvas
        assert is_overlay or overlap_x <= 1e-6 or overlap_y <= 1e-6, "an axes overlaps the dot panel"
    plt.close("all")


def test_dotplot_size_follows_the_number_of_groups_and_genes():
    sizes = {}
    for ng, nf in [(2, 3), (6, 10), (20, 30)]:
        a = _dot_adata(ng, nf)
        fig, _ = bk.pl.dotplot(a, var_names=list(a.var_names), groupby="grp", layer="log1p_cpm", show=False)
        sizes[(ng, nf)] = tuple(fig.get_size_inches())
        plt.close("all")
    (w1, h1), (w2, h2), (w3, h3) = sizes[(2, 3)], sizes[(6, 10)], sizes[(20, 30)]
    assert w1 < w2 < w3 and h1 < h2 < h3
    assert h1 < 3.5, "a 2 x 3 dotplot should not get the old 4.2 in minimum height"


def test_dotplot_cell_size_and_fixed_figsize():
    a = _dot_adata(5, 6)
    kw = dict(var_names=list(a.var_names), groupby="grp", layer="log1p_cpm", show=False)
    fig_small, _ = bk.pl.dotplot(a, cell_size=0.25, **kw)
    fig_big, _ = bk.pl.dotplot(a, cell_size=0.5, **kw)
    assert fig_big.get_size_inches()[0] > fig_small.get_size_inches()[0]
    fig_fixed, _ = bk.pl.dotplot(a, figsize=(7, 4), **kw)
    assert tuple(fig_fixed.get_size_inches()) == (7.0, 4.0)
    assert _inside_figure(fig_fixed)
    plt.close("all")


def test_dotplot_dendrogram_leaves_sit_on_the_dots():
    a = _dot_adata(8, 7)
    fig, ax = bk.pl.dotplot(
        a, var_names=list(a.var_names), groupby="grp", layer="log1p_cpm",
        dendrogram_top=True, dendrogram_rows=True, show=False,
    )
    tops = [x for x in fig.axes if x is not ax and x.get_xlim() == ax.get_xlim() and x.lines]
    rows = [x for x in fig.axes if x is not ax and x.get_ylim() == ax.get_ylim() and x.lines]
    assert tops and rows, "dendrogram axes must share the panel's limits so leaves align with dots"
    plt.close("all")


def test_dotplot_size_legend_is_labelled_per_dot_and_follows_dot_max():
    a = _dot_adata(4, 5)
    kw = dict(var_names=list(a.var_names), groupby="grp", layer="log1p_cpm", show=False)
    fig, _ = bk.pl.dotplot(a, **kw)
    texts = {t.get_text() for x in fig.axes for t in x.texts}
    assert {"20", "40", "60", "80", "100"} <= texts
    plt.close("all")
    fig, _ = bk.pl.dotplot(a, dot_max=0.6, **kw)
    texts = {t.get_text() for x in fig.axes for t in x.texts}
    assert "100" not in texts and "60" in texts, "legend must not advertise sizes above dot_max"
    plt.close("all")


def test_dotplot_small_figsize_grows_and_titles_can_be_omitted(caplog):
    a = _dot_adata(3, 8)
    fig, ax = bk.pl.dotplot(
        a, var_names=list(a.var_names), groupby="grp", layer="log1p_cpm",
        figsize=(4, 2), colorbar_title=None, show=False,
    )
    w, h = fig.get_size_inches()
    assert (w, h) != (4.0, 2.0) and w > 4.0, "a too-small figsize must be enlarged, not squashed"
    assert _inside_figure(fig)
    assert not any("Mean expression" in t.get_text() for t in fig.texts)
    assert not ax.xaxis.get_gridlines()[0].get_visible()
    plt.close("all")
