from __future__ import annotations

from .pca import (
    pca,
    pca_loadings,
)
from .neighbors import (
    neighbors,
)
from .clustering import (
    cluster,
    cluster_metrics,
    categorical_confusion,
    leiden_resolution_scan,
)
from .umap import (
    umap,
)
from .de import (
    de,
)
from .de_glm import (
    de_glm,
)
from .rank_genes_groups import (
    rank_genes_groups,
)
from .adjusted_rand_index import (
    adjusted_rand_index,
)
from .score_genes import (
    score_genes,
    score_genes_dict,
    score_genes_cell_cycle,
)
from .tcga_groups import (
    tcga_define_groups,
)
from .dispersion import (
    sample_dispersion,
    signature_dispersion_by_group,
)
from .metaprograms import (
    score_metaprograms,
    metaprogram_dispersion_metrics,
    metaprogram_heterogeneity,
    summarize_by_group,
    metaprogram_sample_metrics,
    metaprogram_dispersion_by_group,
    metaprogram_ne_enrichment,
    metaprogram_topk_contribution,
)
from .correlations import (
    top_gene_gene_correlations,
    gene_gene_correlations,
    top_gene_obs_correlations,
    top_obs_obs_correlations,
    obs_obs_corr_matrix,
    partial_corr,
)
from .associations import (
    association,
    gene_categorical_association,
    obs_categorical_association,
    rank_genes_groups_fast,
    categorical_association,
    posthoc_per_gene,
    gene_metadata_association_scan,
    signature_axis,
    quantile_groups,
)
from .posthoc import (
    pairwise_posthoc,
)
from .cox import (
    cox_gene_association,
    cox_multivariate,
    cox_penalizer_cv,
    cox_stability_selection,
    panel_size_cindex,
    recommended_cox_panel,
    cox_fit_penalized,
    cox_univariate,
    cox_interaction,
    surv_1d_bins,
    surv_2x2_bins,
    run_cox_per_group,
)
from .signature import (
    filter_genes_var,
    rank_genes_univariate_pr,
    gene_corr_redundancy,
    signature_score,
    validate_signature_external,
    validate_signature_auc_external,
    design_signatures,
    signature_permutation_importance,
    correlation_redundancy_filter,
    redundancy_filter_by_corr,
    meta_rank_genes_from_signatures,
    recommended_auc_panel,
    make_signature_frozen,
    quick_fit_weights,
    add_signature_from_genes,
    benchmark_signatures_pr_auc,
    fit_intercept_only,
    evaluate_signatures_store,
    pick_best_signature,
    score_signature_on_ppv,
    build_signature_bank_from_genes,
    bootstrap_pr_bands,
)
from .gsea_prerank import (
    gsea_preranked,
    list_enrichr_libraries,
)
from .gsea_leading_edge import (
    leading_edge_pathway_clusters,
    leading_edge_cluster_driver_genes,
    export_leading_edge_clusters_cytoscape,
)

__all__ = [
    "pca",
    "neighbors",
    "cluster",
    "umap",
    "de",
    "de_glm",
    "pca_loadings",
    "rank_genes_groups",
    "adjusted_rand_index",
    "cluster_metrics",
    "leiden_resolution_scan",
    "categorical_confusion",
    "score_genes",
    "score_genes_dict",
    "score_genes_cell_cycle",
    "top_gene_gene_correlations",
    "top_gene_obs_correlations",
    "top_obs_obs_correlations",
    "gene_gene_correlations",
    "obs_obs_corr_matrix",
    "partial_corr",
    "association",
    "gene_categorical_association",
    "obs_categorical_association",
    "rank_genes_groups_fast",
    "categorical_association",
    "posthoc_per_gene",
    "pairwise_posthoc",
    "cox_gene_association",
    "cox_multivariate",
    "cox_penalizer_cv",
    "cox_stability_selection",
    "panel_size_cindex",
    "recommended_cox_panel",
    "cox_fit_penalized",
    "cox_univariate",
    "cox_interaction",
    "surv_1d_bins",
    "surv_2x2_bins",
    "filter_genes_var",
    "rank_genes_univariate_pr",
    "gene_corr_redundancy",
    "signature_score",
    "validate_signature_external",
    "validate_signature_auc_external",
    "design_signatures",
    "signature_permutation_importance",
    "correlation_redundancy_filter",
    "redundancy_filter_by_corr",
    "meta_rank_genes_from_signatures",
    "recommended_auc_panel",
    "make_signature_frozen",
    "quick_fit_weights",
    "fit_intercept_only",
    "benchmark_signatures_pr_auc",
    "add_signature_from_genes",
    "evaluate_signatures_store",
    "pick_best_signature",
    "score_signature_on_ppv",
    "build_signature_bank_from_genes",
    "gsea_preranked",
    "list_enrichr_libraries",
    "tcga_define_groups",
    "sample_dispersion",
    "signature_dispersion_by_group",
    "score_metaprograms",
    "metaprogram_dispersion_metrics",
    "metaprogram_heterogeneity",
    "summarize_by_group",
    "signature_axis",
    "quantile_groups",
    "gene_metadata_association_scan",
    "metaprogram_sample_metrics",
    "metaprogram_dispersion_by_group",
    "metaprogram_ne_enrichment",
    "metaprogram_topk_contribution",
    "bootstrap_pr_bands",
    "run_cox_per_group",
    "leading_edge_pathway_clusters",
    "leading_edge_cluster_driver_genes",
    "export_leading_edge_clusters_cytoscape",
]

from .._deprecated import moved_functions

__getattr__ = moved_functions(
    "bullkpy.tl",
    {
        "umap_graph": "bullkpy.tl._legacy.umap_graph",
        "store_gene_categorical_association": "bullkpy.tl._legacy.store_gene_categorical_association",
        "store_rank_genes_groups_fast": "bullkpy.tl._legacy.store_rank_genes_groups_fast",
        "plot_corr_scatter": "bullkpy.pl.plot_corr_scatter",
        "plot_corr_heatmap": "bullkpy.pl.plot_corr_heatmap",
        "mp_dispersion_metrics": "bullkpy.tl.metaprogram_dispersion_metrics",
        "read_metaprograms_xlsx": "bullkpy.io.read_metaprograms_xlsx",
        "metaprogram_scores_get": "bullkpy.get.metaprogram_scores",
        "score_signature_on_PPV": "bullkpy.tl.score_signature_on_ppv",
        "save_signature_bank": "bullkpy.io.save_signature_bank",
        "load_signature_bank": "bullkpy.io.load_signature_bank",
    },
)
