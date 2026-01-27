from ._style import set_style, _savefig, get_palette
 
from .qc import (
    qc_metrics,
    library_size_vs_genes,
    mt_fraction_vs_counts,
    genes_vs_mt_fraction,
    qc_scatter_panel,
    qc_by_group, qc_pairplot,
)

from .pca import pca_scatter, pca_variance_ratio
from .umap import umap
from .de import volcano, rankplot, ma
from .dotplot import dotplot
from .heatmap_de import heatmap_de
from .violin import violin
from .corr_heatmap import corr_heatmap, gene_panel_correlation_heatmap
from .pca_loadings import pca_loadings_bar, pca_loadings_heatmap
from .gene_plot import gene_plot
from .rank_genes_groups import rank_genes_groups, rank_genes_groups_dotplot
from .sample_distances import sample_distances, sample_correlation_clustergram
from ._colors import _get_series, get_categorical_colors, categorical_colors_array
from .clustering_plots import ari_resolution_heatmap, categorical_confusion
from .corrplot_obs import corrplot_obs
from .association import association_heatmap, boxplot_with_stats, categorical_confusion
from .association_rankplots import (
    rankplot_association, dotplot_association, heatmap_association,
    )

from .dispersion import (
    axis_vs_dispersion, dispersion_summary,
    metaprogram_corr, metaprogram_heatmap,
    metaprogram_dispersion_heatmap, metaprogram_metrics_summary,
    metaprogram_ne_scatter, metaprogram_dominance_ridgeplot_like,
    )

from .cox import (
    cox_volcano, cox_forest, km_plot_signature, time_dependent_roc,
    panel_size_cindex_plot, stability_barplot, permutation_importance_barplot,
    stability_freq_bar, 
    cox_forest_from_uns, run_cox_per_group,
    km_univariate, km_2x2_interaction,
    metaprogram_rank1_composition_stackedbar,
    )

from .signature_plots import (
    auc_cv_heatmap, pr_curve_signatures, plot_pr_curves_signatures,
    bootstrap_pr_bands, plot_pr_curves_signatures_bootstrap,
    plot_pub_template_pr, signature_dispersion_heatmap,
    )

from .gene_association import gene_association, gene_association_volcano
from .volcano_categorical import volcano_categorical
from .gsea_leading_edge import (
     gsea_leading_edge_heatmap, leading_edge_jaccard_heatmap, leading_edge_overlap_matrix,
     leading_edge_pathway_clusters, leading_edge_cluster_driver_genes,  
     leading_edge_cluster_bubbles, export_leading_edge_clusters_cytoscape,
     )
from .oncoprint import oncoprint

from .gsea_bubbleplot import gsea_bubbleplot

__all__ = ["qc_metrics", "library_size_vs_genes",
           "mt_fraction_vs_counts", "genes_vs_mt_fraction",
           "qc_scatter_panel",
           "set_style", "_savefig", "get_palette",
           "qc_pairplot", "qc_by_group", "pca_scatter", "umap",
           "volcano", "rankplot", "ma",
           "dotplot",
           "heatmap_de",
           "violin",
           "pca_variance_ratio",
           "corr_heatmap", "gene_panel_correlation_heatmap",
           "pca_loadings_bar", "pca_loadings_heatmap",
           "gene_plot",
           "rank_genes_groups", "rank_genes_groups_dotplot",
           "library_size_vs_genes",
           "sample_distances", "sample_correlation_clustergram",
           "_get_series", "get_categorical_colors", "categorical_colors_array",
           "ari_resolution_heatmap", "categorical_confusion",
           "corrplot_obs",
           "association_heatmap", "boxplot_with_stats", "categorical_confusion",
           "rankplot_association", "dotplot_association", "heatmap_association",
           "gene_association", "volcano_categorical", "gene_association_volcano",
           "cox_volcano", "cox_forest", "km_plot_signature",
           "cox_forest_from_uns", "run_cox_per_group",
           "metaprogram_rank1_composition_stackedbar",
           "km_univariate", "km_2x2_interaction",
           "permutation_importance_barplot",
           "stability_freq_bar", "time_dependent_roc",
           "panel_size_cindex_plot", "stability_barplot",
           "auc_cv_heatmap", 
           "pr_curve_signatures", "plot_pr_curves_signatures",
           "bootstrap_pr_bands", "plot_pr_curves_signatures_bootstrap",
           "plot_pub_template_pr", 
           "gsea_leading_edge_heatmap", "leading_edge_jaccard_heatmap",
           "leading_edge_overlap_matrix",
           "leading_edge_pathway_clusters", "leading_edge_cluster_driver_genes", 
           "leading_edge_cluster_bubbles", "export_leading_edge_clusters_cytoscape",
           "oncoprint", "gsea_bubbleplot",
           "axis_vs_dispersion", "dispersion_summary",
           "metaprogram_corr", "metaprogram_heatmap",
           "signature_dispersion_heatmap",
           "metaprogram_dispersion_heatmap", "metaprogram_metrics_summary",
           "metaprogram_ne_scatter", "metaprogram_dominance_ridgeplot_like",
]
