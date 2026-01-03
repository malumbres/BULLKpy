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
from .corr_heatmap import corr_heatmap
from .pca_loadings import pca_loadings_bar, pca_loadings_heatmap
from .gene_plot import gene_plot
from .rank_genes_groups import rank_genes_groups, rank_genes_groups_dotplot
from .sample_distances import sample_distances, sample_correlation_clustergram
from ._colors import _get_series, get_categorical_colors, categorical_colors_array
from .ari_resolution_heatmap import ari_resolution_heatmap
from .corrplot_obs import corrplot_obs
from .association import association_heatmap, boxplot_with_stats, categorical_confusion
from .association_rankplots import (
    rankplot_association, dotplot_association, heatmap_association,
    )
from .gene_association import gene_association, gene_association_volcano
from .volcano_categorical import volcano_categorical
from .gsea_leading_edge import (gsea_leading_edge_heatmap, leading_edge_jaccard_heatmap,
    leading_edge_overlap_matrix)
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
           "corr_heatmap",
           "pca_loadings_bar", "pca_loadings_heatmap",
           "gene_plot",
           "rank_genes_groups", "rank_genes_groups_dotplot",
           "library_size_vs_genes",
           "sample_distances", "sample_correlation_clustergram",
           "_get_series", "get_categorical_colors", "categorical_colors_array",
           "ari_resolution_heatmap", "corrplot_obs",
           "association_heatmap", "boxplot_with_stats", "categorical_confusion",
           "rankplot_association", "dotplot_association", "heatmap_association",
           "gene_association", "volcano_categorical", "gene_association_volcano",
           "gsea_leading_edge_heatmap", "leading_edge_jaccard_heatmap",
           "leading_edge_overlap_matrix",
           "oncoprint", "gsea_bubbleplot",
]
