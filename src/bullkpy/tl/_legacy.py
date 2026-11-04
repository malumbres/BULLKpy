"""Old entry points kept for one deprecation cycle. Reached through ``bk.tl.<name>``."""
from __future__ import annotations

from .associations import gene_categorical_association, rank_genes_groups_fast
from .umap import umap


def umap_graph(adata, **kwargs):
    """Deprecated: use ``bk.tl.umap(adata, method="graph")``."""
    return umap(adata, method="graph", **kwargs)


def store_gene_categorical_association(adata, *, groupby, key="associations", subkey=None, **kwargs):
    """Deprecated: use ``bk.tl.gene_categorical_association(..., key_added=...)``."""
    return gene_categorical_association(adata, groupby=groupby, key_added=key, subkey=subkey, **kwargs)


def store_rank_genes_groups_fast(
    adata, *, groupby, group, reference=None, key="rank_genes_groups_fast", subkey=None, **kwargs
):
    """Deprecated: use ``bk.tl.rank_genes_groups_fast(..., key_added=...)``."""
    return rank_genes_groups_fast(
        adata, groupby=groupby, group=group, reference=reference, key_added=key, subkey=subkey, **kwargs
    )
