from __future__ import annotations

from typing import Mapping, Sequence
import numpy as np
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def score_genes_dict(
    adata: ad.AnnData,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    layer: str | None = "log1p_cpm",
    prefix: str = "score_",
    scale: bool = False,
    store_params: bool = True,
) -> None:
    """
    Compute gene set (signature) scores and store them in `adata.obs`.

    Parameters
    ----------
    adata
        AnnData object
    gene_sets
        Dict mapping {signature_name: [gene1, gene2, ...]}
    layer
        Layer to use (None -> adata.X)
    prefix
        Prefix for columns added to adata.obs
    scale
        If True, z-score each signature across samples
    store_params
        If True, store minimal metadata in adata.uns["score_genes_dict"]

    Notes
    -----
    - Scores are mean expression of genes per sample
    - Missing genes are ignored (with warning)
    - This follows scanpy.tl.score_genes philosophy
    """
    if not isinstance(gene_sets, Mapping):
        raise TypeError("gene_sets must be a dict-like mapping")

    X = adata.layers[layer] if layer is not None else adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    var_names = np.asarray(adata.var_names, dtype=str)

    stored = {}

    for name, genes in gene_sets.items():
        genes = [str(g) for g in genes]
        idx = np.isin(var_names, genes)

        if idx.sum() == 0:
            warn(f"Gene set '{name}' has no genes present in adata.var_names; skipped.")
            continue

        score = X[:, idx].mean(axis=1)

        if scale:
            sd = score.std(ddof=1)
            if sd > 0:
                score = (score - score.mean()) / sd

        col = f"{prefix}{name}"
        adata.obs[col] = score

        stored[name] = {
            "n_genes_used": int(idx.sum()),
            "genes_present": var_names[idx].tolist(),
        }

        info(f"Stored gene set score '{col}' in adata.obs ({idx.sum()} genes)")

    if store_params:
        adata.uns["score_genes_dict"] = {
            "layer": layer,
            "prefix": prefix,
            "scale": scale,
            "gene_sets": stored,
        }