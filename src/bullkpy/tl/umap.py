from __future__ import annotations

import numpy as np
import anndata as ad

from ..logging import info


def umap(
    adata: ad.AnnData,
    *,
    n_neighbors: int = 15,
    n_pcs: int = 20,
    use_rep: str = "X_pca",
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    metric: str = "euclidean",
    random_state: int = 0,
    init: str = "spectral",
) -> None:
    """
    Compute UMAP embedding from a representation (default: PCA).

    This mirrors Scanpy practice: UMAP is computed from `X_pca` with
    `n_neighbors`/`min_dist`, consistent with the neighbors graph settings.

    Stores:
      - adata.obsm['X_umap']
      - adata.uns['umap']
    """
    if use_rep not in adata.obsm:
        raise KeyError(f"adata.obsm['{use_rep}'] not found. Run bk.tl.pca(adata) first.")

    try:
        import umap  # umap-learn
    except Exception as e:
        raise ImportError(
            "UMAP requires `umap-learn`. Install with: pip install umap-learn\n"
            f"Original error: {e}"
        )

    X = np.asarray(adata.obsm[use_rep], dtype=float)
    if n_pcs is not None:
        n_pcs = int(min(n_pcs, X.shape[1]))
        X = X[:, :n_pcs]

    info(
        "Running UMAP "
        f"(rep={use_rep}, n_pcs={n_pcs}, n_neighbors={n_neighbors}, min_dist={min_dist}, "
        f"n_components={n_components})"
    )

    reducer = umap.UMAP(
        n_neighbors=int(n_neighbors),
        n_components=int(n_components),
        min_dist=float(min_dist),
        spread=float(spread),
        metric=metric,
        random_state=int(random_state),
        init=init,
    )

    emb = reducer.fit_transform(X)

    adata.obsm["X_umap"] = np.asarray(emb, dtype=float)
    adata.uns["umap"] = {
        "params": {
            "use_rep": use_rep,
            "n_neighbors": int(n_neighbors),
            "n_pcs": None if n_pcs is None else int(n_pcs),
            "min_dist": float(min_dist),
            "spread": float(spread),
            "n_components": int(n_components),
            "metric": metric,
            "random_state": int(random_state),
            "init": init,
        }
    }