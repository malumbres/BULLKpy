from __future__ import annotations

import numpy as np
import scipy.sparse as sp
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


def umap_graph(
    adata: ad.AnnData,
    *,
    graph_key: str = "connectivities",
    use_rep: str = "X_pca",
    n_pcs: int = 20,
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
    random_state: int = 0,
    init: str = "spectral",  # "spectral" or "random"
    negative_sample_rate: int = 5,
    n_epochs: int | None = None,
) -> None:
    """
    Compute UMAP embedding strictly from a precomputed neighbor graph.

    Requires:
      - adata.obsp[graph_key]  (typically 'connectivities' from bk.tl.neighbors)

    Uses:
      - adata.obsm[use_rep] only for initialization (spectral/random), not for graph construction.

    Stores:
      - adata.obsm['X_umap_graph']
      - adata.uns['umap_graph']
    """
    if graph_key not in adata.obsp:
        raise KeyError(f"adata.obsp['{graph_key}'] not found. Run bk.tl.neighbors(adata) first.")

    if use_rep not in adata.obsm:
        raise KeyError(f"adata.obsm['{use_rep}'] not found. Run bk.tl.pca(adata) first.")

    try:
        import umap  # umap-learn
    except Exception as e:
        raise ImportError(
            "UMAP requires `umap-learn`. Install with: pip install umap-learn\n"
            f"Original error: {e}"
        )

    # Graph (fixed)
    graph = adata.obsp[graph_key]
    if not sp.issparse(graph):
        graph = sp.csr_matrix(graph)
    graph = graph.tocsr().astype(np.float32)
    graph.eliminate_zeros()

    # Representation only for init
    X = np.asarray(adata.obsm[use_rep], dtype=np.float32)
    if n_pcs is not None:
        n_pcs = int(min(n_pcs, X.shape[1]))
        X = X[:, :n_pcs]

    info(
        "Running graph-based UMAP "
        f"(graph={graph_key}, init={init}, n_pcs={n_pcs}, min_dist={min_dist}, n_components={n_components})"
    )

    # Create a reducer mainly to get (a, b) parameters for min_dist/spread
    reducer = umap.UMAP(
        n_neighbors=15,  # required by umap-learn validator; graph is overridden below
        n_components=int(n_components),
        min_dist=float(min_dist),
        spread=float(spread),
        metric="euclidean",
        random_state=int(random_state),
        init=init,
        negative_sample_rate=int(negative_sample_rate),
        n_epochs=n_epochs,
        verbose=False,
    )

    # Fit once to initialize internal parameters (a, b, etc.). We will override graph_.

    # Fit once to initialize internal parameters (a, b, etc.) and get a valid init embedding
    reducer.fit(X)

    # Use the embedding from the initial fit as a valid initializer (works across versions)
    init_emb = getattr(reducer, "embedding_", None)
    if init_emb is None:
        raise RuntimeError("UMAP did not produce an initial embedding during fit().")

    # Override the graph with *your* connectivities and re-embed.
    reducer.graph_ = graph

    # Call private embedding routine with a version-compatible signature
    try:
        # older versions
        emb = reducer._fit_embed_data(X)
    except TypeError:
        # newer versions require (X, n_epochs, init, random_state)

        n_epochs_eff = n_epochs if n_epochs is not None else getattr(reducer, "n_epochs", None)
        if n_epochs_eff is None:
            n_epochs_eff = 500  # safe default for small cohorts


        rs = np.random.RandomState(int(random_state))
        emb = reducer._fit_embed_data(X, int(n_epochs_eff), init_emb, rs)

    if isinstance(emb, tuple):
        emb = emb[0]

    # umap-learn may return (embedding, aux). Keep only the embedding.
    if isinstance(emb, (tuple, list)) and len(emb) >= 1:
        emb = emb[0]
    emb = np.asarray(emb, dtype=float)
    if emb.ndim != 2 or emb.shape[0] != adata.n_obs:
        raise ValueError(f"UMAP graph embedding has unexpected shape: {emb.shape}")
    adata.obsm["X_umap_graph"] = emb

    adata.uns["umap_graph"] = {
        "params": {
            "mode": "graph",
            "graph_key": graph_key,
            "use_rep_init": use_rep,
            "n_pcs_init": None if n_pcs is None else int(n_pcs),
            "min_dist": float(min_dist),
            "spread": float(spread),
            "n_components": int(n_components),
            "random_state": int(random_state),
            "init": init,
            "negative_sample_rate": int(negative_sample_rate),
            "n_epochs": n_epochs,
        }
    }