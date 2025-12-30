from __future__ import annotations

import numpy as np
import anndata as ad

from ..logging import info, warn


def cluster(
    adata: ad.AnnData,
    *,
    method: str = "leiden",     # "leiden" or "kmeans"
    key_added: str = "clusters",
    resolution: float = 1.0,    # used for leiden
    n_clusters: int = 8,        # used for kmeans fallback
    use_rep: str = "X_pca",
    n_pcs: int = 20,
    random_state: int = 0,
) -> None:
    """
    Cluster samples.

    - Leiden: uses `adata.obsp['connectivities']` (requires igraph + leidenalg).
    - KMeans fallback: clusters in PCA space (no extra deps).
    """
    if method == "leiden":
        if "connectivities" not in adata.obsp:
            raise KeyError("adata.obsp['connectivities'] not found. Run bk.tl.neighbors(adata) first.")

        try:
            import igraph as ig  # type: ignore
            import leidenalg     # type: ignore
        except Exception as e:
            raise ImportError(
                "Leiden clustering requires `igraph` and `leidenalg`.\n"
                "Install with: pip install igraph leidenalg\n"
                f"Original error: {e}"
            )

    if method == "leiden":
        conn = adata.obsp["connectivities"].tocsr()
        info(f"Clustering with Leiden (resolution={resolution})")

        # Build igraph from sparse adjacency
        sources, targets = conn.nonzero()
        weights = np.asarray(conn[sources, targets]).ravel()

        g = ig.Graph(n=adata.n_obs, edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
        g.es["weight"] = weights.tolist()

        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=float(resolution),
            seed=int(random_state),
        )

        labels = np.array(part.membership, dtype=int)
        adata.obs[key_added] = labels.astype(str)
        adata.obs[key_added] = adata.obs[key_added].astype("category")

        adata.uns.setdefault("clusters", {})
        adata.uns["clusters"][key_added] = {
            "method": "leiden",
            "resolution": float(resolution),
            "random_state": int(random_state),
        }
        info(f"Leiden produced {len(np.unique(labels))} clusters stored in adata.obs['{key_added}'].")

    elif method == "kmeans":
        if use_rep not in adata.obsm:
            raise KeyError(f"adata.obsm['{use_rep}'] not found. Run bk.tl.pca(adata) first.")

        X = np.asarray(adata.obsm[use_rep], dtype=float)
        if n_pcs is not None:
            X = X[:, : int(min(n_pcs, X.shape[1]))]

        info(f"Clustering with kmeans (k={n_clusters}, rep={use_rep}, n_pcs={X.shape[1]})")
        labels = _kmeans(X, k=int(n_clusters), random_state=int(random_state))
        adata.obs[key_added] = labels.astype(str)
        adata.obs[key_added] = adata.obs[key_added].astype("category")

        adata.uns.setdefault("clusters", {})
        adata.uns["clusters"][key_added] = {
            "method": "kmeans",
            "n_clusters": int(n_clusters),
            "random_state": int(random_state),
            "use_rep": use_rep,
            "n_pcs": int(X.shape[1]),
        }
        info(f"KMeans produced {len(np.unique(labels))} clusters stored in adata.obs['{key_added}'].")

    else:
        raise ValueError("method must be 'leiden' or 'kmeans'")


def _kmeans(X: np.ndarray, k: int, random_state: int = 0, n_iter: int = 200) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    n = X.shape[0]
    k = max(2, min(k, n))

    # init: pick k random points
    centers = X[rng.choice(n, size=k, replace=False)]

    labels = np.zeros(n, dtype=int)
    for _ in range(n_iter):
        # assign
        d2 = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        new_labels = d2.argmin(axis=1)

        if np.array_equal(new_labels, labels):
            break
        labels = new_labels

        # update
        for j in range(k):
            mask = labels == j
            if mask.any():
                centers[j] = X[mask].mean(axis=0)
            else:
                centers[j] = X[rng.integers(0, n)]
    return labels