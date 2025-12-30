from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
import anndata as ad

from ..logging import info, warn


def neighbors(
    adata: ad.AnnData,
    *,
    n_neighbors: int = 15,
    n_pcs: int = 20,
    use_rep: str = "X_pca",
    metric: str = "euclidean",  # "euclidean" or "cosine"
    key_added: str = "neighbors",
) -> None:
    """
    Compute kNN graph on samples (obs) using PCA representation.

    Stores:
      - adata.obsp["distances"] (CSR sparse)
      - adata.obsp["connectivities"] (CSR sparse; Gaussian kernel)
      - adata.uns[key_added] with parameters
    """
    if use_rep not in adata.obsm:
        raise KeyError(f"adata.obsm['{use_rep}'] not found. Run bk.tl.pca(adata) first.")

    X = np.asarray(adata.obsm[use_rep], dtype=float)
    if X.ndim != 2:
        raise ValueError(f"{use_rep} must be 2D; got shape {X.shape}")

    if n_pcs is not None:
        n_pcs = int(min(n_pcs, X.shape[1]))
        X = X[:, :n_pcs]

    n_obs = X.shape[0]
    k = int(min(n_neighbors, n_obs - 1))
    info(f"Computing neighbors: k={k}, rep={use_rep}, n_pcs={n_pcs}, metric={metric}")

    X_use = X
    if metric == "cosine":
        # cosine distance = 1 - cosine similarity
        # Use L2 normalization then euclidean in that space (monotonic with cosine).
        denom = np.linalg.norm(X_use, axis=1, keepdims=True)
        denom[denom == 0] = 1.0
        X_use = X_use / denom
        # For normalized vectors, euclidean distance relates to cosine similarity.

    elif metric != "euclidean":
        raise ValueError("metric must be 'euclidean' or 'cosine'")

    tree = cKDTree(X_use)
    # query k+1 because the closest neighbor is the point itself (dist=0)
    dists, idxs = tree.query(X_use, k=k + 1)

    # remove self
    dists = dists[:, 1:]
    idxs = idxs[:, 1:]

    # Build sparse distance matrix
    rows = np.repeat(np.arange(n_obs), k)
    cols = idxs.reshape(-1)
    data = dists.reshape(-1)

    distances = sp.csr_matrix((data, (rows, cols)), shape=(n_obs, n_obs))

    # Symmetrize by taking min distance (or keep directed; Scanpy effectively uses symmetric connectivities)
    distances = _symmetrize_min(distances)

    # Convert distances -> connectivities with a local scaling kernel
    # Simple, robust: sigma per node = distance to kth neighbor (after sym)

    # sigma per node = distance to kth neighbor (after sym)
    mx = distances.max(axis=1)  # may be sparse matrix on some SciPy versions
    if sp.issparse(mx):
        kth = mx.toarray().ravel()
    else:
        kth = np.asarray(mx).ravel()

    kth = kth.astype(float, copy=False)
    pos = kth > 0
    kth[~pos] = np.median(kth[pos]) if np.any(pos) else 1.0


    # Gaussian kernel: exp(-d^2 / (2 sigma_i sigma_j))
    connectivities = _gaussian_connectivities(distances, kth)

    adata.obsp["distances"] = distances
    adata.obsp["connectivities"] = connectivities

    adata.uns[key_added] = {
        "params": {
            "n_neighbors": int(n_neighbors),
            "n_pcs": None if n_pcs is None else int(n_pcs),
            "use_rep": use_rep,
            "metric": metric,
        }
    }

    info("Stored neighbors graph in adata.obsp['distances'] and adata.obsp['connectivities'].")


def _symmetrize_min(A: sp.csr_matrix) -> sp.csr_matrix:
    A = A.tocsr()
    AT = A.transpose().tocsr()
    # elementwise min where both present; if only one present keep that
    # Use: min(A, AT) + abs(A-AT) to keep existing edges. Safer: take minimum of nonzeros via sparse ops
    # Simple approach: keep the smaller where both exist by taking (A + AT)/2 on union but that changes values.
    # We'll do: union with minimum by constructing both and taking elementwise minimum on dense mask is too heavy.
    # Practical: take symmetric by: A_sym = A.minimum(AT) + (A - A.minimum(AT)) + (AT - A.minimum(AT))
    M = A.minimum(AT)
    A_sym = M + (A - M) + (AT - M)
    A_sym.eliminate_zeros()
    return A_sym


def _gaussian_connectivities(distances: sp.csr_matrix, sigma: np.ndarray) -> sp.csr_matrix:
    D = distances.tocsr()
    rows, cols = D.nonzero()
    d = np.asarray(D[rows, cols]).ravel()

    sig_i = sigma[rows]
    sig_j = sigma[cols]
    denom = 2.0 * sig_i * sig_j
    denom[denom == 0] = np.median(denom[denom > 0]) if np.any(denom > 0) else 1.0

    w = np.exp(-(d ** 2) / denom)
    W = sp.csr_matrix((w, (rows, cols)), shape=D.shape)
    # ensure symmetry
    W = 0.5 * (W + W.T)
    W.eliminate_zeros()
    return W