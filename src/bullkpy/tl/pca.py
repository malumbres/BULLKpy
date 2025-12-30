from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def pca(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    n_comps: int = 50,
    center: bool = True,
    scale: bool = False,
    use_highly_variable: bool = False,
    key_added: str = "pca",
) -> None:
    X = adata.layers[layer] if layer is not None else adata.X

    # Track which genes are used so we can write loadings back to full var space
    used_mask = np.ones(adata.n_vars, dtype=bool)

    if use_highly_variable:
        if "highly_variable" not in adata.var.columns:
            warn("use_highly_variable=True but adata.var['highly_variable'] not found; using all genes.")
        else:
            used_mask = adata.var["highly_variable"].to_numpy(bool)
            if used_mask.sum() == 0:
                raise ValueError("adata.var['highly_variable'] has 0 True values.")
            X = X[:, used_mask]

    if sp.issparse(X):
        X = X.toarray()

    X = np.asarray(X, dtype=float)  # samples x genes_used
    n_obs, n_vars_used = X.shape
    n_comps = int(min(n_comps, n_obs - 1, n_vars_used))

    info(f"Running PCA on matrix {n_obs} samples × {n_vars_used} genes; n_comps={n_comps}")

    # Center/scale per gene
    mean_ = X.mean(axis=0) if center else np.zeros(n_vars_used, dtype=float)
    Xc = X - mean_

    std_ = None
    if scale:
        std_ = Xc.std(axis=0, ddof=1)
        std_[std_ == 0] = 1.0
        Xc = Xc / std_

    # SVD: Xc = U S Vt
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)

    # Scores (samples x comps)
    X_pca = U[:, :n_comps] * S[:n_comps]

    # Loadings for used genes (genes_used x comps)
    PCs_used = Vt[:n_comps, :].T

    # Explained variance
    eigvals = (S**2) / max(n_obs - 1, 1)
    var_ratio = eigvals / eigvals.sum()

    # Store scores
    adata.obsm["X_pca"] = X_pca

    # Store loadings in FULL gene space to satisfy AnnData alignment
    PCs_full = np.full((adata.n_vars, n_comps), np.nan, dtype=float)
    PCs_full[used_mask, :] = PCs_used
    adata.varm["PCs"] = PCs_full

    # Store metadata
    adata.uns[key_added] = {
        "params": {
            "layer": layer,
            "n_comps": n_comps,
            "center": center,
            "scale": scale,
            "use_highly_variable": use_highly_variable,
            "n_vars_used": int(n_vars_used),
        },
        "variance": eigvals[:n_comps],
        "variance_ratio": var_ratio[:n_comps],
        "used_genes_mask": used_mask,  # helpful for debugging / downstream
        "mean": mean_,
        "std": std_,
    }

    info(f"PCA stored in adata.obsm['X_pca'], adata.varm['PCs'], adata.uns['{key_added}']")