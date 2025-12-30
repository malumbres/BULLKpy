from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import anndata as ad

from ..logging import info, warn
from .._settings import settings
from ..pl._style import set_style, _savefig


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
    """
    Compute PCA on samples x genes matrix (bulk).

    Stores:
      - adata.obsm["X_pca"] : (n_obs, n_comps)
      - adata.varm["PCs"]   : (n_vars, n_comps) loadings (full space, NaN for unused genes)
      - adata.uns[key_added]["variance_ratio"] etc.
    """
    # pick matrix
    if layer is not None:
        if layer not in adata.layers:
            warn(f"layer='{layer}' not found in adata.layers; using adata.X")
            X = adata.X
            layer_used = None
        else:
            X = adata.layers[layer]
            layer_used = layer
    else:
        X = adata.X
        layer_used = None

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
    n_comps_eff = int(min(int(n_comps), max(n_obs - 1, 1), n_vars_used))

    info(f"Running PCA on matrix {n_obs} samples × {n_vars_used} genes; n_comps={n_comps_eff}")

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
    X_pca = U[:, :n_comps_eff] * S[:n_comps_eff]

    # Loadings for used genes (genes_used x comps)
    PCs_used = Vt[:n_comps_eff, :].T

    # Explained variance
    eigvals = (S**2) / max(n_obs - 1, 1)
    var_ratio = eigvals / eigvals.sum() if eigvals.sum() > 0 else np.zeros_like(eigvals)

    # Store scores
    adata.obsm["X_pca"] = X_pca

    # Store loadings in FULL gene space to satisfy AnnData alignment
    PCs_full = np.full((adata.n_vars, n_comps_eff), np.nan, dtype=float)
    PCs_full[used_mask, :] = PCs_used
    adata.varm["PCs"] = PCs_full

    # Store metadata
    adata.uns[key_added] = {
        "params": {
            "layer": layer_used,
            "n_comps": n_comps_eff,
            "center": bool(center),
            "scale": bool(scale),
            "use_highly_variable": bool(use_highly_variable),
            "n_vars_used": int(n_vars_used),
        },
        "variance": eigvals[:n_comps_eff],
        "variance_ratio": var_ratio[:n_comps_eff],
        "used_genes_mask": used_mask,
        "mean": mean_,
        "std": std_,
    }

    info(f"PCA stored in adata.obsm['X_pca'], adata.varm['PCs'], adata.uns['{key_added}']")


def pca_variance_ratio(
    adata: ad.AnnData,
    *,
    key: str = "pca",
    n_comps: int | None = None,
    cumulative: bool = True,
    figsize: tuple[float, float] = (6.5, 4.5),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scree plot of PCA variance ratio (+ optional cumulative curve).
    Requires adata.uns[key]['variance_ratio'].
    """
    set_style()

    if key not in adata.uns or "variance_ratio" not in adata.uns[key]:
        raise KeyError(f"Missing adata.uns['{key}']['variance_ratio']. Run bk.tl.pca() first.")

    vr = np.asarray(adata.uns[key]["variance_ratio"], dtype=float)
    if n_comps is not None:
        vr = vr[: int(n_comps)]

    xs = np.arange(1, len(vr) + 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(xs, vr, marker="o")
    ax.set_xlabel("PC")
    ax.set_ylabel("Explained variance ratio")
    ax.set_title("PCA variance ratio")

    if cumulative:
        ax2 = ax.twinx()
        ax2.plot(xs, np.cumsum(vr), marker="o")
        ax2.set_ylabel("Cumulative variance")

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax