from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


# -------------------------
# Helpers: ComBat EB priors
# -------------------------
def _aprior(delta_hat: np.ndarray) -> float:
    m = np.mean(delta_hat)
    s2 = np.var(delta_hat, ddof=1)
    return (2.0 * s2 + m**2) / (s2 + 1e-12)


def _bprior(delta_hat: np.ndarray) -> float:
    m = np.mean(delta_hat)
    s2 = np.var(delta_hat, ddof=1)
    return (m * s2 + m**3) / (s2 + 1e-12)


def _it_sol(
    sdat: np.ndarray,
    gamma_hat: np.ndarray,
    delta_hat: np.ndarray,
    gamma_bar: float,
    t2: float,
    a_prior: float,
    b_prior: float,
    *,
    conv: float = 1e-4,
    max_iter: int = 100,
) -> tuple[float, float]:
    """
    Iterative solution for posterior gamma* and delta* for one batch and one gene-set.
    This is the standard parametric ComBat update.
    """
    g_old = gamma_hat
    d_old = delta_hat
    n = sdat.size

    for _ in range(max_iter):
        # posterior mean of gamma (Normal)
        g_new = (t2 * n * gamma_hat + d_old * gamma_bar) / (t2 * n + d_old + 1e-12)

        # posterior mode/mean-ish of delta (Inv-Gamma)
        sum2 = np.sum((sdat - g_new) ** 2)
        d_new = (0.5 * sum2 + b_prior) / (0.5 * n + a_prior - 1.0 + 1e-12)

        if np.max(np.abs(g_new - g_old) / (np.abs(g_old) + 1e-12)) < conv and np.max(
            np.abs(d_new - d_old) / (np.abs(d_old) + 1e-12)
        ) < conv:
            return float(g_new), float(d_new)

        g_old = g_new
        d_old = d_new

    return float(g_old), float(d_old)


# -------------------------
# Main: ComBat for bulk
# -------------------------
def batch_correct_combat(
    adata: ad.AnnData,
    *,
    batch_key: str,
    layer: str | None = "log1p_cpm",
    covariates: Sequence[str] | None = None,
    key_added: str = "combat",
    overwrite: bool = False,
    inplace: bool = True,
) -> np.ndarray | None:
    """
    ComBat batch correction (Johnson et al.) for bulk expression.

    Notes
    -----
    - ComBat is intended for approximately Gaussian data:
      use log-transformed normalized expression (e.g., log1p_cpm), not raw counts.
    - Writes corrected matrix to `adata.layers[key_added]` by default.

    Parameters
    ----------
    batch_key
        adata.obs column with batch labels (categorical recommended).
    layer
        Which matrix to correct. If None, uses adata.X.
    covariates
        Optional adata.obs columns to include in the design (biological covariates to preserve).
    key_added
        Layer name to store corrected values (if overwrite=False).
    overwrite
        If True, write corrected values back into the selected layer / X.
    inplace
        If True, store results in adata; if False, return corrected matrix (samples x genes).
    """
    if batch_key not in adata.obs.columns:
        raise KeyError(f"batch_key='{batch_key}' not found in adata.obs")

    covariates = list(covariates) if covariates is not None else []
    for c in covariates:
        if c not in adata.obs.columns:
            raise KeyError(f"covariate '{c}' not found in adata.obs")

    # -------- data matrix --------
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)  # samples x genes
    n_samples, n_genes = X.shape

    # -------- batch variable --------
    batch = adata.obs[batch_key].astype("category")
    batches = list(batch.cat.categories)
    n_batch = len(batches)

    if n_batch < 2:
        warn(f"ComBat: batch_key='{batch_key}' has <2 batches; skipping correction.")
        if inplace:
            return None
        return X

    info(
        f"ComBat: correcting layer='{layer}' with batch_key='{batch_key}' "
        f"({n_batch} batches), covariates={covariates}"
    )

    # -------- design matrix: intercept + batch dummies + covariates --------
    # Intercept
    design = pd.DataFrame({"Intercept": np.ones(n_samples, dtype=float)}, index=adata.obs_names)

    # Covariates (to preserve)
    for c in covariates:
        s = adata.obs[c]
        if pd.api.types.is_numeric_dtype(s):
            design[c] = s.astype(float).values
        else:
            d = pd.get_dummies(s.astype("category"), prefix=c, drop_first=True)
            design = pd.concat([design, d], axis=1)

    # Batch dummies (NO drop_first; we need all batches for estimation)
    batch_dum = pd.get_dummies(batch, prefix="batch", drop_first=False)
    design = pd.concat([design, batch_dum], axis=1)

    design_mat = design.to_numpy(dtype=float)  # n x p
    # columns that correspond to batch indicators
    batch_cols = batch_dum.columns.tolist()
    batch_col_idx = [design.columns.get_loc(c) for c in batch_cols]

    # -------- standardize genes --------
    # Fit linear model: X = design * B + E  (least squares)
    # B_hat = (D'D)^-1 D'X
    DtD_inv = np.linalg.pinv(design_mat.T @ design_mat)
    B_hat = DtD_inv @ (design_mat.T @ X)  # p x genes

    # Grand mean: intercept + covariates (exclude batch effects)
    # Use design without batch columns for "biological part"
    keep_cols = [i for i in range(design_mat.shape[1]) if i not in batch_col_idx]
    design_keep = design_mat[:, keep_cols]
    B_keep = B_hat[keep_cols, :]  # (p_keep x genes)
    grand_mean = design_keep @ B_keep  # samples x genes

    # Pooled variance of residuals
    resid = X - (design_mat @ B_hat)
    var_pooled = np.var(resid, axis=0, ddof=1)
    var_pooled[var_pooled == 0] = 1.0

    sdat = (X - grand_mean) / np.sqrt(var_pooled)  # standardized data

    # -------- estimate batch effects on standardized data --------
    gamma_hat = np.zeros((n_batch, n_genes), dtype=float)
    delta_hat = np.zeros((n_batch, n_genes), dtype=float)

    batch_indices = []
    for i, b in enumerate(batches):
        idx = np.where(batch.to_numpy() == b)[0]
        batch_indices.append(idx)
        gamma_hat[i, :] = np.mean(sdat[idx, :], axis=0)
        delta_hat[i, :] = np.var(sdat[idx, :], axis=0, ddof=1)
        delta_hat[i, delta_hat[i, :] == 0] = 1.0

    # -------- empirical Bayes shrinkage (parametric) --------
    gamma_star = np.zeros_like(gamma_hat)
    delta_star = np.zeros_like(delta_hat)

    for i in range(n_batch):
        g_i = gamma_hat[i, :]
        d_i = delta_hat[i, :]

        gamma_bar = float(np.mean(g_i))
        t2 = float(np.var(g_i, ddof=1) if n_genes > 1 else 1.0)
        if t2 == 0:
            t2 = 1.0

        a_prior = float(_aprior(d_i))
        b_prior = float(_bprior(d_i))

        idx = batch_indices[i]
        for j in range(n_genes):
            g, d = _it_sol(
                sdat[idx, j],
                g_i[j],
                d_i[j],
                gamma_bar,
                t2,
                a_prior,
                b_prior,
            )
            gamma_star[i, j] = g
            delta_star[i, j] = d

    # -------- adjust data --------
    bayesdata = sdat.copy()
    for i in range(n_batch):
        idx = batch_indices[i]
        bayesdata[idx, :] = (bayesdata[idx, :] - gamma_star[i, :]) / np.sqrt(delta_star[i, :])

    # de-standardize
    corrected = bayesdata * np.sqrt(var_pooled) + grand_mean  # samples x genes

    # -------- store / return --------
    if inplace:
        if overwrite:
            if layer is None:
                adata.X = corrected
            else:
                if layer not in adata.layers:
                    warn(f"overwrite=True but layer='{layer}' not in adata.layers; writing to adata.layers['{layer}'].")
                adata.layers[layer] = corrected
        else:
            adata.layers[key_added] = corrected
            adata.uns.setdefault("combat", {})
            adata.uns["combat"] = {
                "params": {
                    "batch_key": batch_key,
                    "layer": layer,
                    "covariates": covariates,
                    "key_added": key_added,
                    "overwrite": overwrite,
                },
                "batches": batches,
            }
        info(f"ComBat: stored corrected matrix in {'adata.X' if (overwrite and layer is None) else f'adata.layers[{key_added!r}]' if not overwrite else f'adata.layers[{layer!r}]'}")
        return None

    return corrected