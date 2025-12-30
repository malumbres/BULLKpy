from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def _as_dense(X):
    if sp.issparse(X):
        return X.toarray()
    return np.asarray(X)


def deseq2_size_factors(counts: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """
    DESeq2 median-of-ratios size factors.
    counts: (n_samples, n_genes) raw counts (non-negative).
    """
    counts = np.asarray(counts, dtype=float)
    # geometric means per gene (ignore zeros)
    with np.errstate(divide="ignore", invalid="ignore"):
        logc = np.log(counts)
    logc[~np.isfinite(logc)] = np.nan
    gmean = np.exp(np.nanmean(logc, axis=0))
    # genes with gmean==0 are uninformative
    valid = gmean > 0

    ratios = counts[:, valid] / gmean[valid][None, :]
    ratios[~np.isfinite(ratios)] = np.nan
    sf = np.nanmedian(ratios, axis=1)
    sf[~np.isfinite(sf)] = 1.0
    sf[sf <= 0] = 1.0

    # normalize to have geometric mean 1
    sf = sf / np.exp(np.mean(np.log(sf + eps)))
    return sf


def estimate_dispersion_mom(norm_counts: np.ndarray, eps: float = 1e-8) -> tuple[np.ndarray, np.ndarray]:
    """
    Method-of-moments NB dispersion on normalized counts.
    Returns (mu, alpha_hat) per gene.

    For NB2 Var(Y) = mu + alpha * mu^2  =>  alpha = max((var - mu) / mu^2, eps)
    """
    Y = np.asarray(norm_counts, dtype=float)
    mu = Y.mean(axis=0)
    var = Y.var(axis=0, ddof=1)
    alpha = (var - mu) / (mu**2 + eps)
    alpha = np.maximum(alpha, eps)
    return mu, alpha


def shrink_dispersion_to_trend(mu: np.ndarray, alpha: np.ndarray, prior_df: float = 10.0, eps: float = 1e-8) -> np.ndarray:
    """
    Simple shrinkage: fit log(alpha) ~ a + b*log(mu) and shrink toward trend.
    prior_df controls strength (bigger = more shrinkage).
    """
    mu = np.asarray(mu, dtype=float)
    alpha = np.asarray(alpha, dtype=float)

    m = mu > 0
    x = np.log(mu[m] + eps)
    y = np.log(alpha[m] + eps)

    # robust-ish linear fit
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a, b = coef
    alpha_trend = np.exp(a + b * np.log(mu + eps))

    # shrink in log space (empirical Bayes-ish)
    # weight ~ prior_df / (prior_df + n_samples) is a simple approximation.
    # We'll let prior_df be the shrinkage strength directly.
    w = prior_df / (prior_df + 1.0)
    log_alpha = np.log(alpha + eps)
    log_trend = np.log(alpha_trend + eps)
    alpha_shrunk = np.exp((1 - w) * log_alpha + w * log_trend)

    return np.maximum(alpha_shrunk, eps)