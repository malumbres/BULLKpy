"""Statistical helpers shared across modules."""
from __future__ import annotations

import numpy as np

__all__ = ["bh_fdr"]


def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values.

    NaN p-values mark tests that could not be run: they stay NaN and do not count
    towards the number of tests.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out

    order = np.argsort(p[ok])
    q = p[ok][order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = np.clip(q, 0.0, 1.0)
    out[ok] = adjusted
    return out
