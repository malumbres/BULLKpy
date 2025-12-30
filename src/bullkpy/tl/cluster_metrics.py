from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
import pandas as pd

from ..logging import info, warn

try:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
except Exception:  # pragma: no cover
    adjusted_rand_score = None
    normalized_mutual_info_score = None
    silhouette_score = None

try:
    from scipy.stats import chi2_contingency
except Exception:  # pragma: no cover
    chi2_contingency = None


def _mask_valid_pair(a: pd.Series, b: pd.Series) -> np.ndarray:
    return a.notna().to_numpy() & b.notna().to_numpy()


def _cramers_v_from_crosstab(ct: pd.DataFrame) -> float:
    """
    Cramér’s V for association between two categorical variables.
    Robust to rectangular contingency tables.
    """
    if chi2_contingency is None:
        raise ImportError("cluster_metrics (Cramér’s V) requires scipy.")

    if ct.size == 0:
        return np.nan

    chi2, _, _, _ = chi2_contingency(ct.to_numpy(), correction=False)
    n = ct.to_numpy().sum()
    if n == 0:
        return np.nan

    r, k = ct.shape
    # bias-corrected Cramér's V (Bergsma 2013 style)
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - ((k - 1) * (r - 1)) / max(n - 1, 1))
    rcorr = r - ((r - 1) ** 2) / max(n - 1, 1)
    kcorr = k - ((k - 1) ** 2) / max(n - 1, 1)
    denom = min(kcorr - 1, rcorr - 1)
    if denom <= 0:
        return np.nan
    return float(np.sqrt(phi2corr / denom))


def cluster_metrics(
    adata,
    *,
    true_key: str,
    cluster_key: str = "leiden",
    use_rep: str = "X_pca",
    layer: str | None = None,
    n_pcs: int | None = None,
    silhouette_on: Literal["rep", "X", "layer"] = "rep",
    metric: str = "euclidean",
    dropna: bool = True,
) -> dict[str, float]:
    """
    Compute agreement + quality metrics between known labels and clustering.

    Returns dict with: n_used, ari, nmi, cramers_v, silhouette

    Notes:
    - ARI/NMI/Cramér’s V compare `true_key` vs `cluster_key`.
    - Silhouette measures cluster separation on chosen data:
        * silhouette_on="rep": adata.obsm[use_rep] (recommended)
        * silhouette_on="X": adata.X
        * silhouette_on="layer": adata.layers[layer]
    """
    if true_key not in adata.obs:
        raise KeyError(f"true_key='{true_key}' not in adata.obs")
    if cluster_key not in adata.obs:
        raise KeyError(f"cluster_key='{cluster_key}' not in adata.obs")

    y_true = adata.obs[true_key]
    y_pred = adata.obs[cluster_key]

    if dropna:
        mask = _mask_valid_pair(y_true, y_pred)
        y_true = y_true.loc[mask].astype(str)
        y_pred = y_pred.loc[mask].astype(str)
    else:
        # still need to avoid NaNs for sklearn metrics
        if y_true.isna().any() or y_pred.isna().any():
            raise ValueError("NaNs present; use dropna=True (recommended).")

    n_used = int(len(y_true))
    if n_used == 0:
        raise ValueError("No samples with both labels present after filtering.")

    out: dict[str, float] = {"n_used": float(n_used)}

    # ARI/NMI
    if adjusted_rand_score is None or normalized_mutual_info_score is None:
        raise ImportError("cluster_metrics requires scikit-learn (sklearn).")

    out["ari"] = float(adjusted_rand_score(y_true, y_pred))
    out["nmi"] = float(normalized_mutual_info_score(y_true, y_pred))

    # Cramér’s V
    ct = pd.crosstab(y_true, y_pred)
    out["cramers_v"] = float(_cramers_v_from_crosstab(ct))

    # Silhouette
    if silhouette_score is None:
        raise ImportError("cluster_metrics (silhouette) requires scikit-learn.")

    # silhouette requires >=2 clusters and less than n samples
    n_clusters = pd.Series(y_pred).nunique(dropna=True)
    if n_clusters < 2 or n_clusters >= n_used:
        out["silhouette"] = np.nan
        warn("Silhouette not defined (need 2..n-1 clusters). Returning NaN.")
    else:
        if silhouette_on == "rep":
            if use_rep not in adata.obsm:
                raise KeyError(f"adata.obsm['{use_rep}'] not found (run PCA/UMAP/etc).")
            X = np.asarray(adata.obsm[use_rep], dtype=float)
            if dropna:
                X = X[mask, :]
            if n_pcs is not None:
                X = X[:, : int(n_pcs)]
        elif silhouette_on == "X":
            X = adata.X
            X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=float)
            if dropna:
                X = X[mask, :]
        else:  # "layer"
            if layer is None:
                raise ValueError("silhouette_on='layer' requires layer=...")
            if layer not in adata.layers:
                raise KeyError(f"layer='{layer}' not in adata.layers")
            X = adata.layers[layer]
            X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=float)
            if dropna:
                X = X[mask, :]

        out["silhouette"] = float(silhouette_score(X, y_pred.to_numpy(), metric=metric))

    info(
        f"cluster_metrics: n_used={n_used} ari={out['ari']:.3f} nmi={out['nmi']:.3f} "
        f"cramers_v={out['cramers_v']:.3f} silhouette={out['silhouette'] if np.isfinite(out['silhouette']) else np.nan}"
    )
    return out