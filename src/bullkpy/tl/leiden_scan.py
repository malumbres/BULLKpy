from __future__ import annotations

from typing import Sequence
import numpy as np
import pandas as pd
import anndata as ad

from ..logging import info, warn
from .association import categorical_association
from .clustering import cluster


def _cramers_v(x: pd.Series, y: pd.Series) -> float:
    """Cramér's V for two categorical vectors."""
    xt = pd.crosstab(x, y)
    if xt.size == 0:
        return np.nan
    n = xt.to_numpy().sum()
    if n == 0:
        return np.nan

    # chi2 without scipy dependency (manual)
    obs = xt.to_numpy(dtype=float)
    row_sum = obs.sum(axis=1, keepdims=True)
    col_sum = obs.sum(axis=0, keepdims=True)
    exp = row_sum @ col_sum / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum((obs - exp) ** 2 / exp)

    r, k = obs.shape
    if min(r, k) <= 1:
        return 0.0
    return float(np.sqrt((chi2 / n) / (min(r - 1, k - 1))))


def leiden_resolution_scan(
    adata: ad.AnnData,
    *,
    true_key: str,
    resolutions: Sequence[float] | None = None,
    base_key: str = "leiden",
    store_key: str = "leiden_scan",
    use_rep: str = "X_pca",
    n_pcs: int = 20,
    n_neighbors: int = 15,
    metric: str = "euclidean",
    recompute_neighbors: bool = False,
    inplace: bool = True,
    random_state: int = 0,
) -> pd.DataFrame:
    """
    Scan multiple Leiden resolutions and score vs a ground-truth annotation.
    Scan many Leiden resolutions and score against a 'true' categorical label
    using ARI / NMI / Cramér’s V.

    Notes:
      - Uses bk.tl.neighbors() to build adata.obsp['connectivities'] if missing.
      - Uses bk.tl.cluster(method='leiden') to compute clusters for each resolution.
      - Computes ARI, NMI (if sklearn available) + Cramér's V (always).
    """
    if true_key not in adata.obs.columns:
        raise KeyError(f"true_key='{true_key}' not in adata.obs")

    if resolutions is None:
        resolutions = np.round(np.linspace(0.2, 2.0, 10), 3)

    # local imports to avoid circular imports
    from .neighbors import neighbors
    from .clustering import cluster

    # Ensure neighbors graph exists once (resolution changes don't require rebuilding)
    if recompute_neighbors or ("connectivities" not in adata.obsp):
        info(
            f"Computing neighbors graph for scan (use_rep={use_rep}, n_pcs={n_pcs}, "
            f"n_neighbors={n_neighbors}, metric={metric})"
        )
        neighbors(
            adata,
            n_neighbors=int(n_neighbors),
            n_pcs=int(n_pcs),
            use_rep=str(use_rep),
            metric=str(metric),
        )

    # Optional sklearn metrics
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
        have_sklearn = True
    except Exception:
        have_sklearn = False
        warn("sklearn not available → ARI/NMI will be NaN (Cramér's V still computed).")

    y_true_all = adata.obs[true_key]

    rows: list[dict] = []
    for r in resolutions:
        key_added = f"{base_key}_{r:g}"
        info(f"Leiden resolution={r:g}")

        cluster(
            adata,
            method="leiden",
            resolution=float(r),
            key_added=key_added,
            use_rep=str(use_rep),
            n_pcs=int(n_pcs),
            random_state=int(random_state),
        )

        y_pred_all = adata.obs[key_added]

        # drop NaNs for metrics
        mask = (~pd.isna(y_true_all)) & (~pd.isna(y_pred_all))
        y_true = y_true_all[mask].astype(str)
        y_pred = y_pred_all[mask].astype(str)

        n_clusters = int(pd.Series(y_pred).nunique(dropna=True))

        ari = np.nan
        nmi = np.nan
        if have_sklearn and len(y_true) > 0:
            ari = float(adjusted_rand_score(y_true, y_pred))
            nmi = float(normalized_mutual_info_score(y_true, y_pred))

        cv = _cramers_v(pd.Series(y_true), pd.Series(y_pred))

        rows.append(
            {
                "resolution": float(r),
                "key": key_added,
                "n_clusters": n_clusters,
                "ARI": ari,
                "NMI": nmi,
                "cramers_v": cv,
                "n_used": int(mask.sum()),
                "n_missing": int((~mask).sum()),
            }
        )

    df = pd.DataFrame(rows).sort_values("resolution").reset_index(drop=True)
    if inplace:
        adata.uns[store_key] = df
    return df