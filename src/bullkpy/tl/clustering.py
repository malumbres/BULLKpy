from __future__ import annotations

from typing import Literal, Optional, Dict, Any, Sequence
import numpy as np
import pandas as pd
import anndata as ad

from dataclasses import dataclass

from ..logging import info, warn
from .associations import categorical_association



##########################################################
### CLUSTERING
##########################################################



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




##########################################################
### CLUSTER METRICS
##########################################################

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

# If categorical_association is in the same module, call it directly.
# Otherwise, import it properly, e.g.:
# from .association import categorical_association

def categorical_confusion(
    adata: ad.AnnData,
    *,
    key1: str,
    key2: str,
    normalize: Literal["none", "row", "col", "all"] = "row",
    dropna: bool = True,
    min_count: int = 1,
) -> Dict[str, Any]:
    """
    Compute a confusion-style table for two categorical obs columns plus association metrics.

    Returns dict with:
      - table: raw contingency table (counts)
      - matrix: normalized matrix (float) according to `normalize`
      - normalize: normalize mode used
      - metrics: dict of association metrics (chi2, cramers_v, ari, nmi, etc. if available)
      - key1, key2: input keys

    Notes:
      - No plotting here (bk.tl is compute-only). Use bk.pl for heatmap rendering.
      - `dropna=True` removes rows where either key is NA before tabulation.
      - `min_count` can be used to drop rare categories after tabulation (optional).
    """
    if key1 not in adata.obs.columns:
        raise KeyError(f"{key1!r} not found in adata.obs")
    if key2 not in adata.obs.columns:
        raise KeyError(f"{key2!r} not found in adata.obs")
    if normalize not in {"none", "row", "col", "all"}:
        raise ValueError("normalize must be one of {'none','row','col','all'}")

    # Pull columns
    s1 = adata.obs[key1]
    s2 = adata.obs[key2]

    # Clean NA handling
    if dropna:
        ok = s1.notna() & s2.notna()
        s1 = s1.loc[ok]
        s2 = s2.loc[ok]

    # Ensure categorical-ish labels (strings are fine)
    s1 = s1.astype(str)
    s2 = s2.astype(str)

    # Raw contingency table
    tab = pd.crosstab(s1, s2, dropna=False)

    # Optionally drop ultra-rare categories
    if int(min_count) > 1:
        row_keep = tab.sum(axis=1) >= int(min_count)
        col_keep = tab.sum(axis=0) >= int(min_count)
        tab = tab.loc[row_keep, col_keep]

    # Normalized matrix
    mat = tab.to_numpy(dtype=float)
    if normalize == "row":
        denom = np.maximum(mat.sum(axis=1, keepdims=True), 1.0)
        mat = mat / denom
    elif normalize == "col":
        denom = np.maximum(mat.sum(axis=0, keepdims=True), 1.0)
        mat = mat / denom
    elif normalize == "all":
        denom = np.maximum(mat.sum(), 1.0)
        mat = mat / denom
    # normalize == "none": keep counts

    # Metrics: prefer your existing categorical_association implementation
    metrics: Dict[str, Any] = {}
    try:
        res = categorical_association(adata, key1=key1, key2=key2)  # noqa: F821
        # pull everything except the table to avoid duplication
        metrics = {k: v for k, v in res.items() if k != "table"}
        # but use our cleaned/filtered tab as authoritative
    except Exception as e:
        metrics = {"_warning": f"categorical_association failed: {e}"}

    return {
        "key1": key1,
        "key2": key2,
        "normalize": normalize,
        "table": tab,
        "matrix": mat,
        "metrics": metrics,
    }

