from __future__ import annotations

from typing import Mapping, Sequence
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn


def _get_X(adata: ad.AnnData, layer: str | None) -> sp.spmatrix | np.ndarray:
    if layer is None:
        return adata.X
    if layer in adata.layers:
        return adata.layers[layer]
    raise KeyError(f"layer='{layer}' not found in adata.layers")


def _as_1d(a) -> np.ndarray:
    return np.asarray(a).ravel()


def _mean_per_gene(X) -> np.ndarray:
    # shape (n_vars,)
    if sp.issparse(X):
        return _as_1d(X.mean(axis=0))
    return np.asarray(X, dtype=float).mean(axis=0)


def _std_per_gene(X, mean: np.ndarray | None = None, ddof: int = 0) -> np.ndarray:
    # shape (n_vars,)
    if mean is None:
        mean = _mean_per_gene(X)
    mean = np.asarray(mean, dtype=float)

    if sp.issparse(X):
        ex2 = _as_1d(X.power(2).mean(axis=0))
        var = np.maximum(ex2 - mean**2, 0.0)
        if ddof != 0:
            n = X.shape[0]
            if n - ddof > 0:
                var = var * (n / (n - ddof))
        return np.sqrt(var)

    Xd = np.asarray(X, dtype=float)
    return Xd.std(axis=0, ddof=ddof)


def _subset_matrix(X, idx: np.ndarray) -> np.ndarray:
    # returns dense (n_obs, len(idx))
    if sp.issparse(X):
        return X[:, idx].toarray()
    return np.asarray(X[:, idx], dtype=float)


def _validate_gene_list(adata: ad.AnnData, genes: Sequence[str], *, name: str = "genes") -> list[str]:
    genes = [str(g) for g in genes]
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        warn(f"{name}: {len(missing)} genes not in adata.var_names (first 10): {missing[:10]}")
    if len(present) == 0:
        raise ValueError(f"None of the provided {name} were found in adata.var_names.")
    return present


def score_genes(
    adata: ad.AnnData,
    genes: Sequence[str],
    *,
    score_name: str = "score",
    layer: str | None = None,
    gene_pool: Sequence[str] | None = None,
    ctrl_size: int = 50,
    n_bins: int = 25,
    random_state: int | None = 0,
    scale: bool = False,
) -> None:
    """
    Scanpy-like signature scoring (bulk-friendly).

    Score per sample:
        mean(expression[genes]) - mean(expression[control_genes])

    Control genes are sampled from gene_pool (default: all genes) by matching
    expression bins of signature genes (similar to scanpy.tl.score_genes).

    Stores result in: adata.obs[score_name]
    """
    X = _get_X(adata, layer)

    sig = _validate_gene_list(adata, genes, name="genes")

    if gene_pool is None:
        pool = list(map(str, adata.var_names))
    else:
        pool = _validate_gene_list(adata, gene_pool, name="gene_pool")

    # remove signature genes from pool
    pool_set = set(pool)
    for g in sig:
        pool_set.discard(g)
    pool = list(pool_set)
    if len(pool) == 0:
        raise ValueError("gene_pool became empty after removing signature genes.")

    # name -> index
    name_to_idx = {str(g): i for i, g in enumerate(adata.var_names)}
    sig_idx = np.array([name_to_idx[g] for g in sig], dtype=int)
    pool_idx = np.array([name_to_idx[g] for g in pool], dtype=int)

    # binning by mean expression
    gene_means = _mean_per_gene(X)
    pool_means = gene_means[pool_idx]
    sig_means = gene_means[sig_idx]

    n_bins_eff = int(max(1, n_bins))
    if n_bins_eff == 1:
        bin_edges = np.array([np.min(pool_means) - 1e-9, np.max(pool_means) + 1e-9], dtype=float)
    else:
        qs = np.linspace(0, 1, n_bins_eff + 1)
        bin_edges = np.quantile(pool_means, qs)
        bin_edges = np.unique(bin_edges)  # handle ties
        if len(bin_edges) < 2:
            bin_edges = np.array([np.min(pool_means) - 1e-9, np.max(pool_means) + 1e-9], dtype=float)

    pool_bins = np.digitize(pool_means, bin_edges[1:-1], right=True)

    rng = np.random.default_rng(random_state)

    ctrl_indices: list[int] = []
    for m in sig_means:
        b = int(np.digitize(m, bin_edges[1:-1], right=True))
        candidates = pool_idx[pool_bins == b]
        if candidates.size == 0:
            candidates = pool_idx  # fallback
        k = int(min(max(1, ctrl_size), candidates.size))
        picked = rng.choice(candidates, size=k, replace=False)
        ctrl_indices.extend(picked.tolist())

    ctrl_idx = np.unique(np.array(ctrl_indices, dtype=int))
    if ctrl_idx.size == 0:
        raise ValueError("Failed to sample any control genes. Try lowering n_bins or ctrl_size.")

    X_sig = _subset_matrix(X, sig_idx)
    X_ctrl = _subset_matrix(X, ctrl_idx)

    if scale:
        mu = gene_means
        sd = _std_per_gene(X, mean=gene_means, ddof=0)
        sd[sd == 0] = 1.0
        X_sig = (X_sig - mu[sig_idx]) / sd[sig_idx]
        X_ctrl = (X_ctrl - mu[ctrl_idx]) / sd[ctrl_idx]

    score = X_sig.mean(axis=1) - X_ctrl.mean(axis=1)
    adata.obs[score_name] = np.asarray(score, dtype=float)

    info(f"score_genes: stored '{score_name}' in adata.obs (n_sig={len(sig)}, n_ctrl={ctrl_idx.size}, layer={layer})")


def score_genes_dict(
    adata: ad.AnnData,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    layer: str | None = None,
    prefix: str = "",
    suffix: str = "_score",
    gene_pool: Sequence[str] | None = None,
    ctrl_size: int = 50,
    n_bins: int = 25,
    random_state: int | None = 0,
    scale: bool = False,
) -> None:
    """
    Compute multiple signature scores and store each in adata.obs.

    Output column name:
        f"{prefix}{name}{suffix}"
    """
    if not isinstance(gene_sets, Mapping) or len(gene_sets) == 0:
        raise ValueError("gene_sets must be a non-empty mapping of name -> genes.")

    for name, genes in gene_sets.items():
        out = f"{prefix}{name}{suffix}"
        score_genes(
            adata,
            genes,
            score_name=out,
            layer=layer,
            gene_pool=gene_pool,
            ctrl_size=ctrl_size,
            n_bins=n_bins,
            random_state=random_state,
            scale=scale,
        )


def score_genes_cell_cycle(
    adata: ad.AnnData,
    *,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    layer: str | None = None,
    gene_pool: Sequence[str] | None = None,
    ctrl_size: int = 50,
    n_bins: int = 25,
    random_state: int | None = 0,
    scale: bool = False,
    s_score: str = "S_score",
    g2m_score: str = "G2M_score",
    phase: str = "phase",
) -> None:
    """
    Score cell cycle phases similar to scanpy.tl.score_genes_cell_cycle.

    Writes to:
      - adata.obs[s_score]   (float)
      - adata.obs[g2m_score] (float)
      - adata.obs[phase]     (category: "S", "G2M", "G1")

    Phase calling rule (Scanpy-like):
      - if S_score > G2M_score and S_score > 0 -> "S"
      - elif G2M_score > S_score and G2M_score > 0 -> "G2M"
      - else -> "G1"
    """
    if s_genes is None or len(s_genes) == 0:
        raise ValueError("Provide non-empty s_genes.")
    if g2m_genes is None or len(g2m_genes) == 0:
        raise ValueError("Provide non-empty g2m_genes.")

    # Compute scores (stored in obs)
    score_genes(
        adata,
        s_genes,
        score_name=s_score,
        layer=layer,
        gene_pool=gene_pool,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        random_state=random_state,
        scale=scale,
    )
    score_genes(
        adata,
        g2m_genes,
        score_name=g2m_score,
        layer=layer,
        gene_pool=gene_pool,
        ctrl_size=ctrl_size,
        n_bins=n_bins,
        random_state=None if random_state is None else (int(random_state) + 1),
        scale=scale,
    )

    s = adata.obs[s_score].to_numpy(dtype=float)
    g = adata.obs[g2m_score].to_numpy(dtype=float)

    # Phase assignment (close to Scanpy; robust for bulk too)
    ph = np.full(adata.n_obs, "G1", dtype=object)
    ph[(s > g) & (s > 0)] = "S"
    ph[(g > s) & (g > 0)] = "G2M"

    adata.obs[phase] = pd.Categorical(ph, categories=["G1", "S", "G2M"], ordered=True)

    info(
        "score_genes_cell_cycle: stored "
        f"obs['{s_score}'], obs['{g2m_score}'], obs['{phase}'] "
        f"(layer={layer}, scale={scale})"
    )

    # Helpful warning if nothing matched (common if genes missing / not expressed)
    if (adata.obs[phase] == "G1").mean() > 0.95:
        warn(
            "Most samples assigned to G1. This can happen if S/G2M genes are missing "
            "from var_names, or if expression values are not comparable. "
            "Consider using a more suitable layer (e.g. log1p CPM/TPM) and/or check gene symbols."
        )