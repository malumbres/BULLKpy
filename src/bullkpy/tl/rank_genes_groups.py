from __future__ import annotations

from typing import Literal, Sequence
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn

try:
    from scipy import stats
except Exception:  # pragma: no cover
    stats = None

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None


def _get_X(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    p = np.clip(p, 0.0, 1.0)
    if multipletests is not None:
        return multipletests(p, method="fdr_bh")[1]
    # fallback BH
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    ranked[order] = p[order] * n / (np.arange(n) + 1.0)
    # enforce monotonicity
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.clip(ranked, 0.0, 1.0)


def rank_genes_groups(
    adata: ad.AnnData,
    *,
    groupby: str,
    groups: Sequence[str] | None = None,
    reference: str | Literal["rest"] = "rest",
    layer: str | None = "log1p_cpm",
    method: Literal["t-test", "wilcoxon"] = "t-test",
    corr_method: Literal["benjamini-hochberg"] = "benjamini-hochberg",
    use_abs: bool = False,
    n_genes: int = 100,
    key_added: str = "rank_genes_groups",
) -> None:
    """
    Scanpy-like ranking of genes per group.

    Statistics
    ----------
    method="t-test"    -> Welch t-test per gene, score = t statistic
    method="wilcoxon"  -> Mann–Whitney U per gene, score = z-like from U (approx)

    Output
    ------
    Stores results in adata.uns[key_added] with:
      - params
      - names[group], scores[group], logfoldchanges[group], pvals[group], pvals_adj[group]
      - mean_group[group], mean_ref[group] (extra, useful for bulk)
    """
    if stats is None:
        raise ImportError("rank_genes_groups requires scipy (scipy.stats).")

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    X = _get_X(adata, layer)  # samples x genes
    obs = adata.obs.copy()

    g = obs[groupby].astype("category")
    all_groups = list(g.cat.categories)

    if groups is None:
        groups = all_groups
    else:
        groups = [str(x) for x in groups]
        missing = [x for x in groups if x not in set(all_groups)]
        if missing:
            raise KeyError(f"Requested groups not found in {groupby}: {missing}")

    if reference != "rest":
        reference = str(reference)
        if reference not in set(all_groups):
            raise KeyError(f"reference='{reference}' not found in {groupby} categories")

    info(
        f"rank_genes_groups: groupby='{groupby}', method='{method}', "
        f"layer='{layer}', reference='{reference}', n_genes={n_genes}"
    )

    genes = adata.var_names.astype(str).to_numpy()
    n_genes_total = adata.n_vars

    # Prepare outputs in Scanpy-like per-group dict-of-arrays
    out_names = {}
    out_scores = {}
    out_logfc = {}
    out_pvals = {}
    out_pvals_adj = {}
    out_mean_g = {}
    out_mean_r = {}

    # Precompute group indices
    group_to_idx = {cat: np.where(g.to_numpy() == cat)[0] for cat in all_groups}

    eps = 1e-9

    for grp in groups:
        idx_g = group_to_idx[grp]
        if idx_g.size == 0:
            warn(f"Group '{grp}' has 0 samples; skipping.")
            continue

        if reference == "rest":
            idx_r = np.setdiff1d(np.arange(adata.n_obs), idx_g)
        else:
            idx_r = group_to_idx[reference]

        if idx_r.size == 0:
            warn(f"Reference for group '{grp}' has 0 samples; skipping.")
            continue

        Xg = X[idx_g, :]  # ng x genes
        Xr = X[idx_r, :]  # nr x genes

        mean_g = Xg.mean(axis=0)
        mean_r = Xr.mean(axis=0)

        # log2FC: since layer may already be log-scale, we still report delta / ln2
        # (Scanpy reports logfoldchanges in the space of input; for log1p this is ok as "logFC")
        log2fc = (mean_g - mean_r) / np.log(2.0)

        if method == "t-test":
            # Welch t-test per gene; vectorized via scipy is limited -> do manual moments
            ng = Xg.shape[0]
            nr = Xr.shape[0]

            vg = Xg.var(axis=0, ddof=1)
            vr = Xr.var(axis=0, ddof=1)

            se = np.sqrt(vg / max(ng, 1) + vr / max(nr, 1)) + eps
            tstat = (mean_g - mean_r) / se

            # Welch-Satterthwaite df
            df = (vg / ng + vr / nr) ** 2 / ((vg**2) / (ng**2 * (ng - 1 + eps)) + (vr**2) / (nr**2 * (nr - 1 + eps)) + eps)
            pvals = 2.0 * stats.t.sf(np.abs(tstat), df=df)

            scores = tstat

        elif method == "wilcoxon":
            # Mann–Whitney U per gene (slow-ish but OK for bulk sizes)
            # We compute U and approximate z for a score.
            pvals = np.ones(n_genes_total, dtype=float)
            scores = np.zeros(n_genes_total, dtype=float)

            # ranks per gene -> compute with scipy per gene
            for j in range(n_genes_total):
                a = Xg[:, j]
                b = Xr[:, j]
                # skip constant vectors
                if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
                    pvals[j] = 1.0
                    scores[j] = 0.0
                    continue
                res = stats.mannwhitneyu(a, b, alternative="two-sided")
                pvals[j] = float(res.pvalue)

                # approximate z-score from U (no tie correction here; good enough v1)
                U = float(res.statistic)
                ng = a.size
                nr = b.size
                mu = ng * nr / 2.0
                sigma = np.sqrt(ng * nr * (ng + nr + 1) / 12.0) + eps
                z = (U - mu) / sigma
                scores[j] = z
        else:
            raise ValueError("method must be 't-test' or 'wilcoxon'")

        # multiple testing
        if corr_method != "benjamini-hochberg":
            raise ValueError("Only corr_method='benjamini-hochberg' supported in v1.")
        pvals_adj = _bh_fdr(pvals)

        # ranking
        if use_abs:
            rank_idx = np.argsort(np.abs(scores))[::-1]
        else:
            rank_idx = np.argsort(scores)[::-1]

        rank_idx = rank_idx[: int(n_genes)]

        out_names[grp] = genes[rank_idx]
        out_scores[grp] = scores[rank_idx].astype(float)
        out_logfc[grp] = log2fc[rank_idx].astype(float)
        out_pvals[grp] = pvals[rank_idx].astype(float)
        out_pvals_adj[grp] = pvals_adj[rank_idx].astype(float)
        out_mean_g[grp] = mean_g[rank_idx].astype(float)
        out_mean_r[grp] = mean_r[rank_idx].astype(float)

    # store Scanpy-like structure
    adata.uns[key_added] = {
        "params": {
            "groupby": groupby,
            "groups": list(groups),
            "reference": reference,
            "layer": layer,
            "method": method,
            "corr_method": corr_method,
            "use_abs": use_abs,
            "n_genes": int(n_genes),
        },
        "names": out_names,
        "scores": out_scores,
        "logfoldchanges": out_logfc,
        "pvals": out_pvals,
        "pvals_adj": out_pvals_adj,
        # bulk-friendly extras
        "mean_group": out_mean_g,
        "mean_ref": out_mean_r,
    }

    info(f"rank_genes_groups: stored results in adata.uns['{key_added}']")