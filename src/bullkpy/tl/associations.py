from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
import anndata as ad

from ..logging import info, warn
from ..get import vector

try:
    from scipy import stats
except Exception:  # pragma: no cover
    stats = None

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None


Method = Literal["mwu", "ttest", "kruskal", "anova"]


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    if multipletests is not None:
        return multipletests(p, method="fdr_bh")[1]
    # fallback BH
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    ranked[order] = np.arange(1, n + 1)
    q = p * n / ranked
    # enforce monotonicity
    q_ordered = q[order]
    q_ordered = np.minimum.accumulate(q_ordered[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q_ordered, 0, 1)
    return out


def _rank_biserial_from_mwu(u_stat: float, n1: int, n2: int) -> float:
    # rank-biserial correlation in [-1, 1]
    # r_rb = 1 - (2U)/(n1*n2) when U is for group1 "wins"
    denom = float(n1 * n2)
    if denom <= 0:
        return np.nan
    return 1.0 - (2.0 * float(u_stat)) / denom


def _log2fc(mean_a: float, mean_b: float, eps: float = 1e-9) -> float:
    return float(np.log2((mean_a + eps) / (mean_b + eps)))


def rank_genes_categorical(
    adata: ad.AnnData,
    *,
    groupby: str,
    group: str | None = None,
    reference: str = "rest",   # "rest" or a specific group name
    layer: str | None = None,
    genes: list[str] | None = None,   # default: all genes
    method: Method = "mwu",
    store_key: str | None = None,     # e.g. "assoc:Project_ID:ACC_vs_rest"
) -> pd.DataFrame:
    """
    Scanpy-like (but bulk-friendly) rank of genes associated with a categorical obs.

    Two-group comparison:
      - MWU (default): effect_size = rank-biserial correlation
      - t-test: effect_size = Cohen's d (approx)
    Multi-group comparison:
      - Kruskal or ANOVA: effect_size = (rough) eta^2

    Returns a DataFrame with:
      gene, pval, qval, effect_size, mean_group, mean_ref, log2FC
    """
    if stats is None:
        raise ImportError("rank_genes_categorical requires scipy (scipy.stats).")

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    g = adata.obs[groupby].astype("category")
    cats = list(g.cat.categories)

    if group is None:
        if len(cats) != 2:
            raise ValueError(
                f"group is None but {groupby} has {len(cats)} categories. "
                f"Provide group='...' (and optionally reference=...)."
            )
        group = str(cats[0])
        reference = str(cats[1])

    group = str(group)
    if group not in cats:
        raise KeyError(f"group='{group}' not in {groupby} categories")

    if reference != "rest" and str(reference) not in cats:
        raise KeyError(f"reference='{reference}' not in categories (or use reference='rest')")

    genes_use = list(adata.var_names) if genes is None else [str(x) for x in genes]
    missing = [x for x in genes_use if x not in adata.var_names]
    if missing:
        raise KeyError(f"Genes not in adata.var_names (first 10): {missing[:10]}")

    # masks
    m_group = (g.astype(str) == group).to_numpy()
    if reference == "rest":
        m_ref = ~m_group
        ref_label = "rest"
    else:
        m_ref = (g.astype(str) == str(reference)).to_numpy()
        ref_label = str(reference)

    n1 = int(m_group.sum())
    n2 = int(m_ref.sum())
    if n1 < 2 or n2 < 2:
        raise ValueError(f"Not enough samples: n_group={n1}, n_ref={n2}")

    pvals = []
    effects = []
    mean_a = []
    mean_b = []

    # do one gene at a time (bulk n_obs can be large but n_vars also large)
    for gene in genes_use:
        y = vector(adata, gene, layer=layer)
        a = y[m_group]
        b = y[m_ref]

        ma = float(np.nanmean(a))
        mb = float(np.nanmean(b))
        mean_a.append(ma)
        mean_b.append(mb)

        if method == "mwu":
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            eff = _rank_biserial_from_mwu(float(u), len(a), len(b))
        elif method == "ttest":
            t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
            # Cohen's d approx
            sa = np.nanstd(a, ddof=1)
            sb = np.nanstd(b, ddof=1)
            s = np.sqrt((sa * sa + sb * sb) / 2.0) if (sa > 0 or sb > 0) else np.nan
            eff = float((ma - mb) / s) if (s and np.isfinite(s) and s > 0) else np.nan
        else:
            # If user asked kruskal/anova for 2 groups: still OK, but use global stat
            if method == "kruskal":
                h, p = stats.kruskal(a, b, nan_policy="omit")
                # eta^2 from H (rough)
                eff = float((h - 1) / (len(y) - 1)) if len(y) > 1 else np.nan
            elif method == "anova":
                f, p = stats.f_oneway(a, b)
                eff = float((f) / (f + (len(y) - 2))) if len(y) > 2 else np.nan
            else:
                raise ValueError(f"Unsupported method='{method}'")

        pvals.append(float(p))
        effects.append(float(eff))

    pvals = np.asarray(pvals, dtype=float)
    qvals = _bh_fdr(pvals)

    df = pd.DataFrame(
        {
            "gene": genes_use,
            "pval": pvals,
            "qval": qvals,
            "effect_size": np.asarray(effects, dtype=float),
            "mean_group": np.asarray(mean_a, dtype=float),
            "mean_ref": np.asarray(mean_b, dtype=float),
        }
    )
    df["log2FC"] = [_log2fc(a, b) for a, b in zip(df["mean_group"], df["mean_ref"])]

    df = df.sort_values(["qval", "pval"], ascending=True).reset_index(drop=True)

    if store_key is not None:
        adata.uns.setdefault("assoc", {})
        adata.uns["assoc"][store_key] = {
            "params": {
                "groupby": groupby,
                "group": group,
                "reference": ref_label,
                "layer": layer,
                "method": method,
                "n_group": n1,
                "n_ref": n2,
            },
            "results": df,
        }
        info(f"Stored categorical gene association in adata.uns['assoc']['{store_key}'].")

    return df


def posthoc_per_gene(
    adata: ad.AnnData,
    *,
    genes: list[str],
    groupby: str,
    layer: str | None = None,
    method: Literal["mwu", "ttest"] = "mwu",
) -> dict[str, pd.DataFrame]:
    """
    For each gene, run pairwise posthoc tests across all categories in `groupby`.
    Returns dict: gene -> DataFrame (pairwise results).

    Requires bk.tl.pairwise_posthoc(df, method=...) to exist.
    """
    from .posthoc import pairwise_posthoc  # your existing helper

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    out: dict[str, pd.DataFrame] = {}
    grp = adata.obs[groupby].astype(str)

    for gene in genes:
        y = vector(adata, gene, layer=layer)
        tmp = pd.DataFrame({"y": y, "grp": grp})
        out[gene] = pairwise_posthoc(tmp, method=method)

    return out