from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import anndata as ad
from scipy import stats

from ..logging import warn
from ._style import set_style, _savefig


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    p0 = p[ok]
    n = p0.size
    order = np.argsort(p0)
    ranked = p0[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out_ok = np.empty_like(p0)
    out_ok[order] = q
    out[ok] = out_ok
    return out


def gene_association_volcano(
    adata: ad.AnnData,
    *,
    groupby: str,
    group: str,
    reference: str | None = None,   # if None -> "rest"
    layer: str | None = "log1p_cpm",
    genes: Sequence[str] | None = None,
    method: Literal["mwu"] = "mwu",
    effect: Literal["delta_mean", "delta_median"] = "delta_mean",
    alpha: float = 0.05,
    top_labels: int = 12,
    figsize: tuple[float, float] = (6.5, 5.5),
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Effect-size aware volcano for categorical association:
      x = effect size (delta_mean or delta_median between group vs reference/rest)
      y = -log10(q)

    This is NOT a full DE method — it's a fast association scan.
    """
    set_style()
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    grp = adata.obs[groupby].astype(str)
    m1 = grp.eq(str(group)).to_numpy()
    if reference is None:
        m2 = ~m1
        ref_name = "rest"
    else:
        m2 = grp.eq(str(reference)).to_numpy()
        ref_name = str(reference)

    if m1.sum() < 2 or m2.sum() < 2:
        raise ValueError(f"Not enough samples in '{group}' and '{ref_name}' for volcano.")

    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    if genes is None:
        genes_use = list(map(str, adata.var_names))
    else:
        genes_use = [str(g) for g in genes]
        miss = [g for g in genes_use if g not in adata.var_names]
        if miss:
            raise KeyError(f"Genes not in adata.var_names (first 10): {miss[:10]}")

    idx = [adata.var_names.get_loc(g) for g in genes_use]
    M = X[:, idx]
    if sp.issparse(M):
        M = M.toarray()
    else:
        M = np.asarray(M, dtype=float)

    pvals = np.full(len(genes_use), np.nan, dtype=float)
    effs = np.full(len(genes_use), np.nan, dtype=float)

    for k in range(len(genes_use)):
        x1 = M[m1, k]
        x2 = M[m2, k]
        # MWU (two-sided)
        try:
            _, p = stats.mannwhitneyu(x1, x2, alternative="two-sided")
        except Exception:
            p = np.nan
        pvals[k] = p

        if effect == "delta_mean":
            effs[k] = np.nanmean(x1) - np.nanmean(x2)
        else:
            effs[k] = np.nanmedian(x1) - np.nanmedian(x2)

    qvals = _bh_fdr(pvals)
    df = pd.DataFrame({"gene": genes_use, "effect": effs, "pval": pvals, "qval": qvals})
    df["neglog10q"] = -np.log10(df["qval"].clip(lower=1e-300))
    df = df.sort_values("qval", na_position="last").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.scatter(df["effect"], df["neglog10q"], s=10, alpha=0.75, edgecolors="none")

    # threshold line
    ax.axhline(-np.log10(alpha), lw=1.0, color="0.3")

    ax.set_xlabel(effect)
    ax.set_ylabel("-log10(q)")
    ax.set_title(title or f"{groupby}: {group} vs {ref_name}")

    # label top genes (by q and magnitude)
    lab = df[np.isfinite(df["qval"])].head(max(top_labels * 3, top_labels))
    lab = lab.reindex(lab["effect"].abs().sort_values(ascending=False).index).head(top_labels)
    for _, r in lab.iterrows():
        ax.text(r["effect"], r["neglog10q"], str(r["gene"]), fontsize=9)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax, df