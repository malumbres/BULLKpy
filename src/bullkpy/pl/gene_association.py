from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
from scipy import stats

from ..logging import warn
from ._style import set_style, _savefig
from ..tl.posthoc import pairwise_posthoc


def _get_gene_vector(adata: ad.AnnData, gene: str, *, layer: str | None = "log1p_cpm") -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"Gene '{gene}' not in adata.var_names")
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    idx = adata.var_names.get_loc(gene)
    v = X[:, idx]
    if sp.issparse(v):
        v = v.toarray().ravel()
    else:
        v = np.asarray(v).ravel()
    return v.astype(float)


def _stars(q: float) -> str:
    if not np.isfinite(q):
        return ""
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    if q < 5e-2:
        return "*"
    return ""


def _add_brackets(
    ax: plt.Axes,
    post: pd.DataFrame,
    *,
    order: list[str],
    alpha: float = 0.05,
    max_brackets: int = 6,
    bracket_height: float = 0.05,  # fraction of y-range
    lw: float = 1.0,
):
    # choose significant pairs
    sig = post[np.isfinite(post["qval"]) & (post["qval"] <= alpha)].copy()
    if sig.shape[0] == 0:
        return
    sig = sig.sort_values(["qval", "pval"], na_position="last").head(max_brackets)

    y0, y1 = ax.get_ylim()
    yr = (y1 - y0) if (y1 - y0) != 0 else 1.0
    step = bracket_height * yr

    # start above current max
    ymax = y1
    used = 0
    for _, r in sig.iterrows():
        g1, g2 = str(r["group1"]), str(r["group2"])
        if g1 not in order or g2 not in order:
            continue
        i1, i2 = order.index(g1), order.index(g2)
        if i1 == i2:
            continue
        x1, x2 = (min(i1, i2), max(i1, i2))
        y = ymax + step * (used + 1)

        ax.plot([x1, x1, x2, x2], [y, y + 0.3 * step, y + 0.3 * step, y], lw=lw, color="0.2", clip_on=False)
        ax.text((x1 + x2) / 2, y + 0.35 * step, _stars(float(r["qval"])), ha="center", va="bottom", color="0.2")
        used += 1

    ax.set_ylim(y0, ymax + step * (used + 2))


def gene_association(
    adata: ad.AnnData,
    *,
    gene: str | Sequence[str],
    groupby: str,
    layer: str | None = "log1p_cpm",
    kind: Literal["violin", "box"] = "violin",
    order: Sequence[str] | None = None,
    rotate_xticklabels: int = 45,
    figsize: tuple[float, float] | None = None,
    panel_size: tuple[float, float] = (4.2, 3.2),
    show_points: bool = True,
    point_size: float = 2.0,
    point_alpha: float = 0.35,
    palette: str = "Set2",
    annotate_posthoc: bool = True,
    posthoc_method: Literal["mwu", "ttest"] = "mwu",
    posthoc_alpha: float = 0.05,
    max_brackets: int = 6,
    bracket_height: float = 0.06,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Gene vs categorical obs association plot (Scanpy-like panels).

    - gene can be a string or list of genes -> row of panels
    - violin/box + optional strip points
    - optional automatic pairwise post-hoc + significance brackets (BH corrected)

    Returns (fig, axes).
    """
    set_style()
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not in adata.obs")

    genes = [gene] if isinstance(gene, str) else [str(g) for g in gene]
    grp = adata.obs[groupby].astype(str)

    if order is None:
        cats = list(pd.Categorical(grp).categories)
    else:
        cats = [str(x) for x in order]

    # build plot df once per gene (keeps memory low)
    n = len(genes)
    if figsize is None:
        figsize = (panel_size[0] * n, panel_size[1])

    fig, axes = plt.subplots(1, n, figsize=figsize, constrained_layout=True)
    if n == 1:
        axes = [axes]

    for ax, g in zip(axes, genes):
        y = _get_gene_vector(adata, g, layer=layer)
        df = pd.DataFrame({"y": y, "grp": grp.values})
        df["grp"] = pd.Categorical(df["grp"], categories=cats, ordered=True)

        if kind == "violin":
            sns.violinplot(data=df, x="grp", y="y", ax=ax, cut=0, inner="quartile", palette=palette)
        else:
            sns.boxplot(data=df, x="grp", y="y", ax=ax, palette=palette)

        if show_points:
            sns.stripplot(
                data=df, x="grp", y="y", ax=ax,
                color="k", size=float(point_size), alpha=float(point_alpha), jitter=0.25
            )

        ax.set_title(g)
        ax.set_xlabel("")
        ax.set_ylabel("Expression" if layer is None else str(layer))
        ax.tick_params(axis="x", rotation=int(rotate_xticklabels))

        # posthoc
        if annotate_posthoc and len(cats) >= 2:
            try:
                post = pairwise_posthoc(df, group_col="grp", value_col="y", method=posthoc_method, correction="bh")
                _add_brackets(
                    ax, post,
                    order=cats,
                    alpha=float(posthoc_alpha),
                    max_brackets=int(max_brackets),
                    bracket_height=float(bracket_height),
                )
            except Exception as e:
                warn(f"posthoc annotation failed for gene '{g}': {e}")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, np.array(axes, dtype=object)

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