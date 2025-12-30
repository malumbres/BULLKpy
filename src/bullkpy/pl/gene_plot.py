from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None

try:
    from scipy import stats
except Exception:  # pragma: no cover
    stats = None

import anndata as ad

from ._style import set_style, _savefig


def _get_vector(adata: ad.AnnData, gene: str, layer: str | None) -> np.ndarray:
    if gene not in adata.var_names:
        raise KeyError(f"gene '{gene}' not found in adata.var_names")
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    j = int(adata.var_names.get_loc(gene))
    v = X[:, j]
    if sp.issparse(v):
        v = np.asarray(v.toarray()).ravel()
    else:
        v = np.asarray(v).ravel()
    return v.astype(float)


def _format_p(p: float) -> str:
    if not np.isfinite(p):
        return "p = NA"
    if p < 1e-4:
        return f"p < 1e-4"
    if p < 1e-3:
        return f"p = {p:.1e}"
    return f"p = {p:.3g}"


def _bracket(ax: plt.Axes, x1: float, x2: float, y: float, h: float, text: str) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.0, c="0.2")
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom")


def gene_plot(
    adata: ad.AnnData,
    *,
    gene: str,
    groupby: str,
    layer: str | None = "log1p_cpm",
    groups: Sequence[str] | None = None,
    order: Sequence[str] | None = None,
    kind: Literal["violin", "box", "both"] = "both",
    show_points: bool = True,
    point_size: float = 3.0,
    point_alpha: float = 0.7,
    test: Literal["auto", "t-test", "wilcoxon", "anova", "kruskal"] = "auto",
    compare: Sequence[str] | None = None,   # e.g. ["LumA","LumB"] for pairwise
    add_pvalue: bool = True,
    show_n: bool = True,
    ylabel: str | None = None,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes, dict]:
    """
    Plot one gene across groups (bulk): violin/box + optional p-value.

    - If compare=['A','B'] -> pairwise test + bracket annotation.
    - Else -> overall test across all groups (ANOVA or Kruskal).

    Returns
    -------
    fig, ax, stats_dict
    """
    set_style()
    if sns is None:
        raise ImportError("gene_plot requires seaborn. Please install seaborn.")
    if stats is None:
        raise ImportError("gene_plot requires scipy. Please install scipy.")

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby '{groupby}' not found in adata.obs")

    y = _get_vector(adata, gene, layer)
    g = adata.obs[groupby].astype(str).copy()

    # subset groups
    mask = np.ones(adata.n_obs, dtype=bool)
    if groups is not None:
        groups = [str(x) for x in groups]
        mask = g.isin(groups).to_numpy()
    y = y[mask]
    g = g.iloc[mask]

    # ordering
    if order is None:
        cats = sorted(pd.unique(g))
    else:
        cats = [str(x) for x in order if str(x) in set(g)]
    g = pd.Categorical(g, categories=cats, ordered=True)

    df = pd.DataFrame({"group": g, "value": y}).dropna()
    if df.shape[0] == 0:
        raise ValueError("No data to plot after filtering/NA removal.")

    # figure
    if figsize is None:
        figsize = (max(5.0, 0.8 * len(cats) + 2.2), 4.2)
    fig, ax = plt.subplots(figsize=figsize)

    # main plot
    if kind in ("violin", "both"):
        sns.violinplot(
            data=df, x="group", y="value",
            inner=None, cut=0, linewidth=0.8,
            ax=ax,
        )
    if kind in ("box", "both"):
        sns.boxplot(
            data=df, x="group", y="value",
            width=0.35, showcaps=True,
            boxprops={"facecolor": "white", "edgecolor": "0.2", "linewidth": 1.0},
            whiskerprops={"color": "0.2", "linewidth": 1.0},
            medianprops={"color": "0.2", "linewidth": 1.0},
            showfliers=False,
            ax=ax,
        )

    if show_points:
        sns.stripplot(
            data=df, x="group", y="value",
            jitter=0.25, size=float(point_size),
            alpha=float(point_alpha),
            edgecolor="none",
            ax=ax,
        )

    ax.set_xlabel("")
    if ylabel is None:
        ylabel = f"{gene} expression" + (f" ({layer})" if layer is not None else "")
    ax.set_ylabel(ylabel)

    if title is None:
        title = gene
    ax.set_title(title)

    # add n per group
    if show_n:
        counts = df.groupby("group", observed=True)["value"].size().to_dict()
        xt = ax.get_xticks()
        for i, cat in enumerate(cats):
            if cat in counts:
                ax.text(xt[i], 0.98, f"n={counts[cat]}", transform=ax.get_xaxis_transform(),
                        ha="center", va="top", fontsize=max(7, plt.rcParams["font.size"] - 2), color="0.35")

    # statistics
    stats_out: dict = {"test": None, "pvalue": None, "compare": None}

    if add_pvalue:
        # group arrays
        arrays = [df.loc[df["group"] == c, "value"].to_numpy() for c in cats if np.any(df["group"] == c)]
        cats_present = [c for c in cats if np.any(df["group"] == c)]

        # choose test
        if compare is not None:
            a, b = map(str, compare)
            if a not in cats_present or b not in cats_present:
                raise ValueError(f"compare={compare} must be present in the plotted groups.")
            x1 = cats_present.index(a)
            x2 = cats_present.index(b)
            va = df.loc[df["group"] == a, "value"].to_numpy()
            vb = df.loc[df["group"] == b, "value"].to_numpy()

            if test == "auto":
                # default: nonparametric for bulk group comparisons
                test_use = "wilcoxon"
            else:
                test_use = test

            if test_use in ("t-test", "ttest"):
                p = stats.ttest_ind(va, vb, equal_var=False, nan_policy="omit").pvalue
                test_name = "t-test (Welch)"
            else:
                # Mann–Whitney U (wilcoxon-style)
                p = stats.mannwhitneyu(va, vb, alternative="two-sided").pvalue
                test_name = "Mann–Whitney U"

            stats_out.update({"test": test_name, "pvalue": float(p), "compare": [a, b]})

            # bracket
            ymax = np.nanmax(df["value"].to_numpy())
            yr = np.nanmax(df["value"].to_numpy()) - np.nanmin(df["value"].to_numpy())
            yr = yr if np.isfinite(yr) and yr > 0 else 1.0
            yb = ymax + 0.08 * yr
            _bracket(ax, x1, x2, yb, 0.03 * yr, _format_p(float(p)))

        else:
            # overall test across groups
            if len(arrays) < 2:
                warn_txt = "Not enough groups for p-value."
                stats_out.update({"test": warn_txt, "pvalue": np.nan})
            else:
                if test == "auto":
                    test_use = "kruskal"
                else:
                    test_use = test

                if test_use == "anova":
                    p = stats.f_oneway(*arrays).pvalue
                    test_name = "ANOVA"
                else:
                    p = stats.kruskal(*arrays).pvalue
                    test_name = "Kruskal–Wallis"

                stats_out.update({"test": test_name, "pvalue": float(p)})

                ax.text(
                    0.99, 0.02,
                    f"{test_name}\n{_format_p(float(p))}",
                    transform=ax.transAxes,
                    ha="right", va="bottom",
                    color="0.25",
                )

    plt.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax, stats_out