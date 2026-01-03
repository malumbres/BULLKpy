from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import anndata as ad
import seaborn as sns

from ._style import set_style, _savefig
from ..logging import info, warn


def qc_metrics(
    adata: ad.AnnData,
    *,
    color: str | None = "pct_counts_mt",
    vars_to_plot: Sequence[str] = (
        "total_counts",
        "n_genes_detected",
        "pct_counts_mt",
        "pct_counts_ribo",
    ),
    log1p_total_counts: bool = True,
    log1p_n_genes: bool = False,
    point_size: float = 20.0,
    alpha: float = 0.8,
    figsize: tuple[float, float] = (10, 7),
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, np.ndarray]:
    """
    Plot bulk RNA-seq QC metrics (robust to missing columns).

    - If total_counts + n_genes_detected exist: scatter (library size vs detected genes)
    - Otherwise: skip scatter and show histograms only
    """
    set_style()

    # --- vars that actually exist ---
    vars_use = [v for v in vars_to_plot if v in adata.obs.columns]
    missing = [v for v in vars_to_plot if v not in adata.obs.columns]
    if missing:
        warn(f"qc_metrics: skipping missing QC columns: {missing}")

    if len(vars_use) == 0:
        raise KeyError(
            "qc_metrics: none of vars_to_plot exist in adata.obs. "
            "Run bk.pp.qc_metrics first (or compute n_genes_detected manually)."
        )

    # --- color handling (optional) ---
    color_vals = None
    if color is not None:
        if color not in adata.obs.columns:
            warn(f"qc_metrics: color='{color}' not found in adata.obs; coloring disabled.")
        else:
            s = adata.obs[color]
            # If categorical, we'll just pass codes (matplotlib handles numeric colors)
            if str(s.dtype) == "category" or s.dtype == object:
                color_vals = pd.Categorical(s.astype(str)).codes
            else:
                color_vals = s.to_numpy()

    # --- do we have scatter ingredients? ---
    has_total = "total_counts" in adata.obs.columns
    has_ngenes = "n_genes_detected" in adata.obs.columns
    do_scatter = has_total and has_ngenes

    # --- layout ---
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    if do_scatter:
        gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1.0])
        ax_scatter = fig.add_subplot(gs[:, 0])
        ax_h1 = fig.add_subplot(gs[0, 1])
        ax_h2 = fig.add_subplot(gs[1, 1])

        # Scatter data
        x = adata.obs["total_counts"].to_numpy(dtype=float)
        y = adata.obs["n_genes_detected"].to_numpy(dtype=float)
        x_plot = np.log1p(x) if log1p_total_counts else x
        y_plot = np.log1p(y) if log1p_n_genes else y

        if color_vals is None:
            ax_scatter.scatter(x_plot, y_plot, s=point_size, alpha=alpha, edgecolors="none")
        else:
            sc = ax_scatter.scatter(
                x_plot, y_plot, s=point_size, alpha=alpha, c=color_vals, edgecolors="none"
            )
            cbar = fig.colorbar(sc, ax=ax_scatter, pad=0.02, fraction=0.05)
            cbar.set_label(color if color is not None else "")

        ax_scatter.set_xlabel("log1p(total_counts)" if log1p_total_counts else "total_counts")
        ax_scatter.set_ylabel("log1p(n_genes_detected)" if log1p_n_genes else "n_genes_detected")
        ax_scatter.set_title("Library size vs detected genes")

    else:
        # hist-only layout
        warn("qc_metrics: 'total_counts' or 'n_genes_detected' missing; scatter skipped (histograms only).")
        gs = fig.add_gridspec(2, 2)
        ax_scatter = None
        ax_h1 = fig.add_subplot(gs[0, :])
        ax_h2 = fig.add_subplot(gs[1, :])

    # --- histograms ---
    # first histogram
    v1 = vars_use[0]
    ax_h1.hist(adata.obs[v1].to_numpy(), bins=30, alpha=0.85)
    ax_h1.set_title(v1)
    ax_h1.set_ylabel("n samples")

    # second histogram
    if len(vars_use) >= 2:
        v2 = vars_use[1]
        ax_h2.hist(adata.obs[v2].to_numpy(), bins=30, alpha=0.85)
        ax_h2.set_title(v2)
        ax_h2.set_ylabel("n samples")
    else:
        ax_h2.axis("off")

    # optional overlays (3rd and 4th) on twin axes
    if len(vars_use) >= 3:
        v3 = vars_use[2]
        ax_h1_t = ax_h1.twinx()
        ax_h1_t.hist(adata.obs[v3].to_numpy(), bins=30, alpha=0.45)
        ax_h1_t.set_ylabel("n samples (overlay)")
        ax_h1.set_title(f"{v1} + {v3}")

    if len(vars_use) >= 4 and ax_h2 is not None and ax_h2.axes is not None:
        v4 = vars_use[3]
        ax_h2_t = ax_h2.twinx()
        ax_h2_t.hist(adata.obs[v4].to_numpy(), bins=30, alpha=0.45)
        ax_h2_t.set_ylabel("n samples (overlay)")
        ax_h2.set_title(f"{v2} + {v4}")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    axes = [a for a in [ax_scatter, ax_h1, ax_h2] if a is not None]
    return fig, np.array(axes, dtype=object)


# -----------------------------------------------------------------------------
# Scanpy-like QC scatter plots with thresholds (bulk)
# -----------------------------------------------------------------------------

def _qc_mask(
    df: pd.DataFrame,
    *,
    x: str,
    y: str,
    min_x: float | None = None,
    max_x: float | None = None,
    min_y: float | None = None,
    max_y: float | None = None,
) -> np.ndarray:
    ok = np.ones(len(df), dtype=bool)
    if min_x is not None:
        ok &= df[x].to_numpy() >= float(min_x)
    if max_x is not None:
        ok &= df[x].to_numpy() <= float(max_x)
    if min_y is not None:
        ok &= df[y].to_numpy() >= float(min_y)
    if max_y is not None:
        ok &= df[y].to_numpy() <= float(max_y)
    return ok


def _scatter_qc(
    adata: ad.AnnData,
    *,
    ax: plt.Axes,
    x: str,
    y: str,
    groupby: str | None = None,
    min_x: float | None = None,
    max_x: float | None = None,
    min_y: float | None = None,
    max_y: float | None = None,
    logx: bool = False,
    logy: bool = False,
    s: float = 18.0,
    alpha: float = 0.85,
    linewidth: float = 0.25,
    edgecolor: str = "0.15",
    show_outliers: bool = True,
    outlier_color: str = "crimson",
    outlier_marker: str = "x",
    outlier_size: float = 28.0,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    legend: bool = True,
    legend_loc: str = "best",
) -> tuple[np.ndarray, pd.DataFrame]:
    for k in (x, y):
        if k not in adata.obs.columns:
            raise KeyError(f"adata.obs['{k}'] not found")

    df = adata.obs[[x, y]].copy()
    df[x] = pd.to_numeric(df[x], errors="coerce")
    df[y] = pd.to_numeric(df[y], errors="coerce")
    df = df.dropna(subset=[x, y])

    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"adata.obs['{groupby}'] not found")
        df[groupby] = adata.obs.loc[df.index, groupby].astype(str)

    ok = _qc_mask(df, x=x, y=y, min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y)

    if groupby is None:
        ax.scatter(
            df.loc[ok, x],
            df.loc[ok, y],
            s=s,
            alpha=alpha,
            linewidths=linewidth,
            edgecolors=edgecolor,
        )
    else:
        cats = pd.Categorical(df[groupby]).categories
        for c in cats:
            m = ok & (df[groupby] == c).to_numpy()
            if m.sum() == 0:
                continue
            ax.scatter(
                df.loc[m, x],
                df.loc[m, y],
                s=s,
                alpha=alpha,
                linewidths=linewidth,
                edgecolors=edgecolor,
                label=str(c),
            )

    if show_outliers and (~ok).any():
        ax.scatter(
            df.loc[~ok, x],
            df.loc[~ok, y],
            s=outlier_size,
            alpha=0.95,
            marker=outlier_marker,
            c=outlier_color,
            linewidths=0.8,
            label="QC fail" if groupby is not None else None,
        )

    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    def _vline(val: float) -> None:
        ax.axvline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    def _hline(val: float) -> None:
        ax.axhline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    if min_x is not None:
        _vline(float(min_x))
    if max_x is not None:
        _vline(float(max_x))
    if min_y is not None:
        _hline(float(min_y))
    if max_y is not None:
        _hline(float(max_y))

    ax.set_xlabel(xlabel if xlabel is not None else x)
    ax.set_ylabel(ylabel if ylabel is not None else y)
    if title is not None:
        ax.set_title(title)

    if legend and groupby is not None:
        ax.legend(loc=legend_loc, frameon=False)

    return ok, df


def library_size_vs_genes(
    adata: ad.AnnData,
    *,
    x: str = "total_counts",
    y: str = "n_genes_detected",
    groupby: str | None = None,
    # thresholds
    min_counts: float | None = None,
    max_counts: float | None = None,
    min_genes: float | None = None,
    max_genes: float | None = None,
    # display
    logx: bool = True,
    logy: bool = True,
    s: float = 18.0,
    alpha: float = 0.85,
    linewidth: float = 0.25,
    edgecolor: str = "0.15",
    # highlight outliers
    show_outliers: bool = True,
    outlier_color: str = "crimson",
    outlier_marker: str = "x",
    outlier_size: float = 28.0,
    # labels
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    # figure
    figsize: tuple[float, float] = (5.5, 4.5),
    legend: bool = True,
    legend_loc: str = "best",
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    QC scatter: library size vs detected genes with optional thresholds.

    Requires columns in adata.obs:
      - x (default: total_counts)
      - y (default: n_genes_detected)

    Typical usage:
      bk.pl.library_size_vs_genes(adata, min_counts=1e6, min_genes=12000, groupby="Subtype")
    """
    set_style()

    if x not in adata.obs.columns:
        raise KeyError(f"adata.obs['{x}'] not found")
    if y not in adata.obs.columns:
        raise KeyError(f"adata.obs['{y}'] not found")

    df = adata.obs[[x, y]].copy()
    df[x] = pd.to_numeric(df[x], errors="coerce")
    df[y] = pd.to_numeric(df[y], errors="coerce")
    df = df.dropna(subset=[x, y])

    # compute pass/fail mask
    ok = np.ones(len(df), dtype=bool)

    if min_counts is not None:
        ok &= df[x].to_numpy() >= float(min_counts)
    if max_counts is not None:
        ok &= df[x].to_numpy() <= float(max_counts)
    if min_genes is not None:
        ok &= df[y].to_numpy() >= float(min_genes)
    if max_genes is not None:
        ok &= df[y].to_numpy() <= float(max_genes)

    # handle grouping/palette via seaborn set_palette in set_style()
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise KeyError(f"adata.obs['{groupby}'] not found")
        df[groupby] = adata.obs.loc[df.index, groupby].astype(str)

    fig, ax = plt.subplots(figsize=figsize)

    # scatter (inliers)
    if groupby is None:
        ax.scatter(
            df.loc[ok, x],
            df.loc[ok, y],
            s=s,
            alpha=alpha,
            linewidths=linewidth,
            edgecolors=edgecolor,
        )
    else:
        # plot per category (keeps legend manageable and stable)
        cats = pd.Categorical(df[groupby]).categories
        for c in cats:
            m = ok & (df[groupby] == c).to_numpy()
            if m.sum() == 0:
                continue
            ax.scatter(
                df.loc[m, x],
                df.loc[m, y],
                s=s,
                alpha=alpha,
                linewidths=linewidth,
                edgecolors=edgecolor,
                label=str(c),
            )

    # outliers on top
    if show_outliers and (~ok).any():
        ax.scatter(
            df.loc[~ok, x],
            df.loc[~ok, y],
            s=outlier_size,
            alpha=0.95,
            marker=outlier_marker,
            c=outlier_color,
            linewidths=0.8,
            label="QC fail" if groupby is not None else None,
        )

    # log scales
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    # threshold lines
    def _vline(val):
        ax.axvline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    def _hline(val):
        ax.axhline(val, color="0.35", lw=1.0, ls="--", zorder=0)

    if min_counts is not None:
        _vline(min_counts)
    if max_counts is not None:
        _vline(max_counts)
    if min_genes is not None:
        _hline(min_genes)
    if max_genes is not None:
        _hline(max_genes)

    # labels
    ax.set_xlabel(xlabel if xlabel is not None else x)
    ax.set_ylabel(ylabel if ylabel is not None else y)

    if title is None:
        n_fail = int((~ok).sum())
        title = f"{y} vs {x} (QC fail: {n_fail})" if (min_counts or max_counts or min_genes or max_genes) else f"{y} vs {x}"
    ax.set_title(title)

    if legend and groupby is not None:
        ax.legend(loc=legend_loc, frameon=False, ncol=1)

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax

def mt_fraction_vs_counts(
    adata: ad.AnnData,
    *,
    x: str = "total_counts",
    y: str = "pct_counts_mt",
    groupby: str | None = None,
    min_counts: float | None = None,
    max_counts: float | None = None,
    min_mt: float | None = None,
    max_mt: float | None = None,
    logx: bool = True,
    logy: bool = False,
    figsize: tuple[float, float] = (5.5, 4.5),
    save: str | Path | None = None,
    show: bool = True,
):
    """Scatter QC: total counts vs mitochondrial fraction (with thresholds)."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)
    ok, _ = _scatter_qc(
        adata,
        ax=ax,
        x=x,
        y=y,
        groupby=groupby,
        min_x=min_counts,
        max_x=max_counts,
        min_y=min_mt,
        max_y=max_mt,
        logx=logx,
        logy=logy,
        title=f"{y} vs {x} (QC fail: {(~ok).sum()})"
        if any(v is not None for v in (min_counts, max_counts, min_mt, max_mt))
        else f"{y} vs {x}",
    )
    fig.tight_layout()
    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


def genes_vs_mt_fraction(
    adata: ad.AnnData,
    *,
    x: str = "pct_counts_mt",
    y: str = "n_genes_detected",
    groupby: str | None = None,
    min_mt: float | None = None,
    max_mt: float | None = None,
    min_genes: float | None = None,
    max_genes: float | None = None,
    logx: bool = False,
    logy: bool = True,
    figsize: tuple[float, float] = (5.5, 4.5),
    save: str | Path | None = None,
    show: bool = True,
):
    """Scatter QC: mt fraction vs detected genes (with thresholds)."""
    set_style()
    fig, ax = plt.subplots(figsize=figsize)
    ok, _ = _scatter_qc(
        adata,
        ax=ax,
        x=x,
        y=y,
        groupby=groupby,
        min_x=min_mt,
        max_x=max_mt,
        min_y=min_genes,
        max_y=max_genes,
        logx=logx,
        logy=logy,
        title=f"{y} vs {x} (QC fail: {(~ok).sum()})"
        if any(v is not None for v in (min_mt, max_mt, min_genes, max_genes))
        else f"{y} vs {x}",
    )
    fig.tight_layout()
    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


def qc_scatter_panel(
    adata: ad.AnnData,
    *,
    groupby: str | None = None,
    min_counts: float | None = None,
    max_counts: float | None = None,
    min_genes: float | None = None,
    max_genes: float | None = None,
    min_mt: float | None = None,
    max_mt: float | None = None,
    total_counts_key: str = "total_counts",
    n_genes_key: str = "n_genes_detected",
    pct_mt_key: str = "pct_counts_mt",
    log_counts: bool = True,
    log_genes: bool = True,
    figsize: tuple[float, float] = (16.0, 4.6),
    save: str | Path | None = None,
    show: bool = True,
    share_legend: bool = True,
):
    """1×3 QC panel: (genes vs counts) | (mt% vs counts) | (genes vs mt%)."""
    set_style()
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    _scatter_qc(
        adata,
        ax=axes[0],
        x=total_counts_key,
        y=n_genes_key,
        groupby=groupby,
        min_x=min_counts,
        max_x=max_counts,
        min_y=min_genes,
        max_y=max_genes,
        logx=log_counts,
        logy=log_genes,
        legend=(groupby is not None and not share_legend),
        title=f"{n_genes_key} vs {total_counts_key}",
    )

    _scatter_qc(
        adata,
        ax=axes[1],
        x=total_counts_key,
        y=pct_mt_key,
        groupby=groupby,
        min_x=min_counts,
        max_x=max_counts,
        min_y=min_mt,
        max_y=max_mt,
        logx=log_counts,
        logy=False,
        legend=False,
        title=f"{pct_mt_key} vs {total_counts_key}",
    )

    _scatter_qc(
        adata,
        ax=axes[2],
        x=pct_mt_key,
        y=n_genes_key,
        groupby=groupby,
        min_x=min_mt,
        max_x=max_mt,
        min_y=min_genes,
        max_y=max_genes,
        logx=False,
        logy=log_genes,
        legend=False,
        title=f"{n_genes_key} vs {pct_mt_key}",
    )

    if share_legend and groupby is not None:
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="center right", frameon=False)
            fig.subplots_adjust(right=0.90)

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes

def qc_by_group(
    adata: ad.AnnData,
    *,
    groupby: str,
    keys: Sequence[str] = ("total_counts", "n_genes_detected", "pct_counts_mt", "pct_counts_ribo"),
    kind: Literal["violin", "box"] = "violin",
    log1p: Sequence[str] = ("total_counts",),
    figsize: tuple[float, float] = (11, 4),
    rotate_xticks: int = 45,
    save: str | Path | None = None,
    show: bool = True,
    show_n: bool = True,
):
    """
    Plot QC metrics grouped by a metadata column in `adata.obs` (e.g. batch, cohort).
    """
    set_style()

    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    for k in keys:
        if k not in adata.obs.columns:
            raise KeyError(f"Missing '{k}' in adata.obs. Run bk.pp.qc_metrics(adata) first.")

    groups = adata.obs[groupby].astype("category")
    cat = groups.cat.categories.tolist()

    counts = groups.value_counts().reindex(cat)

    fig, axes = plt.subplots(1, len(keys), figsize=figsize, constrained_layout=True)
    if len(keys) == 1:
        axes = [axes]

    for ax, k in zip(axes, keys):
        vals = adata.obs[k].to_numpy(dtype=float)
        if k in log1p:
            vals = np.log1p(vals)

        data = [vals[groups == g] for g in cat]

        if kind == "violin":
            parts = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
            # keep default coloring; just clean edges
            for pc in parts.get("bodies", []):
                pc.set_alpha(0.8)
        else:
            ax.boxplot(data, showfliers=False)

        ax.set_title(f"{_label(k, log1p)}")

        ax.set_xticks(range(1, len(cat) + 1))

        if show_n:
            labels = [f"{c} (n={counts[c]})" for c in cat]
        else:
            labels = cat

        ax.set_xticklabels(labels, rotation=rotate_xticks, ha="right")


        ax.set_ylabel("value")

    fig.suptitle(f"QC by {groupby}", y=1.05)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes


def _label(k: str, log1p: Sequence[str]) -> str:
    return f"log1p({k})" if k in log1p else k


def qc_pairplot(
    adata: ad.AnnData,
    *,
    keys: Sequence[str] = ("total_counts", "n_genes_detected", "pct_counts_mt"),
    color: str | None = "pct_counts_mt",
    log1p: Sequence[str] = ("total_counts",),
    point_size: float = 14.0,
    alpha: float = 0.7,
    figsize: tuple[float, float] = (8, 8),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scatter-matrix (pairplot) of QC metrics in `adata.obs`.

    Diagonal: histograms
    Off-diagonal: scatter plots
    """
    set_style()

    for k in keys:
        if k not in adata.obs.columns:
            raise KeyError(f"Missing '{k}' in adata.obs. Run bk.pp.qc_metrics(adata) first.")

    C = None
    if color is not None:
        if color in adata.obs.columns:
            C = adata.obs[color].to_numpy()
        else:
            warn(f"color='{color}' not found in adata.obs; coloring disabled.")

    # prepare transformed columns
    data = {}
    for k in keys:
        v = adata.obs[k].to_numpy(dtype=float)
        if k in log1p:
            v = np.log1p(v)
        data[k] = v

    n = len(keys)
    fig, axes = plt.subplots(n, n, figsize=figsize, constrained_layout=True)

    for i, yi in enumerate(keys):
        for j, xj in enumerate(keys):
            ax = axes[i, j]
            if i == j:
                ax.hist(data[xj], bins=30)
                ax.set_ylabel("")
            else:
                sc = ax.scatter(
                    data[xj], data[yi],
                    c=C, s=point_size, alpha=alpha, edgecolors="none"
                )
            # labels only on left and bottom
            if i < n - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel(_label(xj, log1p))
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(_label(yi, log1p))

    if C is not None:
        cbar = fig.colorbar(sc, ax=axes, shrink=0.75, pad=0.01)
        cbar.set_label(color)

    fig.suptitle("QC pairplot", y=1.02)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, axes


def _label(k: str, log1p: Sequence[str]) -> str:
    return f"log1p({k})" if k in log1p else k
