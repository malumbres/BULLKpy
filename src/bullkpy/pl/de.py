from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

from ..logging import warn
from ._style import set_style, _savefig


def _get_de_table(
    adata: ad.AnnData | None,
    *,
    res: pd.DataFrame | None = None,
    uns_key: str | None = None,
    contrast: str | None = None,
) -> pd.DataFrame:
    if res is not None:
        return res
    if adata is None:
        raise ValueError("Provide either `res=` or `adata=` with `uns_key`/`contrast`.")
    if uns_key is None:
        uns_key = "de"
    if contrast is None:
        raise ValueError("When using `adata`, you must provide `contrast=`.")
    return adata.uns[uns_key][contrast]["results"]



def volcano(
    res: pd.DataFrame,
    *,
    gene_col: str = "gene",
    fc_col: str = "log2FC",
    p_col: str = "pval",
    q_col: str = "qval",
    use_qval: bool = True,
    p_cutoff: float = 0.05,
    fc_cutoff: float = 0.0,  # optional, set >0 to require |log2FC| >= fc_cutoff
    title: str | None = None,
    figsize: tuple[float, float] = (6.5, 5.5),
    alpha: float = 0.85,
    point_size: float = 18.0,
    # colors
    color_ns: str = "#BDBDBD",   # grey
    color_up: str = "#D62728",   # red
    color_down: str = "#1F77B4", # blue
    # labels
    label_genes: Sequence[str] = (),
    top_n_labels: int = 10,
    bottom_n_labels: int = 10,
    label_fontsize: float | None = None,
    label_offset: tuple[float, float] = (0.02, 0.02),  # in data units (x,y); small nudge
    # thresholds lines
    show_thresholds: bool = True,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Volcano plot for DE results.

    Coloring:
      - non-significant: grey
      - significant & log2FC > 0: red
      - significant & log2FC < 0: blue

    Labeling:
      - label_genes: explicit list
      - top_n_labels: top upregulated significant genes
      - bottom_n_labels: top downregulated significant genes (most negative log2FC)
    """
    set_style()

    df = res.copy()
    for c in (gene_col, fc_col, p_col):
        if c not in df.columns:
            raise KeyError(f"'{c}' column not found in results.")
    if use_qval and (q_col not in df.columns):
        raise KeyError(f"use_qval=True but '{q_col}' column not found in results.")

    # numeric coercion
    df[fc_col] = pd.to_numeric(df[fc_col], errors="coerce")
    p_used_col = q_col if use_qval else p_col
    df[p_used_col] = pd.to_numeric(df[p_used_col], errors="coerce")

    # drop NA for plotting
    df = df.dropna(subset=[fc_col, p_used_col, gene_col]).copy()

    # y = -log10(p)
    pvals = df[p_used_col].to_numpy(dtype=float)
    pvals = np.clip(pvals, 1e-300, 1.0)
    df["_neglog10p"] = -np.log10(pvals)

    # significance mask
    sig = df[p_used_col] <= float(p_cutoff)
    if fc_cutoff > 0:
        sig = sig & (df[fc_col].abs() >= float(fc_cutoff))

    up = sig & (df[fc_col] > 0)
    down = sig & (df[fc_col] < 0)
    ns = ~sig

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # plot in layers: ns first, then down, then up (so sig points are on top)
    ax.scatter(
        df.loc[ns, fc_col],
        df.loc[ns, "_neglog10p"],
        s=point_size,
        alpha=alpha,
        c=color_ns,
        edgecolors="none",
        label="NS",
        rasterized=True,
    )
    ax.scatter(
        df.loc[down, fc_col],
        df.loc[down, "_neglog10p"],
        s=point_size,
        alpha=alpha,
        c=color_down,
        edgecolors="none",
        label="Down",
        rasterized=True,
    )
    ax.scatter(
        df.loc[up, fc_col],
        df.loc[up, "_neglog10p"],
        s=point_size,
        alpha=alpha,
        c=color_up,
        edgecolors="none",
        label="Up",
        rasterized=True,
    )

    # threshold lines
    if show_thresholds:
        ax.axhline(-np.log10(max(p_cutoff, 1e-300)), lw=1, ls="--", color="0.4")
        if fc_cutoff > 0:
            ax.axvline(+fc_cutoff, lw=1, ls="--", color="0.4")
            ax.axvline(-fc_cutoff, lw=1, ls="--", color="0.4")

    ax.set_xlabel("log2 fold change")
    ax.set_ylabel(f"-log10({'qval' if use_qval else 'pval'})")
    if title is not None:
        ax.set_title(title)

    # ---------- labels ----------
    label_fontsize = label_fontsize if label_fontsize is not None else 10.0
    labels = set(str(g) for g in label_genes)

    # top/bottom labels among significant only
    if top_n_labels and top_n_labels > 0:
        top_up = (
            df.loc[up]
            .sort_values(fc_col, ascending=False)
            .head(int(top_n_labels))
        )
        labels.update(top_up[gene_col].astype(str).tolist())

    if bottom_n_labels and bottom_n_labels > 0:
        top_down = (
            df.loc[down]
            .sort_values(fc_col, ascending=True)  # most negative first
            .head(int(bottom_n_labels))
        )
        labels.update(top_down[gene_col].astype(str).tolist())

    # annotate (only if present in df)
    if labels:
        df_idx = df.set_index(df[gene_col].astype(str), drop=False)
        for g in labels:
            if g not in df_idx.index:
                continue
            row = df_idx.loc[g]
            # if duplicates exist, pick the first
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            x = float(row[fc_col])
            y = float(row["_neglog10p"])
            ax.text(
                x + label_offset[0],
                y + label_offset[1],
                str(g),
                fontsize=label_fontsize,
                ha="left",
                va="bottom",
            )

    ax.legend(frameon=False, loc="upper right")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax



def rankplot(
    adata: ad.AnnData | None = None,
    *,
    res: pd.DataFrame | None = None,
    uns_key: str | None = None,
    contrast: str | None = None,
    n_genes: int = 20,
    sort_by: str = "qval",           # "qval" | "pval" | "log2FC"
    direction: str = "both",         # "up" | "down" | "both"
    fc_col: str | None = None,       # auto: "log2FC" else "logFC"
    figsize: tuple[float, float] = (7, 6),
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Ranked horizontal barplot of top DE genes (up/down).

    - Upregulated bars: red
    - Downregulated bars: blue
    - sort_by supports: qval, pval, log2FC
    - For direction="both": top shows strongest upregulated; bottom shows strongest downregulated (last).
    """
    set_style()
    df = _get_de_table(adata, res=res, uns_key=uns_key, contrast=contrast).copy()

    gene_col = "gene" if "gene" in df.columns else df.columns[0]

    if fc_col is None:
        fc_col = "log2FC" if "log2FC" in df.columns else ("logFC" if "logFC" in df.columns else None)
    if fc_col is None:
        raise KeyError("No fold-change column found. Provide fc_col='log2FC' or 'logFC'.")

    if sort_by not in df.columns:
        raise KeyError(f"sort_by='{sort_by}' not found in DE table columns.")

    # Clean types / drop NaNs in key columns
    df = df.dropna(subset=[gene_col, fc_col, sort_by]).copy()
    df[fc_col] = pd.to_numeric(df[fc_col], errors="coerce")
    df = df.dropna(subset=[fc_col]).copy()

    n_genes = int(n_genes)
    if n_genes <= 0:
        raise ValueError("n_genes must be > 0")

    # Base selection ranking:
    # - if sort_by is qval/pval: choose most significant first
    # - if sort_by is log2FC: choose strongest effects by sign
    df_sig = df.sort_values(sort_by, ascending=True) if sort_by != fc_col else df

    if direction == "up":
        # pick top up by significance (or all if sort_by=log2FC), then order by effect size
        sub = df_sig[df_sig[fc_col] > 0].copy()
        if sort_by == fc_col:
            sub = sub.sort_values(fc_col, ascending=False).head(n_genes)
        else:
            sub = sub.head(n_genes).sort_values(fc_col, ascending=False)

        # strongest up should be first (top)
        sub = sub.sort_values(fc_col, ascending=False)

    elif direction == "down":
        sub = df_sig[df_sig[fc_col] < 0].copy()
        if sort_by == fc_col:
            sub = sub.sort_values(fc_col, ascending=True).head(n_genes)  # most negative first
        else:
            sub = sub.head(n_genes).sort_values(fc_col, ascending=True)

        # most downregulated should be LAST → ensure order is weak->strong negative?
        # You requested: "most downregulated gene should be the last one"
        # That means values go from closer-to-0 (top) to most negative (bottom).
        sub = sub.sort_values(fc_col, ascending=False)  # -0.1, -0.5, -2.0 (most negative last)

    elif direction == "both":
        n_up = n_genes // 2
        n_down = n_genes - n_up

        up = df_sig[df_sig[fc_col] > 0].copy()
        down = df_sig[df_sig[fc_col] < 0].copy()

        if sort_by == fc_col:
            # strongest effects
            up = up.sort_values(fc_col, ascending=False).head(n_up)
            down = down.sort_values(fc_col, ascending=True).head(n_down)  # most negative
        else:
            # most significant
            up = up.head(n_up)
            down = down.head(n_down)

        # Now enforce final display order:
        # - Up: strongest at TOP
        up = up.sort_values(fc_col, ascending=False)

        # - Down: keep most downregulated LAST (bottom)
        #   So show from closer-to-0 → most negative
        down = down.sort_values(fc_col, ascending=False)  # e.g. -0.1, -0.5, -2.0

        sub = pd.concat([up, down], axis=0)

    else:
        raise ValueError("direction must be 'up', 'down', or 'both'.")

    genes = sub[gene_col].astype(str).tolist()
    vals = sub[fc_col].to_numpy(dtype=float)

    colors = np.where(vals >= 0, "#D62728", "#1F77B4")  # red / blue

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.barh(genes, vals, color=colors)
    ax.axvline(0, linewidth=1)

    ax.set_xlabel(fc_col)
    ax.set_ylabel("gene")

    # Put first item at top (Scanpy-like)
    ax.invert_yaxis()

    ax.set_title(title or (contrast or "Ranked DE genes"))

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax


def ma(
    *,
    result: pd.DataFrame,
    mean_col: str = "mean_norm",
    fc_col: str = "log2FC",
    gene_col: str = "gene",
    pval_col: str | None = "pval",
    qval_col: str | None = "qval",
    alpha: float = 0.05,
    use_qval: bool = True,
    min_abs_fc: float = 0.5,
    show_fc_lines: bool = True,
    fc_line_kwargs: dict | None = None,
    point_size: float = 10.0,
    point_alpha: float = 0.6,
    label_genes: Sequence[str] = (),
    top_n_labels: int = 10,
    bottom_n_labels: int = 10,
    label_only_large_fc: bool = True,
    title: str | None = "MA plot",
    xlabel: str | None = None,
    ylabel: str | None = None,
    figsize: tuple[float, float] = (8.5, 5.5),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    MA plot for DE results.

    Colors:
      - non-significant: grey
      - significant up (fc >= +min_abs_fc): red
      - significant down (fc <= -min_abs_fc): blue

    Adds automatic FC cutoff lines at ±min_abs_fc (optional).

    Parameters
    ----------
    result
        DataFrame with at least [mean_col, fc_col] and (gene_col recommended).
    alpha
        Significance threshold on q-value or p-value (see use_qval).
    min_abs_fc
        Minimum absolute fold-change (in the same units as fc_col, typically log2FC)
        to consider a gene biologically relevant and to color as up/down.
    show_fc_lines
        Draw horizontal lines at ±min_abs_fc.
    label_genes
        Explicit genes to label (if present).
    top_n_labels / bottom_n_labels
        Additionally label top/bottom genes by fold-change among "relevant" hits
        (significant + abs(fc) >= min_abs_fc if label_only_large_fc=True).
    """
    set_style()

    if result is None or not isinstance(result, pd.DataFrame):
        raise TypeError("ma(...): 'result' must be a pandas DataFrame.")

    df = result.copy()

    for c in (mean_col, fc_col):
        if c not in df.columns:
            raise KeyError(f"ma(...): '{c}' not found in result columns.")

    if gene_col not in df.columns:
        # keep plotting, but labeling will be limited
        df[gene_col] = df.index.astype(str)

    # Determine significance column
    stat_col = None
    if use_qval and qval_col is not None and qval_col in df.columns:
        stat_col = qval_col
    elif (not use_qval) and pval_col is not None and pval_col in df.columns:
        stat_col = pval_col
    else:
        # fallback to any available
        if qval_col is not None and qval_col in df.columns:
            stat_col = qval_col
            use_qval = True
        elif pval_col is not None and pval_col in df.columns:
            stat_col = pval_col
            use_qval = False

    if stat_col is None:
        # no significance info; treat all as non-significant
        sig = np.zeros(len(df), dtype=bool)
    else:
        vals = pd.to_numeric(df[stat_col], errors="coerce").to_numpy()
        sig = np.isfinite(vals) & (vals <= float(alpha))

    x = pd.to_numeric(df[mean_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[fc_col], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if not np.all(ok):
        df = df.loc[ok].copy()
        x = x[ok]
        y = y[ok]
        sig = sig[ok]

    # Define classes using BOTH significance and effect size
    min_abs_fc = float(min_abs_fc)
    up = sig & (y >= +min_abs_fc)
    down = sig & (y <= -min_abs_fc)
    nonsig = ~(up | down)

    # Plot
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # non-significant (grey)
    ax.scatter(
        x[nonsig],
        y[nonsig],
        s=point_size,
        alpha=point_alpha,
        edgecolors="none",
        c="#BDBDBD",
        rasterized=True,
        label="ns",
    )
    # down (blue)
    ax.scatter(
        x[down],
        y[down],
        s=point_size,
        alpha=point_alpha,
        edgecolors="none",
        c="#2C7FB8",
        rasterized=True,
        label=f"down (|{fc_col}|≥{min_abs_fc:g})",
    )
    # up (red)
    ax.scatter(
        x[up],
        y[up],
        s=point_size,
        alpha=point_alpha,
        edgecolors="none",
        c="#D7301F",
        rasterized=True,
        label=f"up (|{fc_col}|≥{min_abs_fc:g})",
    )

    # Baseline + FC cut lines
    ax.axhline(0.0, linewidth=1.0, color="black", alpha=0.7)
    if show_fc_lines and min_abs_fc > 0:
        kw = dict(linestyle="--", linewidth=1.0, color="black", alpha=0.6)
        if isinstance(fc_line_kwargs, dict):
            kw.update(fc_line_kwargs)
        ax.axhline(+min_abs_fc, **kw)
        ax.axhline(-min_abs_fc, **kw)

    ax.set_title(title or "MA plot")
    ax.set_xlabel(xlabel or mean_col)
    ax.set_ylabel(ylabel or fc_col)

    # --- labeling ---
    # Define "label pool"
    if label_only_large_fc:
        label_pool = (up | down)
    else:
        label_pool = sig

    df["_x"] = x
    df["_y"] = y
    df["_sig"] = sig
    df["_label_pool"] = label_pool

    # Explicit labels
    to_label = set()
    if label_genes:
        want = set(map(str, label_genes))
        hit = df[df[gene_col].astype(str).isin(want)]
        for g in hit[gene_col].astype(str).tolist():
            to_label.add(g)

    # Top/bottom labels among label pool
    pool = df[df["_label_pool"]].copy()
    if len(pool) > 0:
        pool = pool.sort_values(fc_col, ascending=True)  # most down first
        # bottom_n_labels: most downregulated (negative) -> from the start
        if bottom_n_labels and bottom_n_labels > 0:
            for g in pool.head(int(bottom_n_labels))[gene_col].astype(str).tolist():
                to_label.add(g)
        # top_n_labels: most upregulated (positive) -> from the end
        if top_n_labels and top_n_labels > 0:
            for g in pool.tail(int(top_n_labels))[gene_col].astype(str).tolist():
                to_label.add(g)

    if to_label:
        lab = df[df[gene_col].astype(str).isin(to_label)].copy()
        # small repel-like jitter (deterministic) to reduce exact overlaps a bit
        rng = np.random.default_rng(0)
        for _, r in lab.iterrows():
            ax.text(
                float(r["_x"]),
                float(r["_y"]) + float(rng.normal(0, 0.03)),
                str(r[gene_col]),
                fontsize=max(8, plt.rcParams.get("font.size", 10) - 2),
                ha="center",
                va="bottom",
            )

    # Optional legend (kept minimal)
    # Comment out if you prefer no legend:
    ax.legend(frameon=False, loc="best")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax