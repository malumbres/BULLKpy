from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

from ._style import set_style, _savefig


Metric = Literal["ARI", "NMI", "cramers_v"]


def ari_resolution_heatmap(
    adata: ad.AnnData,
    *,
    df: pd.DataFrame | None = None,
    store_key: str = "leiden_scan",
    metric: Metric = "ARI",
    show_n_clusters: bool = True,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Heatmap-like summary of clustering quality vs Leiden resolution.

    Expects `df` with at least:
      - 'resolution'
      - metric column: 'ARI' or 'NMI' or 'cramers_v'
    Optional:
      - 'n_clusters' (for a second row annotation)

    Example:
        df = bk.tl.leiden_resolution_scan(...)
        bk.pl.ari_resolution_heatmap(adata, df=df, metric="ARI")
    """
    set_style()

    if df is None:
        if store_key not in adata.uns:
            raise KeyError(
                f"adata.uns['{store_key}'] not found. Run bk.tl.leiden_resolution_scan(...) first "
                f"or pass df=..."
            )
        df = adata.uns[store_key]
        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame(df)

    required = {"resolution", metric}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"Missing columns in df: {sorted(missing)}")

    d = df.copy().sort_values("resolution").reset_index(drop=True)

    # columns = resolutions (as strings for tick labels)
    cols = [f"{r:g}" for r in d["resolution"].astype(float).to_numpy()]
    n = len(cols)

    metric_vals = np.asarray(d[metric].to_numpy(dtype=float))[None, :]  # (1, n)

    have_clusters = bool(show_n_clusters and ("n_clusters" in d.columns))
    if have_clusters:
        cluster_vals = np.asarray(d["n_clusters"].to_numpy(dtype=float))[None, :]  # (1, n)

    rows = [metric] + (["n_clusters"] if have_clusters else [])
    n_rows = len(rows)

    # default size
    if figsize is None:
        w = max(6.0, 0.55 * n + 1.8)
        h = 2.3 if n_rows == 1 else 3.2
        figsize = (w, h)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # ---- draw images with explicit extents (no ax.images.pop) ----
    # We use a coordinate system where:
    #   x in [0, n], y in [0, n_rows]
    # and each cell is 1×1 centered at (j+0.5, i+0.5).
    im_metric = ax.imshow(
        metric_vals,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=(0, n, 0, 1),  # row 0
        origin="lower",
    )

    if have_clusters:
        ax.imshow(
            cluster_vals,
            aspect="auto",
            interpolation="nearest",
            cmap="Greys",
            extent=(0, n, 1, 2),  # row 1
            origin="lower",
        )

    # ---- ticks/labels ----
    ax.set_xlim(0, n)
    ax.set_ylim(0, n_rows)

    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels(cols, rotation=45, ha="right")
    ax.set_xlabel("Leiden resolution")

    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(rows)

    # gridlines between cells (scanpy-ish)
    ax.set_xticks(np.arange(n + 1), minor=True)
    ax.set_yticks(np.arange(n_rows + 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # ---- annotate metric cells ----
    for j in range(n):
        v = float(metric_vals[0, j])
        if np.isfinite(v):
            ax.text(j + 0.5, 0.5, f"{v:.2f}", ha="center", va="center", fontsize=9, color="white")

    if have_clusters:
        for j in range(n):
            v = float(cluster_vals[0, j])
            if np.isfinite(v):
                ax.text(j + 0.5, 1.5, f"{int(round(v))}", ha="center", va="center", fontsize=9, color="black")

    ax.set_title(title or f"{metric} vs Leiden resolution")

    # colorbar for metric only
    cbar = fig.colorbar(im_metric, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(metric)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax