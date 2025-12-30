from __future__ import annotations

from pathlib import Path
from typing import Literal, Sequence

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.sparse as sp

import anndata as ad

from ._style import set_style, _savefig


def _get_color_vector(
    adata: ad.AnnData,
    key: str,
    *,
    layer: str | None = None,
) -> tuple[np.ndarray, str, str]:
    """
    Returns (values, kind, label)
      kind in {"numeric", "categorical"}.
    key can be obs column or gene in var_names.
    """
    # obs
    if key in adata.obs.columns:
        s = adata.obs[key]
        if pd.api.types.is_numeric_dtype(s):
            return s.to_numpy(dtype=float), "numeric", key
        return s.astype(str).to_numpy(dtype=object), "categorical", key

    # gene
    if key in adata.var_names:
        X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
        gidx = int(adata.var_names.get_loc(key))
        if sp.issparse(X):
            vals = np.asarray(X[:, gidx].toarray()).ravel()
        else:
            vals = np.asarray(X[:, gidx], dtype=float).ravel()
        return vals.astype(float), "numeric", key

    raise KeyError(f"color/hue key '{key}' not found in adata.obs columns or adata.var_names.")


def _corr_stats(
    xvals: np.ndarray,
    yvals: np.ndarray,
    *,
    method: Literal["pearson", "spearman", "both"] = "both",
) -> dict[str, float | int]:
    out: dict[str, float | int] = {"n": int(len(xvals))}
    if method in ("pearson", "both"):
        try:
            from scipy.stats import pearsonr
            r, p = pearsonr(xvals, yvals)
            out["pearson_r"] = float(r)
            out["pearson_p"] = float(p)
        except Exception:
            out["pearson_r"] = float("nan")
            out["pearson_p"] = float("nan")
    if method in ("spearman", "both"):
        try:
            from scipy.stats import spearmanr
            rho, p = spearmanr(xvals, yvals, nan_policy="omit")
            out["spearman_rho"] = float(rho)
            out["spearman_p"] = float(p)
        except Exception:
            out["spearman_rho"] = float("nan")
            out["spearman_p"] = float("nan")
    return out


def _plot_one(
    ax: plt.Axes,
    df: pd.DataFrame,
    *,
    x: str,
    y: str,
    color_key: str | None,
    layer: str | None,
    palette: str,
    cmap: str,
    point_size: float,
    alpha: float,
    add_regline: bool,
    annotate: bool,
    method: Literal["pearson", "spearman", "both"],
    title: str | None,
    legend: bool,
) -> dict[str, float | int]:
    xvals = df[x].to_numpy(float)
    yvals = df[y].to_numpy(float)

    stats = _corr_stats(xvals, yvals, method=method)

    slope = intercept = None
    if add_regline:
        slope, intercept = np.polyfit(xvals, yvals, 1)
        stats["slope"] = float(slope)
        stats["intercept"] = float(intercept)

    if color_key is None:
        ax.scatter(xvals, yvals, s=point_size, alpha=alpha, edgecolors="none")
    else:
        vals, kind, label = _get_color_vector(adata=df.attrs["_adata_"], key=color_key, layer=layer)
        # align to df index
        vals = pd.Series(vals, index=df.attrs["_adata_"].obs_names).reindex(df.index).to_numpy()

        if kind == "numeric":
            sc = ax.scatter(
                xvals, yvals,
                c=vals.astype(float),
                cmap=cmap,
                s=point_size,
                alpha=alpha,
                edgecolors="none",
            )
            if legend:
                cbar = ax.figure.colorbar(sc, ax=ax, pad=0.02, fraction=0.05)
                cbar.set_label(label)
        else:
            cats = pd.Categorical(pd.Series(vals).astype(str))
            names = list(cats.categories)

            pal = mpl.cm.get_cmap(palette) if palette in plt.colormaps() else None
            if pal is None:
                pal = mpl.cm.get_cmap("tab20")
            cmap_map = {str(n): pal(i % pal.N) for i, n in enumerate(names)}

            for name in names:
                m = (cats == name)
                ax.scatter(
                    xvals[m],
                    yvals[m],
                    s=point_size,
                    alpha=alpha,
                    edgecolors="none",
                    color=cmap_map[str(name)],
                    label=str(name),
                )

            if legend:
                ax.legend(
                    title=label,
                    bbox_to_anchor=(1.02, 1),
                    loc="upper left",
                    frameon=False,
                )

    ax.set_xlabel(x)
    ax.set_ylabel(y)

    if add_regline and slope is not None and intercept is not None:
        xs = np.linspace(np.nanmin(xvals), np.nanmax(xvals), 200)
        ax.plot(xs, slope * xs + intercept)

    if annotate:
        lines = [f"n={stats['n']}"]
        if "pearson_r" in stats:
            lines.append(f"Pearson r={stats['pearson_r']:.3g}, p={stats['pearson_p']:.2g}")
        if "spearman_rho" in stats:
            lines.append(f"Spearman ρ={stats['spearman_rho']:.3g}, p={stats['spearman_p']:.2g}")
        ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, ha="left", va="top")

    ax.set_title(title or (color_key if color_key is not None else f"{x} vs {y}"))
    return stats


def corrplot_obs(
    adata: ad.AnnData,
    *,
    x: str,
    y: str,
    color: str | Sequence[str] | None = None,
    hue: str | Sequence[str] | None = None,  # alias for color
    layer: str | None = None,
    palette: str = "tab20",
    cmap: str = "viridis",
    legend: bool = True,
    method: Literal["pearson", "spearman", "both"] = "both",
    add_regline: bool = True,
    annotate: bool = True,
    dropna: bool = True,
    point_size: float = 18.0,
    alpha: float = 0.75,
    figsize: tuple[float, float] = (5.5, 4.5),
    panel_size: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Scatter + correlations between two quantitative obs columns.

    Multi-panel:
      color=["DLL3","SOX10"] makes one panel per color key in a single row.
    """
    set_style()

    if hue is not None and color is None:
        color = hue

    if x not in adata.obs.columns:
        raise KeyError(f"x='{x}' not in adata.obs")
    if y not in adata.obs.columns:
        raise KeyError(f"y='{y}' not in adata.obs")

    df = adata.obs[[x, y]].copy()
    df[x] = pd.to_numeric(df[x], errors="coerce")
    df[y] = pd.to_numeric(df[y], errors="coerce")

    # store for _plot_one
    df.attrs["_adata_"] = adata

    if dropna:
        df = df.dropna(subset=[x, y])

    if df.shape[0] < 3:
        raise ValueError(f"Not enough valid points for correlation: n={df.shape[0]}")

    # normalize color -> list
    if color is None:
        colors = [None]
    elif isinstance(color, (list, tuple)):
        colors = [str(c) for c in color]
    else:
        colors = [str(color)]

    n = len(colors)

    if panel_size is not None:
        figsize_eff = (panel_size[0] * n, panel_size[1])
    else:
        figsize_eff = (figsize[0] * n, figsize[1]) if n > 1 else figsize

    fig, axes = plt.subplots(1, n, figsize=figsize_eff, constrained_layout=True)
    if n == 1:
        axes = [axes]

    stats_list: list[dict[str, float | int]] = []
    for ax, c in zip(axes, colors):
        st = _plot_one(
            ax,
            df,
            x=x,
            y=y,
            color_key=c,
            layer=layer,
            palette=palette,
            cmap=cmap,
            point_size=point_size,
            alpha=alpha,
            add_regline=add_regline,
            annotate=annotate,
            method=method,
            title=(title if (title is not None and n == 1) else None),
            legend=legend,
        )
        stats_list.append(st)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, np.array(axes, dtype=object), stats_list