from __future__ import annotations

from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from bullkpy._settings import settings


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

def _scaled_fontsize(figsize: tuple[float, float], base: float) -> float:
    """
    Scanpy-like font scaling:
    - fontsize stays ~constant on screen
    - small figures do NOT get huge text
    - large figures do NOT get tiny text
    """
    if not settings.scale_fonts_with_figsize:
        return base

    # reference: ~6x6 inches is "normal"
    ref = 6.0
    scale = math.sqrt((figsize[0] * figsize[1]) / (ref * ref))
    scale = max(0.75, min(scale, 1.25))  # clamp like scanpy
    return base * scale


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def set_style(figsize: tuple[float, float] | None = None) -> None:
    """
    Apply BULLKpy global plotting style (matplotlib + seaborn).

    This is intentionally lightweight and can be called repeatedly.
    """
    # ---- base font size ----
    fs = float(settings.plot_fontsize) if settings.plot_fontsize is not None else 12.0
    dpi = int(settings.plot_dpi) if settings.plot_dpi is not None else 150

    if figsize is not None:
        fs = _scaled_fontsize(figsize, fs)

    # ---- matplotlib rcParams ----
    mpl.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": dpi,
            "font.size": fs,
            "axes.titlesize": fs,
            "axes.labelsize": fs,
            "xtick.labelsize": fs - 1,
            "ytick.labelsize": fs - 1,
            "legend.fontsize": fs - 1,
            "legend.title_fontsize": fs,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    # ---- seaborn theme ----
    if settings.plot_theme == "paper":
        sns.set_theme(
            style="white",
            context="paper",
            font_scale=1.0,
        )
    elif settings.plot_theme == "talk":
        sns.set_theme(
            style="white",
            context="talk",
            font_scale=1.2,
        )
    else:  # default
        sns.set_theme(
            style="whitegrid",
            context="notebook",
            font_scale=1.0,
        )

    # ---- default categorical palette ----
    sns.set_palette(settings.plot_palette)


def _savefig(fig: plt.Figure, path: str | Path) -> None:
    """
    Save figure respecting bk.settings.figdir and dpi.
    """
    path = Path(path)
    if not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        path,
        dpi=settings.plot_dpi,
        bbox_inches="tight",
    )

def get_palette(n: int | None = None, palette: str | None = None):
    """
    Return a seaborn palette (list of RGB tuples).
    If n is None, returns the default palette object.
    """
    pal_name = palette or settings.plot_palette
    if n is None:
        return sns.color_palette(pal_name)
    return sns.color_palette(pal_name, n_colors=int(n))

def _apply_clustergrid_style(cg, *, fontsize: float | None = None):
    """
    Style a seaborn ClusterGrid (clustermap) in a Scanpy-like compact way.
    Safe to call even if some axes are missing.
    """
    fs = float(fontsize) if fontsize is not None else float(settings.plot_fontsize)

    # Heatmap axis
    ax = getattr(cg, "ax_heatmap", None)
    if ax is not None:
        ax.tick_params(axis="x", labelrotation=90, labelsize=max(fs - 2, 6))
        ax.tick_params(axis="y", labelsize=max(fs - 2, 6))
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Dendrograms
    for a in [getattr(cg, "ax_row_dendrogram", None), getattr(cg, "ax_col_dendrogram", None)]:
        if a is not None:
            a.set_xticks([])
            a.set_yticks([])
            for spn in a.spines.values():
                spn.set_visible(False)

    # Colorbar
    cax = getattr(cg, "cax", None)
    if cax is not None:
        cax.tick_params(labelsize=max(fs - 2, 6))

    # Col/row color labels (if present)
    for a in [getattr(cg, "ax_col_colors", None), getattr(cg, "ax_row_colors", None)]:
        if a is not None:
            a.set_xticks([])
            a.set_yticks([])
            for spn in a.spines.values():
                spn.set_visible(False)