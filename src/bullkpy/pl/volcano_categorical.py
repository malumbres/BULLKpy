from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ._style import set_style, _savefig
from .._deprecated import renamed_kwargs


@renamed_kwargs(q_col="qval_col")
def volcano_categorical(
    results: pd.DataFrame,
    *,
    x: str = "effect_size",   # or "log2FC"
    qval_col: str = "qval",
    label_top: int = 12,
    q_thr: float = 0.05,
    effect_thr: float | None = None,
    figsize: tuple[float, float] = (6.8, 5.2),
    title: str | None = None,
    save: str | Path | None = None,
    show: bool = True, ax: plt.Axes | None = None,
):
    """
    Volcano plot for categorical associations.
    Expects columns: gene, qval, and x (effect_size or log2FC).
    """
    set_style()

    df = results.copy()
    if "gene" not in df.columns:
        raise KeyError("results must contain column 'gene'")
    if qval_col not in df.columns:
        raise KeyError(f"results must contain column '{qval_col}'")
    if x not in df.columns:
        raise KeyError(f"results must contain column '{x}'")

    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, qval_col])
    df["mlog10q"] = -np.log10(np.clip(df[qval_col].to_numpy(dtype=float), 1e-300, 1.0))

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    else:
        fig = ax.figure

    sig = df[qval_col].to_numpy(dtype=float) <= float(q_thr)
    ax.scatter(df.loc[~sig, x], df.loc[~sig, "mlog10q"], alpha=0.35, s=12, edgecolors="none")
    ax.scatter(df.loc[sig, x], df.loc[sig, "mlog10q"], alpha=0.85, s=14, edgecolors="none")

    ax.axhline(-np.log10(q_thr), ls="--", lw=1)

    if effect_thr is not None:
        ax.axvline(-abs(effect_thr), ls="--", lw=1)
        ax.axvline(abs(effect_thr), ls="--", lw=1)

    ax.set_xlabel(x)
    ax.set_ylabel(f"-log10({qval_col})")
    if title is not None:
        ax.set_title(title)

    # label top hits
    top = df.sort_values(qval_col, ascending=True).head(int(label_top))
    for _, r in top.iterrows():
        ax.text(float(r[x]), float(r["mlog10q"]), str(r["gene"]), fontsize=9)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax