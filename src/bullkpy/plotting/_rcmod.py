"""Matplotlib rcParams helpers (Scanpy-inspired)."""

from __future__ import annotations

from typing import Optional
import matplotlib as mpl


def set_figure_params(
    *,
    dpi: Optional[int] = None,
    facecolor: Optional[str] = "white",
    transparent: Optional[bool] = False,
    format: Optional[str] = "png",
) -> None:
    if dpi is not None:
        mpl.rcParams["figure.dpi"] = dpi
        mpl.rcParams["savefig.dpi"] = dpi
    if facecolor is not None:
        mpl.rcParams["figure.facecolor"] = facecolor
        mpl.rcParams["axes.facecolor"] = facecolor
        mpl.rcParams["savefig.facecolor"] = facecolor
    if transparent is not None:
        mpl.rcParams["savefig.transparent"] = transparent
    if format is not None:
        mpl.rcParams["savefig.format"] = format
