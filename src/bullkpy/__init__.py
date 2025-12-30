"""BULLKpy: a Scanpy-inspired toolkit for bulk RNA-seq analysis in Python.

This project is inspired by Scanpy (BSD 3-Clause) and reuses some design patterns
(settings/logging/plotting style) with attribution in source headers where applicable.
"""

from __future__ import annotations

from ._settings import settings
from .io import read_counts, add_metadata

# --- auto Jupyter fix (Scanpy-like) ---
try:
    from ._jupyter import in_jupyter, apply_inline_dpi_defaults

    if in_jupyter():
        apply_inline_dpi_defaults()
except Exception:
    pass

from . import pp, pl, tl, get, logging, io
from .logging import info, warn, debug, error  # noqa: F401

__all__ = [
    "settings",
    "read_counts",
    "add_metadata",
    "pp",
    "pl",
    "tl",
    "info",
    "warn",
    "debug",
    "error",
    "get",
    "io",
    "logging",
]

__version__ = "0.0.1"