"""BULLKpy: a Scanpy-inspired toolkit for bulk RNA-seq analysis in Python.

This project is inspired by Scanpy (BSD 3-Clause) and reuses some design patterns
(settings/logging/plotting style) with attribution in source headers where applicable.
"""

from __future__ import annotations

from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("bullkpy")
except PackageNotFoundError:
    __version__ = "unknown"

from ._settings import settings
from .io import read_counts, add_metadata

from . import pp, pl, tl, get, logging, io, datasets
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
    "datasets",
    "logging",
    "__version__",
]