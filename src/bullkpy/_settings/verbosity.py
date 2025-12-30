"""Verbosity configuration (inspired by scanpy._settings.verbosity)."""

from __future__ import annotations

from dataclasses import dataclass
import logging


@dataclass
class Verbosity:
    """Verbosity state."""
    level: int = logging.INFO
