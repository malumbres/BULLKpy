from __future__ import annotations

import logging as _logging
from datetime import datetime, timezone

from ._settings import settings

_LOGGER_NAME = "bullkpy"


def _verbosity_to_level(verbosity) -> int:
    """
    Map verbosity (int or Verbosity dataclass) to python logging level.
    Accepts:
      - int: 0..3 (scanpy-like)
      - Verbosity: object with .level (python logging level) OR .value/.verbosity
    """
    # If Verbosity dataclass/object
    if hasattr(verbosity, "level"):
        v = verbosity.level
        # if already a logging level (10/20/30/40)
        if isinstance(v, int) and v in (10, 20, 30, 40, 50):
            return v
        verbosity = v  # fall back to treating it as an int scale

    # Normalize to int scale (0..3) or a logging level
    if isinstance(verbosity, int) and verbosity in (10, 20, 30, 40, 50):
        return verbosity

    try:
        verbosity = int(verbosity)
    except Exception:
        return _logging.INFO

    if verbosity <= 0:
        return _logging.ERROR
    elif verbosity == 1:
        return _logging.WARNING
    elif verbosity == 2:
        return _logging.INFO
    else:
        return _logging.DEBUG


def _get_logger() -> _logging.Logger:
    """
    Get (or create) the BULLKpy logger.

    This function is intentionally self-contained to avoid circular imports.
    """
    logger = _logging.getLogger(_LOGGER_NAME)

    if not logger.handlers:
        handler = _logging.StreamHandler()
        formatter = _logging.Formatter(
            fmt="%(levelname)s: %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False

    logger.setLevel(_verbosity_to_level(settings.verbosity))
    return logger


def info(msg: str) -> None:
    """Log an info-level message."""
    _get_logger().info(msg)


def warn(msg: str) -> None:
    """Log a warning-level message."""
    _get_logger().warning(msg)


def debug(msg: str) -> None:
    """Log a debug-level message."""
    _get_logger().debug(msg)


def error(msg: str) -> None:
    """Log an error-level message."""
    _get_logger().error(msg)


def time_stamp() -> str:
    """Return a UTC timestamp string."""
    return datetime.now(timezone.utc).isoformat()