"""Keep renamed functions and parameters working, with a warning, for one deprecation cycle."""
from __future__ import annotations

import functools
import importlib
import warnings
from typing import Callable

__all__ = ["renamed_kwargs", "RENAMED_KWARGS", "moved_functions", "MOVED"]

#: "module.function" -> {old parameter name: new parameter name}. Used by the tests and
#: the notebook contract check to recognise deprecated spellings.
RENAMED_KWARGS: dict[str, dict[str, str]] = {}

#: old public name ("namespace.name") -> new public name ("namespace.name")
MOVED: dict[str, str] = {}


def renamed_kwargs(**old_to_new: str) -> Callable:
    """Accept ``old=`` keyword arguments for parameters that were renamed to ``new=``."""

    def decorator(func: Callable) -> Callable:
        RENAMED_KWARGS[f"{func.__module__}.{func.__qualname__}"] = dict(old_to_new)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for old, new in old_to_new.items():
                if old in kwargs:
                    if new in kwargs:
                        raise TypeError(f"{func.__name__}() got both '{old}' (deprecated) and '{new}'")
                    warnings.warn(
                        f"{func.__name__}(): parameter '{old}' was renamed to '{new}' and "
                        f"'{old}' will be removed in a future release.",
                        DeprecationWarning,
                        stacklevel=2,
                    )
                    kwargs[new] = kwargs.pop(old)
            return func(*args, **kwargs)

        wrapper.__deprecated_kwargs__ = dict(old_to_new)
        return wrapper

    return decorator


def moved_functions(namespace: str, moved: dict[str, str]) -> Callable[[str], object]:
    """Build a module ``__getattr__`` that serves functions whose home changed.

    ``moved`` maps the old attribute name in ``namespace`` (e.g. ``"bullkpy.tl"``) to the
    new location as ``"bullkpy.pl.plot_corr_scatter"``. Access still works and warns.
    """
    for old, target in moved.items():
        MOVED[f"{namespace.removeprefix('bullkpy.')}.{old}"] = target.removeprefix("bullkpy.")

    def __getattr__(name: str):
        if name in moved:
            module_name, _, attr = moved[name].rpartition(".")
            new_location = moved[name].removeprefix("bullkpy.")
            warnings.warn(
                f"{namespace.removeprefix('bullkpy.')}.{name} moved to {new_location}; "
                f"the old location will be removed in a future release.",
                DeprecationWarning,
                stacklevel=2,
            )
            return getattr(importlib.import_module(module_name), attr)
        raise AttributeError(f"module {namespace!r} has no attribute {name!r}")

    return __getattr__
