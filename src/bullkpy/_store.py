"""Optionally store a function's result table in ``adata.uns``."""
from __future__ import annotations

import functools
import inspect
from typing import Callable

__all__ = ["stores_result"]

_DOC = """

    Storing the result
    ------------------
    key_added
        If given, also store ``{"results": table, "params": ...}`` in
        ``adata.uns[key_added][subkey]``. The table is returned either way.
    subkey
        Name under ``key_added``. Defaults to ``__DEFAULT__``.
"""


def stores_result(*, params: tuple[str, ...], default_subkey: str, subkey_of: Callable[[dict], str]):
    """Give a function that returns a table the optional ``key_added=`` / ``subkey=`` storage.

    ``params`` lists the scalar arguments kept as provenance; ``subkey_of`` builds the default
    subkey from the bound arguments.
    """

    def decorator(func: Callable) -> Callable:
        sig = inspect.signature(func)
        extra = [
            inspect.Parameter("key_added", inspect.Parameter.KEYWORD_ONLY, default=None, annotation="str | None"),
            inspect.Parameter("subkey", inspect.Parameter.KEYWORD_ONLY, default=None, annotation="str | None"),
        ]

        @functools.wraps(func)
        def wrapper(adata, *args, key_added=None, subkey=None, **kwargs):
            result = func(adata, *args, **kwargs)
            if key_added is not None:
                bound = sig.bind(adata, *args, **kwargs)
                bound.apply_defaults()
                arguments = dict(bound.arguments)
                chosen = subkey or subkey_of(arguments)
                adata.uns.setdefault(key_added, {})[chosen] = {
                    "results": result,
                    "params": {k: arguments[k] for k in params if arguments.get(k) is not None},
                }
            return result

        wrapper.__signature__ = sig.replace(parameters=[*sig.parameters.values(), *extra])
        wrapper.__doc__ = (func.__doc__ or "").rstrip() + _DOC.replace("__DEFAULT__", default_subkey)
        return wrapper

    return decorator
