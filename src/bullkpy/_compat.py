"""Small dtype helpers that behave the same across pandas 2.x and 3.x.

pandas 3 introduced a dedicated ``str`` dtype, so string columns are no longer
``object`` dtype. Code written as ``s.dtype == object`` silently stopped
recognising string columns as categorical, which routed them into numeric code
paths. These helpers express the intent directly instead.
"""
from __future__ import annotations

import pandas as pd

__all__ = ["is_categorical_like", "is_numeric_like", "is_true_categorical"]


def is_numeric_like(s) -> bool:
    """True for integer/float (and other numeric) columns."""
    return bool(pd.api.types.is_numeric_dtype(getattr(s, "dtype", s)))


def is_true_categorical(s) -> bool:
    """True only for pandas Categorical columns."""
    return isinstance(getattr(s, "dtype", s), pd.CategoricalDtype)


def is_categorical_like(s) -> bool:
    """True for columns that should be treated as discrete groups.

    Covers pandas ``category``, legacy ``object`` and pandas 3 ``str`` dtypes,
    plus booleans; numeric and datetime columns are excluded.
    """
    dtype = getattr(s, "dtype", s)
    if isinstance(dtype, pd.CategoricalDtype):
        return True
    if pd.api.types.is_bool_dtype(dtype):
        return True
    if pd.api.types.is_numeric_dtype(dtype):
        return False
    if pd.api.types.is_datetime64_any_dtype(dtype):
        return False
    return bool(dtype == object or pd.api.types.is_string_dtype(dtype))
