"""Refactors must not change numeric results. See golden_cases.py to regenerate."""
from __future__ import annotations

import json

import numpy as np
import pytest

from golden_cases import GOLDEN, compute


def _compare(actual, expected, path=""):
    if isinstance(expected, dict):
        assert set(actual) == set(expected), f"{path}: keys {set(actual)} != {set(expected)}"
        for k in expected:
            _compare(actual[k], expected[k], f"{path}/{k}")
    elif isinstance(expected, list) and expected and isinstance(expected[0], (int, float, type(None))):
        a = np.array([np.nan if v is None else v for v in actual], dtype=float)
        e = np.array([np.nan if v is None else v for v in expected], dtype=float)
        assert a.shape == e.shape, path
        np.testing.assert_allclose(a, e, rtol=1e-6, atol=1e-8, equal_nan=True, err_msg=path)
    else:
        assert actual == expected, path


def test_numeric_results_are_unchanged(adata, genes):
    _compare(compute(adata, genes), json.loads(GOLDEN.read_text()))
