"""Renamed parameters and moved functions keep working, with a DeprecationWarning."""
from __future__ import annotations

import importlib
import inspect
import warnings

import numpy as np
import pytest

import bullkpy as bk
from bullkpy._deprecated import MOVED, RENAMED_KWARGS


def _resolve(qualified: str):
    module, _, name = qualified.rpartition(".")
    return getattr(importlib.import_module(module), name)


@pytest.mark.parametrize("qualified", sorted(RENAMED_KWARGS))
def test_signature_has_new_names_and_not_old(qualified):
    func = _resolve(qualified)
    params = set(inspect.signature(func).parameters)
    for old, new in RENAMED_KWARGS[qualified].items():
        assert new in params, f"{qualified}: new parameter {new!r} missing"
        assert old not in params, f"{qualified}: old parameter {old!r} still in the signature"
        assert func.__deprecated_kwargs__[old] == new


def test_renamed_parameters_cover_the_documented_set():
    # Guards against the decorator silently not being applied.
    assert len(RENAMED_KWARGS) >= 35


def test_old_keyword_works_and_warns(adata):
    new = bk.tl.cox_univariate(
        adata, time_col="OS_time", event_col="OS_event", x_keys=["age"], groupby=None, strata=None
    )
    with pytest.warns(DeprecationWarning, match="time_key.*time_col"):
        old = bk.tl.cox_univariate(
            adata, time_key="OS_time", event_key="OS_event", x_keys=["age"], groupby=None, strata=None
        )
    assert type(old) is type(new)


def test_old_and_new_together_is_an_error(adata):
    with pytest.raises(TypeError, match="both"):
        bk.tl.cox_univariate(
            adata, time_key="OS_time", time_col="OS_time", event_col="OS_event", x_keys=["age"]
        )


def test_random_state_is_now_seed(adata_pca):
    bk.tl.cluster(adata_pca, method="kmeans", n_clusters=3, use_rep="X_pca", n_pcs=5, seed=1)
    a = adata_pca.obs["clusters"].astype(str).to_numpy().copy()
    with pytest.warns(DeprecationWarning, match="random_state"):
        bk.tl.cluster(adata_pca, method="kmeans", n_clusters=3, use_rep="X_pca", n_pcs=5, random_state=1)
    assert (adata_pca.obs["clusters"].astype(str).to_numpy() == a).all()


@pytest.mark.parametrize("old", sorted(MOVED))
def test_moved_functions_resolve_and_warn(old):
    namespace, _, name = old.rpartition(".")
    new_namespace, _, new_name = MOVED[old].rpartition(".")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        via_old = getattr(getattr(bk, namespace), name)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert via_old is getattr(getattr(bk, new_namespace), new_name)
    assert name not in getattr(getattr(bk, namespace), "__all__", [])
