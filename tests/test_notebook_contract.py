"""The published tutorial must stay callable against the shipped API.

The notebook is the only tutorial on Read the Docs and its stored outputs make
stale calls look like they work, so drift here is invisible without a check.
These tests bind every ``bk.*`` call in the notebook against the real
signatures without executing anything.
"""
from __future__ import annotations

import ast
import inspect
import json
import pathlib

import pytest

import bullkpy as bk

ROOT = pathlib.Path(__file__).resolve().parents[1]
NOTEBOOKS = sorted((ROOT / "notebooks").glob("*.ipynb"))

SUBMODULES = {"pp", "tl", "pl", "get", "io"}


def _code_cells(path: pathlib.Path):
    nb = json.loads(path.read_text())
    for idx, cell in enumerate(nb.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        # drop IPython magics and shell escapes so the cell parses as Python
        src = "\n".join(
            line for line in "".join(cell["source"]).splitlines()
            if not line.strip().startswith(("%", "!", "?"))
        )
        yield idx, src


def _bk_calls(src: str):
    """Yield (submodule, name, ast.Call) for every bk.<sub>.<name>(...) call."""
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if (
            isinstance(fn, ast.Attribute)
            and isinstance(fn.value, ast.Attribute)
            and isinstance(fn.value.value, ast.Name)
            and fn.value.value.id == "bk"
            and fn.value.attr in SUBMODULES
        ):
            yield fn.value.attr, fn.attr, node


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_calls_existing_functions(nb_path):
    missing = []
    for idx, src in _code_cells(nb_path):
        for sub, name, _ in _bk_calls(src):
            if not hasattr(getattr(bk, sub), name):
                missing.append(f"cell {idx}: bk.{sub}.{name}")
    assert not missing, "notebook calls functions that do not exist:\n" + "\n".join(missing)


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_calls_match_signatures(nb_path):
    """Catch calls that would raise TypeError before any computation happens."""
    problems = []
    for idx, src in _code_cells(nb_path):
        for sub, name, node in _bk_calls(src):
            obj = getattr(getattr(bk, sub), name, None)
            if obj is None or not callable(obj):
                continue
            try:
                sig = inspect.signature(obj)
            except (ValueError, TypeError):
                continue
            # argument unpacking makes the call unverifiable statically
            if any(k.arg is None for k in node.keywords) or any(
                isinstance(a, ast.Starred) for a in node.args
            ):
                continue
            try:
                sig.bind(*([None] * len(node.args)), **{k.arg: None for k in node.keywords})
            except TypeError as exc:
                problems.append(f"cell {idx}: bk.{sub}.{name}(): {exc}")

    assert not problems, "notebook calls no longer match the API:\n" + "\n".join(problems)


@pytest.mark.parametrize("nb_path", NOTEBOOKS, ids=lambda p: p.name)
def test_notebook_has_no_foreign_absolute_paths(nb_path):
    """Absolute paths must point at this machine's data root, not a previous one."""
    allowed_prefix = "/Users/mmalumbres/Library/CloudStorage/OneDrive-VHIO/"
    offenders = []
    for idx, src in _code_cells(nb_path):
        for node in ast.walk(ast.parse(src)) if _parses(src) else []:
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                v = node.value
                if v.startswith("/Users/") and not v.startswith(allowed_prefix):
                    offenders.append(f"cell {idx}: {v}")
    assert not offenders, "notebook references paths outside the data root:\n" + "\n".join(offenders)


def _parses(src: str) -> bool:
    try:
        ast.parse(src)
        return True
    except SyntaxError:
        return False
