"""Contracts that keep the public API, its exports and its docs in sync.

These are the checks that caught the drift fixed in the 0.1.1 cleanup, so they
are worth keeping green.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

import bullkpy as bk

SUBMODULES = ["pp", "tl", "pl", "get", "io"]
ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "bullkpy"
DOCS_API = ROOT / "docs" / "api"


def _exports(name: str) -> list[str]:
    return list(getattr(getattr(bk, name), "__all__", []))


@pytest.mark.parametrize("sub", SUBMODULES)
def test_all_declares_exports(sub):
    assert _exports(sub), f"bullkpy.{sub} must declare __all__"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_every_export_resolves(sub):
    mod = getattr(bk, sub)
    missing = [n for n in _exports(sub) if not hasattr(mod, n)]
    assert not missing, f"bullkpy.{sub}.__all__ names that do not exist: {missing}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_no_duplicate_exports(sub):
    all_ = _exports(sub)
    dupes = sorted({n for n in all_ if all_.count(n) > 1})
    assert not dupes, f"duplicated entries in bullkpy.{sub}.__all__: {dupes}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_no_private_names_exported(sub):
    private = [n for n in _exports(sub) if n.startswith("_")]
    assert not private, f"private names exported from bullkpy.{sub}: {private}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_every_export_is_documented(sub):
    """Each public function needs a docs/api/<sub>/<name>.md page."""
    docdir = DOCS_API / sub
    pages = {p.stem for p in docdir.glob("*.md")} - {"index"}
    undocumented = sorted(set(_exports(sub)) - pages)
    assert not undocumented, f"exported but undocumented in docs/api/{sub}/: {undocumented}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_no_orphan_doc_pages(sub):
    """No docs page may describe a function that no longer exists."""
    docdir = DOCS_API / sub
    pages = {p.stem for p in docdir.glob("*.md")} - {"index"}
    orphans = sorted(pages - set(_exports(sub)))
    assert not orphans, f"docs/api/{sub}/ pages with no matching export: {orphans}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_toctree_matches_files(sub):
    """Every docs page is reachable from its index.md toctree, and vice versa."""
    idx = DOCS_API / sub / "index.md"
    assert idx.exists(), f"missing docs/api/{sub}/index.md"

    listed, in_toc = set(), False
    for line in idx.read_text().splitlines():
        s = line.strip()
        if s.startswith("```{toctree}"):
            in_toc = True
            continue
        if in_toc and s.startswith("```"):
            in_toc = False
            continue
        if in_toc and s and not s.startswith(":"):
            listed.add(s.split("<")[-1].rstrip(">").strip())

    actual = {p.stem for p in (DOCS_API / sub).glob("*.md")} - {"index"}
    assert not (actual - listed), f"docs/api/{sub}: files missing from toctree: {sorted(actual - listed)}"
    assert not (listed - actual), f"docs/api/{sub}: toctree entries with no file: {sorted(listed - actual)}"


def test_version_is_single_sourced():
    """__version__ must come from package metadata, not a hardcoded literal."""
    tree = ast.parse((SRC / "__init__.py").read_text())
    # The only permitted literal is the "unknown" fallback for a non-installed package.
    literals = [
        node.value.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(getattr(t, "id", None) == "__version__" for t in node.targets)
        and isinstance(node.value, ast.Constant)
        and node.value.value != "unknown"
    ]
    assert not literals, (
        f"__version__ is hardcoded to {literals}; it must come from importlib.metadata "
        "so pyproject.toml stays the single source of truth"
    )
    assert bk.__version__ and bk.__version__ != "unknown"


def test_no_source_file_uses_the_bk_alias():
    """Library code must not reference the `bk` notebook alias (NameError at runtime)."""
    offenders = []
    for f in SRC.rglob("*.py"):
        tree = ast.parse(f.read_text())
        defined = {"bk"} & {
            n.asname or n.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for n in node.names
        }
        if defined:
            continue
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "bk"
            ):
                offenders.append(f"{f.relative_to(ROOT)}:{node.lineno}")
    assert not offenders, f"`bk.` used without importing bullkpy: {offenders}"


def test_source_files_use_lf_line_endings():
    """Guard against the legacy CR-only line endings that broke io.py."""
    bad = [
        str(f.relative_to(ROOT))
        for f in SRC.rglob("*.py")
        if b"\r" in f.read_bytes()
    ]
    assert not bad, f"files with CR/CRLF line endings: {bad}"


@pytest.mark.parametrize("sub", SUBMODULES)
def test_every_export_has_a_docstring(sub):
    """A published API page with no docstring renders as an empty stub."""
    import inspect

    mod = getattr(bk, sub)
    undocumented = [
        n for n in _exports(sub)
        if callable(getattr(mod, n, None))
        and not (inspect.getdoc(getattr(mod, n)) or "").strip()
    ]
    assert not undocumented, f"exported without a docstring in bullkpy.{sub}: {undocumented}"


@pytest.mark.parametrize("sub", ["pp", "tl", "pl", "get"])
def test_no_dataset_specific_defaults(sub):
    """Defaults must not hardcode column names from one particular study."""
    import inspect

    smell = ("Project_ID", "Project ID", "OS.time", "PFS_6m", "NR_6m",
             "Neuroendocrine_score", "mp_heterogeneity_entropy", "NE25_score")
    # tcga_define_groups is TCGA-specific by name, so Project_ID is legitimate there
    allowed = {("tcga_define_groups", "project_key")}

    mod = getattr(bk, sub)
    offenders = []
    for name in _exports(sub):
        obj = getattr(mod, name, None)
        if not callable(obj):
            continue
        try:
            sig = inspect.signature(obj)
        except (ValueError, TypeError):
            continue
        for pname, p in sig.parameters.items():
            if p.default is inspect.Parameter.empty or (name, pname) in allowed:
                continue
            text = repr(p.default)
            if any(s in text for s in smell):
                offenders.append(f"{sub}.{name}({pname}={p.default!r})")

    assert not offenders, (
        "public defaults reference a specific dataset's columns; make them "
        f"required instead:\n  " + "\n  ".join(offenders)
    )
