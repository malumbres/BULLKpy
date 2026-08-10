"""Docs must not reference APIs or image files that do not exist."""
from __future__ import annotations

import pathlib
import re

import bullkpy as bk

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
SUBMODULES = ["pp", "tl", "pl", "get", "io"]

_API_REF = re.compile(r"\bbk\.(pp|tl|pl|get|io)\.([A-Za-z_][A-Za-z0-9_]*)")
_IMG_REF = re.compile(r"_static/([A-Za-z0-9_.\-]+)")


def _markdown_files():
    return [p for p in DOCS.rglob("*.md") if "_build" not in p.parts]


def _public_api() -> set[str]:
    api = set()
    for sub in SUBMODULES:
        for name in getattr(getattr(bk, sub), "__all__", []):
            api.add(f"{sub}.{name}")
    return api


def test_docs_reference_only_existing_functions():
    api = _public_api()
    offenders = []
    for md in _markdown_files():
        for m in _API_REF.finditer(md.read_text()):
            ref = f"{m.group(1)}.{m.group(2)}"
            if ref not in api:
                offenders.append(f"{md.relative_to(ROOT)}: bk.{ref}")
    assert not offenders, "docs reference non-existent API:\n" + "\n".join(sorted(set(offenders)))


def test_docs_reference_only_existing_images():
    available = {p.name for p in (DOCS / "_static").glob("*")}
    offenders = []
    for md in _markdown_files():
        for m in _IMG_REF.finditer(md.read_text()):
            if m.group(1) not in available:
                offenders.append(f"{md.relative_to(ROOT)}: _static/{m.group(1)}")
    assert not offenders, "docs reference missing images:\n" + "\n".join(sorted(set(offenders)))


def test_no_private_helpers_in_doc_examples():
    """Examples must show the public API, not underscore-prefixed internals."""
    offenders = []
    for md in _markdown_files():
        for m in re.finditer(r"\bbk\.(?:pp|tl|pl|get|io)\._[A-Za-z0-9_]*", md.read_text()):
            offenders.append(f"{md.relative_to(ROOT)}: {m.group(0)}")
    assert not offenders, "docs use private helpers:\n" + "\n".join(sorted(set(offenders)))
