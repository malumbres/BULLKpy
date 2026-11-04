"""What is committed for a release: tutorial notebooks, versions and dataset checksums.

The notebook checks read the notebooks as committed (``git show HEAD:...``), not the working copies, because the
working copies hold full-resolution outputs and local paths and are trimmed by ``scripts/prepare_notebooks.py``
when staged for commit.
"""
from __future__ import annotations

import base64
import io
import json
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
TUTORIALS = ["BULLKpy_nb1", "BULLKpy_nb2", "BULLKpy_nb3"]
MAX_NOTEBOOK_MB = 8.0
MAX_PNG_WIDTH = 1400


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True, text=True).stdout


def _committed_notebooks() -> list[str]:
    try:
        out = _git("ls-tree", "-r", "--name-only", "HEAD", "--", "notebooks")
    except (OSError, subprocess.CalledProcessError):
        return []
    return [p for p in out.splitlines() if p.endswith(".ipynb")]


COMMITTED = _committed_notebooks()
needs_git = pytest.mark.skipif(not COMMITTED, reason="needs a git checkout with committed notebooks")


def _load(path: str) -> tuple[str, dict]:
    text = _git("show", f"HEAD:{path}")
    return text, json.loads(text)


@needs_git
def test_the_three_tutorials_are_committed():
    assert {pathlib.Path(p).stem for p in COMMITTED} >= set(TUTORIALS)


@needs_git
@pytest.mark.parametrize("path", COMMITTED, ids=lambda p: pathlib.Path(p).name)
def test_committed_notebook_is_small_clean_and_fully_run(path):
    text, nb = _load(path)
    assert len(text) / 1e6 <= MAX_NOTEBOOK_MB, f"{path} is {len(text) / 1e6:.1f} MB; run scripts/prepare_notebooks.py"
    assert "/Users/" not in text and "\\Users\\" not in text, "local paths in the committed notebook"

    from PIL import Image

    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    for cell in code:
        for out in cell.get("outputs", []):
            assert out["output_type"] != "error", "an error output is committed"
            assert not (out["output_type"] == "stream" and out.get("name") == "stderr"), "stderr (warnings) committed"
            if out["output_type"] in ("display_data", "execute_result") and "image/png" in out["data"]:
                png = out["data"]["image/png"]
                png = "".join(png) if isinstance(png, list) else png
                width = Image.open(io.BytesIO(base64.b64decode(png))).width
                assert width <= MAX_PNG_WIDTH, f"a {width}px figure is committed; run scripts/prepare_notebooks.py"
    unrun = [i for i, c in enumerate(code) if "".join(c["source"]).strip() and c.get("execution_count") is None]
    assert not unrun, f"code cells that were never run: {unrun}"
    assert not any("importlib.reload" in "".join(c["source"]) for c in code), "development reload cell left in"


@needs_git
def test_docs_tutorial_index_lists_exactly_the_committed_notebooks():
    index = (ROOT / "docs" / "notebooks" / "index.md").read_text()
    block = index.split("```{toctree}", 1)[1].split("```", 1)[0]
    listed = {line.strip().split("<")[-1].rstrip(">").strip() for line in block.splitlines() if line.strip() and not line.startswith(":")}
    assert listed == {pathlib.Path(p).stem for p in COMMITTED}


def test_version_is_the_same_in_pyproject_citation_and_changelog():
    import tomllib

    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    cff = re.search(r"^version:\s*(\S+)", (ROOT / "CITATION.cff").read_text(), re.M).group(1)
    first_release = re.search(r"^## \[(\d[^\]]*)\]", (ROOT / "CHANGELOG.md").read_text(), re.M).group(1)
    assert version == cff == first_release, (version, cff, first_release)


def test_dataset_checksums_match_the_local_manifest_when_there_is_one():
    manifest = ROOT / "data" / "zenodo_tcga_lung" / "TCGA_lung_MANIFEST.tsv"
    if not manifest.exists():
        pytest.skip("no local copy of the tutorial data")
    from bullkpy import datasets

    rows = dict(line.split("\t")[::2] for line in manifest.read_text().splitlines()[1:] if line)
    for name, digest in datasets._FILES.items():
        assert rows.get(name) == digest, f"{name}: update bullkpy.datasets._FILES from the manifest"


def test_extras_kept_for_backwards_compatibility_exist():
    import tomllib

    extras = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["optional-dependencies"]
    assert {"umap", "leiden", "network"} <= set(extras)
    assert not extras["umap"] and not extras["leiden"] and not extras["network"], "these extras should install nothing"
