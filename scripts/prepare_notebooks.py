#!/usr/bin/env python
"""Make executed tutorial notebooks small and clean enough to commit.

Read the Docs renders the stored outputs (it does not execute the notebooks), so the committed copies
need their figures, but not at 300 dpi, and not the warnings or local paths of the machine that ran them.

For every notebook given, a trimmed copy is produced:

* PNG figures wider than ``--max-width`` pixels are scaled down (default 1000, a little wider than a
  documentation page) and stored with a 256-colour palette, which is visually identical for plots and
  3-5 times smaller (``--no-quantize`` keeps full colour);
* ``stderr`` streams (warnings, log lines) are dropped;
* editor-specific output formats that no renderer shows (``application/vnd.*``) are dropped, and an
  HTML table larger than 150 kB is replaced by its plain-text form;
* absolute paths (``/Users/...``, ``C:\\Users\\...``) in text outputs are replaced by ``<path>``;
* ``metadata.widgets`` is removed;
* short code cells that only reload modules (``importlib.reload``) are removed: they are development scaffolding.

Code, markdown and execution counts are untouched, and so are error outputs (there should be none).

Usage::

    python scripts/prepare_notebooks.py notebooks/BULLKpy_nb*.ipynb --stage   # trim and stage for commit
    python scripts/prepare_notebooks.py notebooks/BULLKpy_nb1.ipynb --out /tmp/clean

``--stage`` writes the trimmed notebook into the git index only; the file in your working folder keeps its
full-resolution outputs. Run it again after editing a notebook.
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image

KEEP_MIME = {"text/plain", "text/html", "image/png", "image/svg+xml", "text/markdown", "text/latex", "application/json"}
PATH_RE = re.compile(r"(?:/Users/|/home/)[^\s'\"<>)\]]+|[A-Za-z]:\\Users\\[^\s'\"<>)\]]+")
HTML_LIMIT = 150_000


def _join(value) -> str:
    return "".join(value) if isinstance(value, list) else str(value)


def _shrink_png(b64: str, max_width: int, quantize: bool) -> tuple[str, tuple[int, int]]:
    raw = base64.b64decode(b64)
    img = Image.open(io.BytesIO(raw))
    resized = img.width > max_width
    if not resized and not quantize:
        return b64, img.size  # nothing to gain; keep matplotlib's own bytes
    if resized:
        img = img.resize((max_width, round(img.height * max_width / img.width)), Image.LANCZOS)
    if img.mode in ("RGBA", "LA"):
        flat = Image.new("RGB", img.size, "white")
        flat.paste(img, mask=img.getchannel("A"))
        img = flat
    else:
        img = img.convert("RGB")
    if quantize:
        img = img.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    out = io.BytesIO()
    img.save(out, format="PNG", optimize=True)
    return base64.b64encode(out.getvalue()).decode("ascii"), img.size


def _clean_text(value):
    text = PATH_RE.sub("<path>", _join(value))
    return text


def _is_reload_cell(cell: dict) -> bool:
    source = _join(cell.get("source", ""))
    return cell.get("cell_type") == "code" and "importlib.reload" in source and len(source.splitlines()) <= 6


def trim_notebook(nb: dict, *, max_width: int = 1000, quantize: bool = True) -> dict:
    nb.get("metadata", {}).pop("widgets", None)
    nb["cells"] = [c for c in nb.get("cells", []) if not _is_reload_cell(c)]
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        outputs = []
        for out in cell.get("outputs", []):
            kind = out.get("output_type")
            if kind == "stream":
                if out.get("name") == "stderr":
                    continue
                out["text"] = _clean_text(out.get("text", ""))
            elif kind in ("display_data", "execute_result"):
                data = {k: v for k, v in out.get("data", {}).items() if k in KEEP_MIME}
                if "text/html" in data and "text/plain" in data and len(_join(data["text/html"])) > HTML_LIMIT:
                    del data["text/html"]
                for mime in ("text/plain", "text/html"):
                    if mime in data:
                        data[mime] = _clean_text(data[mime])
                if "image/png" in data:
                    data["image/png"], size = _shrink_png(_join(data["image/png"]), max_width, quantize)
                    meta = out.setdefault("metadata", {})
                    meta["image/png"] = {"width": size[0], "height": size[1]}
                if not data:
                    continue
                out["data"] = data
            elif kind == "error":
                out["traceback"] = [_clean_text(t) for t in out.get("traceback", [])]
            outputs.append(out)
        cell["outputs"] = outputs
    return nb


def _dump(nb: dict) -> str:
    return json.dumps(nb, indent=1, ensure_ascii=False) + "\n"


def _stage(path: Path, text: str) -> None:
    blob = subprocess.run(
        ["git", "hash-object", "-w", "--stdin"], input=text.encode(), capture_output=True, check=True
    ).stdout.decode().strip()
    rel = subprocess.run(
        ["git", "ls-files", "--full-name", "--others", "--cached", str(path)], capture_output=True, check=True, text=True
    ).stdout.split("\n")[0] or str(path)
    subprocess.run(["git", "update-index", "--add", "--cacheinfo", f"100644,{blob},{rel}"], check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("notebooks", nargs="+", type=Path)
    parser.add_argument("--max-width", type=int, default=1000)
    parser.add_argument("--no-quantize", action="store_true", help="keep full-colour PNGs")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--stage", action="store_true", help="stage trimmed copies in the git index")
    group.add_argument("--out", type=Path, help="write trimmed copies to this folder")
    args = parser.parse_args()

    for path in args.notebooks:
        original = path.read_text()
        trimmed = _dump(trim_notebook(json.loads(original), max_width=args.max_width, quantize=not args.no_quantize))
        if args.stage:
            _stage(path, trimmed)
            target = "git index"
        else:
            args.out.mkdir(parents=True, exist_ok=True)
            (args.out / path.name).write_text(trimmed)
            target = str(args.out / path.name)
        print(f"{path.name}: {len(original) / 1e6:.1f} MB -> {len(trimmed) / 1e6:.1f} MB ({target})")


if __name__ == "__main__":
    sys.exit(main())
