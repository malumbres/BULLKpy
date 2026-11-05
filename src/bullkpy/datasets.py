"""Small public datasets used by the tutorials."""
from __future__ import annotations

import hashlib
import os
import shutil
import urllib.request
from pathlib import Path

import anndata as ad
import pandas as pd

from .io import add_metadata, read_counts
from .logging import info, warn

__all__ = ["tcga_lung"]

#: Zenodo record holding the tutorial files. Set when the record is published; until then
#: ``tcga_lung()`` works from a local folder (``path=``) but cannot download.
ZENODO_RECORD: str | None = None

#: SHA-256 of the files in the record (written by ``scripts/make_tcga_lung_subset.py``).
_FILES = {
    "TCGA_lung_counts.tsv.gz": "6c2a2da2bd35fe9339beaabfc1b2fa2a777222efcb0b6acfb9516f28b658a68f",
    "TCGA_lung_genes.tsv": "74b16243d257c47ca6665c00add54057c47939d03db92db39383d312f2cead0b",
    "TCGA_lung_metadata.tsv": "a32c03cf525ff867bb5f102ad8745c1ff05dc0a9824a655678dcf8686f41ac49",
}

_SAMPLE_TYPE = {"01": "Primary", "02": "Recurrent", "11": "Normal"}


def _cache_dir() -> Path:
    base = os.environ.get("BULLKPY_DATA") or (
        Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "bullkpy"
    )
    return Path(base) / "tcga_lung"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` through a temporary file, so a failed download leaves nothing."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url) as response, open(tmp, "wb") as out:  # noqa: S310
            shutil.copyfileobj(response, out)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def tcga_lung(path: str | Path | None = None, *, primary_only: bool = False) -> ad.AnnData:
    """
    TCGA lung cancer RNA-seq (LUAD and LUSC) as an AnnData object with integer counts.

    The data are the open-access TCGA samples from the UCSC Xena GDC Pan-Cancer hub. Xena
    distributes them as log2(count + 1); integer counts were recovered as ``2**x - 1`` (see
    ``scripts/make_tcga_lung_subset.py``). All values were integers to within floating-point precision,
    except for 10 samples whose values were half-integers; those were rounded, and their barcodes are in
    ``TCGA_lung_half_integer_samples.txt`` in the data record.

    The object has 1,108 samples (567 LUAD, 541 LUSC) and about 58,000 genes, with Ensembl IDs in
    ``.var["ensembl_id"]``, gene symbols as ``var_names`` (the Ensembl ID where there is no symbol),
    the clinical annotation from the TCGA Pan-Cancer Clinical Data Resource in ``.obs``, and
    ``.obs["Sample_type"]`` (Primary, Recurrent or Normal) derived from the sample barcode.

    Parameters
    ----------
    path
        Folder containing ``TCGA_lung_counts.tsv.gz``, ``TCGA_lung_genes.tsv`` and
        ``TCGA_lung_metadata.tsv``. By default the files are downloaded once from Zenodo into a cache
        folder (``~/.cache/bullkpy``, or ``$BULLKPY_DATA``) and checked against their SHA-256.
    primary_only
        Keep only the primary tumours (drops normal and recurrent samples).
    """
    if path is None:
        folder = _cache_dir()
        for name, expected in _FILES.items():
            target = folder / name
            if target.exists() and _sha256(target) == expected:
                continue
            if ZENODO_RECORD is None:
                raise RuntimeError(
                    "The tutorial data record on Zenodo is not set in this version. Download the three "
                    "TCGA_lung_* files manually and pass the folder as tcga_lung(path=...)."
                )
            url = f"https://zenodo.org/records/{ZENODO_RECORD}/files/{name}?download=1"
            info(f"Downloading {name} from Zenodo (record {ZENODO_RECORD}) ...")
            _download(url, target)
            if _sha256(target) != expected:
                target.unlink(missing_ok=True)
                raise RuntimeError(f"Checksum mismatch for {name}; the download was discarded.")
    else:
        folder = Path(path)
        missing = [n for n in _FILES if not (folder / n).exists()]
        if missing:
            raise FileNotFoundError(f"{folder} is missing: {', '.join(missing)}")
        changed = [n for n, expected in _FILES.items() if _sha256(folder / n) != expected]
        if changed:
            warn(f"{', '.join(changed)} differ from the published files (different checksum).")

    counts, genes, meta = (folder / n for n in _FILES)
    adata = read_counts(counts, sep="\t", orientation="genes_by_samples")

    gene_table = pd.read_csv(genes, sep="\t", index_col="ensembl_id")
    adata.var = adata.var.join(gene_table)
    adata.var["ensembl_id"] = adata.var_names.values
    names = adata.var["gene_name"].fillna(adata.var["ensembl_id"]).astype(str)
    adata.var_names = names.values
    adata.var_names_make_unique()

    add_metadata(adata, meta, sep="\t", index_col="sample", low_memory=False)
    adata.obs["Sample_type"] = adata.obs_names.str[13:15].map(_SAMPLE_TYPE)

    if primary_only:
        adata = adata[adata.obs["Sample_type"] == "Primary"].copy()
    return adata
