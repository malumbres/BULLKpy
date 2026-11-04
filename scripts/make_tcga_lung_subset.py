#!/usr/bin/env python
"""Build the TCGA lung (LUAD + LUSC) starting material for the BULLKpy tutorials.

UCSC Xena distributes the GDC Pan-Cancer table as log2(count + 1) under the name
``GDC-PANCAN.htseq_counts.tsv``. This script recovers integer counts for lung
samples, so the tutorials can start from true counts with ``bk.io.read_counts``.

Outputs in ``--outdir``:
    TCGA_lung_counts.tsv.gz    integer counts, genes x samples, Ensembl IDs
    TCGA_lung_genes.tsv        Ensembl ID -> gene name, gene type, chromosome
    TCGA_lung_metadata.tsv     sample metadata, restricted to the columns in --columns
    TCGA_lung_half_integer_samples.txt   barcodes of samples whose values were half-integers
    TCGA_lung_MANIFEST.tsv     sizes and SHA-256 checksums

The script stops if back-transformed values are not whole or half-integers within
``--max-residual``, because that would mean the input is not log2(count + 1).

A few samples hold exact half-integers throughout (possibly two sequencing files averaged by
Xena; not verified). They are rounded half-to-even, which avoids a directional bias, and
listed in ``TCGA_lung_half_integer_samples.txt`` so they can be excluded.

The source metadata file also holds analysis results saved by an earlier pipeline (clusters,
scores, cell-cycle phase) and constant columns. ``--columns`` names the ones that are kept; with
``--manifest-only`` only the checksum manifest of an existing output folder is rewritten.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

PROJECTS = ["LUAD", "LUSC"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--counts", type=Path, default=Path("data/external/GDC-PANCAN.htseq_counts.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/external/260127_TCGA_metadata.tsv"))
    parser.add_argument("--biomart", type=Path, default=Path("data/external/251224_biomart_human.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("data/external/lung"))
    parser.add_argument("--max-residual", type=float, default=0.01)
    parser.add_argument("--columns", type=Path, default=Path("scripts/tcga_lung_metadata_columns.txt"))
    parser.add_argument("--manifest-only", action="store_true", help="only rewrite the manifest of --outdir")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.manifest_only:
        _write_manifest(args.outdir)
        return

    meta = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    meta = meta[meta["Project_ID"].isin(PROJECTS)]

    header = pd.read_csv(args.counts, sep="\t", nrows=0).columns
    id_col = header[0]
    samples = [s for s in header[1:] if s in set(meta["sample"])]
    print(f"{len(meta)} lung samples in metadata, {len(samples)} with expression")

    log_values = pd.read_csv(
        args.counts,
        sep="\t",
        usecols=[id_col, *samples],
        index_col=id_col,
        dtype={s: "float64" for s in samples},
    )[samples]

    raw = np.power(2.0, log_values.to_numpy()) - 1.0
    rounded = np.rint(raw)
    # Doubling turns whole and half-integers into integers, so this tests for both.
    residual = float(np.abs(2.0 * raw - np.rint(2.0 * raw)).max()) / 2.0
    print(f"max distance from a whole or half-integer: {residual:.2e}")
    if residual > args.max_residual:
        raise SystemExit(f"Values are not log2(count + 1): distance {residual} > {args.max_residual}")
    if rounded.min() < 0:
        raise SystemExit("Negative counts after back-transform")

    has_half = np.abs(raw - rounded).max(axis=0) > args.max_residual
    half_samples = [s for s, flagged in zip(samples, has_half) if flagged]
    print(f"samples with half-integer values: {len(half_samples)} {half_samples}")

    counts = pd.DataFrame(rounded.astype("int64"), index=log_values.index, columns=samples)
    n_all = len(counts)
    # All-zero genes do not affect library sizes, so dropping them loses nothing.
    counts = counts[counts.sum(axis=1) > 0]
    print(f"genes: {n_all} total, {len(counts)} non-zero in lung")
    lib = counts.sum(axis=0)
    print(f"library size: min {lib.min():,}  median {int(lib.median()):,}  max {lib.max():,}")

    counts_path = args.outdir / "TCGA_lung_counts.tsv.gz"
    counts.to_csv(counts_path, sep="\t", compression={"method": "gzip", "mtime": 0})

    biomart = pd.read_csv(
        args.biomart,
        sep="\t",
        usecols=["Gene stable ID", "Gene name", "Gene type", "Chromosome/scaffold name"],
        low_memory=False,
    ).drop_duplicates("Gene stable ID")
    genes = pd.DataFrame({"ensembl_id": counts.index, "stable_id": counts.index.str.split(".").str[0]})
    genes = genes.merge(biomart, how="left", left_on="stable_id", right_on="Gene stable ID")
    genes = genes.drop(columns=["stable_id", "Gene stable ID"]).rename(
        columns={"Gene name": "gene_name", "Gene type": "gene_type", "Chromosome/scaffold name": "chromosome"}
    )
    print(f"genes with a symbol: {genes['gene_name'].notna().sum()} of {len(genes)}")
    genes_path = args.outdir / "TCGA_lung_genes.tsv"
    genes.to_csv(genes_path, sep="\t", index=False)

    meta = meta[meta["sample"].isin(samples)].dropna(axis=1, how="all")
    meta = meta.set_index("sample").loc[samples].reset_index()
    if args.columns.exists():
        wanted = [
            line.strip() for line in args.columns.read_text().splitlines() if line.strip() and not line.startswith("#")
        ]
        absent = [c for c in wanted if c not in meta.columns]
        if absent:
            print(f"WARNING: {len(absent)} columns in {args.columns} are not in the metadata: {absent[:5]}")
        meta = meta[["sample", *[c for c in wanted if c in meta.columns]]]
    (args.outdir / "TCGA_lung_half_integer_samples.txt").write_text("\n".join(half_samples) + "\n")
    print(f"metadata: {meta.shape[0]} samples, {meta.shape[1]} columns")
    print(meta["Project_ID"].value_counts().to_dict())
    meta_path = args.outdir / "TCGA_lung_metadata.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    _write_manifest(args.outdir)


def _write_manifest(outdir: Path) -> None:
    names = [
        "TCGA_lung_counts.tsv.gz",
        "TCGA_lung_genes.tsv",
        "TCGA_lung_metadata.tsv",
        "TCGA_lung_half_integer_samples.txt",
    ]
    rows = [
        {"file": n, "bytes": (outdir / n).stat().st_size, "sha256": sha256(outdir / n)}
        for n in names
        if (outdir / n).exists()
    ]
    manifest = pd.DataFrame(rows)
    manifest.to_csv(outdir / "TCGA_lung_MANIFEST.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
