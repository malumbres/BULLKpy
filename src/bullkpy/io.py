from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import pandas as pd
import numpy as np
import anndata as ad

from .logging import info, warn

__all__ = [
    "read_counts",
    "add_metadata",
    "read_metaprograms_xlsx",
    "save_signature_bank",
    "load_signature_bank",
]


def read_counts(
    filename: str | Path,
    *,
    sep: str = "\t",
    orientation: Literal["genes_by_samples", "samples_by_genes"] = "genes_by_samples",
    dtype: str | None = "int64",
) -> ad.AnnData:
    """
    Read a bulk RNA-seq count matrix and return an AnnData object.

    Parameters
    ----------
    filename
        Path to counts file (tsv/csv). Rows or columns must be gene IDs.
    sep
        Column separator (default: tab).
    orientation
        - "genes_by_samples": genes in rows, samples in columns (common output of HTSeq, featureCounts)
        - "samples_by_genes": samples in rows, genes in columns
    dtype
        Cast count matrix to dtype (default: int64). Set to None to disable casting.

    Returns
    -------
    AnnData
        AnnData object with samples in `.obs` and genes in `.var`.
    """
    filename = Path(filename)
    info(f"Reading count matrix from {filename}")

    df = pd.read_csv(filename, sep=sep, index_col=0)

    info(f"Raw matrix shape: {df.shape}")

    if orientation == "genes_by_samples":
        info("Interpreting rows as genes and columns as samples")
        df = df.T
    elif orientation == "samples_by_genes":
        info("Interpreting rows as samples and columns as genes")
    else:
        raise ValueError(
            "orientation must be 'genes_by_samples' or 'samples_by_genes'"
        )

    # Basic sanity checks
    if df.index.has_duplicates:
        warn("Sample names are duplicated")

    if df.columns.has_duplicates:
        warn("Gene names are duplicated")

    if dtype is not None:
        # Casting float expression values to an integer dtype truncates them
        # (10.77 -> 10, 0.9999 -> 0), which silently destroys the data. Several
        # public resources distribute log2(count+1) matrices under a "counts"
        # filename, so check before casting rather than after.
        is_integer_target = np.issubdtype(np.dtype(dtype), np.integer)
        values = df.to_numpy()
        looks_integral = (
            np.issubdtype(values.dtype, np.integer)
            or np.allclose(values[np.isfinite(values)],
                           np.round(values[np.isfinite(values)]), atol=1e-8)
            if values.size
            else True
        )

        if is_integer_target and not looks_integral:
            warn(
                f"{filename.name} does not contain integer values, so it was NOT cast "
                f"to {dtype}: doing so would truncate every value (e.g. 10.77 -> 10) "
                "and destroy the data. The matrix has been loaded unchanged as float.\n"
                "  If this file holds log2(count+1) values (as UCSC Xena distributes "
                "them), recover counts with: adata.X = 2 ** adata.X - 1\n"
                "  If it is already normalised, pass dtype=None and skip "
                "bk.pp.normalize_cpm()."
            )
        else:
            try:
                df = df.astype(dtype)
            except Exception as e:
                warn(f"Could not cast counts to {dtype}: {e}")

    # Warn if data does not look like counts
    if np.any(df.values < 0):
        warn("Count matrix contains negative values")

    adata = ad.AnnData(X=df)

    adata.obs_names = df.index.astype(str)
    adata.var_names = df.columns.astype(str)

    info(
        f"Created AnnData object with "
        f"{adata.n_obs} samples x {adata.n_vars} genes"
    )

    return adata


def add_metadata(
    adata: ad.AnnData,
    metadata_file: str | Path,
    *,
    index_col: str,
    sep: str = "\t",
    low_memory = False, 
    how: Literal["left", "inner"] = "left",
) -> ad.AnnData:
    """
    Add sample metadata to an AnnData object.

    Parameters
    ----------
    adata
        AnnData object with samples in `.obs`.
    metadata_file
        Path to metadata file (tsv, csv, or xlsx).
    index_col
        Column in metadata that matches `adata.obs_names`.
    sep
        Column separator for tsv/csv files.
    low_memory
        Pandas ``low_memory`` parameter. Only applies to delimited (tsv/csv)
        files; it is ignored for Excel input.
    how
        Merge strategy:
        - "left": keep all samples in adata (default)
        - "inner": keep only samples present in metadata

    Returns
    -------
    AnnData
        The same AnnData object with updated `.obs`.
    """
    metadata_file = Path(metadata_file)
    info(f"Adding metadata from {metadata_file}")

    # Load metadata
    if metadata_file.suffix in {".xls", ".xlsx"}:
        # NB: read_excel has no `low_memory` parameter.
        meta = pd.read_excel(metadata_file)
    else:
        meta = pd.read_csv(metadata_file, sep=sep, low_memory=low_memory)

    if index_col not in meta.columns:
        raise ValueError(
            f"index_col '{index_col}' not found in metadata columns"
        )

    meta = meta.set_index(index_col)

    # Sanity checks
    if meta.index.has_duplicates:
        warn("Metadata index contains duplicated sample IDs")

    missing = adata.obs_names.difference(meta.index)
    if len(missing) > 0:
        warn(
            f"{len(missing)} samples in AnnData are missing metadata "
            f"(showing up to 5): {list(missing[:5])}"
        )

    overlap = adata.obs_names.intersection(meta.index)
    info(f"Found metadata for {len(overlap)} / {adata.n_obs} samples")

    # Merge
    adata.obs = adata.obs.merge(
        meta,
        left_index=True,
        right_index=True,
        how=how,
    )
    # A wide table arrives with one internal block per column; later inserts into it
    # (QC metrics, group labels) would each warn about fragmentation. One copy fixes it.
    adata.obs = adata.obs.copy()

    return adata


def save_signature_bank(sig_bank: dict, outdir: str):
    """
    Write a bank of signatures to disk.

    Each signature becomes a weights table plus a JSON sidecar holding its
    intercept and fit parameters, so a bank can be version-controlled and
    reloaded without refitting.

    Parameters
    ----------
    sig_bank
        Mapping of name to signature dict (``weights_df``, ``intercept``, ``params``).
    outdir
        Destination directory; created if absent.

    Returns
    -------
    None

    See Also
    --------
    load_signature_bank : read a bank back.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = {"signatures": {}}

    for name, sig in sig_bank.items():
        # save weights as CSV (clean + easy diff)
        w_path = outdir / f"{name}.weights.tsv"
        sig["weights_df"].to_csv(w_path, sep="\t", index=False)

        # store everything else in JSON-friendly form
        meta["signatures"][name] = {
            "intercept": float(sig.get("intercept", 0.0)),
            "params": sig.get("params", {}),
            "genes_requested": sig.get("genes_requested", sig.get("genes", [])),
            "genes_used": sig.get("genes_used", sig.get("genes", [])),
            "panel_nonzero": sig.get("panel_nonzero", []),
            "n_nonzero": int(sig.get("n_nonzero", 0)),
            # optional quick metrics snapshots (if present)
            "metrics_train": sig.get("metrics_train", None),
            "metrics_ext": sig.get("metrics_ext", None),
            "weights_tsv": w_path.name,
        }

    with open(outdir / "sig_bank.json", "w") as f:
        json.dump(meta, f, indent=2)

    info(f"Saved {len(sig_bank)} signatures to: {outdir}")
def load_signature_bank(indir: str):
    """
    Load a bank of signatures previously written by :func:`save_signature_bank`.

    Parameters
    ----------
    indir
        Directory containing the saved weights tables and their JSON sidecars.

    Returns
    -------
    dict
        Mapping of name to signature dict, ready to pass to
        :func:`evaluate_signatures_store` or :func:`benchmark_signatures_pr_auc`.
    """
    indir = Path(indir)
    with open(indir / "sig_bank.json", "r") as f:
        meta = json.load(f)

    sig_bank = {}
    for name, rec in meta["signatures"].items():
        w = pd.read_csv(indir / rec["weights_tsv"], sep="\t")
        sig_bank[name] = {
            "name": name,
            "weights_df": w,
            "intercept": float(rec.get("intercept", 0.0)),
            "params": rec.get("params", {}),
            "genes_requested": rec.get("genes_requested", []),
            "genes_used": rec.get("genes_used", []),
            "panel_nonzero": rec.get("panel_nonzero", []),
            "n_nonzero": int(rec.get("n_nonzero", 0)),
            "metrics_train": rec.get("metrics_train", None),
            "metrics_ext": rec.get("metrics_ext", None),
        }
    return sig_bank



def read_metaprograms_xlsx(
    xlsx_path: str,
    *,
    sheet_name: str = "Malignant",
    min_genes: int = 5,
) -> dict[str, list[str]]:
    """
    Excel format expected: each column is an MP, rows are gene symbols.
    Returns dict mp_name -> list of genes (non-empty strings).
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

    mp = {}
    for col in df.columns:
        genes = (
            df[col]
            .dropna()
            .astype(str)
            .str.strip()
            .replace({"": np.nan, "nan": np.nan, "None": np.nan})
            .dropna()
            .tolist()
        )
        genes = [g for g in genes if g and g.lower() != "nan"]
        if len(genes) >= int(min_genes):
            mp[str(col)] = genes

    if len(mp) == 0:
        raise ValueError(f"No metaprograms found in sheet='{sheet_name}' of {xlsx_path}")

    return mp
