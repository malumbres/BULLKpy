from __future__ import annotations

from pathlib import Path
from typing import Literal

import pandas as pd
import numpy as np
import anndata as ad

from .logging import info, warn

__all__ = ["read_counts", "add_metadata"]


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

    return adata