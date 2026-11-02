"""Gene-name lookups that fail clearly when ``adata.var_names`` is not unique.

A pandas ``Index`` with repeated labels breaks name-based lookup in confusing
ways: ``get_indexer`` raises ``Reindexing only valid with uniquely valued Index
objects`` without saying which names are the problem, and ``get_loc`` quietly
returns a slice or boolean mask where callers expect one integer.

Duplicates are common in bulk data (several Ensembl IDs sharing a symbol such as
``Y_RNA``, or IDs with no symbol that become ``NaN``). A lookup only needs to fail
when a *requested* name is ambiguous; every other gene still resolves normally.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["gene_loc", "gene_indexer"]

_HINT = (
    "Make the names unique before analysis, e.g. "
    "`adata.var_names_make_unique()`; drop unnamed genes with "
    "`adata = adata[:, adata.var_names.notna()].copy()`; or use stable IDs "
    "(Ensembl) as var_names."
)


def _preview(names, k: int = 8) -> str:
    names = [str(n) for n in names]
    more = f", ... (+{len(names) - k} more)" if len(names) > k else ""
    return ", ".join(names[:k]) + more


def _summary(vn: pd.Index) -> str:
    n_nan = int(vn.isna().sum())
    dup = vn[vn.notna() & vn.duplicated(keep=False)]
    parts = []
    if n_nan:
        parts.append(f"{n_nan} missing (NaN) names")
    if len(dup):
        parts.append(f"{len(dup)} rows sharing {dup.nunique()} duplicated names")
    return "adata.var_names is not unique: " + " and ".join(parts) + "."


def gene_loc(adata, gene) -> int:
    """Integer position of ``gene`` in ``adata.var_names``.

    Raises ``KeyError`` if absent (as ``Index.get_loc`` does) and ``ValueError``
    if the name occurs more than once.
    """
    vn = adata.var_names
    loc = vn.get_loc(gene)
    if isinstance(loc, (slice, np.ndarray)):
        raise ValueError(
            f"Gene {gene!r} appears more than once in adata.var_names, so it is "
            f"ambiguous which column to use. {_summary(vn)} {_HINT}"
        )
    return loc


def gene_indexer(adata, genes) -> np.ndarray:
    """Positions of ``genes`` in ``adata.var_names``; ``-1`` where absent.

    Same contract as ``Index.get_indexer``, but it also works when ``var_names``
    holds duplicates, as long as none of the *requested* genes is among them.
    """
    vn = adata.var_names
    if vn.is_unique:
        return vn.get_indexer(genes)

    genes = pd.Index(genes)
    ambiguous = vn[vn.notna() & vn.duplicated(keep=False)].unique()
    bad = [g for g in genes.unique() if g in ambiguous]
    if bad:
        raise ValueError(
            f"{len(bad)} requested gene(s) appear more than once in "
            f"adata.var_names: {_preview(bad)}. {_summary(vn)} {_HINT}"
        )

    keep = ~vn.duplicated(keep=False)
    unique_part = pd.Series(np.flatnonzero(keep), index=vn[keep])
    return unique_part.reindex(genes).fillna(-1).astype(int).to_numpy()
