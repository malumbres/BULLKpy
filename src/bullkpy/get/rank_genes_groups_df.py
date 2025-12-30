from __future__ import annotations

from typing import Literal
import numpy as np
import pandas as pd
import anndata as ad


def rank_genes_groups_df(
    adata: ad.AnnData,
    *,
    group: str,
    key: str = "rank_genes_groups",
    sort_by: Literal["scores", "logfoldchanges", "pvals_adj", "pvals"] = "scores",
    ascending: bool | None = None,
) -> pd.DataFrame:
    """
    Return rank_genes_groups results for one group as a tidy DataFrame.

    Expects `bk.tl.rank_genes_groups(...)` results in `adata.uns[key]`.

    Columns:
      gene, scores, log2FC, pval, qval, mean_group, mean_ref
    """
    if key not in adata.uns:
        raise KeyError(f"adata.uns['{key}'] not found. Run bk.tl.rank_genes_groups first.")

    rg = adata.uns[key]
    for needed in ("names", "scores", "logfoldchanges", "pvals", "pvals_adj"):
        if needed not in rg:
            raise KeyError(f"adata.uns['{key}'] missing '{needed}'")

    group = str(group)
    if group not in rg["names"]:
        available = list(rg["names"].keys())
        raise KeyError(f"group='{group}' not found in adata.uns['{key}']['names']. Available: {available}")

    names = np.asarray(rg["names"][group], dtype=str)
    scores = np.asarray(rg["scores"][group], dtype=float)
    logfc = np.asarray(rg["logfoldchanges"][group], dtype=float)
    pvals = np.asarray(rg["pvals"][group], dtype=float)
    qvals = np.asarray(rg["pvals_adj"][group], dtype=float)

    mean_g = None
    mean_r = None
    if "mean_group" in rg and group in rg["mean_group"]:
        mean_g = np.asarray(rg["mean_group"][group], dtype=float)
    if "mean_ref" in rg and group in rg["mean_ref"]:
        mean_r = np.asarray(rg["mean_ref"][group], dtype=float)

    df = pd.DataFrame(
        {
            "gene": names,
            "scores": scores,
            "log2FC": logfc,
            "pval": pvals,
            "qval": qvals,
        }
    )
    if mean_g is not None:
        df["mean_group"] = mean_g
    if mean_r is not None:
        df["mean_ref"] = mean_r

    if ascending is None:
        ascending = sort_by in ("pvals", "pvals_adj")

    sort_map = {
        "scores": "scores",
        "logfoldchanges": "log2FC",
        "pvals": "pval",
        "pvals_adj": "qval",
    }
    df = df.sort_values(sort_map[sort_by], ascending=bool(ascending)).reset_index(drop=True)
    return df