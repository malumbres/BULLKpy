from __future__ import annotations

from typing import Literal, Sequence
import numpy as np
import pandas as pd
import anndata as ad


def rank_genes_groups_df_all(
    adata: ad.AnnData,
    *,
    key: str = "rank_genes_groups",
    groups: Sequence[str] | None = None,
    sort_by: Literal["scores", "logfoldchanges", "pvals_adj", "pvals"] = "scores",
    ascending: bool | None = None,
) -> pd.DataFrame:
    """
    Return rank_genes_groups results for all (or selected) groups as one tidy DataFrame.

    Output columns:
      group, gene, scores, log2FC, pval, qval, mean_group?, mean_ref?
    """
    if key not in adata.uns:
        raise KeyError(f"adata.uns['{key}'] not found. Run bk.tl.rank_genes_groups first.")

    rg = adata.uns[key]
    for needed in ("names", "scores", "logfoldchanges", "pvals", "pvals_adj"):
        if needed not in rg:
            raise KeyError(f"adata.uns['{key}'] missing '{needed}'")

    available = list(rg["names"].keys())
    if groups is None:
        groups = available
    else:
        groups = [str(g) for g in groups]
        missing = [g for g in groups if g not in set(available)]
        if missing:
            raise KeyError(f"Requested groups not in results: {missing}. Available: {available}")

    rows = []
    for g in groups:
        names = np.asarray(rg["names"][g], dtype=str)
        scores = np.asarray(rg["scores"][g], dtype=float)
        logfc = np.asarray(rg["logfoldchanges"][g], dtype=float)
        pvals = np.asarray(rg["pvals"][g], dtype=float)
        qvals = np.asarray(rg["pvals_adj"][g], dtype=float)

        df = pd.DataFrame(
            {
                "group": g,
                "gene": names,
                "scores": scores,
                "log2FC": logfc,
                "pval": pvals,
                "qval": qvals,
            }
        )

        if "mean_group" in rg and g in rg["mean_group"]:
            df["mean_group"] = np.asarray(rg["mean_group"][g], dtype=float)
        if "mean_ref" in rg and g in rg["mean_ref"]:
            df["mean_ref"] = np.asarray(rg["mean_ref"][g], dtype=float)

        rows.append(df)

    out = pd.concat(rows, axis=0, ignore_index=True)

    if ascending is None:
        ascending = sort_by in ("pvals", "pvals_adj")

    sort_map = {
        "scores": "scores",
        "logfoldchanges": "log2FC",
        "pvals": "pval",
        "pvals_adj": "qval",
    }
    out = out.sort_values(["group", sort_map[sort_by]], ascending=[True, bool(ascending)]).reset_index(drop=True)
    return out