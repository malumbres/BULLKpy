from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal, Sequence

import numpy as np
import pandas as pd
import anndata as ad

from ..logging import info, warn


def pca_loadings(
    adata: ad.AnnData,
    *,
    key: str = "pca",
    loadings_key: str = "PCs",
    pcs: Sequence[int] | None = None,
    n_top: int = 50,
    use_abs: bool = False,
    include_negative: bool = True,
    dropna: bool = True,
    gene_col: str = "gene",
    store_key: str = "pca_loadings",
    export: str | Path | None = None,
    export_sep: Literal["\t", ","] = "\t",
    export_gmt: str | Path | None = None,
    min_abs_loading: float | None = None,
) -> dict[str, pd.DataFrame]:
    """
    Rank genes by PCA loadings per PC and (optionally) export for enrichment.

    Assumes PCA results like:
      - adata.varm[loadings_key] == array (n_vars x n_comps)
      - adata.uns[key] contains PCA params (optional)

    Parameters
    ----------
    pcs
        PCs to process, 1-based indexing (e.g. [1,2,3]). If None -> all available.
    n_top
        Number of top genes to return per PC (and per sign if include_negative=True and use_abs=False).
    use_abs
        If True: rank by absolute loading (single list per PC).
        If False: returns positive top list (and negative list if include_negative=True).
    include_negative
        Only used when use_abs=False. If True returns both pos/neg lists.
    export
        If provided, writes a long-format table (PC, sign, gene, loading, rank).
        Also writes a wide-format file alongside: "<stem>.wide<suffix>".
    export_gmt
        If provided, writes GMT gene sets:
          PC1_pos, PC1_neg, PC2_pos, ...
        If use_abs=True: PC1_abs, PC2_abs, ...
    min_abs_loading
        Optional filter: only keep genes with abs(loading) >= threshold.

    Returns
    -------
    dict mapping keys like "PC1_pos", "PC1_neg", "PC1_abs" -> DataFrame.
    Also stores results in adata.uns[store_key].
    """
    if loadings_key not in adata.varm:
        raise KeyError(f"adata.varm['{loadings_key}'] not found. Run bk.tl.pca first.")

    PCs = adata.varm[loadings_key]
    PCs = np.asarray(PCs, dtype=float)  # (n_vars, n_comps)
    n_vars, n_comps = PCs.shape

    if pcs is None:
        pcs = list(range(1, n_comps + 1))
    else:
        pcs = [int(p) for p in pcs]
        for p in pcs:
            if p < 1 or p > n_comps:
                raise ValueError(f"Requested PC={p}, but only 1..{n_comps} are available.")

    genes = adata.var_names.astype(str).to_numpy()

    out: dict[str, pd.DataFrame] = {}
    long_rows: list[pd.DataFrame] = []

    def _rank_one(pc_index0: int) -> dict[str, pd.DataFrame]:
        v = PCs[:, pc_index0].copy()  # length n_vars
        df = pd.DataFrame({gene_col: genes, "loading": v})

        if dropna:
            df = df.dropna(subset=["loading"])

        if min_abs_loading is not None:
            df = df.loc[df["loading"].abs() >= float(min_abs_loading)]

        if df.shape[0] == 0:
            return {}

        pc_name = f"PC{pc_index0 + 1}"

        if use_abs:
            df2 = df.assign(abs_loading=df["loading"].abs()).sort_values(
                "abs_loading", ascending=False
            ).head(int(n_top))
            df2 = df2.drop(columns=["abs_loading"])
            df2["pc"] = pc_name
            df2["sign"] = "abs"
            df2["rank"] = np.arange(1, df2.shape[0] + 1)
            return {f"{pc_name}_abs": df2[[ "pc", "sign", "rank", gene_col, "loading" ]]}

        # positive
        pos = df[df["loading"] > 0].sort_values("loading", ascending=False).head(int(n_top))
        pos = pos.copy()
        pos["pc"] = pc_name
        pos["sign"] = "pos"
        pos["rank"] = np.arange(1, pos.shape[0] + 1)

        res: dict[str, pd.DataFrame] = {f"{pc_name}_pos": pos[[ "pc", "sign", "rank", gene_col, "loading" ]]}

        if include_negative:
            neg = df[df["loading"] < 0].sort_values("loading", ascending=True).head(int(n_top))
            neg = neg.copy()
            neg["pc"] = pc_name
            neg["sign"] = "neg"
            neg["rank"] = np.arange(1, neg.shape[0] + 1)
            res[f"{pc_name}_neg"] = neg[[ "pc", "sign", "rank", gene_col, "loading" ]]

        return res

    for p in pcs:
        res = _rank_one(p - 1)
        for k, dfk in res.items():
            out[k] = dfk
            long_rows.append(dfk)

    if len(out) == 0:
        warn("pca_loadings: no genes returned (all NaN or filtered out).")
        adata.uns[store_key] = {"params": {"pcs": pcs, "n_top": n_top}, "results": {}}
        return out

    # ---- store in AnnData ----
    long_df = pd.concat(long_rows, axis=0, ignore_index=True)

    adata.uns[store_key] = {
        "params": {
            "key": key,
            "loadings_key": loadings_key,
            "pcs": list(pcs),
            "n_top": int(n_top),
            "use_abs": bool(use_abs),
            "include_negative": bool(include_negative),
            "min_abs_loading": min_abs_loading,
        },
        "results": {k: v for k, v in out.items()},
        "table": long_df,
    }

    # ---- export tabular ----
    if export is not None:
        export = Path(export)
        export.parent.mkdir(parents=True, exist_ok=True)
        long_df.to_csv(export, sep=export_sep, index=False)

        # wide format: one column per PC/sign, values are genes
        wide = {}
        for k, dfk in out.items():
            wide[k] = dfk[gene_col].tolist()
        wide_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in wide.items()]))

        wide_path = export.with_name(export.stem + ".wide" + export.suffix)
        wide_df.to_csv(wide_path, sep=export_sep, index=False)

        info(f"pca_loadings: exported {export} and {wide_path}")

    # ---- export GMT ----
    if export_gmt is not None:
        export_gmt = Path(export_gmt)
        export_gmt.parent.mkdir(parents=True, exist_ok=True)

        lines = []
        for k, dfk in out.items():
            # GMT format: set_name  description  gene1 gene2 ...
            set_name = k
            desc = f"{key}:{loadings_key}"
            geneset = dfk[gene_col].astype(str).tolist()
            # avoid empty sets
            if len(geneset) == 0:
                continue
            lines.append("\t".join([set_name, desc] + geneset))

        export_gmt.write_text("\n".join(lines) + ("\n" if lines else ""))
        info(f"pca_loadings: exported GMT to {export_gmt}")

    return out