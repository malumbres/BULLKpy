from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from ..logging import info, warn

def pca(
    adata: ad.AnnData,
    *,
    layer: str | None = "log1p_cpm",
    n_comps: int = 50,
    center: bool = True,
    scale: bool = False,
    use_highly_variable: bool = False,
    key_added: str = "pca",
) -> None:
    X = adata.layers[layer] if layer is not None else adata.X

    # Track which genes are used so we can write loadings back to full var space
    used_mask = np.ones(adata.n_vars, dtype=bool)

    if use_highly_variable:
        if "highly_variable" not in adata.var.columns:
            warn("use_highly_variable=True but adata.var['highly_variable'] not found; using all genes.")
        else:
            used_mask = adata.var["highly_variable"].to_numpy(bool)
            if used_mask.sum() == 0:
                raise ValueError("adata.var['highly_variable'] has 0 True values.")
            X = X[:, used_mask]

    if sp.issparse(X):
        X = X.toarray()

    X = np.asarray(X, dtype=float)  # samples x genes_used
    n_obs, n_vars_used = X.shape
    n_comps = int(min(n_comps, n_obs - 1, n_vars_used))

    info(f"Running PCA on matrix {n_obs} samples × {n_vars_used} genes; n_comps={n_comps}")

    # Center/scale per gene
    mean_ = X.mean(axis=0) if center else np.zeros(n_vars_used, dtype=float)
    Xc = X - mean_

    std_ = None
    if scale:
        std_ = Xc.std(axis=0, ddof=1)
        std_[std_ == 0] = 1.0
        Xc = Xc / std_

    # SVD: Xc = U S Vt
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)

    # Scores (samples x comps)
    X_pca = U[:, :n_comps] * S[:n_comps]

    # Loadings for used genes (genes_used x comps)
    PCs_used = Vt[:n_comps, :].T

    # Explained variance
    eigvals = (S**2) / max(n_obs - 1, 1)
    var_ratio = eigvals / eigvals.sum()

    # Store scores
    adata.obsm["X_pca"] = X_pca

    # Store loadings in FULL gene space to satisfy AnnData alignment
    PCs_full = np.full((adata.n_vars, n_comps), np.nan, dtype=float)
    PCs_full[used_mask, :] = PCs_used
    adata.varm["PCs"] = PCs_full

    # Store metadata
    adata.uns[key_added] = {
        "params": {
            "layer": layer,
            "n_comps": n_comps,
            "center": center,
            "scale": scale,
            "use_highly_variable": use_highly_variable,
            "n_vars_used": int(n_vars_used),
        },
        "variance": eigvals[:n_comps],
        "variance_ratio": var_ratio[:n_comps],
        "used_genes_mask": used_mask,  # helpful for debugging / downstream
        "mean": mean_,
        "std": std_,
    }

    info(f"PCA stored in adata.obsm['X_pca'], adata.varm['PCs'], adata.uns['{key_added}']")


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