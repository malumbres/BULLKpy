from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

import anndata as ad

from ..logging import info, warn

try:
    import gseapy as gp
except Exception:  # pragma: no cover
    gp = None


def _ensure_gseapy() -> None:
    if gp is None:
        raise ImportError("gsea_preranked requires gseapy. Install with: pip install gseapy")


def list_enrichr_libraries() -> list[str]:
    """
    Return available Enrichr libraries as known by gseapy.
    Useful to discover exact names like 'MSigDB_Hallmark_2020'.

    Notes
    -----
    - Requires internet access (Enrichr).
    """
    _ensure_gseapy()
    try:
        libs = gp.get_library_name()
        return list(libs)
    except Exception as e:
        raise RuntimeError(f"Failed to fetch Enrichr library names: {e}") from e


def _normalize_gene_sets_arg(
    gene_sets: str | Sequence[str] | Mapping[str, Sequence[str]],
) -> str | list[str] | dict[str, list[str]]:
    """
    Accept:
      - dict of gene sets: {"SetA": [...], ...}
      - path to .gmt
      - enrichr library name: "MSigDB_Hallmark_2020"
      - convenience aliases: "hallmark", "hallmarks", "c2", "curated"
      - list of the above (mix allowed)
    """
    # dict passthrough
    if isinstance(gene_sets, Mapping):
        return {str(k): [str(x) for x in v] for k, v in gene_sets.items()}

    # turn comma-separated string into list
    if isinstance(gene_sets, str) and ("," in gene_sets):
        gene_sets_list = [x.strip() for x in gene_sets.split(",") if x.strip()]
        return _normalize_gene_sets_arg(gene_sets_list)

    def _alias_one(x: str) -> str:
        xl = x.strip().lower()
        # Convenience aliases → real Enrichr library names (common)
        # If your Enrichr instance has different years, use list_enrichr_libraries()
        if xl in {"hallmark", "hallmarks"}:
            return "MSigDB_Hallmark_2020"
        if xl in {"c2", "curated", "c2_curated"}:
            return "MSigDB_C2_Curated_2021"
        return x

    # list/tuple
    if isinstance(gene_sets, (list, tuple)):
        out: list[str] = []
        for x in gene_sets:
            if not isinstance(x, str):
                raise TypeError("gene_sets list must contain strings (library names or .gmt paths).")
            out.append(_alias_one(x))
        return out

    # plain string
    if isinstance(gene_sets, str):
        return _alias_one(gene_sets)

    raise TypeError("gene_sets must be a str, list[str], or dict[str, list[str]].")


def _gseapy_result_to_table(pre_res: Any) -> pd.DataFrame:
    """
    Convert gseapy Prerank result object to a tidy pandas DataFrame.
    """
    # gseapy stores results in pre_res.res2d (DataFrame)
    if not hasattr(pre_res, "res2d"):
        raise ValueError("Unexpected gseapy result object: missing attribute 'res2d'.")
    df = pre_res.res2d.copy()

    # make column names consistent and flat
    df.columns = [str(c) for c in df.columns]
    df.index = df.index.astype(str)
    df = df.reset_index().rename(columns={"index": "term"})

    # common columns you may want:
    # term, es, nes, pval, fdr, ... (depends on gseapy version)
    return df


def gsea_preranked(
    adata: ad.AnnData | None = None,
    *,
    res: pd.DataFrame,
    gene_sets: str | Sequence[str] | Mapping[str, Sequence[str]] = "MSigDB_Hallmark_2020",
    comparison: str | None = None,
    store_key: str = "gsea",
    results_key: str = "results",
    gene_col: str = "gene",
    score_col: str = "t",
    ascending: bool = False,
    min_size: int = 5,
    max_size: int = 1000,
    permutation_num: int = 1000,
    seed: int = 6,
    processes: int = 4,
    outdir: str | Path = "gsea",
    save_rnk: bool = True,
    return_pre_res: bool = True,
    add_comparison_column: bool = True,
) -> pd.DataFrame | tuple[pd.DataFrame, Any]:
    """
    Run GSEApy prerank and return a tidy results table.

    Key improvements
    ----------------
    1) Writes to outdir/<comparison>/ (folder created with parents=True)
    2) Adds a 'comparison' column to the returned results
    3) gene_sets can be:
         - Enrichr library name, e.g. "MSigDB_Hallmark_2020"
         - .gmt path
         - dict of gene sets
         - convenience aliases: "hallmark(s)", "c2"/"curated"
         - list of libraries/paths (mix allowed)
    4) Optionally returns the raw gseapy 'pre_res' object so you can plot later.

    Notes
    -----
    - The returned 'pre_res' is NOT stored in adata.uns by default, because it is
      not h5ad-serializable. Store tables only.
    """
    _ensure_gseapy()

    if res is None or not isinstance(res, pd.DataFrame):
        raise TypeError("res must be a pandas DataFrame (your DE table).")
    if gene_col not in res.columns:
        raise KeyError(f"gene_col='{gene_col}' not found in res.columns.")
    if score_col not in res.columns:
        raise KeyError(f"score_col='{score_col}' not found in res.columns.")

    # infer comparison if missing (best-effort)
    if comparison is None:
        comparison = "comparison"

    # build output folder: outdir/comparison/
    outdir = Path(outdir) / str(comparison)
    outdir.mkdir(parents=True, exist_ok=True)

    # build rank table
    rnk = res[[gene_col, score_col]].copy()
    rnk = rnk.dropna(subset=[gene_col, score_col])
    rnk[gene_col] = rnk[gene_col].astype(str).str.upper()
    rnk[score_col] = pd.to_numeric(rnk[score_col], errors="coerce")
    rnk = rnk.dropna(subset=[score_col])

    # handle duplicates: keep best score by absolute value (common)
    if rnk[gene_col].duplicated().any():
        frac_dup = float(rnk[gene_col].duplicated().mean()) * 100.0
        warn(f"Duplicated genes found in preranked stats: {frac_dup:.2f}% of genes. Keeping max |score| per gene.")
        rnk = (
            rnk.assign(_abs=np.abs(rnk[score_col].to_numpy()))
            .sort_values("_abs", ascending=False)
            .drop_duplicates(subset=[gene_col], keep="first")
            .drop(columns="_abs")
        )

    # sort for prerank
    rnk = rnk.sort_values(score_col, ascending=bool(ascending))

    # optionally save .rnk
    if save_rnk:
        rnk_path = outdir / f"{comparison}.rnk.tsv"
        rnk.to_csv(rnk_path, sep="\t", index=False)

    # resolve gene_sets aliases / list handling
    gs = _normalize_gene_sets_arg(gene_sets)

    info(f"Running GSEApy prerank: comparison={comparison} n_genes={len(rnk)} outdir={outdir}")

    # Run prerank (no_plot=True keeps it light; you can plot later from pre_res)
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gs,
        min_size=int(min_size),
        max_size=int(max_size),
        permutation_num=int(permutation_num),
        seed=int(seed),
        processes=int(processes),
        outdir=str(outdir),
        format="png",
        no_plot=True,
    )

    df = _gseapy_result_to_table(pre_res)

    if add_comparison_column:
        df.insert(0, "comparison", str(comparison))

    # store table in adata.uns (serializable)
    if adata is not None:
        adata.uns.setdefault(store_key, {})
        adata.uns[store_key].setdefault(str(comparison), {})
        adata.uns[store_key][str(comparison)][results_key] = df

        # helpful provenance
        adata.uns[store_key][str(comparison)]["params"] = {
            "gene_sets": gs,
            "gene_col": gene_col,
            "score_col": score_col,
            "ascending": bool(ascending),
            "min_size": int(min_size),
            "max_size": int(max_size),
            "permutation_num": int(permutation_num),
            "seed": int(seed),
            "processes": int(processes),
            "outdir": str(outdir),
            "save_rnk": bool(save_rnk),
        }

    return (df, pre_res) if return_pre_res else df