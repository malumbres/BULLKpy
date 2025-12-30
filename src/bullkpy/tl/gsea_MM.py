from __future__ import annotations

import os
import gseapy

from ..logging import info, warn


def _organize_folders(dataname, comparison):
    """Organize variables & folders"""
    outdir = os.path.join(DESKTOP, "de", dataname, comparison)
    os.makedirs(outdir, exist_ok=True)
    return outdir

def _res_to_rnk_mouse_to_human(res, outdir, dataname, comparison):
    """Generates a .rnk file from differential expression for GSEA analysis."""

    res["scores"] = res["t"] * np.abs(res["log2FC"])
    rnk = res[["gene", "scores"]].copy()
    rnk.columns = ["#names", "scores"]
    rnk["#names"] = rnk["#names"].astype(str).str.upper()                     # all to uppercase
    rnk.to_csv(outdir + "/" + comparison + ".rnk", sep="\t", index=False)
    
    return rnk

def _pre_res():
    """Empty table for final results"""
    pre_res = pd.DataFrame(columns=["Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes", "Sample"])
    return pre_res

def _rnk_to_geseapy(pre_res, rnk, gset, outdir, comparison):
    """Run GSEApy and merge index results with a previous table."""
    
    gene_set    = gset + ".gmt"
    pre_res_new = gspy.prerank(rnk=rnk, gene_sets= gene_set, processes=4, permutation_num=1000,
                               outdir=outdir, graph_num=60, format='png', seed=6)
    pre_res_new_df = pd.DataFrame(pre_res_new.res2d.sort_index())
    pre_res_new_df = pre_res_new_df.drop(columns="Name")
    pre_res_new_df["Sample"] = comparison
    pre_res_new_df.to_csv(outdir + "/" + "GSEApy_results.tsv", sep="\t")
    pre_res     = pd.concat([pre_res, pre_res_new_df], join="inner", ignore_index=True)

    return pre_res


def pre_res_from_gseapy(dataname_comparison):
    outdir = _organize_folders(dataname, comparison)
    rnk = res_to_rnk_mouse_to_human(res, outdir, dataname, comparison)
    pre_res = _pre_res()
    pre_res = pre_res = _rnk_to_geseapy(pre_res, rnk, gset, outdir, comparison)
    return pre_res

