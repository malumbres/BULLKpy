from __future__ import annotations

from typing import Optional
import pandas as pd
from anndata import AnnData


def tcga_define_groups(
    adata: AnnData,
    *,
    project_key: str = "Project_ID",
    stage_key: str = "clinical_stage",
    subtype_key: str = "Subtype_Selected",
    out_prefix: str = "grp",
    make_lung_nsclc_subtypes: bool = True,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Create standardized grouping columns in `.obs` for TCGA-style analyses.

    Adds (examples):
    - f"{out_prefix}_tumor"           : tumor type (TCGA project)
    - f"{out_prefix}_stage_simple"   : I/II vs III/IV vs NA
    - f"{out_prefix}_subtype"        : subtype selected (if present)
    - f"{out_prefix}_lung_nsclc"     : LUAD vs LUSC (if present)

    Returns AnnData if copy=True else modifies in place.
    """
    ad = adata.copy() if copy else adata
    obs = ad.obs

    if project_key not in obs:
        raise KeyError(f"{project_key!r} not found in adata.obs")

    obs[f"{out_prefix}_tumor"] = pd.Categorical(obs[project_key].astype("string"))

    # Stage simplification
    if stage_key in obs:
        st = obs[stage_key].astype("string")
        st_simple = pd.Series(pd.NA, index=st.index, dtype="string")
        st_simple[st.str.contains("^stage i", case=False, na=False)] = "Stage I"
        st_simple[st.str.contains("^stage ii", case=False, na=False)] = "Stage II"
        st_simple[st.str.contains("^stage iii", case=False, na=False)] = "Stage III"
        st_simple[st.str.contains("^stage iv", case=False, na=False)] = "Stage IV"
        # combine
        st_2 = pd.Series(pd.NA, index=st.index, dtype="string")
        st_2[st_simple.isin(["Stage I", "Stage II"])] = "I/II"
        st_2[st_simple.isin(["Stage III", "Stage IV"])] = "III/IV"
        obs[f"{out_prefix}_stage_simple"] = pd.Categorical(st_2)

    # Selected subtype
    if subtype_key in obs:
        obs[f"{out_prefix}_subtype"] = pd.Categorical(obs[subtype_key].astype("string"))

    # Lung NSCLC: LUAD vs LUSC (TCGA reality)
    if make_lung_nsclc_subtypes:
        lung = obs[project_key].astype("string")
        lung_grp = pd.Series(pd.NA, index=lung.index, dtype="string")
        lung_grp[lung == "LUAD"] = "LUAD"
        lung_grp[lung == "LUSC"] = "LUSC"
        if lung_grp.notna().any():
            obs[f"{out_prefix}_lung_nsclc"] = pd.Categorical(lung_grp)

    ad.obs = obs
    return ad if copy else None