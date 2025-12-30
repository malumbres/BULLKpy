from __future__ import annotations

import pandas as pd
import anndata as ad


def get_de_table(
    adata: ad.AnnData,
    *,
    uns_key: str,
    contrast: str,
) -> pd.DataFrame:
    obj = adata.uns[uns_key][contrast]
    if isinstance(obj, dict) and "results" in obj:
        return obj["results"]
    if isinstance(obj, pd.DataFrame):
        return obj
    raise TypeError(f"Unsupported DE object at adata.uns['{uns_key}']['{contrast}']")