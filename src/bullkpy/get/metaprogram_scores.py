from __future__ import annotations

from typing import Sequence, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData



def _to_dense(X) -> np.ndarray:
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def metaprogram_scores(
    adata: AnnData,
    *,
    obsm_key_candidates: Sequence[str] = ("X_metaprogram", "X_metaprograms", "X_mp", "X_metaprog"),
    uns_df_key_candidates: Sequence[str] = ("metaprogram_scores", "mp_scores", "metaprograms"),
    names_key_candidates: Sequence[str] = ("metaprogram_names", "mp_names"),
    return_df: bool = True,
) -> Union[pd.DataFrame, Tuple[np.ndarray, list[str]]]:
    """
    Retrieve metaprogram score matrix (samples x programs) from AnnData.

    Supported storage:
      - adata.obsm[<key>] -> matrix; names in adata.uns['metaprogram_names'] (or similar)
      - adata.uns[<key>]  -> DataFrame indexed by sample, columns programs

    Parameters
    ----------
    return_df
        If True returns a DataFrame indexed by adata.obs_names.
        If False returns (X, names).

    Returns
    -------
    pd.DataFrame OR (np.ndarray, list[str])
    """
    # 1) Prefer obsm
    for k in obsm_key_candidates:
        if k in adata.obsm:
            X = _to_dense(adata.obsm[k]).astype(float)
            names = None
            for nk in names_key_candidates:
                if nk in adata.uns:
                    names = list(adata.uns[nk])
                    break
            if names is None:
                names = [f"MP{i+1}" for i in range(X.shape[1])]
            if return_df:
                return pd.DataFrame(X, index=adata.obs_names, columns=names)
            return X, names

    # 2) Fall back to uns dataframe
    for k in uns_df_key_candidates:
        if k in adata.uns:
            df = adata.uns[k]
            if not isinstance(df, pd.DataFrame):
                df = pd.DataFrame(df)
            # align to current obs order
            df = df.loc[adata.obs_names]
            if return_df:
                return df
            return df.to_numpy(dtype=float), list(df.columns)

    raise KeyError(
        "Metaprogram scores not found. Expected in adata.obsm (e.g. 'X_metaprogram') "
        "or adata.uns (e.g. 'metaprogram_scores')."
    )
