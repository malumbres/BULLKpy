from __future__ import annotations

import anndata as ad
import pandas as pd

__all__ = ["de"]

_STORES = ("de", "de_glm")


def _norm(key: str) -> str:
    return key.replace(":", "_")


def de(adata: ad.AnnData, contrast: str, *, key: str | None = None) -> pd.DataFrame:
    """
    Differential-expression results for one contrast, as a DataFrame.

    Looks in ``adata.uns["de"]`` (from ``bk.tl.de``, keys like ``"Project_ID_LUAD_vs_LUSC"``)
    and ``adata.uns["de_glm"]`` (from ``bk.tl.de_glm``, keys like ``"Subtype:LumB_vs_LumA"``).
    ``:`` and ``_`` are treated as the same separator, so either spelling finds either store.

    Parameters
    ----------
    contrast
        The contrast name.
    key
        Search only ``adata.uns[key]``. By default both ``"de"`` and ``"de_glm"`` are searched.

    Raises
    ------
    KeyError
        If no match is found; the message lists the contrasts that exist.
    """
    stores = (key,) if key is not None else _STORES
    available: list[str] = []
    for store in stores:
        entries = adata.uns.get(store)
        if not isinstance(entries, dict):
            continue
        for name, obj in entries.items():
            available.append(f"{store}: {name}")
            if name == contrast or _norm(str(name)) == _norm(contrast):
                if isinstance(obj, dict) and "results" in obj:
                    return pd.DataFrame(obj["results"])
                if isinstance(obj, pd.DataFrame):
                    return obj
                raise TypeError(f"Unsupported object at adata.uns['{store}']['{name}']: {type(obj).__name__}")
    shown = "; ".join(available) if available else f"no results in adata.uns{list(stores)}"
    raise KeyError(f"No DE results for contrast {contrast!r}. Available: {shown}")
