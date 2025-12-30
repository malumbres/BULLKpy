from __future__ import annotations

from typing import Literal, Sequence

import pandas as pd

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None

from matplotlib.colors import to_hex

from .._settings import settings

Where = Literal["obs", "var"]


def _get_series(adata, key: str, where: Where) -> pd.Series:
    if where == "obs":
        if key not in adata.obs.columns:
            raise KeyError(f"adata.obs['{key}'] not found")
        s = adata.obs[key]
        s = s.copy()
        s.name = key
        return s
    if where == "var":
        if key not in adata.var.columns:
            raise KeyError(f"adata.var['{key}'] not found")
        s = adata.var[key]
        s = s.copy()
        s.name = key
        return s
    raise ValueError("where must be 'obs' or 'var'")


def get_categorical_colors(
    adata,
    *,
    key: str,
    where: Where = "obs",
    categories: Sequence[str] | None = None,
    palette: str | None = None,
    force: bool = False,
) -> dict[str, str]:
    """
    Stable mapping category -> hex color for a categorical column.
    Cached in `bk.settings.categorical_colors[(where, key)]`.

    If you set a per-column palette override via:
      settings.categorical_palette_overrides[(where, key)] = "Set2"
    it will be used automatically.
    """
    if sns is None:
        raise ImportError("Categorical color mapping requires seaborn.")

    cache_key = (where, key)

    # return cached unless forcing
    if (not force) and cache_key in settings.categorical_colors:
        cmap = settings.categorical_colors[cache_key]
        if categories is None:
            return cmap
        return {str(c): cmap.get(str(c), "#B0B0B0") for c in categories}

    s = _get_series(adata, key, where).astype(str)

    if categories is None:
        cats = list(pd.Categorical(s).categories)
    else:
        cats = [str(c) for c in categories]

    pal_name = palette
    if pal_name is None:
        pal_name = settings.categorical_palette_overrides.get(cache_key, settings.plot_palette)

    pal = sns.color_palette(pal_name, n_colors=len(cats))
    cmap = {cats[i]: to_hex(pal[i]) for i in range(len(cats))}

    settings.categorical_colors[cache_key] = cmap
    return cmap


def categorical_colors_array(
    adata,
    *,
    key: str,
    where: Where = "obs",
    categories: Sequence[str] | None = None,
    palette: str | None = None,
    force: bool = False,
) -> pd.Series:
    """Series of hex colors aligned to obs_names / var_names."""
    cmap = get_categorical_colors(
        adata, key=key, where=where, categories=categories, palette=palette, force=force
    )
    s = _get_series(adata, key, where).astype(str)
    return s.map(lambda x: cmap.get(str(x), "#B0B0B0"))