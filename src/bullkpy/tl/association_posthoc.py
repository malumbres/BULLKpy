from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind

from ..logging import info


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    m = np.isfinite(p)
    if m.sum() == 0:
        return out
    pv = p[m]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out[m] = q[inv]
    return out


def pairwise_posthoc(
    df: pd.DataFrame,
    *,
    y_col: str = "y",
    group_col: str = "grp",
    method: Literal["mwu", "ttest"] = "mwu",
    adjust: Literal["fdr_bh", "none"] = "fdr_bh",
    min_group_size: int = 2,
) -> pd.DataFrame:
    """
    Pairwise post-hoc comparisons between all group pairs.

    Input df must have:
      - y_col: numeric values
      - group_col: categorical group labels

    Returns tidy table with:
      group1, group2, stat, pval, qval, mean1, mean2, diff
    """
    d = df[[y_col, group_col]].copy()
    d[group_col] = d[group_col].astype(str)
    d[y_col] = pd.to_numeric(d[y_col], errors="coerce")
    d = d.dropna()

    cats = list(pd.Categorical(d[group_col]).categories)
    info(f"pairwise_posthoc: {len(cats)} groups, method={method}")

    rows = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            g1, g2 = cats[i], cats[j]
            x = d.loc[d[group_col] == g1, y_col].to_numpy(dtype=float)
            y = d.loc[d[group_col] == g2, y_col].to_numpy(dtype=float)
            if x.size < min_group_size or y.size < min_group_size:
                continue

            if method == "mwu":
                stat, pval = mannwhitneyu(x, y, alternative="two-sided")
            elif method == "ttest":
                stat, pval = ttest_ind(x, y, equal_var=False, nan_policy="omit")
            else:
                raise ValueError(f"Unknown method='{method}'")

            rows.append(
                {
                    "group1": g1,
                    "group2": g2,
                    "statistic": float(stat),
                    "pval": float(pval),
                    "mean1": float(np.mean(x)),
                    "mean2": float(np.mean(y)),
                    "diff": float(np.mean(x) - np.mean(y)),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    if adjust == "fdr_bh":
        out["qval"] = _bh_fdr(out["pval"].to_numpy())
    else:
        out["qval"] = np.nan

    return out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)