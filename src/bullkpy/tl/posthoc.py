from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from scipy import stats

from ..logging import info, warn


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR correction."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out

    p0 = p[ok]
    n = p0.size
    order = np.argsort(p0)
    ranked = p0[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    out_ok = np.empty_like(p0)
    out_ok[order] = q
    out[ok] = out_ok
    return out


def _rank_biserial_from_u(U: float, n1: int, n2: int) -> float:
    # RBC in [-1, 1] (positive means group1 > group2 in ranks)
    denom = float(n1) * float(n2)
    if denom <= 0:
        return np.nan
    return 1.0 - (2.0 * float(U) / denom)


def pairwise_posthoc(
    df: pd.DataFrame,
    *,
    group_col: str = "grp",
    value_col: str = "y",
    method: Literal["mwu", "ttest"] = "mwu",
    correction: Literal["bh"] = "bh",
    dropna: bool = True,
) -> pd.DataFrame:
    """
    Pairwise post-hoc tests between groups.

    Parameters
    ----------
    df : DataFrame with columns [group_col, value_col]
    method : "mwu" (Mann-Whitney U, two-sided) or "ttest" (Welch t-test)
    correction : currently only "bh" (Benjamini-Hochberg)
    dropna : drop rows with NA in group/value

    Returns
    -------
    DataFrame with columns:
      group1, group2, n1, n2, pval, qval, effect, delta_mean, delta_median
    """
    if group_col not in df.columns or value_col not in df.columns:
        raise KeyError(f"df must contain '{group_col}' and '{value_col}'")

    d = df[[group_col, value_col]].copy()
    if dropna:
        d = d.dropna(subset=[group_col, value_col])

    d[group_col] = d[group_col].astype(str)
    groups = list(pd.Categorical(d[group_col]).categories)
    if len(groups) < 2:
        raise ValueError("Need at least 2 groups for posthoc.")

    rows = []
    pvals = []

    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            g1, g2 = groups[i], groups[j]
            x1 = d.loc[d[group_col] == g1, value_col].to_numpy(dtype=float)
            x2 = d.loc[d[group_col] == g2, value_col].to_numpy(dtype=float)
            n1, n2 = x1.size, x2.size
            if n1 < 2 or n2 < 2:
                p = np.nan
                eff = np.nan
            else:
                if method == "mwu":
                    # two-sided MWU
                    U, p = stats.mannwhitneyu(x1, x2, alternative="two-sided")
                    eff = _rank_biserial_from_u(U, n1, n2)
                elif method == "ttest":
                    t, p = stats.ttest_ind(x1, x2, equal_var=False, nan_policy="omit")
                    # Cohen's d (approx with pooled SD for reporting)
                    s1 = np.nanstd(x1, ddof=1)
                    s2 = np.nanstd(x2, ddof=1)
                    sp = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / max(n1 + n2 - 2, 1))
                    eff = (np.nanmean(x1) - np.nanmean(x2)) / (sp if sp > 0 else np.nan)
                else:
                    raise ValueError("method must be 'mwu' or 'ttest'")

            dm = np.nanmean(x1) - np.nanmean(x2)
            dmed = np.nanmedian(x1) - np.nanmedian(x2)

            rows.append(
                dict(
                    group1=g1,
                    group2=g2,
                    n1=int(n1),
                    n2=int(n2),
                    pval=float(p) if np.isfinite(p) else np.nan,
                    effect=float(eff) if np.isfinite(eff) else np.nan,
                    delta_mean=float(dm) if np.isfinite(dm) else np.nan,
                    delta_median=float(dmed) if np.isfinite(dmed) else np.nan,
                )
            )
            pvals.append(p)

    out = pd.DataFrame(rows)
    if correction == "bh":
        out["qval"] = _bh_fdr(out["pval"].to_numpy())
    else:
        out["qval"] = out["pval"]

    out = out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)
    return out