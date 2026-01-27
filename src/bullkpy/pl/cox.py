from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal, Optional, Union, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData

from ._style import set_style, _savefig

from ..logging import warn



def cox_volcano(
    df_cox: pd.DataFrame,
    *,
    x: Literal["logHR", "coef"] = "logHR",
    sig_col: Literal["qval", "pval"] = "qval",
    alpha: float = 0.05,
    min_effect_loghr: float | None = None,   # e.g. np.log(1.5)
    top_labels: int = 12,
    label_by: Literal["qval", "pval", "abs_logHR"] = "qval",
    figsize: tuple[float, float] = (6.6, 5.6),
    title: str | None = None,
    highlight_color: str = "crimson",
    base_color: str = "0.35",
    point_size: float = 12.0,
    point_alpha: float = 0.75,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Cox volcano plot: log(HR) vs -log10(significance).

    Expects df_cox from bk.tl.cox_gene_association with at least:
      - gene
      - HR and/or coef
      - pval
      - qval (if BH correction used)

    Parameters
    ----------
    x
        "logHR" (default) uses log(HR) on x-axis, or "coef" uses Cox coefficient.
    sig_col
        Column used for y-axis: "qval" (default) or "pval".
    alpha
        Significance threshold (horizontal line at -log10(alpha)).
    min_effect_loghr
        Optional effect-size threshold in log(HR); if provided, adds vertical lines at +/- value.
        Example: np.log(1.5) to highlight HR >= 1.5 or <= 1/1.5.
    top_labels
        Number of genes to label (selected by label_by).
    label_by
        How to choose which genes to label.
    """
    set_style()

    if not isinstance(df_cox, pd.DataFrame):
        raise TypeError("df_cox must be a pandas DataFrame.")

    req = {"gene", "pval"}
    missing = req - set(df_cox.columns)
    if missing:
        raise ValueError(f"df_cox missing required columns: {sorted(missing)}")

    df = df_cox.copy()
    df["gene"] = df["gene"].astype(str)

    # ---- x-axis values ----
    if x == "coef":
        if "coef" not in df.columns:
            # try deriving from HR
            if "HR" in df.columns:
                df["coef"] = np.log(pd.to_numeric(df["HR"], errors="coerce"))
            else:
                raise ValueError("Need 'coef' or 'HR' in df_cox for x='coef'.")
        xvals = pd.to_numeric(df["coef"], errors="coerce").to_numpy(dtype=float)
        xlabel = "Cox coefficient"
    else:
        # logHR
        if "HR" not in df.columns:
            if "coef" in df.columns:
                df["HR"] = np.exp(pd.to_numeric(df["coef"], errors="coerce"))
            else:
                raise ValueError("Need 'HR' or 'coef' in df_cox for x='logHR'.")
        hr = pd.to_numeric(df["HR"], errors="coerce").to_numpy(dtype=float)
        xvals = np.log(hr)
        xlabel = "log(HR)"

    # ---- y-axis values ----
    if sig_col not in df.columns:
        raise ValueError(f"sig_col='{sig_col}' not in df_cox. Available: {list(df.columns)}")
    sig = pd.to_numeric(df[sig_col], errors="coerce").to_numpy(dtype=float)
    sig = np.clip(sig, 1e-300, 1.0)
    yvals = -np.log10(sig)

    # ---- mask finite ----
    ok = np.isfinite(xvals) & np.isfinite(yvals)
    df = df.loc[ok].copy()
    xvals = xvals[ok]
    yvals = yvals[ok]
    sig = sig[ok]

    if df.shape[0] == 0:
        raise ValueError("No finite points to plot.")

    # ---- significance mask ----
    is_sig = sig <= float(alpha)
    if min_effect_loghr is not None:
        is_sig = is_sig & (np.abs(xvals) >= float(min_effect_loghr))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # base points
    ax.scatter(
        xvals,
        yvals,
        s=float(point_size),
        alpha=float(point_alpha),
        c=base_color,
        edgecolors="none",
    )
    # highlight
    if np.any(is_sig):
        ax.scatter(
            xvals[is_sig],
            yvals[is_sig],
            s=float(point_size) * 1.1,
            alpha=min(1.0, float(point_alpha) + 0.15),
            c=highlight_color,
            edgecolors="none",
        )

    # threshold lines
    ax.axhline(-np.log10(alpha), lw=1.0, color="0.25")
    if min_effect_loghr is not None:
        ax.axvline(+float(min_effect_loghr), lw=1.0, color="0.25", ls="--")
        ax.axvline(-float(min_effect_loghr), lw=1.0, color="0.25", ls="--")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"-log10({sig_col})")
    ax.set_title(title or f"Cox volcano ({sig_col}, alpha={alpha})")

    # ---- label selection ----
    if int(top_labels) > 0:
        df_lab = df.copy()
        df_lab["x"] = xvals
        df_lab["y"] = yvals
        df_lab["sig"] = sig

        if label_by == "abs_logHR":
            df_lab["ranker"] = np.abs(df_lab["x"])
            df_lab = df_lab.sort_values("ranker", ascending=False)
        elif label_by == "pval":
            if "pval" not in df_lab.columns:
                df_lab["pval"] = df_lab["sig"]
            df_lab = df_lab.sort_values("pval", ascending=True)
        else:  # qval
            if "qval" in df_lab.columns:
                df_lab = df_lab.sort_values("qval", ascending=True)
            else:
                df_lab = df_lab.sort_values("pval", ascending=True)

        df_lab = df_lab.head(int(top_labels))

        for _, r in df_lab.iterrows():
            ax.text(float(r["x"]), float(r["y"]), str(r["gene"]), fontsize=9)

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax, df

def cox_forest(
    df_cox: pd.DataFrame,
    *,
    genes: Sequence[str] | None = None,
    sort_by: Literal["qval", "pval", "abs_logHR"] = "qval",
    top_n: int | None = 20,
    log_scale: bool = False,  # <- default to decimal scale
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    show_ci: bool = True,
    color_risk: bool = True,
    risk_color: str = "crimson",
    protect_color: str = "royalblue",
    neutral_color: str = "0.4",
    save: str | Path | None = None,
    show: bool = True,
    hr_text_x: float = 1.02,          # <- x position in axes-fraction (to the right)
    right_margin: float = 0.78,       # <- makes room for the text column
):
    """
    Forest plot for Cox proportional hazards results.

    Expects df_cox with columns:
      - gene
      - HR
      - CI_lower
      - CI_upper
      - pval and/or qval
    """
    set_style()

    required = {"gene", "HR", "CI_lower", "CI_upper"}
    missing = required - set(df_cox.columns)
    if missing:
        raise ValueError(f"df_cox missing required columns: {sorted(missing)}")

    df = df_cox.copy()
    df["gene"] = df["gene"].astype(str)

    # ---- select genes ----
    if genes is not None:
        genes = [str(g) for g in genes]
        df = df[df["gene"].isin(genes)].copy()
    else:
        if sort_by == "abs_logHR":
            df["_rank"] = np.abs(np.log(df["HR"]))
            df = df.sort_values("_rank", ascending=False)
        else:
            if sort_by not in df.columns:
                raise ValueError(f"sort_by='{sort_by}' not in df_cox.")
            df = df.sort_values(sort_by, ascending=True)

        if top_n is not None:
            df = df.head(int(top_n))

    if df.shape[0] == 0:
        raise ValueError("No genes selected for forest plot.")

    # ---- ordering for plotting (top at top) ----
    df = df.iloc[::-1].reset_index(drop=True)

    y = np.arange(df.shape[0])

    HR = df["HR"].to_numpy(dtype=float)
    lo = df["CI_lower"].to_numpy(dtype=float)
    hi = df["CI_upper"].to_numpy(dtype=float)

    # ---- colors ----
    if color_risk:
        colors = [
            risk_color if hr > 1 else protect_color if hr < 1 else neutral_color
            for hr in HR
        ]
    else:
        colors = [neutral_color] * len(HR)

    # ---- figure size ----
    if figsize is None:
        h = max(4.0, 0.35 * len(HR) + 1.8)
        figsize = (7.6, h)  # a touch wider to accommodate the text column

    # IMPORTANT: use subplots_adjust to reserve space for right-side text
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=False)
    fig.subplots_adjust(right=right_margin)

    # CI bars
    if show_ci:
        ax.hlines(y, lo, hi, color="0.3", lw=1.5)

    # HR points
    ax.scatter(HR, y, c=colors, s=55, zorder=3)

    # reference line
    ax.axvline(1.0, color="0.25", lw=1.2, ls="--")

    # y-axis
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"].tolist())
    ax.set_ylim(-0.5, len(y) - 0.5)

    # x-axis
    if log_scale:
        ax.set_xscale("log")
        ax.set_xlabel("Hazard ratio (log scale)")
    else:
        ax.set_xlabel("Hazard ratio")
        # Optional: sane x-limits with padding on linear scale
        xmin = np.nanmin(lo)
        xmax = np.nanmax(hi)
        if np.isfinite(xmin) and np.isfinite(xmax) and xmax > xmin:
            pad = 0.08 * (xmax - xmin)
            ax.set_xlim(max(0, xmin - pad), xmax + pad)

    # title
    ax.set_title(title or "Cox proportional hazards (forest plot)")

    # ---- annotate HR + CI OUTSIDE the plotting area ----
    # x in axes fraction, y in data coords
    from matplotlib.transforms import blended_transform_factory
    trans = blended_transform_factory(ax.transAxes, ax.transData)

    for i, (hr, l, u) in enumerate(zip(HR, lo, hi)):
        ax.text(
            hr_text_x,
            y[i],
            f"{hr:.2f} [{l:.2f}, {u:.2f}]",
            transform=trans,
            va="center",
            ha="left",
            fontsize=9,
            color="0.25",
            clip_on=False,   # <- allows drawing outside axes
        )

    if save is not None:
        _savefig(fig, save)  # ideally your _savefig uses bbox_inches="tight"
    if show:
        plt.show()

    return fig, ax, df


def km_plot_signature(
    adata: ad.AnnData,
    *,
    time_col: str,
    event_col: str,
    score_col: str,
    split: str = "median",
    title: str | None = None,
    figsize: tuple[float, float] = (6.3, 5.0),
    show_ci: bool = True,
    ci_alpha: float = 0.18,
    show_risk_table: bool = True,
    risk_times: Sequence[float] | None = None,  # if None -> auto from x ticks
    # NEW: colors
    high_color: str = "#d62728",   # red
    low_color: str = "#1f77b4",    # blue
    palette: dict[str, str] | None = None,  # optional override: {"High": "...", "Low": "..."}
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Kaplan–Meier plot split by signature score (default: median).

    Adds:
      - optional 95% CI (log-log Greenwood)
      - optional numbers-at-risk table under the plot

    Coloring:
      - By default: High=red, Low=blue
      - Override with palette={"High": "...", "Low": "..."}
    """
    set_style()

    for c in (time_col, event_col, score_col):
        if c not in adata.obs.columns:
            raise KeyError(f"'{c}' not in adata.obs")

    t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(float)
    e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(float).astype(int)
    s = pd.to_numeric(adata.obs[score_col], errors="coerce").to_numpy(float)

    ok = np.isfinite(t) & np.isfinite(s) & np.isfinite(e)
    t, e, s = t[ok], e[ok], s[ok]

    if t.size < 5:
        raise ValueError("Not enough valid samples for KM plot.")

    thr = np.nanmedian(s) if split == "median" else np.nanmedian(s)
    hi = s >= thr
    lo = ~hi

    # choose colors (default or palette override)
    if palette is None:
        palette = {"High": high_color, "Low": low_color}
    c_hi = palette.get("High", high_color)
    c_lo = palette.get("Low", low_color)

    # layout: add extra bottom space for risk table
    if show_risk_table:
        fig, ax = plt.subplots(figsize=figsize)
        fig.subplots_adjust(bottom=0.25)
    else:
        fig, ax = plt.subplots(figsize=figsize)

    # --- KM curves ---
    def _plot_group(mask, label, color):
        xs, ys, vs, _ = _km_curve(t[mask], e[mask])
        ax.step(xs, ys, where="post", lw=2.0, label=label, color=color)

        if show_ci:
            lo_ci, hi_ci = _km_ci_loglog(ys, vs, z=1.96)
            ax.fill_between(
                xs, lo_ci, hi_ci,
                step="post",
                alpha=float(ci_alpha),
                color=color,
                linewidth=0,
            )
        return xs, ys

    _plot_group(lo, label=f"Low (n={int(lo.sum())})", color=c_lo)
    _plot_group(hi, label=f"High (n={int(hi.sum())})", color=c_hi)

    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.set_title(title or f"KM by {score_col} ({split} split)")
    ax.legend(frameon=False)

    # log-rank
    p = _logrank_pvalue(t[lo], e[lo], t[hi], e[hi])
    ax.text(
        0.99, 0.02, f"log-rank p={p:.2g}",
        transform=ax.transAxes, ha="right", va="bottom", color="0.25"
    )

    # --- numbers at risk table ---
    if show_risk_table:
        if risk_times is None:
            fig.canvas.draw_idle()
            xt = ax.get_xticks()
            xt = xt[np.isfinite(xt)]
            xt = xt[xt >= 0]
            if len(xt) >= 2:
                risk_times = xt
            else:
                tmax = float(np.nanmax(t))
                risk_times = np.linspace(0, tmax, 5)

        risk_times = np.asarray(list(risk_times), dtype=float)
        risk_lo = _n_at_risk(t[lo], risk_times)
        risk_hi = _n_at_risk(t[hi], risk_times)

        cell_text = [
            [str(x) for x in risk_lo],
            [str(x) for x in risk_hi],
        ]
        row_labels = ["Low", "High"]
        col_labels = [f"{x:g}" for x in risk_times]

        tbl = ax.table(
            cellText=cell_text,
            rowLabels=row_labels,
            colLabels=col_labels,
            loc="bottom",
            cellLoc="center",
            rowLoc="center",
            bbox=[0.0, -0.33, 1.0, 0.23],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(max(8, plt.rcParams.get("font.size", 12) - 3))

        # Optional: color row labels to match curves (nice touch)
        # Row-label cells are in column -1 for matplotlib tables.
        try:
            tbl[(1, -1)].get_text().set_color(c_lo)  # "Low"
            tbl[(2, -1)].get_text().set_color(c_hi)  # "High"
        except Exception:
            pass

        ax.text(
            0.0, -0.13, "Number at risk",
            transform=ax.transAxes, ha="left", va="center", color="0.35"
        )

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax, {"threshold": float(thr), "logrank_p": float(p)}


def _km_curve(t, e):
    """
    Compute KM step curve + Greenwood variance components.

    Returns
    -------
    times_step : np.ndarray
        Step x positions (including 0).
    surv_step : np.ndarray
        Step survival values.
    greenwood_var_step : np.ndarray
        Greenwood variance of S(t) at step points.
    uniq_event_times : np.ndarray
        Unique times considered (all unique times, including censor-only times).
    """
    order = np.argsort(t)
    t = t[order]
    e = e[order].astype(int)

    uniq = np.unique(t)
    at_risk = len(t)
    S = 1.0
    gv = 0.0  # Greenwood cumulative sum: Σ d/(n(n-d))

    xs = [0.0]
    ys = [1.0]
    vs = [0.0]

    for ti in uniq:
        di = int(((t == ti) & (e == 1)).sum())
        ci = int(((t == ti) & (e == 0)).sum())

        if at_risk > 0:
            if di > 0:
                # KM update
                S *= (1.0 - di / at_risk)
                # Greenwood increment
                if at_risk - di > 0:
                    gv += di / (at_risk * (at_risk - di))

            xs.extend([ti, ti])
            ys.extend([ys[-1], S])
            vs.extend([vs[-1], (S * S) * gv])  # Var(S) ≈ S^2 * Σ d/(n(n-d))

        at_risk -= (di + ci)

    return np.asarray(xs), np.asarray(ys), np.asarray(vs), np.asarray(uniq)


def _km_ci_loglog(surv, var_s, z=1.96):
    """
    Log-log transformed CI for KM curve (more stable than linear CI).
    surv: survival values
    var_s: Greenwood variance of surv
    """
    surv = np.asarray(surv, dtype=float)
    var_s = np.asarray(var_s, dtype=float)

    lo = np.full_like(surv, np.nan)
    hi = np.full_like(surv, np.nan)

    # only defined for 0 < S < 1
    m = (surv > 0) & (surv < 1) & np.isfinite(var_s)

    se_s = np.zeros_like(surv)
    se_s[m] = np.sqrt(np.maximum(var_s[m], 0.0))

    # log(-log(S)) transform
    ll = np.zeros_like(surv)
    ll[m] = np.log(-np.log(surv[m]))

    # delta method: var(log(-log(S))) ≈ var(S) / (S^2 * (log(S))^2)
    denom = (surv[m] * np.abs(np.log(surv[m])))
    denom[denom == 0] = np.nan
    se_ll = se_s[m] / denom

    lo_ll = ll[m] - z * se_ll
    hi_ll = ll[m] + z * se_ll

    lo[m] = np.exp(-np.exp(hi_ll))  # note swap for proper ordering
    hi[m] = np.exp(-np.exp(lo_ll))

    # endpoints
    lo[surv == 1.0] = 1.0
    hi[surv == 1.0] = 1.0
    lo[surv == 0.0] = 0.0
    hi[surv == 0.0] = 0.0

    return lo, hi


def _n_at_risk(t, times):
    t = np.asarray(t, dtype=float)
    times = np.asarray(times, dtype=float)
    return np.array([int(np.sum(t >= ti)) for ti in times], dtype=int)


def _logrank_pvalue(t1, e1, t2, e2):
    # Simple log-rank using normal approximation; good enough for QC/EDA.
    times = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    O1 = E1 = V1 = 0.0

    for ti in times:
        r1 = np.sum(t1 >= ti)
        r2 = np.sum(t2 >= ti)
        d1 = np.sum((t1 == ti) & (e1 == 1))
        d2 = np.sum((t2 == ti) & (e2 == 1))
        r = r1 + r2
        d = d1 + d2
        if r <= 1 or d == 0:
            continue
        exp1 = d * (r1 / r)
        var1 = (r1 * r2 * d * (r - d)) / (r**2 * (r - 1))
        O1 += d1
        E1 += exp1
        V1 += var1

    if V1 <= 0:
        return 1.0

    z = (O1 - E1) / np.sqrt(V1)

    # two-sided p using Normal distribution
    try:
        from scipy.stats import norm  # scipy is already in bullkpy deps
        p = 2.0 * norm.sf(np.abs(z))
    except Exception:
        # fallback without scipy
        import math
        # sf(z) = 0.5 * erfc(z/sqrt(2))
        p = 2.0 * 0.5 * math.erfc(abs(float(z)) / math.sqrt(2.0))

    return float(p)

def time_dependent_roc(
    adata,
    *,
    time_col: str,
    event_col: str,
    score_col: str,
    times: Sequence[float],
    title: str | None = None,
    figsize: tuple[float, float] = (6.3, 5.0),
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Time-dependent ROC/AUC using scikit-survival (optional dependency).

    `score_col`: higher = higher risk (Cox-style linear predictor).
    """
    set_style()

    try:
        from sksurv.metrics import cumulative_dynamic_auc
        from sksurv.util import Surv
    except Exception as e:
        raise ImportError(
            "time_dependent_roc requires scikit-survival.\n"
            "Install with: pip install scikit-survival\n"
            f"Original error: {e}"
        )

    for c in (time_col, event_col, score_col):
        if c not in adata.obs.columns:
            raise KeyError(f"'{c}' not in adata.obs")

    t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(float)
    e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(float).astype(bool)
    s = pd.to_numeric(adata.obs[score_col], errors="coerce").to_numpy(float)

    ok = np.isfinite(t) & np.isfinite(s)
    t, e, s = t[ok], e[ok], s[ok]

    y = Surv.from_arrays(event=e, time=t)
    times = np.asarray(list(times), dtype=float)

    # Using same cohort for train/test AUC is optimistic;
    # for proper estimation, pass external/test data or use CV.
    auc, mean_auc = cumulative_dynamic_auc(y, y, s, times)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(times, auc, marker="o")
    ax.set_xlabel("Time")
    ax.set_ylabel("AUC(t)")
    ax.set_title(title or f"Time-dependent AUC (score={score_col})")
    ax.set_ylim(0.0, 1.0)
    ax.text(0.99, 0.02, f"mean AUC={mean_auc:.2f}", transform=ax.transAxes, ha="right", va="bottom", color="0.25")

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax, {"times": times, "auc": auc, "mean_auc": float(mean_auc)}


def panel_size_cindex_plot(
    df_curve: pd.DataFrame,
    *,
    title: str | None = "Panel size vs C-index",
    xlabel: str = "Panel size (k)",
    ylabel: str = "CV C-index",
    show_errorbar: bool = True,
    figsize: tuple[float, float] = (6.5, 4.5),
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Plot panel size vs (mean±sd) C-index curve produced by `bk.tl.panel_size_cindex`.

    Parameters
    ----------
    df_curve
        DataFrame with columns: size, c_index_mean, c_index_sd (optional).
    """
    set_style()
    if not isinstance(df_curve, pd.DataFrame):
        raise TypeError("df_curve must be a pandas DataFrame.")
    if "size" not in df_curve.columns or "c_index_mean" not in df_curve.columns:
        raise ValueError("df_curve must contain columns: ['size', 'c_index_mean'] (and optionally 'c_index_sd').")

    df = df_curve.sort_values("size").copy()

    fig, ax = plt.subplots(figsize=figsize)
    x = df["size"].to_numpy(dtype=float)
    y = df["c_index_mean"].to_numpy(dtype=float)

    if show_errorbar and "c_index_sd" in df.columns:
        yerr = df["c_index_sd"].to_numpy(dtype=float)
        ax.errorbar(x, y, yerr=yerr, marker="o", linestyle="-", capsize=3)
    else:
        ax.plot(x, y, marker="o")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    # helpful: show failures if present
    if "n_fails" in df.columns:
        fails = df["n_fails"].to_numpy()
        if np.any(fails > 0):
            ax.text(
                0.99, 0.02,
                f"fails (total): {int(np.nansum(fails))}",
                transform=ax.transAxes,
                ha="right", va="bottom",
                color="0.35",
            )

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax


def stability_barplot(
    df_stab: pd.DataFrame,
    *,
    top_n: int = 30,
    min_freq: float | None = None,
    title: str | None = "Stability selection (gene frequency)",
    figsize: tuple[float, float] | None = None,
    save: str | Path | None = None,
    show: bool = True,
) -> tuple[plt.Figure, plt.Axes]:
    """
    Barplot for stability selection output from `bk.tl.cox_stability_selection`.

    Shows top genes by frequency, optionally filtered by min_freq.
    """
    set_style()
    if not isinstance(df_stab, pd.DataFrame):
        raise TypeError("df_stab must be a pandas DataFrame.")
    if "gene" not in df_stab.columns or "freq" not in df_stab.columns:
        raise ValueError("df_stab must contain columns: ['gene', 'freq'].")

    df = df_stab.copy()
    df["gene"] = df["gene"].astype(str)
    df["freq"] = pd.to_numeric(df["freq"], errors="coerce")

    if min_freq is not None:
        df = df[df["freq"] >= float(min_freq)]
    df = df.sort_values("freq", ascending=False).head(int(top_n))

    if df.shape[0] == 0:
        raise ValueError("No genes to plot after filtering.")

    if figsize is None:
        h = max(3.2, 0.22 * df.shape[0] + 1.2)
        figsize = (7.0, h)

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(df["gene"].iloc[::-1], df["freq"].iloc[::-1])
    ax.set_xlabel("Selection frequency")
    ax.set_ylabel("")
    ax.set_xlim(0, 1.0)

    if title:
        ax.set_title(title)

    fig.tight_layout()
    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()
    return fig, ax


########
## Signatures
#######


def permutation_importance_barplot(
    df_perm: pd.DataFrame,
    *,
    metric: Literal["cindex", "auc"] = "cindex",
    top_n: int = 20,
    sort: bool = True,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    show_values: bool = False,
    value_fontsize: float = 8.0,
    save: str | Path | None = None,
    show: bool = True,
):
    """
    Barplot of permutation importance for signature genes with error bars.

    Expects output from `bk.tl.signature_permutation_importance(...)`.

    For metric="cindex":
        uses columns:
          - gene
          - delta_cindex_mean
          - delta_cindex_sd

    For metric="auc":
        uses columns:
          - gene
          - delta_auc_mean
          - delta_auc_sd
    """
    set_style()

    if not isinstance(df_perm, pd.DataFrame):
        raise TypeError("df_perm must be a pandas DataFrame")

    if metric == "cindex":
        col_mean, col_sd = "delta_cindex_mean", "delta_cindex_sd"
        xlab = "Permutation importance (Δ C-index)"
    elif metric == "auc":
        col_mean, col_sd = "delta_auc_mean", "delta_auc_sd"
        xlab = "Permutation importance (Δ AUC)"
    else:
        raise ValueError("metric must be 'cindex' or 'auc'")

    required = {"gene", col_mean, col_sd}
    missing = required - set(df_perm.columns)
    if missing:
        raise ValueError(f"df_perm missing required columns: {sorted(missing)}")

    df = df_perm.copy()
    df["gene"] = df["gene"].astype(str)
    df[col_mean] = pd.to_numeric(df[col_mean], errors="coerce")
    df[col_sd] = pd.to_numeric(df[col_sd], errors="coerce")

    df = df[np.isfinite(df[col_mean])].copy()
    if df.shape[0] == 0:
        raise ValueError("No finite permutation-importance values to plot.")

    if sort:
        df = df.sort_values(col_mean, ascending=False)

    df = df.head(int(top_n)).copy()

    # plotting order: top at top
    df = df.iloc[::-1]

    if figsize is None:
        h = max(3.2, 0.30 * df.shape[0] + 1.4)
        figsize = (6.8, h)

    y = np.arange(df.shape[0])
    means = df[col_mean].to_numpy(float)
    sds = df[col_sd].to_numpy(float)
    genes = df["gene"].to_numpy(str)

    fig, ax = plt.subplots(figsize=figsize)

    ax.barh(y, means, xerr=sds, capsize=3, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(genes)
    ax.set_xlabel(xlab)

    if title is None:
        title = f"Permutation importance (top {min(int(top_n), int(df_perm.shape[0]))})"
    if title:
        ax.set_title(title)

    # reference line at 0
    ax.axvline(0.0, linewidth=1.0, alpha=0.4)

    if show_values:
        for yi, m in zip(y, means):
            ax.text(
                m,
                yi,
                f" {m:.3g}",
                va="center",
                ha="left" if m >= 0 else "right",
                fontsize=float(value_fontsize),
                alpha=0.85,
            )

    fig.tight_layout()

    if save is not None:
        _savefig(fig, save)
    if show:
        plt.show()

    return fig, ax




def auc_cv_heatmap(
    df_cv: pd.DataFrame,
    *,
    value_col: str = "auc_mean",
    x_col: str = "C",
    y_col: str = "l1_ratio",
    title: str = "AUC CV grid (mean AUC)",
    figsize: tuple[float, float] = (5.0, 3.8),
    annotate: bool = True,
    annot_fontsize: float = 7.0,
    save: str | None = None,
    show: bool = True,
):
    """
    Heatmap of mean AUC across (C, l1_ratio).
    Expects df_cv from recommended_auc_panel (columns: C, l1_ratio, auc_mean).
    """
    req = {x_col, y_col, value_col}
    if not req.issubset(df_cv.columns):
        raise KeyError(f"df_cv must contain {sorted(req)}")

    piv = df_cv.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc="mean")
    piv = piv.sort_index(ascending=True)
    piv = piv.loc[:, sorted(piv.columns, key=float)]

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(piv.to_numpy(), aspect="auto")

    ax.set_title(title, pad=10)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)

    ax.set_xticks(np.arange(piv.shape[1]))
    ax.set_xticklabels([f"{float(c):g}" for c in piv.columns], rotation=90)
    ax.set_yticks(np.arange(piv.shape[0]))
    ax.set_yticklabels([f"{float(r):g}" for r in piv.index])

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(value_col)

    if annotate and piv.size <= 400:
        arr = piv.to_numpy()
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                v = arr[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=float(annot_fontsize))

    fig.tight_layout()
    if save is not None:
        fig.savefig(save, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    return fig, ax, piv


def stability_freq_bar(
    df_stability: pd.DataFrame,
    *,
    top: int = 30,
    freq_col: str = "freq",
    gene_col: str = "gene",
    title: str = "Stability selection frequency",
    figsize: tuple[float, float] = (6.0, 4.5),
    show_values: bool = True,
    value_fontsize: float = 8.0,
    save: str | None = None,
    show: bool = True,
):
    """
    Horizontal barplot of stability frequencies.
    Expects columns: gene, freq (and optionally sign_consistency).
    """
    if gene_col not in df_stability.columns or freq_col not in df_stability.columns:
        raise KeyError(f"df_stability must contain columns '{gene_col}' and '{freq_col}'")

    df = df_stability.copy()
    df = df[np.isfinite(df[freq_col])].sort_values(freq_col, ascending=False).head(int(top))
    df = df.iloc[::-1]  # so best is at top

    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(df.shape[0])
    ax.barh(y, df[freq_col].to_numpy(float))

    ax.set_yticks(y)
    ax.set_yticklabels(df[gene_col].astype(str).tolist())
    ax.set_xlabel("Selection frequency")
    ax.set_title(title, pad=10)
    ax.set_xlim(0, 1.0)

    if show_values:
        for i, v in enumerate(df[freq_col].to_numpy(float)):
            ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=float(value_fontsize))

    fig.tight_layout()
    if save is not None:
        fig.savefig(save, bbox_inches="tight", dpi=300)
    if show:
        plt.show()

    return fig, ax, df



def pr_curve_signatures(
    adata,
    *,
    sigs: dict,
    label_col: str = "PFS_6m",
    positive_label: str = "NR",
    layer: str | None = "log1p_cpm",
    use_proba: bool = True,
    figsize: tuple[float, float] = (6.0, 5.0),
    title: str | None = None,
    lw: float = 2.0,
    alpha: float = 0.95,
    show_baseline: bool = True,
    baseline_ls: str = "--",
    seed: int = 0,
):
    """
    Plot PR curves for multiple signatures on the same cohort.

    Parameters
    ----------
    sigs : dict
        Mapping name -> dict with:
            - "weights": pd.DataFrame with columns ["gene","beta"]
            - optional "intercept": float (default 0.0)

        Example:
        sigs = {
          "S12": {"weights": w12, "intercept": 0.0},
          "NEW12": {"weights": w_new, "intercept": out["intercept"]},
        }

    Returns
    -------
    fig, ax, df_scores
        df_scores has per-signature PR-AUC and prevalence.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    try:
        from sklearn.metrics import precision_recall_curve, average_precision_score
    except Exception as e:
        raise ImportError(f"pr_curve_signatures requires scikit-learn. ({e})")

    # local import to avoid circulars
    from .signature import signature_score

    if label_col not in adata.obs.columns:
        raise KeyError(f"label_col='{label_col}' not found in adata.obs")

    # encode labels
    yraw = adata.obs[label_col].astype("object")
    ok_y = yraw.notna().to_numpy()
    y = (yraw.iloc[np.where(ok_y)[0]].astype(str) == str(positive_label)).astype(int).to_numpy()

    if len(np.unique(y)) < 2:
        raise ValueError("Need both classes present to plot PR curve.")

    prevalence = float(y.mean())

    def _sigmoid(z):
        z = np.clip(np.asarray(z, float), -50, 50)
        return 1.0 / (1.0 + np.exp(-z))

    rows = []
    fig, ax = plt.subplots(figsize=figsize)

    # Plot each signature
    for name, spec in sigs.items():
        if "weights" not in spec:
            raise ValueError(f"sigs['{name}'] must include 'weights'")

        w = spec["weights"]
        intercept = float(spec.get("intercept", 0.0))

        score_key = f"_sigscore_{name}_{seed}"
        signature_score(adata, weights=w, layer=layer, key_added=score_key)

        s = pd.to_numeric(adata.obs[score_key], errors="coerce").to_numpy(float)
        s = s[ok_y]

        ok = np.isfinite(s)
        yy = y[ok]
        ss = s[ok]

        if len(np.unique(yy)) < 2 or ss.size < 10:
            rows.append({"signature": name, "pr_auc": np.nan, "n": int(ss.size), "n_pos": int(yy.sum()), "prevalence": prevalence})
            continue

        pred = _sigmoid(intercept + ss) if use_proba else ss

        pr_auc = float(average_precision_score(yy, pred))
        precision, recall, _ = precision_recall_curve(yy, pred)

        ax.plot(
            recall, precision,
            linewidth=float(lw),
            alpha=float(alpha),
            label=f"{name} (PR-AUC={pr_auc:.3f})",
        )

        rows.append(
            {
                "signature": name,
                "pr_auc": pr_auc,
                "n": int(ss.size),
                "n_pos": int(yy.sum()),
                "prevalence": prevalence,
            }
        )

    # baseline prevalence line
    if show_baseline:
        ax.hlines(prevalence, 0.0, 1.0, linestyles=baseline_ls, linewidth=1.5, label=f"Baseline (prev={prevalence:.3f})")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    if title is None:
        title = f"Precision–Recall curves ({label_col}, positive='{positive_label}')"
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)

    df_scores = pd.DataFrame(rows).sort_values("pr_auc", ascending=False, na_position="last").reset_index(drop=True)
    return fig, ax, df_scores


#################################################################
### Core COX and KM plots
#################################################################

def _set_xscale_and_limits(
    ax,
    *,
    xmin,
    xmax,
    xscale="log",
    xlim=None,
):
    """
    Set x-scale and limits for Cox forest plot.
    Ensures HR=1.0 is always visible.
    Returns x_text position for p-value annotations.
    """
    eps = 1e-6

    # --- base limits from data ---
    xmin_data = max(eps, float(xmin)) if np.isfinite(xmin) else eps
    xmax_data = max(xmin_data * 1.05, float(xmax)) if np.isfinite(xmax) else xmin_data * 2

    # --- user-supplied limits ---
    if xlim is not None:
        xmin_use, xmax_use = map(float, xlim)
    else:
        xmin_use, xmax_use = xmin_data, xmax_data

    # --- ALWAYS include HR = 1.0 ---
    xmin_use = min(xmin_use, 1.0)
    xmax_use = max(xmax_use, 1.0)

    if xscale == "log":
        ax.set_xscale("log")

        # guard again for log scale
        xmin_use = max(eps, xmin_use)
        xmax_use = max(xmin_use * 1.05, xmax_use)

        ax.set_xlim(xmin_use, xmax_use)

        # annotation position
        x_text = xmax_use * 1.05

    else:
        ax.set_xscale("linear")

        ax.set_xlim(xmin_use, xmax_use)
        x_text = xmax_use * 1.02

    return x_text


def cox_forest_from_uns(
    adata,
    *,
    uns_key: str = "cox",
    uns: dict | None = None,                 # : accept payload directly
    title: Optional[str] = None,
    x_keys: Optional[Sequence[str]] = None,
    xscale: str = "log",                     # : applied ("log"|"linear")
    rename: Optional[dict] = None,
    sort_by: str = "HR",                     # "HR" | "p" | "name"
    xlim: Optional[tuple[float, float]] = None,   
    figsize: tuple[float, float] | None = None,  
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Forest plot for Cox results stored in `adata.uns[uns_key]` OR passed via `uns=`.

    Supports two formats:
    1) Multi-x tidy table produced by `bk.tl.cox_univariate(..., x_keys=[...])`
       -> payload["summary"] columns:
          x_key, term, HR, HR_lower_95, HR_upper_95, p, n_used, events

    2) Single-model summary produced by `bk.tl.cox_interaction()` (lifelines style)
       -> payload["summary"] is a DataFrame indexed by covariate

    Returns matplotlib Axes.
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # -----------------------
    # get payload
    # -----------------------
    if uns is None:
        if adata is None:
            raise ValueError("Provide either `adata` (with adata.uns[uns_key]) or `uns=...`.")
        if uns_key not in adata.uns:
            raise KeyError(f"{uns_key!r} not found in adata.uns")
        payload = adata.uns[uns_key]
    else:
        payload = uns

    if not isinstance(payload, dict) or "summary" not in payload:
        raise ValueError("payload/uns must be a dict containing key 'summary'")

    summ = payload["summary"]
    if not isinstance(summ, pd.DataFrame):
        summ = pd.DataFrame(summ)

    xscale = str(xscale).lower()
    if xscale not in ("log", "linear"):
        raise ValueError("xscale must be 'log' or 'linear'")

    def _fmt_p(p):
        try:
            p = float(p)
        except Exception:
            return ""
        if np.isnan(p):
            return ""
        if p < 1e-4:
            return "p<1e-4"
        if p < 1e-3:
            return "p<1e-3"
        return f"p={p:.3g}"

    def _apply_scale_and_limits(ax, lo, hi):
        lo = np.asarray(lo, float)
        hi = np.asarray(hi, float)

        # robust xmin/xmax
        x_text = _set_xscale_and_limits(
            ax,
            xmin=np.nanmin(lo),
            xmax=np.nanmax(hi),
            xscale=xscale,
            xlim=xlim,
        )

        return x_text

    # -----------------------
    # Case A: multi-x tidy table
    # -----------------------
    if {"x_key", "HR", "HR_lower_95", "HR_upper_95", "p"}.issubset(summ.columns):
        df = summ.copy()
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=["HR", "HR_lower_95", "HR_upper_95"])

        # subset
        if x_keys is not None:
            want = [str(k) for k in x_keys]
            have = set(df["x_key"].astype(str).unique())
            missing = [k for k in want if k not in have]
            if missing:
                raise KeyError(f"x_keys not found in summary: {missing}")
            df = df[df["x_key"].astype(str).isin(want)].copy()

        # labels
        labels = df["x_key"].astype(str).tolist()
        if rename:
            labels = [rename.get(x, x) for x in labels]

        # sorting
        if sort_by == "HR":
            order = np.argsort(df["HR"].to_numpy(float))[::-1]
        elif sort_by == "p":
            order = np.argsort(df["p"].to_numpy(float))
        elif sort_by == "name":
            order = np.argsort(np.array(labels, dtype=str))
        else:
            raise ValueError("sort_by must be one of {'HR','p','name'}")

        df = df.iloc[order].reset_index(drop=True)
        labels = [labels[i] for i in order]

        HR = df["HR"].to_numpy(float)
        lo = df["HR_lower_95"].to_numpy(float)
        hi = df["HR_upper_95"].to_numpy(float)
        ptxt = [_fmt_p(p) for p in df["p"].to_numpy(float)]

        y = np.arange(len(df))[::-1]

        # figsize
        if figsize is None:
            fig_h = max(2.2, 0.45 * len(df) + 1.0)
            figsize = (6.8, fig_h)

        fig, ax = plt.subplots(figsize=figsize)
        xerr = np.vstack([HR - lo, hi - HR])
        ax.errorbar(HR, y, xerr=xerr, fmt="o")
        ax.axvline(1.0, linestyle="--")

        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Hazard ratio (95% CI)")

        params = payload.get("params", {})
        if title is None:
            time_key = params.get("time_key", "")
            strata = params.get("strata")
            title = f"Cox univariate comparison ({time_key})" if time_key else "Cox univariate comparison"
            if strata:
                title += f" | stratified by {strata}"
        ax.set_title(title)

        x_text = _apply_scale_and_limits(ax, lo, hi)
        if x_text is not None:
            for yi, t in zip(y, ptxt):
                ax.text(x_text, yi, t, va="center", fontsize=9)

        plt.tight_layout()
        if save:
            fig.savefig(save, bbox_inches="tight")
        if show:
            plt.show()
        return ax

    # -----------------------
    # Case B: single-model (lifelines-like) summary table
    # -----------------------
    if summ.index.name is None and "covariate" in summ.columns:
        summ = summ.set_index("covariate")

    required = {"HR", "HR_lower_95", "HR_upper_95", "p"}
    if not required.issubset(set(summ.columns)):
        raise ValueError(f"single-model summary must contain columns: {required}")

    df = summ.copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["HR", "HR_lower_95", "HR_upper_95"])

    labels = df.index.astype(str).tolist()
    if rename:
        labels = [rename.get(t, t) for t in labels]

    # sorting for single-model too
    if sort_by == "HR":
        order = np.argsort(df["HR"].to_numpy(float))[::-1]
    elif sort_by == "p":
        order = np.argsort(df["p"].to_numpy(float))
    elif sort_by == "name":
        order = np.argsort(np.array(labels, dtype=str))
    else:
        raise ValueError("sort_by must be one of {'HR','p','name'}")

    df = df.iloc[order].copy()
    labels = [labels[i] for i in order]

    HR = df["HR"].to_numpy(float)
    lo = df["HR_lower_95"].to_numpy(float)
    hi = df["HR_upper_95"].to_numpy(float)
    ptxt = [_fmt_p(p) for p in df["p"].to_numpy(float)]
    y = np.arange(len(df))[::-1]

    if figsize is None:
        fig_h = max(2.2, 0.45 * len(df) + 1.0)
        figsize = (6.8, fig_h)

    fig, ax = plt.subplots(figsize=figsize)
    xerr = np.vstack([HR - lo, hi - HR])
    ax.errorbar(HR, y, xerr=xerr, fmt="o")
    ax.axvline(1.0, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Hazard ratio (95% CI)")
    ax.set_title(title or "Cox model")

    x_text = _apply_scale_and_limits(ax, lo, hi)
    if x_text is not None:
        for yi, t in zip(y, ptxt):
            ax.text(x_text, yi, t, va="center", fontsize=9)

    plt.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def run_cox_per_group(
    adata,
    *,
    group_col: str,
    min_n: int = 30,
    dropna_group: bool = True,
    time_key: str = "OS.time",
    event_key: str = "OS",
    x_keys=(
        "mp_heterogeneity_entropy",
        "mp_dispersion_mad_zwithin_Project_ID",
        "Neuroendocrine_score",
        "Proliferation_score",
    ),
    x_mode: str = "continuous",
    strata: str | None = None,
    out_key_prefix: str = "cox_",
    verbose: bool = True,
    **cox_kwargs,
):
    """
    Run tl.cox_univariate() separately within each level of `group_col`.

    Returns dict:
      {group_value: payload_dict}  where payload_dict is the same object stored in adata.uns[out_key].
    """
    import numpy as np
    import pandas as pd

    # IMPORTANT: import the function directly (no 'bk' global)
    # Adjust the relative import to your actual module layout:
    from ..tl.cox import cox_univariate

    if group_col not in adata.obs.columns:
        raise KeyError(f"group_col='{group_col}' not found in adata.obs")

    g = adata.obs[group_col]
    if dropna_group:
        g = g.dropna()

    results = {}

    # group -> list of obs_names
    groups = g.astype(str).groupby(g.astype(str)).groups

    for grp, idx in groups.items():
        idx = list(idx)
        if len(idx) < int(min_n):
            if verbose:
                print(f"[SKIP] {grp}: n={len(idx)} < min_n={min_n}")
            continue

        ad = adata[idx].copy()

        # survival columns present?
        missing_cols = [k for k in (time_key, event_key) if k not in ad.obs.columns]
        if missing_cols:
            if verbose:
                print(f"[WARN] {grp}: missing {missing_cols} -> skipped")
            continue

        # within-group: strata cannot be the same as group_col
        strata_eff = strata
        if strata_eff is not None and str(strata_eff) == str(group_col):
            if verbose:
                print(f"[WARN] {grp}: strata==group_col inside subset -> setting strata=None")
            strata_eff = None

        out_key = f"{out_key_prefix}{grp}"

        try:
            cox_univariate(
                ad,
                time_key=time_key,
                event_key=event_key,
                x_keys=list(x_keys),
                x_mode=x_mode,
                strata=strata_eff,
                out_key=out_key,
                **cox_kwargs,
            )

            payload = ad.uns.get(out_key, None)
            if payload is None:
                if verbose:
                    print(f"[WARN] {grp}: no payload found at ad.uns[{out_key!r}]")
                continue

            results[grp] = payload
            if verbose:
                print(f"[OK] {grp}: n={len(idx)} stored '{out_key}'")

        except Exception as e:
            if verbose:
                print(f"[WARN] {grp}: {type(e).__name__}: {e}")

    return results



def metaprogram_rank1_composition_stackedbar(
    adata: AnnData,
    *,
    topk_key: str = "mp_topk",
    groupby: str = "Project_ID",
    rank: int = 1,
    top_groups: int = 12,
    top_mps: int = 12,
    figsize=(6, 4),
    legend: bool = True,
    legend_bbox_to_anchor=(2.2, 1.0),
    legend_loc: str = "upper right",
    rotation: int = 90,
    title: str = "Rank-1 metaprogram composition across tumor types",
    ylabel: str = "Fraction of samples where MP is rank-1",
    show: bool = True,
    save: Optional[str] = None,
    out_key: Optional[str] = "mp_rank1_composition",
):
    """
    Stacked barplot showing, for each tumor type (groupby), the fraction of samples
    whose rank-1 metaprogram is each MP.

    Uses `adata.uns[topk_key]` with at least columns: 'sample', 'rank', 'metaprogram'.
    Stores the plotted matrix in `adata.uns[out_key]` if provided.
    """
    if topk_key not in adata.uns:
        raise KeyError(f"{topk_key!r} not found in adata.uns")
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    topk = adata.uns[topk_key].copy()
    if "rank" not in topk.columns or "sample" not in topk.columns or "metaprogram" not in topk.columns:
        raise ValueError(f"adata.uns[{topk_key!r}] must contain columns: 'sample', 'rank', 'metaprogram'")

    topk = topk[topk["rank"] == rank].copy()

    # map sample -> groupby
    topk[groupby] = adata.obs.loc[topk["sample"], groupby].astype(str).to_numpy()

    # pick top groups by sample count (from full obs)
    top_group_names = adata.obs[groupby].astype(str).value_counts().head(top_groups).index.astype(str)
    topk = topk[topk[groupby].isin(top_group_names)]

    # compute fractions (rank-1 metaprogram composition)
    mat = (
        topk.groupby(groupby)["metaprogram"]
            .value_counts(normalize=True)
            .unstack(fill_value=0)
    )

    # order groups (largest first, based on full obs)
    order = (
        adata.obs[groupby].astype(str)
            .value_counts()
            .loc[list(top_group_names)]
            .sort_values(ascending=False)
            .index.astype(str)
    )
    mat = mat.loc[order]

    # keep only top MPs globally for readability
    top_mp_names = mat.sum(axis=0).sort_values(ascending=False).head(top_mps).index
    mat_plot = mat[top_mp_names]

    ax = mat_plot.plot(kind="bar", stacked=True, figsize=figsize, legend=legend)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=rotation)

    if legend:
        plt.legend(loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor)

    plt.tight_layout()

    if save:
        plt.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    if out_key:
        adata.uns[out_key] = {
            "params": dict(
                topk_key=topk_key,
                groupby=groupby,
                rank=rank,
                top_groups=top_groups,
                top_mps=top_mps,
            ),
            "mat": mat_plot,
        }

    return ax


###########################################################
### Kaplan-Meier curves
###########################################################

def km_univariate(
    adata: AnnData,
    *,
    x_key: str = "mp_heterogeneity_entropy",
    x_source: Literal["auto", "obs", "gene"] = "auto",
    layer: str | None = None,
    time_key: str = "OS.time",
    event_key: str = "OS",
    group_col: str | None = "surv_group_1d",
    groupby_for_binning: Optional[str] = None,
    binning: str = "global",
    q: Tuple[float, float] = (0.25, 0.75),
    keep: str = "extremes",
    labels: Tuple[str, str] = ("low", "high"),
    min_n: int = 40,
    figsize=(6, 4),
    ylim: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    logrank: bool = False,
    annotate_p: bool = True,
    show: bool = True,
    save: Optional[str] = None,
    out_key: Optional[str] = "km_1d",
    force_rebin: bool = False,
    cleanup_tmp_cols: bool = True,

    # NEW: binary handling + colors
    binary_mode: Literal["auto", "force", "off"] = "auto",
    binary_labels: Tuple[str, str] = ("0", "1"),   # e.g. ("WT","Mut")
    palette: Optional[object] = None,              # tuple/list/dict; see below
):
    """
    Kaplan–Meier plot for one variable (x_key).

    Behavior:
      - If x is binary (two unique values) and binary_mode in {"auto","force"},
        groups are taken directly from x (no quantile binning).
      - Otherwise, uses surv_1d_bins() (quantile-based low/high).

    palette:
      - None -> use matplotlib defaults
      - tuple/list (len>=2) -> mapped to [low, high] (or [0,1] for binary)
      - dict -> explicit mapping {group_label: color}
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    def _as_1d_gene(_adata, gene: str, _layer: str | None):
        if gene not in _adata.var_names:
            raise KeyError(f"Gene '{gene}' not found in adata.var_names")
        X = _adata.layers[_layer] if (_layer is not None and _layer in _adata.layers) else _adata.X
        j = _adata.var_names.get_loc(gene)
        v = X[:, j]
        try:
            import scipy.sparse as sp
            if sp.issparse(v):
                v = v.toarray()
        except Exception:
            pass
        return np.asarray(v).ravel().astype(float)

    def _resolve_source(key: str, source: str):
        if source == "obs":
            if key not in adata.obs.columns:
                raise KeyError(f"'{key}' not found in adata.obs (x_source='obs').")
            return "obs"
        if source == "gene":
            if key not in adata.var_names:
                raise KeyError(f"Gene '{key}' not found in adata.var_names (x_source='gene').")
            return "gene"
        if key in adata.obs.columns:
            return "obs"
        if key in adata.var_names:
            return "gene"
        raise KeyError(f"'{key}' not found in adata.obs.columns or adata.var_names (source='auto').")

    def _is_binary_like(v: np.ndarray) -> bool:
        v = v[np.isfinite(v)]
        if v.size == 0:
            return False
        u = np.unique(v)
        if u.size != 2:
            return False
        # typical cases: {0,1} or {False,True} or {0,2} etc.
        return True

    def _pick_color(group_label: str, idx: int, default=None):
        if palette is None:
            return default
        if isinstance(palette, dict):
            return palette.get(group_label, default)
        # tuple/list-like
        try:
            return list(palette)[idx]
        except Exception:
            return default

    # If group_col None -> make it unique per key
    if group_col is None:
        group_col = f"surv_group_1d__{x_key}"

    # resolve x
    tmp_cols = []
    x_src = _resolve_source(x_key, x_source)
    x_obs_key = x_key
    if x_src == "gene":
        x_obs_key = f"__gene__{x_key}"
        adata.obs[x_obs_key] = _as_1d_gene(adata, x_key, layer)
        tmp_cols.append(x_obs_key)

    # Decide if we do binary direct grouping
    x_raw = pd.to_numeric(adata.obs[x_obs_key], errors="coerce").to_numpy(float)
    binary_direct = False
    if binary_mode == "force":
        binary_direct = True
    elif binary_mode == "auto":
        binary_direct = _is_binary_like(x_raw)

    # Determine if we need to (re)create group_col
    need_rebin = (
        force_rebin
        or (group_col not in adata.obs)
        or adata.obs[group_col].isna().all()
    )

    # also rebin if prior params mismatch
    if (not need_rebin) and out_key and (out_key in adata.uns) and isinstance(adata.uns[out_key], dict):
        prev = adata.uns[out_key].get("params", {})
        # if previously binary vs quantile mode differs, rebin
        if prev.get("binary_direct") != bool(binary_direct):
            need_rebin = True
        if (
            prev.get("x_key") != x_key
            or prev.get("binning") != binning
            or tuple(prev.get("q", (None, None))) != tuple(q)
            or prev.get("keep") != keep
            or tuple(prev.get("labels", (None, None))) != tuple(labels)
            or prev.get("groupby_for_binning") != groupby_for_binning
        ):
            need_rebin = True

    # Create groups
    if need_rebin:
        if binary_direct:
            # direct binary grouping
            s = pd.to_numeric(adata.obs[x_obs_key], errors="coerce")
            # keep only finite
            m = np.isfinite(s.to_numpy(float))
            # use the two observed values in sorted order
            u = np.unique(s.to_numpy(float)[m])
            if u.size < 2:
                # fallback: everything NA or single value
                adata.obs[group_col] = pd.NA
            else:
                # map lower->label0, higher->label1
                lo_val, hi_val = float(np.min(u)), float(np.max(u))
                lab0, lab1 = binary_labels
                grp = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
                grp.loc[m & (s.to_numpy(float) == lo_val)] = str(lab0)
                grp.loc[m & (s.to_numpy(float) == hi_val)] = str(lab1)
                adata.obs[group_col] = grp.astype("category")
        else:
            # quantile binning
            from bullkpy.tl.cox import surv_1d_bins
            surv_1d_bins(
                adata,
                x_key=x_obs_key,
                time_key=time_key,
                event_key=event_key,
                groupby=groupby_for_binning,
                binning=binning,
                q=q,
                keep=keep,
                labels=labels,
                out_col=group_col,
                copy=False,
            )

    # Build plotting df
    obs = adata.obs
    for k in [time_key, event_key, group_col]:
        if k not in obs:
            raise KeyError(f"{k!r} not found in adata.obs")

    df = obs[[time_key, event_key, group_col]].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[time_key, event_key, group_col])
    df = df[df[time_key].astype(float) > 0]
    df[event_key] = df[event_key].astype(float).astype(int)

    # filter small groups
    sizes = df[group_col].value_counts()
    keep_groups = sizes[sizes >= min_n].index.tolist()
    df = df[df[group_col].isin(keep_groups)].copy()

    if len(keep_groups) < 2:
        # warn early; still plot what exists
        print(f"[WARN] km_univariate: only {len(keep_groups)} group(s) remain after min_n={min_n}. sizes={sizes.to_dict()}")

    # Order + colors
    if binary_direct:
        lab0, lab1 = map(str, binary_labels)
        canonical = [lab0, lab1]
    else:
        low_label, high_label = map(str, labels)
        canonical = [low_label, high_label]

    order = [g for g in canonical if g in keep_groups] + [g for g in keep_groups if g not in canonical]

    fig, ax = plt.subplots(figsize=figsize)
    kmf = KaplanMeierFitter()

    for i, g in enumerate(order):
        sub = df[df[group_col].astype(str) == str(g)]
        color = _pick_color(str(g), i, default=None)
        label = f"{g} (n={sub.shape[0]})"
        kmf.fit(durations=sub[time_key], event_observed=sub[event_key], label=label)
        kmf.plot(ax=ax, ci_show=False, color=color)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title or (f"KM: {x_key} (binary)" if binary_direct else f"KM: {x_key} (low vs high)"))
    ax.legend(frameon=False, fontsize=8)

    if ylim is not None:
        ax.set_ylim(*ylim)

    # log-rank (only meaningful when exactly 2 groups survive)
    logrank_res = None
    if logrank and df[group_col].nunique() == 2:
        try:
            from lifelines.statistics import logrank_test
            g1, g2 = df[group_col].astype(str).unique().tolist()[:2]
            d1 = df[df[group_col].astype(str) == str(g1)]
            d2 = df[df[group_col].astype(str) == str(g2)]
            logrank_res = logrank_test(
                d1[time_key],
                d2[time_key],
                event_observed_A=d1[event_key],
                event_observed_B=d2[event_key],
            )
            if annotate_p:
                p = float(logrank_res.p_value)
                txt = f"log-rank p = {p:.2e}" if p < 1e-3 else f"log-rank p = {p:.3f}"
                ax.text(0.98, 0.02, txt, ha="right", va="bottom", transform=ax.transAxes, fontsize=9)
        except Exception:
            logrank_res = None

    plt.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    if out_key:
        adata.uns[out_key] = {
            "params": dict(
                x_key=x_key,
                x_source=x_src,
                layer=layer,
                time_key=time_key,
                event_key=event_key,
                group_col=group_col,
                groupby_for_binning=groupby_for_binning,
                binning=binning,
                q=q,
                keep=keep,
                labels=labels,
                min_n=min_n,
                logrank=logrank,
                force_rebin=force_rebin,
                binary_mode=binary_mode,
                binary_direct=bool(binary_direct),
                binary_labels=binary_labels,
            ),
            "group_sizes": df[group_col].value_counts().to_dict(),
            "logrank_p": None if logrank_res is None else float(logrank_res.p_value),
            "df_used": df,
        }

    if cleanup_tmp_cols and tmp_cols:
        for c in tmp_cols:
            if c in adata.obs.columns:
                del adata.obs[c]

    return ax


def km_2x2_interaction(
    adata,
    *,
    x_key: str = "Neuroendocrine_score",
    z_key: str = "mp_heterogeneity_entropy",
    x_source: Literal["auto", "obs", "gene"] = "auto",
    z_source: Literal["auto", "obs", "gene"] = "auto",
    layer: str | None = None,
    time_key: str = "OS.time",
    event_key: str = "OS",
    group_col: str | None = "surv_group_2x2",
    groupby_for_binning: Optional[str] = None,
    binning: str = "global",  # "global" | "within_group"
    q: Tuple[float, float] = (0.25, 0.75),
    keep: str = "extremes",
    min_n: int = 40,
    figsize=(6, 4),
    ylim: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    logrank: bool = False,
    annotate_p: bool = True,
    show: bool = True,
    save: Optional[str] = None,
    out_key: Optional[str] = "km_2x2",
    force_rebin: bool = False,
    cleanup_tmp_cols: bool = True,

    # NEW: binary handling per axis
    x_binary_mode: Literal["auto", "force", "off"] = "auto",
    z_binary_mode: Literal["auto", "force", "off"] = "auto",
    x_binary_labels: Tuple[str, str] = ("0", "1"),  # e.g. ("WT","Mut")
    z_binary_labels: Tuple[str, str] = ("0", "1"),

    # NEW: color selection
    # - dict: {group_label: color}
    # - list/tuple len>=4: mapped to canonical order
    palette: Optional[object] = None,
):
    """
    Kaplan–Meier plot for a 2x2 interaction: x low/high × z low/high

    Supports:
      - x_key / z_key as obs or genes (var_names)
      - binary variables: directly use their two levels (no quantile binning)
      - palette control for the 4 curves
      - force_rebin to avoid reusing stale bins/group_col
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    def _as_1d_gene(_adata, gene: str, _layer: str | None):
        if gene not in _adata.var_names:
            raise KeyError(f"Gene '{gene}' not found in adata.var_names")
        X = _adata.layers[_layer] if (_layer is not None and _layer in _adata.layers) else _adata.X
        j = _adata.var_names.get_loc(gene)
        v = X[:, j]
        try:
            import scipy.sparse as sp
            if sp.issparse(v):
                v = v.toarray()
        except Exception:
            pass
        return np.asarray(v).ravel().astype(float)

    def _resolve_source(key: str, source: str):
        if source == "obs":
            if key not in adata.obs.columns:
                raise KeyError(f"'{key}' not found in adata.obs (source='obs').")
            return "obs"
        if source == "gene":
            if key not in adata.var_names:
                raise KeyError(f"Gene '{key}' not found in adata.var_names (source='gene').")
            return "gene"
        if key in adata.obs.columns:
            return "obs"
        if key in adata.var_names:
            return "gene"
        raise KeyError(f"'{key}' not found in adata.obs.columns or adata.var_names (source='auto').")

    def _is_binary_like(arr: np.ndarray) -> bool:
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return False
        u = np.unique(arr)
        return u.size == 2

    def _make_binary_bins(values: np.ndarray, labels: Tuple[str, str]) -> pd.Series:
        """Map the two observed numeric values to labels[0]/labels[1] (low->0, high->1)."""
        s = pd.Series(values, index=adata.obs_names, dtype=float)
        m = np.isfinite(s.to_numpy(float))
        u = np.unique(s.to_numpy(float)[m])
        if u.size < 2:
            return pd.Series(pd.NA, index=adata.obs_names, dtype="object")
        lo, hi = float(np.min(u)), float(np.max(u))
        lab0, lab1 = map(str, labels)
        out = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
        out.loc[m & (s.to_numpy(float) == lo)] = lab0
        out.loc[m & (s.to_numpy(float) == hi)] = lab1
        return out

    def _pick_color(group_label: str, idx: int):
        if palette is None:
            return None
        if isinstance(palette, dict):
            return palette.get(group_label, None)
        try:
            pal = list(palette)
            return pal[idx] if idx < len(pal) else None
        except Exception:
            return None

    # If group_col None -> make it unique per pair
    if group_col is None:
        group_col = f"surv_group_2x2__{x_key}__{z_key}"

    # Determine if we need to rebin
    need_rebin = (
        force_rebin
        or (group_col not in adata.obs)
        or adata.obs[group_col].isna().all()
    )

    # Also rebin if prior params mismatch
    if (not need_rebin) and out_key and (out_key in adata.uns) and isinstance(adata.uns[out_key], dict):
        prev = adata.uns[out_key].get("params", {})
        if (
            prev.get("x_key") != x_key
            or prev.get("z_key") != z_key
            or prev.get("binning") != binning
            or tuple(prev.get("q", (None, None))) != tuple(q)
            or prev.get("keep") != keep
            or prev.get("groupby_for_binning") != groupby_for_binning
            or prev.get("x_binary_mode") != x_binary_mode
            or prev.get("z_binary_mode") != z_binary_mode
            or tuple(prev.get("x_binary_labels", (None, None))) != tuple(x_binary_labels)
            or tuple(prev.get("z_binary_labels", (None, None))) != tuple(z_binary_labels)
        ):
            need_rebin = True

    tmp_cols = []

    # Resolve sources and create obs columns if genes
    x_src = _resolve_source(x_key, x_source)
    z_src = _resolve_source(z_key, z_source)

    x_obs_key = x_key
    z_obs_key = z_key

    if x_src == "gene":
        x_obs_key = f"__gene__{x_key}"
        adata.obs[x_obs_key] = _as_1d_gene(adata, x_key, layer)
        tmp_cols.append(x_obs_key)

    if z_src == "gene":
        z_obs_key = f"__gene__{z_key}"
        adata.obs[z_obs_key] = _as_1d_gene(adata, z_key, layer)
        tmp_cols.append(z_obs_key)

    # Decide binary-direct per axis
    x_raw = pd.to_numeric(adata.obs[x_obs_key], errors="coerce").to_numpy(float)
    z_raw = pd.to_numeric(adata.obs[z_obs_key], errors="coerce").to_numpy(float)

    x_binary_direct = (x_binary_mode == "force") or (x_binary_mode == "auto" and _is_binary_like(x_raw))
    z_binary_direct = (z_binary_mode == "force") or (z_binary_mode == "auto" and _is_binary_like(z_raw))

    # Build / rebuild group_col if needed
    if need_rebin:
        # Create temporary bin columns
        x_bin_col = f"__x_bin_2x2__{x_key}"
        z_bin_col = f"__z_bin_2x2__{z_key}"
        tmp_cols.extend([x_bin_col, z_bin_col])

        if x_binary_direct:
            adata.obs[x_bin_col] = _make_binary_bins(x_raw, x_binary_labels).astype("category")
            x_low, x_high = map(str, x_binary_labels)
        else:
            # quantile bins -> "low"/"high" from your labels tuple
            from bullkpy.tl.cox import surv_1d_bins
            surv_1d_bins(
                adata,
                x_key=x_obs_key,
                time_key=time_key,
                event_key=event_key,
                groupby=groupby_for_binning,
                binning=binning,
                q=q,
                keep=keep,
                labels=("low", "high"),
                out_col=x_bin_col,
                copy=False,
            )
            x_low, x_high = "low", "high"

        if z_binary_direct:
            adata.obs[z_bin_col] = _make_binary_bins(z_raw, z_binary_labels).astype("category")
            z_low, z_high = map(str, z_binary_labels)
        else:
            from bullkpy.tl.cox import surv_1d_bins
            surv_1d_bins(
                adata,
                z_key=z_obs_key,
                time_key=time_key,
                event_key=event_key,
                groupby=groupby_for_binning,
                binning=binning,
                q=q,
                keep=keep,
                labels=("low", "high"),
                out_col=z_bin_col,
                copy=False,
            )
            z_low, z_high = "low", "high"

        # Combine into 2x2 group labels (canonical)
        xb = adata.obs[x_bin_col].astype(str)
        zb = adata.obs[z_bin_col].astype(str)

        out = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
        # only define groups where both bins are present
        m = xb.notna() & zb.notna()
        out.loc[m & (xb == x_low)  & (zb == z_low)]  = "X_low & Z_low"
        out.loc[m & (xb == x_low)  & (zb == z_high)] = "X_low & Z_high"
        out.loc[m & (xb == x_high) & (zb == z_low)]  = "X_high & Z_low"
        out.loc[m & (xb == x_high) & (zb == z_high)] = "X_high & Z_high"

        adata.obs[group_col] = out.astype("category")

    # Build plotting df
    obs = adata.obs
    for k in [time_key, event_key, group_col]:
        if k not in obs:
            raise KeyError(f"{k!r} not found in adata.obs")

    df = obs[[time_key, event_key, group_col]].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[time_key, event_key, group_col])
    df = df[df[time_key].astype(float) > 0]
    df[event_key] = df[event_key].astype(float).astype(int)

    # filter small groups
    sizes = df[group_col].value_counts()
    keep_groups = sizes[sizes >= min_n].index.tolist()
    df = df[df[group_col].isin(keep_groups)].copy()

    if len(keep_groups) < 2:
        print(f"[WARN] km_2x2_interaction: only {len(keep_groups)} group(s) remain after min_n={min_n}. sizes={sizes.to_dict()}")

    canonical = ["X_low & Z_low", "X_low & Z_high", "X_high & Z_low", "X_high & Z_high"]
    order = [g for g in canonical if g in keep_groups] + [g for g in keep_groups if g not in canonical]

    fig, ax = plt.subplots(figsize=figsize)
    kmf = KaplanMeierFitter()

    for i, g in enumerate(order):
        sub = df[df[group_col].astype(str) == str(g)]
        color = _pick_color(str(g), i)
        label = f"{g} (n={sub.shape[0]})"
        kmf.fit(durations=sub[time_key], event_observed=sub[event_key], label=label)
        kmf.plot(ax=ax, ci_show=False, color=color)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title or f"KM: {x_key} × {z_key}")
    ax.legend(frameon=False, fontsize=8)

    if ylim is not None:
        ax.set_ylim(*ylim)

    # global log-rank across groups
    logrank_p = None
    if logrank and df[group_col].nunique() >= 2:
        try:
            from lifelines.statistics import multivariate_logrank_test
            res = multivariate_logrank_test(
                df[time_key].astype(float),
                df[group_col].astype(str),
                df[event_key].astype(int),
            )
            logrank_p = float(res.p_value)
            if annotate_p:
                txt = f"global log-rank p = {logrank_p:.2e}" if logrank_p < 1e-3 else f"global log-rank p = {logrank_p:.3f}"
                ax.text(0.98, 0.02, txt, ha="right", va="bottom", transform=ax.transAxes, fontsize=9)
        except Exception:
            logrank_p = None

    plt.tight_layout()
    if save:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()

    if out_key:
        adata.uns[out_key] = {
            "params": dict(
                x_key=x_key,
                z_key=z_key,
                x_source=x_src,
                z_source=z_src,
                layer=layer,
                time_key=time_key,
                event_key=event_key,
                group_col=group_col,
                groupby_for_binning=groupby_for_binning,
                binning=binning,
                q=q,
                keep=keep,
                min_n=min_n,
                logrank=logrank,
                force_rebin=force_rebin,
                x_binary_mode=x_binary_mode,
                z_binary_mode=z_binary_mode,
                x_binary_labels=x_binary_labels,
                z_binary_labels=z_binary_labels,
                x_binary_direct=bool(x_binary_direct),
                z_binary_direct=bool(z_binary_direct),
            ),
            "group_sizes": df[group_col].value_counts().to_dict(),
            "logrank_p": logrank_p,
            "df_used": df,
        }

    if cleanup_tmp_cols and tmp_cols:
        for c in tmp_cols:
            if c in adata.obs.columns:
                del adata.obs[c]

    return ax