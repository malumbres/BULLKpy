from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData

from ._style import set_style, _savefig

from ..logging import warn

from sklearn.metrics import precision_recall_curve, average_precision_score



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

def _sigmoid(z):
    z = np.clip(np.asarray(z, float), -50, 50)
    return 1.0 / (1.0 + np.exp(-z))


def _y_from_adata(adata, label_col, positive_label):
    """
    Return y_full aligned to adata.obs (length = adata.n_obs) and ok_y mask.
    """
    yraw = adata.obs[label_col].astype("object")
    ok_y = yraw.notna().to_numpy()

    # y_full is always length n_obs; undefined labels become 0 (they'll be masked by ok_y anyway)
    y_full = (yraw.astype(str).to_numpy() == str(positive_label)).astype(int)

    return ok_y, y_full


def _score_from_weights(adata, weights_df, layer="log1p_cpm"):
    # uses bullkpy signature_score (recommended)
    import bullkpy as bk
    tmp_key = "__tmp_sig_score__"
    bk.tl.signature_score(adata, weights=weights_df, layer=layer, key_added=tmp_key)
    s_full = pd.to_numeric(adata.obs[tmp_key], errors="coerce").to_numpy(float)
    return s_full


def plot_pr_curves_signatures(
    adata,
    signatures: dict,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    use_proba=True,
    title=None,
    figsize=(6, 5),          # NEW
):
    """
    signatures: dict name -> dict with keys:
        - "weights_df" (pd.DataFrame gene,beta)
        - "intercept" (float) optional if use_proba=True
    """
    import matplotlib.pyplot as plt
    from sklearn.metrics import average_precision_score, precision_recall_curve

    ok_y, y_full = _y_from_adata(adata, label_col, positive_label)

    # Use only labeled samples to check that both classes exist
    y_labeled = y_full[ok_y]
    if len(np.unique(y_labeled)) < 2:
        raise ValueError("Need 2 classes to plot PR curve.")

    prevalence = float(y_labeled.mean())

    plt.figure(figsize=figsize)   # <-- use figsize here

    # prevalence baseline
    plt.plot(
        [0, 1], [prevalence, prevalence],
        linestyle="--",
        color="grey",
        label=f"Baseline (prev={prevalence:.3f})"
    )

    rows = []
    for name, sig in signatures.items():
        w = sig["weights_df"]
        b0 = float(sig.get("intercept", 0.0))

        s_full = _score_from_weights(adata, w, layer=layer)

        # single consistent mask in full adata space
        ok = ok_y & np.isfinite(s_full)
        yy = y_full[ok]
        ss = s_full[ok]

        if yy.size < 2 or len(np.unique(yy)) < 2:
            rows.append({
                "Signature": name,
                "pr_auc": np.nan,
                "n": int(yy.size),
                "pos": int(yy.sum()),
                "neg": int((yy == 0).sum())
            })
            continue

        if use_proba:
            pred = _sigmoid(b0 + ss)
        else:
            pred = ss

        ap = float(average_precision_score(yy, pred))
        prec, rec, _ = precision_recall_curve(yy, pred)

        plt.plot(rec, prec, linewidth=2, label=f"{name} (PR-AUC={ap:.3f})")

        rows.append({
            "Signature": name,
            "pr_auc": ap,
            "n": int(len(yy)),
            "pos": int(yy.sum()),
            "neg": int((yy == 0).sum())
        })

    plt.xlabel("Recall")
    plt.ylabel("Precision (PPV)")
    plt.ylim(0, 1.02)
    plt.xlim(0, 1.0)

    # ---- legend OUTSIDE ----
    plt.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False
    )

    if title is not None:
        plt.title(title)

    # leave room on the right for legend
    plt.tight_layout(rect=[0, 0, 0.78, 1])

    return (
        pd.DataFrame(rows)
        .sort_values("pr_auc", ascending=False)
        .reset_index(drop=True)
    )


#######################################################################
## BOOTSTRAP
#######################################################################



def _pr_curve_interp(y_true, y_score, recall_grid):
    """
    Compute PR curve and interpolate precision onto a fixed recall grid.
    Returns precision_on_grid (shape = len(recall_grid)).
    """
    prec, rec, _ = precision_recall_curve(y_true, y_score)

    # sklearn PR curve can have recall either increasing or decreasing depending on version/edges;
    # make it strictly increasing for interpolation.
    order = np.argsort(rec)
    rec_s = rec[order]
    prec_s = prec[order]

    # Deduplicate recall values (interp needs monotonic x)
    # Keep the max precision for duplicate recall points (conservative-ish).
    rec_u, idx = np.unique(rec_s, return_index=True)
    prec_u = prec_s[idx]

    # Interpolate precision at recall_grid; outside range -> boundary values
    # (recall is in [0,1])
    p_on_grid = np.interp(recall_grid, rec_u, prec_u, left=prec_u[0], right=prec_u[-1])
    return p_on_grid


def bootstrap_pr_bands(
    y_true,
    y_score,
    *,
    recall_grid=None,
    n_boot=500,
    seed=0,
    ci=(2.5, 97.5),
    min_pos=2,
):
    """
    Bootstrap PR curves by resampling samples with replacement.
    Returns dict with:
      - recall_grid
      - prec_med, prec_lo, prec_hi (arrays)
      - ap, ap_lo, ap_hi (scalars) for PR-AUC (Average Precision)
      - n_boot_ok
    """
    y_true = np.asarray(y_true, int)
    y_score = np.asarray(y_score, float)

    if recall_grid is None:
        recall_grid = np.linspace(0, 1, 201)

    rng = np.random.default_rng(int(seed))
    n = y_true.size

    curves = []
    aps = []
    n_ok = 0

    for _ in range(int(n_boot)):
        idx = rng.integers(0, n, size=n)
        yb = y_true[idx]
        if (yb.sum() < int(min_pos)) or (len(np.unique(yb)) < 2):
            continue
        sb = y_score[idx]

        curves.append(_pr_curve_interp(yb, sb, recall_grid))
        aps.append(float(average_precision_score(yb, sb)))
        n_ok += 1

    if n_ok == 0:
        return {
            "recall_grid": recall_grid,
            "prec_med": np.full_like(recall_grid, np.nan, dtype=float),
            "prec_lo": np.full_like(recall_grid, np.nan, dtype=float),
            "prec_hi": np.full_like(recall_grid, np.nan, dtype=float),
            "ap": np.nan,
            "ap_lo": np.nan,
            "ap_hi": np.nan,
            "n_boot_ok": 0,
        }

    curves = np.vstack(curves)  # (n_ok, len(grid))
    aps = np.asarray(aps, float)

    lo_q, hi_q = np.percentile(curves, ci, axis=0)
    med = np.median(curves, axis=0)

    ap_lo, ap_hi = np.percentile(aps, ci)
    ap = float(average_precision_score(y_true, y_score))

    return {
        "recall_grid": recall_grid,
        "prec_med": med,
        "prec_lo": lo_q,
        "prec_hi": hi_q,
        "ap": float(ap),
        "ap_lo": float(ap_lo),
        "ap_hi": float(ap_hi),
        "n_boot_ok": int(n_ok),
    }

## PR plot with bootstrap bands (single axis; one or multiple signatures)
def plot_pr_curves_signatures_bootstrap(
    adata,
    signatures: dict,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    use_proba=True,
    title=None,
    figsize=(6.5, 5.0),
    legend_outside=True,
    n_boot=500,
    seed=0,
    ci=(2.5, 97.5),
    recall_grid=None,
    alpha_band=0.18,
    ax=None,
):
    """
    signatures: dict name -> dict with keys:
        - "weights_df" (pd.DataFrame gene,beta)
        - "intercept" (float) optional if use_proba=True
    Returns: (df_summary, ax)
    """
    ok_y, y_full = _y_from_adata(adata, label_col, positive_label)
    y_labeled = y_full[ok_y]
    if len(np.unique(y_labeled)) < 2:
        raise ValueError("Need 2 classes to plot PR curve.")

    prevalence = float(y_labeled.mean())
    if recall_grid is None:
        recall_grid = np.linspace(0, 1, 201)

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()

    # baseline prevalence line
    ax.plot([0, 1], [prevalence, prevalence], linestyle="--",
            label=f"Baseline (prev={prevalence:.3f})")

    rows = []
    for name, sig in signatures.items():
        w = sig["weights_df"]
        b0 = float(sig.get("intercept", 0.0))

        s_full = _score_from_weights(adata, w, layer=layer)
        ok = ok_y & np.isfinite(s_full)
        yy = y_full[ok]
        ss = s_full[ok]

        if yy.size < 2 or len(np.unique(yy)) < 2:
            rows.append({
                "Signature": name, "pr_auc": np.nan, "pr_auc_ci_low": np.nan, "pr_auc_ci_high": np.nan,
                "n": int(yy.size), "pos": int(yy.sum()), "neg": int((yy == 0).sum()),
                "n_boot_ok": 0
            })
            continue

        pred = _sigmoid(b0 + ss) if use_proba else ss

        # main curve (from full data)
        prec, rec, _ = precision_recall_curve(yy, pred)
        ap = float(average_precision_score(yy, pred))

        # bootstrap bands
        boot = bootstrap_pr_bands(
            yy, pred,
            recall_grid=recall_grid,
            n_boot=int(n_boot),
            seed=int(seed),
            ci=ci,
            min_pos=2,
        )

        # plot main curve + CI band
        ax.plot(rec, prec, linewidth=2.0,
                label=f"{name} (AP={ap:.3f} [{boot['ap_lo']:.3f},{boot['ap_hi']:.3f}])")
        ax.fill_between(
            boot["recall_grid"], boot["prec_lo"], boot["prec_hi"],
            alpha=float(alpha_band),
            linewidth=0,
        )

        rows.append({
            "Signature": name,
            "pr_auc": ap,
            "pr_auc_ci_low": boot["ap_lo"],
            "pr_auc_ci_high": boot["ap_hi"],
            "n": int(len(yy)),
            "pos": int(yy.sum()),
            "neg": int((yy == 0).sum()),
            "n_boot_ok": int(boot["n_boot_ok"]),
        })

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision (PPV)")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.02)

    if title is not None:
        ax.set_title(title)

    if legend_outside:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    else:
        ax.legend(loc="lower left")

    ax.figure.tight_layout()

    df = pd.DataFrame(rows).sort_values("pr_auc", ascending=False).reset_index(drop=True)
    return df, ax


def plot_pub_template_pr(
    *,
    adata,
    signatures,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    use_proba=True,
    lead_signature=None,   # name in signatures
    title=None,
    figsize=(12.5, 6.8),   # taller by default to fit legends cleanly
    n_boot=1000,
    seed=0,
    ci=(2.5, 97.5),
    legend_fontsize=9,
    legend_max_cols=2,
    legend_row_height=0.40,   # fraction of plot-row height (not figure fraction)
    legend_loc="upper left",
    legend_handlelength=2.5,
    legend_columnspacing=1.2,
    axC_ylabel_rotation=0,     # 0 keeps horizontal labels; 30/45 also ok
    axC_ylabel_fontsize=9,
    axC_ylabel_pad=6,
):
    """
    Creates a publication-style multi-panel figure for one cohort.

    Layout:
      Row 0: A, B, C plots
      Row 1: legend for A, legend for B, empty

    Legends are NEVER drawn inside A/B axes; they are placed below and wrapped
    within their own legend axes.

    Returns df_summary used for bars/table (dfB).
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if not isinstance(signatures, dict) or len(signatures) == 0:
        raise ValueError("signatures must be a non-empty dict")

    if lead_signature is None:
        lead_signature = list(signatures.keys())[0]
    if lead_signature not in signatures:
        raise KeyError(f"lead_signature='{lead_signature}' not found in signatures")

    # ----------------------
    # Grid: 2 rows x 3 cols
    # ----------------------
    # IMPORTANT: legend_row_height is a *ratio*, so legends can grow without overlapping plots.
    # e.g. legend_row_height=0.40 means legends get 0.40 units for every 1.0 unit of plot height.
    hr0 = 1.0
    hr1 = float(legend_row_height)
    if hr1 <= 0:
        hr1 = 0.25

    fig = plt.figure(figsize=figsize)

    gs = gridspec.GridSpec(
        2, 3,
        height_ratios=[hr0, hr1],
        width_ratios=[1.15, 1.15, 0.9],
        hspace=0.22,  # enough vertical gap so x-labels don't collide with legend row
        wspace=0.35,
    )

    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2])

    axA_leg = fig.add_subplot(gs[1, 0])
    axB_leg = fig.add_subplot(gs[1, 1])
    ax_empty = fig.add_subplot(gs[1, 2])

    for ax in (axA_leg, axB_leg, ax_empty):
        ax.axis("off")

    # ----------------------
    # Panel A: lead signature
    # ----------------------
    dfA, _ = plot_pr_curves_signatures_bootstrap(
        adata,
        signatures={lead_signature: signatures[lead_signature]},
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        use_proba=use_proba,
        title="A  Lead signature",
        figsize=None,
        legend_outside=False,   # legend handled here (we will remove it)
        n_boot=n_boot,
        seed=seed,
        ci=ci,
        ax=axA,
    )

    # Remove any legend inside panel A
    legA = axA.get_legend()
    if legA is not None:
        legA.remove()

    # Collect handles/labels AFTER drawing
    hA, lA = axA.get_legend_handles_labels()

    # spacing + consistent axes formatting
    axA.set_xlabel("Recall", labelpad=10)
    axA.set_ylabel("Precision (PPV)")
    axA.set_xlim(0, 1.0)
    axA.set_ylim(0, 1.02)

    # Legend for A in its own axis (full width)
    if hA and lA:
        ncolA = min(int(legend_max_cols), max(1, len(lA)))
        axA_leg.legend(
            hA, lA,
            loc=legend_loc,
            bbox_to_anchor=(0.0, 0.0, 1.0, 1.0),  # occupy full legend-axis box
            mode="expand",
            frameon=False,
            fontsize=legend_fontsize,
            ncol=ncolA,
            handlelength=float(legend_handlelength),
            columnspacing=float(legend_columnspacing),
            borderaxespad=0.0,
        )

    # ----------------------
    # Panel B: comparison
    # ----------------------
    dfB, _ = plot_pr_curves_signatures_bootstrap(
        adata,
        signatures=signatures,
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        use_proba=use_proba,
        title="B  Signature comparison",
        figsize=None,
        legend_outside=False,   # legend handled here (we will remove it)
        n_boot=n_boot,
        seed=seed,
        ci=ci,
        ax=axB,
    )

    legB = axB.get_legend()
    if legB is not None:
        legB.remove()

    hB, lB = axB.get_legend_handles_labels()

    axB.set_xlabel("Recall", labelpad=10)
    axB.set_ylabel("Precision (PPV)")
    axB.set_xlim(0, 1.0)
    axB.set_ylim(0, 1.02)

    # Legend for B in its own axis
    # Auto choose more columns if many signatures to reduce overlap height
    if hB and lB:
        # heuristic: more items -> more columns, capped by legend_max_cols
        # (you can increase legend_max_cols to 3–4 for many signatures)
        n_items = len(lB)
        ncolB = min(int(legend_max_cols), max(1, n_items))
        axB_leg.legend(
            hB, lB,
            loc=legend_loc,
            bbox_to_anchor=(0.0, 0.0, 1.0, 1.0),
            mode="expand",
            frameon=False,
            fontsize=legend_fontsize,
            ncol=ncolB,
            handlelength=float(legend_handlelength),
            columnspacing=float(legend_columnspacing),
            borderaxespad=0.0,
        )

    # ----------------------
    # Panel C: PR-AUC bars
    # ----------------------
    df_plot = dfB.copy()
    df_plot = df_plot[np.isfinite(df_plot["pr_auc"])].copy()
    df_plot = df_plot.sort_values("pr_auc", ascending=True)  # bottom->top

    y = np.arange(df_plot.shape[0])
    ap = df_plot["pr_auc"].to_numpy(float)

    # errors if present
    xerr = None
    if ("pr_auc_ci_low" in df_plot.columns) and ("pr_auc_ci_high" in df_plot.columns):
        lo = df_plot["pr_auc_ci_low"].to_numpy(float)
        hi = df_plot["pr_auc_ci_high"].to_numpy(float)
        if np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)):
            xerr = np.vstack([ap - lo, hi - ap])

    axC.barh(y, ap, xerr=xerr, capsize=3)
    axC.set_yticks(y)

    # KEY FIX: rotate labels + anchor them so they stay within panel C
    axC.set_yticklabels(
        df_plot["Signature"].tolist(),
        fontsize=axC_ylabel_fontsize,
        rotation=axC_ylabel_rotation,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    axC.tick_params(axis="y", pad=float(axC_ylabel_pad))

    axC.set_xlabel("PR-AUC (Average Precision)")
    axC.set_title("C  PR-AUC ± bootstrap CI")
    axC.set_xlim(0, 1.0)

    # add small internal margin so labels don't collide with left border
    axC.margins(y=0.08)

    if title is not None:
        fig.suptitle(title, y=1.02, fontsize=12)

    # Reserve enough bottom area; legends live there, but we still avoid clipping
    # (tight_layout can sometimes clip legends, so we give it a rect)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))

    return dfB




def signature_dispersion_heatmap(
    adata: AnnData,
    *,
    uns_key: str = "signature_dispersion",
    title: Optional[str] = None,
    show: bool = True,
    save: Optional[str] = None,
):
    """
    Heatmap of signature dispersion across groups.
    Expects `adata.uns[uns_key]` to be a DataFrame with columns:
      group, signature, dispersion
    (as returned by `bk.tl.signature_dispersion_by_group`).
    """
    if uns_key not in adata.uns:
        raise KeyError(f"{uns_key!r} not found in adata.uns")

    df = adata.uns[uns_key]
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    req = {"group", "signature", "dispersion"}
    if not req.issubset(df.columns):
        raise ValueError(f"{uns_key} must contain columns {req}")

    mat = df.pivot(index="signature", columns="group", values="dispersion")
    arr = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(arr, aspect="auto")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index.tolist(), fontsize=7)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns.tolist(), fontsize=7, rotation=90)
    ax.set_title(title or "Signature dispersion by group")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return ax