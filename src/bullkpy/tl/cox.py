from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Sequence, Optional, Tuple, Union

import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData

try:
    from statsmodels.duration.hazard_regression import PHReg
except Exception:  # pragma: no cover
    PHReg = None

from ..logging import info, warn



def _get_matrix(adata, layer, use="samples"):
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    return X

def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
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

def _get_expr_matrix(
    adata: ad.AnnData,
    genes: list[str],
    layer: str | None = "log1p_cpm",
) -> np.ndarray:
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    idx = adata.var_names.get_indexer(genes)
    if np.any(idx < 0):
        missing = [g for g, i in zip(genes, idx) if i < 0]
        raise KeyError(f"Genes not found in adata.var_names (first 10): {missing[:10]}")
    M = X[:, idx]
    if hasattr(M, "toarray"):
        M = M.toarray()
    else:
        M = np.asarray(M, dtype=float)
    return np.asarray(M, dtype=float)


def cox_gene_association(
    adata: ad.AnnData,
    *,
    time_col: str,
    event_col: str,
    genes: Sequence[str] | None = None,
    layer: str | None = "log1p_cpm",
    covariates: Sequence[str] | None = None,
    standardize: bool = True,
    min_events: int = 5,
    min_var: float = 1e-12,
    penalizer: float = 0.1,
    robust: bool = True,
    adjust: Literal["bh", "none"] = "bh",
    key_added: str | None = None,
) -> pd.DataFrame:
    """
    Univariate Cox association for many genes (bulk survival screen).

    For each gene g:
        CoxPH: time_col ~ g (+ covariates)

    Returns columns:
        gene, coef, HR, CI_lower, CI_upper, se, z, pval, qval, n, n_events
    """
    # ---- dependencies ----
    try:
        from lifelines import CoxPHFitter  # type: ignore
    except Exception as e:
        raise ImportError(
            "cox_gene_association requires `lifelines`.\n"
            "Install with: pip install lifelines\n"
            f"Original error: {e}"
        )

    if time_col not in adata.obs.columns:
        raise KeyError(f"time_col='{time_col}' not found in adata.obs")
    if event_col not in adata.obs.columns:
        raise KeyError(f"event_col='{event_col}' not found in adata.obs")

    # ---- genes ----
    if genes is None:
        genes_use = list(map(str, adata.var_names))
    else:
        genes_use = [str(g) for g in genes]

    # ---- covariates ----
    covs = [str(c) for c in covariates] if covariates is not None else []
    for c in covs:
        if c not in adata.obs.columns:
            raise KeyError(f"covariate '{c}' not found in adata.obs")

    # ---- base clinical df ----
    df_base = adata.obs[[time_col, event_col] + covs].copy()
    df_base[time_col] = pd.to_numeric(df_base[time_col], errors="coerce")
    df_base[event_col] = pd.to_numeric(df_base[event_col], errors="coerce")
    df_base = df_base.dropna(subset=[time_col, event_col])

    # sanity: event coding to 0/1
    ev = df_base[event_col].to_numpy()
    ev = np.where(ev > 0, 1, 0).astype(int)
    df_base[event_col] = ev

    n_events_total = int(ev.sum())
    if n_events_total < int(min_events):
        raise ValueError(
            f"Not enough events for Cox scan: n_events={n_events_total} < min_events={min_events}"
        )

    # coerce covariates to numeric if possible
    for c in covs:
        if not np.issubdtype(df_base[c].dtype, np.number):
            df_base[c] = pd.to_numeric(df_base[c], errors="coerce")

    # ---- expression matrix aligned to df_base ----
    genes_present = [g for g in genes_use if g in adata.var_names]
    if len(genes_present) == 0:
        raise ValueError("No requested genes are present in adata.var_names.")
    if len(genes_present) < len(genes_use):
        miss = [g for g in genes_use if g not in adata.var_names]
        warn(f"cox_gene_association: skipping missing genes (first 10): {miss[:10]}")

    M = _get_expr_matrix(adata, genes_present, layer=layer)  # (n_obs x n_genes)
    keep_idx = adata.obs.index.get_indexer(df_base.index)
    M = M[keep_idx, :]

    # ---- standardize per gene ----
    if standardize:
        mu = np.nanmean(M, axis=0, keepdims=True)
        sd = np.nanstd(M, axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        M = (M - mu) / sd

    info(
        f"Running Cox gene association: genes={M.shape[1]}, n={M.shape[0]}, "
        f"events={int(df_base[event_col].sum())}, covariates={covs}, penalizer={penalizer}"
    )

    out = []
    z_ci = 1.96  # 95% CI

    # counters for debugging
    n_try = 0
    n_ok = 0
    n_lowvar = 0
    n_fail = 0
    n_lowvar_event = 0
    n_too_few_events = 0

    # precompute event mask (based on df_base, before NA drop on covariates)
    base_events = df_base[event_col].astype(int).to_numpy().astype(bool)

    for j, g in enumerate(genes_present):
        x = M[:, j]
        if not np.isfinite(x).any():
            n_fail += 1
            continue

        # overall variance filter
        if float(np.nanvar(x)) <= float(min_var):
            n_lowvar += 1
            continue

        # event-stratified variance filter (prevents complete separation issues)
        v1 = float(np.nanvar(x[base_events])) if base_events.any() else 0.0
        v0 = float(np.nanvar(x[~base_events])) if (~base_events).any() else 0.0
        if (v1 <= float(min_var)) or (v0 <= float(min_var)):
            n_lowvar_event += 1
            continue

        df = df_base.copy()
        df["_gene_"] = x
        df = df.dropna()  # drops NA covariates too

        n = int(df.shape[0])
        n_events = int(df[event_col].sum())
        if n_events < int(min_events):
            n_too_few_events += 1
            continue

        n_try += 1

        try:
            # instantiate per-gene for safety
            cph = CoxPHFitter(penalizer=float(penalizer))

            cph.fit(
                df,
                duration_col=time_col,
                event_col=event_col,
                robust=bool(robust),
                show_progress=False,
            )

            summ = cph.summary.loc["_gene_"]

            coef = float(summ["coef"])
            se = float(summ["se(coef)"]) if "se(coef)" in summ else float(summ.get("se", np.nan))
            z = float(summ["z"]) if "z" in summ else np.nan
            p = float(summ["p"]) if "p" in summ else np.nan

            # HR + CI
            hr = float(np.exp(coef))
            ci_lower = float(np.exp(coef - z_ci * se)) if np.isfinite(se) else np.nan
            ci_upper = float(np.exp(coef + z_ci * se)) if np.isfinite(se) else np.nan

            # guard against inf/nan explosions
            if not (np.isfinite(hr) and np.isfinite(p)):
                n_fail += 1
                continue

            out.append(
                {
                    "gene": g,
                    "coef": coef,
                    "HR": hr,
                    "CI_lower": ci_lower,
                    "CI_upper": ci_upper,
                    "se": se,
                    "z": z,
                    "pval": p,
                    "n": n,
                    "n_events": n_events,
                }
            )
            n_ok += 1

        except Exception:
            n_fail += 1
            continue

    res = pd.DataFrame(out)
    info(
        "Cox summary: "
        f"attempted={n_try}, ok={n_ok}, lowvar={n_lowvar}, lowvar_event={n_lowvar_event}, "
        f"too_few_events={n_too_few_events}, failed={n_fail}"
    )

    if res.shape[0] == 0:
        raise ValueError(
            "No genes produced valid Cox fits (after filtering). "
            "Try increasing penalizer (e.g. 0.5 or 1.0), lowering min_var, or relaxing min_events."
        )

    # adjust p-values
    if adjust == "bh":
        res["qval"] = _bh_fdr(res["pval"].to_numpy(dtype=float))
    else:
        res["qval"] = np.nan

    res = res.sort_values(["qval", "pval"], ascending=[True, True], na_position="last").reset_index(drop=True)

    # store
    if key_added is not None:
        adata.uns.setdefault("cox_gene_association", {})
        adata.uns["cox_gene_association"][key_added] = {
            "params": {
                "time_col": time_col,
                "event_col": event_col,
                "layer": layer,
                "genes_tested": int(len(genes_present)),
                "standardize": bool(standardize),
                "covariates": covs,
                "min_events": int(min_events),
                "min_var": float(min_var),
                "penalizer": float(penalizer),
                "robust": bool(robust),
                "adjust": adjust,
            },
            "results": res,
        }

    return res


def cox_multivariate(
    adata: ad.AnnData,
    *,
    genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    penalizer: float = 0.1,
    l1_ratio: float = 0.0,
    robust: bool = True,
    min_events: int = 5,
    min_var: float = 1e-12,
    key_added: str | None = None,
    backend: Literal["lifelines", "statsmodels"] = "lifelines",
    tie_method: Literal["breslow", "efron"] = "breslow",  # used only for statsmodels
) -> pd.DataFrame:
    """
    Multivariate Cox PH on a gene panel.

    Default backend="lifelines" (recommended): supports penalization via `penalizer` and `l1_ratio`.

    Returns
    -------
    pd.DataFrame
        Columns: gene, coef, HR, CI_lower, CI_upper, se, z, pval
    """
    # ---- input checks ----
    genes = [str(g) for g in genes if str(g) in adata.var_names]
    if len(genes) < 1:
        raise ValueError("No provided genes found in adata.var_names.")

    if time_col not in adata.obs.columns:
        raise KeyError(f"time_col '{time_col}' not found in adata.obs")
    if event_col not in adata.obs.columns:
        raise KeyError(f"event_col '{event_col}' not found in adata.obs")

    # expression (samples x genes)
    X = _get_matrix(adata, layer=layer, use="samples")
    gidx = adata.var_names.get_indexer(genes)
    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # clinical
    t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(dtype=float)
    e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(dtype=float)
    e = np.where(np.isfinite(e) & (e > 0), 1, 0).astype(int)

    # base filtering (finite)
    ok = np.isfinite(t) & np.isfinite(e) & np.all(np.isfinite(M), axis=1)
    if int(ok.sum()) < 5:
        raise ValueError(f"Not enough valid samples for Cox fit: n={int(ok.sum())}")

    # ensure enough events
    if int(e[ok].sum()) < int(min_events):
        raise ValueError(f"Not enough events for Cox fit: events={int(e[ok].sum())} < {min_events}")

    M = M[ok, :]
    t = t[ok]
    e = e[ok]
    n = int(len(t))
    n_events = int(e.sum())

    # drop low-variance genes (helps convergence)
    var = np.nanvar(M, axis=0)
    keep = var > float(min_var)
    if keep.sum() < 1:
        raise ValueError("All genes filtered out by min_var.")
    if keep.sum() < len(genes):
        genes = [g for g, k in zip(genes, keep) if k]
        M = M[:, keep]

    # standardize predictors
    if standardize:
        mu = np.nanmean(M, axis=0, keepdims=True)
        sd = np.nanstd(M, axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        M = (M - mu) / sd

    info(f"Fitting multivariate Cox: genes={len(genes)}, n={n}, events={n_events}, backend={backend}")

    # ---- backend: lifelines (recommended) ----
    if backend == "lifelines":
        try:
            from lifelines import CoxPHFitter  # type: ignore
        except Exception as ex:
            raise ImportError(
                "cox_multivariate(backend='lifelines') requires `lifelines`.\n"
                "Install with: pip install lifelines\n"
                f"Original error: {ex}"
            )

        df = pd.DataFrame(M, columns=[f"g__{g}" for g in genes])
        df[time_col] = t
        df[event_col] = e

        cph = CoxPHFitter(penalizer=float(penalizer), l1_ratio=float(l1_ratio))
        cph.fit(
            df,
            duration_col=time_col,
            event_col=event_col,
            robust=bool(robust),
            show_progress=False,
        )

        s = cph.summary.copy()
        # lifelines uses: coef, se(coef), z, p, exp(coef), exp(coef) lower 95%, exp(coef) upper 95%
        out = pd.DataFrame(
            {
                "gene": [c.replace("g__", "", 1) for c in s.index.astype(str)],
                "coef": s["coef"].to_numpy(dtype=float),
                "se": s["se(coef)"].to_numpy(dtype=float),
                "z": s["z"].to_numpy(dtype=float),
                "pval": s["p"].to_numpy(dtype=float),
                "HR": s["exp(coef)"].to_numpy(dtype=float),
                "CI_lower": s["exp(coef) lower 95%"].to_numpy(dtype=float),
                "CI_upper": s["exp(coef) upper 95%"].to_numpy(dtype=float),
            }
        ).sort_values("pval", ascending=True).reset_index(drop=True)

    # ---- backend: statsmodels PHReg (no penalizer) ----
    else:
        if PHReg is None:
            raise ImportError("cox_multivariate(backend='statsmodels') requires statsmodels with PHReg available.")

        model = PHReg(endog=t, exog=M, status=e.astype(int), ties=tie_method)
        res = model.fit(disp=0)

        beta = np.asarray(res.params, dtype=float)
        se = np.asarray(res.bse, dtype=float)
        z = beta / se
        pval = 2.0 * (1.0 - _norm_cdf(np.abs(z)))

        lo = beta - 1.96 * se
        hi = beta + 1.96 * se

        out = pd.DataFrame(
            {
                "gene": genes,
                "coef": beta,
                "se": se,
                "z": z,
                "pval": pval,
                "HR": np.exp(beta),
                "CI_lower": np.exp(lo),
                "CI_upper": np.exp(hi),
            }
        ).sort_values("pval", ascending=True).reset_index(drop=True)

    # ---- store results (optional) ----
    if key_added is not None:
        adata.uns.setdefault("cox_multivariate", {})
        adata.uns["cox_multivariate"][str(key_added)] = {
            "params": {
                "time_col": time_col,
                "event_col": event_col,
                "layer": layer,
                "genes": list(map(str, genes)),
                "standardize": bool(standardize),
                "backend": backend,
                "penalizer": float(penalizer),
                "l1_ratio": float(l1_ratio),
                "robust": bool(robust),
                "min_events": int(min_events),
                "min_var": float(min_var),
            },
            "results": out,
        }

    return out

def _norm_cdf(x: np.ndarray) -> np.ndarray:
    from scipy.stats import norm
    return norm.cdf(x)



def _require_lifelines():
    try:
        from lifelines import CoxPHFitter  # type: ignore
        from lifelines.utils import concordance_index  # type: ignore
    except Exception as e:
        raise ImportError(
            "This function requires `lifelines`.\n"
            "Install with: pip install lifelines\n"
            f"Original error: {e}"
        )
    return CoxPHFitter, concordance_index


def _as_list(x):
    if x is None:
        return None
    if isinstance(x, (list, tuple, np.ndarray, pd.Index)):
        return list(x)
    return [x]


def _prepare_cox_df(
    adata: ad.AnnData,
    genes: Sequence[str],
    *,
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    min_var: float = 1e-12,
) -> tuple[pd.DataFrame, list[str]]:
    if time_col not in adata.obs.columns:
        raise KeyError(f"time_col='{time_col}' not found in adata.obs")
    if event_col not in adata.obs.columns:
        raise KeyError(f"event_col='{event_col}' not found in adata.obs")

    genes = [str(g) for g in genes]
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) == 0:
        raise ValueError("None of the requested genes are present in adata.var_names.")
    if len(genes_present) < len(genes):
        miss = [g for g in genes if g not in adata.var_names]
        warn(f"cox: skipping missing genes (first 10): {miss[:10]}")

    # expression matrix (samples x genes_present)
    X = _get_matrix(adata, layer=layer, use="samples")  # samples x all_genes
    gidx = adata.var_names.get_indexer(genes_present)
    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # drop low variance genes (pre-standardization)
    var = np.nanvar(M, axis=0)
    keepg = var > float(min_var)
    genes_kept = [g for g, k in zip(genes_present, keepg) if bool(k)]
    if len(genes_kept) == 0:
        raise ValueError("All requested genes were filtered by min_var.")
    M = M[:, keepg]

    # standardize per gene
    if standardize:
        mu = np.nanmean(M, axis=0, keepdims=True)
        sd = np.nanstd(M, axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        M = (M - mu) / sd

    t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(dtype=float)
    e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(dtype=float)
    e = np.where(e > 0, 1, 0).astype(int)

    df = pd.DataFrame(M, index=adata.obs_names, columns=genes_kept)
    df[time_col] = t
    df[event_col] = e

    # drop rows with missing survival or any gene NA
    df = df.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")

    if df.shape[0] < 5:
        raise ValueError(f"Not enough valid samples after filtering: n={df.shape[0]}")

    return df, genes_kept


def _cox_fit_and_risk(
    df: pd.DataFrame,
    *,
    duration_col: str,
    event_col: str,
    penalizer: float,
    l1_ratio: float,
    robust: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    CoxPHFitter, _ = _require_lifelines()
    cph = CoxPHFitter(penalizer=float(penalizer), l1_ratio=float(l1_ratio))
    cph.fit(df, duration_col=duration_col, event_col=event_col, robust=bool(robust), show_progress=False)

    # lifelines stores params as a Series indexed by covariate names
    beta = cph.params_.to_numpy(dtype=float)
    X = df.drop(columns=[duration_col, event_col]).to_numpy(dtype=float)
    risk = X @ beta
    return beta, risk


def _c_index_from_risk(
    times: np.ndarray, events: np.ndarray, risk: np.ndarray
) -> float:
    _, concordance_index = _require_lifelines()
    # higher risk -> shorter time, so use -risk to align with concordance_index definition
    return float(concordance_index(times, -risk, event_observed=events))


def cox_penalizer_cv(
    adata: ad.AnnData,
    *,
    genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    penalizers: Sequence[float] = tuple(np.logspace(-2, 2, 20)),  # stronger defaults
    l1_ratios: Sequence[float] = (0.0, 0.2, 0.5),                 # include ridge for stability
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 0,
    min_events: int = 10,
    min_var: float = 1e-12,
    robust: bool = False,
    return_fits: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """
    Cross-validated selection of Cox penalizer (lambda) and optional l1_ratio.

    Uses lifelines CoxPHFitter(penalizer=..., l1_ratio=...).
    Scores on held-out folds via Harrell's C-index.

    Returns
    -------
    df_cv : pd.DataFrame
        Columns: penalizer, l1_ratio, c_index_mean, c_index_sd, n_selected_median, n_fails
    best_params : dict
        {"penalizer": float, "l1_ratio": float}
    """
    import warnings
    import numpy as np
    import pandas as pd

    CoxPHFitter, _ = _require_lifelines()

    # Build full df once (NO standardization here; do it per-fold)
    df_all, genes_kept = _prepare_cox_df(
        adata,
        genes,
        time_col=time_col,
        event_col=event_col,
        layer=layer,
        standardize=False,   # IMPORTANT: fold-wise standardization below
        min_var=min_var,
    )

    events_total = int(df_all[event_col].sum())
    if events_total < int(min_events):
        raise ValueError(f"Not enough events: n_events={events_total} < min_events={min_events}")

    Xcols = [c for c in df_all.columns if c not in (time_col, event_col)]
    n = df_all.shape[0]

    penalizers = [float(x) for x in penalizers]
    l1_ratios = [float(x) for x in l1_ratios]

    info(
        f"cox_penalizer_cv: genes={len(Xcols)}, n={n}, events={events_total}, "
        f"splits={n_splits} repeats={n_repeats}, grid={len(penalizers)*len(l1_ratios)}"
    )

    # ---- event-stratified folds ----
    def _stratified_folds(event: np.ndarray, k: int, rng: np.random.Generator):
        event = np.asarray(event).astype(int)
        idx1 = np.where(event == 1)[0]
        idx0 = np.where(event == 0)[0]
        rng.shuffle(idx1)
        rng.shuffle(idx0)
        f1 = np.array_split(idx1, k)
        f0 = np.array_split(idx0, k)
        folds = [np.concatenate([f1[i], f0[i]]) for i in range(k)]
        return folds

    rows = []
    fits = {}  # optional: (pen, l1)->list of fitted models

    for pen in penalizers:
        for l1 in l1_ratios:
            scores: list[float] = []
            nsel: list[int] = []
            nfail = 0

            for rep in range(int(n_repeats)):
                rng = np.random.default_rng(int(seed) + 1000 * rep)
                folds = _stratified_folds(df_all[event_col].to_numpy(dtype=int), int(n_splits), rng)

                for k in range(int(n_splits)):
                    te = folds[k]
                    tr = np.concatenate([folds[i] for i in range(int(n_splits)) if i != k])

                    df_tr = df_all.iloc[tr, :].copy()
                    df_te = df_all.iloc[te, :].copy()

                    if int(df_tr[event_col].sum()) < int(min_events):
                        nfail += 1
                        continue

                    # ---- fold-wise standardization (train -> apply to test) ----
                    if standardize:
                        Xtr = df_tr[Xcols].to_numpy(dtype=float)
                        mu = np.nanmean(Xtr, axis=0, keepdims=True)
                        sd = np.nanstd(Xtr, axis=0, keepdims=True)
                        sd[sd == 0] = 1.0
                        df_tr.loc[:, Xcols] = (Xtr - mu) / sd

                        Xte = df_te[Xcols].to_numpy(dtype=float)
                        df_te.loc[:, Xcols] = (Xte - mu) / sd

                    try:
                        cph = CoxPHFitter(penalizer=float(pen), l1_ratio=float(l1))
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")  # silence ConvergenceWarning/RuntimeWarning spam
                            cph.fit(
                                df_tr,
                                duration_col=time_col,
                                event_col=event_col,
                                robust=bool(robust),
                                show_progress=False,
                            )

                        beta = cph.params_.reindex(Xcols).to_numpy(dtype=float)
                        n_selected = int(np.sum(np.abs(beta) > 1e-10))
                        nsel.append(n_selected)

                        Xte = df_te[Xcols].to_numpy(dtype=float)
                        risk_te = Xte @ beta

                        # optional numeric guard (rarely needed, but helps with extreme cases)
                        risk_te = np.clip(risk_te, -50, 50)

                        ci = _c_index_from_risk(
                            df_te[time_col].to_numpy(dtype=float),
                            df_te[event_col].to_numpy(dtype=int),
                            risk_te,
                        )
                        if np.isfinite(ci):
                            scores.append(float(ci))
                        else:
                            nfail += 1

                        if return_fits:
                            fits.setdefault((pen, l1), []).append(cph)

                    except Exception:
                        nfail += 1
                        continue

            rows.append(
                {
                    "penalizer": pen,
                    "l1_ratio": l1,
                    "c_index_mean": float(np.nanmean(scores)) if len(scores) else np.nan,
                    "c_index_sd": float(np.nanstd(scores)) if len(scores) else np.nan,
                    "n_selected_median": float(np.nanmedian(nsel)) if len(nsel) else np.nan,
                    "n_fails": int(nfail),
                    "n_scores": int(len(scores)),
                }
            )

    df_cv = pd.DataFrame(rows)

    # Prefer higher C-index, then smaller median panel, then fewer failures
    df_cv = df_cv.sort_values(
        ["c_index_mean", "n_selected_median", "n_fails"],
        ascending=[False, True, True],
        na_position="last",
    ).reset_index(drop=True)

    df_ok = df_cv[np.isfinite(df_cv["c_index_mean"])].copy()
    if df_ok.shape[0] == 0:
        raise ValueError(
            "cox_penalizer_cv: no valid CV fits. "
            "Try fewer genes (e.g. 100–200), stronger penalizer grid (>=0.1), "
            "and/or filter low-variance / separation-like genes."
        )

    best = df_ok.iloc[0]
    best_params = {"penalizer": float(best["penalizer"]), "l1_ratio": float(best["l1_ratio"])}

    return df_cv, best_params


def cox_stability_selection(
    adata: ad.AnnData,
    *,
    genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    penalizer: float = 0.1,
    l1_ratio: float = 0.8,
    n_boot: int = 200,
    sample_fraction: float = 0.8,
    seed: int = 0,
    min_events: int = 10,
    min_var: float = 1e-12,
    robust: bool = False,
) -> pd.DataFrame:
    """
    Stability selection for penalized Cox: bootstrap/subsample, refit, compute selection frequency.

    Selected if |beta| > 1e-10.

    Returns
    -------
    df_stab : pd.DataFrame
        Columns:
          gene, freq, mean_coef, sign_consistency, n_runs, n_selected_mean, n_selected_median
    """
    _require_lifelines()
    df_all, genes_kept = _prepare_cox_df(
        adata, genes,
        time_col=time_col, event_col=event_col, layer=layer,
        standardize=standardize, min_var=min_var
    )

    n = df_all.shape[0]
    rng = np.random.default_rng(int(seed))

    events_total = int(df_all[event_col].sum())
    if events_total < int(min_events):
        raise ValueError(f"Not enough events: n_events={events_total} < min_events={min_events}")

    Xcols = [c for c in df_all.columns if c not in (time_col, event_col)]
    p = len(Xcols)
    if p < 1:
        raise ValueError("cox_stability_selection: no gene columns found after preparation.")

    # global index for accumulation (fixed across boots)
    col_to_j = {c: j for j, c in enumerate(Xcols)}

    sel_counts = np.zeros(p, dtype=int)      # times selected
    coef_sum = np.zeros(p, dtype=float)      # sum of betas when selected
    sign_pos = np.zeros(p, dtype=int)        # selected runs with beta>0
    sign_neg = np.zeros(p, dtype=int)        # selected runs with beta<0
    n_selected_each = []
    n_ok = 0
    n_fail = 0
    n_lowvar_drop_total = 0

    info(
        f"cox_stability_selection: genes={p}, n={n}, events={events_total}, "
        f"penalizer={penalizer}, l1_ratio={l1_ratio}, n_boot={n_boot}, frac={sample_fraction}"
    )

    import warnings

    for _ in range(int(n_boot)):
        m = int(np.ceil(float(sample_fraction) * n))
        idx = rng.choice(np.arange(n), size=m, replace=True)
        df_b = df_all.iloc[idx, :].copy()

        # enough events?
        if int(df_b[event_col].sum()) < int(min_events):
            n_fail += 1
            continue

        # ---- per-bootstrap low-variance filter ----
        # (this is what removes most convergence warnings)
        Xb = df_b[Xcols].to_numpy(dtype=float)
        v = np.nanvar(Xb, axis=0)
        keep_mask = np.isfinite(v) & (v > float(min_var))
        Xcols_eff = [c for c, keep in zip(Xcols, keep_mask) if keep]

        dropped = int(np.sum(~keep_mask))
        if dropped > 0:
            n_lowvar_drop_total += dropped

        # must have at least 1 gene column
        if len(Xcols_eff) < 1:
            n_fail += 1
            continue

        # Fit only with kept columns
        df_fit = df_b[[time_col, event_col] + Xcols_eff].copy()

        try:
            CoxPHFitter, _ = _require_lifelines()
            cph = CoxPHFitter(penalizer=float(penalizer), l1_ratio=float(l1_ratio))

            # Silence lifelines convergence spam (still counted via try/except)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cph.fit(
                    df_fit,
                    duration_col=time_col,
                    event_col=event_col,
                    robust=bool(robust),
                    show_progress=False,
                )

            # IMPORTANT: align betas by name
            # cph.params_ is a pd.Series indexed by covariate name
            params = cph.params_.astype(float)

            chosen = (params.abs() > 1e-10)
            chosen_cols = params.index[chosen].tolist()

            # update accumulators
            for c in params.index.tolist():
                j = col_to_j.get(c, None)
                if j is None:
                    continue
                b = float(params.loc[c])
                if abs(b) > 1e-10:
                    sel_counts[j] += 1
                    coef_sum[j] += b
                    if b > 0:
                        sign_pos[j] += 1
                    elif b < 0:
                        sign_neg[j] += 1

            n_selected_each.append(int(len(chosen_cols)))
            n_ok += 1

        except Exception:
            n_fail += 1
            continue

    if n_ok == 0:
        raise ValueError("cox_stability_selection: no successful fits. Increase penalizer or reduce genes.")

    freq = sel_counts / float(n_ok)
    mean_coef = np.where(sel_counts > 0, coef_sum / np.maximum(sel_counts, 1), 0.0)

    # sign consistency among selected runs
    sign_cons = np.full(p, np.nan, dtype=float)
    for j in range(p):
        if sel_counts[j] == 0:
            continue
        maj = max(sign_pos[j], sign_neg[j])
        sign_cons[j] = float(maj / sel_counts[j])

    df = pd.DataFrame(
        {
            "gene": Xcols,
            "freq": freq,
            "mean_coef": mean_coef,
            "sign_consistency": sign_cons,
            "n_runs": int(n_ok),
            "n_selected_mean": float(np.mean(n_selected_each)) if n_selected_each else np.nan,
            "n_selected_median": float(np.median(n_selected_each)) if n_selected_each else np.nan,
        }
    ).sort_values(["freq", "sign_consistency", "gene"], ascending=[False, False, True]).reset_index(drop=True)

    if n_fail > 0:
        warn(f"cox_stability_selection: {n_fail} bootstrap fits skipped/failed (events too low or convergence).")
    if n_lowvar_drop_total > 0:
        info(f"cox_stability_selection: dropped low-variance columns (total across boots) = {n_lowvar_drop_total}")

    return df

def panel_size_cindex(
    adata: ad.AnnData,
    *,
    ranked_genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    penalizer: float = 0.1,
    l1_ratio: float = 0.8,
    sizes: Sequence[int] = (1, 2, 3, 5, 8, 10, 12, 15, 20),
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 0,
    min_events: int = 10,
    min_var: float = 1e-12,
    robust: bool = False,
) -> pd.DataFrame:
    """
    Panel size vs CV C-index curve for a ranked gene list.

    For each size k: fit penalized multivariate Cox on top-k genes, evaluate C-index on held-out folds.

    Returns
    -------
    df_curve : pd.DataFrame
        Columns: size, c_index_mean, c_index_sd, n_selected_median, n_fails
    """
    _require_lifelines()

    ranked_genes = [str(g) for g in ranked_genes]
    rng = np.random.default_rng(int(seed))

    # Prepare once for the max size genes (reduces work)
    kmax = int(np.max(list(sizes)))
    genes_use = [g for g in ranked_genes if g in adata.var_names][:kmax]
    if len(genes_use) == 0:
        raise ValueError("panel_size_cindex: none of ranked_genes present in adata.var_names.")
    df_all, genes_kept = _prepare_cox_df(
        adata, genes_use, time_col=time_col, event_col=event_col, layer=layer,
        standardize=standardize, min_var=min_var
    )
    n = df_all.shape[0]
    events_total = int(df_all[event_col].sum())
    if events_total < int(min_events):
        raise ValueError(f"Not enough events: n_events={events_total} < min_events={min_events}")

    def _folds(n: int, k: int, rng: np.random.Generator):
        idx = np.arange(n)
        rng.shuffle(idx)
        return np.array_split(idx, k)

    rows = []
    info(
        f"panel_size_cindex: n={n}, events={events_total}, penalizer={penalizer}, l1_ratio={l1_ratio}, "
        f"sizes={list(sizes)}, splits={n_splits} repeats={n_repeats}"
    )

    for k in [int(x) for x in sizes]:
        genes_k = genes_use[:k]
        cols_k = [c for c in genes_k if c in df_all.columns]

        if len(cols_k) < 1:
            rows.append({"size": k, "c_index_mean": np.nan, "c_index_sd": np.nan,
                         "n_selected_median": np.nan, "n_fails": 0})
            continue

        scores = []
        nsel = []
        nfail = 0

        for rep in range(int(n_repeats)):
            folds = _folds(n, int(n_splits), rng)
            for i in range(int(n_splits)):
                te = folds[i]
                tr = np.concatenate([folds[j] for j in range(int(n_splits)) if j != i])

                df_tr = df_all.iloc[tr, :][cols_k + [time_col, event_col]].copy()
                df_te = df_all.iloc[te, :][cols_k + [time_col, event_col]].copy()

                if int(df_tr[event_col].sum()) < int(min_events):
                    nfail += 1
                    continue

                try:
                    CoxPHFitter, _ = _require_lifelines()
                    cph = CoxPHFitter(penalizer=float(penalizer), l1_ratio=float(l1_ratio))
                    cph.fit(df_tr, duration_col=time_col, event_col=event_col, robust=bool(robust), show_progress=False)

                    beta = cph.params_.to_numpy(dtype=float)
                    nsel.append(int(np.sum(np.abs(beta) > 1e-10)))

                    Xte = df_te[cols_k].to_numpy(dtype=float)
                    risk_te = Xte @ beta
                    ci = _c_index_from_risk(
                        df_te[time_col].to_numpy(dtype=float),
                        df_te[event_col].to_numpy(dtype=int),
                        risk_te,
                    )
                    scores.append(ci)
                except Exception:
                    nfail += 1
                    continue

        rows.append(
            {
                "size": int(k),
                "c_index_mean": float(np.nanmean(scores)) if scores else np.nan,
                "c_index_sd": float(np.nanstd(scores)) if scores else np.nan,
                "n_selected_median": float(np.nanmedian(nsel)) if nsel else np.nan,
                "n_fails": int(nfail),
            }
        )

    return pd.DataFrame(rows).sort_values("size").reset_index(drop=True)


def recommended_cox_panel(
    adata: ad.AnnData,
    *,
    genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    # CV grid
    penalizers: Sequence[float] = tuple(np.logspace(-3, 1, 15)),
    l1_ratios: Sequence[float] = (0.5, 0.8, 1.0),
    n_splits: int = 5,
    n_repeats: int = 3,
    seed: int = 0,
    # stability
    n_boot: int = 200,
    sample_fraction: float = 0.8,
    # panel rule
    freq_threshold: float = 0.6,
    sign_threshold: float = 0.8,
    max_genes: int = 12,
    # safety
    min_events: int = 10,
    min_var: float = 1e-12,
    robust: bool = False,
    # NEW: allow reuse (skip heavy steps)
    best_params: dict | None = None,
    df_stability: pd.DataFrame | None = None,
    df_cv: pd.DataFrame | None = None,
    keep_nonzero: bool = True,
) -> dict:
    """
    One-call helper:
      1) CV select penalizer + l1_ratio            (skipped if best_params provided)
      2) Stability selection with best params      (skipped if df_stability provided)
      3) Choose a compact panel
      4) Refit final penalized Cox on full data for weights

    Returns dict with:
      best_params, df_cv, df_stability, panel, weights_df
    """
    _require_lifelines()

    # ---------- (A) CV (optional) ----------
    if best_params is None:
        df_cv_run, best = cox_penalizer_cv(
            adata,
            genes=genes,
            time_col=time_col,
            event_col=event_col,
            layer=layer,
            standardize=standardize,
            penalizers=penalizers,
            l1_ratios=l1_ratios,
            n_splits=n_splits,
            n_repeats=n_repeats,
            seed=seed,
            min_events=min_events,
            min_var=min_var,
            robust=robust,
        )
        if df_cv is None:
            df_cv = df_cv_run
    else:
        # user provided best params; do not recompute
        best = {"penalizer": float(best_params["penalizer"]), "l1_ratio": float(best_params["l1_ratio"])}

    # ---------- (B) Stability (optional) ----------
    if df_stability is None:
        df_stab = cox_stability_selection(
            adata,
            genes=genes,
            time_col=time_col,
            event_col=event_col,
            layer=layer,
            standardize=standardize,
            penalizer=best["penalizer"],
            l1_ratio=best["l1_ratio"],
            n_boot=n_boot,
            sample_fraction=sample_fraction,
            seed=seed,
            min_events=min_events,
            min_var=min_var,
            robust=robust,
        )
    else:
        df_stab = df_stability.copy()

    # ---------- (C) Panel selection ----------
    df_pick = df_stab.copy()
    df_pick = df_pick[np.isfinite(df_pick["freq"])].copy()
    df_pick = df_pick[df_pick["freq"] >= float(freq_threshold)].copy()
    if "sign_consistency" in df_pick.columns:
        df_pick = df_pick[
            (~np.isfinite(df_pick["sign_consistency"])) | (df_pick["sign_consistency"] >= float(sign_threshold))
        ]

    if df_pick.shape[0] == 0:
        warn("recommended_cox_panel: no genes pass thresholds; falling back to top genes by stability frequency.")
        df_pick = df_stab.copy()

    # rank + take top K
    sort_cols = ["freq"]
    asc = [False]
    if "sign_consistency" in df_pick.columns:
        sort_cols.append("sign_consistency")
        asc.append(False)

    df_pick = df_pick.sort_values(sort_cols, ascending=asc, na_position="last")
    panel = df_pick["gene"].astype(str).head(int(max_genes)).tolist()

    # ---------- (D) Final penalized fit (cheap) ----------
    weights_df = cox_fit_penalized(
        adata,
        genes=panel,
        time_col=time_col,
        event_col=event_col,
        layer=layer,
        standardize=standardize,
        penalizer=float(best["penalizer"]),
        l1_ratio=float(best["l1_ratio"]),
        robust=bool(robust),
        min_events=int(min_events),
        min_var=float(min_var),
        adjust="none",             # multivariate penalized p/q are not super meaningful
        keep_nonzero=bool(keep_nonzero),
    )

    return {
        "best_params": best,
        "df_cv": df_cv,
        "df_stability": df_stab,
        "panel": panel,
        "weights_df": weights_df,
    }


def permutation_importance_cindex(
    adata,
    *,
    weights: pd.DataFrame,
    time_col: str,
    event_col: str,
    layer: str = "log1p_cpm",
    score_key: str = "sig_tmp",
    n_repeats: int = 20,
    seed: int = 0,
):
    from lifelines.utils import concordance_index

    # ---- infer columns ----
    w = weights.copy()
    if "beta" not in w.columns and "coef" in w.columns:
        w = w.rename(columns={"coef": "beta"})
    if not {"gene", "beta"}.issubset(w.columns):
        raise ValueError("weights must contain columns: gene and beta (or coef).")

    genes = [g for g in w["gene"].astype(str).tolist() if g in adata.var_names]
    w = w[w["gene"].astype(str).isin(genes)].copy()
    if len(genes) == 0:
        raise ValueError("None of the signature genes are present in adata.var_names.")

    # ---- baseline score + baseline c-index ----
    bk.tl.signature_score(adata, weights=w, layer=layer, key_added=score_key)
    t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(float)
    e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(float)
    e = np.where(e > 0, 1, 0).astype(int)
    s = pd.to_numeric(adata.obs[score_key], errors="coerce").to_numpy(float)

    ok = np.isfinite(t) & np.isfinite(s) & np.isfinite(e)
    t0, e0, s0 = t[ok], e[ok], s[ok]

    # c-index expects higher score = longer survival; for a “risk score” invert sign if needed
    c_base = float(concordance_index(t0, -s0, event_observed=e0))

    # ---- pull expression matrix for signature genes ----
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    gidx = adata.var_names.get_indexer(genes)
    M = X[:, gidx].toarray() if hasattr(X, "toarray") else np.asarray(X[:, gidx], dtype=float)
    M = M[ok, :]  # align to ok samples

    beta = w.set_index("gene").loc[genes, "beta"].to_numpy(float)

    rng = np.random.default_rng(int(seed))
    rows = []

    for j, g in enumerate(genes):
        drops = []
        for r in range(int(n_repeats)):
            Mp = M.copy()
            Mp[:, j] = rng.permutation(Mp[:, j])  # permute that gene only
            sp = Mp @ beta
            c_perm = float(concordance_index(t0, -sp, event_observed=e0))
            drops.append(c_base - c_perm)

        rows.append(
            {
                "gene": g,
                "c_index_base": c_base,
                "c_index_perm_mean": c_base - float(np.mean(drops)),
                "delta_cindex_mean": float(np.mean(drops)),
                "delta_cindex_sd": float(np.std(drops)),
            }
        )

    df_perm = pd.DataFrame(rows).sort_values("delta_cindex_mean", ascending=False).reset_index(drop=True)
    return df_perm


def cox_fit_penalized(
    adata: ad.AnnData,
    *,
    genes: Sequence[str],
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    penalizer: float = 0.1,
    l1_ratio: float = 0.0,
    robust: bool = False,
    min_events: int = 10,
    min_var: float = 1e-12,
    adjust: Literal["bh", "none"] = "bh",
    keep_nonzero: bool = True,
    key_added: str | None = None,
) -> pd.DataFrame:
    """
    Fit a penalized Cox model on the full dataset for a candidate gene set.

    Uses lifelines CoxPHFitter(penalizer=..., l1_ratio=...).
    Returns a per-gene table suitable as signature weights.

    Returns
    -------
    pd.DataFrame with columns:
      gene, beta, HR, CI_lower, CI_upper, se, z, pval, qval, n, n_events
    """
    CoxPHFitter, _ = _require_lifelines()

    df_all, genes_kept = _prepare_cox_df(
        adata,
        genes,
        time_col=time_col,
        event_col=event_col,
        layer=layer,
        standardize=standardize,
        min_var=min_var,
    )

    n = int(df_all.shape[0])
    n_events = int(df_all[event_col].sum())
    if n_events < int(min_events):
        raise ValueError(f"Not enough events: n_events={n_events} < min_events={min_events}")

    # Fit
    cph = CoxPHFitter(penalizer=float(penalizer), l1_ratio=float(l1_ratio))
    cph.fit(
        df_all,
        duration_col=time_col,
        event_col=event_col,
        robust=bool(robust),
        show_progress=False,
    )

    summ = cph.summary.copy()

    out = pd.DataFrame(
        {
            "gene": summ.index.astype(str),
            "beta": summ["coef"].to_numpy(dtype=float),
            "se": summ["se(coef)"].to_numpy(dtype=float) if "se(coef)" in summ.columns else np.nan,
            "z": summ["z"].to_numpy(dtype=float) if "z" in summ.columns else np.nan,
            "pval": summ["p"].to_numpy(dtype=float) if "p" in summ.columns else np.nan,
        }
    )

    out["HR"] = np.exp(out["beta"].to_numpy(dtype=float))

    if "exp(coef) lower 95%" in summ.columns and "exp(coef) upper 95%" in summ.columns:
        out["CI_lower"] = summ["exp(coef) lower 95%"].to_numpy(dtype=float)
        out["CI_upper"] = summ["exp(coef) upper 95%"].to_numpy(dtype=float)
    else:
        lo = out["beta"] - 1.96 * out["se"]
        hi = out["beta"] + 1.96 * out["se"]
        out["CI_lower"] = np.exp(lo.to_numpy(dtype=float))
        out["CI_upper"] = np.exp(hi.to_numpy(dtype=float))

    out["n"] = n
    out["n_events"] = int(n_events)

    if keep_nonzero:
        out = out.loc[np.abs(out["beta"].to_numpy(dtype=float)) > 1e-12].copy()

    if adjust == "bh" and out.shape[0] > 0:
        out["qval"] = _bh_fdr(out["pval"].to_numpy(dtype=float))
    else:
        out["qval"] = np.nan

    if out.shape[0] > 0:
        out = out.assign(_abs=np.abs(out["beta"])).sort_values(
            ["_abs", "pval"], ascending=[False, True]
        ).drop(columns=["_abs"]).reset_index(drop=True)

    if key_added is not None:
        adata.uns.setdefault("cox_fit_penalized", {})
        adata.uns["cox_fit_penalized"][key_added] = {
            "params": {
                "time_col": time_col,
                "event_col": event_col,
                "layer": layer,
                "genes_requested": int(len(list(genes))),
                "genes_used": int(len(genes_kept)),
                "standardize": bool(standardize),
                "penalizer": float(penalizer),
                "l1_ratio": float(l1_ratio),
                "robust": bool(robust),
                "min_events": int(min_events),
                "min_var": float(min_var),
                "adjust": adjust,
                "keep_nonzero": bool(keep_nonzero),
            },
            "results": out,
        }

    return out




def cox_univariate(
    adata: AnnData,
    *,
    time_key: str = "OS.time",
    event_key: str = "OS",
    x_keys: Sequence[str] = ("mp_heterogeneity_entropy",),
    x_mode: str = "continuous",   # "continuous" | "binary"
    q: Tuple[float, float] = (0.25, 0.75),
    groupby: str = "Project_ID",
    binning: str = "global",      # "global" | "within_group"
    strata: Union[str, Sequence[str], None] = "Project_ID",
    min_per_stratum: int = 10,
    covariates: Optional[Sequence[str]] = None,
    robust: bool = True,
    penalizer: float = 0.0,
    out_key: str = "cox_univariate_multi",
):
    """
    Fit Cox models for multiple predictors (x_keys) under the same settings.
    Stores a combined tidy summary in `adata.uns[out_key]["summary"]`.

    x_mode:
      - "continuous": predictor is used as numeric (recommended for entropy)
      - "binary": predictor is binned by quantiles (q) and only extremes are kept

    Stratification:
      - Use `strata="Project_ID"` for pan-cancer analyses.

    Returns
    -------
    (results_df, models)
        results_df: tidy DataFrame (one row per x_key + term) with HR and p
        models: dict {x_key: CoxPHFitter}
    """
    from lifelines import CoxPHFitter

    obs = adata.obs
    if time_key not in obs or event_key not in obs:
        raise KeyError(f"{time_key!r} and/or {event_key!r} not found in adata.obs")

    x_keys = list(x_keys)
    for x in x_keys:
        if x not in obs:
            raise KeyError(f"x_key {x!r} not found in adata.obs")

    covariates = list(covariates) if covariates is not None else []
    for c in covariates:
        if c not in obs:
            raise KeyError(f"covariate {c!r} not found in adata.obs")

    # strata columns
    strata_cols: list[str] = []
    if strata is not None:
        strata_cols = [strata] if isinstance(strata, str) else list(strata)
        for s in strata_cols:
            if s not in obs:
                raise KeyError(f"strata column {s!r} not found in adata.obs")

    # base frame (shared)
    base_cols = [time_key, event_key] + covariates + strata_cols
    if binning == "within_group":
        if groupby not in obs:
            raise KeyError(f"{groupby!r} not found in adata.obs (required for binning='within_group')")
        base_cols += [groupby]

    base = obs[base_cols].copy()
    base = base.replace([np.inf, -np.inf], np.nan).dropna(subset=[time_key, event_key])
    base = base[base[time_key].astype(float) > 0]
    base[event_key] = base[event_key].astype(float).astype(int)

    def _bin_extremes(s: pd.Series, lo: float, hi: float) -> pd.Series:
        out = pd.Series(np.nan, index=s.index, dtype=float)
        out[s <= lo] = 0.0
        out[s >= hi] = 1.0
        return out

    all_rows = []
    models = {}

    for x_key in x_keys:
        df = base.join(obs[[x_key]], how="inner")
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_key])

        if x_mode == "continuous":
            df["x"] = df[x_key].astype(float)
            df_used = df.dropna(subset=["x"]).copy()
            term = "x"
            formula_terms = ["x"] + covariates

        elif x_mode == "binary":
            x = df[x_key].astype(float)

            if binning == "global":
                lo, hi = x.quantile(q[0]), x.quantile(q[1])
                df["x_bin"] = _bin_extremes(x, lo, hi)

            elif binning == "within_group":
                df["x_bin"] = np.nan
                for g, idx in df.groupby(groupby, observed=False).groups.items():
                    idx = list(idx)
                    sub = df.loc[idx]
                    if sub.shape[0] < 10:
                        continue
                    lo, hi = sub[x_key].astype(float).quantile(q[0]), sub[x_key].astype(float).quantile(q[1])
                    df.loc[idx, "x_bin"] = _bin_extremes(sub[x_key].astype(float), lo, hi)
            else:
                raise ValueError("binning must be one of {'global','within_group'}")

            df_used = df.dropna(subset=["x_bin"]).copy()
            df_used["x_bin"] = df_used["x_bin"].astype(int)
            term = "x_bin"
            formula_terms = ["x_bin"] + covariates

        else:
            raise ValueError("x_mode must be 'continuous' or 'binary'")

        # strata filtering (after binning)
        if strata_cols:
            df_used = df_used.dropna(subset=strata_cols)
            s0 = strata_cols[0]
            counts = df_used[s0].astype("string").value_counts()
            keep_levels = counts[counts >= min_per_stratum].index
            df_used = df_used[df_used[s0].astype("string").isin(keep_levels)]

        # fit
        cph = CoxPHFitter(penalizer=penalizer)
        fit_kwargs = dict(
            duration_col=time_key,
            event_col=event_key,
            robust=robust,
            formula=" + ".join(formula_terms),
        )
        if strata_cols:
            fit_kwargs["strata"] = strata_cols

        # need at least some events
        if df_used.shape[0] < 30 or df_used[event_key].sum() < 5:
            # store NA row and skip model
            all_rows.append({
                "x_key": x_key,
                "term": term,
                "coef": np.nan,
                "HR": np.nan,
                "HR_lower_95": np.nan,
                "HR_upper_95": np.nan,
                "p": np.nan,
                "n_used": int(df_used.shape[0]),
                "events": int(df_used[event_key].sum()),
            })
            continue

        cph.fit(df_used, **fit_kwargs)
        models[x_key] = cph

        summ = cph.summary.copy()
        # pull the predictor term only (x or x_bin)
        if term not in summ.index:
            # should not happen; but be defensive
            coef = hr = lo95 = hi95 = pval = np.nan
        else:
            coef = float(summ.loc[term, "coef"])
            hr = float(np.exp(coef))
            lo95 = float(np.exp(summ.loc[term, "coef lower 95%"]))
            hi95 = float(np.exp(summ.loc[term, "coef upper 95%"]))
            pval = float(summ.loc[term, "p"])

        all_rows.append({
            "x_key": x_key,
            "term": term,
            "coef": coef,
            "HR": hr,
            "HR_lower_95": lo95,
            "HR_upper_95": hi95,
            "p": pval,
            "n_used": int(df_used.shape[0]),
            "events": int(df_used[event_key].sum()),
        })

    results_df = pd.DataFrame(all_rows)

    adata.uns[out_key] = {
        "params": {
            "time_key": time_key,
            "event_key": event_key,
            "x_keys": x_keys,
            "x_mode": x_mode,
            "q": q if x_mode == "binary" else None,
            "binning": binning if x_mode == "binary" else None,
            "groupby": groupby if x_mode == "binary" else None,
            "strata": strata,
            "min_per_stratum": min_per_stratum,
            "covariates": covariates,
            "robust": robust,
            "penalizer": penalizer,
        },
        "summary": results_df,
    }

    return results_df, models


def cox_interaction(
    adata: AnnData,
    *,
    time_key: str = "OS.time",
    event_key: str = "OS",
    x_key: str = "Neuroendocrine_score",
    z_key: str = "mp_heterogeneity_entropy",
    groupby: str = "Project_ID",
    strata: Union[str, Sequence[str], None] = "Project_ID",
    q: Tuple[float, float] = (0.25, 0.75),
    binning: str = "global",  # "global" | "within_group"
    keep: str = "extremes",   # "extremes" only | "all" (keeps middle as NaN unless you handle)
    min_per_stratum: int = 10,
    covariates: Optional[Sequence[str]] = None,
    robust: bool = True,
    penalizer: float = 0.0,
    out_key: str = "cox_interaction",
):
    """
    Fit Cox proportional hazards model with an interaction term:
        x_bin + z_bin + x_bin:z_bin
    optionally stratified (e.g., by tumor type Project_ID).

    Binning (default):
      - Define x_bin and z_bin using quantiles (q[0], q[1]).
      - Keep only extremes (low/high) by default.

    Parameters
    ----------
    adata
        AnnData with survival data and covariates in `.obs`.
    time_key, event_key
        Survival time and event indicator columns in `.obs`.
        event_key should be 0/1 (will be coerced).
    x_key
        Primary variable (e.g., NE score).
    z_key
        Modifier variable (e.g., heterogeneity).
    groupby
        Grouping key used when binning='within_group'. Often same as tumor type.
    strata
        Column name or list of column names used for Cox stratification.
        Use "Project_ID" for tumor-type stratification.
        Set None for no stratification.
    q
        Quantile cutoffs for low/high.
    binning
        "global": cutoffs computed across all samples.
        "within_group": cutoffs computed within each `groupby` group.
    keep
        "extremes": keep only low/high (drop middle).
        "all": keep all rows (middle bins become NaN; you typically don’t want this).
    min_per_stratum
        Minimum samples per stratum retained (after filtering).
    covariates
        Additional covariate column names from `.obs` to include in the model.
    robust
        Use robust variance estimates (recommended for TCGA).
    penalizer
        L2 penalizer for CoxPHFitter (helps with collinearity).
    out_key
        Store results in `adata.uns[out_key]`.

    Returns
    -------
    (cph, df_used)
        lifelines.CoxPHFitter object and the dataframe used for fitting.

    Notes
    -----
    Model formula: x_bin + z_bin + x_bin:z_bin [+ covariates]
    Interpretation:
      - x_bin: effect of x when z_bin==0 (low)
      - z_bin: effect of z when x_bin==0 (low)
      - interaction: additional effect when both high
    """
    # lazy import to keep dependency optional until used
    from lifelines import CoxPHFitter

    obs = adata.obs

    needed = [time_key, event_key, x_key, z_key]
    for k in needed:
        if k not in obs:
            raise KeyError(f"{k!r} not found in adata.obs")

    if binning == "within_group" and groupby not in obs:
        raise KeyError(f"{groupby!r} not found in adata.obs (required for binning='within_group')")

    covariates = list(covariates) if covariates is not None else []
    for c in covariates:
        if c not in obs:
            raise KeyError(f"covariate {c!r} not found in adata.obs")

    # Build dataframe
    cols = [time_key, event_key, x_key, z_key] + covariates
    if isinstance(strata, str):
        if strata is not None:
            cols += [strata]
    elif strata is not None:
        cols += list(strata)
    if binning == "within_group":
        cols += [groupby]

    df = obs[cols].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[time_key, event_key, x_key, z_key])

    # Coerce event to 0/1
    df[event_key] = df[event_key].astype(float).astype(int)

    # Keep positive time
    df = df[df[time_key].astype(float) > 0]

    # --- Binning helpers
    def _bin_extremes(s: pd.Series, lo: float, hi: float) -> pd.Series:
        out = pd.Series(np.nan, index=s.index, dtype=float)
        out[s <= lo] = 0.0
        out[s >= hi] = 1.0
        return out

    if binning == "global":
        x_lo, x_hi = df[x_key].quantile(q[0]), df[x_key].quantile(q[1])
        z_lo, z_hi = df[z_key].quantile(q[0]), df[z_key].quantile(q[1])
        df["x_bin"] = _bin_extremes(df[x_key], x_lo, x_hi)
        df["z_bin"] = _bin_extremes(df[z_key], z_lo, z_hi)

    elif binning == "within_group":
        df["x_bin"] = np.nan
        df["z_bin"] = np.nan
        for g, idx in df.groupby(groupby, observed=False).groups.items():
            idx = list(idx)
            sub = df.loc[idx]
            if sub.shape[0] < 10:
                continue
            x_lo, x_hi = sub[x_key].quantile(q[0]), sub[x_key].quantile(q[1])
            z_lo, z_hi = sub[z_key].quantile(q[0]), sub[z_key].quantile(q[1])
            df.loc[idx, "x_bin"] = _bin_extremes(sub[x_key], x_lo, x_hi)
            df.loc[idx, "z_bin"] = _bin_extremes(sub[z_key], z_lo, z_hi)
    else:
        raise ValueError("binning must be one of {'global','within_group'}")

    if keep == "extremes":
        df = df.dropna(subset=["x_bin", "z_bin"])
    elif keep != "all":
        raise ValueError("keep must be 'extremes' or 'all'")

    df["x_bin"] = df["x_bin"].astype(int)
    df["z_bin"] = df["z_bin"].astype(int)
    df["x_z"] = (df["x_bin"] * df["z_bin"]).astype(int)

    # Filter strata with too few samples (after binning)
    if strata is not None:
        if isinstance(strata, str):
            strata_cols = [strata]
        else:
            strata_cols = list(strata)

        # drop rows with missing strata values
        df = df.dropna(subset=strata_cols)

        # enforce min samples per stratum (based on first stratum col)
        # (Simple rule: keep strata levels with >= min_per_stratum)
        s0 = strata_cols[0]
        counts = df[s0].astype("string").value_counts()
        keep_levels = counts[counts >= min_per_stratum].index
        df = df[df[s0].astype("string").isin(keep_levels)]

    # --- Model formula
    terms = ["x_bin", "z_bin", "x_z"] + covariates
    formula = " + ".join(terms)

    # Fit Cox
    cph = CoxPHFitter(penalizer=penalizer)
    fit_kwargs = dict(
        duration_col=time_key,
        event_col=event_key,
        formula=formula,
        robust=robust,
    )
    if strata is not None:
        fit_kwargs["strata"] = [strata] if isinstance(strata, str) else list(strata)


    # ---- critical lifelines fix: reset index + cast strata col(s) to str
    df = df.reset_index(drop=True)
    if strata is not None:
        strata_cols = fit_kwargs["strata"]
        for s in strata_cols:
            if s in df.columns:
                df[s] = df[s].astype(str)

    cph.fit(df, **fit_kwargs)

    # Store summary in adata.uns
    summ = cph.summary.copy()
    summ["HR"] = np.exp(summ["coef"])
    summ["HR_lower_95"] = np.exp(summ["coef lower 95%"])
    summ["HR_upper_95"] = np.exp(summ["coef upper 95%"])
    summ = summ[["coef", "HR", "HR_lower_95", "HR_upper_95", "p"]]

    adata.uns[out_key] = {
        "params": {
            "time_key": time_key,
            "event_key": event_key,
            "x_key": x_key,
            "z_key": z_key,
            "groupby": groupby,
            "strata": strata,
            "q": q,
            "binning": binning,
            "keep": keep,
            "covariates": covariates,
            "robust": robust,
            "penalizer": penalizer,
            "n_used": int(df.shape[0]),
        },
        "summary": summ,
    }

    return cph, df



def surv_1d_bins(
    adata: AnnData,
    *,
    x_key: str,
    time_key: str,
    event_key: str,
    groupby: Optional[str] = None,
    binning: str = "global",  # "global" | "within_group"
    q: Tuple[float, float] = (0.25, 0.75),
    keep: str = "extremes",   # "extremes" | "all"
    labels: Tuple[str, str] = ("low", "high"),
    out_col: str = "surv_group_1d",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Create survival bins for a single score x_key (e.g., heterogeneity entropy).

    Bins:
      - low:  x <= q[0]
      - high: x >= q[1]
      - middle becomes NA (dropped if keep='extremes')

    Parameters
    ----------
    x_key
        Continuous score in adata.obs.
    time_key, event_key
        Survival time and event columns in adata.obs.
    groupby
        Required if binning="within_group" (e.g. "Project_ID").
    binning
        "global": quantiles computed across all samples.
        "within_group": quantiles computed within each groupby category.
    q
        Quantiles defining low/high.
    keep
        "extremes": keep only low/high (recommended).
        "all": keep all, with middle as NA in out_col.
    labels
        Names for (low, high).
    out_col
        Column in adata.obs storing labels.
    copy
        Return a copy if True.

    Returns
    -------
    If copy=True returns AnnData, else None (modifies adata in place).
    """
    ad = adata.copy() if copy else adata
    obs = ad.obs

    for k in [x_key, time_key, event_key]:
        if k not in obs:
            raise KeyError(f"{k!r} not found in adata.obs")

    if binning == "within_group":
        if groupby is None:
            raise ValueError("groupby must be provided when binning='within_group'")
        if groupby not in obs:
            raise KeyError(f"{groupby!r} not found in adata.obs")

    df = obs[[x_key, time_key, event_key] + ([groupby] if groupby else [])].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_key, time_key, event_key])
    df = df[df[time_key].astype(float) > 0]

    low_label, high_label = labels

    def _bin_extremes(s: pd.Series, lo: float, hi: float) -> pd.Series:
        out = pd.Series(pd.NA, index=s.index, dtype="string")
        out[s <= lo] = low_label
        out[s >= hi] = high_label
        return out

    x_bin = pd.Series(pd.NA, index=df.index, dtype="string")

    if binning == "global":
        lo, hi = df[x_key].quantile(q[0]), df[x_key].quantile(q[1])
        x_bin = _bin_extremes(df[x_key].astype(float), lo, hi)

    elif binning == "within_group":
        for g, idx in df.groupby(groupby, observed=False).groups.items():
            idx = list(idx)
            sub = df.loc[idx]
            if sub.shape[0] < 10:
                continue
            lo, hi = sub[x_key].quantile(q[0]), sub[x_key].quantile(q[1])
            x_bin.loc[idx] = _bin_extremes(sub[x_key].astype(float), lo, hi)
    else:
        raise ValueError("binning must be one of {'global','within_group'}")

    # attach to full obs
    obs[out_col] = pd.NA
    obs.loc[x_bin.index, out_col] = x_bin

    if keep == "extremes":
        valid = obs[out_col].isin([low_label, high_label])
        obs.loc[~valid, out_col] = pd.NA
    elif keep != "all":
        raise ValueError("keep must be 'extremes' or 'all'")

    return ad if copy else None


def surv_2x2_bins(
    adata: AnnData,
    *,
    x_key: str,
    z_key: str,
    time_key: str,
    event_key: str,
    groupby: Optional[str] = None,
    binning: str = "global",  # "global" | "within_group"
    q: Tuple[float, float] = (0.25, 0.75),
    keep: str = "extremes",   # "extremes" | "all"
    out_col: str = "surv_group_2x2",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Create 2x2 survival groups based on low/high bins of x_key and z_key.

    Groups:
      X_low & Z_low
      X_low & Z_high
      X_high & Z_low
      X_high & Z_high

    Parameters
    ----------
    x_key, z_key
        Continuous scores in adata.obs (e.g. NE score, heterogeneity entropy).
    time_key, event_key
        Survival time and event columns in adata.obs.
    groupby
        Required if binning="within_group" (e.g. "Project_ID").
    binning
        "global": quantiles computed across all samples.
        "within_group": quantiles computed within each groupby category.
    q
        Quantiles defining low/high cutoffs.
    keep
        "extremes": drop middle samples (between q[0] and q[1]).
        "all": keep all samples, middle gets NA for bins (not recommended for plotting).
    out_col
        Column in adata.obs storing group labels.
    copy
        Return a copy if True.

    Returns
    -------
    If copy=True returns AnnData, else None (modifies adata in place).
    """
    ad = adata.copy() if copy else adata
    obs = ad.obs

    for k in [x_key, z_key, time_key, event_key]:
        if k not in obs:
            raise KeyError(f"{k!r} not found in adata.obs")

    if binning == "within_group":
        if groupby is None:
            raise ValueError("groupby must be provided when binning='within_group'")
        if groupby not in obs:
            raise KeyError(f"{groupby!r} not found in adata.obs")

    df = obs[[x_key, z_key, time_key, event_key] + ([groupby] if groupby else [])].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_key, z_key, time_key, event_key])
    df = df[df[time_key].astype(float) > 0]

    def _bin_extremes(s: pd.Series, lo: float, hi: float) -> pd.Series:
        out = pd.Series(pd.NA, index=s.index, dtype="string")
        out[s <= lo] = "low"
        out[s >= hi] = "high"
        return out

    x_bin = pd.Series(pd.NA, index=df.index, dtype="string")
    z_bin = pd.Series(pd.NA, index=df.index, dtype="string")

    if binning == "global":
        x_lo, x_hi = df[x_key].quantile(q[0]), df[x_key].quantile(q[1])
        z_lo, z_hi = df[z_key].quantile(q[0]), df[z_key].quantile(q[1])
        x_bin = _bin_extremes(df[x_key].astype(float), x_lo, x_hi)
        z_bin = _bin_extremes(df[z_key].astype(float), z_lo, z_hi)

    elif binning == "within_group":
        for g, idx in df.groupby(groupby, observed=False).groups.items():
            idx = list(idx)
            sub = df.loc[idx]
            if sub.shape[0] < 10:
                continue
            x_lo, x_hi = sub[x_key].quantile(q[0]), sub[x_key].quantile(q[1])
            z_lo, z_hi = sub[z_key].quantile(q[0]), sub[z_key].quantile(q[1])
            x_bin.loc[idx] = _bin_extremes(sub[x_key].astype(float), x_lo, x_hi)
            z_bin.loc[idx] = _bin_extremes(sub[z_key].astype(float), z_lo, z_hi)

    else:
        raise ValueError("binning must be one of {'global','within_group'}")

    # attach bins into obs (align to full obs index)
    obs["x_bin_2x2"] = pd.NA
    obs["z_bin_2x2"] = pd.NA
    obs.loc[x_bin.index, "x_bin_2x2"] = x_bin
    obs.loc[z_bin.index, "z_bin_2x2"] = z_bin

    if keep == "extremes":
        valid = obs["x_bin_2x2"].isin(["low", "high"]) & obs["z_bin_2x2"].isin(["low", "high"])
    elif keep == "all":
        valid = pd.Series(True, index=obs.index)
    else:
        raise ValueError("keep must be 'extremes' or 'all'")

    # build group label
    group = pd.Series(pd.NA, index=obs.index, dtype="string")
    group[valid] = "X_" + obs.loc[valid, "x_bin_2x2"].astype("string") + " & Z_" + obs.loc[valid, "z_bin_2x2"].astype("string")
    obs[out_col] = group

    return ad if copy else None

