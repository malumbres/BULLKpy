from __future__ import annotations

from pathlib import Path
from typing import Sequence, Literal, Optional
import json
import numpy as np
import pandas as pd
import anndata as ad

from dataclasses import dataclass

from ..logging import info, warn

# optional: use your package logger if present
try:
    from ..logging import info, warn
except Exception:  # pragma: no cover
    def info(x): print(f"INFO: {x}")
    def warn(x): print(f"WARNING: {x}")

from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


############################################
## HELPERS
############################################

def _get_matrix(adata, layer, use="samples"):
    X = adata.layers[layer] if (layer is not None and layer in adata.layers) else adata.X
    return X

# Extract expression matrix for genes (train only)
def get_expr_matrix(adata, genes, *, layer="log1p_cpm"):
    genes = [g for g in genes if g in adata.var_names]
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    idx = adata.var_names.get_indexer(genes)
    M = X[:, idx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, float)
    return genes, M

# Apply the label mask (drop unknowns)
def mask_known_labels(adata, label_col, label):
    yraw = adata.obs[label_col].astype("object")
    ok = yraw.notna().to_numpy()
    y = (yraw[ok].astype(str).to_numpy() == label).astype(int)
    return ok, y


############################################
## FUNCTIONS
############################################

def filter_genes_var(
    adata,
    *,
    label_col,
    label="NR",
    layer="log1p_cpm",
    min_var=0.05,
    top_var=None,
):
    """
    Basic variance-based gene filter.

    Returns
    -------
    df_var : pd.DataFrame
        Columns:
            - gene
            - var
        Sorted by descending variance.
    """
    # mask samples with known labels
    ok, y = mask_known_labels(adata, label_col, label)

    genes_all = list(adata.var_names)
    genes, M = get_expr_matrix(adata, genes_all, layer=layer)
    M = M[ok, :]

    # variance per gene
    var = np.nanvar(M, axis=0)

    df = pd.DataFrame({
        "gene": genes,
        "var": var,
    })

    # apply min_var
    df = df[df["var"] >= float(min_var)].copy()

    # sort by variance
    df = df.sort_values("var", ascending=False).reset_index(drop=True)

    # optionally keep top_var
    if top_var is not None and df.shape[0] > int(top_var):
        df = df.head(int(top_var)).copy()

    print(f"Basic filter: kept {df.shape[0]} genes (min_var={min_var}, top_var={top_var})")

    return df

def rank_genes_univariate_pr(
    adata,
    genes,
    *,
    label_col,
    label,
    layer="log1p_cpm",
):
    """
    Univariate PR-AUC ranking that is direction-aware.

    Returns:
      gene, pr_auc, direction
        direction = +1 → high expression predicts NR
        direction = -1 → low expression predicts NR
    """
    ok, y = mask_known_labels(adata, label_col, label)
    genes, M = get_expr_matrix(adata, genes, layer=layer)
    M = M[ok, :]

    rows = []
    for j, g in enumerate(genes):
        x = M[:, j]
        if not np.all(np.isfinite(x)):
            continue

        try:
            ap_pos = average_precision_score(y, x)
            ap_neg = average_precision_score(y, -x)
        except Exception:
            continue

        rows.append({
            "gene": g,
            "pr_auc": max(ap_pos, ap_neg),
            "ap_pos": ap_pos,
            "ap_neg": ap_neg,
            "direction": +1 if ap_pos >= ap_neg else -1,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("pr_auc", ascending=False).reset_index(drop=True)
    return df


def gene_corr_redundancy(
    adata,
    genes_ranked,
    *,
    label_col,
    label,
    layer="log1p_cpm",
    corr_threshold=0.8,
    prefer_genes=None,                 # e.g. {"SOCS2", "CXCL13"}
    prefer_regex=None,                 # e.g. r"^(SOCS2|CXCL13)$"
    penalize_regex=None,               # e.g. r"(-AS\d*$|^RP\d+|^LINC|antisense)"
    # behavior:
    allow_prefer_swap=True,            # if a preferred gene is redundant with a kept gene, swap them
):
    """
    Greedy redundancy filter with audit trail.

    Returns
    -------
    df_audit : pd.DataFrame with columns:
        gene, main_redundant_gene, redundancy_value, keep
    kept_genes : list[str]
    """
    # ---- subset to labeled samples (same as your current helper) ----
    ok, _ = mask_known_labels(adata, label_col, label)
    genes, M = get_expr_matrix(adata, genes_ranked, layer=layer)
    M = M[ok, :]

    # ---- correlation matrix ----
    dfX = pd.DataFrame(M, columns=genes)
    C = dfX.corr().abs().fillna(0.0)

    prefer_set = set([str(g) for g in prefer_genes]) if prefer_genes is not None else set()

    def is_preferred(g: str) -> bool:
        if g in prefer_set:
            return True
        if prefer_regex is not None:
            import re
            return re.search(prefer_regex, g) is not None
        return False

    def is_penalized(g: str) -> bool:
        if penalize_regex is None:
            return False
        import re
        return re.search(penalize_regex, g) is not None

    kept = []
    rows = []

    for g in genes:
        if len(kept) == 0:
            kept.append(g)
            rows.append({
                "gene": g,
                "main_redundant_gene": None,
                "redundancy_value": 0.0,
                "keep": True,
            })
            continue

        # compute redundancy vs currently kept genes
        vals = C.loc[g, kept].to_numpy(float)
        j = int(np.argmax(vals))
        g_main = kept[j]
        rmax = float(vals[j])

        if rmax < float(corr_threshold):
            kept.append(g)
            rows.append({
                "gene": g,
                "main_redundant_gene": g_main,
                "redundancy_value": rmax,
                "keep": True,
            })
            continue

        # redundant: by default filtered out
        kept_flag = False

        # OPTIONAL: swap logic to protect “preferred” genes (e.g., SOCS2)
        if allow_prefer_swap:
            g_pref = is_preferred(g)
            main_pref = is_preferred(g_main)

            # If g is preferred and the kept representative is NOT preferred, swap
            if g_pref and (not main_pref):
                # remove g_main from kept, add g instead
                kept.remove(g_main)
                kept.append(g)

                # mark current gene as kept
                kept_flag = True

                # also retroactively mark g_main as filtered out in audit table
                # (we'll update after building df)
                rows.append({
                    "gene": g,
                    "main_redundant_gene": g_main,
                    "redundancy_value": rmax,
                    "keep": True,
                })
            # If main is penalized (antisense/lnc) and g is not penalized, swap
            elif (is_penalized(g_main) and (not is_penalized(g))):
                kept.remove(g_main)
                kept.append(g)
                kept_flag = True
                rows.append({
                    "gene": g,
                    "main_redundant_gene": g_main,
                    "redundancy_value": rmax,
                    "keep": True,
                })

        if not kept_flag:
            rows.append({
                "gene": g,
                "main_redundant_gene": g_main,
                "redundancy_value": rmax,
                "keep": False,
            })

    df_audit = pd.DataFrame(rows)

    # If swaps happened, the earlier representative (g_main) might still be marked kept=True.
    # Fix by enforcing: kept flag == membership in final kept list
    kept_set = set(kept)
    df_audit["keep"] = df_audit["gene"].astype(str).isin(kept_set)

    # Add a quick summary
    print(f"Redundancy filter (audit): {len(genes_ranked)} -> {len(kept)} genes (thr={corr_threshold})")

    return df_audit, kept


def signature_score(
    adata: ad.AnnData,
    *,
    weights: pd.DataFrame,
    gene_col: str = "gene",
    beta_col: str = "beta",          # can also be "logHR" if you prefer
    layer: str | None = "log1p_cpm",
    standardize: bool = True,
    missing: Literal["drop", "error"] = "drop",
    key_added: str = "signature_score",
) -> np.ndarray:
    """
    Compute a weighted gene-expression signature score per sample.

    Score_i = sum_g beta_g * z(x_{ig})
    where z() is optional standardization per gene across samples.

    Parameters
    ----------
    adata
        AnnData with samples in rows (obs) and genes in columns (var).
    weights
        DataFrame containing at least [gene_col, beta_col].
        Typically beta_col is the Cox coefficient (log(HR)).
    gene_col
        Column in `weights` with gene names.
    beta_col
        Column in `weights` with coefficients (log-HR).
    layer
        Expression layer to use. If None or missing, falls back to adata.X.
    standardize
        If True, z-score each gene across samples before scoring.
        Recommended if genes are on different scales.
    missing
        What to do if some genes in `weights` are not in adata.var_names:
        - "drop": silently drop them
        - "error": raise an error
    key_added
        Store the resulting score in `adata.obs[key_added]`.

    Returns
    -------
    score : np.ndarray
        Signature score per sample (length = adata.n_obs).
    """
    if not isinstance(weights, pd.DataFrame):
        weights = pd.DataFrame(weights)

    for c in (gene_col, beta_col):
        if c not in weights.columns:
            raise KeyError(f"weights is missing required column '{c}'")

    w = weights[[gene_col, beta_col]].copy()
    w[gene_col] = w[gene_col].astype(str)
    w = w.dropna(subset=[gene_col, beta_col])

    # keep only genes present
    present = w[gene_col].isin(adata.var_names)
    if not present.all():
        missing_genes = w.loc[~present, gene_col].tolist()
        if missing == "error":
            raise KeyError(f"{len(missing_genes)} signature genes not in adata.var_names: {missing_genes[:10]}")
        w = w.loc[present].copy()

    if w.shape[0] == 0:
        raise ValueError("No signature genes remain after matching to adata.var_names.")

    genes = w[gene_col].tolist()
    beta = w[beta_col].to_numpy(dtype=float)

    # expression matrix: samples x genes
    X = _get_matrix(adata, layer=layer, use="samples")  # must return dense or sparse
    gidx = adata.var_names.get_indexer(genes)
    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # standardize per gene
    if standardize:
        mu = np.nanmean(M, axis=0, keepdims=True)
        sd = np.nanstd(M, axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        M = (M - mu) / sd

    score = M @ beta

    adata.obs[key_added] = pd.to_numeric(score, errors="coerce")
    info(f"Stored signature score in adata.obs['{key_added}'] (n_genes={len(genes)})")

    return np.asarray(score, dtype=float)



def c_index_harrell(time: np.ndarray, event: np.ndarray, score: np.ndarray) -> float:
    """
    Harrell's C-index (higher score = higher risk).
    """
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    s = np.asarray(score, float)

    ok = np.isfinite(t) & np.isfinite(s) & np.isfinite(e)
    t, e, s = t[ok], e[ok], s[ok]

    n = 0
    n_conc = 0.0
    for i in range(len(t)):
        for j in range(len(t)):
            if t[i] < t[j] and e[i] == 1:
                n += 1
                if s[i] > s[j]:
                    n_conc += 1
                elif s[i] == s[j]:
                    n_conc += 0.5
    return float(n_conc / n) if n > 0 else float("nan")


def validate_signature_external(
    adata_train: ad.AnnData,
    adata_ext: ad.AnnData,
    *,
    weights: pd.DataFrame,          # gene,beta
    time_col: str,
    event_col: str,
    layer: str | None = "log1p_cpm",
    score_key: str = "signature_score",
    use_training_threshold: bool = True,
) -> dict:
    """
    Apply a trained signature to an external cohort and compute C-index + threshold.

    Returns a dict with:
      - c_index_ext
      - threshold_used
      - n_ext
    """
    # compute train score (mainly to get threshold)
    signature_score(adata_train, weights=weights, layer=layer, key_added=score_key)
    signature_score(adata_ext, weights=weights, layer=layer, key_added=score_key)

    t_tr = pd.to_numeric(adata_train.obs[time_col], errors="coerce").to_numpy(float)
    e_tr = pd.to_numeric(adata_train.obs[event_col], errors="coerce").to_numpy(float).astype(int)
    s_tr = pd.to_numeric(adata_train.obs[score_key], errors="coerce").to_numpy(float)

    t_ex = pd.to_numeric(adata_ext.obs[time_col], errors="coerce").to_numpy(float)
    e_ex = pd.to_numeric(adata_ext.obs[event_col], errors="coerce").to_numpy(float).astype(int)
    s_ex = pd.to_numeric(adata_ext.obs[score_key], errors="coerce").to_numpy(float)

    thr = float(np.nanmedian(s_tr)) if use_training_threshold else float(np.nanmedian(s_ex))
    c_ext = c_index_harrell(t_ex, e_ex, s_ex)

    return {"c_index_ext": float(c_ext), "threshold_used": thr, "n_ext": int(np.isfinite(s_ex).sum())}


####################################################################
## IMPROVE SIGNATURES
####################################################################

############################################
## HELPERS
############################################

def _z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - mu) / sd


def _clip_logq(q: np.ndarray, cap: float = 10.0) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    q = np.where(q <= 0, np.nan, q)
    s = -np.log10(q)
    s = np.clip(s, 0.0, float(cap))
    s = np.where(np.isfinite(s), s, 0.0)
    return s


def _sign(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    out[x > 0] = 1.0
    out[x < 0] = -1.0
    return out


@dataclass
class DesignedSignature:
    name: str
    genes: list[str]
    weights: pd.DataFrame  # gene, beta
    meta_table: pd.DataFrame

def _sigmoid(z):
    z = np.asarray(z, dtype=float)
    # numeric stability
    z = np.clip(z, -50, 50)
    return 1.0 / (1.0 + np.exp(-z))

def _encode_binary_from_obs(obs_series, positive_label):
    y = pd.Series(obs_series).astype("object")
    ok = y.notna()
    y = y[ok].astype(str)
    ybin = (y == str(positive_label)).astype(int).to_numpy()
    if len(np.unique(ybin)) < 2:
        raise ValueError("After filtering, label has <2 classes.")
    return ok.to_numpy(), ybin



############################################
## FUNCTIONS
############################################

def design_signatures(
    df_gene: pd.DataFrame,
    *,
    gene_col: str = "gene",
    # Cox
    beta_col: str | None = None,          # e.g. "beta" or "coef"
    hr_col: str | None = "HR",            # optional if beta_col missing
    q_cox: str | None = "qval_cox",
    # Corr
    corr_col: str | None = "corr_r",
    q_corr: str | None = "qval_corr",
    # Categorical
    cat_effect_col: str | None = None,    # e.g. "delta_mean" (optional)
    q_cat: str | None = "qval_cat",
    # Optional: stability/permutation
    df_stab: pd.DataFrame | None = None,  # expects: gene,freq,sign_consistency,mean_coef(optional)
    df_perm: pd.DataFrame | None = None,  # expects: gene,delta_cindex or delta_auc
    perm_col: str | None = "delta_cindex",
    # scoring weights
    cap_logq: float = 10.0,
    w_cox: float = 1.0,
    w_cat: float = 0.6,
    w_corr: float = 0.4,
    w_stab: float = 1.0,
    w_imp: float = 0.8,
    # design choices
    panel_sizes: Sequence[int] = (12, 20, 36),
    direction_anchor: Literal["cox", "corr", "majority"] = "cox",
    min_dir_agreement: float = 0.5,  # below this, downweight heavily
    name_prefix: str = "sig_designed",
    # weight estimation for the final signature
    weight_mode: Literal["meta", "cox_beta"] = "cox_beta",
) -> tuple[pd.DataFrame, dict[int, DesignedSignature]]:
    """
    Build "designed signatures" from pooled candidate genes + evidence tables.

    Returns
    -------
    df_ranked : pd.DataFrame
        gene-level table with meta_score, direction, support components.
    sigs : dict[int, DesignedSignature]
        one signature per requested panel size K.
    """
    df = df_gene.copy()
    if gene_col not in df.columns:
        raise KeyError(f"df_gene missing '{gene_col}'")
    df[gene_col] = df[gene_col].astype(str)

    # --- determine beta (direction & weights) ---
    beta = None
    if beta_col is not None and beta_col in df.columns:
        beta = pd.to_numeric(df[beta_col], errors="coerce").to_numpy(float)
    elif hr_col is not None and hr_col in df.columns:
        hr = pd.to_numeric(df[hr_col], errors="coerce").to_numpy(float)
        beta = np.log(hr)
    else:
        beta = np.full(df.shape[0], np.nan, dtype=float)

    # --- supports (evidence strength) ---
    S_cox = _clip_logq(pd.to_numeric(df[q_cox], errors="coerce").to_numpy(float), cap=cap_logq) if (q_cox and q_cox in df.columns) else 0.0
    S_cat = _clip_logq(pd.to_numeric(df[q_cat], errors="coerce").to_numpy(float), cap=cap_logq) if (q_cat and q_cat in df.columns) else 0.0
    S_corr = _clip_logq(pd.to_numeric(df[q_corr], errors="coerce").to_numpy(float), cap=cap_logq) if (q_corr and q_corr in df.columns) else 0.0

    if not isinstance(S_cox, np.ndarray): S_cox = np.zeros(df.shape[0])
    if not isinstance(S_cat, np.ndarray): S_cat = np.zeros(df.shape[0])
    if not isinstance(S_corr, np.ndarray): S_corr = np.zeros(df.shape[0])

    E = (w_cox * S_cox) + (w_cat * S_cat) + (w_corr * S_corr)
    E_z = _z(E)

    # --- direction signals ---
    dir_cox = _sign(beta) if np.isfinite(beta).any() else np.zeros(df.shape[0])
    dir_corr = _sign(pd.to_numeric(df[corr_col], errors="coerce").to_numpy(float)) if (corr_col and corr_col in df.columns) else np.zeros(df.shape[0])
    dir_cat = _sign(pd.to_numeric(df[cat_effect_col], errors="coerce").to_numpy(float)) if (cat_effect_col and cat_effect_col in df.columns) else np.zeros(df.shape[0])

    # pick anchor sign
    if direction_anchor == "cox":
        anchor = dir_cox
        # if beta missing, fallback to corr
        anchor = np.where(anchor == 0, dir_corr, anchor)
    elif direction_anchor == "corr":
        anchor = dir_corr
        anchor = np.where(anchor == 0, dir_cox, anchor)
    else:  # majority
        anchor = _sign(dir_cox + dir_corr + dir_cat)

    # agreement fraction across available dirs
    dirs = []
    if np.any(dir_cox != 0): dirs.append(dir_cox)
    if np.any(dir_corr != 0): dirs.append(dir_corr)
    if np.any(dir_cat != 0): dirs.append(dir_cat)

    if len(dirs) == 0:
        P = np.ones(df.shape[0], dtype=float)
    else:
        A = np.vstack(dirs)  # k x n
        # treat 0 as missing; only compare where non-zero
        agree = np.zeros(df.shape[0], dtype=float)
        denom = np.zeros(df.shape[0], dtype=float)
        for row in A:
            m = row != 0
            denom += m.astype(float)
            agree += ((row == anchor) & m).astype(float)
        denom = np.maximum(denom, 1.0)
        P = agree / denom
        # downweight poor agreement
        P = np.where(P < float(min_dir_agreement), P * 0.25, P)

    # --- stability ---
    Stab_z = np.zeros(df.shape[0], dtype=float)
    if df_stab is not None and df_stab.shape[0] > 0:
        ds = df_stab.copy()
        if "gene" not in ds.columns:
            raise KeyError("df_stab must contain column 'gene'")
        ds["gene"] = ds["gene"].astype(str)
        ds = ds.set_index("gene")
        freq = ds["freq"] if "freq" in ds.columns else None
        sc = ds["sign_consistency"] if "sign_consistency" in ds.columns else None
        stab = None
        if freq is not None and sc is not None:
            stab = (freq * sc).astype(float)
        elif freq is not None:
            stab = freq.astype(float)
        if stab is not None:
            stab = df[gene_col].map(stab).to_numpy(dtype=float)
            Stab_z = _z(np.nan_to_num(stab, nan=0.0))

    # --- permutation importance ---
    Imp_z = np.zeros(df.shape[0], dtype=float)
    if df_perm is not None and df_perm.shape[0] > 0 and perm_col is not None:
        dp = df_perm.copy()
        if "gene" not in dp.columns:
            raise KeyError("df_perm must contain column 'gene'")
        if perm_col not in dp.columns:
            raise KeyError(f"df_perm must contain column '{perm_col}'")
        dp["gene"] = dp["gene"].astype(str)
        dp = dp.set_index("gene")
        imp = df[gene_col].map(dp[perm_col]).to_numpy(dtype=float)
        Imp_z = _z(np.nan_to_num(imp, nan=0.0))

    # --- meta score ---
    meta = E_z + (w_stab * Stab_z) + (w_imp * Imp_z)
    meta = meta * P

    df_out = df.copy()
    df_out["beta"] = beta
    df_out["dir"] = anchor
    df_out["support_cox"] = S_cox
    df_out["support_cat"] = S_cat
    df_out["support_corr"] = S_corr
    df_out["evidence"] = E
    df_out["dir_agreement"] = P
    df_out["stab_z"] = Stab_z
    df_out["imp_z"] = Imp_z
    df_out["meta_score"] = meta
    df_out = df_out.sort_values("meta_score", ascending=False).reset_index(drop=True)

    # --- build signatures ---
    sigs: dict[int, DesignedSignature] = {}
    for K in panel_sizes:
        K = int(K)
        dfK = df_out.head(K).copy()
        genes = dfK[gene_col].astype(str).tolist()

        if weight_mode == "cox_beta":
            w = dfK[["gene", "beta"]].copy()
            # enforce direction exists
            w["beta"] = np.where(np.isfinite(w["beta"]), w["beta"], 0.0)
            # align sign to anchor direction if needed
            w["beta"] = np.where(w["beta"] == 0.0, dfK["dir"].to_numpy(float), w["beta"])
        else:
            # meta weights: direction * meta_score (then normalized)
            ww = (dfK["dir"].to_numpy(float) * dfK["meta_score"].to_numpy(float))
            # scale for readability
            if np.nanstd(ww) > 0:
                ww = ww / np.nanstd(ww)
            w = pd.DataFrame({"gene": genes, "beta": ww})

        name = f"{name_prefix}_{K}"
        sigs[K] = DesignedSignature(name=name, genes=genes, weights=w, meta_table=dfK)

    return df_out, sigs


def signature_permutation_importance(
    adata,
    *,
    weights: pd.DataFrame,
    layer: str | None = "log1p_cpm",
    metric: str = "cindex",  # "cindex" | "auc"
    # survival:
    time_col: str | None = None,
    event_col: str | None = None,
    # auc:
    label_col: str | None = None,
    positive_label=None,
    intercept: float = 0.0,      # NEW: for AUC-native probability
    use_proba: bool = True,      # NEW: if True -> sigmoid(intercept+linear_score) for AUC
    n_repeats: int = 50,
    seed: int = 0,
    score_key: str | None = None,
    key_added: str | None = None,
) -> pd.DataFrame:
    """
    Permutation importance for a linear gene signature.

    For metric="cindex": delta_cindex = cindex_base - cindex_perm
    For metric="auc":    delta_auc    = auc_base    - auc_perm
      (AUC uses probability by default: sigmoid(intercept + linear_score))

    weights: DataFrame with columns ["gene","beta"] at minimum.
    """
    metric = str(metric).lower()
    if metric not in {"cindex", "auc"}:
        raise ValueError("metric must be 'cindex' or 'auc'")

    if "gene" not in weights.columns or "beta" not in weights.columns:
        raise ValueError("weights must contain columns: gene, beta")

    genes = weights["gene"].astype(str).tolist()
    betas = weights["beta"].to_numpy(dtype=float)

    # compute baseline linear score (in adata.obs)
    if score_key is None:
        score_key = f"_sig_score_tmp_{metric}"
    # signature_score should ignore missing genes internally; if not, filter here.
    from .signature import signature_score  # if same file, remove this import and call directly
    signature_score(adata, weights=weights, layer=layer, key_added=score_key)

    s_all = pd.to_numeric(adata.obs[score_key], errors="coerce").to_numpy(float)

    rng = np.random.default_rng(int(seed))

    # prepare y / baseline metric
    if metric == "auc":
        if label_col is None:
            raise ValueError("label_col is required for metric='auc'")
        if positive_label is None:
            raise ValueError("positive_label is required for metric='auc' (e.g. 'R')")

        ok_y, y = _encode_binary_from_obs(adata.obs[label_col], positive_label=positive_label)

        ok = ok_y & np.isfinite(s_all)
        y0 = y[ok[ok_y]]  # careful alignment: ok_y filters first
        s0 = s_all[ok]

        # simpler: rebuild aligned arrays
        idx = np.where(ok)[0]
        y0 = (pd.Series(adata.obs[label_col].iloc[idx]).astype(str) == str(positive_label)).astype(int).to_numpy()
        s0 = s_all[idx]

        if s0.size < 10:
            raise ValueError("Not enough valid samples for permutation importance (auc).")
        if len(np.unique(y0)) < 2:
            raise ValueError("label_col has <2 classes after filtering; cannot compute AUC.")

        try:
            from sklearn.metrics import roc_auc_score
        except Exception as e:
            raise ImportError(f"AUC requires scikit-learn. ({e})")

        p_base = _sigmoid(float(intercept) + s0) if use_proba else s0
        auc_base = float(roc_auc_score(y0, p_base))

    else:
        # cindex
        if time_col is None or event_col is None:
            raise ValueError("time_col and event_col are required for metric='cindex'")
        t = pd.to_numeric(adata.obs[time_col], errors="coerce").to_numpy(float)
        e = pd.to_numeric(adata.obs[event_col], errors="coerce").to_numpy(float).astype(int)

        ok = np.isfinite(t) & np.isfinite(e) & np.isfinite(s_all)
        t0, e0, s0 = t[ok], e[ok], s_all[ok]
        if s0.size < 10:
            raise ValueError("Not enough valid samples for permutation importance (cindex).")

        # use your internal c-index helper if you have it
        from .cox import _c_index_from_risk  # adjust import if needed
        cindex_base = float(_c_index_from_risk(t0, e0, s0))

    # expression matrix for permuting one gene at a time
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    gidx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    genes_used = [g for g in genes if g in gene_to_idx]
    betas_used = np.array([weights.loc[weights["gene"].astype(str) == g, "beta"].values[0] for g in genes_used], dtype=float)

    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # align to the same subset used for baseline metric
    # baseline index "idx" for auc; "ok" mask for cindex
    if metric == "auc":
        M0 = M[idx, :]
        s_base = s0
        # sanity: s0 should equal linear combo; but we rely on signature_score for baseline
    else:
        M0 = M[ok, :]
        s_base = s0

    # compute each gene’s contribution in the linear score:
    # score = sum_j beta_j * x_j ; permuting x_g means replace that term.
    # Need per-sample x_g for each gene.
    rows = []
    for j, g in enumerate(genes_used):
        xg = M0[:, j].astype(float)
        bg = float(betas_used[j])

        deltas = []
        for _ in range(int(n_repeats)):
            xp = xg.copy()
            rng.shuffle(xp)
            # new score = base_score - beta*xg + beta*xp
            s_perm = s_base - bg * xg + bg * xp

            if metric == "auc":
                p_perm = _sigmoid(float(intercept) + s_perm) if use_proba else s_perm
                auc_perm = float(roc_auc_score(y0, p_perm))
                deltas.append(auc_base - auc_perm)
            else:
                c_perm = float(_c_index_from_risk(t0, e0, s_perm))
                deltas.append(cindex_base - c_perm)

        deltas = np.asarray(deltas, dtype=float)
        rows.append(
            {
                "gene": g,
                "metric": metric,
                "n_repeats": int(n_repeats),
                "delta_mean": float(np.mean(deltas)),
                "delta_sd": float(np.std(deltas)),
                "delta_values": deltas,  # you can drop this if you want lighter output
            }
        )

    df_perm = pd.DataFrame(rows).sort_values("delta_mean", ascending=False).reset_index(drop=True)

    if key_added is not None:
        adata.uns.setdefault("signature_permutation_importance", {})
        adata.uns["signature_permutation_importance"][key_added] = df_perm

    return df_perm



def correlation_redundancy_filter(
    adata,
    *,
    ranked_genes: Sequence[str],
    layer: str | None = "log1p_cpm",
    max_corr: float = 0.8,
    K: int | None = None,
    corr_method: Literal["pearson", "spearman"] = "pearson",
    use_abs: bool = True,
    min_var: float = 1e-12,
    return_report: bool = False,
):
    """
    Greedy redundancy filter for a ranked gene list.

    Algorithm:
      - iterate genes in ranked order
      - keep g if max(|corr(g, s)| for s in selected) < max_corr
      - stop when K genes selected (if K is not None)

    Parameters
    ----------
    adata : AnnData-like
    ranked_genes : list of genes in priority order (best first)
    layer : expression layer to use; if None or missing uses adata.X
    max_corr : correlation threshold (e.g. 0.7 or 0.8)
    K : target number of genes to keep (optional)
    corr_method : "pearson" or "spearman"
    use_abs : if True uses absolute correlation
    min_var : skip genes with variance <= min_var
    return_report : if True returns a DataFrame explaining rejections

    Returns
    -------
    selected : list[str]
    report : pd.DataFrame (optional) with per-gene decision + max_corr_to_selected
    """
    if not (0.0 <= float(max_corr) <= 1.0):
        raise ValueError("max_corr must be within [0, 1].")

    # --- pick expression matrix
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X

    # genes present and in order
    ranked = [str(g) for g in ranked_genes if str(g) in adata.var_names]
    if len(ranked) == 0:
        raise ValueError("None of ranked_genes are in adata.var_names.")

    # index lookup
    idx = adata.var_names.get_indexer(ranked)

    # extract samples x genes
    M = X[:, idx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)

    # optional: drop genes with ~0 variance (stabilizes corr)
    v = np.nanvar(M, axis=0)
    keep = v > float(min_var)
    ranked = [g for g, ok in zip(ranked, keep) if ok]
    M = M[:, keep]

    if M.shape[1] == 0:
        raise ValueError("All genes were filtered out by min_var.")

    # compute correlation matrix for all ranked genes once (fast enough for <= a few thousand genes)
    if corr_method == "spearman":
        # rank-transform each column (ties average)
        M_rank = np.apply_along_axis(lambda x: pd.Series(x).rank(method="average").to_numpy(), 0, M)
        C = np.corrcoef(M_rank, rowvar=False)
    elif corr_method == "pearson":
        C = np.corrcoef(M, rowvar=False)
    else:
        raise ValueError("corr_method must be 'pearson' or 'spearman'.")

    # numerical safety
    C = np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(C, 0.0)  # ignore self-corr

    if use_abs:
        C = np.abs(C)

    # greedy selection
    selected_idx: list[int] = []
    report_rows = []

    for j, g in enumerate(ranked):
        if K is not None and len(selected_idx) >= int(K):
            break

        if len(selected_idx) == 0:
            selected_idx.append(j)
            if return_report:
                report_rows.append(
                    {"gene": g, "keep": True, "max_corr_to_selected": 0.0, "most_correlated_with": None}
                )
            continue

        corrs = C[j, selected_idx]
        mx = float(np.max(corrs)) if corrs.size else 0.0
        if mx < float(max_corr):
            selected_idx.append(j)
            keep_flag = True
            most = None
        else:
            keep_flag = False
            # who caused the rejection?
            kbest = int(selected_idx[int(np.argmax(corrs))])
            most = ranked[kbest]

        if return_report:
            report_rows.append(
                {"gene": g, "keep": keep_flag, "max_corr_to_selected": mx, "most_correlated_with": most}
            )

    selected = [ranked[j] for j in selected_idx]

    if return_report:
        report = pd.DataFrame(report_rows)
        return selected, report
    return selected



def meta_rank_genes_from_signatures(
    df_sig_perf: pd.DataFrame,
    df_bin: pd.DataFrame,
    *,
    sig_col: str = "signature",
    metric_col: str = "auc",
    hr_col: str | None = "hr",
    pval_col: str | None = "pval",
    top_n: int = 10,
    w_auc: float = 1.0,
    w_hr: float = 0.2,
    w_p: float = 0.2,
) -> pd.DataFrame:
    """
    Rank genes by how strongly they are supported by high-performing signatures.

    df_sig_perf: one row per signature with external metrics
    df_bin: gene x signature membership (0/1), columns must match signature ids
    """
    perf = df_sig_perf.copy()
    if sig_col not in perf.columns:
        raise KeyError(f"'{sig_col}' not in df_sig_perf")
    if metric_col not in perf.columns:
        raise KeyError(f"'{metric_col}' not in df_sig_perf")

    # align signatures
    perf = perf.set_index(sig_col)
    common = [s for s in perf.index if s in df_bin.columns]
    if len(common) < 3:
        raise ValueError("Too few overlapping signatures between df_sig_perf and df_bin.")

    perf = perf.loc[common]
    X = df_bin[common].astype(float)  # gene x sig

    # weights from metrics (z-scored for stability)
    def _z(x):
        x = np.asarray(x, float)
        sd = np.nanstd(x)
        return (x - np.nanmean(x)) / (sd if sd > 1e-12 else 1.0)

    w = w_auc * _z(perf[metric_col].to_numpy())

    if hr_col is not None and hr_col in perf.columns:
        # encourage stronger separation, use |log(HR)|
        w = w + w_hr * _z(np.abs(np.log(perf[hr_col].to_numpy(dtype=float))))

    if pval_col is not None and pval_col in perf.columns:
        w = w + w_p * _z(-np.log10(np.clip(perf[pval_col].to_numpy(dtype=float), 1e-300, 1.0)))

    w = pd.Series(w, index=common, name="sig_weight")

    # gene meta score = sum membership * signature weight
    meta = (X * w).sum(axis=1)

    # frequency info
    freq_all = X.mean(axis=1)
    top_sigs = perf[metric_col].sort_values(ascending=False).head(int(top_n)).index
    freq_top = X[top_sigs].mean(axis=1) if len(top_sigs) else np.nan

    out = pd.DataFrame({
        "gene": X.index.astype(str),
        "meta_score": meta.to_numpy(dtype=float),
        "freq_all": freq_all.to_numpy(dtype=float),
        f"freq_top{top_n}": np.asarray(freq_top, float),
        "n_signatures": X.sum(axis=1).to_numpy(dtype=int),
    }).sort_values(["meta_score", f"freq_top{top_n}", "n_signatures"], ascending=False).reset_index(drop=True)

    return out


def redundancy_filter_by_corr(
    adata,
    ranked_genes: list[str],
    *,
    layer: str | None = "log1p_cpm",
    method: str = "pearson",
    max_corr: float = 0.80,
    use_abs: bool = True,
    max_genes: int | None = None,
) -> list[str]:
    """
    Greedy redundancy filter on expression correlation.
    Keeps a gene only if it is not too correlated with already-selected genes.
    """
    # expression getter
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X

    sel: list[str] = []
    for g in ranked_genes:
        if g not in adata.var_names:
            continue
        if len(sel) == 0:
            sel.append(g)
        else:
            # compute corr(g, each selected)
            gidx = int(adata.var_names.get_loc(g))
            sidx = [int(adata.var_names.get_loc(x)) for x in sel if x in adata.var_names]
            v = X[:, gidx].toarray().ravel() if hasattr(X[:, gidx], "toarray") else np.asarray(X[:, gidx]).ravel()
            S = X[:, sidx].toarray() if hasattr(X[:, sidx], "toarray") else np.asarray(X[:, sidx])

            # correlation with each selected gene
            if method == "spearman":
                v = pd.Series(v).rank().to_numpy()
                S = np.apply_along_axis(lambda x: pd.Series(x).rank().to_numpy(), 0, S)

            # corr with columns
            vv = v - np.nanmean(v)
            denom_v = np.sqrt(np.nanmean(vv**2)) + 1e-12
            ok = []
            for j in range(S.shape[1]):
                sj = S[:, j]
                sj = sj - np.nanmean(sj)
                denom_s = np.sqrt(np.nanmean(sj**2)) + 1e-12
                c = float(np.nanmean(vv * sj) / (denom_v * denom_s))
                ok.append(abs(c) if use_abs else c)

            if np.nanmax(ok) < float(max_corr):
                sel.append(g)

        if max_genes is not None and len(sel) >= int(max_genes):
            break

    return sel

#############################################################
## AUC
#############################################################


###########################
## HELPERS
###########################

# If you already have a matrix getter in your codebase, prefer that.
# This one is safe + minimal.
def _get_expr_matrix(adata, genes: Sequence[str], layer: str | None = "log1p_cpm"):
    genes = [str(g) for g in genes if str(g) in adata.var_names]
    if len(genes) == 0:
        raise ValueError("No provided genes found in adata.var_names.")
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    gidx = adata.var_names.get_indexer(genes)
    M = X[:, gidx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, dtype=float)
    return genes, np.asarray(M, dtype=float)

def _standardize_cols(M: np.ndarray) -> np.ndarray:
    mu = np.nanmean(M, axis=0, keepdims=True)
    sd = np.nanstd(M, axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (M - mu) / sd

def _encode_binary_labels(y, positive_label):
    # y can be strings ("R","NR") or numbers
    ys = pd.Series(y).astype("object")
    ok = ys.notna()
    ys = ys[ok]
    if positive_label is None:
        # choose max label as positive for numeric; for strings choose most frequent? -> force user to set it ideally
        uniq = list(pd.unique(ys))
        if len(uniq) != 2:
            raise ValueError("positive_label must be provided when label_col is non-binary or ambiguous.")
        positive_label = uniq[0]
    ybin = (ys.astype(str) == str(positive_label)).astype(int).to_numpy()
    if len(np.unique(ybin)) < 2:
        raise ValueError("Binary labels have <2 classes after filtering.")
    return ok.to_numpy(), ybin

def _corr_redundancy_filter_ranked(
    M: np.ndarray,
    genes: list[str],
    *,
    corr_threshold: float = 0.8,
    method: Literal["pearson"] = "pearson",
    use_abs: bool = True,
):
    """
    Greedy redundancy filter along the provided order of `genes`:
    keep gene i if corr(gene_i, any_kept) < threshold.
    """
    if len(genes) <= 1:
        return genes

    kept_idx: list[int] = []
    for i in range(len(genes)):
        if len(kept_idx) == 0:
            kept_idx.append(i)
            continue

        xi = M[:, i]
        # compute corr to already kept genes
        xk = M[:, kept_idx]  # n x k
        # pearson corr: corr(xi, xk_j)
        xi0 = xi - np.nanmean(xi)
        den_i = np.nanstd(xi0)
        if den_i == 0 or not np.isfinite(den_i):
            # constant -> redundant/noisy -> drop
            continue

        xk0 = xk - np.nanmean(xk, axis=0, keepdims=True)
        den_k = np.nanstd(xk0, axis=0)
        den_k[den_k == 0] = np.nan

        corr = (np.nanmean(xi0[:, None] * xk0, axis=0) / (den_i * den_k))
        corr = np.nan_to_num(corr, nan=0.0)
        if use_abs:
            corr = np.abs(corr)

        if float(np.nanmax(corr)) < float(corr_threshold):
            kept_idx.append(i)

    return [genes[i] for i in kept_idx]

def _to_binary_labels(y: pd.Series, positive_label):
    """
    Convert labels to 0/1.
    positive_label can be 'R', 'NR', 1, True, etc.
    """
    y = pd.Series(y).astype(str) if y.dtype == object else pd.Series(y)
    if positive_label is None:
        raise ValueError("positive_label must be provided for binary tasks.")
    y01 = (y == positive_label).astype(int).to_numpy()
    return y01


def _score_binary_metric(y01: np.ndarray, prob: np.ndarray, metric: str) -> float:
    """
    metric: 'roc_auc' or 'pr_auc'
    prob must be P(y=1)
    """
    if metric == "roc_auc":
        return float(roc_auc_score(y01, prob))
    if metric in {"pr_auc", "ap"}:
        return float(average_precision_score(y01, prob))
    raise ValueError("metric must be 'roc_auc' or 'pr_auc' (or 'ap').")



###########################
## FUNCTIONS
###########################

def recommended_auc_panel(
    adata,
    *,
    genes: Sequence[str],
    label_col: str = "PFS_6m",
    positive_label: str = "R",
    layer: str | None = "log1p_cpm",
    standardize: bool = True,

    metric: Literal["roc_auc", "pr_auc"] = "pr_auc",   # NEW

    Cs: Sequence[float] = tuple(np.logspace(-4, -1, 10)),
    l1_ratios: Sequence[float] = (0.8, 1.0),
    n_splits: int = 3,
    n_repeats: int = 10,
    seed: int = 0,

    n_boot: int = 200,
    sample_fraction: float = 0.8,
    freq_threshold: float = 0.6,
    sign_threshold: float = 0.8,
    max_genes: int = 12,


    preselect_top: int | None = None,     # NEW: reduce p before CV
    max_iter: int = 50000,              # NEW
    tol: float = 1e-3,                  # NEW
    refit_max_tries: int = 4,           # NEW
    refit_C_shrink: float = 10.0,       # NEW

    min_var: float | None = None,   # NEW: if None, skip variance filter entirely
    top_var: int | None = None,     # NEW: if set, keep top_var genes by variance (after min_var)

    redundancy_filter: bool = True,
    corr_threshold: float = 0.8,

    min_per_class: int | None = None,  # CHANGED (was 8)
    coef_eps: float = 1e-10,

    key_added: str | None = None,
) -> dict:
    """
    One-call helper to design an AUC-optimized signature for a binary label (e.g. PFS_6m: R/NR).

    Steps
    -----
    1) (Optional) redundancy filter by expression correlation (greedy, respects gene order)
    2) CV grid search of LogisticRegression(elastic-net) by mean ROC AUC
    3) Stability selection (bootstrap/subsample) using best hyperparams
    4) Pick compact panel by stability freq (+ optional sign consistency)
    5) Fit final model on full data -> weights_df (gene, beta) + intercept

    Returns dict with:
      best_params, df_cv, df_stability, panel, weights_df, intercept, auc_cv_best
    """
    # sklearn dependency
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import roc_auc_score, average_precision_score
    except Exception as e:
        raise ImportError(f"recommended_auc_panel requires scikit-learn. ({e})")

    if label_col not in adata.obs.columns:
        raise KeyError(f"label_col='{label_col}' not found in adata.obs")

    # ---- prepare X / y ----
    genes_in, M = _get_expr_matrix(adata, genes, layer=layer)
    ok_y, ybin = _encode_binary_labels(adata.obs[label_col].to_numpy(), positive_label=positive_label)

    M = M[ok_y, :]
    genes_in = list(genes_in)

    # drop rows with NaNs in X
    ok_x = np.all(np.isfinite(M), axis=1)
    M = M[ok_x, :]
    ybin = ybin[ok_x]

    # class balance
    n0 = int((ybin == 0).sum())
    n1 = int((ybin == 1).sum())
    min_class = min(n0, n1)

    if min_per_class is None:
        min_per_class_eff = 2 if metric == "pr_auc" else 8
    else:
        min_per_class_eff = int(min_per_class)

    if min_class < min_per_class_eff:
        raise ValueError(
            f"Not enough samples per class after filtering: n0={n0}, n1={n1}, "
            f"min_per_class={min_per_class_eff} (metric={metric})."
        )

    # CV folds cannot exceed minority class count
    n_splits_eff = int(min(n_splits, min_class))
    if n_splits_eff < 2:
        raise ValueError(f"Need at least 2 folds, got n_splits_eff={n_splits_eff} (min_class={min_class}).")



    # ---- optional redundancy filter (cheap win for AUC + interpretability) ----
    if redundancy_filter and len(genes_in) > 2:
        genes_nr = _corr_redundancy_filter_ranked(
            M, genes_in, corr_threshold=float(corr_threshold), use_abs=True
        )
        if len(genes_nr) >= 2:
            keep = [genes_in.index(g) for g in genes_nr]
            M = M[:, keep]
            genes_in = genes_nr
            info(f"recommended_auc_panel: redundancy filter kept {len(genes_in)} genes (thr={corr_threshold}).")
        else:
            warn("recommended_auc_panel: redundancy filter removed too many genes; keeping original list.")

    if preselect_top is not None and len(genes_in) > int(preselect_top):
        genes_in = list(genes_in)[: int(preselect_top)]
        M = M[:, : int(preselect_top)]
        info(f"recommended_auc_panel: preselect_top kept {len(genes_in)} genes.")


    # ---- OPTIONAL variance filter (OFF by default) ----
    # Important: PR-AUC ranking is label-aware; variance filtering can remove predictive low-variance markers.
    if (min_var is not None) or (top_var is not None):
        v = np.nanvar(M, axis=0)

        keep = np.ones(M.shape[1], dtype=bool)

        if min_var is not None:
            keep &= (v >= float(min_var))

        if top_var is not None:
            idx = np.where(keep)[0]
            if idx.size > int(top_var):
                top = idx[np.argsort(v[idx])[::-1][:int(top_var)]]
                keep2 = np.zeros_like(keep, dtype=bool)
                keep2[top] = True
                keep = keep2

        # apply
        genes_before = len(genes_in)
        genes_in = [g for g, k in zip(genes_in, keep) if k]
        M = M[:, keep]

        info(
            f"recommended_auc_panel: variance filter kept {len(genes_in)}/{genes_before} genes "
            f"(min_var={min_var}, top_var={top_var})."
        )

        if len(genes_in) < 2:
            raise ValueError(
                "recommended_auc_panel: variance filter removed too many genes. "
                "Disable it (min_var=None, top_var=None) or relax thresholds."
            )


    # ---- CV grid search (AUC) ----
    rng = np.random.default_rng(int(seed))
    Cs = [float(c) for c in Cs]
    l1_ratios = [float(x) for x in l1_ratios]

    info(
        f"recommended_auc_panel: genes={len(genes_in)}, n={M.shape[0]}, "
        f"n0={n0}, n1={n1}, splits={n_splits_eff}, repeats={n_repeats}, grid={len(Cs)*len(l1_ratios)}"
    )

    rows = []
    for C in Cs:
        for l1 in l1_ratios:
            aucs = []
            nfail = 0

            for rep in range(int(n_repeats)):
                # shuffle with new seed per repeat
                skf = StratifiedKFold(
                    n_splits=int(n_splits_eff),
                    shuffle=True,
                    random_state=int(rng.integers(0, 2**31 - 1)),
                )
                for tr, te in skf.split(M, ybin):
                    Xtr, Xte = M[tr], M[te]
                    ytr, yte = ybin[tr], ybin[te]

                    try:
                        clf = LogisticRegression(
                            penalty="elasticnet",
                            solver="saga",
                            l1_ratio=float(l1),
                            C=float(C),
                            max_iter=int(max_iter),
                            tol=float(tol),
                            class_weight="balanced",
                            random_state=int(rng.integers(0, 2**31 - 1)),
                        )

                        clf.fit(Xtr, ytr)
                        if np.max(clf.n_iter_) >= int(0.98 * clf.max_iter):
                            nfail += 1
                            continue

                        p = clf.predict_proba(Xte)[:, 1]

                        # guard: both classes needed for AUC/AP
                        if len(np.unique(yte)) < 2:
                            nfail += 1
                            continue

                        if metric == "roc_auc":
                            score = roc_auc_score(yte, p)
                        elif metric == "pr_auc":
                            score = average_precision_score(yte, p)
                        else:
                            raise ValueError("metric must be 'roc_auc' or 'pr_auc'")

                        aucs.append(float(score))

                    except Exception:
                        nfail += 1
                        continue

            rows.append(
                {
                    "C": float(C),
                    "l1_ratio": float(l1),
                    "score_mean": float(np.nanmean(aucs)) if len(aucs) else np.nan,
                    "score_sd": float(np.nanstd(aucs)) if len(aucs) else np.nan,
                    "n_fails": int(nfail),
                }
            )

    df_cv = pd.DataFrame(rows)

    df_ok = df_cv[np.isfinite(df_cv["score_mean"])].copy()
    if df_ok.shape[0] == 0:
        raise ValueError("recommended_auc_panel: no valid CV fits. Try fewer genes or adjust grid.")

    # best by AUC, tie-break by stronger sparsity (higher l1) then smaller C (more regularization)
    df_ok = df_ok.sort_values(
        ["score_mean", "l1_ratio", "C"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    best = df_ok.iloc[0]
    best_params = {"C": float(best["C"]), "l1_ratio": float(best["l1_ratio"])}
    score_cv_best = float(best["score_mean"])



    # ---- stability selection ----
    p = len(genes_in)
    sel_counts = np.zeros(p, dtype=int)
    coef_sum = np.zeros(p, dtype=float)
    sign_pos = np.zeros(p, dtype=int)
    sign_neg = np.zeros(p, dtype=int)
    n_selected_each = []
    n_ok = 0
    n_fail = 0

    for _ in range(int(n_boot)):


        pos = np.where(ybin == 1)[0]
        neg = np.where(ybin == 0)[0]

        m = int(np.ceil(float(sample_fraction) * M.shape[0]))

        # enforce at least 2 positives per bootstrap (adjust if you want 3)
        min_pos = 2
        n_pos = max(min_pos, int(round(m * (len(pos) / len(ybin)))))
        n_pos = min(n_pos, len(pos))          # cannot exceed available
        n_neg = m - n_pos

        idx_pos = rng.choice(pos, size=n_pos, replace=True)
        idx_neg = rng.choice(neg, size=n_neg, replace=True)
        idx = np.r_[idx_pos, idx_neg]
        rng.shuffle(idx)

        Xb = M[idx, :]
        yb = ybin[idx]


        try:
            clf = LogisticRegression(
                penalty="elasticnet",
                solver="saga",
                l1_ratio=float(best_params["l1_ratio"]),
                C=float(best_params["C"]),
                max_iter=int(max_iter),
                tol=float(tol),
                class_weight="balanced",
                random_state=int(rng.integers(0, 2**31 - 1)),
            )
            clf.fit(Xb, yb)

            if np.max(clf.n_iter_) >= clf.max_iter:
                n_fail += 1
                continue

            beta = clf.coef_.ravel().astype(float)

            chosen = np.abs(beta) > float(coef_eps)
            sel_counts += chosen.astype(int)
            coef_sum += np.where(chosen, beta, 0.0)
            sign_pos += (chosen & (beta > 0)).astype(int)
            sign_neg += (chosen & (beta < 0)).astype(int)

            n_selected_each.append(int(chosen.sum()))
            n_ok += 1
        except Exception:
            n_fail += 1
            continue

    if n_ok == 0:
        raise ValueError("recommended_auc_panel: no successful stability fits. Try fewer genes or different grid.")

    freq = sel_counts / float(n_ok)
    mean_coef = np.where(sel_counts > 0, coef_sum / np.maximum(sel_counts, 1), 0.0)

    sign_cons = np.full(p, np.nan, dtype=float)
    for j in range(p):
        if sel_counts[j] == 0:
            continue
        sign_cons[j] = float(max(sign_pos[j], sign_neg[j]) / sel_counts[j])

    df_stab = pd.DataFrame(
        {
            "gene": genes_in,
            "freq": freq,
            "mean_coef": mean_coef,
            "sign_consistency": sign_cons,
            "n_runs": int(n_ok),
            "n_selected_mean": float(np.mean(n_selected_each)) if n_selected_each else np.nan,
            "n_selected_median": float(np.median(n_selected_each)) if n_selected_each else np.nan,
        }
    ).sort_values(["freq", "sign_consistency", "gene"], ascending=[False, False, True]).reset_index(drop=True)

    if n_fail > 0:
        warn(f"recommended_auc_panel: {n_fail} bootstrap fits skipped/failed (single-class or convergence).")

    # ---- pick panel ----
    df_pick = df_stab.copy()
    df_pick = df_pick[np.isfinite(df_pick["freq"])].copy()
    df_pick = df_pick[df_pick["freq"] >= float(freq_threshold)].copy()
    df_pick = df_pick[df_pick["sign_consistency"].isna() | (df_pick["sign_consistency"] >= float(sign_threshold))].copy()

    if df_pick.shape[0] == 0:
        warn("recommended_auc_panel: no genes pass thresholds; falling back to top genes by stability frequency.")
        df_pick = df_stab.head(int(max_genes)).copy()

    df_pick = df_pick.sort_values(["freq", "sign_consistency"], ascending=[False, False], na_position="last")
    panel = df_pick["gene"].astype(str).head(int(max_genes)).tolist()

    # ---- final fit on full data (panel only) ----
    # ---- final fit on full data (panel only) ----
    df_pick = df_pick.reset_index(drop=True)  # ensure clean indexing
    ranked_candidates = df_pick["gene"].astype(str).tolist()

    panel = ranked_candidates[: int(max_genes)]
    panel = list(dict.fromkeys(panel))  # preserve order, unique


    # ---- final fit on full data (panel only) ----
    df_pick = df_pick.reset_index(drop=True)
    ranked_candidates = df_pick["gene"].astype(str).tolist()

    panel = ranked_candidates[: int(max_genes)]
    panel = list(dict.fromkeys(panel))  # preserve order, unique

    def _fit_final(panel_genes, C_value):
        # panel_genes defines which columns to use
        keep = [genes_in.index(g) for g in panel_genes]
        Xf = M[:, keep]

        clf = LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            l1_ratio=float(best_params["l1_ratio"]),
            C=float(C_value),
            max_iter=int(max_iter),
            tol=float(tol),
            class_weight="balanced",
            random_state=int(seed),
        )
        clf.fit(Xf, ybin)
        return clf, Xf

    # retry convergence by shrinking C
    C_try = float(best_params["C"])
    clf_final = None
    Xf = None

    for attempt in range(int(refit_max_tries)):
        clf_final, Xf = _fit_final(panel, C_try)

        if np.max(clf_final.n_iter_) < clf_final.max_iter:
            break

        warn(
            f"recommended_auc_panel: final fit not converged (attempt {attempt+1}); "
            f"shrinking C ({C_try} -> {C_try/refit_C_shrink})."
        )
        C_try = C_try / float(refit_C_shrink)

    if clf_final is None or (np.max(clf_final.n_iter_) >= clf_final.max_iter):
        raise ValueError(
            "recommended_auc_panel: final fit did not converge after retries. "
            "Try smaller C grid, higher l1_ratio, or preselect_top."
        )

    # store the refit C used
    best_params_final = dict(best_params)
    best_params_final["C_refit"] = float(C_try)


    # retry convergence by shrinking C (your logic kept)
    for attempt in range(int(refit_max_tries)):
        clf_final, Xf = _fit_final(panel, C_try)
        if np.max(clf_final.n_iter_) < clf_final.max_iter:
            break
        warn(
            f"recommended_auc_panel: final fit not converged (attempt {attempt+1}); "
            f"shrinking C ({C_try} -> {C_try/refit_C_shrink})."
        )
        C_try = C_try / float(refit_C_shrink)

    if np.max(clf_final.n_iter_) >= clf_final.max_iter:
        raise ValueError(
            "recommended_auc_panel: final fit did not converge after retries. "
            "Try smaller C grid, higher l1_ratio, or preselect_top."
        )

    # ---- enforce K nonzero genes (refill loop) ----
    max_refill_rounds = 10
    for _ in range(max_refill_rounds):
        beta = clf_final.coef_.ravel().astype(float)
        nz_mask = np.abs(beta) > float(coef_eps)
        nz_genes = [g for g, ok in zip(panel, nz_mask) if ok]

        if len(nz_genes) >= int(max_genes):
            panel = nz_genes[: int(max_genes)]
            clf_final, Xf = _fit_final(panel, C_try)
            beta = clf_final.coef_.ravel().astype(float)
            break

        # refill using ranked candidates not yet in panel
        needed = int(max_genes) - len(nz_genes)
        extra = [g for g in ranked_candidates if g not in nz_genes][:needed]

        # if no extra genes left, stop
        if len(extra) == 0:
            panel = nz_genes
            clf_final, Xf = _fit_final(panel, C_try)
            beta = clf_final.coef_.ravel().astype(float)
            break

        panel = nz_genes + extra
        clf_final, Xf = _fit_final(panel, C_try)

    # still not converged
    if np.max(clf_final.n_iter_) >= clf_final.max_iter:
        raise ValueError(
            "recommended_auc_panel: final fit did not converge after retries. "
            "Try smaller C grid, higher l1_ratio, or preselect_top."
        )

    # if you want, store the refit C used:
    best_params_final = dict(best_params)
    best_params_final["C_refit"] = float(C_try)



    beta = clf_final.coef_.ravel().astype(float)
    intercept = float(clf_final.intercept_.ravel()[0])

    weights_df = pd.DataFrame(
        {
            "gene": panel,
            "beta": beta,
        }
    ).assign(abs_beta=lambda d: np.abs(d["beta"])).sort_values("abs_beta", ascending=False).drop(columns="abs_beta").reset_index(drop=True)


    # ---- effective (non-zero) genes ----
    nz_mask = np.abs(beta) > float(coef_eps)

    n_nonzero = int(nz_mask.sum())
    panel_nonzero = [g for g, keep in zip(panel, nz_mask) if keep]



    out = {
        "best_params": best_params_final,
        "score_cv_best": score_cv_best,
        "df_cv": df_cv.sort_values(["score_mean"], ascending=False).reset_index(drop=True),
        "df_stability": df_stab,

        # panels
        "panel": panel,
        "panel_nonzero": panel_nonzero,         # effective genes
        "n_nonzero": n_nonzero,

        "weights_df": weights_df,
        "intercept": intercept, 

        "label_col": str(label_col),
        "positive_label": str(positive_label),
        "n": int(M.shape[0]),
        "n0": int((ybin == 0).sum()),
        "n1": int((ybin == 1).sum()),
        "redundancy_filter": bool(redundancy_filter),
        "corr_threshold": float(corr_threshold),
    }

    if key_added is not None:
        adata.uns.setdefault("recommended_auc_panel", {})
        adata.uns["recommended_auc_panel"][key_added] = out

    return out



def validate_signature_auc_external(
    adata,
    *,
    weights: pd.DataFrame,
    label_col: str = "PFS_6m",
    positive_label: str = "R",
    layer: str | None = "log1p_cpm",
    intercept: float = 0.0,
    score_key: str | None = None,
    use_proba: bool = True,
    metric: str = "roc_auc",     # "roc_auc" or "pr_auc"
    n_boot: int = 0,
    seed: int = 0,
    ppv_recalls: tuple[float, ...] = (0.3, 0.5, 0.7, 0.8),
    ppv_targets: tuple[float, ...] = (0.3, 0.5, 0.7, 0.8),
    spec_targets: tuple[float, ...] = (0.3, 0.8, 0.9, 0.95),
    key_added: str | None = None,
) -> pd.DataFrame:
    try:
        from sklearn.metrics import (
            roc_auc_score,
            average_precision_score,
            precision_recall_curve,
            roc_curve,
        )
    except Exception as e:
        raise ImportError(f"validate_signature_auc_external requires scikit-learn. ({e})")

    metric = str(metric).lower()
    if metric not in ("roc_auc", "pr_auc"):
        raise ValueError("metric must be 'roc_auc' or 'pr_auc'")

    if "gene" not in weights.columns or "beta" not in weights.columns:
        raise ValueError("weights must contain columns: gene, beta")
    if label_col not in adata.obs.columns:
        raise KeyError(f"label_col='{label_col}' not found in adata.obs")

    if score_key is None:
        score_key = "signature_score"

    # If this is inside signature.py, you can call signature_score directly.
    from .signature import signature_score
    signature_score(adata, weights=weights, layer=layer, key_added=score_key)

    s = pd.to_numeric(adata.obs[score_key], errors="coerce").to_numpy(float)
    yraw = adata.obs[label_col].astype("object")

    ok = np.isfinite(s) & yraw.notna().to_numpy()
    s = s[ok]
    y = (yraw.iloc[np.where(ok)[0]].astype(str) == str(positive_label)).astype(int).to_numpy()

    if s.size < 10 or len(np.unique(y)) < 2:
        raise ValueError("Not enough valid samples or <2 classes to compute metric.")

    def _sigmoid(z):
        z = np.clip(np.asarray(z, float), -50, 50)
        return 1.0 / (1.0 + np.exp(-z))

    pred = _sigmoid(float(intercept) + s) if use_proba else s

    def _confusion_from_threshold(y_true, y_score, thr):
        yhat = (y_score >= float(thr)).astype(int)
        tp = int(((yhat == 1) & (y_true == 1)).sum())
        fp = int(((yhat == 1) & (y_true == 0)).sum())
        tn = int(((yhat == 0) & (y_true == 0)).sum())
        fn = int(((yhat == 0) & (y_true == 1)).sum())
        return tp, fp, tn, fn

    def _ppv_at_recall(y_true, y_score, target_recall: float):
        prec, rec, thr = precision_recall_curve(y_true, y_score)
        idx = np.where(rec >= float(target_recall))[0]
        if idx.size == 0:
            return np.nan, np.nan
        j = idx[np.nanargmax(prec[idx])]
        ppv = float(prec[j])

        # thresholds has length len(rec)-1
        if j == 0:
            thresh = float("inf")
        elif (j - 1) < thr.size:
            thresh = float(thr[j - 1])
        else:
            thresh = np.nan
        return ppv, thresh

    def _recall_at_ppv(y_true, y_score, target_ppv: float):
        prec, rec, thr = precision_recall_curve(y_true, y_score)
        # thr has length len(prec)-1; ignore last point (no threshold)
        prec2, rec2 = prec[:-1], rec[:-1]
        thr2 = thr

        idx = np.where(prec2 >= float(target_ppv))[0]
        if idx.size == 0:
            return np.nan, np.nan

        j = idx[np.nanargmax(rec2[idx])]
        return float(rec2[j]), float(thr2[j])

    def _npv_at_spec(y_true, y_score, target_spec: float):
        fpr, tpr, thr = roc_curve(y_true, y_score)
        spec = 1.0 - fpr

        idx = np.where(spec >= float(target_spec))[0]
        if idx.size == 0:
            return np.nan, np.nan

        # choose threshold with max tpr among feasible
        j = idx[np.nanargmax(tpr[idx])]
        thr_j = float(thr[j])

        tp, fp, tn, fn = _confusion_from_threshold(y_true, y_score, thr_j)
        denom = tn + fn
        npv = float(tn / denom) if denom > 0 else np.nan
        return npv, thr_j

    # -------------------------
    # Main score
    # -------------------------
    if metric == "roc_auc":
        score = float(roc_auc_score(y, pred))
    else:
        score = float(average_precision_score(y, pred))

    out = {
        "metric": metric,
        "score": score,
        "n": int(s.size),
        "n_pos": int(y.sum()),
        "n_neg": int((y == 0).sum()),
        "use_proba": bool(use_proba),
    }

    # Point estimates
    for r in ppv_recalls:
        ppv, thr = _ppv_at_recall(y, pred, r)
        out[f"ppv_at_recall_{int(round(100*r))}"] = ppv
        out[f"thr_at_recall_{int(round(100*r))}"] = thr

    for p in ppv_targets:
        rec, thr = _recall_at_ppv(y, pred, p)
        out[f"recall_at_ppv_{int(round(100*p))}"] = rec
        out[f"thr_at_ppv_{int(round(100*p))}"] = thr

    for sp in spec_targets:
        npv, thr = _npv_at_spec(y, pred, sp)
        out[f"npv_at_spec_{int(round(100*sp))}"] = npv
        out[f"thr_at_spec_{int(round(100*sp))}"] = thr

    # -------------------------
    # Bootstrap CIs
    # -------------------------
    # store only metric distributions (thresholds are usually too unstable)
    ppv_boot = {r: [] for r in ppv_recalls}
    recall_ppv_boot = {p: [] for p in ppv_targets}
    npv_spec_boot = {sp: [] for sp in spec_targets}

    if int(n_boot) > 0:
        rng = np.random.default_rng(int(seed))
        boots = []
        n = y.size

        for _ in range(int(n_boot)):
            idx = rng.integers(0, n, size=n)
            yb = y[idx]
            if len(np.unique(yb)) < 2:
                continue
            pb = pred[idx]

            # metric
            boots.append(
                float(roc_auc_score(yb, pb)) if metric == "roc_auc" else float(average_precision_score(yb, pb))
            )

            # PPV@fixed recall
            for r in ppv_recalls:
                ppv_r, _ = _ppv_at_recall(yb, pb, r)
                if np.isfinite(ppv_r):
                    ppv_boot[r].append(float(ppv_r))

            # Recall@fixed PPV
            for p in ppv_targets:
                rec_p, _ = _recall_at_ppv(yb, pb, p)
                if np.isfinite(rec_p):
                    recall_ppv_boot[p].append(float(rec_p))

            # NPV@fixed specificity
            for sp in spec_targets:
                npv_sp, _ = _npv_at_spec(yb, pb, sp)
                if np.isfinite(npv_sp):
                    npv_spec_boot[sp].append(float(npv_sp))

        out["n_boot_ok"] = int(len(boots))

        if len(boots) >= 10:
            lo, hi = np.quantile(boots, [0.025, 0.975])
            out["score_ci_lower"] = float(lo)
            out["score_ci_upper"] = float(hi)
            out["score_boot_sd"] = float(np.std(boots))
        else:
            out["score_ci_lower"] = np.nan
            out["score_ci_upper"] = np.nan
            out["score_boot_sd"] = np.nan

        # PPV@recall CIs
        for r in ppv_recalls:
            arr = np.asarray(ppv_boot[r], float)
            key = f"ppv_at_recall_{int(round(100*r))}"
            if arr.size >= 20:
                lo, hi = np.quantile(arr, [0.025, 0.975])
                out[f"{key}_ci_lower"] = float(lo)
                out[f"{key}_ci_upper"] = float(hi)
                out[f"{key}_boot_sd"] = float(np.std(arr))
            else:
                out[f"{key}_ci_lower"] = np.nan
                out[f"{key}_ci_upper"] = np.nan
                out[f"{key}_boot_sd"] = np.nan

        # Recall@PPV CIs
        for p in ppv_targets:
            arr = np.asarray(recall_ppv_boot[p], float)
            key = f"recall_at_ppv_{int(round(100*p))}"
            if arr.size >= 20:
                lo, hi = np.quantile(arr, [0.025, 0.975])
                out[f"{key}_ci_lower"] = float(lo)
                out[f"{key}_ci_upper"] = float(hi)
                out[f"{key}_boot_sd"] = float(np.std(arr))
            else:
                out[f"{key}_ci_lower"] = np.nan
                out[f"{key}_ci_upper"] = np.nan
                out[f"{key}_boot_sd"] = np.nan

        # NPV@spec CIs
        for sp in spec_targets:
            arr = np.asarray(npv_spec_boot[sp], float)
            key = f"npv_at_spec_{int(round(100*sp))}"
            if arr.size >= 20:
                lo, hi = np.quantile(arr, [0.025, 0.975])
                out[f"{key}_ci_lower"] = float(lo)
                out[f"{key}_ci_upper"] = float(hi)
                out[f"{key}_boot_sd"] = float(np.std(arr))
            else:
                out[f"{key}_ci_lower"] = np.nan
                out[f"{key}_ci_upper"] = np.nan
                out[f"{key}_boot_sd"] = np.nan

    else:
        out["score_ci_lower"] = np.nan
        out["score_ci_upper"] = np.nan
        out["score_boot_sd"] = np.nan
        out["n_boot_ok"] = 0

        for r in ppv_recalls:
            key = f"ppv_at_recall_{int(round(100*r))}"
            out[f"{key}_ci_lower"] = np.nan
            out[f"{key}_ci_upper"] = np.nan
            out[f"{key}_boot_sd"] = np.nan

        for p in ppv_targets:
            key = f"recall_at_ppv_{int(round(100*p))}"
            out[f"{key}_ci_lower"] = np.nan
            out[f"{key}_ci_upper"] = np.nan
            out[f"{key}_boot_sd"] = np.nan

        for sp in spec_targets:
            key = f"npv_at_spec_{int(round(100*sp))}"
            out[f"{key}_ci_lower"] = np.nan
            out[f"{key}_ci_upper"] = np.nan
            out[f"{key}_boot_sd"] = np.nan

    df = pd.DataFrame([out])

    if key_added is not None:
        adata.uns.setdefault("validate_signature_auc_external", {})
        adata.uns["validate_signature_auc_external"][key_added] = df

    return df


def score_signature_on_PPV(
    adata_dev,
    *,
    weights_df,
    intercept,
    label_col="NR_6m",
    positive_label="NR",
    layer="log1p_cpm",
    recall_target=0.5,
    n_boot=200,
    seed=0,
):
    df = bk.tl.validate_signature_auc_external(
        adata_dev,
        weights=weights_df,
        intercept=float(intercept),
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        use_proba=True,
        metric="pr_auc",
        ppv_recalls=(float(recall_target),),
        n_boot=int(n_boot),
        seed=int(seed),
    )
    rkey = f"ppv_at_recall_{int(round(100*recall_target))}"
    return df.loc[0, rkey], df


def evaluate_signatures_store(
    adata,
    signatures_store: dict,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    use_proba=True,
    metric="pr_auc",
    # clinical metrics
    ppv_recalls=(0.3, 0.5, 0.7, 0.8),
    ppv_targets=(0.3, 0.5, 0.7, 0.8),
    spec_targets=(0.3, 0.8, 0.9, 0.95),
    # bootstrap
    n_boot=200,
    seed=0,
):
    rows = []
    for name, sig in signatures_store.items():
        w = sig["weights_df"]
        b0 = float(sig.get("intercept", 0.0))

        df = bk.tl.validate_signature_auc_external(
            adata,
            weights=w,
            intercept=b0,
            label_col=label_col,
            positive_label=positive_label,
            layer=layer,
            use_proba=use_proba,
            metric=metric,
            ppv_recalls=ppv_recalls,
            ppv_targets=ppv_targets,
            spec_targets=spec_targets,
            n_boot=int(n_boot),
            seed=int(seed),
            key_added=None,
        )

        out = df.iloc[0].to_dict()
        out["Signature"] = name

        # convenience counts
        out["n_genes"] = int(w.shape[0])
        out["n_nonzero"] = int((w["beta"].abs().to_numpy() > 1e-12).sum())

        # carry params (optional)
        params = sig.get("params", {})
        if isinstance(params, dict):
            # flatten a few common ones if present
            for k in ["C", "l1_ratio", "standardize", "train_cohort", "label_def"]:
                if k in params:
                    out[k] = params[k]

        rows.append(out)

    df_all = pd.DataFrame(rows)

    # nicer ordering
    front = ["Signature", "metric", "score", "score_ci_lower", "score_ci_upper", "n", "n_pos", "n_neg", "n_genes", "n_nonzero"]
    keep_front = [c for c in front if c in df_all.columns]
    other = [c for c in df_all.columns if c not in keep_front]
    df_all = df_all[keep_front + other]

    return df_all



def pick_best_signature(df_eval: pd.DataFrame, *, criterion="ppv_at_recall_50", higher_is_better=True):
    df = df_eval.copy()
    if criterion not in df.columns:
        raise KeyError(f"criterion '{criterion}' not found in df_eval columns.")
    df = df[np.isfinite(df[criterion].to_numpy())].copy()
    if df.shape[0] == 0:
        raise ValueError(f"No finite values for criterion '{criterion}'.")

    df = df.sort_values(criterion, ascending=not bool(higher_is_better)).reset_index(drop=True)
    best_name = str(df.loc[0, "Signature"])
    return best_name, df



def _build_Xy(adata, genes, *, label_col, positive_label, layer="log1p_cpm"):
    genes = [g for g in genes if g in adata.var_names]
    if len(genes) < 1:
        raise ValueError("No genes found in adata.var_names")

    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    idx = adata.var_names.get_indexer(genes)
    M = X[:, idx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, float)

    yraw = adata.obs[label_col].astype("object")
    ok = np.all(np.isfinite(M), axis=1) & yraw.notna().to_numpy()
    M = M[ok, :]
    y = (yraw.iloc[np.where(ok)[0]].astype(str) == str(positive_label)).astype(int).to_numpy()

    if len(np.unique(y)) < 2:
        raise ValueError("Need 2 classes to fit/evaluate.")
    return genes, M, y

def _sigmoid(z):
    z = np.clip(np.asarray(z, float), -50, 50)
    return 1.0 / (1.0 + np.exp(-z))

def _eval_with_frozen_weights(adata, weights_df, intercept, *, label_col, positive_label, layer="log1p_cpm", use_sigmoid=True):
    w = weights_df.copy()
    w["gene"] = w["gene"].astype(str)
    genes = [g for g in w["gene"].tolist() if g in adata.var_names]
    if len(genes) < 1:
        return {"roc_auc": np.nan, "pr_auc": np.nan, "n": 0, "pos": 0, "neg": 0}

    w = w.set_index("gene").loc[genes]
    betas = w["beta"].to_numpy(float)

    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    idx = adata.var_names.get_indexer(genes)
    M = X[:, idx]
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M, float)

    yraw = adata.obs[label_col].astype("object")
    ok = np.all(np.isfinite(M), axis=1) & yraw.notna().to_numpy()
    M = M[ok, :]
    y = (yraw.iloc[np.where(ok)[0]].astype(str) == str(positive_label)).astype(int).to_numpy()

    if len(np.unique(y)) < 2:
        return {"roc_auc": np.nan, "pr_auc": np.nan, "n": int(len(y)), "pos": int(y.sum()), "neg": int((y==0).sum())}

    score = M @ betas
    z = float(intercept) + score
    pred = _sigmoid(z) if use_sigmoid else z

    return {
        "roc_auc": float(roc_auc_score(y, pred)),
        "pr_auc": float(average_precision_score(y, pred)),
        "n": int(len(y)),
        "pos": int(y.sum()),
        "neg": int((y == 0).sum()),
    }

def make_signature_frozen(
    name,
    adata_train,
    adata_ext,
    genes,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    standardize=True,
    l1_ratio=0.8,                 # default elastic-net (recommended)
    C_grid=(0.01, 0.03, 0.1, 0.3, 1.0),   # tiny grid; uses TRAIN only
    max_iter=50000,
    tol=1e-3,
    coef_eps=1e-6,
    seed=0,
):

    # --- freeze requested genes ASAP (before any filtering happens) ---
    genes_req = [str(g) for g in list(genes)]  # full list you *intended* to test
    
    # --- build TRAIN X/y ---
    genes_used, Xtr, ytr = _build_Xy(
        adata_train, genes_req,  # <-- pass genes_req, not genes
        label_col=label_col, positive_label=positive_label, layer=layer
    )

    scaler = None
    if standardize:
        scaler = StandardScaler()
        Xtr_s = scaler.fit_transform(Xtr)
    else:
        Xtr_s = Xtr

    # --- pick C on TRAIN only (fast “in-sample” criterion; keep minimal to avoid overfitting) ---
    # With NR=6, proper CV is unstable. This grid is just to avoid pathological regularization.
    bestC = None
    best_obj = -np.inf
    best_clf = None

    rng = np.random.default_rng(int(seed))

    for C in C_grid:
        clf = LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            C=float(C),
            l1_ratio=float(l1_ratio),
            class_weight="balanced",
            max_iter=int(max_iter),
            tol=float(tol),
            random_state=int(rng.integers(0, 2**31 - 1)),
        )
        clf.fit(Xtr_s, ytr)

        # objective on TRAIN: PR-AUC (more relevant for imbalance)
        p = clf.predict_proba(Xtr_s)[:, 1]
        pr = float(average_precision_score(ytr, p))
        # tie-break: prefer sparsity (fewer nonzero)
        nnz = int((np.abs(clf.coef_.ravel()) > float(coef_eps)).sum())
        obj = pr - 0.001 * nnz

        if obj > best_obj:
            best_obj = obj
            bestC = float(C)
            best_clf = clf

    beta = best_clf.coef_.ravel().astype(float)
    intercept = float(best_clf.intercept_.ravel()[0])

    # ---- weights for ALL requested genes (old + new) ----
    # genes_used is the subset found in adata_train.var_names (order matches Xtr columns)
    weights_fit = pd.DataFrame({"gene": genes_used, "beta": beta})

    weights_df = (
        pd.DataFrame({"gene": genes_req})      # <-- full requested list
        .merge(weights_fit, on="gene", how="left")
    )
    weights_df["beta"] = weights_df["beta"].fillna(0.0)

    # sort by abs(beta) for readability (optional)
    weights_df["abs_beta"] = np.abs(weights_df["beta"])
    weights_df = (
        weights_df.sort_values("abs_beta", ascending=False)
        .drop(columns="abs_beta")
        .reset_index(drop=True)
    )

    # convenience: dict gene->beta
    weights_map = dict(zip(weights_df["gene"].tolist(), weights_df["beta"].astype(float).tolist()))

    panel_nonzero = weights_df.loc[np.abs(weights_df["beta"]) > float(coef_eps), "gene"].tolist()
    n_nonzero = int(len(panel_nonzero))
    
    # --- evaluate frozen weights on TRAIN and EXT (no refit) ---
    m_tr  = _eval_with_frozen_weights(adata_train, weights_df, intercept, label_col=label_col, positive_label=positive_label, layer=layer)
    m_ext = _eval_with_frozen_weights(adata_ext,   weights_df, intercept, label_col=label_col, positive_label=positive_label, layer=layer)

    sig_obj = {
        "name": str(name),
        "genes": list(genes_used),
        "weights_df": weights_df,
        "genes_requested": genes_req,      # <- old+new 
        "genes_used": list(genes_used),                  # <- actually in var_names
        "weights_map": weights_map,                      # <- quick lookup
        "intercept": float(intercept),
        "panel_nonzero": panel_nonzero,
        "n_nonzero": n_nonzero,
        "params": {
            "model": "LogReg(elastic-net)",
            "l1_ratio": float(l1_ratio),
            "C": float(bestC),
            "standardize": bool(standardize),
            "label_col": str(label_col),
            "positive_label": str(positive_label),
            "layer": str(layer),
        },
        "metrics_train": m_tr,
        "metrics_ext": m_ext,
    }

    row = {
        "Signature": sig_obj["name"],
        "n_genes_requested": int(len(sig_obj["genes_requested"])),
        "n_genes_used": int(len(sig_obj["genes_used"])),
        "n_genes": int(len(sig_obj["genes"])),
        "n_nonzero": int(sig_obj["n_nonzero"]),
        "roc_auc_train": m_tr["roc_auc"],
        "pr_auc_train": m_tr["pr_auc"],
        "roc_auc_ext": m_ext["roc_auc"],
        "pr_auc_ext": m_ext["pr_auc"],
        "C": float(bestC),
        "l1_ratio": float(l1_ratio),
        "pos_train": m_tr["pos"],
        "neg_train": m_tr["neg"],
        "pos_ext": m_ext["pos"],
        "neg_ext": m_ext["neg"],
    }

    return sig_obj, pd.DataFrame([row])


def quick_fit_weights(
    adata,
    genes,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    standardize=True,
    C=0.1,
    l1_ratio=0.8,
    max_iter=50000,
    tol=1e-3,
    seed=0,
):
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    # build X/y using YOUR helper
    genes_used, X, y = _build_Xy(
        adata, genes,
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
    )

    if standardize:
        X = StandardScaler().fit_transform(X)

    clf = LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        C=float(C),
        l1_ratio=float(l1_ratio),
        class_weight="balanced",
        max_iter=int(max_iter),
        tol=float(tol),
        random_state=int(seed),
    )
    clf.fit(X, y)

    beta = clf.coef_.ravel().astype(float)
    intercept = float(clf.intercept_.ravel()[0])

    weights_df = pd.DataFrame({
        "gene": genes_used,
        "beta": beta,
    }).sort_values("beta", key=np.abs, ascending=False).reset_index(drop=True)

    return weights_df, intercept


def quick_fit_weights_frozen(
    adata_train,
    genes,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    standardize=True,
    l1_ratio=0.8,
    C_grid=(0.01, 0.03, 0.1, 0.3, 1.0),
    max_iter=50000,
    tol=1e-3,
    coef_eps=1e-6,
    seed=0,
):
    # Use the SAME function as your whole pipeline
    sig_obj, _row = bk.tl.make_signature_frozen(
        name="__tmp__",
        adata_train=adata_train,
        adata_ext=adata_train,   # dummy; we won't use metrics_ext
        genes=genes,
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        standardize=standardize,
        l1_ratio=l1_ratio,
        C_grid=C_grid,
        max_iter=max_iter,
        tol=tol,
        coef_eps=coef_eps,
        seed=seed,
    )
    return sig_obj["weights_df"], float(sig_obj["intercept"]), sig_obj

def add_signature_from_genes(
    sigs: dict,
    name: str,
    adata_train,
    genes,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    **fit_kwargs,
):
    w, b0, sig_obj = quick_fit_weights_frozen(
        adata_train,
        genes,
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        **fit_kwargs,
    )
    sigs[name] = {"weights_df": w, "intercept": b0, "sig_obj": sig_obj}
    return sigs

def benchmark_signatures_pr_auc(
    adata,
    sigs: dict,
    *,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    n_boot=200,
    seed=0,
):
    rows = []

    for name, sig in sigs.items():
        df = validate_signature_auc_external(
            adata,
            weights=sig["weights_df"],
            intercept=sig.get("intercept", 0.0),
            label_col=label_col,
            positive_label=positive_label,
            layer=layer,
            metric="pr_auc",
            use_proba=False,     # <-- IMPORTANT
            n_boot=n_boot,
            seed=seed,
        )

        r = df.iloc[0].to_dict()
        r["Signature"] = name
        rows.append(r)

    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)


#compute intercepts from weights
def fit_intercept_only(y, s, max_iter=100, tol=1e-10):
    """
    Fit b0 in sigmoid(b0 + s) by MLE (Newton-Raphson).
    y: {0,1}, s: raw score (fixed)
    """
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    b0 = 0.0
    for _ in range(max_iter):
        p = 1.0 / (1.0 + np.exp(-(b0 + s)))
        g = np.sum(y - p)                  # derivative
        h = -np.sum(p * (1 - p))           # second derivative
        if abs(h) < 1e-12:
            break
        step = g / h
        b0_new = b0 - step
        if abs(b0_new - b0) < tol:
            b0 = b0_new
            break
        b0 = b0_new
    return float(b0)


def add_signature_from_genes(
    signatures_store: dict,
    *,
    name: str,
    genes: list[str],
    adata_train,
    adata_dev=None,
    adata_test=None,
    label_col="NR_6m",
    positive_label="NR",
    layer="log1p_cpm",
    # fitter settings
    standardize=True,
    l1_ratio=0.8,
    C_grid=(0.01, 0.03, 0.1, 0.3, 1.0),
    max_iter=50000,
    tol=1e-3,
    coef_eps=1e-6,
    seed=0,
):
    # if you don't want to evaluate inside, you can pass adata_dev=adata_train (or None and skip)
    adata_ext = adata_dev if adata_dev is not None else adata_train

    sig_obj, df_row = make_signature_frozen(
        name,
        adata_train,
        adata_ext,
        genes,
        label_col=label_col,
        positive_label=positive_label,
        layer=layer,
        standardize=standardize,
        l1_ratio=l1_ratio,
        C_grid=C_grid,
        max_iter=max_iter,
        tol=tol,
        coef_eps=coef_eps,
        seed=seed,
    )

    signatures_store[str(name)] = {
        "weights_df": sig_obj["weights_df"],
        "intercept": float(sig_obj["intercept"]),
        "params": sig_obj.get("params", {}),
    }

    # optional evaluation on test, without refitting:
    if adata_test is not None:
        # uses locked weights; just evaluate
        pass

    return signatures_store, sig_obj, df_row

################################################################
## BUILD Signature Bank from Genes
################################################################

def build_signature_bank_from_genes(
    signature_defs: dict,
    *,
    adata_train,
    adata_dev=None,
    adata_test=None,
    label_col="PFS_6m",
    positive_label="NR",
    layer="log1p_cpm",
    l1_ratio=0.8,
    C_grid=(0.01, 0.03, 0.1, 0.3, 1.0),
    coef_eps=1e-6,
    seed=0,
    min_genes_in_train=2,
    on_missing="skip",   # "skip" or "raise"
):
    """
    Build sig_bank from gene lists.
    Returns:
      sig_bank: dict name -> sig_obj
      df_compare: concatenated rows returned by make_signature_frozen (DEV metrics)
      df_diag: per-signature diagnostics about gene matching
    """
    import numpy as np
    import pandas as pd
    import bullkpy as bk

    train_genes = set(map(str, adata_train.var_names))
    dev_genes = set(map(str, adata_dev.var_names)) if adata_dev is not None else None
    test_genes = set(map(str, adata_test.var_names)) if adata_test is not None else None

    sig_bank = {}
    rows = []
    diag_rows = []

    for name, genes in signature_defs.items():
        genes_req = [str(g) for g in list(genes)]

        genes_in_train = [g for g in genes_req if g in train_genes]
        missing_train = [g for g in genes_req if g not in train_genes]

        genes_in_dev = [g for g in genes_req if (dev_genes is not None and g in dev_genes)] if dev_genes is not None else None
        genes_in_test = [g for g in genes_req if (test_genes is not None and g in test_genes)] if test_genes is not None else None

        diag_rows.append({
            "Signature": name,
            "n_requested": len(genes_req),
            "n_in_train": len(genes_in_train),
            "n_missing_train": len(missing_train),
            "missing_train_preview": ", ".join(missing_train[:10]) + (" ..." if len(missing_train) > 10 else ""),
            "n_in_dev": (len(genes_in_dev) if genes_in_dev is not None else np.nan),
            "n_in_test": (len(genes_in_test) if genes_in_test is not None else np.nan),
        })

        if len(genes_in_train) < int(min_genes_in_train):
            msg = (
                f"[{name}] Only {len(genes_in_train)}/{len(genes_req)} genes found in adata_train.var_names "
                f"(min required={min_genes_in_train}). Example missing: {missing_train[:5]}"
            )
            if on_missing == "raise":
                raise ValueError(msg)
            else:
                print("WARNING:", msg)
                continue

        # IMPORTANT: pass genes_in_train (not genes_req) to avoid _build_Xy crashing
        sig_obj, df_row = make_signature_frozen(
            name=str(name),
            adata_train=adata_train,
            adata_ext=(adata_dev if adata_dev is not None else adata_train),
            genes=genes_in_train,
            label_col=label_col,
            positive_label=positive_label,
            layer=layer,
            standardize=True,
            l1_ratio=l1_ratio,
            C_grid=C_grid,
            coef_eps=coef_eps,
            seed=seed,
        )

        # Optional: evaluate on TEST and store in sig_obj
        if adata_test is not None:
            m_test = validate_signature_auc_external(
                adata_test,
                weights=sig_obj["weights_df"],
                intercept=sig_obj["intercept"],
                label_col=label_col,
                positive_label=positive_label,
                layer=layer,
                metric="pr_auc",
                n_boot=200,
                seed=seed,
            ).iloc[0].to_dict()
            sig_obj["metrics_test"] = m_test

        sig_bank[sig_obj["name"]] = sig_obj
        rows.append(df_row)

    df_compare = pd.concat(rows, ignore_index=True) if len(rows) else pd.DataFrame()
    df_diag = pd.DataFrame(diag_rows).sort_values(["n_in_train", "n_requested"], ascending=[True, False]).reset_index(drop=True)

    return sig_bank, df_compare, df_diag




########################################################
## SAVE and RELOAD Signature Banks
########################################################

def save_signature_bank(sig_bank: dict, outdir: str):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    meta = {"signatures": {}}

    for name, sig in sig_bank.items():
        # save weights as CSV (clean + easy diff)
        w_path = outdir / f"{name}.weights.tsv"
        sig["weights_df"].to_csv(w_path, sep="\t", index=False)

        # store everything else in JSON-friendly form
        meta["signatures"][name] = {
            "intercept": float(sig.get("intercept", 0.0)),
            "params": sig.get("params", {}),
            "genes_requested": sig.get("genes_requested", sig.get("genes", [])),
            "genes_used": sig.get("genes_used", sig.get("genes", [])),
            "panel_nonzero": sig.get("panel_nonzero", []),
            "n_nonzero": int(sig.get("n_nonzero", 0)),
            # optional quick metrics snapshots (if present)
            "metrics_train": sig.get("metrics_train", None),
            "metrics_ext": sig.get("metrics_ext", None),
            "weights_tsv": w_path.name,
        }

    with open(outdir / "sig_bank.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Saved {len(sig_bank)} signatures to: {outdir}")


def load_signature_bank(indir: str):
    indir = Path(indir)
    with open(indir / "sig_bank.json", "r") as f:
        meta = json.load(f)

    sig_bank = {}
    for name, rec in meta["signatures"].items():
        w = pd.read_csv(indir / rec["weights_tsv"], sep="\t")
        sig_bank[name] = {
            "name": name,
            "weights_df": w,
            "intercept": float(rec.get("intercept", 0.0)),
            "params": rec.get("params", {}),
            "genes_requested": rec.get("genes_requested", []),
            "genes_used": rec.get("genes_used", []),
            "panel_nonzero": rec.get("panel_nonzero", []),
            "n_nonzero": int(rec.get("n_nonzero", 0)),
            "metrics_train": rec.get("metrics_train", None),
            "metrics_ext": rec.get("metrics_ext", None),
        }
    return sig_bank
