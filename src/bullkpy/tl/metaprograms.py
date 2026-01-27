from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple, Union
import re
import numpy as np
import pandas as pd
from anndata import AnnData


# -----------------------------
# Internal helpers
# -----------------------------
def _to_dense(X) -> np.ndarray:
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def _softmax_rows(X: np.ndarray) -> np.ndarray:
    X = X - np.nanmax(X, axis=1, keepdims=True)
    E = np.exp(X)
    denom = np.nansum(E, axis=1, keepdims=True)
    denom[denom == 0] = np.nan
    return E / denom


def _mad_1d(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    return float(np.median(np.abs(x - med)))


def _mad_rows(X: np.ndarray) -> np.ndarray:
    med = np.nanmedian(X, axis=1, keepdims=True)
    return np.nanmedian(np.abs(X - med), axis=1)

def _zscore_within_groups(X: np.ndarray, groups: pd.Series) -> np.ndarray:
    """
    Z-score columns within each group.

    Robust to non-integer obs indices by using positional indices.
    """
    # Ensure positional integer index for groups
    # groups is length n_samples; its index may be sample IDs (strings)
    groups = pd.Series(groups.to_numpy(), index=np.arange(len(groups)))

    Z = np.full_like(X, np.nan, dtype=float)

    # observed=False silences future warning + keeps current behavior
    for g, idx in groups.groupby(groups, observed=False).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        Xg = X[idx, :]

        mu = np.nanmean(Xg, axis=0, keepdims=True)
        sd = np.nanstd(Xg, axis=0, keepdims=True)
        sd[sd == 0] = np.nan

        Z[idx, :] = (Xg - mu) / sd

    return Z


# -----------------------------
# Public API
# -----------------------------
def metaprogram_scores_get(
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


def metaprogram_sample_metrics(
    adata: AnnData,
    *,
    groupby: str = "Project_ID",
    mp_source: str = "auto",
    out_prefix: str = "mp",
    heterogeneity: str = "entropy",   # "entropy" | "gini" | "dominance"
    dispersion: str = "mad",          # "mad" | "iqr" | "std"
    zscore_within_group: bool = True,
    store_weights: bool = False,
    weights_key: str = "X_mp_weights",
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Compute per-sample metaprogram heterogeneity and dispersion.

    Heterogeneity (within sample, across programs):
      - entropy: normalized entropy of softmax(program scores) in [0, 1]
      - gini: Gini inequality of softmax weights in [0, 1]
      - dominance: max(weight) or (max - median) style dominance

    Dispersion (within sample, across programs):
      - computed on metaprogram scores; optionally z-scored within `groupby`

    Writes to:
      adata.obs[f"{out_prefix}_heterogeneity_<...>"]
      adata.obs[f"{out_prefix}_dispersion_<...>"]

    Optionally stores metaprogram weights (softmax) to adata.obsm[weights_key].
    """
    ad = adata.copy() if copy else adata
    if groupby not in ad.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    df = metaprogram_scores_get(ad, return_df=True)
    X = df.to_numpy(dtype=float)
    n_samples, n_prog = X.shape

    # Weights for heterogeneity
    W = _softmax_rows(X)
    if store_weights:
        ad.obsm[weights_key] = W
        ad.uns.setdefault("mp_weights_names", list(df.columns))

    # Heterogeneity
    if heterogeneity == "entropy":
        eps = 1e-12
        H = -np.nansum(W * np.log(W + eps), axis=1) / np.log(n_prog)
        het_key = f"{out_prefix}_heterogeneity_entropy"
        ad.obs[het_key] = H
    elif heterogeneity == "gini":
        # Gini on weights per row
        G = np.full(n_samples, np.nan, dtype=float)
        for i in range(n_samples):
            w = W[i, :]
            w = w[np.isfinite(w)]
            if w.size == 0:
                continue
            w = np.sort(w)
            cumw = np.cumsum(w)
            if cumw[-1] == 0:
                continue
            # classic gini
            n = w.size
            G[i] = (n + 1 - 2 * np.sum(cumw) / cumw[-1]) / n
        het_key = f"{out_prefix}_heterogeneity_gini"
        ad.obs[het_key] = G
    elif heterogeneity == "dominance":
        dom = np.nanmax(W, axis=1)
        het_key = f"{out_prefix}_heterogeneity_dominance"
        ad.obs[het_key] = dom
    else:
        raise ValueError("heterogeneity must be one of {'entropy','gini','dominance'}")

    # Dispersion
    groups = ad.obs[groupby].astype("category")
    Z = X
    if zscore_within_group:
        Z = _zscore_within_groups(X, groups)

    if dispersion == "mad":
        D = _mad_rows(Z)
    elif dispersion == "iqr":
        D = np.nanpercentile(Z, 75, axis=1) - np.nanpercentile(Z, 25, axis=1)
    elif dispersion == "std":
        D = np.nanstd(Z, axis=1)
    else:
        raise ValueError("dispersion must be one of {'mad','iqr','std'}")

    disp_key = f"{out_prefix}_dispersion_{dispersion}" + (f"_zwithin_{groupby}" if zscore_within_group else "")
    ad.obs[disp_key] = D

    return ad if copy else None


def metaprogram_dispersion_by_group(
    adata: AnnData,
    *,
    groupby: str = "Project_ID",
    method: str = "mad",             # "mad" | "iqr" | "std"
    out_key: str = "metaprogram_dispersion",
) -> pd.DataFrame:
    """
    Compute per-group dispersion across samples for each metaprogram.

    Output table columns:
      group, metaprogram, dispersion, n

    Stores result in adata.uns[out_key] and returns it.
    """
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} not found in adata.obs")

    df = metaprogram_scores_get(adata, return_df=True)

    rows = []
    # IMPORTANT: use integer positions for indexing X, not string obs_names
    X = df.to_numpy(dtype=float)
    groups = pd.Series(adata.obs[groupby].to_numpy(), index=np.arange(adata.n_obs))

    for grp, idx in groups.groupby(groups, observed=False).groups.items():
        idx = np.asarray(list(idx), dtype=int)
        Xg = X[idx, :]  # samples in this group

        for j, mp in enumerate(df.columns):
            x = Xg[:, j]
            x = x[np.isfinite(x)]
            if x.size == 0:
                d = np.nan
                n = 0
            else:
                n = int(x.size)
                if method == "mad":
                    med = np.median(x)
                    d = float(np.median(np.abs(x - med)))
                elif method == "iqr":
                    d = float(np.percentile(x, 75) - np.percentile(x, 25))
                elif method == "std":
                    d = float(np.std(x))
                else:
                    raise ValueError("method must be one of {'mad','iqr','std'}")

            rows.append({"group": str(grp), "metaprogram": str(mp), "dispersion": d, "n": n})

    out = pd.DataFrame(rows)
    adata.uns[out_key] = out
    return out


def metaprogram_ne_enrichment(
    adata: AnnData,
    *,
    ne_signature_key: str,
    mp_score_threshold: float = 0.0,
    method: str = "spearman",
    out_key: str = "mp_ne_association",
) -> pd.DataFrame:
    """
    Associate metaprogram activities with a neuroendocrine (NE) signature score.

    Parameters
    ----------
    ne_signature_key
        Column in adata.obs with NE score.
    mp_score_threshold
        Optional filter: keep programs whose absolute mean activity exceeds threshold.
    method
        Currently supports 'spearman'.

    Returns
    -------
    DataFrame with columns: metaprogram, rho, pval, qval, n
    Stored in adata.uns[out_key].
    """
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import multipletests

    if ne_signature_key not in adata.obs:
        raise KeyError(f"{ne_signature_key!r} not found in adata.obs")

    df = metaprogram_scores_get(adata, return_df=True)
    y = adata.obs[ne_signature_key].to_numpy(dtype=float)

    # optional filter
    keep = np.ones(df.shape[1], dtype=bool)
    if mp_score_threshold > 0:
        keep = np.abs(np.nanmean(df.to_numpy(dtype=float), axis=0)) >= mp_score_threshold
        df = df.loc[:, keep]

    rows = []
    for mp in df.columns:
        x = df[mp].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        n = int(ok.sum())
        if n < 10:
            rows.append({"metaprogram": mp, "rho": np.nan, "pval": np.nan, "n": n})
            continue
        rho, p = spearmanr(x[ok], y[ok])
        rows.append({"metaprogram": mp, "rho": float(rho), "pval": float(p), "n": n})

    out = pd.DataFrame(rows)
    pvals = out["pval"].to_numpy(dtype=float)
    q = np.full_like(pvals, np.nan, dtype=float)
    okp = np.isfinite(pvals)
    if okp.any():
        q[okp] = multipletests(pvals[okp], method="fdr_bh")[1]
    out["qval"] = q
    out = out.sort_values(["qval", "pval"], na_position="last").reset_index(drop=True)

    adata.uns[out_key] = out
    return out


def metaprogram_topk_contribution(
    adata: AnnData,
    *,
    k: int = 3,
    out_key: str = "mp_topk",
) -> pd.DataFrame:
    """
    Compute, per sample, the top-k metaprograms and their weights (softmax).

    Returns a tidy DataFrame with columns:
      sample, rank, metaprogram, weight, score

    Stores result in adata.uns[out_key].
    """
    df = metaprogram_scores_get(adata, return_df=True)
    X = df.to_numpy(dtype=float)
    W = _softmax_rows(X)

    rows = []
    for i, sample in enumerate(df.index):
        order = np.argsort(W[i, :])[::-1][:k]
        for r, j in enumerate(order, start=1):
            rows.append({
                "sample": sample,
                "rank": r,
                "metaprogram": df.columns[j],
                "weight": float(W[i, j]),
                "score": float(X[i, j]),
            })
    out = pd.DataFrame(rows)
    adata.uns[out_key] = out
    return out


########################################
## METAPROGRAMS
#######################################


def read_metaprograms_xlsx(
    xlsx_path: str,
    *,
    sheet_name: str = "Malignant",
    min_genes: int = 5,
) -> dict[str, list[str]]:
    """
    Excel format expected: each column is an MP, rows are gene symbols.
    Returns dict mp_name -> list of genes (non-empty strings).
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

    mp = {}
    for col in df.columns:
        genes = (
            df[col]
            .dropna()
            .astype(str)
            .str.strip()
            .replace({"": np.nan, "nan": np.nan, "None": np.nan})
            .dropna()
            .tolist()
        )
        genes = [g for g in genes if g and g.lower() != "nan"]
        if len(genes) >= int(min_genes):
            mp[str(col)] = genes

    if len(mp) == 0:
        raise ValueError(f"No metaprograms found in sheet='{sheet_name}' of {xlsx_path}")

    return mp



def _get_X(adata, layer="log1p_cpm"):
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X
    return X

def _dense(M):
    return M.toarray() if hasattr(M, "toarray") else np.asarray(M, float)



def score_metaprograms(
    adata,
    mp_dict: dict[str, list[str]],
    *,
    layer: str | None = "log1p_cpm",
    zscore_genes: bool = True,
    center: bool = True,
    scale: bool = True,
    min_genes_found: int = 5,
    prefix: str = "MP_",
    obsm_key: str = "X_mp",
    uns_key: str = "metaprograms",
    add_obs: bool = True,   # also write each MP to adata.obs[prefix+mp_short]
):
    """
    Scores metaprograms in bulk: score(sample, MP) = mean(zscore(expression genes in MP)).

    Stores:
      - adata.obsm[obsm_key] = (n_samples, n_MPs) array
      - adata.uns[uns_key]   = metadata table about MPs (n_found, genes_found, ...)
      - optionally adata.obs[prefix+...] for each MP score
    Returns (df_scores, df_mp_info)
    """
    # ---- expression ----
    X = adata.layers[layer] if (layer is not None and layer in getattr(adata, "layers", {})) else adata.X

    # ensure dense float
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X, dtype=float)

    var_names = adata.var_names.astype(str).tolist()
    gene_to_idx = {g: i for i, g in enumerate(var_names)}

    mp_names = []
    mp_scores = []
    info_rows = []

    for mp_name, genes in mp_dict.items():
        genes = [str(g) for g in genes]
        idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        genes_found = [var_names[i] for i in idx]

        if len(idx) < int(min_genes_found):
            info_rows.append({
                "mp": mp_name,
                "n_genes_listed": int(len(genes)),
                "n_genes_found": int(len(idx)),
                "genes_found": genes_found,
                "kept": False,
                "reason": f"too_few_genes_found(<{min_genes_found})",
            })
            continue

        M = X[:, idx]  # (n_samples, n_genes_in_mp)

        if zscore_genes:
            # z-score each gene across samples (within this cohort)
            mu = np.nanmean(M, axis=0) if center else 0.0
            sd = np.nanstd(M, axis=0) if scale else 1.0
            sd = np.where(sd == 0, 1.0, sd)
            Mz = (M - mu) / sd
            s = np.nanmean(Mz, axis=1)
        else:
            s = np.nanmean(M, axis=1)

        mp_names.append(mp_name)
        mp_scores.append(s)

        info_rows.append({
            "mp": mp_name,
            "n_genes_listed": int(len(genes)),
            "n_genes_found": int(len(idx)),
            "genes_found": genes_found,
            "kept": True,
            "reason": "",
        })

    if len(mp_names) == 0:
        raise ValueError("No MPs passed min_genes_found; relax threshold or check gene symbols.")

    S = np.vstack(mp_scores).T  # (n_samples, n_kept_mps)
    df_scores = pd.DataFrame(S, index=adata.obs_names, columns=mp_names)

    df_info = pd.DataFrame(info_rows).sort_values(["kept", "n_genes_found"], ascending=[False, False]).reset_index(drop=True)

    # store
    adata.obsm[obsm_key] = S
    adata.uns[uns_key] = {
        "mp_names": mp_names,
        "params": {
            "layer": layer,
            "zscore_genes": bool(zscore_genes),
            "min_genes_found": int(min_genes_found),
            "prefix": str(prefix),
        },
        "df_info": df_info,
    }

    if add_obs:
        # safe column names in obs
        def _safe(name: str) -> str:
            return re.sub(r"[^0-9a-zA-Z_]+", "_", name).strip("_")

        for mp_name in mp_names:
            col = prefix + _safe(mp_name)
            adata.obs[col] = df_scores[mp_name].to_numpy(float)

    return df_scores, df_info


def mp_dispersion_metrics(
    adata,
    *,
    mp_key="MP_scores",
    key_added="MP_dispersion",
    transform="softmax",   # "softmax" or "absnorm"
    temperature=1.0,
    eps=1e-12,
):
    S = np.asarray(adata.obsm[mp_key], float)
    # fill NaN -> 0 for dispersion calc (or you can row-mask)
    S = np.where(np.isfinite(S), S, 0.0)

    if transform == "softmax":
        Z = S / float(temperature)
        Z = Z - np.max(Z, axis=1, keepdims=True)
        P = np.exp(Z)
        P = P / np.maximum(np.sum(P, axis=1, keepdims=True), eps)
    elif transform == "absnorm":
        P = np.abs(S)
        P = P / np.maximum(np.sum(P, axis=1, keepdims=True), eps)
    else:
        raise ValueError("transform must be 'softmax' or 'absnorm'")

    entropy = -np.sum(P * np.log(np.maximum(P, eps)), axis=1)
    entropy_norm = entropy / np.log(P.shape[1])  # 0..1

    top1 = np.max(P, axis=1)
    # top2:
    top2 = np.partition(P, -2, axis=1)[:, -2]
    dominance = top1 - top2

    df = pd.DataFrame({
        "mp_entropy": entropy,
        "mp_entropy_norm": entropy_norm,
        "mp_top1": top1,
        "mp_top1_minus_top2": dominance,
    }, index=adata.obs_names)

    adata.obs[key_added + "_entropy"] = df["mp_entropy"]
    adata.obs[key_added + "_entropy_norm"] = df["mp_entropy_norm"]
    adata.obs[key_added + "_top1"] = df["mp_top1"]
    adata.obs[key_added + "_top1_minus_top2"] = df["mp_top1_minus_top2"]

    adata.uns.setdefault("metaprograms", {})
    adata.uns["metaprograms"][key_added] = {
        "mp_key": mp_key,
        "transform": transform,
        "temperature": float(temperature),
    }

    return df


def metaprogram_heterogeneity(
    adata,
    *,
    groupby,
    obsm_key="X_mp",
    mp_names=None,
    stat="sd",
):
    if groupby not in adata.obs.columns:
        raise KeyError(f"groupby='{groupby}' not found in adata.obs")

    if obsm_key not in adata.obsm:
        raise KeyError(f"obsm_key='{obsm_key}' not found in adata.obsm")

    X = np.asarray(adata.obsm[obsm_key], float)  # (n_samples, n_MPs)

    if mp_names is None:
        # try to recover MP names if stored, else default MP1..MPk
        mp_names = adata.uns.get("mp_names", None)
    if mp_names is None or len(mp_names) != X.shape[1]:
        mp_names = [f"MP{i+1}" for i in range(X.shape[1])]

    g = adata.obs[groupby].astype(str).to_numpy()

    rows = []
    for grp in np.unique(g):
        idx = np.flatnonzero(g == grp)
        if idx.size < 2:
            continue

        M = X[idx, :]  # integer indexing works

        mu = np.nanmean(M, axis=0)

        if stat == "sd":
            het = np.nanstd(M, axis=0, ddof=1)
        elif stat == "mad":
            med = np.nanmedian(M, axis=0)
            het = np.nanmedian(np.abs(M - med), axis=0)
        elif stat == "iqr":
            q75 = np.nanpercentile(M, 75, axis=0)
            q25 = np.nanpercentile(M, 25, axis=0)
            het = q75 - q25
        else:
            raise ValueError("stat must be one of: 'sd','mad','iqr'")

        for j, mp in enumerate(mp_names):
            rows.append(
                {
                    "group": grp,
                    "mp": mp,
                    "mean": float(mu[j]),
                    f"heterogeneity_{stat}": float(het[j]),
                    "n": int(idx.size),
                }
            )

    return pd.DataFrame(rows)



def summarize_by_group(adata, group_col="Project ID", cols=("MP_dispersion_entropy_norm",)):
    out = []
    for g, d in adata.obs.groupby(group_col):
        row = {"group": g, "n": int(d.shape[0])}
        for c in cols:
            x = pd.to_numeric(d[c], errors="coerce")
            row[c + "_mean"] = float(np.nanmean(x))
            row[c + "_sd"]   = float(np.nanstd(x))
        out.append(row)
    return pd.DataFrame(out).sort_values("n", ascending=False)




def get_metaprogram_matrix(adata, key_candidates=("X_metaprogram", "X_metaprograms")):
    # 1) Prefer obsm matrix
    for k in key_candidates:
        if k in adata.obsm:
            X = adata.obsm[k]
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            names = adata.uns.get("metaprogram_names", [f"MP{i+1}" for i in range(X.shape[1])])
            return pd.DataFrame(X, index=adata.obs_names, columns=names)

    # 2) Fall back to uns dataframe
    if "metaprogram_scores" in adata.uns:
        df = adata.uns["metaprogram_scores"]
        return df.loc[adata.obs_names]

    raise KeyError("Metaprogram matrix not found in adata.obsm or adata.uns")



def softmax_rows(X):
    X = X - np.nanmax(X, axis=1, keepdims=True)
    E = np.exp(X)
    return E / np.nansum(E, axis=1, keepdims=True)

def mad_rows(X):
    med = np.nanmedian(X, axis=1, keepdims=True)
    return np.nanmedian(np.abs(X - med), axis=1)

def metaprogram_heterogeneity_and_dispersion(
    adata,
    *,
    mp_key=None,
    groupby="Project_ID",
    out_prefix="mp",
):
    df = get_metaprogram_matrix(adata)
    X = df.to_numpy(dtype=float)

    # (H1) Entropy heterogeneity (0-1)
    W = softmax_rows(X)
    eps = 1e-12
    H = -np.nansum(W * np.log(W + eps), axis=1) / np.log(W.shape[1])

    # (D) Dispersion: z-score within Project_ID, then MAD across programs
    disp = np.full(X.shape[0], np.nan, dtype=float)
    for grp, idx in adata.obs.groupby(groupby).groups.items():
        idx = np.asarray(list(idx))
        Z = X[idx, :]
        mu = np.nanmean(Z, axis=0, keepdims=True)
        sd = np.nanstd(Z, axis=0, keepdims=True)
        sd[sd == 0] = np.nan
        Z = (Z - mu) / sd
        disp[idx] = mad_rows(Z)

    adata.obs[f"{out_prefix}_heterogeneity_entropy"] = H
    adata.obs[f"{out_prefix}_dispersion_mad_zwithin_{groupby}"] = disp

    return adata


