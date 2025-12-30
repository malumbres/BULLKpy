from __future__ import annotations

from typing import Literal
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import anndata as ad
import warnings

from ..logging import info
from ._nb_utils import deseq2_size_factors, estimate_dispersion_mom, shrink_dispersion_to_trend


def de_glm(
    adata: ad.AnnData,
    *,
    formula: str,
    contrast: tuple[str, str, str] | None = None,
    # contrast = ("Subtype", "Basal", "Luminal") meaning Basal - Luminal
    layer_counts: str = "counts",
    shrink_dispersion: bool = True,
    prior_df: float = 10.0,
    key_added: str = "de_glm",
    add_intercept: bool = True,
) -> None:
    """
    NB-GLM DE with design matrix from a formula and a contrast.

    Example:
      bk.tl.de_glm(adata, formula="~ Subtype + Batch", contrast=("Subtype","Basal","Luminal"))

    Stores:
      adata.uns[key_added][contrast_name] = dict(method/params/results)
    """
    try:
        import statsmodels.api as sm
    except Exception as e:
        raise ImportError(
            "de_glm requires statsmodels. Install with:\n  pip install statsmodels\n"
            f"Original error: {e}"
        )

    if layer_counts not in adata.layers:
        raise KeyError(f"layer_counts='{layer_counts}' not found in adata.layers")

    info(f"DE GLM (NB) with formula={formula!r}, layer_counts={layer_counts!r}")

    # Build design matrix
    X, coef_names = _design_matrix_from_formula(adata.obs, formula, add_intercept=add_intercept)

    # Determine which coefficient to test (simple main-effect contrast)
    if contrast is None:
        raise ValueError("contrast is required for now, e.g. ('Subtype','Basal','Luminal').")

    var, level_a, level_b = contrast
    c = _contrast_vector(coef_names, var, level_a, level_b)

    # Counts matrix
    Y = adata.layers[layer_counts]
    if sp.issparse(Y):
        Y = Y.toarray()
    Y = np.asarray(Y, dtype=float)

    # Size factors (DESeq2-style)
    sf = deseq2_size_factors(Y)
    offset = np.log(sf)

    # Dispersion from normalized counts (MoM + optional shrink)
    Yn = Y / sf[:, None]
    mu, alpha_hat = estimate_dispersion_mom(Yn)
    alpha_use = shrink_dispersion_to_trend(mu, alpha_hat, prior_df=prior_df) if shrink_dispersion else alpha_hat

    # Fit gene-wise NB GLM and Wald test for c'beta = 0
    n_genes = Y.shape[1]
    beta = np.zeros((n_genes, X.shape[1]), dtype=float)
    se_c = np.full(n_genes, np.nan, dtype=float)
    wald = np.full(n_genes, np.nan, dtype=float)
    pval = np.ones(n_genes, dtype=float)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Negative binomial dispersion parameter alpha not set.*")

        for j in range(n_genes):
            y = Y[:, j]
            if y.sum() == 0:
                continue

            a = float(alpha_use[j])

            try:
                model = sm.GLM(
                    y,
                    X,
                    family=sm.families.NegativeBinomial(alpha=a),
                    offset=offset,
                )
                fit = model.fit(maxiter=100, disp=0)

                b = fit.params
                cov = fit.cov_params()

                # Wald for contrast: (c^T b) / sqrt(c^T Cov c)
                num = float(c @ b)
                den = float(np.sqrt(c @ cov @ c))
                if den <= 0 or not np.isfinite(den):
                    continue

                z = num / den
                p = 2 * (1 - stats.norm.cdf(abs(z)))

                beta[j, :] = b
                se_c[j] = den
                wald[j] = z
                pval[j] = p

            except Exception:
                continue

    qval = _bh_fdr(pval)

    # Approximate log2FC for the contrast (natural log -> log2)
    log2fc = (beta @ c) / np.log(2)

    contrast_name = f"{var}:{level_a}_vs_{level_b}"
    res = pd.DataFrame(
        {
            "gene": adata.var_names.to_numpy(),
            "log2FC": log2fc,
            "wald_z": wald,
            "pval": pval,
            "qval": qval,
            "dispersion": alpha_use,
            "mean_norm": mu,
        }
    ).sort_values("qval", ascending=True)

    adata.uns.setdefault(key_added, {})
    adata.uns[key_added][contrast_name] = {
        "method": "nb_glm_deseq2_like",
        "formula": formula,
        "contrast": {"var": var, "level_a": level_a, "level_b": level_b},
        "coef_names": coef_names,
        "layer_counts": layer_counts,
        "shrink_dispersion": shrink_dispersion,
        "prior_df": float(prior_df),
        "results": res,
    }

    info(f"Stored DE GLM results in adata.uns['{key_added}']['{contrast_name}'].")


def _design_matrix_from_formula(obs: pd.DataFrame, formula: str, *, add_intercept: bool = True):
    """
    Minimal formula parser for: "~ A + B + C"
    Treats columns as categorical if dtype is category/object, else numeric.
    """
    f = formula.strip()
    if not f.startswith("~"):
        raise ValueError("formula must start with '~', e.g. '~ Subtype + Batch'")

    terms = [t.strip() for t in f[1:].split("+") if t.strip()]
    if not terms:
        raise ValueError("formula has no terms")

    X_parts = []
    names = []

    if add_intercept:
        X_parts.append(np.ones((obs.shape[0], 1), dtype=float))
        names.append("Intercept")

    for t in terms:
        if t not in obs.columns:
            raise KeyError(f"'{t}' not found in adata.obs")

        s = obs[t]
        is_cat = (str(s.dtype) == "category") or (s.dtype == object)

        if is_cat:
            cats = s.astype("category")
            # one-hot encode, drop first level as reference
            levels = list(cats.cat.categories)
            ref = levels[0]
            for lev in levels[1:]:
                col = (cats == lev).to_numpy(dtype=float).reshape(-1, 1)
                X_parts.append(col)
                names.append(f"{t}[{lev}]")
        else:
            col = pd.to_numeric(s).to_numpy(dtype=float).reshape(-1, 1)
            X_parts.append(col)
            names.append(t)

    X = np.concatenate(X_parts, axis=1)
    return X, names


def _contrast_vector(coef_names: list[str], var: str, level_a: str, level_b: str) -> np.ndarray:
    """
    For categorical var with dummy columns var[level], reference is first category (not in coef_names).
    Contrast returns (level_a - level_b) in coefficient space.

    If one of the levels is the reference, it is handled automatically.
    """
    # Detect dummy name format
    def name_for(level: str) -> str:
        return f"{var}[{level}]"

    c = np.zeros(len(coef_names), dtype=float)

    a_name = name_for(level_a)
    b_name = name_for(level_b)

    # If level is reference (dropped), it has no coefficient; effect is 0 in that dimension
    if a_name in coef_names:
        c[coef_names.index(a_name)] += 1.0
    if b_name in coef_names:
        c[coef_names.index(b_name)] -= 1.0

    if not np.any(c):
        # both are reference? or var not categorical as expected
        raise ValueError(
            f"Could not build contrast for {var} {level_a} vs {level_b}. "
            f"Available coefficients: {coef_names}"
        )
    return c


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out