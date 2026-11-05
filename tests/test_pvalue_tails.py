"""Very strong effects must give tiny but non-zero p-values (2*sf(|z|), not 2*(1-cdf(|z|)), which is 0 for |z| > 8)."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

import bullkpy as bk


def _strong_counts(n: int = 120, n_genes: int = 40) -> ad.AnnData:
    rng = np.random.default_rng(0)
    group = np.array(["A"] * (n // 2) + ["B"] * (n // 2))
    mu = np.full((n, n_genes), 50.0)
    mu[group == "B", 0] = 5000.0  # gene 0: a hundred-fold difference, |z| far beyond 8
    X = rng.poisson(mu).astype(np.int64)
    adata = ad.AnnData(
        X.astype(float),
        obs=pd.DataFrame({"grp": group}, index=[f"s{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )
    adata.layers["counts"] = X
    return adata


def test_nb_de_strong_effect_has_positive_pvalue():
    adata = _strong_counts()
    bk.tl.de(adata, groupby="grp", group="B", reference="A", method="nb_glm")
    r = bk.get.de(adata, "grp_B_vs_A").set_index("gene")
    assert 0 < r.loc["g0", "pval"] < 1e-30
    assert 0 < r.loc["g0", "qval"] < 1e-28


def test_de_glm_strong_effect_has_positive_pvalue():
    adata = _strong_counts()
    bk.tl.de_glm(adata, formula="~ grp", contrast=("grp", "B", "A"))
    r = bk.get.de(adata, "grp:B_vs_A", key="de_glm").set_index("gene")
    assert 0 < r.loc["g0", "pval"] < 1e-30
    assert 0 < r.loc["g0", "qval"] < 1e-28


def test_cox_multivariate_statsmodels_strong_effect_has_positive_pvalue():
    rng = np.random.default_rng(1)
    n = 400
    x = rng.normal(size=n)
    obs = pd.DataFrame({"t": rng.exponential(np.exp(-x)), "e": 1}, index=[f"s{i}" for i in range(n)])
    adata = ad.AnnData(x.reshape(-1, 1), obs=obs, var=pd.DataFrame(index=["g"]))
    res = bk.tl.cox_multivariate(
        adata, genes=["g"], time_col="t", event_col="e", layer=None, backend="statsmodels", standardize=False
    )
    assert 0 < float(res["pval"].iloc[0]) < 1e-30  # z is about 14.5


def test_batch_meta_correlation_tail_is_exact():
    from bullkpy.tl.correlations import _corr_batch_aware

    rng = np.random.default_rng(2)
    n = 600
    x = rng.normal(size=n)
    y = x + 1.5 * rng.normal(size=n)  # r about 0.55, Fisher z about 15
    adata = ad.AnnData(np.zeros((n, 1)), obs=pd.DataFrame({"b": ["b1", "b2", "b3"] * (n // 3)}, index=[f"s{i}" for i in range(n)]))
    r, p, _ = _corr_batch_aware(x, y, adata=adata, method="pearson", batch_key="b", batch_mode="meta")
    assert r > 0.4 and 0 < p < 1e-30
