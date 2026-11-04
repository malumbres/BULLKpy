"""Numeric fingerprints of the analysis functions, used to guard refactors.

``tests/golden/golden.json`` holds values produced by the code at the time Phase 2
started. A refactor that is meant to keep results must reproduce them. To accept an
intentional change, regenerate with ``python tests/golden_cases.py`` and explain it in
the CHANGELOG under "Results may change".
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GOLDEN = Path(__file__).parent / "golden" / "golden.json"
LAYER = "log1p_cpm"


def _vec(x, n: int | None = None) -> list[float]:
    a = np.asarray(x, dtype=float).ravel()
    if n:
        a = a[:n]
    return [None if not np.isfinite(v) else round(float(v), 8) for v in a]


def _table(df: pd.DataFrame, cols: list[str], key: str = "gene", n: int = 15) -> dict:
    df = df.sort_values(key).head(n)
    return {c: _vec(df[c]) for c in cols if c in df.columns}


def compute(adata, genes: list[str]) -> dict:
    import bullkpy as bk

    out: dict = {}

    bk.tl.pca(adata, n_comps=10)
    out["pca_variance_ratio"] = _vec(adata.uns["pca"]["variance_ratio"])
    out["pca_abs_scores_pc1"] = _vec(np.abs(adata.obsm["X_pca"][:, 0]), 20)

    bk.tl.neighbors(adata, n_neighbors=8)
    bk.tl.cluster(adata, method="kmeans", n_clusters=3, use_rep="X_pca", n_pcs=5, random_state=0)
    labels = adata.obs["clusters"].astype(str).to_numpy()
    out["kmeans_partition_sizes"] = sorted(int((labels == u).sum()) for u in np.unique(labels))

    bk.tl.de(adata, groupby="Subtype", group="Basal", reference="Luminal", method="welch_ttest")
    res = adata.uns["de"]["Subtype_Basal_vs_Luminal"]["results"]
    out["de_welch"] = _table(res, ["log2FC", "t", "pval", "qval"])

    bk.tl.de_glm(adata, formula="~ Subtype", contrast=("Subtype", "Basal", "Luminal"))
    res = adata.uns["de_glm"]["Subtype:Basal_vs_Luminal"]["results"]
    out["de_glm"] = _table(res, ["log2FC", "wald_z", "pval", "qval"])

    res = bk.tl.rank_genes_groups_fast(adata, groupby="Subtype", group="Basal", layer=LAYER)
    out["rank_genes_fast"] = _table(res, ["pval", "qval"])

    res = bk.tl.gene_categorical_association(adata, groupby="Subtype", genes=genes, layer=LAYER)
    out["gene_categorical_association"] = _table(res, ["pval", "qval"])

    res = bk.tl.gene_gene_correlations(adata, gene=genes[0], genes=genes, layer=LAYER)
    out["gene_gene_correlations"] = _table(res, [c for c in res.columns if c in ("r", "rho", "pval", "qval")])

    bk.tl.score_genes(adata, genes, score_name="golden_score", layer=LAYER)
    out["score_genes"] = _vec(adata.obs["golden_score"], 20)

    res = bk.tl.cox_gene_association(
        adata, time_col="OS_time", event_col="OS_event", genes=genes, layer=LAYER
    )
    out["cox_gene_association"] = _table(res, ["pval", "qval"])

    return out


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    import conftest  # noqa: F401  (registers the fixtures module for reuse)
    import bullkpy as bk

    rng = np.random.default_rng(0)
    from conftest import N_GENES, N_SAMPLES

    genes_all = [f"GENE{i:04d}" for i in range(N_GENES)]
    genes_all[:5] = ["MT-CO1", "MT-ND1", "MT-ATP6", "RPS3", "RPL7"]
    base = rng.lognormal(mean=4.0, sigma=1.0, size=(N_GENES, 1))
    lib = rng.uniform(0.6, 1.6, size=(1, N_SAMPLES))
    X = rng.poisson(base * lib).astype(np.int64)
    half = N_SAMPLES // 2
    X[:40, :half] = (X[:40, :half] * 3.5).astype(np.int64)
    samples = [f"S{i:03d}" for i in range(N_SAMPLES)]
    r2 = np.random.default_rng(1)
    meta = pd.DataFrame(
        {
            "Subtype": ["Basal"] * half + ["Luminal"] * (N_SAMPLES - half),
            "Batch": r2.choice(["b1", "b2"], size=N_SAMPLES),
            "Stage": r2.choice(["I", "II", "III"], size=N_SAMPLES),
            "age": r2.integers(35, 85, size=N_SAMPLES),
            "purity": r2.uniform(0.3, 0.95, size=N_SAMPLES),
            "OS_time": r2.uniform(1.0, 100.0, size=N_SAMPLES),
            "OS_event": r2.integers(0, 2, size=N_SAMPLES),
        },
        index=samples,
    )
    import anndata as ad

    a = ad.AnnData(X.T.astype(float), obs=meta, var=pd.DataFrame(index=genes_all))
    a.X = X.T.astype("int64")
    bk.pp.set_raw_counts(a)
    bk.pp.qc_metrics(a)
    bk.pp.normalize_cpm(a)
    bk.pp.log1p(a)
    GOLDEN.parent.mkdir(exist_ok=True)
    GOLDEN.write_text(json.dumps(compute(a, list(a.var_names[:10])), indent=1))
    print("wrote", GOLDEN)
