"""io + pp: loading, QC, normalisation, filtering, sanitising."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import bullkpy as bk


# --------------------------------------------------------------------- io
def test_read_counts_orientation(counts_file, counts_df):
    adata = bk.io.read_counts(counts_file, orientation="genes_by_samples")
    assert adata.n_obs == counts_df.shape[1]  # samples
    assert adata.n_vars == counts_df.shape[0]  # genes
    assert list(adata.obs_names) == list(counts_df.columns)
    assert list(adata.var_names) == list(counts_df.index)


def test_read_counts_rejects_bad_orientation(counts_file):
    with pytest.raises(ValueError, match="orientation"):
        bk.io.read_counts(counts_file, orientation="nonsense")


def test_add_metadata_merges_on_index(adata_raw):
    for col in ("Subtype", "Batch", "age", "OS_time"):
        assert col in adata_raw.obs.columns
    assert adata_raw.obs["Subtype"].notna().all()


def test_add_metadata_rejects_missing_index_col(adata_raw, metadata_file):
    with pytest.raises(ValueError, match="index_col"):
        bk.io.add_metadata(adata_raw, metadata_file, index_col="not_a_column")


def test_add_metadata_reads_excel(tmp_path, adata_raw, counts_df):
    """Excel input previously crashed: read_excel has no `low_memory` parameter."""
    meta = pd.DataFrame(
        {"sample_id": list(counts_df.columns), "grp": ["a", "b"] * (counts_df.shape[1] // 2)}
    )
    p = tmp_path / "meta.xlsx"
    meta.to_excel(p, index=False)

    bk.io.add_metadata(adata_raw, p, index_col="sample_id")
    assert "grp" in adata_raw.obs.columns


# --------------------------------------------------------------------- pp
def test_set_raw_counts_and_layers(adata):
    assert "counts" in adata.layers
    assert "cpm" in adata.layers
    assert "log1p_cpm" in adata.layers


def test_qc_metrics_adds_expected_columns(adata):
    for col in ("total_counts", "n_genes_detected", "pct_counts_mt", "pct_counts_ribo"):
        assert col in adata.obs.columns, f"{col} missing from adata.obs"


def test_qc_metrics_percentages_are_in_range(adata):
    for col in ("pct_counts_mt", "pct_counts_ribo"):
        v = adata.obs[col].to_numpy(dtype=float)
        assert np.nanmin(v) >= 0.0
        assert np.nanmax(v) <= 100.0


def test_qc_metrics_honours_compute_flags(adata_raw):
    """compute_pct_mt=False used to be ignored."""
    bk.pp.set_raw_counts(adata_raw)
    bk.pp.qc_metrics(adata_raw, compute_pct_mt=False, compute_pct_ribo=False)
    assert "pct_counts_mt" not in adata_raw.obs.columns
    assert "pct_counts_ribo" not in adata_raw.obs.columns
    assert "total_counts" in adata_raw.obs.columns


def test_normalize_cpm_rows_sum_to_1e6(adata):
    cpm = np.asarray(adata.layers["cpm"])
    assert np.allclose(cpm.sum(axis=1), 1e6, rtol=1e-6)


def test_log1p_is_log_of_cpm(adata):
    cpm = np.asarray(adata.layers["cpm"])
    assert np.allclose(np.asarray(adata.layers["log1p_cpm"]), np.log1p(cpm))


def test_highly_variable_genes(adata):
    bk.pp.highly_variable_genes(adata, n_top_genes=50, layer="log1p_cpm")
    assert "highly_variable" in adata.var.columns
    assert adata.var["highly_variable"].sum() == 50


def test_filter_genes_reduces_features(adata):
    out = bk.pp.filter_genes(adata.copy(), min_samples=2)
    assert out.n_vars <= adata.n_vars
    assert out.n_obs == adata.n_obs


def test_filter_samples_keeps_features(adata):
    out = bk.pp.filter_samples(adata.copy())
    assert out.n_vars == adata.n_vars


def test_batch_correct_combat_writes_layer(adata):
    a = adata.copy()
    bk.pp.batch_correct_combat(a, batch_key="Batch", layer="log1p_cpm")
    assert "combat" in a.layers
    assert np.asarray(a.layers["combat"]).shape == a.shape


def test_sanitize_metadata_takes_a_dataframe(adata):
    out = bk.pp.sanitize_metadata(adata.obs.copy(), verbose=False)
    assert isinstance(out, pd.DataFrame)
    assert len(out) == adata.n_obs


def test_obs_map_categories(adata):
    a = adata.copy()
    bk.pp.obs_map_categories(a, key="Subtype", mapping={"Basal": "B"}, out_key="Subtype2")
    assert "B" in set(a.obs["Subtype2"].astype(str))


def test_h5ad_safety_helpers_roundtrip(tmp_path, adata):
    a = adata.copy()
    bk.pp.make_obs_h5ad_safe_strict(a)
    bk.pp.make_var_h5ad_safe_strict(a)
    p = tmp_path / "out.h5ad"
    a.write_h5ad(p)
    assert p.exists() and p.stat().st_size > 0


# ------------------------------------------------- read_counts truncation guard
def test_read_counts_does_not_truncate_float_input(tmp_path):
    """A float matrix must never be silently cast to int.

    UCSC Xena distributes log2(count+1) matrices under a *_counts.tsv name;
    casting those to int64 turned 10.77 into 10 and 0.9999 into 0.
    """
    df = pd.DataFrame(
        [[10.7698, 10.7211], [3.5, 0.9999], [0.0, 15.25]],
        index=["G1", "G2", "G3"],
        columns=["S1", "S2"],
    )
    p = tmp_path / "log2.tsv"
    df.to_csv(p, sep="\t")

    adata = bk.io.read_counts(p, orientation="genes_by_samples")
    X = np.asarray(adata.X)

    assert not np.issubdtype(X.dtype, np.integer), "float input was cast to an integer dtype"
    assert np.allclose(X, df.to_numpy().T), "values changed while reading"


def test_read_counts_still_casts_genuine_counts(tmp_path):
    df = pd.DataFrame([[10, 7], [3, 1], [0, 15]], index=["G1", "G2", "G3"], columns=["S1", "S2"])
    p = tmp_path / "counts.tsv"
    df.to_csv(p, sep="\t")

    adata = bk.io.read_counts(p, orientation="genes_by_samples")
    assert np.issubdtype(np.asarray(adata.X).dtype, np.integer)


def test_read_counts_dtype_none_never_casts(tmp_path):
    df = pd.DataFrame([[1.5, 2.5]], index=["G1"], columns=["S1", "S2"])
    p = tmp_path / "f.tsv"
    df.to_csv(p, sep="\t")

    adata = bk.io.read_counts(p, orientation="genes_by_samples", dtype=None)
    assert np.allclose(np.asarray(adata.X).ravel(), [1.5, 2.5])
