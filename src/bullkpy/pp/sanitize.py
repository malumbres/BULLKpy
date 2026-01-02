from __future__ import annotations

import re
import numpy as np
import pandas as pd
import tempfile
from pathlib import Path
import anndata as ad

def _safe_filename(name: str, maxlen: int = 80) -> str:
    # Replace characters that can break paths in HDF5 / filesystem
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name))
    return s[:maxlen] if len(s) > maxlen else s


def _pick_informative_rows(s: pd.Series, n: int = 200) -> pd.Index:
    """
    Choose rows likely to expose type issues:
    - non-null rows
    - rows with "weird" python objects (list/dict/set/tuple/bytes)
    - rows with different python types
    Fallback to head.
    """
    idx = []

    # prioritize non-null
    s0 = s.dropna()
    if s0.empty:
        return s.index[: min(len(s), 1)]

    # rows with suspicious objects
    def is_suspicious(x):
        return isinstance(x, (list, dict, set, tuple, bytes, bytearray, memoryview))

    suspicious = s0[s0.map(is_suspicious)]
    if len(suspicious) > 0:
        idx.extend(list(suspicious.index[:n]))

    # rows covering distinct python types
    types = s0.map(lambda x: type(x).__name__)
    for t in types.unique()[:50]:
        idx.extend(list(types[types == t].index[:5]))

    # add a few random non-null
    if len(idx) < n and len(s0) > 1:
        remaining = s0.index.difference(pd.Index(idx))
        take = min(n - len(idx), len(remaining))
        if take > 0:
            # deterministic-ish
            idx.extend(list(remaining[:take]))

    # fallback head
    if not idx:
        idx = list(s0.index[: min(len(s0), n)])

    return pd.Index(idx).unique()


def sanitize_metadata(
    df: pd.DataFrame,
    *,
    index_col: str | None = None,
    numeric_min_frac: float = 0.9,     # if >=90% values parse as numeric → numeric
    category_max_unique: int = 50,     # if <=50 uniques → category (for strings)
    category_max_frac: float = 0.2,    # or if uniques <= 20% of rows → category
    datetime_min_frac: float = 0.9,    # if >=90% parse as datetime → datetime
    drop_high_cardinality_strings: bool = False,
    high_cardinality_frac: float = 0.5,  # if >50% uniques consider “ID-like”
    verbose: bool = True,
) -> pd.DataFrame:


    """
    Sanitize a metadata table before adding it to an AnnData object.

    This function inspects each column of a metadata DataFrame and attempts to
    convert it to an appropriate pandas dtype (numeric, datetime, category, or
    string), minimizing problematic `object` columns that can break `.h5ad` writing
    and downstream analyses.

    The goal is to make metadata:
    - safe to store in `adata.obs`
    - easy to analyze (numeric vs categorical)
    - compatible with HDF5 serialization
    """

    df = df.copy()

    # 0) optional index
    if index_col is not None:
        if index_col not in df.columns:
            raise KeyError(f"index_col='{index_col}' not found in metadata columns.")
        df = df.set_index(index_col)

    # 1) normalize missing tokens
    na_tokens = {"", "NA", "N/A", "na", "n/a", "NaN", "nan", "None", "none", "NULL", "null", ".", "-"}
    for c in df.columns:
        if df[c].dtype == object:
            s = df[c].astype(str)
            m = s.isin(na_tokens)
            if m.any():
                df.loc[m, c] = np.nan

    # 2) per-column inference
    report = []
    n = len(df)

    for c in df.columns:
        s = df[c]

        # keep existing numeric/bool/category
        if pd.api.types.is_bool_dtype(s) or pd.api.types.is_numeric_dtype(s):
            report.append((c, str(s.dtype), str(s.dtype), "kept"))
            continue

        # try datetime
        if s.dtype == object or pd.api.types.is_string_dtype(s):
            parsed_dt = pd.to_datetime(s, errors="coerce", infer_datetime_format=True)
            frac_dt = parsed_dt.notna().mean()
            if frac_dt >= datetime_min_frac and parsed_dt.notna().sum() >= 2:
                df[c] = parsed_dt
                report.append((c, "object", "datetime64[ns]", f"datetime ({frac_dt:.2f})"))
                continue

        # try numeric
        if s.dtype == object or pd.api.types.is_string_dtype(s) or pd.api.types.is_categorical_dtype(s):
            parsed_num = pd.to_numeric(s, errors="coerce")
            frac_num = parsed_num.notna().mean()
            if frac_num >= numeric_min_frac and parsed_num.notna().sum() >= 2:
                df[c] = parsed_num.astype(float)
                report.append((c, "object", "float", f"numeric ({frac_num:.2f})"))
                continue

        # categories vs strings
        # (treat everything else as text-ish)
        s_str = s.astype("string")
        n_unique = s_str.nunique(dropna=True)
        frac_unique = (n_unique / n) if n > 0 else 0.0

        # high-cardinality “ID-like”
        if drop_high_cardinality_strings and frac_unique >= high_cardinality_frac and n_unique >= 100:
            df.drop(columns=[c], inplace=True)
            report.append((c, "object", "dropped", f"high-cardinality ({n_unique})"))
            continue

        # otherwise, category if low-cardinality
        if (n_unique <= category_max_unique) or (frac_unique <= category_max_frac):
            df[c] = pd.Categorical(s_str)
            report.append((c, "object", "category", f"category ({n_unique})"))
        else:
            # keep as pandas string (saves better than python object)
            df[c] = s_str
            report.append((c, "object", "string", f"string ({n_unique})"))

    if verbose:
        rep = pd.DataFrame(report, columns=["col", "from", "to", "rule"]).sort_values("to")
        print(rep.to_string(index=False))

    return df



def find_bad_obs_cols_by_write(
    adata: ad.AnnData,
    *,
    n_rows: int = 200,
    include_index_test: bool = True,
) -> tuple[list[tuple[str, str]], str | None]:
    """
    Try writing each obs column into a tiny AnnData to detect which columns break .write().

    Key improvements vs your version:
    - tests many rows, not only the first
    - selects "informative" rows likely to expose mixed / non-string objects
    - uses safe filenames
    - keeps obs_names aligned
    """
    bad: list[tuple[str, str]] = []

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # Small X with enough rows
        n_total = adata.n_obs
        n_use = min(int(n_rows), n_total)
        Xsmall = np.zeros((n_use, 1), dtype=np.float32)

        # --------- index test (obs_names) ----------
        if include_index_test:
            try:
                tmp = ad.AnnData(X=Xsmall.copy())
                tmp.obs = pd.DataFrame(index=adata.obs_names[:n_use].copy())
                tmp.write(td / "ok_obs_index.h5ad")
            except Exception as e:
                return [("<obs_names/index>", str(e))], str(e)

        # --------- per-column tests ----------
        for col in adata.obs.columns:
            s = adata.obs[col]

            # choose rows likely to expose problems
            idx = _pick_informative_rows(s, n=n_use)
            # Ensure we always have at least 1 row
            if len(idx) == 0:
                idx = adata.obs_names[:1]

            try:
                tmp = ad.AnnData(X=np.zeros((len(idx), 1), dtype=np.float32))
                tmp.obs_names = pd.Index(idx.astype(str))  # safe
                tmp.obs = pd.DataFrame({col: s.loc[idx].to_numpy()}, index=tmp.obs_names)

                fname = f"ok_{_safe_filename(col)}.h5ad"
                tmp.write(td / fname)

            except Exception as e:
                bad.append((col, str(e)))

    return bad, None


def make_obs_h5ad_safe(adata, cols=None):
    if cols is None:
        cols = list(adata.obs.columns)

    # always safe
    adata.obs.columns = adata.obs.columns.map(str)
    adata.obs_names = adata.obs_names.map(str)

    for c in cols:
        s = adata.obs[c]

        # categorical -> force categories to strings
        if pd.api.types.is_categorical_dtype(s):
            adata.obs[c] = s.astype(str).astype("category")
            continue

        # object -> stringify elementwise
        if s.dtype == "object":
            adata.obs[c] = s.map(lambda x: "" if x is None or (isinstance(x, float) and np.isnan(x)) else str(x))

    return adata

