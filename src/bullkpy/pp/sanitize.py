from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import re, tempfile
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import tempfile

from dataclasses import dataclass

@dataclass
class ObsColumnIssue:
    column: str
    dtype: str
    issue: str



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
            parsed_dt = pd.to_datetime(s, errors="coerce")
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



def find_bad_obs_columns(adata: AnnData) -> list[ObsColumnIssue]:
    """
    Identify `.obs` columns likely to break `adata.write_h5ad()`,
    especially Pandas nullable StringArray ('string[python]').

    Returns a list of ObsColumnIssue.
    """
    issues: list[ObsColumnIssue] = []
    for col in adata.obs.columns:
        s = adata.obs[col]
        dt = str(s.dtype)

        # Pandas StringArray / nullable string dtype
        if pd.api.types.is_string_dtype(s.dtype) and not pd.api.types.is_object_dtype(s.dtype):
            issues.append(ObsColumnIssue(col, dt, "nullable pandas string dtype (StringArray)"))

        # Categorical with nullable string categories can also trigger issues
        if pd.api.types.is_categorical_dtype(s.dtype):
            cats = s.cat.categories
            if pd.api.types.is_string_dtype(cats.dtype) and not pd.api.types.is_object_dtype(cats.dtype):
                issues.append(ObsColumnIssue(col, f"{dt} (cats: {cats.dtype})", "categorical has StringArray categories"))

        # Object columns with mixed non-serializable stuff
        if pd.api.types.is_object_dtype(s.dtype):
            # quick heuristic: if it contains lists/dicts/sets, h5ad will likely fail
            bad = s.dropna().head(200)
            if bad.apply(lambda x: isinstance(x, (list, dict, set, tuple))).any():
                issues.append(ObsColumnIssue(col, dt, "object column contains list/dict/set/tuple"))
    return issues




#####################################################
### SANITIZE adata.obs
#####################################################


def _safe_filename(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))
    return s[:150] if len(s) > 150 else s

def _pick_informative_rows(s: pd.Series, n: int) -> pd.Index:
    nonmiss = s.dropna()
    if nonmiss.empty:
        return s.index[: min(n, len(s))]
    take = []
    take.extend(nonmiss.index[: min(n // 3, len(nonmiss))].tolist())
    take.extend(nonmiss.index[-min(n // 3, len(nonmiss)):].tolist())
    remain = [i for i in nonmiss.index.tolist() if i not in set(take)]
    if remain and len(take) < n:
        rng = np.random.default_rng(0)
        extra = rng.choice(remain, size=min(n - len(take), len(remain)), replace=False)
        take.extend(extra.tolist())
    return pd.Index(take)

def find_bad_obs_cols_by_write(
    adata,
    *,
    n_rows: int = 3000,
    include_index_test: bool = True,
):
    import anndata as anndata  # local import avoids notebook shadowing

    bad = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        n_use = min(int(n_rows), adata.n_obs)
        Xsmall = np.zeros((n_use, 1), dtype=np.float32)

        if include_index_test:
            try:
                tmp = anndata.AnnData(X=Xsmall.copy())
                tmp.obs = pd.DataFrame(index=pd.Index(adata.obs_names[:n_use].astype(str)))
                tmp.write(td / "ok_obs_index.h5ad")
            except Exception as e:
                return [("<obs_names/index>", str(e))], str(e)

        for col in adata.obs.columns:
            s = adata.obs[col]
            idx = _pick_informative_rows(s, n=n_use)
            if len(idx) == 0:
                idx = adata.obs_names[:1]

            try:
                tmp = anndata.AnnData(X=np.zeros((len(idx), 1), dtype=np.float32))
                tmp.obs_names = pd.Index(idx.astype(str))
                tmp.obs = pd.DataFrame({col: s.loc[idx].to_numpy()}, index=tmp.obs_names)
                tmp.write(td / f"ok_{_safe_filename(col)}.h5ad")
            except Exception as e:
                bad.append((col, str(e)))

    return bad, None


def make_obs_h5ad_safe_strict(
    adata,
    *,
    columns: list[str] | None = None,
    copy: bool = True,
    numeric_coerce_min_frac: float = 0.85,   # if >=85% of non-missing parses -> numeric
    stringify_complex_objects: bool = True,
    convert_nullable_string_to_object: bool = True,
    convert_categorical_string_categories_to_object: bool = True,
):
    """
    Returns (adata_safe, report_dict).
    """
    adata2 = adata.copy() if copy else adata
    obs = adata2.obs.copy()

    cols = columns if columns is not None else list(obs.columns)
    rep = {"numeric_coerced": [], "string_to_object": [], "stringified": [], "category_fixed": [], "errors": []}

    def _is_missing(x) -> bool:
        if x is None or x is pd.NA:
            return True
        if isinstance(x, float) and np.isnan(x):
            return True
        return False

    def _is_complex(x) -> bool:
        return isinstance(x, (list, dict, set, tuple))

    for col in cols:
        if col not in obs.columns:
            continue

        s = obs[col]
        dt = s.dtype

        # 1) pandas nullable string (string[python]/string[pyarrow]) -> object
        if convert_nullable_string_to_object and pd.api.types.is_string_dtype(dt) and not pd.api.types.is_object_dtype(dt):
            try:
                obs[col] = s.astype("object")
                rep["string_to_object"].append(col)
                s = obs[col]
                dt = s.dtype
            except Exception as e:
                rep["errors"].append((col, f"string->object failed: {e}"))
                continue

        # 2) fix categorical categories dtype if they are nullable strings
        if pd.api.types.is_categorical_dtype(dt) and convert_categorical_string_categories_to_object:
            try:
                cats = s.cat.categories
                if pd.api.types.is_string_dtype(cats.dtype) and not pd.api.types.is_object_dtype(cats.dtype):
                    obs[col] = pd.Categorical(
                        s.astype("object"),
                        categories=cats.astype("object"),
                        ordered=s.cat.ordered,
                    )
                    rep["category_fixed"].append(col)
                    s = obs[col]
                    dt = s.dtype
            except Exception as e:
                rep["errors"].append((col, f"categorical fix failed: {e}"))
                continue

        # 3) object dtype: decide numeric vs stringify
        if pd.api.types.is_object_dtype(obs[col].dtype):
            sample = obs[col].dropna()
            if sample.empty:
                continue

            # 3a) stringify complex objects (lists/dicts/tuples/sets)
            if stringify_complex_objects and sample.head(500).map(_is_complex).any():
                obs[col] = obs[col].map(lambda x: str(x) if _is_complex(x) else x).astype("object")
                rep["stringified"].append(col)
                continue

            # 3b) try numeric coercion if most values parse
            # (this fixes your albumin/days_to/etc columns which should be numeric)
            coerced = pd.to_numeric(obs[col], errors="coerce")
            nonmiss = obs[col].map(lambda x: not _is_missing(x)).sum()
            parsed = coerced.notna().sum()
            frac = (parsed / nonmiss) if nonmiss else 1.0

            if frac >= float(numeric_coerce_min_frac):
                obs[col] = coerced
                rep["numeric_coerced"].append(col)
                continue

            # 3c) otherwise stringify all non-missing non-strings
            obs[col] = obs[col].map(lambda x: x if _is_missing(x) or isinstance(x, str) else str(x)).astype("object")
            rep["stringified"].append(col)

    adata2.obs = obs
    return adata2, rep