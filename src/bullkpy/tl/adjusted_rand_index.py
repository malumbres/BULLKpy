from __future__ import annotations

from sklearn.metrics import adjusted_rand_score


def adjusted_rand_index(
    adata,
    *,
    true_key: str,
    pred_key: str,
):
    """Adjusted Rand Index between two categorical ``.obs`` columns.

    Samples with a missing value in either column are ignored.

    Parameters
    ----------
    adata
        Annotated data matrix.
    true_key
        Column in ``adata.obs`` holding the reference labels.
    pred_key
        Column in ``adata.obs`` holding the predicted/compared labels.

    Returns
    -------
    float
        Adjusted Rand Index in ``[-1, 1]``; 1.0 means identical partitions and
        values around 0 mean agreement no better than chance.
    """
    for key in (true_key, pred_key):
        if key not in adata.obs:
            raise KeyError(f"{key!r} not found in adata.obs")

    mask = (
        adata.obs[true_key].notna() &
        adata.obs[pred_key].notna()
    )

    if mask.sum() == 0:
        raise ValueError("No samples with both labels present")

    return float(adjusted_rand_score(
        adata.obs.loc[mask, true_key],
        adata.obs.loc[mask, pred_key],
    ))
