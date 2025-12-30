def adjusted_rand_index(
    adata,
    *,
    true_key: str,
    pred_key: str,
):
    mask = (
        adata.obs[true_key].notna() &
        adata.obs[pred_key].notna()
    )

    if mask.sum() == 0:
        raise ValueError("No samples with both labels present")

    return adjusted_rand_score(
        adata.obs.loc[mask, true_key],
        adata.obs.loc[mask, pred_key],
    )