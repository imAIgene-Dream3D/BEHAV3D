import numpy as np
import pandas as pd
import scanpy as sc

def select_nonbinary_columnnames(df: pd.DataFrame, cols: list[str], tol: float = 1e-9):
    binary, continuous = [], []
    for c in cols:
        s = df[c].dropna().astype(float)
        if s.empty:
            continuous.append(c)
            continue
        u = np.unique(s.values)
        # binary if all unique values are (approximately) 0 or 1
        if np.all(np.isclose(u, 0.0, atol=tol) | np.isclose(u, 1.0, atol=tol)):
            binary.append(c)
        else:
            continuous.append(c)
    return continuous

def subset_full_tracks(
    df: pd.DataFrame,
    fraction: float = 0.1,
    id_cols: list[str] = ["sample_name", "TrackID"],
    random_state: int = 0,
    return_selected_keys: bool = False,
    sampled_keys=None,
):
    """
    Randomly subsample a fraction of unique tracks defined by id_cols,
    then return df restricted to those full tracks.

    If `selected_keys` is provided (a DataFrame containing columns `id_cols`),
    no random sampling is done; instead, df is restricted to the tracks
    matching those keys. This is useful for reproducing an earlier selection.
    """
    track_keys = df[id_cols].drop_duplicates()

    if sampled_keys is None:
        sampled_keys = track_keys.sample(frac=fraction, random_state=random_state)
    else:
        sampled_keys = sampled_keys[id_cols].drop_duplicates()

    mask = (
        df[id_cols]
        .merge(sampled_keys, on=id_cols, how="left", indicator=True)["_merge"]
        .eq("both")
    ).to_numpy()
    df_sub = df.loc[mask]

    if return_selected_keys:
        return df_sub, sampled_keys
    else:
        return df_sub

def subset_windowed_tracks(
    df_windowed, 
    step_size, 
    time_col="position_t",
    id_cols=("sample_name", "TrackID"), 
    keep="first"
    ):
    
    df_windowed = df_windowed.dropna()
    step_size = max(int(step_size), 1)
    if step_size == 1:
        return df_windowed

    # ensure stable order within each track
    dfw = df_windowed.sort_values(list(id_cols) + [time_col], kind="mergesort").copy()

    if keep == "first":
        rank = dfw.groupby(list(id_cols), sort=False).cumcount()
        mask = (rank % step_size) == 0
    elif keep == "last":
        rank_from_end = dfw.groupby(list(id_cols), sort=False).cumcount(ascending=False)
        mask = (rank_from_end % step_size) == 0
    else:
        raise ValueError("keep must be 'first' or 'last'")

    return dfw.loc[mask]

