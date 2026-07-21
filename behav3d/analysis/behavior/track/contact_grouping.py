import numpy as np
import pandas as pd

from behav3d.features.state_descriptive_features import rle_encode

_CONTACT_COL_SUFFIXES = ("_contact", "_contact_on_distance")


def list_available_contact_columns(df_timepoints):
    """Return per-timepoint contact columns (``{celltype}_contact`` / ``..._on_distance``)
    present in a per-timepoint tracks dataframe, sorted for stable dropdown ordering."""
    if df_timepoints is None or not hasattr(df_timepoints, "columns"):
        return []
    cols = [
        str(c) for c in df_timepoints.columns
        if str(c).endswith(_CONTACT_COL_SUFFIXES) and not str(c).startswith("active_")
    ]
    return sorted(cols)


def _contact_group_col_name(contact_col):
    return f"{contact_col}_group"


def _contact_mean_col_name(contact_col):
    return f"{contact_col}_mean_fraction"


def _contact_max_bout_col_name(contact_col):
    return f"{contact_col}_max_bout_length"


def _max_true_run_length(values):
    """Longest contiguous run of truthy values, via the shared rle_encode helper."""
    runs = rle_encode(list(values))
    true_lengths = [length for value, length in runs if bool(value)]
    return int(max(true_lengths)) if true_lengths else 0


def _normalize_id_column(series):
    """Normalize a key column so common round-trip drift (e.g. int in an h5ad-derived obs vs.
    ``"5.0"`` after a CSV re-read) doesn't break equality joins."""
    numeric = pd.to_numeric(series, errors="coerce")
    is_whole = numeric.notna() & (numeric == numeric.round())
    out = series.astype(str)
    out = out.where(~is_whole, numeric.where(is_whole).astype("Int64").astype(str))
    return out


_WINDOW_COL = "trajectory_window_id"


def compute_track_contact_features(
    df_timepoints,
    adata_tracks,
    *,
    contact_col,
    min_bout_length,
    time_col="position_t",
    groupby_cols=("sample_name", "TrackID"),
    group_col=None,
    mean_col=None,
    max_bout_col=None,
    verbose=False,
):
    """Per-track contact summary, restricted to each track's classified time window.

    For every (sample_name, TrackID[, trajectory_window_id]) present in ``adata_tracks.obs``,
    restricts ``df_timepoints`` to that track's ``[position_t_min, position_t_max]`` window (the
    same truncated window used to build the track-classification features — see
    ``extract_descibing_track_state_features`` in
    ``behav3d/features/state_descriptive_features.py``), then computes:

    - ``{contact_col}_group``: "contact" if the longest contiguous run of ``contact_col`` truthy
      values within that window is >= ``min_bout_length``, else "no_contact".
    - ``{contact_col}_mean_fraction``: mean of ``contact_col`` (0-1) over the window.
    - ``{contact_col}_max_bout_length``: longest contiguous truthy run length (timepoints).

    Every classified track (row of ``windows``) is guaranteed exactly one output row — a track
    with no matching contact timepoints resolves to "no_contact" (max_bout=0) rather than being
    dropped, since ``adata_tracks``/``df_timepoints`` are expected to share the same underlying
    per-timepoint source data. With ``verbose=True``, any such zero-match tracks are logged, since
    they can indicate a stale CSV or track-ID drift even though they now resolve to a definite
    label.

    Returns a DataFrame indexed by ``groupby_cols`` (plus ``trajectory_window_id`` when present in
    ``adata_tracks.obs``) with those three columns.
    """
    groupby_cols = [str(c) for c in list(groupby_cols)]
    time_col = str(time_col)
    contact_col = str(contact_col)
    min_bout_length = int(min_bout_length)
    group_col = str(group_col) if group_col is not None else _contact_group_col_name(contact_col)
    mean_col = str(mean_col) if mean_col is not None else _contact_mean_col_name(contact_col)
    max_bout_col = str(max_bout_col) if max_bout_col is not None else _contact_max_bout_col_name(contact_col)

    has_window_col = _WINDOW_COL in adata_tracks.obs.columns and _WINDOW_COL not in groupby_cols
    key_cols = groupby_cols + ([_WINDOW_COL] if has_window_col else [])

    missing = [c for c in groupby_cols + [time_col, contact_col] if c not in df_timepoints.columns]
    if missing:
        raise KeyError(f"Missing required columns in df_timepoints: {missing}")

    required_obs_cols = key_cols + ["position_t_min", "position_t_max"]
    missing_obs = [c for c in required_obs_cols if c not in adata_tracks.obs.columns]
    if missing_obs:
        raise KeyError(f"Missing required columns in adata_tracks.obs: {missing_obs}")

    windows = adata_tracks.obs[required_obs_cols].drop_duplicates(subset=key_cols).copy()
    for col in groupby_cols:
        windows[col] = _normalize_id_column(windows[col])
    if has_window_col:
        windows[_WINDOW_COL] = windows[_WINDOW_COL].astype(str)

    df = df_timepoints[groupby_cols + [time_col, contact_col]].copy()
    for col in groupby_cols:
        df[col] = _normalize_id_column(df[col])
    df[contact_col] = pd.to_numeric(df[contact_col], errors="coerce").fillna(0).astype(bool)

    # Inner-join the classified windows onto the raw contact timepoints to find, per window, which
    # timepoints actually fall inside it; windows with no match simply won't appear here — they're
    # restored below via a left-merge against the full `windows` key set, so no classified track is
    # ever silently dropped (unlike the previous inner-join-only approach).
    matched = windows.merge(df, on=groupby_cols, how="inner")
    matched = matched[
        (matched[time_col] >= matched["position_t_min"]) & (matched[time_col] <= matched["position_t_max"])
    ]

    if len(matched) > 0:
        stats = matched.groupby(key_cols, sort=False, observed=True)[contact_col].agg(
            **{mean_col: "mean", max_bout_col: lambda s: _max_true_run_length(s.tolist())}
        ).reset_index()
    else:
        stats = pd.DataFrame(columns=key_cols + [mean_col, max_bout_col])

    out = windows[key_cols].merge(stats, on=key_cols, how="left")

    unmatched_mask = out[max_bout_col].isna()
    n_unmatched = int(unmatched_mask.sum())
    out[mean_col] = out[mean_col].fillna(0.0)
    out[max_bout_col] = out[max_bout_col].fillna(0).astype(int)
    out[group_col] = np.where(out[max_bout_col] >= min_bout_length, "contact", "no_contact")

    if verbose and n_unmatched:
        unmatched_keys = out.loc[unmatched_mask, key_cols].to_records(index=False).tolist()[:10]
        print(
            f"Contact grouping on '{contact_col}': {n_unmatched} classified track(s) had no "
            f"matching timepoints in df_timepoints and were labeled 'no_contact' by default "
            f"(possible stale CSV or track-ID mismatch), e.g. {unmatched_keys}"
            + (" ..." if n_unmatched > 10 else "")
        )

    out[group_col] = pd.Categorical(out[group_col], categories=["no_contact", "contact"])
    return out.set_index(key_cols)[[mean_col, max_bout_col, group_col]]


def merge_track_contact_features_into_obs(
    adata_tracks,
    contact_features_df,
    *,
    contact_col,
    min_bout_length,
    groupby_cols=("sample_name", "TrackID"),
):
    """Left-merge ``compute_track_contact_features`` output onto ``adata_tracks.obs`` (mutated
    in place) and record provenance under ``adata_tracks.uns["contact_grouping"]``.

    Raises if any ``obs`` row ends up without a group label after the merge — with
    ``compute_track_contact_features`` guaranteeing full coverage of its input windows, that
    should be unreachable; if it fires, it means ``contact_features_df`` wasn't actually computed
    from (a superset of) this ``adata_tracks``.
    """
    groupby_cols = [str(c) for c in list(groupby_cols)]
    # Match compute_track_contact_features's auto-detected window key so split-track windows
    # (multiple obs rows per TrackID, distinguished by trajectory_window_id) aren't collapsed by
    # drop_duplicates below.
    has_window_col = (
        _WINDOW_COL in adata_tracks.obs.columns
        and _WINDOW_COL in contact_features_df.index.names
        and _WINDOW_COL not in groupby_cols
    )
    key_cols = groupby_cols + ([_WINDOW_COL] if has_window_col else [])

    obs = adata_tracks.obs.copy()
    orig_index = obs.index
    for col in groupby_cols:
        obs[col] = obs[col].astype(str)
    if has_window_col:
        obs[_WINDOW_COL] = obs[_WINDOW_COL].astype(str)

    right = contact_features_df.reset_index()
    for col in groupby_cols:
        right[col] = right[col].astype(str)
    if has_window_col:
        right[_WINDOW_COL] = right[_WINDOW_COL].astype(str)
    new_cols = [c for c in right.columns if c not in key_cols]
    right = right[key_cols + new_cols].drop_duplicates(subset=key_cols)

    # Idempotent: drop any columns from a previous call (e.g. re-running with a different
    # threshold) so the merge below can't produce '_x'/'_y'-suffixed duplicates.
    obs = obs.drop(columns=[c for c in new_cols if c in obs.columns])
    merged = obs.merge(right, on=key_cols, how="left")
    merged.index = orig_index

    group_col = _contact_group_col_name(contact_col)
    unlabeled = merged[group_col].isna()
    if unlabeled.any():
        bad_keys = merged.loc[unlabeled, key_cols].drop_duplicates().to_records(index=False).tolist()
        raise ValueError(
            f"{int(unlabeled.sum())} obs row(s) have no '{group_col}' after merging contact "
            f"features — contact_features_df must be computed from a superset of this "
            f"adata_tracks. Unmatched keys (first 10): {bad_keys[:10]}"
        )

    adata_tracks.obs = merged

    contact_grouping_meta = dict(adata_tracks.uns.get("contact_grouping", {}))
    contact_grouping_meta[str(contact_col)] = {
        "min_bout_length": int(min_bout_length),
        "group_col": group_col,
        "mean_col": _contact_mean_col_name(contact_col),
        "max_bout_col": _contact_max_bout_col_name(contact_col),
    }
    adata_tracks.uns["contact_grouping"] = contact_grouping_meta
    return adata_tracks
