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


def _missing_keys(required_df, available_df, cols):
    """Set of unique ``cols`` tuples present in ``required_df`` but absent from ``available_df``."""
    cols = list(cols)
    required_keys = set(map(tuple, required_df[cols].drop_duplicates().to_numpy()))
    available_keys = set(map(tuple, available_df[cols].drop_duplicates().to_numpy()))
    return required_keys - available_keys


def _assert_tracks_covered_by_csv(
    h5ad_df, csv_df, *, groupby_cols, h5ad_label, available_label="the filtered track-features CSV",
):
    """Hard-error if any track key in ``h5ad_df`` (already ``_normalize_id_column``-normalized
    on ``groupby_cols``) is entirely absent from ``csv_df`` (same normalization).

    A track missing from ``available_label`` altogether means the data sources no longer match —
    almost always because one was regenerated (e.g. filtering re-run) after the other was created.
    """
    missing = _missing_keys(h5ad_df, csv_df, groupby_cols)
    if missing:
        example = sorted(missing)[:10]
        raise ValueError(
            f"{len(missing)} track(s) present in {h5ad_label} are missing entirely from "
            f"{available_label} (matched on {list(groupby_cols)}). The underlying data no longer "
            f"matches {h5ad_label} — this almost always means the input data was regenerated "
            f"(e.g. Filtering re-run) after the behavioral analysis was run. Re-run the "
            f"Behavioral Analysis step to regenerate {h5ad_label} from the current filtered data "
            f"before running this analysis. Missing example track key(s) (first 10, "
            f"{list(groupby_cols)}): {example}"
        )


def _assert_tracks_have_window_data(windows_df, matched_df, *, key_cols, h5ad_label):
    """Hard-error if any track/window key in ``windows_df`` (the full required set, e.g.
    ``adata_tracks.obs``) has zero rows in ``matched_df`` (the CSV already restricted to that
    track's ``[position_t_min, position_t_max]`` window).

    ``position_t_min``/``position_t_max`` are always derived as the min/max of real ``position_t``
    values from rows of this same CSV (only ever dropped upstream, never padded or resampled), so
    on genuinely in-sync data a track's own window bounds are guaranteed to have matching rows.
    Zero overlap therefore only happens when the CSV and h5ad are out of sync — never a legitimate
    "no contact recorded" result — so this must never be silently defaulted or dropped.
    """
    missing = _missing_keys(windows_df, matched_df, key_cols)
    if missing:
        example = sorted(missing)[:10]
        raise ValueError(
            f"{len(missing)} track(s)/window(s) in {h5ad_label} have no matching timepoints in "
            f"the filtered track-features CSV within their classified window (position_t_min–"
            f"position_t_max) (matched on {list(key_cols)}). Since {h5ad_label} is built directly "
            f"from this CSV, this should be impossible on in-sync data — it means the filtering/"
            f"tracking data no longer matches {h5ad_label} (most likely it was regenerated, e.g. "
            f"Filtering re-run, after the behavioral analysis was run). Re-run the Behavioral "
            f"Analysis step to regenerate {h5ad_label} from the current filtered data before "
            f"running this analysis. Missing example track key(s) (first 10, {list(key_cols)}): "
            f"{example}"
        )


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

    Every classified track (row of ``windows``) is required to have at least one matching
    timepoint in ``df_timepoints`` within its classified window — raises ``ValueError`` otherwise
    (see ``_assert_tracks_have_window_data``), since that can only happen when ``df_timepoints``
    no longer matches ``adata_tracks`` (never a legitimate zero-contact result).

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
    # timepoints actually fall inside it.
    matched = windows.merge(df, on=groupby_cols, how="inner")
    matched = matched[
        (matched[time_col] >= matched["position_t_min"]) & (matched[time_col] <= matched["position_t_max"])
    ]

    _assert_tracks_have_window_data(
        windows, matched, key_cols=key_cols,
        h5ad_label="adata_tracks.obs (track classification h5ad)",
    )

    stats = matched.groupby(key_cols, sort=False, observed=True)[contact_col].agg(
        **{mean_col: "mean", max_bout_col: lambda s: _max_true_run_length(s.tolist())}
    ).reset_index()

    # The assert above guarantees every window key has >=1 row in `matched`, so this merge can
    # never leave a NaN row.
    out = windows[key_cols].merge(stats, on=key_cols, how="left")
    out[max_bout_col] = out[max_bout_col].astype(int)
    out[group_col] = np.where(out[max_bout_col] >= min_bout_length, "contact", "no_contact")
    out[group_col] = pd.Categorical(out[group_col], categories=["no_contact", "contact"])

    if verbose:
        n_contact = int((out[group_col] == "contact").sum())
        n_no_contact = int((out[group_col] == "no_contact").sum())
        print(f"Contact grouping on '{contact_col}': {n_contact} contact / {n_no_contact} no_contact tracks.")
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
