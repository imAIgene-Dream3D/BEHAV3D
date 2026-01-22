import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

"""
Modules for single timepoint feature extraction based on raw features over time windows.
"""

def is_nonneg_int_like(values):
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return False
    return np.all(vals >= 0) and np.allclose(vals, np.round(vals))


def infer_signal_types(df_tracks, columns, treat_01_as_bounded_if_other_values_exist=True):
    """
    Inspect the entire df_tracks once and decide signal type per column.
    Returns signal_type_map.

    Rules (global across ALL rows, per-column):
      - 'binary'     : finite unique ⊆ {0,1}
      - 'bounded'    : all finite within [0,1] (up to tiny tolerance)
      - 'count'      : all finite are non-negative integer-like
      - 'continuous' : otherwise
    """
    signal_type_map = {}
    tol = 1e-8

    for col in columns:
        data = pd.to_numeric(df_tracks[col], errors="coerce").to_numpy(dtype=float)
        finite = data[np.isfinite(data)]
        unique = np.unique(finite) if finite.size else np.array([])

        only_01 = unique.size > 0 and np.all(np.isin(unique, [0.0, 1.0]))
        within_01 = unique.size > 0 and (np.nanmin(unique) >= 0 - tol) and (np.nanmax(unique) <= 1 + tol)

        if only_01:
            signal_type = "binary"
        elif within_01:
            signal_type = "bounded"
        elif is_nonneg_int_like(finite):
            signal_type = "count"
        else:
            signal_type = "continuous"

        signal_type_map[col] = signal_type

    return signal_type_map


def compute_mean_value(signal_values):
    """Return the arithmetic mean of the signal, ignoring NaN values."""
    return float(np.nanmean(signal_values))


def compute_median_value(signal_values):
    """Return the median (50th percentile) of the signal."""
    return float(np.nanmedian(signal_values))


def compute_standard_deviation(signal_values):
    """Return the standard deviation of the signal, measuring spread around the mean."""
    return float(np.nanstd(signal_values))


def compute_minimum_value(signal_values):
    """Return the smallest value observed in the signal."""
    return float(np.nanmin(signal_values))


def compute_maximum_value(signal_values):
    """Return the largest value observed in the signal."""
    return float(np.nanmax(signal_values))


def compute_value_range(signal_values):
    """Return the range of the signal: (max − min)."""
    return compute_maximum_value(signal_values) - compute_minimum_value(signal_values)


def compute_quantile(signal_values, quantile_level):
    """Return a specific quantile (e.g., 0.1 = 10th percentile)."""
    return float(np.nanpercentile(signal_values, quantile_level * 100))


def compute_interquartile_range(signal_values):
    """Return the interquartile range (IQR = Q75 − Q25)."""
    q25, q75 = np.nanpercentile(signal_values, [25, 75])
    return float(q75 - q25)


def compute_median_absolute_deviation(signal_values):
    """Return the median absolute deviation (MAD) of the signal."""
    median_value = np.nanmedian(signal_values)
    absolute_deviation = np.abs(signal_values - median_value)
    return float(np.nanmedian(absolute_deviation))


def compute_linear_trend_slope(signal_values, time_values=None):
    """
    Fit a least-squares line to signal vs. time and return the slope.
    Positive = increasing trend; negative = decreasing.
    """
    signal_values = np.asarray(signal_values, dtype=float)
    if signal_values.size < 2 or np.all(np.isnan(signal_values)):
        return np.nan

    valid_mask = np.isfinite(signal_values)
    y = signal_values[valid_mask]
    if y.size < 2:
        return np.nan

    if time_values is not None:
        t = np.asarray(time_values, dtype=float)[valid_mask]
    else:
        t = np.arange(y.size, dtype=float)

    t_centered = t - np.nanmean(t)
    y_centered = y - np.nanmean(y)
    denom = np.dot(t_centered, t_centered)
    return float(np.dot(t_centered, y_centered) / denom) if denom > 0 else np.nan


def compute_lag1_autocorrelation(signal_values):
    """
    Compute lag-1 autocorrelation (similarity between consecutive values).
    High → smooth signal, low → noisy, negative → oscillating.
    """
    signal_values = np.asarray(signal_values, dtype=float)
    valid_mask = np.isfinite(signal_values)
    x = signal_values[valid_mask]
    if x.size < 2:
        return np.nan
    centered = x - np.mean(x)
    numerator = np.dot(centered[:-1], centered[1:])
    denominator = np.dot(centered, centered)
    return float(numerator / (denominator + 1e-12))


def compute_mean_absolute_first_difference(signal_values):
    """Return the average absolute change between consecutive values."""
    signal_values = np.asarray(signal_values, dtype=float)
    if signal_values.size < 2:
        return np.nan
    return float(np.nanmean(np.abs(np.diff(signal_values))))


def convert_to_binary(signal_values, threshold=0.0):
    """
    Convert a numeric signal to a binary float array {0.0, 1.0} while preserving NaNs.
    """
    x = np.asarray(signal_values, dtype=float)
    out = (x > threshold).astype(float)
    out[np.isnan(x)] = np.nan
    return out


def compute_binary_transition_rate(binary_signal):
    """Return the proportion of time steps where the binary signal switches between 0 and 1."""
    b = np.asarray(binary_signal, dtype=float)
    valid = np.isfinite(b)
    b = b[valid]
    if b.size < 2:
        return np.nan
    return float(np.mean(np.diff(b) != 0))


def compute_longest_true_run_length(binary_signal):
    """Return the longest consecutive run of 1s in the binary signal."""
    b = np.asarray(binary_signal, dtype=float)
    valid = np.isfinite(b)
    b = b[valid].astype(int)

    if b.size == 0:
        return np.nan

    best = cur = 0
    for v in b:
        if v == 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return float(best)


def compute_average_true_run_length(binary_signal):
    """
    Compute the average length of contiguous True runs in a binary signal.
    Returns 0 if there are no True runs.
    """
    b = np.asarray(binary_signal, dtype=float)
    valid = np.isfinite(b)
    b = b[valid].astype(int)

    if b.size == 0:
        return np.nan

    # find runs of 1s
    diff = np.diff(np.concatenate(([0], b, [0])))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]

    if starts.size == 0:
        return 0.0

    run_lengths = ends - starts
    return float(np.mean(run_lengths)) if run_lengths.size else np.nan


def compute_fraction_near_bounds(signal_values, lower_bound=0.0, upper_bound=1.0, margin_fraction=0.05):
    """Return the fractions of values near the lower and upper bounds of a bounded signal."""
    x = np.asarray(signal_values, dtype=float)
    valid = np.isfinite(x)
    if not np.any(valid):
        return np.nan, np.nan

    lower_cutoff = lower_bound + margin_fraction * (upper_bound - lower_bound)
    upper_cutoff = upper_bound - margin_fraction * (upper_bound - lower_bound)

    xv = x[valid]
    return float(np.mean(xv <= lower_cutoff)), float(np.mean(xv >= upper_cutoff))


ALL_WINDOW_FEATURES = {
    # metadata-ish
    "signal_type",

    # Basic motion
    "summed_displacement",
    "net_displacement",
    "straightness",
    "directional_persistence",
    "median_turning_angle",
    "fraction_reversed_movement",
    "mean_square_displacement",

    # basic stats
    "mean",
    "median",
    "range",
    "std",
    "min",
    "max",
    "iqr",
    "mad",
    "quantiles",

    # dynamics
    "trend",
    "lag1_autocorr",
    "mean_abs_first_diff",

    # bounded
    "fraction_near_bounds",

    # binary
    "binary_runs",

    # count
    "count_dispersion",
}


def _compute_motion_features_from_xyz(px, py, pz, features_to_compute=ALL_WINDOW_FEATURES):
    coords = np.column_stack([px, py, pz]).astype(float, copy=False)
    coords_rel = coords - coords[0]  # coords_rel[0] == (0,0,0)

    sq_disp_from_origin = np.sum(coords_rel**2, axis=1)
    mean_square_displacement = float(np.nanmean(sq_disp_from_origin)) if sq_disp_from_origin.size else 0.0

    steps = np.diff(coords_rel, axis=0)
    if steps.size == 0:
        return {
            "summed_displacement": 0.0,
            "net_displacement": 0.0,
            "straightness": 0.0,
            "directional_persistence": 0.0,
            "median_turning_angle": 0.0,
            "fraction_reversed_movement": 0.0,
            "mean_square_displacement": mean_square_displacement,
        }

    step_len = np.linalg.norm(steps, axis=1)
    path_length = float(np.nansum(step_len))
    net_disp = float(np.linalg.norm(coords_rel[-1] - coords_rel[0]))
    straightness = (net_disp / path_length) if path_length > 0 else 0.0

    norms = step_len[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.divide(steps, norms, out=np.zeros_like(steps), where=norms > 0)

    if u.shape[0] >= 2:
        dots = np.sum(u[1:] * u[:-1], axis=1)
        dots = np.clip(dots, -1.0, 1.0)
        turn_angles = np.arccos(dots)
        mean_persist = float(np.nanmean(dots)) if dots.size else 0.0
        median_turn = float(np.nanmedian(turn_angles)) if turn_angles.size else 0.0
        frac_reversed = float(np.mean(dots < 0)) if dots.size else 0.0
    else:
        mean_persist = 0.0
        median_turn = 0.0
        frac_reversed = 0.0

    motion_features = {
        "summed_displacement": path_length,
        "net_displacement": net_disp,
        "straightness": straightness,
        "directional_persistence": mean_persist,
        "median_turning_angle": median_turn,
        "fraction_reversed_movement": frac_reversed,
        "mean_square_displacement": mean_square_displacement,
    }

    return {k: v for k, v in motion_features.items() if k in features_to_compute}


def compute_window_features(
    window_dataframe,
    column_name,
    time_column="position_t",
    signal_types=None,
    quantiles=(0.1, 0.25, 0.75, 0.9),
    features_to_compute=ALL_WINDOW_FEATURES,
):
    """
    features_to_compute:
        - Iterable[str] => compute only those features
    """
    requested = set(features_to_compute)
    unknown = requested - ALL_WINDOW_FEATURES
    if unknown:
        raise ValueError(f"Unknown features_to_compute: {sorted(unknown)}")

    if time_column in window_dataframe.columns:
        time_values = window_dataframe[time_column].to_numpy(float, copy=False)
    else:
        time_values = None

    signal_values = window_dataframe[column_name].to_numpy(float, copy=False)
    signal_type = signal_types.get(column_name, "continuous") if signal_types is not None else "continuous"

    features = {}

    if "signal_type" in requested:
        features[f"{column_name}_inferred_signal_type"] = signal_type

    if "mean" in requested:
        features[f"{column_name}_mean_value"] = compute_mean_value(signal_values)

    # Non-binary statistics
    if signal_type != "binary":
        if "median" in requested:
            features[f"{column_name}_median_value"] = compute_median_value(signal_values)
        if "range" in requested:
            features[f"{column_name}_value_range"] = compute_value_range(signal_values)
        if "std" in requested:
            features[f"{column_name}_standard_deviation"] = compute_standard_deviation(signal_values)
        if "min" in requested:
            features[f"{column_name}_minimum_value"] = compute_minimum_value(signal_values)
        if "max" in requested:
            features[f"{column_name}_maximum_value"] = compute_maximum_value(signal_values)
        if "iqr" in requested:
            features[f"{column_name}_interquartile_range"] = compute_interquartile_range(signal_values)
        if "mad" in requested:
            features[f"{column_name}_median_absolute_deviation"] = compute_median_absolute_deviation(signal_values)

        if "quantiles" in requested:
            for q in quantiles:
                features[f"{column_name}_quantile_{int(q * 100)}percent"] = compute_quantile(signal_values, q)

    # Dynamics
    if "trend" in requested:
        features[f"{column_name}_linear_trend_slope_per_time_unit"] = compute_linear_trend_slope(
            signal_values, time_values
        )

    if "lag1_autocorr" in requested:
        features[f"{column_name}_lag1_autocorrelation"] = compute_lag1_autocorrelation(signal_values)

    if signal_type != "binary" and "mean_abs_first_diff" in requested:
        features[f"{column_name}_mean_absolute_first_difference"] = compute_mean_absolute_first_difference(signal_values)

    # bounded-only
    if signal_type == "bounded" and "fraction_near_bounds" in requested:
        low_frac, high_frac = compute_fraction_near_bounds(signal_values)
        features[f"{column_name}_fraction_near_lower_bound"] = low_frac
        features[f"{column_name}_fraction_near_upper_bound"] = high_frac

    # binary-only
    if signal_type == "binary" and "binary_runs" in requested:
        binary_signal = convert_to_binary(signal_values, threshold=0.0)
        features[f"{column_name}_transition_rate"] = compute_binary_transition_rate(binary_signal)
        features[f"{column_name}_longest_true_length"] = compute_longest_true_run_length(binary_signal)
        features[f"{column_name}_mean_true_length"] = compute_average_true_run_length(binary_signal)

        inv = 1.0 - binary_signal
        inv[np.isnan(binary_signal)] = np.nan
        features[f"{column_name}_longest_false_length"] = compute_longest_true_run_length(inv)
        features[f"{column_name}_mean_false_length"] = compute_average_true_run_length(inv)

    # count-only
    if signal_type == "count" and "count_dispersion" in requested:
        finite = signal_values[np.isfinite(signal_values)]
        if finite.size:
            mean_val = np.mean(finite)
            var_val = np.var(finite)
            features[f"{column_name}_dispersion_index_variance_over_mean"] = (
                var_val / mean_val if mean_val > 0 else np.nan
            )
            features[f"{column_name}_fraction_of_zeros"] = float(np.mean(finite == 0))
        else:
            features[f"{column_name}_dispersion_index_variance_over_mean"] = np.nan
            features[f"{column_name}_fraction_of_zeros"] = np.nan

    return features


def _create_descriptive_track_worker(
    group_df,
    columns_to_summarize,
    window_size,
    step_size,
    time_col,
    id_cols,
    signal_types,
    features_to_compute=ALL_WINDOW_FEATURES,
):
    group_df = group_df.reset_index(drop=True)
    n = len(group_df)
    if n == 0:
        return []

    sample_val = group_df[id_cols[0]].iloc[0]
    track_id_val = group_df[id_cols[1]].iloc[0]

    has_time = time_col in group_df.columns
    t_all = group_df[time_col].to_numpy(float, copy=False) if has_time else None

    px_all = group_df["position_x"].to_numpy(float, copy=False)
    py_all = group_df["position_y"].to_numpy(float, copy=False)
    pz_all = group_df["position_z"].to_numpy(float, copy=False)

    sig_arrays = {col: group_df[col].to_numpy(float, copy=False) for col in columns_to_summarize}

    out_rows = []

    # ---------------------------
    # full-track mode
    # ---------------------------
    if window_size is None:
        start_t = float(t_all[0]) if t_all is not None else 0.0
        end_t = float(t_all[-1]) if t_all is not None else float(n - 1)

        sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

        base = {
            id_cols[0]: sample_val,
            id_cols[1]: track_id_val,
            "sub_TrackID": sub_track_id,
            time_col: np.nan,
            f"window_start_{time_col}": start_t,
            f"window_end_{time_col}": end_t,
            "window_length_frames": n,
        }

        base.update(
            _compute_motion_features_from_xyz(
                px_all,
                py_all,
                pz_all,
                features_to_compute=features_to_compute,
            )
        )

        for col in columns_to_summarize:
            stype = signal_types.get(col, "continuous")
            df_for_feat = (
                pd.DataFrame({col: sig_arrays[col], time_col: t_all})
                if t_all is not None
                else pd.DataFrame({col: sig_arrays[col]})
            )
            base.update(
                compute_window_features(
                    df_for_feat,
                    col,
                    time_column=time_col,
                    signal_types={col: stype},
                    features_to_compute=features_to_compute,
                )
            )

        out_rows.append(base)
        return out_rows

    # -----------------------------------------
    # trailing-per-timepoint sliding window mode
    # windows end at each selected timepoint
    # -----------------------------------------
    stride = step_size if step_size is not None else 1
    stride = max(int(stride), 1)

    nan_feature_template = None

    for end_idx in range(0, n, stride):
        start_idx = end_idx - window_size + 1
        end_t = float(t_all[end_idx]) if t_all is not None else float(end_idx)

        # Not enough history: emit NaNs for all features
        if start_idx < 0:
            row = {
                id_cols[0]: sample_val,
                id_cols[1]: track_id_val,
                "sub_TrackID": f"{int(track_id_val)}_tNaN-t{int(end_t)}",
                time_col: end_t,
                f"window_start_{time_col}": np.nan,
                f"window_end_{time_col}": end_t,
                "window_length_frames": np.nan,
            }
            if nan_feature_template is not None:
                row.update(nan_feature_template)
            out_rows.append(row)
            continue

        # Valid trailing window
        start_t = float(t_all[start_idx]) if t_all is not None else float(start_idx)
        sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

        row = {
            id_cols[0]: sample_val,
            id_cols[1]: track_id_val,
            "sub_TrackID": sub_track_id,
            time_col: end_t,
            f"window_start_{time_col}": start_t,
            f"window_end_{time_col}": end_t,
            "window_length_frames": window_size,
        }

        # motion features for trailing window
        row.update(
            _compute_motion_features_from_xyz(
                px_all[start_idx : end_idx + 1],
                py_all[start_idx : end_idx + 1],
                pz_all[start_idx : end_idx + 1],
                features_to_compute=features_to_compute,
            )
        )

        # per-column signal features over trailing window
        for col in columns_to_summarize:
            stype = signal_types.get(col, "continuous")
            sig_win = sig_arrays[col][start_idx : end_idx + 1]

            if t_all is not None:
                t_win = t_all[start_idx : end_idx + 1]
                wdf = pd.DataFrame({col: sig_win, time_col: t_win})
            else:
                wdf = pd.DataFrame({col: sig_win})

            row.update(
                compute_window_features(
                    wdf,
                    col,
                    time_column=time_col,
                    signal_types={col: stype},
                    features_to_compute=features_to_compute,
                )
            )

        # Build a NaN template based on first valid window
        if nan_feature_template is None:
            nan_feature_template = {
                k: np.nan
                for k in row.keys()
                if k
                not in [
                    id_cols[0],
                    id_cols[1],
                    "sub_TrackID",
                    time_col,
                    f"window_start_{time_col}",
                    f"window_end_{time_col}",
                    "window_length_frames",
                ]
            }

        out_rows.append(row)

    # Backfill NaN feature keys for any early rows created before template existed
    if nan_feature_template is not None:
        for r in out_rows:
            for k in nan_feature_template.keys():
                if k not in r:
                    r[k] = np.nan

    return out_rows


def create_descriptive_track_dataset(
    df_tracks,
    columns_to_summarize,
    window_size=None,
    step_size=None,
    time_col="position_t",
    id_cols=("sample_name", "TrackID"),
    signal_types=None,
    n_jobs=None,
    chunksize=8,
    features_to_compute=ALL_WINDOW_FEATURES,
):
    df_sorted = df_tracks.sort_values(list(id_cols) + [time_col], kind="mergesort")

    if signal_types is None:
        signal_types = infer_signal_types(df_sorted, columns=columns_to_summarize)

    groups = [g for _, g in df_sorted.groupby(list(id_cols), sort=False)]
    total_groups = len(groups)

    output_rows = []

    if n_jobs is None or n_jobs != 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            futures = ex.map(
                _create_descriptive_track_worker,
                groups,
                [columns_to_summarize] * total_groups,
                [window_size] * total_groups,
                [step_size] * total_groups,
                [time_col] * total_groups,
                [list(id_cols)] * total_groups,
                [signal_types] * total_groups,
                [features_to_compute] * total_groups,
                chunksize=chunksize,
            )
            for rows in tqdm(futures, total=total_groups):
                if rows:
                    output_rows.extend(rows)
    else:
        for g in tqdm(groups, total=total_groups):
            rows = _create_descriptive_track_worker(
                g,
                columns_to_summarize,
                window_size,
                step_size,
                time_col,
                list(id_cols),
                signal_types,
                features_to_compute,
            )
            if rows:
                output_rows.extend(rows)

    result = pd.DataFrame(output_rows)

    if not result.empty:
        front_cols = [
            id_cols[0],
            id_cols[1],
            "sub_TrackID",
            time_col,
            f"window_start_{time_col}",
            f"window_end_{time_col}",
            "window_length_frames",
        ]
        other_cols = [c for c in result.columns if c not in front_cols]
        result = result[front_cols + other_cols]

    return result
