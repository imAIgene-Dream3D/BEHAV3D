import numpy as np
import pandas as pd
import os

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
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


_MOTION_FEATURE_KEYS = {
    "summed_displacement",
    "net_displacement",
    "straightness",
    "directional_persistence",
    "median_turning_angle",
    "fraction_reversed_movement",
    "mean_square_displacement",
}


def _needs_min_two_timepoints(feature_key: str) -> bool:
    """Heuristic for features that are undefined/unstable for 1-frame windows."""
    if feature_key in _MOTION_FEATURE_KEYS:
        return True
    min2_suffixes = (
        "_linear_trend_slope_per_time_unit",
        "_lag1_autocorrelation",
        "_mean_absolute_first_difference",
        "_transition_rate",
        "_longest_true_length",
        "_mean_true_length",
        "_longest_false_length",
        "_mean_false_length",
    )
    return feature_key.endswith(min2_suffixes)


def _compute_motion_features_from_xyz(px, py, pz, features_to_compute=ALL_WINDOW_FEATURES):
    coords = np.column_stack([px, py, pz]).astype(float, copy=False)
    # Drop invalid coordinate rows so motion features can still be computed
    # from the valid part of the window.
    valid = np.all(np.isfinite(coords), axis=1)
    coords = coords[valid]
    if coords.shape[0] == 0:
        motion_features = {
            "summed_displacement": 0.0,
            "net_displacement": 0.0,
            "straightness": 0.0,
            "directional_persistence": 0.0,
            "median_turning_angle": 0.0,
            "fraction_reversed_movement": 0.0,
            "mean_square_displacement": 0.0,
        }
        return {k: v for k, v in motion_features.items() if k in features_to_compute}
    coords_rel = coords - coords[0]  # coords_rel[0] == (0,0,0)

    sq_disp_from_origin = np.sum(coords_rel**2, axis=1)
    mean_square_displacement = float(np.nanmean(sq_disp_from_origin)) if sq_disp_from_origin.size else 0.0

    steps = np.diff(coords_rel, axis=0)
    if steps.size == 0:
        motion_features = {
            "summed_displacement": 0.0,
            "net_displacement": 0.0,
            "straightness": 0.0,
            "directional_persistence": 0.0,
            "median_turning_angle": 0.0,
            "fraction_reversed_movement": 0.0,
            "mean_square_displacement": mean_square_displacement,
        }
        return {k: v for k, v in motion_features.items() if k in features_to_compute}

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
    if signal_type == "binary":
        binary_signal = convert_to_binary(signal_values, threshold=0.0)
        if "transition_rate" in requested:
            features[f"{column_name}_transition_rate"] = compute_binary_transition_rate(binary_signal)
        
        if "binary_runs" in requested:
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
    only_nonbinary=False,
    incomplete_window_policy="drop",
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
    
    # Separate binary and non-binary columns if only_nonbinary is True
    if only_nonbinary:
        binary_cols = [col for col in columns_to_summarize if signal_types.get(col, "continuous") == "binary"]
        nonbinary_cols = [col for col in columns_to_summarize if signal_types.get(col, "continuous") != "binary"]
    else:
        binary_cols = []
        nonbinary_cols = columns_to_summarize

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

        # Compute window features for non-binary columns
        for col in nonbinary_cols:
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
        
        # For binary columns, just add the last value (or mean for full track)
        for col in binary_cols:
            # For full track mode, we'll use the mean as a summary
            base[f"{col}_value"] = float(np.nanmean(sig_arrays[col]))

        out_rows.append(base)
        return out_rows

    # -----------------------------------------
    # trailing-per-timepoint sliding window mode
    # windows end at each selected timepoint
    # -----------------------------------------
    stride = step_size if step_size is not None else 1
    stride = max(int(stride), 1)

    policy = str(incomplete_window_policy).lower()
    if policy not in {"drop", "partial", "as_far_as_possible"}:
        raise ValueError(
            "incomplete_window_policy must be one of: "
            "'drop', 'partial' (alias: 'as_far_as_possible')."
        )

    idx_all = np.arange(n, dtype=int)
    time_is_integer_like = False
    track_min_t = None
    if t_all is not None:
        finite_t = t_all[np.isfinite(t_all)]
        if finite_t.size > 0:
            time_is_integer_like = bool(np.allclose(finite_t, np.round(finite_t)))
            track_min_t = float(np.nanmin(finite_t))

    for end_idx in range(0, n, stride):
        end_t = float(t_all[end_idx]) if t_all is not None else float(end_idx)
        if t_all is None:
            # Fallback for inputs without time column: keep row-based trailing windows.
            start_idx = end_idx - window_size + 1
            if start_idx < 0 and policy == "drop":
                row = {
                    id_cols[0]: sample_val,
                    id_cols[1]: track_id_val,
                    "sub_TrackID": f"{int(track_id_val)}_tNaN-t{int(end_t)}",
                    time_col: end_t,
                    f"window_start_{time_col}": np.nan,
                    f"window_end_{time_col}": end_t,
                    "window_length_frames": np.nan,
                }
                out_rows.append(row)
                continue
            if start_idx < 0 and policy in {"partial", "as_far_as_possible"}:
                start_idx = 0
            win_idx = idx_all[start_idx : end_idx + 1]
        else:
            # Timepoint-based trailing window: [t-window_size+1, t], clipped in partial mode.
            target_start_t = float(end_t - int(window_size) + 1)
            effective_start_t = (
                max(float(track_min_t), float(target_start_t))
                if track_min_t is not None
                else float(target_start_t)
            )
            win_mask = (
                (idx_all <= int(end_idx))
                & np.isfinite(t_all)
                & (t_all >= effective_start_t)
                & (t_all <= float(end_t))
            )
            win_idx = np.where(win_mask)[0]
            if win_idx.size == 0:
                win_idx = np.asarray([int(end_idx)], dtype=int)

            if policy == "drop":
                has_full_history = False
                if track_min_t is not None and np.isfinite(target_start_t) and float(target_start_t) >= float(track_min_t):
                    if time_is_integer_like:
                        expected = np.arange(
                            int(np.round(target_start_t)),
                            int(np.round(end_t)) + 1,
                            dtype=int,
                        )
                        observed = np.unique(np.round(t_all[win_idx]).astype(int))
                        has_full_history = bool(np.all(np.isin(expected, observed)))
                    else:
                        has_full_history = bool(len(win_idx) >= int(window_size))

                if not has_full_history:
                    row = {
                        id_cols[0]: sample_val,
                        id_cols[1]: track_id_val,
                        "sub_TrackID": f"{int(track_id_val)}_tNaN-t{int(end_t)}",
                        time_col: end_t,
                        f"window_start_{time_col}": np.nan,
                        f"window_end_{time_col}": end_t,
                        "window_length_frames": np.nan,
                    }
                    out_rows.append(row)
                    continue

        start_t = float(t_all[win_idx[0]]) if t_all is not None else float(win_idx[0])
        sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

        row = {
            id_cols[0]: sample_val,
            id_cols[1]: track_id_val,
            "sub_TrackID": sub_track_id,
            time_col: end_t,
            f"window_start_{time_col}": start_t,
            f"window_end_{time_col}": end_t,
            "window_length_frames": int(len(win_idx)),
        }

        # motion features for trailing window
        row.update(
            _compute_motion_features_from_xyz(
                px_all[win_idx],
                py_all[win_idx],
                pz_all[win_idx],
                features_to_compute=features_to_compute,
            )
        )

        # per-column signal features over trailing window (non-binary only)
        for col in nonbinary_cols:
            stype = signal_types.get(col, "continuous")
            sig_win = sig_arrays[col][win_idx]

            if t_all is not None:
                t_win = t_all[win_idx]
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
        
        # For binary columns, just add the value at the current timepoint (end_idx)
        for col in binary_cols:
            row[f"{col}_value"] = float(sig_arrays[col][end_idx])

        out_rows.append(row)

    # In partial mode, first row can be single-frame. For features that need at least
    # 2 timepoints, copy values from t+1 row if available, otherwise nearest later
    # row with window_length_frames >= 2.
    if policy in {"partial", "as_far_as_possible"} and len(out_rows) >= 2:
        first = out_rows[0]
        first_window_len = pd.to_numeric([first.get("window_length_frames", np.nan)], errors="coerce")[0]
        if np.isfinite(first_window_len) and int(first_window_len) <= 1:
            source_row = None
            first_t = pd.to_numeric([first.get(time_col, np.nan)], errors="coerce")[0]
            if np.isfinite(first_t):
                target_t = float(first_t) + 1.0
                for candidate in out_rows[1:]:
                    cand_t = pd.to_numeric([candidate.get(time_col, np.nan)], errors="coerce")[0]
                    cand_len = pd.to_numeric([candidate.get("window_length_frames", np.nan)], errors="coerce")[0]
                    if np.isfinite(cand_t) and np.isfinite(cand_len):
                        if int(cand_len) >= 2 and np.isclose(float(cand_t), float(target_t)):
                            source_row = candidate
                            break
            if source_row is None:
                for candidate in out_rows[1:]:
                    cand_len = pd.to_numeric([candidate.get("window_length_frames", np.nan)], errors="coerce")[0]
                    if np.isfinite(cand_len) and int(cand_len) >= 2:
                        source_row = candidate
                        break
            if source_row is None:
                return out_rows

            protected = {
                id_cols[0],
                id_cols[1],
                "sub_TrackID",
                time_col,
                f"window_start_{time_col}",
                f"window_end_{time_col}",
                "window_length_frames",
            }
            for k in list(first.keys()):
                if k in protected:
                    continue
                if _needs_min_two_timepoints(k) and (k in source_row):
                    first[k] = source_row[k]

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
    only_nonbinary=False,
    incomplete_window_policy="drop",
):
    df_sorted = df_tracks.sort_values(list(id_cols) + [time_col], kind="mergesort")

    if signal_types is None:
        signal_types = infer_signal_types(df_sorted, columns=columns_to_summarize)

    groups = [g for _, g in df_sorted.groupby(list(id_cols), sort=False)]
    total_groups = len(groups)

    output_rows = []

    # n_jobs=None means "auto": use all available CPU cores.
    # Keep a robust fallback to serial if the process pool breaks.
    resolved_n_jobs = os.cpu_count() if n_jobs is None else n_jobs
    if resolved_n_jobs is None:
        resolved_n_jobs = 1
    resolved_n_jobs = int(resolved_n_jobs)
    use_parallel = resolved_n_jobs != 1

    if use_parallel:
        try:
            with ProcessPoolExecutor(max_workers=resolved_n_jobs) as ex:
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
                    [only_nonbinary] * total_groups,
                    [incomplete_window_policy] * total_groups,
                    chunksize=chunksize,
                )
                for rows in tqdm(futures, total=total_groups):
                    if rows:
                        output_rows.extend(rows)
        except BrokenProcessPool:
            # Recover from abrupt worker exits by retrying deterministically in-process.
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
                    only_nonbinary,
                    incomplete_window_policy,
                )
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
                only_nonbinary,
                incomplete_window_policy,
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
