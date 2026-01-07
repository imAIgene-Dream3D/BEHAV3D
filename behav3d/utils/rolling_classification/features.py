import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Literal, List, Iterable
from pandas.api.types import is_numeric_dtype

from tqdm import tqdm

import numpy as np

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from sklearn.feature_selection import VarianceThreshold


from behav3d.utils.rolling_classification import *

SignalType = Literal["binary", "count", "bounded", "continuous"]

"""
Modules for single timepoint feature extraction based on raw features over time windows.
"""
def is_nonneg_int_like(values: np.ndarray) -> bool:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return False
    return np.all(vals >= 0) and np.allclose(vals, np.round(vals))

def infer_signal_types(
    df_tracks: pd.DataFrame,
    columns: list,
    *,
    treat_01_as_bounded_if_other_values_exist: bool = True,
) -> Tuple[Dict[str, SignalType], Dict[str, Tuple[float, float]]]:
    """
    Inspect the entire df_tracks once and decide signal type per column.
    Returns (signal_type_map, bounds_map).
    Rules (global across ALL rows):
      - 'binary'  : finite unique ⊆ {0,1}. If treat_01_as_bounded_if_other_values_exist=True
                    and there exist any values strictly between 0 and 1 globally, prefer 'bounded'.
      - 'bounded' : all finite within [0,1] (up to tiny tolerance) or explicitly chosen as above.
      - 'count'   : all finite are non-negative integer-like.
      - 'continuous': otherwise.
    Also returns bounds_map: for bounded columns (or ones within [0,1]) -> (0,1),
    else (global_min, global_max) for reference.
    """
    signal_type_map: Dict[str, SignalType] = {}
    bounds_map: Dict[str, Tuple[float, float]] = {}
    tol = 1e-8

    for col in columns:
        data = pd.to_numeric(df_tracks[col], errors="coerce").to_numpy(float)
        finite = data[np.isfinite(data)]
        global_min = float(np.nanmin(finite)) if finite.size else np.nan
        global_max = float(np.nanmax(finite)) if finite.size else np.nan
        unique = np.unique(finite) if finite.size else np.array([])

        # detect global properties
        only_01 = unique.size > 0 and np.all(np.isin(unique, [0.0, 1.0]))
        any_between_01 = unique.size > 0 and (np.nanmin(unique) < 0 - tol or np.nanmax(unique) > 1 + tol or
                                              np.any((unique > 0 + tol) & (unique < 1 - tol)))
        within_01 = unique.size > 0 and (np.nanmin(unique) >= 0 - tol) and (np.nanmax(unique) <= 1 + tol)

        if only_01 and treat_01_as_bounded_if_other_values_exist and any_between_01:
            # e.g., some windows 0/1, but globally we also see values in (0,1) → call it bounded globally
            signal_type = "bounded"
            bounds = (0.0, 1.0)
        elif only_01:
            signal_type = "binary"
            bounds = (0.0, 1.0)
        elif within_01:
            signal_type = "bounded"
            bounds = (0.0, 1.0)
        elif is_nonneg_int_like(finite):
            signal_type = "count"
            bounds = (global_min, global_max)
        else:
            signal_type = "continuous"
            bounds = (global_min, global_max)

        signal_type_map[col] = signal_type
        bounds_map[col] = bounds

    return signal_type_map#, bounds_map

def compute_mean_value(signal_values: np.ndarray) -> float:
    """Return the arithmetic mean of the signal, ignoring NaN values."""
    return float(np.nanmean(signal_values))

def compute_median_value(signal_values: np.ndarray) -> float:
    """Return the median (50th percentile) of the signal."""
    return float(np.nanmedian(signal_values))

def compute_standard_deviation(signal_values: np.ndarray) -> float:
    """Return the standard deviation of the signal, measuring spread around the mean."""
    return float(np.nanstd(signal_values))

def compute_minimum_value(signal_values: np.ndarray) -> float:
    """Return the smallest value observed in the signal."""
    return float(np.nanmin(signal_values))

def compute_maximum_value(signal_values: np.ndarray) -> float:
    """Return the largest value observed in the signal."""
    return float(np.nanmax(signal_values))

def compute_value_range(signal_values: np.ndarray) -> float:
    """Return the range of the signal: (max − min)."""
    return compute_maximum_value(signal_values) - compute_minimum_value(signal_values)

def compute_quantile(signal_values: np.ndarray, quantile_level: float) -> float:
    """Return a specific quantile (e.g., 0.1 = 10th percentile)."""
    return float(np.nanpercentile(signal_values, quantile_level * 100))

def compute_interquartile_range(signal_values: np.ndarray) -> float:
    """Return the interquartile range (IQR = Q75 − Q25)."""
    q25, q75 = np.nanpercentile(signal_values, [25, 75])
    return float(q75 - q25)

def compute_median_absolute_deviation(signal_values: np.ndarray) -> float:
    """Return the median absolute deviation (MAD) of the signal."""
    median_value = np.nanmedian(signal_values)
    absolute_deviation = np.abs(signal_values - median_value)
    return float(np.nanmedian(absolute_deviation))

def compute_linear_trend_slope(signal_values: np.ndarray, time_values: Optional[np.ndarray]) -> float:
    """
    Fit a least-squares line to signal vs. time and return the slope.
    Positive = increasing trend; negative = decreasing.
    """
    signal_values = np.asarray(signal_values, dtype=float)
    if signal_values.size < 2 or np.all(np.isnan(signal_values)):
        return np.nan

    valid_mask = ~np.isnan(signal_values)
    signal_values = signal_values[valid_mask]
    if time_values is not None:
        time_vector = np.asarray(time_values, dtype=float)[valid_mask]
    else:
        time_vector = np.arange(signal_values.size, dtype=float)

    time_vector_centered = time_vector - np.nanmean(time_vector)
    numerator = np.dot(time_vector_centered, signal_values - np.nanmean(signal_values))
    denominator = np.dot(time_vector_centered, time_vector_centered)
    return float(numerator / denominator) if denominator > 0 else np.nan

def compute_lag1_autocorrelation(signal_values: np.ndarray) -> float:
    """
    Compute lag-1 autocorrelation (similarity between consecutive values).
    High → smooth signal, low → noisy, negative → oscillating.
    """
    signal_values = np.asarray(signal_values, dtype=float)
    valid_mask = ~np.isnan(signal_values)
    signal_values = signal_values[valid_mask]
    if signal_values.size < 2:
        return np.nan
    centered = signal_values - np.mean(signal_values)
    numerator = np.dot(centered[:-1], centered[1:])
    denominator = np.dot(centered, centered)
    return float(numerator / (denominator + 1e-12))

def compute_mean_absolute_first_difference(signal_values: np.ndarray) -> float:
    """Return the average absolute change between consecutive values."""
    if len(signal_values) < 2:
        return np.nan
    return float(np.nanmean(np.abs(np.diff(signal_values))))

def check_if_nonnegative_integer_like(signal_values: np.ndarray) -> bool:
    """Return True if all finite values are non-negative and integer-like."""
    finite_values = np.asarray(signal_values)[np.isfinite(signal_values)]
    return finite_values.size > 0 and np.all(finite_values >= 0) and np.allclose(finite_values, np.round(finite_values))

def convert_to_binary(signal_values: np.ndarray) -> np.ndarray:
    """Convert a numeric signal to binary (0/1) using a threshold."""
    return np.asarray(signal_values, dtype=bool)

def compute_binary_transition_rate(binary_signal: np.ndarray) -> float:
    """Return the proportion of time steps where the binary signal switches between 0 and 1."""
    valid = ~np.isnan(binary_signal)
    binary_signal = binary_signal[valid]
    return np.nan if binary_signal.size < 2 else float(np.mean(np.diff(binary_signal) != 0))

def compute_longest_true_run_length(binary_signal: np.ndarray) -> float:
    """Return the longest consecutive run of 1s in the binary signal."""
    valid = ~np.isnan(binary_signal)
    binary_signal = binary_signal[valid].astype(int)
    best = cur = 0
    for v in binary_signal:
        if v == 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return float(best) if binary_signal.size else np.nan

def compute_average_true_run_length(binary_signal: np.ndarray) -> float:
    """
    Compute the average length of contiguous True runs in a binary signal.
    Returns NaN if there are no True runs.
    """
    if binary_signal.size == 0:
        return np.nan

    binary_signal = np.asarray(binary_signal, dtype=bool)
    diff = np.diff(np.concatenate(([0], binary_signal.view(np.int8), [0])))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]

    if starts.size == 0:
        return 0

    run_lengths = ends - starts
    return float(np.mean(run_lengths)) if run_lengths.size else np.nan

def compute_fraction_near_bounds(signal_values: np.ndarray, lower_bound=0.0, upper_bound=1.0, margin_fraction=0.05) -> Tuple[float, float]:
    """Return the fractions of values near the lower and upper bounds of a bounded signal."""
    valid = np.isfinite(signal_values)
    if not np.any(valid):
        return np.nan, np.nan
    lower_cutoff = lower_bound + margin_fraction * (upper_bound - lower_bound)
    upper_cutoff = upper_bound - margin_fraction * (upper_bound - lower_bound)
    return (
        float(np.mean(signal_values[valid] <= lower_cutoff)),
        float(np.mean(signal_values[valid] >= upper_cutoff))
    )

ALL_WINDOW_FEATURES = {
    # metadata-ish
    "signal_type",
    
    #Basic motion 
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
    "quantiles",  # includes the q list

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
    coords_rel = coords - coords[0]  # now coords_rel[0] == (0,0,0)

    sq_disp_from_origin = np.sum(coords_rel**2, axis=1)  # ||r(t)-r(0)||^2
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

    motion_features={
        "summed_displacement": path_length,
        "net_displacement": net_disp,
        "straightness": straightness,
        "directional_persistence": mean_persist,
        "median_turning_angle": median_turn,
        "fraction_reversed_movement": frac_reversed,
        "mean_square_displacement": mean_square_displacement,
    }
    
    motion_features = {k: v for k, v in motion_features.items() if k in features_to_compute}
    return motion_features

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
        - None => compute everything (backwards compatible)
        - Iterable[str] => compute only those features (e.g. {"mean", "trend"})
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
                features[f"{column_name}_quantile_{int(q*100)}percent"] = compute_quantile(signal_values, q)

    # Dynamics (can apply to binary too, but your old code computed trend + autocorr always)
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
        binary_signal = convert_to_binary(signal_values)
        features[f"{column_name}_transition_rate"] = compute_binary_transition_rate(binary_signal)
        features[f"{column_name}_longest_true_length"] = compute_longest_true_run_length(binary_signal)
        features[f"{column_name}_mean_true_length"] = compute_average_true_run_length(binary_signal)
        features[f"{column_name}_longest_false_length"] = compute_longest_true_run_length(~binary_signal)
        features[f"{column_name}_mean_false_length"] = compute_average_true_run_length(~binary_signal)

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
    features_to_compute=ALL_WINDOW_FEATURES
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
    # full-track mode (unchanged)
    # ---------------------------
    if window_size is None:
        start_t = float(t_all[0]) if t_all is not None else 0.0
        end_t = float(t_all[-1]) if t_all is not None else float(n - 1)

        sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

        base = {
            id_cols[0]: sample_val,
            id_cols[1]: track_id_val,
            "sub_TrackID": sub_track_id,
            time_col: np.nan,  # <-- no single center timepoint in full-track mode
            f"window_start_{time_col}": start_t,
            f"window_end_{time_col}": end_t,
            "window_length_frames": n,
        }

        base.update(_compute_motion_features_from_xyz(
            px_all, 
            py_all, 
            pz_all, 
            features_to_compute=features_to_compute)
            )

        for col in columns_to_summarize:
            stype = signal_types.get(col, "continuous")
            df_for_feat = (
                pd.DataFrame({col: sig_arrays[col], time_col: t_all})
                if t_all is not None else
                pd.DataFrame({col: sig_arrays[col]})
            )
            base.update(
                compute_window_features(
                    df_for_feat,
                    col,
                    time_column=time_col,
                    signal_types={col: stype},
                    features_to_compute=features_to_compute
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

        # current/end timepoint
        end_t = float(t_all[end_idx]) if t_all is not None else float(end_idx)

        # Not enough history: emit NaNs for all features
        if start_idx < 0:
            row = {
                id_cols[0]: sample_val,
                id_cols[1]: track_id_val,
                "sub_TrackID": f"{int(track_id_val)}_tNaN-t{int(end_t)}",
                time_col: end_t,  # <-- keep center timepoint for merging
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
            time_col: end_t,  # <-- center/endpoint timepoint for merging
            f"window_start_{time_col}": start_t,
            f"window_end_{time_col}": end_t,
            "window_length_frames": window_size,
        }

        # motion features for trailing window
        row.update(
            _compute_motion_features_from_xyz(
                px_all[start_idx:end_idx + 1],
                py_all[start_idx:end_idx + 1],
                pz_all[start_idx:end_idx + 1],
                features_to_compute=features_to_compute
            )
        )

        # per-column signal features over trailing window
        for col in columns_to_summarize:
            stype = signal_types.get(col, "continuous")
            sig_win = sig_arrays[col][start_idx:end_idx + 1]

            if t_all is not None:
                t_win = t_all[start_idx:end_idx + 1]
                wdf = pd.DataFrame({col: sig_win, time_col: t_win})
            else:
                wdf = pd.DataFrame({col: sig_win})

            row.update(
                compute_window_features(
                    wdf,
                    col,
                    time_column=time_col,
                    signal_types={col: stype},
                    features_to_compute=features_to_compute
                )
            )

        # Build a NaN template based on first valid window
        if nan_feature_template is None:
            nan_feature_template = {
                k: np.nan for k in row.keys()
                if k not in [
                    id_cols[0], id_cols[1], "sub_TrackID",
                    time_col,
                    f"window_start_{time_col}", f"window_end_{time_col}",
                    "window_length_frames"
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
    id_cols=["sample_name", "TrackID"],
    signal_types=None,
    n_jobs=None,
    chunksize=8,
    features_to_compute=ALL_WINDOW_FEATURES
):
    df_sorted = df_tracks.sort_values(id_cols + [time_col], kind="mergesort")

    if signal_types is None:
        signal_types = infer_signal_types(df_sorted, columns=columns_to_summarize)

    groups = [g for _, g in df_sorted.groupby(id_cols, sort=False)]
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
                [id_cols] * total_groups,
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
                g, columns_to_summarize, window_size, step_size, time_col, id_cols, signal_types, features_to_compute
            )
            if rows:
                output_rows.extend(rows)

    result = pd.DataFrame(output_rows)

    if not result.empty:
        front_cols = [
            id_cols[0], id_cols[1], "sub_TrackID",
            time_col,  # <-- endpoint timepoint for merging
            f"window_start_{time_col}", f"window_end_{time_col}",
            "window_length_frames"
        ]
        other_cols = [c for c in result.columns if c not in front_cols]
        result = result[front_cols + other_cols]

    return result

"""
Full track features based on states
"""
def rle_encode(states):
    """Return list of (state, run_length)."""
    if len(states) == 0:
        return []
    runs = []
    cur = states[0]
    length = 1
    for s in states[1:]:
        if s == cur:
            length += 1
        else:
            runs.append((cur, length))
            cur = s
            length = 1
    runs.append((cur, length))
    return runs

def fractions_from_runs(runs, states):
    total = sum(l for _, l in runs) if runs else 0
    cols = [f"overall_fraction_{s}" for s in states]
    out = {c: 0.0 for c in cols}
    if total == 0:
        return out, cols
    for s, l in runs:
        out[f"overall_fraction_{s}"] += l / total
    return out, cols

def bout_stats_from_runs(runs, states):
    lengths = {s: [] for s in states}
    for s, l in runs:
        lengths[s].append(l)

    total_bouts = len(runs)
    track_length = sum(l for _, l in runs)

    cols = []
    for s in states:
        cols.extend([f"bouts_nr_{s}", f"bouts_mean_length_{s}", f"bouts_max_length_{s}"])

    if total_bouts == 0 or track_length == 0:
        return {c: 0.0 for c in cols}, cols

    out = {}
    for s in states:
        arr = np.asarray(lengths[s], dtype=float)
        out[f"bouts_nr_{s}"] = float(arr.size) / float(total_bouts)
        out[f"bouts_mean_length_{s}"] = float(arr.mean()) / float(track_length) if arr.size else 0.0
        out[f"bouts_max_length_{s}"] = float(arr.max()) / float(track_length) if arr.size else 0.0

    return out, cols

def transition_probs_from_runs(runs, state_to_idx, states):
    K = len(states)
    M = np.zeros((K, K), dtype=float)
    labels = [s for s, _ in runs]
    if len(labels) < 2:
        return M  # all zeros
    for a, b in zip(labels[:-1], labels[1:]):
        if a in state_to_idx and b in state_to_idx:
            M[state_to_idx[a], state_to_idx[b]] += 1.0
    row_sums = M.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return M / row_sums

def ngram_counts_from_runs(runs, n=3, weight="count"):
    labels = [s for s, _ in runs]
    lens = [l for _, l in runs]
    out = {}
    if len(labels) < n:
        return out
    for i in range(len(labels) - n + 1):
        g = tuple(labels[i : i + n])
        w = float(lens[i]) if weight == "duration" else 1.0
        out[g] = out.get(g, 0.0) + w
    return out

def extract_descibing_track_state_features(
    adata,
    group_col=("sample_name", "TrackID"),
    time_col="position_t",
    state_col="ClusterID",
    use_fractions=True,
    use_bout_stats=True,
    use_transitions=True,
    use_bigrams=True,
    use_trigrams=True,
    ngram_weight="count",
):
    """
    Track-level feature extraction -> returns (track_adata, blocks)

    track_adata.X        = features (n_tracks x n_features)
    track_adata.obs      = per-track metadata (sample_name, TrackID, position_t_min/max, n_timepoints)
    track_adata.var_names= feature names
    """

    group_col = list(group_col) if isinstance(group_col, (list, tuple)) else [group_col]

    obs = adata.obs[group_col + [time_col, state_col]].copy()
    obs = obs.sort_values(group_col + [time_col])

    # Stable state universe
    states = pd.Index(obs[state_col].astype("category").cat.categories).tolist()
    state_to_idx = {s: i for i, s in enumerate(states)}

    # Use observed=True to avoid unobserved categorical groups exploding results
    gb = obs.groupby(group_col, sort=False, observed=True)

    # Collect runs + n-gram vocab
    runs_by_track = []
    bigram_set, trigram_set = set(), set()

    for tid, df_t in gb:
        seq = df_t[state_col].tolist()
        runs = rle_encode(seq)  # <-- your RLE function

        runs_by_track.append((tid, runs))

        if use_bigrams:
            bigram_set |= set(ngram_counts_from_runs(runs, n=2, weight=ngram_weight).keys())
        if use_trigrams:
            trigram_set |= set(ngram_counts_from_runs(runs, n=3, weight=ngram_weight).keys())

    bigram_list = sorted(bigram_set)
    trigram_list = sorted(trigram_set)

    def bg_col(g): return f"bigram_{g[0]}>{g[1]}"
    def tg_col(g): return f"trigram_{g[0]}>{g[1]}>{g[2]}"

    # Track block columns explicitly while building features
    fraction_cols = []
    bout_cols = []
    transition_cols = []
    bigram_cols = [bg_col(g) for g in bigram_list] if use_bigrams else []
    trigram_cols = [tg_col(g) for g in trigram_list] if use_trigrams else []

    rows = []
    id_rows = []

    for tid, runs in runs_by_track:
        # tid is scalar if len(group_col)==1, else tuple
        if len(group_col) == 1:
            id_rows.append([tid])
        else:
            id_rows.append(list(tid))

        feats = {}

        if use_fractions:
            f, fcols = fractions_from_runs(runs, states)
            feats.update(f)
            if not fraction_cols:
                fraction_cols = fcols

        if use_bout_stats:
            b, bcols = bout_stats_from_runs(runs, states)
            feats.update(b)
            if not bout_cols:
                bout_cols = bcols

        if use_transitions:
            T = transition_probs_from_runs(runs, state_to_idx, states)
            if not transition_cols:
                transition_cols = [f"transitions_{a}>{b}" for a in states for b in states]
            for a in states:
                ia = state_to_idx[a]
                for b in states:
                    ib = state_to_idx[b]
                    feats[f"transitions_{a}>{b}"] = float(T[ia, ib])

        if use_bigrams:
            bg = ngram_counts_from_runs(runs, n=2, weight=ngram_weight)
            total = sum(bg.values()) or 1.0
            for g in bigram_list:
                feats[bg_col(g)] = bg.get(g, 0.0) / total

        if use_trigrams:
            tg = ngram_counts_from_runs(runs, n=3, weight=ngram_weight)
            total = sum(tg.values()) or 1.0
            for g in trigram_list:
                feats[tg_col(g)] = tg.get(g, 0.0) / total

        rows.append(feats)

    # Features
    df_feat = pd.DataFrame(rows).fillna(0.0)

    # IDs in the exact same order as df_feat rows
    df_ids = pd.DataFrame(id_rows, columns=group_col)

    # Time summaries (also observed=True to avoid unobserved category combos)
    df_time = (
        obs.groupby(group_col, sort=False, observed=True)[time_col]
           .agg(position_t_min="min", position_t_max="max", n_timepoints="size")
           .reset_index()
    )

    # Build obs table in the same row order as features
    df_meta = df_ids.merge(df_time, on=group_col, how="left")

    # Combine for conversion using your helper
    df_out = pd.concat([df_meta.reset_index(drop=True), df_feat.reset_index(drop=True)], axis=1)

    blocks = [fraction_cols, bout_cols, transition_cols, bigram_cols, trigram_cols]

    feature_cols = df_feat.columns.tolist()
    obs_cols = df_meta.columns.tolist()

    track_adata = df_to_adata(df_out, feature_cols=feature_cols, obs_cols=obs_cols)

    # Optional bookkeeping
    track_adata.uns["feature_blocks"] = {
        name: cols for name, cols in zip(
            ["fractions", "bout_stats", "transitions", "bigrams", "trigrams"], blocks
        ) if cols
    }
    track_adata.uns["feature_params"] = {
        "group_col": group_col,
        "time_col": time_col,
        "state_col": state_col,
        "ngram_weight": ngram_weight,
        "use_fractions": use_fractions,
        "use_bout_stats": use_bout_stats,
        "use_transitions": use_transitions,
        "use_bigrams": use_bigrams,
        "use_trigrams": use_trigrams,
    }

    return track_adata, blocks


"""
Accessory feautre selection or scaling methods
"""

def scale_feature_blocks(
    adata,
    blocks,
    layer=None,
    mode="sqrt"  # "sqrt" or "linear"
):
    for block in blocks:
        d = len(block)
        if d == 0:
            continue

        scale = np.sqrt(d) if mode == "sqrt" else float(d)

        idx = adata.var_names.get_indexer(block)

        if np.any(idx < 0):
            raise ValueError("Some block features not found in adata.var_names")

        if layer is not None:
            adata.layers[layer][:, idx] /= scale
        else:
            adata.X[:, idx] /= scale

    return adata

def l2_normalize_block(adata, cols, layer=None):
    """Row-wise L2 normalize selected columns. Leaves all-zero rows unchanged."""
    adata_norm = adata.copy()
    if layer is not None:
        X = adata_norm[:, cols].layers[layer]
    else:
        X = adata_norm[:, cols].X
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata_norm[:, cols] = X / norms
    return adata_norm

def l2_normalize_all(adata):
    """Row-wise L2 normalize all columns. Leaves all-zero rows unchanged."""
   
    X = adata.X.copy()
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata.X = X / norms
    return(adata)

def l2_normalize_features_blocks(
    adata,
    blocks,
    layer=None
    ):
    for block in blocks:
        adata = l2_normalize_block(adata, block, layer=layer)
    return adata
 
def drop_highly_correlated_features(
    data,
    feature_cols: list,
    threshold: float = 0.98,
    prefer_keep: list | None = None,
    also_drop_constant: bool = True,
    redundancy_tolerance_abs: float = 0.02,       # apply prefer_keep only if |Δredundancy| <= this
    redundancy_tolerance_rel: float = 0.05,
) -> tuple[object, list, pd.DataFrame]:
    """
    Remove one feature from any pair with |corr| > threshold.

    Accepts either:
      - pandas DataFrame (features are columns)
      - AnnData (features are var_names; uses .X via .to_df())

    Returns:
      reduced_data : same type as input (DataFrame or AnnData), with dropped features removed
      dropped      : list of dropped feature names
      report       : DataFrame logging each decision: kept, dropped, correlation, reason
    """
    is_adata = hasattr(data, "var_names") and hasattr(data, "obs_names")

    # Build working feature frame X (pandas) and restrict to existing features
    if is_adata:
        existing = [f for f in feature_cols if f in data.var_names]
        X = data[:, existing].to_df() if existing else pd.DataFrame(index=data.obs_names)
        feature_cols_eff = existing
    else:
        existing = [f for f in feature_cols if f in data.columns]
        X = data[existing].copy() if existing else pd.DataFrame(index=data.index)
        feature_cols_eff = existing

    # Your current behavior: ignore passed prefer_keep and auto-prefer median/mean features
    prefer_keep = [f for f in feature_cols_eff if "median" in f]
    prefer_keep += [f for f in feature_cols_eff if "mean" in f]
    prefer_keep = set(prefer_keep)

    # Numeric-only
    cols_ok = [c for c in X.columns if is_numeric_dtype(X[c])]
    X = X[cols_ok]

    # Replace infs and optionally drop constants
    X = X.replace([np.inf, -np.inf], np.nan)
    dropped: list[str] = []
    decisions: list[dict] = []

    if also_drop_constant and X.shape[1] > 0:
        nunique = X.nunique(dropna=True)
        const_cols = nunique[nunique <= 1].index.tolist()
        for c in const_cols:
            decisions.append({
                "kept_feature": None,
                "dropped_feature": c,
                "abs_correlation": np.nan,
                "reason": "constant_or_single_value"
            })
        if const_cols:
            X = X.drop(columns=const_cols)
            dropped.extend(const_cols)

    if X.shape[1] <= 1:
        report = pd.DataFrame(decisions, columns=["kept_feature", "dropped_feature", "abs_correlation", "reason"])
        if is_adata:
            keep_vars = [v for v in data.var_names if v not in set(dropped)]
            return data[:, keep_vars].copy(), dropped, report
        else:
            reduced_df = data.drop(columns=[c for c in dropped if c in data.columns])
            return reduced_df, dropped, report

    # Impute NaNs so corr works
    X_impute = X.fillna(X.mean(numeric_only=True))

    # Abs correlation and upper triangle
    corr = X_impute.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

    to_drop = set()

    for col in upper.columns:
        partners = upper.index[upper[col] > threshold].tolist()
        for partner in partners:
            if col in to_drop or partner in to_drop:
                continue

            corr_val = float(upper.loc[partner, col])

            col_redund = corr[col].drop(index=[col]).mean()
            part_redund = corr[partner].drop(index=[partner]).mean()
            diff = abs(col_redund - part_redund)

            rel_ok = diff <= redundancy_tolerance_rel * max(col_redund, part_redund, 1e-12)
            abs_ok = diff <= redundancy_tolerance_abs

            if abs_ok or rel_ok:
                if (col in prefer_keep) ^ (partner in prefer_keep):
                    kept = col if col in prefer_keep else partner
                    loser = partner if col in prefer_keep else col
                    reason = "prefer_keep_tie"
                else:
                    if col_redund > part_redund:
                        kept, loser = partner, col
                    elif part_redund > col_redund:
                        kept, loser = col, partner
                    else:
                        kept, loser = sorted([col, partner])[0], sorted([col, partner])[1]
                    reason = "more_redundant"
            else:
                if col_redund >= part_redund:
                    kept, loser = partner, col
                else:
                    kept, loser = col, partner
                reason = "more_redundant"

            to_drop.add(loser)
            decisions.append({
                "kept_feature": kept,
                "dropped_feature": loser,
                "abs_correlation": corr_val,
                "reason": reason,
                "col_redundancy": float(col_redund),
                "partner_redundancy": float(part_redund),
                "redundancy_diff": float(diff),
            })

    dropped.extend(sorted(to_drop))

    report = pd.DataFrame(decisions, columns=["kept_feature", "dropped_feature", "abs_correlation", "reason"])
    report = report.sort_values(
        by=["reason", "abs_correlation"],
        ascending=[True, False],
        na_position="last"
    ).reset_index(drop=True)

    if is_adata:
        keep_vars = [v for v in data.var_names if v not in set(dropped)]
        return data[:, keep_vars].copy(), dropped, report
    else:
        reduced_df = data.drop(columns=[c for c in dropped if c in data.columns])
        return reduced_df, dropped, report


def drop_low_variance_features(
    data,
    feature_cols: list[str],
    low_var_threshold= 1e-4,
) -> tuple[object, list[str], list[str]]:
    """
    Drop low-variance features using sklearn VarianceThreshold.

    Accepts either:
      - pandas DataFrame (features are columns)
      - AnnData (features are var_names; uses .X via .to_df())

    Returns:
      reduced_data      : same type as input (DataFrame or AnnData)
      kept_features     : list of kept feature names
      dropped_low_var   : list of dropped feature names
    """
    is_adata = hasattr(data, "var_names") and hasattr(data, "obs_names")

    if is_adata:
        existing = [f for f in feature_cols if f in data.var_names]
        X = data[:, existing].to_df() if existing else pd.DataFrame(index=data.obs_names)
    else:
        existing = [f for f in feature_cols if f in data.columns]
        X = data[existing].copy() if existing else pd.DataFrame(index=data.index)

    if not existing or X.shape[1] == 0:
        return data, [], []

    selector = VarianceThreshold(threshold=low_var_threshold)
    selector.fit(X)

    keep_mask = selector.get_support()
    kept_features = X.columns[keep_mask].tolist()
    dropped_low_var = X.columns[~keep_mask].tolist()

    if is_adata:
        keep_vars = [v for v in data.var_names if v not in set(dropped_low_var)]
        return data[:, keep_vars].copy(), kept_features, dropped_low_var
    else:
        reduced_df = data.drop(columns=dropped_low_var)
        return reduced_df, kept_features, dropped_low_var
    
# def drop_highly_correlated_features(
#     df: pd.DataFrame,
#     feature_cols: list,
#     threshold: float = 0.98,
#     prefer_keep: list | None = None,
#     also_drop_constant: bool = True,
#     redundancy_tolerance_abs: float = 0.02,       # apply prefer_keep only if |Δredundancy| <= this
#     redundancy_tolerance_rel: float = 0.05,
# ) -> tuple[pd.DataFrame, list, pd.DataFrame]:
#     """
#     Remove one feature from any pair with |corr| > threshold.
#     Returns:
#       X_reduced : DataFrame of remaining features
#       dropped   : list of dropped feature names
#       report    : DataFrame logging each decision: kept, dropped, correlation, reason
#     """
#     # prefer_keep = set(prefer_keep or [])
#     prefer_keep = [f for f in feature_cols if "median" in f]
#     prefer_keep += [f for f in feature_cols if "mean" in f]
#     # Work on a numeric-only copy (keep order)
#     X = df[feature_cols].copy()
#     # Cast non-numeric to numeric where possible, otherwise drop
#     cols_ok = [c for c in X.columns if is_numeric_dtype(X[c])]
#     X = X[cols_ok]

#     # Replace infs and optionally drop constants
#     X = X.replace([np.inf, -np.inf], np.nan)
#     dropped: list[str] = []
#     decisions: list[dict] = []

#     if also_drop_constant:
#         nunique = X.nunique(dropna=True)
#         const_cols = nunique[nunique <= 1].index.tolist()
#         for c in const_cols:
#             decisions.append({
#                 "kept_feature": None,
#                 "dropped_feature": c,
#                 "abs_correlation": np.nan,
#                 "reason": "constant_or_single_value"
#             })
#         if const_cols:
#             X = X.drop(columns=const_cols)
#             dropped.extend(const_cols)

#     if X.shape[1] <= 1:
#         report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
#         return X, dropped, report

#     # Impute NaNs (mean) so corr works; PCA later will use its own scaling
#     X_impute = X.fillna(X.mean(numeric_only=True))

#     # Absolute correlation matrix and its upper triangle
#     corr = X_impute.corr().abs()
#     upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

#     to_drop = set()

#     # Greedy pass over pairs above threshold
#     for col in upper.columns:
#         partners = upper.index[upper[col] > threshold].tolist()
#         for partner in partners:
#             if col in to_drop or partner in to_drop:
#                 continue

#             corr_val = float(upper.loc[partner, col])

#             # compute redundancies (mean abs corr to others)
#             col_redund    = corr[col].drop(index=[col]).mean()
#             part_redund   = corr[partner].drop(index=[partner]).mean()
#             diff = abs(col_redund - part_redund)
#             rel_ok = diff <= redundancy_tolerance_rel * max(col_redund, part_redund, 1e-12)
#             abs_ok = diff <= redundancy_tolerance_abs

#             if abs_ok or rel_ok:
#                 # tie → apply prefer_keep if it singles out one; else fall back to redundancy anyway
#                 if (col in prefer_keep) ^ (partner in prefer_keep):
#                     kept  = col if col in prefer_keep else partner
#                     loser = partner if col in prefer_keep else col
#                     reason = "prefer_keep_tie"
#                 else:
#                     # tie but no clear preference → arbitrarily drop the slightly more redundant (or lexical tiebreak)
#                     if col_redund >= part_redund:
#                         kept, loser = partner, col
#                     else:
#                         kept, loser = col, partner
#                     reason = "more_redundant"
#             else:
#                 # not a tie → drop the more redundant one (ignore prefer_keep)
#                 if col_redund >= part_redund:
#                     kept, loser = partner, col
#                 else:
#                     kept, loser = col, partner
#                 reason = "more_redundant"

#             to_drop.add(loser)
#             decisions.append({
#                 "kept_feature": kept,
#                 "dropped_feature": loser,
#                 "abs_correlation": corr_val,
#                 "reason": reason,
#                 "col_redundancy": float(col_redund),
#                 "partner_redundancy": float(part_redund),
#                 "redundancy_diff": float(diff),
#             })

#     if to_drop:
#         df = df.drop(columns=list(to_drop))
#         dropped.extend(list(to_drop))

#     report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
#     # Sort report: constants first, then by correlation desc
#     report = report.sort_values(
#         by=["reason", "abs_correlation"],
#         ascending=[True, False],
#         na_position="last"
#     ).reset_index(drop=True)

#     return df, dropped, report

# def drop_low_variance_features(
#     df: pd.DataFrame,
#     feature_cols: list[str],
#     low_var_threshold: float,
# ) -> tuple[pd.DataFrame, list[str], list[str]]:
#     """
#     Drop low-variance features using sklearn VarianceThreshold.

#     Returns:
#       df_reduced        : df with low-variance features removed
#       kept_features     : list of kept feature names
#       dropped_low_var   : list of dropped feature names
#     """
#     selector = VarianceThreshold(threshold=low_var_threshold)
#     selector.fit(df[feature_cols])

#     keep_mask = selector.get_support()
#     kept_features = df[feature_cols].columns[keep_mask].tolist()
#     dropped_low_var = df[feature_cols].columns[~keep_mask].tolist()

#     df_reduced = df.drop(columns=dropped_low_var)

#     return df_reduced, kept_features, dropped_low_var

