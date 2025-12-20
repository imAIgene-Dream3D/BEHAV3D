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

SignalType = Literal["binary", "count", "bounded", "continuous"]


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

def _compute_motion_features_from_xyz(px, py, pz):
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

    return {
        "summed_displacement": path_length,
        "net_displacement": net_disp,
        "straightness": straightness,
        "directional_persistence": mean_persist,
        "median_turning_angle": median_turn,
        "fraction_reversed_movement": frac_reversed,
        "mean_square_displacement": mean_square_displacement,
    }
    
def compute_window_features(window_dataframe, column_name, time_column="position_t", signal_types=None):
    if time_column in window_dataframe.columns:
        time_values = window_dataframe[time_column].to_numpy(float, copy=False)
    else:
        time_values = None

    signal_values = window_dataframe[column_name].to_numpy(float, copy=False)
    signal_type = signal_types.get(column_name, "continuous") if signal_types is not None else "continuous"

    features = {}
    features[f"{column_name}_inferred_signal_type"] = signal_type

    features[f"{column_name}_mean_value"] = compute_mean_value(signal_values)

    if signal_type != "binary":
        features[f"{column_name}_median_value"] = compute_median_value(signal_values)
        features[f"{column_name}_value_range"] = compute_value_range(signal_values)
        features[f"{column_name}_standard_deviation"] = compute_standard_deviation(signal_values)
        features[f"{column_name}_minimum_value"] = compute_minimum_value(signal_values)
        features[f"{column_name}_maximum_value"] = compute_maximum_value(signal_values)
        features[f"{column_name}_interquartile_range"] = compute_interquartile_range(signal_values)
        features[f"{column_name}_median_absolute_deviation"] = compute_median_absolute_deviation(signal_values)

        for q in [0.1, 0.25, 0.75, 0.9]:
            features[f"{column_name}_quantile_{int(q*100)}percent"] = compute_quantile(signal_values, q)

    # 1 is postive slope, 0 is no slope, -1 is decreasing slope
    features[f"{column_name}_linear_trend_slope_per_time_unit"] = compute_linear_trend_slope(
        signal_values, time_values
    )
    
    # lag1_autocorrelation gives how similar a value is to itself a step later
    # 1 is smooth curve, 0 is pure noise (current value gives no information on next value), -1 is up/down pattern
    features[f"{column_name}_lag1_autocorrelation"] = compute_lag1_autocorrelation(signal_values)

    # This gives the actual mean value of the differences in value between one point to the next
    if signal_type != "binary":
        features[f"{column_name}_mean_absolute_first_difference"] = compute_mean_absolute_first_difference(signal_values)

    if signal_type == "bounded":
        low_frac, high_frac = compute_fraction_near_bounds(signal_values)
        features[f"{column_name}_fraction_near_lower_bound"] = low_frac
        features[f"{column_name}_fraction_near_upper_bound"] = high_frac

    if signal_type == "binary":
        binary_signal = convert_to_binary(signal_values)
        features[f"{column_name}_transition_rate"] = compute_binary_transition_rate(binary_signal)
        features[f"{column_name}_longest_true_length"] = compute_longest_true_run_length(binary_signal)
        features[f"{column_name}_mean_true_length"] = compute_average_true_run_length(binary_signal)
        features[f"{column_name}_longest_false_length"] = compute_longest_true_run_length(~binary_signal)
        features[f"{column_name}_mean_false_length"] = compute_average_true_run_length(~binary_signal)

    if signal_type == "count":
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
    group_df, columns_to_summarize, window_size, step_size, time_col, id_cols, signal_types
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

        base.update(_compute_motion_features_from_xyz(px_all, py_all, pz_all))

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
    n_jobs=None,      # None => os.cpu_count()
    chunksize=8,      # tune if many small tracks
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
                chunksize=chunksize,
            )
            for rows in tqdm(futures, total=total_groups):
                if rows:
                    output_rows.extend(rows)
    else:
        for g in tqdm(groups, total=total_groups):
            rows = _create_descriptive_track_worker(
                g, columns_to_summarize, window_size, step_size, time_col, id_cols, signal_types
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

def drop_highly_correlated_features(
    df: pd.DataFrame,
    feature_cols: list,
    threshold: float = 0.98,
    prefer_keep: list | None = None,
    also_drop_constant: bool = True,
    redundancy_tolerance_abs: float = 0.02,       # apply prefer_keep only if |Δredundancy| <= this
    redundancy_tolerance_rel: float = 0.05,
) -> tuple[pd.DataFrame, list, pd.DataFrame]:
    """
    Remove one feature from any pair with |corr| > threshold.
    Returns:
      X_reduced : DataFrame of remaining features
      dropped   : list of dropped feature names
      report    : DataFrame logging each decision: kept, dropped, correlation, reason
    """
    # prefer_keep = set(prefer_keep or [])
    prefer_keep = [f for f in feature_cols if "median" in f]
    prefer_keep += [f for f in feature_cols if "mean" in f]
    # Work on a numeric-only copy (keep order)
    X = df[feature_cols].copy()
    # Cast non-numeric to numeric where possible, otherwise drop
    cols_ok = [c for c in X.columns if is_numeric_dtype(X[c])]
    X = X[cols_ok]

    # Replace infs and optionally drop constants
    X = X.replace([np.inf, -np.inf], np.nan)
    dropped: list[str] = []
    decisions: list[dict] = []

    if also_drop_constant:
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
        report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
        return X, dropped, report

    # Impute NaNs (mean) so corr works; PCA later will use its own scaling
    X_impute = X.fillna(X.mean(numeric_only=True))

    # Absolute correlation matrix and its upper triangle
    corr = X_impute.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

    to_drop = set()

    # Greedy pass over pairs above threshold
    for col in upper.columns:
        partners = upper.index[upper[col] > threshold].tolist()
        for partner in partners:
            if col in to_drop or partner in to_drop:
                continue

            corr_val = float(upper.loc[partner, col])

            # compute redundancies (mean abs corr to others)
            col_redund    = corr[col].drop(index=[col]).mean()
            part_redund   = corr[partner].drop(index=[partner]).mean()
            diff = abs(col_redund - part_redund)
            rel_ok = diff <= redundancy_tolerance_rel * max(col_redund, part_redund, 1e-12)
            abs_ok = diff <= redundancy_tolerance_abs

            if abs_ok or rel_ok:
                # tie → apply prefer_keep if it singles out one; else fall back to redundancy anyway
                if (col in prefer_keep) ^ (partner in prefer_keep):
                    kept  = col if col in prefer_keep else partner
                    loser = partner if col in prefer_keep else col
                    reason = "prefer_keep_tie"
                else:
                    # tie but no clear preference → arbitrarily drop the slightly more redundant (or lexical tiebreak)
                    if col_redund >= part_redund:
                        kept, loser = partner, col
                    else:
                        kept, loser = col, partner
                    reason = "more_redundant"
            else:
                # not a tie → drop the more redundant one (ignore prefer_keep)
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

    if to_drop:
        df = df.drop(columns=list(to_drop))
        dropped.extend(list(to_drop))

    report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
    # Sort report: constants first, then by correlation desc
    report = report.sort_values(
        by=["reason", "abs_correlation"],
        ascending=[True, False],
        na_position="last"
    ).reset_index(drop=True)

    return df, dropped, report

 ### Remove features with no variance
    selector = VarianceThreshold(threshold=1e-4)
    selector.fit(df_analysis[descriptive_feature_cols])

    keep_mask = selector.get_support()
    kept_features = df_analysis[descriptive_feature_cols].columns[keep_mask].tolist()
    dropped_low_var = df_analysis[descriptive_feature_cols].columns[~keep_mask].tolist()