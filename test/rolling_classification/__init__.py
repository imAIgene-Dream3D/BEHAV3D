import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Literal, List
from pandas.api.types import is_numeric_dtype

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

# ---------------------------------------------------------------------
# ---------- Basic descriptive statistic functions ----------
# ---------------------------------------------------------------------

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

# ---------------------------------------------------------------------
# ---------- Robust and temporal statistics ----------
# ---------------------------------------------------------------------

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

# ---------------------------------------------------------------------
# ---------- Type inference and specialized features ----------
# ---------------------------------------------------------------------

def check_if_nonnegative_integer_like(signal_values: np.ndarray) -> bool:
    """Return True if all finite values are non-negative and integer-like."""
    finite_values = np.asarray(signal_values)[np.isfinite(signal_values)]
    return finite_values.size > 0 and np.all(finite_values >= 0) and np.allclose(finite_values, np.round(finite_values))

def convert_to_binary(signal_values: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """Convert a numeric signal to binary (0/1) using a threshold."""
    signal_values = np.asarray(signal_values, dtype=float)
    return np.where(np.isfinite(signal_values), (signal_values > threshold).astype(float), np.nan)

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

# ---------------------------------------------------------------------
# ---------- Main orchestrator: compute all per-window features ----------
# ---------------------------------------------------------------------

def compute_window_features(
    window_dataframe: pd.DataFrame, 
    column_name: str, 
    time_column: str = "position_t",
    signal_types = None,
    ) -> Dict[str, float]:
    """
    Compute descriptive, robust, and temporal features for one signal column in a windowed dataframe.
    Each feature has a self-descriptive name and units inferred from context.
    """
    if time_column in window_dataframe.columns:
        window_dataframe = window_dataframe.sort_values(time_column, kind="mergesort")
        time_values = window_dataframe[time_column].to_numpy(float)
    else:
        time_values = None

    signal_values = window_dataframe[column_name].to_numpy(float)
    signal_type = signal_types.get(column_name, "continuous")
    features: Dict[str, float] = {}

    # --- Basic statistics ---
    features[f"{column_name}_mean_value"] = compute_mean_value(signal_values)
    
    if signal_type != "binary":
        features[f"{column_name}_value_range"] = compute_value_range(signal_values)
        features[f"{column_name}_standard_deviation"] = compute_standard_deviation(signal_values)
        features[f"{column_name}_minimum_value"] = compute_minimum_value(signal_values)
        features[f"{column_name}_maximum_value"] = compute_maximum_value(signal_values)
        features[f"{column_name}_median_value"] = compute_median_value(signal_values)
        
    # --- Robust statistics ---
    if signal_type != "binary":
        features[f"{column_name}_interquartile_range"] = compute_interquartile_range(signal_values)
        features[f"{column_name}_median_absolute_deviation"] = compute_median_absolute_deviation(signal_values)

    # --- Quantiles ---
    if signal_type != "binary":
        for q in [0.1, 0.25, 0.75, 0.9]:
            features[f"{column_name}_quantile_{int(q*100)}percent"] = compute_quantile(signal_values, q)

    # --- Temporal structure ---
    features[f"{column_name}_linear_trend_slope_per_time_unit"] = compute_linear_trend_slope(signal_values, time_values)
    features[f"{column_name}_lag1_autocorrelation"] = compute_lag1_autocorrelation(signal_values)
    if signal_type != "binary":
        features[f"{column_name}_mean_absolute_first_difference"] = compute_mean_absolute_first_difference(signal_values)

    # --- Type-specific extras ---
    if signal_type == "bounded":
        low_frac, high_frac = compute_fraction_near_bounds(signal_values)
        features[f"{column_name}_fraction_near_lower_bound"] = low_frac
        features[f"{column_name}_fraction_near_upper_bound"] = high_frac

    if signal_type == "binary":
        binary_signal = convert_to_binary(signal_values)
        # features[f"{column_name}_proportion_true"] = float(np.nanmean(binary_signal))
        features[f"{column_name}_transition_rate"] = compute_binary_transition_rate(binary_signal)
        features[f"{column_name}_longest_true_run_length"] = compute_longest_true_run_length(binary_signal)

    if signal_type == "count":
        finite = signal_values[np.isfinite(signal_values)]
        if finite.size:
            mean_val = np.mean(finite)
            var_val = np.var(finite)
            features[f"{column_name}_dispersion_index_variance_over_mean"] = var_val / mean_val if mean_val > 0 else np.nan
            features[f"{column_name}_fraction_of_zeros"] = float(np.mean(finite == 0))
        else:
            features[f"{column_name}_dispersion_index_variance_over_mean"] = np.nan
            features[f"{column_name}_fraction_of_zeros"] = np.nan

    features[f"{column_name}_inferred_signal_type"] = signal_type
    return features


def create_windowed_track_dataset(
    df_tracks: pd.DataFrame,
    columns_to_summarize: List[str],
    window_size: Optional[int] = None,                  # None => full-track window
    step_size: Optional[int] = None,             # ignored when window_size is None
    time_col: str = "position_t",
    id_cols: List[str] = ["sample_name", "TrackID"],
    signal_types: Optional[Dict[str, str]] = None,   # pass a global map if you have one
    ) -> pd.DataFrame:
    """
    Split each track into windows and compute descriptive features per window.
    If window_size is None, compute features over the *entire track* (one window per track).

    Parameters
    ----------
    df_tracks : DataFrame
        Must include id_cols + time_col + the columns in columns_to_summarize.
    columns_to_summarize : list[str]
        Columns for which to compute per-window descriptive features.
    window_size : int or None
        Number of rows (frames) per window; if None, use full track as one window.
    step_size : int or None
        Stride between consecutive windows (ignored when window_size is None).
        If None (and window_size is not None), defaults to window_size (non-overlapping windows).
    time_col : str
        Time column used to order and to label the window span in sub_TrackID.
    id_cols : list[str]
        Identifier columns for grouping (default ['sample_name','TrackID']).
    signal_types : dict[str, str] or None
        Optional pre-inferred global signal types per column. If None, will call infer_signal_types.

    Returns
    -------
    DataFrame
        One row per window with:
        - sample_name, TrackID
        - sub_TrackID like "534_t0-t50" (based on time_col values)
        - window_start_<time_col>, window_end_<time_col>, window_length_frames
        - computed features for each requested column
    """
    # sort globally for consistency, then group
    df_sorted = df_tracks.sort_values(id_cols + [time_col], kind="mergesort")

    # Infer signal types once if not given
    if signal_types is None:
        # assumes you have this helper defined elsewhere
        signal_types = infer_signal_types(df_sorted, columns=columns_to_summarize)

    output_rows = []

    for _, group in df_sorted.groupby(id_cols, sort=False):
        group = group.reset_index(drop=True)
        n = len(group)
        if n == 0:
            continue

        # ----- Full-track mode -----
        if window_size is None:
            window_df = group
            start_t = window_df[time_col].iloc[0]
            end_t   = window_df[time_col].iloc[-1]
            track_id_val = window_df["TrackID"].iloc[0]
            sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

            base = {
                id_cols[0]: window_df[id_cols[0]].iloc[0],
                id_cols[1]: track_id_val,
                "sub_TrackID": sub_track_id,
                f"window_start_{time_col}": float(start_t),
                f"window_end_{time_col}": float(end_t),
                "window_length_frames": int(len(window_df)),
            }
            for col in columns_to_summarize:
                base.update(
                    compute_window_features(
                        window_df, col, time_column=time_col, signal_types=signal_types
                    )
                )
            output_rows.append(base)
            continue  # next track

        # ----- Sliding-window mode -----
        if n < window_size:
            continue  # too short for any window

        stride = step_size if step_size is not None else window_size

        for start_idx in range(0, n - window_size + 1, stride):
            end_idx = start_idx + window_size  # exclusive
            window_df = group.iloc[start_idx:end_idx]

            start_t = window_df[time_col].iloc[0]
            end_t   = window_df[time_col].iloc[-1]
            track_id_val = window_df["TrackID"].iloc[0]
            sub_track_id = f"{int(track_id_val)}_t{int(start_t)}-t{int(end_t)}"

            base = {
                id_cols[0]: window_df[id_cols[0]].iloc[0],
                id_cols[1]: track_id_val,
                "sub_TrackID": sub_track_id,
                f"window_start_{time_col}": float(start_t),
                f"window_end_{time_col}": float(end_t),
                "window_length_frames": int(len(window_df)),
            }
            for col in columns_to_summarize:
                base.update(
                    compute_window_features(
                        window_df, col, time_column=time_col, signal_types=signal_types
                    )
                )
            output_rows.append(base)

    result = pd.DataFrame(output_rows)

    # Optional: nice column ordering
    if not result.empty:
        front_cols = [
            id_cols[0], id_cols[1], "sub_TrackID",
            f"window_start_{time_col}", f"window_end_{time_col}", "window_length_frames"
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
        X = X.drop(columns=list(to_drop))
        dropped.extend(list(to_drop))

    report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
    # Sort report: constants first, then by correlation desc
    report = report.sort_values(
        by=["reason", "abs_correlation"],
        ascending=[True, False],
        na_position="last"
    ).reset_index(drop=True)

    return X, dropped, report