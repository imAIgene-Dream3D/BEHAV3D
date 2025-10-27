import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Literal, List
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
import math

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from dtaidistance import dtw, dtw_ndim
import seaborn as sns
from sklearn.cluster import KMeans, HDBSCAN
import umap


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

def plot_umap_feature_grid(
    df: pd.DataFrame,
    feature_cols: list[str],
    x_col: str = "UMAP1",
    y_col: str = "UMAP2",
    ncols: int = 4,
    max_plots: int | None = None,
    point_size: int = 5,
    alpha: float = 0.5,
    numeric_cmap: str = "viridis",
    categorical_palette: str = "tab20",
    add_colorbar: bool = True,
    page: int = 0,   # for pagination: 0-based page index
    ):
    """
    Creates a multi-row, multi-column grid of UMAP scatterplots colored by each feature in feature_cols.
    Filters out non-scalar or missing features automatically. Supports pagination via `page`.
    """
    # Filter valid features (exist, scalar, not all NaN)
    valid = []
    for c in feature_cols:
        if c in df.columns and df[c].notna().any():
            valid.append(c)
    if max_plots is not None:
        valid = valid[:max_plots]

    if len(valid) == 0:
        raise ValueError("No valid features to plot.")

    n = len(valid)
    nrows = math.ceil(n / ncols)

    # Pagination support: choose a slice of features per page
    per_page = nrows * ncols
    start = page * per_page
    end = min(start + per_page, len(valid))
    feats = valid[start:end]
    if len(feats) == 0:
        raise ValueError(f"No features to plot on page {page} (only {math.ceil(len(valid)/per_page)} page(s) available).")

    # Axes limits shared across panels
    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    # Build grid
    fig, axes = plt.subplots(
        nrows=math.ceil(len(feats)/ncols),
        ncols=ncols,
        figsize=(4*ncols, 3.5*math.ceil(len(feats)/ncols)),
        squeeze=False,
        constrained_layout=True
    )

    for i, feat in enumerate(feats):
        r, c = divmod(i, ncols)
        ax = axes[r, c]

        s = df[feat]
        # Numeric vs categorical handling
        if is_numeric_dtype(s):
            # Numeric: use matplotlib scatter for easy colorbar handling
            sc = ax.scatter(
                df[x_col], df[y_col],
                s=point_size, alpha=alpha,
                c=s, cmap=numeric_cmap, edgecolors="none"
            )
            if add_colorbar:
                cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cb.ax.tick_params(labelsize=8)
        else:
            # Categorical: enforce category dtype and use seaborn palette
            s_cat = s.astype("category")
            tmp = df.copy()
            tmp[feat] = s_cat
            sns.scatterplot(
                data=tmp, x=x_col, y=y_col, hue=feat,
                palette=categorical_palette, s=point_size, alpha=alpha,
                legend=False, ax=ax
            )

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_title(feat, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Hide any unused axes (if grid not full)
    total_cells = axes.size
    for j in range(len(feats), total_cells):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    # Add a common title
    fig.suptitle("UMAP colored by features", fontsize=14)
    plt.show()     
    
def compute_dtw_window_clusters(
    df_tracks: pd.DataFrame,
    df_windows: pd.DataFrame,
    features: list,
    non_binary_features: list = None,
    min_cluster_size: int = 50,
    umap_n_neighbors: int = 15,
    umap_min_dist: float = 0.1,
    random_state: int = 123,
    sample_frac: float = None,     # e.g. 0.25 to use 25% of windows
    max_windows: int = None        # e.g. 5000 to cap windows for DTW
):
    """
    Returns a DataFrame with keys + DTW_UMAP1/2 + cluster_label_dtw_hdbscan
    Keys will be the intersection of:
      ['sample_name','TrackID','sub_TrackID','window_start_position_t','window_end_position_t']
    present in df_windows.
    """
    # ---- 1) choose keys available to join on ----
    possible_keys = ["sample_name","TrackID","sub_TrackID","window_start_position_t","window_end_position_t"]
    join_keys = [k for k in possible_keys if k in df_windows.columns]

    if len(join_keys) < 2:
        raise ValueError("Need at least two window-identifying columns to merge (e.g. sample_name, TrackID, sub_TrackID).")

    # ---- 2) scale only non-binary (non-contact) features globally on the track table ----
    if non_binary_features is None:
        non_binary_features = [c for c in features if "contact" not in c]

    df_tracks_scaled = df_tracks.copy()
    if len(non_binary_features) > 0:
        scaler = StandardScaler()
        df_tracks_scaled[non_binary_features] = scaler.fit_transform(df_tracks_scaled[non_binary_features])

    # ---- 3) (optional) subsample windows to control DTW O(n^2) cost ----
    dfw = df_windows[join_keys].drop_duplicates().copy()
    if sample_frac is not None:
        dfw = dfw.sample(frac=sample_frac, random_state=random_state)
    if (max_windows is not None) and (len(dfw) > max_windows):
        dfw = dfw.sample(n=max_windows, random_state=random_state)

    # ---- 4) build per-window sequences using your window start/end in frames ----
    sequences = []
    meta_rows = []   # same length as sequences, for later merge
    # For faster filtering, pre-index df_tracks_scaled by (sample_name, TrackID)
    grp = df_tracks_scaled.groupby(["sample_name", "TrackID"], sort=False)

    for _, w in dfw.iterrows():
        # filter rows for this window
        sn = w.get("sample_name")
        tid = w.get("TrackID")
        # fall back to full table if group missing (shouldn't happen)
        g = grp.get_group((sn, tid)) if (sn, tid) in grp.groups else df_tracks_scaled

        # time slicing
        start_t = w.get("window_start_position_t", None)
        end_t   = w.get("window_end_position_t", None)
        if start_t is not None and end_t is not None:
            mask = (g["position_t"] >= start_t) & (g["position_t"] <= end_t)
            gwin = g.loc[mask]
        else:
            # if start/end not present, just use all rows of the sub_TrackID (or entire track)
            if "sub_TrackID" in join_keys and "sub_TrackID" in g.columns:
                gwin = g[g["sub_TrackID"] == w["sub_TrackID"]]
            else:
                gwin = g

        seq = gwin[features].to_numpy(dtype=np.float64)
        # need at least 2 timepoints for DTW to be meaningful
        if seq.shape[0] < 2:
            continue

        sequences.append(seq)
        meta_rows.append(w.to_dict())

    if len(sequences) < 2:
        raise ValueError("Not enough valid windows to run DTW (need >= 2).")

    meta = pd.DataFrame(meta_rows)

    # ---- 5) DTW distance matrix (multivariate) ----
    # Note: uses fast C implementation; still O(n^2)
    D = dtw_ndim.distance_matrix_fast(sequences)  # shape (n_windows, n_windows)

    # ---- 6) UMAP on precomputed distances ----
    um = umap.UMAP(
        n_components=2,
        n_neighbors=umap_n_neighbors,
        min_dist=umap_min_dist,
        init="random",
        metric="precomputed",
        random_state=random_state,
    )
    emb = um.fit_transform(D)  # (n_windows, 2)

    # ---- 7) HDBSCAN clustering on the embedding (or on D via metric='precomputed') ----
    # Clustering on the embedding is common for noisy high-dim distances.
    hdb = HDBSCAN(min_cluster_size=min_cluster_size, metric="euclidean")
    labels = hdb.fit_predict(emb)

    out = meta.copy()
    out["DTW_UMAP1"] = emb[:, 0]
    out["DTW_UMAP2"] = emb[:, 1]
    out["cluster_label_dtw_hdbscan"] = labels

    return out


def plot_clustering_feature_heatmap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 7,
    nr_cols = 2,
    figsize = (8.27, 11.69),
    plot_results=True,
    show_points=False,
    point_alpha=0.5,
    point_size=8,
    mean_marker_size=60,
):
    """
    Produce a PDF with:
      • Page 1: full-page min–max scaled heatmap of cluster means.
      • Subsequent pages: per-feature violin plots tiled across pages.
    """

    info_cols   = list(info_cols) if info_cols is not None else []
    sample_cols = list(sample_cols) if sample_cols is not None else []

    # ---- Helper ----
    def _round_legend_ticks(max_val):
        try:
            return round_legend_ticks(max_val)
        except Exception:
            if not np.isfinite(max_val) or max_val <= 0:
                return 1.0
            magnitude = 10.0 ** np.floor(np.log10(max_val))
            return float(np.ceil(max_val / magnitude) * magnitude)

    # ---- Cluster means ----
    df_for_means = (
        df_umap[list(info_cols) + ["ClusterID"]]
        .drop(columns=sample_cols, errors="ignore")
    )
    cluster_means = (
        df_for_means
        .groupby("ClusterID", observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    # ---- Min-max scaling ----
    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]

    X = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    all_nan_cols = X.columns[X.isna().all()].tolist()
    if all_nan_cols:
        X = X.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]

    if len(scale_columns) > 0:
        X_filled = X.copy()
        med = X_filled.median(numeric_only=True)
        X_filled = X_filled.fillna(med)
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])
        df_heatmap_scaled = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    # ---- Prepare violin plot data ----
    value_cols = [c for c in info_cols if c not in set(sample_cols)]
    df_values = df_umap[["ClusterID"] + value_cols].copy()
    for c in value_cols:
        df_values[c] = pd.to_numeric(df_values[c], errors="coerce")
    df_values.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_values.melt(id_vars="ClusterID", var_name="var", value_name="value")

    cluster_order = sorted(df_values["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_plots = len(feat_names)
    rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
    nr_pages = max(1, (rows_for_plots + rows_per_page - 1) // rows_per_page)

    with PdfPages(outpath) as pdf:
        # ---- Page 1: Full-page heatmap ----
        fig, ax = plt.subplots(figsize=figsize)
        if not overall_heatmap_data.empty:
            try:
                heatmap = sns.heatmap(
                    overall_heatmap_data,
                    ax=ax,
                    cmap="viridis",
                    cbar=True,
                    yticklabels=True
                )
                ax.set_title("Min–Max Scaled Cluster Means", fontsize=14, pad=14)
                ax.set_xlabel("ClusterID", fontsize=10)
                ax.set_ylabel("", fontsize=10)
                ax.tick_params(axis="y", labelsize=6)
                ax.tick_params(axis="x", labelsize=8, rotation=0)
                cbar = heatmap.collections[0].colorbar
                cbar.ax.tick_params(labelsize=8)
                fig.tight_layout(pad=2.0)
            except Exception:
                ax.text(0.5, 0.5, "Overview heatmap unavailable", ha="center", va="center")
                ax.axis("off")
        else:
            ax.text(0.5, 0.5, "No features available for overview scaling", ha="center", va="center")
            ax.axis("off")

        pdf.savefig(fig, dpi=600)
        plt.close(fig)

        # ---- Remaining pages: Violin plots ----
        plot_idx = 0
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(rows_per_page)
                for c in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()
                    continue

                feat = feat_names[plot_idx]
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].dropna(subset=["ClusterID", "value"])
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no finite data)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                try:
                    sns.violinplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        inner=None,
                        ax=ax,
                        cut=0
                    )
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                if show_points:
                    sns.stripplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        ax=ax,
                        dodge=False,
                        jitter=0.2,
                        alpha=point_alpha,
                        size=point_size
                    )

                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)),
                    means.values,
                    s=mean_marker_size,
                    edgecolor="black",
                    linewidths=0.8,
                    zorder=3
                )

                ax.set_title(feat, fontsize=9)
                ax.set_xlabel("ClusterID", fontsize=8)
                ax.set_ylabel("Value", fontsize=8)
                ax.tick_params(axis="x", rotation=0, labelsize=7)
                ax.tick_params(axis="y", labelsize=7)
                plot_idx += 1

            fig.subplots_adjust(left=0.20, right=0.98, top=0.95, bottom=0.08)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    if plot_results:
        print(f"Saved PDF to: {outpath}")
        