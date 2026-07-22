import re
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _vinfo(verbose, prefix, message):
    if not bool(verbose):
        return
    print(f"[{prefix}] INFO {str(message)}")


def _vstart(verbose, prefix, step_name):
    if bool(verbose):
        print(f"[{prefix}] START {step_name}")
    return time.perf_counter()


def _vdone(verbose, prefix, step_name, t_start):
    if bool(verbose):
        dt = time.perf_counter() - float(t_start)
        print(f"[{prefix}] DONE {step_name} | took {dt:.2f}s")


def _vsave(verbose, prefix, label, path):
    if not bool(verbose):
        return
    print(f"[{prefix}] SAVED {str(label)}: {path}")


def _resolve_output_dir(output_dir):
    if output_dir is None:
        raise ValueError("output_dir is required.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return output_dir_path


def _mixed_label_sort_key(value):
    text = str(value)
    if re.fullmatch(r"-?\d+", text):
        return (0, int(text))
    return (1, text)


def _sanitize_filename_token(value, fallback="value"):
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    token = token.strip("._-")
    return token if len(token) > 0 else str(fallback)


def _handle_nan_in_distance_matrix(distance_matrix, context="distance matrix"):
    is_dataframe = isinstance(distance_matrix, pd.DataFrame)
    if is_dataframe:
        nan_count = distance_matrix.isna().sum().sum()
        if nan_count == 0:
            return distance_matrix
        print(f"  ⚠️ Warning: {context} contains {nan_count} NaN values")
        print(f"  → Replacing with maximum finite value")
        max_val = distance_matrix.max().max()
        if pd.isna(max_val):
            print(f"  ⚠️ Warning: All values in {context} are NaN.")
            print("     → Falling back to 0 for all distances. Please check input data and DTW configuration.")
            return distance_matrix.fillna(0)
        return distance_matrix.fillna(max_val)
    else:
        if not np.isnan(distance_matrix).any():
            return distance_matrix
        nan_count = np.isnan(distance_matrix).sum()
        print(f"  ⚠️ Warning: {context} contains {nan_count} NaN values")
        print(f"  → Replacing with maximum finite value")
        max_val = np.nanmax(distance_matrix)
        if np.isnan(max_val):
            print(f"  ⚠️ Warning: All values in {context} are NaN.")
            print("     → Falling back to 0 for all distances. Please check input data and DTW configuration.")
            return np.nan_to_num(distance_matrix, nan=0)
        return np.nan_to_num(distance_matrix, nan=max_val)


def _save_adata_obs_csv(adata, h5ad_path):
    """Write adata.obs to a .csv file next to the given h5ad path (same stem)."""
    h5ad_path = Path(h5ad_path)
    csv_path = h5ad_path.with_suffix(".csv")
    adata.obs.to_csv(csv_path, index=True)
    return csv_path


def _to_numpy_2d(X):
    if hasattr(X, "toarray"):
        X = X.toarray()
    arr = np.asarray(X)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape={arr.shape}.")
    return arr
