from __future__ import annotations

import os
import re
import shutil
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import napari
import numpy as np
import pandas as pd
import zarr

from behav3d.analysis.backprojection import backproject_columns
from behav3d.core.metadata import load_behav3d_metadata
from behav3d.io.images import load_image


# ---------------------------------------------------------------------
output_dir = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor")
sample_name = "BHVD_SB1_Exp012_Img004"


# Optional explicit feature table. If None, the script searches the BEHAV3D run folder.
feature_csv_path = None

# Notebook note:
# If you are running this file by pasting it into a Jupyter notebook, re-run the
# full cell after editing helper functions. Otherwise the kernel may keep using
# older in-memory function definitions.

# Use None to export and show every backprojectable feature column.
# Example: ["mean_speed", "nr_pixels", "dead"]
# cell_type = "organoid"
# selected_features = [
#     "nr_dead_mask_pixels",
#     "increase_dead_mask",
#     "percentage_dead_mask",
#     "smoothed_percentage_dead_mask",
#     "smoothed_increase_dead_mask",
#     "smoothed_nr_dead_mask_pixels",
# ]

cell_type = "tcell"
selected_features = [
    "elongation",
    "extent",
    "solidity",
    "sphericity",
    "percentage_dead_mask",
    "speed"
    # "active_killing"
]
# Choose between one summarized value per track or true per-timepoint values.
# Allowed values: "summary", "timepoint"
overlay_value_mode = "timepoint"

# Number of worker processes used for frame-wise backprojection.
# Use "auto" to pick a conservative process count automatically.
n_workers = "auto"

# Optional: use None for the full time range, or set (start_t, end_t).
timepoint_range = None

# If False, existing backprojected feature zarrs are reused when present.
refresh_backprojection = True

# Set to an integer to skip categorical columns with more than this number of labels.
# Use None to encode all categorical columns.
max_categorical_levels = 500

raw_channel_colormaps = [
    "cyan",
    "yellow",
    "red",
    "green",
    "magenta",
    "blue",
    "gray",
    "turbo",
    "viridis",
    "plasma",
    "inferno",
    "twilight",
]


def _require_editable_settings():
    if str(output_dir) == "/path/to/behav3d/run":
        raise ValueError("Set 'output_dir' at the top of the script before running it.")
    if not str(sample_name).strip() or sample_name == "paste_sample_name_here":
        raise ValueError("Set 'sample_name' at the top of the script before running it.")
    if not str(cell_type).strip():
        raise ValueError("Set 'cell_type' at the top of the script before running it.")


def _mixed_sort_key(value):
    text = str(value)
    if re.fullmatch(r"-?\d+", text):
        return (0, int(text))
    return (1, text)


def _sanitize_filename(name):
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("._")
    return token or "feature"


def _normalize_overlay_value_mode(mode):
    mode = str(mode).strip().lower()
    if mode in {"summary", "summarized", "mean", "track"}:
        return "summary"
    if mode in {"timepoint", "per_timepoint", "time", "per-timepoint"}:
        return "timepoint"
    raise ValueError(
        "overlay_value_mode must be either 'summary' or 'timepoint'. "
        f"Got '{mode}'."
    )


def _resolve_n_workers(value):
    if isinstance(value, str):
        token = value.strip().lower()
        if token == "auto":
            detected = os.cpu_count() or 1
            return max(1, min(int(detected), 8))
        value = token

    resolved = int(value)
    return max(1, resolved)


def _resolve_metadata_path(run_dir):
    metadata_path = Path(run_dir) / "metadata.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Could not find metadata.csv in run folder: {metadata_path}")
    return metadata_path


def _get_sample_row(metadata, selected_sample_name):
    selected = metadata[metadata["sample_name"].astype(str) == str(selected_sample_name)]
    if selected.empty:
        available = ", ".join(metadata["sample_name"].astype(str).tolist())
        raise ValueError(
            f"Sample '{selected_sample_name}' was not found in metadata.csv. "
            f"Available samples: {available}"
        )
    return selected.iloc[0]


def _resolve_raw_image_path(row, run_dir):
    raw_path = row.get("raw_image_path")
    if pd.notna(raw_path) and str(raw_path).strip():
        raw_path = Path(str(raw_path).strip()).expanduser()
        if raw_path.exists():
            return raw_path

    fallback = Path(run_dir) / "images" / str(row["sample_name"]) / f"{row['sample_name']}.zarr"
    if fallback.exists():
        return fallback

    zipped = Path(str(fallback) + ".zip")
    if zipped.exists():
        return zipped

    raise FileNotFoundError(
        "Could not resolve the raw image path from metadata['raw_image_path'] "
        f"or fallback path '{fallback}'."
    )


def _resolve_prefixed_metadata_path(row, cell_type, suffix):
    for prefix in ("im", "or", "ot"):
        col = f"{prefix}_{cell_type}_{suffix}"
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            path = Path(str(value).strip()).expanduser()
            if path.exists():
                return path
    return None


def _resolve_tracked_image_path(row, run_dir, cell_type):
    tracked_path = _resolve_prefixed_metadata_path(row, cell_type, "tracks_image_path")
    if tracked_path is not None:
        return tracked_path

    sample_dir = Path(run_dir) / "images" / str(row["sample_name"])
    candidates = [
        sample_dir / f"{row['sample_name']}_{cell_type}_tracked.zarr",
        sample_dir / f"{row['sample_name']}_{cell_type}_tracked.zarr.zip",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    fallback = sorted(list(sample_dir.glob(f"*_{cell_type}_*tracked.zarr")) + list(sample_dir.glob(f"*_{cell_type}_*tracked.zarr.zip")))
    if fallback:
        return fallback[0]

    raise FileNotFoundError(
        f"Could not resolve tracked labels for sample '{row['sample_name']}' and cell_type '{cell_type}'."
    )


def _resolve_tracked_csv_path(row, run_dir, cell_type):
    tracked_path = _resolve_prefixed_metadata_path(row, cell_type, "tracks_csv_path")
    if tracked_path is not None:
        return tracked_path

    sample_dir = Path(run_dir) / "trackdata" / str(row["sample_name"]) / str(cell_type)
    candidates = [
        sample_dir / f"{row['sample_name']}_{cell_type}_tracks.csv",
        sample_dir / f"{row['sample_name']}_{cell_type}_track.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    fallback = sorted(list(sample_dir.glob("*_tracks.csv")) + list(sample_dir.glob("*_track.csv")))
    if fallback:
        return fallback[0]

    return None


def _resolve_feature_csv_path(run_dir, sample_name, cell_type, explicit_path=None):
    if explicit_path not in {None, ""}:
        candidate = Path(explicit_path).expanduser()
        if not candidate.exists():
            raise FileNotFoundError(f"Feature CSV does not exist: {candidate}")
        return candidate

    run_dir = Path(run_dir)
    sample_name = str(sample_name).strip()
    cell_type = str(cell_type).strip()

    preferred = [
        run_dir / "trackdata" / sample_name / cell_type / f"{sample_name}_{cell_type}_track_features.csv",
        run_dir / "trackadata" / sample_name / cell_type / f"{sample_name}_{cell_type}_track_features.csv",
    ]
    for candidate in preferred:
        if candidate.exists():
            return candidate

    search_dirs = [
        run_dir / "trackdata" / sample_name / cell_type,
        run_dir / "trackadata" / sample_name / cell_type,
        run_dir / "trackdata" / cell_type,
        run_dir / "trackadata" / cell_type,
    ]
    for search_dir in search_dirs:
        if search_dir.exists():
            matches = sorted(search_dir.glob("*features.csv"))
            if matches:
                return matches[0]

    analysis_fallbacks = [
        run_dir / "analysis" / cell_type / "track_features" / f"BEHAV3D_{cell_type}_combined_track_features.csv",
        run_dir / "analysis" / cell_type / "track_features" / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
    ]
    for candidate in analysis_fallbacks:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "Could not find a feature CSV automatically. "
        f"Checked trackdata/ and analysis/ under '{run_dir}'."
    )


def _resolve_track_id_column(df):
    for col in ("TrackID", "track_id"):
        if col in df.columns:
            return col
    raise ValueError("Feature CSV must contain either 'TrackID' or 'track_id'.")


def _normalize_feature_table(df, selected_sample_name, track_col, time_col="position_t", sample_col="sample_name"):
    if sample_col in df.columns:
        sample_mask = df[sample_col].astype("string").str.strip() == str(selected_sample_name)
        df = df[sample_mask].copy()
    else:
        df = df.copy()
        df[sample_col] = str(selected_sample_name)

    if df.empty:
        raise ValueError(
            f"No feature rows found for sample '{selected_sample_name}' in the supplied feature CSV."
        )

    if time_col not in df.columns:
        raise ValueError(f"Feature CSV is missing required time column '{time_col}'.")

    df[track_col] = pd.to_numeric(df[track_col], errors="coerce")
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    df = df.dropna(subset=[track_col, time_col]).copy()

    if df.empty:
        raise ValueError(
            "Feature CSV has no rows with valid numeric TrackID and position_t values after filtering."
        )

    df[track_col] = df[track_col].astype(np.int64)
    df[time_col] = df[time_col].astype(np.int64)
    df = df.sort_values([sample_col, time_col, track_col]).drop_duplicates(
        subset=[sample_col, time_col, track_col],
        keep="last",
    )
    return df


def _normalize_selected_features(selected):
    if selected is None:
        return None
    if isinstance(selected, str):
        return [selected]
    return [str(col) for col in selected]


def _read_csv_header(path):
    path = Path(path)
    if not path.exists():
        return []
    header_df = pd.read_csv(path, nrows=0)
    return [str(col).lstrip("\ufeff") for col in header_df.columns]


TRACKED_CSV_FALLBACK_COLUMNS = {
    "SegmentID",
    "position_x",
    "position_y",
    "position_z",
    "pixel_position_x",
    "pixel_position_y",
    "pixel_position_z",
    "source_cell_type",
}


def _select_feature_columns(
    df,
    track_col,
    time_col,
    selected=None,
    metadata_columns=None,
    tracked_csv_columns=None,
):
    selected = _normalize_selected_features(selected)
    excluded = set()
    structural_columns = {"sample_name", "TrackID", "track_id", track_col, time_col}
    excluded.update(col for col in structural_columns if col in df.columns)
    if metadata_columns is not None:
        excluded.update(str(col) for col in metadata_columns)
    if tracked_csv_columns is not None:
        excluded.update(str(col) for col in tracked_csv_columns)
    else:
        excluded.update(col for col in TRACKED_CSV_FALLBACK_COLUMNS if col in df.columns)

    if selected is None:
        return [col for col in df.columns if col not in excluded]

    missing = [col for col in selected if col not in df.columns]
    if missing:
        raise ValueError(f"Selected feature columns are missing from the CSV: {missing}")

    filtered = [col for col in selected if col not in excluded]
    skipped = [col for col in selected if col in excluded]
    if skipped:
        print(f"Skipping excluded non-feature columns: {skipped}")
    return filtered


def _coerce_boolean_series(series):
    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")

    text = series.astype("string").str.strip().str.lower()
    non_empty = text[text.notna() & (text != "")]
    if non_empty.empty:
        return None

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
    }
    if set(non_empty.unique().tolist()).issubset(set(mapping)):
        return text.map(mapping).astype("boolean")
    return None


def _prepare_feature_spec(df, feature_col, track_col, time_col, out_dir, feature_csv_path, max_levels=None):
    series = df[feature_col]
    work = df[[track_col, time_col, feature_col]].copy()

    bool_series = _coerce_boolean_series(series)
    label_map = None
    feature_kind = "numeric"
    layer_type = "image"
    background_value = 0
    value_col = "_value"

    if bool_series is not None:
        bool_code_map = {False: 1, True: 2}
        work[value_col] = bool_series.map(bool_code_map).astype("Int8")
        feature_kind = "boolean"
        layer_type = "label"
        out_dtype = np.uint8
        label_map = {"1": "False", "2": "True"}
    else:
        numeric_series = pd.to_numeric(series, errors="coerce")
        non_empty_mask = series.notna()
        if pd.api.types.is_string_dtype(series) or series.dtype == object:
            non_empty_mask = non_empty_mask & series.astype("string").str.strip().ne("")

        if non_empty_mask.any() and numeric_series[non_empty_mask].notna().all():
            work[value_col] = numeric_series.astype(np.float32)
            out_dtype = np.float32
        else:
            text = series.astype("string").str.strip()
            text = text.where(text != "")
            unique_labels = sorted(text.dropna().unique().tolist(), key=_mixed_sort_key)
            if len(unique_labels) == 0:
                return None
            if max_levels is not None and len(unique_labels) > int(max_levels):
                print(
                    f"Skipping '{feature_col}': categorical column has {len(unique_labels)} unique labels "
                    f"(limit={int(max_levels)})."
                )
                return None

            label_to_code = {str(label): int(idx + 1) for idx, label in enumerate(unique_labels)}
            work[value_col] = text.map(label_to_code)
            feature_kind = "categorical"
            layer_type = "label"
            out_dtype = np.uint16 if len(unique_labels) <= int(np.iinfo(np.uint16).max) else np.uint32
            label_map = {str(code): str(label) for label, code in label_to_code.items()}

    work = work.dropna(subset=[value_col]).copy()
    if work.empty:
        return None

    work = work.sort_values([time_col, track_col]).drop_duplicates(
        subset=[time_col, track_col],
        keep="last",
    )
    work[track_col] = work[track_col].astype(np.int64)
    work[time_col] = work[time_col].astype(np.int64)
    work[value_col] = work[value_col].astype(out_dtype, copy=False)

    output_path = Path(out_dir) / f"{_sanitize_filename(feature_col)}.zarr"
    return {
        "feature_name": str(feature_col),
        "feature_kind": str(feature_kind),
        "layer_type": str(layer_type),
        "background_value": background_value,
        "label_map": label_map,
        "encoded_rows": work[[track_col, time_col, value_col]].rename(columns={value_col: feature_col}),
        "output_dtype": np.dtype(out_dtype),
        "output_path": output_path,
        "n_rows": int(len(work)),
        "source_feature_csv_path": str(Path(feature_csv_path).expanduser()),
    }


def _remove_path(path):
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def _summarize_encoded_rows(encoded_rows, feature_name, feature_kind, track_col, time_col):
    work = encoded_rows[[track_col, time_col, feature_name]].copy()
    if feature_kind == "numeric":
        return (
            work.groupby(track_col, as_index=False)[feature_name]
            .mean()
        )

    if feature_kind == "boolean":
        return (
            work.groupby(track_col, as_index=False)[feature_name]
            .max()
        )

    def _mode_or_last(series):
        series = series.dropna()
        if series.empty:
            return np.nan
        mode_vals = pd.Series(series).mode(dropna=True)
        if not mode_vals.empty:
            return mode_vals.iloc[0]
        return series.iloc[-1]

    return (
        work.sort_values([track_col, time_col])
        .groupby(track_col, as_index=False)[feature_name]
        .agg(_mode_or_last)
    )


def _write_feature_attrs(spec, tracked_img_path, sample_name, cell_type, overlay_mode):
    root = zarr.open(str(spec["output_path"]), mode="a")
    arr = root[spec["feature_name"]] if hasattr(root, "__getitem__") else root
    arr.attrs["feature_name"] = spec["feature_name"]
    arr.attrs["feature_kind"] = spec["feature_kind"]
    arr.attrs["layer_type"] = spec["layer_type"]
    arr.attrs["background_value"] = int(spec["background_value"])
    arr.attrs["overlay_mode"] = str(overlay_mode)
    arr.attrs["source_feature_csv_path"] = spec["source_feature_csv_path"]
    arr.attrs["source_tracked_path"] = str(tracked_img_path)
    arr.attrs["sample_name"] = str(sample_name)
    arr.attrs["cell_type"] = str(cell_type)
    if spec["label_map"] is not None:
        arr.attrs["label_map"] = dict(spec["label_map"])


def _export_feature_backprojection_zarrs(
    tracked_labels,
    tracked_img_path,
    feature_specs,
    sample_name,
    cell_type,
    overlay_mode,
    n_workers,
    refresh_existing=True,
):
    if not feature_specs:
        raise ValueError("No backprojectable features were found in the feature CSV.")

    reused = 0
    computed_specs = []

    for spec in feature_specs:
        out_path = Path(spec["output_path"])
        if refresh_existing:
            _remove_path(out_path)
        elif out_path.exists():
            reused += 1
            continue

        out_path.parent.mkdir(parents=True, exist_ok=True)
        feature_name = spec["feature_name"]
        if overlay_mode == "summary":
            lookup_df = _summarize_encoded_rows(
                encoded_rows=spec["encoded_rows"],
                feature_name=feature_name,
                feature_kind=spec["feature_kind"],
                track_col="TrackID" if "TrackID" in spec["encoded_rows"].columns else "track_id",
                time_col="position_t",
            )
        else:
            lookup_df = spec["encoded_rows"]

        backproject_columns(
            track_img=tracked_labels,
            tracked_img_path=tracked_img_path,
            zarr_outpath=out_path,
            df_tracks_clustered=lookup_df,
            columns=[feature_name],
            mode=overlay_mode,
            track_col="TrackID" if "TrackID" in lookup_df.columns else "track_id",
            time_col="position_t",
            background_value=spec["background_value"],
            n_workers=n_workers,
        )
        computed_specs.append(spec)

    for spec in computed_specs:
        _write_feature_attrs(
            spec=spec,
            tracked_img_path=tracked_img_path,
            sample_name=sample_name,
            cell_type=cell_type,
            overlay_mode=overlay_mode,
        )

    if computed_specs:
        print(
            f"Exported {len(computed_specs)} backprojected feature stack(s) "
            f"for sample '{sample_name}' ({cell_type}, mode={overlay_mode})."
        )
    if reused > 0:
        print(f"Reused {reused} existing backprojected feature stack(s).")


def _normalize_tracked_labels_shape(tracked_labels):
    if tracked_labels.ndim == 4:
        return tracked_labels
    if tracked_labels.ndim == 5 and int(tracked_labels.shape[1]) == 1:
        return tracked_labels[:, 0]
    raise ValueError(
        "Tracked labels must have shape (T, Z, Y, X) or (T, 1, Z, Y, X), "
        f"but got shape {tracked_labels.shape}."
    )


def _slice_time_range(arr, selected_range):
    if selected_range is None:
        return arr, None

    start_t, end_t = selected_range
    start_t = max(0, min(int(start_t), int(arr.shape[0]) - 1))
    end_t = max(start_t, min(int(end_t), int(arr.shape[0]) - 1))
    return arr[start_t:end_t + 1], (start_t, end_t)


def _slice_rows_to_time_range(df, time_col, sliced_range):
    if sliced_range is None:
        return df
    start_t, end_t = sliced_range
    df = df[(df[time_col] >= start_t) & (df[time_col] <= end_t)].copy()
    df[time_col] = df[time_col] - start_t
    return df


def _build_track_coords(df, track_col, time_col):
    required_cols = [track_col, time_col, "position_z", "position_y", "position_x"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        return None
    return df[required_cols].to_numpy()


def _load_feature_backprojection_layers(feature_specs, selected_range):
    layers = []
    for spec in feature_specs:
        arr = load_image(spec["output_path"], group=spec["feature_name"])
        arr, _ = _slice_time_range(arr, selected_range)
        layers.append(
            {
                "feature_name": spec["feature_name"],
                "layer_type": spec["layer_type"],
                "feature_kind": spec["feature_kind"],
                "label_map": spec["label_map"],
                "arr": arr,
            }
        )
    return layers


def build_feature_backprojection_payload(
    run_dir,
    selected_sample_name,
    selected_cell_type,
    explicit_feature_csv_path=None,
    selected_feature_columns=None,
    selected_timepoint_range=None,
    overlay_mode="timepoint",
    n_workers="auto",
    refresh_existing=True,
    categorical_level_limit=500,
):
    run_dir = Path(run_dir).expanduser()
    overlay_mode = _normalize_overlay_value_mode(overlay_mode)
    n_workers = _resolve_n_workers(n_workers)
    metadata = load_behav3d_metadata(_resolve_metadata_path(run_dir))
    row = _get_sample_row(metadata, selected_sample_name)

    raw_path = _resolve_raw_image_path(row, run_dir)
    tracked_img_path = _resolve_tracked_image_path(row, run_dir, selected_cell_type)
    tracked_csv_path = _resolve_tracked_csv_path(row, run_dir, selected_cell_type)
    feature_path = _resolve_feature_csv_path(
        run_dir=run_dir,
        sample_name=selected_sample_name,
        cell_type=selected_cell_type,
        explicit_path=explicit_feature_csv_path,
    )

    print(f"Sample selected: {selected_sample_name}")
    print(f"Cell type: {selected_cell_type}")
    print(f"Raw image: {raw_path}")
    print(f"Tracked labels: {tracked_img_path}")
    if tracked_csv_path is not None:
        print(f"Tracked CSV: {tracked_csv_path}")
    print(f"Feature CSV: {feature_path}")
    print(f"Backprojection workers: {n_workers}")

    raw_image = load_image(raw_path)
    if raw_image.ndim != 5:
        raise ValueError(f"Expected raw image with shape (T, C, Z, Y, X), got shape {raw_image.shape}")

    tracked_labels = _normalize_tracked_labels_shape(load_image(tracked_img_path))
    if int(tracked_labels.shape[0]) != int(raw_image.shape[0]):
        raise ValueError(
            "Tracked labels and raw image must have the same number of timepoints, "
            f"but got tracked={tracked_labels.shape}, raw={raw_image.shape}."
        )

    df_features = pd.read_csv(feature_path, low_memory=False)
    df_features.columns = [str(col).lstrip("\ufeff") for col in df_features.columns]
    track_col = _resolve_track_id_column(df_features)
    time_col = "position_t"
    df_features = _normalize_feature_table(
        df=df_features,
        selected_sample_name=selected_sample_name,
        track_col=track_col,
        time_col=time_col,
    )

    feature_columns = _select_feature_columns(
        df=df_features,
        track_col=track_col,
        time_col=time_col,
        selected=selected_feature_columns,
        metadata_columns=metadata.columns,
        tracked_csv_columns=_read_csv_header(tracked_csv_path) if tracked_csv_path is not None else None,
    )
    if not feature_columns:
        raise ValueError("No feature columns were selected for backprojection.")
    print(f"Selected feature columns: {len(feature_columns)}")
    if refresh_existing:
        print("Backprojection cache mode: rebuild zarr outputs")
    else:
        print("Backprojection cache mode: reuse existing zarr outputs when available")

    backprojection_dir = Path(
        run_dir,
        "analysis",
        selected_cell_type,
        "feature_backprojection",
        overlay_mode,
        selected_sample_name,
    )
    feature_specs = []
    for feature_col in feature_columns:
        spec = _prepare_feature_spec(
            df=df_features,
            feature_col=feature_col,
            track_col=track_col,
            time_col=time_col,
            out_dir=backprojection_dir,
            feature_csv_path=feature_path,
            max_levels=categorical_level_limit,
        )
        if spec is not None:
            feature_specs.append(spec)

    _export_feature_backprojection_zarrs(
        tracked_labels=tracked_labels,
        tracked_img_path=tracked_img_path,
        feature_specs=feature_specs,
        sample_name=selected_sample_name,
        cell_type=selected_cell_type,
        overlay_mode=overlay_mode,
        n_workers=n_workers,
        refresh_existing=bool(refresh_existing),
    )

    raw_view, sliced_range = _slice_time_range(raw_image, selected_timepoint_range)
    tracked_view, _ = _slice_time_range(tracked_labels, sliced_range)
    df_tracks_view = _slice_rows_to_time_range(df_features, time_col=time_col, sliced_range=sliced_range)
    feature_layers = _load_feature_backprojection_layers(feature_specs, selected_range=sliced_range)

    if sliced_range is not None:
        print(f"Showing timepoints {sliced_range[0]} to {sliced_range[1]}")

    return {
        "sample_name": str(selected_sample_name),
        "cell_type": str(selected_cell_type),
        "raw_path": Path(raw_path),
        "tracked_img_path": Path(tracked_img_path),
        "feature_csv_path": Path(feature_path),
        "raw_image": raw_view,
        "tracked_labels": tracked_view,
        "feature_layers": feature_layers,
        "df_tracks": df_tracks_view,
        "track_col": str(track_col),
        "time_col": str(time_col),
        "overlay_mode": str(overlay_mode),
        "n_workers": int(n_workers),
        "scale_4d": (
            1,
            float(row["pixel_distance_z"]),
            float(row["pixel_distance_xy"]),
            float(row["pixel_distance_xy"]),
        ),
    }


def launch_feature_backprojection_viewer(payload, run=True):
    viewer = napari.Viewer()

    raw_image = payload["raw_image"]
    scale_4d = payload["scale_4d"]
    for ch in range(raw_image.shape[1]):
        viewer.add_image(
            raw_image[:, ch],
            name=f"raw_channel_{ch}",
            colormap=raw_channel_colormaps[ch % len(raw_channel_colormaps)],
            scale=scale_4d,
            blending="additive",
        )

    viewer.add_labels(
        payload["tracked_labels"],
        name=f"{payload['cell_type']} tracked labels",
        scale=scale_4d,
        visible=False,
    )

    track_coords = _build_track_coords(
        payload["df_tracks"],
        track_col=payload["track_col"],
        time_col=payload["time_col"],
    )
    if track_coords is not None and len(track_coords) > 0:
        viewer.add_tracks(
            track_coords,
            name=f"{payload['cell_type']} tracks",
            tail_length=20,
        )

    for idx, layer_payload in enumerate(payload["feature_layers"]):
        layer_name = layer_payload["feature_name"]
        metadata = {
            "feature_kind": layer_payload["feature_kind"],
            "label_map": layer_payload["label_map"],
        }
        is_visible = idx == 0
        if layer_payload["layer_type"] == "label":
            viewer.add_labels(
                layer_payload["arr"],
                name=layer_name,
                scale=scale_4d,
                metadata=metadata,
                visible=is_visible,
            )
        else:
            viewer.add_image(
                layer_payload["arr"],
                name=layer_name,
                scale=scale_4d,
                colormap="inferno",
                opacity=0.75,
                blending="translucent",
                metadata=metadata,
                visible=is_visible,
            )

    print(
        f"Opened feature backprojection viewer for sample '{payload['sample_name']}' "
        f"({payload['cell_type']}) with {len(payload['feature_layers'])} feature layer(s)."
    )
    print(
        "Backprojection matching key: "
        f"(sample_name='{payload['sample_name']}', {payload['time_col']}, {payload['track_col']})"
    )

    if run:
        napari.run()
    return viewer


def show_feature_backprojection(
    run_dir,
    selected_sample_name,
    selected_cell_type,
    explicit_feature_csv_path=None,
    selected_feature_columns=None,
    selected_timepoint_range=None,
    overlay_mode="timepoint",
    n_workers="auto",
    refresh_existing=True,
    categorical_level_limit=500,
):
    payload = build_feature_backprojection_payload(
        run_dir=run_dir,
        selected_sample_name=selected_sample_name,
        selected_cell_type=selected_cell_type,
        explicit_feature_csv_path=explicit_feature_csv_path,
        selected_feature_columns=selected_feature_columns,
        selected_timepoint_range=selected_timepoint_range,
        overlay_mode=overlay_mode,
        n_workers=n_workers,
        refresh_existing=refresh_existing,
        categorical_level_limit=categorical_level_limit,
    )
    return launch_feature_backprojection_viewer(payload, run=True)


if __name__ == "__main__":
    _require_editable_settings()
    show_feature_backprojection(
        run_dir=output_dir,
        selected_sample_name=sample_name,
        selected_cell_type=cell_type,
        explicit_feature_csv_path=feature_csv_path,
        selected_feature_columns=selected_features,
        selected_timepoint_range=timepoint_range,
        overlay_mode=overlay_value_mode,
        n_workers=n_workers,
        refresh_existing=refresh_backprojection,
        categorical_level_limit=max_categorical_levels,
    )
