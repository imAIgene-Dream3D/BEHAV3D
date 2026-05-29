from pathlib import Path

import dask.array as da
import napari
import numpy as np
import pandas as pd

from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    load_behav3d_metadata,
)
from behav3d.io.images import load_image


# ---------------------------------------------------------------------
output_dir = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor")
sample_name = "BHVD_SB1_Exp012_Img001"

# Optional: limit the viewer to specific cell types. Use None to load all.
selected_cell_types = ["tcell1"]

# Optional: use None for the full time range, or set (start_t, end_t).
timepoint_range = None

# Choose what value to fill each tracked object with:
# - "nr_pixels": size in voxels
# - "volume": physical size using pixel_distance_z * pixel_distance_xy * pixel_distance_xy
size_value_mode = "nr_pixels"

# Optional size text at centroids.
show_size_text = False
size_text_decimals = 1

size_layer_colormap = "inferno"
size_layer_opacity = 0.65

show_track_length_layer = True
show_track_length_text = False
track_length_layer_colormap = "viridis"
track_length_layer_opacity = 0.55

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

# Override mask and track colors by cell/organoid name only.
name_colors = {

    "tcell1": "cyan",
}


def _require_editable_settings():
    if str(output_dir) == "/path/to/behav3d/run":
        raise ValueError("Set 'output_dir' at the top of the script before running it.")
    if not str(sample_name).strip() or sample_name == "paste_sample_name_here":
        raise ValueError("Set 'sample_name' at the top of the script before running it.")
    if size_value_mode not in {"nr_pixels", "volume"}:
        raise ValueError("size_value_mode must be either 'nr_pixels' or 'volume'.")


def _patch_napari_animation_labels_state():
    try:
        from napari_animation.viewer_state import ViewerState
    except Exception:
        return

    if getattr(ViewerState, "_behav3d_labels_state_patch", False):
        return

    @classmethod
    def _from_viewer_without_label_colormaps(cls, viewer):
        layers = {}
        for layer in viewer.layers:
            layer_state = dict(layer.as_layer_data_tuple()[1])
            layer_state.pop("metadata", None)
            if layer.__class__.__name__ == "Labels":
                layer_state.pop("colormap", None)
            layers[layer.name] = layer_state

        return cls(
            camera=viewer.camera.dict(),
            dims=viewer.dims.dict(),
            layers=layers,
        )

    ViewerState.from_viewer = _from_viewer_without_label_colormaps
    ViewerState._behav3d_labels_state_patch = True


def _resolve_metadata_path(run_dir):
    metadata_path = Path(run_dir) / "metadata.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"Could not find metadata.csv in run folder: {metadata_path}")
    return metadata_path


def _resolve_raw_image_path(row, run_dir):
    raw_path = row.get("raw_image_path")
    if pd.notna(raw_path) and str(raw_path).strip():
        raw_path = Path(str(raw_path).strip())
        if raw_path.exists():
            return raw_path

    fallback = Path(run_dir) / "images" / str(row["sample_name"]) / f"{row['sample_name']}.zarr"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        "Could not resolve the raw image path from metadata['raw_image_path'] "
        f"or fallback path '{fallback}'."
    )


def _get_sample_row(metadata, selected_sample_name):
    selected = metadata[metadata["sample_name"].astype(str) == str(selected_sample_name)]
    if selected.empty:
        available = ", ".join(metadata["sample_name"].astype(str).tolist())
        raise ValueError(
            f"Sample '{selected_sample_name}' was not found in metadata.csv. "
            f"Available samples: {available}"
        )
    return selected.iloc[0]


def _iter_cell_types(metadata):
    for cell_type in detect_organoid_types_from_metadata(metadata):
        yield "or", cell_type
    for cell_type in detect_immune_cell_types_from_metadata(metadata):
        yield "im", cell_type
    for cell_type in detect_other_cell_types_from_metadata(metadata):
        yield "ot", cell_type


def _resolve_name_color(cell_type):
    return name_colors.get(cell_type, "white")


def _resolve_tracked_inputs(row, metadata):
    tracked_inputs = []
    selected = None if selected_cell_types is None else {str(x) for x in selected_cell_types}

    for prefix, cell_type in _iter_cell_types(metadata):
        if selected is not None and str(cell_type) not in selected:
            continue

        track_img_col = f"{prefix}_{cell_type}_tracks_image_path"
        track_csv_col = f"{prefix}_{cell_type}_tracks_csv_path"

        track_img_val = row.get(track_img_col)
        track_csv_val = row.get(track_csv_col)

        if not (pd.notna(track_img_val) and str(track_img_val).strip()):
            print(f"Skipping {cell_type}: missing {track_img_col}")
            continue
        if not (pd.notna(track_csv_val) and str(track_csv_val).strip()):
            print(f"Skipping {cell_type}: missing {track_csv_col}")
            continue

        track_img_path = Path(str(track_img_val).strip())
        track_csv_path = Path(str(track_csv_val).strip())

        if not track_img_path.exists():
            print(f"Skipping {cell_type}: tracked labels file not found: {track_img_path}")
            continue
        if not track_csv_path.exists():
            print(f"Skipping {cell_type}: track CSV file not found: {track_csv_path}")
            continue

        tracked_inputs.append(
            {
                "prefix": prefix,
                "cell_type": cell_type,
                "track_img_path": track_img_path,
                "track_csv_path": track_csv_path,
                "color": _resolve_name_color(cell_type),
            }
        )

    return tracked_inputs


def _slice_time_range(arr, selected_range):
    if selected_range is None:
        return arr, None

    start_t, end_t = selected_range
    start_t = max(0, min(int(start_t), arr.shape[0] - 1))
    end_t = max(start_t, min(int(end_t), arr.shape[0] - 1))
    return arr[start_t:end_t + 1], (start_t, end_t)


def _slice_tracks(df_tracks, sliced_range):
    if sliced_range is None:
        return df_tracks

    start_t, end_t = sliced_range
    df_tracks = df_tracks[
        (df_tracks["position_t"] >= start_t) & (df_tracks["position_t"] <= end_t)
    ].copy()
    df_tracks["position_t"] = df_tracks["position_t"] - start_t
    return df_tracks


def _resolve_track_id_column(df_tracks):
    for col in ("TrackID", "track_id"):
        if col in df_tracks.columns:
            return col
    raise ValueError("Track CSV must contain either 'TrackID' or 'track_id'.")


def _build_track_coords(df_tracks):
    track_id_col = _resolve_track_id_column(df_tracks)
    required_cols = [track_id_col, "position_t", "position_z", "position_y", "position_x"]
    missing_cols = [col for col in required_cols if col not in df_tracks.columns]
    if missing_cols:
        raise ValueError(f"Track CSV is missing required columns: {missing_cols}")
    return df_tracks[required_cols].to_numpy()


def _build_mask_labels(tracked_labels):
    return (tracked_labels > 0).astype(np.uint8)


def _apply_mask_layer_color(layer, color):
    color_map = {
        0: (0, 0, 0, 0),
        1: color,
    }
    if hasattr(layer, "_set_colormap"):
        layer._set_colormap(color_map)
    layer.color = color_map
    layer.color_mode = "direct"
    layer.refresh()
    if hasattr(layer, "update"):
        layer.update()
    if hasattr(layer, "_update_thumbnail"):
        layer._update_thumbnail()


def _ensure_dask_array(arr):
    if isinstance(arr, da.Array):
        return arr
    return da.from_array(arr, chunks=(1,) + tuple(arr.shape[1:]))


def _rechunk_full_timepoints(arr):
    arr = _ensure_dask_array(arr)
    if arr.ndim != 4:
        raise ValueError(f"Expected tracked labels with shape (T, Z, Y, X), got {arr.shape}")
    return arr.rechunk((1,) + tuple(int(v) for v in arr.shape[1:]))


def _labels_to_lookup_block(block, lookup_ids, lookup_values, out_dtype):
    out = np.zeros(block.shape, dtype=out_dtype)
    lookup_ids = np.asarray(lookup_ids, dtype=np.int64)
    lookup_values = np.asarray(lookup_values, dtype=out_dtype)

    for idx in range(block.shape[0]):
        labels = np.asarray(block[idx], dtype=np.int64)
        if labels.size == 0:
            continue
        unique_labels, inverse = np.unique(labels, return_inverse=True)
        if unique_labels.size == 0:
            continue

        mapped_unique = np.zeros(unique_labels.shape[0], dtype=out_dtype)
        positive_mask = unique_labels > 0
        if np.any(positive_mask) and lookup_ids.size > 0:
            positive_labels = unique_labels[positive_mask]
            positions = np.searchsorted(lookup_ids, positive_labels)
            valid = (positions < lookup_ids.size) & (lookup_ids[positions] == positive_labels)
            if np.any(valid):
                mapped_unique[np.flatnonzero(positive_mask)[valid]] = lookup_values[positions[valid]]

        out[idx] = mapped_unique[inverse].reshape(labels.shape)

    return out


def _labels_to_size_block(block, value_mode, voxel_volume):
    out_dtype = np.float32 if value_mode == "volume" else np.uint32
    out = np.zeros(block.shape, dtype=out_dtype)

    for idx in range(block.shape[0]):
        labels = np.asarray(block[idx], dtype=np.int64)
        if labels.size == 0:
            continue
        counts = np.bincount(labels.ravel())
        if counts.size <= 1:
            continue
        values = counts.astype(np.float32 if value_mode == "volume" else np.uint32, copy=True)
        values[0] = 0
        if value_mode == "volume":
            values = values * float(voxel_volume)
        out[idx] = values[labels]

    return out


def _build_lazy_size_image(tracked_labels, value_mode, voxel_volume):
    tracked_labels = _rechunk_full_timepoints(tracked_labels)
    out_dtype = np.float32 if value_mode == "volume" else np.uint32
    return da.map_blocks(
        _labels_to_size_block,
        tracked_labels,
        value_mode=value_mode,
        voxel_volume=voxel_volume,
        dtype=out_dtype,
        meta=np.array((), dtype=out_dtype),
    )


def _build_lazy_lookup_image(tracked_labels, lookup_map, out_dtype):
    tracked_labels = _rechunk_full_timepoints(tracked_labels)
    if not lookup_map:
        lookup_ids = np.array([], dtype=np.int64)
        lookup_values = np.array([], dtype=out_dtype)
    else:
        lookup_ids = np.asarray(sorted(int(k) for k in lookup_map), dtype=np.int64)
        lookup_values = np.asarray([lookup_map[int(k)] for k in lookup_ids], dtype=out_dtype)
    return da.map_blocks(
        _labels_to_lookup_block,
        tracked_labels,
        lookup_ids=lookup_ids,
        lookup_values=lookup_values,
        out_dtype=out_dtype,
        dtype=out_dtype,
        meta=np.array((), dtype=out_dtype),
    )


def _build_size_table(tracked_labels, df_tracks, value_mode, voxel_volume):
    track_id_col = _resolve_track_id_column(df_tracks)
    rows = []

    for t in sorted(df_tracks["position_t"].dropna().astype(int).unique()):
        labels = np.asarray(tracked_labels[t], dtype=np.int64)
        if labels.size == 0:
            continue

        counts = np.bincount(labels.ravel())
        if counts.size <= 1:
            continue

        active_ids = np.flatnonzero(counts[1:] > 0) + 1
        if active_ids.size == 0:
            continue

        size_values = counts[active_ids].astype(np.float32)
        if value_mode == "volume":
            size_values = size_values * float(voxel_volume)

        rows.append(
            pd.DataFrame(
                {
                    track_id_col: active_ids.astype(np.int64),
                    "position_t": np.full(active_ids.shape[0], t, dtype=np.int64),
                    "nr_pixels": counts[active_ids].astype(np.int64),
                    "size_value": size_values,
                }
            )
        )

    if rows:
        return pd.concat(rows, ignore_index=True), track_id_col
    return pd.DataFrame(columns=[track_id_col, "position_t", "nr_pixels", "size_value"]), track_id_col


def _build_track_length_table(df_tracks_full):
    track_id_col = _resolve_track_id_column(df_tracks_full)
    df_times = (
        df_tracks_full[[track_id_col, "position_t"]]
        .dropna()
        .assign(position_t=lambda df: df["position_t"].astype(int))
    )
    df_lengths = (
        df_times
        .groupby(track_id_col, as_index=False)
        .agg(position_t_min=("position_t", "min"), position_t_max=("position_t", "max"))
    )
    df_lengths["track_length_tp"] = (
        df_lengths["position_t_max"] - df_lengths["position_t_min"] + 1
    ).astype(np.uint32)
    df_lengths = df_lengths[[track_id_col, "track_length_tp"]]
    df_lengths[track_id_col] = df_lengths[track_id_col].astype(np.int64)
    track_length_map = dict(zip(df_lengths[track_id_col], df_lengths["track_length_tp"]))
    return df_lengths, track_length_map, track_id_col


def _format_size_text(value):
    if size_value_mode == "nr_pixels":
        return f"{int(round(float(value)))} px"
    return f"{float(value):.{int(size_text_decimals)}f}"


def _format_track_length_text(value):
    return f"{int(round(float(value)))} tp"


def _build_text_table(df_tracks, df_values, track_id_col, value_col, text_col, formatter):
    required_cols = ["position_t", "pixel_position_z", "pixel_position_y", "pixel_position_x"]
    if not all(col in df_tracks.columns for col in required_cols):
        return pd.DataFrame()

    df_points = (
        df_tracks[[track_id_col, "position_t", "pixel_position_z", "pixel_position_y", "pixel_position_x"]]
        .drop_duplicates([track_id_col, "position_t"])
        .merge(df_values, on=[track_id_col, "position_t"], how="inner")
        .copy()
    )
    if df_points.empty:
        return df_points

    df_points[text_col] = df_points[value_col].map(formatter)
    return df_points


def _build_track_length_text_table(df_tracks, df_lengths, track_id_col):
    if df_lengths.empty:
        return pd.DataFrame()

    df_visible_lengths = (
        df_tracks[[track_id_col, "position_t"]]
        .drop_duplicates([track_id_col, "position_t"])
        .merge(df_lengths[[track_id_col, "track_length_tp"]], on=track_id_col, how="inner")
    )
    return _build_text_table(
        df_tracks=df_tracks,
        df_values=df_visible_lengths,
        track_id_col=track_id_col,
        value_col="track_length_tp",
        text_col="track_length_text",
        formatter=_format_track_length_text,
    )


def _build_size_text_table(df_tracks, df_sizes, track_id_col):
    return _build_text_table(
        df_tracks=df_tracks,
        df_values=df_sizes,
        track_id_col=track_id_col,
        value_col="size_value",
        text_col="size_text",
        formatter=_format_size_text,
    )


def _add_text_layer(viewer, df_points, cell_type, layer_suffix, text_col, show_text, scale_4d, color):
    if df_points.empty or not show_text:
        return

    coords = df_points[["position_t", "pixel_position_z", "pixel_position_y", "pixel_position_x"]].to_numpy()
    try:
        viewer.add_points(
            coords,
            name=f"{cell_type} {layer_suffix}",
            scale=scale_4d,
            size=0.01,
            face_color="transparent",
            edge_color="transparent",
            features=df_points,
            text={
                "string": f"{{{text_col}}}",
                "size": 10,
                "color": color,
                "anchor": "upper_left",
            },
        )
    except Exception as exc:
        print(f"Could not add {layer_suffix} for {cell_type}: {exc}")


def _get_positive_contrast_limits(arr):
    if isinstance(arr, da.Array):
        return None
    positive = np.asarray(arr)[np.asarray(arr) > 0]
    if positive.size == 0:
        return (0, 1)
    lo = float(positive.min())
    hi = float(positive.max())
    if lo == hi:
        hi = lo + 1
    return (lo, hi)


def visualize_sample_tracking_sizes(
    run_dir,
    selected_sample_name,
    selected_timepoint_range=None,
    channel_colormaps=None,
):
    _patch_napari_animation_labels_state()

    run_dir = Path(run_dir).expanduser()
    metadata = load_behav3d_metadata(_resolve_metadata_path(run_dir))
    row = _get_sample_row(metadata, selected_sample_name)

    raw_path = _resolve_raw_image_path(row, run_dir)
    tracked_inputs = _resolve_tracked_inputs(row, metadata)

    if not tracked_inputs:
        raise ValueError(
            f"No tracked cell types with both tracked labels and track CSV were found for sample '{selected_sample_name}'."
        )

    print(f"Sample selected: {selected_sample_name}")
    print(f"Raw image: {raw_path}")

    raw_image = load_image(raw_path)
    if raw_image.ndim != 5:
        raise ValueError(f"Expected raw image with shape (T, C, Z, Y, X), got shape {raw_image.shape}")

    raw_image, sliced_range = _slice_time_range(raw_image, selected_timepoint_range)
    if sliced_range is not None:
        print(f"Showing timepoints {sliced_range[0]} to {sliced_range[1]}")

    elsize_xy = float(row["pixel_distance_xy"])
    elsize_z = float(row["pixel_distance_z"])
    voxel_volume = elsize_z * elsize_xy * elsize_xy
    size_unit = "px" if size_value_mode == "nr_pixels" else "volume"
    scale_4d = (1, elsize_z, elsize_xy, elsize_xy)

    viewer = napari.Viewer()

    for ch in range(raw_image.shape[1]):
        viewer.add_image(
            raw_image[:, ch],
            name=f"raw_channel_{ch}",
            colormap=(channel_colormaps or raw_channel_colormaps)[ch % len(channel_colormaps or raw_channel_colormaps)],
            scale=scale_4d,
            blending="additive",
        )

    for tracked_input in tracked_inputs:
        cell_type = tracked_input["cell_type"]
        color = tracked_input["color"]

        tracked_labels = load_image(tracked_input["track_img_path"])
        tracked_labels, _ = _slice_time_range(tracked_labels, sliced_range)
        tracked_labels = _rechunk_full_timepoints(tracked_labels)

        df_tracks_full = pd.read_csv(tracked_input["track_csv_path"])
        df_track_lengths, track_length_map, track_id_col = _build_track_length_table(df_tracks_full)
        df_tracks = _slice_tracks(df_tracks_full, sliced_range)

        size_image = _build_lazy_size_image(
            tracked_labels=tracked_labels,
            value_mode=size_value_mode,
            voxel_volume=voxel_volume,
        )
        size_layer_kwargs = {
            "name": f"{cell_type} size ({size_unit})",
            "scale": scale_4d,
            "colormap": size_layer_colormap,
            "opacity": size_layer_opacity,
            "blending": "translucent",
        }
        contrast_limits = _get_positive_contrast_limits(size_image)
        if contrast_limits is not None:
            size_layer_kwargs["contrast_limits"] = contrast_limits

        track_length_image = None
        track_length_layer_kwargs = {
            "name": f"{cell_type} track length (tp)",
            "scale": scale_4d,
            "colormap": track_length_layer_colormap,
            "opacity": track_length_layer_opacity,
            "blending": "translucent",
        }
        if show_track_length_layer:
            track_length_image = _build_lazy_lookup_image(
                tracked_labels=tracked_labels,
                lookup_map=track_length_map,
                out_dtype=np.uint32,
            )
            contrast_limits = _get_positive_contrast_limits(track_length_image)
            if contrast_limits is not None:
                track_length_layer_kwargs["contrast_limits"] = contrast_limits

        df_sizes = pd.DataFrame()
        df_size_points = pd.DataFrame()
        if show_size_text:
            df_sizes, _ = _build_size_table(
                tracked_labels=tracked_labels,
                df_tracks=df_tracks,
                value_mode=size_value_mode,
                voxel_volume=voxel_volume,
            )
            df_size_points = _build_size_text_table(df_tracks, df_sizes, track_id_col)

        df_track_length_points = pd.DataFrame()
        if show_track_length_text:
            df_track_length_points = _build_track_length_text_table(
                df_tracks=df_tracks,
                df_lengths=df_track_lengths,
                track_id_col=track_id_col,
            )

        viewer.add_image(
            size_image,
            **size_layer_kwargs,
        )

        if show_track_length_layer:
            viewer.add_image(
                track_length_image,
                **track_length_layer_kwargs,
            )

        viewer.add_labels(
            tracked_labels,
            name=f"{cell_type} tracked labels",
            scale=scale_4d,
            visible=False,
        )

        viewer.add_tracks(
            _build_track_coords(df_tracks),
            name=f"{cell_type} tracks",
            colormap=color,
            tail_length=20,
        )

        mask_layer = viewer.add_labels(
            _build_mask_labels(tracked_labels),
            name=f"{cell_type} mask",
            opacity=0.25,
            scale=scale_4d,
        )
        _apply_mask_layer_color(mask_layer, color)
        _add_text_layer(
            viewer=viewer,
            df_points=df_size_points,
            cell_type=cell_type,
            layer_suffix="size text",
            text_col="size_text",
            show_text=show_size_text,
            scale_4d=scale_4d,
            color=color,
        )
        _add_text_layer(
            viewer=viewer,
            df_points=df_track_length_points,
            cell_type=cell_type,
            layer_suffix="track length text",
            text_col="track_length_text",
            show_text=show_track_length_text,
            scale_4d=scale_4d,
            color=color,
        )

        if show_size_text and not df_sizes.empty:
            print(
                f"{cell_type}: n={len(df_sizes)} objects, "
                f"min={df_sizes['size_value'].min():.2f}, "
                f"median={df_sizes['size_value'].median():.2f}, "
                f"max={df_sizes['size_value'].max():.2f}"
            )
        elif show_size_text:
            print(f"{cell_type}: no tracked objects found in the selected range.")
        else:
            print(f"{cell_type}: size overlay is lazy-loaded from {tracked_input['track_img_path']}")

        if show_track_length_text and not df_track_lengths.empty:
            print(
                f"{cell_type}: track lengths min={df_track_lengths['track_length_tp'].min():.0f}, "
                f"median={df_track_lengths['track_length_tp'].median():.0f}, "
                f"max={df_track_lengths['track_length_tp'].max():.0f}"
            )
        elif show_track_length_layer:
            print(f"{cell_type}: track length overlay is lazy-loaded from {tracked_input['track_csv_path']}")

    return viewer


if __name__ == "__main__":
    _require_editable_settings()
    viewer = visualize_sample_tracking_sizes(
        run_dir=output_dir,
        selected_sample_name=sample_name,
        selected_timepoint_range=timepoint_range,
        channel_colormaps=raw_channel_colormaps,
    )
    napari.run()
