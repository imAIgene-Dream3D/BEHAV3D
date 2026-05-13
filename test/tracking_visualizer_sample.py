from pathlib import Path

import napari
import pandas as pd

from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    load_behav3d_metadata,
)
from behav3d.io.images import load_image


# ---------------------------------------------------------------------
# Editable settings
# Paste your BEHAV3D run folder and sample name here.
# ---------------------------------------------------------------------
output_dir = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/SafetyProfiling")
sample_name = "BHVD_SB1_Exp010_Img004_MDOvs27T_HighBME"

# Optional: use None for the full time range, or set (start_t, end_t).
timepoint_range = None

# Raw channel colors are applied in order of channel index.
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
# Example: {"tcell": "cyan", "organoid1": "magenta"}
name_colors = {
    "27T": "magenta",
    "MDO": "yellow",
    "TEG":"cyan"
}


def _require_editable_settings():
    if str(output_dir) == "/path/to/behav3d/run":
        raise ValueError("Set 'output_dir' at the top of the script before running it.")
    if not str(sample_name).strip() or sample_name == "paste_sample_name_here":
        raise ValueError("Set 'sample_name' at the top of the script before running it.")


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
    return name_colors.get(cell_type)


def _resolve_tracked_inputs(row, metadata):
    tracked_inputs = []

    for prefix, cell_type in _iter_cell_types(metadata):
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

        color = _resolve_name_color(cell_type)
        if color is None:
            print(f"Skipping {cell_type}: no color configured in name_colors")
            continue
        tracked_inputs.append(
            {
                "prefix": prefix,
                "cell_type": cell_type,
                "track_img_path": track_img_path,
                "track_csv_path": track_csv_path,
                "color": color,
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
    return (tracked_labels > 0).astype("uint8")


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


def visualize_sample_tracking(
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

        df_tracks = pd.read_csv(tracked_input["track_csv_path"])
        df_tracks = _slice_tracks(df_tracks, sliced_range)

        viewer.add_labels(
            tracked_labels,
            name=f"{cell_type} tracked labels",
            scale=scale_4d,
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
            opacity=0.35,
            scale=scale_4d,
        )
        _apply_mask_layer_color(mask_layer, color)

    return viewer


if __name__ == "__main__":
    _require_editable_settings()
    viewer = visualize_sample_tracking(
        run_dir=output_dir,
        selected_sample_name=sample_name,
        selected_timepoint_range=timepoint_range,
        channel_colormaps=raw_channel_colormaps,
    )
    napari.run()
