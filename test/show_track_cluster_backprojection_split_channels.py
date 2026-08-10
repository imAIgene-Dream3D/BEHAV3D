from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import napari
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from behav3d.analysis.behavior.state.utils import (
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    _add_mapping_dock_widget,
    _apply_state_code_colors_to_layer,
    _build_code_map,
    _build_state_code_color_map,
    _build_state_mapping_text,
    _resolve_raw_image_path,
    _resolve_tracked_image_path,
)
from behav3d.analysis.behavior.track.visualization.backprojection import (
    add_track_cluster_trajectory_layers,
    backproject_track_cluster_at_timepoint,
    build_track_cluster_frame_assignments,
    build_track_state_sequence_lookup,
    prepare_track_cluster_trajectory_data,
    prepare_track_cluster_code_lookup,
    render_track_statebar_image,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
    merge_coordinate_columns_into_obs,
    resolve_exemplar_positions_csv_path,
)
from behav3d.analysis.behavior.track.utils import (
    get_dtaidistance_track_trajectories_filename,
)
from behav3d.io.images import load_image


# %%
# Interactive configuration
# Fill these in when running the script line by line from the IDE.
# Point this to the BEHAV3D results folder that contains `images/`, `analysis/`,
# and typically `metadata.csv`.
INTERACTIVE_BEHAV3D_FOLDER = (
    r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\BRIM_AW1_DMG_microglia\exp91"
)
# Optional explicit override if your BEHAV3D output directory differs.
INTERACTIVE_OUTPUT_DIR = None
# Leave sample None to auto-pick the first sample found in metadata/images.
INTERACTIVE_SAMPLE_NAME = None
# For this exp91 run, "macrophage" matches the microglia channel. Change to
# "tcell" if you want the T-cell trajectory classes instead.
INTERACTIVE_CELL_TYPE = "macrophage"
INTERACTIVE_CLUSTER_COL = "ClusterID"
# State column shown in the clicked-track state-bar dock.
INTERACTIVE_STATE_COL = "behavioral_state"
INTERACTIVE_OUTPUT_COL = "track_behavioral_cluster"
INTERACTIVE_SHOW_TRAJECTORIES = True
INTERACTIVE_RAW_IMAGE_PATH = None
INTERACTIVE_TRACKED_IMG_PATH = None
INTERACTIVE_TRACK_BP_IMG_PATH = None
INTERACTIVE_TRACK_ADATA_PATH = None
INTERACTIVE_STATE_ADATA_PATH = None
INTERACTIVE_TRACK_LAYER_OPACITY = 0.35
INTERACTIVE_CLASS_LAYER_OPACITY = 0.85
INTERACTIVE_VERBOSE = True


_CHANNEL_COLORS = ["cyan", "yellow", "green", "red", "blue", "magenta"]


def _bp_add_raw_channels(viewer, raw_img, sample_name, saved_channels, dim_order="TCZYX"):
    """Add one napari Image layer per raw channel without importing the GUI module."""
    import dask.array as da

    if not isinstance(raw_img, da.Array):
        raw_img = da.from_array(raw_img)
    if dim_order != "TCZYX" and len(dim_order) == 5:
        try:
            axes = [dim_order.index(d) for d in "TCZYX"]
            raw_img = da.transpose(raw_img, axes)
        except ValueError:
            pass

    n_channels = raw_img.shape[1] if raw_img.ndim >= 5 else 1
    added = []
    for channel_index in range(n_channels):
        channel_data = raw_img[:, channel_index, :, :, :] if raw_img.ndim >= 5 else raw_img
        saved = saved_channels.get(channel_index) or saved_channels.get(str(channel_index))
        color = (
            saved["colormap"]
            if saved and "colormap" in saved
            else _CHANNEL_COLORS[channel_index % len(_CHANNEL_COLORS)]
        )
        layer_name = f"{sample_name} - Ch{channel_index}"
        add_kwargs = dict(name=layer_name, colormap=color, blending="additive", visible=True)
        if saved and "contrast_limits" in saved:
            add_kwargs["contrast_limits"] = tuple(saved["contrast_limits"])
        viewer.add_image(channel_data, **add_kwargs)
        added.append(layer_name)
    return added


def _align_labels_to_raw_shape_for_view_lazy(labels_img, raw_img, layer_name, verbose=False):
    """Align a labels image to the raw viewer TZYX shape without forcing a full load."""
    try:
        import dask.array as da
    except Exception:
        da = None

    raw_shape = tuple(int(v) for v in raw_img.shape)
    if len(raw_shape) not in {4, 5}:
        raise ValueError(f"Unsupported raw image shape {raw_shape} for viewer alignment.")

    labels = da.asarray(labels_img) if da is not None else np.asarray(labels_img)
    if labels.ndim == 5:
        labels = labels[:, 0, ...]
    elif labels.ndim != 4:
        raise ValueError(f"Unsupported {layer_name} shape {labels.shape} for viewer alignment.")

    if tuple(int(v) for v in labels.shape[-3:]) != tuple(raw_shape[-3:]):
        raise ValueError(
            f"{layer_name}/raw spatial mismatch: labels_shape={labels.shape}, raw_shape={raw_shape}."
        )

    raw_T = int(raw_shape[0])
    labels_T = int(labels.shape[0])
    if labels_T == raw_T:
        if verbose:
            print(f"{layer_name} view align: shape already TZYX {tuple(int(v) for v in labels.shape)}")
        return labels

    if labels_T < raw_T:
        pad_shape = (raw_T - labels_T,) + tuple(int(v) for v in labels.shape[1:])
        pad = (
            da.zeros(pad_shape, dtype=labels.dtype)
            if da is not None
            else np.zeros(pad_shape, dtype=labels.dtype)
        )
        labels = da.concatenate([labels, pad], axis=0) if da is not None else np.concatenate([labels, pad], axis=0)
    else:
        labels = labels[:raw_T]

    if verbose:
        print(
            f"{layer_name} view align: T adjusted {labels_T}->{raw_T}, "
            f"final shape {tuple(int(v) for v in labels.shape)}"
        )
    return labels


def _current_viewer_frame(viewer) -> int:
    try:
        return int(viewer.dims.current_step[0])
    except Exception:
        return 0


def _normalize_track_id_key(value) -> str:
    if pd.isna(value):
        return ""
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(numeric):
        numeric = float(numeric)
        if numeric.is_integer():
            return str(int(numeric))
    return str(value).strip()


def _build_track_cluster_lookup(adata_tracks, sample_col="sample_name", track_col="TrackID", cluster_col="ClusterID"):
    if adata_tracks is None or not hasattr(adata_tracks, "obs"):
        return {}
    required_cols = [str(sample_col), str(track_col), str(cluster_col)]
    if any(col not in adata_tracks.obs.columns for col in required_cols):
        return {}

    window_cols = [col for col in ["position_t_min", "position_t_max"] if col in adata_tracks.obs.columns]
    obs = adata_tracks.obs[required_cols + window_cols].copy()
    obs["__sample"] = obs[str(sample_col)].astype("string").str.strip()
    obs["__track"] = obs[str(track_col)].map(_normalize_track_id_key)
    obs["__cluster"] = obs[str(cluster_col)].astype("string").str.strip()
    obs = obs.dropna(subset=["__sample", "__track", "__cluster"]).copy()
    obs = obs[(obs["__sample"] != "") & (obs["__track"] != "")].copy()
    lookup = {}
    for (sample, track), group in obs.groupby(["__sample", "__track"], observed=True, sort=False):
        rows = []
        for _, row in group.iterrows():
            tmin = pd.to_numeric(pd.Series([row.get("position_t_min", pd.NA)]), errors="coerce").iloc[0]
            tmax = pd.to_numeric(pd.Series([row.get("position_t_max", pd.NA)]), errors="coerce").iloc[0]
            rows.append(
                {
                    "cluster": str(row["__cluster"]),
                    "tmin": None if pd.isna(tmin) else float(tmin),
                    "tmax": None if pd.isna(tmax) else float(tmax),
                }
            )
        lookup[(str(sample), str(track))] = rows
    return lookup


def _resolve_track_cluster_label(cluster_lookup, sample_key, track_key_value, cursor_time=None):
    rows = cluster_lookup.get((sample_key, track_key_value), [])
    if len(rows) == 0:
        return "unknown"
    if cursor_time is not None:
        matches = []
        for row in rows:
            tmin = row.get("tmin")
            tmax = row.get("tmax")
            if tmin is None or tmax is None:
                matches.append(str(row.get("cluster", "unknown")))
                continue
            lo = min(float(tmin), float(tmax))
            hi = max(float(tmin), float(tmax))
            if lo <= float(cursor_time) <= hi:
                matches.append(str(row.get("cluster", "unknown")))
        matches = sorted({match for match in matches if match != ""})
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            return "multiple"
    labels = sorted(
        {
            str(row.get("cluster", "unknown"))
            for row in rows
            if str(row.get("cluster", "")).strip() != ""
        }
    )
    if len(labels) == 1:
        return labels[0]
    if len(labels) > 1:
        return "multiple"
    return "unknown"


def _layer_value_at_world_position(layer, world_position):
    data_position = layer.world_to_data(world_position) if hasattr(layer, "world_to_data") else world_position
    data_position = tuple(int(round(float(value))) for value in data_position)
    data = layer.data
    shape = tuple(int(value) for value in getattr(data, "shape", ()))
    if len(shape) == 0 or len(data_position) < len(shape):
        return None, data_position
    index = tuple(
        max(0, min(int(pos), int(size) - 1))
        for pos, size in zip(data_position[:len(shape)], shape)
    )
    value = data[index]
    if hasattr(value, "compute"):
        value = value.compute()
    value = pd.Series([value]).iloc[0]
    try:
        value = int(float(value))
    except Exception:
        value = str(value)
    return value, index


def _set_label_pixmap_from_rgb(label, rgb_img):
    from qtpy.QtGui import QImage, QPixmap

    rgb_img = np.ascontiguousarray(rgb_img, dtype=np.uint8)
    height, width, channels = rgb_img.shape
    if channels != 3:
        raise ValueError(f"Expected RGB image with 3 channels, got shape={rgb_img.shape}.")
    fmt_rgb = getattr(QImage, "Format_RGB888", None)
    if fmt_rgb is None:
        fmt_rgb = QImage.Format.Format_RGB888
    qimg = QImage(rgb_img.data, int(width), int(height), int(3 * width), fmt_rgb).copy()
    label.setPixmap(QPixmap.fromImage(qimg))


def _is_track_cluster_trajectory_layer(layer, output_col: str | None = None) -> bool:
    if layer is None:
        return False
    layer_name = str(getattr(layer, "name", "")).strip()
    if output_col is not None and layer_name.startswith(f"{output_col} trajectory: "):
        return True
    return hasattr(layer, "tail_length") and hasattr(layer, "tail_width") and hasattr(layer, "data")


def _crop_uniform_canvas_border(img, tolerance: int = 2):
    arr = np.asarray(img)
    if arr.ndim != 3 or arr.shape[0] == 0 or arr.shape[1] == 0:
        return arr

    corners = [
        tuple(int(v) for v in arr[0, 0]),
        tuple(int(v) for v in arr[0, -1]),
        tuple(int(v) for v in arr[-1, 0]),
        tuple(int(v) for v in arr[-1, -1]),
    ]
    background_tuple = max(set(corners), key=corners.count)
    background = np.asarray(background_tuple, dtype=np.int16)
    work = arr.astype(np.int16, copy=False)
    bg_mask = np.all(np.abs(work - background) <= int(tolerance), axis=-1)
    row_has_data = ~np.all(bg_mask, axis=1)
    col_has_data = ~np.all(bg_mask, axis=0)
    if (not np.any(row_has_data)) or (not np.any(col_has_data)):
        return arr.copy()

    top = int(np.argmax(row_has_data))
    bottom = int(len(row_has_data) - np.argmax(row_has_data[::-1]))
    left = int(np.argmax(col_has_data))
    right = int(len(col_has_data) - np.argmax(col_has_data[::-1]))
    if bottom <= top or right <= left:
        return arr.copy()
    return arr[top:bottom, left:right].copy()


def _default_screenshot_path(output_dir: Path, sample_name: str, cell_type: str, time_index: int) -> Path:
    screenshots_dir = (
        Path(output_dir)
        / "analysis"
        / str(cell_type)
        / "behavorial_trajectories"
        / "screenshots"
    )
    screenshots_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = (
        f"{sample_name}_{cell_type}_track_clusters_t{int(time_index):04d}_{timestamp}.png"
    )
    return screenshots_dir / filename


def _add_track_viewer_tools_dock(
    viewer,
    *,
    output_dir: Path,
    sample_name: str,
    cell_type: str,
    output_col: str,
    trajectory_layers=None,
    title="Trajectory Tools",
):
    from qtpy.QtCore import QSignalBlocker
    from qtpy.QtWidgets import (
        QDoubleSpinBox,
        QFileDialog,
        QFormLayout,
        QHBoxLayout,
        QLabel,
        QPushButton,
        QSpinBox,
        QVBoxLayout,
        QWidget,
    )

    from napari.utils.io import imsave

    trajectory_layers = list(trajectory_layers or [])

    def _active_trajectory_layer():
        active = viewer.layers.selection.active
        if active in trajectory_layers:
            return active
        for layer in trajectory_layers:
            if layer in viewer.layers:
                return layer
        return None

    info_label = QLabel("")
    info_label.setWordWrap(True)
    tail_length_spin = QSpinBox()
    tail_length_spin.setRange(1, 100000)
    tail_length_spin.setSingleStep(1)
    tail_width_spin = QDoubleSpinBox()
    tail_width_spin.setDecimals(1)
    tail_width_spin.setRange(0.5, 50.0)
    tail_width_spin.setSingleStep(0.5)
    read_button = QPushButton("Load Active Layer Values")
    apply_button = QPushButton("Apply To All Trajectory Layers")
    screenshot_button = QPushButton("Save Cropped Current View")
    status_label = QLabel("")
    status_label.setWordWrap(True)

    form_layout = QFormLayout()
    form_layout.addRow("Tail length (timepoints)", tail_length_spin)
    form_layout.addRow("Tail width (px)", tail_width_spin)

    button_row = QHBoxLayout()
    button_row.addWidget(read_button)
    button_row.addWidget(apply_button)

    widget = QWidget()
    layout = QVBoxLayout()
    layout.addWidget(info_label)
    layout.addLayout(form_layout)
    layout.addLayout(button_row)
    layout.addWidget(screenshot_button)
    layout.addWidget(status_label)
    widget.setLayout(layout)

    def _set_controls_enabled(enabled: bool):
        tail_length_spin.setEnabled(bool(enabled))
        tail_width_spin.setEnabled(bool(enabled))
        read_button.setEnabled(bool(enabled))
        apply_button.setEnabled(bool(enabled))

    def _update_layer_info(*_args):
        layer = _active_trajectory_layer()
        if layer is None:
            _set_controls_enabled(False)
            info_label.setText("No trajectory layers are currently available in this viewer.")
            return

        _set_controls_enabled(True)
        with QSignalBlocker(tail_length_spin):
            tail_length_spin.setValue(int(getattr(layer, "tail_length", 1)))
        with QSignalBlocker(tail_width_spin):
            tail_width_spin.setValue(float(getattr(layer, "tail_width", 2.0)))
        info_label.setText(
            f"Active trajectory layer: {layer.name}<br>"
            f"Blending: {getattr(layer, 'blending', 'opaque')}<br>"
            f"Tail length: {int(getattr(layer, 'tail_length', 1))} timepoints | "
            f"Tail width: {float(getattr(layer, 'tail_width', 2.0)):.1f} px"
        )

    def _apply_to_all_layers():
        if len(trajectory_layers) == 0:
            status_label.setText("No trajectory layers to update.")
            return
        new_tail_length = int(tail_length_spin.value())
        new_tail_width = float(tail_width_spin.value())
        for layer in trajectory_layers:
            if layer not in viewer.layers:
                continue
            layer.tail_length = new_tail_length
            layer.tail_width = new_tail_width
            layer.blending = "opaque"
        status_label.setText(
            f"Applied tail length={new_tail_length} timepoints and tail width={new_tail_width:.1f} px "
            f"to {len([layer for layer in trajectory_layers if layer in viewer.layers])} trajectory layers."
        )
        _update_layer_info()

    def _save_screenshot():
        default_path = _default_screenshot_path(
            output_dir=Path(output_dir),
            sample_name=sample_name,
            cell_type=cell_type,
            time_index=_current_viewer_frame(viewer),
        )
        selected_path, _selected_filter = QFileDialog.getSaveFileName(
            widget,
            "Save Cropped Current View",
            str(default_path),
            "PNG Files (*.png);;TIFF Files (*.tif *.tiff);;All Files (*)",
        )
        if selected_path in {None, ""}:
            status_label.setText("Screenshot save cancelled.")
            return

        screenshot = viewer.screenshot(canvas_only=True, flash=False)
        cropped = _crop_uniform_canvas_border(screenshot, tolerance=2)
        imsave(selected_path, cropped)
        status_label.setText(
            f"Saved cropped screenshot for t={_current_viewer_frame(viewer)} to '{selected_path}'."
        )

    read_button.clicked.connect(_update_layer_info)
    apply_button.clicked.connect(_apply_to_all_layers)
    screenshot_button.clicked.connect(_save_screenshot)
    viewer.layers.selection.events.active.connect(_update_layer_info)
    for layer in trajectory_layers:
        if hasattr(layer, "events"):
            if hasattr(layer.events, "tail_length"):
                layer.events.tail_length.connect(_update_layer_info)
            if hasattr(layer.events, "tail_width"):
                layer.events.tail_width.connect(_update_layer_info)

    _update_layer_info()
    viewer.window.add_dock_widget(widget, area="right", name=str(title))
    return widget


def _add_track_statebar_click_dock_compat(
    viewer,
    *,
    sample_name,
    adata_full,
    adata_tracks=None,
    track_layer_name="filtered TrackID",
    clickable_layer_name="behavioral_state_class",
    state_col="behavioral_state",
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    cluster_col="ClusterID",
    title="Clicked Track State Bar",
):
    from qtpy.QtWidgets import QLabel, QVBoxLayout, QWidget

    layer_names = [layer.name for layer in viewer.layers]
    if track_layer_name not in layer_names:
        raise ValueError(f"Napari viewer is missing required layer '{track_layer_name}'.")
    track_layer = viewer.layers[track_layer_name]
    clickable_layers = []
    if clickable_layer_name in layer_names:
        clickable_layers.append(viewer.layers[clickable_layer_name])
    if track_layer not in clickable_layers:
        clickable_layers.append(track_layer)

    from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
        _build_state_color_map,
    )

    _state_values, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=str(state_col),
        cmap_name="tab20",
    )
    state_lookup = build_track_state_sequence_lookup(
        adata_full=adata_full,
        sample_name=sample_name,
        state_col=state_col,
        sample_col=sample_col,
        track_col=track_col,
        time_col=time_col,
    )
    cluster_lookup = _build_track_cluster_lookup(
        adata_tracks=adata_tracks,
        sample_col=sample_col,
        track_col=track_col,
        cluster_col=cluster_col,
    )

    info_label = QLabel("Click a tracked cell to show its full state bar.")
    info_label.setWordWrap(True)
    image_label = QLabel("")
    image_label.setMinimumWidth(480)
    widget = QWidget()
    layout = QVBoxLayout()
    layout.addWidget(info_label)
    layout.addWidget(image_label)
    widget.setLayout(layout)
    viewer.window.add_dock_widget(widget, area="right", name=str(title))

    sample_key = str(sample_name).strip()

    def _show_message(text):
        info_label.setText(str(text))
        image_label.clear()

    def _update_for_track(track_id, time_idx=None):
        track_key_value = _normalize_track_id_key(track_id)
        track_df = state_lookup.get((sample_key, track_key_value), None)
        if track_df is None or len(track_df) == 0:
            _show_message(
                f"No behavioral-state sequence found for sample={sample_key}, TrackID={track_key_value}."
            )
            return

        times = pd.to_numeric(track_df[str(time_col)], errors="coerce")
        cursor_time = None
        if time_idx is not None:
            cursor_time = float(time_idx)
            if times.notna().any():
                cursor_time = float(times.iloc[(times - cursor_time).abs().argmin()])

        cluster_label = _resolve_track_cluster_label(
            cluster_lookup,
            sample_key,
            track_key_value,
            cursor_time=cursor_time,
        )
        title_text = f"TrackID {track_key_value} | {cluster_col}: {cluster_label}"
        rgb_img = render_track_statebar_image(
            track_df=track_df,
            state_color_map=state_color_map,
            state_col=state_col,
            time_col=time_col,
            cursor_time=cursor_time,
            title=title_text,
        )
        info_label.setText(
            f"Sample: {sample_key}<br>"
            f"TrackID: {track_key_value}<br>"
            f"{cluster_col}: {cluster_label}<br>"
            f"Clicked {time_col}: {cursor_time if cursor_time is not None else 'unknown'}"
        )
        _set_label_pixmap_from_rgb(image_label, rgb_img)

    def _on_click(layer, event):
        track_value, index = _layer_value_at_world_position(track_layer, event.position)
        if track_value is None or str(track_value) in {"0", "0.0", ""}:
            _show_message("No tracked cell under click.")
            return
        time_idx = index[0] if len(index) > 0 else None
        _update_for_track(track_value, time_idx=time_idx)

    for clickable_layer in clickable_layers:
        clickable_layer.mouse_drag_callbacks.append(_on_click)
    return widget


def resolve_behav3d_output_dir(
    behav3d_folder: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Resolve the BEHAV3D output folder used by the backprojection helpers."""
    candidate = output_dir if output_dir not in {None, ""} else behav3d_folder
    if candidate in {None, ""}:
        raise ValueError(
            "Set INTERACTIVE_BEHAV3D_FOLDER to your BEHAV3D results folder "
            "or provide INTERACTIVE_OUTPUT_DIR explicitly."
        )

    resolved = Path(candidate).expanduser()
    if not resolved.is_absolute():
        resolved = resolved.resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"BEHAV3D folder does not exist: '{resolved}'")

    expected_images = resolved / "images"
    expected_analysis = resolved / "analysis"
    if not expected_images.exists() or not expected_analysis.exists():
        raise FileNotFoundError(
            "Expected a BEHAV3D output folder containing 'images/' and 'analysis/', "
            f"but got '{resolved}'."
        )
    return resolved


def _resolve_optional_path(path_value: str | Path | None, base_dir: Path) -> Path | None:
    if path_value in {None, ""}:
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    return path


def _load_behav3d_parameters(output_dir: Path) -> dict:
    params_path = output_dir / "behav3d_parameters.yml"
    if not params_path.exists():
        return {}
    try:
        loaded = yaml.safe_load(params_path.read_text(encoding="utf-8")) or {}
    except Exception as exc:
        raise RuntimeError(f"Could not read '{params_path}': {exc}") from exc
    if not isinstance(loaded, dict):
        return {}
    return loaded


def _resolve_metadata_csv_path(output_dir: Path, params: dict) -> Path | None:
    candidate = (
        (params.get("paths") or {}).get("metadata_csv")
        if isinstance(params, dict)
        else None
    )
    if candidate not in {None, ""}:
        path = Path(candidate).expanduser()
        if not path.is_absolute():
            path = output_dir / path
        if path.exists():
            return path
    fallback = output_dir / "metadata.csv"
    if fallback.exists():
        return fallback
    return None


def _load_metadata_table(metadata_csv_path: Path | None) -> pd.DataFrame | None:
    if metadata_csv_path is None:
        return None
    try:
        return pd.read_csv(metadata_csv_path, low_memory=False)
    except Exception as exc:
        raise RuntimeError(f"Could not read metadata CSV '{metadata_csv_path}': {exc}") from exc


def _detect_available_sample_names(output_dir: Path, metadata: pd.DataFrame | None) -> list[str]:
    names: list[str] = []
    if metadata is not None and "sample_name" in metadata.columns:
        names.extend(
            metadata["sample_name"]
            .astype("string")
            .dropna()
            .str.strip()
            .tolist()
        )

    images_dir = output_dir / "images"
    reserved_non_sample_dirs = {
        "PixelClassification",
        "pixel_classifier",
        "pixelclassification",
        "SignalUnmixing",
    }
    if images_dir.exists():
        for path in images_dir.iterdir():
            if path.is_dir() and str(path.name) not in reserved_non_sample_dirs:
                names.append(str(path.name).strip())
    return sorted({str(x).strip() for x in names if str(x).strip() != ""})


def _extract_sample_names_from_obs(obs: pd.DataFrame | None, sample_col: str = "sample_name") -> list[str]:
    if obs is None or str(sample_col) not in getattr(obs, "columns", []):
        return []
    values = (
        obs[str(sample_col)]
        .astype("string")
        .dropna()
        .str.strip()
    )
    values = values[values != ""]
    return sorted(values.unique().tolist())


def _extract_sample_names_from_adata(adata, sample_col: str = "sample_name") -> list[str]:
    if adata is None or not hasattr(adata, "obs"):
        return []
    return _extract_sample_names_from_obs(adata.obs, sample_col=sample_col)


def _read_h5ad_sample_names(h5ad_path: Path, sample_col: str = "sample_name") -> list[str]:
    adata = sc.read_h5ad(str(h5ad_path), backed="r")
    try:
        return _extract_sample_names_from_adata(adata, sample_col=sample_col)
    finally:
        file_obj = getattr(adata, "file", None)
        if file_obj is not None:
            try:
                file_obj.close()
            except Exception:
                pass


def _intersect_sample_lists(*sample_lists: list[str]) -> list[str]:
    cleaned = []
    for sample_list in sample_lists:
        values = [str(v).strip() for v in (sample_list or []) if str(v).strip() != ""]
        if len(values) > 0:
            cleaned.append(set(values))
    if len(cleaned) == 0:
        return []
    shared = set.intersection(*cleaned)
    return sorted(shared)


def _detect_available_cell_types(output_dir: Path) -> list[str]:
    analysis_dir = output_dir / "analysis"
    if not analysis_dir.exists():
        return []

    detected = []
    for path in sorted(analysis_dir.iterdir()):
        if not path.is_dir():
            continue
        name = str(path.name).strip()
        if (name == "") or (" - Copy" in name):
            continue
        state_path = path / "behavioral_states" / f"BEHAV3D_{name}_behavioral_states.h5ad"
        traj_dir = path / "behavorial_trajectories"
        if not state_path.exists():
            continue
        if not traj_dir.exists():
            continue
        if len(list(traj_dir.glob("*.h5ad"))) == 0:
            continue
        detected.append(name)
    return sorted(detected)


def _resolve_sample_dim_order(metadata: pd.DataFrame | None, sample_name: str) -> str:
    if metadata is None or "sample_name" not in metadata.columns:
        return "TCZYX"
    rows = metadata[metadata["sample_name"].astype("string").str.strip() == str(sample_name).strip()]
    if len(rows) == 0:
        return "TCZYX"
    value = str(rows.iloc[0].get("dimension_order", "TCZYX")).strip()
    return value or "TCZYX"


def _resolve_track_adata_path(output_dir: Path, cell_type: str, track_adata_path: str | Path | None = None) -> Path:
    explicit_path = _resolve_optional_path(track_adata_path, base_dir=output_dir)
    if explicit_path is not None:
        if explicit_path.exists():
            return explicit_path
        raise FileNotFoundError(f"Track-classification h5ad not found: '{explicit_path}'")

    trajectories_dir = output_dir / "analysis" / str(cell_type) / "behavorial_trajectories"
    canonical = trajectories_dir / get_dtaidistance_track_trajectories_filename(str(cell_type))
    if canonical.exists():
        return canonical
    if trajectories_dir.exists():
        matches = sorted(trajectories_dir.glob("*.h5ad"))
        if len(matches) > 0:
            return matches[0]
    raise FileNotFoundError(
        "Could not find track-classification h5ad. Expected "
        f"'{canonical}' or another '*.h5ad' in '{trajectories_dir}'."
    )


def _resolve_state_adata_path(output_dir: Path, cell_type: str, state_adata_path: str | Path | None = None) -> Path:
    explicit_path = _resolve_optional_path(state_adata_path, base_dir=output_dir)
    if explicit_path is not None:
        if explicit_path.exists():
            return explicit_path
        raise FileNotFoundError(f"Behavioral-state h5ad not found: '{explicit_path}'")

    path = (
        output_dir
        / "analysis"
        / str(cell_type)
        / "behavioral_states"
        / f"BEHAV3D_{cell_type}_behavioral_states.h5ad"
    )
    if not path.exists():
        raise FileNotFoundError(f"Behavioral-state h5ad not found: '{path}'")
    return path


def _resolve_sample_image_paths(
    *,
    output_dir: Path,
    sample_name: str,
    cell_type: str,
    raw_image_path: str | Path | None = None,
    tracked_img_path: str | Path | None = None,
    metadata_csv_path: Path | None = None,
    verbose: bool = True,
) -> tuple[Path, Path]:
    raw_path_explicit = _resolve_optional_path(raw_image_path, base_dir=output_dir)
    raw_path = raw_path_explicit if raw_path_explicit is not None else _resolve_raw_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        verbose=bool(verbose),
        metadata_csv_path=None if metadata_csv_path is None else str(metadata_csv_path),
    )
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(f"Could not find raw image for sample '{sample_name}'.")

    tracked_path_explicit = _resolve_optional_path(tracked_img_path, base_dir=output_dir)
    tracked_path = tracked_path_explicit if tracked_path_explicit is not None else _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=bool(verbose),
        metadata_csv_path=None if metadata_csv_path is None else str(metadata_csv_path),
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            f"Could not find tracked image for sample '{sample_name}' and cell_type '{cell_type}'."
        )

    return Path(raw_path), Path(tracked_path)


def _detect_supported_samples_for_cell_type(
    output_dir: Path,
    cell_type: str,
    metadata: pd.DataFrame | None = None,
    track_adata_path: str | Path | None = None,
    state_adata_path: str | Path | None = None,
) -> list[str]:
    track_path = _resolve_track_adata_path(
        output_dir=output_dir,
        cell_type=cell_type,
        track_adata_path=track_adata_path,
    )
    state_path = _resolve_state_adata_path(
        output_dir=output_dir,
        cell_type=cell_type,
        state_adata_path=state_adata_path,
    )
    detected_samples = _detect_available_sample_names(output_dir=output_dir, metadata=metadata)
    track_samples = _read_h5ad_sample_names(track_path, sample_col="sample_name")
    state_samples = _read_h5ad_sample_names(state_path, sample_col="sample_name")
    supported = _intersect_sample_lists(detected_samples, track_samples, state_samples)
    if len(supported) > 0:
        return supported
    return _intersect_sample_lists(track_samples, state_samples)


def _ensure_trajectory_pixel_coordinates(adata_full, output_dir: Path, cell_type: str):
    """Populate pixel coordinates for trajectories, allowing either 3D or 2D pixel positions."""
    obs_cols = set(adata_full.obs.columns)
    pixel_triplet_3d = {"pixel_position_z", "pixel_position_y", "pixel_position_x"}
    pixel_pair_2d = {"pixel_position_y", "pixel_position_x"}
    if pixel_triplet_3d.issubset(obs_cols) or pixel_pair_2d.issubset(obs_cols):
        return {
            "enriched": False,
            "csv_path": None,
            "added_columns": [],
            "coordinate_mode": "3d" if pixel_triplet_3d.issubset(obs_cols) else "2d",
        }

    positions_csv_path = resolve_exemplar_positions_csv_path(output_dir=output_dir, cell_type=cell_type)
    merge_info = merge_coordinate_columns_into_obs(
        adata=adata_full,
        positions_csv_path=positions_csv_path,
        sample_col="sample_name",
        track_col="TrackID",
        time_col="position_t",
    )
    obs_cols = set(adata_full.obs.columns)
    if pixel_triplet_3d.issubset(obs_cols):
        coordinate_mode = "3d"
    elif pixel_pair_2d.issubset(obs_cols):
        coordinate_mode = "2d"
    else:
        raise ValueError(
            "Trajectory overlays require pixel coordinates in adata_full.obs. "
            "Expected either ['pixel_position_z', 'pixel_position_y', 'pixel_position_x'] "
            "or ['pixel_position_y', 'pixel_position_x']. "
            f"Enrichment source: '{positions_csv_path}'."
        )
    return {
        "enriched": True,
        "csv_path": str(positions_csv_path),
        "added_columns": list(merge_info.get("added_columns", [])),
        "coordinate_mode": coordinate_mode,
    }


def _slice_track_frame(tracked_img_view, time_index: int):
    if getattr(tracked_img_view, "ndim", 0) == 4:
        max_index = max(int(tracked_img_view.shape[0]) - 1, 0)
        safe_index = max(0, min(int(time_index), max_index))
        return tracked_img_view[safe_index]
    return tracked_img_view


def _build_track_cluster_overlay_frame(payload, time_index: int):
    labels_frame = _slice_track_frame(payload["tracked_img_view"], time_index)
    return backproject_track_cluster_at_timepoint(
        labels_frame=labels_frame,
        cluster_code_lookup=payload["cluster_code_lookup"],
        time_index=time_index,
        track_col="TrackID",
        time_col="position_t",
        background_value=0,
    )


def build_track_cluster_backprojection_payload(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    cluster_col: str = "ClusterID",
    state_col: str = "behavioral_state",
    output_col: str = "track_behavioral_cluster",
    raw_image_path: str | Path | None = None,
    tracked_img_path: str | Path | None = None,
    track_bp_img_path: str | Path | None = None,
    track_adata_path: str | Path | None = None,
    state_adata_path: str | Path | None = None,
    show_trajectories: bool = True,
    verbose: bool = True,
):
    sample_name = "" if sample_name is None else str(sample_name).strip()
    sample_name_was_provided = len(sample_name) > 0
    sample_name_was_auto_selected = False
    cell_type = "" if cell_type is None else str(cell_type).strip()
    cluster_col = str(cluster_col).strip() or "ClusterID"
    state_col = str(state_col).strip() or "behavioral_state"
    output_col = str(output_col).strip() or "track_behavioral_cluster"
    output_dir = Path(output_dir).expanduser()
    _ = track_bp_img_path

    params = _load_behav3d_parameters(output_dir)
    metadata_csv_path = _resolve_metadata_csv_path(output_dir, params)
    metadata = _load_metadata_table(metadata_csv_path)
    detected_samples = _detect_available_sample_names(output_dir=output_dir, metadata=metadata)
    if len(sample_name) == 0:
        if len(detected_samples) == 0:
            raise ValueError(
                "sample_name was not provided and no samples could be detected from metadata/images."
            )
        sample_name = str(detected_samples[0])
        sample_name_was_auto_selected = True
        if bool(verbose):
            print(f"Auto-selected sample_name='{sample_name}'")
    if len(cell_type) == 0:
        available_cell_types = _detect_available_cell_types(output_dir=output_dir)
        if len(available_cell_types) == 0:
            raise ValueError(
                "cell_type was not provided and no trajectory/state cell types were detected under analysis/."
            )
        cell_type = str(available_cell_types[0])
        if bool(verbose):
            print(f"Auto-selected cell_type='{cell_type}'")
    saved_channels = (
        (params.get("viewer_display") or {}).get("channels", {})
        if isinstance(params, dict)
        else {}
    )
    dim_order = _resolve_sample_dim_order(metadata, sample_name=sample_name)
    raw_path, tracked_path = _resolve_sample_image_paths(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        metadata_csv_path=metadata_csv_path,
        verbose=bool(verbose),
    )

    track_adata_resolved = _resolve_track_adata_path(
        output_dir=output_dir,
        cell_type=cell_type,
        track_adata_path=track_adata_path,
    )
    state_adata_resolved = _resolve_state_adata_path(
        output_dir=output_dir,
        cell_type=cell_type,
        state_adata_path=state_adata_path,
    )

    adata_tracks = sc.read_h5ad(str(track_adata_resolved))
    adata_full = sc.read_h5ad(str(state_adata_resolved))
    track_samples = _extract_sample_names_from_adata(adata_tracks, sample_col="sample_name")
    state_samples = _extract_sample_names_from_adata(adata_full, sample_col="sample_name")
    supported_samples = _intersect_sample_lists(detected_samples, track_samples, state_samples)
    if len(supported_samples) == 0:
        supported_samples = _intersect_sample_lists(track_samples, state_samples)

    if sample_name not in supported_samples:
        if sample_name_was_auto_selected and len(supported_samples) > 0:
            auto_selected_before = sample_name
            sample_name = str(supported_samples[0])
            dim_order = _resolve_sample_dim_order(metadata, sample_name=sample_name)
            raw_path, tracked_path = _resolve_sample_image_paths(
                output_dir=output_dir,
                sample_name=sample_name,
                cell_type=cell_type,
                raw_image_path=raw_image_path,
                tracked_img_path=tracked_img_path,
                metadata_csv_path=metadata_csv_path,
                verbose=bool(verbose),
            )
            if bool(verbose):
                print(
                    "Auto-switched sample because the initial detected sample has no "
                    f"'{cell_type}' track classifications: '{auto_selected_before}' -> '{sample_name}'."
                )
        else:
            availability_suffix = ""
            if len(supported_samples) > 0:
                availability_suffix = f" Available classified samples: {supported_samples}."
            raise ValueError(
                f"Sample '{sample_name}' is not present in the '{cell_type}' track-classification h5ads. "
                "Raw/tracked images may exist for this sample, but no track-cluster assignments are available "
                f"for cell_type '{cell_type}' in this run.{availability_suffix}"
            )

    backproj_obs, missing_windows = build_track_cluster_frame_assignments(
        adata_full=adata_full,
        adata_tracks=adata_tracks,
        cluster_col=cluster_col,
        output_col=output_col,
        sample_name=sample_name,
        verbose=bool(verbose),
    )
    state_order = _get_classification_state_order(adata_tracks, cluster_col)
    code_map = _build_code_map(backproj_obs, state_col=output_col, state_order=state_order)
    if len(code_map) == 0:
        raise ValueError(
            f"'{cluster_col}' has no non-empty labels for sample '{sample_name}'."
        )
    state_colors = _normalize_label_color_map(
        code_map.keys(),
        colors=_get_classification_state_colors(adata_tracks, cluster_col),
    )
    code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    label_map = {str(int(code)): str(label) for label, code in code_map.items()}
    label_to_code = {str(label): int(code) for label, code in code_map.items()}
    mapping_text = _build_state_mapping_text(label_map, code_colors=code_colors)
    cluster_code_lookup = prepare_track_cluster_code_lookup(
        backproj_obs=backproj_obs,
        code_map=code_map,
        output_col=output_col,
    )

    raw_img = load_image(raw_path)
    tracked_img = load_image(tracked_path)

    tracked_img_view = _align_labels_to_raw_shape_for_view_lazy(
        tracked_img,
        raw_img,
        layer_name="TrackID",
        verbose=bool(verbose),
    )

    trajectory_data = {}
    coordinate_info = None
    if bool(show_trajectories):
        coordinate_info = _ensure_trajectory_pixel_coordinates(
            adata_full=adata_full,
            output_dir=output_dir,
            cell_type=cell_type,
        )
        trajectory_data = prepare_track_cluster_trajectory_data(
            adata_full=adata_full,
            backproj_obs=backproj_obs,
            output_col=output_col,
            sample_name=sample_name,
        )

    return {
        "sample_name": sample_name,
        "output_dir": output_dir,
        "cell_type": cell_type,
        "cluster_col": cluster_col,
        "state_col": state_col,
        "output_col": output_col,
        "raw_path": raw_path,
        "tracked_path": tracked_path,
        "track_adata_path": track_adata_resolved,
        "state_adata_path": state_adata_resolved,
        "metadata_csv_path": metadata_csv_path,
        "dim_order": dim_order,
        "saved_channels": saved_channels,
        "track_samples": track_samples,
        "state_samples": state_samples,
        "supported_samples": supported_samples,
        "raw_img": raw_img,
        "raw_shape": tuple(int(v) for v in raw_img.shape),
        "tracked_img_view": tracked_img_view,
        "cluster_code_lookup": cluster_code_lookup,
        "code_map": code_map,
        "label_map": label_map,
        "label_to_code": label_to_code,
        "code_colors": code_colors,
        "state_colors": state_colors,
        "mapping_text": mapping_text,
        "adata_full": adata_full,
        "adata_tracks": adata_tracks,
        "trajectory_data": trajectory_data,
        "show_trajectories": bool(show_trajectories),
        "missing_windows": int(missing_windows),
        "coordinate_info": coordinate_info,
        "track_layer_name": "TrackID",
        "label_layer_name": output_col,
        "track_layer_opacity": float(INTERACTIVE_TRACK_LAYER_OPACITY),
        "class_layer_opacity": float(INTERACTIVE_CLASS_LAYER_OPACITY),
    }


def launch_track_cluster_backprojection_viewer(payload, run: bool = True):
    viewer = napari.Viewer()
    _bp_add_raw_channels(
        viewer,
        raw_img=payload["raw_img"],
        sample_name=payload["sample_name"],
        saved_channels=payload["saved_channels"],
        dim_order=payload["dim_order"],
    )
    track_layer = viewer.add_labels(
        payload["tracked_img_view"],
        name=payload["track_layer_name"],
        visible=True,
        opacity=payload.get("track_layer_opacity", 0.35),
    )
    initial_time_index = _current_viewer_frame(viewer)
    initial_state_frame, _ = _build_track_cluster_overlay_frame(payload, time_index=initial_time_index)
    state_layer = viewer.add_labels(
        initial_state_frame,
        name=payload["label_layer_name"],
        visible=True,
        opacity=payload.get("class_layer_opacity", 0.85),
    )
    _apply_state_code_colors_to_layer(state_layer, payload["code_colors"])

    refresh_state = {"time_index": None}

    def _refresh_state_overlay(event=None, *, force=False):
        time_index = _current_viewer_frame(viewer)
        if (not force) and (refresh_state["time_index"] == int(time_index)):
            return
        mapped_frame, ids_with_value = _build_track_cluster_overlay_frame(payload, time_index=time_index)
        state_layer.data = mapped_frame
        state_layer.metadata["ids_with_value"] = np.asarray(ids_with_value).tolist()
        state_layer.metadata["time_index"] = int(time_index)
        refresh_state["time_index"] = int(time_index)

    viewer.dims.events.current_step.connect(_refresh_state_overlay)
    state_layer.metadata["behav3d_refresh_callback"] = _refresh_state_overlay
    _refresh_state_overlay(force=True)

    trajectory_layers = []
    if bool(payload.get("show_trajectories")):
        trajectory_layers = add_track_cluster_trajectory_layers(
            viewer,
            trajectory_data=payload.get("trajectory_data", {}),
            code_colors=payload["code_colors"],
            label_map=payload["label_map"],
            output_col=payload["output_col"],
            visible=True,
        )
        for trajectory_layer in trajectory_layers:
            trajectory_layer.blending = "opaque"

    _add_track_viewer_tools_dock(
        viewer=viewer,
        output_dir=payload["output_dir"],
        sample_name=payload["sample_name"],
        cell_type=payload["cell_type"],
        output_col=payload["output_col"],
        trajectory_layers=trajectory_layers,
        title="Trajectory + Screenshot",
    )

    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=payload["mapping_text"],
        label_map=payload["label_map"],
        code_colors=payload["code_colors"],
        title="Track Cluster Mapping",
    )
    if (not added_dock) and payload.get("mapping_text"):
        print(payload["mapping_text"])

    _add_track_statebar_click_dock_compat(
        viewer=viewer,
        sample_name=payload["sample_name"],
        adata_full=payload["adata_full"],
        adata_tracks=payload["adata_tracks"],
        track_layer_name=payload["track_layer_name"],
        clickable_layer_name=payload["label_layer_name"],
        state_col=payload["state_col"],
        cluster_col=payload["cluster_col"],
        title="Track State Bar",
    )

    coordinate_info = payload.get("coordinate_info")
    if payload.get("sample_name"):
        summary = (
            "Opened track-classification split-channel backprojection viewer for sample "
            f"'{payload['sample_name']}' with raw='{payload['raw_path'].name}' "
            f"shape={payload['raw_shape']}, tracked='{payload['tracked_path'].name}' "
            f"shape={tuple(int(v) for v in payload['tracked_img_view'].shape)}."
        )
        print(summary)
        print(
            "Track class overlay: "
            f"on-the-fly current-frame backprojection | missing_windows={payload['missing_windows']} | "
            f"track_opacity={track_layer.opacity} | class_opacity={state_layer.opacity}"
        )
        if coordinate_info is not None:
            print(
                "Trajectory coordinates: "
                f"mode={coordinate_info.get('coordinate_mode')} | "
                f"enriched={coordinate_info.get('enriched')} | "
                f"added_columns={coordinate_info.get('added_columns')}"
            )

    if bool(run):
        napari.run()
    return viewer


def show_track_cluster_backprojection_split_channels(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    cluster_col: str = "ClusterID",
    state_col: str = "behavioral_state",
    output_col: str = "track_behavioral_cluster",
    raw_image_path: str | Path | None = None,
    tracked_img_path: str | Path | None = None,
    track_bp_img_path: str | Path | None = None,
    track_adata_path: str | Path | None = None,
    state_adata_path: str | Path | None = None,
    show_trajectories: bool = True,
    run: bool = True,
    verbose: bool = True,
):
    payload = build_track_cluster_backprojection_payload(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        cluster_col=cluster_col,
        state_col=state_col,
        output_col=output_col,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        track_bp_img_path=track_bp_img_path,
        track_adata_path=track_adata_path,
        state_adata_path=state_adata_path,
        show_trajectories=show_trajectories,
        verbose=verbose,
    )
    return launch_track_cluster_backprojection_viewer(payload, run=run)


def build_track_cluster_backprojection_ui(
    behav3d_folder: str | Path | None = None,
    output_dir: str | Path | None = None,
    *,
    default_cell_type: str | None = None,
    default_sample_name: str | None = None,
    cluster_col: str = "ClusterID",
    state_col: str = "behavioral_state",
    output_col: str = "track_behavioral_cluster",
    show_trajectories: bool = True,
    verbose: bool = True,
):
    """Create a small notebook UI with a sample dropdown and open button."""
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except Exception as exc:
        raise RuntimeError(
            "Notebook UI requires 'ipywidgets' and IPython display support."
        ) from exc

    resolved_output_dir = resolve_behav3d_output_dir(
        behav3d_folder=behav3d_folder,
        output_dir=output_dir,
    )
    params = _load_behav3d_parameters(resolved_output_dir)
    metadata_csv_path = _resolve_metadata_csv_path(resolved_output_dir, params)
    metadata = _load_metadata_table(metadata_csv_path)
    all_sample_options = _detect_available_sample_names(resolved_output_dir, metadata)
    cell_type_options = _detect_available_cell_types(resolved_output_dir)

    if len(all_sample_options) == 0:
        raise ValueError(f"No samples were detected under '{resolved_output_dir}'.")
    if len(cell_type_options) == 0:
        raise ValueError(f"No compatible cell types were detected under '{resolved_output_dir / 'analysis'}'.")

    default_cell_type = (
        str(default_cell_type).strip()
        if default_cell_type not in {None, ""}
        else cell_type_options[0]
    )
    if default_cell_type not in cell_type_options:
        default_cell_type = cell_type_options[0]

    supported_samples_cache: dict[str, list[str]] = {}

    def _samples_for_cell_type(selected_cell_type: str) -> list[str]:
        key = str(selected_cell_type).strip()
        if key not in supported_samples_cache:
            supported_samples_cache[key] = _detect_supported_samples_for_cell_type(
                output_dir=resolved_output_dir,
                cell_type=key,
                metadata=metadata,
            )
        return list(supported_samples_cache[key])

    initial_sample_options = _samples_for_cell_type(default_cell_type)
    if len(initial_sample_options) == 0:
        initial_sample_options = list(all_sample_options)

    default_sample_name = (
        str(default_sample_name).strip()
        if default_sample_name not in {None, ""}
        else initial_sample_options[0]
    )
    if default_sample_name not in initial_sample_options:
        default_sample_name = initial_sample_options[0]

    title = widgets.HTML("<b>Track Classification Backprojection</b>")
    outdir_label = widgets.HTML(
        f"<span style='color:#666;'>Output dir:</span> "
        f"<code>{resolved_output_dir}</code>"
    )
    sample_dropdown = widgets.Dropdown(
        options=initial_sample_options,
        value=default_sample_name,
        description="Sample",
        layout=widgets.Layout(width="700px"),
        style={"description_width": "90px"},
    )
    cell_type_dropdown = widgets.Dropdown(
        options=cell_type_options,
        value=default_cell_type,
        description="Cell type",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "90px"},
    )
    cluster_col_text = widgets.Text(
        value=str(cluster_col),
        description="Cluster col",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "90px"},
    )
    state_col_text = widgets.Text(
        value=str(state_col),
        description="State col",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "90px"},
    )
    output_col_text = widgets.Text(
        value=str(output_col),
        description="Output col",
        layout=widgets.Layout(width="320px"),
        style={"description_width": "90px"},
    )
    show_traj_checkbox = widgets.Checkbox(
        value=bool(show_trajectories),
        description="Show trajectory layers",
        indent=False,
    )
    verbose_checkbox = widgets.Checkbox(
        value=bool(verbose),
        description="Verbose logging",
        indent=False,
    )
    open_button = widgets.Button(
        description="Open napari backprojection",
        button_style="success",
        tooltip="Open the track backprojection viewer with live class overlay.",
    )
    sample_hint_html = widgets.HTML("")
    status_html = widgets.HTML(
        "<span style='color:#666;'>Select a sample and click the button.</span>"
    )
    output_widget = widgets.Output()

    def _refresh_sample_dropdown(*_args):
        selected_cell_type = str(cell_type_dropdown.value).strip()
        supported_samples = _samples_for_cell_type(selected_cell_type)
        dropdown_options = supported_samples if len(supported_samples) > 0 else list(all_sample_options)
        current_value = sample_dropdown.value
        sample_dropdown.options = dropdown_options
        if current_value not in dropdown_options:
            sample_dropdown.value = dropdown_options[0]
        if len(supported_samples) > 0:
            sample_hint_html.value = (
                f"<span style='color:#666;'>Classified samples for <code>{selected_cell_type}</code>: "
                f"{len(supported_samples)} available.</span>"
            )
        else:
            sample_hint_html.value = (
                f"<span style='color:#a66;'>No classified samples were detected for "
                f"<code>{selected_cell_type}</code>. Showing all image samples instead.</span>"
            )

    cell_type_dropdown.observe(_refresh_sample_dropdown, names="value")
    _refresh_sample_dropdown()

    def _on_open_clicked(_):
        open_button.disabled = True
        status_html.value = (
            "<span style='color:#666;'>Preparing backprojection and opening napari...</span>"
        )
        output_widget.clear_output()
        with output_widget:
            try:
                viewer = show_track_cluster_backprojection_split_channels(
                    sample_name=sample_dropdown.value,
                    output_dir=resolved_output_dir,
                    cell_type=cell_type_dropdown.value,
                    cluster_col=cluster_col_text.value,
                    state_col=state_col_text.value,
                    output_col=output_col_text.value,
                    show_trajectories=bool(show_traj_checkbox.value),
                    run=True,
                    verbose=bool(verbose_checkbox.value),
                )
                status_html.value = (
                    "<span style='color:#080;'>Viewer opened successfully.</span>"
                )
                return viewer
            except Exception as exc:
                import traceback

                traceback.print_exc()
                status_html.value = (
                    f"<span style='color:#a33;'>Failed: {type(exc).__name__}: {exc}</span>"
                )
            finally:
                open_button.disabled = False

    open_button.on_click(_on_open_clicked)

    ui = widgets.VBox(
        [
            title,
            outdir_label,
            sample_dropdown,
            sample_hint_html,
            widgets.HBox([cell_type_dropdown, cluster_col_text]),
            widgets.HBox([state_col_text, output_col_text]),
            widgets.HBox([show_traj_checkbox, verbose_checkbox]),
            open_button,
            status_html,
            output_widget,
        ]
    )
    display(ui)
    return ui


if __name__ == "__main__":
    # %%
    # Interactive notebook usage:
    # 1. Fill the INTERACTIVE_* values above.
    # 2. Run this file/cell to create a small UI with a sample chooser.
    # 3. Click "Open napari backprojection" to export + open the viewer.
    interactive_output_dir = resolve_behav3d_output_dir(
        behav3d_folder=INTERACTIVE_BEHAV3D_FOLDER,
        output_dir=INTERACTIVE_OUTPUT_DIR,
    )
    backprojection_ui = build_track_cluster_backprojection_ui(
        output_dir=interactive_output_dir,
        default_cell_type=INTERACTIVE_CELL_TYPE,
        default_sample_name=INTERACTIVE_SAMPLE_NAME,
        cluster_col=INTERACTIVE_CLUSTER_COL,
        state_col=INTERACTIVE_STATE_COL,
        output_col=INTERACTIVE_OUTPUT_COL,
        show_trajectories=INTERACTIVE_SHOW_TRAJECTORIES,
        verbose=INTERACTIVE_VERBOSE,
    )
