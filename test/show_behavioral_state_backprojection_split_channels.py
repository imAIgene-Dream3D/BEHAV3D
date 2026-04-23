from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import napari

from behav3d.analysis.clustering.state.visualization.backprojection import (
    _add_mapping_dock_widget,
    _align_labels_to_raw_shape_for_view,
    _build_state_mapping_text,
    _ensure_behavioral_state_backprojection_for_sample,
    _extract_state_label_map,
    _resolve_behavioral_state_image_path,
    _resolve_raw_image_path,
    _resolve_tracked_image_path,
)
from behav3d.io.images import load_image


# %%
# Interactive configuration
# Fill these in when running the script line by line from the IDE.
# Point this to the BEHAV3D results folder that contains `images/`, `analysis/`,
# and typically `metadata.csv`.
INTERACTIVE_BEHAV3D_FOLDER = r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity\behav3d"
# Optional explicit override if your BEHAV3D output directory differs.
INTERACTIVE_SAMPLE_NAME = "BHVD_SB1_Exp009_Img001"
INTERACTIVE_OUTPUT_DIR = None
INTERACTIVE_CELL_TYPE = "tcell"
INTERACTIVE_STATE_COL = "full_behavioral_cluster"
INTERACTIVE_RAW_IMAGE_PATH = None
INTERACTIVE_TRACKED_IMG_PATH = None
INTERACTIVE_STATE_IMG_PATH = None
INTERACTIVE_AUTO_CREATE_IF_MISSING = True
INTERACTIVE_REFRESH_IF_STALE = True
INTERACTIVE_VERBOSE = True


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


def _resolve_inputs(
    sample_name: str,
    output_dir: Path,
    cell_type: str,
    state_col: str,
    raw_image_path: str | None,
    tracked_img_path: str | None,
    state_img_path: str | None,
    auto_create_if_missing: bool,
    refresh_if_stale: bool,
    verbose: bool,
) -> tuple[Path, Path, Path]:
    raw_path = Path(raw_image_path) if raw_image_path is not None else _resolve_raw_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        verbose=verbose,
    )
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(
            "Could not find raw image for sample "
            f"'{sample_name}'. Expected '{Path(output_dir, 'images', sample_name, f'{sample_name}.zarr')}' "
            "or '.zarr.zip'."
        )

    tracked_path = Path(tracked_img_path) if tracked_img_path is not None else _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            "Could not find tracked image for sample "
            f"'{sample_name}' and cell_type '{cell_type}'."
        )

    if state_img_path is not None:
        state_path = Path(state_img_path)
        if not state_path.exists():
            if auto_create_if_missing:
                state_path = _ensure_behavioral_state_backprojection_for_sample(
                    sample_name=sample_name,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    state_col=state_col,
                    track_col="TrackID",
                    time_col="position_t",
                    sample_col="sample_name",
                    background_value=0,
                    enforce_time_coverage=True,
                    refresh_if_stale=refresh_if_stale,
                    verbose=verbose,
                )
            else:
                raise FileNotFoundError(
                    "Could not find behavioral-state image for sample "
                    f"'{sample_name}' and cell_type '{cell_type}'."
                )
    elif auto_create_if_missing:
        state_path = _ensure_behavioral_state_backprojection_for_sample(
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            state_col=state_col,
            track_col="TrackID",
            time_col="position_t",
            sample_col="sample_name",
            background_value=0,
            enforce_time_coverage=True,
            refresh_if_stale=refresh_if_stale,
            verbose=verbose,
        )
    else:
        state_path = _resolve_behavioral_state_image_path(
            output_dir=output_dir,
            sample_name=sample_name,
            cell_type=cell_type,
            verbose=verbose,
        )
        if state_path is None or not Path(state_path).exists():
            raise FileNotFoundError(
                "Could not find behavioral-state image for sample "
                f"'{sample_name}' and cell_type '{cell_type}'."
            )

    return Path(raw_path), Path(tracked_path), Path(state_path)


def build_split_channel_backprojection_payload(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    state_col: str = "full_behavioral_cluster",
    raw_image_path: str | None = None,
    tracked_img_path: str | None = None,
    state_img_path: str | None = None,
    auto_create_if_missing: bool = True,
    refresh_if_stale: bool = True,
    verbose: bool = True,
):
    sample_name = str(sample_name).strip()
    output_dir = Path(output_dir)
    cell_type = str(cell_type).strip()

    if len(sample_name) == 0:
        raise ValueError("sample_name is required")
    if len(cell_type) == 0:
        raise ValueError("cell_type is required")

    raw_path, tracked_path, state_path = _resolve_inputs(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=state_col,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        state_img_path=state_img_path,
        auto_create_if_missing=bool(auto_create_if_missing),
        refresh_if_stale=bool(refresh_if_stale),
        verbose=bool(verbose),
    )

    raw_img = load_image(raw_path)
    raw_shape = tuple(int(v) for v in raw_img.shape)
    if len(raw_shape) != 5:
        raise ValueError(
            "Raw image must have shape (T, C, Z, Y, X), "
            f"but got shape {raw_shape} from '{raw_path}'."
        )
    if int(raw_shape[1]) != 3:
        raise ValueError(
            "Raw image must contain exactly 3 channels on axis 1, "
            f"but got shape {raw_shape} from '{raw_path}'."
        )

    raw_channels = [raw_img[:, ch, ...] for ch in range(3)]

    tracked_img = load_image(tracked_path)
    state_img = load_image(state_path)
    raw_reference = raw_channels[0]

    tracked_img_view = _align_labels_to_raw_shape_for_view(
        tracked_img,
        raw_reference,
        layer_name="TrackID",
        verbose=verbose,
    )
    state_img_view = _align_labels_to_raw_shape_for_view(
        state_img,
        raw_reference,
        layer_name="behavioral_state_class",
        verbose=verbose,
    )

    label_map = _extract_state_label_map(state_path)
    mapping_text = _build_state_mapping_text(label_map)

    return {
        "sample_name": sample_name,
        "output_dir": output_dir,
        "cell_type": cell_type,
        "raw_path": raw_path,
        "tracked_path": tracked_path,
        "state_path": state_path,
        "raw_shape": raw_shape,
        "raw_channels": raw_channels,
        "tracked_img_view": tracked_img_view,
        "state_img_view": state_img_view,
        "label_map": label_map,
        "mapping_text": mapping_text,
    }


def launch_split_channel_backprojection_viewer(payload, run: bool = True):
    viewer = napari.Viewer()
    for ch, channel_img in enumerate(payload["raw_channels"]):
        viewer.add_image(channel_img, name=f"raw_ch{ch}")
    viewer.add_labels(payload["tracked_img_view"], name="TrackID", visible=False)
    viewer.add_labels(
        payload["state_img_view"],
        name="behavioral_state_class",
        visible=True,
    )

    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=payload["mapping_text"],
        title="State Class Mapping",
    )
    if (not added_dock) and payload.get("mapping_text"):
        print(payload["mapping_text"])

    if payload.get("sample_name"):
        print(
            "Opened split-channel backprojection viewer for sample "
            f"'{payload['sample_name']}' with raw='{payload['raw_path'].name}' "
            f"shape={payload['raw_shape']}, "
            f"tracked='{payload['tracked_path'].name}' "
            f"shape={tuple(int(v) for v in payload['tracked_img_view'].shape)}, "
            f"states='{payload['state_path'].name}' "
            f"shape={tuple(int(v) for v in payload['state_img_view'].shape)}."
        )

    if bool(run):
        napari.run()
    return viewer


def show_behavioral_state_backprojection_split_channels(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    state_col: str = "full_behavioral_cluster",
    raw_image_path: str | None = None,
    tracked_img_path: str | None = None,
    state_img_path: str | None = None,
    auto_create_if_missing: bool = True,
    refresh_if_stale: bool = True,
    verbose: bool = True,
):
    payload = build_split_channel_backprojection_payload(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=state_col,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        state_img_path=state_img_path,
        auto_create_if_missing=auto_create_if_missing,
        refresh_if_stale=refresh_if_stale,
        verbose=verbose,
    )
    return launch_split_channel_backprojection_viewer(payload, run=True)


# %%
# Interactive usage example:
# 1. Fill the INTERACTIVE_* values above.
# 2. Run this cell to resolve the BEHAV3D folder and prepare the layers.
# 3. Inspect `payload` in the variable explorer if you want.
#
interactive_output_dir = resolve_behav3d_output_dir(
    behav3d_folder=INTERACTIVE_BEHAV3D_FOLDER,
    output_dir=INTERACTIVE_OUTPUT_DIR,
)

payload = build_split_channel_backprojection_payload(
    sample_name=INTERACTIVE_SAMPLE_NAME,
    output_dir=interactive_output_dir,
    cell_type=INTERACTIVE_CELL_TYPE,
    state_col=INTERACTIVE_STATE_COL,
    raw_image_path=INTERACTIVE_RAW_IMAGE_PATH,
    tracked_img_path=INTERACTIVE_TRACKED_IMG_PATH,
    state_img_path=INTERACTIVE_STATE_IMG_PATH,
    auto_create_if_missing=INTERACTIVE_AUTO_CREATE_IF_MISSING,
    refresh_if_stale=INTERACTIVE_REFRESH_IF_STALE,
    verbose=INTERACTIVE_VERBOSE,
)


# %%
# Launch napari from a prepared payload:
#
viewer = launch_split_channel_backprojection_viewer(payload, run=True)
