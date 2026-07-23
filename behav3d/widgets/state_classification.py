import shutil
import traceback
from pathlib import Path

import ipywidgets as widgets
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from behav3d.analysis.backprojection import backproject_columns, view_napari
from behav3d.analysis.behavior.general import relabel_cluster_ids
from behav3d.analysis.behavior.state.classification import (
    BINARY_GROUP_COL,
    FULL_STATE_COL,
    HMM_INTRINSIC_RAW_STATE_COL,
    INTRINSIC_STATE_COL,
    apply_hmm_deployment_artifact_to_full_dataset,
    load_hmm_deployment_artifact,
    save_hmm_deployment_artifact,
    save_hmm_quality_control_outputs,
    _resolve_hmm_deployment_artifact_path,
    run_hmm_state_clustering,
)
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _coerce_hex_color,
    _get_classification_state_colors,
    _get_classification_state_order,
    _mixed_label_sort_key,
    _normalize_label_color_map,
    _rebuild_full_behavioral_cluster_from_intrinsic,
    _resolve_state_paths,
    _set_classification_state_colors,
    _set_classification_state_order,
)

from behav3d.analysis.behavior.state.visualization.backprojection import (
    show_behavioral_state_backprojection,
)
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
    filter_multicolor_inputs,
)
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    save_state_composition_report,
    save_state_condition_comparison_report,
)
from behav3d.analysis.behavior.utils import _sanitize_filename_token
from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
    save_state_transition_report,
)
from behav3d.analysis.behavior.state.legacy_clustering import (
    apply_state_classifiers_to_full_dataset,
    build_identity_cluster_mapping,
    load_state_classifier_artifact,
)
from behav3d.widgets.base_state_classification import (
    BaseStateClassificationPanel,
    _winfo,
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    _add_mapping_dock_widget,
    _apply_state_code_colors_to_layer,
    _behavioral_state_backprojection_path,
    _build_code_map,
    _build_state_code_color_map,
    _build_state_mapping_text,
    _write_state_color_attrs_to_zarr,
    _resolve_raw_image_path,
    _resolve_tracked_image_path,
)
from behav3d.io.images import load_image
from behav3d.widgets.utils import build_plot_box, build_row_move_buttons, spinning_loader


_TIME_UNIT_DISPLAY_NAMES = {"s": "seconds", "m": "minutes", "h": "hours"}


def _format_timepoint_window_as_time(n_timepoints, metadata):
    """Return an HTML snippet showing ``n_timepoints`` converted to real time.

    Uses the first sample row's ``time_interval``/``time_unit`` metadata
    columns (mirroring :func:`behav3d.core.utils.hours_per_frame_from_metadata`).
    Returns an empty string when metadata is missing/unusable.
    """
    try:
        if metadata is None or len(metadata) == 0:
            return ""
        if "time_interval" not in metadata.columns or "time_unit" not in metadata.columns:
            return ""
        interval = float(metadata["time_interval"].iloc[0])
        unit = str(metadata["time_unit"].iloc[0])
        if not np.isfinite(interval) or interval <= 0 or unit not in _TIME_UNIT_DISPLAY_NAMES:
            return ""
        total = float(n_timepoints) * interval
        unit_name = _TIME_UNIT_DISPLAY_NAMES[unit]
        total_str = f"{total:g}"
        return f"<i>= {total_str} {unit_name}</i>"
    except Exception:
        return ""


def _intrinsic_hmm_backprojection_path(output_dir, sample_name, cell_type):
    return Path(
        output_dir,
        "analysis",
        str(cell_type),
        "behavioral_states",
        "backprojection",
        f"{sample_name}_{cell_type}_{INTRINSIC_STATE_COL}.zarr",
    )


def _remove_path_if_exists(path):
    path = Path(path)
    if not path.exists():
        return False
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()
    return True


def _clear_behavioral_state_backprojection_cache(output_dir, sample_name, cell_type, verbose=True):
    state_path = _behavioral_state_backprojection_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
    )
    removed = []
    for candidate in (state_path, Path(f"{state_path}.zip")):
        if _remove_path_if_exists(candidate):
            removed.append(candidate)
    if verbose and len(removed) > 0:
        print(
            "Removed existing behavioral-state backprojection cache: "
            + ", ".join(str(path) for path in removed)
        )
    return removed


def _clear_intrinsic_hmm_backprojection_cache(output_dir, sample_name, cell_type, verbose=True):
    state_path = _intrinsic_hmm_backprojection_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
    )
    removed = []
    for candidate in (state_path, Path(f"{state_path}.zip")):
        if _remove_path_if_exists(candidate):
            removed.append(candidate)
    if verbose and len(removed) > 0:
        print(
            "Removed existing intrinsic HMM backprojection cache: "
            + ", ".join(str(path) for path in removed)
        )
    return removed

def _filter_track_image_to_ids(track_img, track_ids):
    track_ids = np.asarray(list(track_ids), dtype=np.int64)
    if hasattr(track_img, "chunks"):
        import dask.array as da

        return da.where(da.isin(track_img, track_ids), track_img, 0)
    return np.where(np.isin(track_img, track_ids), track_img, 0)


def _sample_voxel_size(metadata, sample_name):
    if not isinstance(metadata, pd.DataFrame) or "sample_name" not in metadata.columns:
        return [1.0, 1.0, 1.0]
    sample_rows = metadata[metadata["sample_name"].astype(str) == str(sample_name)]
    if len(sample_rows) == 0:
        return [1.0, 1.0, 1.0]
    row = sample_rows.iloc[0]
    try:
        return [
            float(row.get("pixel_distance_z", 1.0)),
            float(row.get("pixel_distance_xy", 1.0)),
            float(row.get("pixel_distance_xy", 1.0)),
        ]
    except Exception:
        return [1.0, 1.0, 1.0]


def _report_intrinsic_hmm_backprojection_overlap(
    *,
    sample_name,
    sample_obs,
    df_tracks_full,
    verbose=True,
):
    key_cols = ["TrackID", "position_t"]
    sample_obs = sample_obs.copy()
    df_tracks_full = df_tracks_full.copy()

    hmm_rows = int(len(sample_obs))
    hmm_tracks = int(sample_obs["TrackID"].dropna().nunique()) if "TrackID" in sample_obs.columns else 0
    track_rows = int(len(df_tracks_full))
    track_lookup_full = df_tracks_full[key_cols].copy()
    track_rows_unique = (
        df_tracks_full.drop_duplicates(subset=key_cols)
        if all(col in df_tracks_full.columns for col in key_cols)
        else df_tracks_full.copy()
    )
    track_tracks = int(df_tracks_full["TrackID"].dropna().nunique()) if "TrackID" in df_tracks_full.columns else 0

    state_lookup = sample_obs[key_cols + [INTRINSIC_STATE_COL]].copy()
    state_lookup_unique = state_lookup.drop_duplicates(subset=key_cols).copy()
    track_lookup = track_rows_unique[key_cols].copy()

    overlap = track_lookup_full.merge(
        state_lookup_unique[key_cols],
        on=key_cols,
        how="inner",
        validate="many_to_one",
    )
    track_vs_hmm = track_lookup_full.merge(
        state_lookup_unique[key_cols],
        on=key_cols,
        how="left",
        indicator=True,
        validate="many_to_one",
    )
    hmm_vs_track = state_lookup_unique[key_cols].merge(
        track_lookup,
        on=key_cols,
        how="left",
        indicator=True,
        validate="one_to_one",
    )

    stats = {
        "sample_name": str(sample_name),
        "hmm_rows": hmm_rows,
        "hmm_tracks": hmm_tracks,
        "track_feature_rows": track_rows,
        "track_feature_tracks": track_tracks,
        "overlap_rows": int(len(overlap)),
        "overlap_tracks": int(overlap["TrackID"].dropna().nunique()) if "TrackID" in overlap.columns else 0,
        "hmm_rows_without_track_match": int((hmm_vs_track["_merge"] == "left_only").sum()),
        "track_rows_without_hmm_match": int((track_vs_hmm["_merge"] == "left_only").sum()),
        "preprocessing": {},
    }

    if verbose:
        print(
            "Intrinsic HMM backprojection diagnostics | "
            f"sample='{sample_name}' | "
            f"hmm_rows={stats['hmm_rows']} | hmm_tracks={stats['hmm_tracks']} | "
            f"track_feature_rows={stats['track_feature_rows']} | track_feature_tracks={stats['track_feature_tracks']} | "
            f"overlap_rows={stats['overlap_rows']} | overlap_tracks={stats['overlap_tracks']} | "
            f"hmm_rows_without_track_match={stats['hmm_rows_without_track_match']} | "
            f"track_rows_without_hmm_match={stats['track_rows_without_hmm_match']}"
        )
    return stats


def _show_intrinsic_hmm_backprojection(
    *,
    adata,
    sample_name,
    output_dir,
    cell_type,
    track_features_csv_path,
    metadata=None,
    n_workers=4,
    run=True,
    verbose=True,
):
    if adata is None or not hasattr(adata, "obs"):
        raise ValueError("model_adata with intrinsic HMM states is required for intrinsic backprojection.")
    if INTRINSIC_STATE_COL not in adata.obs.columns:
        raise ValueError(f"model_adata is missing '{INTRINSIC_STATE_COL}'.")

    sample_name = str(sample_name).strip()
    cell_type = str(cell_type).strip()
    output_dir = Path(output_dir)

    sample_obs = adata.obs[adata.obs["sample_name"].astype(str) == sample_name].copy()
    if len(sample_obs) == 0:
        raise ValueError(f"No rows found for sample '{sample_name}' in model_adata.")

    raw_path = _resolve_raw_image_path(output_dir=output_dir, sample_name=sample_name, verbose=verbose)
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(f"Could not find raw image for sample '{sample_name}'.")

    tracked_path = _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            f"Could not find tracked image for sample '{sample_name}' and cell_type '{cell_type}'."
        )

    state_path = _intrinsic_hmm_backprojection_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
    )
    code_map = _build_code_map(
        adata.obs,
        state_col=INTRINSIC_STATE_COL,
        state_order=_get_classification_state_order(adata, INTRINSIC_STATE_COL),
    )
    if track_features_csv_path is None or not Path(track_features_csv_path).exists():
        raise FileNotFoundError("Track-features CSV is required for intrinsic HMM backprojection.")

    df_tracks_full = pd.read_csv(track_features_csv_path, low_memory=False)
    df_tracks_full = df_tracks_full[df_tracks_full["sample_name"].astype(str) == sample_name].copy()
    if len(df_tracks_full) == 0:
        raise ValueError(f"No track-feature rows found for sample '{sample_name}'.")

    required_pos_cols = ["TrackID", "position_t", "position_z", "position_y", "position_x"]
    missing_pos_cols = [c for c in required_pos_cols if c not in df_tracks_full.columns]
    if missing_pos_cols:
        raise ValueError(
            "Track-features CSV is missing required position columns for backprojection: "
            f"{missing_pos_cols}"
        )

    preprocessing_meta = getattr(adata, "uns", {}).get("preprocessing", {})
    diag_stats = _report_intrinsic_hmm_backprojection_overlap(
        sample_name=sample_name,
        sample_obs=sample_obs,
        df_tracks_full=df_tracks_full,
        verbose=verbose,
    )
    if isinstance(preprocessing_meta, dict):
        diag_stats["preprocessing"] = {
            "n_timepoints_input": preprocessing_meta.get("n_timepoints_input"),
            "n_timepoints_kept": preprocessing_meta.get("n_timepoints_kept"),
            "n_timepoints_dropped_nan": preprocessing_meta.get("n_timepoints_dropped_nan"),
            "n_tracks_with_dropped_timepoints": preprocessing_meta.get("n_tracks_with_dropped_timepoints"),
            "dropped_track_preview": preprocessing_meta.get("dropped_track_preview", []),
            "dropped_track_preview_omitted": preprocessing_meta.get("dropped_track_preview_omitted", 0),
        }
        if verbose and any(v is not None for v in diag_stats["preprocessing"].values()):
            print(
                "HMM preprocessing counts | "
                f"n_timepoints_input={diag_stats['preprocessing'].get('n_timepoints_input')} | "
                f"n_timepoints_kept={diag_stats['preprocessing'].get('n_timepoints_kept')} | "
                f"n_timepoints_dropped_nan={diag_stats['preprocessing'].get('n_timepoints_dropped_nan')}"
            )
        if verbose and int(preprocessing_meta.get("n_timepoints_dropped_nan", 0) or 0) > 0:
            dropped_rows = int(preprocessing_meta.get("n_timepoints_dropped_nan", 0) or 0)
            input_rows = int(preprocessing_meta.get("n_timepoints_input", 0) or 0)
            dropped_tracks = int(preprocessing_meta.get("n_tracks_with_dropped_timepoints", 0) or 0)
            preview = list(preprocessing_meta.get("dropped_track_preview", []) or [])
            preview_text = ", ".join(str(v) for v in preview) if len(preview) > 0 else "none"
            omitted = int(preprocessing_meta.get("dropped_track_preview_omitted", 0) or 0)
            if omitted > 0:
                preview_text = f"{preview_text} (+{omitted} more)"
            print(
                f"Warning: {dropped_rows}/{input_rows} rows were removed during HMM apply before "
                "state assignment because required HMM observation features became NaN after "
                "smoothing, optional window-feature construction, or numeric coercion."
            )
            print(f"Affected tracks: {dropped_tracks}. Preview: {preview_text}")
            print(
                "This can lead to missing segments in TrackID/ClusterID volumes even when "
                "track lines still appear."
            )
    if verbose and diag_stats["overlap_rows"] < min(diag_stats["hmm_rows"], diag_stats["track_feature_rows"]):
        print(
            "Warning: track-features CSV and HMM-labeled rows do not overlap on exact "
            "(TrackID, position_t) pairs."
        )

    state_lookup = sample_obs[["TrackID", "position_t", INTRINSIC_STATE_COL]].copy()
    state_lookup["ClusterID"] = (
        state_lookup[INTRINSIC_STATE_COL]
        .astype(str)
        .map({str(label): int(code) for label, code in code_map.items()})
    )
    state_lookup = state_lookup.drop(columns=[INTRINSIC_STATE_COL])

    df_tracks_clustered = df_tracks_full.merge(
        state_lookup,
        on=["TrackID", "position_t"],
        how="inner",
        validate="many_to_one",
    )
    if len(df_tracks_clustered) == 0:
        raise ValueError(
            "No overlapping rows found between track features and HMM states for sample "
            f"'{sample_name}'. "
            f"hmm_rows={diag_stats['hmm_rows']}, "
            f"hmm_tracks={diag_stats['hmm_tracks']}, "
            f"track_feature_rows={diag_stats['track_feature_rows']}, "
            f"track_feature_tracks={diag_stats['track_feature_tracks']}, "
            f"overlap_rows={diag_stats['overlap_rows']}, "
            f"overlap_tracks={diag_stats['overlap_tracks']}."
        )

    raw_img = load_image(raw_path)
    tracked_img = load_image(tracked_path)
    filtered_track_img = _filter_track_image_to_ids(
        tracked_img,
        df_tracks_clustered["TrackID"].dropna().unique(),
    )
    backprojected = backproject_columns(
        track_img=filtered_track_img,
        zarr_outpath=state_path,
        df_tracks_clustered=df_tracks_clustered,
        columns=["ClusterID"],
        mode="timepoint",
        track_col="TrackID",
        time_col="position_t",
        background_value=0,
        n_workers=max(1, int(n_workers)),
        prefer_parallel=True,
        tracked_img_path=tracked_path,
        verbose=verbose,
        overwrite=True,
    )
    backprojected["ClusterID"]["type"] = "label"
    state_colors = _normalize_label_color_map(
        code_map.keys(),
        colors=_get_classification_state_colors(adata, INTRINSIC_STATE_COL),
    )
    _write_state_color_attrs_to_zarr(state_path, code_map, state_colors=state_colors)

    backproject_data = {
        "raw_data": {"img": raw_img, "type": "image"},
        "TrackID": {"img": filtered_track_img, "type": "label"},
        "ClusterID": backprojected["ClusterID"],
    }
    viewer = view_napari(
        backproject_data=backproject_data,
        df_tracks_full=df_tracks_full,
        df_tracks_clustered=df_tracks_clustered,
        cell_type=cell_type,
        elsize=_sample_voxel_size(metadata, sample_name),
        split_raw_channels=True,
        run=False,
    )
    code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    try:
        _apply_state_code_colors_to_layer(viewer.layers["ClusterID"], code_colors)
    except Exception:
        pass

    label_map = {str(code): str(label) for label, code in code_map.items()}
    mapping_text = _build_state_mapping_text(label_map, code_colors=code_colors)
    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=mapping_text,
        label_map=label_map,
        code_colors=code_colors,
        title="Intrinsic HMM State Mapping",
    )
    if (not added_dock) and verbose:
        print(mapping_text)

    if verbose:
        print(
            "Opened intrinsic HMM backprojection for sample "
            f"'{sample_name}' with raw='{Path(raw_path).name}' shape={tuple(int(v) for v in raw_img.shape)}, "
            f"tracked='{Path(tracked_path).name}' shape={tuple(int(v) for v in filtered_track_img.shape)}, "
            f"states='{Path(state_path).name}'."
        )

    if bool(run):
        import napari

        napari.run()
    return viewer


def _show_hmm_state_backprojection(
    *,
    adata,
    sample_name,
    output_dir,
    cell_type,
    state_col,
    track_features_csv_path,
    metadata=None,
    state_colors=None,
    n_workers=4,
    run=True,
    verbose=True,
):
    if adata is None or not hasattr(adata, "obs"):
        raise ValueError("model_adata with HMM states is required for backprojection.")
    if state_col not in adata.obs.columns:
        raise ValueError(f"model_adata is missing '{state_col}'.")

    sample_name = str(sample_name).strip()
    cell_type = str(cell_type).strip()
    output_dir = Path(output_dir)

    sample_obs = adata.obs[adata.obs["sample_name"].astype(str) == sample_name].copy()
    if len(sample_obs) == 0:
        raise ValueError(f"No rows found for sample '{sample_name}' in model_adata.")

    raw_path = _resolve_raw_image_path(output_dir=output_dir, sample_name=sample_name, verbose=verbose)
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(f"Could not find raw image for sample '{sample_name}'.")

    tracked_path = _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            f"Could not find tracked image for sample '{sample_name}' and cell_type '{cell_type}'."
        )

    state_path = _behavioral_state_backprojection_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
    )
    code_map = _build_code_map(
        adata.obs,
        state_col=state_col,
        state_order=_get_classification_state_order(adata, state_col),
    )
    if track_features_csv_path is None or not Path(track_features_csv_path).exists():
        raise FileNotFoundError("Track-features CSV is required for HMM backprojection.")

    df_tracks_full = pd.read_csv(track_features_csv_path, low_memory=False)
    df_tracks_full = df_tracks_full[df_tracks_full["sample_name"].astype(str) == sample_name].copy()
    if len(df_tracks_full) == 0:
        raise ValueError(f"No track-feature rows found for sample '{sample_name}'.")

    required_pos_cols = ["TrackID", "position_t", "position_z", "position_y", "position_x"]
    missing_pos_cols = [c for c in required_pos_cols if c not in df_tracks_full.columns]
    if missing_pos_cols:
        raise ValueError(
            "Track-features CSV is missing required position columns for backprojection: "
            f"{missing_pos_cols}"
        )

    state_lookup = sample_obs[["TrackID", "position_t", state_col]].copy()
    state_lookup["ClusterID"] = (
        state_lookup[state_col]
        .astype(str)
        .map({str(label): int(code) for label, code in code_map.items()})
    )
    state_lookup = state_lookup.drop(columns=[state_col])

    df_tracks_clustered = df_tracks_full.merge(
        state_lookup,
        on=["TrackID", "position_t"],
        how="inner",
        validate="many_to_one",
    )
    if len(df_tracks_clustered) == 0:
        raise ValueError(
            "No overlapping rows found between track features and HMM states for sample "
            f"'{sample_name}'."
        )

    raw_img = load_image(raw_path)
    tracked_img = load_image(tracked_path)
    filtered_track_img = _filter_track_image_to_ids(
        tracked_img,
        df_tracks_clustered["TrackID"].dropna().unique(),
    )
    backprojected = backproject_columns(
        track_img=filtered_track_img,
        zarr_outpath=state_path,
        df_tracks_clustered=df_tracks_clustered,
        columns=["ClusterID"],
        mode="timepoint",
        track_col="TrackID",
        time_col="position_t",
        background_value=0,
        n_workers=max(1, int(n_workers)),
        prefer_parallel=True,
        tracked_img_path=tracked_path,
        verbose=verbose,
        overwrite=True,
    )
    backprojected["ClusterID"]["type"] = "label"
    _write_state_color_attrs_to_zarr(state_path, code_map, state_colors=state_colors)

    backproject_data = {
        "raw_data": {"img": raw_img, "type": "image"},
        "TrackID": {"img": filtered_track_img, "type": "label"},
        "ClusterID": backprojected["ClusterID"],
    }
    viewer = view_napari(
        backproject_data=backproject_data,
        df_tracks_full=df_tracks_full,
        df_tracks_clustered=df_tracks_clustered,
        cell_type=cell_type,
        elsize=_sample_voxel_size(metadata, sample_name),
        split_raw_channels=True,
        run=False,
    )
    code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    try:
        _apply_state_code_colors_to_layer(viewer.layers["ClusterID"], code_colors)
    except Exception:
        pass

    label_map = {str(code): str(label) for label, code in code_map.items()}
    mapping_text = _build_state_mapping_text(label_map, code_colors=code_colors)
    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=mapping_text,
        label_map=label_map,
        code_colors=code_colors,
        title="Behavioral State Mapping",
    )
    if (not added_dock) and verbose:
        print(mapping_text)

    if verbose:
        print(
            "Opened HMM behavioral-state backprojection for sample "
            f"'{sample_name}' with raw='{Path(raw_path).name}' shape={tuple(int(v) for v in raw_img.shape)}, "
            f"tracked='{Path(tracked_path).name}' shape={tuple(int(v) for v in filtered_track_img.shape)}, "
            f"states='{Path(state_path).name}'."
        )

    if bool(run):
        import napari

        napari.run()
    return viewer


def _update_hmm_report_metadata(
    adata,
    *,
    composition_pdf=None,
    composition_auc_csv=None,
    composition_plot_csvs=None,
    composition_error=None,
    transition_dir=None,
    transition_counts_csv=None,
    transition_probs_csv=None,
    transition_error=None,
    pending_reason=None,
):
    if adata is None:
        return
    clustering_meta = adata.uns.get("clustering", {})
    if not isinstance(clustering_meta, dict):
        clustering_meta = {}
    clustering_meta["state_composition_report_pdf"] = composition_pdf
    clustering_meta["state_composition_report_auc_csv"] = composition_auc_csv
    clustering_meta["state_composition_report_plot_csvs"] = list(composition_plot_csvs or [])
    clustering_meta["state_composition_report_error"] = composition_error
    clustering_meta["state_transition_report_dir"] = transition_dir
    clustering_meta["state_transition_matrix_counts_csv"] = transition_counts_csv
    clustering_meta["state_transition_matrix_probs_csv"] = transition_probs_csv
    clustering_meta["state_transition_report_error"] = transition_error
    clustering_meta["state_reports_ready"] = bool(
        (composition_pdf is not None or transition_dir is not None)
        and (composition_error is None)
        and (transition_error is None)
    )
    clustering_meta["state_reports_reason"] = (
        None
        if bool(clustering_meta["state_reports_ready"])
        else (
            str(pending_reason).strip()
            if pending_reason is not None and len(str(pending_reason).strip()) > 0
            else "State reports have not been created yet. Use 'Create analysis plots' in the HMM widget."
        )
    )
    adata.uns["clustering"] = clustering_meta


class StateClassificationHMMPanel(BaseStateClassificationPanel):
    """Widget for HMM state discovery, relabeling, deployment, and backprojection."""

    def __init__(self, metadata_loader, cell_type=None):
        self._hmm_model_load_status = {"status": "not_checked", "path": None, "message": None}
        self._full_color_pickers = {}
        self._build_analysis_plots_section()
        super().__init__(metadata_loader=metadata_loader, cell_type=cell_type)
        if not hasattr(self, "_hmm_deployment_artifact"):
            self._hmm_deployment_artifact = None
        self._build_steps()
        children = list(self.ui.children)
        if len(children) > 0:
            children[-1] = self.steps
            self.ui.children = tuple(children)
        self.hmm_n_states_mode.observe(self._on_hmm_n_states_mode_changed, names="value")
        self.apply_hmm_artifact_picker.text.observe(self._on_apply_path_changed, names="value")
        self._sync_hmm_state_controls()
        self._load_existing_hmm_deployment_artifact_if_available()
        self.apply_existing_state_classification = widgets.Checkbox(
            description="Apply existing behavioral state classification",
            value=False,
            indent=False,
        )
        self.apply_existing_state_classification.observe(
            self._on_apply_existing_changed, names="value"
        )
        self.output_dir_html.layout.display = "none"
        children = list(self.ui.children)
        if len(children) > 0:
            children[0] = widgets.HTML(
                '<div style="font-size:20px;font-weight:700;">Behavioral state classification</div>'
            )
            children.insert(2, self.apply_existing_state_classification)
            self.ui.children = tuple(children)

    def _panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.setdefault("behavioral_state_classification", {})
        cell_type = self._current_cell_type()
        if cell_type not in section:
            section[cell_type] = {}
        return section[cell_type]

    def _effective_panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.get("behavioral_state_classification", {})
        defaults = section.get("defaults", {})
        cell_cfg = section.get(self._current_cell_type(), {})
        return {**defaults, **cell_cfg}


    def _build_steps(self):
        analysis_plots_section = getattr(self, "analysis_plots_section", widgets.VBox([]))
        self._rename_full_with_description = widgets.VBox(
            [
                widgets.HTML(
                    "<span style='color:#555;'>Rename HMM intrinsic clusters combined with binary group "
                    "values (e.g. organoid_contact)</span>"
                ),
                self.rename_full_section,
            ],
            layout=widgets.Layout(gap="6px"),
        )
        self.steps = widgets.Accordion(
            children=[
                self.clustering_section,
                self.rename_intrinsic_section,
                self._rename_full_with_description,
                analysis_plots_section,
                self.backprojection_section,
            ],
            selected_index=None,
        )
        self.steps.set_title(0, "Assign HMM intrinsic behavioral states")
        self.steps.set_title(1, "Rename intrinsic states")
        self.steps.set_title(2, "Rename full states")
        self.steps.set_title(3, "Create plots")
        self.steps.set_title(4, "Backprojection")

    def _build_apply_section(self):
        super()._build_apply_section()
        self.apply_hmm_artifact_picker = self._make_hmm_artifact_picker()
        self.apply_hmm_artifact_picker.filter_pattern = "*.pkl"
        self.apply_hmm_default_paths_html = widgets.HTML("")
        self.btn_apply_hmm_artifact = widgets.Button(
            description="Apply saved HMM artifact",
            button_style="success",
            layout=widgets.Layout(width="220px"),
        )
        self.btn_apply_hmm_artifact.on_click(self._on_apply_hmm_artifact_clicked)
        self.apply_hmm_spinner = widgets.HTML(value=spinning_loader)
        self.apply_hmm_spinner.layout.display = "none"
        self.out_apply_hmm = widgets.Output()
        self.hmm_artifact_section = widgets.VBox(
            [
                widgets.HTML("<hr style='margin:10px 0;'>"),
                widgets.HTML("<b>Apply saved HMM deployment artifact</b>"),
                self.apply_hmm_artifact_picker,
                self.apply_hmm_default_paths_html,
                widgets.HBox([self.btn_apply_hmm_artifact, self.apply_hmm_spinner]),
                self.out_apply_hmm,
            ]
        )
        self.apply_section = widgets.VBox(
            [
                widgets.HTML("<b>Apply saved classifier artifacts</b>"),
                self.apply_full_pkl_picker,
                self.apply_intrinsic_pkl_picker,
                self.apply_default_paths_html,
                widgets.HBox([self.btn_apply, self.apply_spinner]),
                self.out_apply,
            ]
        )

    def _make_hmm_artifact_picker(self):
        from behav3d.widgets.utils import PathPicker

        return PathPicker(
            mode="file",
            start_dir=self.output_dir or ".",
            default="",
            description="HMM PKL",
            placeholder="Path to saved HMM deployment artifact .pkl",
            width="100%",
        )

    def _default_hmm_deployment_artifact_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return _resolve_hmm_deployment_artifact_path(output_dir=self.output_dir, cell_type=ct)

    def _load_existing_hmm_deployment_artifact_if_available(self):
        self._hmm_deployment_artifact = None
        candidates = []
        picker_value = str(getattr(self, "apply_hmm_artifact_picker", widgets.Text(value="")).value).strip()
        if picker_value != "":
            candidates.append(Path(picker_value))
        candidates.append(self._default_hmm_deployment_artifact_path())
        seen = set()
        for candidate in candidates:
            candidate = Path(candidate)
            candidate_key = str(candidate)
            if candidate_key in seen or not candidate.exists():
                continue
            seen.add(candidate_key)
            try:
                self._hmm_deployment_artifact = load_hmm_deployment_artifact(candidate)
                return
            except Exception:
                self._hmm_deployment_artifact = None

    def _is_hmm_model_adata(self, adata):
        if adata is None or not hasattr(adata, "uns"):
            return False
        clustering_meta = adata.uns.get("clustering", {})
        preprocessing_meta = adata.uns.get("preprocessing", {})
        classification_meta = adata.uns.get("classification", {})
        return (
            (isinstance(clustering_meta, dict) and str(clustering_meta.get("clustering_method", "")).lower() == "hmm")
            or (
                isinstance(preprocessing_meta, dict)
                and str(preprocessing_meta.get("observation_mode", "")).lower() == "timepoint_hmm"
            )
            or (
                isinstance(classification_meta, dict)
                and str(classification_meta.get("intrinsic_output_col", "")) == INTRINSIC_STATE_COL
            )
        )

    def _labels_look_like_raw_hmm_states(self, labels):
        values = pd.Series(labels, dtype="string").dropna().astype(str).str.strip()
        values = values[values != ""]
        if len(values) == 0:
            return False
        return bool(values.map(lambda v: v.isdigit()).all())

    def _normalize_loaded_hmm_model_adata(self):
        if self.model_adata is None or not hasattr(self.model_adata, "obs"):
            return {"normalized": False, "reason": "no_model"}
        if not self._is_hmm_model_adata(self.model_adata):
            return {"normalized": False, "reason": "not_hmm"}

        obs = self.model_adata.obs
        changes = []
        if INTRINSIC_STATE_COL not in obs.columns and "intrinsic_behavioral_cluster" in obs.columns:
            obs[INTRINSIC_STATE_COL] = pd.Categorical(
                pd.Series(obs["intrinsic_behavioral_cluster"], index=obs.index, dtype="string")
            )
            changes.append(INTRINSIC_STATE_COL)
        if FULL_STATE_COL not in obs.columns and "full_behavioral_cluster" in obs.columns:
            obs[FULL_STATE_COL] = pd.Categorical(
                pd.Series(obs["full_behavioral_cluster"], index=obs.index, dtype="string")
            )
            changes.append(FULL_STATE_COL)
        if (
            HMM_INTRINSIC_RAW_STATE_COL not in obs.columns
            and INTRINSIC_STATE_COL in obs.columns
            and self._labels_look_like_raw_hmm_states(obs[INTRINSIC_STATE_COL])
        ):
            obs[HMM_INTRINSIC_RAW_STATE_COL] = pd.Categorical(
                pd.Series(obs[INTRINSIC_STATE_COL], index=obs.index, dtype="string")
            )
            changes.append(HMM_INTRINSIC_RAW_STATE_COL)
        return {"normalized": len(changes) > 0, "changes": list(changes)}

    def _load_existing_model_if_available(self):
        self.model_adata = None
        p = self._model_adata_path()
        self._hmm_model_load_status = {"status": "missing", "path": str(p), "message": None}
        if p.exists():
            try:
                self.model_adata = sc.read_h5ad(p)
                normalize_result = self._normalize_loaded_hmm_model_adata()
                if INTRINSIC_STATE_COL in getattr(self.model_adata, "obs", {}).columns:
                    self.model_adata = self._merge_metadata_into_obs(self.model_adata)
                    message = None
                    if bool(normalize_result.get("normalized", False)):
                        message = "Normalized older HMM model columns: " + ", ".join(normalize_result["changes"])
                    self._hmm_model_load_status = {"status": "loaded", "path": str(p), "message": message}
                elif self._is_hmm_model_adata(self.model_adata):
                    self._hmm_model_load_status = {
                        "status": "missing_intrinsic_col",
                        "path": str(p),
                        "message": f"Loaded HMM model is missing '{INTRINSIC_STATE_COL}'.",
                    }
                else:
                    self._hmm_model_load_status = {
                        "status": "not_hmm",
                        "path": str(p),
                        "message": "Loaded model adata does not look like an HMM model.",
                    }
            except Exception as exc:
                self.model_adata = None
                self._hmm_model_load_status = {"status": "failed", "path": str(p), "message": str(exc)}
        self._load_existing_hmm_deployment_artifact_if_available()

    def _build_clustering_section(self):
        self.feature_groups_status = widgets.HTML("<i>No features loaded yet.</i>")
        self.feature_groups_box = widgets.VBox([])
        self.selected_features_box = widgets.Select(
            options=[],
            rows=10,
            layout=widgets.Layout(width="520px", height="180px"),
        )
        self.hmm_log_scale_status = widgets.HTML("<i>Select features to choose optional log1p scaling.</i>")
        self.hmm_log_scale_box = widgets.VBox([])
        self._hmm_log_scale_checkboxes = {}
        # Distribution preview: sits just below the "Log scaling" header and
        # above the per-feature checkboxes, so the user can inspect which
        # features are skewed before applying log scaling.
        self.hmm_log_scale_dist_btn = widgets.Button(
            description="📊 Preview feature distributions",
            layout=widgets.Layout(width="260px"),
        )
        self.hmm_log_scale_dist_btn.on_click(self._on_preview_log_scale_distributions)
        self.hmm_log_scale_dist_out = widgets.Output()
        self.binary_group_status = widgets.HTML("<i>No binary columns detected yet.</i>")
        self.binary_group_checks_box = widgets.VBox([])
        self._binary_group_checkboxes = {}

        # Keep hidden compatibility widgets so the parent class can initialize cleanly.
        self.window_size = widgets.IntText(value=1)
        self.min_spacing = widgets.Text(value="")
        self.max_samples = widgets.Text(value="")
        self.n_neighbors = widgets.IntText(value=30)
        self.min_dist = widgets.FloatText(value=0.1)
        self.resolution = widgets.FloatText(value=0.2)
        self.pca_var_selection = widgets.FloatText(value=0.95)
        self.use_pca = widgets.Checkbox(value=False, indent=False)
        self.clustering_method = widgets.Dropdown(options=["hmm"], value="hmm")
        self.incomplete_window_policy = widgets.Dropdown(options=["timepoint"], value="timepoint")
        self.reuse_prepared_dataset = widgets.Checkbox(value=False, indent=False)
        self.describe_window_feature_cbs = {}
        self.additional_window_feature_cbs = {
            "net_displacement": widgets.Checkbox(description="net_displacement", value=False, indent=False),
            "straightness": widgets.Checkbox(description="straightness", value=False, indent=False),
            "mean_square_displacement": widgets.Checkbox(
                description="mean_square_displacement",
                value=False,
                indent=False,
            ),
        }
        for cb in self.additional_window_feature_cbs.values():
            cb.observe(self._on_feature_checkbox_changed, names="value")
        self.hmm_window_features_window = widgets.IntText(
            description="Window size (timepoints)",
            value=5,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="230px"),
        )
        self.hmm_window_features_window.observe(
            self._on_time_conversion_widget_changed, names="value"
        )
        self.hmm_window_features_help = widgets.Button(
            description="?",
            tooltip=(
                "Window size is in timepoints (frames). The line below shows "
                "the equivalent real time, computed from this sample's "
                "time_interval/time_unit metadata."
            ),
            layout=widgets.Layout(width="24px"),
        )
        self.hmm_window_features_time_html = widgets.HTML(value="")

        self.hmm_n_states_mode = widgets.Dropdown(
            description="State selection",
            options=[("Fixed", "fixed"), ("Auto (BIC)", "auto")],
            value="fixed",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.hmm_advanced = widgets.Checkbox(
            description="Advanced",
            value=False,
            indent=False,
        )
        self.hmm_advanced.observe(self._on_hmm_advanced_changed, names="value")
        self.hmm_n_states = widgets.IntText(
            description="n_states",
            value=4,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="170px"),
        )
        self.hmm_k_min = widgets.IntText(
            description="k_min",
            value=2,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.hmm_k_max = widgets.IntText(
            description="k_max",
            value=8,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.hmm_covariance_type = widgets.Dropdown(
            description="Covariance",
            options=["full", "diag", "spherical", "tied"],
            value="full",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.hmm_n_iter = widgets.IntText(
            description="n_iter",
            value=200,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="160px"),
        )
        self.hmm_tol = widgets.FloatText(
            description="tol",
            value=1e-3,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.hmm_sticky = widgets.Checkbox(
            description="Sticky HMM",
            value=False,
            indent=False,
        )
        self.hmm_sticky.observe(self._on_hmm_sticky_changed, names="value")
        self.hmm_stickiness_kappa = widgets.FloatText(
            description="kappa",
            value=8.0,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.hmm_transmat_alpha = widgets.FloatText(
            description="alpha",
            value=1.0,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.hmm_min_covar = widgets.FloatText(
            description="min_covar",
            value=1e-3,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="170px"),
        )
        self.hmm_feature_smoothing_window = widgets.IntText(
            description="Smooth window (timepoints)",
            value=1,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="250px"),
        )
        self.hmm_feature_smoothing_window.observe(
            self._on_time_conversion_widget_changed, names="value"
        )
        self.hmm_feature_smoothing_help = widgets.Button(
            description="?",
            tooltip=(
                "Smoothing window is in timepoints (frames). The line below "
                "shows the equivalent real time, computed from this sample's "
                "time_interval/time_unit metadata."
            ),
            layout=widgets.Layout(width="24px"),
        )
        self.hmm_feature_smoothing_time_html = widgets.HTML(value="")
        self.hmm_smoothing_min_periods = widgets.IntText(
            description="Min periods",
            value=1,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="180px"),
        )
        self.hmm_start_offset = widgets.BoundedIntText(
            description="Start offset",
            value=0,
            min=0,
            max=100000,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.hmm_start_offset_fill_mode = widgets.Dropdown(
            description="Skipped frames",
            options=[("Backfill", "backfill"), ("Leave unassigned", "leave_unassigned")],
            value="backfill",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="240px"),
        )
        self.feature_quantile_capping_low_percentile = widgets.Text(
            description="Low percentile cap",
            value="",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.feature_quantile_capping_high_percentile = widgets.Text(
            description="High percentile cap",
            value="0.99",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        # Aliases so legacy base-class methods (_apply_cfg_defaults, _persist_current_settings)
        # that still reference the old names continue to work against the same widget objects.
        self.lower_quantile_cap = self.feature_quantile_capping_low_percentile
        self.upper_quantile_cap = self.feature_quantile_capping_high_percentile
        self.random_state = widgets.IntText(
            description="Seed",
            value=123,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )

        self.btn_cluster = widgets.Button(
            description="Run intrinsic HMM clustering",
            button_style="success",
            layout=widgets.Layout(width="240px"),
        )
        self.btn_cluster.on_click(self._on_cluster_clicked)
        self.cluster_spinner = widgets.HTML(value=spinning_loader)
        self.cluster_spinner.layout.display = "none"
        self.out_cluster = widgets.Output()

        feature_select_block = widgets.VBox(
            [
                self.feature_groups_status,
                self.feature_groups_box,
            ]
        )
        window_features_block = widgets.VBox(
            [
                widgets.HTML("<b>window features</b>"),
                widgets.HBox(
                    [self.hmm_window_features_window, self.hmm_window_features_help],
                    layout=widgets.Layout(align_items="center", gap="4px"),
                ),
                self.hmm_window_features_time_html,
                widgets.VBox(list(self.additional_window_feature_cbs.values())),
            ]
        )
        feature_selection_children = widgets.Accordion(
            children=[feature_select_block, window_features_block],
            selected_index=None,
        )
        feature_selection_children.set_title(0, "timepoint features")
        feature_selection_children.set_title(1, "window features")
        feature_selection_block = widgets.VBox(
            [
                feature_selection_children,
                widgets.HTML("<b>selected features</b>"),
                self.selected_features_box,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.hmm_log_scale_accordion = widgets.Accordion(
            children=[
                widgets.VBox(
                    [
                        self.hmm_log_scale_dist_btn,
                        self.hmm_log_scale_dist_out,
                        self.hmm_log_scale_status,
                        self.hmm_log_scale_box,
                    ],
                    layout=widgets.Layout(gap="6px"),
                )
            ],
            selected_index=None,
        )
        self.hmm_log_scale_accordion.set_title(0, "Log scaling")
        self.hmm_quantile_row = widgets.HBox(
            [self.feature_quantile_capping_low_percentile, self.feature_quantile_capping_high_percentile],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        feature_processing_block = widgets.VBox(
            [
                self.hmm_log_scale_accordion,
                widgets.HTML("<b>feature processing</b>"),
                widgets.HBox(
                    [self.hmm_feature_smoothing_window, self.hmm_feature_smoothing_help],
                    layout=widgets.Layout(align_items="center", gap="4px"),
                ),
                self.hmm_feature_smoothing_time_html,
                self.hmm_quantile_row,
            ]
        )
        binary_group_block = widgets.VBox(
            [
                widgets.HTML("<b>binary group selection</b>"),
                self.binary_group_status,
                self.binary_group_checks_box,
            ]
        )
        self.feature_selection_subaccordion = widgets.Accordion(
            children=[feature_selection_block, feature_processing_block, binary_group_block],
            selected_index=None,
        )
        self.feature_selection_subaccordion.set_title(0, "Feature selection")
        self.feature_selection_subaccordion.set_title(1, "Feature processing")
        self.feature_selection_subaccordion.set_title(2, "Binary group selection")

        self.hmm_primary_row = widgets.HBox(
            [self.hmm_n_states, self.random_state, self.hmm_advanced],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.hmm_state_mode_row = widgets.HBox(
            [self.hmm_n_states_mode, self.hmm_k_min, self.hmm_k_max],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.hmm_advanced_row = widgets.HBox(
            [
                self.hmm_start_offset,
                self.hmm_start_offset_fill_mode,
                self.hmm_covariance_type,
                self.hmm_n_iter,
                self.hmm_tol,
                self.hmm_min_covar,
                self.hmm_sticky,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.hmm_sticky_params_row = widgets.HBox(
            [self.hmm_stickiness_kappa, self.hmm_transmat_alpha],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.hmm_controls_block = widgets.VBox(
            [
                self.hmm_primary_row,
                self.hmm_state_mode_row,
                self.hmm_advanced_row,
                self.hmm_sticky_params_row,
            ],
            layout=widgets.Layout(gap="8px"),
        )

        self.clustering_section = widgets.VBox(
            [
                self.feature_selection_subaccordion,
                widgets.HTML("<b>Intrinsic HMM clustering</b>"),
                self.hmm_controls_block,
                widgets.HBox([self.btn_cluster, self.cluster_spinner]),
                self.out_cluster,
            ]
        )

    def _build_intrinsic_rename_section(self):
        self.rename_intrinsic_status = widgets.HTML("<i>Run HMM clustering or load existing model first.</i>")
        self.rename_intrinsic_rows = widgets.VBox([])
        self.rename_intrinsic_rows.layout = widgets.Layout(width="560px")
        self.intrinsic_combine_name = widgets.Text(
            description="Combine to",
            value="",
            placeholder="New intrinsic state name",
            style={"description_width": "90px"},
            layout=widgets.Layout(width="360px"),
        )
        self.btn_combine_intrinsic = widgets.Button(
            description="combine",
            button_style="info",
            layout=widgets.Layout(width="110px"),
        )
        self.btn_combine_intrinsic.on_click(self._on_combine_intrinsic_clicked)
        self.combine_intrinsic_spinner = widgets.HTML(value=spinning_loader)
        self.combine_intrinsic_spinner.layout.display = "none"
        self.intrinsic_combine_box = widgets.VBox(
            [
                widgets.HTML("<b>Combine selected</b>"),
                self.intrinsic_combine_name,
                widgets.HBox([self.btn_combine_intrinsic, self.combine_intrinsic_spinner]),
            ],
            layout=widgets.Layout(width="390px", align_items="flex-start"),
        )
        self.btn_rename_intrinsic = widgets.Button(
            description="Rename intrinsic clusters",
            button_style="warning",
            layout=widgets.Layout(width="240px"),
        )
        self.btn_rename_intrinsic.on_click(self._on_rename_intrinsic_clicked)
        self.rename_intrinsic_spinner = widgets.HTML(value=spinning_loader)
        self.rename_intrinsic_spinner.layout.display = "none"
        self.out_rename_intrinsic = widgets.Output()
        self.rename_intrinsic_section = widgets.VBox(
            [
                self.rename_intrinsic_status,
                widgets.HBox(
                    [self.rename_intrinsic_rows, self.intrinsic_combine_box],
                    layout=widgets.Layout(align_items="flex-start", gap="14px"),
                ),
                widgets.HBox([self.btn_rename_intrinsic, self.rename_intrinsic_spinner]),
                self.out_rename_intrinsic,
            ]
        )

    def _build_backprojection_section(self):
        self.backprojection_status = widgets.HTML("<i>No samples detected yet.</i>")
        self.backproj_sample_dd = widgets.Dropdown(
            description="Sample",
            options=[],
            value=None,
            layout=widgets.Layout(width="360px"),
            style={"description_width": "90px"},
        )
        self.hmm_backprojection_workers = widgets.BoundedIntText(
            description="Workers",
            value=4,
            min=1,
            max=32,
            style={"description_width": "90px"},
            layout=widgets.Layout(width="170px"),
        )
        self.btn_open_intrinsic_backprojection = widgets.Button(
            description="Open intrinsic HMM backprojection",
            button_style="success",
            layout=widgets.Layout(width="270px"),
        )
        self.btn_open_intrinsic_backprojection.on_click(self._on_open_intrinsic_backprojection_clicked)
        self.btn_open_backprojection = widgets.Button(
            description="Open full cluster backprojection",
            button_style="success",
            layout=widgets.Layout(width="250px"),
        )
        self.btn_open_backprojection.on_click(self._on_open_backprojection_clicked)
        self.backprojection_spinner = widgets.HTML(value=spinning_loader)
        self.backprojection_spinner.layout.display = "none"
        self.out_backprojection = widgets.Output()
        self.backprojection_section = widgets.VBox(
            [
                self.backprojection_status,
                widgets.HBox([self.backproj_sample_dd, self.hmm_backprojection_workers]),
                widgets.HBox(
                    [
                        self.btn_open_intrinsic_backprojection,
                        self.btn_open_backprojection,
                        self.backprojection_spinner,
                    ]
                ),
                self.out_backprojection,
            ]
        )

    def _detect_cell_types(self):
        md = getattr(self.metadata_loader, "metadata", None)
        cell_types = []
        if md is not None:
            try:
                cell_types.extend(filter_multicolor_inputs(detect_organoid_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)))
            except Exception:
                pass

        # Filesystem fallback
        out_dir = Path(self.output_dir) if self.output_dir else None
        if out_dir is not None:
            analysis_dir = out_dir / "analysis"
            if analysis_dir.exists():
                for p in analysis_dir.iterdir():
                    if p.is_dir():
                        cell_types.append(p.name)
        return sorted({str(x).strip() for x in cell_types if str(x).strip() != ""})

    def _detect_sample_names(self):
        md = getattr(self.metadata_loader, "metadata", None)
        if isinstance(md, pd.DataFrame) and "sample_name" in md.columns:
            meta_names = sorted(
                {
                    str(x).strip()
                    for x in md["sample_name"].astype(str).dropna().unique().tolist()
                    if str(x).strip() != ""
                }
            )
            if len(meta_names) > 0:
                return meta_names

        sample_names = []
        if self.model_adata is not None and hasattr(self.model_adata, "obs"):
            if "sample_name" in self.model_adata.obs.columns:
                sample_names.extend(
                    self.model_adata.obs["sample_name"].astype(str).dropna().unique().tolist()
                )

        out_dir = Path(self.output_dir) if self.output_dir else None
        images_dir = (out_dir / "images") if out_dir is not None else None
        if images_dir is not None and images_dir.exists():
            for p in images_dir.iterdir():
                if p.is_dir():
                    sample_names.append(str(p.name))

        return sorted({str(x).strip() for x in sample_names if str(x).strip() != ""})

    def _build_analysis_plots_section(self):
        self.analysis_plots_status = widgets.HTML(
            "<i>Run HMM clustering and finish curated renaming first, then create plots manually.</i>"
        )
        self.btn_create_state_composition_plots = widgets.Button(
            description="Create state composition plots",
            button_style="info",
            layout=widgets.Layout(width="260px"),
        )
        self.btn_create_state_composition_plots.on_click(self._on_create_state_composition_plots_clicked)
        self.state_composition_spinner = widgets.HTML(value=spinning_loader)
        self.state_composition_spinner.layout.display = "none"
        self.composition_group_x_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in X:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.composition_group_y_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in Y:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.composition_group_cols_select = widgets.SelectMultiple(
            options=[],
            value=[],
            description="",
            rows=2,
            layout=widgets.Layout(width="360px"),
            disabled=True,
        )
        self.btn_create_state_transition_plots = widgets.Button(
            description="Create state transition plots",
            button_style="info",
            layout=widgets.Layout(width="250px"),
        )
        self.btn_create_state_transition_plots.on_click(self._on_create_state_transition_plots_clicked)
        self.state_transition_spinner = widgets.HTML(value=spinning_loader)
        self.state_transition_spinner.layout.display = "none"
        self.comparison_condition_col_dd = widgets.Dropdown(
            options=[], description="Compare condition:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.comparison_condition_col_dd.observe(self._on_comparison_condition_col_changed, names="value")
        self.comparison_group_x_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in X:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.comparison_group_y_text = widgets.Text(
            value="", description="Group in Y:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.comparison_group_cols_select = widgets.SelectMultiple(
            options=[],
            value=[],
            description="",
            rows=2,
            layout=widgets.Layout(width="360px"),
            disabled=True,
        )
        self.btn_create_condition_comparison = widgets.Button(
            description="Create condition comparison plot",
            button_style="info",
            layout=widgets.Layout(width="260px"),
            disabled=True,
        )
        self.btn_create_condition_comparison.on_click(self._on_create_condition_comparison_clicked)
        self.condition_comparison_spinner = widgets.HTML(value=spinning_loader)
        self.condition_comparison_spinner.layout.display = "none"
        self.out_analysis_plots = widgets.Output()

        state_composition_box = build_plot_box(
            title="State composition plots",
            description=(
                "Summarizes how full behavioral clusters are distributed per sample "
                "and in pooled overviews across time."
            ),
            run_row=widgets.HBox(
                [self.btn_create_state_composition_plots, self.state_composition_spinner],
                layout=widgets.Layout(align_items="center", gap="8px"),
            ),
            settings=[
                widgets.HTML(
                    "<b>Group by condition:</b><br>"
                    "<span style='color:#555;font-size:11px;'>\"Group in X\"/\"Group in Y\" pick a single "
                    "condition each to arrange the grouped-composition grid. \"Group per page\" pools "
                    "additional metadata columns into that grid instead — hold Ctrl/Cmd to select multiple.</span>"
                ),
                self.composition_group_x_dd,
                self.composition_group_y_dd,
                widgets.HTML("<span style='color:#555;font-size:11px;'>Group per page:</span>"),
                self.composition_group_cols_select,
            ],
        )
        state_transition_box = build_plot_box(
            title="State transition plots",
            description=(
                "Builds transition matrices and transition summaries between full "
                "behavioral clusters along tracks."
            ),
            run_row=widgets.HBox(
                [self.btn_create_state_transition_plots, self.state_transition_spinner],
                layout=widgets.Layout(align_items="center", gap="8px"),
            ),
        )
        condition_comparison_box = build_plot_box(
            title="Condition comparison plot",
            description=(
                "Per-cluster overall proportion difference between every pairwise combination of a "
                "condition's levels (Welch's t-test), shown as one row per pairwise comparison."
            ),
            run_row=widgets.HBox(
                [self.btn_create_condition_comparison, self.condition_comparison_spinner],
                layout=widgets.Layout(align_items="center", gap="8px"),
            ),
            settings=[
                self.comparison_condition_col_dd,
                widgets.HTML(
                    "<b>Group by condition:</b><br>"
                    "<span style='color:#555;font-size:11px;'>Each row is one pairwise comparison of "
                    "\"Compare condition\"'s levels (shown in \"Group in Y\"). \"Group in X\" splits the "
                    "comparison into side-by-side columns from another condition. \"Group per page\" pools "
                    "additional metadata columns into that same columns axis instead — hold Ctrl/Cmd to "
                    "select multiple.</span>"
                ),
                self.comparison_group_x_dd,
                self.comparison_group_y_text,
                widgets.HTML("<span style='color:#555;font-size:11px;'>Group per page:</span>"),
                self.comparison_group_cols_select,
            ],
        )

        self.analysis_plots_section = widgets.VBox(
            [
                self.analysis_plots_status,
                state_composition_box,
                state_transition_box,
                condition_comparison_box,
                self.out_analysis_plots,
            ]
        )

    def _selected_descriptive_features(self):
        return ["timepoint_hmm"]

    def _selected_window_feature_columns(self):
        return [
            str(col)
            for col, cb in self.additional_window_feature_cbs.items()
            if bool(cb.value)
        ]

    def _selected_hmm_feature_columns(self):
        out = []
        seen = set()
        for feature_name in list(self._selected_feature_columns()) + list(self._selected_window_feature_columns()):
            if feature_name in seen:
                continue
            out.append(str(feature_name))
            seen.add(feature_name)
        return out

    def _update_selected_features_box(self):
        self.selected_features_box.options = self._selected_hmm_feature_columns()

    def _build_feature_groups(self):
        cell_cfg = self._panel_cfg()
        if isinstance(cell_cfg, dict):
            effective = self._effective_panel_cfg()
            if not cell_cfg.get("selected_features") and effective.get("selected_features"):
                cell_cfg["selected_features"] = list(effective["selected_features"])
            if not cell_cfg.get("binary_features_to_group") and effective.get("binary_features_to_group"):
                cell_cfg["binary_features_to_group"] = list(effective["binary_features_to_group"])
        super()._build_feature_groups()
        grouped = getattr(self, "_feature_groups", {}) or {}
        if len(grouped) > 0:
            children = []
            titles = []
            for group_name, feats in grouped.items():
                row = self._group_rows.get(group_name, {})
                child_cbs = list(row.get("child_cbs", []))
                if len(child_cbs) == 0:
                    continue
                grid = widgets.GridBox(
                    child_cbs,
                    layout=widgets.Layout(
                        grid_template_columns="repeat(3, max-content)",
                        grid_gap="2px 10px",
                    ),
                )
                children.append(grid)
                titles.append(group_name)

            if len(children) > 0:
                fg_acc = widgets.Accordion(children=children, selected_index=None)
                for idx, title in enumerate(titles):
                    fg_acc.set_title(idx, title)
                self.feature_groups_box.children = [fg_acc]
        self._update_selected_features_box()
        self._rebuild_log_scale_feature_controls()

    def _sync_pca_controls(self):
        return

    def _selected_log_scale_features(self):
        return [
            str(col)
            for col, cb in self._hmm_log_scale_checkboxes.items()
            if bool(cb.value)
        ]

    def _on_preview_log_scale_distributions(self, *_):
        """Plot histograms of the currently selected features so the user can
        judge which benefit from log scaling before applying it."""
        from behav3d.widgets.utils import build_feature_distribution_figure
        from IPython.display import display

        out = self.hmm_log_scale_dist_out
        out.clear_output(wait=True)
        ct = self._current_cell_type()
        csv_path = self._resolve_track_features_csv()
        feats = self._selected_hmm_feature_columns()
        with out:
            if csv_path is None:
                print(f"No track-features CSV found for '{ct}'. Run feature extraction/filtering first.")
                return
            if not feats:
                print("Select at least one feature above to preview its distribution.")
                return
            fig, truncated = build_feature_distribution_figure(
                csv_path, feats, title=f"Feature distributions — {ct}"
            )
            if fig is None:
                print(f"Could not read features from: {csv_path}")
                return
            if truncated:
                print("Showing the first 36 selected features.")
            display(fig)

    def _rebuild_log_scale_feature_controls(self):
        selected_features = self._selected_hmm_feature_columns()
        cfg = self._effective_panel_cfg()
        saved_log = set(cfg.get("hmm_log_scale_features", [])) if isinstance(cfg, dict) else set()
        previous_log = {
            str(col)
            for col, cb in getattr(self, "_hmm_log_scale_checkboxes", {}).items()
            if bool(cb.value)
        }
        effective_selected = saved_log.union(previous_log)

        self._hmm_log_scale_checkboxes = {}
        if len(selected_features) == 0:
            self.hmm_log_scale_status.value = "<i>Select features to choose optional log1p scaling.</i>"
            self.hmm_log_scale_box.children = []
            return

        children = []
        for feature_name in selected_features:
            cb = widgets.Checkbox(
                description=str(feature_name),
                value=(str(feature_name) in effective_selected),
                indent=False,
                layout=widgets.Layout(width="360px"),
            )
            self._hmm_log_scale_checkboxes[str(feature_name)] = cb
            children.append(cb)
        self.hmm_log_scale_status.value = (
            f"<b>Selected features available for log1p scaling:</b> {len(selected_features)}"
        )
        self.hmm_log_scale_box.children = children

    def _sync_hmm_state_controls(self):
        advanced = bool(self.hmm_advanced.value)
        auto_mode = str(self.hmm_n_states_mode.value) == "auto"
        self.hmm_n_states.disabled = auto_mode
        self.hmm_k_min.disabled = not auto_mode
        self.hmm_k_max.disabled = not auto_mode
        if hasattr(self, "hmm_state_mode_row"):
            self.hmm_state_mode_row.layout.display = None if advanced else "none"
        if hasattr(self, "hmm_advanced_row"):
            self.hmm_advanced_row.layout.display = None if advanced else "none"
        if hasattr(self, "hmm_sticky_params_row"):
            self.hmm_sticky_params_row.layout.display = (
                None if (advanced and bool(self.hmm_sticky.value)) else "none"
            )
        if hasattr(self, "hmm_quantile_row"):
            self.hmm_quantile_row.layout.display = None if advanced else "none"

    def _on_use_pca_changed(self, _):
        return

    def _on_hmm_n_states_mode_changed(self, _):
        self._sync_hmm_state_controls()

    def _on_hmm_advanced_changed(self, _):
        self._sync_hmm_state_controls()

    def _on_hmm_sticky_changed(self, _):
        self._sync_hmm_state_controls()

    def _on_apply_existing_changed(self, _):
        self._sync_apply_existing_mode()

    def _sync_apply_existing_mode(self):
        if not hasattr(self, "steps"):
            return
        analysis_plots_section = getattr(self, "analysis_plots_section", widgets.VBox([]))
        if self.apply_existing_state_classification.value:
            self.steps.children = [
                self.hmm_artifact_section,
                analysis_plots_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Apply existing classification")
            self.steps.set_title(1, "Create plots")
            self.steps.set_title(2, "Backprojection")
            self.steps.selected_index = 0
        else:
            rename_full = getattr(self, "_rename_full_with_description", self.rename_full_section)
            self.steps.children = [
                self.clustering_section,
                self.rename_intrinsic_section,
                rename_full,
                analysis_plots_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Assign HMM intrinsic behavioral states")
            self.steps.set_title(1, "Rename intrinsic states")
            self.steps.set_title(2, "Rename full states")
            self.steps.set_title(3, "Create plots")
            self.steps.set_title(4, "Backprojection")

    def _collapse_all_steps(self):
        if hasattr(self, "steps") and isinstance(self.steps, widgets.Accordion):
            self.steps.selected_index = None

    def _open_step(self, index):
        if hasattr(self, "steps") and isinstance(self.steps, widgets.Accordion):
            n = len(self.steps.children)
            if index is not None and 0 <= index < n:
                self.steps.selected_index = index

    def _on_feature_checkbox_changed(self, _):
        super()._on_feature_checkbox_changed(_)
        self._rebuild_log_scale_feature_controls()

    def _apply_cfg_defaults(self):
        super()._apply_cfg_defaults()
        cfg = self._effective_panel_cfg()
        if not isinstance(cfg, dict):
            return
        saved_artifact_path = str(cfg.get("apply_hmm_deployment_artifact_path", "")).strip()
        if saved_artifact_path:
            self.apply_hmm_artifact_picker.value = saved_artifact_path
        self.hmm_n_states_mode.value = str(cfg.get("hmm_n_states_mode", self.hmm_n_states_mode.value))
        self.hmm_advanced.value = bool(cfg.get("hmm_advanced", self.hmm_advanced.value))
        self.hmm_n_states.value = int(cfg.get("hmm_n_states", self.hmm_n_states.value))
        self.hmm_k_min.value = int(cfg.get("hmm_k_min", self.hmm_k_min.value))
        self.hmm_k_max.value = int(cfg.get("hmm_k_max", self.hmm_k_max.value))
        self.hmm_covariance_type.value = str(cfg.get("hmm_covariance_type", self.hmm_covariance_type.value))
        self.hmm_n_iter.value = int(cfg.get("hmm_n_iter", self.hmm_n_iter.value))
        self.hmm_tol.value = float(cfg.get("hmm_tol", self.hmm_tol.value))
        self.hmm_sticky.value = bool(cfg.get("hmm_sticky", self.hmm_sticky.value))
        self.hmm_stickiness_kappa.value = float(
            cfg.get("hmm_stickiness_kappa", self.hmm_stickiness_kappa.value)
        )
        self.hmm_transmat_alpha.value = float(cfg.get("hmm_transmat_alpha", self.hmm_transmat_alpha.value))
        self.hmm_min_covar.value = float(cfg.get("hmm_min_covar", self.hmm_min_covar.value))
        self.hmm_feature_smoothing_window.value = int(
            cfg.get("hmm_feature_smoothing_window", self.hmm_feature_smoothing_window.value)
        )
        self.hmm_smoothing_min_periods.value = int(
            cfg.get("hmm_smoothing_min_periods", self.hmm_smoothing_min_periods.value)
        )
        self.hmm_start_offset.value = int(cfg.get("hmm_start_offset", self.hmm_start_offset.value))
        self.hmm_start_offset_fill_mode.value = str(
            cfg.get("hmm_start_offset_fill_mode", self.hmm_start_offset_fill_mode.value)
        )
        self.hmm_backprojection_workers.value = int(
            cfg.get("hmm_backprojection_workers", self.hmm_backprojection_workers.value)
        )
        self.hmm_window_features_window.value = int(
            cfg.get("hmm_window_features_window", self.hmm_window_features_window.value)
        )
        self.feature_quantile_capping_low_percentile.value = str(
            cfg.get("feature_quantile_capping_low_percentile", self.feature_quantile_capping_low_percentile.value)
        )
        self.feature_quantile_capping_high_percentile.value = str(
            cfg.get("feature_quantile_capping_high_percentile", self.feature_quantile_capping_high_percentile.value)
        )
        saved_window_features = set(cfg.get("hmm_window_features", []))
        for feature_name, cb in self.additional_window_feature_cbs.items():
            cb.value = str(feature_name) in saved_window_features
        self._update_selected_features_box()
        self._rebuild_log_scale_feature_controls()
        self._sync_hmm_state_controls()

    def _persist_current_settings(self):
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        cfg.update({
            "selected_features": self._selected_feature_columns(),
            "binary_features_to_group": self._selected_binary_columns(),
            "random_state": int(self.random_state.value),
        })
        cfg.update(
            {
                "apply_hmm_deployment_artifact_path": str(self.apply_hmm_artifact_picker.value),
                "hmm_n_states_mode": str(self.hmm_n_states_mode.value),
                "hmm_advanced": bool(self.hmm_advanced.value),
                "hmm_n_states": int(self.hmm_n_states.value),
                "hmm_k_min": int(self.hmm_k_min.value),
                "hmm_k_max": int(self.hmm_k_max.value),
                "hmm_covariance_type": str(self.hmm_covariance_type.value),
                "hmm_n_iter": int(self.hmm_n_iter.value),
                "hmm_tol": float(self.hmm_tol.value),
                "hmm_sticky": bool(self.hmm_sticky.value),
                "hmm_stickiness_kappa": float(self.hmm_stickiness_kappa.value),
                "hmm_transmat_alpha": float(self.hmm_transmat_alpha.value),
                "hmm_min_covar": float(self.hmm_min_covar.value),
                "hmm_feature_smoothing_window": int(self.hmm_feature_smoothing_window.value),
                "hmm_smoothing_min_periods": int(self.hmm_smoothing_min_periods.value),
                "hmm_start_offset": int(self.hmm_start_offset.value),
                "hmm_start_offset_fill_mode": str(self.hmm_start_offset_fill_mode.value),
                "hmm_backprojection_workers": int(self.hmm_backprojection_workers.value),
                "hmm_window_features": self._selected_window_feature_columns(),
                "hmm_window_features_window": int(self.hmm_window_features_window.value),
                "feature_quantile_capping_low_percentile": str(self.feature_quantile_capping_low_percentile.value),
                "feature_quantile_capping_high_percentile": str(self.feature_quantile_capping_high_percentile.value),
                "hmm_log_scale_features": self._selected_log_scale_features(),
            }
        )
        self._save_panel_cfg()

    def _refresh_apply_default_paths(self):
        super()._refresh_apply_default_paths()
        if self.output_dir is None or str(self.output_dir).strip() == "":
            self.apply_hmm_default_paths_html.value = ""
            return

        hmm_path = self._default_hmm_deployment_artifact_path()
        hmm_exists = hmm_path.exists()
        if str(self.apply_hmm_artifact_picker.value).strip() == "" and hmm_exists:
            self.apply_hmm_artifact_picker.value = str(hmm_path)

        if hmm_exists:
            self.apply_hmm_default_paths_html.value = (
                "<b style='color:#080;'>Default HMM deployment artifact path detected and prefilled.</b>"
            )
        else:
            self.apply_hmm_default_paths_html.value = ""

    def _refresh_context(self):
        super()._refresh_context()
        if hasattr(self.apply_hmm_artifact_picker, "_start_dir"):
            self.apply_hmm_artifact_picker._start_dir = self.output_dir or "."
        self._load_existing_hmm_deployment_artifact_if_available()
        self._update_time_conversion_labels()

    def _on_time_conversion_widget_changed(self, change):
        self._update_time_conversion_labels()

    def _update_time_conversion_labels(self):
        metadata = getattr(self.metadata_loader, "metadata", None)
        self.hmm_window_features_time_html.value = _format_timepoint_window_as_time(
            self.hmm_window_features_window.value, metadata
        )
        self.hmm_feature_smoothing_time_html.value = _format_timepoint_window_as_time(
            self.hmm_feature_smoothing_window.value, metadata
        )

    def _refresh_enablement(self):
        super()._refresh_enablement()
        has_cell_type = self._current_cell_type() != ""
        has_features = len(self._selected_hmm_feature_columns()) > 0
        has_model = self.model_adata is not None
        has_intrinsic = has_model and (INTRINSIC_STATE_COL in self.model_adata.obs.columns)
        has_full = has_model and (FULL_STATE_COL in self.model_adata.obs.columns)
        has_backproj_sample = self.backproj_sample_dd.value is not None and len(str(self.backproj_sample_dd.value)) > 0
        has_hmm_artifact_input = str(self.apply_hmm_artifact_picker.value).strip() != ""
        self.btn_cluster.disabled = not (has_cell_type and has_features)
        self.btn_rename_intrinsic.disabled = not has_intrinsic
        self.btn_combine_intrinsic.disabled = not has_intrinsic
        self.intrinsic_combine_name.disabled = not has_intrinsic
        self.btn_rename_full.disabled = not has_full
        self.btn_combine_full.disabled = not has_full
        self.full_combine_name.disabled = not has_full
        self.btn_open_intrinsic_backprojection.disabled = not (has_cell_type and has_backproj_sample and has_intrinsic)
        self.btn_apply_hmm_artifact.disabled = not (has_cell_type and has_hmm_artifact_input)
        self.btn_create_state_composition_plots.disabled = not has_full
        self.btn_create_state_transition_plots.disabled = not has_full
        self._refresh_analysis_plots_status()

    def _save_current_hmm_deployment_artifact(self, verbose=True):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")
        self._ensure_full_state_colors(write=True)
        state_paths = _resolve_state_paths(self.output_dir, self._current_cell_type())
        artifact_path = self._default_hmm_deployment_artifact_path()
        current_artifact = self._hmm_deployment_artifact
        if current_artifact is None and artifact_path.exists():
            current_artifact = load_hmm_deployment_artifact(artifact_path)
        if current_artifact is None:
            raise ValueError(
                "No fitted HMM deployment artifact is available in memory or on disk for updating. "
                "Run intrinsic HMM clustering first."
            )

        saved_artifact = save_hmm_deployment_artifact(
            output_path=artifact_path,
            model_adata=self.model_adata,
            hmm_model=current_artifact["model"],
            state_paths=state_paths,
            source_model_adata_path=self._model_adata_path(),
            verbose=verbose,
        )
        self._hmm_deployment_artifact = saved_artifact
        if str(self.apply_hmm_artifact_picker.value).strip() == "":
            self.apply_hmm_artifact_picker.value = str(artifact_path)
        self._refresh_apply_default_paths()
        return artifact_path

    def _current_hmm_model_for_quality_control(self):
        artifact = self._hmm_deployment_artifact
        artifact_path = self._default_hmm_deployment_artifact_path()
        if artifact is None and artifact_path.exists():
            artifact = load_hmm_deployment_artifact(artifact_path)
            self._hmm_deployment_artifact = artifact
        if isinstance(artifact, dict):
            return artifact.get("model", None)
        return None

    def _regenerate_hmm_quality_control_plots(self, *, verbose=False):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")
        if INTRINSIC_STATE_COL not in self.model_adata.obs.columns:
            raise ValueError(f"model_adata is missing '{INTRINSIC_STATE_COL}'.")

        state_paths = _resolve_state_paths(self.output_dir, self._current_cell_type())
        qc_dir = Path(state_paths.processing_outdir) / "hmm_behavioral_classification" / "quality_control"
        preprocessing_meta = self.model_adata.uns.get("preprocessing", {})
        if not isinstance(preprocessing_meta, dict):
            preprocessing_meta = {}
        feature_cols = preprocessing_meta.get(
            "continuous_feature_cols",
            preprocessing_meta.get("kept_features", list(self.model_adata.var_names)),
        )
        feature_cols = [str(col) for col in list(feature_cols or []) if str(col) in self.model_adata.var_names]
        scaler_meta = preprocessing_meta.get("scaler", {}) if isinstance(preprocessing_meta, dict) else {}
        hmm_model = self._current_hmm_model_for_quality_control()
        qc_out = save_hmm_quality_control_outputs(
            self.model_adata,
            feature_cols=feature_cols,
            output_dir=qc_dir,
            model=hmm_model,
            selection_df=None,
            cluster_col=INTRINSIC_STATE_COL,
            scaler_mean=scaler_meta.get("mean", None) if isinstance(scaler_meta, dict) else None,
            scaler_scale=scaler_meta.get("scale", None) if isinstance(scaler_meta, dict) else None,
            title=f"all_data | hmm | curated {INTRINSIC_STATE_COL}",
            preprocessing_params=preprocessing_meta,
            verbose=verbose,
        )
        clustering_meta = self.model_adata.uns.get("clustering", {})
        if not isinstance(clustering_meta, dict):
            clustering_meta = {}
        clustering_meta["quality_control_dir"] = str(qc_dir)
        clustering_meta["raw_quality_control_dir"] = str(qc_dir / "raw")
        clustering_meta["diagnostics_pdf"] = qc_out.get("diagnostics_pdf", None)
        clustering_meta["diagnostics_csvs"] = dict(qc_out.get("diagnostics_csvs", {}) or {})
        clustering_meta["feature_distribution_pdf"] = qc_out.get("feature_distribution_pdf", None)
        clustering_meta["state_counts_csv"] = qc_out.get("state_counts_csv", None)
        clustering_meta["transition_matrix_csv"] = qc_out.get("transition_matrix_csv", None)
        clustering_meta["model_selection_csv"] = qc_out.get("model_selection_csv", None)
        self.model_adata.uns["clustering"] = clustering_meta
        return qc_out

    def _rename_mapping_yaml_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type).strip()
        state_paths = _resolve_state_paths(self.output_dir, ct)
        return (
            Path(state_paths.processing_outdir)
            / "hmm_behavioral_classification"
            / f"hmm_cluster_name_mappings_{ct}.yml"
        )

    def _full_state_labels(self):
        if self.model_adata is None or FULL_STATE_COL not in getattr(self.model_adata, "obs", {}).columns:
            return []
        labels = pd.Series(self.model_adata.obs[FULL_STATE_COL]).dropna().astype("string").str.strip()
        labels = labels[labels != ""]
        return sorted([str(v) for v in labels.unique().tolist()], key=_mixed_label_sort_key)

    def _current_full_state_color_mapping(self, labels=None, *, prefer_pickers=True, write=False):
        labels = self._full_state_labels() if labels is None else [str(v) for v in list(labels or [])]
        saved_colors = _get_classification_state_colors(self.model_adata, FULL_STATE_COL)
        if bool(prefer_pickers):
            for label, picker in getattr(self, "_full_color_pickers", {}).items():
                if label in labels:
                    saved_colors[str(label)] = _coerce_hex_color(getattr(picker, "value", None))
        colors = _normalize_label_color_map(labels, colors=saved_colors)
        if bool(write) and len(colors) > 0:
            _set_classification_state_colors(self.model_adata, FULL_STATE_COL, colors)
        return colors

    def _ensure_full_state_colors(self, labels=None, write=False):
        return self._current_full_state_color_mapping(
            labels=labels,
            prefer_pickers=False,
            write=write,
        )

    def _sync_and_save_full_adata(self, verbose=True):
        """Push the current color/order mapping from `self.model_adata.uns` into
        `self.adata_full.uns`, and persist `self.adata_full` to the canonical full-dataset
        path. No-ops gracefully if `self.adata_full` has not been populated yet."""
        adata_full = getattr(self, "adata_full", None)
        if adata_full is None or self.model_adata is None:
            return None

        state_paths = _resolve_state_paths(self.output_dir, self._current_cell_type())
        synced_any = False
        for state_col in (FULL_STATE_COL, INTRINSIC_STATE_COL):
            if state_col not in getattr(adata_full, "obs", {}).columns:
                continue
            if state_col not in getattr(self.model_adata, "obs", {}).columns:
                continue
            colors = _get_classification_state_colors(self.model_adata, state_col)
            order = _get_classification_state_order(self.model_adata, state_col)
            if colors:
                _set_classification_state_colors(adata_full, state_col, colors)
            if order:
                _set_classification_state_order(adata_full, state_col, order)
            synced_any = True

        if not synced_any:
            return None

        state_paths.full_output_adata_path.parent.mkdir(parents=True, exist_ok=True)
        adata_full.write(state_paths.full_output_adata_path, compression="gzip")
        if verbose:
            _winfo(
                "state-hmm-widget",
                f"Synced colors/order and saved full dataset: {state_paths.full_output_adata_path}",
            )
        return state_paths.full_output_adata_path

    def _remap_full_state_colors(self, mapping, existing_labels):
        old_labels = sorted([str(v) for v in list(existing_labels or [])], key=_mixed_label_sort_key)
        old_colors = self._current_full_state_color_mapping(labels=old_labels, prefer_pickers=True)
        new_colors = {}
        for old_label in old_labels:
            new_label = str(mapping.get(old_label, old_label)).strip() or old_label
            if new_label not in new_colors:
                new_colors[new_label] = old_colors.get(old_label)
        return _normalize_label_color_map(
            sorted(new_colors.keys(), key=_mixed_label_sort_key),
            colors=new_colors,
        )

    def _mapping_dict_from_obs(self, source_col, target_col):
        if self.model_adata is None or not hasattr(self.model_adata, "obs"):
            return {}
        if source_col not in self.model_adata.obs.columns or target_col not in self.model_adata.obs.columns:
            return {}

        mapping_obs = self.model_adata.obs[[source_col, target_col]].copy()
        mapping_obs[source_col] = mapping_obs[source_col].astype("string").str.strip().fillna("")
        mapping_obs[target_col] = mapping_obs[target_col].astype("string").str.strip().fillna("")
        mapping_obs = mapping_obs[(mapping_obs[source_col] != "") & (mapping_obs[target_col] != "")]
        if len(mapping_obs) == 0:
            return {}

        grouped = mapping_obs.groupby(source_col, observed=False)[target_col].agg(
            lambda s: sorted({str(v) for v in s if str(v) != ""}, key=_mixed_label_sort_key)
        )
        mapping = {}
        for source_value, targets in grouped.items():
            source_key = str(source_value).strip()
            if source_key == "":
                continue
            mapping[source_key] = targets[0] if len(targets) == 1 else list(targets)
        return mapping

    def _write_cluster_name_mappings_yaml(self, *, latest_intrinsic_mapping=None, latest_full_mapping=None):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")

        current_intrinsic_mapping = self._mapping_dict_from_obs(
            HMM_INTRINSIC_RAW_STATE_COL,
            INTRINSIC_STATE_COL,
        )
        current_full_mapping = {}
        current_full_mapping_from_intrinsic = {}
        if all(
            col in self.model_adata.obs.columns
            for col in [BINARY_GROUP_COL, HMM_INTRINSIC_RAW_STATE_COL, FULL_STATE_COL]
        ):
            full_obs = self.model_adata.obs[
                [BINARY_GROUP_COL, HMM_INTRINSIC_RAW_STATE_COL, FULL_STATE_COL]
            ].copy()
            full_obs[BINARY_GROUP_COL] = full_obs[BINARY_GROUP_COL].astype("string").str.strip().fillna("")
            full_obs[HMM_INTRINSIC_RAW_STATE_COL] = (
                pd.to_numeric(full_obs[HMM_INTRINSIC_RAW_STATE_COL], errors="coerce")
                .astype("Int64")
                .astype("string")
                .str.strip()
                .fillna("")
            )
            full_obs[FULL_STATE_COL] = full_obs[FULL_STATE_COL].astype("string").str.strip().fillna("")
            full_obs = full_obs[
                (full_obs[BINARY_GROUP_COL] != "")
                & (full_obs[HMM_INTRINSIC_RAW_STATE_COL] != "")
                & (full_obs[FULL_STATE_COL] != "")
            ].copy()
            if len(full_obs) > 0:
                full_obs["_generated_full_label"] = (
                    full_obs[BINARY_GROUP_COL].astype(str)
                    + "_"
                    + full_obs[HMM_INTRINSIC_RAW_STATE_COL].astype(str)
                )
                grouped = full_obs.groupby("_generated_full_label", observed=False)[FULL_STATE_COL].agg(
                    lambda s: sorted({str(v) for v in s if str(v) != ""})
                )
                for key, values in grouped.items():
                    current_full_mapping[str(key)] = values[0] if len(values) == 1 else list(values)

        if all(
            col in self.model_adata.obs.columns
            for col in [BINARY_GROUP_COL, INTRINSIC_STATE_COL, FULL_STATE_COL]
        ):
            intrinsic_obs = self.model_adata.obs[
                [BINARY_GROUP_COL, INTRINSIC_STATE_COL, FULL_STATE_COL]
            ].copy()
            intrinsic_obs[BINARY_GROUP_COL] = intrinsic_obs[BINARY_GROUP_COL].astype("string").str.strip().fillna("")
            intrinsic_obs[INTRINSIC_STATE_COL] = (
                intrinsic_obs[INTRINSIC_STATE_COL].astype("string").str.strip().fillna("")
            )
            intrinsic_obs[FULL_STATE_COL] = intrinsic_obs[FULL_STATE_COL].astype("string").str.strip().fillna("")
            intrinsic_obs = intrinsic_obs[
                (intrinsic_obs[BINARY_GROUP_COL] != "")
                & (intrinsic_obs[INTRINSIC_STATE_COL] != "")
                & (intrinsic_obs[FULL_STATE_COL] != "")
            ].copy()
            if len(intrinsic_obs) > 0:
                intrinsic_obs["_generated_full_label"] = (
                    intrinsic_obs[BINARY_GROUP_COL].astype(str)
                    + "_"
                    + intrinsic_obs[INTRINSIC_STATE_COL].astype(str)
                )
                grouped = intrinsic_obs.groupby("_generated_full_label", observed=False)[FULL_STATE_COL].agg(
                    lambda s: sorted({str(v) for v in s if str(v) != ""})
                )
                for key, values in grouped.items():
                    current_full_mapping_from_intrinsic[str(key)] = (
                        values[0] if len(values) == 1 else list(values)
                    )

        payload = {
            "cell_type": str(self._current_cell_type()),
            "model_adata_path": str(self._model_adata_path()),
            "latest_intrinsic_rename_mapping": None if latest_intrinsic_mapping is None else dict(latest_intrinsic_mapping),
            "latest_full_rename_mapping": None if latest_full_mapping is None else dict(latest_full_mapping),
            "current_mappings": {
                f"{HMM_INTRINSIC_RAW_STATE_COL}_to_{INTRINSIC_STATE_COL}": dict(current_intrinsic_mapping),
                f"generated_{FULL_STATE_COL}_to_{FULL_STATE_COL}": dict(current_full_mapping),
                f"generated_{FULL_STATE_COL}_from_intrinsic_to_{FULL_STATE_COL}": dict(
                    current_full_mapping_from_intrinsic
                ),
            },
            "current_colors": {
                FULL_STATE_COL: dict(self._current_full_state_color_mapping(write=True)),
            },
        }
        yaml_path = self._rename_mapping_yaml_path()
        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        with yaml_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True)
        return yaml_path

    def _rebuild_intrinsic_rename_rows(self):
        self._intrinsic_name_boxes = {}
        self._intrinsic_select_boxes = {}
        self._intrinsic_row_widgets = {}
        self._intrinsic_row_order = []
        self.intrinsic_combine_name.value = ""
        if self.model_adata is None or INTRINSIC_STATE_COL not in self.model_adata.obs.columns:
            self.rename_intrinsic_rows.children = []
            load_status = getattr(self, "_hmm_model_load_status", {}) or {}
            status = str(load_status.get("status", "not_checked"))
            path = load_status.get("path", None)
            message = load_status.get("message", None)
            if status == "failed":
                self.rename_intrinsic_status.value = (
                    "<i>Could not load HMM model adata"
                    + (f" from {path}" if path else "")
                    + f": {message}</i>"
                )
            elif status == "missing_intrinsic_col":
                self.rename_intrinsic_status.value = f"<i>{message}</i>"
            elif status == "not_hmm":
                self.rename_intrinsic_status.value = f"<i>{message}</i>"
            elif status == "missing":
                self.rename_intrinsic_status.value = (
                    "<i>No saved HMM model adata found. Run HMM clustering first.</i>"
                )
            else:
                self.rename_intrinsic_status.value = "<i>Run HMM clustering or load existing model first.</i>"
            self.btn_rename_intrinsic.disabled = True
            self.btn_combine_intrinsic.disabled = True
            self.intrinsic_combine_name.disabled = True
            return

        mapping = build_identity_cluster_mapping(
            self.model_adata,
            cluster_col=INTRINSIC_STATE_COL,
        )
        saved_order = _get_classification_state_order(self.model_adata, INTRINSIC_STATE_COL)
        self._intrinsic_row_order = _apply_state_order(list(mapping.keys()), saved_order)
        for old_name in self._intrinsic_row_order:
            old_name = str(old_name)
            sel = widgets.Checkbox(
                value=False,
                description="",
                indent=False,
                layout=widgets.Layout(width="26px"),
            )
            txt = widgets.Text(value=old_name, layout=widgets.Layout(width="280px"))
            self._intrinsic_select_boxes[old_name] = sel
            self._intrinsic_name_boxes[old_name] = txt
            move_btns = build_row_move_buttons(
                on_move_up=lambda n=old_name: self._move_intrinsic_row(n, -1),
                on_move_down=lambda n=old_name: self._move_intrinsic_row(n, 1),
            )
            self._intrinsic_row_widgets[old_name] = widgets.HBox(
                [move_btns, sel, widgets.Label(old_name, layout=widgets.Layout(width="190px")), txt],
                layout=widgets.Layout(align_items="center", gap="8px"),
            )
        self._refresh_intrinsic_row_children()
        self.rename_intrinsic_status.value = f"<b>Intrinsic HMM states:</b> {len(self._intrinsic_row_order)}"
        self.btn_rename_intrinsic.disabled = False
        self.btn_combine_intrinsic.disabled = False
        self.intrinsic_combine_name.disabled = False

    def _refresh_intrinsic_row_children(self):
        self.rename_intrinsic_rows.children = [
            self._intrinsic_row_widgets[n] for n in self._intrinsic_row_order if n in self._intrinsic_row_widgets
        ]

    def _move_intrinsic_row(self, old_name, delta):
        order = self._intrinsic_row_order
        old_name = str(old_name)
        if old_name not in order:
            return
        idx = order.index(old_name)
        new_idx = idx + delta
        if not (0 <= new_idx < len(order)):
            return
        order[idx], order[new_idx] = order[new_idx], order[idx]
        self._refresh_intrinsic_row_children()

    def _rebuild_full_rename_rows(self):
        self._full_select_boxes = {}
        self._full_name_boxes = {}
        self._full_color_pickers = {}
        self._full_row_widgets = {}
        self._full_row_order = []
        self.rename_full_rows.layout = widgets.Layout(width="660px")
        self.full_combine_name.value = ""
        if self.model_adata is None or FULL_STATE_COL not in self.model_adata.obs.columns:
            self.rename_full_rows.children = []
            self.rename_full_status.value = "<i>Run intrinsic HMM clustering and rename intrinsic states first.</i>"
            self.btn_rename_full.disabled = True
            self.btn_combine_full.disabled = True
            self.full_combine_name.disabled = True
            return

        mapping = build_identity_cluster_mapping(
            self.model_adata,
            cluster_col=FULL_STATE_COL,
        )
        state_colors = self._ensure_full_state_colors(labels=list(mapping.keys()), write=True)
        saved_order = _get_classification_state_order(self.model_adata, FULL_STATE_COL)
        self._full_row_order = _apply_state_order(list(mapping.keys()), saved_order)
        for old_name in self._full_row_order:
            old_name = str(old_name)
            sel = widgets.Checkbox(
                value=False,
                description="",
                indent=False,
                layout=widgets.Layout(width="26px"),
            )
            txt = widgets.Text(value=old_name, layout=widgets.Layout(width="280px"))
            color = widgets.ColorPicker(
                value=str(state_colors.get(old_name, "#808080")),
                concise=True,
                layout=widgets.Layout(width="70px"),
            )
            self._full_select_boxes[old_name] = sel
            self._full_name_boxes[old_name] = txt
            self._full_color_pickers[old_name] = color
            move_btns = build_row_move_buttons(
                on_move_up=lambda n=old_name: self._move_full_row(n, -1),
                on_move_down=lambda n=old_name: self._move_full_row(n, 1),
            )
            self._full_row_widgets[old_name] = widgets.HBox(
                [move_btns, sel, widgets.Label(old_name, layout=widgets.Layout(width="190px")), txt, color],
                layout=widgets.Layout(align_items="center", gap="8px"),
            )
        self._refresh_full_row_children()
        self.rename_full_status.value = (
            f"<b>HMM states assigned to binary groups:</b> {len(self._full_row_order)} "
            "(observed combinations from training data)"
        )
        self.btn_rename_full.disabled = False
        self.btn_combine_full.disabled = False
        self.full_combine_name.disabled = False

    def _refresh_full_row_children(self):
        self.rename_full_rows.children = [
            self._full_row_widgets[n] for n in self._full_row_order if n in self._full_row_widgets
        ]

    def _move_full_row(self, old_name, delta):
        order = self._full_row_order
        old_name = str(old_name)
        if old_name not in order:
            return
        idx = order.index(old_name)
        new_idx = idx + delta
        if not (0 <= new_idx < len(order)):
            return
        order[idx], order[new_idx] = order[new_idx], order[idx]
        self._refresh_full_row_children()

    def _apply_intrinsic_rename_mapping(self, mapping, save_compression="lzf"):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")

        normalized_mapping = {}
        for old_name, new_name in mapping.items():
            old_s = str(old_name)
            new_s = str(new_name).strip()
            normalized_mapping[old_s] = new_s if new_s != "" else old_s

        existing_labels = {
            str(x)
            for x in self.model_adata.obs.get(INTRINSIC_STATE_COL, pd.Series(dtype="object")).astype(str)
        }
        has_changes = any(
            str(normalized_mapping.get(label, label)) != str(label)
            for label in existing_labels
        )
        new_intrinsic_order = list(dict.fromkeys(
            normalized_mapping.get(n, n) for n in getattr(self, "_intrinsic_row_order", [])
        ))
        order_changes = list(new_intrinsic_order) != list(
            _get_classification_state_order(self.model_adata, INTRINSIC_STATE_COL)
        )
        if not has_changes:
            if not order_changes:
                return {"changed": False}
            adata_full = getattr(self, "adata_full", None)
            _set_classification_state_order(self.model_adata, INTRINSIC_STATE_COL, new_intrinsic_order)
            if adata_full is not None and INTRINSIC_STATE_COL in getattr(adata_full, "obs", {}).columns:
                _set_classification_state_order(adata_full, INTRINSIC_STATE_COL, new_intrinsic_order)
            self._save_model_adata(compression=save_compression)
            self._sync_and_save_full_adata(verbose=False)
            yaml_path = self._write_cluster_name_mappings_yaml(latest_intrinsic_mapping=normalized_mapping)
            self._rebuild_intrinsic_rename_rows()
            self._rebuild_full_rename_rows()
            return {"changed": True, "order_changed": True, "mapping_yaml_path": str(yaml_path)}

        relabel_cluster_ids(
            adata=self.model_adata,
            mapping=normalized_mapping,
            cluster_key=INTRINSIC_STATE_COL,
            new_key=INTRINSIC_STATE_COL,
            keep_unmapped=True,
            overwrite_original=True,
            categories=new_intrinsic_order,
        )
        _set_classification_state_order(self.model_adata, INTRINSIC_STATE_COL, new_intrinsic_order)
        clustering_meta = self.model_adata.uns.get("clustering", {})
        binary_group_constraints = None
        enforce_binary_group_constraints = False
        if isinstance(clustering_meta, dict):
            binary_group_constraints = clustering_meta.get("binary_group_constraints", None)
            enforce_binary_group_constraints = (
                isinstance(binary_group_constraints, dict)
                and ("forbidden_binary_combinations" in binary_group_constraints)
            )
        binary_cols_to_merge = self._selected_binary_columns()
        if len(binary_cols_to_merge) == 0 and isinstance(clustering_meta, dict):
            binary_cols_to_merge = [str(c) for c in list(clustering_meta.get("binary_cols_to_merge", []) or [])]
        old_full_labels = (
            self.model_adata.obs[FULL_STATE_COL].astype(str).copy()
            if FULL_STATE_COL in self.model_adata.obs.columns
            else None
        )
        _rebuild_full_behavioral_cluster_from_intrinsic(
            adata=self.model_adata,
            binary_cols_to_merge=binary_cols_to_merge,
            intrinsic_col=INTRINSIC_STATE_COL,
            binary_group_constraints=binary_group_constraints,
            enforce_binary_group_constraints=enforce_binary_group_constraints,
        )
        self.model_adata.obs[FULL_STATE_COL] = pd.Categorical(
            pd.Series(
                self.model_adata.obs["full_behavioral_cluster"],
                index=self.model_adata.obs.index,
                dtype="string",
            )
        )
        adata_full = getattr(self, "adata_full", None)
        if adata_full is not None and INTRINSIC_STATE_COL in getattr(adata_full, "obs", {}).columns:
            relabel_cluster_ids(
                adata=adata_full,
                mapping=normalized_mapping,
                cluster_key=INTRINSIC_STATE_COL,
                new_key=INTRINSIC_STATE_COL,
                keep_unmapped=True,
                overwrite_original=True,
                categories=new_intrinsic_order,
            )
        if (
            adata_full is not None
            and old_full_labels is not None
            and FULL_STATE_COL in getattr(adata_full, "obs", {}).columns
        ):
            new_full_labels = self.model_adata.obs[FULL_STATE_COL].astype(str)
            derived_full_mapping = dict(zip(old_full_labels.tolist(), new_full_labels.tolist()))
            relabel_cluster_ids(
                adata=adata_full,
                mapping=derived_full_mapping,
                cluster_key=FULL_STATE_COL,
                overwrite_original=True,
                keep_unmapped=True,
            )
        self._ensure_full_state_colors(write=True)
        self._invalidate_curated_state_reports(
            reason="State reports were cleared after curated HMM renaming. Recreate them from 'Create analysis plots'."
        )
        deployment_warning = None
        try:
            artifact_path = self._save_current_hmm_deployment_artifact(verbose=False)
            _winfo("state-hmm-widget", f"Updated HMM deployment artifact: {artifact_path}")
        except Exception as exc:
            deployment_warning = str(exc)
        qc_warning = None
        qc_out = {}
        try:
            qc_out = self._regenerate_hmm_quality_control_plots(verbose=False)
            _winfo(
                "state-hmm-widget",
                f"Recreated curated HMM quality-control plots: {qc_out.get('diagnostics_pdf', None)}",
            )
        except Exception as exc:
            qc_warning = str(exc)
        self._save_model_adata(compression=save_compression)
        self._sync_and_save_full_adata(verbose=False)
        yaml_path = self._write_cluster_name_mappings_yaml(latest_intrinsic_mapping=normalized_mapping)
        self._rebuild_intrinsic_rename_rows()
        self._rebuild_full_rename_rows()
        return {
            "changed": True,
            "deployment_warning": deployment_warning,
            "quality_control_warning": qc_warning,
            "quality_control_outputs": dict(qc_out),
            "mapping_yaml_path": str(yaml_path),
        }

    def _apply_full_rename_mapping(self, mapping, save_compression="lzf"):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")

        normalized_mapping = {}
        for old_name, new_name in mapping.items():
            old_s = str(old_name)
            new_s = str(new_name).strip()
            normalized_mapping[old_s] = new_s if new_s != "" else old_s

        existing_labels = {
            str(x)
            for x in self.model_adata.obs.get(FULL_STATE_COL, pd.Series(dtype="object")).astype(str)
        }
        existing_labels = {label for label in existing_labels if label.strip() != ""}
        saved_colors = self._current_full_state_color_mapping(
            labels=sorted(existing_labels, key=_mixed_label_sort_key),
            prefer_pickers=False,
        )
        remapped_colors = self._remap_full_state_colors(normalized_mapping, existing_labels)
        new_full_order = list(dict.fromkeys(
            normalized_mapping.get(n, n) for n in getattr(self, "_full_row_order", [])
        ))
        has_changes = any(
            str(normalized_mapping.get(label, label)) != str(label)
            for label in existing_labels
        )
        color_changes = dict(saved_colors) != dict(remapped_colors)
        order_changes = list(new_full_order) != list(_get_classification_state_order(self.model_adata, FULL_STATE_COL))
        adata_full = getattr(self, "adata_full", None)
        adata_full_has_full_state = (
            adata_full is not None and FULL_STATE_COL in getattr(adata_full, "obs", {}).columns
        )
        if not has_changes:
            if color_changes or order_changes:
                _set_classification_state_colors(self.model_adata, FULL_STATE_COL, remapped_colors)
                _set_classification_state_order(self.model_adata, FULL_STATE_COL, new_full_order)
                if adata_full_has_full_state:
                    _set_classification_state_colors(adata_full, FULL_STATE_COL, remapped_colors)
                    _set_classification_state_order(adata_full, FULL_STATE_COL, new_full_order)
                yaml_path = self._write_cluster_name_mappings_yaml(latest_full_mapping=normalized_mapping)
                self._rebuild_full_rename_rows()
                result = {
                    "changed": True,
                    "colors_changed": bool(color_changes),
                    "order_changed": bool(order_changes),
                    "mapping_yaml_path": str(yaml_path),
                }
            else:
                result = {"changed": False}
        else:
            relabel_cluster_ids(
                adata=self.model_adata,
                mapping=normalized_mapping,
                cluster_key=FULL_STATE_COL,
                overwrite_original=True,
                keep_unmapped=True,
                categories=new_full_order,
            )
            if adata_full_has_full_state:
                relabel_cluster_ids(
                    adata=adata_full,
                    mapping=normalized_mapping,
                    cluster_key=FULL_STATE_COL,
                    overwrite_original=True,
                    keep_unmapped=True,
                    categories=new_full_order,
                )
                _set_classification_state_colors(adata_full, FULL_STATE_COL, remapped_colors)
                _set_classification_state_order(adata_full, FULL_STATE_COL, new_full_order)
            _set_classification_state_colors(self.model_adata, FULL_STATE_COL, remapped_colors)
            _set_classification_state_order(self.model_adata, FULL_STATE_COL, new_full_order)
            yaml_path = self._write_cluster_name_mappings_yaml(latest_full_mapping=normalized_mapping)
            self._rebuild_full_rename_rows()
            result = {"changed": True, "colors_changed": bool(color_changes), "mapping_yaml_path": str(yaml_path)}
        if bool(result.get("changed", False)):
            self._invalidate_curated_state_reports(
                reason="State reports were cleared after curated HMM renaming. Recreate them from 'Create analysis plots'."
            )
            deployment_warning = None
            try:
                artifact_path = self._save_current_hmm_deployment_artifact(verbose=False)
                _winfo(
                    "state-hmm-widget",
                    f"Full classification pipeline and settings saved to:\n  {artifact_path}\n"
                    "You can apply this directly to new data by ticking "
                    "'Apply existing behavioral state classification'.",
                )
            except Exception as exc:
                deployment_warning = str(exc)
            qc_warning = None
            qc_out = {}
            try:
                qc_out = self._regenerate_hmm_quality_control_plots(verbose=False)
                _winfo(
                    "state-hmm-widget",
                    f"Recreated curated HMM quality-control plots: {qc_out.get('diagnostics_pdf', None)}",
                )
            except Exception as exc:
                qc_warning = str(exc)
            self._save_model_adata(compression=save_compression)
            self._sync_and_save_full_adata(verbose=False)
            result["deployment_warning"] = deployment_warning
            result["quality_control_warning"] = qc_warning
            result["quality_control_outputs"] = dict(qc_out)
        return result

    def _invalidate_curated_state_reports(self, reason=None):
        if self.model_adata is None:
            return
        _update_hmm_report_metadata(
            self.model_adata,
            composition_pdf=None,
            composition_auc_csv=None,
            composition_plot_csvs=[],
            composition_error=None,
            transition_dir=None,
            transition_counts_csv=None,
            transition_probs_csv=None,
            transition_error=None,
            pending_reason=reason,
        )

    _METADATA_TECHNICAL_PATTERNS = (
        "_path", "_dir", "pixel_distance", "time_interval", "unit_",
        "channel_", "dead_channel", "zarr", "dimension_order",
    )

    def _merge_metadata_into_obs(self, adata):
        if adata is None or not hasattr(adata, "obs"):
            return adata
        md = getattr(self.metadata_loader, "metadata", None)
        if md is None or "sample_name" not in getattr(md, "columns", []):
            return adata
        meta_cols = self._metadata_grouping_columns()
        cols_to_add = [c for c in meta_cols if c not in adata.obs.columns and c in md.columns]
        if not cols_to_add:
            return adata
        meta_subset = md[["sample_name"] + cols_to_add].drop_duplicates(subset=["sample_name"])
        orig_index = adata.obs.index
        merged = adata.obs.merge(meta_subset, on="sample_name", how="left")
        merged.index = orig_index
        adata.obs = merged
        return adata

    def _metadata_grouping_columns(self):
        md = getattr(self.metadata_loader, "metadata", None)
        if md is None or not hasattr(md, "columns"):
            return []
        exclude = self._METADATA_TECHNICAL_PATTERNS
        cols = [
            c for c in md.columns
            if not any(pat in c.lower() for pat in exclude)
        ]
        return cols

    def _sync_comparison_group_y_text(self):
        if hasattr(self, "comparison_group_y_text"):
            self.comparison_group_y_text.value = str(self.comparison_condition_col_dd.value or "")

    def _on_comparison_condition_col_changed(self, change):
        if change.get("name") == "value":
            self._sync_comparison_group_y_text()

    def _refresh_analysis_plots_status(self):
        if not hasattr(self, "analysis_plots_status"):
            return
        if self.model_adata is None or FULL_STATE_COL not in getattr(self.model_adata, "obs", {}).columns:
            self.analysis_plots_status.value = (
                "<i>Run HMM clustering and finish curated renaming first, then create plots manually.</i>"
            )
            return

        clustering_meta = self.model_adata.uns.get("clustering", {})
        if not isinstance(clustering_meta, dict):
            clustering_meta = {}
        has_composition = clustering_meta.get("state_composition_report_pdf", None) is not None
        has_transition = clustering_meta.get("state_transition_report_dir", None) is not None
        pending_reason = clustering_meta.get("state_reports_reason", None)

        if has_composition and has_transition:
            self.analysis_plots_status.value = (
                "<b>Analysis plots ready:</b> state composition and state transition outputs have been created."
            )
        elif has_composition or has_transition:
            ready = []
            if has_composition:
                ready.append("state composition")
            if has_transition:
                ready.append("state transition")
            self.analysis_plots_status.value = (
                "<b>Analysis plots partially ready:</b> created "
                + ", ".join(ready)
                + "."
            )
        else:
            self.analysis_plots_status.value = (
                f"<i>{pending_reason}</i>"
                if pending_reason is not None and len(str(pending_reason).strip()) > 0
                else "<i>Analysis plots have not been created yet.</i>"
            )

        if hasattr(self, "composition_group_cols_select"):
            candidate_cols = self._metadata_grouping_columns()
            if candidate_cols != list(self.composition_group_cols_select.options):
                prev = set(self.composition_group_cols_select.value)
                self.composition_group_cols_select.options = candidate_cols
                self.composition_group_cols_select.value = [c for c in candidate_cols if c in prev]
            self.composition_group_cols_select.rows = max(2, min(len(candidate_cols), 6))
            self.composition_group_cols_select.disabled = len(candidate_cols) == 0

            axis_options = ["(none)"] + candidate_cols
            for dd in (self.composition_group_x_dd, self.composition_group_y_dd):
                if axis_options != list(dd.options):
                    prev_axis = dd.value
                    dd.options = axis_options
                    dd.value = prev_axis if prev_axis in axis_options else "(none)"
                dd.disabled = len(candidate_cols) == 0

        if hasattr(self, "comparison_condition_col_dd"):
            candidate_cols = self._metadata_grouping_columns()
            if candidate_cols != list(self.comparison_condition_col_dd.options):
                prev_cond = self.comparison_condition_col_dd.value
                self.comparison_condition_col_dd.options = candidate_cols
                if prev_cond in candidate_cols:
                    self.comparison_condition_col_dd.value = prev_cond
                elif candidate_cols:
                    self.comparison_condition_col_dd.value = candidate_cols[0]
            self.comparison_condition_col_dd.disabled = len(candidate_cols) == 0

            if candidate_cols != list(self.comparison_group_cols_select.options):
                prev = set(self.comparison_group_cols_select.value)
                self.comparison_group_cols_select.options = candidate_cols
                self.comparison_group_cols_select.value = [c for c in candidate_cols if c in prev]
            self.comparison_group_cols_select.rows = max(2, min(len(candidate_cols), 6))
            self.comparison_group_cols_select.disabled = len(candidate_cols) == 0

            axis_options = ["(none)"] + candidate_cols
            if axis_options != list(self.comparison_group_x_dd.options):
                prev_axis = self.comparison_group_x_dd.value
                self.comparison_group_x_dd.options = axis_options
                self.comparison_group_x_dd.value = prev_axis if prev_axis in axis_options else "(none)"
            self.comparison_group_x_dd.disabled = len(candidate_cols) == 0
            self._sync_comparison_group_y_text()

            self.btn_create_condition_comparison.disabled = len(candidate_cols) == 0

    def _regenerate_curated_state_reports(
        self,
        verbose=False,
        *,
        create_composition=True,
        create_transition=True,
        group_cols=None,
        group_x=None,
        group_y=None,
    ):
        if self.model_adata is None:
            raise ValueError("No model adata loaded.")
        if FULL_STATE_COL not in self.model_adata.obs.columns:
            raise ValueError(f"model_adata is missing '{FULL_STATE_COL}'.")

        # Use the full applied dataset when available — it covers all samples and timepoints.
        # Fall back to the training model adata when no artifact has been applied yet.
        adata_for_plots = (
            self.adata_full
            if (
                getattr(self, "adata_full", None) is not None
                and FULL_STATE_COL in getattr(self.adata_full, "obs", {}).columns
            )
            else self.model_adata
        )

        state_paths = _resolve_state_paths(self.output_dir, self._current_cell_type())
        clustering_meta = self.model_adata.uns.get("clustering", {})
        if not isinstance(clustering_meta, dict):
            clustering_meta = {}
        full_state_colors = self._current_full_state_color_mapping(write=True)
        full_state_order = _get_classification_state_order(adata_for_plots, FULL_STATE_COL)
        if adata_for_plots is not self.model_adata and FULL_STATE_COL in getattr(adata_for_plots, "obs", {}).columns:
            observed_labels = set(
                pd.Series(adata_for_plots.obs[FULL_STATE_COL]).dropna().astype(str).unique().tolist()
            )
            missing_colors = sorted(observed_labels - set(full_state_colors.keys()))
            if missing_colors:
                print(
                    f"  Warning: {len(missing_colors)} state label(s) in adata_full have no matching saved "
                    f"color (falling back to defaults): {missing_colors[:10]}"
                    f"{'...' if len(missing_colors) > 10 else ''}"
                )

        composition_pdf = clustering_meta.get("state_composition_report_pdf", None)
        composition_auc_csv = clustering_meta.get("state_composition_report_auc_csv", None)
        composition_plot_csvs = list(clustering_meta.get("state_composition_report_plot_csvs", []) or [])
        composition_error = clustering_meta.get("state_composition_report_error", None)
        if bool(create_composition):
            composition_pdf = None
            composition_auc_csv = None
            composition_plot_csvs = []
            composition_error = None
            try:
                composition_dir = Path(state_paths.state_composition_outdir)
                composition_dir.mkdir(parents=True, exist_ok=True)
                report_pdf_path = composition_dir / f"state_composition_report_{FULL_STATE_COL}.pdf"
                report_csv_path = composition_dir / f"state_composition_report_{FULL_STATE_COL}.csv"
                composition_out = save_state_composition_report(
                    adata=adata_for_plots,
                    output_pdf_path=report_pdf_path,
                    output_csv_path=report_csv_path,
                    time_col="position_t",
                    state_col=FULL_STATE_COL,
                    sample_col="sample_name",
                    include_pooled_summary=True,
                    state_colors=full_state_colors,
                    state_order=full_state_order,
                    verbose=verbose,
                    group_cols=group_cols if group_cols else None,
                    group_x=group_x,
                    group_y=group_y,
                )
                composition_pdf = str(composition_out.get("pdf_path", report_pdf_path))
                composition_auc_csv = str(composition_out.get("csv_path", report_csv_path))
                plot_csv_paths = composition_out.get("plot_data_csv_paths", {})
                if isinstance(plot_csv_paths, dict):
                    composition_plot_csvs = [str(v) for v in plot_csv_paths.values()]
            except Exception as exc:
                composition_error = str(exc)

        transition_dir = clustering_meta.get("state_transition_report_dir", None)
        transition_counts_csv = clustering_meta.get("state_transition_matrix_counts_csv", None)
        transition_probs_csv = clustering_meta.get("state_transition_matrix_probs_csv", None)
        transition_error = clustering_meta.get("state_transition_report_error", None)
        if bool(create_transition):
            transition_dir = None
            transition_counts_csv = None
            transition_probs_csv = None
            transition_error = None
            try:
                transitions_outdir = Path(state_paths.state_transitions_outdir)
                transitions_outdir.mkdir(parents=True, exist_ok=True)
                transition_out = save_state_transition_report(
                    adata=adata_for_plots,
                    output_dir=transitions_outdir,
                    state_col=FULL_STATE_COL,
                    id_cols=("sample_name", "TrackID"),
                    time_col="position_t",
                    state_colors=full_state_colors,
                    state_order=full_state_order,
                    verbose=verbose,
                )
                transition_dir = str(transition_out.get("output_dir", transitions_outdir))
                transition_counts_csv = transition_out.get("transition_matrix_counts_csv", None)
                transition_probs_csv = transition_out.get("transition_matrix_probs_csv", None)
            except Exception as exc:
                transition_error = str(exc)

        _update_hmm_report_metadata(
            self.model_adata,
            composition_pdf=composition_pdf,
            composition_auc_csv=composition_auc_csv,
            composition_plot_csvs=composition_plot_csvs,
            composition_error=composition_error,
            transition_dir=transition_dir,
            transition_counts_csv=transition_counts_csv,
            transition_probs_csv=transition_probs_csv,
            transition_error=transition_error,
            pending_reason="State reports have not been created yet. Use 'Create analysis plots' in the HMM widget.",
        )
        return {
            "composition_pdf": composition_pdf,
            "composition_auc_csv": composition_auc_csv,
            "composition_plot_csvs": list(composition_plot_csvs),
            "composition_error": composition_error,
            "transition_dir": transition_dir,
            "transition_counts_csv": transition_counts_csv,
            "transition_probs_csv": transition_probs_csv,
            "transition_error": transition_error,
        }

    def _on_create_state_composition_plots_clicked(self, _):
        self._set_busy(self.btn_create_state_composition_plots, self.state_composition_spinner, busy=True)
        self.out_analysis_plots.clear_output()
        with self.out_analysis_plots:
            try:
                group_cols = list(self.composition_group_cols_select.value) or None
                group_x = self.composition_group_x_dd.value
                group_x = None if group_x in (None, "(none)") else group_x
                group_y = self.composition_group_y_dd.value
                group_y = None if group_y in (None, "(none)") else group_y
                report_out = self._regenerate_curated_state_reports(
                    verbose=True,
                    create_composition=True,
                    create_transition=False,
                    group_cols=group_cols,
                    group_x=group_x,
                    group_y=group_y,
                )
                self._save_model_adata(compression="lzf")
                _winfo(
                    "state-hmm-widget",
                    "Created state composition plots: "
                    f"{report_out.get('composition_pdf', None)}",
                )
                if report_out.get("composition_error", None):
                    _winfo(
                        "state-hmm-widget",
                        f"State composition plot generation reported an error: {report_out['composition_error']}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_create_state_composition_plots, self.state_composition_spinner, busy=False)
                self._refresh_enablement()

    def _on_create_state_transition_plots_clicked(self, _):
        self._set_busy(self.btn_create_state_transition_plots, self.state_transition_spinner, busy=True)
        self.out_analysis_plots.clear_output()
        with self.out_analysis_plots:
            try:
                report_out = self._regenerate_curated_state_reports(
                    verbose=True,
                    create_composition=False,
                    create_transition=True,
                )
                self._save_model_adata(compression="lzf")
                _winfo(
                    "state-hmm-widget",
                    "Created state transition plots: "
                    f"{report_out.get('transition_dir', None)}",
                )
                if report_out.get("transition_error", None):
                    _winfo(
                        "state-hmm-widget",
                        f"State transition plot generation reported an error: {report_out['transition_error']}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_create_state_transition_plots, self.state_transition_spinner, busy=False)
                self._refresh_enablement()

    def _on_create_condition_comparison_clicked(self, _):
        self._set_busy(self.btn_create_condition_comparison, self.condition_comparison_spinner, busy=True)
        self.out_analysis_plots.clear_output()
        with self.out_analysis_plots:
            try:
                condition_col = self.comparison_condition_col_dd.value
                group_cols = list(self.comparison_group_cols_select.value) or None
                group_x = self.comparison_group_x_dd.value
                group_x = None if group_x in (None, "(none)") else group_x
                if not condition_col:
                    raise ValueError("Select a condition column to compare.")

                adata_for_plots = (
                    self.adata_full
                    if (
                        getattr(self, "adata_full", None) is not None
                        and FULL_STATE_COL in getattr(self.adata_full, "obs", {}).columns
                    )
                    else self.model_adata
                )
                if adata_for_plots is None or FULL_STATE_COL not in getattr(adata_for_plots, "obs", {}).columns:
                    raise ValueError("No model adata loaded.")

                full_state_colors = self._current_full_state_color_mapping(write=True)
                full_state_order = _get_classification_state_order(adata_for_plots, FULL_STATE_COL)
                state_paths = _resolve_state_paths(self.output_dir, self._current_cell_type())
                out_dir = Path(state_paths.state_composition_outdir) / "behavior_proportions"
                out_dir.mkdir(parents=True, exist_ok=True)
                cond_token = _sanitize_filename_token(condition_col, fallback="condition")
                out_pdf = out_dir / f"condition_comparison_{cond_token}.pdf"
                out_csv = out_pdf.with_suffix(".csv")

                result = save_state_condition_comparison_report(
                    adata=adata_for_plots,
                    output_pdf_path=out_pdf,
                    output_csv_path=out_csv,
                    state_col=FULL_STATE_COL,
                    sample_col="sample_name",
                    condition_col=condition_col,
                    group_cols=group_cols,
                    group_x=group_x,
                    state_colors=full_state_colors,
                    state_order=full_state_order,
                    verbose=True,
                )
                _winfo(
                    "state-hmm-widget",
                    f"Created condition comparison plot: {result.get('pdf_path', out_pdf)}",
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_create_condition_comparison, self.condition_comparison_spinner, busy=False)

    def _on_cluster_clicked(self, _):
        self._set_busy(self.btn_cluster, self.cluster_spinner, busy=True)
        self.out_cluster.clear_output()
        with self.out_cluster:
            try:
                self._persist_current_settings()
                selected_features = self._selected_feature_columns()
                selected_window_features = self._selected_window_feature_columns()
                selected_hmm_features = self._selected_hmm_feature_columns()
                if len(selected_hmm_features) == 0:
                    raise ValueError("Select at least one feature before HMM clustering.")
                n_states_value = (
                    "auto"
                    if str(self.hmm_n_states_mode.value) == "auto"
                    else int(self.hmm_n_states.value)
                )
                _winfo(
                    "state-hmm-widget",
                    (
                        "Running intrinsic HMM clustering with "
                        f"n_states={n_states_value}, features={selected_hmm_features}"
                    ),
                )
                cluster_out = run_hmm_state_clustering(
                    features=selected_features,
                    additional_window_features=selected_window_features,
                    binary_features_to_group=self._selected_binary_columns(),
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    n_states=n_states_value,
                    k_min=int(self.hmm_k_min.value),
                    k_max=int(self.hmm_k_max.value),
                    covariance_type=str(self.hmm_covariance_type.value),
                    n_iter=int(self.hmm_n_iter.value),
                    tol=float(self.hmm_tol.value),
                    sticky=bool(self.hmm_sticky.value),
                    stickiness_kappa=float(self.hmm_stickiness_kappa.value),
                    transmat_alpha=float(self.hmm_transmat_alpha.value),
                    min_covar=float(self.hmm_min_covar.value),
                    feature_smoothing_window=int(self.hmm_feature_smoothing_window.value),
                    smoothing_min_periods=int(self.hmm_smoothing_min_periods.value),
                    start_offset=int(self.hmm_start_offset.value),
                    start_offset_fill_mode=str(self.hmm_start_offset_fill_mode.value),
                    window_features_window=int(self.hmm_window_features_window.value),
                    lower_quantile_cap=self._parse_optional_float(self.feature_quantile_capping_low_percentile.value),
                    upper_quantile_cap=self._parse_optional_float(self.feature_quantile_capping_high_percentile.value),
                    log_scale_features=self._selected_log_scale_features(),
                    random_state=int(self.random_state.value),
                    return_details=True,
                    verbose=True,
                )
                self.model_adata = cluster_out["model_adata"]
                self.model_adata = self._merge_metadata_into_obs(self.model_adata)
                self._hmm_model_load_status = {
                    "status": "loaded",
                    "path": str(self._model_adata_path()),
                    "message": None,
                }
                self._hmm_deployment_artifact = cluster_out.get("hmm_deployment_artifact", None)
                artifact_path = self._save_current_hmm_deployment_artifact(verbose=True)
                self._save_model_adata()
                _winfo("state-hmm-widget", f"Saved model adata: {self._model_adata_path()}")
                _winfo("state-hmm-widget", f"Saved HMM deployment artifact: {artifact_path}")
                self._rebuild_intrinsic_rename_rows()
                self._rebuild_full_rename_rows()
                self._open_step(1)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_cluster, self.cluster_spinner, busy=False)
                self._refresh_enablement()

    def _on_apply_clicked(self, _):
        self._set_busy(self.btn_apply, self.apply_spinner, busy=True)
        self.out_apply.clear_output()
        with self.out_apply:
            try:
                self._persist_current_settings()
                full_pkl_path = str(self.apply_full_pkl_picker.value).strip()
                intrinsic_pkl_path = str(self.apply_intrinsic_pkl_picker.value).strip()
                if full_pkl_path == "":
                    raise ValueError("Please provide a full classification .pkl path.")

                full_artifact = self._coerce_classifier_artifact(full_pkl_path, "Full")
                intrinsic_artifact = None
                if intrinsic_pkl_path != "":
                    intrinsic_artifact = self._coerce_classifier_artifact(intrinsic_pkl_path, "Intrinsic")

                self.adata_full = self._apply_classifiers_to_full_dataset(
                    intrinsic_artifact=intrinsic_artifact,
                    full_artifact=full_artifact,
                )

                n_rows = int(self.adata_full.n_obs)
                n_intrinsic = (
                    int(self.adata_full.obs[INTRINSIC_STATE_COL].astype(str).nunique())
                    if INTRINSIC_STATE_COL in self.adata_full.obs.columns
                    else 0
                )
                n_full = (
                    int(self.adata_full.obs[FULL_STATE_COL].astype(str).nunique())
                    if FULL_STATE_COL in self.adata_full.obs.columns
                    else 0
                )
                _winfo(
                    "state-hmm-widget",
                    (
                        "Apply finished: "
                        f"rows={n_rows}, intrinsic_clusters={n_intrinsic}, full_clusters={n_full}"
                    ),
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_apply, self.apply_spinner, busy=False)
                self._refresh_enablement()

    def _on_apply_hmm_artifact_clicked(self, _):
        self._set_busy(self.btn_apply_hmm_artifact, self.apply_hmm_spinner, busy=True)
        self.out_apply_hmm.clear_output()
        with self.out_apply_hmm:
            try:
                self._persist_current_settings()
                hmm_pkl_path = str(self.apply_hmm_artifact_picker.value).strip()
                if hmm_pkl_path == "":
                    raise ValueError("Please provide an HMM deployment artifact .pkl path.")

                hmm_artifact = self._coerce_hmm_deployment_artifact(hmm_pkl_path)
                self.adata_full = self._apply_hmm_deployment_artifact_to_full_dataset(
                    hmm_artifact=hmm_artifact,
                )
                self.adata_full = self._merge_metadata_into_obs(self.adata_full)
                self._hmm_deployment_artifact = hmm_artifact

                n_rows = int(self.adata_full.n_obs)
                n_intrinsic = (
                    int(self.adata_full.obs[INTRINSIC_STATE_COL].astype(str).nunique())
                    if INTRINSIC_STATE_COL in self.adata_full.obs.columns
                    else 0
                )
                n_full = (
                    int(self.adata_full.obs[FULL_STATE_COL].astype(str).nunique())
                    if FULL_STATE_COL in self.adata_full.obs.columns
                    else 0
                )
                _winfo(
                    "state-hmm-widget",
                    (
                        "HMM deployment apply finished: "
                        f"rows={n_rows}, intrinsic_clusters={n_intrinsic}, full_clusters={n_full}"
                    ),
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_apply_hmm_artifact, self.apply_hmm_spinner, busy=False)
                self._refresh_enablement()

    def _coerce_classifier_artifact(self, path, label):
        if path == "":
            raise ValueError(f"Please provide a {label.lower()} classification .pkl path.")
        if not Path(path).exists():
            raise FileNotFoundError(f"{label} classifier file not found: {path}")
        return load_state_classifier_artifact(path)

    def _coerce_hmm_deployment_artifact(self, path):
        if path == "":
            raise ValueError("Please provide an HMM deployment artifact .pkl path.")
        if not Path(path).exists():
            raise FileNotFoundError(f"HMM deployment artifact file not found: {path}")
        return load_hmm_deployment_artifact(path)

    def _apply_classifiers_to_full_dataset(self, *, intrinsic_artifact, full_artifact):
        return apply_state_classifiers_to_full_dataset(
            output_dir=self.output_dir,
            cell_type=self._current_cell_type(),
            label_classifier_artifact=intrinsic_artifact,
            full_label_classifier_artifact=full_artifact,
            combine_binary_with_continuous=False,
            verbose=True,
        )

    def _apply_hmm_deployment_artifact_to_full_dataset(self, *, hmm_artifact):
        return apply_hmm_deployment_artifact_to_full_dataset(
            output_dir=self.output_dir,
            cell_type=self._current_cell_type(),
            hmm_deployment_artifact=hmm_artifact,
            start_offset=int(self.hmm_start_offset.value),
            start_offset_fill_mode=str(self.hmm_start_offset_fill_mode.value),
            verbose=True,
        )

    def _on_rename_intrinsic_clicked(self, _):
        self._set_busy(self.btn_rename_intrinsic, self.rename_intrinsic_spinner, busy=True)
        self.out_rename_intrinsic.clear_output()
        with self.out_rename_intrinsic:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                mapping = {}
                for old_name, txt in self._intrinsic_name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)

                result = self._apply_intrinsic_rename_mapping(mapping, save_compression="lzf")
                if not bool(result.get("changed", False)):
                    _winfo("state-hmm-widget", "No intrinsic-mapping changes to apply.")
                else:
                    _winfo(
                        "state-hmm-widget",
                        f"Renamed intrinsic states and saved: {self._model_adata_path()}",
                    )
                    if result.get("mapping_yaml_path", None):
                        _winfo(
                            "state-hmm-widget",
                            f"Saved cluster name mappings YAML: {result['mapping_yaml_path']}",
                        )
                    if result.get("deployment_warning", None):
                        _winfo(
                            "state-hmm-widget",
                            "Intrinsic rename succeeded, but the HMM deployment artifact was not updated: "
                            f"{result['deployment_warning']}",
                        )
                    if result.get("quality_control_warning", None):
                        _winfo(
                            "state-hmm-widget",
                            "Intrinsic rename succeeded, but curated quality-control plots were not recreated: "
                            f"{result['quality_control_warning']}",
                        )
                self._open_step(2)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_rename_intrinsic, self.rename_intrinsic_spinner, busy=False)
                self._refresh_enablement()

    def _on_combine_intrinsic_clicked(self, _):
        self._set_busy(
            self.btn_combine_intrinsic,
            self.combine_intrinsic_spinner,
            busy=True,
            disable_buttons=[self.btn_rename_intrinsic],
        )
        self.out_rename_intrinsic.clear_output()
        with self.out_rename_intrinsic:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                if len(self._intrinsic_name_boxes) == 0:
                    raise ValueError("No intrinsic-state rows are available to combine.")

                target = str(self.intrinsic_combine_name.value).strip()
                if target == "":
                    raise ValueError("Cannot rename to empty.")

                selected = [name for name, cb in self._intrinsic_select_boxes.items() if bool(cb.value)]
                if len(selected) == 0:
                    raise ValueError("Select at least one intrinsic state to combine.")

                for old_name in selected:
                    if old_name in self._intrinsic_name_boxes:
                        self._intrinsic_name_boxes[old_name].value = target

                mapping = {}
                for old_name, txt in self._intrinsic_name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)

                result = self._apply_intrinsic_rename_mapping(mapping, save_compression="lzf")
                if not bool(result.get("changed", False)):
                    _winfo("state-hmm-widget", "No intrinsic-state combine changes to apply.")
                else:
                    _winfo(
                        "state-hmm-widget",
                        f"Combined {len(selected)} intrinsic states into '{target}' and saved.",
                    )
                    if result.get("mapping_yaml_path", None):
                        _winfo(
                            "state-hmm-widget",
                            f"Saved cluster name mappings YAML: {result['mapping_yaml_path']}",
                        )
                    if result.get("deployment_warning", None):
                        _winfo(
                            "state-hmm-widget",
                            "Intrinsic combine succeeded, but the HMM deployment artifact was not updated: "
                            f"{result['deployment_warning']}",
                        )
                    if result.get("quality_control_warning", None):
                        _winfo(
                            "state-hmm-widget",
                            "Intrinsic combine succeeded, but curated quality-control plots were not recreated: "
                            f"{result['quality_control_warning']}",
                        )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(
                    self.btn_combine_intrinsic,
                    self.combine_intrinsic_spinner,
                    busy=False,
                    disable_buttons=[self.btn_rename_intrinsic],
                )
                self._refresh_enablement()

    def _open_backprojection_for_state_col(self, state_col):
        active_button = (
            self.btn_open_intrinsic_backprojection
            if str(state_col) == INTRINSIC_STATE_COL
            else self.btn_open_backprojection
        )
        other_button = (
            self.btn_open_backprojection
            if active_button is self.btn_open_intrinsic_backprojection
            else self.btn_open_intrinsic_backprojection
        )
        self._set_busy(
            active_button,
            self.backprojection_spinner,
            busy=True,
            disable_buttons=[other_button],
        )
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                sample_name = self.backproj_sample_dd.value
                if sample_name is None or len(str(sample_name).strip()) == 0:
                    raise ValueError("Please select a sample name.")
                _winfo(
                    "state-hmm-widget",
                    (
                        f"Opening backprojection for state_col='{state_col}', "
                        f"sample '{sample_name}', cell_type '{self._current_cell_type()}'"
                    ),
                )
                if str(state_col) == INTRINSIC_STATE_COL:
                    track_features_csv = self._resolve_track_features_csv()
                    _clear_intrinsic_hmm_backprojection_cache(
                        output_dir=self.output_dir,
                        sample_name=str(sample_name),
                        cell_type=self._current_cell_type(),
                        verbose=True,
                    )
                    _show_intrinsic_hmm_backprojection(
                        adata=self.model_adata,
                        sample_name=str(sample_name),
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        track_features_csv_path=track_features_csv,
                        metadata=getattr(self.metadata_loader, "metadata", None),
                        n_workers=int(self.hmm_backprojection_workers.value),
                        run=True,
                        verbose=True,
                    )
                else:
                    track_features_csv = self._resolve_track_features_csv()
                    _clear_behavioral_state_backprojection_cache(
                        output_dir=self.output_dir,
                        sample_name=str(sample_name),
                        cell_type=self._current_cell_type(),
                        verbose=True,
                    )
                    _show_hmm_state_backprojection(
                        adata=self.model_adata,
                        sample_name=str(sample_name),
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        state_col=str(state_col),
                        track_features_csv_path=track_features_csv,
                        metadata=getattr(self.metadata_loader, "metadata", None),
                        state_colors=self._current_full_state_color_mapping(write=True),
                        n_workers=int(self.hmm_backprojection_workers.value),
                        run=True,
                        verbose=True,
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(
                    active_button,
                    self.backprojection_spinner,
                    busy=False,
                    disable_buttons=[other_button],
                )
                self._refresh_enablement()

    def _on_open_backprojection_clicked(self, _):
        self._open_backprojection_for_state_col(FULL_STATE_COL)

    def _on_open_intrinsic_backprojection_clicked(self, _):
        self._open_backprojection_for_state_col(INTRINSIC_STATE_COL)


class StateClassificationHMMDeploymentPanel(StateClassificationHMMPanel):
    """HMM-only deployment panel that saves and applies a single deployment artifact."""

    def __init__(self, metadata_loader, cell_type=None):
        super().__init__(metadata_loader=metadata_loader, cell_type=cell_type)

    def _panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.setdefault("behavioral_state_classification", {})
        cell_type = self._current_cell_type()
        if cell_type not in section:
            section[cell_type] = {}
        return section[cell_type]

    def _effective_panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.get("behavioral_state_classification", {})
        defaults = section.get("defaults", {})
        cell_cfg = section.get(self._current_cell_type(), {})
        return {**defaults, **cell_cfg}


    def _build_steps(self):
        analysis_plots_section = getattr(self, "analysis_plots_section", widgets.VBox([]))
        self._rename_full_with_description = widgets.VBox(
            [
                widgets.HTML(
                    "<span style='color:#555;'>Rename HMM intrinsic clusters combined with binary group "
                    "values (e.g. organoid_contact)</span>"
                ),
                self.rename_full_section,
            ],
            layout=widgets.Layout(gap="6px"),
        )
        self.steps = widgets.Accordion(
            children=[
                self.clustering_section,
                self.rename_intrinsic_section,
                self._rename_full_with_description,
                analysis_plots_section,
                self.backprojection_section,
            ],
            selected_index=None,
        )
        self.steps.set_title(0, "Assign HMM intrinsic behavioral states")
        self.steps.set_title(1, "Rename intrinsic states")
        self.steps.set_title(2, "Rename full states")
        self.steps.set_title(3, "Create plots")
        self.steps.set_title(4, "Backprojection")

    def _sync_apply_existing_mode(self):
        if not hasattr(self, "steps"):
            return
        analysis_plots_section = getattr(self, "analysis_plots_section", widgets.VBox([]))
        if self.apply_existing_state_classification.value:
            self.steps.children = [
                self.apply_section,
                analysis_plots_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Apply existing classification")
            self.steps.set_title(1, "Create plots")
            self.steps.set_title(2, "Backprojection")
            self.steps.selected_index = 0
        else:
            rename_full = getattr(self, "_rename_full_with_description", self.rename_full_section)
            self.steps.children = [
                self.clustering_section,
                self.rename_intrinsic_section,
                rename_full,
                analysis_plots_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Assign HMM intrinsic behavioral states")
            self.steps.set_title(1, "Rename intrinsic states")
            self.steps.set_title(2, "Rename full states")
            self.steps.set_title(3, "Create plots")
            self.steps.set_title(4, "Backprojection")
            self.steps.selected_index = None

    def _build_apply_section(self):
        BaseStateClassificationPanel._build_apply_section(self)
        self.apply_hmm_artifact_picker = self._make_hmm_artifact_picker()
        self.apply_hmm_artifact_picker.filter_pattern = "*.pkl"
        self.apply_hmm_default_paths_html = widgets.HTML("")
        self.btn_apply_hmm_artifact = widgets.Button(
            description="Apply saved HMM artifact",
            button_style="success",
            layout=widgets.Layout(width="220px"),
        )
        self.btn_apply_hmm_artifact.on_click(self._on_apply_hmm_artifact_clicked)
        self.apply_hmm_spinner = widgets.HTML(value=spinning_loader)
        self.apply_hmm_spinner.layout.display = "none"
        self.out_apply_hmm = widgets.Output()
        self.apply_section = widgets.VBox(
            [
                widgets.HTML(
                    "<i>This deployment workflow skips downstream classifier training and uses only the saved HMM deployment artifact.</i>"
                ),
                self.apply_hmm_artifact_picker,
                self.apply_hmm_default_paths_html,
                widgets.HBox([self.btn_apply_hmm_artifact, self.apply_hmm_spinner]),
                self.out_apply_hmm,
            ]
        )

    def _refresh_enablement(self):
        super()._refresh_enablement()
        has_cell_type = self._current_cell_type() != ""
        has_hmm_artifact_input = str(self.apply_hmm_artifact_picker.value).strip() != ""
        if hasattr(self, "btn_apply"):
            self.btn_apply.disabled = True
        self.btn_apply_hmm_artifact.disabled = not (has_cell_type and has_hmm_artifact_input)


StateClassificationPanel = StateClassificationHMMDeploymentPanel
