import time
import random
import pickle
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_hex
import seaborn as sns
import umap
import scanpy as sc

from sklearn.cluster import KMeans, HDBSCAN, AgglomerativeClustering
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    confusion_matrix,
    classification_report,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold


import re
import anndata as ad
import yaml
try:
    from dtaidistance import dtw_ndim
except Exception:
    dtw_ndim = None
from sklearn.metrics import silhouette_score
from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.feature_dtw import run_tcell_analysis
from behav3d.analysis.behavior.track.visualization.plots.feature_dtw import (
    plot_cluster_percentage_bars,
    plot_clustering_feature_heatmap,
    plot_feature_umap,
)
from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata
from behav3d.analysis.behavior.general import relabel_cluster_ids
from behav3d.features.state_descriptive_features import (
    extract_descibing_track_state_features, 
    scale_feature_blocks, 
    l2_normalize_features_blocks,
    drop_highly_correlated_features,
    drop_low_variance_features
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.analysis.behavior.general.leiden import (
    run_pca, 
    run_leiden_clustering
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    export_behavioral_state_backprojection_zarrs,
    show_behavioral_state_backprojection,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    select_exemplar_tracks_by_cluster,
    _build_state_color_map,
    _compute_state_bar_segments,
    _plot_statebar_segments_on_ax,
    _prepare_exemplar_backprojection_data,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
    ensure_exemplar_coordinate_columns as _shared_ensure_exemplar_coordinate_columns,
    merge_coordinate_columns_into_obs as _shared_merge_coordinate_columns_into_obs,
    resolve_exemplar_positions_csv_path as _shared_resolve_exemplar_positions_csv_path,
)
from behav3d.analysis.behavior.general.visualization.plots import plot_top_ranking_features
from behav3d.analysis.behavior.utils import (
    _mixed_label_sort_key,
    _resolve_output_dir,
    _sanitize_filename_token,
    _to_numpy_2d,
    _vdone,
    _vinfo,
    _vsave,
    _vstart,
)


def assign_track_clusters_to_full_dataset(
    adata_full,
    adata_tracks,
    cluster_col="ClusterID",
    output_col="track_behavioral_cluster",
    id_cols=("sample_name", "TrackID"),
    unassigned_label="unassigned",
    inplace=False,
):
    """
    Broadcast track-level cluster labels to all rows in the full per-timepoint dataset.
    """
    if adata_full is None or not hasattr(adata_full, "obs"):
        raise ValueError("adata_full with .obs is required.")
    if adata_tracks is None or not hasattr(adata_tracks, "obs"):
        raise ValueError("adata_tracks with .obs is required.")

    required_full = [str(c) for c in id_cols]
    missing_full = [c for c in required_full if c not in adata_full.obs.columns]
    if len(missing_full) > 0:
        raise ValueError(f"adata_full.obs missing required id columns: {missing_full}")

    required_tracks = [str(c) for c in id_cols] + [str(cluster_col)]
    missing_tracks = [c for c in required_tracks if c not in adata_tracks.obs.columns]
    if len(missing_tracks) > 0:
        raise ValueError(
            "adata_tracks.obs missing required columns for track-label assignment: "
            f"{missing_tracks}"
        )

    adata_out = adata_full if bool(inplace) else adata_full.copy()
    id_cols = [str(c) for c in id_cols]

    track_labels = (
        adata_tracks.obs[id_cols + [cluster_col]]
        .copy()
        .dropna(subset=[cluster_col])
        .groupby(id_cols, observed=True, as_index=False)[cluster_col]
        .first()
    )

    full_obs = adata_out.obs[id_cols].copy()
    merged = full_obs.merge(track_labels, on=id_cols, how="left")
    labels = pd.Series(merged[cluster_col], index=adata_out.obs.index).astype("string")
    labels = labels.fillna(str(unassigned_label)).astype(str)
    adata_out.obs[output_col] = pd.Categorical(labels)
    return adata_out


def export_track_cluster_backprojection(
    adata_full,
    adata_tracks,
    output_dir,
    cell_type,
    cluster_col="ClusterID",
    output_col="track_behavioral_cluster",
    sample_name=None,
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    n_workers=1,
    verbose=True,
):
    """
    Export track-cluster backprojections by assigning each predicted track label to the
    selected per-track window (typically the last N timepoints used for feature extraction),
    then matching by (sample_name, TrackID, position_t) during voxel mapping.
    """
    required_full = [str(sample_col), str(track_col), str(time_col)]
    missing_full = [c for c in required_full if c not in adata_full.obs.columns]
    if len(missing_full) > 0:
        raise ValueError(
            "Cannot export track-cluster backprojection: adata_full.obs missing required columns "
            f"{missing_full}"
        )

    required_tracks = [str(sample_col), str(track_col), str(cluster_col)]
    missing_tracks = [c for c in required_tracks if c not in adata_tracks.obs.columns]
    if len(missing_tracks) > 0:
        raise ValueError(
            "Cannot export track-cluster backprojection: adata_tracks.obs missing required columns "
            f"{missing_tracks}"
        )

    tracks_obs = adata_tracks.obs.copy()
    full_obs = adata_full.obs[[str(sample_col), str(track_col), str(time_col)]].copy()

    tracks_obs["__sample"] = tracks_obs[str(sample_col)].astype("string").str.strip()
    tracks_obs["__track"] = pd.to_numeric(tracks_obs[str(track_col)], errors="coerce")
    tracks_obs["__cluster"] = tracks_obs[str(cluster_col)].astype("string").str.strip()
    tracks_obs = tracks_obs.dropna(subset=["__sample", "__track", "__cluster"]).copy()
    tracks_obs = tracks_obs[(tracks_obs["__sample"] != "") & (tracks_obs["__cluster"] != "")].copy()

    if sample_name is not None and len(str(sample_name).strip()) > 0:
        selected = str(sample_name).strip()
        tracks_obs = tracks_obs[tracks_obs["__sample"] == selected].copy()
        if len(tracks_obs) == 0:
            raise ValueError(
                f"Could not export track-cluster backprojection: sample '{selected}' not found in adata_tracks."
            )

    full_obs["__sample"] = full_obs[str(sample_col)].astype("string").str.strip()
    full_obs["__track"] = pd.to_numeric(full_obs[str(track_col)], errors="coerce")
    full_obs["__time"] = pd.to_numeric(full_obs[str(time_col)], errors="coerce")
    full_obs = full_obs.dropna(subset=["__sample", "__track", "__time"]).copy()
    full_obs = full_obs[full_obs["__sample"] != ""].copy()
    full_obs["__track"] = full_obs["__track"].astype(np.int64)
    full_obs["__time"] = full_obs["__time"].astype(np.int64)
    full_obs = full_obs.sort_values(["__sample", "__track", "__time"], kind="mergesort")

    times_by_track = {
        (str(sample), int(track)): np.asarray(g["__time"].to_numpy(dtype=np.int64), dtype=np.int64)
        for (sample, track), g in full_obs.groupby(["__sample", "__track"], observed=True, sort=False)
    }

    has_tmin = "position_t_min" in tracks_obs.columns
    has_tmax = "position_t_max" in tracks_obs.columns
    has_n_timepoints = "n_timepoints" in tracks_obs.columns
    expanded_parts = []
    missing_windows = 0

    for _, row in tracks_obs.iterrows():
        s = str(row["__sample"])
        tr = int(row["__track"])
        cluster_label = str(row["__cluster"])
        times = times_by_track.get((s, tr), None)
        if times is None or times.size == 0:
            missing_windows += 1
            continue

        chosen = times
        if has_tmin and has_tmax:
            tmin = pd.to_numeric(pd.Series([row.get("position_t_min", np.nan)]), errors="coerce").iloc[0]
            tmax = pd.to_numeric(pd.Series([row.get("position_t_max", np.nan)]), errors="coerce").iloc[0]
            if pd.notna(tmin) and pd.notna(tmax):
                lo = int(min(tmin, tmax))
                hi = int(max(tmin, tmax))
                chosen = chosen[(chosen >= lo) & (chosen <= hi)]

        if has_n_timepoints:
            n_pos = pd.to_numeric(pd.Series([row.get("n_timepoints", np.nan)]), errors="coerce").iloc[0]
            if pd.notna(n_pos):
                n_pos = max(0, int(n_pos))
                if n_pos > 0 and chosen.size > n_pos:
                    chosen = chosen[-n_pos:]

        if chosen.size == 0:
            missing_windows += 1
            continue

        expanded_parts.append(
            pd.DataFrame(
                {
                    str(sample_col): np.repeat(s, chosen.size),
                    str(track_col): np.repeat(tr, chosen.size),
                    str(time_col): chosen.astype(np.int64, copy=False),
                    str(output_col): np.repeat(cluster_label, chosen.size),
                }
            )
        )

    if len(expanded_parts) == 0:
        raise ValueError(
            "Could not export track-cluster backprojection: no matching (sample, TrackID, position_t) "
            "rows were resolved from adata_tracks windows."
        )

    backproj_obs = pd.concat(expanded_parts, ignore_index=True)
    backproj_obs = backproj_obs.sort_values(
        [str(sample_col), str(track_col), str(time_col)],
        kind="mergesort",
    ).drop_duplicates(
        subset=[str(sample_col), str(track_col), str(time_col)],
        keep="last",
    )
    _vinfo(
        verbose,
        "trajectory-apply",
        "expanded track windows to frame assignments | "
        f"tracks={int(len(tracks_obs))} | rows={int(len(backproj_obs))} | "
        f"missing_tracks={int(missing_windows)}",
    )

    adata_export = SimpleNamespace(
        obs=backproj_obs,
        uns={},
        filename=getattr(adata_full, "filename", None),
    )

    manifest = export_behavioral_state_backprojection_zarrs(
        adata=adata_export,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=output_col,
        sample_col=sample_col,
        track_col=track_col,
        time_col=time_col,
        enforce_time_coverage=False,
        n_workers=n_workers,
        verbose=verbose,
    )
    manifest["window_backprojection_rows"] = int(len(backproj_obs))
    manifest["window_backprojection_missing_tracks"] = int(missing_windows)
    return manifest


def build_track_state_sequence_lookup(
    adata_full,
    sample_name=None,
    state_col="behavioral_state",
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
):
    """Build a lookup from (sample_name, TrackID) to the full ordered state sequence."""
    required_cols = [str(sample_col), str(track_col), str(time_col), str(state_col)]
    missing_cols = [col for col in required_cols if col not in adata_full.obs.columns]
    if missing_cols:
        raise ValueError(
            "Cannot build track state sequence lookup: adata_full.obs missing required columns "
            f"{missing_cols}"
        )

    obs = adata_full.obs[required_cols].copy()
    obs["__sample"] = obs[str(sample_col)].astype("string").str.strip()
    obs["__track"] = obs[str(track_col)].map(_normalize_track_id_key)
    obs["__time"] = pd.to_numeric(obs[str(time_col)], errors="coerce")
    obs = obs.dropna(subset=["__sample", "__track", "__time"]).copy()
    obs = obs[(obs["__sample"] != "") & (obs["__track"] != "")].copy()
    if sample_name is not None and len(str(sample_name).strip()) > 0:
        obs = obs[obs["__sample"] == str(sample_name).strip()].copy()

    lookup = {}
    for (sample, track), group in obs.groupby(["__sample", "__track"], observed=True, sort=False):
        lookup[(str(sample), str(track))] = (
            group.sort_values("__time", kind="mergesort")
            [[str(sample_col), str(track_col), str(time_col), str(state_col)]]
            .copy()
        )
    return lookup


def render_track_statebar_image(
    track_df,
    state_color_map,
    state_col="behavioral_state",
    time_col="position_t",
    cursor_time=None,
    title=None,
    figsize=(6.0, 1.2),
    dpi=160,
):
    """Render a single full-track state bar as an RGB image array."""
    segments, xlim = _compute_state_bar_segments(
        track_df=track_df,
        state_key=str(state_col),
        time_key=str(time_col),
        state_color_map=state_color_map,
    )
    fig, ax = plt.subplots(figsize=figsize, dpi=int(dpi))
    _plot_statebar_segments_on_ax(
        ax=ax,
        segments=segments,
        xlim=xlim,
        title=title,
        cursor_x=cursor_time,
    )
    fig.tight_layout(pad=0.3)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    rgb = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return rgb


def _build_track_cluster_lookup(adata_tracks, sample_col="sample_name", track_col="TrackID", cluster_col="ClusterID"):
    if adata_tracks is None or not hasattr(adata_tracks, "obs"):
        return {}
    required_cols = [str(sample_col), str(track_col), str(cluster_col)]
    if any(col not in adata_tracks.obs.columns for col in required_cols):
        return {}

    obs = adata_tracks.obs[required_cols].copy()
    obs["__sample"] = obs[str(sample_col)].astype("string").str.strip()
    obs["__track"] = obs[str(track_col)].map(_normalize_track_id_key)
    obs["__cluster"] = obs[str(cluster_col)].astype("string").str.strip()
    obs = obs.dropna(subset=["__sample", "__track", "__cluster"]).copy()
    obs = obs[(obs["__sample"] != "") & (obs["__track"] != "")].copy()
    return {
        (str(row["__sample"]), str(row["__track"])): str(row["__cluster"])
        for _, row in obs.drop_duplicates(["__sample", "__track"], keep="last").iterrows()
    }


def _normalize_track_id_key(value):
    if pd.isna(value):
        return ""
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(numeric):
        numeric = float(numeric)
        if np.isfinite(numeric) and numeric.is_integer():
            return str(int(numeric))
    return str(value).strip()


def _layer_value_at_world_position(layer, world_position):
    data_position = layer.world_to_data(world_position) if hasattr(layer, "world_to_data") else world_position
    data_position = tuple(int(round(float(v))) for v in data_position)
    data = layer.data
    shape = tuple(int(v) for v in getattr(data, "shape", ()))
    if len(shape) == 0 or len(data_position) < len(shape):
        return None, data_position
    idx = tuple(max(0, min(int(pos), int(size) - 1)) for pos, size in zip(data_position[:len(shape)], shape))
    value = data[idx]
    if hasattr(value, "compute"):
        value = value.compute()
    value = np.asarray(value).squeeze()
    if value.size == 0:
        return None, idx
    try:
        value = int(value.item())
    except Exception:
        value = str(value.item())
    return value, idx


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


def _add_track_statebar_click_dock(
    viewer,
    *,
    sample_name,
    adata_full,
    adata_tracks=None,
    track_layer_name="TrackID",
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
    clickable_layer = viewer.layers[clickable_layer_name] if clickable_layer_name in layer_names else viewer.layers[track_layer_name]
    track_layer = viewer.layers[track_layer_name]

    state_values, state_color_map = _build_state_color_map(
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
    dock = viewer.window.add_dock_widget(widget, area="right", name=str(title))

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

        cluster_label = cluster_lookup.get((sample_key, track_key_value), "unknown")
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
        track_value, idx = _layer_value_at_world_position(track_layer, event.position)
        if track_value is None or str(track_value) in {"0", "0.0", ""}:
            _show_message("No tracked cell under click.")
            return
        time_idx = idx[0] if len(idx) > 0 else None
        _update_for_track(track_value, time_idx=time_idx)

    clickable_layer.mouse_drag_callbacks.append(_on_click)
    viewer._behav3d_track_statebar_dock = {
        "dock": dock,
        "widget": widget,
        "info_label": info_label,
        "image_label": image_label,
        "callback": _on_click,
        "state_values": state_values,
    }
    return widget


def show_track_cluster_backprojection(
    sample_name,
    output_dir,
    cell_type,
    *,
    adata_full,
    adata_tracks=None,
    cluster_col="ClusterID",
    state_col="behavioral_state",
    state_img_path=None,
    output_col="track_behavioral_cluster",
    run=True,
    verbose=True,
):
    """
    Open the track-cluster backprojection viewer and show a full-track state bar on click.
    """
    viewer = show_behavioral_state_backprojection(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=output_col,
        state_img_path=state_img_path,
        auto_create_if_missing=False,
        refresh_if_stale=False,
        run=False,
        verbose=verbose,
    )
    _add_track_statebar_click_dock(
        viewer=viewer,
        sample_name=sample_name,
        adata_full=adata_full,
        adata_tracks=adata_tracks,
        state_col=state_col,
        cluster_col=cluster_col,
        title="Track State Bar",
    )
    if bool(run):
        import napari

        napari.run()
    return viewer
