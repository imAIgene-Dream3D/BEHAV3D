import ipywidgets as widgets
from pathlib import Path
import fnmatch 
import os
import shutil
import yaml
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from io import BytesIO
from PIL import Image as PILImage, ImageDraw, ImageFont
from copy import deepcopy
from .utils import (
    _cfg_get,
    spinning_loader,
    detect_cell_type_category,
    _DEFAULT_CONFIG,
    behav3d_calculated_features,
    Debouncer,
)
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
    filter_multicolor_inputs,
    has_dead_channel,
    is_multicolor_celltype,
    list_grouping_candidate_columns,
)
from behav3d.core.utils import expand_column_patterns
from behav3d.io.formats.zarr import load_zarr
from behav3d.io.images import load_image
from behav3d.features.timepoint_features import run_feature_extraction
from behav3d.analysis.backprojection import backproject_columns, view_napari, backproject_feature_at_timepoint
from behav3d.analysis.behavior.track.feature_dtw import run_tcell_analysis
from behav3d.analysis.behavior.track.visualization.plots.feature_dtw import (
    plot_cluster_percentage_bars,
    plot_clustering_feature_heatmap,
    plot_feature_umap,
)
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    hash_stable_label_color_map,
)
from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
)
from behav3d.analysis.organoid_analysis import (
    run_organoid_analysis,
    plot_multi_organoid_death_dynamics,
    run_organoid_morphology_dead_analysis
)
from behav3d.analysis.filtering import filter_tracks, preview_track_length_before_filtering
from behav3d.analysis.behavior.state.visualization.backprojection import (
    _resolve_raw_image_path,
    _resolve_tracked_image_path,
)

from behav3d.io.images import load_zarr
from behav3d.analysis import summarize_track_features
from behav3d.features.advanced_timepoint_features import run_active_killing_analysis
from behav3d.analysis.interaction_analysis import (
    run_interaction_analysis,
    run_multi_organoid_interaction_comparison,
)
from behav3d.analysis.invasiveness_analysis import (
    run_invasiveness_analysis,
    detect_invasiveness_targets,
    SUMMARY_STATS,
)

# Fallback for default features: behav3d_calculated_features is now imported from .utils
'''try:
    from behav3d.defaults import behav3d_calculated_features
except ImportError:
    behav3d_calculated_features = {
        "movement": ["speed", "msd", "track_duration", "dist_to_origin"],
        "intensity": ["mean_intensity_*"],
        "morphology": ["volume", "sphericity", "surface_area"],
        "contact": ["*_contact", "*_distance"],
        "death": ["mean_dead_dye", "percentage_dead_mask"],
        "active_killing": ["is_active_killing", "killing_efficiency"]
    }'''

# Settle time for death-threshold widgets before recalculating/persisting, so
# retyping a value (e.g. 0.035 -> 0.040) or holding the spinner arrows doesn't
# trigger a recompute/YAML write on every intermediate value.
THRESHOLD_DEBOUNCE_SECONDS = 1.2

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

def _feature_output_csv_path(output_dir, cell_type):
    return Path(output_dir, "analysis", cell_type, "track_features", f"BEHAV3D_{cell_type}_combined_track_features.csv")


def _filtered_feature_output_csv_path(output_dir, cell_type):
    return Path(
        output_dir, "analysis", cell_type, "track_features",
        f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
    )


def _preferred_feature_csv_path(output_dir, cell_type):
    """Filtered combined-features CSV if it exists, else the unfiltered one.

    Preferring the filtered CSV (when present) means ``original_TrackID``
    handling in ``backproject_feature_at_timepoint`` is meaningful: that
    column is only added by filtering (to snapshot pre-split TrackIDs), so
    always reading the unfiltered CSV would make that preference dead code.
    """
    filtered_path = _filtered_feature_output_csv_path(output_dir, cell_type)
    if filtered_path.exists():
        return filtered_path
    return _feature_output_csv_path(output_dir, cell_type)


def _percentile_contrast_limits(img, mask=None):
    """1st/99th percentile contrast limits over ``img[mask]`` (or ``img > 0`` if mask is None)."""
    vals = img[mask] if mask is not None else img[img > 0]
    if vals.size == 0:
        return 0, max(1, float(img.max()) if img.size else 1)
    vmin, vmax = float(np.percentile(vals, 1)), float(np.percentile(vals, 99))
    if vmin == vmax:
        vmax = vmin + 1e-6
    return vmin, vmax


def _shared_zarr_time_max(*arrays):
    lengths = [int(a.shape[0]) for a in arrays if a is not None and hasattr(a, "shape") and len(a.shape) > 0]
    return max(0, min(lengths) - 1) if lengths else 0


def _project_label_frame(label_vol):
    """Collapse a (Z,Y,X) label volume to 2D via first-nonzero-along-Z; passthrough if already 2D."""
    if label_vol.ndim == 3:
        first_nz_idx = np.argmax(label_vol > 0, axis=0)
        return np.take_along_axis(label_vol, first_nz_idx[np.newaxis], axis=0)[0]
    return label_vol


def _csv_has_columns(csv_path, required_columns):
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return False
    try:
        columns = set(pd.read_csv(csv_path, nrows=0).columns)
    except Exception:
        return False
    return set(required_columns).issubset(columns)


def _parse_analysis_period_min(start_widget, end_widget):
    """Return ``analysis_period_min`` from optional minute fields; blank = whole movie."""

    def _val(w):
        v = w.value
        if v is None:
            return None
        if isinstance(v, str):
            s = v.strip()
            if not s:
                return None
            return float(s)
        if v == "":
            return None
        return float(v)

    start, end = _val(start_widget), _val(end_widget)
    if start is None and end is None:
        return None
    return (start, end)


def _detect_downstream_cell_types(metadata):
    """Return (organoid, immune, other) excluding per-channel multicolor and keeping merged/grouped outputs."""
    organoid_types = filter_multicolor_inputs(detect_organoid_types_from_metadata(metadata))
    immune_types = filter_multicolor_inputs(detect_immune_cell_types_from_metadata(metadata))
    other_types = filter_multicolor_inputs(detect_other_cell_types_from_metadata(metadata))

    merged_or_grouped_types = detect_merged_cell_types_from_metadata(metadata)
    for ct in merged_or_grouped_types:
        try:
            category = detect_cell_type_category(ct, metadata)
        except Exception:
            category = "immune"
        if category == "organoid":
            if ct not in organoid_types:
                organoid_types.append(ct)
        elif category == "other":
            if ct not in other_types:
                other_types.append(ct)
        else:
            if ct not in immune_types:
                immune_types.append(ct)

    return organoid_types, immune_types, other_types


def _get_channel_labels_for_sample(metadata_loader, sample_name=None):
    cellpose_cfg = _cfg_get(metadata_loader.behav3d_parameters, "cellpose", {}) or {}
    mode = str(cellpose_cfg.get("labels_mode", "same_for_all"))
    if sample_name is not None and mode == "per_sample":
        per_sample = cellpose_cfg.get("per_sample_channel_labels", {}) or {}
        labels = per_sample.get(sample_name)
        if labels:
            return labels
    return cellpose_cfg.get("channel_labels", {}) or {}


def _guess_channel_index(metadata_loader, sample_name, cell_type):
    labels = _get_channel_labels_for_sample(metadata_loader, sample_name=sample_name)
    for idx, label in labels.items():
        if str(label).strip() == str(cell_type).strip():
            try:
                return int(idx)
            except Exception:
                continue
    return None


def _persist_threshold_for_cell_types(metadata_loader, cell_types, threshold):
    params = dict(metadata_loader.behav3d_parameters or {})
    params.setdefault("features", {})
    threshold = float(threshold)

    for ct in cell_types:
        params["features"].setdefault(ct, {})["dead_mask_percentage_threshold"] = threshold

    if getattr(metadata_loader, "metadata_csv_path", None):
        params.setdefault("paths", {})["metadata_csv"] = str(Path(metadata_loader.metadata_csv_path).expanduser())
    if getattr(metadata_loader, "output_dir", None):
        params.setdefault("paths", {})["output_dir"] = str(Path(metadata_loader.output_dir).expanduser())

    metadata_loader.behav3d_parameters = params
    with metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(params, f, sort_keys=False)


def _prepare_threshold_preview_features(
    metadata_loader,
    cell_type,
    threshold,
    sample_name,
    *,
    contact_threshold=None,
    n_workers=None,
    overwrite=False,
):
    output_dir = Path(metadata_loader.output_dir).expanduser()
    csv_path = _feature_output_csv_path(output_dir, cell_type)
    required_columns = {"sample_name", "TrackID", "position_t", "percentage_dead_mask"}

    if _csv_has_columns(csv_path, required_columns):
        try:
            existing = pd.read_csv(csv_path, usecols=["sample_name"], dtype=str)
            if sample_name in existing["sample_name"].values:
                return csv_path
        except Exception:
            pass

    fcfg = _cfg_get(metadata_loader.behav3d_parameters, f"features.{cell_type}", {}) or {}
    if n_workers is None:
        n_workers = int(fcfg.get("n_workers", max(8, (os.cpu_count() or 8))))
    if contact_threshold is None:
        contact_threshold = float(fcfg.get("contact_threshold", 0.0))

    sample_metadata = metadata_loader.metadata[
        metadata_loader.metadata["sample_name"] == sample_name
    ]

    run_feature_extraction(
        metadata=sample_metadata,
        output_dir=str(output_dir),
        cell_type=cell_type,
        features_choice=["intensity", "death"],
        dead_mask_percentage_threshold=float(threshold),
        contact_threshold=float(contact_threshold),
        overwrite=bool(overwrite),
        n_workers=int(n_workers),
    )
    return csv_path


def _reset_feature_intermediate_csvs(metadata_loader, output_dir, cell_type, *, has_dead):
    sample_names = metadata_loader.metadata.get("sample_name", pd.Series(dtype=str))
    for sample_name in sample_names.dropna().astype(str).unique():
        intermediate_dir = Path(output_dir, "trackdata", sample_name, cell_type, "intermediate_results")
        stale_csvs = [
            intermediate_dir / f"{sample_name}_{cell_type}_intensity.csv",
            intermediate_dir / f"{sample_name}_{cell_type}_contact.csv",
        ]
        if has_dead:
            stale_csvs.append(intermediate_dir / f"{sample_name}_{cell_type}_dead_mask.csv")
        for csv_path in stale_csvs:
            if csv_path.exists():
                csv_path.unlink()


class _CellTypeFeatureExtractionPanel:
    """
    Single cell-type feature extraction row.
    """
    def __init__(
        self,
        metadata_loader,
        cell_type,
        *,
        shared_threshold_widget=None,
        threshold_targets=None,
        enable_preview=True,
        show_title=True,
        shared_n_workers_widget=None,
    ):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = Path(self.metadata_loader.output_dir).expanduser()
        self.category = detect_cell_type_category(self.cell_type, self.metadata_loader.metadata)
        self.has_dead = has_dead_channel(self.metadata_loader.metadata)
        self.shared_threshold_widget = shared_threshold_widget
        self.threshold_targets = list(threshold_targets or [self.cell_type])
        self._owns_threshold_widget = self.shared_threshold_widget is None
        self._show_preview = bool(enable_preview and self.has_dead and self._owns_threshold_widget)

        fcfg = _cfg_get(self.metadata_loader.behav3d_parameters, f"features.{self.cell_type}", {}) or {}

        if self._owns_threshold_widget:
            self.dead_mask_threshold = widgets.BoundedFloatText(
                description="Dead % threshold",
                value=float(fcfg.get("dead_mask_percentage_threshold", 0.05)),
                min=0.0,
                max=1.0,
                step=0.001,
                style={'description_width': '160px'},
                layout=widgets.Layout(width="220px")
            )
            if self.has_dead:
                self._threshold_debouncer = Debouncer(self._on_threshold_changed, wait=THRESHOLD_DEBOUNCE_SECONDS)
                self.dead_mask_threshold.observe(self._threshold_debouncer, names="value")
            else:
                self.dead_mask_threshold.layout.display = "none"
        else:
            self.dead_mask_threshold = self.shared_threshold_widget

        self._all_features = ["movement", "intensity", "morphology", "contact", "invasiveness", "death"]
        self._mandatory_features = {"intensity", "contact"}
        if self.category in {"immune", "other"}:
            self._mandatory_features.add("movement")
        if self.has_dead:
            self._mandatory_features.add("death")
        self._optional_features = {"movement", "morphology"} if self.category == "organoid" else {"morphology"}
        # Add invasiveness as optional feature for immune cells only
        if self.category == "immune":
            self._optional_features.add("invasiveness")

        default_feats = fcfg.get("features_choice", self._all_features)
        if not isinstance(default_feats, (list, tuple)):
            default_feats = self._all_features
        default_feats = [f for f in default_feats if f in self._all_features] or self._all_features
        for feature_name in self._mandatory_features:
            if feature_name not in default_feats:
                default_feats.append(feature_name)

        self.feature_checks = {
            f: widgets.Checkbox(description=f.capitalize(), value=(f in default_feats), indent=False)
            for f in self._all_features
        }
        for feature_name in self._mandatory_features:
            if feature_name in self.feature_checks:
                self.feature_checks[feature_name].value = True
                self.feature_checks[feature_name].disabled = True

        self.contact_threshold = widgets.FloatText(
            description="Contact Threshold (µm)",
            value=float(fcfg.get("contact_threshold", 0.0)),
            layout=widgets.Layout(width="240px"),
            style={'description_width': '180px'}
        )

        self.contact_threshold.description = "Contact Threshold (\u00B5m)"
        
        # Build contact_row with contact checkbox and threshold
        contact_row_items = [self.feature_checks["contact"], self.contact_threshold]
        
        # Add invasiveness checkbox inline for immune cells only
        if self.category == "immune":
            invasiveness_label = "Organoid Invasiveness (Advanced)"
            self.feature_checks["invasiveness"].description = invasiveness_label
            contact_row_items.append(self.feature_checks["invasiveness"])
        else:
            # Hide invasiveness checkbox for non-immune cell types
            if "invasiveness" in self.feature_checks:
                self.feature_checks["invasiveness"].layout.display = "none"
        
        contact_row = widgets.HBox(contact_row_items, layout=widgets.Layout(align_items="center", gap="12px"))

        def _toggle_contact_threshold(change=None):
            self.contact_threshold.layout.display = None if self.feature_checks["contact"].value else "none"
            # Also toggle invasiveness visibility with contact
            if "invasiveness" in self.feature_checks and self.category == "immune":
                self.feature_checks["invasiveness"].layout.display = None if self.feature_checks["contact"].value else "none"

        _toggle_contact_threshold()
        self.feature_checks["contact"].observe(_toggle_contact_threshold, names="value")

        feat_rows = []
        for f in self._all_features:
            if f == "contact":
                feat_rows.append(contact_row)
            elif f == "invasiveness":
                # Skip invasiveness here - it's already in contact_row for immune cells
                continue
            elif f == "death" and not self.has_dead:
                continue
            else:
                feat_rows.append(self.feature_checks[f])

        self.features_box = widgets.VBox([
            widgets.HTML("<b>Features</b>"),
            widgets.HTML("<div style='font-size:12px;color:#666;'>Checked disabled features are always computed. Only the category-specific optional features stay editable.</div>"),
        ] + feat_rows)

        if shared_n_workers_widget is not None:
            self.n_workers = shared_n_workers_widget
            self._owns_n_workers_widget = False
        else:
            self.n_workers = widgets.IntText(
                description="Workers",
                value=int(fcfg.get("n_workers", max(8, (os.cpu_count() or 8)))),
                max=max(8, (os.cpu_count() or 8)),
                style={'description_width': '160px'}
            )
            self._owns_n_workers_widget = True

        self.overwrite = widgets.Checkbox(description="Overwrite existing", value=bool(fcfg.get("overwrite", False)))

        self.btn_run = widgets.Button(description=f"Run {self.cell_type} feature extraction", button_style="success", layout=widgets.Layout(width="fit-content", flex="0 0 auto"))
        self.btn_run.on_click(self._on_run_clicked)

        self.btn_preview = None
        self.preview_spinner_html = widgets.HTML(value=spinning_loader)
        self.preview_spinner_html.layout.display = "none"
        self.preview_status = widgets.Output()
        self.preview_container = widgets.VBox([])
        self.sample_dd = None
        if self._show_preview:
            samples = self.metadata_loader.metadata['sample_name'].tolist()
            self.sample_dd = widgets.Dropdown(
                options=samples,
                value=samples[0] if samples else None,
                description="Sample:",
                layout=widgets.Layout(width="250px"),
            )
            self.btn_preview = widgets.ToggleButton(
                description="Preview Threshold",
                button_style="info",
                layout=widgets.Layout(width="180px")
            )
            self.btn_preview.observe(self._on_preview_toggled, names="value")

        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"

        self.run_row = widgets.HBox([self.btn_run, self.spinner_html], layout=widgets.Layout(align_items="center", gap="10px"))
        self.out = widgets.Output()

        title = widgets.HTML(f"<div style='font-size:14px;font-weight:700;'>{self.cell_type} Feature Extraction</div>")

        controls = []
        if self.has_dead:
            if self._show_preview:
                controls.append(
                    widgets.HTML(
                        "<div style='font-size:12px;color:#666;'>Select a sample, then preview. The final extraction replaces the preview CSV with the complete feature table for all samples.</div>"
                    )
                )
                controls.append(self.sample_dd)
                controls.append(
                    widgets.HBox(
                        [self.dead_mask_threshold, self.btn_preview, self.preview_spinner_html],
                        layout=widgets.Layout(align_items="center", gap="10px")
                    )
                )
                controls.append(self.preview_container)
            elif self._owns_threshold_widget:
                controls.append(self.dead_mask_threshold)
            else:
                controls.append(
                    widgets.HTML("<div style='font-size:12px;color:#666;'>Uses the shared organoid dead % threshold above.</div>")
                )

        ui_children = []
        if show_title:
            ui_children.append(title)
        ui_children.extend(controls)
        n_workers_row = []
        if self._owns_n_workers_widget:
            n_workers_row.append(widgets.HBox([self.n_workers], layout=widgets.Layout(align_items="center")))
        ui_children.extend(n_workers_row)
        ui_children.extend([
            self.features_box,
            self.overwrite,
            self.run_row,
            self.out,
            widgets.HTML("<hr>")
        ])
        self.ui = widgets.VBox(ui_children, layout=widgets.Layout(gap="6px"))

    def _selected_features(self):
        selected = []
        for feature_name in self._all_features:
            if feature_name == "death" and not self.has_dead:
                continue
            if feature_name in self._mandatory_features or self.feature_checks[feature_name].value:
                selected.append(feature_name)
        return selected

    def _current_dead_threshold(self):
        return float(self.dead_mask_threshold.value)

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        prof = params.setdefault("features", {}).setdefault(self.cell_type, {})
        prof.update({
            "dead_mask_percentage_threshold": self._current_dead_threshold(),
            "features_choice": self._selected_features(),
            "n_workers": int(self.n_workers.value),
            "overwrite": bool(self.overwrite.value),
            "contact_threshold": float(self.contact_threshold.value)
        })
        if getattr(self.metadata_loader, "metadata_csv_path", None):
            params.setdefault("paths", {})["metadata_csv"] = str(Path(self.metadata_loader.metadata_csv_path).expanduser())
        if getattr(self.metadata_loader, "output_dir", None):
            params.setdefault("paths", {})["output_dir"] = str(Path(self.metadata_loader.output_dir).expanduser())
        if self.has_dead:
            for ct in self.threshold_targets:
                params.setdefault("features", {}).setdefault(ct, {})["dead_mask_percentage_threshold"] = self._current_dead_threshold()

        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    def _lock(self, state: bool):
        widgets_to_lock = [self.n_workers, self.overwrite, self.btn_run, self.contact_threshold]
        if self._owns_threshold_widget and self.has_dead:
            widgets_to_lock.append(self.dead_mask_threshold)
        if self.btn_preview is not None:
            widgets_to_lock.append(self.btn_preview)
        for w in widgets_to_lock:
            if hasattr(w, "disabled"):
                w.disabled = state
        for feature_name, checkbox in self.feature_checks.items():
            if feature_name == "death" and not self.has_dead:
                continue
            checkbox.disabled = bool(state or feature_name in self._mandatory_features)

    def _on_threshold_changed(self, _):
        if self.has_dead:
            _persist_threshold_for_cell_types(
                self.metadata_loader,
                self.threshold_targets,
                self._current_dead_threshold(),
            )

    def _prepare_preview_data(self, sample_name):
        _persist_threshold_for_cell_types(
            self.metadata_loader,
            self.threshold_targets,
            self._current_dead_threshold(),
        )
        return _prepare_threshold_preview_features(
            self.metadata_loader,
            self.cell_type,
            self._current_dead_threshold(),
            sample_name=sample_name,
            contact_threshold=float(self.contact_threshold.value),
            n_workers=int(self.n_workers.value),
            overwrite=bool(self.overwrite.value),
        )

    def _reset_required_intermediates(self):
        _reset_feature_intermediate_csvs(
            self.metadata_loader,
            self.output_dir,
            self.cell_type,
            has_dead=self.has_dead,
        )

    def _on_preview_toggled(self, change):
        if not change["new"]:
            self.preview_container.children = []
            return

        selected_sample = self.sample_dd.value if self.sample_dd is not None else None
        if selected_sample is None:
            with self.preview_status:
                print("Please select a sample first.")
            self.btn_preview.value = False
            return

        self.preview_spinner_html.layout.display = None
        self.preview_status.clear_output()
        try:
            with self.preview_status:
                print(f"Preparing preview data for {self.cell_type} – sample {selected_sample}...")
                self._prepare_preview_data(selected_sample)
                print("Preview data ready.")

            preview = DeathThresholdPreview(
                self.metadata_loader,
                cell_type=self.cell_type,
                threshold_widget=self.dead_mask_threshold,
                sample_name=selected_sample,
            )
            self.preview_container.children = [self.preview_status, preview.ui]
        except Exception:
            with self.preview_status:
                traceback.print_exc()
            self.preview_container.children = [self.preview_status]
        finally:
            self.preview_spinner_html.layout.display = "none"

    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                self.output_dir.mkdir(parents=True, exist_ok=True)
                if bool(self.overwrite.value):
                    self._reset_required_intermediates()
                run_feature_extraction(
                    dead_mask_percentage_threshold=(self._current_dead_threshold() if self.has_dead else None),
                    contact_threshold=float(self.contact_threshold.value),
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(self.output_dir),
                    features_choice=self._selected_features(),
                    cell_type=self.cell_type,
                    n_workers=int(self.n_workers.value),
                    overwrite=bool(self.overwrite.value)
                )
                print("Feature extraction finished.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)


class _SharedOrganoidThresholdPanel:
    """
    Shared preview/threshold row for all detected organoid types.
    """
    def __init__(self, metadata_loader, cell_types, *, shared_n_workers_widget=None):
        self.metadata_loader = metadata_loader
        self.cell_types = [str(ct).strip() for ct in cell_types]
        self.shared_n_workers_widget = shared_n_workers_widget

        threshold_default = None
        for ct in self.cell_types:
            fcfg = _cfg_get(self.metadata_loader.behav3d_parameters, f"features.{ct}", {}) or {}
            if "dead_mask_percentage_threshold" in fcfg:
                threshold_default = float(fcfg["dead_mask_percentage_threshold"])
                break
        if threshold_default is None:
            threshold_default = 0.05

        samples = self.metadata_loader.metadata['sample_name'].tolist()
        self.sample_dd = widgets.Dropdown(
            options=samples,
            value=samples[0] if samples else None,
            description="Sample:",
            layout=widgets.Layout(width="250px"),
        )

        self.dead_mask_threshold = widgets.BoundedFloatText(
            description="Dead % threshold",
            value=float(threshold_default),
            min=0.0,
            max=1.0,
            step=0.001,
            style={'description_width': '160px'},
            layout=widgets.Layout(width="220px")
        )
        self._threshold_debouncer = Debouncer(self._on_threshold_changed, wait=THRESHOLD_DEBOUNCE_SECONDS)
        self.dead_mask_threshold.observe(self._threshold_debouncer, names="value")

        self.btn_preview = widgets.ToggleButton(
            description="Preview Threshold",
            button_style="info",
            layout=widgets.Layout(width="180px")
        )
        self.btn_preview.observe(self._on_preview_toggled, names="value")

        self.preview_spinner_html = widgets.HTML(value=spinning_loader)
        self.preview_spinner_html.layout.display = "none"
        self.preview_status = widgets.Output()
        self.preview_container = widgets.VBox([])

        self.ui = widgets.VBox([
            widgets.HTML("<div style='font-size:12px;color:#666;'>Shared organoid threshold preview. Select a sample, then preview. The final extraction replaces the preview CSV with the complete feature table for all samples.</div>"),
            self.sample_dd,
            widgets.HBox(
                [self.dead_mask_threshold, self.btn_preview, self.preview_spinner_html],
                layout=widgets.Layout(align_items="center", gap="10px")
            ),
            self.preview_container,
        ], layout=widgets.Layout(gap="6px"))

        _persist_threshold_for_cell_types(self.metadata_loader, self.cell_types, self.dead_mask_threshold.value)

    def _on_threshold_changed(self, _):
        _persist_threshold_for_cell_types(
            self.metadata_loader,
            self.cell_types,
            float(self.dead_mask_threshold.value),
        )

    def _prepare_preview_data(self, sample_name):
        _persist_threshold_for_cell_types(
            self.metadata_loader,
            self.cell_types,
            float(self.dead_mask_threshold.value),
        )
        n_workers_val = None
        if self.shared_n_workers_widget is not None:
            n_workers_val = int(self.shared_n_workers_widget.value)
        for ct in self.cell_types:
            fcfg = _cfg_get(self.metadata_loader.behav3d_parameters, f"features.{ct}", {}) or {}
            _prepare_threshold_preview_features(
                self.metadata_loader,
                ct,
                float(self.dead_mask_threshold.value),
                sample_name=sample_name,
                contact_threshold=float(fcfg.get("contact_threshold", 0.0)),
                n_workers=n_workers_val if n_workers_val is not None else int(fcfg.get("n_workers", max(8, (os.cpu_count() or 8)))),
                overwrite=bool(fcfg.get("overwrite", False)),
            )

    def _on_preview_toggled(self, change):
        if not change["new"]:
            self.preview_container.children = []
            return

        selected_sample = self.sample_dd.value
        if selected_sample is None:
            with self.preview_status:
                print("Please select a sample first.")
            self.btn_preview.value = False
            return

        self.preview_spinner_html.layout.display = None
        self.preview_status.clear_output()
        try:
            with self.preview_status:
                print(f"Preparing shared organoid preview data for sample {selected_sample}...")
                self._prepare_preview_data(selected_sample)
                print("Preview data ready.")

            preview = DeathThresholdPreview(
                self.metadata_loader,
                cell_types=self.cell_types,
                threshold_widget=self.dead_mask_threshold,
                sample_name=selected_sample,
            )
            self.preview_container.children = [self.preview_status, preview.ui]
        except Exception:
            with self.preview_status:
                traceback.print_exc()
            self.preview_container.children = [self.preview_status]
        finally:
            self.preview_spinner_html.layout.display = "none"


class FeatureExtractionPanel:
    """
    Feature extraction UI.

    When ``cell_type`` is supplied, renders a single per-cell-type panel.
    When ``cell_type`` is omitted, renders the grouped dynamic feature
    extraction UI used by the main notebook.
    """
    def __init__(self, metadata_loader, cell_type=None):
        self.metadata_loader = metadata_loader
        self.output_dir = Path(self.metadata_loader.output_dir).expanduser()
        self.metadata = getattr(self.metadata_loader, "metadata", None)
        self.has_dead = has_dead_channel(self.metadata)
        self.cell_type = str(cell_type).strip() if cell_type else None

        if self.cell_type:
            self.single_panel = _CellTypeFeatureExtractionPanel(
                metadata_loader=self.metadata_loader,
                cell_type=self.cell_type,
                show_title=True,
                enable_preview=self.has_dead,
            )
            self.ui = self.single_panel.ui
            self.btn_run = self.single_panel.btn_run
            self.run_buttons = [self.btn_run]
            return

        organoid_types, immune_types, other_types = _detect_downstream_cell_types(self.metadata)

        self.cell_panels = {}
        self.run_buttons = []
        blocks = [widgets.HTML("<h3>Feature Extraction</h3>")]

        if not (organoid_types or immune_types or other_types):
            blocks.append(widgets.HTML("<i>No cell types detected in metadata.</i>"))
            self.ui = widgets.VBox(blocks)
            return

        # ── Shared n_workers widget at the top of the panel ──
        # Read default from the first cell type that has a config, or fall back to CPU count.
        all_cell_types = list(organoid_types) + list(immune_types) + list(other_types)
        default_n_workers = max(8, (os.cpu_count() or 8))
        for ct in all_cell_types:
            fcfg = _cfg_get(self.metadata_loader.behav3d_parameters, f"features.{ct}", {}) or {}
            if "n_workers" in fcfg:
                default_n_workers = int(fcfg["n_workers"])
                break

        self.n_workers = widgets.IntText(
            description="Workers",
            value=default_n_workers,
            max=max(8, (os.cpu_count() or 8)),
            style={'description_width': '160px'}
        )
        blocks.append(
            widgets.HBox(
                [self.n_workers],
                layout=widgets.Layout(align_items="center", margin="0 0 6px 0")
            )
        )

        self.organoid_threshold_panel = None
        if organoid_types:
            blocks.append(widgets.HTML("<h4>Organoids</h4>"))
            if self.has_dead:
                self.organoid_threshold_panel = _SharedOrganoidThresholdPanel(
                    metadata_loader=self.metadata_loader,
                    cell_types=organoid_types,
                    shared_n_workers_widget=self.n_workers,
                )
                blocks.append(self.organoid_threshold_panel.ui)
            for ct in organoid_types:
                panel = _CellTypeFeatureExtractionPanel(
                    metadata_loader=self.metadata_loader,
                    cell_type=ct,
                    shared_threshold_widget=(
                        self.organoid_threshold_panel.dead_mask_threshold
                        if self.organoid_threshold_panel is not None else None
                    ),
                    threshold_targets=organoid_types,
                    enable_preview=False,
                    show_title=True,
                    shared_n_workers_widget=self.n_workers,
                )
                self.cell_panels[ct] = panel
                self.run_buttons.append(panel.btn_run)
                blocks.append(panel.ui)

        if immune_types:
            blocks.append(widgets.HTML("<h4>Immune Cells</h4>"))
            for ct in immune_types:
                panel = _CellTypeFeatureExtractionPanel(
                    metadata_loader=self.metadata_loader,
                    cell_type=ct,
                    show_title=True,
                    enable_preview=self.has_dead,
                    shared_n_workers_widget=self.n_workers,
                )
                self.cell_panels[ct] = panel
                self.run_buttons.append(panel.btn_run)
                blocks.append(panel.ui)

        if other_types:
            blocks.append(widgets.HTML("<h4>Other Cells</h4>"))
            for ct in other_types:
                panel = _CellTypeFeatureExtractionPanel(
                    metadata_loader=self.metadata_loader,
                    cell_type=ct,
                    show_title=True,
                    enable_preview=self.has_dead,
                    shared_n_workers_widget=self.n_workers,
                )
                self.cell_panels[ct] = panel
                self.run_buttons.append(panel.btn_run)
                blocks.append(panel.ui)

        self.ui = widgets.VBox(blocks, layout=widgets.Layout(gap="10px"))

class TrackFilterPanel:
    """
    Generic track filtering panel that works for ANY cell type.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.has_dead_channel_in_metadata = has_dead_channel(self.metadata_loader.metadata)
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        params = self.metadata_loader.behav3d_parameters
        cfg = params.setdefault("track_filtering", {}).setdefault(self.cell_type, {})
        
        self.exp_duration = widgets.IntText(description=f"Max timepoints ({cfg.get('time_type','frames')})", value=int(cfg.get("exp_duration", 350)), style={'description_width': '180px'})
        self.en_exp_duration = widgets.Checkbox(description="Trim down full time series", value=bool(cfg.get("exp_duration_enabled", True)), indent=False)
        self.row_exp = widgets.HBox([self.exp_duration], layout=widgets.Layout(display=(None if self.en_exp_duration.value else "none")))
        self.en_exp_duration.observe(lambda c: setattr(self.row_exp.layout, 'display', None if c['new'] else 'none'), names="value")

        self.min_track_length = widgets.IntText(description=f"Minimal length ({cfg.get('time_type','frames')})", value=int(cfg.get("min_track_length", 30)), style={'description_width': '180px'})
        self.en_min_length = widgets.Checkbox(description="Select only tracks with minimal length", value=bool(cfg.get("min_length_enabled", True)), indent=False)
        self.row_min = widgets.HBox([self.min_track_length], layout=widgets.Layout(display=(None if self.en_min_length.value else "none")))
        self.en_min_length.observe(lambda c: setattr(self.row_min.layout, 'display', None if c['new'] else 'none'), names="value")

        self.max_track_length = widgets.IntText(description=f"Maximal length ({cfg.get('time_type','frames')})", value=int(cfg.get("max_track_length", 30)), style={'description_width': '180px'})
        self.en_max_length = widgets.Checkbox(description="Trim down tracks to supplied length", value=bool(cfg.get("max_length_enabled", True)), indent=False)
        self.split_long_tracks = widgets.Checkbox(
            description="Split long tracks into chunks (new TrackIDs)",
            value=bool(cfg.get("split_long_tracks", True)),
            indent=False,
        )
        self.row_max = widgets.HBox([self.max_track_length, self.split_long_tracks], layout=widgets.Layout(display=(None if self.en_max_length.value else "none")))
        self.en_max_length.observe(lambda c: setattr(self.row_max.layout, 'display', None if c['new'] else 'none'), names="value")

        self.has_dead = has_dead_channel(self.metadata_loader.metadata)
        
        # --- Min size at first timepoint (available for ALL cell types) ---
        self.en_min_size = widgets.Checkbox(description="Filter by minimal size at first timepoint", value=bool(cfg.get("filter_min_size", False)), indent=False)
        self.min_size_val = widgets.IntText(description="Minimal size (px)", value=int(cfg.get("min_size", 1000)), style={'description_width': '180px'})
        self.row_min_size = widgets.HBox([self.min_size_val], layout=widgets.Layout(display=(None if self.en_min_size.value else "none")))
        self.en_min_size.observe(lambda c: setattr(self.row_min_size.layout, 'display', None if c['new'] else 'none'), names="value")

        # --- Dead at t=0 (available for ALL cell types, hidden if no dead channel) ---
        self.filter_t0_dead = widgets.Checkbox(description="Filter tracks that are dead at first timepoint", value=bool(cfg.get("filter_t0_dead", False)), indent=False)
        if not self.has_dead: self.filter_t0_dead.layout.display = "none"
        
        # --- Time unit toggle (available for ALL cell types) ---
        self.time_type = widgets.ToggleButtons(options=["frames", "hours"], value=cfg.get("time_type", "frames"))
        def _on_time_unit_change(change):
            if change['name'] == 'value':
                u = change['new']
                self.exp_duration.description = f"Max timepoints ({u})"
                self.min_track_length.description = f"Minimal length ({u})"
                self.max_track_length.description = f"Maximal length ({u})"
        self.time_type.observe(_on_time_unit_change, names='value')

        self.btn_run = widgets.Button(description=f"Filter {self.cell_type} tracks & summarize", button_style="success", layout=widgets.Layout(width="fit-content"))
        self.btn_run.on_click(self._on_run_clicked)
        self.btn_preview_lengths = widgets.Button(
            description="Show tracks length distribution (before filter)",
            button_style="info",
            layout=widgets.Layout(width="fit-content"),
            tooltip="Histogram per sample from combined_track_features CSV (no filtering)",
        )
        self.btn_preview_lengths.on_click(self._on_preview_lengths_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out_run = widgets.Output()
        self.out_preview = widgets.Output()

        _qc_rel = Path("analysis", self.cell_type, "quality_control", "BEHAV3D_filter_counts.pdf")
        self._qc_hint = widgets.HTML(
            f"<div style=\"font-size:12px;color:#444;margin:6px 0;\">"
            f"<b>QC:</b> After filtering, before/after track-length histograms and filter-step counts are saved in "
            f"<code>{_qc_rel.as_posix()}</code> (under your output folder). "
            f"Use the blue button to preview lengths from the <i>current</i> combined features file without filtering."
            f"</div>"
        )

        ui_elements = [
            widgets.HTML(f'<div style="font-size:14px;font-weight:700;">{self.cell_type.capitalize()} Track Filtering</div>'),
            self._qc_hint,
            widgets.HBox([self.btn_preview_lengths, self.spinner_html]),
            self.out_preview,
            self.en_exp_duration, self.row_exp,
            self.en_min_length, self.row_min,
            self.en_max_length, self.row_max,
            self.en_min_size, self.row_min_size,
            self.filter_t0_dead,
            widgets.HTML('<b>Unit for time-based filters:</b>'), self.time_type,
            self.btn_run,
            self.out_run
        ]
        self.ui = widgets.VBox(ui_elements)

    def _lock(self, locked):
        w_list = [self.en_exp_duration, self.exp_duration, self.en_min_length, self.min_track_length, self.en_max_length, self.max_track_length, self.en_min_size, self.min_size_val, self.filter_t0_dead, self.time_type, self.btn_run, self.btn_preview_lengths]
        for w in w_list: w.disabled = locked

    def _get_advanced_features_path(self):
        from behav3d.features.advanced_timepoint_features import find_advanced_features_csv
        return find_advanced_features_csv(self.output_dir, self.cell_type)

    def _on_preview_lengths_clicked(self, *_):
        self._lock(True)
        self.out_preview.clear_output()
        self.spinner_html.layout.display = None
        with self.out_preview:
            try:
                adv_path = self._get_advanced_features_path()
                df_input_path = str(adv_path) if adv_path is not None else None
                preview_track_length_before_filtering(
                    metadata=self.metadata_loader.metadata,
                    output_dir=self.output_dir,
                    cell_type=self.cell_type,
                    df_input_path=df_input_path,
                )
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

    def _on_run_clicked(self, *_):
        self._lock(True); self.out_run.clear_output(); self.spinner_html.layout.display = None
        with self.out_run:
            try:
                cfg = self.metadata_loader.behav3d_parameters["track_filtering"][self.cell_type]
                cfg.update({
                    "exp_duration_enabled": self.en_exp_duration.value,
                    "exp_duration": int(self.exp_duration.value),
                    "min_length_enabled": self.en_min_length.value,
                    "min_track_length": int(self.min_track_length.value),
                    "max_length_enabled": self.en_max_length.value,
                    "max_track_length": int(self.max_track_length.value),
                    "split_long_tracks": bool(self.split_long_tracks.value),
                    "filter_min_size": self.en_min_size.value,
                    "min_size": int(self.min_size_val.value),
                    "filter_t0_dead": bool(self.filter_t0_dead.value),
                    "time_type": str(self.time_type.value)
                })
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)

                adv_path = self._get_advanced_features_path()
                df_input_path = str(adv_path) if adv_path is not None else None
                filter_tracks(
                    metadata=self.metadata_loader.metadata,
                    output_dir=self.output_dir,
                    exp_duration=(int(self.exp_duration.value) if self.en_exp_duration.value else None),
                    min_track_length=(int(self.min_track_length.value) if self.en_min_length.value else None),
                    max_track_length=(int(self.max_track_length.value) if self.en_max_length.value else None),
                    split_long_tracks=bool(self.split_long_tracks.value),
                    filter_t0_dead=(bool(self.filter_t0_dead.value) if self.has_dead else False),
                    min_size=(int(self.min_size_val.value) if self.en_min_size.value else None),
                    cell_type=self.cell_type,
                    time_type=str(self.time_type.value),
                    plot_results=False,
                    df_input_path=df_input_path
                )
                print("✅ Filtering done. Summarizing tracks...")
                summarize_track_features(output_dir=self.output_dir, cell_type=self.cell_type)
                print(f"✅ {self.cell_type} filtering complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self._lock(False)

class ActiveKillingPanel:
    """
    Advanced feature extraction panel for Active Killing Analysis.
    """
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.killing_results = None
        
        md = self.metadata_loader.metadata
        if md is None: raise RuntimeError("Metadata not loaded.")

        self.organoid_types, self.immune_types, self.other_types = _detect_downstream_cell_types(md)
        self.potential_immune = self.immune_types + self.other_types
        self.target_types = self.organoid_types + self.other_types
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        self._cfg = params.setdefault("active_killing", deepcopy(_DEFAULT_CONFIG.get("active_killing", {})))
        
        # UI
        immune_options = self.potential_immune if self.potential_immune else ["(none detected)"]
        self.immune_dd = widgets.Dropdown(options=immune_options, value=immune_options[0] if immune_options else None, description="Immune cell:", style={'description_width': '120px'}, layout=widgets.Layout(width="280px"))
        
        target_options = self.target_types if self.target_types else ["(none detected)"]
        self.target_dd = widgets.SelectMultiple(
            options=target_options,
            value=tuple(target_options) if self.target_types else tuple(["(none detected)"]),
            description="Target cell(s):",
            style={'description_width': '120px'},
            layout=widgets.Layout(width="280px")
        )
        
        self.observation_window = widgets.IntText(description="Observation window:", value=int(self._cfg.get("observation_window", 5)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        self.death_signal_dd = widgets.Dropdown(options=["percentage_dead_mask", "mean_dead_dye", "nr_dead_mask_pixels"], value=self._cfg.get("death_signal_column", "percentage_dead_mask"), description="Death signal:", style={'description_width': '150px'}, layout=widgets.Layout(width="300px"))
        self.killing_threshold = widgets.FloatText(description="Killing threshold:", value=float(self._cfg.get("killing_threshold_multiplier", 1.5)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        self.min_contact_duration = widgets.IntText(description="Min contact duration:", value=int(self._cfg.get("min_contact_duration", 1)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        
        # Absolute threshold option
        self.use_absolute_threshold = widgets.Checkbox(description="Use absolute threshold", value=bool(self._cfg.get("use_absolute_threshold", False)), indent=False, layout=widgets.Layout(width="200px"))
        self.absolute_threshold = widgets.FloatText(description="Absolute threshold:", value=float(self._cfg.get("absolute_killing_threshold", 0.0) or 0.0), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"), disabled=not self._cfg.get("use_absolute_threshold", False))
        self.use_absolute_threshold.observe(self._on_absolute_threshold_toggle, names="value")
        self.absolute_hint_html = widgets.HTML("")
        self.absolute_hint_html.layout.display = "none"
        self.death_signal_dd.observe(lambda _: self._update_absolute_hint(), names="value")

        self.save_results = widgets.Checkbox(description="Save results to CSV", value=bool(self._cfg.get("save_results", True)), indent=False)
        self.gallery_item_count = widgets.IntText(description="Immune cells in gallery:", value=5, style={'description_width': '180px'}, layout=widgets.Layout(width="280px"))
        
        self.btn_run = widgets.Button(description="Run Active Killing Analysis", button_style="danger", icon="bolt", layout=widgets.Layout(width="260px"))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.btn_show_results = widgets.Button(description="Display Existing Results", button_style="success", icon="eye", layout=widgets.Layout(width="260px"))
        self.btn_show_results.on_click(lambda _: self._update_gallery_from_file(self._existing_results_path))
        self.btn_show_results.layout.display = "none"
        self._existing_results_path = None
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.out = widgets.Output()
        self.insights_out = widgets.Output()
        self.gallery_out = widgets.Output()
        
        self.insights_accordion = widgets.Accordion(children=[widgets.VBox([self.insights_out, self.gallery_out])], selected_index=None)
        self.insights_accordion.set_title(0, "Active Killing Insights & Gallery")
        self.insights_accordion.layout.display = "none"
        self.validation_html = widgets.HTML("")
        self._validate_inputs()
        self._update_absolute_hint()
        self.immune_dd.observe(lambda _: self._validate_inputs(), names="value")
        self.target_dd.observe(lambda _: self._validate_inputs(), names="value")
        self.target_dd.observe(lambda _: self._check_existing_results(), names="value")
        self.immune_dd.observe(lambda _: self._check_existing_results(), names="value")
        
        self.ui = widgets.VBox([
            widgets.HTML('<div style="font-size:22px;font-weight:700;">Active Killing Analysis</div>'),
            widgets.HTML('<div style="color:#555;font-size:13px;margin-bottom:10px;">Detects functional killing events. <b>Targets:</b> Select the target cell types to analyze (organoids and/or other cell types).</div>'),
            self.immune_dd, self.target_dd, self.validation_html, widgets.HTML("<hr>"),
            widgets.HBox([self.observation_window, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints after contact</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.killing_threshold, widgets.HTML('<span style="color:#666;font-size:12px;">× organoid\'s own signal at contact start</span>')], layout=widgets.Layout(align_items="center")),
            widgets.HBox([self.min_contact_duration, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.death_signal_dd], layout=widgets.Layout(align_items="center")),
            widgets.HBox([self.use_absolute_threshold, self.absolute_threshold, widgets.HTML('<span style="color:#666;font-size:12px;">death signal increase (bypasses multiplier)</span>')], layout=widgets.Layout(align_items="center")),
            self.absolute_hint_html,
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html, self.save_results, self.gallery_item_count, widgets.HTML('<span style="color:#666;font-size:12px;">x sample</span>')], layout=widgets.Layout(align_items="center", gap="15px")),
            self.btn_show_results,
            self.insights_accordion,
            self.out
        ])
        
        # Check if results already exist to enable visualize button
        self._check_existing_results()
    
    def _validate_inputs(self):
        immune = self.immune_dd.value
        messages = []
        valid = True
        if immune == "(none detected)": messages.append("⚠️ No immune cell types detected."); valid = False
        if not self.target_dd.value or self.target_dd.value[0] == "(none detected)": messages.append("⚠️ No target cell types selected."); valid = False
        if valid:
            p = Path(self.output_dir, "analysis", immune, "track_features", f"BEHAV3D_{immune}_combined_track_features.csv")
            if not p.with_name(f"BEHAV3D_{immune}_combined_track_features_filtered.csv").exists() and not p.exists():
                messages.append(f"⚠️ {immune} tracks not found. Run feature extraction first."); valid = False
        self.validation_html.value = '<span style="color:green;">✓ Ready</span>' if valid else '<br>'.join([f'<span style="color:#c00;">{m}</span>' for m in messages])
        self.btn_run.disabled = not valid

    def _on_absolute_threshold_toggle(self, change):
        """Enable/disable absolute threshold input based on checkbox."""
        self.absolute_threshold.disabled = not change["new"]
        # Disable multiplier display hint when using absolute
        if change["new"]:
            self.killing_threshold.disabled = True
        else:
            self.killing_threshold.disabled = False
        self._update_absolute_hint()

    def _update_absolute_hint(self):
        """Show a non-blocking recommendation to use nr_dead_mask_pixels with an absolute threshold."""
        if self.use_absolute_threshold.value and self.death_signal_dd.value != "nr_dead_mask_pixels":
            self.absolute_hint_html.value = (
                '<div style="color:#8a6d00;font-size:12px;">'
                '💡 Recommended: use <b>nr_dead_mask_pixels</b> as the death signal when using an '
                'absolute threshold — a flat pixel-count cutoff is easier to reason about than one '
                'expressed in a fraction or intensity scale.</div>'
            )
            self.absolute_hint_html.layout.display = None
        else:
            self.absolute_hint_html.layout.display = "none"

    def _check_existing_results(self):
        """Check if analysis results already exist to enable the display button."""
        immune = self.immune_dd.value
        if immune == "(none detected)": return
        
        selected_targets = list(self.target_dd.value)
        if not selected_targets or selected_targets[0] == "(none detected)": return
        
        subfolder = "combined" if len(selected_targets) > 1 else selected_targets[0]
        results_dir = Path(self.output_dir, "analysis", immune, "active_killing", subfolder)
        advanced_path = results_dir / f"BEHAV3D_{immune}_advanced_track_features.csv"
        if advanced_path.exists():
            self._existing_results_path = advanced_path
            self.btn_show_results.layout.display = None
        else:
            self.btn_show_results.layout.display = "none"

    def _update_gallery_from_file(self, advanced_path):
        """Load results from CSV and display gallery."""
        self.gallery_out.clear_output()
        with self.gallery_out:
            display(widgets.HTML('<b style="color:#007bff;">⏳ Loading gallery insights...</b>'))
            try:
                df_killing = pd.read_csv(advanced_path)
                # Normalize columns immediately
                if "TrackID" in df_killing.columns and "immune_track_id" not in df_killing.columns:
                    df_killing["immune_track_id"] = df_killing["TrackID"]
                
                self.killing_results = df_killing # Store for fallback
                
                # Coordinate normalization
                coord_map = {"centroid-0": "position_z", "centroid-1": "position_y", "centroid-2": "position_x"}
                for old_col, new_col in coord_map.items():
                    if old_col in df_killing.columns and new_col not in df_killing.columns:
                        df_killing[new_col] = df_killing[old_col]
                
                # Ensure visibility
                self.insights_accordion.layout.display = None
                self.insights_accordion.selected_index = 0
                
                self._display_insights(df_killing)
                self._display_gallery(df_killing)
            except Exception:
                traceback.print_exc()

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        self.insights_out.clear_output()
        self.gallery_out.clear_output()
        self.insights_accordion.selected_index = None # Collapse accordion
        self.insights_accordion.layout.display = "none"
        
        immune = self.immune_dd.value
        if immune == "(none detected)":
            with self.out: print("Please select an immune cell type.")
            self.btn_run.disabled = False
            self.spinner_html.layout.display = "none"
            return
            
        selected_targets = list(self.target_dd.value)
        if not selected_targets or selected_targets[0] == "(none detected)":
            with self.out: print("Please select at least one target cell type.")
            self.btn_run.disabled = False
            self.spinner_html.layout.display = "none"
            return

        with self.out:
            try:
                print(f"Starting Active Killing Analysis for {immune}...")
                
                # Save current parameters
                self._cfg["observation_window"] = self.observation_window.value
                self._cfg["death_signal_column"] = self.death_signal_dd.value
                self._cfg["killing_threshold_multiplier"] = self.killing_threshold.value
                self._cfg["min_contact_duration"] = self.min_contact_duration.value
                self._cfg["use_absolute_threshold"] = self.use_absolute_threshold.value
                self._cfg["absolute_killing_threshold"] = self.absolute_threshold.value
                self._cfg["save_results"] = self.save_results.value
                
                # Ensure the 'active_killing' key exists in behav3d_parameters
                if "active_killing" not in self.metadata_loader.behav3d_parameters:
                    self.metadata_loader.behav3d_parameters["active_killing"] = {}
                self.metadata_loader.behav3d_parameters["active_killing"].update(self._cfg)
                
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)

                abs_thresh = (
                    float(self.absolute_threshold.value)
                    if self.use_absolute_threshold.value
                    else None
                )

                # Run analysis
                def run_for_targets(targets, subfolder):
                    return run_active_killing_analysis(
                        metadata=self.metadata_loader.metadata,
                        output_dir=self.output_dir,
                        immune_cell_type=immune,
                        target_cell_types=targets,
                        observation_window=self.observation_window.value,
                        death_signal_column=self.death_signal_dd.value,
                        killing_threshold_multiplier=float(self.killing_threshold.value),
                        min_contact_duration=int(self.min_contact_duration.value),
                        absolute_killing_threshold=abs_thresh,
                        save_results=bool(self.save_results.value),
                        output_subfolder=subfolder
                    )
                
                df_killing = None
                stats = {}
                stale_filtered = set()
                for t in selected_targets:
                    print(f"--- Running independent analysis for target: {t} ---")
                    df_killing, _, stats = run_for_targets([t], t)
                    stale_filtered.update(stats.get("filtering_needs_rerun_for") or [])

                if len(selected_targets) > 1:
                    print("--- Running combined analysis for all selected targets ---")
                    combined_subfolder = "combined"
                    df_killing, _, stats = run_for_targets(selected_targets, combined_subfolder)
                    stale_filtered.update(stats.get("filtering_needs_rerun_for") or [])
                else:
                    combined_subfolder = selected_targets[0]

                print("✅ Active Killing Analysis complete! Loading gallery...")
                if stale_filtered:
                    display(widgets.HTML(
                        '<div style="color:#8a6d00;background:#fff8e1;border:1px solid #ffe082;'
                        'padding:6px 10px;border-radius:4px;margin:4px 0;">'
                        f'⚠️ Re-run <b>Filtering</b> for <b>{", ".join(sorted(stale_filtered))}</b> — '
                        'its filtered CSV was built before this run and doesn\'t include these '
                        'Active Killing results yet.</div>'
                    ))
                # Use the automated loader to ensure enriched columns (coordinates) are present
                results_dir = Path(self.output_dir, "analysis", immune, "active_killing", combined_subfolder)
                advanced_path = results_dir / f"BEHAV3D_{immune}_advanced_track_features.csv"
                
                if advanced_path.exists():
                    self._update_gallery_from_file(advanced_path)
                else:
                    # Fallback if file wasn't created for some reason
                    self.killing_results = df_killing
                    self.insights_accordion.layout.display = None
                    self.insights_accordion.selected_index = 0
                    self._display_insights(df_killing)
                    self._display_gallery(df_killing)

            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False

    def _display_insights(self, df_killing):
        """Display graphical insights in the accordion."""
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        self.insights_out.clear_output()
        if df_killing.empty: return
        
        # Normalize columns: enriched CSV uses TrackID, analysis uses immune_track_id
        if "TrackID" in df_killing.columns and "immune_track_id" not in df_killing.columns:
            df_killing["immune_track_id"] = df_killing["TrackID"]
        
        with self.insights_out:
            # 1. Per-sample kinetics (4 plots each)
            samples = df_killing["sample_name"].unique()
            for sample_name in samples:
                df_sample_all = df_killing[df_killing["sample_name"] == sample_name]
                df_sample_active = df_sample_all[df_sample_all["is_active_killing"]].copy()
                
                display(widgets.HTML(f"<h3 style='margin-bottom:5px;'>Sample: {sample_name}</h3>"))
                
                fig, axes = plt.subplots(2, 2, figsize=(16, 10))
                axes = axes.flatten()
                
                # Plot 1: Efficiency Distribution (Per Sample)
                if not df_sample_active.empty:
                    sns.histplot(df_sample_active["killing_efficiency"], kde=True, ax=axes[0], color='red')
                axes[0].set_title("1. Killing Efficiency Distribution", fontsize=14, fontweight='bold')
                axes[0].set_xlabel("Efficiency Score (signal increase / expected background)")
                axes[0].set_ylabel("Active Killing Events")
                
                # Plot 2: Smoothed Kinetics
                temp_counts = df_sample_active.groupby("position_t").size()
                if not temp_counts.empty:
                    max_t = int(df_sample_all["position_t"].max())
                    full_t = pd.Series(0, index=range(max_t + 1))
                    full_t.update(temp_counts)
                    window = max(5, max_t // 20)
                    smoothed = full_t.rolling(window=window, center=True).mean()
                    
                    axes[1].plot(full_t.index, full_t.values, color='darkred', alpha=0.2, label='Raw Counts')
                    axes[1].plot(smoothed.index, smoothed.values, color='red', linewidth=2, label='Kinetics Trend')
                    axes[1].fill_between(smoothed.index, 0, smoothed.values, color='red', alpha=0.1)
                    axes[1].set_title("2. Killing Intensity (Smoothed)", fontsize=14, fontweight='bold')
                    axes[1].set_xlabel("Timepoint")
                    axes[1].set_ylabel("Events / Timepoint")
                    axes[1].legend()

                # Plot 3: Cumulative Progress
                if not temp_counts.empty:
                    cumulative = full_t.cumsum()
                    axes[2].plot(cumulative.index, cumulative.values, color='darkblue', linewidth=3)
                    axes[2].fill_between(cumulative.index, 0, cumulative.values, color='blue', alpha=0.1)
                    axes[2].set_title("3. Cumulative Killing Progress", fontsize=14, fontweight='bold')
                    axes[2].set_xlabel("Timepoint")
                    axes[2].set_ylabel("Cumulative Sum of Active Killing Events")
                    axes[2].grid(True, linestyle='--', alpha=0.6)

                # Plot 4: Distribution of active killing events per cell
                if not df_sample_active.empty:
                    events_per_cell = df_sample_active.groupby("immune_track_id").size()
                    max_events = int(events_per_cell.max())
                    counts_per_bin = events_per_cell.value_counts().reindex(range(1, max_events + 1), fill_value=0)
                    axes[3].bar(counts_per_bin.index, counts_per_bin.values, color='red', edgecolor='black', alpha=0.8)
                    axes[3].set_xticks(range(1, max_events + 1))
                    axes[3].set_title("4. Distribution of Killing Events per Cell", fontsize=14, fontweight='bold')
                    axes[3].set_xlabel("Number of Active Killing Events")
                    axes[3].set_ylabel("Number of Cells")
                    axes[3].grid(True, axis='y', linestyle='--', alpha=0.6)
                else:
                    axes[3].set_title("4. Distribution of Killing Events per Cell", fontsize=14, fontweight='bold')
                    axes[3].set_xlabel("Number of Active Killing Events")
                    axes[3].set_ylabel("Number of Cells")

                plt.tight_layout()
                
                # Save per-sample plot
                
                selected_targets = list(self.target_dd.value)
                subfolder = "combined" if len(selected_targets) > 1 else selected_targets[0]
                sample_plot_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", subfolder, "plots", sample_name)
                sample_plot_dir.mkdir(parents=True, exist_ok=True)
                plot_path = sample_plot_dir / f"killing_kinetics_summary_{sample_name}.png"
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.show()
                display(widgets.HTML(f"<div style='color:green;font-size:11px;margin-top:-10px;margin-bottom:20px;'>📊 <b>Kinetics saved to:</b> {plot_path.relative_to(Path(self.output_dir))}</div>"))

            # 2. Combined Killing Efficiency Distribution (ONLY that plot)
            df_active = df_killing[df_killing["is_active_killing"]].copy()
            if not df_active.empty:
                display(widgets.HTML("<hr><h3 style='margin-bottom:5px;'>Combined Samples</h3>"))
                fig, ax = plt.subplots(figsize=(10, 6))
                sns.histplot(df_active["killing_efficiency"], kde=True, ax=ax, color='purple')
                ax.set_title("Combined Killing Efficiency Distribution", fontsize=16, fontweight='bold')
                ax.set_xlabel("Efficiency Score (signal increase / expected background)")
                ax.set_ylabel("Active Killing Events")
                
                # Save combined plot
                selected_targets = list(self.target_dd.value)
                subfolder = "combined" if len(selected_targets) > 1 else selected_targets[0]
                combined_plot_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", subfolder, "plots")
                combined_plot_dir.mkdir(parents=True, exist_ok=True)
                combined_path = combined_plot_dir / "combined_killing_efficiency_distribution.png"
                plt.savefig(combined_path, dpi=150, bbox_inches='tight')
                plt.show()
                display(widgets.HTML(f"<div style='color:green;font-size:11px;margin-top:-10px;'>📊 <b>Combined distribution saved to:</b> {combined_path.relative_to(Path(self.output_dir))}</div>"))

            # Summary Table
            n_items = int(self.gallery_item_count.value)
            top_hitters = df_killing[df_killing["is_active_killing"]].groupby(["sample_name", "immune_track_id"]).size().sort_values(ascending=False).head(n_items)
            if not top_hitters.empty:
                display(widgets.HTML("<br><b>Top Killing T cells across all samples (by events):</b>"))
                display(top_hitters.to_frame("Count"))

    def _display_gallery(self, df_killing):
        """Display the Active Killing Events Gallery organized by sample."""
        n_items_per_sample = int(self.gallery_item_count.value)
        
        display(widgets.HTML(f'<hr><b style="font-size:18px;">Active Killing Events Gallery (Top {n_items_per_sample} per sample)</b>'))
        display(widgets.HTML('<p style="font-size:12px;color:#666;">Showing Maximum Intensity Projections (MIP) using raw image data.<br>'
                             '<b>Legend:</b> <span style="color:#666;">Grayscale (Raw Structure)</span>, '
                             '<span style="color:red;">Translucent Red (Dead Mask)</span>, '
                             '<span style="color:purple;">Translucent Purple (ACTIVE T-cell)</span></p>'))
        
        df_active = df_killing[df_killing["is_active_killing"]].copy()
        if not df_active.empty:
            samples = df_active["sample_name"].unique()
            
            for sample_name in samples:
                display(widgets.HTML(f"<h3 style='margin-top:20px;background:#f0faf0;padding:5px;border-left:5px solid #28a745;'>Sample: {sample_name}</h3>"))
                
                df_sample = df_active[df_active["sample_name"] == sample_name].copy()
                all_killers = df_sample.groupby("immune_track_id").size().sort_values(ascending=False)
                max_available = len(all_killers)
                
                current_n = min(n_items_per_sample, max_available)
                if n_items_per_sample > max_available:
                     display(widgets.HTML(f'<small style="color:orange;">Only {max_available} active killers found in this sample.</small>'))
                
                top_hitters_idx = all_killers.head(current_n).index
                
                gallery_items = []
                for tid in top_hitters_idx:
                    df_hitter = df_sample[df_sample["immune_track_id"] == tid]
                    best_event = df_hitter.sort_values("killing_efficiency", ascending=False).iloc[0]
                    
                    display(widgets.HTML(f"<i>📸 Processing T-cell {tid}...</i>"))
                    item = self._generate_killing_gallery_item(best_event, df_killing)
                    if item: gallery_items.append(item)
                
                if gallery_items:
                    # No easy way to clear previous "Processing" messages without clearing everything, 
                    # but we can wrap them in a VBox to show them all at once at the end of sample processing.
                    # Or just display them. Let's just display them.
                    display(widgets.VBox(gallery_items))
                else:
                    display(widgets.HTML("<i>Could not generate gallery items for this sample.</i>"))
        else:
            display(widgets.HTML("<i>No active killing events found to display in gallery.</i>"))

    def _generate_killing_gallery_item(self, event_row, df_killing=None):
        """Generate a widget representing one killing event with images."""
        try:
            # Redundant imports removed as they are at the top
            sample_name = event_row["sample_name"]
            t_id = int(event_row["immune_track_id"])
            o_id = int(event_row["targeted_track_id"]) if pd.notna(event_row["targeted_track_id"]) else -1
            t_start = int(event_row["position_t"])
            observation_window = int(self.observation_window.value)
            t_end = t_start + (2 * observation_window) # Extend window for context
            
            # Metadata lookup
            md_match = self.metadata_loader.metadata[self.metadata_loader.metadata["sample_name"] == sample_name]
            if md_match.empty:
                print(f"  ❌ Sample '{sample_name}' not found in metadata.")
                return None
            md_row = md_match.iloc[0]
            
            # Paths - Robust column lookup (handles im_ and or_ prefixes)
            def get_col_path(suffix):
                cols = [c for c in self.metadata_loader.metadata.columns if c.endswith(suffix)]
                if not cols: return None
                val = md_row[cols[0]]
                if pd.isna(val) or str(val).strip() == "": return None
                return Path(val)

            immune_tracks_path = get_col_path(f"{self.immune_dd.value}_tracks_image_path")
            if immune_tracks_path is None:
                print(f"  ❌ Metadata column not found for: {self.immune_dd.value}_tracks_image_path")
                return None
            if not immune_tracks_path.exists():
                print(f"  ❌ Immune tracks file not found at: {immune_tracks_path}")
                return None
            
            # Auto-detect organoid type
            org_type = "organoid"
            for ct in self.metadata_loader.metadata.columns:
                if "_tracks_image_path" in ct and f"im_{self.immune_dd.value}" not in ct and f"_{self.immune_dd.value}" not in ct:
                    org_type = ct.replace("_tracks_image_path", "").replace("or_", "")
                    break
            
            raw_img_path = get_col_path("raw_image_path")
            if raw_img_path is None:
                 raw_img_path = get_col_path("image_path") # Typical fallback

            org_tracks_path = get_col_path(f"{org_type}_tracks_image_path")
            if org_tracks_path is None:
                print(f"  ❌ Metadata column not found for: {org_type}_tracks_image_path")
                return None
            if not org_tracks_path.exists():
                print(f"  ❌ Target tracks file not found at: {org_tracks_path}")
                return None
                
            dead_mask_path = Path(self.output_dir, "images", sample_name, f"{sample_name}_mask_dead.zarr")
            
            # Helper to get MIP crop
            def get_mip_crop(t, track_id, organoid_id, center_coords=None):
                try:
                    res_xy = float(md_row.get("pixel_distance_xy", 1.0))
                    res_z = float(md_row.get("pixel_distance_z", 1.0))
                    
                    # Determine coordinates: Use provided center or look up from data
                    if center_coords is not None:
                        cz, cy, cx = center_coords
                    else:
                        # Try to get coordinates from provided dataframe first
                        row_t = pd.DataFrame()
                        if df_killing is not None:
                            row_t = df_killing[(df_killing["sample_name"] == sample_name) & 
                                               (df_killing["immune_track_id"] == track_id) & 
                                               (df_killing["position_t"] == t)]
                        if row_t.empty:
                            # Fallback to self.killing_results
                            if hasattr(self, "killing_results") and self.killing_results is not None:
                                 row_t = self.killing_results[(self.killing_results["sample_name"] == sample_name) & 
                                                        (self.killing_results["immune_track_id"] == track_id) & 
                                                        (self.killing_results["position_t"] == t)]
                        
                        if row_t.empty: return None
                        cz = int(row_t["position_z"].iloc[0] / res_z)
                        cy = int(row_t["position_y"].iloc[0] / res_xy)
                        cx = int(row_t["position_x"].iloc[0] / res_xy)
                    
                    win = 60 # Slightly larger window for raw data
                    
                    def load_tp(path, timepoint):
                        if not path or not path.exists(): return None
                        try:
                            if str(path).endswith(".zarr") or str(path).endswith(".zarr.zip"):
                                img = load_zarr(path)
                                if timepoint >= img.shape[0]: return None
                                return img[timepoint]
                            img = load_image(path)
                            if timepoint >= img.shape[0]: return None
                            return img[timepoint]
                        except: return None

                    def get_crop(img):
                        if img is None: return None
                        try:
                            if img.ndim == 4: # C,Z,Y,X
                                c, z, y, x = img.shape
                                y1, y2 = max(0, cy-win), min(y, cy+win)
                                x1, x2 = max(0, cx-win), min(x, cx+win)
                                return img[:, :, y1:y2, x1:x2]
                            else: # Z,Y,X
                                z, y, x = img.shape
                                y1, y2 = max(0, cy-win), min(y, cy+win)
                                x1, x2 = max(0, cx-win), min(x, cx+win)
                                return img[:, y1:y2, x1:x2]
                        except: return None

                    # Load raw data and masks
                    raw_tp = load_tp(raw_img_path, t)
                    m_t_full = load_tp(immune_tracks_path, t)
                    m_d_full = load_tp(dead_mask_path, t)

                    c_raw = get_crop(raw_tp)
                    c_mt = get_crop(m_t_full)
                    c_md = get_crop(m_d_full)

                    if c_raw is None: return None
                    
                    # Ensure numpy arrays for computation
                    if c_mt is not None:
                        if hasattr(c_mt, 'compute'): c_mt = c_mt.compute()
                        c_mt = np.asarray(c_mt)
                    if c_raw is not None:
                        if hasattr(c_raw, 'compute'): c_raw = c_raw.compute()
                        c_raw = np.asarray(c_raw)
                    if c_md is not None:
                        if hasattr(c_md, 'compute'): c_md = c_md.compute()
                        c_md = np.asarray(c_md)

                    # 1. Grayscale Base (all channels)
                    if c_raw is not None and c_raw.ndim == 4:
                        raw_gray = np.max(c_raw, axis=0) # Z,Y,X
                        mip_gray = np.max(raw_gray, axis=0).astype(np.float32)
                        p1, p99 = np.percentile(mip_gray, [1, 99])
                        mip_gray = np.clip((mip_gray - p1) / (p99 - p1 + 1e-10), 0, 1)
                    else:
                        mip_gray = np.zeros(c_raw.shape[1:] if c_raw is not None else (120, 120), dtype=np.float32)

                    # 2. Tracks & Color Blending
                    mip_dead = np.max(c_md > 0, axis=0) if c_md is not None else np.zeros_like(mip_gray)
                    mip_active = np.max(c_mt == track_id, axis=0) if c_mt is not None else np.zeros_like(mip_gray)

                    rgb = np.zeros((mip_gray.shape[0], mip_gray.shape[1], 3), dtype=np.float32)
                    rgb[:, :, 0] = mip_gray; rgb[:, :, 1] = mip_gray; rgb[:, :, 2] = mip_gray
                    
                    alpha_red = 0.3
                    alpha_purple = 0.5
                    
                    rgb[mip_dead > 0, 0] = rgb[mip_dead > 0, 0] * (1-alpha_red) + 1.0 * alpha_red
                    rgb[mip_dead > 0, 1] = rgb[mip_dead > 0, 1] * (1-alpha_red)
                    rgb[mip_dead > 0, 2] = rgb[mip_dead > 0, 2] * (1-alpha_red)
                    
                    rgb[mip_active > 0, 0] = rgb[mip_active > 0, 0] * (1-alpha_purple) + 1.0 * alpha_purple
                    rgb[mip_active > 0, 1] = rgb[mip_active > 0, 1] * (1-alpha_purple)
                    rgb[mip_active > 0, 2] = rgb[mip_active > 0, 2] * (1-alpha_purple) + 1.0 * alpha_purple
                    
                    return np.clip(rgb * 255, 0, 255).astype(np.uint8)
                except: return None

            # Max time calculation
            if "timelapse_duration" in md_row:
                max_t = int(md_row["timelapse_duration"]) - 1
            elif df_killing is not None:
                max_t = int(df_killing[df_killing["sample_name"] == sample_name]["position_t"].max())
            else:
                max_t = t_end  # Fallback
            
            # Determine fixed center coordinates at event start
            res_xy = float(md_row.get("pixel_distance_xy", 1.0))
            res_z = float(md_row.get("pixel_distance_z", 1.0))
            fixed_coords = (
                int(event_row["position_z"] / res_z),
                int(event_row["position_y"] / res_xy),
                int(event_row["position_x"] / res_xy)
            )

            frames = []
            for t in range(t_start, min(t_end + 1, max_t + 1)):
                img = get_mip_crop(t, t_id, o_id, center_coords=fixed_coords)
                if img is not None:
                    pil_img = PILImage.fromarray(img)
                    draw = ImageDraw.Draw(pil_img)
                    text = f"T={t}"
                    # Simple text overlay at top-right
                    # Using default font and drawing a small shadow for visibility
                    w, h = pil_img.size
                    pos = (w - 50, 5) 
                    draw.text((pos[0]+1, pos[1]+1), text, fill="black") # shadow
                    draw.text(pos, text, fill="white")
                    frames.append(pil_img)
            
            if not frames:
                return None
            
            # Save GIF to output directory
            selected_targets = list(self.target_dd.value)
            subfolder = "combined" if len(selected_targets) > 1 else selected_targets[0]
            gallery_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", subfolder, "gallery", sample_name)
            gallery_dir.mkdir(parents=True, exist_ok=True)
            gif_name = f"killing_event_{sample_name}_T{t_id}_start{t_start}.gif"
            gif_path = gallery_dir / gif_name
            
            # Save with Pillow
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=200, # 200ms per frame
                loop=0
            )
            
            with open(gif_path, "rb") as f:
                gif_data = f.read()
                
            img_widget = widgets.Image(value=gif_data, format='gif')
            
            killing_window_end = min(t_start + observation_window, max_t)
            viz_window_end = min(t_end, max_t)
            
            info_html = widgets.HTML(f"""
                <div style='border: 1px solid #ddd; padding: 10px; margin-bottom: 5px; background: #f9f9f9; width: 600px;'>
                    <b>Event Details:</b> Sample: {sample_name} | T-cell ID: {t_id} | Target ID: {o_id if o_id != -1 else 'N/A'}<br>
                    <b>Efficiency:</b> {event_row['killing_efficiency']:.2f}<br>
                    <div style='font-size: 11px; color: #444; margin-top: 5px;'>
                        <b>🧪 Killing Window:</b> T={t_start} to T={killing_window_end} (Analysis period)<br>
                        <b>🎬 Visualization Window:</b> T={t_start} to T={viz_window_end} (Forensic context)
                    </div>
                    <small style='color: #666;'>Saved to: {gif_path.relative_to(Path(self.output_dir))}</small>
                </div>
            """)
            return widgets.VBox([info_html, img_widget])
            
        except Exception as e:
            with self.gallery_out:
                print(f"  ❌ Error generating gallery item for {sample_name}:")
                traceback.print_exc()
            return None


class DeathThresholdPreview:
    """
    Interactive threshold preview based on unfiltered track features.

    Supports either a single fixed cell type or a dropdown of selectable cell
    types. Shared organoid previews automatically expose the segment selector,
    while single-cell previews do not.
    """
    def __init__(
        self,
        metadata_loader,
        cell_type=None,
        *,
        cell_types=None,
        threshold_widget=None,
        sample_name=None,
    ):
        self.metadata_loader = metadata_loader
        self.output_dir = Path(self.metadata_loader.output_dir).expanduser()
        self._syncing_segment_options = False
        self._tracks_cache = {}

        if cell_types is None:
            if cell_type is None:
                raise ValueError("Either cell_type or cell_types must be provided.")
            self.cell_types = [str(cell_type).strip()]
        else:
            self.cell_types = [str(ct).strip() for ct in cell_types if str(ct).strip()]
            if not self.cell_types:
                raise ValueError("cell_types cannot be empty.")
        self.sample_name = str(sample_name).strip() if sample_name else None

        self._owns_threshold_widget = threshold_widget is None
        if threshold_widget is None:
            self.threshold_widget = widgets.BoundedFloatText(
                value=0.05,
                min=0,
                max=1,
                step=0.0005,
                description="Threshold:",
                style={'description_width': 'initial'},
                layout={'width': '130px'}
            )
        else:
            self.threshold_widget = threshold_widget

        if len(self.cell_types) > 1:
            self.cell_type_dd = widgets.Dropdown(
                options=self.cell_types,
                value=self.cell_types[0],
                description="Cell type:",
                layout=widgets.Layout(width="250px"),
            )
        else:
            self.cell_type_dd = None

        self.sample_label = widgets.HTML(
            f"<b>Sample:</b> {self.sample_name or '(none)'}",
        )
        self.time_slider = widgets.IntSlider(
            description="Time:",
            continuous_update=False,
            layout={'width': '400px'}
        )
        self.org_ch_input = widgets.IntText(
            value=0,
            description="Cell Channel:",
            style={'description_width': 'initial'},
            layout={'width': '190px'}
        )
        self.dead_ch_input = widgets.IntText(
            value=1,
            description="Dead Channel:",
            style={'description_width': 'initial'},
            layout={'width': '190px'}
        )
        self.out = widgets.Output(layout={'overflow': 'visible', 'height': 'auto'})

        self._threshold_debouncer = Debouncer(self._update, wait=THRESHOLD_DEBOUNCE_SECONDS)

        if self.cell_type_dd is not None:
            self.cell_type_dd.observe(self._on_context_change, names='value')
        self.time_slider.observe(self._update, names='value')
        self.threshold_widget.observe(self._threshold_debouncer, names='value')
        self.org_ch_input.observe(self._update, names='value')
        self.dead_ch_input.observe(self._update, names='value')

        header_row = [self.sample_label]
        if self.cell_type_dd is not None:
            header_row.insert(0, self.cell_type_dd)
        if self._owns_threshold_widget:
            header_row.append(self.threshold_widget)

        self.ui = widgets.VBox([
            widgets.HTML("<b>Visual Threshold Preview</b> <i>(Max Intensity Projection)</i>"),
            widgets.HBox(header_row, layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HBox([
                widgets.HTML("<b>Channels: </b>"),
                self.org_ch_input,
                self.dead_ch_input,
                widgets.HTML("<i>(Indices starting at 0)</i>")
            ]),
            self.time_slider,
            self.out,
        ])

        self.raw_path = None
        self.label_path = None
        self.dead_mask_path = None
        self.tracks_csv_path = None
        self._preview_block_reason = None

        if self.sample_name:
            self._on_context_change(None)

    def _current_cell_type(self):
        return self.cell_type_dd.value if self.cell_type_dd is not None else self.cell_types[0]

    def _tracks_dataframe(self, cell_type):
        if cell_type in self._tracks_cache:
            return self._tracks_cache[cell_type]

        csv_path = _feature_output_csv_path(self.output_dir, cell_type)
        if not csv_path.exists():
            self._tracks_cache[cell_type] = None
            return None

        try:
            df_tracks = pd.read_csv(csv_path)
            if "TrackID" in df_tracks.columns:
                df_tracks["TrackID"] = df_tracks["TrackID"].astype(str)
        except Exception:
            df_tracks = None
        self._tracks_cache[cell_type] = df_tracks
        return df_tracks

    def _resolve_label_path(self, row, sample, cell_type):
        for prefix in ['or', 'im', 'ot']:
            col = f"{prefix}_{cell_type}_tracks_image_path"
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                p = Path(row[col])
                if p.exists():
                    return p
                return None
        return None

    def _resolve_tracks_csv_path(self, row, cell_type):
        for prefix in ['or', 'im', 'ot']:
            col = f"{prefix}_{cell_type}_tracks_csv_path"
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                p = Path(row[col])
                if p.exists():
                    return p
                return None
        legacy_col = f"{cell_type}_tracks_csv_path"
        if legacy_col in row and pd.notna(row[legacy_col]) and str(row[legacy_col]).strip():
            p = Path(row[legacy_col])
            if p.exists():
                return p
        return None

    def _tracks_csv_has_track_rows(self, csv_path):
        try:
            df_tracks = pd.read_csv(csv_path, usecols=["TrackID"])
        except Exception:
            return False
        if df_tracks.empty or "TrackID" not in df_tracks.columns:
            return False
        return pd.to_numeric(df_tracks["TrackID"], errors="coerce").notna().any()

    def _on_context_change(self, _):
        sample = self.sample_name
        cell_type = self._current_cell_type()
        if sample is None or cell_type is None:
            return

        self._preview_block_reason = None

        md_match = self.metadata_loader.metadata[self.metadata_loader.metadata['sample_name'] == sample]
        if md_match.empty:
            self._preview_block_reason = f"Sample {sample} not found in metadata."
            return
        row = md_match.iloc[0]

        self.raw_path = None
        if 'raw_image_path' in row and pd.notna(row['raw_image_path']) and str(row['raw_image_path']).strip():
            raw_path = Path(row['raw_image_path'])
            if raw_path.exists():
                self.raw_path = raw_path
            else:
                self._preview_block_reason = (
                    f"{cell_type} preview for sample {sample} is blocked "
                    f"(raw_image_path file not found on disk: {raw_path})."
                )
        else:
            self._preview_block_reason = (
                f"{cell_type} preview for sample {sample} is blocked "
                "(raw_image_path is missing in metadata)."
            )

        self.label_path = self._resolve_label_path(row, sample, cell_type)
        if self.label_path is None:
            self._preview_block_reason = (
                f"{cell_type} channel is empty for sample {sample} "
                "(no valid tracks image path in metadata or file missing on disk)."
            )

        self.tracks_csv_path = self._resolve_tracks_csv_path(row, cell_type)
        if self.tracks_csv_path is None:
            self._preview_block_reason = (
                f"{cell_type} channel is empty for sample {sample} "
                "(no valid tracks CSV path in metadata or file missing on disk)."
            )
        elif not self._tracks_csv_has_track_rows(self.tracks_csv_path):
            self._preview_block_reason = (
                f"{cell_type} channel is empty for sample {sample} "
                "(tracks CSV has no valid TrackID rows)."
            )

        self.dead_mask_path = None
        if 'dead_mask_path' in row and pd.notna(row['dead_mask_path']) and str(row['dead_mask_path']).strip():
            dp = Path(row['dead_mask_path'])
            if dp.exists():
                self.dead_mask_path = dp
            else:
                self._preview_block_reason = (
                    f"{cell_type} preview for sample {sample} is blocked "
                    f"(dead_mask_path file not found on disk: {dp})."
                )
        else:
            self._preview_block_reason = (
                f"{cell_type} preview for sample {sample} is blocked "
                "(dead_mask_path is missing in metadata)."
            )

        # Resolve immune tracked image paths so we can subtract immune voxels
        # from the displayed dead mask (mirrors the feature-extraction-time
        # cleaning in `calculate_image_based_track_features`). Only relevant
        # when the previewed cell type is NOT itself an immune type.
        self.immune_label_paths = {}
        try:
            all_immune_types = detect_immune_cell_types_from_metadata(self.metadata_loader.metadata) or []
        except Exception:
            all_immune_types = []
        if cell_type not in all_immune_types:
            for im_type in all_immune_types:
                col = f"im_{im_type}_tracks_image_path"
                if col in row and pd.notna(row[col]) and str(row[col]).strip():
                    p = Path(row[col])
                    if p.exists():
                        self.immune_label_paths[im_type] = p

        guessed_channel = _guess_channel_index(self.metadata_loader, sample, cell_type)
        if guessed_channel is not None:
            self.org_ch_input.value = int(guessed_channel)
        if 'dead_channel' in row and pd.notna(row['dead_channel']):
            self.dead_ch_input.value = int(row['dead_channel'])

        try:
            time_lengths = []
            raw = load_zarr(self.raw_path)
            time_lengths.append(int(raw.shape[0]))
            if self.label_path is not None and Path(self.label_path).exists():
                label = load_zarr(self.label_path)
                time_lengths.append(int(label.shape[0]))
            if self.dead_mask_path is not None and Path(self.dead_mask_path).exists():
                dead_mask = load_zarr(self.dead_mask_path)
                time_lengths.append(int(dead_mask.shape[0]))

            shared_max = max(0, min(time_lengths) - 1)
            self.time_slider.max = shared_max
            if self.time_slider.value > shared_max:
                self.time_slider.value = shared_max
            elif self.time_slider.value == 0 and shared_max > 0:
                self.time_slider.value = shared_max // 2
        except Exception:
            pass

        self._update(None)

    @staticmethod
    def _contrast_limits(img):
        return _percentile_contrast_limits(img)

    @staticmethod
    def _shared_time_max(*arrays):
        return _shared_zarr_time_max(*arrays)

    def _update(self, _):
        sample = self.sample_name
        cell_type = self._current_cell_type()
        t = self.time_slider.value
        thr = float(self.threshold_widget.value)

        self.out.clear_output(wait=False)
        with self.out:
            try:
                if sample is None or cell_type is None:
                    print("Select a sample and cell type to preview.")
                    return
                if getattr(self, "_preview_block_reason", None):
                    print(self._preview_block_reason)
                    return
                if self.raw_path is None or not Path(self.raw_path).exists():
                    print("Preview raw data is not available yet. Use the preview button to prepare preview data first.")
                    return
                if self.label_path is None or not Path(self.label_path).exists():
                    print(f"Tracks image not found for {cell_type} in sample {sample}.")
                    return
                if self.dead_mask_path is None or not Path(self.dead_mask_path).exists():
                    print(f"Dead mask not found for sample {sample}.")
                    return

                org_ch = int(self.org_ch_input.value)
                dead_ch = int(self.dead_ch_input.value)

                raw_zarr = load_zarr(self.raw_path)
                label_zarr = load_zarr(self.label_path)
                dead_mask_zarr = load_zarr(self.dead_mask_path)
                shared_max = self._shared_time_max(raw_zarr, label_zarr, dead_mask_zarr)
                if t > shared_max:
                    print(
                        f"Preview timepoint {t} is outside the shared range for raw image, tracks, and dead mask. "
                        f"Use a timepoint between 0 and {shared_max}."
                    )
                    if self.time_slider.value != shared_max:
                        self.time_slider.value = shared_max
                    return
                has_channels = (raw_zarr.ndim == 5)
                n_channels = raw_zarr.shape[1] if has_channels else 1
                if has_channels:
                    if org_ch >= n_channels:
                        print(f"{cell_type} Channel {org_ch} out of range (0-{n_channels - 1}).")
                        return
                    if dead_ch >= n_channels:
                        print(f"Dead Channel {dead_ch} out of range (0-{n_channels - 1}).")
                        return
                    raw_img = np.asarray(raw_zarr[t, org_ch]).max(axis=0)
                    dead_raw = np.asarray(raw_zarr[t, dead_ch]).max(axis=0)
                else:
                    raw_img = np.asarray(raw_zarr[t]).max(axis=0)
                    dead_raw = raw_img

                label_vol = np.asarray(label_zarr[t])
                label_img = _project_label_frame(label_vol)

                # No label restriction: preview always uses all segments.
                label_view = label_img

                dead_vol = np.asarray(dead_mask_zarr[t])
                dead_mask = dead_vol.max(axis=0) if dead_vol.ndim == 3 else dead_vol

                # Mirror the feature-extraction-time cleaning: zero dead-mask
                # voxels that fall inside any immune segment at this timepoint
                # when the previewed cell type is non-immune.
                immune_label_paths = getattr(self, "immune_label_paths", {}) or {}
                if immune_label_paths:
                    immune_overlap = np.zeros_like(dead_mask, dtype=bool)
                    for im_type, im_path in immune_label_paths.items():
                        try:
                            im_zarr = load_zarr(im_path)
                            if t >= int(im_zarr.shape[0]):
                                continue
                            im_vol = np.asarray(im_zarr[t])
                            im_proj = im_vol.max(axis=0) if im_vol.ndim == 3 else im_vol
                            if im_proj.shape != dead_mask.shape:
                                continue
                            immune_overlap |= (im_proj > 0)
                        except Exception:
                            continue
                    if immune_overlap.any():
                        dead_mask = np.where(immune_overlap, 0, dead_mask)

                df_tracks = self._tracks_dataframe(cell_type)
                class_img = np.zeros((*label_view.shape, 3), dtype=float)
                if df_tracks is not None and "percentage_dead_mask" in df_tracks.columns:
                    time_mask = pd.to_numeric(df_tracks['position_t'], errors='coerce').fillna(-1).astype(int) == int(t)
                    df_t = df_tracks[(df_tracks['sample_name'] == sample) & time_mask]
                    for label_id in np.unique(label_view):
                        if label_id == 0:
                            continue
                        row = df_t[df_t['TrackID'] == str(label_id)]
                        if not row.empty:
                            perc = float(row['percentage_dead_mask'].iloc[0])
                            class_img[label_view == label_id] = [1, 0, 0] if perc >= thr else [0, 1, 0]
                        else:
                            class_img[label_view == label_id] = [0.5, 0.5, 0.5]

                vmin_r, vmax_r = self._contrast_limits(raw_img)
                vmin_d, vmax_d = self._contrast_limits(dead_raw)

                seg_display = np.ma.masked_where(label_view == 0, label_view)
                overlay_source = dead_mask
                dead_mask_overlay = np.ma.masked_where(overlay_source == 0, overlay_source)
                red_cmap = ListedColormap(['#FF2020'])

                fig, axes = plt.subplots(1, 5, figsize=(30, 6), dpi=150)
                axes[0].imshow(raw_img, cmap='gray', vmin=vmin_r, vmax=vmax_r)
                axes[0].set_title(f"Raw {cell_type} (Ch {org_ch})")
                axes[1].imshow(np.zeros_like(label_view), cmap='gray', vmin=0, vmax=1)
                axes[1].imshow(seg_display, cmap='nipy_spectral', interpolation='nearest')
                axes[1].set_title(f"Segments ({Path(self.label_path).stem})")
                axes[2].imshow(dead_raw, cmap='gray', vmin=vmin_d, vmax=vmax_d)
                axes[2].set_title(f"Raw Dead (Ch {dead_ch})")
                axes[3].imshow(raw_img, cmap='gray', vmin=vmin_r, vmax=vmax_r)
                axes[3].imshow(dead_mask_overlay, cmap=red_cmap, alpha=0.55)
                axes[3].set_title(f"Dead Mask on {cell_type}")
                axes[4].imshow(class_img)
                axes[4].set_title(f"Classification (Thr: {thr:.3f})\nRed=Dead, Green=Alive")
                for ax in axes:
                    ax.axis('off')
                plt.tight_layout()
                plt.show()

                legend_html = """
                <div style='padding: 5px; border: 1px solid #ccc; background-color: #f9f9f9; font-size: 12px;'>
                    <b>Legend:</b>
                    <span style='color: green;'>●</span> Alive |
                    <span style='color: red;'>●</span> Dead |
                    <span style='color: gray;'>●</span> Track missing in preview CSV
                </div>
                """
                display(widgets.HTML(legend_html))
                if df_tracks is None:
                    print("Preview CSV is missing. Use the preview button in Feature Extraction to prepare intensity/death features first.")
            except Exception as e:
                traceback.print_exc()
                print(f"Preview unavailable: {e}")


_FEATURE_BP_EXCLUDED_COLUMNS = {
    "TrackID", "original_TrackID", "sample_name", "position_t",
    "position_x", "position_y", "position_z", "ClusterID", "UMAP1", "UMAP2",
}


class FeatureBackprojectionPreviewPanel:
    """
    Quick, current-timepoint-only preview: pick a cell type + numeric
    feature, backproject it onto the segmentation labels for the
    currently viewed timepoint. Mirrors ``DeathThresholdPreview``'s lazy
    single-frame zarr reads (fast, low memory — no full-experiment zarr
    is written, no napari viewer opened).
    """

    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.output_dir = Path(self.metadata_loader.output_dir).expanduser()
        self._tracks_cache = {}
        self._preview_block_reason = None
        self.raw_path = None
        self.label_path = None

        md = self.metadata_loader.metadata
        organoid_types, immune_types, other_types = _detect_downstream_cell_types(md)
        self.cell_types = organoid_types + immune_types + other_types
        sample_list = sorted(map(str, md["sample_name"].unique().tolist())) if md is not None else []

        self.cell_type_dd = widgets.Dropdown(
            options=self.cell_types or ["(none)"],
            value=(self.cell_types[0] if self.cell_types else "(none)"),
            description="Cell type:", layout=widgets.Layout(width="220px"),
        )
        self.sample_dd = widgets.Dropdown(
            options=sample_list or ["(none)"],
            value=(sample_list[0] if sample_list else "(none)"),
            description="Sample:", layout=widgets.Layout(width="260px"),
        )
        self.feature_dd = widgets.Dropdown(
            options=["(none)"], value="(none)",
            description="Feature:", layout=widgets.Layout(width="260px"),
        )
        self.org_ch_input = widgets.IntText(
            value=0, description="Cell Channel:",
            style={'description_width': 'initial'}, layout={'width': '190px'},
        )
        self.time_slider = widgets.IntSlider(
            description="Time:", continuous_update=False, layout={'width': '400px'},
        )
        self.out = widgets.Output(layout={'overflow': 'visible', 'height': 'auto'})

        self._channel_debouncer = Debouncer(self._update, wait=THRESHOLD_DEBOUNCE_SECONDS)

        self.cell_type_dd.observe(self._on_context_change, names='value')
        self.sample_dd.observe(self._on_context_change, names='value')
        self.feature_dd.observe(self._update, names='value')
        self.time_slider.observe(self._update, names='value')
        self.org_ch_input.observe(self._channel_debouncer, names='value')

        self.ui = widgets.VBox([
            widgets.HTML("<b>Feature Backprojection Preview</b> <i>(current timepoint only)</i>"),
            widgets.HBox([self.cell_type_dd, self.sample_dd, self.feature_dd],
                         layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HBox([self.org_ch_input, widgets.HTML("<i>(Index starting at 0)</i>")]),
            self.time_slider,
            self.out,
        ])

        if self.cell_types and sample_list:
            self._on_context_change(None)

    # --- path / CSV resolution -------------------------------------------------
    def _resolve_label_path(self, row, cell_type):
        for prefix in ['or', 'im', 'ot']:
            col = f"{prefix}_{cell_type}_tracks_image_path"
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                p = Path(row[col])
                return p if p.exists() else None
        return None

    def _tracks_dataframe(self, cell_type):
        if cell_type in self._tracks_cache:
            return self._tracks_cache[cell_type]
        csv_path = _preferred_feature_csv_path(self.output_dir, cell_type)
        if not csv_path.exists():
            self._tracks_cache[cell_type] = (None, csv_path)
            return None, csv_path
        try:
            df = pd.read_csv(csv_path, low_memory=False)
        except Exception:
            df = None
        self._tracks_cache[cell_type] = (df, csv_path)
        return df, csv_path

    # --- context / options -------------------------------------------------
    def _update_feature_options(self, cell_type):
        df, _csv_path = self._tracks_dataframe(cell_type)
        if df is None or df.empty:
            self.feature_dd.options = ["(none)"]
            self.feature_dd.value = "(none)"
            return
        numeric_cols = [
            c for c in df.columns
            if c not in _FEATURE_BP_EXCLUDED_COLUMNS and pd.api.types.is_numeric_dtype(df[c])
        ]
        options = ["(none)"] + sorted(numeric_cols)
        self.feature_dd.options = options
        self.feature_dd.value = options[1] if len(options) > 1 else "(none)"

    def _on_context_change(self, _):
        cell_type = self.cell_type_dd.value
        sample = self.sample_dd.value
        self._preview_block_reason = None
        self.raw_path = None
        self.label_path = None

        if cell_type in (None, "(none)") or sample in (None, "(none)"):
            return

        self._update_feature_options(cell_type)

        md = self.metadata_loader.metadata
        md_match = md[md['sample_name'].astype(str) == str(sample)]
        if md_match.empty:
            self._preview_block_reason = f"Sample {sample} not found in metadata."
            return
        row = md_match.iloc[0]

        if 'raw_image_path' in row and pd.notna(row['raw_image_path']) and str(row['raw_image_path']).strip():
            p = Path(row['raw_image_path'])
            self.raw_path = p if p.exists() else None
            if self.raw_path is None:
                self._preview_block_reason = f"raw_image_path file not found on disk: {p}"
        else:
            self._preview_block_reason = "raw_image_path is missing in metadata."

        self.label_path = self._resolve_label_path(row, cell_type)
        if self.label_path is None and self._preview_block_reason is None:
            self._preview_block_reason = (
                f"{cell_type} channel is empty for sample {sample} "
                "(no valid tracks image path in metadata or file missing on disk)."
            )

        guessed_channel = _guess_channel_index(self.metadata_loader, sample, cell_type)
        if guessed_channel is not None:
            self.org_ch_input.value = int(guessed_channel)

        if self.raw_path is not None and self.label_path is not None:
            try:
                raw_zarr = load_zarr(self.raw_path)
                label_zarr = load_zarr(self.label_path)
                shared_max = _shared_zarr_time_max(raw_zarr, label_zarr)
                self.time_slider.max = shared_max
                if self.time_slider.value > shared_max:
                    self.time_slider.value = shared_max
                elif self.time_slider.value == 0 and shared_max > 0:
                    self.time_slider.value = shared_max // 2
            except Exception:
                pass

        self._update(None)

    # --- render -------------------------------------------------
    def _update(self, _):
        cell_type = self.cell_type_dd.value
        sample = self.sample_dd.value
        feature_col = self.feature_dd.value
        t = self.time_slider.value

        self.out.clear_output(wait=False)
        with self.out:
            try:
                if cell_type in (None, "(none)") or sample in (None, "(none)"):
                    print("Select a cell type and sample to preview.")
                    return
                if self._preview_block_reason:
                    print(self._preview_block_reason)
                    return
                if self.raw_path is None or self.label_path is None:
                    print("Preview data not available yet.")
                    return
                if feature_col in (None, "(none)"):
                    print("Select a feature to backproject.")
                    return

                df_tracks, csv_path = self._tracks_dataframe(cell_type)
                if df_tracks is None:
                    print(
                        f"No track-features CSV found for {cell_type} "
                        f"(expected at {csv_path}). Run Feature Extraction first."
                    )
                    return
                if feature_col not in df_tracks.columns:
                    print(f"Feature '{feature_col}' not found in {csv_path.name}.")
                    return

                org_ch = int(self.org_ch_input.value)
                raw_zarr = load_zarr(self.raw_path)
                label_zarr = load_zarr(self.label_path)
                shared_max = _shared_zarr_time_max(raw_zarr, label_zarr)
                if t > shared_max:
                    print(f"Timepoint {t} outside shared range (0-{shared_max}).")
                    if self.time_slider.value != shared_max:
                        self.time_slider.value = shared_max
                    return

                has_channels = (raw_zarr.ndim == 5)
                if has_channels:
                    n_channels = raw_zarr.shape[1]
                    if org_ch >= n_channels:
                        print(f"Channel {org_ch} out of range (0-{n_channels - 1}).")
                        return
                    raw_img = np.asarray(raw_zarr[t, org_ch]).max(axis=0)
                else:
                    raw_img = np.asarray(raw_zarr[t]).max(axis=0)

                label_vol = np.asarray(label_zarr[t])
                label_img = _project_label_frame(label_vol)

                pos_mask = pd.to_numeric(df_tracks['position_t'], errors='coerce').fillna(-1).astype(int) == int(t)
                sample_mask = df_tracks['sample_name'].astype(str) == str(sample)
                df_t = df_tracks[sample_mask & pos_mask]

                if df_t.empty:
                    print(f"No feature rows for sample {sample} at timepoint {t}.")
                    mapped = np.zeros_like(label_img)
                    ids_with_value = np.array([], dtype=np.int64)
                else:
                    try:
                        mapped, ids_with_value = backproject_feature_at_timepoint(
                            labels_frame=label_img,
                            df_features=df_t,
                            feature_col=feature_col,
                            track_col="TrackID",
                            background_value=0,
                        )
                    except ValueError as e:
                        print(f"Cannot backproject '{feature_col}': {e}")
                        return

                present_mask = (label_img > 0) & np.isin(label_img, ids_with_value)
                missing_mask = (label_img > 0) & ~present_mask

                vmin_r, vmax_r = _percentile_contrast_limits(raw_img)
                vmin_f, vmax_f = _percentile_contrast_limits(mapped.astype(float), mask=present_mask)

                seg_display = np.ma.masked_where(label_img == 0, label_img)
                feature_display = np.ma.masked_where(~present_mask, mapped.astype(float))

                cmap = plt.get_cmap('inferno').copy()
                cmap.set_bad(color='#2a2a2a')

                fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=150)
                axes[0].imshow(raw_img, cmap='gray', vmin=vmin_r, vmax=vmax_r)
                axes[0].set_title(f"Raw {cell_type} (Ch {org_ch})")
                axes[1].imshow(np.zeros_like(label_img), cmap='gray', vmin=0, vmax=1)
                axes[1].imshow(seg_display, cmap='nipy_spectral', interpolation='nearest')
                axes[1].set_title(f"Segments ({Path(self.label_path).stem})")
                im3 = axes[2].imshow(feature_display, cmap=cmap, vmin=vmin_f, vmax=vmax_f)
                if missing_mask.any():
                    hatch_overlay = np.ma.masked_where(~missing_mask, missing_mask)
                    axes[2].imshow(hatch_overlay, cmap=ListedColormap(['#888888']), alpha=0.6)
                axes[2].set_title(f"{feature_col} (t={t})")
                fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)
                for ax in axes:
                    ax.axis('off')
                plt.tight_layout()
                plt.show()

                if missing_mask.any():
                    display(widgets.HTML(
                        "<div style='padding:5px;border:1px solid #ccc;background:#f9f9f9;font-size:12px;'>"
                        "<span style='color:#888;'>&#9632;</span> Gray = segment present but no matching "
                        f"'{feature_col}' row in {csv_path.name} at this timepoint.</div>"
                    ))
            except Exception as e:
                traceback.print_exc()
                print(f"Preview unavailable: {e}")


class DeathDynamicsPanel:
    """
    Death dynamics analysis panel (organoid-specific).
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.has_dead_channel_in_metadata = has_dead_channel(self.metadata_loader.metadata)

        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        # Check the *filtered* CSV first — this is what run_organoid_analysis
        # actually reads (df_tracks_path=None resolves to it), and it's the
        # only file that exists for a cell-type *group* (see
        # behav3d.analysis.grouping.create_cell_type_group), which merges
        # already-filtered populations and never produces an unfiltered CSV.
        filtered_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        unfiltered_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        p = filtered_p if filtered_p.exists() else unfiltered_p
        self.has_death_features = False
        if p.exists():
            try:
                cols = set(pd.read_csv(p, nrows=0).columns)
                self.has_death_features = bool({"dead", "percentage_dead_mask", "nr_dead_mask_pixels"}.issubset(cols))
            except Exception:
                pass

        self.btn_run = widgets.Button(
            description=f"Run {cell_type} death dynamics",
            button_style="warning",
            layout=widgets.Layout(width="300px")
        )
        self.btn_run.on_click(self._on_run_clicked)

        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()

        # Condition columns to group/color the mean +/- SEM plots by (e.g.
        # organoid_line + macrophage_line), mirroring the state-composition
        # "group by" multi-select pattern.
        md = self.metadata_loader.metadata
        group_col_options = []
        if md is not None and hasattr(md, "columns"):
            group_col_options = [
                c for c in md.columns
                if c in ("exp_nr", "well") or c.endswith("_line_condition")
            ]
        self.select_group_cols = widgets.SelectMultiple(
            options=group_col_options,
            value=[],
            layout=widgets.Layout(width="360px", height="80px"),
            disabled=len(group_col_options) == 0,
        )

        if not self.has_dead_channel_in_metadata:
            self.ui = widgets.VBox([
                widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'),
                widgets.HTML('<div style="color:#b00;">No dead channel is configured in metadata for this dataset.</div>')
            ])
        elif not self.has_death_features:
            self.ui = widgets.VBox([
                widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'),
                widgets.HTML('<div style="color:#b00;">No final death classification found. Run feature extraction and track filtering first.</div>')
            ])
        else:
            self.ui = widgets.VBox([
                widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'),
                widgets.HTML('Group plots by condition column(s) (hold Ctrl/Cmd to select multiple):'),
                self.select_group_cols,
                widgets.HBox([self.btn_run, self.spinner_html]),
                self.out
            ])

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        with self.out:
            try:
                run_organoid_analysis(
                    output_dir=self.output_dir,
                    df_tracks_path=None,
                    org_type=self.cell_type,
                    metadata=self.metadata_loader.metadata,
                    group_cols=list(self.select_group_cols.value) or None,
                )
                print(f"{self.cell_type} death dynamics complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False


class MultiOrganoidDeathDynamicsPanel:
    """
    Combined death dynamics comparison panel for one or more organoid types.
    Shows all organoid types together on the same plot for easy comparison.
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types  # List of organoid type names
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Check which organoid types have death dynamics data available
        self.available_data = {}
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "results", f"combined_general_{org_type}_dynamics_analysis.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
        
        self.btn_refresh = widgets.Button(description="🔄 Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Death Comparison Plot", button_style="warning", layout=widgets.Layout(width="250px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.out = widgets.Output()
        
        self._refresh_status()
        
        self.ui = widgets.VBox([
            widgets.HTML('<b>Death Dynamics Comparison</b>'),
            widgets.HTML(
                '<div style="font-size:12px;color:#666;">'
                'Compare death curves across organoid types on the same plot. '
                'Works with one or more organoid types.'
                '</div>'
            ),
            widgets.HBox(
                [self.status_html, self.btn_refresh],
                layout=widgets.Layout(align_items="center", gap="10px"),
            ),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out,
        ])
    
    def _refresh_status(self):
        """Check which organoid types have death dynamics data available."""
        self.available_data = {}
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "results", f"combined_general_{org_type}_dynamics_analysis.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
        
        if len(self.available_data) < 1:
            missing = [ot for ot in self.organoid_types if ot not in self.available_data]
            self.status_html.value = f'<div style="color:#b00;">Waiting for death dynamics data from: {", ".join(missing)}</div>'
            self.btn_run.disabled = True
        else:
            self.status_html.value = f'<div style="color:#080;">Ready: {", ".join(self.available_data.keys())}</div>'
            self.btn_run.disabled = False
    
    def _on_refresh_clicked(self, *_):
        self._refresh_status()
    
    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        
        with self.out:
            try:
                params = getattr(self.metadata_loader, "behav3d_parameters", {}) or {}
                feat_cfg = params.get("features", {}) or {}
                default_thr = float(
                    _DEFAULT_CONFIG.get("death_dynamics", {})
                    .get("organoid", {})
                    .get("dead_perc_threshold", 0.02)
                )
                dead_perc_threshold_map = {}
                for org_type in self.organoid_types:
                    raw_thr = feat_cfg.get(org_type, {}).get("dead_mask_percentage_threshold", default_thr)
                    try:
                        dead_perc_threshold_map[str(org_type)] = float(raw_thr)
                    except Exception:
                        dead_perc_threshold_map[str(org_type)] = default_thr

                # Call the analysis function from organoid_analysis module
                plot_multi_organoid_death_dynamics(
                    output_dir=self.output_dir,
                    organoid_types=self.organoid_types,
                    dead_perc_threshold_map=dead_perc_threshold_map,
                )
                print("✅ Multi-organoid comparison complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False


class MultiOrganoidInteractionComparisonPanel:
    """
    Interaction Overview: violin, before-death curves, and active-killing
    dashboard across one or more organoid types.
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        md = self.metadata_loader.metadata
        _, immune_types, other_types = _detect_downstream_cell_types(md)
        self.immune_types = immune_types + other_types

        # --- widgets ---
        group_options = [
            ("By organoid type", "organoid_type"),
            ("By treatment (immune cell)", "treatment"),
        ]
        self.group_by = widgets.Dropdown(
            options=group_options,
            value="organoid_type",
            description="Group by:",
            style={"description_width": "80px"},
            layout=widgets.Layout(width="300px"),
        )

        self.annotate_lc = widgets.Checkbox(
            value=False,
            description="Annotate by immune line condition",
            indent=False,
        )

        self.time_window = widgets.IntText(
            value=60,
            description="Before-death window (min):",
            style={"description_width": "160px"},
            layout=widgets.Layout(width="260px"),
        )

        # Temporal Range: radio + Start T / End T
        self.temporal_radio = widgets.RadioButtons(
            options=["All timepoints", "Custom time range"],
            value="All timepoints",
            layout=widgets.Layout(width="200px"),
        )
        self.temporal_radio.observe(self._on_temporal_radio_changed, names="value")

        self.start_t = widgets.IntText(
            value=0,
            description="Start T:",
            style={"description_width": "60px"},
            layout=widgets.Layout(width="160px"),
        )
        self.end_t = widgets.IntText(
            value=100,
            description="End T:",
            style={"description_width": "60px"},
            layout=widgets.Layout(width="160px"),
        )
        self.custom_t_box = widgets.HBox(
            [self.start_t, self.end_t],
            layout=widgets.Layout(display="none"),
        )

        self.immune_checkboxes = {}
        for ct in self.immune_types:
            self.immune_checkboxes[ct] = widgets.Checkbox(
                value=True, description=ct, indent=False,
            )
        self.immune_box = widgets.VBox(list(self.immune_checkboxes.values()))

        self.btn_refresh = widgets.Button(
            description="🔄 Refresh", button_style="info",
            layout=widgets.Layout(width="100px"),
        )
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(
            description="Run Interaction Overview", button_style="warning",
            layout=widgets.Layout(width="280px"),
        )
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.out = widgets.Output()

        self._refresh_status()

        self.ui = widgets.VBox([
            widgets.HTML('<b>Interaction Overview</b>'),
            widgets.HTML(
                '<div style="font-size:12px;color:#666;">'
                'Produces: violin (cumulative contacts per organoid), '
                'cumulative-to-death curves, and active-killing dashboard. '
                'Works with one or more organoid types.'
                '</div>'
            ),
            widgets.HBox(
                [self.status_html, self.btn_refresh],
                layout=widgets.Layout(align_items="center", gap="10px"),
            ),
            widgets.HBox([self.group_by, self.annotate_lc]),
            widgets.HTML(
                '<div style="font-size:11px;color:#888;margin:4px 0;">'
                '<b>Temporal Range</b> — applies to: '
                '<i>cumulative interactions per organoid</i> and '
                '<i>active killing plots</i>'
                '</div>'
            ),
            self.temporal_radio,
            self.custom_t_box,
            widgets.HTML(
                '<div style="font-size:11px;color:#888;margin:4px 0;">'
                '<b>Before-death window</b> — applies to: '
                '<i>cumulative interactions before death plot</i> only'
                '</div>'
            ),
            self.time_window,
            widgets.HTML("<b>Immune / interacting cell types:</b>"),
            self.immune_box,
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out,
        ])

    def _on_temporal_radio_changed(self, change):
        if change["new"] == "Custom time range":
            self.custom_t_box.layout.display = None
        else:
            self.custom_t_box.layout.display = "none"

    # ---- status ----
    def _refresh_status(self):
        ready_types = []
        for org_type in self.organoid_types:
            csv = Path(
                self.output_dir, "analysis", org_type, "track_features",
                f"BEHAV3D_{org_type}_combined_track_features_filtered.csv",
            )
            if csv.exists():
                cols = pd.read_csv(csv, nrows=0).columns
                if any(f"{ct}_contact" in cols for ct in self.immune_types):
                    ready_types.append(org_type)

        if not ready_types:
            self.status_html.value = (
                '<div style="color:#b00;">Waiting for filtered track data '
                'with contact columns.</div>'
            )
            self.btn_run.disabled = True
        else:
            self.status_html.value = (
                f'<div style="color:#080;">Ready: {", ".join(ready_types)}</div>'
            )
            self.btn_run.disabled = False

    def _on_refresh_clicked(self, *_):
        self._refresh_status()

    def _parse_period(self):
        """Return ``analysis_period_t`` as ``(start_t, end_t)`` or ``None``."""
        if self.temporal_radio.value == "All timepoints":
            return None
        return (int(self.start_t.value), int(self.end_t.value))

    # ---- run ----
    def _on_run_clicked(self, *_):
        selected_im = [
            ct for ct, cb in self.immune_checkboxes.items() if cb.value
        ]
        if not selected_im:
            return
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()

        with self.out:
            try:
                run_multi_organoid_interaction_comparison(
                    output_dir=self.output_dir,
                    organoid_types=self.organoid_types,
                    immune_types=selected_im,
                    metadata=self.metadata_loader.metadata,
                    group_by=self.group_by.value,
                    time_window_min=float(self.time_window.value),
                    analysis_period_t=self._parse_period(),
                    annotate_line_condition=bool(self.annotate_lc.value),
                    show_plots=True,
                )
                print("✅ Interaction Overview complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False


class MultiOrganoidAnalysisPanel:
    """
    Wrapper that groups cross-target comparison panels:
      1. Death Dynamics Comparison  (works with ≥1 organoid type)
      2. Interaction Overview       (works with ≥1 organoid type)
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types

        # --- Sub-panel 1: Death Dynamics ---
        self.death_panel = MultiOrganoidDeathDynamicsPanel(
            metadata_loader=metadata_loader,
            organoid_types=organoid_types,
        )

        # --- Sub-panel 2: Interaction Comparison ---
        self.interaction_panel = MultiOrganoidInteractionComparisonPanel(
            metadata_loader=metadata_loader,
            organoid_types=organoid_types,
        )

        self.ui = widgets.VBox([
            widgets.HTML('<h3>Cross-Target Analysis</h3>'),
            widgets.HTML('<h4>1. Death Dynamics Comparison</h4>'),
            self.death_panel.ui,
            widgets.HTML('<h4>2. Interaction Overview</h4>'),
            self.interaction_panel.ui,
        ])

    @property
    def btn_run(self):
        """Expose first run button for register_buttons_in compatibility."""
        return self.death_panel.btn_run


class MultiOrganoidMorphologyDeathPanel:
    """
    Morphology vs Death analysis for all organoid types combined.
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        params = dict(self.metadata_loader.behav3d_parameters or {})
        cfg = params.setdefault("morpho_dead", {})

        self.window_start = widgets.IntText(
            description="Window Start (t >= N)",
            value=int(cfg.get("window_start", 0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.window_end = widgets.IntText(
            description="Window End (t <= N)",
            value=int(cfg.get("window_end", 5)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_n_neighbors = widgets.IntText(
            description="UMAP neighbors",
            value=int(cfg.get("umap_n_neighbors", 15)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_n_neighbors.tooltip = "Number of local neighbors used for UMAP manifold construction. Controls local vs global structure."
        
        self.umap_min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=float(cfg.get("umap_min_dist", 0.1)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_min_dist.tooltip = "Minimum distance between points in UMAP embedding. Controls how tightly points are packed."
        
        self.umap_spread = widgets.FloatText(
            description="UMAP spread",
            value=float(cfg.get("umap_spread", 1.0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_spread.tooltip = "The effective scale of embedded points."
        
        self.distance_metric = widgets.Text(
            value=cfg.get("distance_metric", "euclidean"),
            description="UMAP metric",
            placeholder="euclidean",
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.distance_metric.tooltip = "Distance metric for UMAP. See umap-learn docs for full options."
                
        self.leiden_resolution = widgets.FloatText(
            description="Leiden resolution",
            value=float(cfg.get("leiden_resolution", 1.0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.leiden_resolution.tooltip = "Controls cluster granularity. Higher values lead to more clusters."

        self.status_html = widgets.HTML("")
        self.btn_refresh = widgets.Button(description="Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Initial Morphology vs Death", button_style="warning", layout=widgets.Layout(width="280px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()

        self._refresh_status()

        self.ui = widgets.VBox([
            widgets.HTML('<b>Initial Morphology Related to Death</b>'),
            widgets.HTML('<div style="font-size:12px;color:#666;">Uses baseline morphology averaged over the selected window, excludes organoids dead in that window, and labels dead if any death occurs after.</div>'),
            widgets.HTML('<div style="font-size:12px;color:#666;">Features: 4 key morphology descriptors (volume, sphericity, solidity, surface-to-volume ratio).</div>'),
            widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.window_start, self.window_end], layout=widgets.Layout(gap="12px")),
            widgets.HBox([self.umap_n_neighbors, self.umap_min_dist], layout=widgets.Layout(gap="12px")),
            widgets.HBox([self.distance_metric, self.leiden_resolution], layout=widgets.Layout(gap="12px")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out
        ])

        # Check for organoid presence to show/hide the entire panel
        md = self.metadata_loader.metadata
        has_organoids = any(col.startswith('or_') for col in md.columns)
        if not has_organoids:
            self.ui.layout.display = "none"

    def _refresh_status(self):
        self.available_data = {}
        missing = []
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "track_features", f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
            else:
                missing.append(org_type)

        if not self.available_data:
            self.status_html.value = '<div style="color:#b00;">No organoid track features found. Run track filtering first.</div>'
            self.btn_run.disabled = True
        elif missing:
            self.status_html.value = (
                f'<div style="color:#b80;">Missing data for: {", ".join(missing)}. '
                f'Will run with available: {", ".join(self.available_data.keys())}</div>'
            )
            self.btn_run.disabled = False
        else:
            self.status_html.value = f'<div style="color:#080;">Ready: {", ".join(self.available_data.keys())}</div>'
            self.btn_run.disabled = False

    def _on_refresh_clicked(self, *_):
        self._refresh_status()

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        cfg = params.setdefault("morpho_dead", {})
        cfg.update({
            "window_start": int(self.window_start.value),
            "window_end": int(self.window_end.value),
            "umap_n_neighbors": int(self.umap_n_neighbors.value),
            "umap_min_dist": float(self.umap_min_dist.value),
            #"umap_spread": float(self.umap_spread.value),
            "distance_metric": str(self.distance_metric.value).strip() or "euclidean",
            "leiden_resolution": float(self.leiden_resolution.value),
        })
        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                params,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        with self.out:
            try:
                self._persist_params()
                run_organoid_morphology_dead_analysis(
                    output_dir=self.output_dir,
                    organoid_types=self.organoid_types,
                    initial_window=(int(self.window_start.value), int(self.window_end.value)),
                    umap_n_neighbors=int(self.umap_n_neighbors.value),
                    umap_min_dist=float(self.umap_min_dist.value),
                    #umap_spread=float(self.umap_spread.value),
                    distance_metric=str(self.distance_metric.value).strip() or "euclidean",
                    leiden_resolution=float(self.leiden_resolution.value),
                    metadata=self.metadata_loader.metadata,
                )
                print("Initial morphology vs death analysis complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False



class InteractionAnalysisPanel:
    """
    Interaction analysis panel for organoids.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        md = self.metadata_loader.metadata
        _, immune_types, other_types = _detect_downstream_cell_types(md)
        self.available_types = immune_types + other_types
        self.df_tracks_path = Path(self.output_dir, "analysis", self.cell_type, "track_features", f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        
        self.status_html = widgets.HTML("")
        self.cell_type_box = widgets.VBox([])
        self.cell_type_checkboxes = {}
        self.group_by_line_condition = widgets.Checkbox(
            value=False,
            description="Group overall plots by immune line condition",
            indent=False,
        )
        self.btn_refresh = widgets.Button(description="🔄 Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Interaction Analysis", button_style="warning", layout=widgets.Layout(width="250px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        self.ui = widgets.VBox([
            widgets.HTML(f'<b>{self.cell_type} Interaction Analysis</b>'),
            widgets.HTML(
                '<div style="font-size:12px;color:#666;">'
                'Produces: cumulative contact over time and alive vs dead bar '
                'plots (overall + per-sample). One PDF per immune type.'
                '</div>'
            ),
            widgets.HBox(
                [self.status_html, self.btn_refresh],
                layout=widgets.Layout(align_items="center", gap="10px"),
            ),
            widgets.HTML('<b>Select cell types:</b>'),
            self.cell_type_box,
            self.group_by_line_condition,
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out,
        ])
        self._refresh_data_status()

    def _refresh_data_status(self):
        self.has_data = self.df_tracks_path.exists()
        contact_types = []
        if self.has_data:
            cols = pd.read_csv(self.df_tracks_path, nrows=0).columns
            contact_types = [ct for ct in self.available_types if f"{ct}_contact" in cols]
        
        if not self.has_data or not contact_types:
            self.status_html.value = '<div style="color:#b00;">⚠️ Waiting for data/contacts.</div>'
            self.btn_run.disabled = True
        else:
            self.status_html.value = f'<div style="color:#080;">✅ Ready: {", ".join(contact_types)}</div>'
            self.btn_run.disabled = False
            self.cell_type_checkboxes = {ct: widgets.Checkbox(value=True, description=ct, indent=False) for ct in contact_types}
            self.cell_type_box.children = list(self.cell_type_checkboxes.values())

    def _on_refresh_clicked(self, *_): self._refresh_data_status()

    def _on_run_clicked(self, *_):
        sel = [ct for ct, cb in self.cell_type_checkboxes.items() if cb.value]
        if not sel: return
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                run_interaction_analysis(
                    output_dir=self.output_dir,
                    cell_type=self.cell_type,
                    interacting_cell_types=sel,
                    df_tracks_path=str(self.df_tracks_path),
                    show_plots=True,
                    group_by_line_condition=bool(self.group_by_line_condition.value),
                )
                print("✅ Interaction Analysis complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False

class InvasivenessAnalysisPanel:
    """
    Invasiveness analysis panel (immune-cell perspective).

    Measures how invasive a selected immune cell type is against one or more
    organoid targets (e.g. "27T" vs "MDO"), over time and summarized per movie.
    Reads the immune cell type's filtered track features CSV, which must contain
    ``{target}_invasiveness`` / ``{target}_invasiveness_perc`` columns produced
    during feature extraction (invasiveness feature enabled for immune cells).
    """
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        md = self.metadata_loader.metadata
        if md is None:
            raise RuntimeError("Metadata not loaded.")

        organoid_types, immune_types, other_types = _detect_downstream_cell_types(md)
        self.potential_immune = immune_types + other_types
        self.organoid_types = organoid_types  # Needed for violin-by-fate plots

        # Multi-select immune cell types (one checkbox each). Targets shown are
        # the union across all selected immune types.
        self.immune_checkboxes = {
            im: widgets.Checkbox(
                value=(i == 0), description=im, indent=False,
                layout=widgets.Layout(width="150px"),
            )
            for i, im in enumerate(self.potential_immune)
        }
        for cb in self.immune_checkboxes.values():
            cb.observe(lambda _c: self._refresh_data_status(), names="value")
        self.immune_box = widgets.HBox(
            list(self.immune_checkboxes.values()) or [widgets.HTML("(none detected)")]
        )
        self.target_box = widgets.VBox([])
        self.target_checkboxes = {}
        self.summary_stat_dd = widgets.Dropdown(
            options=list(SUMMARY_STATS),
            value="mean",
            description="Per-movie stat:",
            style={"description_width": "120px"},
            layout=widgets.Layout(width="260px"),
        )
        self.group_by_line_condition = widgets.Checkbox(
            value=False,
            description="Separate by immune line condition",
            indent=False,
        )

        # Temporal Range: radio + Start T / End T
        self.temporal_radio = widgets.RadioButtons(
            options=["All timepoints", "Custom time range"],
            value="All timepoints",
            layout=widgets.Layout(width="200px"),
        )
        self.temporal_radio.observe(self._on_temporal_radio_changed, names="value")

        self.start_t = widgets.IntText(
            value=0,
            description="Start T:",
            style={"description_width": "60px"},
            layout=widgets.Layout(width="160px"),
        )
        self.end_t = widgets.IntText(
            value=100,
            description="End T:",
            style={"description_width": "60px"},
            layout=widgets.Layout(width="160px"),
        )
        self.custom_t_box = widgets.HBox(
            [self.start_t, self.end_t],
            layout=widgets.Layout(display="none"),
        )

        self.status_html = widgets.HTML("")
        self.btn_refresh = widgets.Button(description="🔄 Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Invasiveness Analysis", button_style="warning", layout=widgets.Layout(width="260px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()

        self.ui = widgets.VBox([
            widgets.HTML("<b>Invasiveness Analysis</b>"),
            widgets.HTML(
                '<div style="color:#555;font-size:12px;">'
                'Immune-cell perspective: compare invasion against one or more '
                'organoid targets over time and per movie.<br>'
                '<b>Over-time plots:</b> fraction (% cells invasive), mean % '
                '(average contact across all cells), median % (typical cell).<br>'
                '<b>Per-movie stat</b> collapses the temporal window into one dot '
                'per movie (does not change the median-over-time curve).<br>'
                'Fate violins: one point per organoid. Enable '
                '<i>invasiveness</i> during feature extraction.'
                '</div>'
            ),
            widgets.HTML("<b>Select immune cell type(s):</b>"),
            self.immune_box,
            widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HTML("<b>Select target(s) to compare:</b>"),
            self.target_box,
            self.summary_stat_dd,
            self.group_by_line_condition,
            widgets.HTML(
                '<div style="font-size:11px;color:#888;margin:4px 0;">'
                '<b>Temporal Range</b> (in timepoints)'
                '</div>'
            ),
            self.temporal_radio,
            self.custom_t_box,
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out,
        ])

        self._refresh_data_status()

    def _on_temporal_radio_changed(self, change):
        if change["new"] == "Custom time range":
            self.custom_t_box.layout.display = None
        else:
            self.custom_t_box.layout.display = "none"

    def _filtered_csv(self, immune):
        return Path(
            self.output_dir, "analysis", immune, "track_features",
            f"BEHAV3D_{immune}_combined_track_features_filtered.csv",
        )

    def _selected_immune(self):
        return [im for im, cb in self.immune_checkboxes.items() if cb.value]

    def _refresh_data_status(self):
        selected = self._selected_immune()
        # Preserve which targets were checked across refreshes.
        prev_checked = {t for t, cb in getattr(self, "target_checkboxes", {}).items() if cb.value}
        self.target_checkboxes = {}
        self.target_box.children = []
        if not selected:
            self.status_html.value = '<div style="color:#b00;">⚠️ Select at least one immune cell type.</div>'
            self.btn_run.disabled = True
            return

        # Union of targets across all selected immune types; note any immune
        # types whose filtered CSV is missing.
        targets = set()
        missing = []
        for im in selected:
            csv_path = self._filtered_csv(im)
            if not csv_path.exists():
                missing.append(im)
                continue
            try:
                cols_df = pd.read_csv(csv_path, nrows=0)
            except Exception:
                missing.append(im)
                continue
            targets.update(detect_invasiveness_targets(cols_df))

        if missing:
            self.status_html.value = (
                f'<div style="color:#b00;">⚠️ Filtered tracks not found for: '
                f'{", ".join(missing)}. Run feature extraction (with invasiveness) '
                "+ filtering first.</div>"
            )
            self.btn_run.disabled = True
            return
        if not targets:
            self.status_html.value = (
                '<div style="color:#b00;">⚠️ No invasiveness columns found. '
                "Enable 'invasiveness' during feature extraction for these immune cells.</div>"
            )
            self.btn_run.disabled = True
            return

        targets = sorted(targets)
        self.status_html.value = f'<div style="color:#080;">✅ Ready: {", ".join(targets)}</div>'
        self.btn_run.disabled = False
        # No any_org — only real organoid targets. Default all checked (or keep
        # the user's prior selection when refreshing).
        default_all = not prev_checked
        self.target_checkboxes = {
            t: widgets.Checkbox(
                value=(default_all or t in prev_checked),
                description=t, indent=False,
            )
            for t in targets
        }
        self.target_box.children = list(self.target_checkboxes.values())

    def _on_refresh_clicked(self, *_):
        self._refresh_data_status()

    def _parse_period(self):
        """Return ``analysis_period_t`` as ``(start_t, end_t)`` or ``None``."""
        if self.temporal_radio.value == "All timepoints":
            return None
        return (int(self.start_t.value), int(self.end_t.value))

    def _on_run_clicked(self, *_):
        immune_list = self._selected_immune()
        selected = [t for t, cb in self.target_checkboxes.items() if cb.value]
        if not immune_list or not selected:
            with self.out:
                self.out.clear_output()
                print("Please select at least one immune cell type and one target.")
            return
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        with self.out:
            try:
                run_invasiveness_analysis(
                    output_dir=self.output_dir,
                    immune_cell_types=immune_list,
                    targets=selected,
                    summary_stat=self.summary_stat_dd.value,
                    show_plots=True,
                    group_by_line_condition=bool(self.group_by_line_condition.value),
                    analysis_period_t=self._parse_period(),
                    organoid_types=self.organoid_types,
                )
                print("✅ Invasiveness Analysis complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False


class MotileCellAnalysisPanel:
    """
    Generic behavioral analysis panel (DTW/UMAP/clustering) for ANY cell type.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        self._panel_cfg = params.setdefault("analysis", {}).setdefault(self.cell_type, deepcopy(_DEFAULT_CONFIG["analysis"].get(self.category, {})))
        self._params = params
        
        groups = deepcopy(behav3d_calculated_features)
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        df_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        if not df_p.exists(): df_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        self.df_tracks_path = df_p
        
        self._no_csv = False
        valid_cols = []
        if df_p.exists():
            try:
                df_head = pd.read_csv(df_p, nrows=10)
                # Only keep numeric/bool columns; always omit touching_* (string TrackID lists)
                valid_cols = [
                    c for c in df_head.columns
                    if (pd.api.types.is_numeric_dtype(df_head[c]) or pd.api.types.is_bool_dtype(df_head[c]))
                    and not c.startswith("touching_")
                ]
            except Exception:
                valid_cols = []
        
        # Filter groups: only show features that actually exist in the CSV and are numeric/bool
        if valid_cols:
            expanded = {}
            for g, patterns in groups.items():
                matches = []
                for p in patterns:
                    if '*' in p or '?' in p:
                        # Only add columns that actually match — drop pattern silently if no match
                        matches.extend([c for c in valid_cols if fnmatch.fnmatchcase(c, p)])
                    elif p in valid_cols:
                        matches.append(p)
                if matches:
                    expanded[g] = sorted(set(matches))
            groups = expanded
        else:
            # No CSV found or no valid columns
            groups = {}
            self._no_csv = True

        self.seed_widget = widgets.IntText(description="Seed", value=int(self._panel_cfg.get("seed", 42)), style={"description_width": "80px"})
        saved_scaling_preset = self._panel_cfg.get("feature_scaling_preset", None)
        scaling_mode = "original_behav3d" if saved_scaling_preset == "original_behav3d" else "custom"
        self.feature_scaling_preset = widgets.Dropdown(
            description="Feature scaling",
            options=[
                ("Custom selected features", "custom"),
                ("Original BEHAV3D", "original_behav3d"),
            ],
            value=scaling_mode,
            style={"description_width": "110px"},
            layout=widgets.Layout(width="360px"),
        )
        
        self._group_rows = {}
        sel_boxes = []
        preset_set = set(self._panel_cfg.get("dtw_features_input", []))
        if not preset_set and self.category == "immune":
            preset_set = {
                "mean_square_displacement",
                "speed",
                "mean_dead_dye",
                f"{self.cell_type}_contact",
                "organoid_contact",
            }
        
        for g, feats in groups.items():
            child_cbs = [widgets.Checkbox(value=(f in preset_set), description=f, indent=True) for f in feats]
            gcb = widgets.Checkbox(value=any(cb.value for cb in child_cbs), indent=False)
            header = widgets.HBox([gcb, widgets.HTML(f"<b>{g}</b>")])
            grid = widgets.GridBox(child_cbs, layout=widgets.Layout(grid_template_columns="repeat(3, max-content)", grid_gap="4px 12px", margin="0 0 0 24px"))
            self._group_rows[g] = {"group_cb": gcb, "child_cbs": child_cbs, "container": widgets.VBox([header, grid])}
            sel_boxes.append(self._group_rows[g]["container"])
            
            def make_h(gn):
                def _h(change):
                    for cb in self._group_rows[gn]["child_cbs"]: cb.value = change["new"]
                return _h
            gcb.observe(make_h(g), names="value")

        self.feature_scaling_preset.observe(self._on_feature_scaling_preset_changed, names="value")

        self.umap_dist = widgets.FloatText(description="UMAP min_dist", value=float(self._panel_cfg.get("umap_min_dist", 0.1)))
        self.umap_neigh = widgets.IntText(description="UMAP n_neighbors", value=int(self._panel_cfg.get("umap_n_neighbors", 15)))
        self.clusters = widgets.IntText(description="# clusters", value=int(self._panel_cfg.get("nr_of_clusters", 5)))
        self.en_min_track_length = widgets.Checkbox(
            description="Min track length",
            value=bool(self._panel_cfg.get("min_track_length_enabled", False)),
            indent=False,
            layout=widgets.Layout(width="160px"),
        )
        self.min_track_length = widgets.BoundedIntText(
            description="timepoints",
            value=max(1, int(self._panel_cfg.get("min_track_length", 1))),
            min=1,
            max=999999,
            style={"description_width": "80px"},
            layout=widgets.Layout(width="190px"),
        )
        self.en_max_track_length = widgets.Checkbox(
            description="Max track length",
            value=bool(self._panel_cfg.get("max_track_length_enabled", False)),
            indent=False,
            layout=widgets.Layout(width="160px"),
        )
        self.max_track_length = widgets.BoundedIntText(
            description="timepoints",
            value=max(1, int(self._panel_cfg.get("max_track_length", 999999))),
            min=1,
            max=999999,
            style={"description_width": "80px"},
            layout=widgets.Layout(width="190px"),
        )
        
        self.btn_run = widgets.Button(description=f"Run {cell_type} analysis", button_style="success", layout=widgets.Layout(width="300px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        group_col_options = list_grouping_candidate_columns(getattr(self.metadata_loader, "metadata", None))
        self.group_cols_select = widgets.SelectMultiple(
            options=group_col_options,
            value=[],
            description="",
            layout=widgets.Layout(width="360px", height="80px"),
            disabled=len(group_col_options) == 0,
        )

        self.backproj_pdf = widgets.Checkbox(
            description="Backprojection PDFs",
            value=bool(self._panel_cfg.get("feature_dtw_backprojection_pdf", True)),
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.backproj_mp4 = widgets.Checkbox(
            description="Backprojection MP4",
            value=bool(self._panel_cfg.get("feature_dtw_backprojection_mp4", False)),
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.backproj_examples_per_cluster = widgets.BoundedIntText(
            description="Examples / cluster",
            value=int(self._panel_cfg.get("feature_dtw_backprojection_examples_per_cluster", 5)),
            min=1,
            max=100,
            style={"description_width": "120px"},
            layout=widgets.Layout(width="230px"),
        )
        self.backproj_num_example_ranks = widgets.BoundedIntText(
            description="Example ranks",
            value=int(self._panel_cfg.get("feature_dtw_backprojection_num_example_ranks", 5)),
            min=1,
            max=100,
            style={"description_width": "100px"},
            layout=widgets.Layout(width="200px"),
        )
        self.btn_backproject = widgets.Button(
            description="Generate track backprojection",
            button_style="success",
            layout=widgets.Layout(width="260px"),
        )
        self.btn_backproject.on_click(self._on_backprojection_clicked)
        self.backprojection_spinner = widgets.HTML(value=spinning_loader)
        self.backprojection_spinner.layout.display = "none"
        self.out_backprojection = widgets.Output()
        self.rename_status = widgets.HTML("<i>Run Feature DTW first to rename clusters.</i>")
        self.rename_rows = widgets.VBox([])
        self.btn_refresh_rename = widgets.Button(
            description="Refresh clusters",
            layout=widgets.Layout(width="150px"),
        )
        self.btn_refresh_rename.on_click(self._on_refresh_feature_dtw_rename_clicked)
        self.btn_rename_feature_dtw = widgets.Button(
            description="Apply cluster names",
            button_style="success",
            layout=widgets.Layout(width="190px"),
        )
        self.btn_rename_feature_dtw.on_click(self._on_rename_feature_dtw_clicked)
        self.rename_spinner = widgets.HTML(value=spinning_loader)
        self.rename_spinner.layout.display = "none"
        self.out_rename_feature_dtw = widgets.Output()
        self.btn_track_proportions = widgets.Button(
            description="Create track proportion plots",
            button_style="info",
            layout=widgets.Layout(width="260px"),
        )
        self.btn_track_proportions.on_click(self._on_track_proportions_clicked)
        self.track_proportions_spinner = widgets.HTML(value=spinning_loader)
        self.track_proportions_spinner.layout.display = "none"
        self.out_track_proportions = widgets.Output()
        self.napari_sample_dd = widgets.Dropdown(
            description="Sample",
            options=[],
            value=None,
            layout=widgets.Layout(width="360px"),
            style={"description_width": "80px"},
        )
        self.napari_workers = widgets.BoundedIntText(
            description="Workers",
            value=int(self._panel_cfg.get("feature_dtw_napari_backprojection_workers", 4)),
            min=1,
            max=32,
            style={"description_width": "80px"},
            layout=widgets.Layout(width="170px"),
        )
        self.btn_refresh_napari_samples = widgets.Button(
            description="Refresh samples",
            layout=widgets.Layout(width="150px"),
        )
        self.btn_refresh_napari_samples.on_click(self._on_refresh_napari_samples_clicked)
        self.btn_open_napari_backprojection = widgets.Button(
            description="Open napari backprojection",
            button_style="success",
            layout=widgets.Layout(width="240px"),
        )
        self.btn_open_napari_backprojection.on_click(self._on_open_napari_backprojection_clicked)
        self.napari_spinner = widgets.HTML(value=spinning_loader)
        self.napari_spinner.layout.display = "none"
        self.out_napari_backprojection = widgets.Output()
        
        # Build feature selector or warning
        if self._no_csv:
            feature_section = widgets.HTML(
                '<div style="color:#856404;background:#fff3cd;border:1px solid #ffc107;'
                'padding:10px;border-radius:4px;margin:6px 0;">'
                '<b>⚠️ No feature CSV found.</b> Run feature extraction first so that '
                'available features can be detected from your data.</div>'
            )
        else:
            feature_section = widgets.VBox(sel_boxes) if sel_boxes else widgets.HTML(
                '<div style="color:#856404;background:#fff3cd;border:1px solid #ffc107;'
                'padding:10px;border-radius:4px;margin:6px 0;">'
                '<b>⚠️ No matching features found in CSV.</b></div>'
            )

        self.ui = widgets.VBox([
            widgets.HTML(f'<div style="font-size:22px;font-weight:700;">{self.cell_type} behavioral analysis</div>'),
            self.seed_widget,
            self.feature_scaling_preset,
            widgets.HTML('<b>Select features for DTW:</b>'),
            feature_section,
            widgets.HBox([self.umap_dist, self.umap_neigh, self.clusters]),
            widgets.HBox([self.en_min_track_length, self.min_track_length, self.en_max_track_length, self.max_track_length]),
            widgets.HTML(
                "<b>Group classification plots by:</b><br>"
                "<span style='color:#555;'>Select 1-2 metadata columns to arrange the cluster "
                "percentage grid by condition. Hold Ctrl/Cmd to select multiple.</span>"
            ),
            self.group_cols_select,
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out,
            widgets.HTML("<b>Track backprojection:</b>"),
            widgets.HBox([self.backproj_pdf, self.backproj_mp4]),
            widgets.HBox([self.backproj_examples_per_cluster, self.backproj_num_example_ranks]),
            widgets.HBox([self.btn_backproject, self.backprojection_spinner]),
            self.out_backprojection,
            widgets.HTML("<b>Rename Feature DTW clusters:</b>"),
            self.rename_status,
            self.rename_rows,
            widgets.HBox([self.btn_refresh_rename, self.btn_rename_feature_dtw, self.rename_spinner]),
            self.out_rename_feature_dtw,
            widgets.HTML(
                "<b>Track-class proportions:</b><br>"
                "<span style='color:#555;'>Regenerate the per-sample cluster-percentage plot (and the "
                "grouped grid) using the current renamed cluster names.</span>"
            ),
            widgets.HTML(
                "<b>Group classification plots by:</b><br>"
                "<span style='color:#555;'>Select 1-2 metadata columns to arrange the cluster "
                "percentage grid by condition. Hold Ctrl/Cmd to select multiple.</span>"
            ),
            self.group_cols_select,
            widgets.HBox([self.btn_track_proportions, self.track_proportions_spinner]),
            self.out_track_proportions,
            widgets.HTML("<b>Napari backprojection:</b>"),
            widgets.HBox([self.napari_sample_dd, self.napari_workers, self.btn_refresh_napari_samples]),
            widgets.HBox([self.btn_open_napari_backprojection, self.napari_spinner]),
            self.out_napari_backprojection,
        ])
        self._on_feature_scaling_preset_changed()
        self._rebuild_feature_dtw_rename_rows()
        self._refresh_feature_dtw_napari_samples()

    def _on_feature_scaling_preset_changed(self, *_):
        enabled = self.feature_scaling_preset.value == "custom"
        for row in self._group_rows.values():
            row["group_cb"].disabled = not enabled
            for cb in row["child_cbs"]:
                cb.disabled = not enabled

    def _feature_dtw_clustered_csv_path(self):
        return Path(
            self.output_dir,
            "analysis",
            self.cell_type,
            "timepoint_feature_dtw",
            f"BEHAV3D_{self.cell_type}_combined_track_features_clustered.csv",
        )

    def _feature_dtw_backprojection_path(self, sample_name):
        return Path(
            self.output_dir,
            "analysis",
            self.cell_type,
            "timepoint_feature_dtw",
            "backprojection",
            f"{sample_name}_{self.cell_type}_feature_dtw_clusters.zarr",
        )

    def _refresh_feature_dtw_napari_samples(self):
        previous = self.napari_sample_dd.value
        samples = []
        clustered_path = self._feature_dtw_clustered_csv_path()
        if clustered_path.exists():
            try:
                df = pd.read_csv(clustered_path, usecols=["sample_name"], low_memory=False)
                samples = sorted(df["sample_name"].dropna().astype(str).unique().tolist())
            except Exception:
                samples = []
        self.napari_sample_dd.options = samples
        if previous in samples:
            self.napari_sample_dd.value = previous
        elif samples:
            self.napari_sample_dd.value = samples[0]
        else:
            self.napari_sample_dd.value = None
        self.btn_open_napari_backprojection.disabled = self.napari_sample_dd.value is None

    def _on_refresh_napari_samples_clicked(self, *_):
        self._refresh_feature_dtw_napari_samples()

    def _feature_dtw_rename_mapping_path(self):
        return Path(
            self.output_dir,
            "analysis",
            self.cell_type,
            "timepoint_feature_dtw",
            f"feature_dtw_cluster_names_{self.cell_type}.yml",
        )

    def _load_feature_dtw_name_mapping(self):
        mapping_path = self._feature_dtw_rename_mapping_path()
        if not mapping_path.exists():
            return {}
        try:
            with mapping_path.open("r", encoding="utf-8") as f:
                data = yaml.safe_load(f) or {}
            mapping = data.get("cluster_names", data)
            if not isinstance(mapping, dict):
                return {}
            return {str(k): str(v) for k, v in mapping.items()}
        except Exception:
            return {}

    def _load_feature_dtw_cluster_colors(self):
        mapping_path = self._feature_dtw_rename_mapping_path()
        if not mapping_path.exists():
            return {}
        try:
            with mapping_path.open("r", encoding="utf-8") as f:
                data = yaml.safe_load(f) or {}
            colors = data.get("cluster_colors", {})
            if not isinstance(colors, dict):
                return {}
            return {str(k): str(v) for k, v in colors.items()}
        except Exception:
            return {}

    def _feature_dtw_output_csv_paths(self):
        outdir = Path(self.output_dir, "analysis", self.cell_type, "timepoint_feature_dtw")
        return [
            outdir / f"BEHAV3D_{self.cell_type}_UMAP_clusters.csv",
            outdir / f"BEHAV3D_{self.cell_type}_combined_track_features_clustered.csv",
            outdir / f"BEHAV3D_{self.cell_type}_UMAP_cluster_percentages.csv",
        ]

    def _feature_dtw_outdir(self):
        return Path(self.output_dir, "analysis", self.cell_type, "timepoint_feature_dtw")

    def _feature_dtw_behavior_proportions_outdir(self):
        return _resolve_dtaidistance_paths(self.output_dir, self.cell_type)["behavior_proportions_outfolder"]

    def _snapshot_feature_dtw_before_renaming(self):
        outdir = self._feature_dtw_outdir()
        before_dir = outdir / "before_renaming"
        before_dir.mkdir(parents=True, exist_ok=True)
        copied = []
        patterns = [
            f"BEHAV3D_{self.cell_type}_UMAP_clusters*",
            f"BEHAV3D_{self.cell_type}_UMAP_cluster_feature_heatmap*",
            f"BEHAV3D_{self.cell_type}_UMAP_cluster_percentages*",
            f"BEHAV3D_{self.cell_type}_combined_track_features_clustered*",
        ]
        for pattern in patterns:
            for src in outdir.glob(pattern):
                if not src.is_file():
                    continue
                dst = before_dir / src.name
                if dst.exists():
                    continue
                shutil.copy2(src, dst)
                copied.append(dst)
        return copied

    def _feature_dtw_plot_info_cols(self, df_umap):
        excluded = {"TrackID", "UMAP1", "UMAP2", "ClusterID", "ClusterName"}
        metadata_cols = [
            c for c in df_umap.columns
            if c.endswith("_line_condition") or c in {"sample_name", "well", "exp_nr"}
        ]
        feature_cols = [
            c for c in df_umap.columns
            if c not in excluded and c not in metadata_cols
        ]
        return list(dict.fromkeys(metadata_cols + feature_cols))

    def _feature_dtw_sample_cols(self, df_umap):
        sample_cols = [c for c in df_umap.columns if c.endswith("_line_condition")]
        for col in ["sample_name", "well", "exp_nr"]:
            if col in df_umap.columns:
                sample_cols.append(col)
        return list(dict.fromkeys(sample_cols))

    def _feature_dtw_cluster_percentage_groups(self, df_umap):
        line_cols = [c for c in df_umap.columns if c.endswith("_line_condition")]
        group_cols_sample = line_cols + ["sample_name", "ClusterID"]
        group_cols_sample_only = line_cols + ["sample_name"]
        missing = [c for c in group_cols_sample if c not in df_umap.columns]
        if missing:
            return None

        df_clust_perc = (
            df_umap
            .groupby(group_cols_sample, observed=True)
            .size()
            .reset_index(name="count")
        )
        sample_totals = (
            df_clust_perc
            .groupby(group_cols_sample_only, observed=True)["count"]
            .sum()
            .reset_index(name="sample_total")
        )
        df_clust_perc = df_clust_perc.merge(sample_totals, how="left", on=group_cols_sample_only)
        if line_cols:
            combo_totals = (
                df_clust_perc
                .groupby(line_cols, observed=True)["count"]
                .sum()
                .reset_index(name="combo_total")
            )
            df_clust_perc = df_clust_perc.merge(combo_totals, how="left", on=line_cols)
        df_clust_perc["percentage"] = df_clust_perc["count"] / df_clust_perc["sample_total"]
        return df_clust_perc

    def _regenerate_feature_dtw_renamed_plots(self):
        outdir = self._feature_dtw_outdir()
        umap_csv = outdir / f"BEHAV3D_{self.cell_type}_UMAP_clusters.csv"
        if not umap_csv.exists():
            raise FileNotFoundError(f"Feature DTW UMAP CSV not found: {umap_csv}")

        df_umap = pd.read_csv(umap_csv, low_memory=False)
        if "ClusterName" not in df_umap.columns:
            raise ValueError("Cannot regenerate renamed plots without a ClusterName column.")
        df_plot = df_umap.copy()
        df_plot["ClusterID"] = pd.Categorical(df_plot["ClusterName"].astype(str))

        sample_cols = self._feature_dtw_sample_cols(df_plot)
        info_cols = self._feature_dtw_plot_info_cols(df_plot)
        cluster_colors = self._load_feature_dtw_cluster_colors()

        cluster_UMAP_path = outdir / f"BEHAV3D_{self.cell_type}_UMAP_clusters.pdf"
        plot_feature_umap(
            df_umap=df_plot,
            info_cols=info_cols,
            sample_cols=sample_cols,
            outpath=cluster_UMAP_path,
            rows_per_page=4,
            nr_cols=2,
            rows_first_img=2,
            figsize=(8.27, 11.69),
            plot_results=True,
            cluster_colors=cluster_colors,
        )

        heatmap_path = outdir / f"BEHAV3D_{self.cell_type}_UMAP_cluster_feature_heatmap.pdf"
        plot_clustering_feature_heatmap(
            df_plot,
            info_cols,
            sample_cols,
            heatmap_path,
            rows_per_page=7,
            nr_cols=2,
            figsize=(8.27, 11.69),
            plot_results=True,
        )

        df_clust_perc = self._feature_dtw_cluster_percentage_groups(df_plot)
        prop_outdir = self._feature_dtw_behavior_proportions_outdir()
        prop_outdir.mkdir(parents=True, exist_ok=True)
        perc_prefix = prop_outdir / f"BEHAV3D_{self.cell_type}_UMAP_cluster_percentages"
        if df_clust_perc is not None:
            perc_csv = perc_prefix.with_suffix(".csv")
            df_clust_perc.to_csv(perc_csv, index=False)
            plot_cluster_percentage_bars(
                df_clust_perc,
                perc_prefix,
                group_cols=list(self.group_cols_select.value) or None,
                plot_results=True,
                cluster_colors=cluster_colors,
            )

        return {
            "umap_pdf": cluster_UMAP_path,
            "heatmap_pdf": heatmap_path,
            "percentages_pdf": perc_prefix.with_suffix(".pdf"),
        }

    def _rebuild_feature_dtw_rename_rows(self):
        self._feature_dtw_name_boxes = {}
        clustered_path = self._feature_dtw_clustered_csv_path()
        if not clustered_path.exists():
            self.rename_rows.children = []
            self.rename_status.value = "<i>Run Feature DTW first to rename clusters.</i>"
            self.btn_rename_feature_dtw.disabled = True
            self.btn_track_proportions.disabled = True
            return

        try:
            df = pd.read_csv(clustered_path, usecols=["ClusterID"], low_memory=False)
        except Exception as exc:
            self.rename_rows.children = []
            self.rename_status.value = f"<i>Could not load Feature DTW clusters: {exc}</i>"
            self.btn_rename_feature_dtw.disabled = True
            self.btn_track_proportions.disabled = True
            return

        clusters = sorted(df["ClusterID"].dropna().astype(str).unique().tolist())
        saved_mapping = self._load_feature_dtw_name_mapping()
        rows = []
        for cluster_id in clusters:
            current_name = saved_mapping.get(str(cluster_id), str(cluster_id))
            txt = widgets.Text(value=current_name, layout=widgets.Layout(width="300px"))
            self._feature_dtw_name_boxes[str(cluster_id)] = txt
            rows.append(
                widgets.HBox(
                    [widgets.Label(str(cluster_id), layout=widgets.Layout(width="110px")), txt],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                )
            )

        self.rename_rows.children = rows
        self.rename_status.value = f"<b>Feature DTW clusters:</b> {len(rows)}"
        self.btn_rename_feature_dtw.disabled = len(rows) == 0
        self.btn_track_proportions.disabled = len(rows) == 0

    def _on_refresh_feature_dtw_rename_clicked(self, *_):
        self._rebuild_feature_dtw_rename_rows()

    def _on_rename_feature_dtw_clicked(self, *_):
        self.btn_rename_feature_dtw.disabled = True
        self.rename_spinner.layout.display = None
        self.out_rename_feature_dtw.clear_output()
        with self.out_rename_feature_dtw:
            try:
                mapping = {}
                for cluster_id, txt in self._feature_dtw_name_boxes.items():
                    name = str(txt.value).strip()
                    mapping[str(cluster_id)] = name if name != "" else str(cluster_id)

                copied = self._snapshot_feature_dtw_before_renaming()
                mapping_path = self._feature_dtw_rename_mapping_path()
                mapping_path.parent.mkdir(parents=True, exist_ok=True)
                existing_colors = self._load_feature_dtw_cluster_colors()
                new_names = sorted(set(str(v) for v in mapping.values()))
                cluster_colors = hash_stable_label_color_map(new_names, colors=existing_colors)
                with mapping_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(
                        {
                            "cell_type": self.cell_type,
                            "cluster_id_column": "ClusterID",
                            "cluster_name_column": "ClusterName",
                            "cluster_names": mapping,
                            "cluster_colors": cluster_colors,
                        },
                        f,
                        sort_keys=False,
                    )

                updated = []
                for csv_path in self._feature_dtw_output_csv_paths():
                    if not csv_path.exists():
                        continue
                    df = pd.read_csv(csv_path, low_memory=False)
                    if "ClusterID" not in df.columns:
                        continue
                    df["ClusterName"] = df["ClusterID"].astype(str).map(mapping).fillna(df["ClusterID"].astype(str))
                    df.to_csv(csv_path, index=False)
                    updated.append(csv_path)

                regenerated = self._regenerate_feature_dtw_renamed_plots()
                self._panel_cfg.update({"feature_dtw_cluster_names": dict(mapping)})
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self._params, f, sort_keys=False)
                before_dir = self._feature_dtw_outdir() / "before_renaming"
                if copied:
                    print(f"Saved original Feature DTW outputs to: {before_dir}")
                else:
                    print(f"Original Feature DTW outputs already present in: {before_dir}")
                print(f"Saved Feature DTW cluster names: {mapping_path}")
                if updated:
                    print("Updated CSVs:")
                    for path in updated:
                        print(f"  - {path}")
                print("Regenerated renamed plots:")
                for path in regenerated.values():
                    print(f"  - {path}")
                self._rebuild_feature_dtw_rename_rows()
            except Exception:
                traceback.print_exc()
            finally:
                self.rename_spinner.layout.display = "none"
                self.btn_rename_feature_dtw.disabled = False

    def _on_track_proportions_clicked(self, *_):
        self.btn_track_proportions.disabled = True
        self.track_proportions_spinner.layout.display = None
        self.out_track_proportions.clear_output()
        with self.out_track_proportions:
            try:
                regenerated = self._regenerate_feature_dtw_renamed_plots()
                print("Regenerated track-class proportion plots:")
                for path in regenerated.values():
                    print(f"  - {path}")
            except Exception:
                traceback.print_exc()
            finally:
                self.track_proportions_spinner.layout.display = "none"
                self.btn_track_proportions.disabled = False

    @staticmethod
    def _adata_like_from_obs(obs):
        class _AdataLike:
            pass

        adata_like = _AdataLike()
        adata_like.obs = obs
        return adata_like

    def _load_feature_dtw_backprojection_data(self):
        clustered_path = self._feature_dtw_clustered_csv_path()
        if not clustered_path.exists():
            raise FileNotFoundError(
                "Feature-DTW clustered tracks were not found. Run the timepoint feature DTW analysis first: "
                f"{clustered_path}"
            )

        df = pd.read_csv(clustered_path, low_memory=False)
        required = [
            "sample_name",
            "TrackID",
            "position_t",
            "ClusterID",
            "pixel_position_x",
            "pixel_position_y",
            "pixel_position_z",
        ]
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(
                "Feature-DTW track backprojection requires missing columns in clustered tracks CSV: "
                + ", ".join(missing)
            )

        df = df.dropna(subset=["sample_name", "TrackID", "position_t", "ClusterID"]).copy()
        if len(df) == 0:
            raise ValueError("Feature-DTW clustered tracks CSV has no usable rows.")
        df["position_t"] = pd.to_numeric(df["position_t"], errors="coerce")
        df = df.dropna(subset=["position_t"])
        df["position_t"] = df["position_t"].astype(int)

        chosen_df = (
            df.groupby(["sample_name", "TrackID"], observed=True)
            .agg(
                ClusterID=("ClusterID", "first"),
                position_t_min=("position_t", "min"),
                position_t_max=("position_t", "max"),
            )
            .reset_index()
        )
        chosen_df = chosen_df.dropna(subset=["ClusterID"]).reset_index(drop=True)
        if len(chosen_df) == 0:
            raise ValueError("No Feature-DTW clustered tracks were available for backprojection.")
        return self._adata_like_from_obs(df), chosen_df

    def _on_backprojection_clicked(self, *_):
        if not bool(self.backproj_pdf.value) and not bool(self.backproj_mp4.value):
            return
        self.btn_backproject.disabled = True
        self.backprojection_spinner.layout.display = None
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                self._panel_cfg.update({
                    "feature_dtw_backprojection_pdf": bool(self.backproj_pdf.value),
                    "feature_dtw_backprojection_mp4": bool(self.backproj_mp4.value),
                    "feature_dtw_backprojection_examples_per_cluster": int(self.backproj_examples_per_cluster.value),
                    "feature_dtw_backprojection_num_example_ranks": int(self.backproj_num_example_ranks.value),
                })
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self._params, f, sort_keys=False)

                adata_like, chosen_df = self._load_feature_dtw_backprojection_data()
                out_dir = Path(
                    self.output_dir,
                    "analysis",
                    self.cell_type,
                    "timepoint_feature_dtw",
                    "example_tracks",
                    "backprojection",
                )
                n_per_cluster = int(self.backproj_examples_per_cluster.value)
                num_example_ranks = int(self.backproj_num_example_ranks.value)

                if bool(self.backproj_pdf.value):
                    out_pdf = save_exemplar_statebar_backprojection_pdf(
                        adata_full=adata_like,
                        output_dir=self.output_dir,
                        cell_type=self.cell_type,
                        out_dir=out_dir,
                        chosen_df=chosen_df,
                        n_per_cluster=n_per_cluster,
                        sample_key="sample_name",
                        track_key="TrackID",
                        time_key="position_t",
                        state_key="ClusterID",
                        cluster_key="ClusterID",
                        tmin_key="position_t_min",
                        tmax_key="position_t_max",
                        seed=int(self.seed_widget.value),
                        examples_per_cluster=n_per_cluster,
                        num_example_ranks=num_example_ranks,
                        verbose=True,
                    )
                    print(
                        "Saved backprojection PDFs: "
                        f"clusters={len(out_pdf.get('pdf_paths_by_cluster', {}))}, "
                        f"examples={len(out_pdf.get('pdf_paths_by_example_rank', {}))}"
                    )

                if bool(self.backproj_mp4.value):
                    out_mp4 = save_exemplar_statebar_backprojection_video_per_cluster(
                        adata_full=adata_like,
                        output_dir=self.output_dir,
                        cell_type=self.cell_type,
                        out_dir=out_dir,
                        chosen_df=chosen_df,
                        n_per_cluster=n_per_cluster,
                        sample_key="sample_name",
                        track_key="TrackID",
                        time_key="position_t",
                        state_key="ClusterID",
                        cluster_key="ClusterID",
                        tmin_key="position_t_min",
                        tmax_key="position_t_max",
                        seed=int(self.seed_widget.value),
                        examples_per_cluster=n_per_cluster,
                        num_example_ranks=num_example_ranks,
                        verbose=True,
                    )
                    print(
                        "Saved backprojection MP4: "
                        f"clusters={len(out_mp4.get('video_paths_by_cluster', {}))}, "
                        f"examples={len(out_mp4.get('video_paths_by_example_rank', {}))}"
                    )
                print(f"Output: {out_dir}")
            except Exception:
                traceback.print_exc()
            finally:
                self.backprojection_spinner.layout.display = "none"
                self.btn_backproject.disabled = False

    def _show_feature_dtw_napari_backprojection(self, sample_name):
        sample_name = str(sample_name).strip()
        if sample_name == "":
            raise ValueError("Please select a sample name.")

        clustered_path = self._feature_dtw_clustered_csv_path()
        if not clustered_path.exists():
            raise FileNotFoundError(
                "Feature-DTW clustered tracks were not found. Run Feature DTW first: "
                f"{clustered_path}"
            )

        df_tracks = pd.read_csv(clustered_path, low_memory=False)
        df_tracks = df_tracks[df_tracks["sample_name"].astype(str) == sample_name].copy()
        if len(df_tracks) == 0:
            raise ValueError(f"No Feature-DTW clustered rows found for sample '{sample_name}'.")

        required = [
            "sample_name",
            "TrackID",
            "position_t",
            "position_z",
            "position_y",
            "position_x",
            "ClusterID",
        ]
        missing = [col for col in required if col not in df_tracks.columns]
        if missing:
            raise ValueError(
                "Feature-DTW napari backprojection requires missing columns: "
                + ", ".join(missing)
            )

        raw_path = _resolve_raw_image_path(
            output_dir=self.output_dir,
            sample_name=sample_name,
            verbose=True,
        )
        if raw_path is None or not Path(raw_path).exists():
            raise FileNotFoundError(f"Could not find raw image for sample '{sample_name}'.")

        tracked_path = _resolve_tracked_image_path(
            output_dir=self.output_dir,
            sample_name=sample_name,
            cell_type=self.cell_type,
            verbose=True,
        )
        if tracked_path is None or not Path(tracked_path).exists():
            raise FileNotFoundError(
                f"Could not find tracked image for sample '{sample_name}' and cell_type '{self.cell_type}'."
            )

        raw_img = load_image(raw_path)
        tracked_img = load_image(tracked_path)
        df_tracks = df_tracks.dropna(subset=["TrackID", "position_t", "ClusterID"]).copy()
        df_tracks["TrackID"] = pd.to_numeric(df_tracks["TrackID"], errors="coerce")
        df_tracks["position_t"] = pd.to_numeric(df_tracks["position_t"], errors="coerce")
        df_tracks["ClusterID"] = pd.to_numeric(df_tracks["ClusterID"], errors="coerce")
        df_tracks = df_tracks.dropna(subset=["TrackID", "position_t", "ClusterID"]).copy()
        if len(df_tracks) == 0:
            raise ValueError(f"No usable Feature-DTW clustered rows found for sample '{sample_name}'.")
        df_tracks["TrackID"] = df_tracks["TrackID"].astype(np.int64)
        df_tracks["position_t"] = df_tracks["position_t"].astype(np.int64)
        df_tracks["ClusterID"] = df_tracks["ClusterID"].astype(np.int64)

        filtered_track_img = _filter_track_image_to_ids(
            tracked_img,
            df_tracks["TrackID"].dropna().unique(),
        )
        backproj_path = self._feature_dtw_backprojection_path(sample_name)
        backproj_path.parent.mkdir(parents=True, exist_ok=True)
        backprojected = backproject_columns(
            track_img=filtered_track_img,
            zarr_outpath=backproj_path,
            df_tracks_clustered=df_tracks,
            columns=["ClusterID"],
            mode="summary",
            track_col="TrackID",
            time_col="position_t",
            background_value=0,
            n_workers=max(1, int(self.napari_workers.value)),
            prefer_parallel=True,
            tracked_img_path=tracked_path,
            verbose=True,
            overwrite=True,
        )
        backprojected["ClusterID"]["type"] = "label"

        backproject_data = {
            "raw_data": {"img": raw_img, "type": "image"},
            "TrackID": {"img": filtered_track_img, "type": "label"},
            "ClusterID": backprojected["ClusterID"],
        }
        viewer = view_napari(
            backproject_data=backproject_data,
            df_tracks_full=df_tracks,
            df_tracks_clustered=df_tracks,
            cell_type=self.cell_type,
            elsize=_sample_voxel_size(getattr(self.metadata_loader, "metadata", None), sample_name),
            split_raw_channels=True,
            run=False,
        )
        print(
            "Opened Feature-DTW napari backprojection for sample "
            f"'{sample_name}' with raw='{Path(raw_path).name}', tracked='{Path(tracked_path).name}', "
            f"clusters='{backproj_path.name}'."
        )
        import napari

        napari.run()
        return viewer

    def _on_open_napari_backprojection_clicked(self, *_):
        self.btn_open_napari_backprojection.disabled = True
        self.napari_spinner.layout.display = None
        self.out_napari_backprojection.clear_output()
        with self.out_napari_backprojection:
            try:
                self._panel_cfg.update({
                    "feature_dtw_napari_backprojection_workers": int(self.napari_workers.value),
                })
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self._params, f, sort_keys=False)
                self._show_feature_dtw_napari_backprojection(self.napari_sample_dd.value)
            except Exception:
                traceback.print_exc()
            finally:
                self.napari_spinner.layout.display = "none"
                self.btn_open_napari_backprojection.disabled = self.napari_sample_dd.value is None
        
    def _on_run_clicked(self, *_):
        sel = [cb.description for r in self._group_rows.values() for cb in r["child_cbs"] if cb.value]
        scaling_preset = None if self.feature_scaling_preset.value == "custom" else self.feature_scaling_preset.value
        if scaling_preset is None and not sel: return
        min_track_length = int(self.min_track_length.value) if self.en_min_track_length.value else None
        max_track_length = int(self.max_track_length.value) if self.en_max_track_length.value else None
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                self._panel_cfg.update({"seed": int(self.seed_widget.value), "umap_min_dist": float(self.umap_dist.value), "umap_n_neighbors": int(self.umap_neigh.value), "nr_of_clusters": int(self.clusters.value), "dtw_features_input": sel, "feature_scaling_preset": scaling_preset, "min_track_length_enabled": bool(self.en_min_track_length.value), "min_track_length": int(self.min_track_length.value), "max_track_length_enabled": bool(self.en_max_track_length.value), "max_track_length": int(self.max_track_length.value)})
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self._params, f, sort_keys=False)
                run_tcell_analysis(cell_type=self.cell_type, output_dir=self.output_dir, df_tracks_path=str(self.df_tracks_path), columns_to_use=sel, columns_to_normalize=sel, umap_minimal_distance=float(self.umap_dist.value), umap_n_neighbors=int(self.umap_neigh.value), nr_of_clusters=int(self.clusters.value), plot_results=True, seed=int(self.seed_widget.value), output_subdir_name="timepoint_feature_dtw", feature_scaling_preset=scaling_preset, min_track_length=min_track_length, max_track_length=max_track_length, cluster_percentage_group_by=list(self.group_cols_select.value) or None)
                self._rebuild_feature_dtw_rename_rows()
                self._refresh_feature_dtw_napari_samples()
                print(f"✅ {self.cell_type} analysis complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False
