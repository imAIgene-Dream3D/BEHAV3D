import traceback
from copy import deepcopy
import math
from pathlib import Path
from time import perf_counter

import ipywidgets as widgets
import pandas as pd
import scanpy as sc
import yaml

from behav3d.analysis.behavior.general import relabel_cluster_ids
from behav3d.analysis.behavior.state.legacy_clustering import (
    apply_state_classifiers_to_full_dataset,
    build_identity_cluster_mapping,
    load_state_classifier_artifact,
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    show_behavioral_state_backprojection,
)
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.core.utils import expand_column_patterns
from behav3d.widgets.utils import PathPicker, behav3d_calculated_features, spinning_loader


def normalize_binary_value(value, tol=1e-9):
    """Return 0/1 if ``value`` normalizes to a boolean, else ``None``.

    Accepts native bools, the strings true/false/t/f, and numeric 0/1.
    Shared by the notebook and napari binary-column detectors so both use
    identical semantics.
    """
    if isinstance(value, (bool, pd.BooleanDtype)):
        return 1 if bool(value) else 0
    if isinstance(value, str):
        sval = value.strip().lower()
        if sval in {"true", "t"}:
            return 1
        if sval in {"false", "f"}:
            return 0
        try:
            value = float(sval)
        except Exception:
            return None
    try:
        fval = float(value)
    except Exception:
        return None
    if not math.isfinite(fval):
        return None
    if abs(fval - 0.0) <= tol:
        return 0
    if abs(fval - 1.0) <= tol:
        return 1
    return None


def detect_binary_columns_from_csv(csv_path, cols, chunksize=50000):
    """Value-based binary-column detection over the *full* CSV.

    A column is binary only if every non-NA value normalizes (via
    :func:`normalize_binary_value`) to ``{0, 1}``. This replaces fragile
    dtype sniffing over a small row sample, which mis-classifies numeric
    columns as binary whenever the sampled rows happen to be NaN/blank.
    """
    if csv_path is None or len(cols) == 0:
        return []

    states = {str(c): {"seen": set(), "invalid": False} for c in cols}
    try:
        for chunk in pd.read_csv(csv_path, usecols=cols, chunksize=chunksize, low_memory=False):
            for col in cols:
                st = states[str(col)]
                if st["invalid"]:
                    continue
                series = chunk[col].dropna()
                if len(series) == 0:
                    continue
                unique_vals = pd.unique(series)
                for raw in unique_vals:
                    norm = normalize_binary_value(raw)
                    if norm is None:
                        st["invalid"] = True
                        break
                    st["seen"].add(int(norm))
                    if len(st["seen"]) > 2:
                        st["invalid"] = True
                        break
    except Exception:
        return []

    return sorted(
        [
            c
            for c in cols
            if (not states[str(c)]["invalid"]) and (len(states[str(c)]["seen"]) > 0)
        ]
    )


def _winfo(prefix, message):
    print(f"[{prefix}] INFO {message}")


def _ordered_unique(values):
    out = []
    seen = set()
    for value in list(values or []):
        text = str(value).strip()
        if text == "" or text in seen:
            continue
        out.append(text)
        seen.add(text)
    return out


class BaseStateClassificationPanel:
    """
    Compact state-classification panel for notebook usage.

    Workflow (accordion sections):
    1) Apply
    2) Clustering
    3) Rename primary dynamic state clusters
    4) Rename clusters assigned to binary groups
    5) Train
    """

    def __init__(self, metadata_loader, cell_type=None):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(getattr(metadata_loader, "output_dir", "")).expanduser())
        self.model_adata = None
        self.adata_full = None
        self._available_columns = []
        self._feature_groups = {}
        self._group_rows = {}
        self._intrinsic_name_boxes = {}
        self._full_name_boxes = {}
        self._full_select_boxes = {}
        self._binary_detection_cache = None

        self._descriptive_feature_options = (
            "mean",
            "median",
            "std",
            "net_displacement",
            "straightness",
            "mean_square_displacement",
        )

        self._cell_types = self._detect_cell_types()
        if len(self._cell_types) == 0:
            self._cell_types = ["tcell"]
        initial_cell_type = (
            str(cell_type)
            if cell_type is not None
            else (self._cell_types[0] if len(self._cell_types) > 0 else "tcell")
        )
        if initial_cell_type not in self._cell_types:
            self._cell_types.append(initial_cell_type)

        self.cell_type_dd = widgets.Dropdown(
            description="Cell type",
            options=list(self._cell_types),
            value=initial_cell_type,
            layout=widgets.Layout(width="260px"),
            style={"description_width": "90px"},
        )
        self.refresh_btn = widgets.Button(
            description="Refresh",
            icon="refresh",
            layout=widgets.Layout(width="110px"),
        )
        self.refresh_spinner = widgets.HTML(value=spinning_loader)
        self.refresh_spinner.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.output_dir_html = widgets.HTML("")

        self._build_apply_section()
        self._build_clustering_section()
        self._build_intrinsic_rename_section()
        self._build_full_rename_section()
        self._build_backprojection_section()

        self._build_steps()

        self.ui = widgets.VBox(
            [
                widgets.HTML('<div style="font-size:20px;font-weight:700;">State Classification</div>'),
                widgets.HBox(
                    [self.cell_type_dd, self.refresh_btn, self.refresh_spinner],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                ),
                self.output_dir_html,
                self.status_html,
                self.steps,
            ]
        )

        self.cell_type_dd.observe(self._on_cell_type_changed, names="value")
        self.refresh_btn.on_click(self._on_refresh_clicked)
        self.apply_full_pkl_picker.text.observe(self._on_apply_path_changed, names="value")
        self.apply_intrinsic_pkl_picker.text.observe(self._on_apply_path_changed, names="value")
        self._refresh_context()

    def _build_steps(self):
        step_defs = [
            ("Primary dynamic state clustering (based on continuous features)", self.clustering_section),
            ("Rename primary dynamic state clusters", self.rename_intrinsic_section),
            ("Rename clusters assigned to binary groups", self.rename_full_section),
            ("Apply classification", self.apply_section),
            ("Backprojection", self.backprojection_section),
        ]
        self._step_accordions = []
        for title, section in step_defs:
            acc = widgets.Accordion(children=[section], selected_index=None)
            acc.set_title(0, title)
            self._step_accordions.append(acc)
        self.steps = widgets.VBox(self._step_accordions)

    def _collapse_all_steps(self):
        for acc in self._step_accordions:
            acc.selected_index = None

    def _open_step(self, index):
        if index is None:
            return
        if index < 0 or index >= len(self._step_accordions):
            return
        self._step_accordions[index].selected_index = 0

    def _build_apply_section(self):
        self.apply_full_pkl_picker = PathPicker(
            mode="file",
            start_dir=self.output_dir or ".",
            default="",
            description="Full PKL",
            placeholder="Path to full classification .pkl (required)",
            width="100%",
        )
        self.apply_intrinsic_pkl_picker = PathPicker(
            mode="file",
            start_dir=self.output_dir or ".",
            default="",
            description="Intrinsic PKL",
            placeholder="Path to intrinsic classification .pkl (optional)",
            width="100%",
        )
        self.apply_full_pkl_picker.filter_pattern = "*.pkl"
        self.apply_intrinsic_pkl_picker.filter_pattern = "*.pkl"
        self.apply_default_paths_html = widgets.HTML("")
        self.btn_apply = widgets.Button(
            description="Apply classification",
            button_style="success",
            layout=widgets.Layout(width="170px"),
        )
        self.btn_apply.on_click(self._on_apply_clicked)
        self.apply_spinner = widgets.HTML(value=spinning_loader)
        self.apply_spinner.layout.display = "none"
        self.out_apply = widgets.Output()
        self.apply_section = widgets.VBox(
            [
                self.apply_full_pkl_picker,
                self.apply_intrinsic_pkl_picker,
                self.apply_default_paths_html,
                widgets.HBox([self.btn_apply, self.apply_spinner]),
                self.out_apply,
            ]
        )

    def _build_clustering_section(self):
        self.feature_groups_status = widgets.HTML("<i>No features loaded yet.</i>")
        self.feature_groups_box = widgets.VBox([])
        self.selected_features_box = widgets.Select(
            options=[],
            rows=10,
            layout=widgets.Layout(width="520px", height="180px"),
        )
        self.binary_group_status = widgets.HTML("<i>No binary columns detected yet.</i>")
        self.binary_group_checks_box = widgets.VBox([])

        self._binary_group_checkboxes = {}

        self.window_size = widgets.IntText(
            description="Trailing window",
            value=5,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.min_spacing = widgets.Text(
            description="Min spacing between choses timepoints",
            value="",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="360px"),
        )
        self.max_samples = widgets.Text(
            description="Max samples",
            value="",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.n_neighbors = widgets.IntText(
            description="Neighbors",
            value=60,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=0.1,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.resolution = widgets.FloatText(
            description="Resolution",
            value=0.2,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.pca_var_selection = widgets.FloatText(
            description="PCA var",
            value=0.95,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="180px"),
        )
        self.clustering_method = widgets.Dropdown(
            description="Method",
            options=["leiden", "kmeans"],
            value="leiden",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.lower_quantile_cap = widgets.Text(
            description="Lower quantile",
            value="",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.upper_quantile_cap = widgets.Text(
            description="Upper quantile",
            value="0.99",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="220px"),
        )
        self.incomplete_window_policy = widgets.Dropdown(
            description="Partial window policy",
            options=["partial", "drop"],
            value="partial",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="320px"),
        )
        self.random_state = widgets.IntText(
            description="Seed",
            value=123,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="150px"),
        )
        self.reuse_prepared_dataset = widgets.Checkbox(
            description="Reuse prepared dataset",
            value=True,
            indent=False,
        )
        self.describe_window_feature_cbs = {
            "mean": widgets.Checkbox(description="mean", value=True, indent=False),
            "median": widgets.Checkbox(description="median", value=True, indent=False),
            "std": widgets.Checkbox(description="std", value=True, indent=False),
        }
        self.additional_window_feature_cbs = {
            "net_displacement": widgets.Checkbox(description="net_displacement", value=True, indent=False),
            "straightness": widgets.Checkbox(description="straightness", value=True, indent=False),
            "mean_square_displacement": widgets.Checkbox(
                description="mean_square_displacement",
                value=True,
                indent=False,
            ),
        }

        self.btn_cluster = widgets.Button(
            description="Run primary dynamic state clustering",
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
                widgets.HTML("<b>selected features</b>"),
                self.selected_features_box,
            ]
        )
        binary_group_block = widgets.VBox(
            [
                widgets.HTML("<b>assign binary group labels</b>"),
                self.binary_group_status,
                self.binary_group_checks_box,
            ]
        )
        self.feature_selection_subaccordion = widgets.Accordion(
            children=[feature_select_block, binary_group_block],
            selected_index=None,
        )
        self.feature_selection_subaccordion.set_title(0, "available features")
        self.feature_selection_subaccordion.set_title(1, "assign binary group labels")
        feature_selection_main_block = widgets.VBox([self.feature_selection_subaccordion])

        describe_window_features_block = widgets.VBox(
            [
                widgets.HTML("<b>window descriptors</b>"),
                widgets.VBox(list(self.describe_window_feature_cbs.values())),
            ]
        )
        additional_window_features_block = widgets.VBox(
            [
                widgets.HTML("<b>window_based additional features</b>"),
                widgets.VBox(list(self.additional_window_feature_cbs.values())),
            ]
        )

        windowed_timepoint_processing_block = widgets.VBox(
            [
                describe_window_features_block,
                additional_window_features_block,
                widgets.HTML("<b>window processing parameters</b>"),
                widgets.HBox(
                    [self.window_size, self.incomplete_window_policy],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                widgets.HTML("<b>feature value capping</b>"),
                widgets.HBox(
                    [self.lower_quantile_cap, self.upper_quantile_cap],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
            ]
        )

        clustering_block = widgets.VBox(
            [
                widgets.HBox(
                    [self.min_spacing, self.max_samples, self.n_neighbors, self.min_dist],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                widgets.HBox(
                    [self.resolution, self.pca_var_selection, self.clustering_method],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                widgets.HBox([self.random_state], layout=widgets.Layout(flex_flow="row wrap", gap="8px")),
                self.reuse_prepared_dataset,
            ]
        )

        self.clustering_subaccordion = widgets.Accordion(
            children=[windowed_timepoint_processing_block, clustering_block],
            selected_index=None,
        )
        self.clustering_subaccordion.set_title(0, "Windowed timepoint processing")
        self.clustering_subaccordion.set_title(1, "Clustering")
        clustering_main_block = widgets.VBox([self.clustering_subaccordion])

        self.clustering_main_accordion = widgets.Accordion(
            children=[feature_selection_main_block, clustering_main_block],
            selected_index=None,
        )
        self.clustering_main_accordion.set_title(0, "Feature selection")
        self.clustering_main_accordion.set_title(1, "Clustering")

        self.clustering_section = widgets.VBox(
            [
                self.clustering_main_accordion,
                widgets.HBox([self.btn_cluster, self.cluster_spinner]),
                self.out_cluster,
            ]
        )

    def _selected_descriptive_features(self):
        selected = []
        for feature_name, cb in self.describe_window_feature_cbs.items():
            if bool(cb.value):
                selected.append(str(feature_name))
        for feature_name, cb in self.additional_window_feature_cbs.items():
            if bool(cb.value):
                selected.append(str(feature_name))
        return selected

    def _build_intrinsic_rename_section(self):
        self.rename_intrinsic_status = widgets.HTML("<i>Run clustering or load existing model first.</i>")
        self.rename_intrinsic_rows = widgets.VBox([])
        self.btn_rename_intrinsic = widgets.Button(
            description="Rename primary dynamic state clusters",
            button_style="warning",
            layout=widgets.Layout(width="320px"),
        )
        self.btn_rename_intrinsic.on_click(self._on_rename_intrinsic_clicked)
        self.rename_intrinsic_spinner = widgets.HTML(value=spinning_loader)
        self.rename_intrinsic_spinner.layout.display = "none"
        self.out_rename_intrinsic = widgets.Output()
        self.rename_intrinsic_section = widgets.VBox(
            [
                self.rename_intrinsic_status,
                self.rename_intrinsic_rows,
                widgets.HBox([self.btn_rename_intrinsic, self.rename_intrinsic_spinner]),
                self.out_rename_intrinsic,
            ]
        )

    def _build_full_rename_section(self):
        self.rename_full_status = widgets.HTML("<i>Rename primary dynamic state clusters first.</i>")
        self.rename_full_rows = widgets.VBox([])
        self.rename_full_rows.layout = widgets.Layout(width="560px")
        self.full_combine_name = widgets.Text(
            description="Combine to",
            value="",
            placeholder="New combined name",
            style={"description_width": "90px"},
            layout=widgets.Layout(width="360px"),
        )
        self.btn_combine_full = widgets.Button(
            description="combine",
            button_style="info",
            layout=widgets.Layout(width="110px"),
        )
        self.btn_combine_full.on_click(self._on_combine_full_clicked)
        self.combine_full_spinner = widgets.HTML(value=spinning_loader)
        self.combine_full_spinner.layout.display = "none"
        self.full_combine_box = widgets.VBox(
            [
                widgets.HTML("<b>Combine selected</b>"),
                self.full_combine_name,
                widgets.HBox([self.btn_combine_full, self.combine_full_spinner]),
            ],
            layout=widgets.Layout(width="390px", align_items="flex-start"),
        )
        self.btn_rename_full = widgets.Button(
            description="Rename clusters assigned to binary groups",
            button_style="warning",
            layout=widgets.Layout(width="320px"),
        )
        self.btn_rename_full.on_click(self._on_rename_full_clicked)
        self.rename_full_spinner = widgets.HTML(value=spinning_loader)
        self.rename_full_spinner.layout.display = "none"
        self.out_rename_full = widgets.Output()
        self.rename_full_section = widgets.VBox(
            [
                self.rename_full_status,
                widgets.HBox(
                    [self.rename_full_rows, self.full_combine_box],
                    layout=widgets.Layout(align_items="flex-start", gap="14px"),
                ),
                widgets.HBox([self.btn_rename_full, self.rename_full_spinner]),
                self.out_rename_full,
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
                self.backproj_sample_dd,
                widgets.HBox([self.btn_open_backprojection, self.backprojection_spinner]),
                self.out_backprojection,
            ]
        )

    def _detect_cell_types(self):
        md = getattr(self.metadata_loader, "metadata", None)
        cell_types = []
        if md is not None:
            try:
                cell_types.extend(detect_immune_cell_types_from_metadata(md))
                cell_types.extend(detect_organoid_types_from_metadata(md))
                cell_types.extend(detect_other_cell_types_from_metadata(md))
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
        return _ordered_unique(cell_types)

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

    def _refresh_backprojection_samples(self):
        current = self.backproj_sample_dd.value
        sample_names = self._detect_sample_names()
        if len(sample_names) == 0:
            self.backproj_sample_dd.options = []
            self.backproj_sample_dd.value = None
            self.backprojection_status.value = "<i>No samples detected for backprojection.</i>"
            return

        self.backproj_sample_dd.options = sample_names
        if current in sample_names:
            self.backproj_sample_dd.value = current
        else:
            self.backproj_sample_dd.value = sample_names[0]
        self.backprojection_status.value = f"<b>Available samples:</b> {len(sample_names)}"

    def _current_cell_type(self):
        return str(self.cell_type_dd.value).strip()

    def _state_root(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return Path(self.output_dir, "analysis", ct, "behavioral_states")

    def _model_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return Path(
            self.output_dir,
            "analysis",
            ct,
            "behavioral_states",
            "processing",
            f"BEHAV3D_{ct}_behavioral_states_modeldata.h5ad",
        )

    def _default_classifier_paths(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        base = Path(self.output_dir, "analysis", ct, "behavioral_states", "processing")
        return {
            "intrinsic": base
            / "intrinsic_behavioral_classification"
            / f"intrinsic_state_classification_random_forest_{ct}.pkl",
            "full": base
            / "full_behavioral_classification"
            / f"state_classification_random_forest_{ct}.pkl",
        }

    def _refresh_apply_default_paths(self):
        if self.output_dir is None or str(self.output_dir).strip() == "":
            self.apply_default_paths_html.value = ""
            return

        paths = self._default_classifier_paths()
        full_path = paths["full"]
        intrinsic_path = paths["intrinsic"]
        full_exists = full_path.exists()
        intrinsic_exists = intrinsic_path.exists()

        if str(self.apply_full_pkl_picker.value).strip() == "" and full_exists:
            self.apply_full_pkl_picker.value = str(full_path)
        if str(self.apply_intrinsic_pkl_picker.value).strip() == "" and intrinsic_exists:
            self.apply_intrinsic_pkl_picker.value = str(intrinsic_path)

        if full_exists or intrinsic_exists:
            self.apply_default_paths_html.value = (
                "<b style='color:#080;'>Default classifier path(s) detected and prefilled.</b>"
            )
        else:
            self.apply_default_paths_html.value = ""

    def _resolve_track_features_csv(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        base = Path(self.output_dir, "analysis", ct, "track_features")
        filtered_csv = base / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        combined_csv = base / f"BEHAV3D_{ct}_combined_track_features.csv"
        if filtered_csv.exists():
            return filtered_csv
        if combined_csv.exists():
            return combined_csv
        return None

    def _set_widgets_disabled(self, widgets_list, state):
        for w in widgets_list:
            if hasattr(w, "disabled"):
                w.disabled = bool(state)

    def _parse_optional_int(self, text_value):
        value = str(text_value).strip()
        if value == "":
            return None
        return int(value)

    def _parse_optional_float(self, text_value):
        value = str(text_value).strip()
        if value == "":
            return None
        return float(value)

    def _normalize_binary_value(self, value, tol=1e-9):
        return normalize_binary_value(value, tol=tol)

    def _detect_binary_columns_from_csv(self, csv_path, cols, chunksize=50000):
        if csv_path is None or len(cols) == 0:
            return []

        cache_key = (
            str(csv_path),
            float(csv_path.stat().st_mtime) if Path(csv_path).exists() else -1.0,
            tuple(cols),
        )
        if isinstance(self._binary_detection_cache, dict):
            if self._binary_detection_cache.get("key", None) == cache_key:
                return list(self._binary_detection_cache.get("binary_cols", []))

        binary_cols = detect_binary_columns_from_csv(csv_path, cols, chunksize=chunksize)
        self._binary_detection_cache = {
            "key": cache_key,
            "binary_cols": list(binary_cols),
        }
        return binary_cols

    def _excluded_non_behavior_columns(self, cols):
        col_set = {str(c) for c in cols}
        excluded = {
            "sample_name",
            "TrackID",
            "sub_TrackID",
            "position_t",
            "position_x",
            "position_y",
            "position_z",
            "segment_id",
            "lineage_id",
            "frame",
            "t",
            "interpolated",
            "exp_nr",
            "well",
            "origin_cell_type",
            "origin_TrackID",
        }

        md = getattr(self.metadata_loader, "metadata", None)
        if isinstance(md, pd.DataFrame):
            excluded.update({str(c) for c in md.columns})

        for c in col_set:
            lc = str(c).lower()
            if lc.endswith("_line_condition"):
                excluded.add(c)
            if lc.endswith("_tracks_csv_path"):
                excluded.add(c)
            if lc.endswith("_segments_image_path"):
                excluded.add(c)
            if lc.endswith("_tracks_image_path"):
                excluded.add(c)
            if lc.endswith("_raw_image_path"):
                excluded.add(c)
            if lc.endswith("_dead_mask_path"):
                excluded.add(c)
            if lc.startswith("channel_") and lc.endswith("_label"):
                excluded.add(c)

        return excluded.intersection(col_set)

    def _panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        panel = params.setdefault("state_classification", {})
        return panel.setdefault(self._current_cell_type(), {})

    def _save_panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        cfg_path = getattr(self.metadata_loader, "behav3d_parameters_path", None)
        if not isinstance(params, dict) or cfg_path is None:
            return
        try:
            with Path(cfg_path).open("w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        except Exception:
            pass

    def _apply_cfg_defaults(self):
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        self.window_size.value = int(cfg.get("window_size", self.window_size.value))
        self.min_spacing.value = str(cfg.get("min_spacing", self.min_spacing.value or ""))
        self.max_samples.value = str(cfg.get("max_samples", self.max_samples.value or ""))
        self.n_neighbors.value = int(cfg.get("n_neighbors", self.n_neighbors.value))
        self.min_dist.value = float(cfg.get("min_dist", self.min_dist.value))
        self.resolution.value = float(cfg.get("resolution", self.resolution.value))
        self.pca_var_selection.value = float(cfg.get("pca_var_selection", self.pca_var_selection.value))
        self.clustering_method.value = str(cfg.get("clustering_method", self.clustering_method.value))
        self.lower_quantile_cap.value = str(cfg.get("lower_quantile_cap", self.lower_quantile_cap.value or ""))
        self.upper_quantile_cap.value = str(cfg.get("upper_quantile_cap", self.upper_quantile_cap.value or "0.99"))
        self.incomplete_window_policy.value = str(
            cfg.get("incomplete_window_policy", self.incomplete_window_policy.value)
        )
        self.random_state.value = int(cfg.get("random_state", self.random_state.value))
        self.reuse_prepared_dataset.value = bool(cfg.get("reuse_prepared_dataset", self.reuse_prepared_dataset.value))
        self.apply_full_pkl_picker.value = str(cfg.get("apply_full_classifier_path", self.apply_full_pkl_picker.value))
        self.apply_intrinsic_pkl_picker.value = str(
            cfg.get("apply_intrinsic_classifier_path", self.apply_intrinsic_pkl_picker.value)
        )
        saved_desc = cfg.get("descriptive_features", None)
        if isinstance(saved_desc, (list, tuple)) and len(saved_desc) > 0:
            for cb in self.describe_window_feature_cbs.values():
                cb.value = False
            for cb in self.additional_window_feature_cbs.values():
                cb.value = False
            for feat in saved_desc:
                if feat in self.describe_window_feature_cbs:
                    self.describe_window_feature_cbs[str(feat)].value = True
                if feat in self.additional_window_feature_cbs:
                    self.additional_window_feature_cbs[str(feat)].value = True



    def _persist_current_settings(self):
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        cfg.update(
            {
                "window_size": int(self.window_size.value),
                "min_spacing": str(self.min_spacing.value),
                "max_samples": str(self.max_samples.value),
                "n_neighbors": int(self.n_neighbors.value),
                "min_dist": float(self.min_dist.value),
                "resolution": float(self.resolution.value),
                "pca_var_selection": float(self.pca_var_selection.value),
                "clustering_method": str(self.clustering_method.value),
                "lower_quantile_cap": str(self.lower_quantile_cap.value),
                "upper_quantile_cap": str(self.upper_quantile_cap.value),
                "incomplete_window_policy": str(self.incomplete_window_policy.value),
                "random_state": int(self.random_state.value),
                "reuse_prepared_dataset": bool(self.reuse_prepared_dataset.value),
                "apply_full_classifier_path": str(self.apply_full_pkl_picker.value),
                "apply_intrinsic_classifier_path": str(self.apply_intrinsic_pkl_picker.value),
                "descriptive_features": list(self._selected_descriptive_features()),
                "selected_features": self._selected_feature_columns(),
                "binary_features_to_group": self._selected_binary_columns(),
            }
        )
        self._save_panel_cfg()

    def _build_feature_groups(self):
        self._group_rows = {}
        self._feature_groups = {}
        self._binary_group_checkboxes = {}
        cols = list(self._available_columns)
        if len(cols) == 0:
            self.feature_groups_status.value = (
                "<i>No track-features CSV found for this cell type. "
                "Run feature extraction first.</i>"
            )
            self.feature_groups_box.children = []
            self.selected_features_box.options = []
            self.binary_group_status.value = "<i>No binary columns detected yet.</i>"
            self.binary_group_checks_box.children = []
            return

        excluded = self._excluded_non_behavior_columns(cols)
        usable_cols = [c for c in cols if c not in excluded]
        base_groups = deepcopy(behav3d_calculated_features)
        matched = set()
        grouped = {}
        for group_name, patterns in base_groups.items():
            vals = []
            for pat in patterns:
                vals.extend(expand_column_patterns(pat, usable_cols))
            clean_vals = sorted({x for x in vals if x in usable_cols})
            if len(clean_vals) > 0:
                grouped[group_name] = clean_vals
                matched.update(clean_vals)

        other = sorted([c for c in usable_cols if c not in matched])
        if len(other) > 0:
            grouped["other"] = other

        self._feature_groups = grouped
        cfg = self._panel_cfg()
        preselected = set(cfg.get("selected_features", [])) if isinstance(cfg, dict) else set()
        if len(preselected) == 0:
            # fallback to basic sensible defaults when present
            for f in ["speed", "elongation", "sphericity", "extent", "solidity"]:
                if f in usable_cols:
                    preselected.add(f)

        children = []
        titles = []
        for group_name, feats in grouped.items():
            child_cbs = [
                widgets.Checkbox(
                    description=f,
                    value=(f in preselected),
                    indent=False,
                    layout=widgets.Layout(width="230px"),
                )
                for f in feats
            ]
            for cb in child_cbs:
                cb.observe(self._on_feature_checkbox_changed, names="value")
            group_cb = widgets.Checkbox(
                description=f"Select {group_name}",
                value=all(cb.value for cb in child_cbs) if len(child_cbs) > 0 else False,
                indent=False,
            )

            def _make_group_handler(cbs):
                def _handler(change):
                    new_val = bool(change["new"])
                    for cb in cbs:
                        cb.value = new_val

                return _handler

            group_cb.observe(_make_group_handler(child_cbs), names="value")
            grid = widgets.GridBox(
                child_cbs,
                layout=widgets.Layout(
                    grid_template_columns="repeat(3, max-content)",
                    grid_gap="2px 10px",
                ),
            )
            box = widgets.VBox([group_cb, grid])
            children.append(box)
            titles.append(group_name)
            self._group_rows[group_name] = {"group_cb": group_cb, "child_cbs": child_cbs}

        if len(children) == 0:
            self.feature_groups_status.value = "<i>No usable feature columns detected.</i>"
            self.feature_groups_box.children = []
        else:
            fg_acc = widgets.Accordion(children=children, selected_index=None)
            for idx, title in enumerate(titles):
                fg_acc.set_title(idx, title)
            self.feature_groups_status.value = (
                f"<b>Available features:</b> {len(usable_cols)} usable columns "
                f"(excluded metadata/technical: {len(excluded)})"
            )
            self.feature_groups_box.children = [fg_acc]
        self._update_selected_features_box()

        csv_path = self._resolve_track_features_csv()
        binary_input_cols = [c for c in usable_cols if c not in {"interpolated"}]
        binary_candidates = self._detect_binary_columns_from_csv(
            csv_path=csv_path,
            cols=binary_input_cols,
        )
        saved_binary = tuple(cfg.get("binary_features_to_group", [])) if isinstance(cfg, dict) else tuple()
        safe_binary = tuple([x for x in saved_binary if x in binary_candidates])
        if len(safe_binary) == 0:
            safe_binary = tuple(binary_candidates)
        binary_checks = []
        for col in binary_candidates:
            cb = widgets.Checkbox(
                description=str(col),
                value=(str(col) in safe_binary),
                indent=False,
                layout=widgets.Layout(width="500px"),
            )
            cb.observe(self._on_binary_checkbox_changed, names="value")
            self._binary_group_checkboxes[str(col)] = cb
            binary_checks.append(cb)

        if len(binary_checks) == 0:
            self.binary_group_status.value = "<i>No binary columns detected from values (0/1 or true/false).</i>"
            self.binary_group_checks_box.children = []
        else:
            self.binary_group_status.value = (
                f"<b>Detected binary columns:</b> {len(binary_checks)} "
                "(each selected column is treated as an independent binary group; "
                "full labels show observed training combinations only)"
            )
            self.binary_group_checks_box.children = binary_checks

    def _selected_feature_columns(self):
        out = []
        for row in self._group_rows.values():
            for cb in row["child_cbs"]:
                if cb.value:
                    out.append(str(cb.description))
        # preserve order + de-dup
        seen = set()
        uniq = []
        for x in out:
            if x not in seen:
                uniq.append(x)
                seen.add(x)
        return uniq

    def _selected_binary_columns(self):
        return [
            str(col)
            for col, cb in self._binary_group_checkboxes.items()
            if bool(cb.value)
        ]

    def _update_selected_features_box(self):
        self.selected_features_box.options = self._selected_feature_columns()

    def _on_feature_checkbox_changed(self, _):
        self._update_selected_features_box()
        self._refresh_enablement()

    def _on_binary_checkbox_changed(self, _):
        self._refresh_enablement()

    def _load_columns(self):
        csv_path = self._resolve_track_features_csv()
        if csv_path is None:
            self._available_columns = []
            return
        try:
            self._available_columns = pd.read_csv(csv_path, nrows=0).columns.tolist()
        except Exception:
            self._available_columns = []

    def _load_existing_model_if_available(self):
        self.model_adata = None
        p = self._model_adata_path()
        if p.exists():
            try:
                self.model_adata = sc.read_h5ad(p)
            except Exception:
                self.model_adata = None

    def _save_model_adata(self, compression="gzip"):
        if self.model_adata is None:
            return
        p = self._model_adata_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        self.model_adata.write(p, compression=compression)

    def _set_busy(self, button, spinner, busy=True, disable_buttons=None):
        state = bool(busy)
        if button is not None and hasattr(button, "disabled"):
            button.disabled = state
        if spinner is not None:
            spinner.layout.display = None if state else "none"
        if disable_buttons is not None:
            for b in disable_buttons:
                if b is not None and hasattr(b, "disabled"):
                    b.disabled = state

    def _rebuild_intrinsic_rename_rows(self):
        self._intrinsic_name_boxes = {}
        if self.model_adata is None or "intrinsic_behavioral_cluster" not in self.model_adata.obs.columns:
            self.rename_intrinsic_rows.children = []
            self.rename_intrinsic_status.value = "<i>Run clustering or load existing model first.</i>"
            self.btn_rename_intrinsic.disabled = True
            return

        mapping = build_identity_cluster_mapping(
            self.model_adata,
            cluster_col="intrinsic_behavioral_cluster",
        )
        rows = []
        for old_name in mapping.keys():
            txt = widgets.Text(value=str(old_name), layout=widgets.Layout(width="280px"))
            self._intrinsic_name_boxes[str(old_name)] = txt
            rows.append(
                widgets.HBox(
                    [widgets.Label(str(old_name), layout=widgets.Layout(width="230px")), txt],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                )
            )
        self.rename_intrinsic_rows.children = rows
        self.rename_intrinsic_status.value = f"<b>Primary dynamic state clusters:</b> {len(rows)}"
        self.btn_rename_intrinsic.disabled = False

    def _rebuild_full_rename_rows(self):
        self._full_select_boxes = {}
        self._full_name_boxes = {}
        self.full_combine_name.value = ""
        if self.model_adata is None or "full_behavioral_cluster" not in self.model_adata.obs.columns:
            self.rename_full_rows.children = []
            self.rename_full_status.value = "<i>Run primary dynamic state cluster rename first.</i>"
            self.btn_rename_full.disabled = True
            self.btn_combine_full.disabled = True
            self.full_combine_name.disabled = True
            return

        mapping = build_identity_cluster_mapping(
            self.model_adata,
            cluster_col="full_behavioral_cluster",
        )
        rows = []
        for old_name in mapping.keys():
            sel = widgets.Checkbox(
                value=False,
                description="",
                indent=False,
                layout=widgets.Layout(width="26px"),
            )
            txt = widgets.Text(value=str(old_name), layout=widgets.Layout(width="280px"))
            self._full_select_boxes[str(old_name)] = sel
            self._full_name_boxes[str(old_name)] = txt
            rows.append(
                widgets.HBox(
                    [sel, widgets.Label(str(old_name), layout=widgets.Layout(width="200px")), txt],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                )
            )
        self.rename_full_rows.children = rows
        self.rename_full_status.value = (
            f"<b>Clusters assigned to binary groups:</b> {len(rows)} "
            "(observed combinations from training data)"
        )
        self.btn_rename_full.disabled = False

    def _refresh_enablement(self):
        has_cell_type = self._current_cell_type() != ""
        has_features = len(self._selected_feature_columns()) > 0
        has_model = self.model_adata is not None
        has_intrinsic = has_model and ("intrinsic_behavioral_cluster" in self.model_adata.obs.columns)
        has_full = has_model and ("full_behavioral_cluster" in self.model_adata.obs.columns)
        has_full_classifier_input = str(self.apply_full_pkl_picker.value).strip() != ""
        has_backproj_sample = self.backproj_sample_dd.value is not None and len(str(self.backproj_sample_dd.value)) > 0

        self.btn_apply.disabled = not (has_cell_type and has_full_classifier_input)
        self.btn_cluster.disabled = not (has_cell_type and has_features)
        self.btn_rename_intrinsic.disabled = not has_intrinsic
        self.btn_rename_full.disabled = not has_full
        self.btn_combine_full.disabled = not has_full
        self.full_combine_name.disabled = not has_full
        self.btn_open_backprojection.disabled = not (has_cell_type and has_backproj_sample)

        if has_cell_type:
            self.status_html.value = f"<b>Ready:</b> cell_type='{self._current_cell_type()}'"
        else:
            self.status_html.value = "<b>Select a cell type to continue.</b>"

    def _refresh_context(self):
        self.output_dir = str(Path(getattr(self.metadata_loader, "output_dir", "")).expanduser())
        self.output_dir_html.value = f"<b>Output dir:</b> {self.output_dir}"
        if hasattr(self.apply_full_pkl_picker, "_start_dir"):
            self.apply_full_pkl_picker._start_dir = self.output_dir or "."
        if hasattr(self.apply_intrinsic_pkl_picker, "_start_dir"):
            self.apply_intrinsic_pkl_picker._start_dir = self.output_dir or "."

        # Refresh cell-type options from metadata/filesystem while preserving current selection.
        current_value = self._current_cell_type()
        refreshed_types = self._detect_cell_types()
        if len(refreshed_types) == 0:
            refreshed_types = [current_value or "tcell"]
        if current_value not in refreshed_types:
            refreshed_types.append(current_value)
        new_options = _ordered_unique(refreshed_types)
        if len(new_options) == 0:
            new_options = ["tcell"]
        if list(self.cell_type_dd.options) != list(new_options):
            self.cell_type_dd.options = new_options
        if self.cell_type_dd.value not in self.cell_type_dd.options:
            self.cell_type_dd.value = self.cell_type_dd.options[0]

        self._refresh_apply_default_paths()
        self._apply_cfg_defaults()
        self._load_columns()
        self._build_feature_groups()
        self._load_existing_model_if_available()
        self._refresh_backprojection_samples()
        self._rebuild_intrinsic_rename_rows()
        self._rebuild_full_rename_rows()
        self._collapse_all_steps()
        self._refresh_enablement()

    def _on_refresh_clicked(self, _):
        self._set_busy(self.refresh_btn, self.refresh_spinner, busy=True)
        try:
            self._refresh_context()
        finally:
            self._set_busy(self.refresh_btn, self.refresh_spinner, busy=False)

    def _on_cell_type_changed(self, _):
        self._refresh_context()

    def _on_apply_path_changed(self, _):
        self._refresh_enablement()

    def _on_open_backprojection_clicked(self, _):
        self._set_busy(self.btn_open_backprojection, self.backprojection_spinner, busy=True)
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                sample_name = self.backproj_sample_dd.value
                if sample_name is None or len(str(sample_name).strip()) == 0:
                    raise ValueError("Please select a sample name.")
                _winfo(
                    "state-widget",
                    (
                        "Opening full-cluster backprojection for "
                        f"sample '{sample_name}' and cell_type '{self._current_cell_type()}'"
                    ),
                )
                show_behavioral_state_backprojection(
                    sample_name=str(sample_name),
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    state_col="full_behavioral_cluster",
                    auto_create_if_missing=True,
                    refresh_if_stale=True,
                    run=True,
                    verbose=True,
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_open_backprojection, self.backprojection_spinner, busy=False)
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
                if not Path(full_pkl_path).exists():
                    raise FileNotFoundError(f"Full classifier file not found: {full_pkl_path}")

                full_artifact = load_state_classifier_artifact(full_pkl_path)
                if not (isinstance(full_artifact, dict) and ("continuous_feature_cols" in full_artifact)):
                    raise ValueError(
                        "The supplied full classifier path is not a full-classifier artifact. "
                        "Please provide the full classifier .pkl."
                    )
                intrinsic_artifact = None
                if intrinsic_pkl_path != "":
                    if not Path(intrinsic_pkl_path).exists():
                        raise FileNotFoundError(f"Intrinsic classifier file not found: {intrinsic_pkl_path}")
                    intrinsic_artifact = load_state_classifier_artifact(intrinsic_pkl_path)
                    if isinstance(intrinsic_artifact, dict) and ("continuous_feature_cols" in intrinsic_artifact):
                        raise ValueError(
                            "The supplied intrinsic classifier path appears to be a full-classifier artifact. "
                            "Provide an intrinsic classifier .pkl or leave intrinsic empty."
                        )

                ct = self._current_cell_type()
                self.adata_full = apply_state_classifiers_to_full_dataset(
                    output_dir=self.output_dir,
                    cell_type=ct,
                    label_classifier_artifact=intrinsic_artifact,
                    full_label_classifier_artifact=full_artifact,
                    combine_binary_with_continuous=False,
                    verbose=True,
                )

                n_rows = int(self.adata_full.n_obs)
                n_intrinsic = (
                    int(self.adata_full.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
                    if "intrinsic_behavioral_cluster" in self.adata_full.obs.columns
                    else 0
                )
                n_full = (
                    int(self.adata_full.obs["full_behavioral_cluster"].astype(str).nunique())
                    if "full_behavioral_cluster" in self.adata_full.obs.columns
                    else 0
                )
                _winfo(
                    "state-widget",
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



    def _on_rename_full_clicked(self, _):
        self._set_busy(
            self.btn_rename_full,
            self.rename_full_spinner,
            busy=True,
            disable_buttons=[self.btn_combine_full],
        )
        self.out_rename_full.clear_output()
        with self.out_rename_full:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                mapping = {}
                for old_name, txt in self._full_name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)

                result = self._apply_full_rename_mapping(mapping, save_compression="lzf")
                if not bool(result.get("changed", False)):
                    _winfo("state-widget", "No full-mapping changes to apply.")
                else:
                    _winfo("state-widget", f"Renamed full-mapping clusters and saved: {self._model_adata_path()}")
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(
                    self.btn_rename_full,
                    self.rename_full_spinner,
                    busy=False,
                    disable_buttons=[self.btn_combine_full],
                )
                self._refresh_enablement()

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
            for x in self.model_adata.obs.get("full_behavioral_cluster", pd.Series(dtype="object")).astype(str)
        }
        has_changes = any(
            str(normalized_mapping.get(label, label)) != str(label)
            for label in existing_labels
        )
        if not has_changes:
            return {"changed": False, "relabel_s": 0.0, "save_s": 0.0, "rebuild_s": 0.0, "total_s": 0.0}

        t0 = perf_counter()
        t_relabel0 = perf_counter()
        relabel_cluster_ids(
            adata=self.model_adata,
            mapping=normalized_mapping,
            cluster_key="full_behavioral_cluster",
            overwrite_original=True,
            keep_unmapped=True,
        )
        t_relabel = perf_counter() - t_relabel0

        t_save0 = perf_counter()
        self._save_model_adata(compression=save_compression)
        t_save = perf_counter() - t_save0

        t_rebuild0 = perf_counter()
        self._rebuild_full_rename_rows()
        t_rebuild = perf_counter() - t_rebuild0
        t_total = perf_counter() - t0
        return {
            "changed": True,
            "relabel_s": float(t_relabel),
            "save_s": float(t_save),
            "rebuild_s": float(t_rebuild),
            "total_s": float(t_total),
        }

    def _on_combine_full_clicked(self, _):
        self._set_busy(
            self.btn_combine_full,
            self.combine_full_spinner,
            busy=True,
            disable_buttons=[self.btn_rename_full],
        )
        self.out_rename_full.clear_output()
        with self.out_rename_full:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                if len(self._full_name_boxes) == 0:
                    raise ValueError("No full-state rows are available to combine.")

                target = str(self.full_combine_name.value).strip()
                if target == "":
                    raise ValueError("Cannot rename to empty.")

                selected = [name for name, cb in self._full_select_boxes.items() if bool(cb.value)]
                if len(selected) == 0:
                    raise ValueError("Select at least one full state to combine.")

                # Apply combined name to selected rows first.
                for old_name in selected:
                    if old_name in self._full_name_boxes:
                        self._full_name_boxes[old_name].value = target

                mapping = {}
                for old_name, txt in self._full_name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)

                result = self._apply_full_rename_mapping(mapping, save_compression="lzf")
                if not bool(result.get("changed", False)):
                    _winfo("state-widget", "No full-state combine changes to apply.")
                else:
                    _winfo("state-widget", f"Combined {len(selected)} full states into '{target}'.")
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(
                    self.btn_combine_full,
                    self.combine_full_spinner,
                    busy=False,
                    disable_buttons=[self.btn_rename_full],
                )
                self._refresh_enablement()

