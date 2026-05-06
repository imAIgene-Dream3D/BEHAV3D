import traceback
import re
from pathlib import Path

import ipywidgets as widgets
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import scanpy as sc
import yaml

from behav3d.analysis.clustering.state.visualization.backprojection import (
    show_behavioral_state_backprojection,
)
from behav3d.analysis.clustering.track.classification import (
    apply_track_classifier_to_subtracks,
    build_identity_cluster_mapping,
    export_track_cluster_backprojection,
    generate_track_clustering_report_pdfs,
    get_track_classifier_filename,
    get_track_trajectories_filename,
    rename_track_clusters,
    run_state_based_analysis,
    train_track_classifier,
)
from behav3d.analysis.clustering.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
    filter_multicolor_inputs,
)
from behav3d.widgets.utils import spinning_loader


def _winfo(prefix, message):
    print(f"[{prefix}] INFO {message}")


class TrackClassificationPanel:
    """
    Widget for behavioral state trajectory clustering/classification.

    Workflow:
    1) Run behavioral state trajectory clustering
    2) Rename discovered behavioral state trajectory clusters
    3) Train/apply behavioral state trajectory classifier on (sub)trajectories
    4) Backproject renamed behavioral state trajectory clusters
    """

    def __init__(self, metadata_loader, cell_type=None, show_fixed_io=False, show_ngram_weight=False):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(getattr(metadata_loader, "output_dir", "")).expanduser())
        self.model_adata = None
        self._loaded_model_path = None
        self._name_boxes = {}
        self.show_fixed_io = bool(show_fixed_io)
        self.show_ngram_weight = bool(show_ngram_weight)

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
            options=sorted(set(self._cell_types)),
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

        self._build_clustering_section()
        self._build_rename_section()
        self._build_classifier_section()
        self._build_backprojection_section()

        self._build_steps()

        self.ui = widgets.VBox(
            [
                widgets.HTML('<div style="font-size:20px;font-weight:700;">Behavioral State Trajectory Classification</div>'),
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
        self.classifier_path.observe(self._on_classifier_path_changed, names="value")
        self._refresh_context()

    def _build_steps(self):
        step_defs = [
            ("Behavioral state trajectory clustering", self.clustering_section),
            ("Rename behavioral state trajectory clusters", self.rename_section),
            ("Train behavioral state trajectory classifier", self.train_classifier_section),
            ("Apply behavioral state trajectory classifier", self.apply_classifier_section),
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

    def _build_clustering_section(self):
        self.fixed_io_html = widgets.HTML("")

        self.behavioral_trajectory_size = widgets.IntText(
            description="Behavioral trajectory size",
            value=100,
            layout=widgets.Layout(width="300px"),
            style={"description_width": "170px"},
        )
        self.use_fractions = widgets.Checkbox(
            description="state_proportion",
            value=True,
            indent=False,
            layout=widgets.Layout(width="220px"),
        )
        self.use_bouts_mean_length = widgets.Checkbox(
            description="bouts_mean_length",
            value=True,
            indent=False,
            layout=widgets.Layout(width="220px"),
        )
        self.use_bouts_nr = widgets.Checkbox(
            description="bouts_nr",
            value=True,
            indent=False,
            layout=widgets.Layout(width="220px"),
        )
        self.use_bouts_max_length = widgets.Checkbox(
            description="bouts_max_length",
            value=True,
            indent=False,
            layout=widgets.Layout(width="220px"),
        )
        self.use_transitions = widgets.Checkbox(
            description="transition probabilities",
            value=True,
            indent=False,
            layout=widgets.Layout(width="220px"),
        )
        self.use_bigrams = widgets.Checkbox(
            description="bigrams",
            value=True,
            indent=False,
            layout=widgets.Layout(width="180px"),
        )
        self.use_trigrams = widgets.Checkbox(
            description="trigrams",
            value=True,
            indent=False,
            layout=widgets.Layout(width="180px"),
        )
        self.ngram_weight = widgets.Dropdown(
            description="N-gram weight",
            options=["count", "duration"],
            value="count",
            layout=widgets.Layout(width="230px"),
            style={"description_width": "100px"},
        )
        self.use_bigrams.observe(self._sync_track_feature_group_controls, names="value")
        self.use_trigrams.observe(self._sync_track_feature_group_controls, names="value")
        self._sync_track_feature_group_controls()

        self.n_neighbors = widgets.IntText(
            description="Neighbors",
            value=30,
            layout=widgets.Layout(width="210px"),
            style={"description_width": "90px"},
        )
        self.leiden_resolution = widgets.FloatText(
            description="Resolution",
            value=0.2,
            layout=widgets.Layout(width="220px"),
            style={"description_width": "90px"},
        )
        self.leiden_metric = widgets.Dropdown(
            description="Metric",
            options=["euclidean", "cosine"],
            value="euclidean",
            layout=widgets.Layout(width="220px"),
            style={"description_width": "90px"},
        )
        self.umap_min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=0.1,
            layout=widgets.Layout(width="220px"),
            style={"description_width": "90px"},
        )
        self.random_state = widgets.IntText(
            description="Seed",
            value=123,
            layout=widgets.Layout(width="170px"),
            style={"description_width": "70px"},
        )
        self.n_per_cluster = widgets.IntText(
            description="Exemplars/cluster",
            value=10,
            layout=widgets.Layout(width="220px"),
            style={"description_width": "120px"},
        )

        self.btn_cluster = widgets.Button(
            description="Run behavioral state trajectory clustering",
            button_style="success",
            layout=widgets.Layout(width="360px"),
        )
        self.btn_cluster.on_click(self._on_cluster_clicked)
        self.cluster_spinner = widgets.HTML(value=spinning_loader)
        self.cluster_spinner.layout.display = "none"
        self.out_cluster = widgets.Output()

        self.clustering_section = widgets.VBox(
            self._build_clustering_children()
        )
        if not self.show_fixed_io:
            self.fixed_io_html.layout.display = "none"
        if not self.show_ngram_weight:
            self.ngram_weight.layout.display = "none"

    def _build_clustering_children(self):
        children = []
        if self.show_fixed_io:
            children.append(self.fixed_io_html)
        children.extend(
            [
                widgets.HTML("<h4 style='margin:0;'>Filter and cut trajectories</h4>"),
                widgets.HBox([self.behavioral_trajectory_size]),
                widgets.HTML("<h4 style='margin:6px 0 0 0;'>Feature selection</h4>"),
                widgets.HTML("<i>Select which state-trajectory feature families to include in clustering.</i>"),
                widgets.HBox([self.use_fractions]),
                widgets.HBox([self.use_bouts_mean_length, self.use_bouts_nr, self.use_bouts_max_length]),
                widgets.HBox([self.use_transitions]),
            ]
        )
        if self.show_ngram_weight:
            children.append(widgets.HBox([self.use_bigrams, self.use_trigrams, self.ngram_weight]))
        else:
            children.append(widgets.HBox([self.use_bigrams, self.use_trigrams]))
        children.extend(
            [
                widgets.HTML("<h4 style='margin:6px 0 0 0;'>Clustering</h4>"),
                widgets.HBox([self.n_neighbors, self.leiden_resolution, self.leiden_metric, self.umap_min_dist]),
                widgets.HBox([self.btn_cluster, self.cluster_spinner]),
                self.out_cluster,
            ]
        )
        return children

    def _build_rename_section(self):
        self.rename_status = widgets.HTML("<i>Run clustering or load existing behavioral state trajectory model first.</i>")
        self.rename_rows = widgets.VBox([])
        self.btn_rename = widgets.Button(
            description="Rename behavioral state trajectory clusters",
            button_style="warning",
            layout=widgets.Layout(width="360px"),
        )
        self.btn_rename.on_click(self._on_rename_clicked)
        self.rename_spinner = widgets.HTML(value=spinning_loader)
        self.rename_spinner.layout.display = "none"
        self.out_rename = widgets.Output()

        self.rename_section = widgets.VBox(
            [
                self.rename_status,
                self.rename_rows,
                widgets.HBox([self.btn_rename, self.rename_spinner]),
                self.out_rename,
            ]
        )

    def _build_classifier_section(self):
        self.classifier_status = widgets.HTML("<i>Run clustering or load existing behavioral state trajectory model first.</i>")
        self.classifier_path = widgets.Text(
            description="Classifier .pkl",
            value="",
            placeholder="Auto-detected default path if present",
            layout=widgets.Layout(width="760px"),
            style={"description_width": "120px"},
        )
        self.cls_n_estimators = widgets.IntText(
            description="Trees",
            value=300,
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_min_samples_leaf = widgets.IntText(
            description="Min leaf",
            value=2,
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_n_jobs = widgets.IntText(
            description="n_jobs",
            value=-1,
            layout=widgets.Layout(width="180px"),
            style={"description_width": "80px"},
        )
        self.cls_max_depth = widgets.Text(
            description="Max depth",
            value="",
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_min_samples_split = widgets.IntText(
            description="Min split",
            value=2,
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_max_features = widgets.Text(
            description="Max feat",
            value="sqrt",
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_class_weight = widgets.Text(
            description="Class weight",
            value="",
            layout=widgets.Layout(width="220px"),
            style={"description_width": "110px"},
        )
        self.cls_validation_test_size = widgets.FloatText(
            description="Val size",
            value=0.05,
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_validation_random_state = widgets.Text(
            description="Val seed",
            value="",
            layout=widgets.Layout(width="200px"),
            style={"description_width": "100px"},
        )
        self.cls_validation_stratify = widgets.Checkbox(
            description="Stratify validation",
            value=True,
            indent=False,
        )

        self.btn_train_classifier = widgets.Button(
            description="Train behavioral state trajectory classifier",
            button_style="success",
            layout=widgets.Layout(width="330px"),
        )
        self.btn_train_classifier.on_click(self._on_train_classifier_clicked)
        self.train_classifier_spinner = widgets.HTML(value=spinning_loader)
        self.train_classifier_spinner.layout.display = "none"
        self.out_train_classifier = widgets.Output()

        self.btn_apply_classifier = widgets.Button(
            description="Apply to (sub)trajectories",
            button_style="success",
            layout=widgets.Layout(width="230px"),
        )
        self.btn_apply_classifier.on_click(self._on_apply_classifier_clicked)
        self.apply_classifier_spinner = widgets.HTML(value=spinning_loader)
        self.apply_classifier_spinner.layout.display = "none"
        self.out_apply_classifier = widgets.Output()
        self.apply_default_paths_html = widgets.HTML("")

        classifier_params_box = widgets.VBox(
            [
                widgets.HBox([self.cls_n_estimators, self.cls_min_samples_leaf, self.cls_n_jobs]),
                widgets.HBox([self.cls_max_depth, self.cls_min_samples_split, self.cls_max_features]),
                widgets.HBox([self.cls_class_weight, self.cls_validation_test_size, self.cls_validation_random_state]),
                self.cls_validation_stratify,
            ]
        )
        self.classifier_params_accordion = widgets.Accordion(
            children=[classifier_params_box],
            selected_index=None,
        )
        self.classifier_params_accordion.set_title(0, "classifier parameters")

        self.train_classifier_section = widgets.VBox(
            [
                self.classifier_status,
                self.classifier_params_accordion,
                widgets.HBox([self.btn_train_classifier, self.train_classifier_spinner]),
                self.out_train_classifier,
            ]
        )

        self.apply_classifier_section = widgets.VBox(
            [
                self.classifier_path,
                self.apply_default_paths_html,
                widgets.HBox([self.btn_apply_classifier, self.apply_classifier_spinner]),
                self.out_apply_classifier,
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
        self.btn_backproject = widgets.Button(
            description="Backproject selected sample",
            button_style="success",
            layout=widgets.Layout(width="240px"),
        )
        self.btn_backproject.on_click(self._on_backproject_clicked)
        self.backprojection_spinner = widgets.HTML(value=spinning_loader)
        self.backprojection_spinner.layout.display = "none"
        self.out_backprojection = widgets.Output()

        self.backprojection_section = widgets.VBox(
            [
                self.backprojection_status,
                self.backproj_sample_dd,
                widgets.HBox([self.btn_backproject, self.backprojection_spinner]),
                self.out_backprojection,
            ]
        )

    def _sync_track_feature_group_controls(self, _=None):
        self.ngram_weight.disabled = not (bool(self.use_bigrams.value) or bool(self.use_trigrams.value))

    def _detect_cell_types(self):
        md = getattr(self.metadata_loader, "metadata", None)
        cell_types = []
        if md is not None:
            try:
                cell_types.extend(filter_multicolor_inputs(detect_organoid_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)))
                cell_types.extend(detect_merged_cell_types_from_metadata(md))
            except Exception:
                pass

        out_dir = Path(self.output_dir) if self.output_dir else None
        analysis_dir = (out_dir / "analysis") if out_dir is not None else None
        if analysis_dir is not None and analysis_dir.exists():
            for p in analysis_dir.iterdir():
                if p.is_dir():
                    cell_types.append(p.name)
        return sorted({str(x).strip() for x in cell_types if str(x).strip() != ""})

    def _detect_sample_names(self):
        names = []
        md = getattr(self.metadata_loader, "metadata", None)
        if isinstance(md, pd.DataFrame):
            metadata_candidate_cols = [
                "sample_name",
                "sample",
                "sample_id",
                "sampleid",
                "name",
            ]
            for col in metadata_candidate_cols:
                if col in md.columns:
                    names.extend(
                        md[col]
                        .astype("string")
                        .dropna()
                        .str.strip()
                        .tolist()
                    )
                    break
            if len(names) > 0:
                return sorted({str(x).strip() for x in names if str(x).strip() != ""})

        if self.model_adata is not None and "sample_name" in self.model_adata.obs.columns:
            names.extend(
                self.model_adata.obs["sample_name"].astype("string").dropna().str.strip().tolist()
            )
        out_dir = Path(self.output_dir) if self.output_dir else None
        images_dir = (out_dir / "images") if out_dir is not None else None
        if images_dir is not None and images_dir.exists():
            reserved_non_sample_dirs = {
                "SignalUnmixing",
                "PixelClassification",
                "pixel_classifier",
                "pixelclassification",
            }
            for p in images_dir.iterdir():
                if p.is_dir() and str(p.name) not in reserved_non_sample_dirs:
                    names.append(str(p.name))
        return sorted({str(x).strip() for x in names if str(x).strip() != ""})

    def _current_cell_type(self):
        return str(self.cell_type_dd.value).strip()

    @staticmethod
    def _fixed_state_col():
        return "full_behavioral_cluster"

    @staticmethod
    def _fixed_cluster_key():
        return "ClusterID"

    @staticmethod
    def _fixed_time_col():
        return "position_t"

    @staticmethod
    def _fixed_backprojection_output_col():
        return "track_behavioral_cluster"

    @staticmethod
    def _parse_optional_int(text_value):
        value = str(text_value).strip()
        if value == "":
            return None
        return int(value)

    @staticmethod
    def _sanitize_filename_token(value, fallback="plot"):
        token = str(value).strip()
        if token == "":
            token = str(fallback)
        token = re.sub(r"[^A-Za-z0-9._-]+", "_", token)
        token = token.strip("._-")
        return token if token != "" else str(fallback)

    def _track_outdir(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return Path(self.output_dir, "analysis", ct, "behavioral_state_trajectories")

    def _clustering_plot_outdir(self, cell_type=None):
        return self._track_outdir(cell_type=cell_type) / "clustering"

    def _default_classifier_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return self._track_outdir(cell_type=ct) / "classification" / get_track_classifier_filename(ct)

    def _model_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return self._track_outdir(cell_type=ct) / get_track_trajectories_filename(ct)

    def _legacy_model_adata_path(self, cell_type=None):
        return self._track_outdir(cell_type=cell_type) / "adata_state_features_clustered.h5ad"

    def _default_adata_full_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return Path(self.output_dir, "analysis", ct, "behavioral_states", f"BEHAV3D_{ct}_behavioral_states.h5ad")

    def _panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        panel = params.setdefault("track_classification", {})
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
        self.use_fractions.value = bool(cfg.get("use_fractions", self.use_fractions.value))
        self.use_bouts_mean_length.value = bool(
            cfg.get("use_bouts_mean_length", self.use_bouts_mean_length.value)
        )
        self.use_bouts_nr.value = bool(cfg.get("use_bouts_nr", self.use_bouts_nr.value))
        self.use_bouts_max_length.value = bool(
            cfg.get("use_bouts_max_length", self.use_bouts_max_length.value)
        )
        self.use_transitions.value = bool(cfg.get("use_transitions", self.use_transitions.value))
        self.use_bigrams.value = bool(cfg.get("use_bigrams", self.use_bigrams.value))
        self.behavioral_trajectory_size.value = int(
            cfg.get("behavioral_trajectory_size", self.behavioral_trajectory_size.value)
        )
        self.use_trigrams.value = bool(cfg.get("use_trigrams", self.use_trigrams.value))
        self.ngram_weight.value = str(cfg.get("ngram_weight", self.ngram_weight.value))
        self.n_neighbors.value = int(cfg.get("n_neighbors", self.n_neighbors.value))
        self.leiden_resolution.value = float(cfg.get("leiden_resolution", self.leiden_resolution.value))
        self.leiden_metric.value = str(cfg.get("leiden_metric", self.leiden_metric.value))
        self.umap_min_dist.value = float(cfg.get("umap_min_dist", self.umap_min_dist.value))
        self.n_per_cluster.value = int(cfg.get("n_per_cluster", self.n_per_cluster.value))
        self.random_state.value = int(cfg.get("random_state", self.random_state.value))
        self.classifier_path.value = str(cfg.get("classifier_path", self.classifier_path.value))
        self.cls_n_estimators.value = int(
            cfg.get("classifier_n_estimators", self.cls_n_estimators.value)
        )
        self.cls_min_samples_leaf.value = int(
            cfg.get("classifier_min_samples_leaf", self.cls_min_samples_leaf.value)
        )
        self.cls_n_jobs.value = int(cfg.get("classifier_n_jobs", self.cls_n_jobs.value))
        self.cls_max_depth.value = str(cfg.get("classifier_max_depth", self.cls_max_depth.value))
        self.cls_min_samples_split.value = int(
            cfg.get("classifier_min_samples_split", self.cls_min_samples_split.value)
        )
        self.cls_max_features.value = str(
            cfg.get("classifier_max_features", self.cls_max_features.value)
        )
        self.cls_class_weight.value = str(
            cfg.get("classifier_class_weight", self.cls_class_weight.value)
        )
        self.cls_validation_test_size.value = float(
            cfg.get("validation_test_size", self.cls_validation_test_size.value)
        )
        self.cls_validation_random_state.value = str(
            cfg.get("validation_random_state", self.cls_validation_random_state.value)
        )
        self.cls_validation_stratify.value = bool(
            cfg.get("validation_stratify", self.cls_validation_stratify.value)
        )
        self._sync_track_feature_group_controls()

    def _persist_current_settings(self):
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        cfg.update(
            {
                "use_fractions": bool(self.use_fractions.value),
                "use_bouts_mean_length": bool(self.use_bouts_mean_length.value),
                "use_bouts_nr": bool(self.use_bouts_nr.value),
                "use_bouts_max_length": bool(self.use_bouts_max_length.value),
                "use_transitions": bool(self.use_transitions.value),
                "use_bigrams": bool(self.use_bigrams.value),
                "behavioral_trajectory_size": int(self.behavioral_trajectory_size.value),
                "use_trigrams": bool(self.use_trigrams.value),
                "ngram_weight": str(self.ngram_weight.value),
                "n_neighbors": int(self.n_neighbors.value),
                "leiden_resolution": float(self.leiden_resolution.value),
                "leiden_metric": str(self.leiden_metric.value),
                "umap_min_dist": float(self.umap_min_dist.value),
                "n_per_cluster": int(self.n_per_cluster.value),
                "random_state": int(self.random_state.value),
                "classifier_path": str(self.classifier_path.value),
                "classifier_n_estimators": int(self.cls_n_estimators.value),
                "classifier_min_samples_leaf": int(self.cls_min_samples_leaf.value),
                "classifier_n_jobs": int(self.cls_n_jobs.value),
                "classifier_max_depth": str(self.cls_max_depth.value),
                "classifier_min_samples_split": int(self.cls_min_samples_split.value),
                "classifier_max_features": str(self.cls_max_features.value),
                "classifier_class_weight": str(self.cls_class_weight.value),
                "validation_test_size": float(self.cls_validation_test_size.value),
                "validation_random_state": str(self.cls_validation_random_state.value),
                "validation_stratify": bool(self.cls_validation_stratify.value),
            }
        )
        self._save_panel_cfg()

    def _set_busy(self, button, spinner, busy=True):
        state = bool(busy)
        if button is not None and hasattr(button, "disabled"):
            button.disabled = state
        if spinner is not None:
            spinner.layout.display = None if state else "none"

    def _load_existing_model_if_available(self):
        self.model_adata = None
        self._loaded_model_path = None
        candidate_paths = [self._model_adata_path(), self._legacy_model_adata_path()]
        for p in candidate_paths:
            if not p.exists():
                continue
            try:
                self.model_adata = sc.read_h5ad(p)
                self._loaded_model_path = p
                return
            except Exception:
                self.model_adata = None
                self._loaded_model_path = None

    def _save_model_adata(self, compression="gzip"):
        if self.model_adata is None:
            return
        p = self._model_adata_path()
        p.parent.mkdir(parents=True, exist_ok=True)
        self.model_adata.write(p, compression=compression)
        self._loaded_model_path = p

    def _refresh_apply_default_paths(self):
        if self.output_dir is None or str(self.output_dir).strip() == "":
            self.apply_default_paths_html.value = ""
            return
        default_path = self._default_classifier_path()
        default_exists = default_path.exists()
        if str(self.classifier_path.value).strip() == "" and default_exists:
            self.classifier_path.value = str(default_path)
        if default_exists:
            self.apply_default_paths_html.value = (
                "<b style='color:#080;'>Default classifier path detected and prefilled.</b>"
            )
        else:
            self.apply_default_paths_html.value = ""

    def _refresh_classifier_default_path(self):
        self._refresh_apply_default_paths()

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

    def _rebuild_rename_rows(self):
        self._name_boxes = {}
        cluster_col = self._fixed_cluster_key()

        if self.model_adata is None or cluster_col not in self.model_adata.obs.columns:
            self.rename_rows.children = []
            self.rename_status.value = (
                f"<i>Run clustering first. Missing '{cluster_col}' in model adata.</i>"
            )
            self.btn_rename.disabled = True
            return

        mapping = build_identity_cluster_mapping(self.model_adata, cluster_col=cluster_col)
        rows = []
        for old_name in mapping.keys():
            txt = widgets.Text(value=str(old_name), layout=widgets.Layout(width="300px"))
            self._name_boxes[str(old_name)] = txt
            rows.append(
                widgets.HBox(
                    [widgets.Label(str(old_name), layout=widgets.Layout(width="210px")), txt],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                )
            )
        self.rename_rows.children = rows
        self.rename_status.value = (
            f"<b>Clusters in '{cluster_col}':</b> {len(rows)}"
        )
        self.btn_rename.disabled = len(rows) == 0

    def _refresh_enablement(self):
        has_cell_type = self._current_cell_type() != ""
        has_model = self.model_adata is not None
        has_sample = (
            self.backproj_sample_dd.value is not None
            and len(str(self.backproj_sample_dd.value).strip()) > 0
        )
        self.btn_cluster.disabled = not has_cell_type
        self.btn_rename.disabled = not has_model
        self.btn_train_classifier.disabled = not has_model
        has_classifier_path = str(self.classifier_path.value).strip() != "" or self._default_classifier_path().exists()
        self.btn_apply_classifier.disabled = not (has_cell_type and has_classifier_path)
        self.btn_backproject.disabled = not (has_model and has_sample and has_cell_type)
        if has_model:
            loaded_name = self._model_adata_path().name
            if self._loaded_model_path is not None:
                loaded_name = Path(self._loaded_model_path).name
            self.status_html.value = (
                f"<b>Model loaded:</b> {loaded_name} "
                f"({self.model_adata.n_obs} rows)"
            )
            self.classifier_status.value = "<b>Model loaded.</b> Train or apply behavioral state trajectory classifier."
        else:
            self.status_html.value = "<i>Run behavioral state trajectory clustering or load existing model.</i>"
            self.classifier_status.value = "<i>Run clustering or load existing behavioral state trajectory model first.</i>"

    def _refresh_context(self):
        self.output_dir = str(Path(getattr(self.metadata_loader, "output_dir", "")).expanduser())
        self.output_dir_html.value = f"<b>Output dir:</b> {self.output_dir}"
        self.fixed_io_html.value = (
            f"<b>Input h5ad (auto):</b> {self._default_adata_full_path()}<br>"
            f"<b>Output folder (fixed):</b> {self._track_outdir()}<br>"
            "<b>Fixed columns:</b> "
            f"state_col='{self._fixed_state_col()}', "
            f"cluster_key='{self._fixed_cluster_key()}', "
            f"time_col='{self._fixed_time_col()}'"
        )

        current_value = self._current_cell_type()
        refreshed_types = self._detect_cell_types()
        if len(refreshed_types) == 0:
            refreshed_types = [current_value or "tcell"]
        if current_value not in refreshed_types:
            refreshed_types.append(current_value)
        new_options = sorted({str(x).strip() for x in refreshed_types if str(x).strip() != ""})
        if len(new_options) == 0:
            new_options = ["tcell"]
        if list(self.cell_type_dd.options) != list(new_options):
            self.cell_type_dd.options = new_options
        if self.cell_type_dd.value not in self.cell_type_dd.options:
            self.cell_type_dd.value = self.cell_type_dd.options[0]

        self._apply_cfg_defaults()
        self._load_existing_model_if_available()
        self._refresh_apply_default_paths()
        self._refresh_backprojection_samples()
        self._rebuild_rename_rows()
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

    def _on_classifier_path_changed(self, _):
        self._refresh_enablement()

    def _on_cluster_clicked(self, _):
        self._set_busy(self.btn_cluster, self.cluster_spinner, busy=True)
        self.out_cluster.clear_output()
        with self.out_cluster:
            try:
                self._persist_current_settings()
                use_bout_stats = any(
                    [
                        bool(self.use_bouts_mean_length.value),
                        bool(self.use_bouts_nr.value),
                        bool(self.use_bouts_max_length.value),
                    ]
                )
                selected_track_feature_count = sum(
                    [
                        bool(self.use_fractions.value),
                        bool(use_bout_stats),
                        bool(self.use_transitions.value),
                        bool(self.use_bigrams.value),
                        bool(self.use_trigrams.value),
                    ]
                )
                if selected_track_feature_count == 0:
                    raise ValueError(
                        "At least one trajectory-state feature option must be enabled "
                        "(fractions, bout stats, transitions, bigrams, or trigrams)."
                    )
                ct = self._current_cell_type()
                cluster_key = self._fixed_cluster_key()
                self.model_adata = run_state_based_analysis(
                    output_dir=self.output_dir,
                    cell_type=ct,
                    adata_full_path=None,
                    state_col=self._fixed_state_col(),
                    groupby_cols=("sample_name", "TrackID"),
                    time_col=self._fixed_time_col(),
                    behavioral_trajectory_size=int(self.behavioral_trajectory_size.value),
                    use_fractions=bool(self.use_fractions.value),
                    use_bout_stats=bool(use_bout_stats),
                    include_bouts_mean_length=bool(self.use_bouts_mean_length.value),
                    include_bouts_nr=bool(self.use_bouts_nr.value),
                    include_bouts_max_length=bool(self.use_bouts_max_length.value),
                    use_transitions=bool(self.use_transitions.value),
                    use_bigrams=bool(self.use_bigrams.value),
                    use_trigrams=bool(self.use_trigrams.value),
                    ngram_weight=str(self.ngram_weight.value).strip() or "count",
                    do_block_scaling=False,
                    do_l2_normalization=False,
                    drop_highly_correlated=False,
                    corr_threshold=0.95,
                    drop_low_variance=False,
                    low_var_threshold=1e-4,
                    do_pca=True,
                    pca_var_selection=0.95,
                    n_neighbors=int(self.n_neighbors.value),
                    leiden_resolution=float(self.leiden_resolution.value),
                    leiden_metric=str(self.leiden_metric.value),
                    leiden_use_rep="X",
                    cluster_key=cluster_key,
                    umap_min_dist=float(self.umap_min_dist.value),
                    plot_results=True,
                    plot_exemplars=True,
                    n_per_cluster=int(self.n_per_cluster.value),
                    autosave_plots=True,
                    save_outputs=True,
                    output_subdir_name="behavioral_state_trajectories",
                    random_state=int(self.random_state.value),
                )
                self._save_model_adata()
                _winfo("trajectory-widget", f"Saved model adata: {self._model_adata_path()}")
                self._rebuild_rename_rows()
                self._refresh_backprojection_samples()
                self._open_step(1)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_cluster, self.cluster_spinner, busy=False)
                self._refresh_enablement()

    def _on_rename_clicked(self, _):
        self._set_busy(self.btn_rename, self.rename_spinner, busy=True)
        self.out_rename.clear_output()
        with self.out_rename:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                cluster_col = self._fixed_cluster_key()
                mapping = {}
                for old_name, txt in self._name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)
                rename_track_clusters(
                    adata=self.model_adata,
                    mapping=mapping,
                    cluster_col=cluster_col,
                    keep_unmapped=True,
                )
                self._save_model_adata(compression="lzf")
                self._regenerate_clustering_plots(cluster_col=cluster_col)
                self._rebuild_rename_rows()
                _winfo("trajectory-widget", f"Renamed clusters and saved: {self._model_adata_path()}")
                self._open_step(2)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_rename, self.rename_spinner, busy=False)
                self._refresh_enablement()

    def _regenerate_clustering_plots(self, cluster_col):
        if self.model_adata is None:
            return

        plot_outdir = self._clustering_plot_outdir()
        plot_outdir.mkdir(parents=True, exist_ok=True)
        report_paths = generate_track_clustering_report_pdfs(
            adata_tracks=self.model_adata,
            outfolder=plot_outdir,
            cluster_key=cluster_col,
            heatmap_figsize=(25, 20),
            matrixplot_figsize=(20, 40),
            umap_size=1,
            umap_alpha=0.5,
            plot_dpi=300,
            filename_suffix="_renamed",
        )
        _winfo("trajectory-widget", f"Regenerated diagnostics PDF: {report_paths['diagnostics_pdf']}")

        adata_full_path = self._default_adata_full_path()
        if not adata_full_path.exists():
            _winfo("trajectory-widget", f"Skipped exemplar replot (missing full dataset): {adata_full_path}")
            return
        adata_full = sc.read_h5ad(adata_full_path)
        filter_cfg = {}
        if hasattr(self.model_adata, "uns"):
            raw_cfg = self.model_adata.uns.get("track_filtering", {})
            if isinstance(raw_cfg, dict):
                filter_cfg = dict(raw_cfg)
        groupby_cols = filter_cfg.get("groupby_cols", ("sample_name", "TrackID"))
        if not isinstance(groupby_cols, (list, tuple)) or len(groupby_cols) == 0:
            groupby_cols = ("sample_name", "TrackID")
        time_col = str(filter_cfg.get("time_col", self._fixed_time_col()))
        if "behavioral_trajectory_size" not in filter_cfg:
            raise ValueError(
                "Model adata is missing required 'track_filtering.behavioral_trajectory_size'. "
                "Regenerate the behavioral state trajectory model."
            )
        behavioral_trajectory_size = int(filter_cfg["behavioral_trajectory_size"])
        adata_plot = filter_and_truncate_tracks_anndata(
            adata_full,
            groupby_cols=[str(c) for c in groupby_cols],
            time_col=time_col,
            min_length=int(behavioral_trajectory_size),
            max_length=int(behavioral_trajectory_size),
        )
        state_key = self._fixed_state_col()
        fig, _, _ = plot_exemplar_tracks_by_cluster(
            adata_full=adata_plot,
            adata_tracks=self.model_adata,
            n_per_cluster=int(self.n_per_cluster.value),
            state_key=state_key,
            cluster_key=str(cluster_col),
            seed=int(self.random_state.value),
        )
        plot_name = "example_tracks_overview_renamed.pdf"
        exemplar_outdir = (
            self._track_outdir()
            / "clustering"
            / "example_tracks"
            / "per_cluster"
            / "pdf"
        )
        exemplar_outdir.mkdir(parents=True, exist_ok=True)
        plot_path = exemplar_outdir / plot_name
        width, height = fig.get_size_inches()
        if abs(width - height) < 0.05:
            fig.set_size_inches(max(width, height * 1.25), height, forward=True)
        with PdfPages(plot_path) as pdf:
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
        _winfo("trajectory-widget", f"Regenerated exemplar plot: {plot_path}")
        plt.close(fig)

    def _plot_classifier_exemplars(self, cluster_col):
        if self.model_adata is None:
            return
        if cluster_col not in self.model_adata.obs.columns:
            raise ValueError(
                f"Cannot plot classifier exemplars: missing '{cluster_col}' in model adata."
            )

        adata_full_path = self._default_adata_full_path()
        if not adata_full_path.exists():
            raise FileNotFoundError(f"Full dataset h5ad not found: {adata_full_path}")
        adata_full = sc.read_h5ad(adata_full_path)
        filter_cfg = {}
        if hasattr(self.model_adata, "uns"):
            raw_cfg = self.model_adata.uns.get("track_filtering", {})
            if isinstance(raw_cfg, dict):
                filter_cfg = dict(raw_cfg)
        groupby_cols = filter_cfg.get("groupby_cols", ("sample_name", "TrackID"))
        if not isinstance(groupby_cols, (list, tuple)) or len(groupby_cols) == 0:
            groupby_cols = ("sample_name", "TrackID")
        time_col = str(filter_cfg.get("time_col", self._fixed_time_col()))
        if "behavioral_trajectory_size" not in filter_cfg:
            raise ValueError(
                "Model adata is missing required 'track_filtering.behavioral_trajectory_size'. "
                "Regenerate the behavioral state trajectory model."
            )
        behavioral_trajectory_size = int(filter_cfg["behavioral_trajectory_size"])

        adata_plot = filter_and_truncate_tracks_anndata(
            adata_full,
            groupby_cols=[str(c) for c in groupby_cols],
            time_col=time_col,
            min_length=int(behavioral_trajectory_size),
            max_length=int(behavioral_trajectory_size),
        )
        state_key = self._fixed_state_col()
        fig, _, _ = plot_exemplar_tracks_by_cluster(
            adata_full=adata_plot,
            adata_tracks=self.model_adata,
            n_per_cluster=int(self.n_per_cluster.value),
            state_key=state_key,
            cluster_key=str(cluster_col),
            seed=int(self.random_state.value),
        )
        class_outdir = (
            self._track_outdir()
            / "classification"
            / "example_tracks"
            / "per_cluster"
            / "pdf"
        )
        class_outdir.mkdir(parents=True, exist_ok=True)
        plot_name = "example_tracks_overview_classified.pdf"
        plot_path = class_outdir / plot_name
        width, height = fig.get_size_inches()
        if abs(width - height) < 0.05:
            fig.set_size_inches(max(width, height * 1.25), height, forward=True)
        with PdfPages(plot_path) as pdf:
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
        plt.close(fig)
        _winfo("trajectory-widget", f"Saved classifier exemplar plot: {plot_path}")

    def _on_train_classifier_clicked(self, _):
        self._set_busy(self.btn_train_classifier, self.train_classifier_spinner, busy=True)
        self.out_train_classifier.clear_output()
        with self.out_train_classifier:
            try:
                self._persist_current_settings()
                if self.model_adata is None:
                    self._load_existing_model_if_available()
                if self.model_adata is None:
                    raise ValueError(
                        f"No model adata found at: {self._model_adata_path()}. "
                        "Run clustering first."
                    )

                max_depth = self._parse_optional_int(self.cls_max_depth.value)
                class_weight_raw = str(self.cls_class_weight.value).strip()
                class_weight = None if class_weight_raw == "" else class_weight_raw
                val_seed = self._parse_optional_int(self.cls_validation_random_state.value)
                out = train_track_classifier(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    model_adata=self.model_adata,
                    cluster_col=self._fixed_cluster_key(),
                    classifier_n_estimators=int(self.cls_n_estimators.value),
                    classifier_min_samples_leaf=int(self.cls_min_samples_leaf.value),
                    classifier_n_jobs=int(self.cls_n_jobs.value),
                    classifier_max_depth=max_depth,
                    classifier_min_samples_split=int(self.cls_min_samples_split.value),
                    classifier_max_features=str(self.cls_max_features.value).strip(),
                    classifier_class_weight=class_weight,
                    output_subdir_name="behavioral_state_trajectories",
                    save_classifier=True,
                    classifier_path=(str(self.classifier_path.value).strip() or None),
                    random_state=int(self.random_state.value),
                    validation_test_size=float(self.cls_validation_test_size.value),
                    validation_random_state=val_seed,
                    validation_stratify=bool(self.cls_validation_stratify.value),
                    verbose=True,
                )
                self.model_adata = out.get("model_adata", self.model_adata)
                self._save_model_adata()
                classifier_path = out.get("classifier_path", None)
                if classifier_path is not None:
                    self.classifier_path.value = str(classifier_path)
                self._refresh_apply_default_paths()
                qc_metrics_csv = out.get("classifier_metrics_csv", None)
                _winfo("trajectory-widget", "Classifier training finished.")
                _winfo("trajectory-widget", f"Classifier path: {classifier_path}")
                if qc_metrics_csv is not None:
                    _winfo("trajectory-widget", f"Classifier metrics CSV: {qc_metrics_csv}")
                self._open_step(3)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_train_classifier, self.train_classifier_spinner, busy=False)
                self._refresh_enablement()

    def _on_apply_classifier_clicked(self, _):
        self._set_busy(self.btn_apply_classifier, self.apply_classifier_spinner, busy=True)
        self.out_apply_classifier.clear_output()
        with self.out_apply_classifier:
            try:
                self._persist_current_settings()
                classifier_path = str(self.classifier_path.value).strip()
                if classifier_path == "":
                    classifier_path = str(self._default_classifier_path())
                if not Path(classifier_path).exists():
                    raise FileNotFoundError(f"Classifier file not found: {classifier_path}")

                out = apply_track_classifier_to_subtracks(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    classifier_artifact_or_path=classifier_path,
                    adata_full_path=None,
                    state_col=self._fixed_state_col(),
                    output_col=self._fixed_cluster_key(),
                    confidence_col=f"{self._fixed_cluster_key()}_confidence",
                    output_subdir_name="behavioral_state_trajectories",
                    n_per_cluster=int(self.n_per_cluster.value),
                    save_outputs=True,
                    save_as_model=True,
                    random_state=int(self.random_state.value),
                    verbose=True,
                )
                self.model_adata = out.get("adata_tracks", None)
                if self.model_adata is None:
                    raise RuntimeError("Classifier apply returned no adata_tracks.")
                self._save_model_adata()
                self._rebuild_rename_rows()
                self._refresh_backprojection_samples()
                self._plot_classifier_exemplars(cluster_col=self._fixed_cluster_key())
                _winfo("trajectory-widget", "Classifier apply finished.")
                _winfo("trajectory-widget", f"Saved model adata: {self._model_adata_path()}")
                self._open_step(4)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_apply_classifier, self.apply_classifier_spinner, busy=False)
                self._refresh_enablement()

    def _on_backproject_clicked(self, _):
        self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=True)
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                if self.model_adata is None:
                    raise ValueError("No model adata loaded.")
                sample_name = self.backproj_sample_dd.value
                if sample_name is None or len(str(sample_name).strip()) == 0:
                    raise ValueError("Please select a sample.")

                adata_full_path = self._default_adata_full_path()
                if not adata_full_path.exists():
                    raise FileNotFoundError(f"Full dataset h5ad not found: {adata_full_path}")

                adata_full = sc.read_h5ad(adata_full_path)
                cluster_col = self._fixed_cluster_key()
                output_col = self._fixed_backprojection_output_col()
                manifest = export_track_cluster_backprojection(
                    adata_full=adata_full,
                    adata_tracks=self.model_adata,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    cluster_col=cluster_col,
                    output_col=output_col,
                    sample_name=str(sample_name),
                    verbose=True,
                )
                sample_key = str(sample_name).strip()
                state_img_path = manifest.get("output_paths", {}).get(sample_key, None)
                if state_img_path is None:
                    raise RuntimeError(
                        "Backprojection export finished but no state image was written for sample "
                        f"'{sample_key}'. manifest={manifest}"
                    )

                _winfo(
                    "trajectory-widget",
                    f"Opening backprojection for sample '{sample_key}' | state_image={state_img_path}",
                )
                show_behavioral_state_backprojection(
                    sample_name=sample_key,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    state_col=output_col,
                    state_img_path=state_img_path,
                    auto_create_if_missing=False,
                    refresh_if_stale=False,
                    run=True,
                    verbose=True,
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=False)
                self._refresh_enablement()
