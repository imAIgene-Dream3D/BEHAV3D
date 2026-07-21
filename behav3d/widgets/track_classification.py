from pathlib import Path
import traceback

import anndata as ad
import ipywidgets as widgets
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _set_classification_state_colors,
    _set_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    hash_stable_label_color_map,
)
from behav3d.analysis.behavior.track.bouts import (
    apply_track_classifier_to_subtracks,
    rename_track_clusters,
)
from behav3d.analysis.behavior.track.visualization.backprojection import (
    export_track_cluster_backprojection,
    show_track_cluster_backprojection,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    select_exemplar_tracks_by_cluster,
)
from behav3d.analysis.behavior.track.feature_dtw import (
    _create_original_behav3d_adata,
    _feature_dtw_clustered_csv_path,
    _feature_dtw_outdir,
    _feature_dtw_output_csv_paths,
    _feature_dtw_rename_mapping_path,
    _feature_dtw_umap_csv_path,
    _load_feature_dtw_cluster_order,
    _load_feature_dtw_name_mapping,
    _save_feature_dtw_quality_control,
    run_tcell_analysis,
)
from behav3d.analysis.behavior.track.state_dtw import (
    run_categorical_dtaidistance_trajectory_clustering,
    save_dtaidistance_diagnostics,
    save_dtaidistance_exemplar_plots,
    train_dtaidistance_trajectory_classifier,
)
from behav3d.analysis.behavior.track.visualization.plots.reports import (
    save_track_class_proportions_by_sample_plot,
    save_track_condition_comparison_report,
    save_track_contact_group_analysis,
)
from behav3d.analysis.behavior.track.contact_grouping import list_available_contact_columns
from behav3d.analysis.behavior.track.utils import (
    _build_identity_cluster_mapping_from_obs,
    _default_behavioral_states_path,
    _filter_tracks_for_dtaidistance,
    _ordered_unique,
    _resolve_dtaidistance_paths,
    _resolve_optional_int,
    _winfo,
    get_dtaidistance_track_trajectories_filename,
)
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
    filter_multicolor_inputs,
    list_grouping_candidate_columns,
    merge_condition_columns_into_obs,
)
from behav3d.core.utils import rmtree_ignore_missing
from behav3d.widgets.utils import build_plot_box, build_row_move_buttons, spinning_loader


class TrackClassificationPanel:
    """Notebook widget for one-hot categorical track clustering."""

    def __init__(self, metadata_loader, cell_type=None, **kwargs):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(getattr(metadata_loader, "output_dir", "")).expanduser())
        self.model_adata = None

        cell_types = self._detect_cell_types()
        if len(cell_types) == 0:
            cell_types = ["tcell"]
        initial_cell_type = str(cell_type) if cell_type is not None else cell_types[0]
        if initial_cell_type not in cell_types:
            cell_types.append(initial_cell_type)

        self.cell_type_dd = widgets.Dropdown(
            description="Cell type",
            options=list(cell_types),
            value=initial_cell_type,
            layout=widgets.Layout(width="260px"),
            style={"description_width": "90px"},
        )
        self.refresh_btn = widgets.Button(description="Refresh", icon="refresh", layout=widgets.Layout(width="110px"))
        self.refresh_spinner = widgets.HTML(value=spinning_loader)
        self.refresh_spinner.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.output_dir_html = widgets.HTML("")
        self.use_original_behav3d = widgets.Checkbox(
            description="Use original feature-based BEHAV3D DTW clustering",
            value=False,
            indent=False,
            layout=widgets.Layout(width="420px"),
        )
        self.use_original_behav3d.observe(self._on_use_original_behav3d_changed, names="value")
        self.apply_pretrained_classifier = widgets.Checkbox(
            description="Apply pretrained trajectory classifier",
            value=False,
            indent=False,
            layout=widgets.Layout(width="380px"),
        )
        self.apply_pretrained_classifier.observe(self._on_apply_pretrained_classifier_changed, names="value")

        self.state_col_html = widgets.HTML("")
        self.behavioral_trajectory_size = widgets.Text(
            description="Trajectory size",
            value="100",
            placeholder="blank = variable-length DTW",
            layout=widgets.Layout(width="230px"),
            style={"description_width": "110px"},
        )
        self.behavioral_trajectory_size.observe(self._on_trajectory_size_changed, names="value")
        self.n_clusters = widgets.IntText(
            description="Clusters",
            value=6,
            layout=widgets.Layout(width="170px"),
            style={"description_width": "80px"},
        )
        self.n_per_cluster = widgets.IntText(
            description="Exemplars/cluster",
            value=10,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "125px"},
        )
        self.random_state = widgets.IntText(
            description="Seed",
            value=123,
            layout=widgets.Layout(width="150px"),
            style={"description_width": "60px"},
        )
        self.advanced = widgets.Checkbox(description="Advanced", value=False, indent=False)
        self.advanced.observe(self._on_advanced_changed, names="value")

        self.window = widgets.Text(
            description="Window",
            value="",
            placeholder="blank = unconstrained",
            layout=widgets.Layout(width="220px"),
            style={"description_width": "80px"},
        )
        self.max_dist = widgets.Text(
            description="Max dist",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="180px"),
            style={"description_width": "80px"},
        )
        self.max_length_diff = widgets.Text(
            description="Max len diff",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="200px"),
            style={"description_width": "95px"},
        )
        self.penalty = widgets.Text(
            description="Penalty",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="170px"),
            style={"description_width": "70px"},
        )
        self.psi = widgets.Text(
            description="Psi",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="150px"),
            style={"description_width": "45px"},
        )
        self.parallel = widgets.Checkbox(description="Parallel", value=True, indent=False)
        self.linkage = widgets.Dropdown(
            description="Linkage",
            options=["average", "complete", "single"],
            value="average",
            layout=widgets.Layout(width="190px"),
            style={"description_width": "80px"},
        )
        self.missing_policy = widgets.Dropdown(
            description="Missing",
            options=[("Keep as category", "keep"), ("Drop missing timepoints", "drop")],
            value="keep",
            layout=widgets.Layout(width="260px"),
            style={"description_width": "80px"},
        )
        self.save_distance_matrix = widgets.Checkbox(description="Save distance matrix CSV", value=False, indent=False)
        self.trajectory_trim_mode = widgets.Dropdown(
            description="Trim",
            options=[("Last timepoints", "last"), ("First timepoints", "first")],
            value="last",
            layout=widgets.Layout(width="210px"),
            style={"description_width": "60px"},
        )
        self.split_long_tracks = widgets.Checkbox(
            description="Divide long tracks",
            value=False,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )

        self.btn_run = widgets.Button(
            description="Run one-hot dtaidistance clustering",
            button_style="success",
            layout=widgets.Layout(width="330px"),
        )
        self.run_spinner = widgets.HTML(value=spinning_loader)
        self.run_spinner.layout.display = "none"
        self.out_run = widgets.Output()
        self.plot_status_html = widgets.HTML("<i>Run clustering first, then create plots as needed.</i>")
        self.btn_diagnostics = widgets.Button(
            description="Create diagnostics plots",
            button_style="info",
            layout=widgets.Layout(width="230px"),
        )
        self.diagnostics_spinner = widgets.HTML(value=spinning_loader)
        self.diagnostics_spinner.layout.display = "none"
        self.btn_track_proportions = widgets.Button(
            description="Create track proportion plots",
            button_style="info",
            layout=widgets.Layout(width="260px"),
        )
        self.track_proportions_spinner = widgets.HTML(value=spinning_loader)
        self.track_proportions_spinner.layout.display = "none"
        self.btn_exemplars = widgets.Button(
            description="Create exemplar PDFs",
            button_style="info",
            layout=widgets.Layout(width="210px"),
        )
        self.exemplar_spinner = widgets.HTML(value=spinning_loader)
        self.exemplar_spinner.layout.display = "none"
        self.make_overview_statebars = widgets.Checkbox(
            description="Overview statebars",
            value=True,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.make_backprojection_pdf = widgets.Checkbox(
            description="Backprojection PDFs",
            value=False,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.make_backprojection_mp4 = widgets.Checkbox(
            description="Backprojection MP4",
            value=False,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.state_outputs_warning = widgets.HTML("")
        self.out_plots = widgets.Output()
        self.rename_status = widgets.HTML("<i>Run clustering first to rename clusters.</i>")
        self.rename_rows = widgets.VBox([])
        self.btn_refresh_rename = widgets.Button(
            description="Refresh clusters",
            icon="refresh",
            layout=widgets.Layout(width="170px"),
        )
        self.btn_rename = widgets.Button(
            description="Apply cluster names",
            button_style="warning",
            layout=widgets.Layout(width="210px"),
        )
        self.rename_spinner = widgets.HTML(value=spinning_loader)
        self.rename_spinner.layout.display = "none"
        self.out_rename = widgets.Output()
        self._name_boxes = {}
        self.backprojection_status = widgets.HTML("<i>No samples detected yet.</i>")
        self.backproj_sample_dd = widgets.Dropdown(
            description="Sample",
            options=[],
            value=None,
            layout=widgets.Layout(width="360px"),
            style={"description_width": "90px"},
        )
        self.backprojection_workers = widgets.BoundedIntText(
            description="Workers",
            value=4,
            min=1,
            max=32,
            style={"description_width": "90px"},
            layout=widgets.Layout(width="170px"),
        )
        self.btn_backproject = widgets.Button(
            description="Open backprojection",
            button_style="success",
            layout=widgets.Layout(width="220px"),
        )
        self.backprojection_spinner = widgets.HTML(value=spinning_loader)
        self.backprojection_spinner.layout.display = "none"
        self.out_backprojection = widgets.Output()

        self.advanced_row_1 = widgets.HBox(
            [self.window, self.max_dist, self.max_length_diff, self.penalty, self.psi],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.advanced_row_2 = widgets.HBox(
            [
                self.trajectory_trim_mode,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
                self.use_original_behav3d,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.backend_summary_html = widgets.HTML(
            "<b>dtaidistance backend:</b> one-hot vectors with <code>inner_dist='squared euclidean'</code>; "
            "the DTW window is the main speed constraint."
        )

        self.apply_group_x_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in X:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.apply_group_y_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in Y:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.apply_group_cols_select = widgets.SelectMultiple(
            options=[],
            value=[],
            description="",
            rows=2,
            layout=widgets.Layout(width="360px"),
            disabled=True,
        )
        self.apply_group_cols_html = widgets.HTML(
            "<b>Group classification plots by:</b><br>"
            "<span style='color:#555;'>Select metadata columns to group by. "
            "Hold Ctrl/Cmd to select multiple.</span>"
        )
        self.track_proportions_group_html = widgets.HTML(
            "<b>Group by condition:</b><br>"
            "<span style='color:#555;'>\"Group in X\"/\"Group in Y\" pick a single condition each to "
            "arrange the track-class proportion grid. \"Group per page\" pools additional metadata "
            "columns into that grid instead — hold Ctrl/Cmd to select multiple.</span>"
        )

        self.track_comparison_condition_col_dd = widgets.Dropdown(
            options=[], description="Compare condition:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.track_comparison_condition_col_dd.observe(
            self._on_track_comparison_condition_col_changed, names="value",
        )
        self.track_comparison_group_x_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in X:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.track_comparison_group_y_text = widgets.Text(
            value="", description="Group in Y:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.track_comparison_group_cols_select = widgets.SelectMultiple(
            options=[],
            value=[],
            description="",
            rows=2,
            layout=widgets.Layout(width="360px"),
            disabled=True,
        )
        self.btn_track_condition_comparison = widgets.Button(
            description="Create condition comparison plot",
            button_style="info",
            layout=widgets.Layout(width="260px"),
            disabled=True,
        )
        self.track_condition_comparison_spinner = widgets.HTML(value=spinning_loader)
        self.track_condition_comparison_spinner.layout.display = "none"
        self.btn_track_condition_comparison.on_click(self._on_track_condition_comparison_clicked)

        self.contact_col_dd = widgets.Dropdown(
            options=[], description="Contact column:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.contact_min_bout_length = widgets.IntText(
            description="Min. contiguous contact bout (timepoints):",
            value=5,
            style={"description_width": "260px"},
            layout=widgets.Layout(width="360px"),
        )
        self.contact_group_x_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in X:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.contact_group_y_dd = widgets.Dropdown(
            options=["(none)"], value="(none)", description="Group in Y:",
            style={"description_width": "130px"}, layout=widgets.Layout(width="360px"), disabled=True,
        )
        self.contact_group_cols_select = widgets.SelectMultiple(
            options=[],
            value=[],
            description="",
            rows=2,
            layout=widgets.Layout(width="360px"),
            disabled=True,
        )
        self.contact_group_cols_html = widgets.HTML(
            "<b>Also split by condition:</b><br>"
            "<span style='color:#555;'>\"Group in X\"/\"Group in Y\" pick a single condition each to "
            "arrange the grid. \"Group per page\" pools additional metadata columns into it instead — "
            "hold Ctrl/Cmd to select multiple.</span>"
        )
        self.btn_contact_analysis = widgets.Button(
            description="Run contact-vs-no-contact analysis",
            button_style="info",
            layout=widgets.Layout(width="280px"),
            disabled=True,
        )
        self.contact_analysis_spinner = widgets.HTML(value=spinning_loader)
        self.contact_analysis_spinner.layout.display = "none"

        self.run_section = widgets.VBox(
            [
                self.state_col_html,
                widgets.HBox(
                    [
                        self.behavioral_trajectory_size,
                        self.n_clusters,
                        self.random_state,
                        self.advanced,
                    ],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                self.advanced_row_1,
                self.advanced_row_2,
                self.backend_summary_html,
                widgets.HBox([self.btn_run, self.run_spinner]),
                self.out_run,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        diagnostics_box = build_plot_box(
            title="Diagnostics",
            description="Writes cluster counts, medoids, and a distance-matrix heatmap.",
            run_row=widgets.HBox([self.btn_diagnostics, self.diagnostics_spinner], layout=widgets.Layout(gap="8px")),
        )
        track_proportions_box = build_plot_box(
            title="Track-class proportions",
            description=(
                "Writes a stacked bar of track-class proportions per sample, "
                "plus optional condition-grouped grid pages."
            ),
            run_row=widgets.HBox(
                [self.btn_track_proportions, self.track_proportions_spinner],
                layout=widgets.Layout(gap="8px"),
            ),
            settings=[
                self.track_proportions_group_html,
                self.apply_group_x_dd,
                self.apply_group_y_dd,
                widgets.HTML("<span style='color:#555;font-size:11px;'>Group per page:</span>"),
                self.apply_group_cols_select,
            ],
        )
        track_condition_comparison_box = build_plot_box(
            title="Condition comparison plot",
            description=(
                "Per-cluster overall track-class proportion difference between every pairwise "
                "combination of a condition's levels (Welch's t-test), shown as one row per "
                "pairwise comparison."
            ),
            run_row=widgets.HBox(
                [self.btn_track_condition_comparison, self.track_condition_comparison_spinner],
                layout=widgets.Layout(gap="8px"),
            ),
            settings=[
                self.track_comparison_condition_col_dd,
                widgets.HTML(
                    "<b>Group by condition:</b><br>"
                    "<span style='color:#555;'>Each row is one pairwise comparison of \"Compare "
                    "condition\"'s levels (shown in \"Group in Y\"). \"Group in X\" splits the "
                    "comparison into side-by-side columns from another condition. \"Group per page\" "
                    "pools additional metadata columns into that same columns axis instead — hold "
                    "Ctrl/Cmd to select multiple.</span>"
                ),
                self.track_comparison_group_x_dd,
                self.track_comparison_group_y_text,
                widgets.HTML("<span style='color:#555;font-size:11px;'>Group per page:</span>"),
                self.track_comparison_group_cols_select,
            ],
        )
        contact_group_box = build_plot_box(
            title="Contact-based grouping",
            description=(
                "Groups tracks by whether they had a sufficiently long contiguous bout of contact "
                "with another cell type ('contact' vs 'no_contact'), then writes cluster proportions "
                "for each group (and the reverse view), a condition-comparison report between the two "
                "groups, and violin plots of mean contact fraction / max contact-bout length per cluster."
            ),
            run_row=widgets.HBox(
                [self.btn_contact_analysis, self.contact_analysis_spinner],
                layout=widgets.Layout(gap="8px"),
            ),
            settings=[
                self.contact_col_dd,
                self.contact_min_bout_length,
                self.contact_group_cols_html,
                self.contact_group_x_dd,
                self.contact_group_y_dd,
                widgets.HTML("<span style='color:#555;font-size:11px;'>Group per page:</span>"),
                self.contact_group_cols_select,
            ],
        )
        exemplar_box = build_plot_box(
            title="Exemplar tracks",
            description="Writes exemplar overview and per-cluster statebar PDFs.",
            run_row=widgets.HBox([self.btn_exemplars, self.exemplar_spinner], layout=widgets.Layout(gap="8px")),
            settings=[
                widgets.HBox(
                    [
                        self.n_per_cluster,
                        self.make_overview_statebars,
                        self.make_backprojection_pdf,
                        self.make_backprojection_mp4,
                    ],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
            ],
            extra=[self.state_outputs_warning],
        )
        self.plot_section = widgets.VBox(
            [
                self.plot_status_html,
                diagnostics_box,
                track_proportions_box,
                track_condition_comparison_box,
                contact_group_box,
                exemplar_box,
                self.out_plots,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.rename_section = widgets.VBox(
            [
                self.rename_status,
                self.rename_rows,
                widgets.HBox(
                    [self.btn_refresh_rename, self.btn_rename, self.rename_spinner],
                    layout=widgets.Layout(gap="8px"),
                ),
                self.out_rename,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.backprojection_section = widgets.VBox(
            [
                self.backprojection_status,
                widgets.HBox([self.backproj_sample_dd, self.backprojection_workers], layout=widgets.Layout(gap="8px")),
                widgets.HBox([self.btn_backproject, self.backprojection_spinner], layout=widgets.Layout(gap="8px")),
                self.out_backprojection,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        # --- Classify tracks section (dtaidistance mode only) ---
        self.classifier_n_estimators = widgets.IntText(
            description="Number of trees",
            value=100,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.classifier_min_samples_leaf = widgets.IntText(
            description="Min samples leaf",
            value=1,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="200px"),
        )
        self.classifier_test_size = widgets.Text(
            description="Test holdout %",
            value="0.2",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="180px"),
        )
        self.classifier_advanced = widgets.Checkbox(
            description="Advanced",
            value=False,
            indent=False,
        )
        self.classifier_advanced.observe(self._on_classifier_advanced_changed, names="value")
        self.classifier_max_depth = widgets.Text(
            description="Max depth",
            value="",
            placeholder="None",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="170px"),
        )
        self.classifier_min_samples_split = widgets.IntText(
            description="Min samples split",
            value=2,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="210px"),
        )
        self.classifier_max_features = widgets.Dropdown(
            description="Max features",
            options=["sqrt", "log2", "None"],
            value="sqrt",
            style={"description_width": "initial"},
            layout=widgets.Layout(width="190px"),
        )
        self.classifier_n_jobs = widgets.IntText(
            description="n_jobs",
            value=-1,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="140px"),
        )
        self.classifier_advanced_row = widgets.HBox(
            [
                self.classifier_min_samples_leaf,
                self.classifier_min_samples_split,
                self.classifier_max_features,
                self.classifier_n_jobs,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.classifier_advanced_row.layout.display = "none"

        self.btn_train_classifier = widgets.Button(
            description="Train RF classifier",
            button_style="success",
            layout=widgets.Layout(width="200px"),
        )
        self.train_classifier_spinner = widgets.HTML(value=spinning_loader)
        self.train_classifier_spinner.layout.display = "none"
        self.out_train_classifier = widgets.Output()

        self.classifier_artifact_path = widgets.Text(
            description="Classifier path",
            value="",
            placeholder="auto-filled after training",
            style={"description_width": "120px"},
            layout=widgets.Layout(width="560px"),
        )
        self.classifier_apply_states_path = widgets.Text(
            description="States h5ad path",
            value="",
            placeholder="defaults to project behavioral states",
            style={"description_width": "120px"},
            layout=widgets.Layout(width="560px"),
        )
        self.btn_apply_classifier = widgets.Button(
            description="Apply classifier",
            button_style="info",
            layout=widgets.Layout(width="200px"),
        )
        self.apply_classifier_spinner = widgets.HTML(value=spinning_loader)
        self.apply_classifier_spinner.layout.display = "none"
        self.out_apply_classifier = widgets.Output()

        self.classifier_section = widgets.VBox(
            [
                widgets.HTML("<b>Train RF classifier on named clusters</b>"),
                widgets.HBox(
                    [self.classifier_n_estimators, self.classifier_max_depth, self.classifier_test_size, self.classifier_advanced],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                self.classifier_advanced_row,
                widgets.HBox([self.btn_train_classifier, self.train_classifier_spinner], layout=widgets.Layout(gap="8px")),
                self.out_train_classifier,
                widgets.HTML("<b>Apply classifier to new data</b>"),
                self.classifier_artifact_path,
                self.classifier_apply_states_path,
                self.apply_group_cols_html,
                self.apply_group_cols_select,
                widgets.HBox([self.btn_apply_classifier, self.apply_classifier_spinner], layout=widgets.Layout(gap="8px")),
                self.out_apply_classifier,
            ],
            layout=widgets.Layout(gap="8px"),
        )

        self.apply_pretrained_section = widgets.VBox(
            [
                widgets.HTML("<b>Apply pretrained classifier</b>"),
                self.classifier_artifact_path,
                self.classifier_apply_states_path,
                self.apply_group_cols_html,
                self.apply_group_cols_select,
                widgets.HBox([self.btn_apply_classifier, self.apply_classifier_spinner], layout=widgets.Layout(gap="8px")),
                self.out_apply_classifier,
            ],
            layout=widgets.Layout(gap="8px"),
        )

        self.steps = widgets.Accordion(
            children=[self.run_section, self.rename_section, self.plot_section, self.backprojection_section],
            selected_index=0,
        )
        self.steps.set_title(0, "Run clustering")
        self.steps.set_title(1, "Rename clusters")
        self.steps.set_title(2, "Create plots")
        self.steps.set_title(3, "Backprojection")

        self.original_trajectory_size = widgets.IntText(
            description="Trajectory size",
            value=100,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "110px"},
        )
        self.original_n_clusters = widgets.IntText(
            description="Clusters",
            value=5,
            layout=widgets.Layout(width="170px"),
            style={"description_width": "80px"},
        )
        self.original_umap_n_neighbors = widgets.IntText(
            description="UMAP n_neighbors",
            value=15,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "130px"},
        )
        self.original_umap_min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=0.1,
            layout=widgets.Layout(width="220px"),
            style={"description_width": "110px"},
        )
        self.btn_run_original = widgets.Button(
            description="Run original BEHAV3D",
            button_style="success",
            layout=widgets.Layout(width="230px"),
        )
        self.original_spinner = widgets.HTML(value=spinning_loader)
        self.original_spinner.layout.display = "none"
        self.out_original = widgets.Output()
        self.original_mode_switch_row = widgets.HBox(
            [self.use_original_behav3d],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.original_description_html = widgets.HTML(
            "<span style='color:#555;'>Original BEHAV3D clusters tracks from continuous timepoint features. "
            "For each track, selected features are normalized, DTW distances are calculated across feature "
            "trajectories, UMAP projects the distance structure, and K-means assigns the requested clusters.</span>"
        )
        self.original_controls_row = widgets.HBox(
            [
                self.original_trajectory_size,
                self.original_n_clusters,
                self.original_umap_n_neighbors,
                self.original_umap_min_dist,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.original_run_row = widgets.HBox(
            [self.btn_run_original, self.original_spinner],
            layout=widgets.Layout(gap="8px"),
        )
        self.original_run_section = widgets.VBox(
            [
                self.original_description_html,
                self.original_controls_row,
                self.apply_group_cols_html,
                self.apply_group_cols_select,
                self.original_run_row,
                self.out_original,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.mode_body = widgets.VBox([self.steps], layout=widgets.Layout(gap="8px"))

        self.ui = widgets.VBox(
            [
                widgets.HTML('<div style="font-size:20px;font-weight:700;">Behavioral trajectory classification</div>'),
                widgets.HTML(
                    "<span style='color:#555;'>Classifies track super-behaviors using dynamic time warping "
                    "based on state transitions and proportions</span>"
                ),
                widgets.HBox([self.cell_type_dd, self.refresh_btn, self.refresh_spinner], layout=widgets.Layout(gap="8px")),
                self.apply_pretrained_classifier,
                self.status_html,
                self.mode_body,
            ],
            layout=widgets.Layout(gap="8px"),
        )

        self.cell_type_dd.observe(self._on_cell_type_changed, names="value")
        self.refresh_btn.on_click(self._on_refresh_clicked)
        self.btn_run.on_click(self._on_run_clicked)
        self.btn_run_original.on_click(self._on_run_original_clicked)
        self.btn_refresh_rename.on_click(self._on_refresh_rename_clicked)
        self.btn_rename.on_click(self._on_rename_clicked)
        self.btn_diagnostics.on_click(self._on_diagnostics_clicked)
        self.btn_track_proportions.on_click(self._on_track_proportions_clicked)
        self.btn_contact_analysis.on_click(self._on_contact_analysis_clicked)
        self.btn_exemplars.on_click(self._on_exemplars_clicked)
        self.make_overview_statebars.observe(self._on_exemplar_output_changed, names="value")
        self.make_backprojection_pdf.observe(self._on_exemplar_output_changed, names="value")
        self.make_backprojection_mp4.observe(self._on_exemplar_output_changed, names="value")
        self.btn_backproject.on_click(self._on_backproject_clicked)
        self.btn_train_classifier.on_click(self._on_train_classifier_clicked)
        self.btn_apply_classifier.on_click(self._on_apply_classifier_clicked)
        self._refresh_context()
        self._sync_advanced_visibility()
        self._sync_mode()

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

    def _panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.setdefault("behavioral_track_classification", {})
        cell_type = self._current_cell_type()
        if cell_type not in section:
            section[cell_type] = {}
        return section[cell_type]

    def _effective_panel_cfg(self):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return {}
        section = params.get("behavioral_track_classification", {})
        defaults = section.get("defaults", {})
        cell_cfg = section.get(self._current_cell_type(), {})
        return {**defaults, **cell_cfg}

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
        cfg = self._effective_panel_cfg()
        if not isinstance(cfg, dict):
            return
        self.behavioral_trajectory_size.value = str(cfg.get("behavioral_trajectory_size", self.behavioral_trajectory_size.value))
        self.n_clusters.value = int(cfg.get("n_clusters", self.n_clusters.value))
        self.n_per_cluster.value = int(cfg.get("n_per_cluster", self.n_per_cluster.value))
        self.random_state.value = int(cfg.get("random_state", self.random_state.value))
        self.trajectory_trim_mode.value = str(cfg.get("trajectory_trim_mode", self.trajectory_trim_mode.value))
        self.split_long_tracks.value = bool(cfg.get("split_long_tracks", self.split_long_tracks.value))
        self.linkage.value = str(cfg.get("linkage", self.linkage.value))
        self.missing_policy.value = str(cfg.get("missing_policy", self.missing_policy.value))
        self.parallel.value = bool(cfg.get("parallel", self.parallel.value))
        self.save_distance_matrix.value = bool(cfg.get("save_distance_matrix", self.save_distance_matrix.value))
        self.window.value = str(cfg.get("window", self.window.value))
        self.max_dist.value = str(cfg.get("max_dist", self.max_dist.value))
        self.max_length_diff.value = str(cfg.get("max_length_diff", self.max_length_diff.value))
        self.penalty.value = str(cfg.get("penalty", self.penalty.value))
        self.psi.value = str(cfg.get("psi", self.psi.value))
        self.make_overview_statebars.value = bool(cfg.get("make_overview_statebars", self.make_overview_statebars.value))
        self.make_backprojection_pdf.value = bool(cfg.get("make_backprojection_pdf", self.make_backprojection_pdf.value))
        self.make_backprojection_mp4.value = bool(cfg.get("make_backprojection_mp4", self.make_backprojection_mp4.value))
        self.backprojection_workers.value = int(cfg.get("backprojection_workers", self.backprojection_workers.value))
        self.classifier_n_estimators.value = int(cfg.get("classifier_n_estimators", self.classifier_n_estimators.value))
        self.classifier_min_samples_leaf.value = int(cfg.get("classifier_min_samples_leaf", self.classifier_min_samples_leaf.value))
        self.classifier_min_samples_split.value = int(cfg.get("classifier_min_samples_split", self.classifier_min_samples_split.value))
        self.classifier_max_features.value = str(cfg.get("classifier_max_features", self.classifier_max_features.value))
        self.classifier_max_depth.value = str(cfg.get("classifier_max_depth", self.classifier_max_depth.value))
        self.classifier_test_size.value = str(cfg.get("classifier_test_size", self.classifier_test_size.value))
        self.classifier_n_jobs.value = int(cfg.get("classifier_n_jobs", self.classifier_n_jobs.value))
        saved_artifact = str(cfg.get("classifier_artifact_path", "")).strip()
        if saved_artifact and not str(self.classifier_artifact_path.value).strip():
            self.classifier_artifact_path.value = saved_artifact
        self._sync_trim_mode_visibility()

    def _persist_current_settings(self):
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        cfg.update({
            "behavioral_trajectory_size": str(self.behavioral_trajectory_size.value),
            "n_clusters": int(self.n_clusters.value),
            "n_per_cluster": int(self.n_per_cluster.value),
            "random_state": int(self.random_state.value),
            "trajectory_trim_mode": str(self.trajectory_trim_mode.value),
            "split_long_tracks": bool(self.split_long_tracks.value),
            "linkage": str(self.linkage.value),
            "missing_policy": str(self.missing_policy.value),
            "parallel": bool(self.parallel.value),
            "save_distance_matrix": bool(self.save_distance_matrix.value),
            "window": str(self.window.value),
            "max_dist": str(self.max_dist.value),
            "max_length_diff": str(self.max_length_diff.value),
            "penalty": str(self.penalty.value),
            "psi": str(self.psi.value),
            "make_overview_statebars": bool(self.make_overview_statebars.value),
            "make_backprojection_pdf": bool(self.make_backprojection_pdf.value),
            "make_backprojection_mp4": bool(self.make_backprojection_mp4.value),
            "backprojection_workers": int(self.backprojection_workers.value),
            "classifier_n_estimators": int(self.classifier_n_estimators.value),
            "classifier_min_samples_leaf": int(self.classifier_min_samples_leaf.value),
            "classifier_min_samples_split": int(self.classifier_min_samples_split.value),
            "classifier_max_features": str(self.classifier_max_features.value),
            "classifier_max_depth": str(self.classifier_max_depth.value),
            "classifier_test_size": str(self.classifier_test_size.value),
            "classifier_n_jobs": int(self.classifier_n_jobs.value),
            "classifier_artifact_path": str(self.classifier_artifact_path.value),
        })
        self._save_panel_cfg()

    def _detect_cell_types(self):
        md = getattr(self.metadata_loader, "metadata", None)
        cell_types = []
        if md is not None:
            try:
                cell_types.extend(filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_organoid_types_from_metadata(md)))
                cell_types.extend(filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)))
            except Exception:
                pass
        analysis_dir = Path(self.output_dir) / "analysis" if self.output_dir else None
        if analysis_dir is not None and analysis_dir.exists():
            for path in analysis_dir.iterdir():
                if path.is_dir():
                    cell_types.append(path.name)
        return _ordered_unique(cell_types)

    def _current_cell_type(self):
        return str(self.cell_type_dd.value).strip()

    def _model_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return (
            Path(self.output_dir)
            / "analysis"
            / ct
            / "behavorial_trajectories"
            / get_dtaidistance_track_trajectories_filename(ct)
        )

    def _original_track_features_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        feature_outdir = Path(self.output_dir) / "analysis" / ct / "track_features"
        filtered = feature_outdir / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        if filtered.exists():
            return filtered
        return feature_outdir / f"BEHAV3D_{ct}_combined_track_features.csv"

    def _state_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return _default_behavioral_states_path(self.output_dir, ct)

    def _has_behavioral_states(self):
        return Path(self._state_adata_path()).exists()

    def _any_exemplar_outputs_selected(self):
        return any(
            [
                bool(self.make_overview_statebars.value),
                bool(self.make_backprojection_pdf.value),
                bool(self.make_backprojection_mp4.value),
            ]
        )

    def _sync_exemplar_state_controls(self):
        has_states = self._has_behavioral_states()
        for widget in [
            self.make_overview_statebars,
            self.make_backprojection_pdf,
            self.make_backprojection_mp4,
        ]:
            widget.disabled = not has_states

        if not has_states:
            self.state_outputs_warning.value = (
                "<b style='color:#a66;'>No states have been defined.</b> "
                "Run behavioral state classification before creating overview statebars, "
                "backprojection PDFs, or backprojection MP4s."
            )
            self.btn_exemplars.disabled = True
            return

        self.state_outputs_warning.value = ""
        if bool(self.use_original_behav3d.value):
            has_clusters = _feature_dtw_clustered_csv_path(self.output_dir, self._current_cell_type()).exists()
        else:
            has_clusters = self.model_adata is not None or self._model_adata_path().exists()
        self.btn_exemplars.disabled = not (has_clusters and self._any_exemplar_outputs_selected())

    def _on_exemplar_output_changed(self, _):
        self._sync_exemplar_state_controls()

    def _sync_trim_mode_visibility(self):
        has_size = str(self.behavioral_trajectory_size.value).strip() != ""
        self.trajectory_trim_mode.layout.display = None if has_size else "none"
        self.split_long_tracks.layout.display = None if has_size else "none"

    def _on_trajectory_size_changed(self, _):
        self._sync_trim_mode_visibility()

    def _sync_advanced_visibility(self):
        display = None if bool(self.advanced.value) else "none"
        self.advanced_row_1.layout.display = "none"
        self.advanced_row_2.layout.display = display
        self.backend_summary_html.layout.display = "none"
        self._sync_trim_mode_visibility()

    def _on_advanced_changed(self, _):
        self._sync_advanced_visibility()

    def _on_classifier_advanced_changed(self, _):
        display = None if bool(self.classifier_advanced.value) else "none"
        self.classifier_advanced_row.layout.display = display

    def _on_apply_pretrained_classifier_changed(self, _):
        self._sync_mode()

    def _sync_mode(self):
        if bool(self.apply_pretrained_classifier.value):
            self.steps.children = [
                self.apply_pretrained_section,
                self.plot_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Apply classifier")
            self.steps.set_title(1, "Create plots")
            self.steps.set_title(2, "Backprojection")
            self.steps.selected_index = 0
            return
        if bool(self.use_original_behav3d.value):
            self.advanced_row_2.children = [
                self.trajectory_trim_mode,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
            ]
            self.original_run_section.children = [
                self.original_mode_switch_row,
                self.original_description_html,
                self.original_controls_row,
                self.apply_group_cols_html,
                self.apply_group_cols_select,
                self.original_run_row,
                self.out_original,
            ]
            self.steps.children = [
                self.original_run_section,
                self.rename_section,
                self.plot_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Run original BEHAV3D")
            self.steps.set_title(1, "Rename clusters")
            self.steps.set_title(2, "Create plots")
            self.steps.set_title(3, "Backprojection")
            self.status_html.value = (
                "<b>Original BEHAV3D mode:</b> running feature-based DTW on continuous track features."
            )
            self.btn_diagnostics.disabled = not _feature_dtw_umap_csv_path(
                self.output_dir,
                self._current_cell_type(),
            ).exists()
            self.btn_track_proportions.disabled = self.btn_diagnostics.disabled
            self.btn_exemplars.disabled = True
            self.plot_status_html.value = "<i>Run original BEHAV3D clustering first, then create QC or rename clusters.</i>"
            self._rebuild_rename_rows()
            self._sync_exemplar_state_controls()
        else:
            self.advanced_row_2.children = [
                self.trajectory_trim_mode,
                self.split_long_tracks,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
                self.use_original_behav3d,
            ]
            self.original_run_section.children = [
                self.original_description_html,
                self.original_controls_row,
                self.apply_group_cols_html,
                self.apply_group_cols_select,
                self.original_run_row,
                self.out_original,
            ]
            self.steps.children = [
                self.run_section,
                self.rename_section,
                self.classifier_section,
                self.plot_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Run clustering")
            self.steps.set_title(1, "Rename clusters")
            self.steps.set_title(2, "Classify tracks")
            self.steps.set_title(3, "Create plots")
            self.steps.set_title(4, "Backprojection")
            self._refresh_context()

    def _on_use_original_behav3d_changed(self, _):
        self._sync_mode()

    def _set_busy(self, button, spinner, busy=True):
        button.disabled = bool(busy)
        spinner.layout.display = None if bool(busy) else "none"

    def _refresh_context(self):
        self.output_dir = str(Path(getattr(self.metadata_loader, "output_dir", "")).expanduser())
        self.output_dir_html.value = ""
        current = self._current_cell_type()
        cell_types = self._detect_cell_types()
        if len(cell_types) == 0:
            cell_types = [current or "tcell"]
        if current and current not in cell_types:
            cell_types.append(current)
        self.cell_type_dd.options = _ordered_unique(cell_types)
        if self.cell_type_dd.value not in self.cell_type_dd.options:
            self.cell_type_dd.value = self.cell_type_dd.options[0]
        model_path = self._model_adata_path()
        if model_path.exists():
            if not bool(self.use_original_behav3d.value):
                self.status_html.value = "<b style='color:#080;'>Ready for plots:</b> an existing model is available."
            self.plot_status_html.value = "<b>Ready for plots:</b> an existing DTAI model is available."
            self.btn_diagnostics.disabled = False
            self.btn_track_proportions.disabled = False
            self.btn_exemplars.disabled = False
        else:
            if not bool(self.use_original_behav3d.value):
                self.status_html.value = (
                    "<i>No one-hot dtaidistance model detected yet. Press Run to one-hot encode the "
                    f"<code>{FULL_STATE_COL}</code> trajectories and create it.</i>"
                )
            self.plot_status_html.value = "<i>Run clustering first, then create plots as needed.</i>"
            self.btn_diagnostics.disabled = True
            self.btn_track_proportions.disabled = True
            self.btn_exemplars.disabled = True
        self._refresh_backprojection_samples()
        self._rebuild_rename_rows()
        candidate_group_cols = list_grouping_candidate_columns(getattr(self.metadata_loader, "metadata", None))
        if candidate_group_cols != list(self.apply_group_cols_select.options):
            prev_selection = set(self.apply_group_cols_select.value)
            self.apply_group_cols_select.options = candidate_group_cols
            self.apply_group_cols_select.value = [c for c in candidate_group_cols if c in prev_selection]
        self.apply_group_cols_select.rows = max(2, min(len(candidate_group_cols), 6))
        self.apply_group_cols_select.disabled = len(candidate_group_cols) == 0

        if candidate_group_cols != list(self.contact_group_cols_select.options):
            prev_contact_group_selection = set(self.contact_group_cols_select.value)
            self.contact_group_cols_select.options = candidate_group_cols
            self.contact_group_cols_select.value = [
                c for c in candidate_group_cols if c in prev_contact_group_selection
            ]
        self.contact_group_cols_select.rows = max(2, min(len(candidate_group_cols), 6))
        self.contact_group_cols_select.disabled = len(candidate_group_cols) == 0

        if candidate_group_cols != list(self.track_comparison_group_cols_select.options):
            prev_comparison_group_selection = set(self.track_comparison_group_cols_select.value)
            self.track_comparison_group_cols_select.options = candidate_group_cols
            self.track_comparison_group_cols_select.value = [
                c for c in candidate_group_cols if c in prev_comparison_group_selection
            ]
        self.track_comparison_group_cols_select.rows = max(2, min(len(candidate_group_cols), 6))
        self.track_comparison_group_cols_select.disabled = len(candidate_group_cols) == 0

        axis_options = ["(none)"] + candidate_group_cols
        for dd in (
            self.apply_group_x_dd, self.apply_group_y_dd,
            self.contact_group_x_dd, self.contact_group_y_dd,
            self.track_comparison_group_x_dd,
        ):
            if axis_options != list(dd.options):
                prev_axis = dd.value
                dd.options = axis_options
                dd.value = prev_axis if prev_axis in axis_options else "(none)"
            dd.disabled = len(candidate_group_cols) == 0

        if candidate_group_cols != list(self.track_comparison_condition_col_dd.options):
            prev_cond = self.track_comparison_condition_col_dd.value
            self.track_comparison_condition_col_dd.options = candidate_group_cols
            if prev_cond in candidate_group_cols:
                self.track_comparison_condition_col_dd.value = prev_cond
            elif candidate_group_cols:
                self.track_comparison_condition_col_dd.value = candidate_group_cols[0]
        self.track_comparison_condition_col_dd.disabled = len(candidate_group_cols) == 0
        self._sync_track_comparison_group_y_text()
        self.btn_track_condition_comparison.disabled = len(candidate_group_cols) == 0

        contact_columns = []
        features_path = self._original_track_features_path()
        if features_path.exists():
            try:
                contact_columns = list_available_contact_columns(pd.read_csv(features_path, nrows=0))
            except Exception:
                contact_columns = []
        if contact_columns != list(self.contact_col_dd.options):
            prev_contact_col = self.contact_col_dd.value
            self.contact_col_dd.options = contact_columns
            if prev_contact_col in contact_columns:
                self.contact_col_dd.value = prev_contact_col
            elif contact_columns:
                self.contact_col_dd.value = contact_columns[0]
        self.contact_col_dd.disabled = len(contact_columns) == 0
        self.btn_contact_analysis.disabled = len(contact_columns) == 0 or not model_path.exists()
        if bool(self.use_original_behav3d.value):
            self.btn_diagnostics.disabled = not _feature_dtw_umap_csv_path(
                self.output_dir,
                self._current_cell_type(),
            ).exists()
            self.btn_track_proportions.disabled = self.btn_diagnostics.disabled
            self._sync_exemplar_state_controls()
        else:
            self._sync_exemplar_state_controls()
        self._apply_cfg_defaults()
        # Pre-fill classifier paths with defaults if the fields are empty
        if not str(self.classifier_apply_states_path.value).strip():
            self.classifier_apply_states_path.value = str(self._state_adata_path())
        if not str(self.classifier_artifact_path.value).strip():
            clf_path = self._classifier_artifact_path()
            if clf_path.exists():
                self.classifier_artifact_path.value = str(clf_path)

    def _on_refresh_clicked(self, _):
        self._set_busy(self.refresh_btn, self.refresh_spinner, busy=True)
        try:
            self._refresh_context()
        finally:
            self._set_busy(self.refresh_btn, self.refresh_spinner, busy=False)

    def _on_cell_type_changed(self, _):
        self._refresh_context()

    def _load_model_adata(self):
        if self.model_adata is not None:
            return self.model_adata
        model_path = self._model_adata_path()
        if not model_path.exists():
            raise FileNotFoundError(f"No one-hot dtaidistance model found at: {model_path}")
        self.model_adata = sc.read_h5ad(model_path)
        return self.model_adata

    def _detect_backprojection_sample_names(self):
        names = []
        adata_tracks = self.model_adata
        model_path = self._model_adata_path()
        if adata_tracks is None and model_path.exists():
            try:
                adata_tracks = sc.read_h5ad(model_path)
            except Exception:
                adata_tracks = None
        if adata_tracks is not None and "sample_name" in adata_tracks.obs.columns:
            names.extend(
                adata_tracks.obs["sample_name"].astype("string").dropna().str.strip().tolist()
            )
        return sorted({str(name).strip() for name in names if str(name).strip() != ""})

    def _refresh_backprojection_samples(self):
        sample_names = self._detect_backprojection_sample_names()
        self.backproj_sample_dd.options = sample_names
        if len(sample_names) == 0:
            self.backproj_sample_dd.value = None
            self.backprojection_status.value = "<i>No samples detected for backprojection.</i>"
            self.btn_backproject.disabled = True
            return
        if self.backproj_sample_dd.value not in sample_names:
            self.backproj_sample_dd.value = sample_names[0]
        self.backprojection_status.value = f"<b>Available samples:</b> {len(sample_names)}"
        self.btn_backproject.disabled = False

    def _rebuild_rename_rows(self):
        self._name_boxes = {}
        self._track_row_widgets = {}
        self._track_row_order = []
        if bool(self.use_original_behav3d.value):
            clustered_path = _feature_dtw_clustered_csv_path(self.output_dir, self._current_cell_type())
            if not clustered_path.exists():
                self.rename_rows.children = []
                self.rename_status.value = "<i>Run original BEHAV3D clustering first to rename clusters.</i>"
                self.btn_rename.disabled = True
                return
            try:
                df = pd.read_csv(clustered_path, usecols=["ClusterID"], low_memory=False)
            except Exception as exc:
                self.rename_rows.children = []
                self.rename_status.value = f"<i>Could not load original BEHAV3D clusters: {exc}</i>"
                self.btn_rename.disabled = True
                return
            mapping = _build_identity_cluster_mapping_from_obs(df, cluster_col="ClusterID")
            saved_mapping = _load_feature_dtw_name_mapping(self.output_dir, self._current_cell_type())
            for key in list(mapping.keys()):
                mapping[key] = saved_mapping.get(key, key)
            saved_order = _load_feature_dtw_cluster_order(self.output_dir, self._current_cell_type())
            label = "Original BEHAV3D clusters"
        else:
            try:
                adata_tracks = self.model_adata
                if adata_tracks is None and self._model_adata_path().exists():
                    adata_tracks = sc.read_h5ad(self._model_adata_path())
            except Exception:
                adata_tracks = None
            if adata_tracks is None or "ClusterID" not in adata_tracks.obs.columns:
                self.rename_rows.children = []
                self.rename_status.value = "<i>Run DTAI clustering first to rename clusters.</i>"
                self.btn_rename.disabled = True
                return
            mapping = _build_identity_cluster_mapping_from_obs(adata_tracks.obs, cluster_col="ClusterID")
            saved_order = _get_classification_state_order(adata_tracks, "ClusterID")
            label = "DTAI clusters"

        # `saved_order` lists display names (post-rename); `mapping` maps raw keys -> display
        # names, and for the DTAI branch keys already equal display names. Reorder by display
        # name, then translate back to the raw keys the rows/apply-mapping logic is keyed on.
        raw_keys = list(mapping.keys())
        name_to_key = {}
        for k in raw_keys:
            name_to_key.setdefault(str(mapping[k]), k)
        ordered_names = _apply_state_order(list(name_to_key.keys()), saved_order)
        self._track_row_order = [name_to_key[n] for n in ordered_names]
        for old_name in self._track_row_order:
            old_name = str(old_name)
            current_name = mapping[old_name]
            txt = widgets.Text(value=str(current_name), layout=widgets.Layout(width="300px"))
            self._name_boxes[old_name] = txt
            move_btns = build_row_move_buttons(
                on_move_up=lambda n=old_name: self._move_track_row(n, -1),
                on_move_down=lambda n=old_name: self._move_track_row(n, 1),
            )
            self._track_row_widgets[old_name] = widgets.HBox(
                [move_btns, widgets.Label(old_name, layout=widgets.Layout(width="110px")), txt],
                layout=widgets.Layout(align_items="center", gap="8px"),
            )
        self._refresh_track_row_children()
        self.rename_status.value = f"<b>{label}:</b> {len(self._track_row_order)}"
        self.btn_rename.disabled = len(self._track_row_order) == 0

    def _refresh_track_row_children(self):
        self.rename_rows.children = [
            self._track_row_widgets[n] for n in self._track_row_order if n in self._track_row_widgets
        ]

    def _move_track_row(self, old_name, delta):
        order = self._track_row_order
        old_name = str(old_name)
        if old_name not in order:
            return
        idx = order.index(old_name)
        new_idx = idx + delta
        if not (0 <= new_idx < len(order)):
            return
        order[idx], order[new_idx] = order[new_idx], order[idx]
        self._refresh_track_row_children()

    def _on_refresh_rename_clicked(self, _):
        self._rebuild_rename_rows()

    def _on_rename_clicked(self, _):
        self._set_busy(self.btn_rename, self.rename_spinner, busy=True)
        self.out_rename.clear_output()
        with self.out_rename:
            try:
                mapping = {}
                for old_name, txt in self._name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)
                if bool(self.use_original_behav3d.value):
                    self._rename_original_behav3d_clusters(mapping)
                else:
                    self._rename_dtaidistance_clusters(mapping)
                self._rebuild_rename_rows()
                self._persist_current_settings()
                self.steps.selected_index = 2
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_rename, self.rename_spinner, busy=False)

    def _rename_dtaidistance_clusters(self, mapping):
        adata_tracks = self._load_model_adata()
        old_colors = _get_classification_state_colors(adata_tracks, "ClusterID")
        new_colors = {}
        for old_label, new_label in mapping.items():
            new_colors.setdefault(str(new_label), old_colors.get(str(old_label)))
        new_colors = _normalize_label_color_map(list(new_colors.keys()), colors=new_colors)
        rename_track_clusters(
            adata=adata_tracks,
            mapping=mapping,
            cluster_col="ClusterID",
            keep_unmapped=True,
        )
        _set_classification_state_colors(adata_tracks, "ClusterID", new_colors)
        new_order = list(dict.fromkeys(
            mapping.get(n, n) for n in getattr(self, "_track_row_order", [])
        ))
        _set_classification_state_order(adata_tracks, "ClusterID", new_order)
        adata_tracks.write(self._model_adata_path(), compression="lzf")
        qc_dir = _resolve_dtaidistance_paths(self.output_dir, self._current_cell_type())[
            "quality_control_outfolder"
        ] / "after_renaming"
        plot_paths = save_dtaidistance_diagnostics(
            adata_tracks,
            output_dir=self.output_dir,
            cell_type=self._current_cell_type(),
            cluster_key="ClusterID",
            outfolder=qc_dir,
            random_state=int(self.random_state.value),
            save_outputs=True,
            verbose=True,
        )
        self.plot_status_html.value = "<b>Renamed QC ready:</b> diagnostics were written after renaming."
        _winfo("trajectory-dtai-widget", f"Renamed clusters and wrote QC: {plot_paths.get('diagnostics_pdf')}")

    def _rename_original_behav3d_clusters(self, mapping):
        ct = self._current_cell_type()
        mapping_path = _feature_dtw_rename_mapping_path(self.output_dir, ct)
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        existing_colors = {}
        if mapping_path.exists():
            try:
                existing_colors = (yaml.safe_load(mapping_path.read_text()) or {}).get("cluster_colors", {}) or {}
            except Exception:
                existing_colors = {}
        new_names = sorted(set(str(v) for v in mapping.values()))
        cluster_colors = hash_stable_label_color_map(new_names, colors=existing_colors)
        new_order = list(dict.fromkeys(
            mapping.get(n, n) for n in getattr(self, "_track_row_order", [])
        ))
        with mapping_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(
                {
                    "cell_type": ct,
                    "cluster_id_column": "ClusterID",
                    "cluster_name_column": "ClusterName",
                    "cluster_names": dict(mapping),
                    "cluster_colors": cluster_colors,
                    "cluster_order": new_order,
                },
                f,
                sort_keys=False,
            )

        updated = []
        for csv_path in _feature_dtw_output_csv_paths(self.output_dir, ct):
            if not csv_path.exists():
                continue
            df = pd.read_csv(csv_path, low_memory=False)
            if "ClusterID" not in df.columns:
                continue
            df["ClusterName"] = df["ClusterID"].astype(str).map(mapping).fillna(df["ClusterID"].astype(str))
            df.to_csv(csv_path, index=False)
            updated.append(csv_path)

        plot_paths = _save_feature_dtw_quality_control(
            self.output_dir,
            ct,
            cluster_percentage_group_by=list(self.apply_group_cols_select.value) or None,
            proportions_outfolder=_resolve_dtaidistance_paths(self.output_dir, ct)["behavior_proportions_outfolder"],
        )
        self.plot_status_html.value = "<b>Renamed QC ready:</b> original BEHAV3D diagnostics were written after renaming."
        _winfo("trajectory-dtai-widget", f"Saved original BEHAV3D cluster names: {mapping_path}")
        for path in updated:
            _winfo("trajectory-dtai-widget", f"Updated renamed CSV: {path}")
        for path in plot_paths.values():
            _winfo("trajectory-dtai-widget", f"Created renamed QC: {path}")

    def _load_original_exemplar_data(self):
        ct = self._current_cell_type()
        clustered_path = _feature_dtw_clustered_csv_path(self.output_dir, ct)
        if not clustered_path.exists():
            raise FileNotFoundError(f"Original BEHAV3D clustered CSV not found: {clustered_path}")
        state_path = self._state_adata_path(ct)
        if not Path(state_path).exists():
            raise FileNotFoundError(f"Behavioral-state h5ad not found: {state_path}")

        adata_full = sc.read_h5ad(state_path)
        adata_filt = filter_and_truncate_tracks_anndata(
            adata_full,
            groupby_cols=["sample_name", "TrackID"],
            time_col="position_t",
            min_length=int(self.original_trajectory_size.value),
            max_length=int(self.original_trajectory_size.value),
        )
        adata_filt = adata_filt.copy()
        adata_filt.obs["sample_name"] = adata_filt.obs["sample_name"].astype(str)
        adata_filt.obs["TrackID"] = adata_filt.obs["TrackID"].astype(str)
        clusters = pd.read_csv(clustered_path, low_memory=False)
        required = ["sample_name", "TrackID", "ClusterID"]
        missing = [col for col in required if col not in clusters.columns]
        if missing:
            raise ValueError(f"Original BEHAV3D clustered CSV missing required columns: {missing}")
        clusters = (
            clusters[required]
            .dropna(subset=["sample_name", "TrackID", "ClusterID"])
            .assign(
                sample_name=lambda df: df["sample_name"].astype(str),
                TrackID=lambda df: df["TrackID"].astype(str),
                ClusterID=lambda df: df["ClusterID"].astype(str),
            )
            .drop_duplicates(subset=["sample_name", "TrackID"])
        )
        track_obs = (
            adata_filt.obs[["sample_name", "TrackID", "position_t"]]
            .copy()
            .assign(
                sample_name=lambda df: df["sample_name"].astype(str),
                TrackID=lambda df: df["TrackID"].astype(str),
            )
            .groupby(["sample_name", "TrackID"], observed=True, as_index=False)
            .agg(position_t_min=("position_t", "min"), position_t_max=("position_t", "max"))
            .merge(clusters, on=["sample_name", "TrackID"], how="inner")
        )
        if len(track_obs) == 0:
            raise ValueError("No original BEHAV3D clustered tracks matched the behavioral-state h5ad.")
        adata_tracks = ad.AnnData(
            X=np.zeros((len(track_obs), 1), dtype=float),
            obs=track_obs,
            var=pd.DataFrame(index=["feature_dtw_cluster"]),
        )
        return adata_filt, adata_tracks

    def _save_original_exemplar_plots(self):
        adata_filt, adata_tracks = self._load_original_exemplar_data()
        ct = self._current_cell_type()
        outdir = _feature_dtw_outdir(self.output_dir, ct) / "example_tracks"
        outdir.mkdir(parents=True, exist_ok=True)
        n_per_cluster = int(self.n_per_cluster.value)
        num_example_ranks = 5
        chosen_exemplars, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=n_per_cluster,
            sample_key="sample_name",
            track_key="TrackID",
            cluster_key="ClusterID",
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            seed=int(self.random_state.value),
        )
        selection_csv = outdir / "example_track_selection_cluster_ClusterID_state_behavioral_state_original_behav3d.csv"
        chosen_exemplars.to_csv(selection_csv, index=False)

        plot_paths = {"exemplar_selection_csv": str(selection_csv)}
        if bool(self.make_overview_statebars.value):
            fig_exemplar, _, _ = plot_exemplar_tracks_by_cluster(
                adata_filt,
                adata_tracks,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                seed=int(self.random_state.value),
            )
            overview_pdf = outdir / "example_tracks_overview.pdf"
            with PdfPages(overview_pdf) as pdf:
                pdf.savefig(fig_exemplar, bbox_inches="tight", dpi=300)
            plt.close(fig_exemplar)
            statebar_out = save_exemplar_statebar_track_pdf_per_cluster(
                adata_full=adata_filt,
                out_dir=outdir,
                chosen_df=chosen_exemplars,
                adata_tracks=None,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                rows_per_page=6,
                plot_dpi=300,
                seed=int(self.random_state.value),
                cmap_name="tab20",
                layout_mode="both",
                num_example_ranks=num_example_ranks,
            )
            plot_paths["exemplar_tracks_overview_pdf"] = str(overview_pdf)
            plot_paths["exemplar_statebar_track_pdf_by_cluster"] = dict(
                statebar_out.get("pdf_paths_by_cluster", {})
            )

        if bool(self.make_backprojection_pdf.value):
            pdf_out = save_exemplar_statebar_backprojection_pdf(
                adata_full=adata_filt,
                output_dir=self.output_dir,
                cell_type=ct,
                out_dir=outdir / "backprojection",
                chosen_df=chosen_exemplars,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                plot_dpi=220,
                seed=int(self.random_state.value),
                examples_per_cluster=n_per_cluster,
                num_example_ranks=num_example_ranks,
                verbose=True,
            )
            plot_paths["exemplar_backprojection_pdf_by_cluster"] = dict(pdf_out.get("pdf_paths_by_cluster", {}))

        if bool(self.make_backprojection_mp4.value):
            mp4_out = save_exemplar_statebar_backprojection_video_per_cluster(
                adata_full=adata_filt,
                output_dir=self.output_dir,
                cell_type=ct,
                out_dir=outdir / "backprojection",
                chosen_df=chosen_exemplars,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                dpi=180,
                seed=int(self.random_state.value),
                examples_per_cluster=n_per_cluster,
                num_example_ranks=num_example_ranks,
                verbose=True,
            )
            plot_paths["exemplar_backprojection_mp4_by_cluster"] = dict(mp4_out.get("video_paths_by_cluster", {}))

        return plot_paths

    def _original_behav3d_features(self, csv_path):
        df_head = pd.read_csv(csv_path, nrows=5)
        preferred = [
            "mean_square_displacement",
            "speed",
            "mean_dead_dye",
            f"{self._current_cell_type()}_contact",
            "organoid_contact",
        ]
        return [col for col in preferred if col in df_head.columns]

    def _on_run_clicked(self, _):
        self._set_busy(self.btn_run, self.run_spinner, busy=True)
        self.out_run.clear_output()
        with self.out_run:
            try:
                import shutil
                _traj_dir = Path(self.output_dir) / "analysis" / self._current_cell_type() / "behavorial_trajectories"
                if _traj_dir.exists():
                    rmtree_ignore_missing(_traj_dir)
                trajectory_size = _resolve_optional_int(self.behavioral_trajectory_size.value)
                _winfo(
                    "trajectory-dtai-widget",
                    f"Running one-hot dtaidistance on fixed state column='{FULL_STATE_COL}'",
                )
                self.model_adata = run_categorical_dtaidistance_trajectory_clustering(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    behavioral_trajectory_size=trajectory_size,
                    min_track_length=trajectory_size,
                    trajectory_trim_mode=str(self.trajectory_trim_mode.value),
                    split_long_tracks=bool(self.split_long_tracks.value),
                    max_tracks=None,
                    n_clusters=int(self.n_clusters.value),
                    window=None,
                    max_dist=None,
                    max_length_diff=None,
                    penalty=None,
                    psi=None,
                    parallel=bool(self.parallel.value),
                    linkage=str(self.linkage.value),
                    missing_policy="keep",
                    save_distance_matrix=bool(self.save_distance_matrix.value),
                    plot_results=True,
                    plot_exemplars=False,
                    n_per_cluster=int(self.n_per_cluster.value),
                    random_state=int(self.random_state.value),
                    verbose=True,
                )
                self.status_html.value = (
                    f"<b style='color:#080;'>Saved one-hot dtaidistance model:</b> {self._model_adata_path().name} "
                    f"({self.model_adata.n_obs} tracks)"
                )
                umap_error = (self.model_adata.uns.get("visualization", {}) or {}).get("umap_error")
                if umap_error:
                    self.plot_status_html.value = (
                        f"<b style='color:#a60;'>⚠ UMAP was skipped in diagnostics:</b> {umap_error}"
                    )
                else:
                    self.plot_status_html.value = "<b>Ready for plots:</b> clustering finished."
                self.btn_diagnostics.disabled = False
                self.btn_track_proportions.disabled = False
                self.btn_exemplars.disabled = False
                self._rebuild_rename_rows()
                self._refresh_backprojection_samples()
                self._sync_exemplar_state_controls()
                self._persist_current_settings()
                self.steps.selected_index = 1
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_run, self.run_spinner, busy=False)

    def _on_run_original_clicked(self, _):
        self._set_busy(self.btn_run_original, self.original_spinner, busy=True)
        self.out_original.clear_output()
        with self.out_original:
            try:
                ct = self._current_cell_type()
                import shutil
                _traj_dir = Path(self.output_dir) / "analysis" / ct / "behavorial_trajectories"
                if _traj_dir.exists():
                    rmtree_ignore_missing(_traj_dir)
                csv_path = self._original_track_features_path(ct)
                if not csv_path.exists():
                    raise FileNotFoundError(f"Original BEHAV3D track-features CSV not found: {csv_path}")
                features = self._original_behav3d_features(csv_path)
                if len(features) == 0:
                    raise ValueError(
                        "None of the default original BEHAV3D features were found in the track-features CSV."
                    )
                normalize = [
                    col
                    for col in ["mean_square_displacement", "speed", "mean_dead_dye"]
                    if col in features
                ]
                _winfo(
                    "trajectory-dtai-widget",
                    f"Running original BEHAV3D feature DTW with features={features}",
                )
                run_tcell_analysis(
                    cell_type=ct,
                    output_dir=self.output_dir,
                    df_tracks_path=str(csv_path),
                    columns_to_use=features,
                    columns_to_normalize=normalize,
                    umap_minimal_distance=float(self.original_umap_min_dist.value),
                    umap_n_neighbors=int(self.original_umap_n_neighbors.value),
                    nr_of_clusters=int(self.original_n_clusters.value),
                    plot_results=False,
                    seed=int(self.random_state.value),
                    output_subdir_name="behavorial_trajectories/original_behav3d/raw",
                    feature_scaling_preset="original_behav3d",
                    min_track_length=int(self.original_trajectory_size.value),
                    max_track_length=int(self.original_trajectory_size.value),
                    cluster_percentage_group_by=list(self.apply_group_cols_select.value) or None,
                )
                _create_original_behav3d_adata(self.output_dir, ct)
                _winfo("trajectory-dtai-widget", "Original BEHAV3D feature-DTW analysis complete.")
                self.btn_diagnostics.disabled = False
                self.btn_track_proportions.disabled = False
                self.plot_status_html.value = (
                    "<b>Ready for renaming:</b> raw clusters written. Click 'Create diagnostics plots' "
                    "to generate the original_behav3d QC (with any renamed cluster names)."
                )
                self._rebuild_rename_rows()
                self._sync_exemplar_state_controls()
                self.steps.selected_index = 1
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_run_original, self.original_spinner, busy=False)

    def _on_diagnostics_clicked(self, _):
        self._set_busy(self.btn_diagnostics, self.diagnostics_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if bool(self.use_original_behav3d.value):
                    plot_paths = _save_feature_dtw_quality_control(
                        self.output_dir,
                        self._current_cell_type(),
                        cluster_percentage_group_by=list(self.apply_group_cols_select.value) or None,
                        proportions_outfolder=_resolve_dtaidistance_paths(
                            self.output_dir, self._current_cell_type()
                        )["behavior_proportions_outfolder"],
                    )
                    self.plot_status_html.value = "<b>Diagnostics ready:</b> original BEHAV3D QC was written."
                    for path in plot_paths.values():
                        _winfo("trajectory-dtai-widget", f"Created original BEHAV3D QC: {path}")
                else:
                    adata_tracks = self._load_model_adata()
                    plot_paths = save_dtaidistance_diagnostics(
                        adata_tracks,
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        random_state=int(self.random_state.value),
                        save_outputs=True,
                        verbose=True,
                    )
                    umap_error = plot_paths.get("umap_error")
                    if umap_error:
                        self.plot_status_html.value = (
                            f"<b style='color:#a60;'>⚠ UMAP was skipped in diagnostics:</b> {umap_error}"
                        )
                    else:
                        self.plot_status_html.value = "<b>Diagnostics ready:</b> cluster diagnostics were written."
                    _winfo(
                        "trajectory-dtai-widget",
                        f"Created diagnostics plots: {plot_paths.get('diagnostics_pdf')}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_diagnostics, self.diagnostics_spinner, busy=False)

    def _on_track_proportions_clicked(self, _):
        self._set_busy(self.btn_track_proportions, self.track_proportions_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                group_cols = list(self.apply_group_cols_select.value) or None
                group_x = self.apply_group_x_dd.value
                group_x = None if group_x in (None, "(none)") else group_x
                group_y = self.apply_group_y_dd.value
                group_y = None if group_y in (None, "(none)") else group_y
                if bool(self.use_original_behav3d.value):
                    plot_paths = _save_feature_dtw_quality_control(
                        self.output_dir,
                        self._current_cell_type(),
                        cluster_percentage_group_by=group_cols,
                        proportions_outfolder=_resolve_dtaidistance_paths(
                            self.output_dir, self._current_cell_type()
                        )["behavior_proportions_outfolder"],
                    )
                    self.plot_status_html.value = (
                        "<b>Track proportions ready:</b> original BEHAV3D per-sample proportions were written."
                    )
                    for path in plot_paths.values():
                        _winfo("trajectory-dtai-widget", f"Created original BEHAV3D proportions QC: {path}")
                else:
                    adata_tracks = self._load_model_adata()
                    md = getattr(self.metadata_loader, "metadata", None)
                    all_cols = (group_cols or []) + [c for c in (group_x, group_y) if c]
                    cols_to_merge = [c for c in all_cols if c not in adata_tracks.obs.columns]
                    if cols_to_merge and md is not None:
                        merge_condition_columns_into_obs(adata_tracks, md, cols_to_merge)
                    paths = _resolve_dtaidistance_paths(self.output_dir, self._current_cell_type())
                    plot_paths = save_track_class_proportions_by_sample_plot(
                        adata_tracks,
                        paths["behavior_proportions_outfolder"],
                        sample_col="sample_name",
                        class_col="ClusterID",
                        group_cols=group_cols,
                        group_x=group_x,
                        group_y=group_y,
                        verbose=True,
                    )
                    self.plot_status_html.value = (
                        "<b>Track proportions ready:</b> per-sample class proportions were written."
                    )
                    _winfo(
                        "trajectory-dtai-widget",
                        f"Created track-class proportions plot: {plot_paths.get('pdf_path')}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_track_proportions, self.track_proportions_spinner, busy=False)

    def _sync_track_comparison_group_y_text(self):
        if hasattr(self, "track_comparison_group_y_text"):
            self.track_comparison_group_y_text.value = str(self.track_comparison_condition_col_dd.value or "")

    def _on_track_comparison_condition_col_changed(self, change):
        if change.get("name") == "value":
            self._sync_track_comparison_group_y_text()

    def _on_track_condition_comparison_clicked(self, _):
        self._set_busy(self.btn_track_condition_comparison, self.track_condition_comparison_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if bool(self.use_original_behav3d.value):
                    raise ValueError(
                        "Condition comparison plots are only available for the one-hot dtaidistance method."
                    )
                condition_col = self.track_comparison_condition_col_dd.value
                group_cols = list(self.track_comparison_group_cols_select.value) or None
                group_x = self.track_comparison_group_x_dd.value
                group_x = None if group_x in (None, "(none)") else group_x
                if not condition_col:
                    raise ValueError("Select a condition column to compare.")

                adata_tracks = self._load_model_adata()
                md = getattr(self.metadata_loader, "metadata", None)
                all_cols = [condition_col] + (group_cols or []) + [c for c in (group_x,) if c]
                cols_to_merge = [c for c in all_cols if c not in adata_tracks.obs.columns]
                if cols_to_merge and md is not None:
                    merge_condition_columns_into_obs(adata_tracks, md, cols_to_merge)
                paths = _resolve_dtaidistance_paths(self.output_dir, self._current_cell_type())
                result = save_track_condition_comparison_report(
                    adata_tracks,
                    paths["behavior_comparisons_outfolder"],
                    sample_col="sample_name",
                    class_col="ClusterID",
                    condition_col=condition_col,
                    group_cols=group_cols,
                    group_x=group_x,
                    verbose=True,
                )
                self.plot_status_html.value = "<b>Condition comparison ready:</b> plot was written."
                _winfo(
                    "trajectory-dtai-widget",
                    f"Created track condition comparison plot: {result.get('pdf_path')}",
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_track_condition_comparison, self.track_condition_comparison_spinner, busy=False)

    def _on_contact_analysis_clicked(self, _):
        self._set_busy(self.btn_contact_analysis, self.contact_analysis_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if bool(self.use_original_behav3d.value):
                    raise ValueError(
                        "Contact-based grouping is only available for the one-hot dtaidistance method."
                    )
                contact_col = self.contact_col_dd.value
                if not contact_col:
                    raise ValueError("Select a contact column to group tracks by.")
                min_bout_length = int(self.contact_min_bout_length.value)
                if min_bout_length < 1:
                    raise ValueError("Min. contiguous contact bout must be at least 1 timepoint.")

                adata_tracks = self._load_model_adata()
                df_timepoints = pd.read_csv(self._original_track_features_path())
                extra_group_cols = list(self.contact_group_cols_select.value) or None
                group_x = self.contact_group_x_dd.value
                group_x = None if group_x in (None, "(none)") else group_x
                group_y = self.contact_group_y_dd.value
                group_y = None if group_y in (None, "(none)") else group_y
                md = getattr(self.metadata_loader, "metadata", None)
                all_extra_cols = (extra_group_cols or []) + [c for c in (group_x, group_y) if c]
                cols_to_merge = [c for c in all_extra_cols if c not in adata_tracks.obs.columns]
                if cols_to_merge and md is not None:
                    merge_condition_columns_into_obs(adata_tracks, md, cols_to_merge)
                paths = _resolve_dtaidistance_paths(self.output_dir, self._current_cell_type())
                result = save_track_contact_group_analysis(
                    adata_tracks,
                    df_timepoints,
                    paths["outfolder"],
                    contact_col=contact_col,
                    min_bout_length=min_bout_length,
                    sample_col="sample_name",
                    class_col="ClusterID",
                    extra_group_cols=extra_group_cols,
                    group_x=group_x,
                    group_y=group_y,
                    verbose=True,
                )
                self.plot_status_html.value = (
                    "<b>Contact-based grouping ready:</b> proportion, condition-comparison, "
                    "and violin plots were written."
                )
                _winfo(
                    "trajectory-dtai-widget",
                    f"Created contact-group analysis bundle: {result.get('condition_comparison', {}).get('pdf_path')}",
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_contact_analysis, self.contact_analysis_spinner, busy=False)

    def _on_exemplars_clicked(self, _):
        self._set_busy(self.btn_exemplars, self.exemplar_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if not self._any_exemplar_outputs_selected():
                    raise ValueError("Select at least one exemplar output.")
                if not self._has_behavioral_states():
                    raise ValueError(
                        "No states have been defined. Run behavioral state classification before creating exemplar outputs."
                    )
                if bool(self.use_original_behav3d.value):
                    plot_paths = self._save_original_exemplar_plots()
                else:
                    adata_tracks = self._load_model_adata()
                    plot_paths = save_dtaidistance_exemplar_plots(
                        adata_tracks,
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        n_per_cluster=int(self.n_per_cluster.value),
                        make_overview_statebars=bool(self.make_overview_statebars.value),
                        make_backprojection_pdf=bool(self.make_backprojection_pdf.value),
                        make_backprojection_mp4=bool(self.make_backprojection_mp4.value),
                        random_state=int(self.random_state.value),
                        save_outputs=True,
                        verbose=True,
                    )
                self.plot_status_html.value = "<b>Exemplar PDFs ready:</b> exemplar track outputs were written."
                _winfo(
                    "trajectory-dtai-widget",
                    f"Created exemplar overview: {plot_paths.get('exemplar_tracks_overview_pdf')}",
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_exemplars, self.exemplar_spinner, busy=False)
                self._sync_exemplar_state_controls()

    def _on_backproject_clicked(self, _):
        self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=True)
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                adata_tracks = self._load_model_adata()
                sample_name = self.backproj_sample_dd.value
                if sample_name is None or len(str(sample_name).strip()) == 0:
                    raise ValueError("Please select a sample.")
                adata_full_path = self._state_adata_path()
                if not Path(adata_full_path).exists():
                    raise FileNotFoundError(f"Behavioral-state h5ad not found: {adata_full_path}")

                adata_full = sc.read_h5ad(adata_full_path)
                output_col = "dtai_track_behavioral_cluster"
                manifest = export_track_cluster_backprojection(
                    adata_full=adata_full,
                    adata_tracks=adata_tracks,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    cluster_col="ClusterID",
                    output_col=output_col,
                    sample_name=str(sample_name),
                    n_workers=max(1, int(self.backprojection_workers.value)),
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
                    "trajectory-dtai-widget",
                    f"Opening backprojection for sample '{sample_key}' | state_image={state_img_path}",
                )
                show_track_cluster_backprojection(
                    sample_name=sample_key,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    adata_full=adata_full,
                    adata_tracks=adata_tracks,
                    cluster_col="ClusterID",
                    state_col=FULL_STATE_COL,
                    state_img_path=state_img_path,
                    output_col=output_col,
                    run=True,
                    verbose=True,
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=False)
                self._refresh_backprojection_samples()

    def _classifier_artifact_path(self):
        from behav3d.analysis.behavior.track.bouts import (
            _resolve_track_classifier_path,
        )
        return _resolve_track_classifier_path(
            self.output_dir,
            self._current_cell_type(),
            output_subdir_name="behavorial_trajectories",
        )

    def _on_train_classifier_clicked(self, _):
        self._set_busy(self.btn_train_classifier, self.train_classifier_spinner, busy=True)
        self.out_train_classifier.clear_output()
        with self.out_train_classifier:
            try:
                if self.model_adata is None:
                    self.model_adata = self._load_model_adata()
                test_size_text = str(self.classifier_test_size.value).strip()
                test_size = float(test_size_text) if test_size_text else 0.2
                max_depth_text = str(self.classifier_max_depth.value).strip()
                max_depth = None if max_depth_text in ("", "None") else int(max_depth_text)
                max_features_val = str(self.classifier_max_features.value)
                max_features = None if max_features_val == "None" else max_features_val
                result = train_dtaidistance_trajectory_classifier(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    model_adata=self.model_adata,
                    classifier_n_estimators=int(self.classifier_n_estimators.value),
                    classifier_min_samples_leaf=int(self.classifier_min_samples_leaf.value),
                    classifier_min_samples_split=int(self.classifier_min_samples_split.value),
                    classifier_max_features=max_features,
                    classifier_max_depth=max_depth,
                    classifier_n_jobs=int(self.classifier_n_jobs.value),
                    validation_test_size=test_size,
                    random_state=int(self.random_state.value),
                    verbose=True,
                )
                saved_path = result.get("classifier_path", "")
                self.classifier_artifact_path.value = str(saved_path)
                _winfo("trajectory-dtai-widget", f"RF classifier saved: {saved_path}")
                fit_info = result.get("fit_info", {})
                tr_acc = fit_info.get("train_accuracy", "n/a")
                tr_acc_str = f"{tr_acc:.3f}" if isinstance(tr_acc, (int, float)) else str(tr_acc)
                _winfo(
                    "trajectory-dtai-widget",
                    f"Train accuracy: {tr_acc_str} | "
                    f"n_train: {fit_info.get('n_train_rows', 'n/a')} | "
                    f"classes: {fit_info.get('classes', [])}",
                )
                val = result.get("classifier_validation", {})
                test_metrics = (val.get("test_metrics") or {}) if val else {}
                if test_metrics:
                    acc = test_metrics.get("accuracy", "n/a")
                    bal = test_metrics.get("balanced_accuracy", "n/a")
                    acc_str = f"{acc:.3f}" if isinstance(acc, (int, float)) else str(acc)
                    bal_str = f"{bal:.3f}" if isinstance(bal, (int, float)) else str(bal)
                    _winfo(
                        "trajectory-dtai-widget",
                        f"Holdout accuracy: {acc_str} | balanced: {bal_str}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_train_classifier, self.train_classifier_spinner, busy=False)

    def _on_apply_classifier_clicked(self, _):
        self._set_busy(self.btn_apply_classifier, self.apply_classifier_spinner, busy=True)
        self.out_apply_classifier.clear_output()
        with self.out_apply_classifier:
            try:
                clf_path = str(self.classifier_artifact_path.value).strip() or str(self._classifier_artifact_path())
                states_path_text = str(self.classifier_apply_states_path.value).strip()
                states_path = states_path_text if states_path_text else str(self._state_adata_path())
                if not Path(clf_path).exists():
                    raise FileNotFoundError(f"Classifier artifact not found: {clf_path}")
                if not Path(states_path).exists():
                    raise FileNotFoundError(f"Behavioral states file not found: {states_path}")
                ct = self._current_cell_type()
                _winfo("trajectory-dtai-widget", f"Applying classifier: {clf_path}")
                apply_result = apply_track_classifier_to_subtracks(
                    output_dir=self.output_dir,
                    cell_type=ct,
                    classifier_artifact_or_path=clf_path,
                    adata_full_path=states_path,
                    output_subdir_name="behavorial_trajectories",
                    group_cols=list(self.apply_group_cols_select.value) or None,
                    metadata=getattr(self.metadata_loader, "metadata", None),
                    verbose=True,
                )
                out_path = apply_result.get("output_path", "")
                _winfo("trajectory-dtai-widget", f"Applied classifier, output: {out_path}")
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_apply_classifier, self.apply_classifier_spinner, busy=False)
