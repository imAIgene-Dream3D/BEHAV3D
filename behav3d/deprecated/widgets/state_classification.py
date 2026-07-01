import ipywidgets as widgets
from pathlib import Path
from behav3d.widgets.base_state_classification import BaseStateClassificationPanel, _winfo
from behav3d.deprecated.state_classification_clustering import (
    run_state_clustering,
    train_state_classifiers,
    rename_intrinsic_behavioral_clusters,
)

class StateClassificationPanel(BaseStateClassificationPanel):
    def __init__(self, metadata_loader, cell_type=None):
        self._build_train_section()
        super().__init__(metadata_loader=metadata_loader, cell_type=cell_type)

    def _build_steps(self):
        step_defs = [
            ("Primary dynamic state clustering (based on continuous features)", self.clustering_section),
            ("Rename primary dynamic state clusters", self.rename_intrinsic_section),
            ("Rename clusters assigned to binary groups", self.rename_full_section),
            ("Train classification", self.train_section),
            ("Apply classification", self.apply_section),
            ("Backprojection", self.backprojection_section),
        ]
        self._step_accordions = []
        for title, section in step_defs:
            acc = widgets.Accordion(children=[section], selected_index=None)
            acc.set_title(0, title)
            self._step_accordions.append(acc)
        self.steps = widgets.VBox(self._step_accordions)

    def _apply_cfg_defaults(self):
        super()._apply_cfg_defaults()
        cfg = self._panel_cfg()
        if not isinstance(cfg, dict):
            return
        
        self.n_estimators.value = int(cfg.get("classifier_n_estimators", self.n_estimators.value))
        self.min_samples_leaf.value = int(cfg.get("classifier_min_samples_leaf", self.min_samples_leaf.value))
        self.n_jobs.value = int(cfg.get("classifier_n_jobs", self.n_jobs.value))
        self.max_depth.value = str(cfg.get("classifier_max_depth", self.max_depth.value))
        self.min_samples_split.value = int(cfg.get("classifier_min_samples_split", self.min_samples_split.value))
        self.max_features.value = str(cfg.get("classifier_max_features", self.max_features.value))
        self.class_weight.value = str(cfg.get("classifier_class_weight", self.class_weight.value))
        self.validation_test_size.value = float(cfg.get("validation_test_size", self.validation_test_size.value))
        self.validation_random_state.value = str(
            cfg.get("validation_random_state", self.validation_random_state.value)
        )
        self.validation_stratify.value = bool(cfg.get("validation_stratify", self.validation_stratify.value))

    def _persist_current_settings(self):
        cfg = self._panel_cfg()
        if isinstance(cfg, dict):
            cfg.update(
                {
                    "classifier_n_estimators": int(self.n_estimators.value),
                    "classifier_min_samples_leaf": int(self.min_samples_leaf.value),
                    "classifier_n_jobs": int(self.n_jobs.value),
                    "classifier_max_depth": str(self.max_depth.value),
                    "classifier_min_samples_split": int(self.min_samples_split.value),
                    "classifier_max_features": str(self.max_features.value),
                    "classifier_class_weight": str(self.class_weight.value),
                    "validation_test_size": float(self.validation_test_size.value),
                    "validation_random_state": str(self.validation_random_state.value),
                    "validation_stratify": bool(self.validation_stratify.value),
                    "save_label_classifier": True,
                    "save_full_label_classifier": True,
                    "train_continuous_classifier": True,
                    "train_full_classifier": True,
                }
            )
        super()._persist_current_settings()
    def _refresh_enablement(self):
        super()._refresh_enablement()
        has_model = self.model_adata is not None
        self.btn_train.disabled = not has_model

        if has_model:
            self.train_status.value = (
                f"<b>Model loaded:</b> {self._model_adata_path().name} "
                f"({self.model_adata.n_obs} rows)"
            )
        else:
            self.train_status.value = "<i>Load or create model adata first.</i>"

    def _build_train_section(self):
        self.train_status = widgets.HTML("<i>Load or create model adata first.</i>")
        self.n_estimators = widgets.IntText(description="Trees", value=300, style={"description_width": "130px"})
        self.min_samples_leaf = widgets.IntText(description="Min leaf", value=2, style={"description_width": "130px"})
        self.n_jobs = widgets.IntText(description="n_jobs", value=-1, style={"description_width": "130px"})
        self.max_depth = widgets.Text(description="Max depth", value="", style={"description_width": "130px"})
        self.min_samples_split = widgets.IntText(
            description="Min split",
            value=2,
            style={"description_width": "130px"},
        )
        self.max_features = widgets.Text(description="Max feat", value="sqrt", style={"description_width": "130px"})
        self.class_weight = widgets.Text(description="Class weight", value="", style={"description_width": "130px"})
        self.validation_test_size = widgets.FloatText(
            description="Val size",
            value=0.05,
            style={"description_width": "130px"},
        )
        self.validation_random_state = widgets.Text(
            description="Val seed",
            value="",
            style={"description_width": "130px"},
        )
        self.validation_stratify = widgets.Checkbox(
            description="Stratify validation",
            value=True,
            indent=False,
        )

        self.btn_train = widgets.Button(
            description="Train classification",
            button_style="success",
            layout=widgets.Layout(width="170px"),
        )
        self.btn_train.on_click(self._on_train_clicked)
        self.train_spinner = widgets.HTML(value=spinning_loader)
        self.train_spinner.layout.display = "none"
        self.out_train = widgets.Output()

        classifier_params_box = widgets.VBox(
            [
                widgets.HBox([self.n_estimators, self.min_samples_leaf, self.n_jobs]),
                widgets.HBox([self.max_depth, self.min_samples_split, self.max_features]),
                widgets.HBox([self.class_weight, self.validation_test_size, self.validation_random_state]),
                self.validation_stratify,
            ]
        )
        self.classifier_params_accordion = widgets.Accordion(
            children=[classifier_params_box],
            selected_index=None,
        )
        self.classifier_params_accordion.set_title(0, "classifier parameters")

        self.train_section = widgets.VBox(
            [
                self.train_status,
                self.classifier_params_accordion,
                widgets.HBox([self.btn_train, self.train_spinner]),
                self.out_train,
            ]
        )

    def _on_cluster_clicked(self, _):
        self._set_busy(self.btn_cluster, self.cluster_spinner, busy=True)
        self.out_cluster.clear_output()
        with self.out_cluster:
            try:
                self._persist_current_settings()
                selected_features = self._selected_feature_columns()
                if len(selected_features) == 0:
                    raise ValueError("Select at least one feature before clustering.")
                selected_descriptive_features = self._selected_descriptive_features()
                if len(selected_descriptive_features) == 0:
                    raise ValueError("Select at least one window descriptive feature before clustering.")
                ct = self._current_cell_type()
                self.model_adata = run_state_clustering(
                    features=selected_features,
                    binary_features_to_group=self._selected_binary_columns(),
                    output_dir=self.output_dir,
                    cell_type=ct,
                    window_size=int(self.window_size.value),
                    min_spacing=self._parse_optional_int(self.min_spacing.value),
                    max_samples=self._parse_optional_int(self.max_samples.value),
                    n_neighbors=int(self.n_neighbors.value),
                    min_dist=float(self.min_dist.value),
                    resolution=float(self.resolution.value),
                    descriptive_features=selected_descriptive_features,
                    pca_var_selection=float(self.pca_var_selection.value),
                    clustering_method=str(self.clustering_method.value),
                    lower_quantile_cap=self._parse_optional_float(self.lower_quantile_cap.value),
                    upper_quantile_cap=self._parse_optional_float(self.upper_quantile_cap.value),
                    incomplete_window_policy=str(self.incomplete_window_policy.value),
                    random_state=int(self.random_state.value),
                    reuse_prepared_dataset=bool(self.reuse_prepared_dataset.value),
                    verbose=True,
                )
                self._save_model_adata()
                _winfo("state-widget", f"Saved model adata: {self._model_adata_path()}")
                self._rebuild_intrinsic_rename_rows()
                self._rebuild_full_rename_rows()
                # Keep compact behavior: open rename-intrinsic section next.
                self._open_step(1)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_cluster, self.cluster_spinner, busy=False)
                self._refresh_enablement()

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

                rename_intrinsic_behavioral_clusters(
                    adata=self.model_adata,
                    mapping=mapping,
                    binary_cols_to_merge=self._selected_binary_columns(),
                )
                self._save_model_adata(compression="lzf")
                self._rebuild_intrinsic_rename_rows()
                self._rebuild_full_rename_rows()
                _winfo("state-widget", f"Renamed intrinsic clusters and saved: {self._model_adata_path()}")
                # Open full mapping next.
                self._open_step(2)
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_rename_intrinsic, self.rename_intrinsic_spinner, busy=False)
                self._refresh_enablement()

    def _on_train_clicked(self, _):
        self._set_busy(self.btn_train, self.train_spinner, busy=True)
        self.out_train.clear_output()
        with self.out_train:
            try:
                self._persist_current_settings()
                if self.model_adata is None:
                    self._load_existing_model_if_available()
                if self.model_adata is None:
                    raise ValueError(
                        f"No model adata found at: {self._model_adata_path()}. "
                        "Run clustering first."
                    )

                max_depth = self._parse_optional_int(self.max_depth.value)
                class_weight_raw = str(self.class_weight.value).strip()
                class_weight = None if class_weight_raw == "" else class_weight_raw
                val_seed = self._parse_optional_int(self.validation_random_state.value)

                out = train_state_classifiers(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    model_adata=self.model_adata,
                    label_transfer_method="classifier",
                    classifier_backend="random_forest",
                    classifier_n_estimators=int(self.n_estimators.value),
                    classifier_min_samples_leaf=int(self.min_samples_leaf.value),
                    classifier_n_jobs=int(self.n_jobs.value),
                    classifier_max_depth=max_depth,
                    classifier_min_samples_split=int(self.min_samples_split.value),
                    classifier_max_features=str(self.max_features.value).strip(),
                    classifier_class_weight=class_weight,
                    save_label_classifier=True,
                    train_continuous_classifier=True,
                    train_full_classifier=True,
                    save_full_label_classifier=True,
                    validation_test_size=float(self.validation_test_size.value),
                    validation_random_state=val_seed,
                    validation_stratify=bool(self.validation_stratify.value),
                    verbose=True,
                )
                self._save_model_adata()
                _winfo("state-widget", "Training finished.")
                _winfo("state-widget", f"Intrinsic classifier: {out.get('partial_classifier_path', None)}")
                _winfo("state-widget", f"Full classifier: {out.get('full_classifier_path', None)}")
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_train, self.train_spinner, busy=False)
                self._refresh_enablement()

