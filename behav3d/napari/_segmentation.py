
import napari
from behav3d.napari._widgets import (
    HelpButton,
    make_help_row,
    browse_file_or_zarr,
    prompt_axis_order,
    resolve_external_path,
)
import yaml
from magicgui.widgets import create_widget
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QStackedWidget, QPushButton, QGroupBox, QFormLayout,
    QSpinBox, QDoubleSpinBox, QCheckBox, QFileDialog, QScrollArea,
    QLineEdit, QPlainTextEdit, QTextEdit, QMessageBox, QFrame, QGridLayout,
    QDialog, QDialogButtonBox, QTableWidget, QTableWidgetItem, QHeaderView,
    QSizePolicy,
)
from qtpy.QtCore import Qt
import pandas as pd
from pathlib import Path
import traceback

from behav3d.preprocessing.segmentation.convpaint_train import AnnotationLegendTab
import numpy as np
import dask.array as da
import shutil
import gc
import os
from functools import partial
# Scikit-learn imports for training
from sklearn.ensemble import RandomForestClassifier
import joblib


from behav3d.core.utils import convert_distance
from behav3d.napari._units import UnitGroupManager
from behav3d.io.images import load_image, get_image_shape, save_as_zarr, append_to_zarr, load_zarr, convert_label_file_to_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.preprocessing.segmentation.napari_pixelclassifier import (
    train_pixel_classifier,
    run_pixel_classifier_segmentation,
    segment_mask,
    postprocess_mask,
    make_features,
    _set_dead_mask_layer_color,
)
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
    fire_extra_callback,
)


# make_features is imported from napari_pixelclassifier and built per-sample
# inside _on_load_training_clicked using actual pixel sizes from metadata,
# so training and batch inference use the same physical-unit sigma values.

class SegmentationTab(QWidget):
    def __init__(self, viewer: napari.Viewer, metadata_loader):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._init_ui()

    def _init_ui(self):
        self.main_layout = QVBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)

        self.stack = QStackedWidget()
        self.main_layout.addWidget(self.stack)

        # -- Page 0: Placeholder --
        self.placeholder_page = QWidget()
        place_lay = QVBoxLayout(self.placeholder_page)
        self.placeholder_label = QLabel("Load metadata in the Data Preparation tab to see segmentation options.")
        self.placeholder_label.setAlignment(Qt.AlignCenter)
        self.placeholder_label.setStyleSheet("color: #888; font-style: italic; font-size: 14px; padding: 20px;")
        place_lay.addWidget(self.placeholder_label)
        self.stack.addWidget(self.placeholder_page)

        # -- Page 1: Main Content --
        self.main_content = QWidget()
        layout = QVBoxLayout(self.main_content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # Shared progress row used by every sub-page's Run button.  Each
        # sub-widget pulls this from its parent ``SegmentationTab`` via
        # the ``tab_progress_row`` constructor argument so all four
        # methods drive the same visual widget.  Created up-front so
        # sub-widgets can capture a reference during construction.
        self.progress_row = ProgressBarRow()

        # Method Selection
        method_group = QGroupBox("Segmentation Method")
        method_layout = QHBoxLayout()
        self.method_combo = QComboBox()

        methods = [
            "APOC (GPU)",
            "ConvPaint (DL pixel classifier)",
            "Pixel Classifier (Random Forest)",
            "Cellpose (Deep Learning)",
            "Import segmentation",
        ]
        self.method_combo.addItems(methods)
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.method_combo)
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # Stacked Widget for Method Parameters
        self.param_stack = QStackedWidget()
        
        # 0. APOC (GPU) Page  ← matches combo index 0
        self.apoc_page = APOCWidget(
            self.viewer,
            self.metadata_loader,
            log_callback=self._log,
            switch_to_visualization_callback=self._switch_to_visualization,
            tab_progress_row=self.progress_row,
        )
        self.param_stack.addWidget(self.apoc_page)

        # 1. ConvPaint Page  ← matches combo index 1
        self.convpaint_page = ConvPaintWidget(
            self.viewer,
            self.metadata_loader,
            log_callback=self._log,
            switch_to_visualization_callback=self._switch_to_visualization,
            tab_progress_row=self.progress_row,
        )
        self.param_stack.addWidget(self.convpaint_page)

        # 2. Pixel Classifier Page
        self.pixel_classifier_page = PixelClassifierWidget(
            self.viewer,
            self.metadata_loader,
            log_callback=self._log,
            tab_progress_row=self.progress_row,
        )
        self.param_stack.addWidget(self.pixel_classifier_page)

        # 2. Cellpose Page  ← matches combo index 2
        self.cellpose_page = CellposeWidget(
            self.viewer,
            self.metadata_loader,
            log_callback=self._log,
            tab_progress_row=self.progress_row,
        )
        self.param_stack.addWidget(self.cellpose_page)
        
        # 3. Import Page ← matches combo index 3
        self.import_page = ImportWidget(
            self.viewer,
            self.metadata_loader,
            log_callback=self._log,
            switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
        )
        self.param_stack.addWidget(self.import_page)


        layout.addWidget(self.param_stack)
        
        # Add stretch to keep things at the top
        # layout.addStretch()

        # Add the progress row underneath the parameter stack so it sits
        # right above the log.
        layout.addWidget(self.progress_row)

        # Log Window
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(120)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log)

        self.stack.addWidget(self.main_content)
        
        # Connect to metadata signal
        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_loaded)
            
        # Initial state
        if self.metadata_loader.metadata is not None:
            self.stack.setCurrentIndex(1)
        else:
            self.stack.setCurrentIndex(0)

    def _on_metadata_loaded(self, metadata):
        """Metadata signal received: show options and refresh sub-widgets."""
        self.stack.setCurrentIndex(1)
        self._on_metadata_updated()

    def _on_metadata_updated(self):
        """Trigger updates in sub-widgets."""
        if hasattr(self, 'pixel_classifier_page'):
            self.pixel_classifier_page._on_metadata_updated()
        self.cellpose_page._on_metadata_updated()
        self.import_page._on_metadata_updated()
        if hasattr(self, 'apoc_page'):
            self.apoc_page._on_metadata_updated()
        if hasattr(self, 'convpaint_page'):
            self.convpaint_page._on_metadata_updated()

    def _log(self, msg):
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        self.log.append(f"[{timestamp}] {msg}")
        self.log.verticalScrollBar().setValue(
            self.log.verticalScrollBar().maximum()
        )

    def _switch_to_data_prep_edit(self):
        """Switch the main window to the Data Preparation tab with the
        Metadata Builder already open in edit mode.

        Used by ``ImportWidget`` (the "Add a new sample or cell type"
        shortcut) so users don't need to duplicate the sample/cell-type
        editing UI here — they jump straight to it instead.
        """
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(0)
            if hasattr(parent, 'data_prep_tab'):
                parent.data_prep_tab.enter_metadata_edit_mode()

    def _switch_to_visualization(self):
        """Switch the main window to the Visualization tab and load the first sample.

        Used by sub-widgets (e.g. APOCWidget) to ask the host window to
        navigate without walking the parent chain.
        """
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(1)
            if hasattr(parent, 'visualization_tab'):
                self._log("  Loading first sample in Visualization Tab...")
                parent.visualization_tab.sample_combo.setCurrentIndex(0)
                parent.visualization_tab._on_load_dataset()
                for layer in self.viewer.layers:
                    if "Segments" in layer.name:
                        layer.visible = True

    def _on_method_changed(self, index):
        current_idx = self.param_stack.currentIndex()

        # If switching away from Pixel Classifier (index 2), check session
        if current_idx == 2 and self.pixel_classifier_page.is_session_active:
            from qtpy.QtWidgets import QMessageBox
            res = QMessageBox.question(
                self,
                "Switch Method?",
                "A segmentation session is active. Switching methods will clear the viewer "
                "and reset your unsaved labeling session.\n\nDo you want to proceed?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                self.method_combo.blockSignals(True)
                self.method_combo.setCurrentIndex(current_idx)
                self.method_combo.blockSignals(False)
                return
            else:
                self.cleanup_session()

        # If switching away from APOC tab (index 0), remove training layers
        elif current_idx == 0 and self.apoc_page._is_session_active:
            from qtpy.QtWidgets import QMessageBox
            res = QMessageBox.question(
                self,
                "Switch Method?",
                "APOC training layers are loaded in the viewer. Switching methods will remove them.\n\nDo you want to proceed?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                self.method_combo.blockSignals(True)
                self.method_combo.setCurrentIndex(current_idx)
                self.method_combo.blockSignals(False)
                return
            else:
                self.apoc_page.cleanup_session()

        # If switching away from ConvPaint tab (index 1), remove training layers
        elif current_idx == 1 and getattr(self.convpaint_page, "_is_session_active", False):
            from qtpy.QtWidgets import QMessageBox
            res = QMessageBox.question(
                self,
                "Switch Method?",
                "ConvPaint training layers are loaded in the viewer. Switching methods will remove them.\n\nDo you want to proceed?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                self.method_combo.blockSignals(True)
                self.method_combo.setCurrentIndex(current_idx)
                self.method_combo.blockSignals(False)
                return
            else:
                self.convpaint_page.cleanup_session()

        self.param_stack.setCurrentIndex(index)

        # If switching to APOC tab (index 0), refresh it
        if index == 0:
            self.apoc_page._on_metadata_updated()
        elif index == 1:
            self.convpaint_page._on_metadata_updated()
        # If switching to Import tab (index 4), refresh it to avoid outdated info
        elif index == 4:
            self.import_page._on_metadata_updated()

    def request_tab_exit(self):
        """Called by the main widget when the user tries to leave this tab."""
        from qtpy.QtWidgets import QMessageBox

        # 0. Block tab switching while any background run is in flight.
        for page_name in (
            "pixel_classifier_page",
            "cellpose_page",
            "apoc_page",
            "convpaint_page",
        ):
            page = getattr(self, page_name, None)
            bg = getattr(page, "_bg", None) if page is not None else None
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "A segmentation run is still in progress. Please wait "
                    "for it to finish before switching tabs.",
                )
                return False

        # 1. Pixel Classifier Page
        if self.pixel_classifier_page.is_session_active:
            res = QMessageBox.question(
                self,
                "Leave Segmentation Tab?",
                "A segmentation session is active. Switching tabs will clear the viewer "
                "and reset your unsaved labeling session.\n\nDo you want to leave?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                return False
            save_res = QMessageBox.question(
                self, "Save Labels?",
                "Do you want to save your user labels before leaving?",
                QMessageBox.Yes | QMessageBox.No
            )
            if save_res == QMessageBox.Yes:
                self.pixel_classifier_page._save_user_labels()
            self.cleanup_session()
            return True

        # 2. APOC Page
        if self.apoc_page._is_session_active:
            res = QMessageBox.question(
                self,
                "Leave Segmentation Tab?",
                "APOC training layers are loaded in the viewer. Leaving this tab will remove them.\n\nDo you want to leave?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                return False
            save_res = QMessageBox.question(
                self, "Save Labels?",
                "Do you want to save your user labels before leaving?",
                QMessageBox.Yes | QMessageBox.No
            )
            if save_res == QMessageBox.Yes:
                self.apoc_page._save_apoc_user_labels()
            self.apoc_page.cleanup_session()
            return True

        # 3. ConvPaint Page
        if getattr(self.convpaint_page, "_is_session_active", False):
            res = QMessageBox.question(
                self,
                "Leave Segmentation Tab?",
                "ConvPaint training layers are loaded in the viewer. Leaving this tab will remove them.\n\nDo you want to leave?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.No:
                return False
            save_res = QMessageBox.question(
                self, "Save Labels?",
                "Do you want to save your user labels before leaving?",
                QMessageBox.Yes | QMessageBox.No
            )
            if save_res == QMessageBox.Yes:
                tw = getattr(self.convpaint_page, "_training_widget", None)
                if tw is not None:
                    tw.save_user_labels()
            self.convpaint_page.cleanup_session()
            return True

        return True

    def cleanup_session(self):
        """Clear viewer and reset tab state (Pixel Classifier, APOC, ConvPaint)."""
        self.pixel_classifier_page.reset_ui()
        self.apoc_page.cleanup_session()
        if hasattr(self, "convpaint_page"):
            self.convpaint_page.cleanup_session()

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    """Get a value from a nested dict using dotted key notation."""
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


# Layer-name patterns used by all three segmentation methods (CPU Pixel
# Classifier, APOC, ConvPaint).  Any layer whose name matches one of these
# is considered "owned" by the segmentation session and must be removed
# when the session ends, so the viewer is clean for the next tab/method.
_SEGMENTATION_LAYER_PATTERNS = (
    "Channel ",                # input channel images: "Channel 0", "Channel 1", ...
    "User Provided Labels",    # user-painted labels
    "Pixel Classification",    # CPU/APOC prediction-mask labels
    "Probability Map",         # APOC / ConvPaint probability maps
    " Segments",               # instance segmentation labels: "<CellType> Segments"
)


def _remove_segmentation_layers(viewer) -> int:
    """Remove every layer matching the shared segmentation patterns.

    Returns the number of layers removed (mostly for logging).  Safe to
    call repeatedly; missing layers are silently ignored.
    """
    if viewer is None:
        return 0
    to_remove = []
    for layer in list(viewer.layers):
        name = getattr(layer, "name", "") or ""
        for pat in _SEGMENTATION_LAYER_PATTERNS:
            if pat in name:
                to_remove.append(layer)
                break
    for layer in to_remove:
        try:
            viewer.layers.remove(layer)
        except Exception:
            pass
    return len(to_remove)


class PerClassLegendWidget(QWidget):
    """Clickable legend for Random Forest and APOC that use one layer per cell type."""
    def __init__(self, viewer, all_cell_types, has_death=False, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.all_cell_types = all_cell_types
        self.has_death = has_death

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        intro = QLabel(
            "<b>Annotation legend.</b> Click a button to select the target layer and set the active brush (Eraser, Background, or Foreground)."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #888; font-size: 11px;")
        layout.addWidget(intro)

        def create_button(text, color_hex, layer_name, label_idx):
            btn = QPushButton(text)
            text_color = "white" if color_hex not in ("none", "white", "#ffffff", "transparent") else "black"
            bg_color = "transparent" if color_hex == "none" else color_hex
            btn.setStyleSheet(f"background-color: {bg_color}; color: {text_color}; font-size: 11px; padding: 4px; border: 1px solid #555; border-radius: 2px;")
            btn.clicked.connect(lambda _, ln=layer_name, idx=label_idx: self._select_brush(ln, idx))
            return btn

        if self.has_death:
            row_layout = QHBoxLayout()
            row_layout.addWidget(QLabel("<b>Dead</b>"), stretch=1)
            layer_name = 'User Provided Labels (Dead)'
            row_layout.addWidget(create_button("0: Eraser", "none", layer_name, 0))
            row_layout.addWidget(create_button("1: Background", "#8b3a26", layer_name, 1))
            row_layout.addWidget(create_button("2: Foreground", "#00ffff", layer_name, 2))
            layout.addLayout(row_layout)

        for ct in self.all_cell_types:
            row_layout = QHBoxLayout()
            row_layout.addWidget(QLabel(f"<b>{ct.capitalize()}</b>"), stretch=1)
            layer_name = f'User Provided Labels ({ct.capitalize()})'
            row_layout.addWidget(create_button("0: Eraser", "none", layer_name, 0))
            row_layout.addWidget(create_button("1: Background", "#8b3a26", layer_name, 1))
            row_layout.addWidget(create_button("2: Foreground", "#00ffff", layer_name, 2))
            layout.addLayout(row_layout)

        layout.addStretch()
        self.setLayout(layout)

    def _select_brush(self, layer_name, label_idx):
        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == layer_name:
                layer = lyr
                break
        if layer is None:
            return
        
        layer.selected_label = label_idx
        try:
            layer.mode = "paint"
        except Exception:
            pass
        try:
            self.viewer.layers.selection.active = layer
        except Exception:
            try:
                self.viewer.layers.selection = {layer}
            except Exception:
                pass


def _run_multicolor_cleanup_if_needed(metadata_loader, log_fn, skip_unified=False):
    """Run multicolor redundancy cleanup after segmentation."""
    if skip_unified:
        log_fn(
            "Skipping multicolor cleanup: unified ConvPaint already "
            "produces mutually-exclusive per-cell-type masks."
        )
        return metadata_loader.metadata

    metadata = metadata_loader.metadata
    if metadata is None or metadata.empty:
        return metadata

    from behav3d.core.metadata import (
        is_multicolor_celltype, 
        multicolor_base_name,
        detect_organoid_types_from_metadata,
        detect_immune_cell_types_from_metadata,
        detect_other_cell_types_from_metadata
    )
    all_cell_types = (
        detect_organoid_types_from_metadata(metadata) +
        detect_immune_cell_types_from_metadata(metadata) +
        detect_other_cell_types_from_metadata(metadata)
    )

    families = {}
    for cell_type in all_cell_types:
        if not is_multicolor_celltype(cell_type):
            continue
        base_name = multicolor_base_name(cell_type)
        families.setdefault(base_name, []).append(cell_type)

    if not families:
        return metadata

    from behav3d.preprocessing.segmentation.multicolor_segment_processing import (
        apply_multicolor_segment_correction_for_base,
    )

    cleaned_metadata = metadata
    for base_name, family_cell_types in sorted(families.items()):
        if len(family_cell_types) < 2:
            continue
        log_fn(
            f"Running multicolor cleanup for '{base_name}' after batch segmentation: "
            f"{sorted(family_cell_types)}"
        )
        cleaned_metadata = apply_multicolor_segment_correction_for_base(
            metadata=cleaned_metadata,
            output_dir=str(metadata_loader.output_dir),
            base_cell_type=base_name,
            n_channels=len(family_cell_types),
            overwrite=True,
            n_workers=1,
        )

    return cleaned_metadata


class PixelClassifierWidget(QWidget):
    def __init__(self, viewer, metadata_loader, log_callback=None,
                 tab_progress_row=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback if log_callback else print
        # Background-execution infrastructure (shared progress row sits
        # on the parent SegmentationTab).
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        self._detect_cell_types()
        self._init_ui()
        self.all_images = None
        self.all_features = None
        self.is_session_active = False

    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
            is_combined_multicolor_celltype,
        )
        metadata = self.metadata_loader.metadata
        
        if metadata is None:
            self.organoid_types = []
            self.immune_types = []
            self.other_types = []
            self.has_death = False
            self.all_cell_types = []
            return

        self.organoid_types = [ct for ct in detect_organoid_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.immune_types = [ct for ct in detect_immune_cell_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.other_types = [ct for ct in detect_other_cell_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        content = QWidget()
        content_layout = QVBoxLayout()
        content_layout.setContentsMargins(4, 4, 4, 4)
        content_layout.setSpacing(4)
        content.setLayout(content_layout)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        # Read saved config values
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        # Training Configuration
        train_group = QGroupBox("1. Labeling & Training")
        train_form = QFormLayout()
        train_form.setContentsMargins(6, 6, 6, 6)
        train_form.setSpacing(4)
        train_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        
        self.spin_examples = QSpinBox()
        self.spin_examples.setValue(int(pc.get("examples_per_sample", 3)))
        self.spin_examples.setRange(1, 10)
        self.spin_examples.setMaximumWidth(70)
        train_form.addRow("Examples/sample:", make_help_row(
            self.spin_examples,
            "Examples per Sample",
            "Number of timepoints selected from each dataset for labeling and "
            "training the pixel classifier.\n\n"
            "Timepoints are chosen evenly spaced across the time series "
            "(including the first and last), not at random.\n\n"
            "More examples → better generalization but more labeling work "
            "and slower training."
        ))
        

        
        self.btn_load_training = QPushButton("Generate Training Data")
        self.btn_load_training.setToolTip("Clears viewer and loads selected timepoints for labeling")
        self.btn_load_training.setStyleSheet("background-color: #007bff; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_load_training.clicked.connect(lambda _: self._on_load_training_clicked(interactive=True))
        
        self.btn_save_labels = QPushButton("Save Labels")
        self.btn_save_labels.setToolTip("Save all user-provided label layers to disk")
        self.btn_save_labels.setStyleSheet("background-color: #fd7e14; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_save_labels.clicked.connect(self._save_user_labels)
        
        self.btn_clear_layer = QPushButton("Clear Layer")
        self.btn_clear_layer.setToolTip("Erase painted labels in the currently selected user label layer")
        self.btn_clear_layer.setStyleSheet("background-color: #dc3545; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_clear_layer.clicked.connect(self._clear_selected_user_labels)
        
        self.btn_clear_all = QPushButton("Clear All")
        self.btn_clear_all.setToolTip("Erase painted labels from ALL user label layers")
        self.btn_clear_all.setStyleSheet("background-color: #dc3545; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_clear_all.clicked.connect(self._clear_all_user_labels)
        
        self.btn_train_update = QPushButton("Train Classifier")
        self.btn_train_update.setStyleSheet("background-color: #28a745; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_train_update.clicked.connect(self._on_train_clicked)
        
        train_btn_row = QHBoxLayout()
        train_btn_row.setSpacing(4)
        train_btn_row.addWidget(self.btn_save_labels)
        train_btn_row.addWidget(self.btn_clear_layer)
        train_btn_row.addWidget(self.btn_clear_all)
        train_btn_row.addWidget(self.btn_train_update)

        # +🛒 button for Train
        self.btn_queue_train = QPushButton("+🛒")
        self.btn_queue_train.setFixedSize(36, 28)
        self.btn_queue_train.setToolTip("Add Train Classifier to Processing Queue")
        self.btn_queue_train.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )
        self.btn_queue_train.setVisible(False)
        train_btn_row.addWidget(self.btn_queue_train)
        
        train_layout = QVBoxLayout()
        train_layout.addLayout(train_form)
        train_layout.addWidget(self.btn_load_training)

        # Instructional label for labeling
        self.instruction_label = QLabel("Select a 'User Provided Labels ...' layer and annotate your ground truth data in the left panel")
        self.instruction_label.setWordWrap(True)
        self.instruction_label.setStyleSheet("color: white; font-style: italic; margin-top: 4px;")
        self.instruction_label.setVisible(False)
        train_layout.addWidget(self.instruction_label)

        # Add training related buttons
        train_layout.addLayout(train_btn_row)

        # Placeholder for dynamic legend
        self.legend_container = QWidget()
        self.legend_container.setLayout(QVBoxLayout())
        self.legend_container.layout().setContentsMargins(0, 0, 0, 0)
        self.legend_container.setVisible(False)
        train_layout.addWidget(self.legend_container)

        # Hide action buttons until training data is loaded
        self.btn_save_labels.setVisible(False)
        self.btn_clear_layer.setVisible(False)
        self.btn_clear_all.setVisible(False)
        self.btn_train_update.setVisible(False)
        self.instruction_label.setVisible(False)
        
        train_group.setLayout(train_layout)
        
        content_layout.addWidget(train_group)

        # Segmentation Parameters (Dynamic)
        self.param_group = QGroupBox("2. Fine-tune Segmentation Parameters (Optional)")
        self.param_layout = QVBoxLayout()
        self.param_group.setLayout(self.param_layout)

        # Workers spinner (applies to batch segmentation)
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        self.spin_workers = QSpinBox()
        default_workers = min(int(pc.get("workers", n_cores)), max_allowed)
        self.spin_workers.setValue(default_workers)
        self.spin_workers.setRange(1, max_allowed)
        self.spin_workers.setMaximumWidth(70)
        self.spin_workers.valueChanged.connect(self._on_workers_changed)
        self.workers_form = QFormLayout()
        self.workers_form.setContentsMargins(6, 6, 6, 6)
        self.workers_form.setSpacing(4)
        n_cores = os.cpu_count() or 4
        self.workers_form.addRow("Workers:", make_help_row(
            self.spin_workers,
            "Workers",
            f"Number of parallel CPU threads to use.\n\n"
            f"Your machine has {n_cores} cores.\n"
            f"Recommendation: Use at most {max(1, n_cores - 1)} cores to keep the system responsive."
        ))
        self.param_layout.addLayout(self.workers_form)

        self.btn_resegment = QPushButton("Test Segmentation Parameters")
        self.btn_resegment.setToolTip("Re-runs segmentation on the current slice/time point using current parameters.")
        self.btn_resegment.setStyleSheet("background-color: #17a2b8; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_resegment.clicked.connect(self._on_resegment_clicked)
        
        self.param_widgets = {} # {cell_type: {param_name: widget}}
        self._refresh_params_ui() # Populate initially
        
        content_layout.addWidget(self.param_group)

        # Execution Configuration
        seg_group = QGroupBox("3. Batch Segmentation")
        seg_form = QFormLayout()
        seg_form.setContentsMargins(6, 6, 6, 6)
        seg_form.setSpacing(4)
        seg_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        
        self.check_process_all = QCheckBox("All timepoints")
        self.check_process_all.setChecked(bool(pc.get("use_all_timepoints", True)))
        self.check_process_all.toggled.connect(self._on_process_all_toggled)
        seg_form.addRow(self.check_process_all)
        
        use_all = bool(pc.get("use_all_timepoints", True))
        self.spin_t_start = QSpinBox()
        self.spin_t_start.setRange(0, 9999)
        self.spin_t_start.setValue(int(pc.get("tp_start", 0)))
        self.spin_t_start.setEnabled(not use_all)
        self.spin_t_start.setMaximumWidth(60)
        self.spin_t_end = QSpinBox()
        self.spin_t_end.setRange(0, 9999)
        self.spin_t_end.setValue(int(pc.get("tp_end", 0)))
        self.spin_t_end.setEnabled(not use_all)
        self.spin_t_end.setMaximumWidth(60)
        
        t_range_layout = QHBoxLayout()
        t_range_layout.setSpacing(4)
        t_range_layout.addWidget(QLabel("T:"))
        t_range_layout.addWidget(self.spin_t_start)
        t_range_layout.addWidget(QLabel("–"))
        t_range_layout.addWidget(self.spin_t_end)
        t_range_layout.addStretch()
        seg_form.addRow("Range:", t_range_layout)
        
        self.btn_run_segmentation = QPushButton("Run Batch Segmentation")
        self.btn_run_segmentation.setStyleSheet("background-color: #28a745; color: white; font-weight: bold; border-radius: 4px; padding: 6px; font-size: 12px;")
        self.btn_run_segmentation.clicked.connect(self._on_run_segmentation_clicked)

        # +🛒 button for Batch Segmentation
        self.btn_queue_segment = QPushButton("+🛒")
        self.btn_queue_segment.setFixedSize(36, 28)
        self.btn_queue_segment.setToolTip("Add Batch Segmentation to Processing Queue")
        self.btn_queue_segment.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )
        self.btn_queue_segment.setVisible(False)

        seg_btn_row = QHBoxLayout()
        seg_btn_row.setSpacing(4)
        seg_btn_row.addWidget(self.btn_run_segmentation, stretch=1)
        seg_btn_row.addWidget(self.btn_queue_segment)
        
        seg_layout = QVBoxLayout()
        seg_layout.addLayout(seg_form)
        seg_layout.addLayout(seg_btn_row)
        seg_group.setLayout(seg_layout)
        
        content_layout.addWidget(seg_group)
        content_layout.addStretch()

    def _on_metadata_updated(self):
        self.log("Metadata updated, refreshing Segmentation Tab...")
        self._detect_cell_types()
        self._refresh_params_ui()
        # Show segmentation queue button once metadata is loaded
        self.btn_queue_segment.setVisible(True)

    def _refresh_params_ui(self):
        # Collect all items first, then process
        items = []
        while self.param_layout.count():
            items.append(self.param_layout.takeAt(0))
        
        # Delete dynamic items, keep workers_form and btn_resegment
        for item in items:
            widget = item.widget()
            layout = item.layout()
            if widget and widget is self.btn_resegment:
                pass  # keep it
            elif layout and layout is self.workers_form:
                pass  # keep it
            elif widget:
                widget.deleteLater()
            elif layout:
                while layout.count():
                    child = layout.takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
        
        # Re-add preserved items (workers first, resegment re-added at end)
        self.param_layout.addLayout(self.workers_form)
        
        self.param_widgets = {}

        if not self.all_cell_types:
            self.param_layout.addWidget(QLabel("No cell types detected or no metadata loaded."))
            return

        # Read saved per-cell-type parameters from config
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        for cell_type in self.all_cell_types:
            ct_group = QGroupBox(f"{cell_type.capitalize()}")
            ct_form = QFormLayout()
            ct_form.setContentsMargins(6, 4, 6, 4)
            ct_form.setSpacing(3)
            ct_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

            # Per-group physical(µm)/pixel unit toggle. EDT & min-size are
            # entered in the chosen unit and converted to native px/voxels
            # before running/persisting; opening stays fixed in pixels.
            unit_mgr = UnitGroupManager(self.metadata_loader.metadata, default_physical=True)
            ct_form.addRow("Units:", unit_mgr.header_row(label=""))

            # Determine defaults based on category
            if cell_type in self.organoid_types:
                def_edt, def_size, def_open, def_fill = 12.0, 1000, 3, True
            elif cell_type in self.immune_types:
                def_edt, def_size, def_open, def_fill = 2.5, 10, 0, True
            else:
                def_edt, def_size, def_open, def_fill = 1.0, 10, 0, True
            
            # Override defaults with saved config values if present
            saved_edt = pc.get(f"{cell_type}_edt_threshold", def_edt)
            saved_size = pc.get(f"{cell_type}_segment_size_min", def_size)
            saved_open = pc.get(f"{cell_type}_opening_nr_pixels", def_open)
            saved_fill = pc.get(f"{cell_type}_fill_holes", def_fill)

            w_edt = QDoubleSpinBox()
            # Range widened so physical-unit display never clamps.
            w_edt.setRange(0.0, 1000.0)
            w_edt.setSingleStep(0.5)
            w_edt.setValue(float(saved_edt))
            w_edt.setMaximumWidth(90)
            ct_form.addRow("EDT:", make_help_row(
                w_edt,
                "EDT (Distance Transform Threshold)",
                "Controls sensitivity for separating touching objects.\n\n"
                "Lower values → more sensitive, splits objects more aggressively.\n"
                "Higher values → less splitting, objects stay merged.\n\n"
                "Shown in the unit selected above (µm by default; converted "
                "using the image's pixel size). Typical values in pixels:\n"
                "  • Organoids: 8–15 px\n"
                "  • Immune cells: 1.5–4.0 px\n"
                "  • Disable: 0.0 (not recommended, leads to under-segmentation)"
            ))

            w_size = QSpinBox()
            w_size.setRange(1, 1000000000)
            w_size.setValue(int(saved_size))
            w_size.setMaximumWidth(90)
            ct_form.addRow("Min size:", make_help_row(
                w_size,
                "Minimum Object Size",
                "Minimum volume for a segmented object to be kept.\n\n"
                "Objects smaller than this are removed as noise.\n\n"
                "Shown in the unit selected above (µm³ by default; converted "
                "using the image's pixel size). Typical values in voxels:\n"
                "  • Organoids: 500–2000 voxels\n"
                "  • Immune cells: 5–50 voxels"
            ))

            w_open = QSpinBox()
            w_open.setRange(0, 10)
            w_open.setValue(int(saved_open))
            w_open.setMaximumWidth(60)
            ct_form.addRow("Opening:", make_help_row(
                w_open,
                "Morphological Opening (pixels)",
                "Radius of morphological opening applied to the mask.\n\n"
                "Smooths edges and removes small protrusions/noise.\n"
                "0 = no opening applied.\n\n"
                "Increase for smoother boundaries; keep low for small objects."
            ))

            w_fill = QCheckBox("Fill holes")
            w_fill.setChecked(bool(saved_fill))
            fill_row = QHBoxLayout()
            fill_row.addWidget(w_fill)
            fill_row.addWidget(HelpButton(
                "Fill Holes",
                "When checked, internal gaps or holes within segmented "
                "objects are automatically filled.\n\n"
                "Recommended ON for most use cases."
            ))
            fill_row.addStretch()
            ct_form.addRow("", fill_row)

            # Register distance/volume spinboxes with the unit toggle so they
            # display in the chosen unit; opening/fill-holes stay native.
            unit_mgr.register(w_edt, "distance", float(saved_edt))
            unit_mgr.register(w_size, "volume", int(saved_size))

            ct_group.setLayout(ct_form)
            self.param_layout.addWidget(ct_group)

            self.param_widgets[cell_type] = {
                'edt': w_edt,
                'min_size': w_size,
                'opening': w_open,
                'fill_holes': w_fill,
                'unit_mgr': unit_mgr,
            }
        
        # Add Test Button at the end of params list
        self.param_layout.addWidget(self.btn_resegment)

    def _persist_params(self):
        """Save all segmentation widget values to behav3d_parameters.yml."""
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["examples_per_sample"] = int(self.spin_examples.value())

        pc["workers"] = int(self.spin_workers.value())
        pc["use_all_timepoints"] = bool(self.check_process_all.isChecked())
        pc["tp_start"] = int(self.spin_t_start.value())
        pc["tp_end"] = int(self.spin_t_end.value())

        # Save per-cell-type parameters
        for cell_type in self.all_cell_types:
            if cell_type in self.param_widgets:
                w = self.param_widgets[cell_type]
                mgr = w.get('unit_mgr')
                if mgr is not None:
                    # Persist canonical native (pixel/voxel) values regardless
                    # of the current display unit.
                    edt_native = mgr.get_native(w['edt'])
                    size_native = mgr.get_native(w['min_size'])
                else:
                    edt_native = float(w['edt'].value())
                    size_native = int(w['min_size'].value())
                pc[f"{cell_type}_edt_threshold"] = float(edt_native)
                pc[f"{cell_type}_segment_size_min"] = int(round(size_native))
                pc[f"{cell_type}_opening_nr_pixels"] = int(w['opening'].value())
                pc[f"{cell_type}_fill_holes"] = bool(w['fill_holes'].isChecked())

        # Write to YAML file
        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    def _save_user_labels(self):
        """Save all user-provided label layers to zarr on disk."""
        if not self.metadata_loader.output_dir:
            self.log("Cannot save labels: no output directory set.")
            return

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        saved_any = False

        if self.has_death:
            layer_name = 'User Provided Labels (Dead)'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved Dead labels to {outpath.name}")
                saved_any = True

        for cell_type in self.all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved {cell_type} labels to {outpath.name}")
                saved_any = True

        if not saved_any:
            self.log("No user label layers found in viewer to save.")
        else:
            self.log("All user labels saved.")

    def reset_ui(self):
        """Hide action buttons and clear session data."""
        self.btn_save_labels.setVisible(False)
        self.btn_clear_layer.setVisible(False)
        self.btn_clear_all.setVisible(False)
        self.btn_train_update.setVisible(False)
        # Delete elements inside legend_container
        while self.legend_container.layout().count():
            child = self.legend_container.layout().takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.legend_container.setVisible(False)
        self.instruction_label.setVisible(False)
        self.btn_queue_train.setVisible(False)
        
        self.all_images = None
        self.all_features = None
        self.is_session_active = False

        n_removed = _remove_segmentation_layers(self.viewer)
        self.log(f"Segmentation session reset ({n_removed} layers removed).")

    def _configure_user_label_layer(self, layer):
        """Configure a user label layer for simplified pixel classifier usage.
        
        - Additive blending
        - Labels restricted to 0 (Eraser), 1 (Background), 2 (Foreground)
        - Default selected label = 1 (Background)
        - Custom colors: 1=blue (background), 2=red (foreground)
        """
        import numpy as np

        # Additive blending
        layer.blending = "additive"

        # Custom color map: 0=transparent, 1=red (bg), 2=cyan (fg)
        layer.color = {
            0: (0, 0, 0, 0),       # Eraser / transparent
            1: (1.0, 0.2, 0.2, 1), # Background – red
            2: (0.0, 1.0, 1.0, 1), # Foreground – cyan
        }

        # Default selected label
        layer.selected_label = 1

        # Clamp selected_label to 0–2 whenever user changes it
        def _clamp_label(event):
            if layer.selected_label > 2:
                layer.selected_label = 2
        layer.events.selected_label.connect(_clamp_label)

    def _clear_selected_user_labels(self):
        """Zero out the currently selected user label layer."""
        import numpy as np
        import napari

        active = self.viewer.layers.selection.active
        if active is None or not isinstance(active, napari.layers.Labels):
            self.log("Please select a User Provided Labels layer first.")
            return
        if not active.name.startswith("User Provided Labels"):
            self.log(f"'{active.name}' is not a user label layer — select a 'User Provided Labels' layer.")
            return

        active.data = np.zeros_like(active.data)
        active.refresh()
        self.log(f"Cleared labels in '{active.name}'.")

    def _clear_all_user_labels(self):
        """Zero out every 'User Provided Labels' layer (keeps layers alive)."""
        import numpy as np
        from qtpy.QtWidgets import QMessageBox

        reply = QMessageBox.question(
            self, "Clear All Labels",
            "This will erase ALL painted labels in every user label layer.\n\nContinue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply != QMessageBox.Yes:
            return

        cleared = 0
        for layer in self.viewer.layers:
            if layer.name.startswith("User Provided Labels"):
                layer.data = np.zeros_like(layer.data)
                layer.refresh()
                cleared += 1
        self.log(f"Cleared labels in {cleared} layer(s).")

    def _on_process_all_toggled(self, checked):
        self.spin_t_start.setEnabled(not checked)
        self.spin_t_end.setEnabled(not checked)

    def _on_workers_changed(self, value):
        import os
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        if value > max_allowed:
            self.spin_workers.setValue(max_allowed)
            self.log(
                f"⚠️ Workers capped to {max_allowed} (system has {n_cores} cores). "
                f"Using all cores can freeze the system."
            )

    def _on_run_segmentation_clicked(self):
        """User-triggered batch — runs asynchronously."""
        self.run_batch_segmentation(interactive=True, block=False)

    def run_batch_segmentation(self, interactive=True, skip_existing=False,
                               block=True, extra_callbacks=None):
        """Run batch segmentation.

        ``block=True`` (default, queue path) runs synchronously.
        ``block=False`` (GUI button) executes the backend in a background
        worker with a sample-level determinate progress bar.

        When ``interactive=False``, skips overwrite + completion dialogs.

        ``extra_callbacks`` is an optional dict accepted by every public
        batch method so the processing queue can chain
        ``{"on_done": cb, "on_failed": cb, "progress": cb}`` without
        owning the underlying :class:`BackgroundOperation`.
        """
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run segmentation: No metadata loaded.")
            fire_extra_callback(extra_callbacks, "on_failed", "no metadata loaded")
            return

        if not block and self._bg.is_running():
            self.log("⚠️ A pixel-classifier segmentation run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        self.log("Starting batch segmentation...")
        self._persist_params()
        if interactive:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
            choice = prompt_overwrite_batch(
                self,
                "Existing Segmentation Results",
                ["existing segmentation results may already be present"],
                body_prefix=(
                    "Segmentation results may already exist for some timepoints."
                ),
            )
            if choice == "cancel":
                self.log("Segmentation cancelled.")
                return
            overwrite = (choice == "overwrite")
        else:
            overwrite = not skip_existing

        # Timepoint range (Qt widget reads — must stay on Qt thread).
        if self.check_process_all.isChecked():
            timepoint_range = None
            self.log("Processing ALL timepoints.")
        else:
            t_start = self.spin_t_start.value()
            t_end = self.spin_t_end.value()
            if t_start > t_end:
                self.log("Error: Start timepoint must be <= End timepoint.")
                return
            timepoint_range = (t_start, t_end)
            self.log(f"Processing timepoints: {t_start} - {t_end}")

        # Build per-category parameter dicts matching run_pixel_classifier_segmentation signature
        org_edt, org_size, org_open, org_fill = {}, {}, {}, {}
        imm_edt, imm_size, imm_open, imm_fill = {}, {}, {}, {}
        oth_edt, oth_size, oth_open, oth_fill = {}, {}, {}, {}

        for cell_type in self.all_cell_types:
            ct_key = cell_type.lower()
            if ct_key in self.param_widgets:
                w = self.param_widgets[ct_key]
                mgr = w.get('unit_mgr')
                edt = float(mgr.get_native(w['edt'])) if mgr is not None else float(w['edt'].value())
                sz  = int(round(mgr.get_native(w['min_size']))) if mgr is not None else int(w['min_size'].value())
                opn = int(w['opening'].value())
                fh  = bool(w['fill_holes'].isChecked())
            else:
                if cell_type in self.organoid_types:
                    edt, sz, opn, fh = 12.0, 1000, 3, True
                elif cell_type in self.immune_types:
                    edt, sz, opn, fh = 2.5, 10, 0, True
                else:
                    edt, sz, opn, fh = 1.0, 10, 0, True

            if cell_type in self.organoid_types:
                org_edt[cell_type] = edt; org_size[cell_type] = sz
                org_open[cell_type] = opn; org_fill[cell_type] = fh
            elif cell_type in self.immune_types:
                imm_edt[cell_type] = edt; imm_size[cell_type] = sz
                imm_open[cell_type] = opn; imm_fill[cell_type] = fh
            else:
                oth_edt[cell_type] = edt; oth_size[cell_type] = sz
                oth_open[cell_type] = opn; oth_fill[cell_type] = fh

        # Collect classifier paths
        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"

        clf_organoid_paths = {}
        clf_immune_paths = {}
        clf_other_paths = {}
        clf_death_path = None

        def get_clf_path(target):
             p = pixel_class_outdir / f'PixelClassifier_{target.capitalize()}.joblib'
             return p if p.exists() else None

        for ct in self.organoid_types:
            clf_organoid_paths[ct] = get_clf_path(ct)
        for ct in self.immune_types:
            clf_immune_paths[ct] = get_clf_path(ct)
        for ct in self.other_types:
            clf_other_paths[ct] = get_clf_path(ct)
        if self.has_death:
            clf_death_path = get_clf_path('Death')

        n_workers_val = int(self.spin_workers.value())
        log_bridge = ThreadSafeLogger(self.log)

        def _do_segmentation(progress_cb=None):
            updated_metadata = run_pixel_classifier_segmentation(
                metadata=self.metadata_loader.metadata,
                output_dir=self.metadata_loader.output_dir,
                organoid_edt_thresholds=org_edt or None,
                immune_edt_thresholds=imm_edt or None,
                other_edt_thresholds=oth_edt or None,
                organoid_segment_size_mins=org_size or None,
                immune_segment_size_mins=imm_size or None,
                other_segment_size_mins=oth_size or None,
                organoid_opening_nr_pixels=org_open or None,
                immune_opening_nr_pixels=imm_open or None,
                other_opening_nr_pixels=oth_open or None,
                organoid_fill_holes=org_fill or None,
                immune_fill_holes=imm_fill or None,
                other_fill_holes=oth_fill or None,
                only_segment=False,
                n_workers=n_workers_val,
                log_callback=log_bridge,
                overwrite_existing=overwrite,
                timepoint_range=timepoint_range,
                clf_organoid_paths=clf_organoid_paths,
                clf_immune_paths=clf_immune_paths,
                clf_other_paths=clf_other_paths,
                clf_death_path=clf_death_path,
                progress_cb=progress_cb,
            )
            return updated_metadata

        def _apply_result(updated_metadata):
            if updated_metadata is not None:
                self.metadata_loader.metadata = updated_metadata
                updated_metadata = _run_multicolor_cleanup_if_needed(self.metadata_loader, self.log, skip_unified=False)
                self.metadata_loader.metadata = updated_metadata
            self.metadata_loader.metadata = updated_metadata
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                try:
                    updated_metadata.to_csv(csv_path, index=False)
                except Exception as e:
                    self.log(f"  Warning: Could not save metadata CSV: {e}")
            self.log("Batch segmentation finished successfully!")

        if block:
            try:
                updated_metadata = _do_segmentation(progress_cb=None)
                _apply_result(updated_metadata)
                if interactive:
                    self._prompt_switch_to_viz_after_seg()
                fire_extra_callback(extra_callbacks, "on_done", updated_metadata)
            except Exception as e:
                traceback.print_exc()
                self.log(f"Error during batch segmentation: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
            return

        def _on_done(updated_metadata):
            _apply_result(updated_metadata)
            if interactive:
                self._prompt_switch_to_viz_after_seg()
            fire_extra_callback(extra_callbacks, "on_done", updated_metadata)

        def _on_failed(err: str):
            self.log(f"Error during batch segmentation: {err}")
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_segmentation,
            desc="Pixel classifier segmentation\u2026",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run_segmentation] if hasattr(self, 'btn_run_segmentation') else [],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _prompt_switch_to_viz_after_seg(self):
        """Post-run dialog: offer to jump to the Visualization tab."""
        from qtpy.QtWidgets import QMessageBox
        res = QMessageBox.question(
            self,
            "Segmentation Finished",
            "Batch segmentation finished successfully! \n\n"
            "Do you want to switch to the Visualization Tab and see the segments?",
            QMessageBox.Yes | QMessageBox.No,
        )
        if res == QMessageBox.Yes:
            parent = self.parent()
            while parent and not hasattr(parent, 'tabs'):
                parent = parent.parent()
            if parent and hasattr(parent, 'tabs'):
                parent.tabs.setCurrentIndex(1)
                if hasattr(parent, 'visualization_tab'):
                    self.log("  Loading first sample in Visualization Tab...")
                    parent.visualization_tab.sample_combo.setCurrentIndex(0)
                    parent.visualization_tab._on_load_dataset()
                    for layer in self.viewer.layers:
                        if "Segments" in layer.name:
                            layer.visible = True

    def _on_load_training_clicked(self, interactive=True):
        try:
            if self.metadata_loader.metadata is None:
                self.log("⚠️ Cannot generate training data: No metadata loaded.")
                return

            self.log("Loading training data...")
            self._persist_params()
            
            if not self.metadata_loader.metadata_loaded:
                self.log("Metadata not loaded!")
                return

            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            pixel_class_outdir.mkdir(parents=True, exist_ok=True)
            
            features_outpath = pixel_class_outdir / 'PixelClassifier_Features.zarr'
            image_outpath = pixel_class_outdir / 'PixelClassifier_Images.zarr'
            
            metadata = self.metadata_loader.metadata
            examples_per_sample = self.spin_examples.value()
            
            # --- CACHING: Check if cached data exists ---
            pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
            saved_examples = pc.get("examples_per_sample", None)
            exists = features_outpath.exists() and image_outpath.exists()

            # Invalidate cache if it was built with the old multiscale_basic_features
            # extractor (no version marker) so we never load incompatible features.
            FEATURES_VERSION = "make_features_v1"
            if exists and pc.get("features_extractor_version") != FEATURES_VERSION:
                self.log(
                    "⚠️ Cached features were built with an older extractor — "
                    "regenerating to ensure training/inference compatibility."
                )
                exists = False

            load_existing = False
            if exists:
                if interactive:
                    # Read the actual number of stored training images from the zarr header
                    # Pixel Classifier saves as (T, C, Z, Y, X) → shape[0] = image count
                    try:
                        import zarr as _zarr
                        _cached = _zarr.open(str(image_outpath), mode="r")
                        saved_image_count = _cached.shape[0]
                        del _cached
                    except Exception:
                        saved_image_count = None

                    n_samples = len(metadata['sample_name'].unique())
                    total_new_images = n_samples * examples_per_sample
                    if saved_image_count is not None:
                        msg = (
                            f"Training data already exists with {saved_image_count} images.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            f"Overwrite with new training data, or load the existing data?"
                        )
                    else:
                        msg = (
                            f"Training data already exists.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            f"Overwrite with new training data, or load the existing data?"
                        )

                    from qtpy.QtWidgets import QMessageBox
                    box = QMessageBox(self)
                    box.setWindowTitle("Training Data Detected")
                    box.setText(msg)
                    btn_generate = box.addButton("Generate New", QMessageBox.AcceptRole)
                    btn_load = box.addButton("Load Existing", QMessageBox.YesRole)
                    btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                    box.exec_()

                    if box.clickedButton() == btn_cancel:
                        self.log("Action cancelled.")
                        return
                    elif box.clickedButton() == btn_load:
                        load_existing = True
                else:
                    # Non-interactive (Queue): default to loading existing data
                    load_existing = True

            self.viewer.layers.clear() # Clear after prompting, to avoid removing if cancel

            if load_existing:
                # --- FAST PATH: load from cached zarr files ---
                self.log(f"Loading cached images and features from disk ({saved_examples} examples/sample)...")
                stacked_images = np.asarray(load_image(image_outpath))
                stacked_features = np.asarray(load_image(features_outpath))
                self.log(f"  Cached images shape: {stacked_images.shape}")
                self.log(f"  Cached features shape: {stacked_features.shape}")
            else:
                # --- SLOW PATH: regenerate images and features ---
                if exists and not load_existing:
                     self.log("Generating new training data (overwriting)...")
                elif saved_examples is not None and saved_examples != examples_per_sample:
                    self.log(f"examples_per_sample changed ({saved_examples} → {examples_per_sample}), regenerating...")
                else:
                    self.log("No cached data found, generating images and features...")
                
                all_images = []
                
                self.log(f"Checking {len(metadata)} samples...")
                
                for idx, sample in metadata.iterrows():
                    sample_name = sample['sample_name']
                    raw_image_path = Path(sample['raw_image_path'])
                    
                    if not raw_image_path.exists():
                        self.log(f"Error: Raw image path not found: {raw_image_path}")
                        continue
                    
                    raw_image_zarr = output_dir / "images" / sample_name / f"{sample_name}.zarr"
                    
                    if not raw_image_zarr.exists():
                        self.log(f"- Converting {sample_name} to .zarr...")
                        try:
                            images = load_image(raw_image_path)
                            save_as_zarr(img=images, path=raw_image_zarr)
                        except Exception as e:
                            self.log(f"Error loading/converting {sample_name}: {e}")
                            traceback.print_exc()
                            continue
                    
                    try:
                        images = load_image(raw_image_zarr)
                        if len(images.shape) == 5:
                             max_t = images.shape[0]-1
                        elif len(images.shape) == 4:
                            self.log(f"Warning: Image {sample_name} has shape {images.shape}. Assuming T,C,Z,Y,X logic might fail.")
                            max_t = images.shape[0]-1
                        else:
                            self.log(f"Warning: unexpected shape {images.shape}")
                            max_t = 0
                            
                        idc = np.linspace(0, max_t, examples_per_sample, dtype=int)
                        idc = np.unique(idc)
                        
                        self.log(f"  {sample_name}: Taking timepoints {idc}")
                        
                        sample_images = [images[t] for t in idc]
                        all_images += sample_images
                        
                    except Exception as e:
                         self.log(f"Error reading {raw_image_zarr}: {e}")
                         traceback.print_exc()

                if not all_images:
                    self.log("No images loaded.")
                    return

                max_shape = [0] * all_images[0].ndim
                for img in all_images:
                    for i in range(img.ndim):
                        max_shape[i] = max(max_shape[i], img.shape[i])
                max_shape = tuple(max_shape)
                
                self.log(f"Max shape for padding: {max_shape}")
                
                # Remove old cached files, labels, and predictions before regenerating
                self.log("Cleaning up old training files, labels and predictions...")
                to_delete = [
                    image_outpath, 
                    features_outpath,
                    pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr',
                    pixel_class_outdir / 'PixelClassifier_Death_PredictedLabels.zarr'
                ]
                for ct in self.all_cell_types:
                    to_delete.append(pixel_class_outdir / f'PixelClassifier_User{ct.capitalize()}Labels.zarr')
                    to_delete.append(pixel_class_outdir / f'PixelClassifier_{ct.capitalize()}_PredictedLabels.zarr')
                
                for p in to_delete:
                    if p.exists():
                        if p.is_dir(): shutil.rmtree(p)
                        else: p.unlink()
                
                self.log("Calculating features (this may take a while)...")

                # Build extractor from the first sample's pixel sizes so
                # sigma values are in physical units — matching batch inference.
                pixel_xy = float(metadata['pixel_distance_xy'].iloc[0])
                pixel_z  = float(metadata['pixel_distance_z'].iloc[0])
                extract_features = make_features(pixel_xy, pixel_z, verbose=True)
                self.log(
                    f"Feature extractor: make_features "
                    f"(pixel_xy={pixel_xy:.4f} µm, pixel_z={pixel_z:.4f} µm)"
                )

                processed_images = []
                processed_features = []
                
                for i, img in enumerate(all_images):
                    self.log(f"Processing image {i+1}/{len(all_images)}")
                    try:
                        feats = extract_features(img)
                        img_padded = zeropad_image_to_match_shape(img, max_shape)
                        processed_images.append(img_padded)
                        processed_features.append(feats)
                    except Exception as e:
                        self.log(f"Error calculating features: {e}")
                        traceback.print_exc()

                if not processed_images:
                     return

                stacked_images = np.stack(processed_images)
                
                feat_max_shape = [0] * processed_features[0].ndim
                for f in processed_features:
                    for i in range(f.ndim):
                         feat_max_shape[i] = max(feat_max_shape[i], f.shape[i])
                feat_max_shape = tuple(feat_max_shape)
                
                processed_features_padded = []
                for f in processed_features:
                    processed_features_padded.append(zeropad_image_to_match_shape(f, feat_max_shape))
                
                stacked_features = np.stack(processed_features_padded)
                
                save_as_zarr(stacked_images, image_outpath)
                save_as_zarr(stacked_features, features_outpath)

                # Persist version marker and pixel sizes used so cache
                # invalidation and the pre-flight check work correctly.
                pc_cfg = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
                pc_cfg["features_extractor_version"] = "make_features_v1"
                pc_cfg["pixel_size_xy"] = pixel_xy
                pc_cfg["pixel_size_z"] = pixel_z

                self.log("Images and features saved to disk.")
            
            # --- Common path: store in memory and update config ---
            self.all_images = stacked_images
            self.all_features = stacked_features
            
            # Save the current examples_per_sample to config so we can detect changes next time
            self._persist_params()
            
            # --- Add Layers to Viewer ---
            self.log("Adding layers to viewer...")
            channel_colors = ['cyan', 'yellow', 'red', 'green', 'magenta', 'blue', 'gray']
            
            n_channels = stacked_images.shape[1]
            for ch in range(n_channels):
                channel_data = stacked_images[:, ch, :, :, :]
                
                layer = self.viewer.add_image(
                    channel_data,
                    name=f"Channel {ch+1}",
                    colormap=channel_colors[ch % len(channel_colors)],
                    blending='additive',
                    visible=True,
                )
                try:
                    flat = channel_data.flatten()
                    flat = flat[flat > 0]
                    if flat.size > 0:
                        clim = np.percentile(flat, 99.8)
                        layer.contrast_limits = (0, clim)
                except:
                    pass
            
            # Create label layer dimensions
            dims = stacked_images.shape  # (T, C, Z, Y, X)
            self.log(f"  Images shape: {dims}")
            
            # label_dims should be (T, Z, Y, X)
            # If dims is 5D: (T, C, Z, Y, X)
            if len(dims) == 5:
                label_dims = (dims[0], dims[2], dims[3], dims[4])
            elif len(dims) == 4:
                # Assuming (T, C, Y, X)
                label_dims = (dims[0], 1, dims[2], dims[3])
                self.log("  Warning: 4D images detected, assuming (T, C, Y, X)")
            else:
                label_dims = dims # fallback
            
            self.log(f"  Labels shape: {label_dims}")
            
            # --- Reload existing user labels and predictions from disk ---
            # 1. User Provided Labels
            if self.has_death:
                dead_labels_path = pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr'
                if dead_labels_path.exists():
                    self.log("  Loading existing dead labels from disk...")
                    dead_labels = np.asarray(load_zarr(dead_labels_path))
                    if dead_labels.shape != label_dims:
                        self.log(f"  Warning: Dead labels shape mismatch {dead_labels.shape} vs {label_dims}. Resizing...")
                        # Simple zero-pad or crop could be here, but we'll just use zeros if wrong
                        dead_labels = np.zeros(label_dims, dtype=np.uint16)
                else:
                    dead_labels = np.zeros(label_dims, dtype=np.uint16)
                layer = self.viewer.add_labels(dead_labels, name='User Provided Labels (Dead)', opacity=0.3)
                self._configure_user_label_layer(layer)

            for cell_type in self.all_cell_types:
                user_label_path = pixel_class_outdir / f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr'
                if user_label_path.exists():
                    self.log(f"  Loading existing {cell_type} user labels from disk...")
                    user_labels = np.asarray(load_zarr(user_label_path))
                    if user_labels.shape != label_dims:
                         self.log(f"  Warning: {cell_type} user labels shape mismatch {user_labels.shape} vs {label_dims}. Using zeros.")
                         user_labels = np.zeros(label_dims, dtype=np.uint16)
                else:
                    user_labels = np.zeros(label_dims, dtype=np.uint16)
                layer = self.viewer.add_labels(
                    user_labels,
                    name=f'User Provided Labels ({cell_type.capitalize()})',
                    opacity=0.3
                )
                self._configure_user_label_layer(layer)

            # 2. Pixel Classification predictions
            if self.has_death:
                pred_death_path = pixel_class_outdir / 'PixelClassifier_Death_PredictedLabels.zarr'
                if pred_death_path.exists():
                    self.log("  Loading existing dead predictions from disk...")
                    pred_dead = np.asarray(load_zarr(pred_death_path))
                    if pred_dead.shape != label_dims:
                         self.log(f"  Warning: Dead prediction shape mismatch {pred_dead.shape} vs {label_dims}. Using zeros.")
                         pred_dead = np.zeros(label_dims, dtype=np.uint16)
                else:
                    pred_dead = np.zeros(label_dims, dtype=np.uint16)
                dead_pred_layer = self.viewer.add_labels(pred_dead, name='Pixel Classification (Dead)', opacity=0.3, visible=False)
                _set_dead_mask_layer_color(dead_pred_layer)

            for cell_type in self.all_cell_types:
                pred_label_path = pixel_class_outdir / f'PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr'
                if pred_label_path.exists():
                    self.log(f"  Loading existing {cell_type} predictions from disk...")
                    pred_data = np.asarray(load_zarr(pred_label_path))
                    if pred_data.shape != label_dims:
                        self.log(f"  Warning: {cell_type} prediction shape mismatch {pred_data.shape} vs {label_dims}. Using zeros.")
                        pred_data = np.zeros(label_dims, dtype=np.uint16)
                else:
                    pred_data = np.zeros(label_dims, dtype=np.uint16)
                self.viewer.add_labels(
                    pred_data,
                    name=f'Pixel Classification ({cell_type.capitalize()})',
                    opacity=0.3,
                    visible=False
                )

            # 3. Segments
            for cell_type in self.all_cell_types:
                 self.viewer.add_labels(
                    np.zeros(label_dims, dtype=np.uint16),
                    name=f'{cell_type.capitalize()} Segments',
                    opacity=0.7,
                    visible=False
                )

            self.log("Training data loaded successfully!")

            # Show action buttons now that data is loaded
            self.btn_save_labels.setVisible(True)
            self.btn_clear_layer.setVisible(True)
            self.btn_clear_all.setVisible(True)
            self.btn_train_update.setVisible(True)
            
            # Update and show legend
            while self.legend_container.layout().count():
                child = self.legend_container.layout().takeAt(0)
                if child.widget():
                    child.widget().deleteLater()
            legend = PerClassLegendWidget(self.viewer, self.all_cell_types, self.has_death)
            self.legend_container.layout().addWidget(legend)
            self.legend_container.setVisible(True)
            
            self.instruction_label.setVisible(True)
            self.btn_queue_train.setVisible(True)  # Show queue button now that training data is loaded
            self.is_session_active = True
            
        except Exception as e:
            traceback.print_exc()
            self.log(f"Failed to load training data: {e}")

    def run_train(self, interactive=True, extra_callbacks=None):
        """Run classifier training. Called by queue (interactive=False) or button (interactive=True).

        Training stays synchronous on the Qt thread (deeply interleaved
        with napari viewer state).  ``extra_callbacks`` is the queue's
        chaining hook fired once the synchronous training finishes.
        """
        if not self.is_session_active:
            self.log("Session not active. Loading training data...")
            self._on_load_training_clicked(interactive=interactive)

        if self.is_session_active:
            try:
                self._on_train_clicked()
                fire_extra_callback(extra_callbacks, "on_done", None)
            except Exception as e:
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
                raise
        else:
            self.log("Error: Failed to load training data for classifier.")
            fire_extra_callback(extra_callbacks, "on_failed", "no training data")
            if not interactive:
                raise RuntimeError("Training data not loaded. Please load training data manually first.")

    def _on_train_clicked(self):
        """Train classifiers, predict pixels, and segment — runs on a background thread.

        The work is split into three phases:
          1. Qt thread  — read viewer layers & widget values into plain NumPy arrays
                          via :meth:`_gather_train_inputs`.
          2. Worker     — RandomForestClassifier fit, batch predict, instance segment.
          3. Qt thread  — write result arrays back to viewer layers via
                          :meth:`_apply_train_results` (the ``on_done`` callback).

        The in-tab progress bar shows live per-target labels, e.g.
        ``"Training Dead… — 1/3 — ETA 0:22"``.
        """
        if self._bg.is_running():
            self.log("⚠️ Training is already in progress.")
            return

        # ── Phase 1 (Qt thread): gather all data we need ────────────────
        ctx = self._gather_train_inputs()
        if ctx is None:
            return  # validation failed; error already logged

        all_features       = ctx["all_features"]
        targets            = ctx["targets"]
        label_data         = ctx["label_data"]
        params_by_target   = ctx["params_by_target"]
        seg_params         = ctx["seg_params"]
        all_cell_types     = ctx["all_cell_types"]
        pixel_class_outdir = ctx["pixel_class_outdir"]

        from behav3d.napari._background_runner import ThreadSafeLogger
        safe_log = ThreadSafeLogger(self.log)

        # ── Phase 2 (worker thread) ──────────────────────────────────────
        def _do_train(progress_cb=None):
            classifiers = {}
            valid_targets = [t for t in targets if t in label_data]
            n_valid = len(valid_targets)

            # 1. Train one RF per target ─────────────────────────────────
            for i, target in enumerate(valid_targets):
                if progress_cb:
                    progress_cb(i, n_valid * 2, f"Training {target}…")
                safe_log(f"  Training {target}…")

                labels  = label_data[target]
                opening = params_by_target[target]["opening"]
                fill    = params_by_target[target]["fill"]

                # Postprocess user labels
                fg_mask = (labels == 2).astype(bool)
                if np.any(fg_mask):
                    processed_fg = postprocess_mask(
                        fg_mask, fill_holes=fill, opening_nr_pixels=opening
                    )
                    labels = np.where(
                        processed_fg, 2, np.where(labels == 1, 1, 0)
                    ).astype(labels.dtype)

                # Save postprocessed labels to zarr
                labels_outpath = (
                    pixel_class_outdir / f"PixelClassifier_User{target}Labels.zarr"
                )
                if labels_outpath.exists():
                    shutil.rmtree(labels_outpath)
                save_as_zarr(labels, labels_outpath)
                safe_log(f"  Saved {target} user labels (postprocessed)")

                flat_labels   = labels.ravel()
                flat_features = all_features.reshape(-1, all_features.shape[-1])
                label_indices = np.flatnonzero(flat_labels > 0)
                if label_indices.size == 0:
                    safe_log(f"  No labeled pixels for {target}. Skipping.")
                    continue
                selected_features = flat_features[label_indices]
                selected_labels   = flat_labels[label_indices]

                nr_bg = int(np.sum(selected_labels == 1))
                nr_fg = int(np.sum(selected_labels == 2))
                total = nr_bg + nr_fg
                safe_log(f"    BG pixels: {nr_bg}, FG pixels: {nr_fg}")
                if total == 0:
                    safe_log(f"  No labeled pixels for {target}. Skipping.")
                    continue

                class_weights = {1: nr_bg / total, 2: nr_fg / total}
                clf = RandomForestClassifier(
                    n_estimators=50, n_jobs=-1, max_depth=20,
                    class_weight=class_weights,
                )
                from skimage import future
                clf = future.fit_segmenter(selected_labels, selected_features, clf)

                clf_path = pixel_class_outdir / f"PixelClassifier_{target}.joblib"
                joblib.dump(clf, clf_path)
                safe_log(f"  Saved classifier: {clf_path.name}")
                classifiers[target] = clf

            # 2. Predict all pixels ──────────────────────────────────────
            pred_masks   = {}
            n_timepoints = all_features.shape[0]

            for j, (target, clf) in enumerate(classifiers.items()):
                if progress_cb:
                    progress_cb(
                        n_valid + j, n_valid * 2,
                        f"Predicting {target} pixels…",
                    )
                safe_log(f"  Predicting {target} pixels…")
                prediction_stack = []
                for t in range(n_timepoints):
                    feat_t        = all_features[t]
                    spatial_shape = feat_t.shape[:-1]
                    X_t           = feat_t.reshape(-1, feat_t.shape[-1])
                    batch_size    = 100_000
                    y_pred_t      = np.zeros(X_t.shape[0], dtype=np.uint8)
                    for start in range(0, X_t.shape[0], batch_size):
                        y_pred_t[start : start + batch_size] = clf.predict(
                            X_t[start : start + batch_size]
                        )
                    prediction_stack.append(y_pred_t.reshape(spatial_shape))
                full_prediction = np.stack(prediction_stack)

                opening = params_by_target.get(target, {}).get("opening", 0)
                fill    = params_by_target.get(target, {}).get("fill", True)
                fg_mask = (full_prediction == 2).astype(bool)
                processed_fg = postprocess_mask(
                    fg_mask, fill_holes=fill, opening_nr_pixels=opening
                )
                full_prediction = np.where(
                    processed_fg, 2, np.where(full_prediction == 1, 1, 0)
                ).astype(np.uint8)
                pred_masks[target] = full_prediction

                pred_outpath = (
                    pixel_class_outdir
                    / f"PixelClassifier_{target}_PredictedLabels.zarr"
                )
                if pred_outpath.exists():
                    shutil.rmtree(pred_outpath)
                save_as_zarr(full_prediction, pred_outpath)

            # 3. Segment instances ───────────────────────────────────────
            seg_results = {}
            safe_log("Segmenting cell instances…")
            for cell_type in all_cell_types:
                target = cell_type.capitalize()
                if target not in pred_masks:
                    continue
                sp               = seg_params.get(cell_type, {})
                edt_threshold    = sp.get("edt_threshold", 1.0)
                segment_size_min = sp.get("segment_size_min", 10)
                safe_log(
                    f"  Segmenting {target} (EDT={edt_threshold}, "
                    f"min_size={segment_size_min})…"
                )
                pred_mask = pred_masks[target]
                segmented_timepoints = []
                for t_idx in range(pred_mask.shape[0]):
                    mask_t = (pred_mask[t_idx] == 2)
                    if int(mask_t.sum()) == 0:
                        segmented_timepoints.append(
                            np.zeros_like(mask_t, dtype=np.uint16)
                        )
                    else:
                        segmented_timepoints.append(
                            segment_mask(
                                mask=mask_t,
                                edt_thr=edt_threshold,
                                edt_thr_refined=None,
                                segment_size_min=segment_size_min,
                                use_dims=2,
                                n_workers=1,
                            ).astype(np.uint16, copy=False)
                        )
                full_seg = np.stack(segmented_timepoints, axis=0)
                safe_log(f"  {target}: max label = {int(full_seg.max())}")
                seg_results[cell_type] = full_seg

            return {"predictions": pred_masks, "segments": seg_results}

        # ── Phase 3 (Qt thread): apply results to viewer ─────────────────
        def _on_done(results):
            self._apply_train_results(results)

        def _on_failed(err):
            self.log(f"❌ Error during training/segmentation: {err}")

        self._bg.run(
            fn=_do_train,
            desc="Training pixel classifier…",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_train_update],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
            inject_progress=True,
            indeterminate=False,
        )

    # ── Training helpers ─────────────────────────────────────────────────

    def _gather_train_inputs(self):
        """Phase 1 (Qt thread): read viewer layers & widget values.

        Persists params and labels, loads features, reads label layer data
        and all widget values into plain Python/NumPy objects so the worker
        thread never touches Qt widgets or napari layers directly.

        Returns a context dict on success, or ``None`` when validation fails.
        """
        self._persist_params()
        self._save_user_labels()

        # Load features
        if self.all_features is None:
            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            features_outpath = pixel_class_outdir / "PixelClassifier_Features.zarr"
            if features_outpath.exists():
                self.log("Loading features from disk…")
                self.all_features = np.asarray(load_zarr(features_outpath))
            else:
                self.log("Error: No features found. Please load training data first.")
                return None

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        all_features = np.asarray(self.all_features)
        self.log(f"Features shape: {all_features.shape}")

        targets = []
        if self.has_death:
            targets.append("Dead")
        targets.extend([t.capitalize() for t in self.all_cell_types])

        # Read label layers — must be on the Qt thread
        label_data = {}
        for target in targets:
            layer_name = f"User Provided Labels ({target})"
            if layer_name not in self.viewer.layers:
                self.log(f"  Layer {layer_name} not found. Skipping.")
                continue
            labels = self.viewer.layers[layer_name].data.copy()
            if labels.max() == 0:
                self.log(f"  No labels for {target}. Skipping.")
                continue
            label_data[target] = labels

        # Read postprocessing params from widgets — must be on the Qt thread
        params_by_target = {}
        for target in targets:
            cell_type_key = target.lower()
            if cell_type_key in self.param_widgets:
                w       = self.param_widgets[cell_type_key]
                opening = int(w["opening"].value())
                fill    = bool(w["fill_holes"].isChecked())
            elif cell_type_key == "dead":
                opening, fill = 0, True
            else:
                opening, fill = 0, True
            params_by_target[target] = {"opening": opening, "fill": fill}

        # Read segmentation params from widgets — must be on the Qt thread
        seg_params = {}
        for cell_type in self.all_cell_types:
            cell_type_key = cell_type.lower()
            if cell_type_key in self.param_widgets:
                w                = self.param_widgets[cell_type_key]
                mgr              = w.get('unit_mgr')
                edt_threshold    = float(mgr.get_native(w["edt"])) if mgr is not None else float(w["edt"].value())
                segment_size_min = int(round(mgr.get_native(w["min_size"]))) if mgr is not None else int(w["min_size"].value())
            else:
                if cell_type in self.organoid_types:
                    edt_threshold, segment_size_min = 12.0, 1000
                elif cell_type in self.immune_types:
                    edt_threshold, segment_size_min = 2.5, 10
                else:
                    edt_threshold, segment_size_min = 1.0, 10
            seg_params[cell_type] = {
                "edt_threshold": edt_threshold,
                "segment_size_min": segment_size_min,
            }

        return {
            "all_features": all_features,
            "targets": targets,
            "label_data": label_data,
            "params_by_target": params_by_target,
            "seg_params": seg_params,
            "all_cell_types": list(self.all_cell_types),
            "pixel_class_outdir": pixel_class_outdir,
        }

    def _apply_train_results(self, results):
        """Phase 3 (Qt thread): write training results back to viewer layers."""
        if results is None:
            return
        pred_masks  = results.get("predictions", {})
        seg_results = results.get("segments", {})

        for target, full_prediction in pred_masks.items():
            pred_layer_name = f"Pixel Classification ({target})"
            if pred_layer_name in self.viewer.layers:
                pred_layer = self.viewer.layers[pred_layer_name]
                pred_layer.data = full_prediction
                pred_layer.refresh()
            else:
                pred_layer = self.viewer.add_labels(
                    full_prediction, name=pred_layer_name, opacity=0.3
                )
            if target == "Dead":
                _set_dead_mask_layer_color(pred_layer)
            self.log(f"  {target} prediction updated.")

        for cell_type, full_seg in seg_results.items():
            target         = cell_type.capitalize()
            seg_layer_name = f"{target} Segments"
            if seg_layer_name in self.viewer.layers:
                self.viewer.layers[seg_layer_name].data    = full_seg
                self.viewer.layers[seg_layer_name].visible = True
                self.viewer.layers[seg_layer_name].refresh()
            else:
                self.viewer.add_labels(full_seg, name=seg_layer_name, opacity=0.7)

        self.log("Training, prediction, and segmentation complete!")


    def _on_resegment_clicked(self):
        self.log("Testing segmentation parameters on current view...")
        try:
            # 1. Identify active timepoint
            current_step = self.viewer.dims.current_step
            # Assuming (T, Z, Y, X) or (T, C, Z, Y, X)
            # We need to find T index. 
            # If dims order is standard BEHAV3D: T is 0
            t_idx = current_step[0]
            self.log(f"Current timepoint index: {t_idx}")

            # 2. Iterate over cell types
            targets = []
            if self.has_death:
                targets.append('Dead')
            targets.extend([t.capitalize() for t in self.all_cell_types])

            for target in targets:
                # Get Pixel Classification layer (Prediction)
                pred_layer_name = f'Pixel Classification ({target})'
                if pred_layer_name not in self.viewer.layers:
                    self.log(f"  Warning: Layer '{pred_layer_name}' not found. skipping {target}.")
                    continue
                
                pred_layer = self.viewer.layers[pred_layer_name]
                pred_data = pred_layer.data # (T, Z, Y, X)
                
                if t_idx >= pred_data.shape[0]:
                    self.log(f"  Error: Current timepoint {t_idx} out of range for {target}.")
                    continue

                # Get prediction for current timepoint
                # Expected shape (Z, Y, X)
                mask_t = pred_data[t_idx]
                
                # Check for content
                # Background=1, Foreground=2 (usually)
                # But let's check values
                
                # Get parameters for this target
                # If 'Dead', use defaults or specific widget if added (currently Dead doesn't have params widget in UI logic above, need check)
                # The UI logic in _refresh_params_ui iterates over self.all_cell_types only.
                # So 'Dead' might not have adjustable params in UI. 
                # Let's skip 'Dead' if not in param_widgets or use defaults.
                
                cell_type_key = target.lower()
                if cell_type_key == 'dead':
                     # Use defaults for dead
                     edt, size, opening, fill = 2.5, 10, 0, True
                else:
                    if cell_type_key in self.param_widgets:
                        widgets = self.param_widgets[cell_type_key]
                        mgr = widgets.get('unit_mgr')
                        edt = float(mgr.get_native(widgets['edt'])) if mgr is not None else widgets['edt'].value()
                        size = int(round(mgr.get_native(widgets['min_size']))) if mgr is not None else widgets['min_size'].value()
                        opening = widgets['opening'].value()
                        fill = widgets['fill_holes'].isChecked()
                    else:
                        self.log(f"  Warning: No parameters found for {target}. Using defaults.")
                        edt, size, opening, fill = 2.5, 10, 0, True

                self.log(f"  Segmenting {target} (EDT={edt}, Size={size}, Open={opening})...")
                
                # Preprocess mask (Postprocess prediction)
                # logic from napari_pixelclassifier.segment_and_update:
                # mask_t has labels: 0=background, 1=background_label, 2=foreground
                # extract fg (2)
                
                # Postprocess the prediction mask
                fg_mask = (mask_t == 2).astype(bool)
                
                if np.any(fg_mask):
                    fg_mask = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                    
                # Now Segment
                if np.sum(fg_mask) == 0:
                    segmented = np.zeros_like(fg_mask, dtype=np.uint16)
                    self.log(f"  {target}: Empty mask.")
                else:
                    segmented = segment_mask(
                        fg_mask, 
                        edt_thr=edt, 
                        segment_size_min=size
                    )
                    self.log(f"  {target}: {segmented.max()} objects found.")

                # Update 'Segments' layer
                seg_layer_name = f'{target} Segments' # e.g. "Organoid1 Segments"
                # "Dead" usually doesn't have a segments layer in the original code snippet 
                # (only 'User Provided' and 'Pixel Classification' for Dead).
                # But let's check create_layers logic. 
                # In _on_load_training_clicked:
                # self.viewer.add_labels(..., name=f'{cell_type.capitalize()} Segments', ...) for cell types
                # Not for Dead.
                
                if target == 'Dead':
                     continue 
                
                if seg_layer_name in self.viewer.layers:
                    layer = self.viewer.layers[seg_layer_name]
                    self.viewer.layers[seg_layer_name].visible = True

                    # Check shape
                    if layer.data.shape != pred_data.shape:
                        self.log(f"  Warning: Shape mismatch for {seg_layer_name}. Re-initializing layer data.")
                        layer.data = np.zeros_like(pred_data, dtype=np.uint16)
                    
                    # Update only current timepoint
                    # We need to write to the layer.data
                    # Accessing layer.data returns numpy-like (or tensorstore/dask)
                    # If it's in-memory numpy:
                    try:
                        layer.data[t_idx] = segmented
                        layer.refresh()
                    except Exception as e:
                         self.log(f"  Error updating layer {seg_layer_name}: {e}")
                else:
                    self.log(f"  Layer {seg_layer_name} not found. Creating...")
                    # Create full empty array
                    new_data = np.zeros_like(pred_data, dtype=np.uint16)
                    new_data[t_idx] = segmented
                    self.viewer.add_labels(new_data, name=seg_layer_name, opacity=0.7)

            self.log("Segmentation preview updated.")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error resegmenting: {e}")




class CellposeWidget(QWidget):
    """Napari widget for Cellpose channel configuration, model loading, and inference."""

    def __init__(self, viewer, metadata_loader, log_callback=None,
                 tab_progress_row=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.pretrained_model_dir = None
        # Background-execution infrastructure (shared progress row comes
        # from the parent SegmentationTab).
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        self._detect_cell_types()
        _cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {}) or {}
        self._manual_n_channels = _cellpose_cfg.get("manual_n_channels")
        self._init_ui()

    # ── cell-type detection ─────────────────────────────────────────────
    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
            is_combined_multicolor_celltype,
        )
        metadata = self.metadata_loader.metadata
        if metadata is None:
            self.organoid_types = []
            self.immune_types = []
            self.other_types = []
            self.has_death = False
            self.all_cell_types = []
            return
        self.organoid_types = [ct for ct in detect_organoid_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.immune_types = [ct for ct in detect_immune_cell_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.other_types = [ct for ct in detect_other_cell_types_from_metadata(metadata) if not is_combined_multicolor_celltype(ct)]
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _calculate_number_of_channels(self):
        n = len(self.organoid_types) + len(self.immune_types) + len(self.other_types)
        if self.has_death:
            n += 1
        return n

    def _channel_detail_text(self):
        parts = [
            f"Organoids: {len(self.organoid_types)}",
            f"Immune: {len(self.immune_types)}",
            f"Other: {len(self.other_types)}",
        ]
        if self.has_death:
            parts.append("Dead: 1")
        if self._manual_n_channels is not None:
            parts.append(f"Manual total: {self._manual_n_channels}")
        return "  ·  ".join(parts)

    def _get_channel_options(self):
        opts = list(self.organoid_types) + list(self.immune_types) + list(self.other_types)
        if self.has_death:
            opts.append("dead")
        # Add a single "none" option when manual override adds extra channels
        if self._manual_n_channels is not None and self._manual_n_channels > len(opts):
            opts.append("none")
        return opts

    def _get_sample_names(self):
        if self.metadata_loader.metadata is None:
            return []
        return list(self.metadata_loader.metadata["sample_name"].unique())

    # ── helpers ─────────────────────────────────────────────────────────
    def _get_prefix(self, cell_type):
        if cell_type in self.organoid_types:
            return "or"
        if cell_type in self.immune_types:
            return "im"
        return "ot"

    # ── UI ──────────────────────────────────────────────────────────────
    def _init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        content = QWidget()
        content_layout = QVBoxLayout()
        content_layout.setContentsMargins(4, 4, 4, 4)
        content_layout.setSpacing(4)
        content.setLayout(content_layout)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        # ── Section 1: Channel Configuration ─────────────────────────
        self.channel_group = QGroupBox("1. Channel Configuration")
        ch_layout = QVBoxLayout()
        self.channel_group.setLayout(ch_layout)

        self.n_channels = self._calculate_number_of_channels()
        # Primary: compact channel count
        self.n_channels_label = QLabel(f"<b>Detected:</b> {self.n_channels} channel(s)")
        self.n_channels_label.setWordWrap(True)
        ch_layout.addWidget(self.n_channels_label)

        # Secondary: per-type breakdown
        self.n_channels_detail_label = QLabel(self._channel_detail_text())
        self.n_channels_detail_label.setWordWrap(True)
        self.n_channels_detail_label.setStyleSheet("color: #aaa; font-size: 11px;")
        ch_layout.addWidget(self.n_channels_detail_label)

        # Manual override row
        override_row = QWidget()
        override_hl = QHBoxLayout()
        override_hl.setContentsMargins(0, 2, 0, 2)
        override_row.setLayout(override_hl)
        self.btn_override_channels = QPushButton("⚙ Set Manually")
        self.btn_override_channels.setToolTip("Manually set a higher channel count for extra unlabelled channels")
        self.btn_override_channels.setFixedWidth(110)
        self.btn_override_channels.clicked.connect(self._on_show_override)
        override_hl.addWidget(self.btn_override_channels)
        self.spin_manual_channels = QSpinBox()
        self.spin_manual_channels.setRange(1, 99)
        self.spin_manual_channels.setValue(max(1, self.n_channels))
        self._lbl_manual_ch = QLabel("Channels:")
        override_hl.addWidget(self._lbl_manual_ch)
        override_hl.addWidget(self.spin_manual_channels)
        self.btn_apply_override = QPushButton("Apply")
        self.btn_apply_override.setFixedWidth(55)
        self.btn_apply_override.clicked.connect(self._on_apply_override)
        override_hl.addWidget(self.btn_apply_override)
        self.btn_clear_override = QPushButton("↺")
        self.btn_clear_override.setFixedWidth(32)
        self.btn_clear_override.setToolTip("Clear manual override, revert to auto-detected count")
        self.btn_clear_override.clicked.connect(self._on_clear_override)
        override_hl.addWidget(self.btn_clear_override)
        override_hl.addStretch()
        ch_layout.addWidget(override_row)
        # Hide spinbox controls initially; show if a saved override exists
        _override_active = self._manual_n_channels is not None
        self.btn_override_channels.setVisible(not _override_active)
        self.spin_manual_channels.setVisible(_override_active)
        self._lbl_manual_ch.setVisible(_override_active)
        self.btn_apply_override.setVisible(_override_active)
        self.btn_clear_override.setVisible(_override_active)
        if _override_active:
            self.spin_manual_channels.setValue(self._manual_n_channels)

        # Mode selection
        from qtpy.QtWidgets import QRadioButton, QButtonGroup
        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Mode:"))
        self.radio_same = QRadioButton("Same for all")
        self.radio_per_sample = QRadioButton("Per sample")
        self.mode_group = QButtonGroup(self)
        self.mode_group.addButton(self.radio_same, 0)
        self.mode_group.addButton(self.radio_per_sample, 1)
        mode_row.addWidget(self.radio_same)
        mode_row.addWidget(self.radio_per_sample)
        ch_layout.addLayout(mode_row)

        # Container for dynamic dropdowns
        self.channel_container = QWidget()
        self.channel_container_layout = QVBoxLayout()
        self.channel_container_layout.setContentsMargins(0, 0, 0, 0)
        self.channel_container.setLayout(self.channel_container_layout)
        ch_layout.addWidget(self.channel_container)

        # Save button
        self.btn_save_channels = QPushButton("💾  Save Channel Configuration")
        self.btn_save_channels.setStyleSheet(
            "QPushButton { background: #28a745; color: white; border-radius: 4px; "
            "padding: 4px 8px; font-weight: bold; }"
            "QPushButton:hover { background: #218838; }"
        )
        self.btn_save_channels.clicked.connect(self._on_save_channels)
        ch_layout.addWidget(self.btn_save_channels)

        content_layout.addWidget(self.channel_group)

        # Load saved mode and build dropdowns
        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        saved_mode = cellpose_cfg.get("labels_mode", "same_for_all") if cellpose_cfg else "same_for_all"
        if saved_mode == "per_sample":
            self.radio_per_sample.setChecked(True)
        else:
            self.radio_same.setChecked(True)

        self.channel_dropdowns = {}
        self._build_channel_dropdowns()
        self.mode_group.buttonClicked.connect(lambda _btn: self._build_channel_dropdowns())

        # ── Section 2: Model & Cell Type ─────────────────────────────
        model_group = QGroupBox("2. Model & Cell Type")
        model_layout = QVBoxLayout()
        model_group.setLayout(model_layout)

        # Model path
        path_row = QHBoxLayout()
        path_row.addWidget(QLabel("Model:"))
        self.model_path_edit = QLineEdit()
        self.model_path_edit.setPlaceholderText("Path to Cellpose model file…")
        path_row.addWidget(self.model_path_edit, stretch=1)
        self.btn_browse = QPushButton("Browse…")
        self.btn_browse.clicked.connect(self._on_browse)
        path_row.addWidget(self.btn_browse)
        model_layout.addLayout(path_row)

        # Load model button
        self.btn_load_model = QPushButton("📂  Load Model")
        self.btn_load_model.setStyleSheet(
            "QPushButton { background: #17a2b8; color: white; border-radius: 4px; "
            "padding: 4px 8px; font-weight: bold; }"
            "QPushButton:hover { background: #138496; }"
        )
        self.btn_load_model.clicked.connect(self._on_load_model)
        model_layout.addWidget(self.btn_load_model)

        self.model_info_label = QLabel("")
        self.model_info_label.setStyleSheet("color: #aaa; font-size: 11px;")
        model_layout.addWidget(self.model_info_label)

        # Cell type selector
        ct_row = QHBoxLayout()
        ct_row.addWidget(QLabel("Cell type to segment:"))
        self.cell_type_combo = QComboBox()
        if self.all_cell_types:
            self.cell_type_combo.addItems(self.all_cell_types)
        else:
            self.cell_type_combo.addItem("(load metadata)")
        ct_row.addWidget(self.cell_type_combo, stretch=1)
        model_layout.addLayout(ct_row)

        # Timepoint range — checkbox on its own row, spinboxes on a second row
        tp_check_row = QHBoxLayout()
        self.check_process_all = QCheckBox("Process all timepoints")
        self.check_process_all.setChecked(True)
        tp_check_row.addWidget(self.check_process_all)
        tp_check_row.addStretch()
        model_layout.addLayout(tp_check_row)

        tp_range_row = QHBoxLayout()
        tp_range_row.addWidget(QLabel("  Start:"))
        self.spin_t_start = QSpinBox()
        self.spin_t_start.setRange(0, 9999)
        self.spin_t_start.setValue(0)
        self.spin_t_start.setMaximumWidth(70)
        tp_range_row.addWidget(self.spin_t_start)
        tp_range_row.addWidget(QLabel("End:"))
        self.spin_t_end = QSpinBox()
        self.spin_t_end.setRange(0, 9999)
        self.spin_t_end.setValue(0)
        self.spin_t_end.setMaximumWidth(70)
        tp_range_row.addWidget(self.spin_t_end)
        tp_range_row.addStretch()
        model_layout.addLayout(tp_range_row)

        self.check_process_all.toggled.connect(lambda on: (
            self.spin_t_start.setEnabled(not on),
            self.spin_t_end.setEnabled(not on),
        ))
        self.spin_t_start.setEnabled(False)
        self.spin_t_end.setEnabled(False)

        content_layout.addWidget(model_group)

        # ── Section 3: Run & Queue ───────────────────────────────────
        run_group = QGroupBox("3. Run")
        run_layout = QVBoxLayout()
        run_group.setLayout(run_layout)

        # Run Cellpose button (text updates with cell type)
        run_btn_row = QHBoxLayout()
        run_btn_row.setSpacing(4)
        self.btn_run_cellpose = QPushButton("▶  Run Cellpose")
        self.btn_run_cellpose.setStyleSheet(
            "QPushButton { background: #28a745; color: white; border-radius: 4px; "
            "padding: 6px; font-weight: bold; font-size: 12px; }"
            "QPushButton:hover { background: #218838; }"
        )
        self.btn_run_cellpose.clicked.connect(lambda: self.run_batch_cellpose(interactive=True, block=False))
        run_btn_row.addWidget(self.btn_run_cellpose, stretch=1)

        # +🛒 queue button
        self.btn_queue_cellpose = QPushButton("+🛒")
        self.btn_queue_cellpose.setFixedSize(36, 28)
        self.btn_queue_cellpose.setToolTip("Add Cellpose Segmentation to Processing Queue")
        self.btn_queue_cellpose.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )
        self.btn_queue_cellpose.setVisible(False)  # shown after metadata loaded
        run_btn_row.addWidget(self.btn_queue_cellpose)
        run_layout.addLayout(run_btn_row)

        # Dead mask (Otsu) button + queue button
        otsu_row = QHBoxLayout()
        otsu_row.setSpacing(4)
        self.btn_run_otsu = QPushButton("☠  Run Dead Mask (Otsu)")
        self.btn_run_otsu.setStyleSheet(
            "QPushButton { background: #6c757d; color: white; border-radius: 4px; "
            "padding: 4px; font-weight: bold; }"
            "QPushButton:hover { background: #5a6268; }"
        )
        self.btn_run_otsu.clicked.connect(lambda: self.run_otsu_threshold(interactive=True, block=False))
        self.btn_run_otsu.setVisible(self.has_death)
        otsu_row.addWidget(self.btn_run_otsu, stretch=1)
        self.btn_queue_otsu = QPushButton("+🛒")
        self.btn_queue_otsu.setFixedSize(36, 28)
        self.btn_queue_otsu.setToolTip("Add Dead Mask (Otsu) to Processing Queue")
        self.btn_queue_otsu.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )
        self.btn_queue_otsu.setVisible(False)  # shown after metadata loaded
        otsu_row.addWidget(self.btn_queue_otsu)
        run_layout.addLayout(otsu_row)

        content_layout.addWidget(run_group)
        content_layout.addStretch()

        # Update run button label when cell type changes
        self.cell_type_combo.currentTextChanged.connect(self._update_run_button_label)

    # ── dynamic run button label ────────────────────────────────────────
    def _update_run_button_label(self, text):
        if text and text != "(load metadata)":
            self.btn_run_cellpose.setText(f"▶  Run Cellpose — {text}")
        else:
            self.btn_run_cellpose.setText("▶  Run Cellpose")

    # ── channel dropdown builder ────────────────────────────────────────
    def _build_channel_dropdowns(self):
        """Rebuild the channel label dropdowns based on current mode."""
        # Clear existing
        while self.channel_container_layout.count():
            item = self.channel_container_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
            elif item.layout():
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()

        self.channel_dropdowns = {}
        n_channels = self._manual_n_channels if self._manual_n_channels is not None else self._calculate_number_of_channels()
        channel_options = self._get_channel_options()
        auto_n = self._calculate_number_of_channels()  # channels with a real label
        if n_channels == 0:
            self.channel_container_layout.addWidget(
                QLabel('<span style="color: orange;">No channels detected. Load metadata first.</span>')
            )
            return

        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {}) or {}
        is_per_sample = self.radio_per_sample.isChecked()

        if not is_per_sample:
            # ── Same for all ──
            saved_labels = cellpose_cfg.get("channel_labels", {})
            self.channel_dropdowns["global"] = {}
            for i in range(n_channels):
                row = QHBoxLayout()
                row.addWidget(QLabel(f"Channel {i}:"))
                combo = QComboBox()
                combo.addItems(channel_options)
                saved_val = saved_labels.get(str(i), saved_labels.get(i))
                if saved_val and saved_val in channel_options:
                    combo.setCurrentText(saved_val)
                elif i < auto_n:
                    combo.setCurrentIndex(i)  # auto channel — assign by position
                elif "none" in channel_options:
                    combo.setCurrentText("none")  # extra channel — default to none
                row.addWidget(combo, stretch=1)
                self.channel_dropdowns["global"][i] = combo
                self.channel_container_layout.addLayout(row)
        else:
            # ── Per sample — same vertical layout as "same for all"
            saved_per_sample = cellpose_cfg.get("per_sample_channel_labels", {})
            for sample_name in self._get_sample_names():
                header = QLabel(f"<b>{sample_name}</b>")
                self.channel_container_layout.addWidget(header)
                self.channel_dropdowns[sample_name] = {}
                sample_saved = saved_per_sample.get(sample_name, {})
                for i in range(n_channels):
                    row = QHBoxLayout()
                    row.addWidget(QLabel(f"Channel {i}:"))
                    combo = QComboBox()
                    combo.addItems(channel_options)
                    saved_val = sample_saved.get(str(i), sample_saved.get(i))
                    if saved_val and saved_val in channel_options:
                        combo.setCurrentText(saved_val)
                    elif i < auto_n:
                        combo.setCurrentIndex(i)  # auto channel — assign by position
                    elif "none" in channel_options:
                        combo.setCurrentText("none")  # extra channel — default to none
                    row.addWidget(combo, stretch=1)
                    self.channel_dropdowns[sample_name][i] = combo
                    self.channel_container_layout.addLayout(row)

    # ── save / persist ──────────────────────────────────────────────────
    def _persist_channel_config(self):
        """Save channel configuration to behav3d_parameters.yml."""
        cellpose_cfg = self.metadata_loader.behav3d_parameters.setdefault("cellpose", {})
        cellpose_cfg["number_of_channels"] = self._manual_n_channels if self._manual_n_channels is not None else self._calculate_number_of_channels()
        if self._manual_n_channels is not None:
            cellpose_cfg["manual_n_channels"] = self._manual_n_channels
        else:
            cellpose_cfg.pop("manual_n_channels", None)

        is_per_sample = self.radio_per_sample.isChecked()
        cellpose_cfg["labels_mode"] = "per_sample" if is_per_sample else "same_for_all"

        if not is_per_sample:
            channel_labels = {}
            if "global" in self.channel_dropdowns:
                for i, combo in self.channel_dropdowns["global"].items():
                    channel_labels[i] = combo.currentText()
            cellpose_cfg["channel_labels"] = channel_labels
            cellpose_cfg["per_sample_channel_labels"] = {}
        else:
            per_sample_labels = {}
            for sample_name, combos in self.channel_dropdowns.items():
                if sample_name != "global":
                    per_sample_labels[sample_name] = {i: c.currentText() for i, c in combos.items()}
            cellpose_cfg["per_sample_channel_labels"] = per_sample_labels
            cellpose_cfg["channel_labels"] = {}

        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False,
            )

    def _on_save_channels(self):
        self._persist_channel_config()
        self.log("✅ Cellpose channel configuration saved.")

    # ── override channel count ──────────────────────────────────────────
    def _on_show_override(self):
        calc = self._calculate_number_of_channels()
        self.spin_manual_channels.setRange(max(1, calc), 99)
        self.spin_manual_channels.setValue(self._manual_n_channels or max(1, calc))
        self.btn_override_channels.setVisible(False)
        self._lbl_manual_ch.setVisible(True)
        self.spin_manual_channels.setVisible(True)
        self.btn_apply_override.setVisible(True)
        self.btn_clear_override.setVisible(True)

    def _on_apply_override(self):
        self._manual_n_channels = self.spin_manual_channels.value()
        self.n_channels_detail_label.setText(self._channel_detail_text())
        self._build_channel_dropdowns()
        self.log(f"✅ Channel count manually set to {self._manual_n_channels}.")

    def _on_clear_override(self):
        self._manual_n_channels = None
        self._lbl_manual_ch.setVisible(False)
        self.spin_manual_channels.setVisible(False)
        self.btn_apply_override.setVisible(False)
        self.btn_clear_override.setVisible(False)
        self.btn_override_channels.setVisible(True)
        self.n_channels_detail_label.setText(self._channel_detail_text())
        self._build_channel_dropdowns()
        self.log("↺ Channel count reset to auto-detected.")

    # ── browse / load model ─────────────────────────────────────────────
    def _on_browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Cellpose Model File")
        if path:
            self.model_path_edit.setText(path)

    def _on_load_model(self):
        path = self.model_path_edit.text().strip()
        if not path:
            self.log("⚠️ Please enter or browse for a model path first.")
            return
        self.pretrained_model_dir = path
        # Parse channel order from model dirname
        try:
            channel_order = path.replace("\\", "/").split("/")[-1].split("__")[-1].split("_")
            channel_order = [ch.split("-")[-1] for ch in channel_order]
        except Exception:
            channel_order = ["Unknown"]
        self.model_info_label.setText(f"Model loaded  ·  Channels: {', '.join(channel_order)}")
        self.log(f"✅ Loaded Cellpose model from: {path}")
        self.log(f"   Channels detected: {channel_order}")

    # ── metadata update ─────────────────────────────────────────────────
    def _on_metadata_updated(self, _metadata=None):
        self._detect_cell_types()
        self.n_channels = self._calculate_number_of_channels()
        self.n_channels_label.setText(f"<b>Detected:</b> {self.n_channels} channel(s)")
        self.n_channels_detail_label.setText(self._channel_detail_text())
        self._build_channel_dropdowns()
        # Refresh cell type combo
        current = self.cell_type_combo.currentText()
        self.cell_type_combo.blockSignals(True)
        self.cell_type_combo.clear()
        if self.all_cell_types:
            self.cell_type_combo.addItems(self.all_cell_types)
            if current in self.all_cell_types:
                self.cell_type_combo.setCurrentText(current)
        else:
            self.cell_type_combo.addItem("(load metadata)")
        self.cell_type_combo.blockSignals(False)
        self._update_run_button_label(self.cell_type_combo.currentText())
        # Show queue & otsu buttons
        self.btn_queue_cellpose.setVisible(True)
        self.btn_run_otsu.setVisible(self.has_death)
        self.btn_queue_otsu.setVisible(self.has_death)

    # ── queue validation ────────────────────────────────────────────────
    def validate_for_queue(self):
        """Return (ok, error_msg) — checks channel config + model + cell type."""
        # 1. Channel config saved?
        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        if not cellpose_cfg or not cellpose_cfg.get("channel_labels") and not cellpose_cfg.get("per_sample_channel_labels"):
            return False, "Channel configuration not saved. Please configure and save channel labels first (Section 1)."
        # 2. Model loaded?
        if not self.pretrained_model_dir:
            return False, "No Cellpose model loaded. Please load a model first (Section 2)."
        # 3. Cell type selected?
        ct = self.cell_type_combo.currentText()
        if not ct or ct == "(load metadata)":
            return False, "No cell type selected. Please select a cell type to segment (Section 2)."
        return True, ""

    def get_queue_params(self):
        """Return the params dict for a queue step."""
        return {
            "cell_type": self.cell_type_combo.currentText(),
            "model_path": self.pretrained_model_dir,
        }

    # ── run cellpose ────────────────────────────────────────────────────
    def run_batch_cellpose(self, interactive=True, skip_existing=False,
                           cell_type_override=None, model_path_override=None,
                           block=True, extra_callbacks=None):
        """Run Cellpose segmentation for the selected cell type.

        ``block=True`` (default, queue) runs synchronously; ``block=False``
        (GUI) executes in a background worker with sample-level progress.

        ``extra_callbacks`` is the queue's hook for chaining
        ``{"on_done": cb, "on_failed": cb}``; see
        :meth:`PixelClassifierWidget.run_batch_segmentation` for details.
        """
        from behav3d.preprocessing.segmentation.cellpose_prediction import (
            run_cellpose_and_sync_metadata,
        )
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run Cellpose: No metadata loaded.")
            fire_extra_callback(extra_callbacks, "on_failed", "no metadata loaded")
            return

        if not block and self._bg.is_running():
            self.log("⚠️ A cellpose run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        cell_type = cell_type_override or self.cell_type_combo.currentText()
        model_path = model_path_override or self.pretrained_model_dir

        if not model_path:
            self.log("⚠️ Cannot run Cellpose: No model loaded.")
            return
        if not cell_type or cell_type == "(load metadata)":
            self.log("⚠️ Cannot run Cellpose: No cell type selected.")
            return

        self._persist_channel_config()

        if interactive:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
            md = self.metadata_loader.metadata
            existing = []
            if md is not None:
                out_dir = Path(self.metadata_loader.output_dir)
                for sn in md["sample_name"].unique():
                    seg_path = out_dir / "images" / sn / f"{sn}_{cell_type}_segments.zarr"
                    if seg_path.exists():
                        existing.append(f"{cell_type} segments for {sn}")
            if not existing:
                existing = [f"existing {cell_type} Cellpose results"]
            choice = prompt_overwrite_batch(
                self,
                "Existing Cellpose Results",
                existing,
            )
            if choice == "cancel":
                self.log("Cellpose segmentation cancelled.")
                fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                return
            overwrite = (choice == "overwrite")
        else:
            overwrite = not skip_existing

        # Timepoint range (Qt widget reads on caller thread).
        if self.check_process_all.isChecked():
            timepoint_range = None
        else:
            t_start = self.spin_t_start.value()
            t_end = self.spin_t_end.value()
            if t_start > t_end:
                self.log("Error: Start timepoint must be <= End timepoint.")
                fire_extra_callback(extra_callbacks, "on_failed", "bad timepoint range")
                return
            timepoint_range = (t_start, t_end)

        self.log(f"Starting Cellpose segmentation for '{cell_type}'...")

        def _do_cellpose(progress_cb=None):
            return run_cellpose_and_sync_metadata(
                output_dir=self.metadata_loader.output_dir,
                metadata_loader=self.metadata_loader,
                pretrained_model_dir=model_path,
                input_channels=[cell_type],
                label_name=cell_type,
                timepoint_range=timepoint_range,
                progress_cb=progress_cb,
            )

        def _report_summary(summary):
            n_proc = len(summary["processed"])
            n_skip = len(summary["skipped"])
            if n_skip > 0:
                self.log(f"Cellpose finished: {n_proc} processed, {n_skip} skipped.")
            else:
                self.log(f"Cellpose finished for all {n_proc} samples.")
            self.log("Metadata updated.")

        if block:
            try:
                result = _do_cellpose(progress_cb=None)
                _, summary = result
                _report_summary(summary)
                fire_extra_callback(extra_callbacks, "on_done", result)
            except Exception as e:
                traceback.print_exc()
                self.log(f"❌ Cellpose error: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
                if not interactive:
                    raise
            return

        def _on_done(result):
            _, summary = result
            _report_summary(summary)
            fire_extra_callback(extra_callbacks, "on_done", result)

        def _on_failed(err: str):
            self.log(f"❌ Cellpose error: {err}")
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_cellpose,
            desc=f"Cellpose \u2014 {cell_type}\u2026",
            progress_row=self.tab_progress_row,
            buttons=[],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    # ── run Otsu dead mask ──────────────────────────────────────────────
    def run_otsu_threshold(self, interactive=True, block=True,
                           extra_callbacks=None):
        """Run Otsu thresholding on the dead channel.

        ``block=True`` (default, queue) runs synchronously; ``block=False``
        (GUI) executes in a background worker with sample-level progress.
        ``extra_callbacks`` is the queue's chaining hook.
        """
        from behav3d.preprocessing.segmentation.cellpose_prediction import (
            run_otsu_threshold_segmentation_from_zarr,
        )
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run Otsu: No metadata loaded.")
            fire_extra_callback(extra_callbacks, "on_failed", "no metadata loaded")
            return

        if not block and self._bg.is_running():
            self.log("⚠️ A cellpose/otsu run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        self._persist_channel_config()

        # ---- Overwrite check ------------------------------------------------
        out_dir = Path(self.metadata_loader.output_dir)
        md = self.metadata_loader.metadata
        existing = []
        for sn in md["sample_name"].unique():
            mask_path = out_dir / "images" / sn / f"{sn}_mask_dead.zarr"
            if mask_path.exists():
                existing.append(f"{sn} dead mask ({mask_path.name})")

        overwrite_existing = True
        if existing and interactive:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
            choice = prompt_overwrite_batch(
                self,
                "Overwrite Existing Dead Masks?",
                existing,
            )
            if choice == "cancel":
                self.log("Otsu dead mask cancelled.")
                fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                return
            overwrite_existing = (choice == "overwrite")

        self.log("Starting Otsu dead mask segmentation...")

        # Timepoint range (Qt widget reads on caller thread).
        if self.check_process_all.isChecked():
            timepoint_range = None
        else:
            timepoint_range = (self.spin_t_start.value(), self.spin_t_end.value())

        if not overwrite_existing:
            self.log(
                "\u26a0\ufe0f Backend does not yet support partial skipping; "
                "existing dead masks will be overwritten."
            )

        def _do_otsu(progress_cb=None):
            return run_otsu_threshold_segmentation_from_zarr(
                output_dir=self.metadata_loader.output_dir,
                metadata=self.metadata_loader.metadata,
                timepoint_range=timepoint_range,
                progress_cb=progress_cb,
            )

        def _apply_otsu(result):
            updated_metadata, summary = result
            self.metadata_loader.metadata = updated_metadata
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                updated_metadata.to_csv(csv_path, index=False)
            n_proc = len(summary["processed"])
            self.log(f"Otsu dead mask finished: {n_proc} samples processed.")

        if block:
            try:
                result = _do_otsu(progress_cb=None)
                _apply_otsu(result)
                fire_extra_callback(extra_callbacks, "on_done", result)
            except Exception as e:
                traceback.print_exc()
                self.log(f"❌ Otsu error: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
            return

        def _on_done_async(result):
            _apply_otsu(result)
            fire_extra_callback(extra_callbacks, "on_done", result)

        def _on_failed(err: str):
            self.log(f"❌ Otsu error: {err}")
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_otsu,
            desc="Otsu dead mask\u2026",
            progress_row=self.tab_progress_row,
            buttons=[],
            viewer=self.viewer,
            on_done=_on_done_async,
            on_failed=_on_failed,
        )


class ImportWidget(QWidget):
    """Widget to browse-in, convert, and validate pre-existing segmentation
    (and dead-mask) files.

    Shows a per-sample, per-cell-type list of editable path rows (prefilled
    from metadata when a path is already set) with live status/conversion
    actions. Newly browsed/typed paths are staged in the row's own widget —
    metadata.csv is only updated once the corresponding Convert/Save action
    actually runs (or a batch Convert/Save-All).
    """

    def __init__(self, viewer, metadata_loader, log_callback=None,
                 switch_to_data_prep_edit_callback=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self._switch_to_data_prep_edit = switch_to_data_prep_edit_callback
        self.organoid_types = []
        self.immune_types = []
        self.other_types = []
        self.all_cell_types = []
        # (sample_name, cell_type) -> {"path_edit", "browse_btn",
        #  "status_container", "status_layout", "last_value", "row_idx"}
        self._rows = {}
        # sample_name -> same widget dict, for the dead-mask row
        self._dead_rows = {}
        self._init_ui()

    # ── helpers ──────────────────────────────────────────────────────────
    def _get_prefix(self, cell_type):
        """Return the metadata column prefix for a cell type."""
        if cell_type in self.organoid_types:
            return "or"
        if cell_type in self.immune_types:
            return "im"
        return "ot"

    def _seg_col(self, cell_type):
        """Full metadata column name for a cell type's segmentation path."""
        return f"{self._get_prefix(cell_type)}_{cell_type}_segments_image_path"

    def _expected_outpath(self, sample_name, cell_type):
        """Where a converted zarr would be saved (matches batch segmentation)."""
        out_dir = Path(self.metadata_loader.output_dir)
        return out_dir / "images" / sample_name / f"{sample_name}_{cell_type}_segments.zarr"

    def _expected_dead_mask_outpath(self, sample_name):
        """Where a converted dead-mask zarr would be saved."""
        out_dir = Path(self.metadata_loader.output_dir)
        return out_dir / "images" / sample_name / f"{sample_name}_dead_mask.zarr"

    def _metadata_csv_path(self):
        path = getattr(self.metadata_loader, "_loaded_csv_path", None)
        if path:
            return path
        return self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")

    def _resolve_path(self, path_str):
        """Resolve a path string, trying the metadata CSV directory if not absolute/found."""
        return resolve_external_path(path_str, self._metadata_csv_path())

    @staticmethod
    def _clear_layout(layout):
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    # ── zarr validation ─────────────────────────────────────────────────
    @staticmethod
    def _check_zarr_structure(path):
        """Check that a .zarr file matches the expected save_as_zarr structure.

        Expected: single root array (no sub-groups), chunks starting with (1, ...).
        Returns (ok: bool, reason: str).
        """
        import zarr
        try:
            store = zarr.storage.LocalStore(str(path))
            root = zarr.open(store, mode="r")
            # Must be an array at root, not a group with children
            if isinstance(root, zarr.Group):
                # Check if group contains a single array
                members = list(root.arrays())
                if len(members) == 0:
                    return False, "Zarr group contains no arrays"
                return False, "Zarr contains sub-groups instead of a root array"
            # root is zarr.Array
            if root.chunks[0] != 1:
                return False, f"First chunk dim is {root.chunks[0]}, expected 1"
            return True, "OK"
        except Exception as e:
            return False, str(e)

    # ── UI ──────────────────────────────────────────────────────────────
    def _init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        if self._switch_to_data_prep_edit is not None:
            add_row = QHBoxLayout()
            btn_add = QPushButton("➕  Add a new sample or cell type…")
            btn_add.setToolTip(
                "Jumps to the Data Preparation tab's Metadata Builder "
                "(already in edit mode) to add samples/cell types that "
                "don't exist in metadata yet."
            )
            btn_add.setStyleSheet(
                "QPushButton{background:#455A64;color:white;padding:6px 12px;"
                "border-radius:3px}"
                "QPushButton:hover{background:#546E7A}"
            )
            btn_add.clicked.connect(self._switch_to_data_prep_edit)
            add_row.addWidget(btn_add)
            add_row.addStretch()
            layout.addLayout(add_row)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scroll_content = QWidget()
        self.scroll_layout = QVBoxLayout(self.scroll_content)
        self.scroll_layout.setAlignment(Qt.AlignTop)
        scroll.setWidget(self.scroll_content)
        layout.addWidget(scroll)

        # placeholder
        self.placeholder = QLabel("Load metadata to see segmentation import status.")
        self.placeholder.setStyleSheet("color:#888; font-style:italic; padding:20px;")
        self.placeholder.setAlignment(Qt.AlignCenter)
        self.scroll_layout.addWidget(self.placeholder)

    def _on_metadata_updated(self, _metadata=None):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            is_combined_multicolor_celltype,
        )
        md = self.metadata_loader.metadata
        if md is None:
            return
        self.organoid_types = [ct for ct in detect_organoid_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        self.immune_types = [ct for ct in detect_immune_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        self.other_types = [ct for ct in detect_other_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types
        self._rebuild_table()

    # ── table builder ───────────────────────────────────────────────────
    def _rebuild_table(self):
        """Clear and rebuild the entire status table from metadata."""
        self._rows = {}
        self._dead_rows = {}

        # clear old widgets
        while self.scroll_layout.count():
            item = self.scroll_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        md = self.metadata_loader.metadata
        if md is None or md.empty:
            lbl = QLabel("No metadata loaded.")
            lbl.setStyleSheet("color:#888; font-style:italic; padding:20px;")
            self.scroll_layout.addWidget(lbl)
            return

        for idx, row in md.iterrows():
            sample_name = str(row.get("sample_name", f"Row {idx+1}"))
            self._add_sample_section(sample_name, idx, row)

        # Check if ANY button was added to the scroll layout
        any_fixes_needed = False
        for i in range(self.scroll_layout.count()):
            w = self.scroll_layout.itemAt(i).widget()
            if isinstance(w, QWidget) and w.findChildren(QPushButton):
                any_fixes_needed = True
                break

        if any_fixes_needed:
            # ── global Convert All button ──
            btn_all = QPushButton("⚡  Convert / Save All Samples")
            btn_all.setStyleSheet(
                "QPushButton{background:#1565C0;color:white;padding:8px 16px;"
                "border-radius:4px;font-weight:bold;font-size:13px}"
                "QPushButton:hover{background:#1976D2}"
            )
            btn_all.clicked.connect(self._convert_all)
            self.scroll_layout.addWidget(btn_all)
        else:
            # All correct -> Show success message
            success_msg = QLabel("✨ All available segmentation files are in the correct format.")
            success_msg.setWordWrap(True)
            success_msg.setAlignment(Qt.AlignCenter)
            success_msg.setStyleSheet("color:#2E7D32; background:#E8F5E9; padding:12px; border-radius:8px; margin:5px; border:1px solid #C8E6C9; font-weight:bold;")
            self.scroll_layout.addWidget(success_msg)

        # ── Save-timing note (ALWAYS visible at the bottom) ──
        instr = QLabel(
            "Paths typed or browsed above are only written to metadata.csv "
            "once you click Convert/Save for that row (or Convert/Save All)."
        )
        instr.setWordWrap(True)
        instr.setAlignment(Qt.AlignCenter)
        instr.setStyleSheet("color:#555; background:#f9f9f9; padding:15px; border-radius:8px; margin:10px; border:1px solid #ddd;")
        self.scroll_layout.addWidget(instr)

        self.scroll_layout.addStretch()

    def _add_sample_section(self, sample_name, row_idx, row):
        """Add a grouped section for one sample: one editable row per cell
        type, plus a dead-mask row when applicable."""
        from behav3d.core.metadata import has_dead_channel

        header = QLabel(f"📁  {sample_name}")
        header.setStyleSheet("font-weight:bold; font-size:13px; padding:6px 0 2px 0;")
        self.scroll_layout.addWidget(header)

        for ct in self.all_cell_types:
            self._add_cell_type_row(sample_name, row_idx, row, ct)

        md = self.metadata_loader.metadata
        if has_dead_channel(md):
            self._add_dead_mask_row(sample_name, row_idx, row)

        # Check if any "action" buttons were added for this sample
        fixes_needed = False
        found_header = False
        for i in range(self.scroll_layout.count()):
            w = self.scroll_layout.itemAt(i).widget()
            if not found_header:
                if w is header:
                    found_header = True
                continue
            if isinstance(w, QWidget) and w.findChildren(QPushButton):
                fixes_needed = True
                break

        if fixes_needed:
            btn = QPushButton(f"⚡  Convert / Save all for {sample_name}")
            btn.setStyleSheet(
                "QPushButton{background:#2E7D32;color:white;padding:5px 12px;"
                "border-radius:3px;font-size:12px}"
                "QPushButton:hover{background:#388E3C}"
            )
            btn.clicked.connect(
                lambda _checked=False, s=sample_name, r=row_idx: self._convert_sample(s, r)
            )
            wrap = QWidget()
            wrap_lay = QHBoxLayout(wrap)
            wrap_lay.setContentsMargins(16, 2, 4, 6)
            wrap_lay.addWidget(btn)
            wrap_lay.addStretch()
            self.scroll_layout.addWidget(wrap)

        # separator
        sep = QWidget()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background:#ddd;")
        self.scroll_layout.addWidget(sep)

    def _add_cell_type_row(self, sample_name, row_idx, row, cell_type):
        """Build+register one editable segmentation-path row for one
        (sample, cell_type) pair."""
        col = self._seg_col(cell_type)
        raw_val = row.get(col) if col in row.index else None
        has_value = raw_val is not None and pd.notna(raw_val) and str(raw_val).strip() not in ("", "nan")
        initial_value = str(raw_val).strip().strip('"').strip("'") if has_value else ""

        row_widget, path_edit, browse_btn, status_layout = self._build_path_row(f"{cell_type}:", initial_value)

        key = (sample_name, cell_type)
        self._rows[key] = {
            "path_edit": path_edit,
            "browse_btn": browse_btn,
            "status_layout": status_layout,
            "last_value": initial_value,
            "row_idx": row_idx,
        }
        path_edit.editingFinished.connect(partial(self._on_row_path_edited, sample_name, cell_type))
        browse_btn.clicked.connect(partial(self._on_browse_clicked, sample_name, cell_type))

        self.scroll_layout.addWidget(row_widget)
        self._refresh_row_status(sample_name, cell_type)

    def _add_dead_mask_row(self, sample_name, row_idx, row):
        """Build+register one editable dead-mask-path row for one sample."""
        col = "dead_mask_path"
        raw_val = row.get(col) if col in row.index else None
        has_value = raw_val is not None and pd.notna(raw_val) and str(raw_val).strip() not in ("", "nan")
        initial_value = str(raw_val).strip().strip('"').strip("'") if has_value else ""

        row_widget, path_edit, browse_btn, status_layout = self._build_path_row("dead mask:", initial_value)

        self._dead_rows[sample_name] = {
            "path_edit": path_edit,
            "browse_btn": browse_btn,
            "status_layout": status_layout,
            "last_value": initial_value,
            "row_idx": row_idx,
        }
        path_edit.editingFinished.connect(partial(self._on_dead_mask_path_edited, sample_name))
        browse_btn.clicked.connect(partial(self._on_dead_mask_browse_clicked, sample_name))

        self.scroll_layout.addWidget(row_widget)
        self._refresh_dead_mask_status(sample_name)

    @staticmethod
    def _build_path_row(label_text, initial_value):
        """Build the shared visual shape of an editable path row (label,
        path field, Browse button, status container) without wiring any
        signals or registering it anywhere — callers do that."""
        row_widget = QWidget()
        row_lay = QHBoxLayout(row_widget)
        row_lay.setContentsMargins(16, 2, 4, 2)

        lbl = QLabel(label_text)
        lbl.setFixedWidth(120)
        row_lay.addWidget(lbl)

        path_edit = QLineEdit(initial_value)
        path_edit.setPlaceholderText("Path to .tif/.tiff file or .zarr directory")
        path_edit.setMinimumWidth(220)
        row_lay.addWidget(path_edit, stretch=1)

        browse_btn = QPushButton("Browse…")
        row_lay.addWidget(browse_btn)

        status_container = QWidget()
        status_layout = QHBoxLayout(status_container)
        status_layout.setContentsMargins(8, 0, 0, 0)
        row_lay.addWidget(status_container)
        row_lay.addStretch()

        return row_widget, path_edit, browse_btn, status_layout

    # ── live status refresh (reads from the widget, not metadata) ───────
    def _refresh_row_status(self, sample_name, cell_type):
        info = self._rows.get((sample_name, cell_type))
        if info is None:
            return
        self._clear_layout(info["status_layout"])

        path_str = info["path_edit"].text().strip().strip('"').strip("'")
        row_idx = info["row_idx"]

        if not path_str:
            status = QLabel("No segmentation available")
            status.setStyleSheet("color:#999; font-style:italic;")
            info["status_layout"].addWidget(status)
            return

        file_path = self._resolve_path(path_str)
        if file_path is None or not file_path.exists():
            status = QLabel("⚠️  File not found")
            status.setToolTip(str(file_path))
            status.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(status)
            return

        if file_path.suffix == ".zarr" or file_path.is_dir():
            ok, reason = self._check_zarr_structure(file_path)
            if not ok:
                btn = QPushButton("🔄  Fix zarr format")
                btn.setToolTip(f"Issue: {reason}")
                btn.setStyleSheet(
                    "QPushButton{background:#F57C00;color:white;padding:4px 10px;"
                    "border-radius:3px}"
                    "QPushButton:hover{background:#FB8C00}"
                )
                btn.clicked.connect(
                    lambda _checked=False, s=sample_name, c=cell_type, r=row_idx: self._convert_single(s, c, r)
                )
                info["status_layout"].addWidget(btn)
                return

            # --- Dimension check against the sample's raw image ---
            dims_match = True
            try:
                md = self.metadata_loader.metadata
                raw_path_str = md.at[row_idx, "raw_image_path"] if "raw_image_path" in md.columns else ""
                raw_path = Path(str(raw_path_str)) if raw_path_str else None
                if raw_path is not None and raw_path.exists():
                    import zarr
                    raw_store = zarr.open(str(raw_path), mode='r')
                    seg_store = zarr.open(str(file_path), mode='r')
                    raw_shape = raw_store.shape
                    seg_shape = seg_store.shape
                    if len(seg_shape) < 3 or raw_shape[-3:] != seg_shape[-3:]:
                        dims_match = False
                        reason = f"Spatial mismatch: Raw {raw_shape[-3:]} vs Seg {seg_shape[-3:]}"
                    elif raw_shape[0] != seg_shape[0]:
                        dims_match = False
                        reason = f"Time mismatch: Raw {raw_shape[0]}T vs Seg {seg_shape[0]}T"
            except Exception:
                pass

            if not dims_match:
                status = QLabel("⚠️  Dimension mismatch")
                status.setToolTip(reason)
                status.setStyleSheet("color:#E65100; font-weight:bold;")
                info["status_layout"].addWidget(status)
            elif str(file_path) == info["last_value"]:
                status = QLabel("✅  Ready for tracking")
                status.setStyleSheet("color:#2E7D32; font-weight:bold;")
                info["status_layout"].addWidget(status)
            else:
                # Already a valid zarr, just not yet persisted to metadata.
                btn = QPushButton("💾  Save path to metadata")
                btn.setStyleSheet(
                    "QPushButton{background:#2E7D32;color:white;padding:4px 10px;"
                    "border-radius:3px}"
                    "QPushButton:hover{background:#388E3C}"
                )
                btn.clicked.connect(
                    lambda _checked=False, s=sample_name, c=cell_type, r=row_idx: self._convert_single(s, c, r)
                )
                info["status_layout"].addWidget(btn)

        elif file_path.suffix.lower() in (".tif", ".tiff"):
            btn = QPushButton("🔄  Convert TIFF → zarr")
            btn.setStyleSheet(
                "QPushButton{background:#1565C0;color:white;padding:4px 10px;"
                "border-radius:3px}"
                "QPushButton:hover{background:#1976D2}"
            )
            btn.clicked.connect(
                lambda _checked=False, s=sample_name, c=cell_type, r=row_idx: self._convert_single(s, c, r)
            )
            info["status_layout"].addWidget(btn)
        else:
            status = QLabel(f"⚠️  Format not supported ({file_path.suffix})")
            status.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(status)

    def _refresh_dead_mask_status(self, sample_name):
        info = self._dead_rows.get(sample_name)
        if info is None:
            return
        self._clear_layout(info["status_layout"])

        path_str = info["path_edit"].text().strip().strip('"').strip("'")
        row_idx = info["row_idx"]

        if not path_str:
            status = QLabel("No dead mask available")
            status.setStyleSheet("color:#999; font-style:italic;")
            info["status_layout"].addWidget(status)
            return

        file_path = self._resolve_path(path_str)
        if file_path is None or not file_path.exists():
            status = QLabel("⚠️  File not found")
            status.setToolTip(str(file_path))
            status.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(status)
            return

        if file_path.suffix == ".zarr" or file_path.is_dir():
            ok, reason = self._check_zarr_structure(file_path)
            if not ok:
                btn = QPushButton("🔄  Fix zarr format")
                btn.setToolTip(f"Issue: {reason}")
                btn.setStyleSheet(
                    "QPushButton{background:#F57C00;color:white;padding:4px 10px;"
                    "border-radius:3px}"
                    "QPushButton:hover{background:#FB8C00}"
                )
                btn.clicked.connect(
                    lambda _checked=False, s=sample_name, r=row_idx: self._convert_dead_mask(s, r)
                )
                info["status_layout"].addWidget(btn)
            elif str(file_path) == info["last_value"]:
                status = QLabel("✅  Ready")
                status.setStyleSheet("color:#2E7D32; font-weight:bold;")
                info["status_layout"].addWidget(status)
            else:
                btn = QPushButton("💾  Save path to metadata")
                btn.setStyleSheet(
                    "QPushButton{background:#2E7D32;color:white;padding:4px 10px;"
                    "border-radius:3px}"
                    "QPushButton:hover{background:#388E3C}"
                )
                btn.clicked.connect(
                    lambda _checked=False, s=sample_name, r=row_idx: self._convert_dead_mask(s, r)
                )
                info["status_layout"].addWidget(btn)
        elif file_path.suffix.lower() in (".tif", ".tiff"):
            btn = QPushButton("🔄  Convert TIFF → zarr")
            btn.setStyleSheet(
                "QPushButton{background:#1565C0;color:white;padding:4px 10px;"
                "border-radius:3px}"
                "QPushButton:hover{background:#1976D2}"
            )
            btn.clicked.connect(
                lambda _checked=False, s=sample_name, r=row_idx: self._convert_dead_mask(s, r)
            )
            info["status_layout"].addWidget(btn)
        else:
            status = QLabel(f"⚠️  Format not supported ({file_path.suffix})")
            status.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(status)

    # ── path-field change handlers ───────────────────────────────────────
    def _on_browse_clicked(self, sample_name, cell_type):
        info = self._rows.get((sample_name, cell_type))
        if info is None:
            return
        new_path = browse_file_or_zarr(
            self, f"Select {cell_type} segmentation for {sample_name}",
            "Image files (*.tif *.tiff *.zarr);; All Files (*)",
            allow_zarr=True,
        )
        if not new_path:
            return
        self._maybe_accept_new_value(
            info, new_path,
            label=f"{cell_type} / {sample_name}",
            on_accept=lambda: self._refresh_row_status(sample_name, cell_type),
        )

    def _on_row_path_edited(self, sample_name, cell_type):
        info = self._rows.get((sample_name, cell_type))
        if info is None:
            return
        new_value = info["path_edit"].text().strip()
        self._maybe_accept_new_value(
            info, new_value,
            label=f"{cell_type} / {sample_name}",
            on_accept=lambda: self._refresh_row_status(sample_name, cell_type),
            already_in_field=True,
        )

    def _on_dead_mask_browse_clicked(self, sample_name):
        info = self._dead_rows.get(sample_name)
        if info is None:
            return
        new_path = browse_file_or_zarr(
            self, f"Select dead mask for {sample_name}",
            "Image files (*.tif *.tiff *.zarr);; All Files (*)",
            allow_zarr=True,
        )
        if not new_path:
            return
        self._maybe_accept_new_value(
            info, new_path,
            label=f"dead mask / {sample_name}",
            on_accept=lambda: self._refresh_dead_mask_status(sample_name),
        )

    def _on_dead_mask_path_edited(self, sample_name):
        info = self._dead_rows.get(sample_name)
        if info is None:
            return
        new_value = info["path_edit"].text().strip()
        self._maybe_accept_new_value(
            info, new_value,
            label=f"dead mask / {sample_name}",
            on_accept=lambda: self._refresh_dead_mask_status(sample_name),
            already_in_field=True,
        )

    def _maybe_accept_new_value(self, info, new_value, *, label, on_accept, already_in_field=False):
        """Shared change-acceptance logic for both cell-type and dead-mask
        rows: confirms overwriting an already-set metadata value, then
        (re)writes the field and refreshes that row's status."""
        new_value = str(new_value).strip()
        old_value = info["last_value"]

        if new_value != old_value and old_value:
            res = QMessageBox.question(
                self, "Replace existing path?",
                f"This will replace the existing path for {label}:\n{old_value}\n→\n{new_value}\n\nContinue?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                info["path_edit"].blockSignals(True)
                info["path_edit"].setText(old_value)
                info["path_edit"].blockSignals(False)
                return

        if not already_in_field:
            info["path_edit"].blockSignals(True)
            info["path_edit"].setText(new_value)
            info["path_edit"].blockSignals(False)

        on_accept()

    # ── conversion logic ────────────────────────────────────────────────
    def _row_action(self, info, col, dest_path, row_idx, *, save_metadata, refresh_ui, label):
        """Shared conversion engine for one row (cell-type or dead-mask):
        saves an already-valid zarr path as-is, repairs a malformed zarr, or
        converts a TIFF (with axis-order confirmation) — then writes the
        result into metadata. Returns True if metadata changed."""
        path_str = info["path_edit"].text().strip().strip('"').strip("'")
        if not path_str:
            return False
        src = self._resolve_path(path_str)
        if src is None or not src.exists():
            return False

        md = self.metadata_loader.metadata

        if src.suffix == ".zarr" or src.is_dir():
            ok, _ = self._check_zarr_structure(src)
            if ok:
                if str(src) == info["last_value"]:
                    return False  # already saved, nothing to do
                if col not in md.columns:
                    md[col] = pd.NA
                md[col] = md[col].astype("object")
                md.at[row_idx, col] = str(src)
                info["last_value"] = str(src)
                self.log(f"  {col} for {label} -> {src}")
                if save_metadata:
                    self._save_metadata(refresh_ui=refresh_ui)
                return True

            # bad zarr structure -> repair via raw save_as_zarr (unchanged behavior)
            if dest_path.exists():
                res = QMessageBox.question(
                    self, "Overwrite?",
                    f"There is a pre-existing file:\n{dest_path}\n\nDo you want to overwrite it?",
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
                )
                if res != QMessageBox.Yes:
                    self.log(f"Skipped {label} (user declined overwrite)")
                    return False
            try:
                dest_path.parent.mkdir(parents=True, exist_ok=True)
                self.log(f"Repairing zarr format: {src.name} → {dest_path.name} ...")
                img = load_image(src)
                if dest_path.exists():
                    shutil.rmtree(dest_path)
                save_as_zarr(img, dest_path)
            except Exception as e:
                traceback.print_exc()
                self.log(f"❌  Error repairing {label}: {e}")
                if refresh_ui:
                    QMessageBox.warning(self, "Conversion Error", str(e))
                return False

        elif src.suffix.lower() in (".tif", ".tiff"):
            if dest_path.exists():
                res = QMessageBox.question(
                    self, "Overwrite?",
                    f"There is a pre-existing file:\n{dest_path}\n\nDo you want to overwrite it?",
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
                )
                if res != QMessageBox.Yes:
                    self.log(f"Skipped {label} (user declined overwrite)")
                    return False
            axis_order = prompt_axis_order(self, src, log=self.log)
            if axis_order is None:
                self.log(f"Cancelled: no axis order selected for {label}")
                return False
            try:
                self.log(f"Converting {src.name} → {dest_path.name} ...")
                convert_label_file_to_zarr(path=src, outpath=dest_path, axis_order=axis_order, overwrite=True)
            except Exception as e:
                traceback.print_exc()
                self.log(f"❌  Error converting {label}: {e}")
                if refresh_ui:
                    QMessageBox.warning(self, "Conversion Error", str(e))
                return False
        else:
            self.log(f"⚠️  Unsupported format for {label}: {src.suffix}")
            return False

        self.log(f"✅  Saved: {dest_path}")
        if col not in md.columns:
            md[col] = pd.NA
        md[col] = md[col].astype("object")
        md.at[row_idx, col] = str(dest_path)
        info["last_value"] = str(dest_path)
        if save_metadata:
            self._save_metadata(refresh_ui=refresh_ui)
        return True

    def _convert_single(self, sample_name, cell_type, row_idx, *, save_metadata=True, refresh_ui=True):
        """Convert/save one cell type's segmentation row for one sample."""
        info = self._rows.get((sample_name, cell_type))
        if info is None:
            return False
        dest = self._expected_outpath(sample_name, cell_type)
        col = self._seg_col(cell_type)
        return self._row_action(
            info, col, dest, row_idx,
            save_metadata=save_metadata, refresh_ui=refresh_ui,
            label=f"{cell_type} / {sample_name}",
        )

    def _convert_dead_mask(self, sample_name, row_idx, *, save_metadata=True, refresh_ui=True):
        """Convert/save the dead-mask row for one sample."""
        info = self._dead_rows.get(sample_name)
        if info is None:
            return False
        dest = self._expected_dead_mask_outpath(sample_name)
        return self._row_action(
            info, "dead_mask_path", dest, row_idx,
            save_metadata=save_metadata, refresh_ui=refresh_ui,
            label=f"dead mask / {sample_name}",
        )

    def _convert_sample(self, sample_name, row_idx, *, save_metadata=True, refresh_ui=True):
        """Convert/save all actionable rows (cell types + dead mask) for one sample."""
        converted_any = False

        for ct in self.all_cell_types:
            if (sample_name, ct) in self._rows:
                # refresh_ui=False to avoid rebuilding the table mid-loop, which crashes.
                if self._convert_single(sample_name, ct, row_idx, save_metadata=False, refresh_ui=False):
                    converted_any = True

        if sample_name in self._dead_rows:
            if self._convert_dead_mask(sample_name, row_idx, save_metadata=False, refresh_ui=False):
                converted_any = True

        if converted_any:
            if save_metadata:
                self._save_metadata(refresh_ui=refresh_ui)
        elif refresh_ui:
            self.log(f"No conversions needed for {sample_name}")

        return converted_any

    def _convert_all(self):
        """Convert/save all actionable rows across all samples."""
        md = self.metadata_loader.metadata
        if md is None:
            return

        converted_any = False
        for idx, row in md.iterrows():
            sample_name = str(row.get("sample_name", f"Row {idx+1}"))
            # refresh_ui=False to avoid crash
            if self._convert_sample(sample_name, idx, save_metadata=False, refresh_ui=False):
                converted_any = True

        if converted_any:
            self._save_metadata(refresh_ui=True)
        else:
            self.log("No conversions needed.")

    def _save_metadata(self, show_popup=True, refresh_ui=True):
        """Save metadata CSV to disk and show info popup."""
        md = self.metadata_loader.metadata
        csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv", "")
        if not csv_path:
            self.log("⚠️  Cannot save metadata: CSV path unknown")
            return
        md.to_csv(csv_path, index=False)
        self.log(f"Metadata updated: {csv_path}")
        if show_popup:
            QMessageBox.information(self, "Metadata Updated",
                                    f"The metadata CSV has been updated:\n{csv_path}")
        if refresh_ui:
            self._rebuild_table()


# ═══════════════════════════════════════════════════════════════════════════
# APOC Segmentation Widget — GPU-accelerated pixel classification
# ═══════════════════════════════════════════════════════════════════════════
class APOCWidget(QWidget):
    """APOC (GPU) segmentation page with embedded training and batch inference.

    All training UI (per-cell-type tabs, strategy combo, instance controls,
    per-tab strategy combos) is provided natively by
    :class:`behav3d.preprocessing.segmentation.apoc_train.APOCTrainingWidget`
    through its ``per_cell_type_strategy=True`` / ``instance_controls_mode="inline"``
    constructor flags. The plugin only wires GPU/device selection, batch run
    options, queueing, and metadata interactions on top of that widget.
    """

    # The training widget exposes its canonical list as ``APOCTrainingWidget.STRATEGIES``
    # plus the ``ADVANCED_STRATEGY`` constant. We compose both here for convenience
    # so other modules (e.g. queue helpers) keep working unchanged.
    @classmethod
    def all_strategies(cls):
        from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget as _ATW
        return list(_ATW.STRATEGIES) + [_ATW.ADVANCED_STRATEGY]

    @classmethod
    def per_tab_strategies(cls):
        from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget as _ATW
        return list(_ATW.STRATEGIES)

    def __init__(
        self,
        viewer,
        metadata_loader,
        log_callback=None,
        parent=None,
        switch_to_visualization_callback=None,
        tab_progress_row=None,
    ):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self._switch_to_visualization = switch_to_visualization_callback
        self._training_widget = None
        self._is_session_active = False
        self._organoid_types = []          # set when metadata loaded
        self.all_cell_types = []           # set when metadata loaded
        self._unify_organoids = False      # state of the unified-organoid checkbox
        self._syncing_organoids = False    # recursion guard for continuous organoid sync
        # Background-execution infrastructure (shared progress row owned
        # by the parent SegmentationTab).
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        self._init_ui()

    def _init_ui(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        # ── Info banner ─────────────────────────────────────────
        info = QLabel(
            "<b>APOC GPU Segmentation</b><br>"
            "Uses <tt>apoc</tt> + <tt>pyclesperanto</tt> for GPU-accelerated<br>"
            "pixel classification and instance segmentation."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color: #888; font-size: 11px; padding: 4px 0 6px 0;")
        layout.addWidget(info)

        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {}) or {}

        # Discover available OpenCL devices and restore last saved choice.
        self._apoc_gpu_devices = self._discover_apoc_gpu_devices()

        gpu_row = QHBoxLayout()
        gpu_row.addWidget(QLabel("GPU Device:"))
        self.combo_gpu_device = QComboBox()
        self.combo_gpu_device.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        if self._apoc_gpu_devices:
            self.combo_gpu_device.addItems(self._apoc_gpu_devices)
            saved_gpu = str(pc.get("gpu_device_name", "")).strip()
            if saved_gpu and saved_gpu in self._apoc_gpu_devices:
                self.combo_gpu_device.setCurrentText(saved_gpu)
        else:
            self.combo_gpu_device.addItem("Auto (default OpenCL device)")
            self.combo_gpu_device.setEnabled(False)
        gpu_row.addWidget(self.combo_gpu_device, stretch=1)
        gpu_row.addStretch()
        layout.addLayout(gpu_row)

        force_cpu_row = QHBoxLayout()
        self.btn_force_cpu = QCheckBox("Force CPU-only processing")
        self.btn_force_cpu.setToolTip("Override GPU selection and run pyclesperanto on the CPU")
        is_force_cpu = bool(pc.get("force_cpu", False))
        self.btn_force_cpu.setChecked(is_force_cpu)
        if is_force_cpu:
            self.combo_gpu_device.setEnabled(False)
        force_cpu_row.addWidget(self.btn_force_cpu)
        force_cpu_row.addWidget(HelpButton(
            "Force CPU-only Processing",
            "Overrides the GPU Device selection above and asks pyclesperanto "
            "to run APOC (pixel classification + EDT/watershed) on the CPU "
            "instead of the GPU.\n\n"
            "Requires an OpenCL CPU runtime to be installed (e.g. the Intel "
            "CPU Runtime for OpenCL). If none is found, pyclesperanto silently "
            "falls back to whatever OpenCL device it can find and a warning "
            "is logged.\n\n"
            "Use this if you don't have a compatible GPU, or if GPU "
            "processing is unstable/unavailable on this machine."
        ))
        force_cpu_row.addStretch()
        layout.addLayout(force_cpu_row)

        # ── Training section (embedded APOCTrainingWidget) ──────
        self.training_group = QGroupBox("🎯 APOC Classifier Training")
        self.training_layout = QVBoxLayout(self.training_group)
        self.training_layout.setContentsMargins(4, 4, 4, 4)

        # Strategy config is now inside APOCTrainingWidget
        
        train_ctrl_lay = QHBoxLayout()
        train_ctrl_lay.addWidget(QLabel("Examples/sample:"))
        self.spin_examples = QSpinBox()
        self.spin_examples.setValue(int(pc.get("examples_per_sample", 3)))
        self.spin_examples.setRange(1, 10)
        self.spin_examples.setMaximumWidth(70)
        train_ctrl_lay.addWidget(self.spin_examples)
        
        self.btn_load_training = QPushButton("Generate Training Data")
        self.btn_load_training.setToolTip("Clears viewer and loads selected timepoints for labeling")
        self.btn_load_training.setStyleSheet("background-color: #007bff; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_load_training.clicked.connect(lambda _: self._on_load_training_clicked(interactive=True))
        train_ctrl_lay.addWidget(self.btn_load_training)
        train_ctrl_lay.addStretch()
        
        self.training_layout.addLayout(train_ctrl_lay)

        # ── Import Training Data row ──────────────────────────────────
        import_ctrl_lay = QHBoxLayout()
        self.btn_import_training = QPushButton("Import Training Data")
        self.btn_import_training.setToolTip(
            "Select a training_metadata.yml file from a previous BEHAV3D experiment.\n"
            "Features and sigma values will be locked to match the imported classifier."
        )
        self.btn_import_training.setStyleSheet(
            "background-color: #6f42c1; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        self.btn_import_training.clicked.connect(self._on_import_training_data_clicked)
        self.btn_import_training.setEnabled(False)   # enabled once metadata is loaded
        self.btn_import_training.setVisible(False)   # hidden until training data is generated
        import_ctrl_lay.addWidget(self.btn_import_training)
        self._btn_import_help = HelpButton(
            "Import Training Data",
            "Select a training_metadata.yml file from a previous BEHAV3D experiment to reuse its "
            "labeled training data (pixel features + labels).\n\n"
            "⚠️ Features and sigma values will be locked to match the imported classifier. "
            "This is required because the classifier was trained on a specific set of features "
            "computed at specific scales (sigma values) — changing them would make the classifier "
            "incompatible with the imported labels.\n\n"
            "⚠️ Even when importing, you must first click 'Generate Training Data' to load this "
            "experiment's images into the viewer. Importing replaces the need to manually label "
            "pixels — but the viewer session must be active first."
        )
        self._btn_import_help.setToolTip(
            "Requires: training_metadata.yml from a previous BEHAV3D experiment.\n"
            "Locks features and sigma values to match the imported classifier — because the "
            "classifier was trained on those exact features and cannot be used with different ones.\n"
            "Note: click 'Generate Training Data' first — importing replaces labeling, not image loading."
        )
        self._btn_import_help.setVisible(False)
        import_ctrl_lay.addWidget(self._btn_import_help)
        import_ctrl_lay.addStretch()
        self.training_layout.addLayout(import_ctrl_lay)

        # ── Imported data info panel (hidden until import is applied) ─
        self.import_info_group = QGroupBox("📦 Imported Training Data")
        self.import_info_group.setVisible(False)
        import_info_lay = QVBoxLayout(self.import_info_group)
        import_info_lay.setContentsMargins(6, 6, 6, 6)
        import_info_lay.setSpacing(4)

        self.import_source_label = QLabel("")
        self.import_source_label.setWordWrap(True)
        self.import_source_label.setStyleSheet("color: #666; font-size: 10px;")
        import_info_lay.addWidget(self.import_source_label)

        self.import_details_label = QLabel("")
        self.import_details_label.setWordWrap(True)
        self.import_details_label.setStyleSheet(
            "color: #444; font-size: 10px; font-family: monospace;"
        )
        import_info_lay.addWidget(self.import_details_label)

        btn_clear_import = QPushButton("Clear Import")
        btn_clear_import.setToolTip("Remove the imported training data and unlock feature controls.")
        btn_clear_import.clicked.connect(self._on_clear_import_clicked)
        import_info_lay.addWidget(btn_clear_import)

        self.training_layout.addWidget(self.import_info_group)

        # ── Export Training Bundle button ─────────────────────────────
        self.btn_export_bundle = QPushButton("Export Training Bundle…")
        self.btn_export_bundle.setToolTip(
            "Copy PixelClassifier_TrainingData.zip + .cl classifier files into a shareable archive."
        )
        self.btn_export_bundle.setEnabled(False)
        self.btn_export_bundle.clicked.connect(self._on_export_training_bundle)
        self.training_layout.addWidget(self.btn_export_bundle)

        # Legend container for APOC
        self.legend_container = QWidget()
        self.legend_container.setLayout(QVBoxLayout())
        self.legend_container.layout().setContentsMargins(0, 0, 0, 0)
        self.legend_container.setVisible(False)
        self.training_layout.addWidget(self.legend_container)

        self.training_placeholder = QLabel(
            "Load metadata to enable APOC classifier training."
        )
        self.training_placeholder.setWordWrap(True)
        self.training_placeholder.setStyleSheet("color:#888; font-style:italic; padding:10px;")
        self.training_layout.addWidget(self.training_placeholder)
        layout.addWidget(self.training_group)

        # (Strategy selector is now inside training_group above)

        # ── Batch controls ──────────────────────────────────────
        batch_group = QGroupBox("🚀 Batch Segmentation")
        batch_lay = QVBoxLayout(batch_group)
        batch_lay.setSpacing(4)

        # Timepoint range
        self.check_process_all = QCheckBox("Process ALL timepoints")
        self.check_process_all.setChecked(True)
        batch_lay.addWidget(self.check_process_all)

        tp_row = QHBoxLayout()
        tp_row.addWidget(QLabel("From t:"))
        self.spin_t_start = QSpinBox()
        self.spin_t_start.setRange(0, 99999)
        self.spin_t_start.setValue(0)
        self.spin_t_start.setMaximumWidth(70)
        tp_row.addWidget(self.spin_t_start)
        tp_row.addWidget(QLabel("to t:"))
        self.spin_t_end = QSpinBox()
        self.spin_t_end.setRange(0, 99999)
        self.spin_t_end.setValue(100)
        self.spin_t_end.setMaximumWidth(70)
        tp_row.addWidget(self.spin_t_end)
        tp_row.addStretch()
        self.tp_widget = QWidget()
        self.tp_widget.setLayout(tp_row)
        batch_lay.addWidget(self.tp_widget)

        def _toggle_tp(state):
            self.tp_widget.setVisible(not self.check_process_all.isChecked())
        self.check_process_all.stateChanged.connect(_toggle_tp)
        _toggle_tp(None)

        # Workers
        n_cores = os.cpu_count() or 4
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, max(1, n_cores - 1))
        self.spin_workers.setValue(1)
        self.spin_workers.setMaximumWidth(60)
        w_form = QFormLayout()
        w_form.setContentsMargins(0, 0, 0, 0)
        w_form.addRow("Workers:", self.spin_workers)
        batch_lay.addLayout(w_form)

        # Overwrite checkbox
        self.check_overwrite = QCheckBox("Overwrite existing results")
        self.check_overwrite.setChecked(False)
        batch_lay.addWidget(self.check_overwrite)

        # Run button
        self.btn_run_segmentation = QPushButton("▶ Run APOC Batch Segmentation")
        self.btn_run_segmentation.setStyleSheet(
            "QPushButton { background: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 13px; } "
            "QPushButton:hover { background: #0069d9; }"
        )
        self.btn_run_segmentation.clicked.connect(lambda _: self._on_run_segmentation(interactive=True, block=False))

        # Add-to-queue button
        self.btn_queue_apoc_segment = QPushButton("+🛒")
        self.btn_queue_apoc_segment.setFixedSize(36, 36)
        self.btn_queue_apoc_segment.setToolTip("Add APOC Segmentation to Queue")
        self.btn_queue_apoc_segment.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_btn_row = QHBoxLayout()
        batch_btn_row.setSpacing(4)
        batch_btn_row.addWidget(self.btn_run_segmentation, stretch=1)
        batch_btn_row.addWidget(self.btn_queue_apoc_segment)
        batch_lay.addLayout(batch_btn_row)

        layout.addWidget(batch_group)

        layout.addStretch()
        scroll.setWidget(content)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)

        # ── Auto-save bindings ─────────────────────────────────────
        self.spin_examples.valueChanged.connect(lambda _: self._save_apoc_params_to_yaml())
        self.spin_workers.valueChanged.connect(lambda _: self._save_apoc_params_to_yaml())
        self.combo_gpu_device.currentTextChanged.connect(self._on_gpu_device_changed)
        self.btn_force_cpu.toggled.connect(self._on_force_cpu_toggled)
        self.check_process_all.stateChanged.connect(lambda _: self._save_apoc_params_to_yaml())
        self.spin_t_start.valueChanged.connect(lambda _: self._save_apoc_params_to_yaml())
        self.spin_t_end.valueChanged.connect(lambda _: self._save_apoc_params_to_yaml())

        # Apply saved device immediately so training/preview runs on the intended GPU.
        self._apply_apoc_gpu_selection(log_message=False)

    def _discover_apoc_gpu_devices(self):
        try:
            import pyclesperanto_prototype as cle
            names = list(cle.available_device_names() or [])
        except Exception:
            names = []

        # Keep order stable and remove duplicates.
        seen = set()
        unique = []
        for name in names:
            key = str(name)
            if key in seen:
                continue
            seen.add(key)
            unique.append(key)
        return unique

    def _on_gpu_device_changed(self, _text):
        self._save_apoc_params_to_yaml()
        self._apply_apoc_gpu_selection(log_message=True)

    def _on_force_cpu_toggled(self, checked):
        self.combo_gpu_device.setEnabled(not checked)
        self._save_apoc_params_to_yaml()
        self._apply_apoc_gpu_selection(log_message=True)

    def _selected_gpu_device_name(self):
        if not hasattr(self, 'combo_gpu_device'):
            return ""
        if not self.combo_gpu_device.isEnabled():
            return ""
        return str(self.combo_gpu_device.currentText() or "").strip()

    def _apply_apoc_gpu_selection(self, log_message=True):
        if hasattr(self, 'btn_force_cpu') and self.btn_force_cpu.isChecked():
            try:
                import pyclesperanto_prototype as cle
                import warnings
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    device = cle.select_device("CPU")
                    
                    device_name = device.name if hasattr(device, 'name') else str(device)
                    warning_msg = str(w[-1].message) if len(w) > 0 else ""
                    
                    if log_message:
                        if "CPU" not in device_name and "No OpenCL device found" in warning_msg:
                            self.log(f"⚠️ OpenCL CPU driver missing. Fell back to: {device_name}")
                            self.log("To run APOC on the CPU, install an OpenCL CPU runtime (e.g. Intel CPU Runtime for OpenCL).")
                        else:
                            self.log(f"APOC: Forced OpenCL device '{device_name}'")
                return True
            except Exception as e:
                if log_message:
                    self.log(f"⚠️ Could not force CPU device: {e}")
                return False

        gpu_device = self._selected_gpu_device_name()
        if not gpu_device:
            return False
        try:
            import pyclesperanto_prototype as cle
            cle.select_device(gpu_device)
            if log_message:
                self.log(f"APOC: using OpenCL device '{gpu_device}'")
            return True
        except Exception as e:
            if log_message:
                self.log(f"⚠️ Could not select OpenCL device '{gpu_device}': {e}")
            return False

    def _update_training_controls_state(self):
        """Lock/unlock APOC training controls based on whether training data is loaded.

        Delegates to the training widget's native
        :meth:`APOCTrainingWidget.set_training_enabled`. Per-tab tuning
        widgets (feature presets, RF params, etc.) are also toggled here
        because the training widget only flips the high-impact controls.
        """
        if self._training_widget is None:
            return

        can_train = bool(self._is_session_active)
        tw = self._training_widget
        tw.set_training_enabled(can_train)

        # Lock the tuning controls that are not part of the training widget's
        # public "training-enabled" set so users can't tweak presets/RF params
        # before data is loaded but can still tune post-processing spinners.
        for tab in getattr(tw, "tabs", {}).values():
            for name in (
                "feature_combo",
                "tune_group",
                "consider_original_cb",
                "stats_btn",
                "max_depth_spin",
                "num_ensembles_spin",
            ):
                w = getattr(tab, name, None)
                if w is not None:
                    w.setEnabled(can_train)

        # Re-apply import lock for tabs whose features were locked by a previous import.
        # The enable loop above would otherwise undo locks set by apply_import().
        if can_train:
            for tab in getattr(tw, "tabs", {}).values():
                fs = getattr(tab, "_locked_feature_string", None)
                if fs:
                    tab.lock_features(fs)

        # Show/hide training parameters and import button based on session state
        if hasattr(tw, "set_classifier_params_visible"):
            tw.set_classifier_params_visible(can_train)
        else:
            tw.setVisible(can_train)
        self.btn_import_training.setVisible(can_train)
        if hasattr(self, '_btn_import_help'):
            self._btn_import_help.setVisible(can_train)

        if can_train:
            self.log("APOC training controls unlocked (training data loaded).")
        else:
            self.log("APOC training controls locked. Load training data to enable classifier training.")

    def _channel_index_from_name(self, channel_name):
        text = str(channel_name or "")
        import re
        match = re.search(r"\((\d+)\)", text)
        if match:
            return int(match.group(1))
        match = re.search(r"(?:^|\b)(?:channel|ch)\s*(\d+)(?:\b|$)", text, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
        nums = re.findall(r"\d+", text)
        if nums:
            return int(nums[-1])
        return None

    def _restore_tab_channels_fallback(self, tab, configured_channels):
        configured = [str(name) for name in (configured_channels or []) if str(name).strip()]
        exact = set(configured)
        configured_idx = {
            idx for idx in (self._channel_index_from_name(name) for name in configured)
            if idx is not None
        }

        any_checked = False
        for cb in getattr(tab, "channel_checkboxes", []):
            checked = cb.text() in exact
            if not checked and configured_idx:
                current_idx = self._channel_index_from_name(cb.text())
                checked = current_idx in configured_idx if current_idx is not None else False
            cb.setChecked(checked)
            any_checked = any_checked or checked

        if hasattr(tab, "_default_channel_names"):
            if any_checked:
                tab._default_channel_names = [cb.text() for cb in getattr(tab, "channel_checkboxes", []) if cb.isChecked()]
            else:
                tab._default_channel_names = list(configured)

    def _extract_tab_config(self, apoc_params, ct):
        prefix = f"apoc_{ct}_"
        cfg = {}
        for key, value in (apoc_params or {}).items():
            if key.startswith(prefix):
                cfg[key[len(prefix):]] = value
        return cfg

    def _on_training_channels_refreshed(self, ct):
        """Slot: APOCTrainingWidget rebuilt the channel checkboxes for *ct*.

        Restores saved channel selections from the persisted tab config and
        re-wires the organoid-mirroring signal handlers. Used instead of
        monkey-patching ``tab.refresh_channel_checkboxes``.
        """
        if self._training_widget is None:
            return
        tab = self._training_widget.tabs.get(ct)
        if tab is None:
            return
        # Keep channel inputs locked until training data is loaded.
        can_train = bool(self._is_session_active)
        for cb in getattr(tab, "channel_checkboxes", []):
            cb.setEnabled(can_train)
        cfg = getattr(tab, "_saved_config_for_restore", None)
        if not (isinstance(cfg, dict) and cfg.get("channels")):
            # Stash is empty or has no channel selection — read from live YAML.
            live_params = self.metadata_loader.behav3d_parameters.get("apoc", {})
            cfg = self._extract_tab_config(live_params, ct) or {}
        if isinstance(cfg, dict) and cfg.get("channels"):
            self._restore_tab_channels_fallback(tab, cfg.get("channels", []))
        # Channel checkboxes are recreated by refresh; re-wire sync handlers.
        self._wire_organoid_tab_sync_signals(ct, tab, include_channels=True, include_rf=False)
        self._wire_non_organoid_tab_channel_signals(ct, tab)

    def _collect_apoc_tab_config(self):
        """Collect flattened per-cell-type APOC config from training tabs."""
        apoc_config = {}
        tw = self._training_widget
        if tw is None:
            return apoc_config

        for ct, tab in tw.tabs.items():
            cfg = tab.get_config()
            if not tab.channel_selection_is_complete():
                # Channel layers are still loading in — don't overwrite the
                # saved selection with this partial checkbox state.
                cfg.pop("channels", None)
            for k, v in cfg.items():
                apoc_config[f"apoc_{ct}_{k}"] = v
        return apoc_config

    def _on_organoid_tab_input_changed(self, source_ct, *_args):
        """Persist tab params and mirror organoid settings from the edited source tab."""
        if self._syncing_organoids:
            return
        if source_ct not in self._organoid_types:
            return

        apoc_config = self._collect_apoc_tab_config()
        self._save_apoc_params_to_yaml(
            updated_apoc_params=apoc_config,
            sync_source_ct=source_ct,
        )

    def _wire_organoid_tab_sync_signals(self, ct, tab, include_channels=True, include_rf=True):
        """Wire organoid tab widgets so edits trigger immediate source-aware mirroring."""
        if ct not in self._organoid_types:
            return

        def _connect_once(widget, signal_name):
            if widget is None:
                return
            prop_name = f"_behav3d_sync_connected_{signal_name}"
            try:
                if bool(widget.property(prop_name)):
                    return
                signal = getattr(widget, signal_name, None)
                if signal is None:
                    return
                signal.connect(lambda *_args, _ct=ct: self._on_organoid_tab_input_changed(_ct))
                widget.setProperty(prop_name, True)
            except Exception:
                return

        if include_rf:
            _connect_once(getattr(tab, "max_depth_spin", None), "valueChanged")
            _connect_once(getattr(tab, "num_ensembles_spin", None), "valueChanged")

        if include_channels:
            for cb in getattr(tab, "channel_checkboxes", []):
                _connect_once(cb, "stateChanged")

    def _wire_non_organoid_tab_channel_signals(self, ct, tab):
        """Wire channel checkboxes for non-organoid tabs to immediately persist
        the selection to YAML.  Organoid tabs already do this via
        ``_wire_organoid_tab_sync_signals``; this covers immune / other tabs."""
        if ct in self._organoid_types:
            return  # already handled by organoid sync

        def _connect_once(widget, signal_name):
            if widget is None:
                return
            prop_name = f"_behav3d_persist_connected_{signal_name}"
            try:
                if bool(widget.property(prop_name)):
                    return
                signal = getattr(widget, signal_name, None)
                if signal is None:
                    return
                signal.connect(lambda *_: self._on_non_organoid_channel_changed())
                widget.setProperty(prop_name, True)
            except Exception:
                pass

        for cb in getattr(tab, "channel_checkboxes", []):
            _connect_once(cb, "stateChanged")

    def _on_non_organoid_channel_changed(self):
        """Persist all tab configs to YAML when a non-organoid channel changes."""
        apoc_config = self._collect_apoc_tab_config()
        self._save_apoc_params_to_yaml(updated_apoc_params=apoc_config)

    def _restore_training_tab_configs(self, apoc_params):
        if self._training_widget is None:
            return

        for ct, tab in self._training_widget.tabs.items():
            cfg = self._extract_tab_config(apoc_params, ct)
            if not cfg:
                continue
            # Stash so :meth:`_on_training_channels_refreshed` can re-apply
            # channel selections every time the checkboxes get rebuilt.
            tab._saved_config_for_restore = dict(cfg)
            try:
                tab.apply_config(cfg)
            except Exception as e:
                self.log(f"⚠️ Could not restore APOC tab config for '{ct}': {e}")
            if cfg.get("channels"):
                self._restore_tab_channels_fallback(tab, cfg.get("channels", []))

    def _reorder_apoc_training_layers(self):
        """Enforce grouped training layer order: channels, labels, probabilities, segments."""
        try:
            layer_names = [layer.name for layer in self.viewer.layers]

            channel_names = [name for name in layer_names if name.startswith("Channel ")]
            label_names = [name for name in layer_names if name.startswith("User Provided Labels")]
            prob_names = [name for name in layer_names if name.startswith("Probability Map")]
            seg_names = [
                name for name in layer_names
                if name.endswith(" Segments") or name == "Pixel Classification (Dead)"
            ]

            grouped = channel_names + label_names + prob_names + seg_names
            grouped_unique = []
            seen = set()
            for name in grouped:
                if name in seen:
                    continue
                seen.add(name)
                grouped_unique.append(name)

            trailing = [name for name in layer_names if name not in seen]
            final_order = grouped_unique + trailing

            for target_idx, name in enumerate(final_order):
                if name not in [layer.name for layer in self.viewer.layers]:
                    continue
                current_idx = self.viewer.layers.index(name)
                if current_idx != target_idx:
                    self.viewer.layers.move(current_idx, target_idx)
        except Exception as e:
            self.log(f"⚠️ Could not reorder APOC training layers: {e}")

    def _prompt_visualize_after_apoc_segmentation(self):
        """Prompt user to switch to Visualization tab and load first sample.

        Delegates the actual tab switch to the host-supplied callback
        (``switch_to_visualization_callback``) so we don't walk the parent
        chain at runtime.
        """
        res = QMessageBox.question(
            self,
            "Segmentation Finished",
            "Batch segmentation finished successfully! \n\nDo you want to switch to the Visualization Tab and see the segments?",
            QMessageBox.Yes | QMessageBox.No,
        )

        if res != QMessageBox.Yes:
            return

        if callable(self._switch_to_visualization):
            self._switch_to_visualization()
        else:
            self.log("  Visualization tab callback unavailable.")

    # ── Load Training Data ──────────────────────────────────────
    def _on_load_training_clicked(self, interactive=True):
        try:
            if self.metadata_loader.metadata is None:
                self.log("⚠️ Cannot generate training data: No metadata loaded.")
                return

            self._apply_apoc_gpu_selection(log_message=False)
            gpu_name = self._selected_gpu_device_name() or "default OpenCL device"
            self.log(f"APOC training pipeline using GPU device: {gpu_name}")

            self.log("Loading training data...")
            
            # --- Save params if needed ---
            pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
            pc["examples_per_sample"] = int(self.spin_examples.value())
            
            output_dir = Path(self.metadata_loader.output_dir)
            if output_dir:
                params_path = Path(output_dir) / "behav3d_parameters.yml"
                try:
                    with open(params_path, "w") as f:
                        yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
                except Exception as e:
                    self.log(f"Warning: Could not save parameters: {e}")

            md = self.metadata_loader.metadata
            from behav3d.core.metadata import (
                detect_organoid_types_from_metadata,
                detect_immune_cell_types_from_metadata,
                detect_other_cell_types_from_metadata,
                is_combined_multicolor_celltype,
            )
            organoid_types = [ct for ct in detect_organoid_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
            immune_types = [ct for ct in detect_immune_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
            other_types = [ct for ct in detect_other_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
            
            # Use apoc_train helper to fetch and process images
            from behav3d.preprocessing.segmentation.apoc_train import (
                _load_training_images, _predicted_labels_path, _probability_map_path,
            )
            from behav3d.preprocessing import zeropad_image_to_match_shape
            
            # Determine the cached data path and whether it already exists
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            image_outpath = pixel_class_outdir / 'PixelClassifier_Images.zarr'
            saved_examples = pc.get("examples_per_sample", None)

            # Default: generate if no data exists; prompt only if data already exists
            load_existing = False
            if image_outpath.exists():
                if interactive:
                    # Read the actual number of stored training images from the zarr header
                    # APOC saves as (C, T, Z, Y, X) → shape[1] = image count
                    try:
                        import zarr as _zarr
                        _cached = _zarr.open(str(image_outpath), mode="r")
                        saved_image_count = _cached.shape[1] if len(_cached.shape) > 1 else _cached.shape[0]
                        del _cached
                    except Exception:
                        saved_image_count = None

                    n_samples = len(md['sample_name'].unique())
                    examples_per_sample = self.spin_examples.value()
                    total_new_images = n_samples * examples_per_sample
                    if saved_image_count is not None:
                        msg = (
                            f"Training data already exists with {saved_image_count} images.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            f"Overwrite with new training data, or load the existing data?"
                        )
                    else:
                        msg = (
                            f"Training data already exists.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            f"Overwrite with new training data, or load the existing data?"
                        )

                    from qtpy.QtWidgets import QMessageBox
                    box = QMessageBox(self)
                    box.setWindowTitle("Training Data Detected")
                    box.setText(msg)
                    btn_generate = box.addButton("Generate New", QMessageBox.AcceptRole)
                    btn_load = box.addButton("Load Existing", QMessageBox.YesRole)
                    btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                    box.exec_()
                    if box.clickedButton() == btn_cancel:
                        self.log("Action cancelled.")
                        return
                    elif box.clickedButton() == btn_load:
                        load_existing = True
                    # else: Generate New clicked → load_existing stays False
                else:
                    # Non-interactive: load existing to avoid discarding data
                    load_existing = True
            # else: no data exists → generate immediately without any dialog

            # Flush current channel selections and tab config to YAML BEFORE
            # clearing layers, so the subsequent checkbox rebuild can restore
            # the live UI state rather than a stale one-time stash.
            current_tab_cfg = self._collect_apoc_tab_config()
            existing_apoc = self.metadata_loader.behav3d_parameters.get("apoc", {})
            for key in list(current_tab_cfg.keys()):
                if key.endswith("_channels") and not current_tab_cfg[key]:
                    existing_val = existing_apoc.get(key)
                    if existing_val:
                        current_tab_cfg[key] = existing_val
            self._save_apoc_params_to_yaml(updated_apoc_params=current_tab_cfg)
            # Refresh the per-tab stash from the freshly written YAML so that
            # _on_training_channels_refreshed re-applies the correct channels.
            self._restore_training_tab_configs(
                self.metadata_loader.behav3d_parameters.get("apoc", {})
            )

            # Clear all viewer layers entirely
            self.viewer.layers.clear()

            # Load images
            self.log("Running _load_training_images...")
            image_list, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
                metadata=md,
                output_dir=str(output_dir),
                examples_per_sample=self.spin_examples.value(),
                organoid_types=organoid_types,
                immune_types=immune_types,
                other_types=other_types,
                overwrite_images=not load_existing,
            )
            
            if not image_list:
                self.log("⚠️ No training images found!")
                return

            self.log("Padding and stacking images...")
            max_shape = list(image_list[0].shape)
            for img in image_list[1:]:
                for i in range(len(max_shape)):
                    max_shape[i] = max(max_shape[i], img.shape[i])
            image_list = [zeropad_image_to_match_shape(img, max_shape) for img in image_list]
            stacked = np.stack(image_list, axis=0) # (T_total, C, Z, Y, X)

            T_total = stacked.shape[0]
            n_channels = stacked.shape[1]
            channel_colors = [
                "cyan", "yellow", "red", "green", "magenta", "blue",
                "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
            ]

            # Add Image layers. Channels are added one at a time, and each
            # `add_image` call fires its own `layers.events.inserted` event —
            # pause the training widget's own listener for the duration of
            # the loop and refresh the training tabs' channel checkboxes
            # once afterwards, so per-tab channel selections aren't rebuilt
            # against a still-partial set of channel layers (which would
            # look like no channels matching a tab's saved selection).
            # Note: this must only block the training widget's own callback
            # (via `pause_channel_refresh`), not the whole event — napari's
            # internal layer-controls widget registration listens on the
            # same event, and blocking it desyncs the viewer UI.
            def _add_channel_layers():
                for ch in range(n_channels):
                    channel_data = stacked[:, ch, :, :, :]
                    nonzero = channel_data[channel_data > 0]
                    clim = (0, float(np.percentile(nonzero, 99.8))) if nonzero.size > 0 else (0, 1e-3)
                    img_layer = self.viewer.add_image(
                        channel_data,
                        name=f"Channel {ch}",
                        contrast_limits=clim,
                        colormap=channel_colors[ch % len(channel_colors)],
                        blending="additive",
                        opacity=0.8,
                    )
                    img_layer.contrast_limits_range = (0, max(float(channel_data.max()), 1e-3))

            if self._training_widget is not None:
                with self._training_widget.pause_channel_refresh():
                    _add_channel_layers()
                self._training_widget._refresh_all_channels()
            else:
                _add_channel_layers()

            # Add Label layers
            label_shape = (T_total,) + stacked.shape[2:] 

            for cell_type in all_cell_types:
                saved_path = Path(pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
                if saved_path.exists():
                    existing = np.asarray(load_zarr(saved_path))
                    if existing.shape == label_shape:
                        user_labels = existing
                        self.log(f"  ↩ Restored saved labels for '{cell_type}'")
                    else:
                        self.log(f"  ⚠️ Saved labels shape {existing.shape} ≠ expected {label_shape} — starting fresh")
                        user_labels = np.zeros(label_shape, dtype=np.int16)
                else:
                    user_labels = np.zeros(label_shape, dtype=np.int16)

                # Use original pixel classifier setup function for consistency
                new_layer = self.viewer.add_labels(
                    user_labels,
                    name=f"User Provided Labels ({cell_type.capitalize()})",
                    opacity=0.5,
                )
                self._configure_user_label_layer(new_layer)

            if has_death:
                dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
                if dead_path.exists():
                    dead_labels = np.asarray(load_zarr(dead_path))
                    if dead_labels.shape == label_shape:
                        self.log(f"  ↩ Restored saved labels for 'dead'")
                    else:
                        dead_labels = np.zeros(label_shape, dtype=np.int16)
                else:
                    dead_labels = np.zeros(label_shape, dtype=np.int16)
                    
                dead_layer = self.viewer.add_labels(dead_labels, name="User Provided Labels (Dead)", opacity=0.5)
                self._configure_user_label_layer(dead_layer)

            # Restore previously generated probability maps and segmentations
            for cell_type in all_cell_types:
                seg_p = _predicted_labels_path(pixel_class_outdir, cell_type)
                if seg_p and seg_p.exists():
                    try:
                        pred = np.asarray(load_zarr(seg_p))
                        if pred.shape == label_shape:
                            self.viewer.add_labels(
                                pred,
                                name=f"{cell_type.capitalize()} Segments",
                                opacity=0.8, visible=False,
                            )
                            self.log(f"  ↩ Restored predicted labels for '{cell_type}'")
                    except Exception:
                        pass
                prob_p = _probability_map_path(pixel_class_outdir, cell_type)
                if prob_p and prob_p.exists():
                    try:
                        prob = np.asarray(load_zarr(prob_p))
                        if prob.shape == label_shape:
                            self.viewer.add_image(
                                prob,
                                name=f"Probability Map ({cell_type.capitalize()})",
                                opacity=0.6, blending="additive", colormap="magma",
                                contrast_limits=(0.0, 1.0), visible=False,
                            )
                            self.log(f"  ↩ Restored probability map for '{cell_type}'")
                    except Exception:
                        pass

            if has_death:
                seg_p = _predicted_labels_path(pixel_class_outdir, "dead")
                if seg_p and seg_p.exists():
                    try:
                        pred = np.asarray(load_zarr(seg_p))
                        if pred.shape == label_shape:
                            dead_pred_layer = self.viewer.add_labels(
                                pred, name="Pixel Classification (Dead)",
                                opacity=0.8, visible=False,
                            )
                            _set_dead_mask_layer_color(dead_pred_layer)
                            self.log("  ↩ Restored predicted labels for 'dead'")
                    except Exception:
                        pass
                prob_p = _probability_map_path(pixel_class_outdir, "dead")
                if prob_p and prob_p.exists():
                    try:
                        prob = np.asarray(load_zarr(prob_p))
                        if prob.shape == label_shape:
                            self.viewer.add_image(
                                prob, name="Probability Map (Dead)",
                                opacity=0.6, blending="additive", colormap="magma",
                                contrast_limits=(0.0, 1.0), visible=False,
                            )
                            self.log("  ↩ Restored probability map for 'dead'")
                    except Exception:
                        pass

            self._reorder_apoc_training_layers()

            self.log("✅ Training data generated/loaded in viewer!")
            self._is_session_active = True
            self._update_training_controls_state()
            
            # Update and show legend
            while self.legend_container.layout().count():
                child = self.legend_container.layout().takeAt(0)
                if child.widget():
                    child.widget().deleteLater()
            legend = PerClassLegendWidget(self.viewer, all_cell_types, has_death)
            self.legend_container.layout().addWidget(legend)
            self.legend_container.setVisible(True)

        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Error during training data generation: {e}")

    def _configure_user_label_layer(self, layer):
        """Configure a user label layer for simplified APOC usage."""
        import numpy as np
        layer.blending = "additive"
        layer.color = {0: (0, 0, 0, 0), 1: (1.0, 0.2, 0.2, 1), 2: (0.0, 1.0, 1.0, 1)}
        layer.selected_label = 1
        def _clamp_label(event):
            if layer.selected_label > 2:
                layer.selected_label = 2
        layer.events.selected_label.connect(_clamp_label)

    def cleanup_session(self):
        """Remove all APOC layers from the viewer and reset session state.

        Removes input channels, user-painted labels, probability maps and
        instance segmentations (everything created during a training/preview
        session) so the viewer is clean when the user switches tab or
        method.  Mirrors PixelClassifierWidget.reset_ui() and
        ConvPaintWidget.cleanup_session().
        """
        n_removed = _remove_segmentation_layers(self.viewer)
        self._is_session_active = False
        self.log(f"APOC training session cleared ({n_removed} layers removed).")
        self._update_training_controls_state()
        
        # Hide legend
        while self.legend_container.layout().count():
            child = self.legend_container.layout().takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.legend_container.setVisible(False)

    def _save_apoc_user_labels(self):
        """Save APOC user labels to disk via the training widget."""
        if self._training_widget is not None:
            self._training_widget.save_user_labels(log=self.log)

    # ── Strategy changed ───────────────────────────────────────────
    def _on_training_widget_strategy_changed(self, strategy):
        """Persist the new APOC strategy when the training widget signals a change."""
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["apoc_strategy"] = str(strategy)
        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
            except Exception:
                pass

    # ── Training widget signal slots ───────────────────────────────
    def _on_training_started(self, cell_types_to_train):
        """Log + GPU select when the training widget starts a run."""
        self._apply_apoc_gpu_selection(log_message=False)
        gpu_name = self._selected_gpu_device_name() or "default OpenCL device"
        tw = self._training_widget
        strategy_name = (
            tw.strategy_combo.currentText()
            if tw is not None and tw.strategy_combo is not None
            else ""
        )
        self.log(f"\U0001f5a5\ufe0f APOC training device: {gpu_name}")
        self.log(f"\u25b6 Starting APOC training for: {cell_types_to_train}")
        self._apoc_cli_info(
            f"training start | strategy={strategy_name} | gpu={gpu_name} | cell_types={cell_types_to_train}"
        )

    def _on_training_finished(self, cell_types_to_train, successes):
        """Reorder layers + normalize classifier headers after training."""
        self._reorder_apoc_training_layers()
        try:
            pixel_class_outdir = Path(self.metadata_loader.output_dir) / "images" / "PixelClassification"
            for ct in cell_types_to_train:
                ok, msg = self._apoc_try_migrate_from_tab_config(pixel_class_outdir, ct)
                if not ok:
                    self.log(f"  \u26a0\ufe0f {ct}: could not normalize classifier header after training ({msg})")
                    self._apoc_cli_info(f"warning | {ct}: could not normalize classifier header ({msg})")
        except Exception as e:
            self.log(f"  \u26a0\ufe0f Could not normalize classifier headers after training: {e}")
            self._apoc_cli_info(f"warning | classifier-header normalization failed: {e}")
        self.log("\u2705 APOC training process completed.")
        self._apoc_cli_info("training complete")

    def _on_instance_preview_started(self, ct):
        """Clear stale preview layers + select GPU when the preview starts."""
        self._apply_apoc_gpu_selection(log_message=False)
        tw = self._training_widget
        strategy = tw._resolve_strategy(ct) if tw is not None else ""
        self._clear_apoc_preview_layers_for_cell_type(ct)
        if strategy == "APOC (Direct Instance Segmentation)":
            self._apoc_cli_info(
                f"instance preview skipped (direct strategy) | stale layers cleared | cell_type={ct}"
            )
        self.log(f"\U0001f50d Running instance segmentation preview for '{ct}' (strategy: {strategy})...")
        self._apoc_cli_info(f"instance preview start | cell_type={ct} | strategy={strategy}")

    def _on_instance_preview_finished(self, ct):
        """Reorder layers after a preview completes."""
        self._reorder_apoc_training_layers()
        self._apoc_cli_info(f"instance preview complete | cell_type={ct}")

    def _clear_apoc_preview_layers_for_cell_type(self, ct: str):
        """Remove stale APOC preview output layers for one cell type (GUI-side only)."""
        names = []
        if str(ct).strip().lower() == "dead":
            names = ["Pixel Classification (Dead)", "Probability Map (Dead)"]
        else:
            ctc = str(ct).capitalize()
            names = [f"{ctc} Segments", f"Probability Map ({ctc})"]

        for lname in names:
            if lname in [layer.name for layer in self.viewer.layers]:
                try:
                    self.viewer.layers.remove(self.viewer.layers[lname])
                except Exception:
                    pass

    def _apoc_cli_info(self, message: str):
        """Emit concise APOC progress info to stdout (best effort)."""
        try:
            print(f"APOC | {message}", flush=True)
        except Exception:
            pass

    # ── Metadata update ─────────────────────────────────────────
    def _on_metadata_updated(self):
        """Called when metadata is loaded/updated. Rebuild training widget."""
        md = self.metadata_loader.metadata
        if md is None:
            return

        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
            is_combined_multicolor_celltype,
        )

        org = [ct for ct in detect_organoid_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        imm = [ct for ct in detect_immune_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        oth = [ct for ct in detect_other_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)]
        all_types = org + imm + oth
        has_death = has_dead_channel(md)
        self._organoid_types = list(org)   # keep for unified-organoid checkbox
        self.all_cell_types = list(all_types)

        if not all_types:
            return

        # Remove the previous training widget cleanly: ``cleanup()`` is the
        # native disconnect method exposed by APOCTrainingWidget so napari
        # layer event handlers never outlive this object.
        if self._training_widget is not None:
            try:
                self._training_widget.cleanup()
            except Exception:
                pass
            self.training_layout.removeWidget(self._training_widget)
            self._training_widget.hide()
            self._training_widget.deleteLater()
            self._training_widget = None

        if self.training_placeholder.parent() is not None:
            self.training_layout.removeWidget(self.training_placeholder)
            self.training_placeholder.hide()

        # Create and embed APOCTrainingWidget.
        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget
        params = self.metadata_loader.behav3d_parameters
        apoc_params = dict(params.get("apoc", {}) or {})
        pc = params.get("pixel_classifier", {}) or {}

        all_strats = self.all_strategies()
        current_strategy = pc.get("apoc_strategy", "APOC Probability Map + Watershed")
        if current_strategy not in all_strats:
            current_strategy = all_strats[0]
        apoc_params["apoc_strategy"] = current_strategy

        # Optional "Unify all organoids" toolbar checkbox shown above the tab
        # widget. Built first so we can pass it via ``extra_toolbar_widgets``.
        toolbar_widgets = []
        self._check_unify_organoids = None
        if len(org) >= 2:
            unify_row = QHBoxLayout()
            self._check_unify_organoids = QCheckBox(
                "\u26d3 Apply same settings to all organoids"
            )
            self._check_unify_organoids.setChecked(self._unify_organoids)
            self._check_unify_organoids.setToolTip(
                "When checked, clicking will immediately mirror the first organoid tab's\n"
                "config (features, strategy, params) to all other organoid tabs.\n"
                "Toggle again to re-sync after making further changes."
            )
            self._check_unify_organoids.stateChanged.connect(self._on_unify_organoids_changed)
            unify_row.addWidget(self._check_unify_organoids)
            unify_row.addWidget(HelpButton(
                "Apply Same Settings to All Organoids",
                "When there are multiple organoid cell types, checking this box "
                "immediately copies the first organoid tab's configuration "
                "(selected channels/features, strategy, and EDT/opening/min "
                "size/fill-holes parameters) to all other organoid tabs.\n\n"
                "It only syncs once at the moment you check it \u2014 further edits "
                "to any organoid tab are NOT automatically propagated. Toggle "
                "the checkbox off and back on to re-sync after making more "
                "changes.\n\n"
                "Immune and other non-organoid cell types are not affected."
            ))
            unify_widget = QWidget()
            unify_widget.setLayout(unify_row)
            toolbar_widgets.append(unify_widget)

        _pc_params = params.get("pixel_classifier", {}) or {}
        _md = self.metadata_loader.metadata
        if _md is not None and "pixel_distance_xy" in _md.columns and "distance_unit" in _md.columns:
            _unit = str(_md["distance_unit"].iloc[0])
            _xy_from_md = convert_distance(float(_md["pixel_distance_xy"].iloc[0]), _unit)
            _z_from_md  = convert_distance(float(_md["pixel_distance_z"].iloc[0]),  _unit)
        else:
            _xy_from_md = None
            _z_from_md  = None
        _pixel_sizes = {
            "xy_um": _pc_params.get("pixel_size_xy") or _pc_params.get("pixel_size_xy_um") or _xy_from_md,
            "z_um":  _pc_params.get("pixel_size_z")  or _pc_params.get("pixel_size_z_um")  or _z_from_md,
        }

        self._training_widget = APOCTrainingWidget(
            viewer=self.viewer,
            pixel_class_outdir=str(pixel_class_outdir),
            all_cell_types=all_types,
            has_death=has_death,
            initial_params=apoc_params,
            on_params_changed=self._save_apoc_params_to_yaml,
            instance_controls_mode="inline",
            per_cell_type_strategy=True,
            extra_toolbar_widgets=toolbar_widgets,
            pixel_sizes=_pixel_sizes,
        )
        # Pre-select the saved global strategy via the widget's combo (this
        # also propagates per-tab visibility for ``Advanced`` mode).
        if self._training_widget.strategy_combo is not None:
            self._training_widget.strategy_combo.setCurrentText(current_strategy)
        # Redirect backend internal logs to our GUI log.
        self._training_widget.set_log_fn(self.log)
        self._restore_training_tab_configs(apoc_params)

        # Wire signals: channel rebuilds restore saved channels & re-wire
        # organoid mirroring; training / preview hooks add logging, GPU
        # selection and layer reordering without monkey-patching backend
        # methods.
        tw = self._training_widget
        tw.channels_refreshed.connect(self._on_training_channels_refreshed)
        tw.training_started.connect(self._on_training_started)
        tw.training_finished.connect(self._on_training_finished)
        tw.instance_preview_started.connect(self._on_instance_preview_started)
        tw.instance_preview_finished.connect(self._on_instance_preview_finished)
        tw.strategy_changed.connect(self._on_training_widget_strategy_changed)
        tw.training_finished.connect(self._on_training_finished_update_export_btn)
        tw.import_applied.connect(self._on_import_applied)
        tw.import_cleared.connect(self._on_import_cleared)

        # Add global "▶ Run [CellType] Segmentation" button below the
        # train buttons. Lives in the same APOC layout as the Train buttons.
        # Per-cell type time range limits for the global Run button
        ct_tp_row = QHBoxLayout()
        self.check_limit_timerange_ct = QCheckBox("Process all timepoints")
        self.check_limit_timerange_ct.setChecked(True)
        ct_tp_row.addWidget(self.check_limit_timerange_ct)
        
        self.spin_t_start_ct = QSpinBox()
        self.spin_t_start_ct.setRange(0, 9999)
        self.spin_t_start_ct.setValue(0)
        self.spin_t_start_ct.setMaximumWidth(70)
        self.spin_t_end_ct = QSpinBox()
        self.spin_t_end_ct.setRange(0, 9999)
        self.spin_t_end_ct.setValue(0)
        self.spin_t_end_ct.setMaximumWidth(70)

        ct_tp_range_row = QHBoxLayout()
        ct_tp_range_row.addWidget(QLabel("  From t:"))
        ct_tp_range_row.addWidget(self.spin_t_start_ct)
        ct_tp_range_row.addWidget(QLabel("to t:"))
        ct_tp_range_row.addWidget(self.spin_t_end_ct)
        ct_tp_range_row.addStretch()

        self.spin_t_start_ct.setVisible(False)
        self.spin_t_end_ct.setVisible(False)
        ct_tp_range_row_widget = QWidget()
        ct_tp_range_row_widget.setLayout(ct_tp_range_row)
        ct_tp_range_row_widget.setVisible(False)

        def _toggle_ct_tp(_state):
            ct_tp_range_row_widget.setVisible(not self.check_limit_timerange_ct.isChecked())
            self.spin_t_start_ct.setVisible(not self.check_limit_timerange_ct.isChecked())
            self.spin_t_end_ct.setVisible(not self.check_limit_timerange_ct.isChecked())

        self.check_limit_timerange_ct.stateChanged.connect(_toggle_ct_tp)

        tw._main_layout.addLayout(ct_tp_row)
        tw._main_layout.addWidget(ct_tp_range_row_widget)

        first_ct = tw._tab_cell_types[0] if tw._tab_cell_types else "?"
        self._global_run_instance_btn = QPushButton(f"▶ Run {first_ct.capitalize()} Segmentation")
        self._global_run_instance_btn.setStyleSheet(
            "background-color:#1976D2; color:white; font-weight:bold; padding:6px 12px;"
        )
        self._global_run_instance_btn.setToolTip(
            "Run instance segmentation preview for the currently active cell-type tab.\n"
            "Uses the strategy and parameters configured in that tab."
        )
        self._global_run_instance_btn.clicked.connect(self._on_global_run_instance)
        tw._main_layout.addWidget(self._global_run_instance_btn)
        tw.tab_widget.currentChanged.connect(self._update_global_run_btn_label)

        # Wire reactive organoid mirroring for channels + RF controls.
        for ct, tab in tw.tabs.items():
            self._wire_organoid_tab_sync_signals(ct, tab)
            self._wire_non_organoid_tab_channel_signals(ct, tab)

        self.training_layout.addWidget(self._training_widget)
        self._training_widget.setVisible(True)
        self._update_training_controls_state()

        # Enable import button now that cell types are known
        self.btn_import_training.setEnabled(True)

        # Re-apply persisted import if available
        self._restore_import_on_session_load()

        # Enable export bundle if training bundle already exists
        self._refresh_export_btn_state()

    # ── Import training data ────────────────────────────────────

    def _on_import_training_data_clicked(self):
        """Open file dialog, validate and apply imported training data."""
        from behav3d.preprocessing.segmentation.apoc_train import (
            _load_training_metadata,
            _load_training_data_for_celltype,
            _load_training_bundle,
        )

        if self._training_widget is None:
            QMessageBox.warning(self, "No Metadata", "Load metadata first before importing training data.")
            return

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select training data bundle or metadata from source experiment",
            "",
            "BEHAV3D Training Bundle (*.zip);;YAML files (training_metadata.yml);;All files (*)",
        )
        if not path:
            return

        if path.endswith(".zip"):
            try:
                meta, data_by_celltype = _load_training_bundle(path)
            except Exception as exc:
                QMessageBox.critical(self, "Import Error", f"Could not read training bundle from:\n{path}\n\n{exc}")
                return
            if meta is None:
                QMessageBox.critical(self, "Import Error", f"Could not read training metadata from:\n{path}")
                return
            imported_cts = list(meta.get("cell_types", []))
            if meta.get("has_death"):
                imported_cts.append("dead")
            missing = [ct for ct in imported_cts if ct not in data_by_celltype]
        else:
            meta = _load_training_metadata(path)
            if meta is None:
                QMessageBox.critical(self, "Import Error", f"Could not read training metadata from:\n{path}")
                return
            source_td = Path(path).parent
            imported_cts = list(meta.get("cell_types", []))
            if meta.get("has_death"):
                imported_cts.append("dead")
            data_by_celltype = {}
            missing = []
            for ct in imported_cts:
                X, y = _load_training_data_for_celltype(source_td, ct)
                if X is None:
                    missing.append(ct)
                else:
                    data_by_celltype[ct] = (X, y)

        if not data_by_celltype:
            QMessageBox.critical(
                self, "Import Error",
                "No training data arrays found in the selected file.\n"
                "Ensure the training bundle is complete."
            )
            return

        if missing:
            QMessageBox.warning(
                self, "Partial Import",
                f"Training data not found for: {', '.join(missing)}\n"
                "The remaining cell types will still be imported."
            )

        local_cts = list(self.all_cell_types)
        if self._training_widget.has_death:
            local_cts.append("dead")

        # Build cell type mapping (imported → local)
        cell_type_mapping = self._build_celltype_mapping(imported_cts, local_cts, meta)
        if cell_type_mapping is None:
            return  # user cancelled the remap dialog

        # Check pixel size compatibility
        if not self._check_pixel_size_compatibility(meta):
            return  # user cancelled

        # Show validation preview (pixel counts) before committing
        cell_type_mapping = self._show_import_validation_dialog(meta, data_by_celltype, cell_type_mapping)
        if cell_type_mapping is None:
            return  # user cancelled

        # Apply import
        self._training_widget.apply_import(meta, data_by_celltype, cell_type_mapping, source_path=path)

        # Persist import path to YAML
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["apoc_imported_training_path"] = str(path)
        self._flush_params_to_yaml()

    def _build_celltype_mapping(self, imported_cts, local_cts, meta):
        """Return {imported_ct: local_ct | None} mapping, or None if cancelled.

        When all names match exactly, returns the identity mapping immediately.
        Otherwise shows the remap dialog.
        """
        exact = {ct: ct for ct in imported_cts if ct in local_cts}
        unmatched = [ct for ct in imported_cts if ct not in local_cts]
        if not unmatched:
            return exact

        remap = self._show_celltype_remap_dialog(imported_cts, local_cts, exact)
        return remap

    def _show_celltype_remap_dialog(self, imported_cts, local_cts, initial_mapping):
        """Show a dialog letting the user map imported cell types to local ones.

        Returns the mapping dict {imported_ct: local_ct | None}, or None if cancelled.
        """
        dlg = QDialog(self)
        dlg.setWindowTitle("Cell Type Mapping")
        dlg.setMinimumWidth(420)
        lay = QVBoxLayout(dlg)

        info = QLabel(
            "Some cell type names in the imported data do not match the current experiment.\n"
            "Please map each imported type to a local type (or '(skip)' to exclude it)."
        )
        info.setWordWrap(True)
        lay.addWidget(info)

        combos = {}
        for imp_ct in imported_cts:
            row = QHBoxLayout()
            row.addWidget(QLabel(f"<b>{imp_ct}</b>  →"))
            combo = QComboBox()
            options = ["(skip)"] + local_cts
            combo.addItems(options)
            default = initial_mapping.get(imp_ct, None)
            if default and default in local_cts:
                combo.setCurrentText(default)
            elif imp_ct in local_cts:
                combo.setCurrentText(imp_ct)
            row.addWidget(combo, stretch=1)
            lay.addLayout(row)
            combos[imp_ct] = combo

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        lay.addWidget(btns)

        if dlg.exec_() != QDialog.Accepted:
            return None

        result = {}
        for imp_ct, combo in combos.items():
            sel = combo.currentText()
            result[imp_ct] = None if sel == "(skip)" else sel
        return result

    def _check_pixel_size_compatibility(self, imported_meta):
        """Warn if pixel sizes differ significantly. Returns True to proceed, False to cancel."""
        img_meta = imported_meta.get("image_metadata", {})
        imp_xy = img_meta.get("pixel_size_xy_um")
        imp_z  = img_meta.get("pixel_size_z_um")
        if imp_xy is None and imp_z is None:
            return True   # no pixel size info in imported data — skip check

        pc = self.metadata_loader.behav3d_parameters.get("pixel_classifier", {}) or {}
        cur_xy = pc.get("pixel_size_xy") or pc.get("pixel_size_xy_um")
        cur_z  = pc.get("pixel_size_z")  or pc.get("pixel_size_z_um")

        mismatches = []
        if imp_xy is not None and cur_xy is not None:
            if abs(float(imp_xy) - float(cur_xy)) / max(float(imp_xy), 1e-9) > 0.10:
                mismatches.append(f"XY: imported {imp_xy} µm/px vs current {cur_xy} µm/px")
        if imp_z is not None and cur_z is not None:
            if abs(float(imp_z) - float(cur_z)) / max(float(imp_z), 1e-9) > 0.10:
                mismatches.append(f"Z: imported {imp_z} µm/px vs current {cur_z} µm/px")

        if not mismatches:
            return True

        msg = (
            "Pixel sizes differ between the imported and current experiment:\n\n"
            + "\n".join(f"  • {m}" for m in mismatches)
            + "\n\nFeatures are computed in pixel space, so the same sigma value covers "
            "a different physical scale. This may reduce classifier performance.\n\n"
            "Proceed anyway?"
        )
        reply = QMessageBox.warning(
            self, "Pixel Size Mismatch", msg,
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        return reply == QMessageBox.Yes

    def _show_import_validation_dialog(self, imported_meta, data_by_celltype, cell_type_mapping):
        """Show a summary of imported pixel counts before committing, with a per-row checkbox
        to opt out of importing specific cell types.

        Returns the (possibly filtered) cell_type_mapping dict to proceed, with unchecked
        rows set to None, or None if the user cancelled.
        """
        td_info = imported_meta.get("training_data", {})

        dlg = QDialog(self)
        dlg.setWindowTitle("Imported Training Data — Summary")
        dlg.setMinimumWidth(480)
        lay = QVBoxLayout(dlg)

        img_meta = imported_meta.get("image_metadata", {})
        header_lines = []
        xy = img_meta.get("pixel_size_xy_um")
        z  = img_meta.get("pixel_size_z_um")
        if xy:
            header_lines.append(f"Pixel size: xy={xy} µm, z={z} µm")
        ch_names = img_meta.get("channel_names", [])
        if ch_names:
            header_lines.append(f"Channels: {', '.join(ch_names)}")
        if header_lines:
            lbl = QLabel("\n".join(header_lines))
            lbl.setStyleSheet("color:#555; font-size:10px;")
            lay.addWidget(lbl)

        active_rows = [
            (imp_ct, local_ct)
            for imp_ct, local_ct in cell_type_mapping.items()
            if local_ct is not None and imp_ct in data_by_celltype
        ]

        table = QTableWidget(len(active_rows), 5)
        table.setHorizontalHeaderLabels(["Import", "Imported type", "→ Local type", "Total px", "Fg / Bg"])
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        row_checkboxes = {}
        for row_idx, (imp_ct, local_ct) in enumerate(active_rows):
            counts = td_info.get(imp_ct, {})
            n = counts.get("n_pixels", len(data_by_celltype[imp_ct][1]))
            pos = counts.get("n_positive", int(np.sum(data_by_celltype[imp_ct][1] == 2)))
            neg = counts.get("n_negative", int(np.sum(data_by_celltype[imp_ct][1] == 1)))

            checkbox = QCheckBox()
            checkbox.setChecked(True)
            cb_container = QWidget()
            cb_layout = QHBoxLayout(cb_container)
            cb_layout.addWidget(checkbox)
            cb_layout.setAlignment(Qt.AlignCenter)
            cb_layout.setContentsMargins(0, 0, 0, 0)
            table.setCellWidget(row_idx, 0, cb_container)
            row_checkboxes[imp_ct] = checkbox

            table.setItem(row_idx, 1, QTableWidgetItem(str(imp_ct)))
            table.setItem(row_idx, 2, QTableWidgetItem(str(local_ct)))
            table.setItem(row_idx, 3, QTableWidgetItem(str(n)))
            table.setItem(row_idx, 4, QTableWidgetItem(f"{pos} / {neg}"))
        lay.addWidget(table)

        note = QLabel(
            "Uncheck a row to exclude that cell type from the import.\n"
            "Checked pixels will be prepended to your new labels when you click Train.\n"
            "The combined dataset is saved so this experiment can itself be imported later."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color: #555; font-style: italic; font-size: 10px; padding-top: 4px;")
        lay.addWidget(note)

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.button(QDialogButtonBox.Ok).setText("Proceed with Import")
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        lay.addWidget(btns)

        if dlg.exec_() != QDialog.Accepted:
            return None

        filtered_mapping = dict(cell_type_mapping)
        for imp_ct, checkbox in row_checkboxes.items():
            if not checkbox.isChecked():
                filtered_mapping[imp_ct] = None
        return filtered_mapping

    def _update_import_panel(self):
        """Refresh the import info panel from the training widget's active import."""
        if self._training_widget is None:
            return
        summary = self._training_widget.get_import_summary()
        if not summary:
            self.import_info_group.setVisible(False)
            return

        meta_path = self._training_widget._imported_metadata_path
        short_path = str(meta_path)
        if len(short_path) > 60:
            short_path = "…" + short_path[-57:]
        self.import_source_label.setText(f"Source: {short_path}")
        self.import_source_label.setToolTip(str(meta_path))
        self.import_details_label.setText(summary)
        self.import_info_group.setVisible(True)

    def _on_import_applied(self):
        """Slot: update panel after import is applied."""
        self._update_import_panel()

    def _on_import_cleared(self):
        """Slot: hide panel after import is cleared."""
        self.import_info_group.setVisible(False)

    def _on_clear_import_clicked(self):
        """Clear imported training data and remove it from YAML."""
        if self._training_widget is not None:
            self._training_widget.clear_import()
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc.pop("apoc_imported_training_path", None)
        self._flush_params_to_yaml()

    def _restore_import_on_session_load(self):
        """Re-apply a persisted import if the YAML holds an import path."""
        from behav3d.preprocessing.segmentation.apoc_train import (
            _load_training_metadata,
            _load_training_data_for_celltype,
            _load_training_bundle,
        )
        if self._training_widget is None:
            return
        pc = self.metadata_loader.behav3d_parameters.get("pixel_classifier", {}) or {}
        saved_path = pc.get("apoc_imported_training_path", "")
        if not saved_path or not Path(saved_path).exists():
            return

        try:
            if saved_path.endswith(".zip"):
                meta, data_by_celltype = _load_training_bundle(saved_path)
                if meta is None or not data_by_celltype:
                    return
                imported_cts = list(meta.get("cell_types", []))
                if meta.get("has_death"):
                    imported_cts.append("dead")
            else:
                meta = _load_training_metadata(saved_path)
                if meta is None:
                    return
                source_td = Path(saved_path).parent
                imported_cts = list(meta.get("cell_types", []))
                if meta.get("has_death"):
                    imported_cts.append("dead")
                data_by_celltype = {}
                for ct in imported_cts:
                    X, y = _load_training_data_for_celltype(source_td, ct)
                    if X is not None:
                        data_by_celltype[ct] = (X, y)
            if not data_by_celltype:
                return

            local_cts = list(self.all_cell_types)
            if self._training_widget.has_death:
                local_cts.append("dead")
            cell_type_mapping = {ct: ct for ct in imported_cts if ct in local_cts}

            # Run pixel size check silently on restore (just warn, no cancel)
            self._check_pixel_size_compatibility(meta)

            self._training_widget.apply_import(
                meta, data_by_celltype, cell_type_mapping, source_path=saved_path
            )
            self._update_import_panel()
            self.log(f"↩ Restored import from {saved_path}")
        except Exception as exc:
            self.log(f"⚠️ Could not restore import from {saved_path}: {exc}")

    def _on_training_finished_update_export_btn(self, *_):
        """Enable the export button once training data has been saved."""
        self._refresh_export_btn_state()

    def _refresh_export_btn_state(self):
        """Enable the export button if the training bundle exists for this experiment."""
        try:
            pixel_class_outdir = (
                Path(self.metadata_loader.output_dir) / "images" / "PixelClassification"
            )
            bundle = pixel_class_outdir / "PixelClassifier_TrainingData.zip"
            self.btn_export_bundle.setEnabled(
                bool(self.metadata_loader.output_dir) and bundle.is_file()
            )
        except Exception:
            pass

    # ── Queue param snapshot ────────────────────────────────────
    def _current_global_strategy(self):
        """Return the current global strategy from the training widget's combo."""
        if (
            self._training_widget is not None
            and self._training_widget.strategy_combo is not None
        ):
            return self._training_widget.strategy_combo.currentText()
        return self.all_strategies()[0]

    def validate_for_queue(self):
        """Return ``(ok, error_msg)`` for the queue-add path.

        Mirrors :meth:`CellposeWidget.validate_for_queue`: refuses to enqueue
        a step that would immediately fail at run-time because metadata is
        missing, no output directory is set, or no trained APOC classifier
        ``.cl`` files exist yet under ``output_dir/images/PixelClassification``.
        In *Advanced* mode every expected cell type must have its classifier;
        otherwise at least one is enough since a single classifier still runs.
        """
        md = self.metadata_loader.metadata
        if md is None:
            return False, "No metadata loaded. Please load metadata first (Data Preparation tab)."

        out_dir = self.metadata_loader.output_dir
        if not out_dir:
            return False, "No output directory established. Please set it in the Data Preparation tab."

        pixel_class_outdir = Path(out_dir) / "images" / "PixelClassification"
        expected_cts = self._apoc_expected_cell_types(md)
        if not expected_cts:
            return False, "No cell types detected from metadata. APOC needs at least one cell type to segment."

        existing = [
            ct for ct in expected_cts
            if self._apoc_classifier_path(pixel_class_outdir, ct).exists()
        ]
        missing = [ct for ct in expected_cts if ct not in existing]

        if not existing:
            return False, (
                "No trained APOC classifiers found in\n"
                f"{pixel_class_outdir}\n\n"
                "Please train the classifier(s) in the legend tab before queueing."
            )

        # Advanced strategy runs a separate classifier per cell type.
        from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget as _ATW
        strategy = self._current_global_strategy()
        if strategy == _ATW.ADVANCED_STRATEGY and missing:
            return False, (
                "Advanced strategy requires a trained classifier for every cell type.\n\n"
                f"Missing classifier(s) for: {', '.join(missing)}\n\n"
                "Please train the missing classifier(s) or switch to a non-Advanced strategy."
            )
        return True, ""

    def get_queue_params(self) -> dict:
        """Snapshot current widget state for use by the processing queue."""
        strategy = self._current_global_strategy()
        all_strats = self.all_strategies()
        params = {
            "strategy_index": all_strats.index(strategy) if strategy in all_strats else 0,
            "strategy_name": strategy,
            "gpu_device_name": self._selected_gpu_device_name(),
            "force_cpu": self.btn_force_cpu.isChecked() if hasattr(self, 'btn_force_cpu') else False,
            "overwrite": self.check_overwrite.isChecked(),
            "workers": self.spin_workers.value(),
            "process_all": self.check_process_all.isChecked(),
            "t_start": self.spin_t_start.value(),
            "t_end": self.spin_t_end.value(),
        }
        # Collect per-cell-type params from training widget tabs
        if self._training_widget is not None:
            per_ct_params = {}
            per_ct_strategies = {}
            for ct, tab in self._training_widget.tabs.items():
                cfg = tab.get_config()
                per_ct_params[ct] = cfg
                per_tab_combo = getattr(tab, "_per_tab_strategy_combo", None)
                if per_tab_combo is not None:
                    per_ct_strategies[ct] = per_tab_combo.currentText()
            params["per_ct_params"] = per_ct_params
            from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget as _ATW
            if strategy == _ATW.ADVANCED_STRATEGY:
                params["per_ct_strategies"] = per_ct_strategies
        return params

    # ── Parameters Saving ───────────────────────────────────────
    def _save_apoc_params_to_yaml(self, updated_apoc_params=None, sync_source_ct=None, skip_sync=False):
        out_dir = self.metadata_loader.output_dir
        if not out_dir:
            return
        
        # Update apoc dict
        apoc = self.metadata_loader.behav3d_parameters.setdefault("apoc", {})
        if updated_apoc_params is not None:
            # Signals like stateChanged / valueChanged may pass an int or str;
            # only merge if it is actually a dict.
            if not isinstance(updated_apoc_params, dict):
                pass  # ignore non-dict args (int from checkbox, str from combo, etc.)
            else:
                apoc.update(updated_apoc_params)
            
        # Collect top-level parameters from the GUI
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        if self._training_widget is not None and self._training_widget.strategy_combo is not None:
            pc["apoc_strategy"] = self._training_widget.strategy_combo.currentText()
        if hasattr(self, 'spin_examples') and self.spin_examples is not None:
            pc["examples_per_sample"] = self.spin_examples.value()
        if hasattr(self, 'combo_gpu_device') and self.combo_gpu_device is not None:
            pc["gpu_device_name"] = self._selected_gpu_device_name()
        if hasattr(self, 'btn_force_cpu') and self.btn_force_cpu is not None:
            pc["force_cpu"] = self.btn_force_cpu.isChecked()
            
        if hasattr(self, 'spin_workers') and hasattr(self, 'check_process_all'):
            pc["workers"] = self.spin_workers.value()
            pc["process_all"] = self.check_process_all.isChecked()
            pc["t_start"] = self.spin_t_start.value()
            pc["t_end"] = self.spin_t_end.value()

        params_path = Path(out_dir) / "behav3d_parameters.yml"
        import yaml
        try:
            with open(params_path, "w") as f:
                yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
        except Exception as e:
            self.log(f"Warning: Could not save parameters: {e}")

        # If unified mode is active, mirror organoid settings across tabs.
        if not skip_sync:
            if sync_source_ct is not None:
                self._maybe_sync_organoid_tabs(source_ct=sync_source_ct)
            else:
                self._maybe_sync_organoid_tabs()

    def _flush_params_to_yaml(self):
        """Persist current behav3d_parameters to the YAML file on disk."""
        self._save_apoc_params_to_yaml(skip_sync=True)

    # ── Export training bundle ──────────────────────────────────

    def _on_export_training_bundle(self):
        """Copy PixelClassifier_TrainingData.zip + .cl classifiers into a shareable archive."""
        import zipfile

        if not self.metadata_loader.output_dir:
            QMessageBox.warning(self, "No Output Directory", "No output directory is set.")
            return

        pixel_class_outdir = (
            Path(self.metadata_loader.output_dir) / "images" / "PixelClassification"
        )
        bundle = pixel_class_outdir / "PixelClassifier_TrainingData.zip"
        if not bundle.is_file():
            QMessageBox.warning(
                self, "No Training Data",
                "PixelClassifier_TrainingData.zip not found. Train at least one classifier first."
            )
            return

        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Training Bundle As",
            str(Path(self.metadata_loader.output_dir) / "behav3d_training_bundle.zip"),
            "ZIP archives (*.zip)",
        )
        if not save_path:
            return

        try:
            with zipfile.ZipFile(save_path, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(bundle, arcname=bundle.name)
                for cl_file in sorted(pixel_class_outdir.glob("*.cl")):
                    zf.write(cl_file, arcname=cl_file.name)

            QMessageBox.information(
                self, "Export Complete",
                f"Training bundle saved to:\n{save_path}"
            )
            self.log(f"📦 Training bundle exported to {save_path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export Error", f"Could not create bundle:\n{exc}")
            self.log(f"❌ Training bundle export failed: {exc}")

    # ── APOC GUI compatibility helpers (frontend-only) ──────────
    def _apoc_normalize_feature_spec(self, feature_spec):
        spec = str(feature_spec or "").replace(",", " ").replace("\t", " ").strip()
        while "  " in spec:
            spec = spec.replace("  ", " ")
        return spec

    def _apoc_parse_feature_tokens(self, feature_spec):
        normalized = self._apoc_normalize_feature_spec(feature_spec)
        return [tok for tok in normalized.split(" ") if tok]

    def _apoc_validate_feature_spec(self, feature_spec):
        tokens = self._apoc_parse_feature_tokens(feature_spec)
        if not tokens:
            return False, "empty feature_specification"

        for tok in tokens:
            lower_tok = tok.lower()
            if "_channel" in lower_tok:
                return False, "legacy expanded token format with _channel suffix"
            if lower_tok == "original":
                continue
            if "=" not in tok:
                return False, f"invalid token without '=': {tok}"
            name, value = tok.split("=", 1)
            if not name.strip() or not value.strip():
                return False, f"invalid token format: {tok}"
            try:
                float(value)
            except ValueError:
                return False, f"non-numeric feature parameter: {tok}"
        return True, "ok"

    def _apoc_read_classifier_header_value(self, clf_path, key, default=None):
        path = Path(clf_path)
        if not path.exists():
            return default

        prefix = f"{key} = "
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as f:
                line_count = 0
                for line in f:
                    line_count += 1
                    if line_count > 120:
                        break
                    stripped = line.strip()
                    if stripped == "*/":
                        break
                    if stripped.startswith(prefix):
                        return stripped[len(prefix):].strip()
        except Exception:
            return default
        return default

    def _apoc_rewrite_classifier_feature_spec(self, clf_path, new_feature_spec):
        path = Path(clf_path)
        if not path.exists():
            return False, f"classifier file not found: {path.name}"

        normalized = self._apoc_normalize_feature_spec(new_feature_spec)
        is_valid, reason = self._apoc_validate_feature_spec(normalized)
        if not is_valid:
            return False, f"cannot rewrite invalid feature specification ({reason})"

        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
            lines = text.splitlines(keepends=True)
            in_header = False
            replaced = False

            for i, line in enumerate(lines):
                stripped = line.strip()
                if stripped.startswith("/*"):
                    in_header = True
                if in_header and stripped.startswith("feature_specification = "):
                    ending = "\r\n" if line.endswith("\r\n") else "\n"
                    lines[i] = f"feature_specification = {normalized}{ending}"
                    replaced = True
                    break
                if in_header and stripped == "*/":
                    break

            if not replaced:
                return False, "feature_specification header line not found"

            path.write_text("".join(lines), encoding="utf-8")
            return True, "ok"
        except Exception as e:
            return False, str(e)

    def _apoc_classifier_path(self, pixel_class_outdir, ct):
        if ct == "dead":
            return Path(pixel_class_outdir) / "PixelClassifier_Death.cl"
        return Path(pixel_class_outdir) / f"PixelClassifier_{ct.capitalize()}.cl"

    def _apoc_expected_cell_types(self, metadata):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
        )
        organoid_types = detect_organoid_types_from_metadata(metadata)
        immune_types = detect_immune_cell_types_from_metadata(metadata)
        other_types = detect_other_cell_types_from_metadata(metadata)
        cell_types = organoid_types + immune_types + other_types
        if has_dead_channel(metadata):
            cell_types.append("dead")
        return cell_types

    def _apoc_try_migrate_from_tab_config(self, pixel_class_outdir, ct):
        if self._training_widget is None:
            return False, "training widget unavailable"
        if ct not in self._training_widget.tabs:
            return False, f"no GUI tab config for '{ct}'"

        cfg = self._training_widget.tabs[ct].get_config()
        canonical_spec = self._apoc_normalize_feature_spec(cfg.get("feature_string", ""))
        ok, msg = self._apoc_rewrite_classifier_feature_spec(
            self._apoc_classifier_path(pixel_class_outdir, ct),
            canonical_spec,
        )
        if ok:
            self.log(f"  🔧 {ct}: updated classifier header feature_specification from GUI tab config")
        return ok, msg

    def _apoc_preflight_classifier_headers(self, output_dir, metadata):
        pixel_class_outdir = Path(output_dir) / "images" / "PixelClassification"
        incompatible = []

        for ct in self._apoc_expected_cell_types(metadata):
            clf_path = self._apoc_classifier_path(pixel_class_outdir, ct)
            if not clf_path.exists():
                continue

            header_spec = self._apoc_read_classifier_header_value(
                clf_path,
                "feature_specification",
                "",
            )
            is_valid, reason = self._apoc_validate_feature_spec(header_spec)

            if not is_valid:
                mig_ok, mig_msg = self._apoc_try_migrate_from_tab_config(pixel_class_outdir, ct)
                if mig_ok:
                    header_spec = self._apoc_read_classifier_header_value(clf_path, "feature_specification", "")
                    is_valid, reason = self._apoc_validate_feature_spec(header_spec)
                else:
                    reason = f"{reason}; migration failed: {mig_msg}"

            if not is_valid:
                incompatible.append((ct, clf_path.name, reason, str(header_spec or "")))
            else:
                self.log(f"  ✅ {ct}: feature_specification OK -> {header_spec}")

        return incompatible

    # ── Global Run [CellType] button helpers ────────────────────
    def _on_global_run_instance(self):
        """Run batch segmentation for the currently active cell type with overwrite prompt."""
        tw = self._training_widget
        if tw is None:
            return
        idx = tw.tab_widget.currentIndex()
        ct = tw._tab_cell_types[idx] if idx < len(tw._tab_cell_types) else None
        if ct is None:
            return

        md = self.metadata_loader.metadata
        if md is None:
            self.log("⚠️ No metadata loaded.")
            return

        output_dir = Path(self.metadata_loader.output_dir)

        # Check for pre-existing segments and prompt
        sample_names = md['sample_name'].unique()
        existing = [
            f"{ct} segments for {sn}"
            for sn in sample_names
            if (output_dir / "images" / sn / f"{sn}_{ct}_segments.zarr").exists()
        ]

        overwrite = False
        if existing:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_single
            choice = prompt_overwrite_single(
                self,
                f"Run {ct.capitalize()} Segmentation",
                existing,
            )
            if choice == "cancel":
                return
            overwrite = (choice == "overwrite")

        # Run segmentation for this cell type only.  Use the background
        # path (``block=False``) so napari remains interactive.
        orig_overwrite = self.check_overwrite.isChecked()
        self.check_overwrite.setChecked(overwrite)
        try:
            self._on_run_segmentation(
                interactive=True, only_cell_types=[ct], block=False,
            )
        finally:
            self.check_overwrite.setChecked(orig_overwrite)

    def _update_global_run_btn_label(self, tab_index):
        """Update the global Run button label to reflect the active cell-type tab."""
        tw = self._training_widget
        btn = getattr(self, "_global_run_instance_btn", None)
        if tw is None or btn is None:
            return
        ct = tw._tab_cell_types[tab_index] if tab_index < len(tw._tab_cell_types) else "?"
        btn.setText(f"▶ Run {ct.capitalize()} Segmentation")

    # ── Helper: unify organoid settings ────────────────────────
    def _on_unify_organoids_changed(self, state):
        """When toggled ON, do an initial full sync from the first organoid tab."""
        from qtpy.QtCore import Qt
        self._unify_organoids = (state == Qt.Checked)
        if not self._unify_organoids or self._training_widget is None:
            return
        # Perform initial full sync from the first organoid tab
        if self._organoid_types:
            self._maybe_sync_organoid_tabs(source_ct=self._organoid_types[0])

    def _maybe_sync_organoid_tabs(self, source_ct=None):
        """If unified mode is on, mirror the source organoid tab to all other organoid tabs.

        If *source_ct* is not given, the currently visible tab is used.
        Does nothing if it's not an organoid tab or the guard flag is set.
        """
        if not self._unify_organoids or self._syncing_organoids:
            return
        tw = self._training_widget
        if tw is None:
            return

        # Determine source cell type
        if source_ct is None:
            idx = tw.tab_widget.currentIndex()
            source_ct = tw._tab_cell_types[idx] if idx < len(tw._tab_cell_types) else None
        if source_ct is None or source_ct not in self._organoid_types:
            return  # only sync when an organoid tab is the source

        source_tab = tw.tabs.get(source_ct)
        if source_tab is None:
            return

        other_org_tabs = [
            tw.tabs.get(ct) for ct in self._organoid_types
            if ct != source_ct and ct in tw.tabs
        ]
        if not other_org_tabs:
            return

        self._syncing_organoids = True
        try:
            src_cfg = source_tab.get_config()
            src_strat = getattr(
                getattr(source_tab, "_per_tab_strategy_combo", None),
                "currentText", lambda: None
            )()
            per_tab_strategies = self.per_tab_strategies()
            for tab in other_org_tabs:
                if tab is None:
                    continue
                tab.apply_config(src_cfg)
                per_combo = getattr(tab, "_per_tab_strategy_combo", None)
                if per_combo is not None and src_strat and src_strat in per_tab_strategies:
                    per_combo.blockSignals(True)
                    per_combo.setCurrentText(src_strat)
                    per_combo.blockSignals(False)
                    params_snap = self.metadata_loader.behav3d_parameters.get("apoc", {})
                    tab.rebuild_instance_controls(strategy=src_strat, initial_params=params_snap)
                # Mirror instance segmentation parameter values from source
                for attr in (
                    "prob_mask_threshold_spin", "prob_seed_threshold_spin",
                    "edt_threshold_spin", "segment_size_min_spin", "opening_nr_pixels_spin",
                ):
                    src_w = getattr(source_tab, attr, None)
                    dst_w = getattr(tab, attr, None)
                    if src_w is not None and dst_w is not None:
                        dst_w.blockSignals(True)
                        dst_w.setValue(src_w.value())
                        dst_w.blockSignals(False)
                src_fh = getattr(source_tab, "fill_holes_cb", None)
                dst_fh = getattr(tab, "fill_holes_cb", None)
                if src_fh is not None and dst_fh is not None:
                    dst_fh.blockSignals(True)
                    dst_fh.setChecked(src_fh.isChecked())
                    dst_fh.blockSignals(False)
                # Refresh preview after all updates to ensure it reflects all changes
                # (apply_config updates all config but then spinners are changed after, so we need to update preview again)
                tab._update_preview()
            # Persist mirrored destination tab values without triggering another sync cycle.
            self._save_apoc_params_to_yaml(
                updated_apoc_params=self._collect_apoc_tab_config(),
                skip_sync=True,
            )
        finally:
            self._syncing_organoids = False

    # ── Run batch segmentation ──────────────────────────────────
    def _on_run_segmentation(self, interactive=True, only_cell_types=None,
                             block=True, extra_callbacks=None):
        """Run APOC batch segmentation.

        ``block=True`` (default) keeps the call synchronous so the
        processing queue and other internal callers see the original
        behaviour.  ``block=False`` (GUI button) executes the work in a
        background worker with an indeterminate busy spinner (the APOC
        backend does not currently expose a per-sample progress hook).
        ``extra_callbacks`` is the queue's chaining hook.
        """
        try:
            md = self.metadata_loader.metadata
            if md is None:
                self.log("\u26a0\ufe0f No metadata loaded.")
                fire_extra_callback(extra_callbacks, "on_failed", "no metadata loaded")
                return

            if not block and self._bg.is_running():
                self.log("\u26a0\ufe0f An APOC run is already in progress.")
                fire_extra_callback(extra_callbacks, "on_failed", "already running")
                return

            from behav3d.preprocessing.segmentation.apoc_segment import run_apoc_segmentation

            output_dir = Path(self.metadata_loader.output_dir)
            strategy = self._current_global_strategy()

            # ---- Overwrite check (batch only — per-cell-type uses its own
            # prompt in ``_on_global_run_instance``) ----------------------
            if interactive and only_cell_types is None:
                check_cts = list(self.all_cell_types or [])
                existing = []
                for sn in md["sample_name"].unique():
                    for ct in check_cts:
                        seg_path = output_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                        if seg_path.exists():
                            existing.append(f"{ct} segments for {sn}")
                if existing:
                    from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
                    choice = prompt_overwrite_batch(
                        self,
                        "Overwrite Existing APOC Segmentations?",
                        existing,
                    )
                    if choice == "cancel":
                        self.log("APOC segmentation cancelled.")
                        fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                        return
                    self.check_overwrite.setChecked(choice == "overwrite")

            # Timepoint range
            if only_cell_types is not None and hasattr(self, 'check_limit_timerange_ct'):
                use_all = self.check_limit_timerange_ct.isChecked()
                t_start = self.spin_t_start_ct.value()
                t_end = self.spin_t_end_ct.value()
            else:
                use_all = self.check_process_all.isChecked()
                t_start = self.spin_t_start.value()
                t_end = self.spin_t_end.value()

            if use_all:
                timepoint_range = None
            else:
                if t_start > t_end:
                    self.log("Error: Start timepoint must be <= End timepoint.")
                    fire_extra_callback(extra_callbacks, "on_failed", "bad timepoint range")
                    return
                timepoint_range = (t_start, t_end)

            # Collect APOC config from training widget tabs (Qt thread).
            apoc_config = {}
            if self._training_widget is not None:
                for ct, tab in self._training_widget.tabs.items():
                    cfg = tab.get_config()
                    for k, v in cfg.items():
                        apoc_config[f"apoc_{ct}_{k}"] = v

            # Save all parameters safely
            self._save_apoc_params_to_yaml(updated_apoc_params=apoc_config)

            # Fix channel names for apoc_segment.py parser (backend expects 'Ch 0' format)
            import re
            for key in list(apoc_config.keys()):
                if key.endswith("_channels") and isinstance(apoc_config[key], list):
                    fixed_chans = []
                    for ch_name in apoc_config[key]:
                        match = re.search(r'\((\d+)\)', str(ch_name))
                        if match:
                            fixed_chans.append(f"Ch {match.group(1)}")
                        else:
                            fixed_chans.append(ch_name)
                    apoc_config[key] = fixed_chans

            # GUI-side compatibility preflight: ensure classifier headers can drive predict(image=...).
            incompatible = self._apoc_preflight_classifier_headers(output_dir, md)
            if incompatible:
                self.log("❌ APOC classifier preflight failed. Incompatible classifier header(s):")
                for ct, fname, reason, header_spec in incompatible:
                    shown = header_spec if header_spec else "<empty>"
                    self.log(f"   - {ct} ({fname}): {reason}. header feature_specification='{shown}'")
                self.log("Please retrain the listed classifier(s) from the GUI and run batch segmentation again.")
                return

            self.log(f"Starting APOC batch segmentation ({strategy})...")

            self._apply_apoc_gpu_selection(log_message=False)

            metadata_csv = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv", "")

            per_ct_strategies = None
            from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget as _ATW
            if strategy == _ATW.ADVANCED_STRATEGY and self._training_widget is not None:
                per_ct_strategies = {}
                for ct, tab in self._training_widget.tabs.items():
                    per_tab_combo = getattr(tab, "_per_tab_strategy_combo", None)
                    per_ct_strategies[ct] = (
                        per_tab_combo.currentText() if per_tab_combo is not None
                        else self.per_tab_strategies()[0]
                    )
                apoc_dict = self.metadata_loader.behav3d_parameters.setdefault("apoc", {})
                apoc_dict["per_ct_strategies"] = per_ct_strategies
                self._save_apoc_params_to_yaml()
                strat_summary = ", ".join(f"{k}:{v[:4]}" for k, v in per_ct_strategies.items())
                self.log(f"  Advanced strategies: {strat_summary}")

            overwrite_val = self.check_overwrite.isChecked()
            n_workers_val = self.spin_workers.value()
            gpu_device_val = self._selected_gpu_device_name()

            def _do_apoc(progress_cb=None):
                return run_apoc_segmentation(
                    output_dir=str(output_dir),
                    metadata=md,
                    metadata_csv_path=metadata_csv,
                    timepoint_range=timepoint_range,
                    apoc_config=apoc_config,
                    overwrite_existing=overwrite_val,
                    n_workers=n_workers_val,
                    gpu_device=gpu_device_val,
                    apoc_strategy=strategy,
                    per_ct_strategies=per_ct_strategies,
                    only_cell_types=only_cell_types,
                )

            def _apply_apoc_result(updated_metadata):
                if updated_metadata is not None:
                    self.metadata_loader.metadata = updated_metadata
                    updated_metadata = _run_multicolor_cleanup_if_needed(self.metadata_loader, self.log, skip_unified=False)
                    self.metadata_loader.metadata = updated_metadata
                    if metadata_csv:
                        try:
                            updated_metadata.to_csv(metadata_csv, index=False)
                        except Exception as e:
                            self.log(f"Warning: Could not save metadata: {e}")
                self.log("✅ APOC batch segmentation finished!")
                if interactive:
                    self._prompt_visualize_after_apoc_segmentation()

            if block:
                try:
                    result = _do_apoc(progress_cb=None)
                    _apply_apoc_result(result)
                    fire_extra_callback(extra_callbacks, "on_done", result)
                except Exception as e:
                    traceback.print_exc()
                    self.log(f"❌ Error during APOC segmentation: {e}")
                    fire_extra_callback(extra_callbacks, "on_failed", str(e))
                return

            def _on_done_async(result):
                _apply_apoc_result(result)
                fire_extra_callback(extra_callbacks, "on_done", result)

            def _on_failed(err: str):
                self.log(f"❌ Error during APOC segmentation: {err}")
                fire_extra_callback(extra_callbacks, "on_failed", err)

            self._bg.run(
                fn=_do_apoc,
                desc=f"APOC segmentation ({strategy})\u2026",
                progress_row=self.tab_progress_row,
                buttons=[],
                viewer=self.viewer,
                on_done=_on_done_async,
                on_failed=_on_failed,
                inject_progress=False,
                indeterminate=True,
            )

        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Error during APOC segmentation: {e}")
            fire_extra_callback(extra_callbacks, "on_failed", str(e))


class ConvPaintWidget(QWidget):
    """ConvPaint (DL pixel classifier) segmentation page.

    Mirrors :class:`APOCWidget` for the unified ConvPaint backbone:
    a single multi-class model annotated through the legend tab inside the
    embedded :class:`ConvPaintTrainingWidget` (which is built lazily after
    the user generates training data because it needs ``all_images``).

    The plugin wraps device selection, batch run options, per-cell-type
    strategy snapshots and queueing on top of the training widget — the
    training widget itself handles all label/strategy UI through its
    ``per_cell_type_strategy=True`` constructor flag.
    """

    @classmethod
    def all_strategies(cls):
        from behav3d.preprocessing.segmentation.convpaint_train import (
            ConvPaintTrainingWidget as _CTW,
        )
        return list(_CTW.STRATEGIES) + [_CTW.ADVANCED_STRATEGY]

    @classmethod
    def per_tab_strategies(cls):
        from behav3d.preprocessing.segmentation.convpaint_train import (
            ConvPaintTrainingWidget as _CTW,
        )
        return list(_CTW.STRATEGIES)

    def __init__(
        self,
        viewer,
        metadata_loader,
        log_callback=None,
        parent=None,
        switch_to_visualization_callback=None,
        tab_progress_row=None,
    ):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self._switch_to_visualization = switch_to_visualization_callback
        self._training_widget = None
        self._is_session_active = False
        self._all_images_cache = None
        self.all_cell_types = []           # set when metadata loaded
        # Background-execution infrastructure (shared progress row).
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        self._init_ui()

    # ── UI construction ──────────────────────────────────────────
    def _init_ui(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setMinimumWidth(0)

        content = QWidget()
        content.setMinimumWidth(0)
        layout = QVBoxLayout(content)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        info = QLabel(
            "<b>ConvPaint Pixel Classifier.</b> "
            "Uses <tt>napari-convpaint</tt> + a deep feature extractor "
            "(VGG/DINOv2/...) to train a single multi-class classifier "
            "for all cell types at once."
        )
        info.setWordWrap(True)
        info.setMinimumWidth(0)
        info.setStyleSheet("color: #888; font-size: 11px; padding: 4px 0 6px 0;")
        layout.addWidget(info)

        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {}) or {}

        # ── Device selection (PyTorch / CUDA) ─────────────────────
        from behav3d.preprocessing.segmentation.convpaint_train import (
            _detect_torch_devices, _default_device,
        )
        self._convpaint_devices = _detect_torch_devices()  # list of (label, value)
        gpu_row = QHBoxLayout()
        gpu_row.addWidget(QLabel("Device:"))
        self.combo_gpu_device = QComboBox()
        self.combo_gpu_device.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        for label, value in self._convpaint_devices:
            self.combo_gpu_device.addItem(label, value)

        saved_device = str(pc.get("convpaint_device", "")).strip() or _default_device()
        idx = self.combo_gpu_device.findData(saved_device)
        if idx >= 0:
            self.combo_gpu_device.setCurrentIndex(idx)
        gpu_row.addWidget(self.combo_gpu_device, stretch=1)
        gpu_row.addStretch()
        layout.addLayout(gpu_row)

        # ── Training section (embedded ConvPaintTrainingWidget) ──
        self.training_group = QGroupBox("🎯 ConvPaint Classifier Training")
        self.training_layout = QVBoxLayout(self.training_group)
        self.training_layout.setContentsMargins(4, 4, 4, 4)

        train_ctrl_lay = QHBoxLayout()
        train_ctrl_lay.addWidget(QLabel("Examples/sample:"))
        self.spin_examples = QSpinBox()
        self.spin_examples.setValue(int(pc.get("examples_per_sample", 3)))
        self.spin_examples.setRange(1, 10)
        self.spin_examples.setMaximumWidth(70)
        train_ctrl_lay.addWidget(self.spin_examples)

        self.btn_load_training = QPushButton("Generate Training Data")
        self.btn_load_training.setToolTip(
            "Clears viewer and loads selected timepoints for unified ConvPaint labeling"
        )
        self.btn_load_training.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        self.btn_load_training.clicked.connect(
            lambda _: self._on_load_training_clicked(interactive=True)
        )
        train_ctrl_lay.addWidget(self.btn_load_training)
        train_ctrl_lay.addStretch()
        self.training_layout.addLayout(train_ctrl_lay)

        # Legend container for ConvPaint
        self.legend_container = QWidget()
        self.legend_container.setLayout(QVBoxLayout())
        self.legend_container.layout().setContentsMargins(0, 0, 0, 0)
        self.legend_container.setVisible(False)
        self.training_layout.addWidget(self.legend_container)

        self.training_placeholder = QLabel(
            "Load metadata and click <b>Generate Training Data</b> to enable "
            "ConvPaint classifier training."
        )
        self.training_placeholder.setWordWrap(True)
        self.training_placeholder.setStyleSheet("color:#888; font-style:italic; padding:10px;")
        self.training_layout.addWidget(self.training_placeholder)

        layout.addWidget(self.training_group)

        # ── Batch controls ───────────────────────────────────────
        batch_group = QGroupBox("🚀 Batch Segmentation")
        batch_lay = QVBoxLayout(batch_group)
        batch_lay.setSpacing(4)

        self.check_process_all = QCheckBox("Process ALL timepoints")
        self.check_process_all.setChecked(True)
        batch_lay.addWidget(self.check_process_all)

        tp_row = QHBoxLayout()
        tp_row.addWidget(QLabel("From t:"))
        self.spin_t_start = QSpinBox()
        self.spin_t_start.setRange(0, 99999)
        self.spin_t_start.setValue(0)
        self.spin_t_start.setMaximumWidth(70)
        tp_row.addWidget(self.spin_t_start)
        tp_row.addWidget(QLabel("to t:"))
        self.spin_t_end = QSpinBox()
        self.spin_t_end.setRange(0, 99999)
        self.spin_t_end.setValue(100)
        self.spin_t_end.setMaximumWidth(70)
        tp_row.addWidget(self.spin_t_end)
        tp_row.addStretch()
        self.tp_widget = QWidget()
        self.tp_widget.setLayout(tp_row)
        batch_lay.addWidget(self.tp_widget)

        def _toggle_tp(_state):
            self.tp_widget.setVisible(not self.check_process_all.isChecked())
        self.check_process_all.stateChanged.connect(_toggle_tp)
        _toggle_tp(None)

        n_cores = os.cpu_count() or 4
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, max(1, n_cores - 1))
        self.spin_workers.setValue(1)
        self.spin_workers.setMaximumWidth(60)
        w_form = QFormLayout()
        w_form.setContentsMargins(0, 0, 0, 0)
        w_form.addRow("Workers:", self.spin_workers)
        batch_lay.addLayout(w_form)

        self.check_overwrite = QCheckBox("Overwrite existing results")
        self.check_overwrite.setChecked(False)
        batch_lay.addWidget(self.check_overwrite)

        self.btn_run_segmentation = QPushButton("▶ Run ConvPaint Batch Segmentation")
        self.btn_run_segmentation.setStyleSheet(
            "QPushButton { background: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px; } "
            "QPushButton:hover { background: #0069d9; }"
        )
        self.btn_run_segmentation.clicked.connect(lambda _: self._on_run_segmentation(interactive=True, block=False))

        self.btn_queue_segment = QPushButton("+🛒")
        self.btn_queue_segment.setFixedSize(36, 36)
        self.btn_queue_segment.setToolTip("Add ConvPaint Segmentation to Queue")
        self.btn_queue_segment.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_btn_row = QHBoxLayout()
        batch_btn_row.setSpacing(4)
        batch_btn_row.addWidget(self.btn_run_segmentation, stretch=1)
        batch_btn_row.addWidget(self.btn_queue_segment)
        batch_lay.addLayout(batch_btn_row)

        layout.addWidget(batch_group)
        layout.addStretch()
        scroll.setWidget(content)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)

        # ── Auto-save bindings ───────────────────────────────────
        self.spin_examples.valueChanged.connect(lambda _: self._save_params_to_yaml())
        self.spin_workers.valueChanged.connect(lambda _: self._save_params_to_yaml())
        self.combo_gpu_device.currentTextChanged.connect(lambda _: self._save_params_to_yaml())
        self.check_process_all.stateChanged.connect(lambda _: self._save_params_to_yaml())
        self.spin_t_start.valueChanged.connect(lambda _: self._save_params_to_yaml())
        self.spin_t_end.valueChanged.connect(lambda _: self._save_params_to_yaml())

    # ── Helpers ──────────────────────────────────────────────────
    def _selected_device(self):
        return str(self.combo_gpu_device.currentData() or "auto")

    def _prompt_visualize_after_segmentation(self):
        res = QMessageBox.question(
            self,
            "Segmentation Finished",
            "Batch segmentation finished successfully! \n\n"
            "Do you want to switch to the Visualization Tab and see the segments?",
            QMessageBox.Yes | QMessageBox.No,
        )
        if res != QMessageBox.Yes:
            return
        if callable(self._switch_to_visualization):
            self._switch_to_visualization()
        else:
            self.log("  Visualization tab callback unavailable.")

    # ── Training widget signal slots ─────────────────────────────
    def _on_training_started(self, cell_types_to_train):
        device = self._selected_device()
        self.log(f"🖥️ ConvPaint training device: {device}")
        self.log(f"▶ Starting ConvPaint training for: {cell_types_to_train}")

    def _on_training_finished(self, cell_types_to_train, _result):
        self._reorder_training_layers()
        self.log("✅ ConvPaint training process completed.")

    def _on_instance_preview_started(self, ct):
        tw = self._training_widget
        strategy = tw._resolve_strategy(ct) if tw is not None else ""
        self.log(f"🔍 Running instance segmentation preview for '{ct}' (strategy: {strategy})...")

    def _on_instance_preview_finished(self, ct):
        self._reorder_training_layers()

    def _on_training_widget_strategy_changed(self, strategy):
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["convpaint_strategy"] = str(strategy)
        self._save_params_to_yaml()

    def _reorder_training_layers(self):
        if self._training_widget is None:
            return
        try:
            from behav3d.preprocessing.segmentation.convpaint_train import (
                _reorder_convpaint_layers,
            )
            tw = self._training_widget
            _reorder_convpaint_layers(
                self.viewer,
                [ct for ct in tw._tab_cell_types if ct != "dead"],
                has_death=tw.has_death,
            )
        except Exception as e:
            self.log(f"⚠️ Could not reorder ConvPaint training layers: {e}")

    # ── Metadata update (rebuilds nothing until training data is loaded) ──
    def _on_metadata_updated(self):
        """Called when metadata is loaded/updated.

        Mirrors the APOC behaviour: build a parameter-only
        :class:`ConvPaintTrainingWidget` with ``all_images=[]`` immediately so
        the user can tune every segmentation parameter (strategy, EDT
        threshold, min segment size, …) even before training data is loaded.
        Training and preview buttons are disabled in this mode; they become
        fully active once the user clicks *Generate Training Data*.
        """
        # Tear down any existing widget (stale session or previous metadata).
        if self._training_widget is not None:
            try:
                self._training_widget.cleanup()
            except Exception:
                pass
            self.training_layout.removeWidget(self._training_widget)
            self._training_widget.hide()
            self._training_widget.deleteLater()
            self._training_widget = None

        self._is_session_active = False
        self._all_images_cache = None

        md = self.metadata_loader.metadata
        if md is None:
            self.all_cell_types = []
            if self.training_placeholder.parent() is None:
                self.training_layout.addWidget(self.training_placeholder)
            self.training_placeholder.show()
            return

        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
        )
        from behav3d.preprocessing.segmentation.convpaint_train import (
            _filter_merge_types,
            ConvPaintTrainingWidget,
            DEFAULT_STRATEGY,
            _normalize_strategy,
        )

        organoid_types = _filter_merge_types(detect_organoid_types_from_metadata(md))
        immune_types   = _filter_merge_types(detect_immune_cell_types_from_metadata(md))
        other_types    = _filter_merge_types(detect_other_cell_types_from_metadata(md))
        all_cell_types = organoid_types + immune_types + other_types
        has_death      = has_dead_channel(md)
        self.all_cell_types = list(all_cell_types)

        if not all_cell_types:
            if self.training_placeholder.parent() is None:
                self.training_layout.addWidget(self.training_placeholder)
            self.training_placeholder.show()
            return

        # Hide the placeholder now that we will embed the real widget.
        if self.training_placeholder.parent() is not None:
            self.training_layout.removeWidget(self.training_placeholder)
            self.training_placeholder.hide()

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        params      = self.metadata_loader.behav3d_parameters
        ip          = dict(params.get("pixel_classifier", {}) or {})
        ip["convpaint_device"] = self._selected_device()
        cp_strategy = _normalize_strategy(ip.get("convpaint_strategy", DEFAULT_STRATEGY))

        if "pixel_distance_xy" in md.columns and "distance_unit" in md.columns:
            _unit = str(md["distance_unit"].iloc[0])
            _xy_from_md = convert_distance(float(md["pixel_distance_xy"].iloc[0]), _unit)
            _z_from_md  = convert_distance(float(md["pixel_distance_z"].iloc[0]),  _unit)
        else:
            _xy_from_md = None
            _z_from_md  = None
        _pixel_sizes = {
            "xy_um": ip.get("pixel_size_xy") or ip.get("pixel_size_xy_um") or _xy_from_md,
            "z_um":  ip.get("pixel_size_z")  or ip.get("pixel_size_z_um")  or _z_from_md,
        }

        self._training_widget = ConvPaintTrainingWidget(
            viewer=self.viewer,
            all_images=[],          # no training data yet — parameter-only mode
            pixel_class_outdir=str(pixel_class_outdir),
            all_cell_types=all_cell_types,
            has_death=has_death,
            initial_params=ip,
            on_params_changed=self._save_params_to_yaml,
            convpaint_strategy=cp_strategy,
            unified_input_channels=[],
            death_input_channels=[],
            pixel_sizes=_pixel_sizes,
            per_cell_type_strategy=True,
            show_device=False,
            external_log=self.log,
        )
        if self._training_widget.strategy_combo is not None:
            self._training_widget.strategy_combo.setCurrentText(cp_strategy)

        # Disable training/preview controls — images are not loaded yet.
        _no_data_tip = (
            "Training data not loaded. Click 'Generate Training Data' above to enable."
        )
        self._training_widget.btn_train.setEnabled(False)
        self._training_widget.btn_train.setToolTip(_no_data_tip)
        self._training_widget.btn_rerun_preview.setEnabled(False)
        self._training_widget.btn_rerun_preview.setToolTip(_no_data_tip)
        for tab in self._training_widget.tabs.values():
            if hasattr(tab, "btn_preview") and tab.btn_preview is not None:
                tab.btn_preview.setEnabled(False)
                tab.btn_preview.setToolTip(_no_data_tip)

        # Disable training-specific parameter groups (Feature Extractor,
        # Classifier, tiling) — the user only needs the segmentation
        # post-processing parameters (EDT, probability, watershed) here.
        self._training_widget.disable_training_params()

        # Wire only the parameter-persistence signal (training signals are
        # connected later in _on_load_training_clicked when images are loaded).
        self._training_widget.strategy_changed.connect(
            self._on_training_widget_strategy_changed
        )

        self.training_layout.addWidget(self._training_widget)

    # ── Generate Training Data ───────────────────────────────────
    def _on_load_training_clicked(self, interactive=True):
        try:
            md = self.metadata_loader.metadata
            if md is None:
                self.log("⚠️ Cannot generate training data: No metadata loaded.")
                return

            self.log("Loading ConvPaint training data...")

            from behav3d.core.metadata import (
                detect_organoid_types_from_metadata,
                detect_immune_cell_types_from_metadata,
                detect_other_cell_types_from_metadata,
                has_dead_channel,
            )
            from behav3d.preprocessing.segmentation.convpaint_train import (
                _load_training_images,
                _filter_merge_types,
                _derive_convpaint_input_channels,
                _reorder_convpaint_layers,
                UNIFIED_LABELS_LAYER_NAME,
                DEAD_LABELS_LAYER_NAME,
                UNIFIED_PREDICTED_LAYER_NAME,
                DEAD_PREDICTED_LAYER_NAME,
                _predicted_labels_path,
                _probability_map_path,
                _segments_layer_name,
                _probability_layer_name,
                CHANNEL_COLORS,
                ConvPaintTrainingWidget,
                DEFAULT_STRATEGY,
                _normalize_strategy,
            )
            from behav3d.preprocessing.segmentation.convpaint_label_map import (
                build_label_map, label_map_matches, load_label_map,
                unified_label_map_path, unified_user_labels_path,
                unified_predicted_labels_path,
            )

            output_dir = Path(self.metadata_loader.output_dir)

            organoid_types = _filter_merge_types(
                detect_organoid_types_from_metadata(md)
            )
            immune_types = _filter_merge_types(
                detect_immune_cell_types_from_metadata(md)
            )
            other_types = _filter_merge_types(
                detect_other_cell_types_from_metadata(md)
            )

            # Persist examples/sample before loading.
            pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
            pc["examples_per_sample"] = int(self.spin_examples.value())
            self._save_params_to_yaml()

            # Decide load-existing vs regenerate via interactive prompt.
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            image_outpath = pixel_class_outdir / "PixelClassifier_Images.zarr"

            overwrite_images = False
            if image_outpath.exists():
                if interactive:
                    try:
                        import zarr as _zarr
                        _cached = _zarr.open(str(image_outpath), mode="r")
                        # ConvPaint saves as (C, T, Z, Y, X) → shape[1] = image count
                        saved_image_count = _cached.shape[1] if len(_cached.shape) > 1 else _cached.shape[0]
                        del _cached
                    except Exception:
                        saved_image_count = None

                    n_samples = len(md["sample_name"].unique())
                    examples_per_sample = self.spin_examples.value()
                    total_new_images = n_samples * examples_per_sample
                    if saved_image_count is not None:
                        msg = (
                            f"Training data already exists with {saved_image_count} images.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            "Overwrite with new training data, or load the existing data?"
                        )
                    else:
                        msg = (
                            "Training data already exists.\n\n"
                            f"Currently selected: {n_samples} sample(s) × {examples_per_sample} examples/sample "
                            f"= {total_new_images} images.\n\n"
                            "Overwrite with new training data, or load the existing data?"
                        )
                    box = QMessageBox(self)
                    box.setWindowTitle("Training Data Detected")
                    box.setText(msg)
                    btn_generate = box.addButton("Generate New", QMessageBox.AcceptRole)
                    btn_load = box.addButton("Load Existing", QMessageBox.YesRole)
                    btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                    box.exec_()
                    clicked = box.clickedButton()
                    if clicked is btn_cancel:
                        self.log("Action cancelled.")
                        return
                    overwrite_images = (clicked is btn_generate)
                else:
                    overwrite_images = False

            self.viewer.layers.clear()

            all_images, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
                metadata=md,
                output_dir=str(output_dir),
                examples_per_sample=int(self.spin_examples.value()),
                organoid_types=organoid_types,
                immune_types=immune_types,
                other_types=other_types,
                overwrite_images=overwrite_images,
            )
            all_cell_types = _filter_merge_types(all_cell_types)
            self.all_cell_types = list(all_cell_types)

            if not all_images:
                self.log("⚠️ No training images loaded.")
                return
            if not all_cell_types:
                self.log("⚠️ No cell types detected — nothing to train.")
                return

            n_channels = all_images[0].shape[0]
            unified_input_channels, death_input_channels = _derive_convpaint_input_channels(
                md, n_channels
            )
            self.log(
                f"Loaded {len(all_images)} training images with {n_channels} channels."
            )
            self.log(f"  Unified input channels: {unified_input_channels}")
            if has_death:
                self.log(f"  Death input channels: {death_input_channels}")

            stacked = np.stack(all_images, axis=0)  # (T_total, C, Z, Y, X)
            T_total = stacked.shape[0]
            label_shape = (T_total,) + stacked.shape[2:]

            for ch in range(n_channels):
                ch_data = stacked[:, ch, :, :, :]
                nonzero = ch_data[ch_data > 0]
                clim = (
                    (0, float(np.percentile(nonzero, 99.8)))
                    if nonzero.size > 0 else (0, 1e-3)
                )
                img_layer = self.viewer.add_image(
                    ch_data,
                    name=f"Channel {ch}",
                    contrast_limits=clim,
                    colormap=CHANNEL_COLORS[ch % len(CHANNEL_COLORS)],
                    blending="additive",
                    opacity=0.8,
                )
                img_layer.contrast_limits_range = (0, float(ch_data.max()))

            # Build (or rebuild) the active label map.
            label_map = build_label_map(all_cell_types)

            # Restore/create unified User Provided Labels.
            unified_path = unified_user_labels_path(pixel_class_outdir)
            saved_map_path = unified_label_map_path(pixel_class_outdir)
            unified_data = np.zeros(label_shape, dtype=np.int16)
            if unified_path.exists():
                try:
                    existing = np.asarray(load_zarr(unified_path))
                    if existing.shape == label_shape:
                        if saved_map_path.exists():
                            saved_map = load_label_map(saved_map_path)
                            if label_map_matches(saved_map, all_cell_types):
                                unified_data = existing
                                self.log(
                                    f"  ↩ Restored unified labels ({unified_path.name})"
                                )
                            else:
                                self.log(
                                    "  ⚠️ Cached unified labels were trained with a "
                                    "different cell-type set — starting fresh."
                                )
                        else:
                            unified_data = existing
                            self.log(
                                f"  ↩ Restored unified labels (no map sidecar)."
                            )
                    else:
                        self.log(
                            f"  ⚠️ Unified labels shape mismatch ({existing.shape} vs "
                            f"{label_shape}); starting fresh."
                        )
                except Exception as exc:
                    self.log(f"  ⚠️ Could not restore unified labels: {exc}")

            self.viewer.add_labels(
                unified_data, name=UNIFIED_LABELS_LAYER_NAME, opacity=0.5,
            )

            if has_death:
                dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
                dead_data = np.zeros(label_shape, dtype=np.int16)
                if dead_path.exists():
                    try:
                        existing = np.asarray(load_zarr(dead_path))
                        if existing.shape == label_shape:
                            dead_data = existing
                            self.log(f"  ↩ Restored death labels ({dead_path.name})")
                    except Exception as exc:
                        self.log(f"  ⚠️ Could not restore death labels: {exc}")
                self.viewer.add_labels(dead_data, name=DEAD_LABELS_LAYER_NAME, opacity=0.5)

            # Restore cached previews (unified predicted + per-cell-type seg/prob).
            pred_path = unified_predicted_labels_path(pixel_class_outdir)
            if pred_path.exists():
                try:
                    pred = np.asarray(load_zarr(pred_path))
                    if pred.shape == label_shape:
                        self.viewer.add_labels(
                            pred, name=UNIFIED_PREDICTED_LAYER_NAME, opacity=0.6, visible=False,
                        )
                except Exception as exc:
                    self.log(f"  ⚠️ Could not restore unified predicted labels: {exc}")

            for ct in all_cell_types:
                seg_p = _predicted_labels_path(pixel_class_outdir, ct)
                if seg_p and seg_p.exists():
                    try:
                        pred = np.asarray(load_zarr(seg_p))
                        if pred.shape == label_shape:
                            self.viewer.add_labels(
                                pred, name=_segments_layer_name(ct),
                                opacity=0.8, visible=False,
                            )
                    except Exception:
                        pass
                prob_p = _probability_map_path(pixel_class_outdir, ct)
                if prob_p and prob_p.exists():
                    try:
                        prob = np.asarray(load_zarr(prob_p))
                        if prob.shape == label_shape:
                            self.viewer.add_image(
                                prob, name=_probability_layer_name(ct),
                                opacity=0.6, blending="additive", colormap="magma",
                                contrast_limits=(0.0, 1.0), visible=False,
                            )
                    except Exception:
                        pass

            if has_death:
                seg_p = _predicted_labels_path(pixel_class_outdir, "dead")
                if seg_p and seg_p.exists():
                    try:
                        pred = np.asarray(load_zarr(seg_p))
                        if pred.shape == label_shape:
                            dead_pred_layer = self.viewer.add_labels(
                                pred, name=DEAD_PREDICTED_LAYER_NAME,
                                opacity=0.8, visible=False,
                            )
                            _set_dead_mask_layer_color(dead_pred_layer)
                    except Exception:
                        pass

            # Build (or rebuild) the training widget itself.
            if self._training_widget is not None:
                try:
                    self._training_widget.cleanup()
                except Exception:
                    pass
                self.training_layout.removeWidget(self._training_widget)
                self._training_widget.hide()
                self._training_widget.deleteLater()
                self._training_widget = None

            if self.training_placeholder.parent() is not None:
                self.training_layout.removeWidget(self.training_placeholder)
                self.training_placeholder.hide()

            params = self.metadata_loader.behav3d_parameters
            ip = dict(params.get("pixel_classifier", {}) or {})
            # Keep the user-selected device in the initial params so the
            # training widget picks the same one.
            ip["convpaint_device"] = self._selected_device()
            cp_strategy = _normalize_strategy(
                ip.get("convpaint_strategy", DEFAULT_STRATEGY)
            )

            if md is not None and "pixel_distance_xy" in md.columns and "distance_unit" in md.columns:
                _unit = str(md["distance_unit"].iloc[0])
                _xy_from_md = convert_distance(float(md["pixel_distance_xy"].iloc[0]), _unit)
                _z_from_md  = convert_distance(float(md["pixel_distance_z"].iloc[0]),  _unit)
            else:
                _xy_from_md = None
                _z_from_md  = None
            _pixel_sizes = {
                "xy_um": ip.get("pixel_size_xy") or ip.get("pixel_size_xy_um") or _xy_from_md,
                "z_um":  ip.get("pixel_size_z")  or ip.get("pixel_size_z_um")  or _z_from_md,
            }

            self._training_widget = ConvPaintTrainingWidget(
                viewer=self.viewer,
                all_images=all_images,
                pixel_class_outdir=str(pixel_class_outdir),
                all_cell_types=all_cell_types,
                has_death=has_death,
                initial_params=ip,
                on_params_changed=self._save_params_to_yaml,
                convpaint_strategy=cp_strategy,
                unified_input_channels=unified_input_channels,
                death_input_channels=death_input_channels,
                pixel_sizes=_pixel_sizes,
                per_cell_type_strategy=True,
                show_device=False,
                external_log=self.log,
            )
            if self._training_widget.strategy_combo is not None:
                self._training_widget.strategy_combo.setCurrentText(cp_strategy)

            # Wire signals.
            tw = self._training_widget
            tw.training_started.connect(self._on_training_started)
            tw.training_finished.connect(self._on_training_finished)
            tw.instance_preview_started.connect(self._on_instance_preview_started)
            tw.instance_preview_finished.connect(self._on_instance_preview_finished)
            tw.strategy_changed.connect(self._on_training_widget_strategy_changed)

            self.training_layout.addWidget(self._training_widget)

            self._all_images_cache = all_images
            self._is_session_active = True

            _reorder_convpaint_layers(self.viewer, all_cell_types, has_death=has_death)
            
            # Show AnnotationLegendTab
            while self.legend_container.layout().count():
                child = self.legend_container.layout().takeAt(0)
                if child.widget():
                    child.widget().deleteLater()
            legend = AnnotationLegendTab(self.viewer, label_map, has_death=has_death)
            self.legend_container.layout().addWidget(legend)
            self.legend_container.setVisible(True)
            
            self.log("✅ ConvPaint training data generated/loaded in viewer!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Error during ConvPaint training data generation: {e}")

    def cleanup_session(self):
        """Remove ConvPaint training layers and tear down the training widget.

        After clearing the viewer layers the parameter-only widget is restored
        (via :meth:`_on_metadata_updated`) so the user can still tune
        segmentation parameters without reloading training data.
        """
        n_removed = _remove_segmentation_layers(self.viewer)

        self._all_images_cache = None
        self._is_session_active = False

        # Delete elements inside legend_container and hide it
        while self.legend_container.layout().count():
            child = self.legend_container.layout().takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.legend_container.setVisible(False)

        # Rebuild the parameter-only widget so parameters remain tunable.
        self._on_metadata_updated()
        self.log(f"ConvPaint training session cleared ({n_removed} layers removed).")

    # ── Queue parameter snapshot ─────────────────────────────────
    def _current_global_strategy(self):
        if (
            self._training_widget is not None
            and self._training_widget.strategy_combo is not None
        ):
            return self._training_widget.strategy_combo.currentText()
        return self.all_strategies()[0]

    def validate_for_queue(self):
        """Return ``(ok, error_msg)`` for the queue-add path.

        Mirrors :meth:`APOCWidget.validate_for_queue` /
        :meth:`CellposeWidget.validate_for_queue`: refuses to enqueue a
        ConvPaint segmentation step that would immediately fail at run-time
        because metadata is missing, no output directory is set, or the
        unified ConvPaint model has not been trained yet.
        """
        md = self.metadata_loader.metadata
        if md is None:
            return False, "No metadata loaded. Please load metadata first (Data Preparation tab)."

        out_dir = self.metadata_loader.output_dir
        if not out_dir:
            return False, "No output directory established. Please set it in the Data Preparation tab."

        from behav3d.preprocessing.segmentation.convpaint_label_map import (
            unified_model_path,
        )
        pixel_class_outdir = Path(out_dir) / "images" / "PixelClassification"
        model_path = unified_model_path(pixel_class_outdir)
        if not model_path.exists():
            return False, (
                "No trained unified ConvPaint model found at\n"
                f"{model_path}\n\n"
                "Please train the classifier in the legend tab before queueing."
            )
        return True, ""

    def get_queue_params(self) -> dict:
        """Snapshot current widget state for use by the processing queue."""
        strategy = self._current_global_strategy()
        all_strats = self.all_strategies()
        params = {
            "strategy_index": all_strats.index(strategy) if strategy in all_strats else 0,
            "strategy_name": strategy,
            "device": self._selected_device(),
            "overwrite": self.check_overwrite.isChecked(),
            "workers": self.spin_workers.value(),
            "process_all": self.check_process_all.isChecked(),
            "t_start": self.spin_t_start.value(),
            "t_end": self.spin_t_end.value(),
        }
        if self._training_widget is not None:
            per_ct_params = {}
            per_ct_strategies = {}
            for ct, tab in self._training_widget.tabs.items():
                try:
                    per_ct_params[ct] = tab.collect_params()
                except Exception:
                    pass
                per_combo = getattr(tab, "_per_tab_strategy_combo", None)
                if per_combo is not None:
                    per_ct_strategies[ct] = per_combo.currentText()
            params["per_ct_params"] = per_ct_params
            from behav3d.preprocessing.segmentation.convpaint_train import (
                ConvPaintTrainingWidget as _CTW,
            )
            if strategy == _CTW.ADVANCED_STRATEGY:
                params["per_ct_strategies"] = per_ct_strategies
        return params

    # ── Parameter persistence ────────────────────────────────────
    def _save_params_to_yaml(self, updated_params=None):
        """Save ConvPaint-related parameters into ``behav3d_parameters.yml``.

        Accepts an optional ``updated_params`` dict from the training widget
        (``on_params_changed`` callback) — those are merged into the
        ``pixel_classifier`` section since ConvPaint uses ``convpaint_*``
        keys at that level.
        """
        out_dir = self.metadata_loader.output_dir
        if not out_dir:
            return

        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        if isinstance(updated_params, dict):
            pc.update(updated_params)

        if hasattr(self, "spin_examples") and self.spin_examples is not None:
            pc["examples_per_sample"] = self.spin_examples.value()
        if hasattr(self, "combo_gpu_device") and self.combo_gpu_device is not None:
            pc["convpaint_device"] = self._selected_device()
        if (
            self._training_widget is not None
            and self._training_widget.strategy_combo is not None
        ):
            pc["convpaint_strategy"] = self._training_widget.strategy_combo.currentText()
        if hasattr(self, "spin_workers") and hasattr(self, "check_process_all"):
            pc["convpaint_workers"] = self.spin_workers.value()
            pc["convpaint_process_all"] = self.check_process_all.isChecked()
            pc["convpaint_t_start"] = self.spin_t_start.value()
            pc["convpaint_t_end"] = self.spin_t_end.value()

        params_path = Path(out_dir) / "behav3d_parameters.yml"
        try:
            with open(params_path, "w") as f:
                yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
        except Exception as e:
            self.log(f"Warning: Could not save parameters: {e}")

    # ── Run batch segmentation ───────────────────────────────────
    def _on_run_segmentation(self, interactive=True, only_cell_types=None,
                             block=True, extra_callbacks=None):
        """Run ConvPaint batch segmentation.

        ``block=True`` (default) is synchronous so the processing queue
        sees unchanged behaviour.  ``block=False`` (GUI button) routes
        the work through :class:`BackgroundOperation` with an
        indeterminate busy spinner (the ConvPaint backend does not
        currently expose a per-sample progress hook).
        ``extra_callbacks`` is the queue's chaining hook.
        """
        try:
            md = self.metadata_loader.metadata
            if md is None:
                self.log("⚠️ No metadata loaded.")
                fire_extra_callback(extra_callbacks, "on_failed", "no metadata loaded")
                return

            if not block and self._bg.is_running():
                self.log("⚠️ A ConvPaint run is already in progress.")
                fire_extra_callback(extra_callbacks, "on_failed", "already running")
                return

            from behav3d.preprocessing.segmentation.convpaint_segment import (
                run_convpaint_segmentation,
            )

            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"

            # ---- Overwrite check (batch only) ---------------------------
            if interactive and only_cell_types is None:
                check_cts = list(self.all_cell_types or [])
                existing = []
                for sn in md["sample_name"].unique():
                    for ct in check_cts:
                        seg_path = output_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                        if seg_path.exists():
                            existing.append(f"{ct} segments for {sn}")
                if existing:
                    from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
                    choice = prompt_overwrite_batch(
                        self,
                        "Overwrite Existing ConvPaint Segmentations?",
                        existing,
                    )
                    if choice == "cancel":
                        self.log("ConvPaint segmentation cancelled.")
                        fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                        return
                    self.check_overwrite.setChecked(choice == "overwrite")

            # Quick preflight: a unified ConvPaint model must exist.
            from behav3d.preprocessing.segmentation.convpaint_label_map import (
                unified_model_path,
            )
            if not unified_model_path(pixel_class_outdir).exists():
                self.log(
                    "❌ No trained unified ConvPaint model found at "
                    f"{unified_model_path(pixel_class_outdir)}. Train the classifier first."
                )
                fire_extra_callback(extra_callbacks, "on_failed", "no trained model")
                return

            strategy = self._current_global_strategy()

            if self.check_process_all.isChecked():
                timepoint_range = None
            else:
                t_start = self.spin_t_start.value()
                t_end = self.spin_t_end.value()
                if t_start > t_end:
                    self.log("Error: Start timepoint must be <= End timepoint.")
                    fire_extra_callback(extra_callbacks, "on_failed", "bad timepoint range")
                    return
                timepoint_range = (t_start, t_end)

            # Collect convpaint_config from training widget tabs (Qt thread).
            convpaint_config = {}
            if self._training_widget is not None:
                for tab in self._training_widget.tabs.values():
                    try:
                        convpaint_config.update(tab.collect_params() or {})
                    except Exception:
                        pass

            self._save_params_to_yaml()

            per_ct_convpaint_strategies = None
            from behav3d.preprocessing.segmentation.convpaint_train import (
                ConvPaintTrainingWidget as _CTW,
            )
            if strategy == _CTW.ADVANCED_STRATEGY and self._training_widget is not None:
                per_ct_convpaint_strategies = {}
                for ct, tab in self._training_widget.tabs.items():
                    per_combo = getattr(tab, "_per_tab_strategy_combo", None)
                    per_ct_convpaint_strategies[ct] = (
                        per_combo.currentText() if per_combo is not None
                        else self.per_tab_strategies()[0]
                    )
                pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
                pc["per_ct_convpaint_strategies"] = per_ct_convpaint_strategies
                self._save_params_to_yaml()
                self.log(
                    "  Advanced strategies: "
                    + ", ".join(
                        f"{k}:{v[:4]}" for k, v in per_ct_convpaint_strategies.items()
                    )
                )

            self.log(f"Starting ConvPaint batch segmentation ({strategy})...")

            metadata_csv = self.metadata_loader.behav3d_parameters.get("paths", {}).get(
                "metadata_csv", ""
            )

            convpaint_config["convpaint_device"] = self._selected_device()

            overwrite_val = self.check_overwrite.isChecked()
            n_workers_val = self.spin_workers.value()

            def _do_convpaint(progress_cb=None):
                return run_convpaint_segmentation(
                    output_dir=str(output_dir),
                    metadata=md,
                    metadata_csv_path=metadata_csv,
                    timepoint_range=timepoint_range,
                    convpaint_config=convpaint_config,
                    overwrite_existing=overwrite_val,
                    n_workers=n_workers_val,
                    convpaint_strategy=strategy,
                    per_ct_convpaint_strategies=per_ct_convpaint_strategies,
                    only_cell_types=only_cell_types,
                )

            def _apply_result(updated_metadata):
                if updated_metadata is not None:
                    self.metadata_loader.metadata = updated_metadata
                    updated_metadata = _run_multicolor_cleanup_if_needed(self.metadata_loader, self.log, skip_unified=True)
                    self.metadata_loader.metadata = updated_metadata
                    if metadata_csv:
                        try:
                            updated_metadata.to_csv(metadata_csv, index=False)
                        except Exception as e:
                            self.log(f"Warning: Could not save metadata: {e}")
                self.log("✅ ConvPaint batch segmentation finished!")
                if interactive:
                    self._prompt_visualize_after_segmentation()

            if block:
                try:
                    result = _do_convpaint(progress_cb=None)
                    _apply_result(result)
                    fire_extra_callback(extra_callbacks, "on_done", result)
                except Exception as e:
                    traceback.print_exc()
                    self.log(f"❌ Error during ConvPaint segmentation: {e}")
                    fire_extra_callback(extra_callbacks, "on_failed", str(e))
                return

            def _on_done_async(result):
                _apply_result(result)
                fire_extra_callback(extra_callbacks, "on_done", result)

            def _on_failed(err: str):
                self.log(f"❌ Error during ConvPaint segmentation: {err}")
                fire_extra_callback(extra_callbacks, "on_failed", err)

            self._bg.run(
                fn=_do_convpaint,
                desc=f"ConvPaint segmentation ({strategy})\u2026",
                progress_row=self.tab_progress_row,
                buttons=[],
                viewer=self.viewer,
                on_done=_on_done_async,
                on_failed=_on_failed,
                inject_progress=False,
                indeterminate=True,
            )

        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Error during ConvPaint segmentation: {e}")
            fire_extra_callback(extra_callbacks, "on_failed", str(e))
