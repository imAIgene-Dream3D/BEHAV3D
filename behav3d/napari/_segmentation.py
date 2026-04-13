
import napari
from behav3d.napari._widgets import HelpButton, make_help_row
import yaml
from magicgui.widgets import create_widget
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, 
    QStackedWidget, QPushButton, QGroupBox, QFormLayout, 
    QSpinBox, QDoubleSpinBox, QCheckBox, QFileDialog, QScrollArea,
    QLineEdit, QPlainTextEdit, QTextEdit, QMessageBox
)
from qtpy.QtCore import Qt
import pandas as pd
from pathlib import Path
import traceback
import numpy as np
import dask.array as da
import shutil
import gc
import os
from functools import partial
from skimage import feature
# Scikit-learn imports for training
from sklearn.ensemble import RandomForestClassifier
import joblib


from behav3d.io.images import load_image, get_image_shape, save_as_zarr, append_to_zarr, load_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.preprocessing.segmentation.napari_pixelclassifier import (
    train_pixel_classifier, 
    run_pixel_classifier_segmentation,
    segment_mask,
    postprocess_mask
)


# Feature calculation function (same as in napari_pixelclassifier.py)
sigma_min = 1
sigma_max = 16
features_func = partial(
        feature.multiscale_basic_features,
        intensity=True,
        edges=True,
        texture=False,
        sigma_max=sigma_max,
        channel_axis=0,
    )

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

        # Method Selection
        method_group = QGroupBox("Segmentation Method")
        method_layout = QHBoxLayout()
        self.method_combo = QComboBox()
        self.method_combo.addItems([
            "APOC (GPU)",
            "Pixel Classifier (Random Forest)",
            "Cellpose (Deep Learning)",
            "Import segmentation",
        ])
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.method_combo)
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # Stacked Widget for Method Parameters
        self.param_stack = QStackedWidget()
        
        # 0. APOC (GPU) Page  ← matches combo index 0
        self.apoc_page = APOCWidget(self.viewer, self.metadata_loader, log_callback=self._log)
        self.param_stack.addWidget(self.apoc_page)

        # 1. Pixel Classifier Page  ← matches combo index 1
        self.pixel_classifier_page = PixelClassifierWidget(self.viewer, self.metadata_loader, log_callback=self._log)
        self.param_stack.addWidget(self.pixel_classifier_page)

        # 2. Cellpose Page  ← matches combo index 2
        self.cellpose_page = CellposeWidget(self.viewer, self.metadata_loader, log_callback=self._log)
        self.param_stack.addWidget(self.cellpose_page)

        # 3. Import Page  ← matches combo index 3
        self.import_page = ImportWidget(self.viewer, self.metadata_loader, log_callback=self._log)
        self.param_stack.addWidget(self.import_page)


        layout.addWidget(self.param_stack)
        
        # Add stretch to keep things at the top
        # layout.addStretch()


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

    def _log(self, msg):
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        self.log.append(f"[{timestamp}] {msg}")
        self.log.verticalScrollBar().setValue(
            self.log.verticalScrollBar().maximum()
        )

    def _on_method_changed(self, index):
        current_idx = self.param_stack.currentIndex()

        # If switching away from Pixel Classifier (index 1), check session
        if current_idx == 1 and self.pixel_classifier_page.is_session_active:
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

        self.param_stack.setCurrentIndex(index)

        # If switching to APOC tab (index 0), refresh it
        if index == 0:
            self.apoc_page._on_metadata_updated()
        # If switching to Import tab (index 3), refresh it to avoid outdated info
        elif index == 3:
            self.import_page._on_metadata_updated()

    def request_tab_exit(self):
        """Called by the main widget when the user tries to leave this tab."""
        from qtpy.QtWidgets import QMessageBox

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
            else:
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
            else:
                self.apoc_page.cleanup_session()
                return True

        return True

    def cleanup_session(self):
        """Clear viewer and reset tab state (both Pixel Classifier and APOC)."""
        self.pixel_classifier_page.reset_ui()
        self.apoc_page.cleanup_session()

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    """Get a value from a nested dict using dotted key notation."""
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


class PixelClassifierWidget(QWidget):
    def __init__(self, viewer, metadata_loader, log_callback=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback if log_callback else print
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
            has_dead_channel
        )
        metadata = self.metadata_loader.metadata
        
        if metadata is None:
            self.organoid_types = []
            self.immune_types = []
            self.other_types = []
            self.has_death = False
            self.all_cell_types = []
            return

        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _init_ui(self):
        # Connect to metadata updates
        if hasattr(self.metadata_loader, 'metadata_loaded'):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

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
            "Number of timepoints randomly sampled from each dataset "
            "for training the pixel classifier.\n\n"
            "More examples → better generalization but slower training."
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

        # Legend for labels
        legend_layout = QHBoxLayout()
        legend_layout.setContentsMargins(6, 4, 6, 4)
        
        def create_legend_item(text, color_hex):
            item = QHBoxLayout()
            box = QLabel()
            box.setFixedSize(12, 12)
            box.setStyleSheet(f"background-color: {color_hex}; border: 1px solid #555;")
            label = QLabel(text)
            label.setStyleSheet("font-size: 10px; font-weight: bold;")
            item.addWidget(box)
            item.addWidget(label)
            return item

        legend_layout.addLayout(create_legend_item("0: Eraser", "none"))
        legend_layout.addLayout(create_legend_item("1: Background", "#8b3a26")) # Red
        legend_layout.addLayout(create_legend_item("2: Foreground", "#00ffff")) # Cyan
        legend_layout.addStretch()
        
        self.legend_widget = QWidget()
        self.legend_widget.setLayout(legend_layout)
        self.legend_widget.setVisible(False) # Hide until data loaded
        train_layout.addWidget(self.legend_widget)

        # Add training related buttons
        train_layout.addLayout(train_btn_row)

        # Hide action buttons until training data is loaded
        self.btn_save_labels.setVisible(False)
        self.btn_clear_layer.setVisible(False)
        self.btn_clear_all.setVisible(False)
        self.btn_train_update.setVisible(False)
        self.legend_widget.setVisible(False)
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
            w_edt.setRange(0.0, 50.0)
            w_edt.setSingleStep(0.5)
            w_edt.setValue(float(saved_edt))
            w_edt.setMaximumWidth(70)
            ct_form.addRow("EDT:", make_help_row(
                w_edt,
                "EDT (Distance Transform Threshold)",
                "Controls sensitivity for separating touching objects.\n\n"
                "Lower values → more sensitive, splits objects more aggressively.\n"
                "Higher values → less splitting, objects stay merged.\n\n"
                "Typical values:\n"
                "  • Organoids: 8–15\n"
                "  • Immune cells: 1.5–4.0"
                "  • Disable: 0.0 (not recommended, leads to under-segmentation)"
            ))

            w_size = QSpinBox()
            w_size.setRange(1, 100000)
            w_size.setValue(int(saved_size))
            w_size.setMaximumWidth(80)
            ct_form.addRow("Min size:", make_help_row(
                w_size,
                "Minimum Object Size (pixels)",
                "Minimum volume (in pixels) for a segmented object to be kept.\n\n"
                "Objects smaller than this are removed as noise.\n\n"
                "Typical values:\n"
                "  • Organoids: 500–2000\n"
                "  • Immune cells: 5–50"
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
            ct_form.addRow(make_help_row(
                w_fill,
                "Fill Holes",
                "When checked, internal gaps or holes within segmented "
                "objects are automatically filled.\n\n"
                "Recommended ON for most use cases."
            ))

            ct_group.setLayout(ct_form)
            self.param_layout.addWidget(ct_group)

            self.param_widgets[cell_type] = {
                'edt': w_edt,
                'min_size': w_size,
                'opening': w_open,
                'fill_holes': w_fill
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
                pc[f"{cell_type}_edt_threshold"] = float(w['edt'].value())
                pc[f"{cell_type}_segment_size_min"] = int(w['min_size'].value())
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
        self.legend_widget.setVisible(False)
        self.instruction_label.setVisible(False)
        self.btn_queue_train.setVisible(False)
        
        self.all_images = None
        self.all_features = None
        self.is_session_active = False
        
        # Clear specific layers related to segmentation/labeling if they exist
        layers_to_remove = [l for l in self.viewer.layers if 'User Provided Labels' in l.name or 'Pixel Classification' in l.name or 'Channel' in l.name or 'Segments' in l.name]
        for l in layers_to_remove:
            self.viewer.layers.remove(l)
        
        self.log("Segmentation session reset.")

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
        self.run_batch_segmentation(interactive=True)

    def run_batch_segmentation(self, interactive=True, skip_existing=False):
        """Run batch segmentation. When interactive=False, skips dialogs."""
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run segmentation: No metadata loaded.")
            return

        self.log("Starting batch segmentation...")
        self._persist_params()
        try:
            if interactive:
                from qtpy.QtWidgets import QMessageBox
                box = QMessageBox(self)
                box.setWindowTitle("Existing Segmentation Results")
                box.setText("Segmentation results may already exist for some timepoints.\n\nWhat do you want to do?")
                btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
                btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
                btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                box.setDefaultButton(btn_cancel)
                box.exec_()
                clicked = box.clickedButton()
                if clicked == btn_cancel:
                    self.log("Segmentation cancelled.")
                    return
                overwrite = (clicked == btn_overwrite)
            else:
                overwrite = not skip_existing
            
            # Timepoint range
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
                    edt = float(w['edt'].value())
                    sz  = int(w['min_size'].value())
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

            # Call the segmentation function
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
                n_workers=self.spin_workers.value(),
                log_callback=self.log,
                overwrite_existing=overwrite,
                timepoint_range=timepoint_range,
                clf_organoid_paths=clf_organoid_paths,
                clf_immune_paths=clf_immune_paths,
                clf_other_paths=clf_other_paths,
                clf_death_path=clf_death_path
            )
            # Update metadata in loader
            self.metadata_loader.metadata = updated_metadata
            
            # Save updated metadata CSV to disk
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                try:
                    updated_metadata.to_csv(csv_path, index=False)
                except Exception as e:
                    self.log(f"  Warning: Could not save metadata CSV: {e}")

            self.log("Batch segmentation finished successfully!")

            # Prompt the user if they want to visualize results (only in interactive mode)
            if interactive:
                from qtpy.QtWidgets import QMessageBox
                res = QMessageBox.question(
                    self,
                    "Segmentation Finished",
                    "Batch segmentation finished successfully! \n\nDo you want to switch to the Visualization Tab and see the segments?",
                    QMessageBox.Yes | QMessageBox.No
                )
                
                if res == QMessageBox.Yes:
                    # Trigger switch to visualization tab and load first sample
                    parent = self.parent()
                    while parent and not hasattr(parent, 'tabs'):
                        parent = parent.parent()
                    
                    if parent and hasattr(parent, 'tabs'):
                        # Switch to visualization tab (index 1)
                        parent.tabs.setCurrentIndex(1)
                        # Load first sample
                        if hasattr(parent, 'visualization_tab'):
                            self.log("  Loading first sample in Visualization Tab...")
                            parent.visualization_tab.sample_combo.setCurrentIndex(0)
                            parent.visualization_tab._on_load_dataset()
                            
                            # Ensure all 'Segments' layers are visible
                            for layer in self.viewer.layers:
                                if "Segments" in layer.name:
                                    layer.visible = True

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during batch segmentation: {e}")

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
                
                processed_images = []
                processed_features = []
                
                for i, img in enumerate(all_images):
                    self.log(f"Processing image {i+1}/{len(all_images)}")
                    try:
                        feats = features_func(img)
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
                self.viewer.add_labels(pred_dead, name='Pixel Classification (Dead)', opacity=0.3, visible=False)

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
            self.legend_widget.setVisible(True)
            self.instruction_label.setVisible(True)
            self.btn_queue_train.setVisible(True)  # Show queue button now that training data is loaded
            self.is_session_active = True
            
        except Exception as e:
            traceback.print_exc()
            self.log(f"Failed to load training data: {e}")

    def run_train(self, interactive=True):
        """Run classifier training. Called by queue (interactive=False) or button (interactive=True)."""
        if not self.is_session_active:
            self.log("Session not active. Loading training data...")
            self._on_load_training_clicked(interactive=interactive)
        
        if self.is_session_active:
            self._on_train_clicked()
        else:
            self.log("Error: Failed to load training data for classifier.")
            if not interactive:
                 raise RuntimeError("Training data not loaded. Please load training data manually first.")

    def _on_train_clicked(self):
        """Train classifiers, predict pixels, and segment — mirrors original segment_and_update()."""
        self.log("Training classifier & segmenting...")
        self._persist_params()
        self._save_user_labels()  # Always persist labels before training
        
        # ── Load features ───────────────────────────────────────────
        if self.all_features is None:
            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            features_outpath = pixel_class_outdir / 'PixelClassifier_Features.zarr'
            if features_outpath.exists():
                self.log("Loading features from disk...")
                self.all_features = np.asarray(load_zarr(features_outpath))
            else:
                self.log("Error: No features found. Please load training data first.")
                return

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)
        n_workers = self.spin_workers.value()

        # all_features shape: (T, Z, Y, X, C_feat) — channel LAST
        all_features = np.asarray(self.all_features)
        self.log(f"Features shape: {all_features.shape}")

        # ── Collect targets ─────────────────────────────────────────
        targets = []
        if self.has_death:
            targets.append('Dead')
        targets.extend([t.capitalize() for t in self.all_cell_types])

        try:
            # ── 1. Save & postprocess user labels, then train ───────
            classifiers = {}
            cell_type_labels = {}

            for target in targets:
                label_layer_name = f'User Provided Labels ({target})'
                if label_layer_name not in self.viewer.layers:
                    self.log(f"  Layer {label_layer_name} not found. Skipping.")
                    continue
                
                labels = self.viewer.layers[label_layer_name].data.copy()

                if labels.max() == 0:
                    self.log(f"  No labels for {target}. Skipping.")
                    continue

                # Get postprocessing params for this target
                cell_type_key = target.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    opening = int(w['opening'].value())
                    fill = bool(w['fill_holes'].isChecked())
                elif cell_type_key == 'dead':
                    opening, fill = 0, True
                else:
                    opening, fill = 0, True

                # Postprocess user labels before training (match original lines 594-618)
                fg_mask = (labels == 2).astype(bool)
                if np.any(fg_mask):
                    processed_fg = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                    labels = np.where(processed_fg, 2,
                                      np.where(labels == 1, 1, 0)).astype(labels.dtype)

                # Save user labels to zarr
                labels_outpath = pixel_class_outdir / f'PixelClassifier_User{target}Labels.zarr'
                if labels_outpath.exists():
                    shutil.rmtree(labels_outpath)
                save_as_zarr(labels, labels_outpath)
                cell_type_labels[target] = labels
                self.log(f"  Saved {target} user labels (postprocessed)")

                # Train classifier (match original train_classifier)
                self.log(f"  Training {target}...")
                flat_labels = labels.ravel()
                flat_features = all_features.reshape(-1, all_features.shape[-1])

                label_indices = np.flatnonzero(flat_labels > 0)
                selected_features = flat_features[label_indices]
                selected_labels = flat_labels[label_indices]

                nr_bg = int(np.sum(selected_labels == 1))
                nr_fg = int(np.sum(selected_labels == 2))
                total = nr_bg + nr_fg
                self.log(f"    BG pixels: {nr_bg}, FG pixels: {nr_fg}")

                if total == 0:
                    self.log(f"  No labeled pixels for {target}. Skipping.")
                    continue

                class_weights = {1: nr_bg / total, 2: nr_fg / total}

                clf = RandomForestClassifier(
                    n_estimators=50,
                    n_jobs=-1,
                    max_depth=20,
                    class_weight=class_weights
                )
                from skimage import future
                clf = future.fit_segmenter(selected_labels, selected_features, clf)

                clf_path = pixel_class_outdir / f'PixelClassifier_{target}.joblib'
                joblib.dump(clf, clf_path)
                self.log(f"  Saved classifier: {clf_path.name}")
                classifiers[target] = clf

            # ── 2. Predict all pixels ───────────────────────────────
            pred_masks = {}
            n_timepoints = all_features.shape[0]

            for target, clf in classifiers.items():
                self.log(f"  Predicting {target} pixels...")
                prediction_stack = []

                for t in range(n_timepoints):
                    # all_features[t] shape: (Z, Y, X, C_feat) — already channel-last
                    feat_t = all_features[t]
                    spatial_shape = feat_t.shape[:-1]  # (Z, Y, X)
                    X_t = feat_t.reshape(-1, feat_t.shape[-1])

                    # Predict in batches
                    batch_size = 100_000
                    y_pred_t = np.zeros(X_t.shape[0], dtype=np.uint8)
                    for i in range(0, X_t.shape[0], batch_size):
                        y_pred_t[i:i+batch_size] = clf.predict(X_t[i:i+batch_size])

                    y_pred_t = y_pred_t.reshape(spatial_shape)
                    prediction_stack.append(y_pred_t)

                full_prediction = np.stack(prediction_stack)  # (T, Z, Y, X)

                # Postprocess predictions (match original lines 728-738)
                cell_type_key = target.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    opening = int(w['opening'].value())
                    fill = bool(w['fill_holes'].isChecked())
                elif cell_type_key == 'dead':
                    opening, fill = 0, True
                else:
                    opening, fill = 0, True

                fg_mask = (full_prediction == 2).astype(bool)
                processed_fg = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                full_prediction = np.where(processed_fg, 2,
                                           np.where(full_prediction == 1, 1, 0)).astype(np.uint8)

                pred_masks[target] = full_prediction

                # Save predictions to zarr
                pred_outpath = pixel_class_outdir / f'PixelClassifier_{target}_PredictedLabels.zarr'
                if pred_outpath.exists():
                    shutil.rmtree(pred_outpath)
                save_as_zarr(full_prediction, pred_outpath)

                # Update viewer layer
                pred_layer_name = f'Pixel Classification ({target})'
                if pred_layer_name in self.viewer.layers:
                    self.viewer.layers[pred_layer_name].data = full_prediction
                    # self.viewer.layers[pred_layer_name].visible = True
                    self.viewer.layers[pred_layer_name].refresh()
                else:
                    self.viewer.add_labels(full_prediction, name=pred_layer_name, opacity=0.3)
                self.log(f"  {target} prediction updated.")

            # ── 3. Segment instances (match original lines 757-819) ──
            self.log("Segmenting cell instances...")
            for cell_type in self.all_cell_types:
                target = cell_type.capitalize()
                if target not in pred_masks:
                    continue

                cell_type_key = cell_type.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    edt_threshold = float(w['edt'].value())
                    segment_size_min = int(w['min_size'].value())
                else:
                    if cell_type in self.organoid_types:
                        edt_threshold, segment_size_min = 12.0, 1000
                    elif cell_type in self.immune_types:
                        edt_threshold, segment_size_min = 2.5, 10
                    else:
                        edt_threshold, segment_size_min = 1.0, 10

                self.log(f"  Segmenting {target} (EDT={edt_threshold}, min_size={segment_size_min})...")
                pred_mask = pred_masks[target]

                segmented_timepoints = []
                for t_idx in range(pred_mask.shape[0]):
                    mask_t = (pred_mask[t_idx] == 2)
                    fg_pixels = int(mask_t.sum())

                    if fg_pixels == 0:
                        segmented = np.zeros_like(mask_t, dtype=np.uint16)
                    else:
                        segmented = segment_mask(
                            mask=mask_t,
                            edt_thr=edt_threshold,
                            edt_thr_refined=None,
                            segment_size_min=segment_size_min,
                            use_dims=2,
                            n_workers=1,
                        ).astype(np.uint16, copy=False)

                    segmented_timepoints.append(segmented)

                full_seg = np.stack(segmented_timepoints, axis=0)
                self.log(f"  {target}: max label = {int(full_seg.max())}")

                seg_layer_name = f'{target} Segments'
                if seg_layer_name in self.viewer.layers:
                    self.viewer.layers[seg_layer_name].data = full_seg
                    self.viewer.layers[seg_layer_name].visible = True
                    self.viewer.layers[seg_layer_name].refresh()
                else:
                    self.viewer.add_labels(full_seg, name=seg_layer_name, opacity=0.7)

            self.log("Training, prediction, and segmentation complete!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during training/segmentation: {e}")

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
                        edt = widgets['edt'].value()
                        size = widgets['min_size'].value()
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

    def __init__(self, viewer, metadata_loader, log_callback=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.pretrained_model_dir = None
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
        )
        metadata = self.metadata_loader.metadata
        if metadata is None:
            self.organoid_types = []
            self.immune_types = []
            self.other_types = []
            self.has_death = False
            self.all_cell_types = []
            return
        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
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
        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

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
        self.btn_run_cellpose.clicked.connect(lambda: self.run_batch_cellpose(interactive=True))
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
        self.btn_run_otsu.clicked.connect(lambda: self.run_otsu_threshold(interactive=True))
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
                           cell_type_override=None, model_path_override=None):
        """Run Cellpose segmentation for the selected cell type."""
        from behav3d.preprocessing.segmentation.cellpose_prediction import (
            run_cellpose_and_sync_metadata,
        )
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run Cellpose: No metadata loaded.")
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
            box = QMessageBox(self)
            box.setWindowTitle("Existing Cellpose Results")
            box.setText(
                f"Cellpose segmentation for '{cell_type}' may already exist.\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == btn_cancel:
                self.log("Cellpose segmentation cancelled.")
                return
            overwrite = clicked == btn_overwrite
        else:
            overwrite = not skip_existing

        # Timepoint range
        if self.check_process_all.isChecked():
            timepoint_range = None
        else:
            t_start = self.spin_t_start.value()
            t_end = self.spin_t_end.value()
            if t_start > t_end:
                self.log("Error: Start timepoint must be <= End timepoint.")
                return
            timepoint_range = (t_start, t_end)

        self.log(f"Starting Cellpose segmentation for '{cell_type}'...")

        try:
            _, summary = run_cellpose_and_sync_metadata(
                output_dir=self.metadata_loader.output_dir,
                metadata_loader=self.metadata_loader,
                pretrained_model_dir=model_path,
                input_channels=[cell_type],
                label_name=cell_type,
                timepoint_range=timepoint_range,
            )
            n_proc = len(summary["processed"])
            n_skip = len(summary["skipped"])
            if n_skip > 0:
                self.log(f"Cellpose finished: {n_proc} processed, {n_skip} skipped.")
            else:
                self.log(f"Cellpose finished for all {n_proc} samples.")
            self.log("Metadata updated.")
        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Cellpose error: {e}")
            if not interactive:
                raise

    # ── run Otsu dead mask ──────────────────────────────────────────────
    def run_otsu_threshold(self, interactive=True):
        """Run Otsu thresholding on the dead channel."""
        from behav3d.preprocessing.segmentation.cellpose_prediction import (
            run_otsu_threshold_segmentation_from_zarr,
        )
        if self.metadata_loader.metadata is None:
            self.log("⚠️ Cannot run Otsu: No metadata loaded.")
            return

        self._persist_channel_config()
        self.log("Starting Otsu dead mask segmentation...")

        # Timepoint range
        if self.check_process_all.isChecked():
            timepoint_range = None
        else:
            timepoint_range = (self.spin_t_start.value(), self.spin_t_end.value())

        try:
            updated_metadata, summary = run_otsu_threshold_segmentation_from_zarr(
                output_dir=self.metadata_loader.output_dir,
                metadata=self.metadata_loader.metadata,
                timepoint_range=timepoint_range,
            )
            self.metadata_loader.metadata = updated_metadata
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                updated_metadata.to_csv(csv_path, index=False)
            n_proc = len(summary["processed"])
            self.log(f"Otsu dead mask finished: {n_proc} samples processed.")
        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Otsu error: {e}")


class ImportWidget(QWidget):
    """Widget to validate and import pre-existing segmentation files.
    
    Shows a per-sample, per-cell-type status table with conversion actions.
    """

    def __init__(self, viewer, metadata_loader, log_callback=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.organoid_types = []
        self.immune_types = []
        self.other_types = []
        self.all_cell_types = []
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

    def _resolve_path(self, path_str):
        """Resolve a path string, trying the metadata CSV directory if not absolute/found."""
        if not path_str:
            return None
        p = Path(path_str)
        if p.exists():
            return p
        
        # Try relative to metadata CSV
        md_path = getattr(self.metadata_loader, "_loaded_csv_path", None)
        if md_path:
            p_rel = Path(md_path).parent / path_str
            if p_rel.exists():
                return p_rel
        return p # Return original if still not found (will trigger 'File not found' UI)

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
        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

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
        )
        md = self.metadata_loader.metadata
        if md is None:
            return
        self.organoid_types = detect_organoid_types_from_metadata(md)
        self.immune_types = detect_immune_cell_types_from_metadata(md)
        self.other_types = detect_other_cell_types_from_metadata(md)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types
        self._rebuild_table()

    # ── table builder ───────────────────────────────────────────────────
    def _rebuild_table(self):
        """Clear and rebuild the entire status table."""
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
            if isinstance(w, QWidget):
                # Look for buttons with specific text/tooltips
                btns = w.findChildren(QPushButton)
                for b in btns:
                    if "Convert" in b.text() or "Fix" in b.text():
                        any_fixes_needed = True
                        break
            if any_fixes_needed:
                break

        if any_fixes_needed:
            # ── global Convert All button ──
            btn_all = QPushButton("⚡  Convert / Fix All Samples")
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

        # ── General Instruction (ALWAYS visible at the bottom) ──
        instr = QLabel(
            "To add or change segmentation paths, go to the <b>Data Preparation</b> tab "
            "and check the <b>Metadata Builder</b>."
        )
        instr.setWordWrap(True)
        instr.setAlignment(Qt.AlignCenter)
        instr.setStyleSheet("color:#555; background:#f9f9f9; padding:15px; border-radius:8px; margin:10px; border:1px solid #ddd;")
        self.scroll_layout.addWidget(instr)

        self.scroll_layout.addStretch()

    def _add_sample_section(self, sample_name, row_idx, row):
        """Add a grouped section for one sample."""
        from functools import partial

        # Sample header
        header = QLabel(f"📁  {sample_name}")
        header.setStyleSheet("font-weight:bold; font-size:13px; padding:6px 0 2px 0;")
        self.scroll_layout.addWidget(header)

        # Track whether ALL cell types have empty paths
        all_empty = True

        for ct in self.all_cell_types:
            col = self._seg_col(ct)
            raw_val = row.get(col) if col in row.index else None
            has_value = raw_val is not None and pd.notna(raw_val) and str(raw_val).strip() not in ("", "nan")

            if has_value:
                all_empty = False

            row_widget = QWidget()
            row_lay = QHBoxLayout(row_widget)
            row_lay.setContentsMargins(16, 2, 4, 2)
            lbl = QLabel(f"{ct}:")
            lbl.setFixedWidth(120)
            row_lay.addWidget(lbl)

            if not has_value:
                status = QLabel("No segmentation available")
                status.setStyleSheet("color:#999; font-style:italic;")
                row_lay.addWidget(status)
            else:
                path_str = str(raw_val).strip().strip('"').strip("'")
                file_path = self._resolve_path(path_str)

                if not file_path.exists():
                    status = QLabel(f"⚠️  File not found")
                    status.setToolTip(str(file_path))
                    status.setStyleSheet("color:#E65100;")
                    row_lay.addWidget(status)

                elif file_path.suffix == ".zarr":
                    ok, reason = self._check_zarr_structure(file_path)
                    if ok:
                        # --- Dimension check ---
                        dims_match = True
                        try:
                            raw_path = Path(row.get("raw_image_path", ""))
                            if raw_path.exists():
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
                        except:
                            pass

                        if dims_match:
                            status = QLabel("✅  Ready for tracking")
                            status.setStyleSheet("color:#2E7D32; font-weight:bold;")
                            row_lay.addWidget(status)
                        else:
                            status = QLabel(f"⚠️  Dimension mismatch")
                            status.setToolTip(reason)
                            status.setStyleSheet("color:#E65100; font-weight:bold;")
                            row_lay.addWidget(status)
                    else:
                        btn = QPushButton("🔄  Fix zarr format")
                        btn.setToolTip(f"Issue: {reason}")
                        btn.setStyleSheet(
                            "QPushButton{background:#F57C00;color:white;padding:4px 10px;"
                            "border-radius:3px}"
                            "QPushButton:hover{background:#FB8C00}"
                        )
                        btn.clicked.connect(partial(self._convert_single, path_str, ct, sample_name, row_idx))
                        row_lay.addWidget(btn)

                elif file_path.suffix.lower() in (".tif", ".tiff"):
                    # Check TIFF dims if possible
                    warning = ""
                    try:
                        raw_path = Path(row.get("raw_image_path", ""))
                        if raw_path.exists():
                            import zarr
                            from behav3d.core.io import load_image
                            raw_store = zarr.open(str(raw_path), mode='r')
                            # For TIFF we'd have to load it to check shape, which is slow for many files.
                            # Let's just do it and log if mismatch.
                            # Actually, maybe just keep it simple and check AFTER conversion? 
                            # No, user wants check during import.
                    except:
                        pass

                    btn = QPushButton("🔄  Convert TIFF → zarr")
                    btn.setStyleSheet(
                        "QPushButton{background:#1565C0;color:white;padding:4px 10px;"
                        "border-radius:3px}"
                        "QPushButton:hover{background:#1976D2}"
                    )
                    btn.clicked.connect(partial(self._convert_single, path_str, ct, sample_name, row_idx))
                    row_lay.addWidget(btn)
                else:
                    status = QLabel(f"⚠️  Format not supported ({file_path.suffix})")
                    status.setStyleSheet("color:#E65100;")
                    row_lay.addWidget(status)

            row_lay.addStretch()
            self.scroll_layout.addWidget(row_widget)

        # If ALL cell types empty → show combined message instead
        if all_empty:
            # remove the per-cell-type rows we just added (they all say "No segmentation")
            # keep the header, remove the rest
            items_to_remove = []
            for i in range(self.scroll_layout.count() - 1, -1, -1):
                item = self.scroll_layout.itemAt(i)
                w = item.widget()
                if w is header:
                    break
                items_to_remove.append(i)
            for i in items_to_remove:
                item = self.scroll_layout.takeAt(i)
                if item.widget():
                    item.widget().deleteLater()

            msg = QLabel("  No segmentation data found for this sample.<br>  (Check the <b>Metadata Builder</b> in Data Prep if you have files to import)")
            msg.setStyleSheet("color:#888; font-style:italic; padding:4px 16px;")
            msg.setWordWrap(True)
            self.scroll_layout.addWidget(msg)
        else:
            # Check if any "action" buttons were added for this sample
            fixes_needed = False
            # We look at the widgets added after the header
            found_header = False
            for i in range(self.scroll_layout.count()):
                w = self.scroll_layout.itemAt(i).widget()
                if not found_header:
                    if w is header:
                        found_header = True
                    continue
                # Now we are past the header
                if isinstance(w, QWidget):
                    if w.findChildren(QPushButton):
                        fixes_needed = True
                        break
            
            if fixes_needed:
                # Per-sample Convert All button
                btn = QPushButton(f"⚡  Convert / Fix all for {sample_name}")
                btn.setStyleSheet(
                    "QPushButton{background:#2E7D32;color:white;padding:5px 12px;"
                    "border-radius:3px;font-size:12px}"
                    "QPushButton:hover{background:#388E3C}"
                )
                btn.clicked.connect(partial(self._convert_sample, sample_name, row_idx))
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

    # ── conversion logic ────────────────────────────────────────────────
    def _convert_single(self, src_path_str, cell_type, sample_name, row_idx, 
                        save_metadata=True, refresh_ui=True):
        """Convert a single file (TIFF or bad zarr) and update metadata."""
        from qtpy.QtWidgets import QMessageBox
        src = Path(src_path_str)
        dest = self._expected_outpath(sample_name, cell_type)

        # Overwrite check
        if dest.exists():
            # If we are in a batch operation (save_metadata=False), we might want to skip prompt
            # but for safety let's prompt or assume user already confirmed global action?
            # Actually, per-file confirmation is safer unless we add "Apply to all".
            res = QMessageBox.question(
                self, "Overwrite?",
                f"There is a pre-existing segmentation file:\n{dest}\n\nDo you want to overwrite it?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                self.log(f"Skipped {cell_type} for {sample_name} (user declined overwrite)")
                return False

        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            self.log(f"Converting {src.name} → {dest.name} ...")
            img = load_image(src)
            if dest.exists():
                import shutil
                shutil.rmtree(dest)
            save_as_zarr(img, dest)
            self.log(f"✅  Saved: {dest}")

            # Update metadata
            md = self.metadata_loader.metadata
            col = self._seg_col(cell_type)
            if col not in md.columns:
                md[col] = pd.NA
            md[col] = md[col].astype("object")
            md.at[row_idx, col] = str(dest)
            
            if save_metadata:
                self._save_metadata(refresh_ui=refresh_ui)
            return True

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.log(f"❌  Error converting {cell_type} for {sample_name}: {e}")
            if refresh_ui: # Only show error dialog if not in background bulk operation
                QMessageBox.warning(self, "Conversion Error", str(e))
            return False

    def _convert_sample(self, sample_name, row_idx, save_metadata=True, refresh_ui=True):
        """Convert all actionable files for one sample."""
        md = self.metadata_loader.metadata
        row = md.iloc[row_idx]
        converted_any = False

        for ct in self.all_cell_types:
            col = self._seg_col(ct)
            raw_val = row.get(col) if col in row.index else None
            has_value = raw_val is not None and pd.notna(raw_val) and str(raw_val).strip() not in ("", "nan")
            if not has_value:
                continue

            path_str = str(raw_val).strip().strip('"').strip("'")
            file_path = self._resolve_path(path_str)

            if not file_path.exists():
                continue

            needs_convert = False
            if file_path.suffix.lower() in (".tif", ".tiff"):
                needs_convert = True
            elif file_path.suffix == ".zarr":
                ok, _ = self._check_zarr_structure(file_path)
                if not ok:
                    needs_convert = True

            if needs_convert:
                # We pass refresh_ui=False here to prevent the table from being rebuilt 
                # inside the loop, which causes the crash!
                if self._convert_single(path_str, ct, sample_name, row_idx, 
                                        save_metadata=False, refresh_ui=False):
                    converted_any = True

        if converted_any:
            if save_metadata:
                self._save_metadata(refresh_ui=refresh_ui)
        elif refresh_ui:
            self.log(f"No conversions needed for {sample_name}")
        
        return converted_any

    def _convert_all(self):
        """Convert all actionable files across all samples."""
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
        from qtpy.QtWidgets import QMessageBox
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
    """APOC (GPU) segmentation page with embedded training and batch inference."""

    STRATEGIES = [
        "APOC (Direct Instance Segmentation)",
        "APOC + EDT/Watershed",
        "APOC Probability Map + Watershed",
    ]

    def __init__(self, viewer, metadata_loader, log_callback=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self._training_widget = None
        self._is_session_active = False
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

        # ── Training section (embedded APOCTrainingWidget) ──────
        self.training_group = QGroupBox("🎯 APOC Classifier Training")
        self.training_layout = QVBoxLayout(self.training_group)
        self.training_layout.setContentsMargins(4, 4, 4, 4)

        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
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

        self.training_placeholder = QLabel(
            "Load metadata to enable APOC classifier training."
        )
        self.training_placeholder.setWordWrap(True)
        self.training_placeholder.setStyleSheet("color:#888; font-style:italic; padding:10px;")
        self.training_layout.addWidget(self.training_placeholder)
        layout.addWidget(self.training_group)

        # ── Strategy selector (global) ──────────────────────────
        strat_group = QGroupBox("⚙️ Segmentation Strategy")
        strat_lay = QVBoxLayout(strat_group)
        strat_lay.setSpacing(4)

        strat_desc = QLabel(
            "Strategy determines how the APOC classifier output is<br>"
            "converted to instance segmentation labels."
        )
        strat_desc.setWordWrap(True)
        strat_desc.setStyleSheet("color:#888; font-size:10px;")
        strat_lay.addWidget(strat_desc)

        self.combo_strategy = QComboBox()
        self.combo_strategy.addItems(self.STRATEGIES)
        strat_lay.addWidget(self.combo_strategy)

        # Strategy-specific parameters (stacked)
        self.strategy_stack = QStackedWidget()

        # Page 0: Direct — only size filter
        direct_page = QWidget()
        direct_lay = QVBoxLayout(direct_page)
        direct_lay.setContentsMargins(0, 0, 0, 0)
        
        # Per-cell-type spinboxes — rebuilt in _rebuild_size_filter_form()
        self._size_filter_spins: dict = {}   # ct → QSpinBox
        self._sf_form_widget = QWidget()
        self._sf_form_layout = QFormLayout(self._sf_form_widget)
        self._sf_form_layout.setContentsMargins(0, 0, 0, 0)
        self._sf_form_layout.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
        self._sf_placeholder = QLabel("Load metadata to configure per-cell-type thresholds.")
        self._sf_placeholder.setStyleSheet("color:#888; font-style:italic; font-size:10px;")
        self._sf_form_layout.addRow(self._sf_placeholder)
        
        sf_desc = QLabel(
            "Remove segments smaller than a minimum voxel count.\n"
            "Applied per-timepoint, per cell type."
        )
        sf_desc.setStyleSheet("color:#888; font-size:10px;")
        direct_lay.addWidget(sf_desc)
        direct_lay.addWidget(self._sf_form_widget)
        direct_lay.addStretch()
        self.strategy_stack.addWidget(direct_page)

        # Page 1: EDT/Watershed
        edt_page = QWidget()
        edt_form = QFormLayout(edt_page)
        edt_form.setContentsMargins(0, 0, 0, 0)
        edt_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.spin_edt_threshold = QDoubleSpinBox()
        self.spin_edt_threshold.setRange(0.1, 100.0)
        self.spin_edt_threshold.setSingleStep(0.5)
        self.spin_edt_threshold.setValue(2.5)
        self.spin_edt_threshold.setMaximumWidth(90)
        edt_form.addRow("EDT threshold:", make_help_row(
            self.spin_edt_threshold,
            "EDT Threshold",
            "Euclidean Distance Transform threshold for watershed\n"
            "seed detection. Higher = fewer, larger segments."
        ))

        self.spin_edt_min_size = QSpinBox()
        self.spin_edt_min_size.setRange(1, 999999)
        self.spin_edt_min_size.setValue(100)
        self.spin_edt_min_size.setMaximumWidth(90)
        edt_form.addRow("Min segment size:", make_help_row(
            self.spin_edt_min_size,
            "Minimum Segment Size",
            "Segments with fewer pixels than this are removed."
        ))
        self.strategy_stack.addWidget(edt_page)

        # Page 2: Probability Map + Watershed
        prob_page = QWidget()
        prob_form = QFormLayout(prob_page)
        prob_form.setContentsMargins(0, 0, 0, 0)
        prob_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.spin_prob_threshold = QDoubleSpinBox()
        self.spin_prob_threshold.setRange(0.0, 1.0)
        self.spin_prob_threshold.setSingleStep(0.05)
        self.spin_prob_threshold.setDecimals(2)
        self.spin_prob_threshold.setValue(0.50)
        self.spin_prob_threshold.setMaximumWidth(90)
        prob_form.addRow("Probability threshold:", make_help_row(
            self.spin_prob_threshold,
            "Probability Threshold",
            "Minimum probability for a pixel to be considered foreground.\n"
            "Range 0-1. Higher = stricter, fewer false positives."
        ))

        self.spin_prob_min_size = QSpinBox()
        self.spin_prob_min_size.setRange(1, 999999)
        self.spin_prob_min_size.setValue(100)
        self.spin_prob_min_size.setMaximumWidth(90)
        prob_form.addRow("Min segment size:", make_help_row(
            self.spin_prob_min_size,
            "Minimum Segment Size",
            "Segments with fewer pixels than this are removed."
        ))
        self.strategy_stack.addWidget(prob_page)

        self.combo_strategy.currentIndexChanged.connect(self.strategy_stack.setCurrentIndex)
        strat_lay.addWidget(self.strategy_stack)
        layout.addWidget(strat_group)

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
        self.btn_run_segmentation.clicked.connect(self._on_run_segmentation)

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

    # ── Load Training Data ──────────────────────────────────────
    def _on_load_training_clicked(self, interactive=True):
        try:
            if self.metadata_loader.metadata is None:
                self.log("⚠️ Cannot generate training data: No metadata loaded.")
                return

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
            )
            organoid_types = detect_organoid_types_from_metadata(md)
            immune_types = detect_immune_cell_types_from_metadata(md)
            other_types = detect_other_cell_types_from_metadata(md)
            
            # Use apoc_train helper to fetch and process images
            from behav3d.preprocessing.segmentation.apoc_train import _load_training_images
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

            # Clear viewer layers related to this
            layers_to_remove = [l for l in self.viewer.layers if 'Channel' in l.name or 'User Provided Labels' in l.name]
            for l in layers_to_remove:
                self.viewer.layers.remove(l)

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

            # Add Image layers
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
                img_layer.contrast_limits_range = (0, float(channel_data.max()))

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

            self.log("✅ Training data generated/loaded in viewer!")
            self._is_session_active = True

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
        """Remove all APOC training layers from the viewer and reset session state.
        Mirrors PixelClassifierWidget.reset_ui() — called when switching tabs or methods."""
        layers_to_remove = [
            l for l in list(self.viewer.layers)
            if 'Channel' in l.name or 'User Provided Labels' in l.name
        ]
        for l in layers_to_remove:
            try:
                self.viewer.layers.remove(l)
            except Exception:
                pass
        self._is_session_active = False
        self.log("APOC training session cleared.")

    # ── Per-CT size filter form ─────────────────────────────────
    def _rebuild_size_filter_form(self, all_types: list):
        """Rebuild per-cell-type min-size spinboxes for the size filter group."""
        # Clear old widgets
        while self._sf_form_layout.rowCount():
            self._sf_form_layout.removeRow(0)
        self._size_filter_spins.clear()

        if not all_types:
            self._sf_form_layout.addRow(self._sf_placeholder)
            return

        for ct in all_types:
            spin = QSpinBox()
            spin.setRange(1, 9_999_999)
            spin.setValue(500)
            spin.setMaximumWidth(100)
            spin.setSuffix(" vx")
            self._sf_form_layout.addRow(f"{ct}:", spin)
            self._size_filter_spins[ct] = spin

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
        )

        org = detect_organoid_types_from_metadata(md)
        imm = detect_immune_cell_types_from_metadata(md)
        oth = detect_other_cell_types_from_metadata(md)
        all_types = org + imm + oth
        has_death = has_dead_channel(md)

        if not all_types:
            return

        # Rebuild per-CT size filter spinboxes
        self._rebuild_size_filter_form(all_types)

        # Remove old training widget if present
        if self._training_widget is not None:
            self.training_layout.removeWidget(self._training_widget)
            self._training_widget.deleteLater()
            self._training_widget = None
        if self.training_placeholder.parent() is not None:
            self.training_layout.removeWidget(self.training_placeholder)
            self.training_placeholder.hide()

        # Create and embed APOCTrainingWidget
        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        from behav3d.preprocessing.segmentation.apoc_train import APOCTrainingWidget
        params = self.metadata_loader.behav3d_parameters
        apoc_params = params.get("apoc", {})

        self._training_widget = APOCTrainingWidget(
            viewer=self.viewer,
            pixel_class_outdir=str(pixel_class_outdir),
            all_cell_types=all_types,
            has_death=has_death,
            initial_params=apoc_params,
        )
        self.training_layout.addWidget(self._training_widget)

    # ── Queue param snapshot ────────────────────────────────────
    def get_queue_params(self) -> dict:
        """Snapshot current widget state for use by the processing queue."""
        idx = self.combo_strategy.currentIndex()
        params = {
            "strategy_index": idx,
            "strategy_name": self.STRATEGIES[idx],
            "overwrite": self.check_overwrite.isChecked(),
            "workers": self.spin_workers.value(),
            "process_all": self.check_process_all.isChecked(),
            "t_start": self.spin_t_start.value(),
            "t_end": self.spin_t_end.value(),
        }
        if idx == 0:  # Direct APOC
            params["min_sizes"] = {ct: spin.value() for ct, spin in self._size_filter_spins.items()}
        elif idx == 1:  # EDT/Watershed
            params["edt_threshold"] = self.spin_edt_threshold.value()
            params["edt_min_size"] = self.spin_edt_min_size.value()
        elif idx == 2:  # Probability Map + Watershed
            params["prob_threshold"] = self.spin_prob_threshold.value()
            params["prob_min_size"] = self.spin_prob_min_size.value()
        return params



    # ── Run batch segmentation ──────────────────────────────────
    def _on_run_segmentation(self):
        try:
            md = self.metadata_loader.metadata
            if md is None:
                self.log("⚠️ No metadata loaded.")
                return

            from behav3d.preprocessing.segmentation.apoc_segment import run_apoc_segmentation

            output_dir = Path(self.metadata_loader.output_dir)
            strategy = self.STRATEGIES[self.combo_strategy.currentIndex()]

            # Timepoint range
            if self.check_process_all.isChecked():
                timepoint_range = None
            else:
                t_start = self.spin_t_start.value()
                t_end = self.spin_t_end.value()
                if t_start > t_end:
                    self.log("Error: Start timepoint must be <= End timepoint.")
                    return
                timepoint_range = (t_start, t_end)

            # Collect APOC config from training widget if available
            apoc_config = {}
            if self._training_widget is not None:
                for ct, tab in self._training_widget.tabs.items():
                    if hasattr(tab, 'get_params'):
                        p = tab.get_params()
                        for k, v in p.items():
                            apoc_config[f"apoc_{ct}_{k}"] = v

            # Embed the per CT size filters into apoc_config so direct APOC sees them
            for ct, spin in self._size_filter_spins.items():
                apoc_config[f"{ct}_segment_size_min"] = spin.value()

            # Classifier paths: scan PixelClassification dir for .cl files
            pixel_class_outdir = output_dir / "images" / "PixelClassification"

            self.log(f"Starting APOC batch segmentation ({strategy})...")

            metadata_csv = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv", "")

            updated_metadata = run_apoc_segmentation(
                output_dir=str(output_dir),
                metadata=md,
                metadata_csv_path=metadata_csv,
                timepoint_range=timepoint_range,
                apoc_config=apoc_config,
                overwrite_existing=self.check_overwrite.isChecked(),
                n_workers=self.spin_workers.value(),
                apoc_strategy=strategy,
            )

            # Apply size filtering immediately if strategy is Direct APOC
            if self.combo_strategy.currentIndex() == 0:
                self.log("Applying size filtering for Direct APOC strategy...")
                try:
                    from behav3d.preprocessing.segmentation.size_filter import filter_segments_by_size
                    for ct, spin in self._size_filter_spins.items():
                        min_size = spin.value()
                        for sample_name in md['sample_name'].unique():
                            seg_path = output_dir / "images" / sample_name / f"{sample_name}_{ct}_segments.zarr"
                            if seg_path.exists():
                                self.log(f"  [{ct}] {sample_name} — min={min_size} vx")
                                filter_segments_by_size(
                                    segments_zarr_path=str(seg_path),
                                    min_size_voxels=min_size,
                                )
                except Exception as ex:
                    self.log(f"Warning: size filtering failed: {ex}")

            # Update metadata
            if updated_metadata is not None:
                self.metadata_loader.metadata = updated_metadata
                if metadata_csv:
                    try:
                        updated_metadata.to_csv(metadata_csv, index=False)
                    except Exception as e:
                        self.log(f"Warning: Could not save metadata: {e}")

            self.log("✅ APOC batch segmentation finished!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"❌ Error during APOC segmentation: {e}")


