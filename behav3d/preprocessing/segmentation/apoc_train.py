"""
Custom APOC (Accelerated Pixel and Object Classification) widget for BEHAV3D.

This file is self-contained: it handles image loading, the Qt training widget
(with per-cell-type tabs), and APOC train/predict calls.
The original pixel classifier in napari_pixelclassifier.py is left completely untouched.
"""

import os
import gc
import json
import time
import shutil
from pathlib import Path

import numpy as np
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox,
    QSpinBox, QLineEdit, QPushButton as QtPushButton,
    QCheckBox, QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
    QSizePolicy, QTabWidget, QFrame, QMessageBox,
    QDialog, QTableWidget, QTableWidgetItem, QHeaderView,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

import dask.array as da
import zarr
import pyclesperanto_prototype as cle

from behav3d.io.images import load_image, load_zarr, save_as_zarr, append_to_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape

# ---------------------------------------------------------------------------
# Feature grid constants (matching the official napari-apoc widget)
# ---------------------------------------------------------------------------

# Sigma columns (these match the default column headers in the official widget)
APOC_SIGMAS = [0.3, 0.5, 1, 2, 3, 4, 5, 10, 15, 25]

# Feature rows shown in presets
APOC_FEATURES = [
    ("gaussian_blur",               "Gauss"),
    ("difference_of_gaussian",      "DoG"),
    ("laplace_box_of_gaussian_blur","LoG"),
    ("sobel_of_gaussian_blur",      "SoG"),
]
# For custom preset every feature is available (same list here; extend if needed)
APOC_ALL_FEATURES = APOC_FEATURES

# Format sigma values nicely (drop trailing .0)
def _fmt_sigma(s):
    return str(int(s)) if float(s) == int(s) else str(s)


# Each preset is a set of (feature_key, sigma_str) pairs that should be
# pre-checked.  The sigma value is stored as a string matching _fmt_sigma().
FEATURE_PRESETS = {
    "small_quick": {
        "label": "Small / Quick",
        "checked": {
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
        },
    },
    "medium_quick": {
        "label": "Medium / Quick",
        "checked": {
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("difference_of_gaussian",      "5"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
            ("laplace_box_of_gaussian_blur","5"),
            ("sobel_of_gaussian_blur",      "1"),
            ("sobel_of_gaussian_blur",      "2"),
            ("sobel_of_gaussian_blur",      "5"),
        },
    },
    "large_quick": {
        "label": "Large / Quick",
        "checked": {
            ("gaussian_blur",               "1"),
            ("gaussian_blur",               "2"),
            ("gaussian_blur",               "5"),
            ("gaussian_blur",               "10"),
            ("gaussian_blur",               "25"),
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("difference_of_gaussian",      "5"),
            ("difference_of_gaussian",      "10"),
            ("difference_of_gaussian",      "25"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
            ("laplace_box_of_gaussian_blur","5"),
            ("laplace_box_of_gaussian_blur","10"),
            ("laplace_box_of_gaussian_blur","25"),
            ("sobel_of_gaussian_blur",      "1"),
            ("sobel_of_gaussian_blur",      "2"),
            ("sobel_of_gaussian_blur",      "5"),
            ("sobel_of_gaussian_blur",      "10"),
            ("sobel_of_gaussian_blur",      "25"),
        },
    },
    "custom": {
        "label": "Custom",
        "checked": set(),  # user starts with all unchecked
    },
}


def _checked_set_for_preset(preset_name):
    """Return the set of (feature_key, sigma_str) pairs that should be checked for a preset."""
    return set(FEATURE_PRESETS.get(preset_name, FEATURE_PRESETS["medium_quick"])["checked"])


def _build_feature_string_from_checked(checked_set, consider_original=False, current_sigmas=None):
    """Build an APOC feature string from a set of (feature_key, sigma_str) pairs."""
    parts = []
    if consider_original:
        parts.append("original")
    if current_sigmas is None:
        current_sigmas = APOC_SIGMAS
    sigma_strs = [_fmt_sigma(s) for s in current_sigmas]
    # Keep deterministic order: rows in APOC_ALL_FEATURES order, sigmas in current_sigmas order
    for feat_key, _label in APOC_ALL_FEATURES:
        for s_str in sigma_strs:
            if (feat_key, s_str) in checked_set:
                parts.append(f"{feat_key}={s_str}")
    return " ".join(parts)


# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Image loading (same logic as CPU pipeline for consistency)
# ---------------------------------------------------------------------------

def _load_training_images(
    metadata, output_dir, examples_per_sample, organoid_types, immune_types, other_types, overwrite_images=False
):
    """
    Load training images from metadata.
    Returns a list of (C, Z, Y, X) arrays — one per selected timepoint — plus
    the pixel_class_outdir path, has_death flag, and all cell types list.
    """
    from behav3d.core.metadata import has_dead_channel

    has_death = has_dead_channel(metadata)
    all_cell_types = organoid_types + immune_types + other_types

    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)

    image_outpath = Path(pixel_class_outdir, "PixelClassifier_Images.zarr")

    # --- Check cached dataset ---
    if image_outpath.exists():
        if overwrite_images:
            print("Overwrite sample images requested — deleting cached data...")
            shutil.rmtree(image_outpath)
        else:
            print("Loading cached training images from previous session...")
            # We assume the cache matches the metadata order of the experimental setup
            # Find the channel axis based on metadata of the first sample
            first_sample = metadata.iloc[0]
            axis_order = first_sample.get("dimension_order", "TCZYX")
            if not isinstance(axis_order, str) or not axis_order:
                axis_order = "TCZYX"
            
            cached = load_zarr(image_outpath)
            # The cached data is (C, T, Z, Y, X)
            # We need to return a list of (C, Z, Y, X) arrays, one per timepoint
            all_images = [np.asarray(cached[:, t, :, :, :]) for t in range(cached.shape[1])]
            return all_images, pixel_class_outdir, has_death, all_cell_types

    # --- Load fresh from raw files ---
    all_images = []
    max_shape = None
    n_samples = 0

    for _, sample in metadata.iterrows():
        sample_name = sample.get("sample_name", "unknown")
        raw_image_path = sample.get("raw_image_path", "")

        if not raw_image_path or not Path(raw_image_path).exists():
            print(f"⚠️ Skipping {sample_name}: raw image not found")
            continue

        n_samples += 1
        img_dir = Path(output_dir, "images", sample_name)
        img_dir.mkdir(parents=True, exist_ok=True)

        # Get dimension order from metadata if available
        axis_order = sample.get("dimension_order", "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"

        # Use raw_image_path directly as it is the source of truth in metadata
        # (and is updated to the .zarr path if conversion was run)
        # load_image already normalizes to TCZYX (default_axis_order in images.py)
        img = load_image(raw_image_path, axis_order=axis_order)

        # Since load_image returned TCZYX, T is always axis 0
        t_axis = 0
            
        n_timepoints = img.shape[t_axis]
        n_to_select = min(examples_per_sample, n_timepoints)

        # Select equally spaced timepoints (first, last, and middle)
        if n_to_select <= 1:
            t_indices = [0]
        else:
            t_indices = np.round(np.linspace(0, n_timepoints - 1, n_to_select)).astype(int)
            t_indices = sorted(list(set(t_indices)))
            
        t_indices_list = [int(t) for t in t_indices]
        print(f"  {sample_name}: selected {len(t_indices)} equidistant timepoints {t_indices_list}")

        for t_idx in t_indices:
            # Fetch only the specific timepoint slice into memory
            frame = np.asarray(np.take(img, t_idx, axis=t_axis)) 
            all_images.append(frame)
            frame_shape = frame.shape
            if max_shape is None:
                max_shape = list(frame_shape)
            else:
                for i in range(len(max_shape)):
                    max_shape[i] = max(max_shape[i], frame_shape[i])

    if max_shape is not None:
        all_images = [zeropad_image_to_match_shape(img, max_shape) for img in all_images]

    # Stash on disk for next time
    import dask.array as da
    import gc
    all_images_stack = da.stack(all_images)
    # Transpose to (C, T, Z, Y, X) to match original classifier
    all_images_stack = all_images_stack.transpose(1, 0, 2, 3, 4)
    save_as_zarr(all_images_stack, image_outpath)
    del all_images_stack
    gc.collect()

    # Re-load to ensure we have the right shape
    cached = load_zarr(image_outpath)
    # The cached data is (C, T, Z, Y, X)
    # We need to return a list of (C, Z, Y, X) arrays, one per timepoint
    all_images = [np.asarray(cached[:, t, :, :, :]) for t in range(cached.shape[1])]

    return all_images, pixel_class_outdir, has_death, all_cell_types


# ---------------------------------------------------------------------------
# Per-cell-type tab panel
# ---------------------------------------------------------------------------

class CellTypeTab(QWidget):
    """
    Widget for a single cell type tab.
    Contains:
      - channel selection group
      - preset dropdown (kept)
      - collapsible "Tune Features" QGroupBox with a sigma × feature checkbox grid
      - "Consider original image as well" standalone checkbox
      - "Show classifier statistics" button
      - max_depth / num_ensembles RF parameters
    """

    def __init__(self, cell_type, viewer, initial_params=None, parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.viewer = viewer
        ip = initial_params or {}

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # ── Channel selection ────────────────────────────────────────────────
        chan_group = QGroupBox("Image Channel Inputs")
        chan_layout = QVBoxLayout()
        chan_layout.setSpacing(2)
        self.channel_checkboxes = []
        self.chan_checkbox_container = QWidget()
        self.chan_checkbox_layout = QVBoxLayout()
        self.chan_checkbox_layout.setContentsMargins(0, 0, 0, 0)
        self.chan_checkbox_container.setLayout(self.chan_checkbox_layout)
        chan_layout.addWidget(self.chan_checkbox_container)
        chan_group.setLayout(chan_layout)
        layout.addWidget(chan_group)

        # Labeling hint
        hint = QLabel("Labels: <b>1</b> = background&nbsp;&nbsp; <b>2</b> = foreground")
        hint.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(hint)

        # ── Preset dropdown ──────────────────────────────────────────────────
        feat_row = QHBoxLayout()
        feat_row.addWidget(QLabel("Preset:"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(FEATURE_PRESETS.keys()))
        saved_preset = ip.get(f"apoc_{cell_type}_feature_preset", "medium_quick")
        if saved_preset not in FEATURE_PRESETS:
            saved_preset = "medium_quick"
        self.feature_combo.setCurrentText(saved_preset)
        feat_row.addWidget(self.feature_combo)
        feat_row.addStretch()
        layout.addLayout(feat_row)

        # ── Tune Features collapsible group ─────────────────────────────────
        self.tune_group = QGroupBox("Tune Features")
        self.tune_group.setCheckable(True)
        self.tune_group.setChecked(False)  # collapsed by default
        self.tune_layout = QVBoxLayout()
        self.tune_layout.setContentsMargins(4, 4, 4, 4)
        self.tune_layout.setSpacing(2)

        # ── Sigmas Input Row ──
        sigma_row = QHBoxLayout()
        sigma_row.addWidget(QLabel("Custom Sigmas:"))
        self.sigma_input = QLineEdit(", ".join(_fmt_sigma(s) for s in APOC_SIGMAS))
        self.sigma_input.setToolTip("Comma or space-separated list of sigmas to use for APOC filtering.")
        sigma_row.addWidget(self.sigma_input)
        self.update_grid_btn = QtPushButton("Update Grid")
        self.update_grid_btn.clicked.connect(self._on_update_grid)
        sigma_row.addWidget(self.update_grid_btn)
        self.tune_layout.addLayout(sigma_row)

        self.current_sigmas = list(APOC_SIGMAS)
        self._feat_sigma_checks = {}
        self._grid_widget = None

        self.tune_group.setLayout(self.tune_layout)
        layout.addWidget(self.tune_group)

        # ── Consider original image checkbox ─────────────────────────────────
        self.consider_original_cb = QCheckBox("Consider original image as well")
        saved_orig = bool(ip.get(f"apoc_{cell_type}_consider_original", False))
        self.consider_original_cb.setChecked(saved_orig)
        self.consider_original_cb.stateChanged.connect(self._update_preview)
        layout.addWidget(self.consider_original_cb)

        # ── Feature string preview ───────────────────────────────────────────
        self.preview_label = QLabel("")
        self.preview_label.setWordWrap(True)
        self.preview_label.setStyleSheet(
            "color: #888; font-size: 10px; padding: 2px 4px; "
            "background: rgba(0,0,0,0.05); border-radius: 3px;"
        )
        layout.addWidget(self.preview_label)

        # ── Show classifier statistics button ────────────────────────────────
        self.stats_btn = QtPushButton("Show classifier statistics")
        self.stats_btn.setToolTip(
            "Display the feature importance distribution of the trained classifier."
        )
        self.stats_btn.clicked.connect(self._on_show_statistics)
        layout.addWidget(self.stats_btn)

        # ── RF parameters ────────────────────────────────────────────────────
        rf_row = QHBoxLayout()
        rf_row.addWidget(QLabel("Max depth:"))
        self.max_depth_spin = QSpinBox()
        self.max_depth_spin.setRange(1, 20)
        self.max_depth_spin.setValue(int(ip.get(f"apoc_{cell_type}_max_depth", 2)))
        rf_row.addWidget(self.max_depth_spin)
        rf_row.addWidget(QLabel("Trees:"))
        self.num_ensembles_spin = QSpinBox()
        self.num_ensembles_spin.setRange(10, 1000)
        self.num_ensembles_spin.setSingleStep(10)
        self.num_ensembles_spin.setValue(int(ip.get(f"apoc_{cell_type}_num_ensembles", 100)))
        rf_row.addWidget(self.num_ensembles_spin)
        rf_row.addStretch()
        layout.addLayout(rf_row)

        layout.addStretch()
        self.setLayout(layout)

        # Wire up preset changes LAST so initial_params can be applied first
        # (we do not auto-reset when the user changes the preset manually —
        #  we only reset if the preset is changed *and* the user hasn't
        #  customised the grid yet; simplest UX: always reset on preset change)
        self.feature_combo.currentTextChanged.connect(self._on_preset_changed)
        self.tune_group.toggled.connect(self._on_tune_toggled)

        # Restore saved grid configuration
        saved_sigmas = ip.get(f"apoc_{cell_type}_grid_sigmas")
        if saved_sigmas:
            self.sigma_input.setText(saved_sigmas)
            self._on_update_grid() # This parses sigmas, builds the grid, etc.
        else:
            self._build_grid()

        # Restore saved checkbox state from config (if available)
        saved_checked = ip.get(f"apoc_{cell_type}_checked_features")
        if saved_checked:
            # saved_checked is a list of [feat_key, sigma_str] pairs
            saved_set = {tuple(pair) for pair in saved_checked}
            self._apply_checked_set(saved_set)
        else:
            # Apply preset defaults
            self._apply_preset_defaults(saved_preset)

        self._update_preview()

    # ────────────────────────────────────────────────────────────────────────
    # Internal helpers
    # ────────────────────────────────────────────────────────────────────────

    def _apply_checked_set(self, checked_set):
        """Check/uncheck grid cells according to the given set of (feat_key, sigma_str)."""
        for key, cb in self._feat_sigma_checks.items():
            cb.blockSignals(True)
            cb.setChecked(key in checked_set)
            cb.blockSignals(False)

    def _apply_preset_defaults(self, preset_name):
        """Reset the grid to the preset's default checked cells."""
        self._apply_checked_set(_checked_set_for_preset(preset_name))

    def _get_current_checked_set(self):
        """Return the set of (feat_key, sigma_str) currently checked in the grid."""
        return {key for key, cb in self._feat_sigma_checks.items() if cb.isChecked()}

    def _on_preset_changed(self, preset_name):
        """Reset grid checkboxes to the new preset's defaults."""
        if preset_name != "custom":
            self.sigma_input.setText(", ".join(_fmt_sigma(s) for s in APOC_SIGMAS))
            self.current_sigmas = list(APOC_SIGMAS)
            self._build_grid()
            self._apply_preset_defaults(preset_name)
        self._update_preview()

    def _on_tune_toggled(self, checked):
        """Show/hide the inner grid widget when the group box is toggled."""
        if self._grid_widget:
            self._grid_widget.setVisible(checked)

    def _on_update_grid(self):
        """Parse custom sigmas from input and rebuild the feature grid."""
        text = self.sigma_input.text()
        try:
            parts = text.replace(",", " ").split()
            new_sigmas = [float(p) for p in parts if p.strip()]
            if not new_sigmas:
                raise ValueError("Empty list")
            self.current_sigmas = new_sigmas
            self.feature_combo.setCurrentText("custom")
            self._build_grid()
        except ValueError:
            QMessageBox.warning(self, "Invalid Sigmas", "Could not parse sigma values. Please enter numbers separated by spaces or commas.")

    def _build_grid(self):
        """Rebuild the checkbox grid based on current_sigmas."""
        old_checked = set()
        if self._feat_sigma_checks:
            old_checked = self._get_current_checked_set()

        if self._grid_widget is not None:
            self.tune_layout.removeWidget(self._grid_widget)
            self._grid_widget.deleteLater()

        self._grid_widget = QWidget()
        grid = QGridLayout()
        grid.setSpacing(2)
        grid.setContentsMargins(0, 0, 0, 0)
        self._feat_sigma_checks = {}

        sigma_header = QLabel("sigma")
        sigma_header.setStyleSheet("font-weight: bold; font-size: 11px;")
        grid.addWidget(sigma_header, 0, 0)

        for col_idx, s in enumerate(self.current_sigmas):
            s_str = _fmt_sigma(s)
            lbl = QLabel(s_str)
            lbl.setStyleSheet("font-size: 10px;")
            lbl.setAlignment(Qt.AlignCenter)
            grid.addWidget(lbl, 0, col_idx + 1)

        for row_idx, (feat_key, feat_label) in enumerate(APOC_ALL_FEATURES):
            feat_lbl = QLabel(feat_label)
            feat_lbl.setStyleSheet("font-size: 11px;")
            grid.addWidget(feat_lbl, row_idx + 1, 0)
            for col_idx, s in enumerate(self.current_sigmas):
                s_str = _fmt_sigma(s)
                cb = QCheckBox()
                cb.setChecked(False)
                cb.setStyleSheet("QCheckBox { margin: 0px; padding: 0px; }")
                cb.stateChanged.connect(self._update_preview)
                grid.addWidget(cb, row_idx + 1, col_idx + 1, alignment=Qt.AlignCenter)
                self._feat_sigma_checks[(feat_key, s_str)] = cb

        self._grid_widget.setLayout(grid)
        self.tune_layout.insertWidget(1, self._grid_widget)

        # Retain checked state for checkboxes that still exist in the new grid
        self._apply_checked_set(old_checked)
        self._update_preview()

    def _update_preview(self):
        """Rebuild the feature string preview from the current grid state."""
        feat_str = self.get_feature_string()
        # Show a truncated preview if very long
        if len(feat_str) > 120:
            display = feat_str[:117] + "…"
        else:
            display = feat_str
        self.preview_label.setText(f"<b>Features:</b> {display}")
        self.preview_label.setToolTip(feat_str)

    def _get_clf_path(self):
        """Return the expected .cl file path for this cell type."""
        if not hasattr(self, '_pixel_class_outdir') or not self._pixel_class_outdir:
            return None
        ct = self.cell_type
        if ct == "dead":
            fname = "PixelClassifier_Death.cl"
        else:
            fname = f"PixelClassifier_{ct.capitalize()}.cl"
        return Path(self._pixel_class_outdir) / fname

    def _on_show_statistics(self):
        """Show a window with the classifier feature statistics."""
        clf_path = self._get_clf_path()
        if clf_path is None or not Path(clf_path).exists():
            QMessageBox.information(
                self,
                "Classifier Statistics",
                f"No trained classifier found for '{self.cell_type}'.\n"
                "Train the classifier first, then click this button."
            )
            return
        try:
            import apoc
            clf = apoc.ObjectSegmenter(opencl_filename=str(clf_path))
            feature_importances = clf.feature_importances()  # dict feature -> importance

            dlg = QDialog(self)
            dlg.setWindowTitle(f"Classifier Statistics — {self.cell_type.capitalize()}")
            dlg_layout = QVBoxLayout(dlg)

            table = QTableWidget(len(feature_importances), 2)
            table.setHorizontalHeaderLabels(["Feature", "Importance"])
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
            for i, (feat, imp) in enumerate(sorted(
                feature_importances.items(), key=lambda x: -x[1]
            )):
                table.setItem(i, 0, QTableWidgetItem(str(feat)))
                item = QTableWidgetItem(f"{imp:.4f}")
                item.setTextAlignment(Qt.AlignCenter)
                # Colour code: green for high importance
                r = int(max(0, 255 - imp * 1200))
                g = int(min(255, 80 + imp * 1200))
                item.setBackground(QColor(r, g, 80))
                table.setItem(i, 1, item)

            dlg_layout.addWidget(table)
            dlg.setMinimumSize(500, 400)
            dlg.exec_()
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Classifier Statistics Error",
                f"Could not load statistics for '{self.cell_type}':\n{exc}"
            )

    # ────────────────────────────────────────────────────────────────────────
    # Public API
    # ────────────────────────────────────────────────────────────────────────

    def get_feature_string(self):
        """Return the full APOC feature string based on current grid state."""
        checked = self._get_current_checked_set()
        consider_orig = self.consider_original_cb.isChecked()
        return _build_feature_string_from_checked(checked, consider_original=consider_orig, current_sigmas=self.current_sigmas)

    def refresh_channel_checkboxes(self):
        """Rebuild channel checkboxes from current Napari image layers."""
        for cb in self.channel_checkboxes:
            self.chan_checkbox_layout.removeWidget(cb)
            cb.deleteLater()
        self.channel_checkboxes = []

        for layer in self.viewer.layers:
            if (
                isinstance(layer, napari.layers.Image)
                and not layer.name.startswith("Pixel Classification")
            ):
                cb = QCheckBox(layer.name)
                cb.setChecked(True)
                self.chan_checkbox_layout.addWidget(cb)
                self.channel_checkboxes.append(cb)

    def get_config(self):
        """Return a dict with all current widget values."""
        checked_set = self._get_current_checked_set()
        return {
            "feature_preset":     self.feature_combo.currentText(),
            "grid_sigmas":        self.sigma_input.text().strip(),
            # Legacy keys kept for backward compat:
            "sigmas":             ",".join(
                sorted({s for _, s in checked_set},
                       key=lambda x: APOC_SIGMAS.index(float(x)) if float(x) in APOC_SIGMAS else 999)
            ),
            "custom_feature_string": "",
            "feature_string":     self.get_feature_string(),
            "consider_original":  self.consider_original_cb.isChecked(),
            # Rich state for round-tripping
            "checked_features":   [list(pair) for pair in checked_set],
            "max_depth":          self.max_depth_spin.value(),
            "num_ensembles":      self.num_ensembles_spin.value(),
            "channels":           [
                cb.text() for cb in self.channel_checkboxes if cb.isChecked()
            ],
        }

    def apply_config(self, cfg):
        """Restore widget values from a config dict (used by 'Apply to all tabs')."""
        if "feature_preset" in cfg:
            # Block the auto-reset signal while we apply config
            self.feature_combo.blockSignals(True)
            self.feature_combo.setCurrentText(cfg["feature_preset"])
            self.feature_combo.blockSignals(False)
        if "grid_sigmas" in cfg:
            self.sigma_input.setText(cfg["grid_sigmas"])
            self._on_update_grid()
        if "checked_features" in cfg and cfg["checked_features"]:
            saved_set = {tuple(pair) for pair in cfg["checked_features"]}
            self._apply_checked_set(saved_set)
        elif "feature_preset" in cfg:
            self._apply_preset_defaults(cfg["feature_preset"])
        if "consider_original" in cfg:
            self.consider_original_cb.setChecked(bool(cfg["consider_original"]))
        if "max_depth" in cfg:
            self.max_depth_spin.setValue(int(cfg["max_depth"]))
        if "num_ensembles" in cfg:
            self.num_ensembles_spin.setValue(int(cfg["num_ensembles"]))
        if "channels" in cfg:
            for cb in self.channel_checkboxes:
                cb.setChecked(cb.text() in cfg["channels"])
        self._update_preview()


# ---------------------------------------------------------------------------
# Main APOC training widget
# ---------------------------------------------------------------------------

class APOCTrainingWidget(QWidget):
    """
    Tabbed APOC training widget docked inside the Napari viewer.
    One tab per cell type (+ Death if applicable), each with:
      - channel selection
      - feature preset + sigma values + auto-preview
      - max_depth / num_ensembles
    Global controls: Apply-to-All, Train Current, Train-All, Save Labels.
    """

    def __init__(
        self,
        viewer,
        pixel_class_outdir,
        all_cell_types,
        has_death,
        initial_params=None,
        on_params_changed=None,
        parent=None,
    ):
        super().__init__(parent)
        self.viewer = viewer
        self.pixel_class_outdir = pixel_class_outdir
        self.all_cell_types = all_cell_types
        self.has_death = has_death
        self._initial_params = initial_params or {}
        self._on_params_changed = on_params_changed

        # All tab labels = cell types + optional Death
        self._tab_cell_types = list(all_cell_types)
        if has_death:
            self._tab_cell_types.append("dead")

        self._build_ui()
        self._connect_signals()
        self._refresh_all_channels()

        # Listen for layer changes
        self.viewer.layers.events.inserted.connect(lambda _: self._refresh_all_channels())
        self.viewer.layers.events.removed.connect(lambda _: self._refresh_all_channels())

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        title = QLabel("<h3>APOC Pixel Classification</h3>")
        layout.addWidget(title)

        # --- Tab widget ---
        self.tab_widget = QTabWidget()
        self.tabs = {}  # cell_type -> CellTypeTab

        for ct in self._tab_cell_types:
            tab = CellTypeTab(ct, self.viewer, initial_params=self._initial_params)
            tab._pixel_class_outdir = self.pixel_class_outdir  # for statistics button
            self.tabs[ct] = tab
            self.tab_widget.addTab(tab, ct.capitalize())

        layout.addWidget(self.tab_widget)

        # --- Separator ---
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line)

        # --- Global buttons ---
        global_row1 = QHBoxLayout()

        self.apply_all_btn = QtPushButton("⬇ Apply config to all tabs")
        self.apply_all_btn.setStyleSheet("padding: 5px 10px;")
        self.apply_all_btn.setToolTip("Copy the current tab's preset, sigmas, depth and trees to ALL other tabs")
        global_row1.addWidget(self.apply_all_btn)
        layout.addLayout(global_row1)

        global_row2 = QHBoxLayout()

        self.train_current_btn = QtPushButton("▶ Train current tab")
        self.train_current_btn.setStyleSheet(
            "background-color: #1976D2; color: white; font-weight: bold; padding: 6px 12px;"
        )
        global_row2.addWidget(self.train_current_btn)

        self.train_all_btn = QtPushButton("▶▶ Train ALL classifiers")
        self.train_all_btn.setStyleSheet(
            "background-color: #2e7d32; color: white; font-weight: bold; padding: 6px 12px;"
        )
        global_row2.addWidget(self.train_all_btn)
        layout.addLayout(global_row2)

        # --- Status label ---
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        layout.addStretch()
        self.setLayout(layout)

    def _connect_signals(self):
        self.apply_all_btn.clicked.connect(self._on_apply_to_all)
        self.train_current_btn.clicked.connect(self._on_train_current)
        self.train_all_btn.clicked.connect(self._on_train_all)

    def _refresh_all_channels(self):
        """Refresh channel checkboxes in all tabs."""
        for tab in self.tabs.values():
            tab.refresh_channel_checkboxes()

    # ------------------------------------------------------------------
    # Button handlers
    # ------------------------------------------------------------------

    def _on_apply_to_all(self):
        """Copy current tab's config to all other tabs."""
        current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
        cfg = self.tabs[current_ct].get_config()
        for ct, tab in self.tabs.items():
            if ct != current_ct:
                tab.apply_config(cfg)
        self.status_label.setText(f"↪ Config applied from '{current_ct}' to all other tabs.")

    def _on_train_current(self):
        """Train only the classifier for the currently visible tab."""
        current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
        self._run_training([current_ct])

    def _on_train_all(self):
        """Train classifiers for ALL cell types using each tab's individual config."""
        self._run_training(self._tab_cell_types)

    # ------------------------------------------------------------------
    # Core training logic
    # ------------------------------------------------------------------

    def _get_images_for_tab(self, ct):
        """Return numpy arrays for the image layers checked in a tab."""
        tab = self.tabs[ct]
        images = []
        for cb in tab.channel_checkboxes:
            if cb.isChecked():
                try:
                    layer = self.viewer.layers[cb.text()]
                    # Keep data lazy (dask/zarr) to save memory
                    images.append(layer.data)
                except KeyError:
                    pass
        return images

    def save_user_labels(self, log=print):
        """Save all user-provided labels for all cell types and Dead (if present)."""
        import shutil
        for cell_type in self.all_cell_types:
            lname = f"User Provided Labels ({cell_type.capitalize()})"
            if lname not in [l.name for l in self.viewer.layers]:
                continue
            label_layer = self.viewer.layers[lname]
            outpath = Path(self.pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
            if outpath.exists():
                shutil.rmtree(outpath)
            save_as_zarr(label_layer.data, outpath)
            log(f"Saved {cell_type} labels → {outpath}")

        if self.has_death and "User Provided Labels (Dead)" in [l.name for l in self.viewer.layers]:
            dead_layer = self.viewer.layers["User Provided Labels (Dead)"]
            dead_outpath = Path(self.pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
            if dead_outpath.exists():
                shutil.rmtree(dead_outpath)
            save_as_zarr(dead_layer.data, dead_outpath)
            log(f"Saved Death labels → {dead_outpath}")

        log("✅ All user labels saved!")

    def _run_training(self, cell_types_to_train):
        """Train (and apply) APOC classifiers for the given list of cell types."""
        import apoc

        # Auto-save labels before training
        print("Auto-saving user labels before training...")
        self.save_user_labels(log=print)

        successes = []
        try:
            for ct in cell_types_to_train:
                tab = self.tabs[ct]
                cfg = tab.get_config()
                feature_string = cfg["feature_string"]
                max_depth = cfg["max_depth"]
                num_ensembles = cfg["num_ensembles"]

                self.status_label.setText(f"Processing {ct}...")
                QApplication.processEvents()

                images = self._get_images_for_tab(ct)
                if not images:
                    self.status_label.setText(f"⚠️ No image layers selected for '{ct}'!")
                    continue

                # Get annotation layer (keep lazy if possible)
                if ct == "dead":
                    layer_name = "User Provided Labels (Dead)"
                else:
                    layer_name = f"User Provided Labels ({ct.capitalize()})"

                try:
                    annotation = self.viewer.layers[layer_name].data
                except KeyError:
                    print(f"Skipping '{ct}': annotation layer not found")
                    continue

                if not np.any(annotation):
                    print(f"Skipping '{ct}': no labels drawn")
                    continue

                if ct == "dead":
                    clf_name = "PixelClassifier_Death.cl"
                else:
                    clf_name = f"PixelClassifier_{ct.capitalize()}.cl"

                clf_path = str(Path(self.pixel_class_outdir, clf_name))

                # Erase existing classifier
                if Path(clf_path).exists():
                    apoc.erase_classifier(clf_path)

                # Train with ObjectSegmenter
                clf = apoc.ObjectSegmenter(
                    opencl_filename=clf_path,
                    max_depth=max_depth,
                    num_ensembles=num_ensembles,
                    positive_class_identifier=2,
                )

                has_trained = False
                n_timepoints = annotation.shape[0] if annotation.ndim == 4 else 1

                if annotation.ndim == 4:
                    for t in range(n_timepoints):
                        # Load only the current timepoint slice into memory
                        ann_t = np.asarray(annotation[t])
                        if not np.any(ann_t):
                            continue
                        imgs_t = [np.asarray(img[t]) for img in images]
                        imgs_to_pass = imgs_t[0] if len(imgs_t) == 1 else imgs_t
                        clf.train(feature_string, ann_t, imgs_to_pass, continue_training=has_trained)
                        has_trained = True
                else:
                    # Single timepoint
                    ann_np = np.asarray(annotation)
                    imgs_np = [np.asarray(img) for img in images]
                    imgs_to_pass = imgs_np[0] if len(imgs_np) == 1 else imgs_np
                    clf.train(feature_string, ann_np, imgs_to_pass)
                    has_trained = True

                if not has_trained:
                    continue

                # Predict (visual confirmation)
                if images[0].ndim == 4:
                    results = []
                    for t in range(n_timepoints):
                        # Load only current timepoint for prediction
                        imgs_t = [np.asarray(img[t]) for img in images]
                        imgs_to_pass = imgs_t[0] if len(imgs_t) == 1 else imgs_t
                        res_t = clf.predict(image=imgs_to_pass)
                        results.append(np.asarray(res_t).astype(np.int16))
                    result = np.stack(results, axis=0)
                else:
                    imgs = images[0] if len(images) == 1 else images
                    result = np.asarray(clf.predict(image=imgs)).astype(np.int16)

                # Show in viewer
                seg_layer_name = (
                    "Pixel Classification (Dead)" if ct == "dead"
                    else f"{ct.capitalize()} Segments"
                )
                if seg_layer_name in [l.name for l in self.viewer.layers]:
                    self.viewer.layers[seg_layer_name].data = result
                    self.viewer.layers[seg_layer_name].visible = True
                else:
                    self.viewer.add_labels(result, name=seg_layer_name, opacity=0.8, visible=True)

                # Save to disk
                arr_path = Path(self.pixel_class_outdir, f"{Path(clf_path).stem}_PredictedLabels.zarr")
                if arr_path.exists():
                    shutil.rmtree(arr_path)
                save_as_zarr(result, arr_path)

                successes.append(ct)

            if successes:
                self.status_label.setText(f"✅ Trained: {', '.join(successes)}")
                # Auto-show statistics for each successfully trained classifier
                for ct in successes:
                    if ct in self.tabs:
                        self.tabs[ct]._on_show_statistics()
            else:
                self.status_label.setText("⚠️ No cell types were trained (check labels).")

            # Persist all params back to caller
            self._persist_params()

        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            import traceback
            traceback.print_exc()

    # ------------------------------------------------------------------
    # Config persistence
    # ------------------------------------------------------------------

    def _collect_all_params(self):
        """Return a flat param dict for all tabs, suitable for storing in BEHAV3D config."""
        params = {}
        for ct, tab in self.tabs.items():
            cfg = tab.get_config()
            params[f"apoc_{ct}_feature_preset"]        = cfg["feature_preset"]
            params[f"apoc_{ct}_sigmas"]                = cfg["sigmas"]
            params[f"apoc_{ct}_custom_feature_string"] = cfg["custom_feature_string"]
            params[f"apoc_{ct}_feature_string"]        = cfg["feature_string"]
            params[f"apoc_{ct}_consider_original"]     = cfg["consider_original"]
            params[f"apoc_{ct}_checked_features"]      = cfg["checked_features"]
            params[f"apoc_{ct}_max_depth"]             = cfg["max_depth"]
            params[f"apoc_{ct}_num_ensembles"]         = cfg["num_ensembles"]
            params[f"apoc_{ct}_channels"]              = cfg["channels"]
        return params

    def _persist_params(self):
        """Push current params back to the notebook via the on_params_changed callback."""
        if callable(self._on_params_changed):
            try:
                self._on_params_changed(self._collect_all_params())
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Entry-point called from segmentation.py
# ---------------------------------------------------------------------------

def train_pixel_classifier_apoc(
    output_dir,
    metadata,
    examples_per_sample=3,
    overwrite_images=False,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    initial_params=None,
    on_params_changed=None,
):
    """
    Open a Napari viewer with the tabbed APOC classification widget.
    The original napari_pixelclassifier.py is NOT called.
    """
    ip = initial_params or {}
    gpu_device = ip.get("gpu_device_name")
    if gpu_device:
        print(f"APOC Training: Selecting device {gpu_device}")
        cle.select_device(gpu_device)

    organoid_types = organoid_types or []
    immune_types   = immune_types   or []
    other_types    = other_types    or []
    all_cell_types = organoid_types + immune_types + other_types

    # --- Load images ---
    print("Loading training images for APOC...")
    image_list, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
        metadata=metadata,
        output_dir=output_dir,
        examples_per_sample=examples_per_sample,
        organoid_types=organoid_types,
        immune_types=immune_types,
        other_types=other_types,
        overwrite_images=overwrite_images,
    )

    if not image_list:
        print("⚠️ No training images found!")
        return None

    # --- Pad images to the same shape and stack along T ---
    max_shape = list(image_list[0].shape)
    for img in image_list[1:]:
        for i in range(len(max_shape)):
            max_shape[i] = max(max_shape[i], img.shape[i])
    image_list = [zeropad_image_to_match_shape(img, max_shape) for img in image_list]

    # image_list is a list of timepoints: [t1(C,Z,Y,X), t2(C,Z,Y,X), ...]
    # Stacking along axis 0 gives (T_total, C, Z, Y, X)
    stacked = np.stack(image_list, axis=0)

    # Since stacked was created by stacking (C,Z,Y,X) frames along axis 0,
    # it is ALWAYS (T_total, C, Z, Y, X).
    c_axis_in_stacked = 1

    # --- Create Napari viewer ---
    viewer = napari.Viewer()

    # stacked is (T_total, C, Z, Y, X).
    # Axis 0 = timepoints, axis 1 = channels.
    T_total = stacked.shape[0]
    n_channels = stacked.shape[1]   # ← was incorrectly stacked.shape[0]
    channel_colors = [
        "cyan", "yellow", "red", "green", "magenta", "blue",
        "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
    ]

    for ch in range(n_channels):
        channel_data = stacked[:, ch, :, :, :]  # (T_total, Z, Y, X)
        nonzero = channel_data[channel_data > 0]
        clim = (0, float(np.percentile(nonzero, 99.8))) if nonzero.size > 0 else (0, 1e-3)
        img_layer = viewer.add_image(
            channel_data,
            name=f"Channel {ch}",
            contrast_limits=clim,
            colormap=channel_colors[ch % len(channel_colors)],
            blending="additive",
            opacity=0.8,
        )
        img_layer.contrast_limits_range = (0, float(channel_data.max()))

    # --- Annotation label layers per cell type ---
    # Labels must be (T_total, Z, Y, X) — NOT (C, Z, Y, X)
    label_shape = (T_total,) + stacked.shape[2:]  # (T_total, Z, Y, X)

    ip = initial_params or {}

    # 1. User Provided Labels (all cell types)
    for cell_type in all_cell_types:
        saved_path = Path(pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
        if saved_path.exists():
            existing = np.asarray(load_zarr(saved_path))
            if existing.shape == label_shape:
                user_labels = existing
                print(f"  ↩ Restored saved labels for '{cell_type}' ({label_shape})")
            else:
                print(f"  ⚠️ Saved labels shape {existing.shape} ≠ expected {label_shape} — starting fresh")
                user_labels = np.zeros(label_shape, dtype=np.int16)
        else:
            user_labels = np.zeros(label_shape, dtype=np.int16)

        viewer.add_labels(
            user_labels,
            name=f"User Provided Labels ({cell_type.capitalize()})",
            opacity=0.5,
        )

    if has_death:
        dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
        if dead_path.exists():
            dead_labels = np.asarray(load_zarr(dead_path))
            if dead_labels.shape == label_shape:
                print(f"  ↩ Restored saved labels for 'dead' ({label_shape})")
            else:
                print(f"  ⚠️ Saved Death labels shape {dead_labels.shape} ≠ {label_shape} — starting fresh")
                dead_labels = np.zeros(label_shape, dtype=np.int16)
        else:
            dead_labels = np.zeros(label_shape, dtype=np.int16)
        viewer.add_labels(dead_labels, name="User Provided Labels (Dead)", opacity=0.5)

    # 2. Results (Pixel Classification / Segments) on top
    # Dead result
    if has_death:
        pred_death_path = Path(pixel_class_outdir, "PixelClassifier_Death_PredictedLabels.zarr")
        if pred_death_path.exists():
            pred_death = np.asarray(load_zarr(pred_death_path))
            if pred_death.shape != label_shape:
                pred_death = np.zeros(label_shape, dtype=np.int16)
        else:
            pred_death = np.zeros(label_shape, dtype=np.int16)
        viewer.add_labels(pred_death, name="Pixel Classification (Dead)", opacity=0.8, visible=False)

    # Cell type segments
    for cell_type in all_cell_types:
        pred_path = Path(pixel_class_outdir, f"PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr")
        if pred_path.exists():
            pred_data = np.asarray(load_zarr(pred_path))
            if pred_data.shape != label_shape:
                pred_data = np.zeros(label_shape, dtype=np.int16)
        else:
            pred_data = np.zeros(label_shape, dtype=np.int16)

        viewer.add_labels(
            pred_data,
            name=f"{cell_type.capitalize()} Segments",
            opacity=0.8,
            visible=False,
        )

    # --- Dock the APOC widget ---
    apoc_widget = APOCTrainingWidget(
        viewer=viewer,
        pixel_class_outdir=pixel_class_outdir,
        all_cell_types=all_cell_types,
        has_death=has_death,
        initial_params=ip,
        on_params_changed=on_params_changed,
    )
    viewer.window.add_dock_widget(apoc_widget, area="right", name="APOC Pixel Classification")

    log_output = QPlainTextEdit()
    log_output.setReadOnly(True)
    log_output.setMaximumHeight(120)
    log_widget = QWidget()
    log_layout = QVBoxLayout()
    log_layout.addWidget(log_output)
    log_widget.setLayout(log_layout)
    viewer.window.add_dock_widget(log_widget, area="right", name="Log Output")

    save_button = QtPushButton("💾 Save User Labels")
    save_button.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 6px;")
    save_button.clicked.connect(lambda: apoc_widget.save_user_labels(log=log_output.appendPlainText))
    apoc_widget.layout().addWidget(save_button)

    return viewer
