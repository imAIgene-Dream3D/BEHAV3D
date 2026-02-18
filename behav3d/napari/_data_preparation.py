"""
Data Preparation tab for BEHAV3D napari plugin.

Combines:
  1. Output-directory picker
  2. Metadata Builder  (create / edit metadata CSV)
  3. Metadata Loader   (load + validate existing CSV)
  4. Dimension Order table
  5. Convert-to-Zarr button

All backed by the existing behav3d.core / behav3d.io / behav3d.preprocessing
functions — no logic duplication.
"""
from __future__ import annotations

import traceback
from copy import deepcopy
from pathlib import Path

import pandas as pd
import yaml
from qtpy.QtCore import Qt, Signal, QThread
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QCheckBox,
    QFileDialog,
    QScrollArea,
    QMessageBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QSplitter,
    QTabWidget,
)

# ---------------------------------------------------------------------------
# Core imports (these are the existing, untouched BEHAV3D modules)
# ---------------------------------------------------------------------------
from behav3d.core.metadata import (
    load_behav3d_metadata,
    check_behav3d_metadata,
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)

# These are imported lazily by the conversion worker to keep startup fast.
# from behav3d.preprocessing import convert_input_files_to_zarr
# from behav3d.io.images import load_image, get_image_shape, get_image_dimension_order

# Default config (mirrors widgets/utils.py _DEFAULT_CONFIG)
_DEFAULT_CONFIG: dict = {
    "seed": 42,
    "paths": {"metadata_csv": "", "output_dir": ""},
    "dim_order": {"default_apply_all": "TCZYX"},
    "pixel_classifier": {
        "examples_per_sample": 3,
        "sample_specific_classifier": False,
        "workers": 8,
        "use_all_timepoints": True,
        "tp_start": 0,
        "tp_end": 0,
        "overwrite_existing": False,
    },
}

DIM_ORDER_OPTIONS = ["TCZYX", "TZCYX", "ZCTYX", "ZTCYX", "CZTYX", "CTZYX"]


# ═══════════════════════════════════════════════════════════════════════════
# Helper: load / save YAML config
# ═══════════════════════════════════════════════════════════════════════════
def _load_config(path: Path) -> dict:
    if path.exists():
        with open(path) as f:
            return yaml.safe_load(f) or {}
    return {}


def _deep_merge(base: dict, override: dict) -> dict:
    out = base.copy()
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


# ═══════════════════════════════════════════════════════════════════════════
# Background worker for zarr conversion (so the GUI stays responsive)
# ═══════════════════════════════════════════════════════════════════════════
class _ZarrWorker(QThread):
    progress = Signal(str)
    finished = Signal(bool, str)  # (success, message)

    def __init__(self, output_dir: str, metadata: pd.DataFrame, parent=None):
        super().__init__(parent)
        self.output_dir = output_dir
        self.metadata = metadata

    def run(self):
        try:
            from behav3d.preprocessing import convert_input_files_to_zarr

            self.progress.emit("Starting zarr conversion…")
            convert_input_files_to_zarr(
                output_dir=self.output_dir,
                metadata=self.metadata,
            )
            self.finished.emit(True, "✅ Zarr conversion complete!")
        except Exception:
            self.finished.emit(False, f"❌ Conversion failed:\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
# DataPreparationTab
# ═══════════════════════════════════════════════════════════════════════════
class DataPreparationTab(QWidget):
    """Combined metadata + dim-order + zarr tab for the BEHAV3D napari plugin."""

    # Signal emitted when metadata is successfully loaded
    metadata_loaded = Signal(object)  # pd.DataFrame

    def __init__(self, parent=None):
        super().__init__(parent)

        # State -----------------------------------------------------------
        self.metadata: pd.DataFrame | None = None
        self.output_dir: str = ""
        self.behav3d_parameters: dict = deepcopy(_DEFAULT_CONFIG)
        self._zarr_worker: _ZarrWorker | None = None

        # Build UI --------------------------------------------------------
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        self._layout = QVBoxLayout(inner)
        self._layout.setContentsMargins(4, 4, 4, 4)
        self._layout.setSpacing(4)
        self._layout.setAlignment(Qt.AlignTop)
        scroll.setWidget(inner)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)

        self._build_output_dir_section()
        self._build_metadata_builder_section()
        self._build_metadata_loader_section()
        self._build_dim_order_section()
        self._build_zarr_section()

        # Status log at the bottom ----------------------------------------
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(120)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        self._layout.addWidget(QLabel("Log"))
        self._layout.addWidget(self.log)

    # ------------------------------------------------------------------
    # Logging helper
    # ------------------------------------------------------------------
    def _log(self, msg: str):
        self.log.append(msg)
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    # ══════════════════════════════════════════════════════════════════════
    # Section 1 – Output Directory
    # ══════════════════════════════════════════════════════════════════════
    def _build_output_dir_section(self):
        grp = QGroupBox("1 · Output Directory")
        lay = QHBoxLayout(grp)

        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText("Select output directory…")
        btn = QPushButton("Browse…")
        btn.clicked.connect(self._browse_output_dir)
        lay.addWidget(self.output_dir_edit, stretch=1)
        lay.addWidget(btn)
        self._layout.addWidget(grp)

    def _browse_output_dir(self):
        d = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if d:
            self.output_dir_edit.setText(d)
            self.output_dir = d

    # ══════════════════════════════════════════════════════════════════════
    # Section 2 – Metadata Builder
    # ══════════════════════════════════════════════════════════════════════
    def _build_metadata_builder_section(self):
        grp = QGroupBox("2 · Metadata Builder")
        grp.setCheckable(True)
        grp.setChecked(False)  # collapsed by default
        lay = QVBoxLayout(grp)

        # --- Number of samples ---
        row = QHBoxLayout()
        row.addWidget(QLabel("Number of samples:"))
        self.n_samples_spin = QSpinBox()
        self.n_samples_spin.setMinimum(1)
        self.n_samples_spin.setMaximum(999)
        self.n_samples_spin.setValue(1)
        row.addWidget(self.n_samples_spin)
        lay.addLayout(row)

        # --- Population counts ---
        pop_form = QFormLayout()
        pop_form.setContentsMargins(4, 2, 4, 2)
        pop_form.setSpacing(3)
        pop_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.n_organoid_spin = QSpinBox(); self.n_organoid_spin.setMinimum(0); self.n_organoid_spin.setValue(0); self.n_organoid_spin.setMaximumWidth(60)
        self.n_immune_spin = QSpinBox(); self.n_immune_spin.setMinimum(0); self.n_immune_spin.setValue(0); self.n_immune_spin.setMaximumWidth(60)
        self.n_other_spin = QSpinBox(); self.n_other_spin.setMinimum(0); self.n_other_spin.setValue(0); self.n_other_spin.setMaximumWidth(60)
        pop_form.addRow("Organoid types:", self.n_organoid_spin)
        pop_form.addRow("Immune types:", self.n_immune_spin)
        pop_form.addRow("Other types:", self.n_other_spin)
        self.include_dead_cb = QCheckBox("Include dead channel")
        pop_form.addRow(self.include_dead_cb)
        lay.addLayout(pop_form)

        # --- Cell type naming inputs (dynamically rebuilt) ---
        self.cell_type_naming_container = QVBoxLayout()
        lay.addLayout(self.cell_type_naming_container)

        btn_configure = QPushButton("Configure Cell Types")
        btn_configure.clicked.connect(self._on_configure_cell_types)
        lay.addWidget(btn_configure)

        # --- Sample form area (dynamically rebuilt) ---
        self.sample_form_container = QVBoxLayout()
        lay.addLayout(self.sample_form_container)

        # --- Fill-down + Save ---
        btns = QHBoxLayout()
        btn_fill = QPushButton("Fill All from Sample 1")
        btn_fill.clicked.connect(self._on_fill_down)
        btn_save = QPushButton("Save Metadata CSV")
        btn_save.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        btn_save.clicked.connect(self._on_save_metadata)
        btns.addWidget(btn_fill)
        btns.addWidget(btn_save)
        lay.addLayout(btns)

        self._layout.addWidget(grp)

        # Internal storage
        self._organoid_name_edits: list[QLineEdit] = []
        self._immune_name_edits: list[QLineEdit] = []
        self._other_name_edits: list[QLineEdit] = []
        self._sample_forms: list[dict] = []

    # --- configure cell-type naming fields --------------------------------
    def _on_configure_cell_types(self):
        # Clear old
        self._clear_layout(self.cell_type_naming_container)

        self._organoid_name_edits = []
        self._immune_name_edits = []
        self._other_name_edits = []

        for i in range(self.n_organoid_spin.value()):
            e = QLineEdit(f"organoid{i + 1}")
            e.setPlaceholderText(f"Organoid type {i + 1} name")
            self.cell_type_naming_container.addWidget(QLabel(f"Organoid {i + 1}:"))
            self.cell_type_naming_container.addWidget(e)
            self._organoid_name_edits.append(e)

        for i in range(self.n_immune_spin.value()):
            default = "tcell" if i == 0 else f"immune{i + 1}"
            e = QLineEdit(default)
            e.setPlaceholderText(f"Immune type {i + 1} name")
            self.cell_type_naming_container.addWidget(QLabel(f"Immune {i + 1}:"))
            self.cell_type_naming_container.addWidget(e)
            self._immune_name_edits.append(e)

        for i in range(self.n_other_spin.value()):
            e = QLineEdit(f"other{i + 1}")
            e.setPlaceholderText(f"Other type {i + 1} name")
            self.cell_type_naming_container.addWidget(QLabel(f"Other {i + 1}:"))
            self.cell_type_naming_container.addWidget(e)
            self._other_name_edits.append(e)

        # After naming is set, build sample forms
        btn = QPushButton("Create Sample Forms")
        btn.clicked.connect(self._build_sample_forms)
        self.cell_type_naming_container.addWidget(btn)

    def _build_sample_forms(self):
        self._clear_layout(self.sample_form_container)
        self._sample_forms = []

        org_names = [e.text().strip() for e in self._organoid_name_edits]
        imm_names = [e.text().strip() for e in self._immune_name_edits]
        oth_names = [e.text().strip() for e in self._other_name_edits]
        include_dead = self.include_dead_cb.isChecked()

        for idx in range(self.n_samples_spin.value()):
            form = self._create_sample_form(idx, org_names, imm_names, oth_names, include_dead)
            self._sample_forms.append(form)
            self.sample_form_container.addWidget(form["group"])

    def _create_sample_form(self, idx: int, org_names, imm_names, oth_names, include_dead) -> dict:
        grp = QGroupBox(f"Sample {idx + 1}")
        layout = QFormLayout(grp)
        layout.setContentsMargins(6, 4, 6, 4)
        layout.setSpacing(3)
        layout.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

        fields: dict[str, QWidget] = {}

        # Basic fields
        fields["sample_name"] = QLineEdit()
        fields["sample_name"].setPlaceholderText("e.g. Sample001")
        layout.addRow("Name*:", fields["sample_name"])

        fields["exp_nr"] = QSpinBox(); fields["exp_nr"].setMaximum(9999); fields["exp_nr"].setValue(1); fields["exp_nr"].setMaximumWidth(70)
        layout.addRow("Exp #*:", fields["exp_nr"])

        fields["well"] = QLineEdit(); fields["well"].setPlaceholderText("e.g. A1")
        layout.addRow("Well*:", fields["well"])

        fields["raw_image_path"] = QLineEdit()
        fields["raw_image_path"].setPlaceholderText("/path/to/image.czi")
        layout.addRow("Image path*:", fields["raw_image_path"])

        fields["dimension_order"] = QLineEdit()
        fields["dimension_order"].setPlaceholderText("e.g. TCZYX")
        layout.addRow("Dim order:", fields["dimension_order"])

        fields["pixel_distance_xy"] = QDoubleSpinBox(); fields["pixel_distance_xy"].setDecimals(4)
        fields["pixel_distance_xy"].setMaximum(9999); fields["pixel_distance_xy"].setValue(0.5); fields["pixel_distance_xy"].setMaximumWidth(90)
        layout.addRow("Pixel XY (μm)*:", fields["pixel_distance_xy"])

        fields["pixel_distance_z"] = QDoubleSpinBox(); fields["pixel_distance_z"].setDecimals(4)
        fields["pixel_distance_z"].setMaximum(9999); fields["pixel_distance_z"].setValue(2.0); fields["pixel_distance_z"].setMaximumWidth(90)
        layout.addRow("Pixel Z (μm)*:", fields["pixel_distance_z"])

        fields["distance_unit"] = QLineEdit("μm")
        fields["distance_unit"].setMaximumWidth(60)
        layout.addRow("Dist unit*:", fields["distance_unit"])

        fields["time_interval"] = QDoubleSpinBox(); fields["time_interval"].setDecimals(4)
        fields["time_interval"].setMaximum(99999); fields["time_interval"].setValue(1.0); fields["time_interval"].setMaximumWidth(90)
        layout.addRow("Time interval*:", fields["time_interval"])

        fields["time_unit"] = QLineEdit("s")
        fields["time_unit"].setMaximumWidth(60)
        layout.addRow("Time unit*:", fields["time_unit"])

        # Dead channel
        dead_fields: dict[str, QWidget] = {}
        if include_dead:
            dead_fields["number"] = QSpinBox(); dead_fields["number"].setMinimum(0); dead_fields["number"].setMaximumWidth(60)
            layout.addRow("Dead ch #:", dead_fields["number"])
            dead_fields["mask_path"] = QLineEdit()
            dead_fields["mask_path"].setPlaceholderText("Optional — dead mask path")
            layout.addRow("Dead mask:", dead_fields["mask_path"])

        # Cell-type fields
        cell_type_fields: dict[str, dict[str, QLineEdit]] = {}

        def _add_ct(name, prefix):
            layout.addRow(QLabel(f"── {name} ──"), QLabel(""))
            ct: dict[str, QLineEdit] = {}
            ct["line"] = QLineEdit(); ct["line"].setPlaceholderText("Line")
            layout.addRow(f"  Line:", ct["line"])
            ct["condition"] = QLineEdit(); ct["condition"].setPlaceholderText("Condition")
            layout.addRow(f"  Condition:", ct["condition"])
            ct["segments_image_path"] = QLineEdit(); ct["segments_image_path"].setPlaceholderText("Optional")
            layout.addRow(f"  Segments:", ct["segments_image_path"])
            ct["tracks_image_path"] = QLineEdit(); ct["tracks_image_path"].setPlaceholderText("Optional")
            layout.addRow(f"  Tracks img:", ct["tracks_image_path"])
            ct["tracks_csv_path"] = QLineEdit(); ct["tracks_csv_path"].setPlaceholderText("Optional")
            layout.addRow(f"  Tracks csv:", ct["tracks_csv_path"])
            cell_type_fields[name] = ct

        for n in org_names:
            _add_ct(n, "or")
        for n in imm_names:
            _add_ct(n, "im")
        for n in oth_names:
            _add_ct(n, "ot")

        return {
            "group": grp,
            "basic": fields,
            "dead_channel": dead_fields,
            "cell_types": cell_type_fields,
            "org_names": org_names,
            "imm_names": imm_names,
            "oth_names": oth_names,
        }

    # --- Fill-down --------------------------------------------------------
    def _on_fill_down(self):
        if len(self._sample_forms) < 2:
            return
        src = self._sample_forms[0]
        for tgt in self._sample_forms[1:]:
            for k, w in src["basic"].items():
                if k == "sample_name":
                    continue
                tw = tgt["basic"][k]
                if isinstance(w, QLineEdit):
                    tw.setText(w.text())
                elif isinstance(w, (QSpinBox, QDoubleSpinBox)):
                    tw.setValue(w.value())
            for ct, fields in src["cell_types"].items():
                if ct in tgt["cell_types"]:
                    for fk, fw in fields.items():
                        tgt["cell_types"][ct][fk].setText(fw.text())
            for dk, dw in src["dead_channel"].items():
                if dk in tgt["dead_channel"]:
                    tw = tgt["dead_channel"][dk]
                    if isinstance(dw, QSpinBox):
                        tw.setValue(dw.value())
                    elif isinstance(dw, QLineEdit):
                        tw.setText(dw.text())
        self._log("✅ Filled all samples from Sample 1")

    # --- Save metadata CSV ------------------------------------------------
    def _on_save_metadata(self):
        out_dir = self.output_dir_edit.text().strip()
        if not out_dir or not Path(out_dir).exists():
            QMessageBox.warning(self, "Error", "Please set a valid output directory first.")
            return

        org_names = [e.text().strip() for e in self._organoid_name_edits]
        imm_names = [e.text().strip() for e in self._immune_name_edits]
        oth_names = [e.text().strip() for e in self._other_name_edits]
        include_dead = self.include_dead_cb.isChecked()

        rows = []
        for form in self._sample_forms:
            row: dict = {}
            for k, w in form["basic"].items():
                if isinstance(w, QLineEdit):
                    row[k] = w.text().strip().strip('"').strip("'")
                elif isinstance(w, QSpinBox):
                    row[k] = w.value()
                elif isinstance(w, QDoubleSpinBox):
                    row[k] = w.value()

            if include_dead:
                if "number" in form["dead_channel"]:
                    row["dead_channel"] = form["dead_channel"]["number"].value()
                if "mask_path" in form["dead_channel"]:
                    row["dead_mask_path"] = form["dead_channel"]["mask_path"].text().strip()

            for ct, fields in form["cell_types"].items():
                if ct in org_names:
                    pfx = "or"
                elif ct in imm_names:
                    pfx = "im"
                elif ct in oth_names:
                    pfx = "ot"
                else:
                    continue
                line = fields["line"].text().strip()
                cond = fields["condition"].text().strip()
                if line and cond:
                    row[f"{pfx}_{ct}_line_condition"] = f"{line}_{cond}"
                elif line:
                    row[f"{pfx}_{ct}_line_condition"] = line
                else:
                    row[f"{pfx}_{ct}_line_condition"] = ""
                row[f"{pfx}_{ct}_segments_image_path"] = fields["segments_image_path"].text().strip().strip('"').strip("'")
                row[f"{pfx}_{ct}_tracks_image_path"] = fields["tracks_image_path"].text().strip().strip('"').strip("'")
                row[f"{pfx}_{ct}_tracks_csv_path"] = fields["tracks_csv_path"].text().strip().strip('"').strip("'")

            rows.append(row)

        df = pd.DataFrame(rows)
        csv_path = Path(out_dir) / "metadata.csv"
        df.to_csv(csv_path, index=False)
        self._log(f"✅ Metadata saved to {csv_path}  ({len(df)} samples, {len(df.columns)} columns)")

    # ══════════════════════════════════════════════════════════════════════
    # Section 3 – Metadata Loader
    # ══════════════════════════════════════════════════════════════════════
    def _build_metadata_loader_section(self):
        grp = QGroupBox("3 · Metadata Loader")
        lay = QVBoxLayout(grp)

        row = QHBoxLayout()
        self.csv_path_edit = QLineEdit()
        self.csv_path_edit.setPlaceholderText("Path to metadata.csv")
        btn_browse = QPushButton("Browse…")
        btn_browse.clicked.connect(self._browse_csv)
        row.addWidget(self.csv_path_edit, stretch=1)
        row.addWidget(btn_browse)
        lay.addLayout(row)

        btn_load = QPushButton("Load Metadata")
        btn_load.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold;")
        btn_load.clicked.connect(self._on_load_metadata)
        lay.addWidget(btn_load)

        self.metadata_info_label = QLabel("")
        self.metadata_info_label.setWordWrap(True)
        lay.addWidget(self.metadata_info_label)

        self._layout.addWidget(grp)

    def _browse_csv(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select metadata CSV", "", "CSV Files (*.csv)")
        if path:
            self.csv_path_edit.setText(path)

    def _on_load_metadata(self):
        csv_path = self.csv_path_edit.text().strip()
        out_dir = self.output_dir_edit.text().strip()

        if not csv_path or not csv_path.lower().endswith(".csv"):
            QMessageBox.warning(self, "Error", "Please select a valid .csv file.")
            return
        if not Path(csv_path).exists():
            QMessageBox.warning(self, "Error", f"File not found: {csv_path}")
            return

        try:
            self.metadata = load_behav3d_metadata(csv_path)
            self.output_dir = out_dir
            check_behav3d_metadata(self.metadata, func=False)

            # Persist parameters
            params_path = Path(out_dir) / "behav3d_parameters.yml" if out_dir else None
            if params_path:
                if params_path.exists():
                    self.behav3d_parameters = _deep_merge(deepcopy(_DEFAULT_CONFIG), _load_config(params_path))
                else:
                    self.behav3d_parameters = deepcopy(_DEFAULT_CONFIG)
                self.behav3d_parameters["paths"]["metadata_csv"] = csv_path
                self.behav3d_parameters["paths"]["output_dir"] = out_dir
                with open(params_path, "w") as f:
                    yaml.safe_dump(self.behav3d_parameters, f, sort_keys=False)

            n = len(self.metadata)
            cols = len(self.metadata.columns)
            org = detect_organoid_types_from_metadata(self.metadata)
            imm = detect_immune_cell_types_from_metadata(self.metadata)
            oth = detect_other_cell_types_from_metadata(self.metadata)
            info_parts = [f"✅ Loaded {n} samples, {cols} columns."]
            if org:
                info_parts.append(f"Organoid: {', '.join(org)}")
            if imm:
                info_parts.append(f"Immune: {', '.join(imm)}")
            if oth:
                info_parts.append(f"Other: {', '.join(oth)}")
            self.metadata_info_label.setText("  |  ".join(info_parts))
            self._log(f"✅ Metadata loaded from {csv_path}")

            # Populate dim-order table
            self._populate_dim_order_table()

            # Emit signal for other tabs
            self.metadata_loaded.emit(self.metadata)

        except Exception as e:
            self.metadata_info_label.setText(f"❌ {e}")
            self._log(f"❌ Error loading metadata: {e}")

    # ══════════════════════════════════════════════════════════════════════
    # Section 4 – Dimension Order Table
    # ══════════════════════════════════════════════════════════════════════
    def _build_dim_order_section(self):
        grp = QGroupBox("4 · Dimension Order")
        lay = QVBoxLayout(grp)

        self.dim_table = QTableWidget(0, 3)
        self.dim_table.setHorizontalHeaderLabels(["Sample", "Shape", "Dim Order"])
        self.dim_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.dim_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.dim_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.dim_table.setMaximumHeight(200)
        lay.addWidget(self.dim_table)

        row = QHBoxLayout()
        self.dim_apply_all_combo = QComboBox()
        self.dim_apply_all_combo.addItems(DIM_ORDER_OPTIONS)
        btn_apply_all = QPushButton("Apply to All")
        btn_apply_all.clicked.connect(self._apply_dim_order_all)
        btn_save_dim = QPushButton("Save to Metadata CSV")
        btn_save_dim.clicked.connect(self._save_dim_order)
        row.addWidget(QLabel("Set all:"))
        row.addWidget(self.dim_apply_all_combo)
        row.addWidget(btn_apply_all)
        row.addWidget(btn_save_dim)
        lay.addLayout(row)

        self._layout.addWidget(grp)

    def _populate_dim_order_table(self):
        if self.metadata is None:
            return

        from behav3d.io.images import get_image_shape, get_image_dimension_order

        samples = self.metadata["sample_name"].tolist()
        self.dim_table.setRowCount(len(samples))

        for i, sample in enumerate(samples):
            # Sample name (read-only)
            item = QTableWidgetItem(str(sample))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.dim_table.setItem(i, 0, item)

            # Shape
            raw_path = self.metadata.loc[self.metadata["sample_name"] == sample, "raw_image_path"]
            shape_str = ""
            detected_order = "TCZYX"
            if not raw_path.empty:
                p = Path(str(raw_path.iloc[0]).strip())
                if p.exists():
                    try:
                        shape_str = str(get_image_shape(p))
                    except Exception:
                        shape_str = "?"
                    try:
                        detected_order = get_image_dimension_order(p)
                    except Exception:
                        pass

            shape_item = QTableWidgetItem(shape_str)
            shape_item.setFlags(shape_item.flags() & ~Qt.ItemIsEditable)
            self.dim_table.setItem(i, 1, shape_item)

            # Dim order combo
            combo = QComboBox()
            combo.addItems(DIM_ORDER_OPTIONS)
            # Use value from metadata if present
            existing = self.metadata.loc[
                self.metadata["sample_name"] == sample, "dimension_order"
            ]
            if not existing.empty and pd.notna(existing.iloc[0]) and str(existing.iloc[0]).strip():
                val = str(existing.iloc[0]).strip()
                if val in DIM_ORDER_OPTIONS:
                    combo.setCurrentText(val)
                else:
                    combo.setCurrentText(detected_order if detected_order in DIM_ORDER_OPTIONS else "TCZYX")
            else:
                combo.setCurrentText(detected_order if detected_order in DIM_ORDER_OPTIONS else "TCZYX")
            self.dim_table.setCellWidget(i, 2, combo)

    def _apply_dim_order_all(self):
        order = self.dim_apply_all_combo.currentText()
        for i in range(self.dim_table.rowCount()):
            combo = self.dim_table.cellWidget(i, 2)
            if combo:
                combo.setCurrentText(order)

    def _save_dim_order(self):
        if self.metadata is None:
            return
        csv_path = self.csv_path_edit.text().strip()
        if not csv_path:
            return

        for i in range(self.dim_table.rowCount()):
            sample = self.dim_table.item(i, 0).text()
            combo = self.dim_table.cellWidget(i, 2)
            if combo:
                self.metadata.loc[
                    self.metadata["sample_name"] == sample, "dimension_order"
                ] = combo.currentText()

        self.metadata.to_csv(csv_path, index=False)
        self._log(f"✅ Dimension orders saved to {csv_path}")

    # ══════════════════════════════════════════════════════════════════════
    # Section 5 – Convert to Zarr
    # ══════════════════════════════════════════════════════════════════════
    def _build_zarr_section(self):
        grp = QGroupBox("5 · Convert to Zarr")
        lay = QVBoxLayout(grp)

        self.zarr_btn = QPushButton("Convert All Images to Zarr")
        self.zarr_btn.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold;")
        self.zarr_btn.clicked.connect(self._on_convert_zarr)
        lay.addWidget(self.zarr_btn)

        self.zarr_status = QLabel("")
        self.zarr_status.setWordWrap(True)
        lay.addWidget(self.zarr_status)

        self._layout.addWidget(grp)

    def _on_convert_zarr(self):
        out_dir = self.output_dir_edit.text().strip()
        if self.metadata is None:
            QMessageBox.warning(self, "Error", "Please load metadata first.")
            return
        if not out_dir:
            QMessageBox.warning(self, "Error", "Please set output directory first.")
            return

        # Save dim orders first
        self._save_dim_order()

        self.zarr_btn.setEnabled(False)
        self.zarr_status.setText("⏳ Converting…")
        self._log("Starting zarr conversion…")

        self._zarr_worker = _ZarrWorker(out_dir, self.metadata.copy(), parent=self)
        self._zarr_worker.progress.connect(lambda msg: self._log(msg))
        self._zarr_worker.finished.connect(self._on_zarr_done)
        self._zarr_worker.start()

    def _on_zarr_done(self, success: bool, message: str):
        self.zarr_btn.setEnabled(True)
        self.zarr_status.setText(message.split("\n")[0])
        self._log(message)

    # ══════════════════════════════════════════════════════════════════════
    # Utility
    # ══════════════════════════════════════════════════════════════════════
    @staticmethod
    def _clear_layout(layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
            elif child.layout():
                DataPreparationTab._clear_layout(child.layout())
