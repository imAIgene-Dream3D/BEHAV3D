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

import os
import shutil
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
    # (success, message, updated_metadata_df_or_None, originals_map_or_None)
    # originals_map is a list of dicts: {"sample": str, "original": str, "new_zarr": str}
    finished = Signal(bool, str, object, object)

    def __init__(self, output_dir: str, metadata: pd.DataFrame,
                 t_start: int = None, t_end: int = None,
                 n_workers: int = 1, parent=None):
        super().__init__(parent)
        self.output_dir = output_dir
        self.metadata = metadata
        self.t_start = t_start
        self.t_end = t_end
        self.n_workers = n_workers

    def run(self):
        try:
            from behav3d.preprocessing import convert_input_files_to_zarr

            # Snapshot of (sample_name, original raw_image_path) BEFORE the call,
            # because convert_input_files_to_zarr rewrites raw_image_path in place.
            originals_snapshot = [
                {
                    "sample": str(row.get("sample_name", "")),
                    "original": str(row.get("raw_image_path", "")),
                }
                for _, row in self.metadata.iterrows()
            ]

            self.progress.emit("Starting zarr conversion…")
            updated_md = convert_input_files_to_zarr(
                output_dir=self.output_dir,
                metadata=self.metadata,
                t_start=self.t_start,
                t_end=self.t_end,
                n_workers=self.n_workers,
            )

            # Pair each original with the new zarr path now stored in raw_image_path.
            originals_map = []
            new_paths_by_sample = {
                str(row.get("sample_name", "")): str(row.get("raw_image_path", ""))
                for _, row in updated_md.iterrows()
            }
            for entry in originals_snapshot:
                new_zarr = new_paths_by_sample.get(entry["sample"], "")
                originals_map.append({
                    "sample": entry["sample"],
                    "original": entry["original"],
                    "new_zarr": new_zarr,
                })

            self.finished.emit(True, "✅ Zarr conversion complete!", updated_md, originals_map)
        except Exception:
            self.finished.emit(
                False,
                f"❌ Conversion failed:\n{traceback.format_exc()}",
                None,
                None,
            )


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
        self._edit_mode: bool = False          # True when editing loaded metadata
        self._loaded_csv_path: str = ""        # CSV path of last loaded metadata

        # Build UI --------------------------------------------------------
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
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
        self._build_metadata_overview_section()
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
        self.output_dir_edit.editingFinished.connect(self._on_output_dir_edited)
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

    def _on_output_dir_edited(self):
        self.output_dir = self.output_dir_edit.text().strip()

    # ══════════════════════════════════════════════════════════════════════
    # Section 2 – Metadata Builder
    # ══════════════════════════════════════════════════════════════════════
    def _build_metadata_builder_section(self):
        self.builder_grp = QGroupBox("2 · Metadata Builder")
        self.builder_grp.setCheckable(True)
        self.builder_grp.setChecked(False)  # collapsed by default
        self.builder_grp.toggled.connect(self._on_builder_toggled)
        lay = QVBoxLayout(self.builder_grp)

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
        btn_configure.clicked.connect(lambda: self._on_configure_cell_types(force=True))
        lay.addWidget(btn_configure)

        # --- Sample form area (dynamically rebuilt) ---
        self.sample_form_container = QVBoxLayout()
        lay.addLayout(self.sample_form_container)

        # --- Fill-down + Save (Horizontal Row) ---
        btns = QHBoxLayout()
        btn_fill = QPushButton("Fill All from Sample 1")
        btn_fill.clicked.connect(self._on_fill_down)
        
        self.btn_save_metadata = QPushButton("Save Metadata CSV")
        self.btn_save_metadata.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        self.btn_save_metadata.clicked.connect(self._on_save_metadata)
        
        btns.addWidget(btn_fill)
        btns.addWidget(self.btn_save_metadata)
        lay.addLayout(btns)

        self._layout.addWidget(self.builder_grp)

        # Internal storage
        self._organoid_name_edits: list[QLineEdit] = []
        self._immune_name_edits: list[QLineEdit] = []
        self._other_name_edits: list[QLineEdit] = []
        self._sample_forms: list[dict] = []

    # --- configure cell-type naming fields --------------------------------
    def _on_configure_cell_types(self, force: bool = False):
        # Idempotent: when the existing naming fields already match the configured
        # counts, skip the destructive rebuild so entered names survive. The
        # assistant may re-trigger this structurally; user button clicks pass
        # force=True to always rebuild.
        if (not force
                and (self._organoid_name_edits or self._immune_name_edits
                     or self._other_name_edits)
                and len(self._organoid_name_edits) == self.n_organoid_spin.value()
                and len(self._immune_name_edits) == self.n_immune_spin.value()
                and len(self._other_name_edits) == self.n_other_spin.value()):
            return
        # Clear old
        self._clear_layout(self.cell_type_naming_container)

        self._organoid_name_edits = []
        self._immune_name_edits = []
        # Multicolor controls parallel to ``_immune_name_edits``.
        self._immune_multicolor_flags: list[QCheckBox] = []
        self._immune_multicolor_counts: list[QSpinBox] = []
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

            # ── Multicolor controls (one row per immune type) ─────────
            multi_row = QHBoxLayout()
            multi_row.setContentsMargins(0, 0, 0, 0)
            multi_chk = QCheckBox("Multicolor")
            multi_chk.setToolTip(
                "Tick when this immune cell type is split across multiple fluorescent\n"
                "channels (e.g. red/green T cells). Each channel becomes its own\n"
                "intermediate cell type named '{base}_{n}_multicolor', and feature\n"
                "extraction operates on the aggregated '{base}_merged' output."
            )
            multi_n = QSpinBox()
            multi_n.setRange(2, 10)
            multi_n.setValue(2)
            multi_n.setMaximumWidth(70)
            multi_n.setSuffix(" ch")
            multi_n.setEnabled(False)

            def _toggle_multi(state, spin=multi_n):
                spin.setEnabled(bool(state))

            multi_chk.stateChanged.connect(_toggle_multi)
            multi_row.addWidget(multi_chk)
            multi_row.addWidget(multi_n)
            multi_row.addStretch()
            self.cell_type_naming_container.addLayout(multi_row)
            self._immune_multicolor_flags.append(multi_chk)
            self._immune_multicolor_counts.append(multi_n)

        for i in range(self.n_other_spin.value()):
            e = QLineEdit(f"other{i + 1}")
            e.setPlaceholderText(f"Other type {i + 1} name")
            self.cell_type_naming_container.addWidget(QLabel(f"Other {i + 1}:"))
            self.cell_type_naming_container.addWidget(e)
            self._other_name_edits.append(e)

        # After naming is set, build sample forms
        btn = QPushButton("Create Sample Forms")
        btn.clicked.connect(lambda: self._build_sample_forms(force=True))
        self.cell_type_naming_container.addWidget(btn)

    def _sample_forms_match_current(self) -> bool:
        """True if the existing sample forms already match the configured sample
        count, cell-type names and dead-channel flag — so a rebuild would only
        destroy values the user (or assistant) already entered."""
        forms = getattr(self, "_sample_forms", [])
        if not forms or len(forms) != self.n_samples_spin.value():
            return False
        org = [e.text().strip() for e in self._organoid_name_edits]
        imm = self._expanded_immune_names()
        oth = [e.text().strip() for e in self._other_name_edits]
        include_dead = self.include_dead_cb.isChecked()
        for form in forms:
            if (form.get("org_names") != org or form.get("imm_names") != imm
                    or form.get("oth_names") != oth):
                return False
            if include_dead != bool(form.get("dead_channel")):
                return False
        return True

    def _build_sample_forms(self, force: bool = False):
        # Idempotent: skip the destructive rebuild when the forms already match the
        # configured structure, so entered values survive structural re-triggers
        # from the assistant. User button clicks pass force=True to always rebuild.
        if not force and self._sample_forms_match_current():
            return
        self._clear_layout(self.sample_form_container)
        self._sample_forms = []

        org_names = [e.text().strip() for e in self._organoid_name_edits]
        imm_names = self._expanded_immune_names()
        oth_names = [e.text().strip() for e in self._other_name_edits]
        include_dead = self.include_dead_cb.isChecked()

        for idx in range(self.n_samples_spin.value()):
            form = self._create_sample_form(idx, org_names, imm_names, oth_names, include_dead)
            self._sample_forms.append(form)
            self.sample_form_container.addWidget(form["group"])

    def _expanded_immune_names(self) -> list[str]:
        """Return immune cell type names with multicolor expansion applied.

        Mirrors :meth:`behav3d.widgets.metadata._on_names_confirmed`: when the
        Multicolor checkbox is ticked next to an immune name, that base name is
        expanded into N entries of the form ``{base}_{n}_multicolor`` (one per
        configured channel). Untouched names pass through unchanged.
        """
        flags = getattr(self, "_immune_multicolor_flags", [])
        counts = getattr(self, "_immune_multicolor_counts", [])
        out: list[str] = []
        for idx, edit in enumerate(self._immune_name_edits):
            base = edit.text().strip()
            if not base:
                continue
            is_multi = (
                idx < len(flags)
                and flags[idx] is not None
                and flags[idx].isChecked()
            )
            if is_multi:
                try:
                    n = int(counts[idx].value()) if idx < len(counts) else 2
                except Exception:
                    n = 2
                n = max(2, n)
                for j in range(n):
                    out.append(f"{base}_{j+1}_multicolor")
            else:
                out.append(base)
        return out

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

        fields["time_unit"] = QComboBox()
        fields["time_unit"].addItems(["s", "m", "h"])
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
                elif isinstance(w, QComboBox):
                    tw.setCurrentText(w.currentText())
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
    def _on_builder_toggled(self, checked):
        """Handle builder toggle — prompt if metadata is already loaded."""
        if not checked:
            self._reset_builder()
            return

        if self.metadata is not None:
            msg = QMessageBox(self)
            msg.setWindowTitle("Metadata Builder")
            msg.setText("Metadata is already loaded. Do you want to edit the current metadata or generate a new one?")
            
            # Custom buttons
            btn_edit = msg.addButton("Edit Current", QMessageBox.AcceptRole)
            btn_new = msg.addButton("Generate New", QMessageBox.DestructiveRole)
            btn_cancel = msg.addButton("Cancel", QMessageBox.RejectRole)
            
            msg.exec_()
            
            clicked = msg.clickedButton()
            if clicked == btn_cancel:
                # Revert toggle
                self.builder_grp.blockSignals(True)
                self.builder_grp.setChecked(False)
                self.builder_grp.blockSignals(False)
                return
            
            if clicked == btn_edit:
                self._enter_edit_mode()
            else: # Generate New
                self._reset_builder()

    def _reset_builder(self):
        """Reset the builder UI and internal state to blank."""
        self._edit_mode = False
        
        # Reset spins without triggering configuration
        for s in [self.n_samples_spin, self.n_organoid_spin, self.n_immune_spin, self.n_other_spin]:
            s.blockSignals(True)
            if s == self.n_samples_spin:
                s.setValue(1)
            else:
                s.setValue(0)
            s.blockSignals(False)

        self.include_dead_cb.blockSignals(True)
        self.include_dead_cb.setChecked(False)
        self.include_dead_cb.blockSignals(False)

        # Clear UI sections
        self._clear_layout(self.cell_type_naming_container)
        self._clear_layout(self.sample_form_container)

        # Clear lists
        self._organoid_name_edits = []
        self._immune_name_edits = []
        self._other_name_edits = []
        self._sample_forms = []

        # Reset button text/style based on current state
        is_checked = self.builder_grp.isChecked()
        if not is_checked and self.metadata is not None:
            # Collapsed + Metadata loaded -> Yellow Edit
            # NOTE: To keep the button active while the group is "collapsed", 
            # we would technically need the group to be checked. 
            # In standard Qt, an unchecked QGroupBox disables all children.
            self.btn_save_metadata.setText("Edit Metadata CSV")
            self.btn_save_metadata.setStyleSheet(
                "background-color: #FFA726; color: white; font-weight: bold;"
            )
        elif is_checked and not self._edit_mode:
            # Open + Not editing -> Green Save New
            self.btn_save_metadata.setText("Save Metadata CSV")
            self.btn_save_metadata.setStyleSheet(
                "background-color: #4CAF50; color: white; font-weight: bold;"
            )
        elif self.metadata is not None and not is_checked:
            # Fallback for yellow state
            self.btn_save_metadata.setText("Edit Metadata CSV")
            self.btn_save_metadata.setStyleSheet(
                "background-color: #FFA726; color: white; font-weight: bold;"
            )
        else:
            # Default Green Save
            self.btn_save_metadata.setText("Save Metadata CSV")
            self.btn_save_metadata.setStyleSheet(
                "background-color: #4CAF50; color: white; font-weight: bold;"
            )

    def _enter_edit_mode(self):
        """Populate builder from loaded metadata and switch to edit mode."""
        self._edit_mode = True
        self._populate_builder_from_metadata()
        self.btn_save_metadata.setText("Save Changes")
        self.btn_save_metadata.setStyleSheet(
            "background-color: #4CAF50; color: white; font-weight: bold;"
        )

    def _populate_builder_from_metadata(self):
        """Fill the metadata builder forms from self.metadata."""
        md = self.metadata
        if md is None or md.empty:
            return

        org_types = detect_organoid_types_from_metadata(md)
        imm_types = detect_immune_cell_types_from_metadata(md)
        oth_types = detect_other_cell_types_from_metadata(md)
        include_dead = "dead_channel" in md.columns

        # Collapse multicolor channel entries (``{base}_{n}_multicolor``) into
        # a base name + channel count so the UI can re-show the checkbox+spin
        # row instead of N separate naming edits.
        imm_bases, imm_multicolor = self._collapse_multicolor_immune_names(imm_types)

        # Set spinners
        self.n_samples_spin.setValue(len(md))
        self.n_organoid_spin.setValue(len(org_types))
        self.n_immune_spin.setValue(len(imm_bases))
        self.n_other_spin.setValue(len(oth_types))
        self.include_dead_cb.setChecked(include_dead)

        # Build cell type naming fields (force: loading replaces any prior structure)
        self._on_configure_cell_types(force=True)

        # Fill naming edits with detected names
        for i, name in enumerate(org_types):
            if i < len(self._organoid_name_edits):
                self._organoid_name_edits[i].setText(name)
        for i, name in enumerate(imm_bases):
            if i < len(self._immune_name_edits):
                self._immune_name_edits[i].setText(name)
        for i, n_channels in enumerate(imm_multicolor):
            if n_channels and n_channels > 1 and i < len(self._immune_multicolor_flags):
                self._immune_multicolor_flags[i].setChecked(True)
                if i < len(self._immune_multicolor_counts):
                    self._immune_multicolor_counts[i].setValue(int(n_channels))
        for i, name in enumerate(oth_types):
            if i < len(self._other_name_edits):
                self._other_name_edits[i].setText(name)

        # Build sample forms (force: loading replaces any prior structure)
        self._build_sample_forms(force=True)

    @staticmethod
    def _collapse_multicolor_immune_names(imm_types):
        """Collapse ``base_N_multicolor`` siblings into ``(bases, n_channels)``.

        Mirrors the notebook builder's expansion logic in reverse so that when
        the user opens an existing CSV containing per-channel multicolor names
        (e.g. ``tcell_1_multicolor`` + ``tcell_2_multicolor``), the builder
        re-shows a single ``tcell`` row with the Multicolor checkbox ticked.

        Non-multicolor names pass through unchanged with ``n_channels=0``.
        """
        import re
        pattern = re.compile(r"^(?P<base>.+)_(?P<idx>\d+)_multicolor$")
        bases: list[str] = []
        counts: dict[str, int] = {}
        order: list[str] = []
        for name in imm_types:
            m = pattern.match(str(name))
            if m:
                base = m.group("base")
                if base not in counts:
                    counts[base] = 0
                    order.append(base)
                counts[base] += 1
            else:
                if name not in counts:
                    counts[name] = 0
                    order.append(name)

        n_channels_per_base = [counts[name] for name in order]
        return order, n_channels_per_base

        # Populate each sample form from the DataFrame row
        for row_idx, (_, row) in enumerate(md.iterrows()):
            if row_idx >= len(self._sample_forms):
                break
            form = self._sample_forms[row_idx]

            # Basic fields
            for key, widget in form["basic"].items():
                val = row.get(key)
                if val is None or pd.isna(val):
                    continue
                if isinstance(widget, QLineEdit):
                    widget.setText(str(val))
                elif isinstance(widget, QComboBox):
                    widget.setCurrentText(str(val))
                elif isinstance(widget, QSpinBox):
                    widget.setValue(int(val))
                elif isinstance(widget, QDoubleSpinBox):
                    widget.setValue(float(val))

            # Dead channel fields
            if include_dead:
                if "number" in form["dead_channel"] and "dead_channel" in row.index:
                    dc_val = row.get("dead_channel")
                    if pd.notna(dc_val):
                        form["dead_channel"]["number"].setValue(int(dc_val))
                if "mask_path" in form["dead_channel"] and "dead_mask_path" in row.index:
                    dm_val = row.get("dead_mask_path")
                    if pd.notna(dm_val):
                        form["dead_channel"]["mask_path"].setText(str(dm_val))

            # Cell type fields
            all_ct = org_types + imm_types + oth_types
            for ct_name in all_ct:
                if ct_name not in form["cell_types"]:
                    continue
                if ct_name in org_types:
                    pfx = "or"
                elif ct_name in imm_types:
                    pfx = "im"
                else:
                    pfx = "ot"

                ct_fields = form["cell_types"][ct_name]

                # line_condition → split into line + condition
                lc_col = f"{pfx}_{ct_name}_line_condition"
                lc_val = row.get(lc_col)
                if pd.notna(lc_val) and str(lc_val).strip():
                    parts = str(lc_val).split("_", 1)
                    ct_fields["line"].setText(parts[0])
                    if len(parts) > 1:
                        ct_fields["condition"].setText(parts[1])

                for fld in ("segments_image_path", "tracks_image_path", "tracks_csv_path"):
                    col = f"{pfx}_{ct_name}_{fld}"
                    v = row.get(col)
                    if pd.notna(v) and str(v).strip():
                        ct_fields[fld].setText(str(v))

        self._log("📝 Metadata loaded into builder for editing.")

    def _on_save_metadata(self):
        # If button is in "Edit" state, open builder in edit mode instead of saving
        if self.btn_save_metadata.text() == "Edit Metadata CSV":
            self.builder_grp.blockSignals(True)
            self.builder_grp.setChecked(True)
            self.builder_grp.blockSignals(False)
            self._enter_edit_mode()
            return

        # Determine save target
        if self._edit_mode and self._loaded_csv_path:
            csv_path = Path(self._loaded_csv_path)
        else:
            out_dir = self.output_dir_edit.text().strip()
            if not out_dir or not Path(out_dir).exists():
                QMessageBox.warning(self, "Error", "Please set a valid output directory first.")
                return
            csv_path = Path(out_dir) / "metadata.csv"

        org_names = [e.text().strip() for e in self._organoid_name_edits]
        imm_names = self._expanded_immune_names()
        oth_names = [e.text().strip() for e in self._other_name_edits]
        include_dead = self.include_dead_cb.isChecked()

        rows = []
        for form in self._sample_forms:
            row: dict = {}
            for k, w in form["basic"].items():
                if isinstance(w, QLineEdit):
                    row[k] = w.text().strip().strip('"').strip("'")
                elif isinstance(w, QComboBox):
                    row[k] = w.currentText()
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

        # Overwrite prompt
        if csv_path.exists():
            res = QMessageBox.question(
                self, "Overwrite?",
                f"A metadata file already exists:\n{csv_path}\n\nDo you want to overwrite it?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                self._log("Save cancelled by user.")
                return

        df.to_csv(csv_path, index=False)
        self.metadata = df
        self._loaded_csv_path = str(csv_path)
        self._populate_metadata_overview()
        self._log(f"✅ Metadata saved to {csv_path}  ({len(df)} samples, {len(df.columns)} columns)")

        # Update behav3d_parameters
        self.behav3d_parameters["paths"]["metadata_csv"] = str(csv_path)

        # If we were in edit mode, exit it and reset button
        if self._edit_mode:
            self._edit_mode = False
            self.builder_grp.blockSignals(True)
            self.builder_grp.setChecked(False)
            self.builder_grp.blockSignals(False)
            self._reset_builder()

        # Reset button to "Edit" state since metadata is now loaded
        self.btn_save_metadata.setText("Edit Metadata CSV")
        self.btn_save_metadata.setStyleSheet(
            "background-color: #FFA726; color: white; font-weight: bold;"
        )

        # Reload + emit so other tabs update
        self.csv_path_edit.setText(str(csv_path))
        self.metadata_loaded.emit(self.metadata)

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
            self._loaded_csv_path = csv_path

            # Switch save button to "Edit Metadata CSV" (yellow)
            self.btn_save_metadata.setText("Edit Metadata CSV")
            self.btn_save_metadata.setStyleSheet(
                "background-color: #FFA726; color: white; font-weight: bold;"
            )

            # Populate overview and dim-order table
            self._populate_metadata_overview()
            self._populate_dim_order_table()

            # Emit signal for other tabs
            self.metadata_loaded.emit(self.metadata)

        except Exception as e:
            self.metadata_info_label.setText(f"❌ {e}")
            self._log(f"❌ Error loading metadata: {e}")

    # ══════════════════════════════════════════════════════════════════════
    # Section 4 – Metadata Overview
    # ══════════════════════════════════════════════════════════════════════
    def _build_metadata_overview_section(self):
        grp = QGroupBox("4 · Metadata Overview")
        lay = QVBoxLayout(grp)

        self.metadata_table = QTableWidget(0, 0)
        self.metadata_table.setAlternatingRowColors(True)
        # Fix for white-on-white text in some themes
        self.metadata_table.setStyleSheet("QTableWidget { alternate-background-color: #333; }")
        self.metadata_table.setMaximumHeight(250)
        lay.addWidget(self.metadata_table)

        self._layout.addWidget(grp)

    def _populate_metadata_overview(self):
        if self.metadata is None or self.metadata.empty:
            self.metadata_table.setRowCount(0)
            self.metadata_table.setColumnCount(0)
            return

        # Show first 10 rows for overview
        df_preview = self.metadata.head(10)
        self.metadata_table.setRowCount(len(df_preview))
        self.metadata_table.setColumnCount(len(df_preview.columns))
        self.metadata_table.setHorizontalHeaderLabels(df_preview.columns)

        for i in range(len(df_preview)):
            for j in range(len(df_preview.columns)):
                val = df_preview.iloc[i, j]
                item = QTableWidgetItem(str(val) if pd.notna(val) else "")
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.metadata_table.setItem(i, j, item)

        # self.metadata_table.resizeColumnsToContents()  # Removed to prevent excessive width

    # ══════════════════════════════════════════════════════════════════════
    # Section 5 – Dimension Order Table
    # ══════════════════════════════════════════════════════════════════════
    def _build_dim_order_section(self):
        grp = QGroupBox("5 · Dimension Order")
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

        # Always apply UI choices back to the in-memory DataFrame first,
        # regardless of whether a CSV path is available.
        for i in range(self.dim_table.rowCount()):
            sample = self.dim_table.item(i, 0).text()
            combo = self.dim_table.cellWidget(i, 2)
            if combo:
                self.metadata.loc[
                    self.metadata["sample_name"] == sample, "dimension_order"
                ] = combo.currentText()

        # Persist to CSV only when a path is known.
        csv_path = self.csv_path_edit.text().strip()
        if csv_path:
            self.metadata.to_csv(csv_path, index=False)
            self._log(f"✅ Dimension orders updated in memory and saved to {csv_path}")
        else:
            self._log("✅ Dimension orders updated in memory (no CSV path set — save metadata first to persist)")

        # Refresh the metadata overview table so the user sees the changes.
        self._populate_metadata_overview()

    # ══════════════════════════════════════════════════════════════════════
    # Section 6 – Convert to Zarr
    # ══════════════════════════════════════════════════════════════════════
    def _build_zarr_section(self):
        grp = QGroupBox("6 · Convert to Zarr")
        lay = QVBoxLayout(grp)

        # Timepoint clipping
        self.zarr_clip_check = QCheckBox("Clip timepoints (select range)")
        self.zarr_clip_check.setChecked(False)
        lay.addWidget(self.zarr_clip_check)

        clip_row = QHBoxLayout()
        clip_row.setContentsMargins(20, 0, 0, 0)
        clip_row.setSpacing(4)
        clip_row.addWidget(QLabel("Start:"))
        self.zarr_t_start = QSpinBox()
        self.zarr_t_start.setRange(0, 99999)
        self.zarr_t_start.setValue(0)
        self.zarr_t_start.setMaximumWidth(80)
        self.zarr_t_start.setEnabled(False)
        clip_row.addWidget(self.zarr_t_start)
        clip_row.addWidget(QLabel("End:"))
        self.zarr_t_end = QSpinBox()
        self.zarr_t_end.setRange(0, 99999)
        self.zarr_t_end.setValue(0)
        self.zarr_t_end.setMaximumWidth(80)
        self.zarr_t_end.setEnabled(False)
        clip_row.addWidget(self.zarr_t_end)
        clip_row.addStretch()
        lay.addLayout(clip_row)

        def _toggle_clip(state):
            enabled = self.zarr_clip_check.isChecked()
            self.zarr_t_start.setEnabled(enabled)
            self.zarr_t_end.setEnabled(enabled)
        self.zarr_clip_check.stateChanged.connect(_toggle_clip)

        # Workers
        n_cores = os.cpu_count() or 2
        default_workers = max(1, n_cores // 2)
        workers_row = QHBoxLayout()
        workers_row.addWidget(QLabel("Workers:"))
        self.zarr_workers_spin = QSpinBox()
        self.zarr_workers_spin.setRange(1, max(1, n_cores - 1))
        self.zarr_workers_spin.setValue(default_workers)
        self.zarr_workers_spin.setMaximumWidth(70)
        self.zarr_workers_spin.setToolTip(
            "Number of parallel worker processes for zarr conversion.\n"
            f"Defaults to half of available CPU cores ({default_workers}/{n_cores})."
        )
        workers_row.addWidget(self.zarr_workers_spin)
        workers_row.addStretch()
        lay.addLayout(workers_row)

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

        # Gather clipping params
        t_start = None
        t_end = None
        if self.zarr_clip_check.isChecked():
            t_start = int(self.zarr_t_start.value())
            t_end = int(self.zarr_t_end.value())
            if t_end <= t_start:
                QMessageBox.warning(self, "Error",
                    f"End timepoint ({t_end}) must be greater than start ({t_start}).")
                return

        self.zarr_btn.setEnabled(False)
        self.zarr_status.setText("⏳ Converting…")
        self._log("Starting zarr conversion…")

        self._zarr_worker = _ZarrWorker(
            out_dir, self.metadata.copy(),
            t_start=t_start, t_end=t_end,
            n_workers=self.zarr_workers_spin.value(), parent=self
        )
        self._zarr_worker.progress.connect(lambda msg: self._log(msg))
        self._zarr_worker.finished.connect(self._on_zarr_done)
        self._zarr_worker.start()

    def _on_zarr_done(self, success: bool, message: str,
                       updated_metadata, originals_map):
        self.zarr_btn.setEnabled(True)
        self.zarr_status.setText(message.split("\n")[0])
        self._log(message)

        if not success or updated_metadata is None:
            return

        # 1) Persist new zarr paths into in-memory metadata
        self.metadata = updated_metadata
        self._log("Updated in-memory metadata with new .zarr paths.")

        # 2) Persist to CSV if a path is known
        csv_path = ""
        if hasattr(self, "csv_path_edit"):
            csv_path = self.csv_path_edit.text().strip()
        if not csv_path and self._loaded_csv_path:
            csv_path = self._loaded_csv_path
        if csv_path:
            try:
                self.metadata.to_csv(csv_path, index=False)
                self._loaded_csv_path = csv_path
                self._log(f"✅ Metadata CSV updated: {csv_path}")
            except Exception as e:
                self._log(f"⚠ Failed to write metadata CSV: {e}")
        else:
            self._log("⚠ No CSV path set — metadata changes are in-memory only.")

        # 3) Refresh metadata overview so the table reflects the new paths
        try:
            self._populate_metadata_overview()
        except Exception as e:
            self._log(f"⚠ Could not refresh metadata overview: {e}")

        # 4) Offer deletion of the now-redundant original files
        self._offer_delete_originals(originals_map or [])

    # ══════════════════════════════════════════════════════════════════════
    # Section 6b – Optional cleanup of original raw files post-conversion
    # ══════════════════════════════════════════════════════════════════════
    @staticmethod
    def _path_size_bytes(p: Path) -> int:
        """Return total size in bytes for a file or a directory tree."""
        try:
            if p.is_file():
                return p.stat().st_size
            if p.is_dir():
                total = 0
                for f in p.rglob("*"):
                    try:
                        if f.is_file():
                            total += f.stat().st_size
                    except OSError:
                        continue
                return total
        except OSError:
            return 0
        return 0

    @staticmethod
    def _human_bytes(n: int) -> str:
        """Format a byte count as a human-readable string."""
        if n is None:
            return "0 B"
        n = float(n)
        if n >= 1024 ** 3:
            return f"{n / (1024 ** 3):.2f} GB"
        if n >= 1024 ** 2:
            return f"{n / (1024 ** 2):.1f} MB"
        if n >= 1024:
            return f"{n / 1024:.1f} KB"
        return f"{int(n)} B"

    def _collect_deletion_candidates(self, originals_map: list) -> list:
        """Filter originals_map down to entries safe to delete.

        Skips:
          * entries where the original is already a .zarr (those are not converted)
          * entries where original and new zarr resolve to the same path
          * entries where the new zarr does not exist or fails to open
          * entries where the original no longer exists on disk
        """
        candidates: list = []
        for entry in originals_map:
            sample = entry.get("sample", "")
            original = entry.get("original", "")
            new_zarr = entry.get("new_zarr", "")
            if not original or not new_zarr:
                continue

            orig_path = Path(original)
            zarr_path = Path(new_zarr)

            # Never offer to delete a .zarr — those are not reconverted.
            if orig_path.suffix == ".zarr" or str(orig_path).endswith(".zarr.zip"):
                continue

            # Same-path safety: do not delete if original IS the new zarr.
            try:
                if orig_path.resolve() == zarr_path.resolve():
                    continue
            except OSError:
                continue

            # Verify the original still exists.
            if not orig_path.exists():
                self._log(f"⚠ Original for '{sample}' not found, skipping: {orig_path}")
                continue

            # Verify the new zarr exists AND can be opened.
            if not zarr_path.exists():
                self._log(
                    f"⚠ New zarr for '{sample}' not found, keeping original: {zarr_path}"
                )
                continue
            try:
                from behav3d.io.formats.zarr import load_zarr
                _ = load_zarr(zarr_path)
            except Exception as e:
                self._log(
                    f"⚠ New zarr for '{sample}' could not be opened, keeping original: {e}"
                )
                continue

            size_bytes = self._path_size_bytes(orig_path)
            candidates.append({
                "sample": sample,
                "path": orig_path,
                "size_bytes": size_bytes,
            })
        return candidates

    def _offer_delete_originals(self, originals_map: list) -> None:
        """Run the two-step pop-up flow to offer deletion of original raw files."""
        candidates = self._collect_deletion_candidates(originals_map)
        if not candidates:
            self._log("Nothing to clean up: no deletable originals found.")
            return

        total_bytes = sum(c["size_bytes"] for c in candidates)
        human_total = self._human_bytes(total_bytes)
        n = len(candidates)

        # ───── Pop-up #1: informational offer ─────
        info = QMessageBox(self)
        info.setIcon(QMessageBox.Information)
        info.setWindowTitle("Free up disk space?")
        info.setTextFormat(Qt.RichText)
        info.setText(
            f"<b>Zarr conversion completed successfully.</b><br><br>"
            f"All {n} image{'s' if n != 1 else ''} have been duplicated as Zarr "
            f"and the metadata now points to the new <code>.zarr</code> files.<br><br>"
            f"The original raw files are still on disk and are now redundant. "
            f"Removing them would free approximately "
            f"<b>~ {human_total}</b>."
        )
        info.setInformativeText(
            "Would you like to delete the original files?\n"
            "(You can also keep them and remove them manually later.)"
        )
        info.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        info.setDefaultButton(QMessageBox.No)
        choice = info.exec_()
        if choice != QMessageBox.Yes:
            self._log("Originals kept on disk (user declined cleanup).")
            return

        # ───── Pop-up #2: final irreversible warning ─────
        warn = QMessageBox(self)
        warn.setIcon(QMessageBox.Warning)
        warn.setWindowTitle("Permanently delete original files?")
        warn.setTextFormat(Qt.RichText)
        warn.setText(
            f"You are about to permanently delete <b>{n}</b> original image "
            f"file{'s' if n != 1 else ''}, freeing approximately "
            f"<b>~ {human_total}</b>."
        )
        warn.setInformativeText(
            "⚠  This action CANNOT be undone. Files will NOT go to the Recycle Bin.\n"
            "Make sure the new .zarr files load correctly before continuing."
        )
        warn.setDetailedText(
            "Files to be removed:\n"
            + "\n".join(
                f"  • [{c['sample']}] {c['path']}   ({self._human_bytes(c['size_bytes'])})"
                for c in candidates
            )
        )
        cancel_btn = warn.addButton("Cancel", QMessageBox.RejectRole)
        delete_btn = warn.addButton("Delete permanently", QMessageBox.DestructiveRole)
        warn.setDefaultButton(cancel_btn)
        delete_btn.setStyleSheet(
            "background-color: #c62828; color: white; font-weight: bold;"
        )
        warn.exec_()
        if warn.clickedButton() is not delete_btn:
            self._log("Originals kept on disk (deletion cancelled).")
            return

        # ───── Perform deletion ─────
        self._delete_originals(candidates)

    def _delete_originals(self, candidates: list) -> None:
        """Delete each candidate's original file or directory and report results."""
        freed = 0
        deleted: list = []
        errors: list = []
        for c in candidates:
            p: Path = c["path"]
            try:
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
                freed += int(c["size_bytes"])
                deleted.append(c)
                self._log(f"🗑 Deleted: {p}")
            except Exception as e:
                errors.append((p, str(e)))
                self._log(f"⚠ Failed to delete {p}: {e}")

        n_ok = len(deleted)
        n_err = len(errors)
        human_freed = self._human_bytes(freed)

        summary = QMessageBox(self)
        summary.setTextFormat(Qt.RichText)
        if n_err == 0:
            summary.setIcon(QMessageBox.Information)
            summary.setWindowTitle("Originals removed")
            summary.setText(
                f"Successfully deleted <b>{n_ok}</b> file"
                f"{'s' if n_ok != 1 else ''}.<br>"
                f"Disk space freed: <b>{human_freed}</b>."
            )
        else:
            summary.setIcon(QMessageBox.Warning)
            summary.setWindowTitle("Cleanup completed with errors")
            summary.setText(
                f"Deleted <b>{n_ok}</b> file{'s' if n_ok != 1 else ''} "
                f"(<b>{human_freed}</b> freed).<br>"
                f"<b>{n_err}</b> file{'s' if n_err != 1 else ''} could not be removed."
            )
            summary.setDetailedText(
                "\n".join(f"{p}: {msg}" for p, msg in errors)
            )
        summary.setStandardButtons(QMessageBox.Ok)
        summary.exec_()

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
