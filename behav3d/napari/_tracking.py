"""
BEHAV3D napari plugin – Tracking Tab.

Provides per-cell-type sub-tabs with method selection
(LapTrack / TrackPy / Fragmentation Propagation),
method-specific parameters, and batch-tracking options.

Run buttons execute backend tracking in a background ``QThread`` (via
:class:`behav3d.napari._background_runner.BackgroundOperation`) so napari
stays responsive.  Progress is forwarded to the per-tab
:class:`ProgressBarRow` and the napari activity dock.  The queue keeps
its existing synchronous behaviour by calling
``run_batch_tracking(block=True)``.
"""
import os
import sys
import traceback
from functools import partial
from pathlib import Path

import pandas as pd
import yaml
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QTabWidget, QGroupBox, QLabel, QPushButton,
    QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox,
    QScrollArea, QTextEdit, QStackedWidget,
    QLineEdit, QFileDialog, QMessageBox,
)
from qtpy.QtCore import Qt, Signal

from behav3d.napari._widgets import (
    make_help_row,
    HelpButton,
    HelpSection,
    IllustratedHelpButton,
    browse_file_or_zarr,
    prompt_axis_order,
    resolve_external_path,
)
from behav3d.napari._units import UnitGroupManager
from behav3d.core.qt_help import reset_scroll_on_page_change
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
    fire_extra_callback,
)


# Canonical order + display labels for the tracking-method selector. The
# internal keys (left column) are load-bearing — they are what gets persisted
# to behav3d_parameters.yml and dispatched to the backend — while the labels
# and their order here are purely presentational.
_METHOD_ORDER = [
    "btrack", "propagation", "bounded_propagation",
    "reporter_propagation", "trackpy", "lap", "import",
]
_METHOD_LABELS = {
    "btrack": "btrack (Bayesian)",
    "propagation": "Fragmentation Propagation",
    "bounded_propagation": "Bounded Propagation",
    "reporter_propagation": "Reporter Propagation",
    "trackpy": "TrackPy",
    "lap": "LapTrack",
    "import": "Import tracking",
}
_METHOD_IDX = {k: i for i, k in enumerate(_METHOD_ORDER)}


# ═══════════════════════════════════════════════════════════════════════════
# _ImportTrackingPage — Page 4 of CellTypeTrackingPanel
# ═══════════════════════════════════════════════════════════════════════════
class _ImportTrackingPage(QWidget):
    """Per-cell-type widget for importing pre-tracked zarr/tiff files.

    Shows an editable per-sample path row (prefilled from the metadata
    column ``{prefix}_{cell_type}_tracks_image_path`` when already set),
    validates/converts it, and writes outputs to the standard BEHAV3D
    output locations — exactly the same paths used by every other tracking
    algorithm. Newly browsed/typed paths are staged in the row's own widget;
    metadata.csv is only updated once Convert/Import/Re-process actually runs
    (or a batch "Process All").
    """

    def __init__(self, cell_type: str, category: str, metadata_loader, parent=None,
                 switch_to_data_prep_edit_callback=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category
        self.metadata_loader = metadata_loader
        self._prefix = {"organoid": "or", "immune": "im"}.get(category, "ot")
        self._switch_to_data_prep_edit = switch_to_data_prep_edit_callback
        # sample_name -> {"path_edit", "browse_btn", "status_layout", "last_value", "row_idx"}
        self._rows = {}
        self._init_ui()
        if hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)
        if getattr(metadata_loader, "metadata", None) is not None:
            self._rebuild()

    # ── column helpers ──────────────────────────────────────────────────
    def _tracks_img_col(self):
        return f"{self._prefix}_{self.cell_type}_tracks_image_path"

    def _tracks_csv_col(self):
        return f"{self._prefix}_{self.cell_type}_tracks_csv_path"

    # ── standard output paths (same as all tracking algorithms) ────────
    def _dest_zarr(self, sample_name: str) -> Path:
        return (Path(self.metadata_loader.output_dir)
                / "images" / sample_name
                / f"{sample_name}_{self.cell_type}_tracked.zarr")

    def _dest_csv(self, sample_name: str) -> Path:
        return (Path(self.metadata_loader.output_dir)
                / "trackdata" / sample_name / self.cell_type
                / f"{sample_name}_{self.cell_type}_tracks.csv")

    # ── path helpers ────────────────────────────────────────────────────
    def _resolve_path(self, path_str: str):
        md_csv = (self.metadata_loader.behav3d_parameters
                  .get("paths", {}).get("metadata_csv"))
        return resolve_external_path(path_str, md_csv)

    @staticmethod
    def _clear_layout(layout):
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

    @staticmethod
    def _check_zarr_structure(path) -> tuple:
        import zarr
        try:
            root = zarr.open(str(path), mode="r")
            if isinstance(root, zarr.Group):
                return False, "Zarr contains sub-groups instead of a root array"
            if root.chunks[0] != 1:
                return False, f"First chunk dim is {root.chunks[0]}, expected 1"
            return True, "OK"
        except Exception as exc:
            return False, str(exc)

    # ── UI ──────────────────────────────────────────────────────────────
    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

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

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # Initial placeholder content
        content = QWidget()
        self.scroll_layout = QVBoxLayout(content)
        self.scroll_layout.setAlignment(Qt.AlignTop)
        self._scroll.setWidget(content)
        layout.addWidget(self._scroll)

    def _on_metadata_updated(self, _=None):
        self._rebuild()

    def _rebuild(self):
        # Detach the old content widget without destroying it — this prevents
        # Qt from synchronously deleting the clicked button while its signal
        # emission is still on the C++ stack.
        old_content = self._scroll.takeWidget()
        self._rows = {}

        content = QWidget()
        scroll_layout = QVBoxLayout(content)
        scroll_layout.setAlignment(Qt.AlignTop)
        self.scroll_layout = scroll_layout

        md = self.metadata_loader.metadata
        if md is None or md.empty:
            lbl = QLabel("No metadata loaded.")
            lbl.setStyleSheet("color:#888; font-style:italic; padding:20px;")
            scroll_layout.addWidget(lbl)
            self._scroll.setWidget(content)
            if old_content is not None:
                old_content.deleteLater()
            return

        info = QLabel(
            f"<b>Import pre-tracked files for <code>{self.cell_type}</code></b>"
        )
        info.setStyleSheet("padding:6px 4px 2px 4px;")
        scroll_layout.addWidget(info)
        desc = QLabel(
            f"Type or browse a pre-tracked source (.zarr or .tif/.tiff) below "
            f"per sample. Column: {self._tracks_img_col()} — only written to "
            f"metadata.csv once you run Convert/Import for that row."
        )
        desc.setWordWrap(True)
        desc.setStyleSheet("color:#888; font-size:10px; padding:0 4px 10px 4px;")
        scroll_layout.addWidget(desc)

        any_action = False
        for idx, row in md.iterrows():
            sample_name = str(row.get("sample_name", f"Row {idx + 1}"))
            if self._add_sample_row(sample_name, int(idx), row):
                any_action = True

        if any_action:
            btn_all = QPushButton("⚡  Process All Samples")
            btn_all.setStyleSheet(
                "QPushButton{background:#1565C0;color:white;padding:7px 16px;"
                "border-radius:4px;font-weight:bold;font-size:12px}"
                "QPushButton:hover{background:#1976D2}"
            )
            btn_all.clicked.connect(self._process_all)
            scroll_layout.addWidget(btn_all)

        scroll_layout.addStretch()
        self._scroll.setWidget(content)
        if old_content is not None:
            old_content.deleteLater()

    def _add_sample_row(self, sample_name: str, row_idx: int, row) -> bool:
        """Build+register one editable tracking-source row for one sample.

        Returns True if a genuinely new conversion is available (matches the
        old "needs_action" semantics — "Re-process" doesn't count).
        """
        header = QLabel(f"📁  {sample_name}")
        header.setStyleSheet("font-weight:bold; font-size:12px; padding:6px 0 2px 0;")
        self.scroll_layout.addWidget(header)

        col = self._tracks_img_col()
        raw_val = row.get(col) if col in row.index else None
        has_value = (
            raw_val is not None
            and pd.notna(raw_val)
            and str(raw_val).strip() not in ("", "nan")
        )
        initial_value = str(raw_val).strip().strip('"').strip("'") if has_value else ""

        row_w = QWidget()
        row_lay = QHBoxLayout(row_w)
        row_lay.setContentsMargins(16, 2, 4, 4)

        path_edit = QLineEdit(initial_value)
        path_edit.setPlaceholderText("Path to a pre-tracked .tif/.tiff file or .zarr directory")
        path_edit.setMinimumWidth(220)
        row_lay.addWidget(path_edit, stretch=1)

        browse_btn = QPushButton("Browse…")
        row_lay.addWidget(browse_btn)

        status_container = QWidget()
        status_layout = QHBoxLayout(status_container)
        status_layout.setContentsMargins(8, 0, 0, 0)
        row_lay.addWidget(status_container)
        row_lay.addStretch()

        self._rows[sample_name] = {
            "path_edit": path_edit,
            "browse_btn": browse_btn,
            "status_layout": status_layout,
            "last_value": initial_value,
            "row_idx": row_idx,
        }
        path_edit.editingFinished.connect(partial(self._on_row_path_edited, sample_name))
        browse_btn.clicked.connect(partial(self._on_browse_clicked, sample_name))

        self.scroll_layout.addWidget(row_w)
        self._refresh_row_status(sample_name)

        needs_action = False
        for i in range(status_layout.count()):
            w = status_layout.itemAt(i).widget()
            if isinstance(w, QPushButton) and "Re-process" not in w.text():
                needs_action = True
                break

        sep = QWidget()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background:#ddd; margin:2px 0;")
        self.scroll_layout.addWidget(sep)
        return needs_action

    # ── live status refresh (reads from the widget, not metadata) ───────
    def _refresh_row_status(self, sample_name: str):
        info = self._rows.get(sample_name)
        if info is None:
            return
        self._clear_layout(info["status_layout"])

        path_str = info["path_edit"].text().strip().strip('"').strip("'")
        if not path_str:
            lbl = QLabel("No source path set")
            lbl.setStyleSheet("color:#999; font-style:italic;")
            info["status_layout"].addWidget(lbl)
            return

        file_path = self._resolve_path(path_str)
        if file_path is None or not file_path.exists():
            lbl = QLabel("⚠️  File not found")
            lbl.setToolTip(str(file_path))
            lbl.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(lbl)
            return

        row_idx = info["row_idx"]
        dest_z = self._dest_zarr(sample_name)
        dest_c = self._dest_csv(sample_name)

        if file_path.suffix.lower() in (".tif", ".tiff"):
            already = dest_z.exists() and dest_c.exists()
            btn_text = "🔄  Re-process" if already else "🔄  Convert TIFF → zarr"
            btn_style = (
                "QPushButton{background:#546E7A;color:white;padding:4px 8px;"
                "border-radius:3px;font-size:10px}"
                "QPushButton:hover{background:#607D8B}"
                if already else
                "QPushButton{background:#1565C0;color:white;padding:4px 10px;border-radius:3px}"
                "QPushButton:hover{background:#1976D2}"
            )
            if already:
                lbl = QLabel("✅  Already processed")
                lbl.setStyleSheet("color:#2E7D32; font-weight:bold;")
                info["status_layout"].addWidget(lbl)
            btn = QPushButton(btn_text)
            btn.setToolTip(f"Will write:\n  {dest_z}\n  {dest_c}")
            btn.setStyleSheet(btn_style)
            btn.clicked.connect(partial(self._process_single, sample_name, row_idx))
            info["status_layout"].addWidget(btn)

        elif file_path.suffix == ".zarr" or file_path.is_dir():
            already = dest_z.exists() and dest_c.exists()
            btn_text = "🔄  Re-process" if already else "📄  Import zarr"
            btn_style = (
                "QPushButton{background:#546E7A;color:white;padding:4px 8px;"
                "border-radius:3px;font-size:10px}"
                "QPushButton:hover{background:#607D8B}"
                if already else
                "QPushButton{background:#2E7D32;color:white;padding:4px 10px;border-radius:3px}"
                "QPushButton:hover{background:#388E3C}"
            )
            if already:
                lbl = QLabel("✅  Already processed")
                lbl.setStyleSheet("color:#2E7D32; font-weight:bold;")
                info["status_layout"].addWidget(lbl)
            btn = QPushButton(btn_text)
            btn.setToolTip(f"Will write:\n  {dest_z}\n  {dest_c}")
            btn.setStyleSheet(btn_style)
            btn.clicked.connect(partial(self._process_single, sample_name, row_idx))
            info["status_layout"].addWidget(btn)
        else:
            lbl = QLabel(f"⚠️  Unsupported format ({file_path.suffix})")
            lbl.setStyleSheet("color:#E65100;")
            info["status_layout"].addWidget(lbl)

    # ── path-field change handlers ───────────────────────────────────────
    def _on_browse_clicked(self, sample_name):
        info = self._rows.get(sample_name)
        if info is None:
            return
        new_path = browse_file_or_zarr(
            self, f"Select {self.cell_type} tracking source for {sample_name}",
            "Image files (*.tif *.tiff *.zarr);; All Files (*)",
            allow_zarr=True,
        )
        if not new_path:
            return
        self._maybe_accept_new_value(sample_name, new_path)

    def _on_row_path_edited(self, sample_name):
        info = self._rows.get(sample_name)
        if info is None:
            return
        new_value = info["path_edit"].text().strip()
        self._maybe_accept_new_value(sample_name, new_value, already_in_field=True)

    def _maybe_accept_new_value(self, sample_name, new_value, already_in_field=False):
        info = self._rows.get(sample_name)
        if info is None:
            return
        new_value = str(new_value).strip()
        old_value = info["last_value"]

        if new_value != old_value and old_value:
            res = QMessageBox.question(
                self, "Replace existing path?",
                f"This will replace the existing tracking source path for "
                f"{self.cell_type} / {sample_name}:\n{old_value}\n→\n{new_value}\n\nContinue?",
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

        self._refresh_row_status(sample_name)

    # ── processing ──────────────────────────────────────────────────────
    def _process_single(self, sample_name: str, row_idx: int, _=None, save: bool = True):
        import shutil
        from behav3d.io.images import load_image, save_as_zarr, convert_label_file_to_zarr
        from behav3d.preprocessing.tracking import convert_tracked_image_to_csv

        info = self._rows.get(sample_name)
        if info is None:
            return
        src_path_str = info["path_edit"].text().strip().strip('"').strip("'")
        src = self._resolve_path(src_path_str)
        if src is None or not src.exists():
            return

        dest_z = self._dest_zarr(sample_name)
        dest_c = self._dest_csv(sample_name)
        dest_z.parent.mkdir(parents=True, exist_ok=True)
        dest_c.parent.mkdir(parents=True, exist_ok=True)

        # Warn only if the destination zarr already exists
        if dest_z.exists():
            res = QMessageBox.question(
                self,
                "Overwrite existing tracking?",
                f"A tracked zarr already exists at:\n{dest_z}\n\nOverwrite it?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                return

        try:
            # Step 1 — ensure valid zarr at dest_z
            zarr_ok = False
            if src.suffix == ".zarr" or src.is_dir():
                zarr_ok, _ = self._check_zarr_structure(src)

            if src.suffix.lower() in (".tif", ".tiff"):
                axis_order = prompt_axis_order(self, src, log=print)
                if axis_order is None:
                    print(f"  Cancelled: no axis order selected for {sample_name}", flush=True)
                    return
                print(f"  Converting {src.name} → {dest_z.name} …", flush=True)
                if dest_z.exists():
                    shutil.rmtree(dest_z)
                convert_label_file_to_zarr(path=src, outpath=dest_z, axis_order=axis_order, overwrite=True)
            elif not zarr_ok:
                print(f"  Repairing zarr format: {src.name} → {dest_z.name} …", flush=True)
                if dest_z.exists():
                    shutil.rmtree(dest_z)
                save_as_zarr(load_image(src), dest_z)
            else:
                # Valid zarr — copy to standard output location if needed
                if src.resolve() != dest_z.resolve():
                    print(f"  Copying {src.name} → {dest_z.name} …", flush=True)
                    if dest_z.exists():
                        shutil.rmtree(dest_z)
                    shutil.copytree(str(src), str(dest_z))

            zarr_path = dest_z

            # Step 2 — generate tracks CSV (unchanged: already-correct logic)
            md_row = self.metadata_loader.metadata.iloc[row_idx]
            element_size_x = float(md_row.get("pixel_distance_xy") or 1)
            element_size_y = float(md_row.get("pixel_distance_xy") or 1)
            element_size_z = float(md_row.get("pixel_distance_z") or 1)

            print(f"  Generating tracks CSV from {zarr_path.name} …", flush=True)
            convert_tracked_image_to_csv(
                img_path=zarr_path,
                outpath=dest_c,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
            )
            print(f"  ✅ {dest_c.name}", flush=True)

            # Step 3 — update metadata (same pattern as _run_tracking_for)
            img_col = self._tracks_img_col()
            csv_col = self._tracks_csv_col()
            md = self.metadata_loader.metadata.copy()
            if img_col not in md.columns:
                md[img_col] = pd.Series(dtype=object)
            if csv_col not in md.columns:
                md[csv_col] = pd.Series(dtype=object)
            md.at[row_idx, img_col] = str(zarr_path)
            md.at[row_idx, csv_col] = str(dest_c)
            self.metadata_loader.metadata = md

            info["last_value"] = str(zarr_path)
            info["path_edit"].blockSignals(True)
            info["path_edit"].setText(str(zarr_path))
            info["path_edit"].blockSignals(False)

            if save:
                self._save_metadata(md)

        except Exception as exc:
            traceback.print_exc()
            print(f"  ❌ Error processing {sample_name}: {exc}", flush=True)

    def _process_all(self, _=None):
        md = self.metadata_loader.metadata
        if md is None:
            return
        for idx, row in md.iterrows():
            sample_name = str(row.get("sample_name", f"Row {idx + 1}"))
            info = self._rows.get(sample_name)
            if info is None:
                continue
            path_str = info["path_edit"].text().strip().strip('"').strip("'")
            if not path_str:
                continue
            fp = self._resolve_path(path_str)
            if fp is None or not fp.exists():
                continue
            self._process_single(sample_name, int(idx), save=False)
        self._save_metadata(self.metadata_loader.metadata)

    def _save_metadata(self, md):
        csv_path = (self.metadata_loader.behav3d_parameters
                    .get("paths", {}).get("metadata_csv"))
        if csv_path:
            md.to_csv(csv_path, sep=",", index=False)
            print(f"  Metadata saved → {csv_path}", flush=True)
        # Defer _rebuild() to the next event-loop cycle so that any button
        # click handler that triggered this call has fully returned before Qt
        # destroys the button widget (avoids 0xC0000005 access violation).
        from qtpy.QtCore import QTimer
        QTimer.singleShot(0, self._rebuild)


def _cfg_get(cfg: dict, dotted_key: str, default=None):
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


# ═══════════════════════════════════════════════════════════════════════════
# CellTypeTrackingPanel — one per cell type
# ═══════════════════════════════════════════════════════════════════════════
class CellTypeTrackingPanel(QWidget):
    """Parameter panel for tracking a single cell type."""

    def __init__(self, cell_type: str, category: str, metadata_loader,
                 all_cell_types: list, category_types: list,
                 log_callback=None, viewer=None, parent=None,
                 tab_progress_row=None, on_tracking_complete_callback=None,
                 switch_to_data_prep_edit_callback=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category          # "organoid" / "immune" / "other"
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types   # all types in same category
        self.log = log_callback or print
        self.viewer = viewer
        self.on_tracking_complete = on_tracking_complete_callback  # callback when tracking finishes
        self._toggle_all_organoids_callback = None  # set by TrackingTab for organoid panels
        self._switch_to_data_prep_edit = switch_to_data_prep_edit_callback

        # Background-execution infrastructure (mirrors editing-mode pattern).
        # The tab supplies a shared ProgressBarRow so all panels feed the
        # same visual widget, but each panel owns its own BackgroundOperation
        # so an in-progress run is tracked per-panel.
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)

        # Determine defaults
        if category == "organoid":
            def_method = "propagation"
        elif category == "immune":
            def_method = "btrack"
        else:
            def_method = "lap"

        # Load saved config
        tcfg = _cfg_get(self.metadata_loader.behav3d_parameters,
                        f"tracking.{cell_type}", {})
        if not isinstance(tcfg, dict):
            tcfg = {}

        self._build_ui(tcfg, def_method)

    # ------------------------------------------------------------------
    def _build_ui(self, tcfg, def_method):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # ── Method selector ──────────────────────────────────────────
        method_group = QGroupBox("Tracking Method")
        method_layout = QHBoxLayout()
        method_layout.setContentsMargins(6, 4, 6, 4)
        self.combo_method = QComboBox()
        self.combo_method.addItems([_METHOD_LABELS[k] for k in _METHOD_ORDER])
        saved_method = tcfg.get("method", def_method)
        self.combo_method.setCurrentIndex(_METHOD_IDX.get(saved_method, 0))
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.combo_method)
        # The three propagation variants differ in ways a paragraph describes
        # poorly, so each carries a schematic instead (behav3d/resources/tracking/).
        method_layout.addWidget(IllustratedHelpButton(
            "Tracking Method",
            [
                HelpSection(
                    "btrack (Bayesian)",
                    "Kalman-filter based tracker with an optional global "
                    "hypothesis optimizer for resolving merges, splits and "
                    "false positives.",
                ),
                HelpSection(
                    "Fragmentation Propagation",
                    "No tunable parameters. Objects are identified by spatial "
                    "overlap between frames rather than a linking cost, so one "
                    "identity survives fragmentation and fusion — which is what "
                    "lets death-dye signal be quantified while an organoid "
                    "breaks apart or merges.",
                    "tracking/fragmentation_tracking.png",
                ),
                HelpSection(
                    "Bounded Propagation",
                    "The same overlap propagation, but a track ID can never span "
                    "more than one physically disconnected region: each frame's "
                    "mask is split into its connected regions first and every "
                    "existing track claims whichever region its previous footprint "
                    "overlaps most. Identities that start out merged are recovered "
                    "as distinct tracks once the objects separate.",
                    "tracking/bounded_propagation.png",
                ),
                HelpSection(
                    "Reporter Propagation",
                    "For near-static objects whose segmentation flickers on and "
                    "off with a fluctuating reporter. Every segment detected at "
                    "any timepoint is pooled, spatially overlapping ones are "
                    "grouped regardless of how far apart in time they are, and "
                    "the largest instance of each group is stamped onto every "
                    "timepoint. No linking or gap parameters, since grouping is "
                    "purely spatial.",
                    "tracking/reporter_propagation.png",
                ),
                HelpSection(
                    "TrackPy",
                    "Crocker-Grier style nearest-neighbour linker with an "
                    "adaptive search range; simple and fast, no merge/split "
                    "support.",
                ),
                HelpSection(
                    "LapTrack",
                    "Links detections frame-to-frame by solving a Linear "
                    "Assignment Problem on distance costs. Supports gap "
                    "closing, merging and splitting.",
                ),
                HelpSection(
                    "Import tracking",
                    "Load pre-computed track IDs from an existing tracks file "
                    "instead of running a tracker.",
                ),
            ],
        ))
        method_layout.addStretch()
        method_group.setLayout(method_layout)

        layout.addWidget(method_group)

        # ── Stacked method params ────────────────────────────────────
        self.param_stack = QStackedWidget()
        # Param pages are built below in their own (historical) source order and
        # registered by method key; the combo order is independent of it.
        self._method_pages = {}

        def _show_params_for_combo_idx(idx):
            key = _METHOD_ORDER[idx] if 0 <= idx < len(_METHOD_ORDER) else None
            page = self._method_pages.get(key)
            if page is not None:
                self.param_stack.setCurrentIndex(page)
        self.combo_method.currentIndexChanged.connect(_show_params_for_combo_idx)

        # Page 0 — LAP params
        lap_cfg = tcfg.get("lap", {})
        lap_page = QWidget()
        lap_form = QFormLayout(lap_page)
        lap_form.setContentsMargins(6, 4, 6, 4)
        lap_form.setSpacing(3)
        lap_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        # Per-method physical(µm)/pixel unit toggle. LAP links on
        # micron-scaled coordinates, so its distance costs are natively in
        # µm even though the config keys are historically named "_px".
        self._lap_unit_mgr = UnitGroupManager(
            self.metadata_loader.metadata, default_physical=True
        )
        lap_form.addRow("Distance units:", self._lap_unit_mgr.header_row(label=""))

        self.lap_track_cost = QSpinBox()
        self.lap_track_cost.setRange(1, 999999)
        self.lap_track_cost.setValue(int(lap_cfg.get("track_cost_px", 45)))
        self.lap_track_cost.setMaximumWidth(90)
        lap_form.addRow("Track cost:", make_help_row(
            self.lap_track_cost,
            "Track Cost",
            "Maximum distance a cell can travel between two consecutive "
            "frames to be linked as the same track.\n\n"
            "Entered in the unit selected above (µm by default; LapTrack links "
            "in µm natively).\n\n"
            "Increase if cells move fast; decrease to avoid false links."
        ))
        self._lap_unit_mgr.register(
            self.lap_track_cost, "distance",
            int(lap_cfg.get("track_cost_px", 45)), native_unit="physical",
        )

        self.lap_gap_cost = QSpinBox()
        self.lap_gap_cost.setRange(1, 999999)
        self.lap_gap_cost.setValue(int(lap_cfg.get("gap_close_cost_px", 60)))
        self.lap_gap_cost.setMaximumWidth(90)
        lap_form.addRow("Gap close cost:", make_help_row(
            self.lap_gap_cost,
            "Gap Closing Cost",
            "Maximum distance allowed when reconnecting a track that was "
            "temporarily lost for one or more frames.\n\n"
            "Entered in the unit selected above (µm by default).\n\n"
            "Should be >= Track cost. Increase if cells disappear "
            "briefly due to segmentation gaps."
        ))
        self._lap_unit_mgr.register(
            self.lap_gap_cost, "distance",
            int(lap_cfg.get("gap_close_cost_px", 60)), native_unit="physical",
        )

        self.lap_gap_frames = QSpinBox()
        self.lap_gap_frames.setRange(0, 100)
        self.lap_gap_frames.setValue(int(lap_cfg.get("gap_close_max_frames", 3)))
        self.lap_gap_frames.setMaximumWidth(60)
        lap_form.addRow("Gap close max frames:", make_help_row(
            self.lap_gap_frames,
            "Gap Closing Max Frames",
            "Maximum number of consecutive frames a cell can be missing "
            "before the gap is too large to close.\n\n"
            "Higher values recover longer disappearances but may "
            "introduce false links."
        ))

        self.lap_merge_cost = QSpinBox()
        self.lap_merge_cost.setRange(0, 999999)
        self.lap_merge_cost.setValue(int(lap_cfg.get("merging_cost_px", 0)))
        self.lap_merge_cost.setMaximumWidth(90)
        lap_form.addRow("Merging cost:", make_help_row(
            self.lap_merge_cost,
            "Merging Cost",
            "Maximum distance for detecting merge events, where two "
            "tracks combine into one object.\n\n"
            "Entered in the unit selected above (µm by default).\n\n"
            "Set to 0 to disable merging detection.\n"
            "Useful when cells fuse or cluster together."
        ))
        self._lap_unit_mgr.register(
            self.lap_merge_cost, "distance",
            int(lap_cfg.get("merging_cost_px", 0)), native_unit="physical",
        )

        self.lap_split_cost = QSpinBox()
        self.lap_split_cost.setRange(0, 999999)
        self.lap_split_cost.setValue(int(lap_cfg.get("splitting_cost_px", 0)))
        self.lap_split_cost.setMaximumWidth(90)
        lap_form.addRow("Splitting cost:", make_help_row(
            self.lap_split_cost,
            "Splitting Cost",
            "Maximum distance for detecting split events, where one "
            "object divides into two tracks.\n\n"
            "Entered in the unit selected above (µm by default).\n\n"
            "Set to 0 to disable splitting detection.\n"
            "Useful for cell division or organoid fragmentation."
        ))
        self._lap_unit_mgr.register(
            self.lap_split_cost, "distance",
            int(lap_cfg.get("splitting_cost_px", 0)), native_unit="physical",
        )

        self._method_pages["lap"] = self.param_stack.addWidget(lap_page)

        # Page 1 — TrackPy params
        tp_cfg = tcfg.get("trackpy", {})
        tp_page = QWidget()
        tp_form = QFormLayout(tp_page)
        tp_form.setContentsMargins(6, 4, 6, 4)
        tp_form.setSpacing(3)
        tp_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        # TrackPy links on micron-scaled coordinates too, so its search
        # range is natively in µm despite the historical "_px" config key.
        self._tp_unit_mgr = UnitGroupManager(
            self.metadata_loader.metadata, default_physical=True
        )
        tp_form.addRow("Distance units:", self._tp_unit_mgr.header_row(label=""))

        self.tp_search_range = QSpinBox()
        self.tp_search_range.setRange(1, 999999)
        self.tp_search_range.setValue(int(tp_cfg.get("search_range_px", 31)))
        self.tp_search_range.setMaximumWidth(90)
        tp_form.addRow("Search range:", make_help_row(
            self.tp_search_range,
            "Search Range",
            "Maximum distance to look for a cell in the next frame.\n\n"
            "Entered in the unit selected above (µm by default; TrackPy "
            "links in µm natively).\n\n"
            "Should be large enough to cover the fastest-moving cells."
        ))
        self._tp_unit_mgr.register(
            self.tp_search_range, "distance",
            int(tp_cfg.get("search_range_px", 31)), native_unit="physical",
        )

        self.tp_memory = QSpinBox()
        self.tp_memory.setRange(0, 100)
        self.tp_memory.setValue(int(tp_cfg.get("memory_frames", 2)))
        self.tp_memory.setMaximumWidth(60)
        tp_form.addRow("Memory (frames):", make_help_row(
            self.tp_memory,
            "Memory (frames)",
            "Number of frames a cell can disappear and still be "
            "reconnected to its previous track.\n\n"
            "Similar to 'gap closing' in LapTrack."
        ))

        self.tp_adaptive_stop = QDoubleSpinBox()
        self.tp_adaptive_stop.setRange(0.0, 100.0)
        self.tp_adaptive_stop.setSingleStep(0.5)
        self.tp_adaptive_stop.setValue(float(tp_cfg.get("adaptive_stop", 10.0)))
        self.tp_adaptive_stop.setMaximumWidth(80)
        tp_form.addRow("Adaptive stop:", make_help_row(
            self.tp_adaptive_stop,
            "Adaptive Stop",
            "When a frame has too many nearby candidates to link "
            "unambiguously (an oversized 'subnet'), TrackPy retries by "
            "shrinking the search range (see Adaptive step) until "
            "linking succeeds.\n\n"
            "This value is the lower bound (in the unit selected above) "
            "the search range is allowed to shrink to before giving up "
            "and raising an error.\n\n"
            "Lower = more retries allowed (more robust, slower on dense "
            "data). Higher = gives up sooner."
        ))

        self.tp_adaptive_step = QDoubleSpinBox()
        self.tp_adaptive_step.setRange(0.01, 1.0)
        self.tp_adaptive_step.setSingleStep(0.01)
        self.tp_adaptive_step.setDecimals(3)
        self.tp_adaptive_step.setValue(float(tp_cfg.get("adaptive_step", 0.95)))
        self.tp_adaptive_step.setMaximumWidth(80)
        tp_form.addRow("Adaptive step:", make_help_row(
            self.tp_adaptive_step,
            "Adaptive Step",
            "Multiplier applied to the search range at each iteration "
            "of the adaptive search.\n\n"
            "Values close to 1.0 = slow reduction.\n"
            "Values close to 0.5 = aggressive reduction.\n\n"
            "Default 0.95 works well in most cases."
        ))

        self._method_pages["trackpy"] = self.param_stack.addWidget(tp_page)

        # Page 2 — Propagation params
        prop_page = QWidget()
        prop_lay = QVBoxLayout(prop_page)
        prop_lay.setContentsMargins(6, 6, 6, 6)
        prop_lay.setSpacing(6)

        # "Track all organoids together" checkbox (organoid panels only)
        self.check_all_together_prop = None
        if self.category == "organoid":
            self.check_all_together_prop = QCheckBox("Track all organoids together")
            self.check_all_together_prop.setChecked(False)  # False: we are in individual mode
            self.check_all_together_prop.setToolTip(
                "When checked, all organoid types are tracked simultaneously\n"
                "using Fragmentation Propagation, collapsing all organoid tabs into one.\n"
                "This is the default behaviour when multiple organoid types exist."
            )
            def _on_check_all(checked, _self=self):
                if _self._toggle_all_organoids_callback:
                    _self._toggle_all_organoids_callback(checked)
            self.check_all_together_prop.toggled.connect(_on_check_all)
            prop_lay.addWidget(self.check_all_together_prop)

        prop_notice = QLabel("No tunable parameters for\nFragmentation Propagation.")
        prop_notice.setWordWrap(True)
        prop_notice.setAlignment(Qt.AlignCenter)
        prop_notice.setStyleSheet("color: #666; font-style: italic; padding: 10px;")
        prop_lay.addWidget(prop_notice)
        prop_lay.addStretch()

        self._method_pages["propagation"] = self.param_stack.addWidget(prop_page)

        # Page 3 — Bounded Propagation params
        bp_cfg = tcfg.get("bounded_propagation", {})
        bp_page = QWidget()
        bp_form = QFormLayout(bp_page)
        bp_form.setContentsMargins(6, 4, 6, 4)
        bp_form.setSpacing(3)
        bp_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        bp_info = QLabel(
            "Same as Fragmentation Propagation, but a track ID can never span more than "
            "one disconnected region. Each region of the current frame's "
            "mask is claimed by whichever existing track overlaps it most, "
            "before watershed runs; a region no track claims becomes a new "
            "track."
        )
        bp_info.setWordWrap(True)
        bp_info.setStyleSheet("color: #888; font-size: 10px; padding: 2px 0 6px 0;")
        bp_form.addRow(bp_info)

        self.bp_min_overlap_fraction = QDoubleSpinBox()
        self.bp_min_overlap_fraction.setRange(0.0, 1.0)
        self.bp_min_overlap_fraction.setSingleStep(0.05)
        self.bp_min_overlap_fraction.setDecimals(2)
        self.bp_min_overlap_fraction.setValue(float(bp_cfg.get("min_overlap_fraction", 0.0)))
        self.bp_min_overlap_fraction.setMaximumWidth(90)
        bp_form.addRow("Min overlap fraction:", make_help_row(
            self.bp_min_overlap_fraction,
            "Min Overlap Fraction",
            "Minimum fraction of a current-frame region's area that a "
            "track's previous-frame footprint must fill for that track "
            "to be allowed to claim the region as its own.\n\n"
            "0 = any shared pixel qualifies (default). Raise this to stop "
            "a track from claiming a region it barely touches."
        ))

        self.bp_segment_size_min = QSpinBox()
        self.bp_segment_size_min.setRange(1, 999999)
        self.bp_segment_size_min.setValue(int(bp_cfg.get("segment_size_min", 20)))
        self.bp_segment_size_min.setMaximumWidth(90)
        bp_form.addRow("Min segment size:", make_help_row(
            self.bp_segment_size_min,
            "Min Segment Size",
            "Minimum segment size in voxels, applied after watershed.\n\n"
            "This is NOT the same as a segmentation-stage size filter — "
            "the input mask is normally already filtered for whole-object "
            "size before it reaches tracking. This instead cleans up small "
            "fragments watershed splitting itself can introduce (a merge "
            "region divided between two tracks, or a small leftover sliver "
            "spun off into a new track), so it's usually set much lower "
            "than a segmentation-stage minimum.\n\n"
            "Segments that are 2D/flat (only 1 voxel thick along any axis) "
            "are always removed as well, regardless of this value."
        ))

        self._method_pages["bounded_propagation"] = self.param_stack.addWidget(bp_page)

        # Page 4 — Reporter Propagation params
        rp_cfg = tcfg.get("reporter_propagation", {})
        rp_page = QWidget()
        rp_form = QFormLayout(rp_page)
        rp_form.setContentsMargins(6, 4, 6, 4)
        rp_form.setSpacing(3)
        rp_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        rp_info = QLabel(
            "For near-static objects with unreliable/intermittent "
            "segmentation. Groups segments that spatially overlap each "
            "other across the whole video, then uses the single largest "
            "instance per group as that object's mask for every timepoint."
        )
        rp_info.setWordWrap(True)
        rp_info.setStyleSheet("color: #888; font-size: 10px; padding: 2px 0 6px 0;")
        rp_form.addRow(rp_info)

        self.rp_min_overlap_fraction = QDoubleSpinBox()
        self.rp_min_overlap_fraction.setRange(0.0, 1.0)
        self.rp_min_overlap_fraction.setSingleStep(0.05)
        self.rp_min_overlap_fraction.setDecimals(2)
        self.rp_min_overlap_fraction.setValue(float(rp_cfg.get("min_overlap_fraction", 0.1)))
        self.rp_min_overlap_fraction.setMaximumWidth(90)
        rp_form.addRow("Min overlap fraction:", make_help_row(
            self.rp_min_overlap_fraction,
            "Min Overlap Fraction",
            "Minimum spatial overlap (as a fraction of the smaller "
            "segment's area/volume) required for two segments — from any "
            "two timepoints, however far apart — to be grouped as the same "
            "object.\n\n"
            "0 = any shared pixel counts. Raise this if unrelated nearby "
            "objects are being incorrectly merged together."
        ))

        self.rp_segment_size_min = QSpinBox()
        self.rp_segment_size_min.setRange(1, 999999)
        self.rp_segment_size_min.setValue(int(rp_cfg.get("segment_size_min", 100)))
        self.rp_segment_size_min.setMaximumWidth(90)
        rp_form.addRow("Min segment size:", make_help_row(
            self.rp_segment_size_min,
            "Min Segment Size",
            "Minimum segment size in voxels. Segments smaller than this "
            "are ignored entirely (treated as noise), both as candidates "
            "for the 'largest instance' and when grouping."
        ))

        self._method_pages["reporter_propagation"] = self.param_stack.addWidget(rp_page)

        # Page 5 — btrack (Bayesian tracking)
        bt_cfg = tcfg.get("btrack", {})
        btrack_page = QWidget()
        btrack_lay = QVBoxLayout(btrack_page)
        btrack_lay.setContentsMargins(6, 4, 6, 4)
        btrack_lay.setSpacing(4)

        # Info banner
        _bt_info = QLabel(
            "<b>btrack</b> \u2014 Bayesian multi-object tracker.<br>"
            "<b>Step 1</b> (Kalman filter): always runs \u2014 tune and validate first.<br>"
            "<b>Step 2</b> (optimizer): opt-in \u2014 enable only once Step 1 looks correct.<br>"
            "See <tt>models/README_btrack.md</tt> for full parameter reference."
        )
        _bt_info.setWordWrap(True)
        _bt_info.setStyleSheet("color: #888; font-size: 10px; padding: 2px 0 6px 0;")
        btrack_lay.addWidget(_bt_info)

        # Per-method physical(µm)/pixel unit toggle. btrack links on
        # micron-scaled coordinates, so its distance thresholds are natively
        # in µm; the toggle lets the user enter them in pixels if preferred.
        self._bt_unit_mgr = UnitGroupManager(
            self.metadata_loader.metadata, default_physical=True
        )
        _bt_unit_row = QHBoxLayout()
        _bt_unit_row.addWidget(QLabel("Distance units:"))
        _bt_unit_row.addWidget(self._bt_unit_mgr.header_row(label=""))
        _bt_unit_wrap = QWidget()
        _bt_unit_wrap.setLayout(_bt_unit_row)
        btrack_lay.addWidget(_bt_unit_wrap)

        # ── Sub-group A: Core Tracking (Step 1) ─────────────────
        step1_group = QGroupBox("Step 1 \u2014 Kalman Filter Tracking")
        step1_form = QFormLayout(step1_group)
        step1_form.setContentsMargins(6, 4, 6, 4)
        step1_form.setSpacing(3)
        step1_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.bt_config_preset = QComboBox()
        self.bt_config_preset.addItems([
            "Cell (bundled)", "Particle (bundled)", "Custom JSON…",
        ])
        preset_val = bt_cfg.get("config_preset", "cell")
        preset_idx = {"cell": 0, "particle": 1}.get(str(preset_val).lower(), 2)
        self.bt_config_preset.setCurrentIndex(preset_idx)
        step1_form.addRow("Config preset:", make_help_row(
            self.bt_config_preset,
            "Configuration Preset",
            "Select a bundled motion/hypothesis model, or provide\n"
            "your own JSON configuration file.\n\n"
            "• Cell — tuned for biological cell tracking\n"
            "• Particle — tuned for small particle tracking\n"
            "• Custom JSON — browse to your own config file"
        ))

        config_path_row = QHBoxLayout()
        self.bt_config_path = QLineEdit()
        self.bt_config_path.setPlaceholderText("Path to custom btrack JSON…")
        self.bt_config_path.setText(bt_cfg.get("config_path", ""))
        self.bt_config_path.setEnabled(preset_idx == 2)
        self.bt_browse_btn = QPushButton("Browse")
        self.bt_browse_btn.setFixedWidth(60)
        self.bt_browse_btn.setEnabled(preset_idx == 2)
        self.bt_browse_btn.clicked.connect(self._bt_browse_config)
        config_path_row.addWidget(self.bt_config_path)
        config_path_row.addWidget(self.bt_browse_btn)
        config_path_widget = QWidget()
        config_path_widget.setLayout(config_path_row)
        step1_form.addRow("", config_path_widget)

        def _on_preset_changed(idx):
            is_custom = idx == 2
            self.bt_config_path.setEnabled(is_custom)
            self.bt_browse_btn.setEnabled(is_custom)
        self.bt_config_preset.currentIndexChanged.connect(_on_preset_changed)

        self.bt_max_search_radius = QSpinBox()
        self.bt_max_search_radius.setRange(1, 999999)
        self.bt_max_search_radius.setValue(int(bt_cfg.get("max_search_radius", 100)))
        self.bt_max_search_radius.setMaximumWidth(90)
        step1_form.addRow("Max search radius:", make_help_row(
            self.bt_max_search_radius,
            "Max Search Radius",
            "Maximum isotropic distance to search for\n"
            "linking objects between frames.\n\n"
            "Entered in the unit selected by the toggle above\n"
            "(µm by default; btrack links in µm natively).\n\n"
            "Increase for fast-moving cells; decrease to\n"
            "prevent long-range false links."
        ))
        self._bt_unit_mgr.register(
            self.bt_max_search_radius, "distance",
            int(bt_cfg.get("max_search_radius", 100)), native_unit="physical",
        )

        self.bt_update_method = QComboBox()
        self.bt_update_method.addItems(["EXACT", "APPROXIMATE"])
        self.bt_update_method.setToolTip(
            "EXACT — full Bayesian belief matrix (accurate, slower).\n"
            "APPROXIMATE — local spatial search (faster, for >1000 cells/frame)."
        )
        self.bt_update_method.setCurrentIndex(
            1 if bt_cfg.get("update_method", "EXACT").upper() == "APPROXIMATE" else 0
        )
        self.bt_update_method.setMaximumWidth(130)
        step1_form.addRow("Update method:", make_help_row(
            self.bt_update_method,
            "Update Method",
            "EXACT — full Bayesian belief matrix (accurate, slower).\n"
            "APPROXIMATE — local spatial search (faster, for >1000\n"
            "cells per frame).\n\n"
            "Use EXACT unless tracking is too slow."
        ))

        self.bt_step_size = QSpinBox()
        self.bt_step_size.setRange(1, 10000)
        self.bt_step_size.setValue(int(bt_cfg.get("step_size", 100)))
        self.bt_step_size.setMaximumWidth(80)
        step1_form.addRow("Step size (frames):", make_help_row(
            self.bt_step_size,
            "Step Size",
            "Number of frames processed per tracking iteration.\n\n"
            "Lower values use less RAM but may be slightly slower."
        ))

        # ── Visual features checkbox ─────────────────────────
        self.bt_use_visual_features = QCheckBox("Use visual features")
        self.bt_use_visual_features.setChecked(bool(bt_cfg.get("use_visual_features", False)))
        visual_features_row = QHBoxLayout()
        visual_features_row.addWidget(self.bt_use_visual_features)
        visual_features_row.addWidget(HelpButton(
            "Visual Features",
            "When enabled, raw image intensity statistics (mean, std per channel)\n"
            "are computed alongside centroids and used by the Kalman filter for\n"
            "more accurate linking.\n\n"
            "Requires raw image data (raw_image_path) in metadata."
        ))
        visual_features_row.addStretch()
        step1_form.addRow("", visual_features_row)

        # ── Workers spinbox ──────────────────────────────────
        n_cores = os.cpu_count() or 4
        max_allowed_bt = max(1, n_cores - 1)
        self.bt_n_workers = QSpinBox()
        self.bt_n_workers.setRange(1, max_allowed_bt)
        self.bt_n_workers.setValue(min(int(bt_cfg.get("n_workers", 1)), max_allowed_bt))
        self.bt_n_workers.setMaximumWidth(60)
        step1_form.addRow("Workers:", make_help_row(
            self.bt_n_workers,
            "Number of Workers",
            f"Number of CPU cores for parallel regionprops extraction\n"
            f"and zarr writing during btrack.\n\n"
            f"Your machine has {n_cores} cores.\n"
            f"Recommendation: Use at most {max_allowed_bt} to keep the system responsive."
        ))

        btrack_lay.addWidget(step1_group)

        # ── Sub-group B: Global Optimizer (Step 2 — opt-in) ─────
        step2_group = QGroupBox("Step 2 — Global Hypothesis Optimizer")
        step2_lay = QVBoxLayout(step2_group)
        step2_lay.setContentsMargins(6, 4, 6, 4)
        step2_lay.setSpacing(3)

        self.bt_use_optimize = QCheckBox("Enable global track optimization")
        self.bt_use_optimize.setChecked(bt_cfg.get("use_optimize", False))
        self.bt_use_optimize.setToolTip(
            "Step 2 runs a global integer-programming optimizer\n"
            "to resolve track merges, splits, and false positives.\n\n"
            "Recommended: first tune Step 1 without optimization,\n"
            "then enable this for refinement."
        )
        step2_lay.addWidget(self.bt_use_optimize)

        # Hypotheses checkboxes
        hyp_group = QGroupBox("Hypotheses")
        hyp_outer_lay = QVBoxLayout(hyp_group)
        hyp_outer_lay.setContentsMargins(4, 2, 4, 2)
        hyp_outer_lay.setSpacing(1)
        hyp_header = QHBoxLayout()
        hyp_header.addWidget(HelpButton(
            "Hypotheses",
            "Each checked hypothesis becomes a candidate explanation the\n"
            "optimizer can assign to a tracklet (short track segment left\n"
            "by Step 1) when resolving ambiguities.\n\n"
            "Distance/time-limited hypotheses (P_link, P_branch, P_dead,\n"
            "P_merge) are constrained by the Distance threshold and Time\n"
            "threshold below.\n\n"
            "P_FP is always required by the optimizer and cannot be "
            "disabled."
        ))
        hyp_header.addStretch()
        hyp_outer_lay.addLayout(hyp_header)
        hyp_lay = QVBoxLayout()
        hyp_lay.setContentsMargins(0, 0, 0, 0)
        hyp_lay.setSpacing(1)
        saved_hyps = bt_cfg.get("hypotheses",
                                ["P_FP", "P_init", "P_term", "P_link"])
        self.bt_hyp_checks = {}
        for hyp_name, hyp_desc, hyp_tooltip, default_on in [
            ("P_FP", "False positive", "Tracklet is a spurious detection "
             "(e.g. segmentation noise) and should be discarded.\n"
             "Always enabled — required by the optimizer.", True),
            ("P_init", "Track initialization", "Tracklet legitimately "
             "starts partway through the movie (a cell entering the "
             "field of view), rather than only at frame 0.", True),
            ("P_term", "Track termination", "Tracklet legitimately ends "
             "partway through the movie (a cell leaving the field of "
             "view), rather than only at the last frame.", True),
            ("P_link", "Track linking", "Two tracklets in different "
             "frames belong to the same object and should be joined "
             "into one track. Limited by Distance/Time threshold.", True),
            ("P_branch", "Track branching", "One tracklet splits into "
             "two (e.g. cell division/mitosis). Limited by "
             "Distance/Time threshold.", False),
            ("P_dead", "Cell death", "Tracklet ends because the cell "
             "died (apoptosis), rather than leaving the field of view "
             "or being a tracking gap.", False),
            ("P_merge", "Track merging", "Two tracklets converge into "
             "one (e.g. cells overlapping/occluding each other). "
             "Limited by Distance/Time threshold.", False),
        ]:
            cb = QCheckBox(f"{hyp_name} — {hyp_desc}")
            cb.setToolTip(hyp_tooltip)
            is_on = hyp_name in saved_hyps if saved_hyps else default_on
            cb.setChecked(is_on)
            if hyp_name == "P_FP":
                cb.setEnabled(False)  # P_FP is always required
            hyp_lay.addWidget(cb)
            self.bt_hyp_checks[hyp_name] = cb
        hyp_outer_lay.addLayout(hyp_lay)
        step2_lay.addWidget(hyp_group)

        step2_form = QFormLayout()
        step2_form.setContentsMargins(0, 0, 0, 0)
        step2_form.setSpacing(3)
        step2_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.bt_dist_thresh = QSpinBox()
        self.bt_dist_thresh.setRange(1, 999999)
        self.bt_dist_thresh.setValue(int(bt_cfg.get("dist_thresh", 60)))
        self.bt_dist_thresh.setMaximumWidth(90)
        step2_form.addRow("Distance threshold:", make_help_row(
            self.bt_dist_thresh,
            "Distance Threshold",
            "Maximum distance for generating link/branch\n"
            "hypotheses in the optimizer.\n\n"
            "Entered in the unit selected by the toggle above\n"
            "(µm by default)."
        ))
        self._bt_unit_mgr.register(
            self.bt_dist_thresh, "distance",
            int(bt_cfg.get("dist_thresh", 60)), native_unit="physical",
        )

        self.bt_time_thresh = QSpinBox()
        self.bt_time_thresh.setRange(1, 999)
        self.bt_time_thresh.setValue(int(bt_cfg.get("time_thresh", 3)))
        self.bt_time_thresh.setMaximumWidth(60)
        step2_form.addRow("Time threshold:", make_help_row(
            self.bt_time_thresh,
            "Time Threshold",
            "Maximum frame gap for generating link\n"
            "hypotheses in the optimizer."
        ))
        step2_lay.addLayout(step2_form)

        btrack_lay.addWidget(step2_group)

        # Toggle Step 2 widgets enabled state
        self._bt_step2_widgets = [
            hyp_group, self.bt_dist_thresh, self.bt_time_thresh,
        ]
        def _on_optimize_toggled(checked):
            for w in self._bt_step2_widgets:
                w.setEnabled(checked)
        self.bt_use_optimize.toggled.connect(_on_optimize_toggled)
        _on_optimize_toggled(self.bt_use_optimize.isChecked())

        btrack_lay.addStretch()
        self._method_pages["btrack"] = self.param_stack.addWidget(btrack_page)

        # Page 6 — Import tracking
        import_page = _ImportTrackingPage(
            cell_type=self.cell_type,
            category=self.category,
            metadata_loader=self.metadata_loader,
            switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
        )
        self._method_pages["import"] = self.param_stack.addWidget(import_page)
        reset_scroll_on_page_change(self.param_stack)

        # Set active page
        _show_params_for_combo_idx(self.combo_method.currentIndex())
        layout.addWidget(self.param_stack)

        # ── Apply settings buttons ──────────────────────────────────
        cat_label = self.category.capitalize() + "s" if self.category != "other" else "Other types"
        sync_row = QHBoxLayout()
        self.btn_apply_cat = QPushButton(f"Apply to all {cat_label}")
        self.btn_apply_cat.clicked.connect(lambda: self._apply_to_others(category_only=True))
        self.btn_apply_all = QPushButton("Apply to all")
        self.btn_apply_all.clicked.connect(lambda: self._apply_to_others(category_only=False))
        
        sync_row.addWidget(self.btn_apply_cat)
        sync_row.addWidget(self.btn_apply_all)
        layout.addLayout(sync_row)

        # ── Run button ───────────────────────────────────────────────
        self.btn_run = QPushButton(f"Run {self.cell_type.capitalize()} Tracking")
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        self.btn_run.clicked.connect(self._on_run_clicked)
        layout.addWidget(self.btn_run)

        # Disable run button for coming-soon methods (only Import tracking)
        def _on_method_idx_changed(idx):
            key = _METHOD_ORDER[idx] if 0 <= idx < len(_METHOD_ORDER) else None
            is_coming_soon = key == "import"
            self.btn_run.setEnabled(not is_coming_soon)
            if is_coming_soon:
                self.btn_run.setToolTip("This tracking method is not yet available.")
            else:
                self.btn_run.setToolTip("")
        self.combo_method.currentIndexChanged.connect(_on_method_idx_changed)
        _on_method_idx_changed(self.combo_method.currentIndex())

        layout.addStretch()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------
    def _get_method_key(self):
        idx = self.combo_method.currentIndex()
        return _METHOD_ORDER[idx] if 0 <= idx < len(_METHOD_ORDER) else _METHOD_ORDER[0]

    def _bt_browse_config(self):
        """Open file dialog for custom btrack JSON config."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Select btrack config JSON", "",
            "JSON files (*.json);;All files (*)",
        )
        if path:
            self.bt_config_path.setText(path)

    def _bt_get_config_preset(self):
        """Return the resolved config preset value."""
        idx = self.bt_config_preset.currentIndex()
        if idx == 0:
            return "cell"
        elif idx == 1:
            return "particle"
        else:
            return self.bt_config_path.text().strip()

    def _bt_get_hypotheses(self):
        """Return list of checked hypothesis names."""
        return [name for name, cb in self.bt_hyp_checks.items() if cb.isChecked()]

    def _collect_params(self) -> dict:
        """Collect current widget values into a dict."""
        return {
            "method": self._get_method_key(),
            "lap": {
                "track_cost_px": int(round(self._lap_unit_mgr.get_native(self.lap_track_cost))),
                "gap_close_cost_px": int(round(self._lap_unit_mgr.get_native(self.lap_gap_cost))),
                "gap_close_max_frames": int(self.lap_gap_frames.value()),
                "merging_cost_px": int(round(self._lap_unit_mgr.get_native(self.lap_merge_cost))),
                "splitting_cost_px": int(round(self._lap_unit_mgr.get_native(self.lap_split_cost))),
            },
            "trackpy": {
                "search_range_px": int(round(self._tp_unit_mgr.get_native(self.tp_search_range))),
                "memory_frames": int(self.tp_memory.value()),
                "adaptive_stop": float(self.tp_adaptive_stop.value()),
                "adaptive_step": float(self.tp_adaptive_step.value()),
            },
            "propagation": {
                # Notice: no tunable params currently exposed
            },
            "bounded_propagation": {
                "min_overlap_fraction": float(self.bp_min_overlap_fraction.value()),
                "segment_size_min": int(self.bp_segment_size_min.value()),
            },
            "reporter_propagation": {
                "min_overlap_fraction": float(self.rp_min_overlap_fraction.value()),
                "segment_size_min": int(self.rp_segment_size_min.value()),
            },
            "btrack": {
                "config_preset": self._bt_get_config_preset(),
                "config_path": self.bt_config_path.text().strip(),
                "use_visual_features": self.bt_use_visual_features.isChecked(),
                "max_search_radius": int(round(self._bt_unit_mgr.get_native(self.bt_max_search_radius))),
                "update_method": "APPROXIMATE" if self.bt_update_method.currentIndex() == 1 else "EXACT",
                "step_size": int(self.bt_step_size.value()),
                "n_workers": max(1, int(self.bt_n_workers.value())),
                "use_optimize": self.bt_use_optimize.isChecked(),
                "hypotheses": self._bt_get_hypotheses(),
                "dist_thresh": int(round(self._bt_unit_mgr.get_native(self.bt_dist_thresh))),
                "time_thresh": int(self.bt_time_thresh.value()),
            },
        }

    def _apply_to_others(self, category_only=False):
        """Copy current settings to other panels."""
        parent_tab = self.parent()
        while parent_tab and not hasattr(parent_tab, 'panels'):
            parent_tab = parent_tab.parent()
        
        if not parent_tab:
            return

        settings = self._collect_params()
        targets = self.category_types if category_only else self.all_cell_types
        
        count = 0
        for ct in targets:
            if ct == self.cell_type:
                continue
            if ct in parent_tab.panels:
                panel = parent_tab.panels[ct]
                # Apply settings
                panel.combo_method.setCurrentIndex(_METHOD_IDX.get(settings["method"], 0))

                # LAP
                panel._lap_unit_mgr.set_native(panel.lap_track_cost, settings["lap"]["track_cost_px"])
                panel._lap_unit_mgr.set_native(panel.lap_gap_cost, settings["lap"]["gap_close_cost_px"])
                panel.lap_gap_frames.setValue(settings["lap"]["gap_close_max_frames"])
                panel._lap_unit_mgr.set_native(panel.lap_merge_cost, settings["lap"]["merging_cost_px"])
                panel._lap_unit_mgr.set_native(panel.lap_split_cost, settings["lap"]["splitting_cost_px"])

                # TrackPy
                panel._tp_unit_mgr.set_native(panel.tp_search_range, settings["trackpy"]["search_range_px"])
                panel.tp_memory.setValue(settings["trackpy"]["memory_frames"])
                panel.tp_adaptive_stop.setValue(settings["trackpy"]["adaptive_stop"])
                panel.tp_adaptive_step.setValue(settings["trackpy"]["adaptive_step"])

                # Bounded Propagation
                bp = settings.get("bounded_propagation", {})
                panel.bp_min_overlap_fraction.setValue(bp.get("min_overlap_fraction", 0.0))
                panel.bp_segment_size_min.setValue(bp.get("segment_size_min", 20))

                # Reporter Propagation
                rp = settings.get("reporter_propagation", {})
                panel.rp_min_overlap_fraction.setValue(rp.get("min_overlap_fraction", 0.1))
                panel.rp_segment_size_min.setValue(rp.get("segment_size_min", 100))

                # btrack
                bt = settings.get("btrack", {})
                preset = bt.get("config_preset", "cell")
                preset_idx = {"cell": 0, "particle": 1}.get(
                    str(preset).lower(), 2
                )
                panel.bt_config_preset.setCurrentIndex(preset_idx)
                panel.bt_config_path.setText(bt.get("config_path", ""))
                panel.bt_use_visual_features.setChecked(bt.get("use_visual_features", False))
                panel._bt_unit_mgr.set_native(panel.bt_max_search_radius, bt.get("max_search_radius", 100))
                panel.bt_update_method.setCurrentIndex(
                    1 if bt.get("update_method", "EXACT").upper() == "APPROXIMATE" else 0
                )
                panel.bt_step_size.setValue(bt.get("step_size", 100))
                panel.bt_n_workers.setValue(bt.get("n_workers", 1))
                panel.bt_use_optimize.setChecked(bt.get("use_optimize", False))
                for hyp_name, cb in panel.bt_hyp_checks.items():
                    if hyp_name == "P_FP":
                        continue
                    cb.setChecked(hyp_name in bt.get("hypotheses", []))
                panel._bt_unit_mgr.set_native(panel.bt_dist_thresh, bt.get("dist_thresh", 60))
                panel.bt_time_thresh.setValue(bt.get("time_thresh", 3))
                
                count += 1
        
        scope = "category" if category_only else "all"
        self.log(f"Applied settings to {count} other cell types ({scope}).")

    def _persist(self):
        """Save this cell type's tracking params to YAML."""
        params = self.metadata_loader.behav3d_parameters
        tracking = params.setdefault("tracking", {})
        tracking[self.cell_type] = self._collect_params()

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    # ------------------------------------------------------------------
    # Running
    # ------------------------------------------------------------------
    def collect_runtime_params(self) -> dict:
        """Snapshot all per-method widget values into a thread-safe dict.

        Called on the Qt thread before kicking off a background worker so
        the worker never needs to touch a Qt widget.  Returns a plain
        ``dict`` of primitive types that fully describes the tracking
        configuration.
        """
        return {
            "method": self._get_method_key(),
            "lap": {
                "track_cost": int(round(self._lap_unit_mgr.get_native(self.lap_track_cost))),
                "gap_cost": int(round(self._lap_unit_mgr.get_native(self.lap_gap_cost))),
                "gap_frames": int(self.lap_gap_frames.value()),
                "merge_cost": int(round(self._lap_unit_mgr.get_native(self.lap_merge_cost))),
                "split_cost": int(round(self._lap_unit_mgr.get_native(self.lap_split_cost))),
            },
            "trackpy": {
                "search_range": int(round(self._tp_unit_mgr.get_native(self.tp_search_range))),
                "memory": int(self.tp_memory.value()),
                "adaptive_stop": float(self.tp_adaptive_stop.value()),
                "adaptive_step": float(self.tp_adaptive_step.value()),
            },
            "bounded_propagation": {
                "min_overlap_fraction": float(self.bp_min_overlap_fraction.value()),
                "segment_size_min": int(self.bp_segment_size_min.value()),
            },
            "reporter_propagation": {
                "min_overlap_fraction": float(self.rp_min_overlap_fraction.value()),
                "segment_size_min": int(self.rp_segment_size_min.value()),
            },
            "btrack": {
                "config_preset": self._bt_get_config_preset(),
                "update_method_idx": int(self.bt_update_method.currentIndex()),
                "use_visual_features": bool(self.bt_use_visual_features.isChecked()),
                "max_search_radius": int(round(self._bt_unit_mgr.get_native(self.bt_max_search_radius))),
                "step_size": int(self.bt_step_size.value()),
                "n_workers": max(1, int(self.bt_n_workers.value())),
                "use_optimize": bool(self.bt_use_optimize.isChecked()),
                "hypotheses": list(self._bt_get_hypotheses()),
                "dist_thresh": int(round(self._bt_unit_mgr.get_native(self.bt_dist_thresh))),
                "time_thresh": int(self.bt_time_thresh.value()),
            },
        }

    def _run_tracking_for(self, cell_type: str, overwrite: bool = False,
                          params: dict = None, progress_cb=None):
        """Run tracking for a single cell type.

        Parameters
        ----------
        cell_type:
            Cell type to track.
        overwrite:
            Whether to overwrite existing tracking outputs.
        params:
            Optional pre-collected snapshot from
            :meth:`collect_runtime_params`.  When ``None`` (queue /
            legacy path) the values are read from the Qt widgets — this
            requires the call to happen on the Qt thread.  Pass a
            snapshot when invoking from a background worker.
        progress_cb:
            Optional ``progress_cb(current, total, label)`` forwarded
            to the backend so a GUI can drive a progress bar.

        Returns
        -------
        pandas.DataFrame
            The updated metadata frame.  Also assigned to
            ``self.metadata_loader.metadata`` for compatibility with
            existing synchronous callers (the queue path).
        """
        if params is None:
            params = self.collect_runtime_params()
        method = params["method"]
        metadata = self.metadata_loader.metadata
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        print(f"\n{'='*50}", file=sys.stderr)
        print(f"  Tracking: {cell_type} ({method.upper()})", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        if method == "lap":
            from behav3d.preprocessing.tracking.laptracking import run_laptracking
            lap = params["lap"]
            tc = int(lap["track_cost"]) ** 2
            gc = int(lap["gap_cost"]) ** 2
            mc = int(lap["merge_cost"])
            sc = int(lap["split_cost"])
            new_md = run_laptracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                track_cost_cutoff=tc,
                gap_closing_cost_cutoff=gc,
                gap_closing_max_frame_count=int(lap["gap_frames"]),
                merging_cost_cutoff=(mc ** 2 if mc > 0 else False),
                splitting_cost_cutoff=(sc ** 2 if sc > 0 else False),
                overwrite=overwrite,
                progress_cb=progress_cb,
            )

        elif method == "trackpy":
            from behav3d.preprocessing.tracking.trackpy_tracking import run_trackpy_tracking_generic
            tp = params["trackpy"]
            new_md = run_trackpy_tracking_generic(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                search_range=int(tp["search_range"]),
                memory=int(tp["memory"]),
                adaptive_stop=float(tp["adaptive_stop"]),
                adaptive_step=float(tp["adaptive_step"]),
                log_callback=ThreadSafeLogger(self.log),
                progress_cb=progress_cb,
            )

        elif method == "btrack":
            from behav3d.preprocessing.tracking.btrack_tracking import run_btracking
            bt = params["btrack"]
            update_method = "APPROXIMATE" if int(bt["update_method_idx"]) == 1 else "EXACT"
            use_opt = bool(bt["use_optimize"])
            new_md = run_btracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                config_preset=bt["config_preset"],
                use_visual_features=bool(bt["use_visual_features"]),
                max_search_radius=int(bt["max_search_radius"]),
                update_method=update_method,
                step_size=int(bt["step_size"]),
                n_workers=int(bt["n_workers"]),
                use_optimize=use_opt,
                hypotheses=list(bt["hypotheses"]) if use_opt else None,
                dist_thresh=int(bt["dist_thresh"]) if use_opt else None,
                time_thresh=int(bt["time_thresh"]) if use_opt else None,
                log_callback=ThreadSafeLogger(self.log),
                progress_cb=progress_cb,
            )

        elif method == "bounded_propagation":
            from behav3d.preprocessing.tracking.bounded_propagation_tracking import run_bounded_propagation_tracking
            bp = params["bounded_propagation"]
            new_md = run_bounded_propagation_tracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                min_overlap_fraction=float(bp["min_overlap_fraction"]),
                segment_size_min=int(bp.get("segment_size_min", 20)),
                progress_cb=progress_cb,
            )

        elif method == "reporter_propagation":
            from behav3d.preprocessing.tracking.reporter_propagation_tracking import run_reporter_propagation_tracking
            rp = params["reporter_propagation"]
            new_md = run_reporter_propagation_tracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                min_overlap_fraction=float(rp["min_overlap_fraction"]),
                segment_size_min=int(rp["segment_size_min"]),
                progress_cb=progress_cb,
            )

        else:  # propagation
            from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
            new_md = run_propagation_tracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                progress_cb=progress_cb,
            )

        print(f"\u2705 {cell_type} tracking finished.", file=sys.stderr)

        # Update shared metadata.  Plain attribute assignment on a
        # DataFrame holder is GIL-protected and the only writer is the
        # active tracking call, so this is safe from the worker.  The
        # caller's ``on_done`` may still re-assign on the Qt thread for
        # parity with downstream code that observes the change there.
        self.metadata_loader.metadata = new_md
        # Save updated metadata CSV (file I/O — safe off the Qt thread).
        csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
        if csv_path:
            new_md.to_csv(csv_path, sep=",", index=False)
        return new_md

    def _check_existing_tracking(self, cell_types: list) -> list:
        """Return descriptions of existing tracking data that would be overwritten."""
        warnings = []
        md = self.metadata_loader.metadata
        if md is None:
            return warnings
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in cell_types:
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                zarr_path = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                csv_dir = out_dir / "trackdata" / sn / ct
                if zarr_path.exists() or (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                    warnings.append(f"{ct} tracking data for {sn}")
        return warnings

    def _on_run_clicked(self):
        """Run tracking for this cell type in the background.

        The Qt-thread phase persists params, prompts for overwrite, and
        snapshots the runtime params; the worker then runs the backend
        with sample-level progress feedback; the ``on_done`` callback
        fires the "switch to visualisation" prompt back on the Qt thread.
        """
        if self._bg.is_running():
            self.log("⚠️ A tracking run is already in progress for this panel.")
            return

        self._persist()
        method = self._get_method_key()
        self.log(f"Running {method.upper()} tracking for: {self.cell_type}")

        # Check for existing data across all samples
        existing = self._check_existing_tracking([self.cell_type])

        if existing:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_single
            choice = prompt_overwrite_single(
                self,
                "Overwrite Existing Tracking?",
                existing,
            )
            if choice != "overwrite":
                self.log(f"Tracking for {self.cell_type} cancelled.")
                return
        overwrite = True

        # Snapshot all Qt widget values now — the worker must not touch them.
        params = self.collect_runtime_params()
        cell_type = self.cell_type

        def _do_tracking(progress_cb=None):
            return self._run_tracking_for(
                cell_type, overwrite=overwrite,
                params=params, progress_cb=progress_cb,
            )

        def _on_done(new_md):
            # Metadata already updated inside _run_tracking_for; this
            # callback is responsible for the post-run dialog only.
            self.log(f"\u2705 {self.cell_type} tracking finished.")
            res = QMessageBox.question(
                self, "Tracking Finished",
                f"Tracking for {self.cell_type} finished! \n\n"
                "Do you want to switch to the Visualization Tab and see the tracks?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res == QMessageBox.Yes:
                self._switch_to_viz_and_show_tracks()
            # Refresh metadata last: this rebuilds the tracking tabs and
            # discards this panel, so nothing may touch ``self`` after it.
            if self.on_tracking_complete is not None:
                self.on_tracking_complete()

        def _on_failed(err: str):
            self.log(f"Error during tracking: {err}")

        self._bg.run(
            fn=_do_tracking,
            desc=f"Tracking {self.cell_type}…",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _switch_to_viz_and_show_tracks(self):
        """Utility to switch to visualization tab and enable track layers."""
        # Stop any running napari animation to avoid stale slider references
        # (prevents RuntimeError: wrapped C/C++ object has been deleted)
        try:
            qt_dims = self.viewer.window._qt_viewer.dims
            if hasattr(qt_dims, '_animation_thread') and qt_dims._animation_thread is not None:
                qt_dims._animation_thread.quit()
                qt_dims._animation_thread.wait()
        except Exception:
            pass

        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(1) # Visualization Tab
            if hasattr(parent, 'visualization_tab'):
                parent.visualization_tab.sample_combo.setCurrentIndex(0)
                parent.visualization_tab._on_load_dataset()
                
                # Make 'tracks' layers visible
                for layer in self.viewer.layers:
                    if "tracks" in layer.name:
                        layer.visible = True


# ═══════════════════════════════════════════════════════════════════════════
# AllOrganoidsPropagationPanel — combined organoid tracking (propagation)
# ═══════════════════════════════════════════════════════════════════════════
class AllOrganoidsPropagationPanel(QWidget):
    """Combined tab used when 'Track all organoids together' is active.

    Locks tracking to Propagation for all organoid types, shows a warning,
    and provides a single run button.  Unchecking the toggle calls back to
    TrackingTab to rebuild individual per-organoid tabs.
    """

    def __init__(self, organoid_types, metadata_loader, log_callback,
                 viewer, toggle_callback, parent=None, tab_progress_row=None,
                 on_tracking_complete_callback=None):
        super().__init__(parent)
        self.organoid_types = list(organoid_types)
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.viewer = viewer
        self._toggle_callback = toggle_callback
        self.tab_progress_row = tab_progress_row
        self.on_tracking_complete = on_tracking_complete_callback
        self._bg = BackgroundOperation(self)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # ── Method group (locked to Propagation) ─────────────────────
        method_group = QGroupBox("Tracking Method")
        method_layout = QHBoxLayout()
        method_layout.setContentsMargins(6, 4, 6, 4)
        method_layout.addWidget(QLabel("Method:"))
        method_label = QLabel("Fragmentation Propagation")
        method_label.setStyleSheet("font-weight: bold;")
        method_layout.addWidget(method_label)
        method_layout.addStretch()
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # ── Warning banner ────────────────────────────────────────────
        warning = QLabel(
            "⚠️  Fragmentation Propagation will be applied to ALL organoid types simultaneously. "
            "To run individual tracking, uncheck 'Track all organoids together'"
        )
        warning.setWordWrap(True)
        warning.setStyleSheet(
            "color: #e65100; background: #fff3e0; border: 1px solid #ffcc80; "
            "border-radius: 4px; padding: 8px; font-size: 11px; margin: 2px 0;"
        )
        layout.addWidget(warning)

        # ── Propagation parameters group ──────────────────────────────
        prop_group = QGroupBox("Fragmentation Propagation")
        prop_lay = QVBoxLayout(prop_group)
        prop_lay.setSpacing(6)

        self.check_all_together = QCheckBox("Track all organoids together")
        self.check_all_together.setChecked(True)
        self.check_all_together.setToolTip(
            "When checked, all organoid types are tracked simultaneously using Fragmentation Propagation.\n"
            "Uncheck to configure and run tracking independently per organoid type."
        )
        self.check_all_together.toggled.connect(self._on_toggled)
        prop_lay.addWidget(self.check_all_together)

        notice = QLabel("No additional tunable parameters for Fragmentation Propagation.")
        notice.setWordWrap(True)
        notice.setStyleSheet("color: #666; font-style: italic; padding: 2px 0;")
        prop_lay.addWidget(notice)
        layout.addWidget(prop_group)

        layout.addStretch()

        # ── Run All Organoids Tracking button ─────────────────────────
        types_str = ", ".join(self.organoid_types)
        self.btn_run = QPushButton(f"▶  Run All Organoids Tracking  ({types_str})")
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 8px; font-size: 13px;"
        )
        self.btn_run.clicked.connect(self._on_run_clicked)
        layout.addWidget(self.btn_run)

    def _on_toggled(self, checked):
        """Relay toggle to TrackingTab so it can rebuild tabs."""
        if self._toggle_callback:
            self._toggle_callback(checked)

    def _on_run_clicked(self):
        """Run propagation tracking for all organoid types in the background."""
        if self._bg.is_running():
            self.log("⚠️ All-organoids tracking is already in progress.")
            return

        md = self.metadata_loader.metadata
        out_dir = self.metadata_loader.output_dir
        if md is None or not out_dir:
            self.log("⚠️ No metadata loaded.")
            return

        out_path = Path(out_dir)
        existing = []
        for ct in self.organoid_types:
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                zarr_path = out_path / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                csv_dir = out_path / "trackdata" / sn / ct
                if zarr_path.exists() or (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                    existing.append(f"{ct} tracking data for {sn}")

        overwrite = True
        if existing:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
            choice = prompt_overwrite_batch(
                self,
                "Overwrite Existing Tracking?",
                existing,
            )
            if choice == "cancel":
                self.log("All-organoids tracking cancelled.")
                return
            overwrite = (choice == "overwrite")

        types_str = ", ".join(self.organoid_types)
        self.btn_run.setText("\u23f3 Running\u2026")
        self.log(f"\u25b6 Fragmentation Propagation \u2014 all organoids: {types_str}\u2026")

        out_dir_str = str(out_path.expanduser())
        first_ct = self.organoid_types[0]

        def _do_all_organoids(progress_cb=None):
            from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
            new_md = run_propagation_tracking(
                metadata=md,
                output_dir=out_dir_str,
                cell_type=first_ct,
                overwrite=overwrite,
                all_organoids=True,
                progress_cb=progress_cb,
            )
            # Persist to CSV from the worker — file I/O is safe off the
            # Qt thread.  metadata_loader.metadata is reassigned in
            # _on_done so the change is observed on the Qt thread.
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                new_md.to_csv(csv_path, sep=",", index=False)
            return new_md

        def _on_done(new_md):
            self.metadata_loader.metadata = new_md
            self.log("\u2705 All organoids propagation tracking finished.")
            self.btn_run.setText(f"\u25b6  Run All Organoids Tracking  ({types_str})")
            res = QMessageBox.question(
                self, "Tracking Finished",
                "All organoid types tracked successfully!\n\n"
                "Do you want to switch to the Visualization Tab to see the tracks?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res == QMessageBox.Yes:
                self._switch_to_viz()
            # Refresh metadata last: this rebuilds the tracking tabs and
            # discards this panel, so nothing may touch ``self`` after it.
            if self.on_tracking_complete is not None:
                self.on_tracking_complete()

        def _on_failed(err: str):
            self.log(f"\u274c Error during all-organoids tracking: {err}")
            self.btn_run.setText(f"\u25b6  Run All Organoids Tracking  ({types_str})")

        self._bg.run(
            fn=_do_all_organoids,
            desc="Tracking all organoids\u2026",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _switch_to_viz(self):
        parent = self.parent()
        while parent and not hasattr(parent, "tabs"):
            parent = parent.parent()
        if parent and hasattr(parent, "tabs"):
            parent.tabs.setCurrentIndex(1)
            if hasattr(parent, "visualization_tab"):
                parent.visualization_tab.sample_combo.setCurrentIndex(0)
                parent.visualization_tab._on_load_dataset()
                for layer in self.viewer.layers:
                    if "tracks" in layer.name:
                        layer.visible = True


# ═══════════════════════════════════════════════════════════════════════════
# MulticolorTrackingPanel — single combined tab for multicolor cell types
# ═══════════════════════════════════════════════════════════════════════════
class MulticolorTrackingPanel(QWidget):
    """Combined tracking tab for a multicolor cell type group.

    Groups per-channel types (e.g. tcell_1_multicolor, tcell_2_multicolor)
    under a single UI.  Parameters are shared across all channels; after
    tracking each channel individually the outputs are automatically merged
    into ``<base_name>_merged``.
    """

    def __init__(self, base_name: str, channel_types: list, category: str,
                 metadata_loader, log_callback=None, viewer=None, parent=None,
                 tab_progress_row=None, on_tracking_complete_callback=None,
                 switch_to_data_prep_edit_callback=None):
        super().__init__(parent)
        self.base_name = base_name
        self.channel_types = sorted(channel_types)
        self.category = category
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.viewer = viewer
        self.tab_progress_row = tab_progress_row
        self.on_tracking_complete = on_tracking_complete_callback

        # Inner panel drives the shared parameter UI; keyed to first channel
        # so saved config is stored/loaded under that channel name.
        first_channel = self.channel_types[0]
        self._inner_panel = CellTypeTrackingPanel(
            cell_type=first_channel,
            category=category,
            metadata_loader=metadata_loader,
            all_cell_types=list(channel_types),
            category_types=list(channel_types),
            log_callback=log_callback,
            viewer=viewer,
            tab_progress_row=tab_progress_row,
            switch_to_data_prep_edit_callback=switch_to_data_prep_edit_callback,
        )
        # A multicolor run executes on the inner panel's runner; alias it here
        # so TrackingTab.request_tab_exit — which looks for ``_bg`` on each
        # panel it holds — can see the run and block the tab switch.
        self._bg = self._inner_panel._bg

        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # ── Warning banner ────────────────────────────────────────────
        channels_str = ", ".join(self.channel_types)
        merged_name = f"{self.base_name}_merged"
        warning = QLabel(
            f"⚠️  Multicolor cell type detected.  Tracking will be applied "
            f"independently to each channel ({channels_str}) using the shared "
            f"parameters below, then the outputs will be automatically merged "
            f"into <b>{merged_name}</b>."
        )
        warning.setWordWrap(True)
        warning.setTextFormat(Qt.RichText)
        warning.setStyleSheet(
            "color: #e65100; background: #fff3e0; border: 1px solid #ffcc80; "
            "border-radius: 4px; padding: 8px; font-size: 11px; margin: 2px 0;"
        )
        layout.addWidget(warning)

        # ── Channel summary ───────────────────────────────────────────
        channels_group = QGroupBox("Channels to track (same parameters applied to all)")
        ch_lay = QVBoxLayout(channels_group)
        ch_lay.setContentsMargins(6, 4, 6, 4)
        ch_lay.setSpacing(2)
        for i, ch in enumerate(self.channel_types, 1):
            lbl = QLabel(f"  {i}. {ch}  →  {self.base_name}_merged")
            lbl.setStyleSheet("font-size: 11px; color: #555;")
            ch_lay.addWidget(lbl)
        layout.addWidget(channels_group)

        # ── Shared parameters (inner CellTypeTrackingPanel) ───────────
        layout.addWidget(self._inner_panel)

        # The inner panel ships its own single-channel Run button.  Running
        # only the first channel would leave the remaining channels untracked
        # and never produce the merged output, so the button is retargeted at
        # the whole group (all channels + merge).
        btn = self._inner_panel.btn_run
        btn.setText(
            f"▶  Run {self.base_name.capitalize()} Multicolor Tracking  "
            f"({len(self.channel_types)} channels + merge)"
        )
        group_tip = (
            "Tracks " + ", ".join(self.channel_types)
            + f" with the shared parameters above, then merges them into {merged_name}."
        )
        btn.setToolTip(group_tip)
        try:
            btn.clicked.disconnect(self._inner_panel._on_run_clicked)
        except (TypeError, RuntimeError):
            pass
        btn.clicked.connect(self._on_run_clicked)

        # The inner panel resets the tooltip whenever the method changes;
        # restore the group wording afterwards (it stays cleared for the
        # not-yet-available methods, which disable the button).
        def _restore_group_tooltip(_idx=None):
            if btn.isEnabled():
                btn.setToolTip(group_tip)
        self._restore_group_tooltip = _restore_group_tooltip
        self._inner_panel.combo_method.currentIndexChanged.connect(_restore_group_tooltip)

    # ------------------------------------------------------------------
    def _persist(self):
        """Persist shared parameters for all channels."""
        params = self.metadata_loader.behav3d_parameters
        tracking = params.setdefault("tracking", {})
        shared_cfg = self._inner_panel._collect_params()
        for ch in self.channel_types:
            tracking[ch] = shared_cfg

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    def _on_run_clicked(self):
        """Run tracking for every channel of this multicolor group + merge.

        Mirrors :meth:`CellTypeTrackingPanel._on_run_clicked` but operates on
        the whole channel family: the shared parameters are persisted for all
        channels, existing outputs for every channel *and* the merged type are
        checked in a single overwrite prompt, and the background worker runs
        :meth:`run_tracking` (per-channel tracking followed by the merge).
        """
        inner = self._inner_panel
        if inner._bg.is_running():
            self.log("⚠️ A tracking run is already in progress for this panel.")
            return

        self._persist()
        merged_name = f"{self.base_name}_merged"
        method = inner._get_method_key()
        self.log(
            f"Running {method.upper()} tracking for multicolor group "
            f"{self.base_name}: {', '.join(self.channel_types)} → {merged_name}"
        )

        # One prompt covering every channel plus the merged output.
        existing = inner._check_existing_tracking(
            list(self.channel_types) + [merged_name]
        )
        if existing:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_single
            choice = prompt_overwrite_single(
                self,
                "Overwrite Existing Tracking?",
                existing,
            )
            if choice != "overwrite":
                self.log(f"Multicolor tracking for {self.base_name} cancelled.")
                return
        overwrite = True

        # Snapshot Qt widget values now — the worker must not touch them.
        params = inner.collect_runtime_params()

        def _do_tracking(progress_cb=None):
            return self.run_tracking(
                overwrite=overwrite, params=params, progress_cb=progress_cb,
            )

        def _on_done(_new_md):
            self.log(
                f"✅ {self.base_name} multicolor tracking finished "
                f"({len(self.channel_types)} channels merged → {merged_name})."
            )
            res = QMessageBox.question(
                self, "Tracking Finished",
                f"Multicolor tracking for {self.base_name} finished!\n\n"
                f"Channels tracked: {', '.join(self.channel_types)}\n"
                f"Merged output: {merged_name}\n\n"
                "Do you want to switch to the Visualization Tab and see the tracks?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res == QMessageBox.Yes:
                inner._switch_to_viz_and_show_tracks()
            # Refresh metadata last: this rebuilds the tracking tabs and
            # discards this panel, so nothing may touch ``self`` after it.
            if self.on_tracking_complete is not None:
                self.on_tracking_complete()

        def _on_failed(err: str):
            self.log(f"Error during multicolor tracking: {err}")

        inner._bg.run(
            fn=_do_tracking,
            desc=f"Tracking {self.base_name} (multicolor)…",
            progress_row=self.tab_progress_row,
            buttons=[inner.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    # ------------------------------------------------------------------
    def run_tracking(self, overwrite: bool = False, params: dict = None,
                     progress_cb=None):
        """Track each channel then merge into ``<base_name>_merged``.

        ``params`` is an optional Qt-thread-safe snapshot from
        :meth:`CellTypeTrackingPanel.collect_runtime_params`.  When
        ``None`` the inner panel's widget values are read directly (only
        safe on the Qt thread).  ``progress_cb`` drives a global
        per-channel-per-sample progress bar.
        """
        n = len(self.channel_types)
        merged_name = f"{self.base_name}_merged"

        # Wrap per-channel progress into a global counter so the bar runs
        # linearly across all channels rather than restarting.
        per_channel_total = [0] * n  # filled in lazily by the first tick
        completed = [0]

        def _make_channel_cb(channel_idx, channel_name):
            def _cb(curr, total, label):
                per_channel_total[channel_idx] = max(per_channel_total[channel_idx], int(total))
                # Estimated grand total assumes equal samples per channel;
                # we refine as new channels report their totals.
                if all(t > 0 for t in per_channel_total[: channel_idx + 1]):
                    avg_per_ch = sum(per_channel_total[: channel_idx + 1]) / (channel_idx + 1)
                    grand_total = int(round(avg_per_ch * n))
                else:
                    grand_total = max(int(total) * n, 1)
                global_curr = completed[0] + int(curr)
                if progress_cb is not None:
                    try:
                        progress_cb(global_curr, grand_total,
                                    f"{channel_name} ({channel_idx + 1}/{n}): {label}")
                    except Exception:
                        pass
            return _cb

        for i, ch in enumerate(self.channel_types):
            self.log(f"  Multicolor [{i + 1}/{n}] Tracking channel: {ch}…")
            print(f"\n  Multicolor [{i + 1}/{n}] Tracking: {ch}", file=sys.stderr)
            self._inner_panel._run_tracking_for(
                ch,
                overwrite=overwrite,
                params=params,
                progress_cb=_make_channel_cb(i, ch),
            )
            completed[0] += per_channel_total[i] if per_channel_total[i] else 0
            self.log(f"  Done: {ch}")

        self.log(f"  Merging {n} channels → {merged_name}…")
        print(f"\n  Merging multicolor channels → {merged_name}", file=sys.stderr)
        try:
            from behav3d.preprocessing.tracking.multicolor_tracking_processing import (
                combine_multicolor_tracked_outputs,
            )
            new_md = combine_multicolor_tracked_outputs(
                metadata=self.metadata_loader.metadata,
                output_dir=str(Path(self.metadata_loader.output_dir).expanduser()),
                source_cell_types=self.channel_types,
                combined_cell_type=merged_name,
                overwrite=overwrite,
            )
            self.metadata_loader.metadata = new_md
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                new_md.to_csv(csv_path, sep=",", index=False)
            self.log(f"  ✅ Merged → {merged_name}")
            print(f"  ✅ Merged → {merged_name}", file=sys.stderr)
        except Exception as e:
            traceback.print_exc()
            self.log(f"  ❌ Merge failed: {e}")
            raise
        return self.metadata_loader.metadata


# ═══════════════════════════════════════════════════════════════════════════
# TrackingTab — main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class TrackingTab(QWidget):
    # Signal emitted when tracking is completed
    tracking_completed = Signal()

    def __init__(self, viewer, metadata_loader, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeTrackingPanel] = {}
        self._multicolor_panels: dict[str, MulticolorTrackingPanel] = {}

        # Background-execution infrastructure for the batch-run button.
        self._bg = BackgroundOperation(self)

        self._init_ui()

        # Listen to metadata changes
        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    # ------------------------------------------------------------------
    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        scroll.setWidget(content)

        # Per-rebuild organoid state
        self._organoid_types = []
        self._all_organoids_panel = None

        # Sub-tab widget (West position = left tabs)
        self.cell_tabs = QTabWidget()
        self.cell_tabs.setTabPosition(QTabWidget.West)
        layout.addWidget(self.cell_tabs)
        reset_scroll_on_page_change(self.cell_tabs)

        # Global Run Button + Queue button
        self.btn_run_batch = QPushButton("Run Batch Tracking (All Cell Types)")
        self.btn_run_batch.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_batch.clicked.connect(self._on_run_batch_clicked)

        self.btn_queue_track = QPushButton("+🛒")
        self.btn_queue_track.setFixedSize(36, 32)
        self.btn_queue_track.setToolTip("Add Batch Tracking to Processing Queue")
        self.btn_queue_track.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_btn_row = QHBoxLayout()
        batch_btn_row.setSpacing(4)
        batch_btn_row.addWidget(self.btn_run_batch, stretch=1)
        batch_btn_row.addWidget(self.btn_queue_track)
        layout.addLayout(batch_btn_row)

        # Hidden until metadata is loaded
        self.btn_run_batch.setVisible(False)
        self.btn_queue_track.setVisible(False)

        # Shared progress row — every Run button on this tab feeds this
        # single widget.  Hidden until a run is active.
        self.progress_row = ProgressBarRow()
        layout.addWidget(self.progress_row)

        # Log window
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(120)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        # Placeholder when no metadata
        self._placeholder = QLabel("Load metadata in the Data Preparation tab to see tracking options.")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setWordWrap(True)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.cell_tabs.addTab(self._placeholder, "—")

    def _log(self, msg):
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    def _switch_to_data_prep_edit(self):
        """Switch the main window to the Data Preparation tab with the
        Metadata Builder already open in edit mode.

        Used by ``_ImportTrackingPage`` (the "Add a new sample or cell type"
        shortcut) — mirrors ``SegmentationTab._switch_to_data_prep_edit``.
        """
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(0)
            if hasattr(parent, 'data_prep_tab'):
                parent.data_prep_tab.enter_metadata_edit_mode()

    # ------------------------------------------------------------------
    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing tracking tabs…")
        self._rebuild_tabs()

    # ------------------------------------------------------------------
    def request_tab_exit(self) -> bool:
        """Block tab switching while a background tracking run is in flight."""
        from qtpy.QtWidgets import QMessageBox

        if self._bg.is_running():
            QMessageBox.information(
                self,
                "Operation in progress",
                "A tracking run is still in progress. Please wait for it "
                "to finish before switching tabs.",
            )
            return False
        for panel in list(self.panels.values()) + list(self._multicolor_panels.values()):
            bg = getattr(panel, "_bg", None)
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "A tracking run is still in progress. Please wait for "
                    "it to finish before switching tabs.",
                )
                return False
        all_org = getattr(self, "_all_organoids_panel", None)
        if all_org is not None:
            bg = getattr(all_org, "_bg", None)
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "A tracking run is still in progress. Please wait for "
                    "it to finish before switching tabs.",
                )
                return False
        return True

    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            is_combined_multicolor_celltype,
        )
        md = self.metadata_loader.metadata
        if md is None:
            return [], [], []
        # Filter out *_merged / *_grouped types — these are tracking outputs,
        # not inputs.  Per-channel *_N_multicolor types are kept and later
        # grouped into a single MulticolorTrackingPanel per base name.
        return (
            [ct for ct in detect_organoid_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)],
            [ct for ct in detect_immune_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)],
            [ct for ct in detect_other_cell_types_from_metadata(md) if not is_combined_multicolor_celltype(ct)],
        )

    def _rebuild_tabs(self):
        # Clear old tabs
        self.cell_tabs.clear()
        self.panels.clear()
        self._multicolor_panels.clear()
        self._all_organoids_panel = None

        org, imm, oth = self._detect_cell_types()

        # ── Partition each category into standalone and multicolor groups ─
        from behav3d.core.metadata import is_multicolor_celltype, multicolor_base_name

        def _partition(types):
            standalone, mc_groups = [], {}
            for ct in types:
                if is_multicolor_celltype(ct):
                    base = multicolor_base_name(ct)
                    mc_groups.setdefault(base, []).append(ct)
                else:
                    standalone.append(ct)
            return standalone, mc_groups

        org_standalone, org_mc = _partition(org)
        imm_standalone, imm_mc = _partition(imm)
        oth_standalone, oth_mc = _partition(oth)

        all_standalone = org_standalone + imm_standalone + oth_standalone
        all_types = all_standalone  # context list for CellTypeTrackingPanel

        has_anything = all_standalone or org_mc or imm_mc or oth_mc
        if not has_anything:
            self.cell_tabs.addTab(self._placeholder, "—")
            self.btn_run_batch.setVisible(False)
            self.btn_queue_track.setVisible(False)
            return

        self.btn_run_batch.setVisible(True)
        self.btn_queue_track.setVisible(True)

        # ── Organoid standalone ───────────────────────────────────────
        self._organoid_types = list(org_standalone)

        track_together = False
        if org_standalone:
            tcfg = _cfg_get(self.metadata_loader.behav3d_parameters, "tracking", {}) or {}
            track_together = tcfg.get("track_organoids_together", True)

        if track_together and org_standalone:
            self._all_organoids_panel = AllOrganoidsPropagationPanel(
                organoid_types=org_standalone,
                metadata_loader=self.metadata_loader,
                log_callback=self._log,
                viewer=self.viewer,
                toggle_callback=self._on_all_organoids_toggled,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
            )
            self.cell_tabs.addTab(self._all_organoids_panel, "🟣 All Organoids")
        else:
            for ct in org_standalone:
                panel = CellTypeTrackingPanel(
                    cell_type=ct, category="organoid",
                    metadata_loader=self.metadata_loader,
                    all_cell_types=all_types, category_types=org_standalone,
                    log_callback=self._log,
                    viewer=self.viewer,
                    tab_progress_row=self.progress_row,
                    on_tracking_complete_callback=self.tracking_completed.emit,
                    switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
                )
                panel._toggle_all_organoids_callback = self._on_all_organoids_toggled
                self.panels[ct] = panel
                self.cell_tabs.addTab(panel, f"🟣 {ct.capitalize()}")

        # ── Organoid multicolor groups ────────────────────────────────
        for base, channels in org_mc.items():
            panel = MulticolorTrackingPanel(
                base_name=base, channel_types=channels, category="organoid",
                metadata_loader=self.metadata_loader,
                log_callback=self._log, viewer=self.viewer,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
                switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
            )
            self._multicolor_panels[base] = panel
            self.cell_tabs.addTab(panel, f"🟣 {base.capitalize()} (multicolor)")

        # ── Immune standalone ─────────────────────────────────────────
        for ct in imm_standalone:
            panel = CellTypeTrackingPanel(
                cell_type=ct, category="immune",
                metadata_loader=self.metadata_loader,
                all_cell_types=all_types, category_types=imm_standalone,
                log_callback=self._log,
                viewer=self.viewer,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
                switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🔵 {ct.capitalize()}")

        # ── Immune multicolor groups ──────────────────────────────────
        for base, channels in imm_mc.items():
            panel = MulticolorTrackingPanel(
                base_name=base, channel_types=channels, category="immune",
                metadata_loader=self.metadata_loader,
                log_callback=self._log, viewer=self.viewer,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
                switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
            )
            self._multicolor_panels[base] = panel
            self.cell_tabs.addTab(panel, f"🔵 {base.capitalize()} (multicolor)")

        # ── Other standalone ──────────────────────────────────────────
        for ct in oth_standalone:
            panel = CellTypeTrackingPanel(
                cell_type=ct, category="other",
                metadata_loader=self.metadata_loader,
                all_cell_types=all_types, category_types=oth_standalone,
                log_callback=self._log,
                viewer=self.viewer,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
                switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🟡 {ct.capitalize()}")

        # ── Other multicolor groups ───────────────────────────────────
        for base, channels in oth_mc.items():
            panel = MulticolorTrackingPanel(
                base_name=base, channel_types=channels, category="other",
                metadata_loader=self.metadata_loader,
                log_callback=self._log, viewer=self.viewer,
                tab_progress_row=self.progress_row,
                on_tracking_complete_callback=self.tracking_completed.emit,
                switch_to_data_prep_edit_callback=self._switch_to_data_prep_edit,
            )
            self._multicolor_panels[base] = panel
            self.cell_tabs.addTab(panel, f"🟡 {base.capitalize()} (multicolor)")

    # ------------------------------------------------------------------
    def _on_all_organoids_toggled(self, checked):
        """Persist 'track_organoids_together' setting and rebuild tabs."""
        params = self.metadata_loader.behav3d_parameters
        tracking = params.setdefault("tracking", {})
        tracking["track_organoids_together"] = bool(checked)
        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception:
                pass
        self._rebuild_tabs()

    def _on_run_batch_clicked(self):
        """User-triggered batch run — runs asynchronously."""
        self.run_batch_tracking(interactive=True, block=False)

    def run_batch_tracking(self, interactive=True, skip_existing=False, block=True,
                           extra_callbacks=None):
        """Sequential run for all configured cell type panels.

        When ``interactive=False``, skips overwrite and visualization
        dialogs (this is what the processing queue uses).

        When ``block=True`` (default) everything runs synchronously on
        the caller's thread — required by the processing queue which
        iterates steps sequentially.  When ``block=False`` the heavy
        work is moved to a background worker (used by the GUI batch
        button) so napari stays interactive and a progress bar is shown.

        ``extra_callbacks`` is the queue's chaining hook
        (``{"on_done": cb, "on_failed": cb}``).
        """
        if not self.panels and not self._multicolor_panels and not self._all_organoids_panel:
            self._log("No cell type panels to track.")
            fire_extra_callback(extra_callbacks, "on_failed", "no cell type panels")
            return

        if not block and self._bg.is_running():
            self._log("⚠️ A batch tracking run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        all_organoids_mode = self._all_organoids_panel is not None
        total = (
            len(self.panels)
            + (len(self._organoid_types) if all_organoids_mode else 0)
            + len(self._multicolor_panels)
        )
        self._log(f"Starting batch tracking for {total} cell types...")
        if all_organoids_mode:
            self._log("  Mode: All organoids tracked together (propagation)")
        if self._multicolor_panels:
            mc_names = ", ".join(
                f"{b} ({len(p.channel_types)} channels)"
                for b, p in self._multicolor_panels.items()
            )
            self._log(f"  Multicolor groups: {mc_names}")

        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Running batch tracking for {total} cell types", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        # Check for existing tracking data across all standalone cell types
        all_cts = list(self.panels.keys())
        existing = []
        existing_cts = set()
        out_dir = Path(self.metadata_loader.output_dir)
        md = self.metadata_loader.metadata
        for ct in all_cts:
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                zarr_path = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                csv_dir = out_dir / "trackdata" / sn / ct
                if zarr_path.exists() or (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                    existing.append(f"{ct} tracking data for {sn}")
                    existing_cts.add(ct)

        skip_existing_flag = skip_existing
        overwrite = not skip_existing
        if existing and interactive:
            from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
            choice = prompt_overwrite_batch(
                self,
                "Overwrite Existing Tracking?",
                existing,
            )
            if choice == "cancel":
                self._log("Batch tracking cancelled by user.")
                fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                return
            skip_existing_flag = (choice == "skip")
            overwrite = not skip_existing_flag

        # Persist all panel params before starting so config is saved even if an error occurs
        for ct, panel in self.panels.items():
            panel._persist()
        for base, mc_panel in self._multicolor_panels.items():
            mc_panel._persist()

        # Snapshot every panel's Qt widget state so the worker thread
        # never has to read from Qt (only safe to read widgets here).
        panel_params: dict[str, dict] = {
            ct: panel.collect_runtime_params() for ct, panel in self.panels.items()
        }
        mc_params: dict[str, dict] = {
            base: mc_panel._inner_panel.collect_runtime_params()
            for base, mc_panel in self._multicolor_panels.items()
        }

        # Cell-type-level progress: one "step" per cell-type panel.  We
        # scale each step into 100 sub-units so the activity dock can
        # estimate ETA from sample-level ticks within a step.
        steps_total = total  # noqa: F841 (kept for readability / future use)
        steps_done = [0]

        def _do_batch(progress_cb=None):
            """Worker function.  ``progress_cb`` is None in block=True mode."""

            def _step_progress(step_idx, label_prefix):
                """Make a per-step progress wrapper that scales to the
                global bar.  ``step_idx`` is 0-based for the current
                step that's executing."""
                def _cb(curr, sub_total, label):
                    if progress_cb is None:
                        return
                    full_label = f"{label_prefix}: {label}" if label else label_prefix
                    # Per-step sub-progress contributes a fractional
                    # advance on top of completed steps.
                    if sub_total and sub_total > 0:
                        frac = min(max(curr / sub_total, 0.0), 1.0)
                    else:
                        frac = 0.0
                    # Scale: each step is one "unit"; show progress in
                    # 100ths to surface ETA via the activity dock.
                    SCALE = 100
                    grand = steps_total * SCALE
                    val = step_idx * SCALE + int(frac * SCALE)
                    try:
                        progress_cb(val, grand, full_label)
                    except Exception:
                        pass
                return _cb

            # ── All-organoids propagation (single step) ─────────────
            organoids_done = set()
            if all_organoids_mode and self._organoid_types:
                org_cts_to_run = [
                    ct for ct in self._organoid_types
                    if not (skip_existing_flag and ct in existing_cts)
                ]
                if org_cts_to_run:
                    self._log(f"--- Tracking ALL organoids together ({', '.join(self._organoid_types)}) ---")
                    print(f"\n▶ Tracking ALL organoids together...", file=sys.stderr)

                    first_org_ct = self._organoid_types[0]
                    from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
                    new_md = run_propagation_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(Path(self.metadata_loader.output_dir).expanduser()),
                        cell_type=first_org_ct,
                        overwrite=overwrite,
                        all_organoids=True,
                        progress_cb=_step_progress(
                            steps_done[0], "all organoids"
                        ) if progress_cb is not None else None,
                    )
                    self.metadata_loader.metadata = new_md
                    csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
                    if csv_path:
                        new_md.to_csv(csv_path, sep=",", index=False)

                    self._log(f"Done: all organoids propagation tracking")
                    organoids_done = set(self._organoid_types)
                    steps_done[0] += len(self._organoid_types)
                else:
                    self._log("All organoid tracking data already exists — skipping.")
                    organoids_done = set(self._organoid_types)
                    steps_done[0] += len(self._organoid_types)

            # ── Per-cell-type tracking ──────────────────────────────
            for i, (ct, panel) in enumerate(self.panels.items(), 1):
                if ct in organoids_done:
                    continue
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- Skipping {ct} (existing data) ---")
                    steps_done[0] += 1
                    continue
                print(f"\n▶ [{steps_done[0] + 1}/{total}] Tracking {ct}...", file=sys.stderr)
                self._log(f"--- [{steps_done[0] + 1}/{total}] Tracking {ct} ---")
                panel._run_tracking_for(
                    ct, overwrite=overwrite,
                    params=panel_params.get(ct),
                    progress_cb=_step_progress(steps_done[0], ct)
                    if progress_cb is not None else None,
                )

                # Verify tracked outputs exist for this cell type after tracking.
                missing_outputs = []
                for _, sample in self.metadata_loader.metadata.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    zarr_path = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                    csv_dir = out_dir / "trackdata" / sn / ct
                    has_csv = csv_dir.exists() and bool(list(csv_dir.glob("*.csv")))
                    if not zarr_path.exists() or not has_csv:
                        missing_outputs.append(str(sn))
                if missing_outputs:
                    self._log(
                        f"⚠️ {ct}: tracked outputs missing for {len(missing_outputs)} sample(s): "
                        + ", ".join(missing_outputs)
                    )
                self._log(f"Done: {ct}")
                steps_done[0] += 1

            # ── Multicolor groups (track channels + auto-merge) ───────
            for base, mc_panel in self._multicolor_panels.items():
                merged_name = f"{base}_merged"
                print(f"\n▶ [{steps_done[0] + 1}/{total}] Multicolor tracking: {base} "
                      f"({len(mc_panel.channel_types)} channels + merge)…",
                      file=sys.stderr)
                self._log(f"--- [{steps_done[0] + 1}/{total}] Multicolor tracking: {base} "
                          f"({', '.join(mc_panel.channel_types)}) → {merged_name} ---")
                mc_panel.run_tracking(
                    overwrite=overwrite,
                    params=mc_params.get(base),
                    progress_cb=_step_progress(steps_done[0], f"multicolor {base}")
                    if progress_cb is not None else None,
                )
                self._log(f"Done: {base} (channels tracked + merged → {merged_name})")
                steps_done[0] += 1

            # Persist final metadata state after full batch execution.
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                try:
                    self.metadata_loader.metadata.to_csv(csv_path, sep=",", index=False)
                    self._log("Saved metadata CSV after batch tracking.")
                except Exception as e:
                    self._log(f"⚠️ Could not save metadata CSV after batch tracking: {e}")

            self._log("\u2705 Batch tracking finished.")
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"  \u2705 Batch tracking complete", file=sys.stderr)
            print(f"{'='*60}\n", file=sys.stderr)
            return self.metadata_loader.metadata

        # ── Dispatch: synchronous (queue) or async (GUI button) ─────
        if block:
            try:
                result = _do_batch(progress_cb=None)
                # Emit tracking completion signal for metadata refresh — must
                # happen before the Visualization-tab prompt so that tab loads
                # the refreshed metadata rather than its stale copy.
                self.tracking_completed.emit()
                if interactive:
                    self._prompt_switch_to_viz_after_batch()
                fire_extra_callback(extra_callbacks, "on_done", result)
            except Exception as e:
                traceback.print_exc()
                self._log(f"\u274c Batch tracking error: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
            return

        def _on_done(result):
            # Emit tracking completion signal for metadata refresh before the
            # Visualization-tab prompt (see the blocking branch above).
            self.tracking_completed.emit()
            if interactive:
                self._prompt_switch_to_viz_after_batch()
            fire_extra_callback(extra_callbacks, "on_done", result)

        def _on_failed(err: str):
            self._log(f"\u274c Batch tracking error: {err}")
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_batch,
            desc=f"Batch tracking ({total} cell types)\u2026",
            progress_row=self.progress_row,
            buttons=[self.btn_run_batch],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _prompt_switch_to_viz_after_batch(self):
        """Post-batch dialog: offer to jump to the Visualization tab."""
        res = QMessageBox.question(
            self, "Batch Tracking Finished",
            "Successfully tracked all cell types! \n\n"
            "Do you want to switch to the Visualization Tab and see the tracks?",
            QMessageBox.Yes | QMessageBox.No,
        )
        if res == QMessageBox.Yes:
            first_panel = (
                next(iter(self.panels.values()), None)
                or next(
                    (p._inner_panel for p in self._multicolor_panels.values()),
                    None,
                )
            )
            if first_panel is not None:
                first_panel._switch_to_viz_and_show_tracks()
