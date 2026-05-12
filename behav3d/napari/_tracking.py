"""
BEHAV3D napari plugin – Tracking Tab.

Provides per-cell-type sub-tabs with method selection (LAP / TrackPy / Propagation),
method-specific parameters, and batch-tracking options.
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
    QLineEdit, QFileDialog,
)
from qtpy.QtCore import Qt

from behav3d.napari._widgets import make_help_row


# ═══════════════════════════════════════════════════════════════════════════
# _ImportTrackingPage — Page 4 of CellTypeTrackingPanel
# ═══════════════════════════════════════════════════════════════════════════
class _ImportTrackingPage(QWidget):
    """Per-cell-type widget for importing pre-tracked zarr/tiff files.

    Reads the source path from the metadata column
    ``{prefix}_{cell_type}_tracks_image_path``, validates/converts it,
    and writes outputs to the standard BEHAV3D output locations — exactly
    the same paths used by every other tracking algorithm.
    """

    def __init__(self, cell_type: str, category: str, metadata_loader, parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category
        self.metadata_loader = metadata_loader
        self._prefix = {"organoid": "or", "immune": "im"}.get(category, "ot")
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
        if not path_str:
            return None
        p = Path(path_str)
        if p.exists():
            return p
        md_csv = (self.metadata_loader.behav3d_parameters
                  .get("paths", {}).get("metadata_csv"))
        if md_csv:
            p_rel = Path(md_csv).parent / path_str
            if p_rel.exists():
                return p_rel
        return p  # return non-existent path so caller can report it

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
            f"Set the tracking image path in your metadata CSV "
            f"to a pre-tracked image (.zarr or .tif/.tiff). "
            f"Column: {self._tracks_img_col()}"
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
        """Add one sample row. Returns True if an action button was added."""
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

        row_w = QWidget()
        row_lay = QHBoxLayout(row_w)
        row_lay.setContentsMargins(16, 2, 4, 4)
        needs_action = False

        if not has_value:
            lbl = QLabel("No source path set in metadata")
            lbl.setStyleSheet("color:#999; font-style:italic;")
            row_lay.addWidget(lbl)

        else:
            path_str = str(raw_val).strip().strip('"').strip("'")
            file_path = self._resolve_path(path_str)

            if not file_path.exists():
                lbl = QLabel(f"⚠️  File not found: {file_path}")
                lbl.setWordWrap(True)
                lbl.setStyleSheet("color:#E65100;")
                row_lay.addWidget(lbl)

            elif file_path.suffix.lower() in (".tif", ".tiff"):
                dest_z = self._dest_zarr(sample_name)
                dest_c = self._dest_csv(sample_name)
                if dest_z.exists() and dest_c.exists():
                    lbl = QLabel("✅  Already processed")
                    lbl.setStyleSheet("color:#2E7D32; font-weight:bold;")
                    row_lay.addWidget(lbl)
                    btn_regen = QPushButton("🔄  Re-process")
                    btn_regen.setToolTip(f"Will overwrite:\n  {dest_z}\n  {dest_c}")
                    btn_regen.setStyleSheet(
                        "QPushButton{background:#546E7A;color:white;padding:4px 8px;"
                        "border-radius:3px;font-size:10px}"
                        "QPushButton:hover{background:#607D8B}"
                    )
                    btn_regen.clicked.connect(
                        partial(self._process_single, path_str, sample_name, row_idx)
                    )
                    row_lay.addWidget(btn_regen)
                else:
                    btn = QPushButton("🔄  Convert TIFF → zarr")
                    btn.setToolTip(f"Convert TIFF to zarr + generate CSV tracking data.\nWill write:\n  {dest_z}\n  {dest_c}")
                    btn.setStyleSheet(
                        "QPushButton{background:#1565C0;color:white;padding:4px 10px;border-radius:3px}"
                        "QPushButton:hover{background:#1976D2}"
                    )
                    btn.clicked.connect(partial(self._process_single, path_str, sample_name, row_idx))
                    row_lay.addWidget(btn)
                    needs_action = True

            elif file_path.suffix == ".zarr" or file_path.is_dir():
                dest_z = self._dest_zarr(sample_name)
                dest_c = self._dest_csv(sample_name)

                if dest_z.exists() and dest_c.exists():
                    lbl = QLabel("✅  Already processed")
                    lbl.setStyleSheet("color:#2E7D32; font-weight:bold;")
                    row_lay.addWidget(lbl)
                    btn_regen = QPushButton("🔄  Re-process")
                    btn_regen.setToolTip(
                        f"Will overwrite:\n  {dest_z}\n  {dest_c}"
                    )
                    btn_regen.setStyleSheet(
                        "QPushButton{background:#546E7A;color:white;padding:4px 8px;"
                        "border-radius:3px;font-size:10px}"
                        "QPushButton:hover{background:#607D8B}"
                    )
                    btn_regen.clicked.connect(
                        partial(self._process_single, path_str, sample_name, row_idx)
                    )
                    row_lay.addWidget(btn_regen)
                    # Re-process is available but doesn't count for "Process All"
                else:
                    btn = QPushButton("📄  Import zarr")
                    btn.setToolTip(
                        f"Import zarr + generate CSV tracking data.\nWill write:\n  {dest_z}\n  {dest_c}"
                    )
                    btn.setStyleSheet(
                        "QPushButton{background:#2E7D32;color:white;padding:4px 10px;border-radius:3px}"
                        "QPushButton:hover{background:#388E3C}"
                    )
                    btn.clicked.connect(
                        partial(self._process_single, path_str, sample_name, row_idx)
                    )
                    row_lay.addWidget(btn)
                    needs_action = True
            else:
                lbl = QLabel(f"⚠️  Unsupported format ({file_path.suffix})")
                lbl.setStyleSheet("color:#E65100;")
                row_lay.addWidget(lbl)

        row_lay.addStretch()
        self.scroll_layout.addWidget(row_w)

        sep = QWidget()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background:#ddd; margin:2px 0;")
        self.scroll_layout.addWidget(sep)
        return needs_action

    # ── processing ──────────────────────────────────────────────────────
    def _process_single(self, src_path_str: str, sample_name: str,
                        row_idx: int, _=None, save: bool = True):
        import shutil
        from qtpy.QtWidgets import QMessageBox
        from behav3d.io.images import load_image, save_as_zarr
        from behav3d.preprocessing.tracking import convert_tracked_image_to_csv

        src = self._resolve_path(src_path_str)
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

            if src.suffix.lower() in (".tif", ".tiff") or not zarr_ok:
                print(f"  Converting {src.name} → {dest_z.name} …", flush=True)
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

            # Step 2 — generate tracks CSV
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

            if save:
                self._save_metadata(md)

        except Exception as exc:
            traceback.print_exc()
            print(f"  ❌ Error processing {sample_name}: {exc}", flush=True)

    def _process_all(self, _=None):
        col = self._tracks_img_col()
        md = self.metadata_loader.metadata
        if md is None:
            return
        for idx, row in md.iterrows():
            raw_val = row.get(col) if col in row.index else None
            if raw_val is None or not pd.notna(raw_val) or str(raw_val).strip() in ("", "nan"):
                continue
            path_str = str(raw_val).strip().strip('"').strip("'")
            fp = self._resolve_path(path_str)
            if fp is None or not fp.exists():
                continue
            sample_name = str(row.get("sample_name", f"Row {idx + 1}"))
            self._process_single(path_str, sample_name, int(idx), save=False)
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
                 log_callback=None, viewer=None, parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category          # "organoid" / "immune" / "other"
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types   # all types in same category
        self.log = log_callback or print
        self.viewer = viewer
        self._toggle_all_organoids_callback = None  # set by TrackingTab for organoid panels

        # Determine defaults
        if category == "organoid":
            def_method = "propagation"
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
        self.combo_method.addItems([
            "LAP (laptrack)", "TrackPy", "Propagation",
            "btrack (Bayesian)", "Import tracking",
        ])
        saved_method = tcfg.get("method", def_method)
        idx_map = {"lap": 0, "trackpy": 1, "propagation": 2, "btrack": 3, "import": 4}
        self.combo_method.setCurrentIndex(idx_map.get(saved_method, 0))
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.combo_method)
        method_group.setLayout(method_layout)

        layout.addWidget(method_group)

        # ── Stacked method params ────────────────────────────────────
        self.param_stack = QStackedWidget()
        self.combo_method.currentIndexChanged.connect(self.param_stack.setCurrentIndex)

        # Page 0 — LAP params
        lap_cfg = tcfg.get("lap", {})
        lap_page = QWidget()
        lap_form = QFormLayout(lap_page)
        lap_form.setContentsMargins(6, 4, 6, 4)
        lap_form.setSpacing(3)
        lap_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.lap_track_cost = QSpinBox()
        self.lap_track_cost.setRange(1, 999)
        self.lap_track_cost.setValue(int(lap_cfg.get("track_cost_px", 45)))
        self.lap_track_cost.setMaximumWidth(80)
        lap_form.addRow("Track cost (px):", make_help_row(
            self.lap_track_cost,
            "Track Cost (pixels)",
            "Maximum pixel distance a cell can travel between two "
            "consecutive frames to be linked as the same track.\n\n"
            "Increase if cells move fast; decrease to avoid false links."
        ))

        self.lap_gap_cost = QSpinBox()
        self.lap_gap_cost.setRange(1, 999)
        self.lap_gap_cost.setValue(int(lap_cfg.get("gap_close_cost_px", 60)))
        self.lap_gap_cost.setMaximumWidth(80)
        lap_form.addRow("Gap close cost (px):", make_help_row(
            self.lap_gap_cost,
            "Gap Closing Cost (pixels)",
            "Maximum distance (in pixels) allowed when reconnecting a "
            "track that was temporarily lost for one or more frames.\n\n"
            "Should be >= Track cost. Increase if cells disappear "
            "briefly due to segmentation gaps."
        ))

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
        self.lap_merge_cost.setRange(0, 999)
        self.lap_merge_cost.setValue(int(lap_cfg.get("merging_cost_px", 0)))
        self.lap_merge_cost.setMaximumWidth(80)
        lap_form.addRow("Merging cost (px):", make_help_row(
            self.lap_merge_cost,
            "Merging Cost (pixels)",
            "Maximum distance for detecting merge events, where two "
            "tracks combine into one object.\n\n"
            "Set to 0 to disable merging detection.\n"
            "Useful when cells fuse or cluster together."
        ))

        self.lap_split_cost = QSpinBox()
        self.lap_split_cost.setRange(0, 999)
        self.lap_split_cost.setValue(int(lap_cfg.get("splitting_cost_px", 0)))
        self.lap_split_cost.setMaximumWidth(80)
        lap_form.addRow("Splitting cost (px):", make_help_row(
            self.lap_split_cost,
            "Splitting Cost (pixels)",
            "Maximum distance for detecting split events, where one "
            "object divides into two tracks.\n\n"
            "Set to 0 to disable splitting detection.\n"
            "Useful for cell division or organoid fragmentation."
        ))

        self.param_stack.addWidget(lap_page)

        # Page 1 — TrackPy params
        tp_cfg = tcfg.get("trackpy", {})
        tp_page = QWidget()
        tp_form = QFormLayout(tp_page)
        tp_form.setContentsMargins(6, 4, 6, 4)
        tp_form.setSpacing(3)
        tp_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.tp_search_range = QSpinBox()
        self.tp_search_range.setRange(1, 999)
        self.tp_search_range.setValue(int(tp_cfg.get("search_range_px", 31)))
        self.tp_search_range.setMaximumWidth(80)
        tp_form.addRow("Search range (px):", make_help_row(
            self.tp_search_range,
            "Search Range (pixels)",
            "Maximum pixel distance to look for a cell in the next frame.\n\n"
            "Should be large enough to cover the fastest-moving cells."
        ))

        self.tp_memory = QSpinBox()
        self.tp_memory.setRange(0, 100)
        self.tp_memory.setValue(int(tp_cfg.get("memory_frames", 2)))
        self.tp_memory.setMaximumWidth(60)
        tp_form.addRow("Memory (frames):", make_help_row(
            self.tp_memory,
            "Memory (frames)",
            "Number of frames a cell can disappear and still be "
            "reconnected to its previous track.\n\n"
            "Similar to 'gap closing' in LAP."
        ))

        self.tp_adaptive_stop = QDoubleSpinBox()
        self.tp_adaptive_stop.setRange(0.0, 100.0)
        self.tp_adaptive_stop.setSingleStep(0.5)
        self.tp_adaptive_stop.setValue(float(tp_cfg.get("adaptive_stop", 10.0)))
        self.tp_adaptive_stop.setMaximumWidth(80)
        tp_form.addRow("Adaptive stop:", make_help_row(
            self.tp_adaptive_stop,
            "Adaptive Stop",
            "Factor that limits how much the search range can shrink "
            "adaptively.\n\n"
            "Higher = more conservative shrinking.\n"
            "Leave at default unless tracking quality is poor."
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

        self.param_stack.addWidget(tp_page)

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
                "using Propagation, collapsing all organoid tabs into one.\n"
                "This is the default behaviour when multiple organoid types exist."
            )
            def _on_check_all(checked, _self=self):
                if _self._toggle_all_organoids_callback:
                    _self._toggle_all_organoids_callback(checked)
            self.check_all_together_prop.toggled.connect(_on_check_all)
            prop_lay.addWidget(self.check_all_together_prop)

        prop_notice = QLabel("No tunable parameters for\nPropagation tracking method.")
        prop_notice.setWordWrap(True)
        prop_notice.setAlignment(Qt.AlignCenter)
        prop_notice.setStyleSheet("color: #666; font-style: italic; padding: 10px;")
        prop_lay.addWidget(prop_notice)
        prop_lay.addStretch()

        self.param_stack.addWidget(prop_page)

        # Page 3 — btrack (Bayesian tracking)
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
        self.bt_max_search_radius.setRange(1, 9999)
        self.bt_max_search_radius.setValue(int(bt_cfg.get("max_search_radius", 100)))
        self.bt_max_search_radius.setMaximumWidth(80)
        step1_form.addRow("Max search radius (px):", make_help_row(
            self.bt_max_search_radius,
            "Max Search Radius (pixels)",
            "Maximum isotropic distance (pixels) to search for\n"
            "linking objects between frames.\n\n"
            "Increase for fast-moving cells; decrease to\n"
            "prevent long-range false links."
        ))

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
        step1_form.addRow("", make_help_row(
            self.bt_use_visual_features,
            "Visual Features",
            "When enabled, raw image intensity statistics (mean, std per channel)\n"
            "are computed alongside centroids and used by the Kalman filter for\n"
            "more accurate linking.\n\n"
            "Requires raw image data (raw_image_path) in metadata."
        ))

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
        hyp_lay = QVBoxLayout(hyp_group)
        hyp_lay.setContentsMargins(4, 2, 4, 2)
        hyp_lay.setSpacing(1)
        saved_hyps = bt_cfg.get("hypotheses",
                                ["P_FP", "P_init", "P_term", "P_link"])
        self.bt_hyp_checks = {}
        for hyp_name, hyp_desc, default_on in [
            ("P_FP",     "False positive",    True),
            ("P_init",   "Track initialization", True),
            ("P_term",   "Track termination",  True),
            ("P_link",   "Track linking",      True),
            ("P_branch", "Track branching",    False),
            ("P_dead",   "Cell death",         False),
            ("P_merge",  "Track merging",      False),
        ]:
            cb = QCheckBox(f"{hyp_name} — {hyp_desc}")
            is_on = hyp_name in saved_hyps if saved_hyps else default_on
            cb.setChecked(is_on)
            if hyp_name == "P_FP":
                cb.setEnabled(False)  # P_FP is always required
            hyp_lay.addWidget(cb)
            self.bt_hyp_checks[hyp_name] = cb
        step2_lay.addWidget(hyp_group)

        step2_form = QFormLayout()
        step2_form.setContentsMargins(0, 0, 0, 0)
        step2_form.setSpacing(3)
        step2_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.bt_dist_thresh = QSpinBox()
        self.bt_dist_thresh.setRange(1, 9999)
        self.bt_dist_thresh.setValue(int(bt_cfg.get("dist_thresh", 60)))
        self.bt_dist_thresh.setMaximumWidth(80)
        step2_form.addRow("Distance threshold:", make_help_row(
            self.bt_dist_thresh,
            "Distance Threshold",
            "Maximum distance (pixels) for generating\n"
            "link/branch hypotheses in the optimizer."
        ))

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
        self.param_stack.addWidget(btrack_page)

        # Page 4 — Import tracking
        import_page = _ImportTrackingPage(
            cell_type=self.cell_type,
            category=self.category,
            metadata_loader=self.metadata_loader,
        )
        self.param_stack.addWidget(import_page)

        # Set active page
        self.param_stack.setCurrentIndex(self.combo_method.currentIndex())
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

        # Disable run button for coming-soon methods (only Import at index 4)
        def _on_method_idx_changed(idx):
            is_coming_soon = idx >= 4
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
        return ["lap", "trackpy", "propagation", "btrack", "import"][
            self.combo_method.currentIndex()
        ]

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
                "track_cost_px": int(self.lap_track_cost.value()),
                "gap_close_cost_px": int(self.lap_gap_cost.value()),
                "gap_close_max_frames": int(self.lap_gap_frames.value()),
                "merging_cost_px": int(self.lap_merge_cost.value()),
                "splitting_cost_px": int(self.lap_split_cost.value()),
            },
            "trackpy": {
                "search_range_px": int(self.tp_search_range.value()),
                "memory_frames": int(self.tp_memory.value()),
                "adaptive_stop": float(self.tp_adaptive_stop.value()),
                "adaptive_step": float(self.tp_adaptive_step.value()),
            },
            "propagation": {
                # Notice: no tunable params currently exposed
            },
            "btrack": {
                "config_preset": self._bt_get_config_preset(),
                "config_path": self.bt_config_path.text().strip(),
                "use_visual_features": self.bt_use_visual_features.isChecked(),
                "max_search_radius": int(self.bt_max_search_radius.value()),
                "update_method": "APPROXIMATE" if self.bt_update_method.currentIndex() == 1 else "EXACT",
                "step_size": int(self.bt_step_size.value()),
                "n_workers": max(1, int(self.bt_n_workers.value())),
                "use_optimize": self.bt_use_optimize.isChecked(),
                "hypotheses": self._bt_get_hypotheses(),
                "dist_thresh": int(self.bt_dist_thresh.value()),
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
                idx_map = {"lap": 0, "trackpy": 1, "propagation": 2, "btrack": 3}
                panel.combo_method.setCurrentIndex(idx_map.get(settings["method"], 0))
                
                # LAP
                panel.lap_track_cost.setValue(settings["lap"]["track_cost_px"])
                panel.lap_gap_cost.setValue(settings["lap"]["gap_close_cost_px"])
                panel.lap_gap_frames.setValue(settings["lap"]["gap_close_max_frames"])
                panel.lap_merge_cost.setValue(settings["lap"]["merging_cost_px"])
                panel.lap_split_cost.setValue(settings["lap"]["splitting_cost_px"])
                
                # TrackPy
                panel.tp_search_range.setValue(settings["trackpy"]["search_range_px"])
                panel.tp_memory.setValue(settings["trackpy"]["memory_frames"])
                panel.tp_adaptive_stop.setValue(settings["trackpy"]["adaptive_stop"])
                panel.tp_adaptive_step.setValue(settings["trackpy"]["adaptive_step"])
                
                # btrack
                bt = settings.get("btrack", {})
                preset = bt.get("config_preset", "cell")
                preset_idx = {"cell": 0, "particle": 1}.get(
                    str(preset).lower(), 2
                )
                panel.bt_config_preset.setCurrentIndex(preset_idx)
                panel.bt_config_path.setText(bt.get("config_path", ""))
                panel.bt_use_visual_features.setChecked(bt.get("use_visual_features", False))
                panel.bt_max_search_radius.setValue(bt.get("max_search_radius", 100))
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
                panel.bt_dist_thresh.setValue(bt.get("dist_thresh", 60))
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
    def _determine_targets(self) -> list:
        """Which cell types to track based on batch checkboxes."""
        if self.check_batch_all.isChecked():
            return list(self.all_cell_types)
        if self.check_batch_category.isChecked():
            return list(self.category_types)
        return [self.cell_type]

    def _run_tracking_for(self, cell_type: str, overwrite: bool = False):
        """Run tracking for a single cell type using current panel settings."""
        method = self._get_method_key()
        metadata = self.metadata_loader.metadata
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        print(f"\n{'='*50}", file=sys.stderr)
        print(f"  Tracking: {cell_type} ({method.upper()})", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        if method == "lap":
            from behav3d.preprocessing.tracking.laptracking import run_laptracking
            tc = int(self.lap_track_cost.value()) ** 2
            gc = int(self.lap_gap_cost.value()) ** 2
            mc = int(self.lap_merge_cost.value())
            sc = int(self.lap_split_cost.value())
            new_md = run_laptracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                track_cost_cutoff=tc,
                gap_closing_cost_cutoff=gc,
                gap_closing_max_frame_count=int(self.lap_gap_frames.value()),
                merging_cost_cutoff=(mc ** 2 if mc > 0 else False),
                splitting_cost_cutoff=(sc ** 2 if sc > 0 else False),
                overwrite=overwrite,
            )

        elif method == "trackpy":
            from behav3d.preprocessing.tracking.trackpy_tracking import run_trackpy_tracking_generic
            new_md = run_trackpy_tracking_generic(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                search_range=int(self.tp_search_range.value()),
                memory=int(self.tp_memory.value()),
                adaptive_stop=float(self.tp_adaptive_stop.value()),
                adaptive_step=float(self.tp_adaptive_step.value()),
                log_callback=self.log,
            )

        elif method == "btrack":
            from behav3d.preprocessing.tracking.btrack_tracking import run_btracking
            config_preset = self._bt_get_config_preset()
            update_method = "APPROXIMATE" if self.bt_update_method.currentIndex() == 1 else "EXACT"
            use_opt = self.bt_use_optimize.isChecked()
            new_md = run_btracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                config_preset=config_preset,
                use_visual_features=bool(self.bt_use_visual_features.isChecked()),
                max_search_radius=int(self.bt_max_search_radius.value()),
                update_method=update_method,
                step_size=int(self.bt_step_size.value()),
                n_workers=max(1, int(self.bt_n_workers.value())),
                use_optimize=use_opt,
                hypotheses=self._bt_get_hypotheses() if use_opt else None,
                dist_thresh=int(self.bt_dist_thresh.value()) if use_opt else None,
                time_thresh=int(self.bt_time_thresh.value()) if use_opt else None,
                log_callback=self.log,
            )

        else:  # propagation
            from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
            new_md = run_propagation_tracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite
            )

        print(f"\u2705 {cell_type} tracking finished.", file=sys.stderr)

        # Update shared metadata
        self.metadata_loader.metadata = new_md
        # Save updated metadata CSV
        csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
        if csv_path:
            new_md.to_csv(csv_path, sep=",", index=False)

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
        self._persist()
        method = self._get_method_key()
        self.log(f"Running {method.upper()} tracking for: {self.cell_type}")

        # Check for existing data across all samples
        existing = self._check_existing_tracking([self.cell_type])

        overwrite = True
        if existing:
            from qtpy.QtWidgets import QMessageBox
            details = "\n".join(f"  \u2022 {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Tracking?")
            box.setText(
                f"The following tracking data already exists:\n\n{details}\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked != btn_overwrite:
                self.log(f"Tracking for {self.cell_type} cancelled.")
                return
            overwrite = True

        try:
            self._run_tracking_for(self.cell_type, overwrite=overwrite)
            self.log(f"\u2705 {self.cell_type} tracking finished.")
            
            # Show visualizer prompt
            from qtpy.QtWidgets import QMessageBox
            res = QMessageBox.question(
                self, "Tracking Finished",
                f"Tracking for {self.cell_type} finished! \n\nDo you want to switch to the Visualization Tab and see the tracks?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.Yes:
                self._switch_to_viz_and_show_tracks()

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during tracking: {e}")

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
                
                # Make 'Tracks' layers visible
                for layer in self.viewer.layers:
                    if "Tracks" in layer.name:
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
                 viewer, toggle_callback, parent=None):
        super().__init__(parent)
        self.organoid_types = list(organoid_types)
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.viewer = viewer
        self._toggle_callback = toggle_callback
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
        method_label = QLabel("Propagation")
        method_label.setStyleSheet("font-weight: bold;")
        method_layout.addWidget(method_label)
        method_layout.addStretch()
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # ── Warning banner ────────────────────────────────────────────
        warning = QLabel(
            "⚠️  Propagation tracking will be applied to ALL organoid types simultaneously. "
            "To run individual tracking, uncheck 'Track all organoids together'"
        )
        warning.setWordWrap(True)
        warning.setStyleSheet(
            "color: #e65100; background: #fff3e0; border: 1px solid #ffcc80; "
            "border-radius: 4px; padding: 8px; font-size: 11px; margin: 2px 0;"
        )
        layout.addWidget(warning)

        # ── Propagation parameters group ──────────────────────────────
        prop_group = QGroupBox("Propagation Parameters")
        prop_lay = QVBoxLayout(prop_group)
        prop_lay.setSpacing(6)

        self.check_all_together = QCheckBox("Track all organoids together")
        self.check_all_together.setChecked(True)
        self.check_all_together.setToolTip(
            "When checked, all organoid types are tracked simultaneously using Propagation.\n"
            "Uncheck to configure and run tracking independently per organoid type."
        )
        self.check_all_together.toggled.connect(self._on_toggled)
        prop_lay.addWidget(self.check_all_together)

        notice = QLabel("No additional tunable parameters for Propagation tracking.")
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
        """Run propagation tracking for all organoid types at once."""
        from qtpy.QtWidgets import QMessageBox
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
            details = "\n".join(f"  \u2022 {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Tracking?")
            box.setText(
                f"The following tracking data already exists:\n\n{details}\n\nWhat do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite", QMessageBox.DestructiveRole)
            box.addButton("Skip", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == btn_cancel:
                self.log("All-organoids tracking cancelled.")
                return
            overwrite = (clicked == btn_overwrite)

        types_str = ", ".join(self.organoid_types)
        self.btn_run.setEnabled(False)
        self.btn_run.setText("\u23f3 Running\u2026")
        try:
            self.log(f"\u25b6 Propagation tracking \u2014 all organoids: {types_str}\u2026")
            from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
            new_md = run_propagation_tracking(
                metadata=md,
                output_dir=str(out_path.expanduser()),
                cell_type=self.organoid_types[0],
                overwrite=overwrite,
                all_organoids=True,
            )
            self.metadata_loader.metadata = new_md
            csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
            if csv_path:
                new_md.to_csv(csv_path, sep=",", index=False)
            self.log("\u2705 All organoids propagation tracking finished.")

            res = QMessageBox.question(
                self, "Tracking Finished",
                "All organoid types tracked successfully!\n\n"
                "Do you want to switch to the Visualization Tab to see the tracks?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res == QMessageBox.Yes:
                self._switch_to_viz()
        except Exception as e:
            traceback.print_exc()
            self.log(f"\u274c Error during all-organoids tracking: {e}")
        finally:
            self.btn_run.setEnabled(True)
            self.btn_run.setText(f"\u25b6  Run All Organoids Tracking  ({types_str})")

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
                    if "Tracks" in layer.name:
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
                 metadata_loader, log_callback=None, viewer=None, parent=None):
        super().__init__(parent)
        self.base_name = base_name
        self.channel_types = sorted(channel_types)
        self.category = category
        self.metadata_loader = metadata_loader
        self.log = log_callback or print
        self.viewer = viewer

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
        )

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

    def run_tracking(self, overwrite: bool = False):
        """Track each channel then merge into <base_name>_merged."""
        n = len(self.channel_types)
        merged_name = f"{self.base_name}_merged"

        for i, ch in enumerate(self.channel_types, 1):
            self.log(f"  Multicolor [{i}/{n}] Tracking channel: {ch}…")
            print(f"\n  Multicolor [{i}/{n}] Tracking: {ch}", file=sys.stderr)
            self._inner_panel._run_tracking_for(ch, overwrite=overwrite)
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


# ═══════════════════════════════════════════════════════════════════════════
# TrackingTab — main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class TrackingTab(QWidget):
    def __init__(self, viewer, metadata_loader, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeTrackingPanel] = {}
        self._multicolor_panels: dict[str, MulticolorTrackingPanel] = {}

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

    # ------------------------------------------------------------------
    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing tracking tabs…")
        self._rebuild_tabs()

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
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🔵 {ct.capitalize()}")

        # ── Immune multicolor groups ──────────────────────────────────
        for base, channels in imm_mc.items():
            panel = MulticolorTrackingPanel(
                base_name=base, channel_types=channels, category="immune",
                metadata_loader=self.metadata_loader,
                log_callback=self._log, viewer=self.viewer,
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
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🟡 {ct.capitalize()}")

        # ── Other multicolor groups ───────────────────────────────────
        for base, channels in oth_mc.items():
            panel = MulticolorTrackingPanel(
                base_name=base, channel_types=channels, category="other",
                metadata_loader=self.metadata_loader,
                log_callback=self._log, viewer=self.viewer,
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
        self.run_batch_tracking(interactive=True)

    def run_batch_tracking(self, interactive=True, skip_existing=False):
        """Sequential run for all configured cell type panels.
        When interactive=False, skips overwrite and visualization dialogs."""
        if not self.panels and not self._multicolor_panels and not self._all_organoids_panel:
            self._log("No cell type panels to track.")
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
            from qtpy.QtWidgets import QMessageBox
            details = "\n".join(f"  \u2022 {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Tracking?")
            box.setText(
                f"The following tracking data already exists:\n\n{details}\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == btn_cancel:
                self._log("Batch tracking cancelled by user.")
                return
            skip_existing_flag = (clicked == btn_skip)
            overwrite = not skip_existing_flag

        # Persist all panel params before starting so config is saved even if an error occurs
        for ct, panel in self.panels.items():
            panel._persist()
        for base, mc_panel in self._multicolor_panels.items():
            mc_panel._persist()

        try:
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
                    )
                    self.metadata_loader.metadata = new_md
                    csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
                    if csv_path:
                        new_md.to_csv(csv_path, sep=",", index=False)

                    self._log(f"Done: all organoids propagation tracking")
                    organoids_done = set(self._organoid_types)
                else:
                    self._log("All organoid tracking data already exists — skipping.")
                    organoids_done = set(self._organoid_types)

            # ── Per-cell-type tracking (non-organoids or individual organoids) ──
            step = len(organoids_done)
            for i, (ct, panel) in enumerate(self.panels.items(), 1):
                if ct in organoids_done:
                    continue
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- [{i}/{total}] Skipping {ct} (existing data) ---")
                    continue
                step += 1
                print(f"\n▶ [{step}/{total}] Tracking {ct}...", file=sys.stderr)
                self._log(f"--- [{step}/{total}] Tracking {ct} ---")
                panel._run_tracking_for(ct, overwrite=overwrite)

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

            # ── Multicolor groups (track channels + auto-merge) ───────
            for base, mc_panel in self._multicolor_panels.items():
                step += 1
                merged_name = f"{base}_merged"
                print(f"\n▶ [{step}/{total}] Multicolor tracking: {base} "
                      f"({len(mc_panel.channel_types)} channels + merge)…",
                      file=sys.stderr)
                self._log(f"--- [{step}/{total}] Multicolor tracking: {base} "
                          f"({', '.join(mc_panel.channel_types)}) → {merged_name} ---")
                mc_panel.run_tracking(overwrite=overwrite)
                self._log(f"Done: {base} (channels tracked + merged → {merged_name})")

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

            if interactive:
                from qtpy.QtWidgets import QMessageBox
                res = QMessageBox.question(
                    self, "Batch Tracking Finished",
                    "Successfully tracked all cell types! \n\nDo you want to switch to the Visualization Tab and see the tracks?",
                    QMessageBox.Yes | QMessageBox.No
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

        except Exception as e:
            traceback.print_exc()
            self._log(f"\u274c Batch tracking error: {e}")
