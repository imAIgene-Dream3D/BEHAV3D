"""
BEHAV3D napari plugin – Backprojection Tab.

A top-level main tab that lets the user:

- Select a cell type (immune / other, non-multicolor)
- Select a sample (from loaded metadata)
- View backprojected behavioral states as napari image layers
- View backprojected track clusters as napari image layers
- Export backprojection PDFs / MP4s

Works on existing clustering outputs: the adata from state classification
(``behavioral_states/BEHAV3D_<ct>_behavioral_states.h5ad``) and track
clustering (``behavorial_trajectories/BEHAV3D_<ct>_track_trajectories.h5ad``)
must already be present on disk.  If not, the buttons will be disabled with a
helpful tooltip.
"""
from __future__ import annotations

import datetime
import traceback
from pathlib import Path
from typing import Optional

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
    QStackedWidget,
)

from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
)
from behav3d.core.qt_help import make_help_row


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _detect_sc_cell_types(metadata_loader) -> list[str]:
    if metadata_loader is None or getattr(metadata_loader, "metadata", None) is None:
        return []
    try:
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            filter_multicolor_inputs,
        )
        md = metadata_loader.metadata
        imm = list(filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)))
        oth = list(filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)))
        return imm + oth
    except Exception:
        traceback.print_exc()
        return []


def _style_primary(btn: QPushButton):
    btn.setStyleSheet(
        "QPushButton { background: #28a745; color: white; font-weight: bold; "
        "border-radius: 4px; padding: 6px 10px; }"
        "QPushButton:hover { background: #218838; }"
        "QPushButton:disabled { background: #3a3a3a; color: #666; }"
    )


def _style_secondary(btn: QPushButton):
    btn.setStyleSheet(
        "QPushButton { background: #2d5a8e; color: #ddd; font-weight: bold; "
        "border-radius: 4px; padding: 6px 10px; }"
        "QPushButton:hover { background: #3a72b3; }"
        "QPushButton:disabled { background: #3a3a3a; color: #666; }"
    )


# ═══════════════════════════════════════════════════════════════════════════
# State Backprojection inner tab
# ═══════════════════════════════════════════════════════════════════════════

class StateBackprojectionSubTab(QWidget):
    """Visualise behavioral states overlaid on raw images."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None,
                 sample_getter=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._get_cell_type = cell_type_getter
        self._get_sample = sample_getter
        self._bg = BackgroundOperation(self)
        self._init_ui()

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        scroll.setWidget(content)

        # ── Napari layer viewer ──────────────────────────────────────
        grp_view = QGroupBox("Live Napari Layer Backprojection")
        g_view = QVBoxLayout(grp_view)
        g_view.setSpacing(4)

        info = QLabel(
            "Loads the state adata and overlays colored state labels onto "
            "the selected sample image in napari."
        )
        info.setStyleSheet("color: #999; font-size: 10px;")
        info.setWordWrap(True)
        g_view.addWidget(info)

        form = QFormLayout()
        form.setSpacing(3)

        self.combo_color_by = QComboBox()
        self.combo_color_by.addItems(
            ["full_behavioral_cluster", "intrinsic_behavioral_cluster"]
        )
        self.combo_color_by.setMinimumWidth(220)
        form.addRow("Color by:", make_help_row(
            self.combo_color_by, "Color by", "Which clustering label to use for coloring the backprojection."
        ))

        self.spin_opacity = QSpinBox()
        self.spin_opacity.setRange(10, 100)
        self.spin_opacity.setValue(80)
        self.spin_opacity.setSuffix(" %")
        self.spin_opacity.setMaximumWidth(90)
        form.addRow("Opacity:", make_help_row(
            self.spin_opacity, "Opacity", "Opacity of the colored state overlay."
        ))

        g_view.addLayout(form)

        view_row = QHBoxLayout()
        self.btn_show_state = QPushButton("▶ Show State Backprojection in Napari")
        _style_primary(self.btn_show_state)
        self.btn_show_state.clicked.connect(self._on_show_state)
        view_row.addWidget(self.btn_show_state, stretch=1)
        g_view.addLayout(view_row)

        lay.addWidget(grp_view)

        # ── PDF / MP4 Export ────────────────────────────────────────
        grp_export = QGroupBox("Export Backprojection (PDF / MP4)")
        g_export = QVBoxLayout(grp_export)
        g_export.setSpacing(4)

        export_form = QFormLayout()
        export_form.setSpacing(3)

        self.spin_export_dpi = QSpinBox()
        self.spin_export_dpi.setRange(50, 600)
        self.spin_export_dpi.setValue(150)
        self.spin_export_dpi.setSuffix(" dpi")
        self.spin_export_dpi.setMaximumWidth(110)
        export_form.addRow("DPI (PDF):", self.spin_export_dpi)
        g_export.addLayout(export_form)

        export_opts = QHBoxLayout()
        self.chk_export_pdf = QCheckBox("PDF")
        self.chk_export_pdf.setChecked(True)
        export_opts.addWidget(self.chk_export_pdf)
        self.chk_export_mp4 = QCheckBox("MP4")
        self.chk_export_mp4.setChecked(False)
        export_opts.addWidget(self.chk_export_mp4)
        export_opts.addStretch()
        g_export.addLayout(export_opts)

        export_row = QHBoxLayout()
        self.btn_export_state = QPushButton("▶ Export State Backprojection")
        _style_secondary(self.btn_export_state)
        self.btn_export_state.clicked.connect(self._on_export_state)
        export_row.addWidget(self.btn_export_state, stretch=1)
        g_export.addLayout(export_row)

        lay.addWidget(grp_export)

        # ── Progress + Log ───────────────────────────────────────────
        self.progress_row = ProgressBarRow()
        lay.addWidget(self.progress_row)
        lay.addWidget(QLabel("Log"))
        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumHeight(120)
        self.log_edit.setStyleSheet(
            "background: #111; color: #ccc; font-family: monospace; font-size: 10px;"
        )
        lay.addWidget(self.log_edit)
        lay.addStretch(1)

    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_edit.append(f"[{ts}] {msg}")

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

    def _sample(self) -> str:
        return self._get_sample() if self._get_sample else ""

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _state_adata_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return out / "analysis" / ct / "behavioral_states" / f"BEHAV3D_{ct}_behavioral_states.h5ad"

    def on_context_updated(self):
        """Called when cell type or sample changes."""
        ct = self._cell_type()
        state_path = self._state_adata_path(ct) if ct else None
        has_state = bool(state_path and state_path.exists())
        self.btn_show_state.setEnabled(has_state)
        if has_state:
            self.btn_show_state.setToolTip(f"Source: {state_path}")
        else:
            tip = "Run State Classification first to generate the adata file."
            if ct:
                expected = str(state_path) if state_path else "—"
                tip += f"\n\nExpected: {expected}"
            self.btn_show_state.setToolTip(tip)
        self.btn_export_state.setEnabled(has_state)

    def _on_show_state(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        state_path = self._state_adata_path(ct)
        if not state_path or not state_path.exists():
            QMessageBox.warning(
                self, "No state adata",
                f"State adata not found:\n{state_path}\n\n"
                "Run State Classification first."
            )
            return
        color_by = self.combo_color_by.currentText()
        opacity = self.spin_opacity.value() / 100.0
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Loading state backprojection for '{ct}' / sample '{sample}'…")

        try:
            from behav3d.napari.backprojection import (
                show_behavioral_state_backprojection,
            )
            show_behavioral_state_backprojection(
                viewer=self.viewer,
                adata_path=str(state_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                opacity=opacity,
                metadata=getattr(self.metadata_loader, "metadata", None),
                output_dir=str(self._out_dir()) if self._out_dir() else "",
            )
            self._log(f"✅ State backprojection loaded.")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Backprojection failed: {e}")

    def _on_export_state(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        state_path = self._state_adata_path(ct)
        if not state_path or not state_path.exists():
            QMessageBox.warning(
                self, "No state adata",
                "Run State Classification first."
            )
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return

        color_by = self.combo_color_by.currentText()
        dpi = int(self.spin_export_dpi.value())
        make_pdf = self.chk_export_pdf.isChecked()
        make_mp4 = self.chk_export_mp4.isChecked()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting state backprojection for '{ct}'…")

        def _run(**kw):
            from behav3d.napari.backprojection import (
                export_behavioral_state_backprojection,
            )
            return export_behavioral_state_backprojection(
                adata_path=str(state_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                output_dir=str(out) if out else "",
                metadata=getattr(self.metadata_loader, "metadata", None),
                dpi=dpi,
                make_pdf=make_pdf,
                make_mp4=make_mp4,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting state backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_state, self.btn_show_state],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: self._log(f"✅ State backprojection export done."),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )


# ═══════════════════════════════════════════════════════════════════════════
# Track Backprojection inner tab
# ═══════════════════════════════════════════════════════════════════════════

class TrackBackprojectionSubTab(QWidget):
    """Visualise track clusters overlaid on raw images."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None,
                 sample_getter=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._get_cell_type = cell_type_getter
        self._get_sample = sample_getter
        self._bg = BackgroundOperation(self)
        self._init_ui()

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        scroll.setWidget(content)

        # ── Napari layer viewer ──────────────────────────────────────
        grp_view = QGroupBox("Live Napari Layer Backprojection (Tracks)")
        g_view = QVBoxLayout(grp_view)
        g_view.setSpacing(4)

        info = QLabel(
            "Loads the track cluster adata and overlays colored trajectory "
            "cluster labels onto the selected sample image in napari."
        )
        info.setStyleSheet("color: #999; font-size: 10px;")
        info.setWordWrap(True)
        g_view.addWidget(info)

        form = QFormLayout()
        form.setSpacing(3)
        self.combo_track_color_by = QComboBox()
        self.combo_track_color_by.addItems(
            ["behavioral_trajectory_cluster", "dtaidistance_cluster", "track_cluster"]
        )
        self.combo_track_color_by.setMinimumWidth(220)
        form.addRow("Color by:", make_help_row(
            self.combo_track_color_by, "Color by", "Which tracking cluster label to use for coloring."
        ))

        self.spin_track_opacity = QSpinBox()
        self.spin_track_opacity.setRange(10, 100)
        self.spin_track_opacity.setValue(80)
        self.spin_track_opacity.setSuffix(" %")
        self.spin_track_opacity.setMaximumWidth(90)
        form.addRow("Opacity:", make_help_row(
            self.spin_track_opacity, "Opacity", "Opacity of the colored track overlay."
        ))

        g_view.addLayout(form)

        view_row = QHBoxLayout()
        self.btn_show_track = QPushButton("▶ Show Track Backprojection in Napari")
        _style_primary(self.btn_show_track)
        self.btn_show_track.clicked.connect(self._on_show_track)
        view_row.addWidget(self.btn_show_track, stretch=1)
        g_view.addLayout(view_row)

        lay.addWidget(grp_view)

        # ── PDF / MP4 Export ────────────────────────────────────────
        grp_export = QGroupBox("Export Track Backprojection (PDF / MP4)")
        g_export = QVBoxLayout(grp_export)
        g_export.setSpacing(4)

        export_form = QFormLayout()
        self.spin_track_export_dpi = QSpinBox()
        self.spin_track_export_dpi.setRange(50, 600)
        self.spin_track_export_dpi.setValue(150)
        self.spin_track_export_dpi.setSuffix(" dpi")
        self.spin_track_export_dpi.setMaximumWidth(110)
        export_form.addRow("DPI (PDF):", self.spin_track_export_dpi)
        g_export.addLayout(export_form)

        opts = QHBoxLayout()
        self.chk_track_pdf = QCheckBox("PDF")
        self.chk_track_pdf.setChecked(True)
        opts.addWidget(self.chk_track_pdf)
        self.chk_track_mp4 = QCheckBox("MP4")
        self.chk_track_mp4.setChecked(False)
        opts.addWidget(self.chk_track_mp4)
        opts.addStretch()
        g_export.addLayout(opts)

        export_row = QHBoxLayout()
        self.btn_export_track = QPushButton("▶ Export Track Backprojection")
        _style_secondary(self.btn_export_track)
        self.btn_export_track.clicked.connect(self._on_export_track)
        export_row.addWidget(self.btn_export_track, stretch=1)
        g_export.addLayout(export_row)

        lay.addWidget(grp_export)

        # Progress + Log
        self.progress_row = ProgressBarRow()
        lay.addWidget(self.progress_row)
        lay.addWidget(QLabel("Log"))
        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumHeight(120)
        self.log_edit.setStyleSheet(
            "background: #111; color: #ccc; font-family: monospace; font-size: 10px;"
        )
        lay.addWidget(self.log_edit)
        lay.addStretch(1)

    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_edit.append(f"[{ts}] {msg}")

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

    def _sample(self) -> str:
        return self._get_sample() if self._get_sample else ""

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _track_adata_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        d = out / "analysis" / ct / "behavorial_trajectories"
        if d.exists():
            for f in sorted(d.glob("*.h5ad")):
                return f
        return d / f"BEHAV3D_{ct}_track_trajectories.h5ad"

    def on_context_updated(self):
        ct = self._cell_type()
        track_path = self._track_adata_path(ct) if ct else None
        has_track = bool(track_path and track_path.exists())
        self.btn_show_track.setEnabled(has_track)
        if has_track:
            self.btn_show_track.setToolTip(f"Source: {track_path}")
        else:
            tip = "Run Track Clustering first to generate the adata file."
            self.btn_show_track.setToolTip(tip)
        self.btn_export_track.setEnabled(has_track)

    def _on_show_track(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        track_path = self._track_adata_path(ct)
        if not track_path or not track_path.exists():
            QMessageBox.warning(
                self, "No track adata",
                "Run Track Clustering first."
            )
            return
        color_by = self.combo_track_color_by.currentText()
        opacity = self.spin_track_opacity.value() / 100.0
        self._log(f"▶ Loading track backprojection for '{ct}' / sample '{sample}'…")

        try:
            from behav3d.napari.backprojection import (
                show_track_cluster_backprojection,
            )
            show_track_cluster_backprojection(
                viewer=self.viewer,
                adata_path=str(track_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                opacity=opacity,
                metadata=getattr(self.metadata_loader, "metadata", None),
                output_dir=str(self._out_dir()) if self._out_dir() else "",
            )
            self._log(f"✅ Track backprojection loaded.")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Backprojection failed: {e}")

    def _on_export_track(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        track_path = self._track_adata_path(ct)
        if not track_path or not track_path.exists():
            QMessageBox.warning(self, "No track adata", "Run Track Clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return

        color_by = self.combo_track_color_by.currentText()
        dpi = int(self.spin_track_export_dpi.value())
        make_pdf = self.chk_track_pdf.isChecked()
        make_mp4 = self.chk_track_mp4.isChecked()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting track backprojection for '{ct}'…")

        def _run(**kw):
            from behav3d.napari.backprojection import (
                export_track_cluster_backprojection,
            )
            return export_track_cluster_backprojection(
                adata_path=str(track_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                output_dir=str(out) if out else "",
                metadata=getattr(self.metadata_loader, "metadata", None),
                dpi=dpi,
                make_pdf=make_pdf,
                make_mp4=make_mp4,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting track backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_track, self.btn_show_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: self._log(f"✅ Track backprojection export done."),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )


# ═══════════════════════════════════════════════════════════════════════════
# Outer BackprojectionTab
# ═══════════════════════════════════════════════════════════════════════════

class BackprojectionTab(QWidget):
    """Top-level main tab for viewing behavioral state / track backprojections."""

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        self._init_ui()

        if metadata_loader is not None and hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

        if (
            metadata_loader is not None
            and getattr(metadata_loader, "metadata", None) is not None
        ):
            self._on_metadata_updated()

    # ── UI ──────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(6)

        # Title / description
        title = QLabel("🔭 Backprojection")
        title.setStyleSheet("font-size: 14px; font-weight: bold; color: #ddd;")
        outer.addWidget(title)

        desc = QLabel(
            "Overlay behavioral state or track cluster labels on raw image data.  "
            "Requires state / track clustering to have been run first in the "
            "Analysis tab → Single Cell sub-tab."
        )
        desc.setStyleSheet("color: #888; font-size: 11px;")
        desc.setWordWrap(True)
        outer.addWidget(desc)

        self.stack = QStackedWidget()
        outer.addWidget(self.stack, stretch=1)

        # -- Page 0: Placeholder --
        self.placeholder_page = QWidget()
        place_lay = QVBoxLayout(self.placeholder_page)
        self.placeholder_label = QLabel("Load metadata in the Data Preparation tab to begin.")
        self.placeholder_label.setAlignment(Qt.AlignCenter)
        self.placeholder_label.setStyleSheet("color: #888; font-style: italic; font-size: 14px; padding: 20px;")
        place_lay.addWidget(self.placeholder_label)
        self.stack.addWidget(self.placeholder_page)

        # -- Page 1: Main Content --
        self.main_content = QWidget()
        main_lay = QVBoxLayout(self.main_content)
        main_lay.setContentsMargins(0, 0, 0, 0)
        main_lay.setSpacing(6)

        # ── Shared selectors ────────────────────────────────────────
        sel_lay = QHBoxLayout()
        sel_lay.setSpacing(8)

        sel_lay.addWidget(QLabel("Cell type:"))
        self.cell_type_combo = QComboBox()
        self.cell_type_combo.setMinimumWidth(70)
        self.cell_type_combo.currentTextChanged.connect(self._on_context_changed)
        sel_lay.addWidget(self.cell_type_combo)

        sel_lay.addWidget(QLabel("Sample:"))
        self.sample_combo = QComboBox()
        self.sample_combo.addItem("— All samples —")
        self.sample_combo.setMinimumWidth(70)
        self.sample_combo.currentTextChanged.connect(self._on_context_changed)
        sel_lay.addWidget(self.sample_combo)
        sel_lay.addStretch()

        self.avail_lbl = QLabel("Run analysis to begin.")
        self.avail_lbl.setStyleSheet("color: #888; font-size: 11px;")
        sel_lay.addWidget(self.avail_lbl)
        main_lay.addLayout(sel_lay)

        # ── Inner sub-tabs ───────────────────────────────────────────
        self.inner_tabs = QTabWidget()
        main_lay.addWidget(self.inner_tabs, stretch=1)

        self.state_bp = StateBackprojectionSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            sample_getter=self._current_sample,
            parent=self,
        )
        self.inner_tabs.addTab(self.state_bp, "🔬 State Backprojection")

        self.track_bp = TrackBackprojectionSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            sample_getter=self._current_sample,
            parent=self,
        )
        self.inner_tabs.addTab(self.track_bp, "🛤️ Track Backprojection")

        self.stack.addWidget(self.main_content)
        
        # Initial state
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        if md is not None:
            self.stack.setCurrentIndex(1)
        else:
            self.stack.setCurrentIndex(0)

    # ── Accessors ───────────────────────────────────────────────────────

    def _current_cell_type(self) -> str:
        return self.cell_type_combo.currentText()

    def _current_sample(self) -> str:
        txt = self.sample_combo.currentText()
        return "" if txt.startswith("—") else txt

    # ── Metadata update ─────────────────────────────────────────────────

    def _on_metadata_updated(self, *_):
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        if md is None:
            self.stack.setCurrentIndex(0)
            return
        
        self.stack.setCurrentIndex(1)
        cell_types = _detect_sc_cell_types(self.metadata_loader)
        cur_ct = self.cell_type_combo.currentText()
        self.cell_type_combo.blockSignals(True)
        self.cell_type_combo.clear()
        if cell_types:
            self.cell_type_combo.addItems(cell_types)
            if cur_ct in cell_types:
                self.cell_type_combo.setCurrentText(cur_ct)
        self.cell_type_combo.blockSignals(False)

        # Samples
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        cur_s = self.sample_combo.currentText()
        self.sample_combo.blockSignals(True)
        self.sample_combo.clear()
        self.sample_combo.addItem("— All samples —")
        if md is not None and "sample_name" in md.columns:
            for s in sorted(md["sample_name"].astype(str).unique()):
                self.sample_combo.addItem(s)
        if cur_s in [self.sample_combo.itemText(i) for i in range(self.sample_combo.count())]:
            self.sample_combo.setCurrentText(cur_s)
        self.sample_combo.blockSignals(False)

        n_ct = len(cell_types)
        n_s = max(0, self.sample_combo.count() - 1)
        if n_ct:
            self.avail_lbl.setText(f"{n_ct} cell type(s), {n_s} sample(s) available.")
        else:
            self.avail_lbl.setText("No immune/other cell types detected.")

        self._on_context_changed()

    def _on_context_changed(self, *_):
        self.state_bp.on_context_updated()
        self.track_bp.on_context_updated()
