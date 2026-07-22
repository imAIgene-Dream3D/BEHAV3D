"""
BEHAV3D napari plugin – Analysis Tab.

Top-level tab wrapping two sub-tabs:

- **Death Dynamics** — multi-select target cell-type picker, per-target /
  combined death dynamics, interaction analysis with its own interaction
  cell-type picker, per-step Advanced Settings, and processing-queue
  integration. Disabled when no dead channel is configured in metadata
  (interaction analysis stays partially enabled).
- **Single Cell** — placeholder for later phases.

Mirrors the notebook flow (DeathDynamicsPanel / InteractionAnalysisPanel /
MultiOrganoidAnalysisPanel in ``behav3d.widgets.analysis``) but exposes a
single combined picker so users can run any subset of cell types from one
location.
"""
from __future__ import annotations

import datetime
import sys
import traceback
from pathlib import Path
from typing import Iterable

import yaml
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QTabWidget, QTextEdit, QCheckBox, QGroupBox,
    QScrollArea, QSpinBox, QDoubleSpinBox, QComboBox, QFrame, QSizePolicy,
    QToolButton, QMessageBox, QSplitter, QStackedWidget,
    QTreeWidget, QTreeWidgetItem, QMenu, QHeaderView,
    QListWidget, QAbstractItemView,
)
from qtpy.QtCore import Qt, QUrl
from qtpy.QtGui import QDesktopServices

from behav3d.napari._pdf_view import open_pdf_in_napari
from behav3d.napari._results_panel import ResultsPanel
from behav3d.core.qt_help import make_help_row, HelpButton


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════
def _period_from_spins(start_spin, end_spin):
    """Build ``analysis_period_min`` from optional start/end minute spinboxes.

    Spin value ``-1`` means an open bound (from start / to end). Both open
    returns ``None`` (whole movie).
    """
    start = None if start_spin.value() < 0 else float(start_spin.value())
    end = None if end_spin.value() < 0 else float(end_spin.value())
    if start is None and end is None:
        return None
    return (start, end)


def _period_from_t_radios(radio_group, start_spin, end_spin):
    """Build ``analysis_period_t`` from a radio group + integer spinboxes.

    Returns ``None`` when "All timepoints" is selected, otherwise
    ``(start_t, end_t)`` from the spin values.
    """
    # Qt radio: checked text; ipywidgets: value string
    if hasattr(radio_group, 'checkedButton'):
        # Qt path
        btn = radio_group.checkedButton()
        if btn is None or btn.text() == "All timepoints":
            return None
    else:
        if getattr(radio_group, 'value', 'All timepoints') == 'All timepoints':
            return None
    return (int(start_spin.value()), int(end_spin.value()))


def _detect_cell_types(metadata_loader):
    """Return (organoid_types, immune_types, other_types) lists.

    Includes both legacy metadata-registered merged types and the current
    ``cell_type_groups`` (yml-based, post-filtering) groups, routed to their
    detected category. Returns three empty lists when metadata is missing.
    """
    if metadata_loader is None or metadata_loader.metadata is None:
        return [], [], []
    from behav3d.core.metadata import (
        detect_organoid_types_from_metadata,
        detect_immune_cell_types_from_metadata,
        detect_other_cell_types_from_metadata,
        detect_merged_cell_types_from_metadata,
        filter_multicolor_inputs,
    )
    from behav3d.widgets.utils import detect_cell_type_category
    from behav3d.analysis.grouping import list_cell_type_groups, group_category

    md = metadata_loader.metadata
    org = filter_multicolor_inputs(detect_organoid_types_from_metadata(md))
    imm = filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md))
    oth = filter_multicolor_inputs(detect_other_cell_types_from_metadata(md))

    for ct in detect_merged_cell_types_from_metadata(md):
        try:
            category = detect_cell_type_category(ct, md)
        except Exception:
            category = "immune"
        if category == "organoid" and ct not in org:
            org.append(ct)
        elif category == "other" and ct not in oth:
            oth.append(ct)
        elif category == "immune" and ct not in imm:
            imm.append(ct)

    params = getattr(metadata_loader, "behav3d_parameters", {}) or {}
    for group_id in list_cell_type_groups(params):
        category = group_category(params, md, group_id)
        if category == "organoid" and group_id not in org:
            org.append(group_id)
        elif category == "other" and group_id not in oth:
            oth.append(group_id)
        elif group_id not in imm:
            imm.append(group_id)
    return org, imm, oth


def _has_dead_channel(metadata) -> bool:
    return (
        metadata is not None
        and "dead_channel" in metadata.columns
        and metadata["dead_channel"].notna().any()
    )


def _filtered_csv(out_dir: Path, ct: str) -> Path:
    return (
        out_dir
        / "analysis"
        / ct
        / "track_features"
        / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
    )


def _contact_cols_in_csv(csv_path: Path) -> set[str]:
    if not csv_path.exists():
        return set()
    try:
        import pandas as pd

        cols = pd.read_csv(csv_path, nrows=0).columns
        return {c[:-len("_contact")] for c in cols if c.endswith("_contact")}
    except Exception:
        return set()


# Summary statistics offered by the invasiveness per-movie plot. Kept as a
# local constant to avoid importing the (matplotlib-heavy) invasiveness module
# at plugin import time; the actual analysis is imported lazily when run.
_INV_SUMMARY_STATS = ("mean", "median", "max", "auc")


def _invasiveness_targets_in_csv(csv_path: Path) -> list:
    """Return available invasiveness targets in the filtered CSV at *csv_path*."""
    if not csv_path.exists():
        return []
    try:
        import pandas as pd
        from behav3d.analysis.invasiveness_analysis import (
            detect_invasiveness_targets,
        )

        cols = pd.read_csv(csv_path, nrows=0)
        return detect_invasiveness_targets(cols)
    except Exception:
        return []


def _csv_has_column(csv_path: Path, col: str) -> bool:
    if not csv_path.exists():
        return False
    try:
        import pandas as pd

        cols = pd.read_csv(csv_path, nrows=0).columns
        return col in cols
    except Exception:
        return False


# ═══════════════════════════════════════════════════════════════════════════
# Collapsible section helper
# ═══════════════════════════════════════════════════════════════════════════
class CollapsibleSection(QWidget):
    """A collapsible container with an arrow-toggle header — no checkbox.

    Unlike ``QGroupBox.setCheckable(True)`` (which adds a checkbox to the
    title and only enables/disables the contents), this widget genuinely
    hides the body when collapsed and exposes a simple ``contentLayout()``
    for callers to populate.
    """

    def __init__(self, title: str, *, expanded: bool = False, parent=None):
        super().__init__(parent)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(2)

        self._toggle = QToolButton()
        self._toggle.setStyleSheet(
            "QToolButton { border: none; font-weight: bold; color: #ddd; "
            "background: transparent; padding: 2px; text-align: left; }"
            "QToolButton:hover { color: #fff; }"
        )
        self._toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._toggle.setCheckable(True)
        self._toggle.setChecked(expanded)
        self._toggle.setText(title)
        self._toggle.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._toggle.clicked.connect(self._on_toggled)
        outer.addWidget(self._toggle)

        self._content = QFrame()
        self._content.setFrameShape(QFrame.NoFrame)
        self._content_layout = QVBoxLayout(self._content)
        self._content_layout.setContentsMargins(12, 2, 4, 4)
        self._content_layout.setSpacing(4)
        outer.addWidget(self._content)
        self._content.setVisible(expanded)

        self._update_arrow()

    def _on_toggled(self, _checked: bool):
        self._content.setVisible(self._toggle.isChecked())
        self._update_arrow()

    def _update_arrow(self):
        self._toggle.setArrowType(
            Qt.DownArrow if self._toggle.isChecked() else Qt.RightArrow
        )

    def contentLayout(self) -> QVBoxLayout:
        return self._content_layout

    def addWidget(self, w):
        self._content_layout.addWidget(w)

    def addLayout(self, lay):
        self._content_layout.addLayout(lay)


# ═══════════════════════════════════════════════════════════════════════════
# Death Dynamics sub-tab
# ═══════════════════════════════════════════════════════════════════════════
class DeathDynamicsTab(QWidget):
    """Death dynamics + interaction analysis controls.

    UI structure
    ------------
    - Target cells multi-select (organoid + other categories)
    - Step 1 — Death Dynamics: Run (per target) + Run Combined + advanced
    - Step 2 — Interaction Analysis: Interaction cells multi-select +
      Run (per target) + Run Combined + advanced
    - Run All Available + log
    """

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        self._target_checks: dict[str, QCheckBox] = {}
        self._interaction_checks: dict[str, QCheckBox] = {}
        self._interaction_labels: dict[str, QLabel] = {}
        self._target_meta: dict[str, dict] = {}

        # Queue button placeholders (wired by ``_widget``)
        self.btn_queue_dd_single = QPushButton("+🛒")
        self.btn_queue_dd_combined = QPushButton("+🛒")
        self.btn_queue_ia_single = QPushButton("+🛒")
        self.btn_queue_ia_combined = QPushButton("+🛒")
        self.btn_queue_all = QPushButton("+🛒")
        for btn in (
            self.btn_queue_dd_single,
            self.btn_queue_dd_combined,
            self.btn_queue_ia_single,
            self.btn_queue_ia_combined,
            self.btn_queue_all,
        ):
            btn.setFixedSize(36, 28)
            btn.setToolTip("Add to Processing Queue")
            btn.setStyleSheet(
                "QPushButton { background: #1a1a2e; color: #ffc107; "
                "border: 1px solid #ffc107; border-radius: 4px; "
                "font-size: 11px; font-weight: bold; }"
                "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
            )
        # Slightly taller queue button for the Run-All row so it visually
        # matches the larger primary "Run All Available" button below.
        self.btn_queue_all.setFixedSize(44, 38)
        self.btn_queue_all.setToolTip(
            "Add every applicable Death Dynamics / Interaction step for the "
            "current selection to the Processing Queue."
        )

        # Per-step "View result" buttons. Enabled when the corresponding
        # result PDF already exists on disk.
        self.btn_view_dd_single = QPushButton("👁")
        self.btn_view_dd_combined = QPushButton("👁")
        self.btn_view_ia_single = QPushButton("👁")
        self.btn_view_ia_combined = QPushButton("👁")
        for btn in (
            self.btn_view_dd_single,
            self.btn_view_dd_combined,
            self.btn_view_ia_single,
            self.btn_view_ia_combined,
        ):
            btn.setFixedSize(36, 28)
            btn.setToolTip("View the corresponding result PDF in napari")
            btn.setStyleSheet(
                "QPushButton { background: #1a2e1a; color: #6dd56d; "
                "border: 1px solid #6dd56d; border-radius: 4px; "
                "font-size: 12px; font-weight: bold; }"
                "QPushButton:hover { background: #6dd56d; color: #1a2e1a; }"
                "QPushButton:disabled { background: #1a1a1a; color: #555; "
                "border: 1px solid #333; }"
            )
            btn.setEnabled(False)
        self.btn_view_dd_single.clicked.connect(
            lambda: self._on_view_clicked("dd_single")
        )
        self.btn_view_dd_combined.clicked.connect(
            lambda: self._on_view_clicked("dd_combined")
        )
        self.btn_view_ia_single.clicked.connect(
            lambda: self._on_view_clicked("ia_single")
        )
        self.btn_view_ia_combined.clicked.connect(
            lambda: self._on_view_clicked("ia_combined")
        )

        self._init_ui()

        if metadata_loader is not None and hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    # ── UI ──────────────────────────────────────────────────────────────
    def _init_ui(self):
        from behav3d.napari._guided import GuidedPanel, make_back_header
        from behav3d.napari.analysis_guided_copy import (
            DEATH_DYNAMICS_ANALYSES, GUIDED_INTRO,
        )

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Stacked pages: Guided overview (0), isolated settings form (1).
        self._stack = QStackedWidget()
        outer.addWidget(self._stack)

        # Page 0 — Guided explainers.
        self._guided_panel = GuidedPanel(
            DEATH_DYNAMICS_ANALYSES,
            start_cb=self._on_guided_start,
            intro=GUIDED_INTRO,
        )
        self._stack.addWidget(self._guided_panel)

        # Page 1 — the existing settings form, wrapped verbatim, with a Back
        # header above it. `_focus_step` shows only the group(s) relevant to
        # whichever analysis Start was pressed for.
        form_page = QWidget()
        form_outer = QVBoxLayout(form_page)
        form_outer.setContentsMargins(0, 0, 0, 0)
        form_outer.setSpacing(0)
        back_bar, self._focus_title = make_back_header(
            on_back=lambda: self._stack.setCurrentIndex(0)
        )
        form_outer.addWidget(back_bar)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        form_outer.addWidget(scroll)
        self._form_scroll = scroll
        self._stack.addWidget(form_page)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        scroll.setWidget(content)

        # ── Target cells ────────────────────────────────────────────
        self.targets_group = QGroupBox("🎯 Target cells  (organoid / other)")
        self.targets_lay = QVBoxLayout(self.targets_group)
        self.targets_lay.setSpacing(2)
        self._targets_placeholder = QLabel(
            "Load metadata in the Data Preparation tab to see available targets."
        )
        self._targets_placeholder.setStyleSheet(
            "color: #888; font-style: italic;"
        )
        self.targets_lay.addWidget(self._targets_placeholder)
        layout.addWidget(self.targets_group)

        # ── Step 1: Death Dynamics ──────────────────────────────────
        self.dd_group = QGroupBox("Step 1 — Death Dynamics")
        dd_lay = QVBoxLayout(self.dd_group)
        dd_lay.setSpacing(4)

        self.dd_warning_label = QLabel(
            "⚠ No dead channel configured in metadata. "
            "Death Dynamics is unavailable until 'dead_channel' is set."
        )
        self.dd_warning_label.setWordWrap(True)
        self.dd_warning_label.setStyleSheet(
            "background: #4a2222; color: #ffdada; padding: 6px; "
            "border-radius: 3px; font-size: 11px;"
        )
        self.dd_warning_label.setVisible(False)
        dd_lay.addWidget(self.dd_warning_label)

        # ── Condition-column grouping (per-target Death Dynamics only) ──
        # Lets the user group/color the mean +/- SEM plots by one or more
        # metadata condition columns (e.g. organoid_line + macrophage_line),
        # mirroring the state-composition "group by" multi-select.
        dd_lay.addWidget(QLabel("Group plots by condition column(s):"))
        self.list_dd_group_cols = QListWidget()
        self.list_dd_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list_dd_group_cols.setMaximumHeight(70)
        dd_lay.addWidget(self.list_dd_group_cols)

        dd_btn_row = QHBoxLayout()
        self.btn_run_dd_single = QPushButton("▶ Run Death Dynamics (per target)")
        self._style_primary(self.btn_run_dd_single)
        self.btn_run_dd_single.clicked.connect(self._on_run_dd_single)
        dd_btn_row.addWidget(self.btn_run_dd_single, stretch=1)
        dd_btn_row.addWidget(self.btn_queue_dd_single)
        dd_btn_row.addWidget(self.btn_view_dd_single)
        dd_lay.addLayout(dd_btn_row)

        dd_combined_row = QHBoxLayout()
        self.btn_run_dd_combined = QPushButton(
            "▶ Run Combined Death Dynamics (≥2 targets)"
        )
        self._style_primary(self.btn_run_dd_combined)
        self.btn_run_dd_combined.clicked.connect(self._on_run_dd_combined)
        dd_combined_row.addWidget(self.btn_run_dd_combined, stretch=1)
        dd_combined_row.addWidget(self.btn_queue_dd_combined)
        dd_combined_row.addWidget(self.btn_view_dd_combined)
        dd_lay.addLayout(dd_combined_row)

        # ── Missing-death-classification disclaimer ──────────────
        # Shown when targets are selected but the filtered CSV has no
        # ``dead`` column (i.e. the user has not yet generated death
        # classification through Feature Extraction). The DD buttons are
        # disabled while this disclaimer is visible.
        self.dd_disclaimer_label = QLabel()
        self.dd_disclaimer_label.setWordWrap(True)
        self.dd_disclaimer_label.setStyleSheet(
            "background: #4a4222; color: #fff3cd; padding: 6px; "
            "border-radius: 3px; font-size: 11px;"
        )
        self.dd_disclaimer_label.setVisible(False)
        dd_lay.addWidget(self.dd_disclaimer_label)

        # ── Death thresholds info (inline, read-only) ────────────
        # Replaces the previous "Advanced settings" group — there is
        # nothing to tune here, so we just inform the user about the
        # configured thresholds and where to change them.
        self.dd_thresholds_frame = QFrame()
        self.dd_thresholds_frame.setFrameShape(QFrame.StyledPanel)
        self.dd_thresholds_frame.setStyleSheet(
            "QFrame { background: #232323; border: 1px solid #333; "
            "border-radius: 3px; }"
        )
        dd_thr_lay = QVBoxLayout(self.dd_thresholds_frame)
        dd_thr_lay.setContentsMargins(8, 6, 8, 6)
        dd_thr_lay.setSpacing(4)
        dd_thr_header = QLabel("Death thresholds (read-only)")
        dd_thr_header.setStyleSheet(
            "color: #ddd; font-weight: bold; font-size: 11px;"
        )
        dd_thr_lay.addWidget(dd_thr_header)
        self.dd_thresholds_form_host = QWidget()
        self.dd_thresholds_form_layout = QFormLayout(self.dd_thresholds_form_host)
        self.dd_thresholds_form_layout.setContentsMargins(0, 0, 0, 0)
        self.dd_thresholds_form_layout.setSpacing(2)
        dd_thr_lay.addWidget(self.dd_thresholds_form_host)
        dd_thr_info = QLabel(
            "ⓘ Dead-mask % threshold is owned by Feature Extraction. To "
            "change it, re-run Feature Extraction (single or batch) with "
            "the new value."
        )
        dd_thr_info.setWordWrap(True)
        dd_thr_info.setStyleSheet("color: #aaa; font-size: 10px;")
        dd_thr_lay.addWidget(dd_thr_info)
        dd_lay.addWidget(self.dd_thresholds_frame)

        layout.addWidget(self.dd_group)

        # ── Step 2: Interaction Analysis ──────────────────────────
        self.ia_group = QGroupBox("Step 2 — Interaction Analysis")
        ia_lay = QVBoxLayout(self.ia_group)
        ia_lay.setSpacing(4)

        self.ia_warning_label = QLabel(
            "ⓘ No dead channel configured — Interaction Analysis can still "
            "run, but the following will be skipped:\n"
            "  • Alive vs Dead overall / per-sample comparison plots\n"
            "  • Fate-based statistics (n_survives / n_dies, "
            "mean_contacts_survives / mean_contacts_dies)\n"
            "  • Cumulative-to-death curves (Interaction Overview)\n"
            "  • Active Killing dashboard (Interaction Overview)"
        )
        self.ia_warning_label.setWordWrap(True)
        self.ia_warning_label.setStyleSheet(
            "background: #4a4222; color: #fff3cd; padding: 6px; "
            "border-radius: 3px; font-size: 11px;"
        )
        self.ia_warning_label.setVisible(False)
        ia_lay.addWidget(self.ia_warning_label)

        # Per-target disclaimer shown when the metadata has a dead channel
        # but the selected target(s) lack a ``dead`` column in their filtered
        # CSV (e.g. Feature Extraction not yet re-run with a dead-mask % > 0).
        # Combined IA stays enabled in that case; the alive/dead panels are
        # simply omitted by the analysis code.
        self.ia_disclaimer_label = QLabel()
        self.ia_disclaimer_label.setWordWrap(True)
        self.ia_disclaimer_label.setStyleSheet(
            "background: #4a4222; color: #fff3cd; padding: 6px; "
            "border-radius: 3px; font-size: 11px;"
        )
        self.ia_disclaimer_label.setVisible(False)
        ia_lay.addWidget(self.ia_disclaimer_label)

        # Interaction cells selector (separate group inside Step 2)
        self.interaction_group = QGroupBox("🤝 Interaction cells  (immune / other)")
        self.interaction_lay = QVBoxLayout(self.interaction_group)
        self.interaction_lay.setSpacing(2)
        self._interaction_placeholder = QLabel(
            "No interaction cell types detected yet."
        )
        self._interaction_placeholder.setStyleSheet(
            "color: #888; font-style: italic;"
        )
        self.interaction_lay.addWidget(self._interaction_placeholder)
        ia_lay.addWidget(self.interaction_group)

        # ── 2a Per-target Interaction Analysis ───────────────────
        ia_single_info = QLabel(
            "2a — Per-target analysis: cumulative contact over time and "
            "alive vs dead bar plots (one PDF per target × immune type)."
        )
        ia_single_info.setWordWrap(True)
        ia_single_info.setStyleSheet("color: #aaa; font-size: 11px;")
        ia_lay.addWidget(ia_single_info)

        ia_btn_row = QHBoxLayout()
        self.btn_run_ia_single = QPushButton("▶ Run Interaction Analysis (per target)")
        self._style_primary(self.btn_run_ia_single)
        self.btn_run_ia_single.clicked.connect(self._on_run_ia_single)
        ia_btn_row.addWidget(self.btn_run_ia_single, stretch=1)
        ia_btn_row.addWidget(self.btn_queue_ia_single)
        ia_btn_row.addWidget(self.btn_view_ia_single)
        ia_lay.addLayout(ia_btn_row)

        self.ia_single_adv = CollapsibleSection(
            "Per-target settings", expanded=False
        )
        ia_single_form = QFormLayout()
        ia_single_form.setContentsMargins(0, 0, 0, 0)
        ia_single_form.setSpacing(4)
        self.check_group_by_lc = QCheckBox(
            "Group overall plots by immune line condition"
        )
        gbl_row = QHBoxLayout()
        gbl_row.addWidget(self.check_group_by_lc)
        gbl_row.addWidget(HelpButton(
            "Group overall plots by immune line condition",
            "Splits the two per-target overall plots by the interacting "
            "immune cell's line condition (im_{type}_line_condition): the "
            "cumulative-contact curve draws one line per condition, and the "
            "alive-vs-dead bar plot colours bars by condition. Only affects "
            "the 2a per-target overall plots.",
        ))
        gbl_row.addStretch()
        ia_single_form.addRow("", gbl_row)
        ia_single_adv_host = QWidget()
        ia_single_adv_host.setLayout(ia_single_form)
        self.ia_single_adv.contentLayout().addWidget(ia_single_adv_host)
        ia_lay.addWidget(self.ia_single_adv)

        # ── 2b Interaction Overview ──────────────────────────────
        ia_overview_info = QLabel(
            "2b — Interaction Overview: violin (cumulative contacts per "
            "organoid), cumulative-to-death curves, and active-killing "
            "dashboard. Works with one or more targets."
        )
        ia_overview_info.setWordWrap(True)
        ia_overview_info.setStyleSheet("color: #aaa; font-size: 11px;")
        ia_lay.addWidget(ia_overview_info)

        ia_combined_row = QHBoxLayout()
        self.btn_run_ia_combined = QPushButton(
            "▶ Run Interaction Overview"
        )
        self._style_primary(self.btn_run_ia_combined)
        self.btn_run_ia_combined.clicked.connect(self._on_run_ia_combined)
        ia_combined_row.addWidget(self.btn_run_ia_combined, stretch=1)
        ia_combined_row.addWidget(self.btn_queue_ia_combined)
        ia_combined_row.addWidget(self.btn_view_ia_combined)
        ia_lay.addLayout(ia_combined_row)

        self.ia_overview_settings = CollapsibleSection(
            "Interaction Overview settings", expanded=False
        )
        ia_overview_form = QFormLayout()
        ia_overview_form.setContentsMargins(0, 0, 0, 0)
        ia_overview_form.setSpacing(4)

        self.spin_time_window = QSpinBox()
        self.spin_time_window.setRange(1, 99999)
        self.spin_time_window.setValue(60)
        self.spin_time_window.setSuffix(" min")
        self.spin_time_window.setMaximumWidth(100)
        ia_overview_form.addRow("Before-death window:", make_help_row(
            self.spin_time_window, "Before-death window",
            "Lookback window (minutes before Time of Death) for the "
            "cumulative-to-death curve only. Requires death classification.",
        ))

        # Temporal Range: radio buttons + Start T / End T spinboxes
        from qtpy.QtWidgets import QRadioButton, QButtonGroup
        self.radio_period_all = QRadioButton("All timepoints")
        self.radio_period_custom = QRadioButton("Custom time range")
        self.radio_period_all.setChecked(True)
        self.period_radio_group = QButtonGroup()
        self.period_radio_group.addButton(self.radio_period_all)
        self.period_radio_group.addButton(self.radio_period_custom)

        self.spin_period_start_t = QSpinBox()
        self.spin_period_start_t.setRange(0, 999999)
        self.spin_period_start_t.setValue(0)
        self.spin_period_start_t.setPrefix("Start T: ")
        self.spin_period_start_t.setMaximumWidth(130)
        self.spin_period_start_t.setEnabled(False)
        self.spin_period_end_t = QSpinBox()
        self.spin_period_end_t.setRange(0, 999999)
        self.spin_period_end_t.setValue(100)
        self.spin_period_end_t.setPrefix("End T: ")
        self.spin_period_end_t.setMaximumWidth(130)
        self.spin_period_end_t.setEnabled(False)

        def _on_period_radio_toggled(checked):
            custom = self.radio_period_custom.isChecked()
            self.spin_period_start_t.setEnabled(custom)
            self.spin_period_end_t.setEnabled(custom)
        self.radio_period_custom.toggled.connect(_on_period_radio_toggled)

        period_radio_row = QHBoxLayout()
        period_radio_row.addWidget(self.radio_period_all)
        period_radio_row.addWidget(self.radio_period_custom)
        period_radio_row.addStretch()
        period_t_row = QHBoxLayout()
        period_t_row.addWidget(self.spin_period_start_t)
        period_t_row.addWidget(QLabel("to"))
        period_t_row.addWidget(self.spin_period_end_t)
        period_t_row.addStretch()

        temporal_info = QLabel(
            "Temporal Range — applies to: cumulative interactions "
            "per organoid and active killing plots"
        )
        temporal_info.setWordWrap(True)
        temporal_info.setStyleSheet("color: #aaa; font-size: 10px;")
        ia_overview_form.addRow(temporal_info)
        ia_overview_form.addRow("", period_radio_row)
        ia_overview_form.addRow("", period_t_row)

        self.check_annotate_lc = QCheckBox(
            "Annotate by immune line condition"
        )
        annotate_lc_row = QHBoxLayout()
        annotate_lc_row.addWidget(self.check_annotate_lc)
        annotate_lc_row.addWidget(HelpButton(
            "Annotate by immune line condition",
            "When enabled, the violin (cumulative contacts per organoid) and "
            "the active-killing dashboard add the immune line condition as a "
            "hue / line-style annotation within the selected 'Group by' "
            "grouping, instead of pooling all line conditions together. Does "
            "not affect the cumulative-to-death curve.",
        ))
        annotate_lc_row.addStretch()
        ia_overview_form.addRow("", annotate_lc_row)

        self.combo_group_by = QComboBox()
        self.combo_group_by.addItem("By target (organoid) type", userData="organoid_type")
        self.combo_group_by.addItem("By treatment (immune cell)", userData="treatment")
        self.combo_group_by.setMaximumWidth(280)
        ia_overview_form.addRow("Group by:", make_help_row(
            self.combo_group_by, "Group by",
            "How results are grouped on the x-axis in Interaction Overview "
            "(violin, cumulative-to-death curve, active-killing dashboard). "
            "'By target (organoid) type' keeps each organoid type as its own "
            "x-axis group. 'By treatment (immune cell)' pools all organoid "
            "types together and groups by the interacting immune cell type "
            "instead.",
        ))

        ia_overview_host = QWidget()
        ia_overview_host.setLayout(ia_overview_form)
        self.ia_overview_settings.contentLayout().addWidget(ia_overview_host)
        ia_lay.addWidget(self.ia_overview_settings)

        # Overview widgets whose enabled state is managed in
        # ``_refresh_button_states``.
        self._ia_overview_setting_widgets = [
            self.spin_time_window,
            self.radio_period_all,
            self.radio_period_custom,
            self.spin_period_start_t,
            self.spin_period_end_t,
            self.check_annotate_lc,
            self.combo_group_by,
        ]

        layout.addWidget(self.ia_group)

        # ── Step 3: Invasiveness Analysis (immune-cell perspective) ──
        # Self-contained: has its own immune-cell picker + target checkboxes,
        # independent of the Death Dynamics / Interaction selectors above.
        self.inv_group = QGroupBox("Step 3 — Invasiveness Analysis")
        inv_lay = QVBoxLayout(self.inv_group)
        inv_lay.setSpacing(4)

        inv_info_row = QHBoxLayout()
        inv_info_row.addStretch()
        inv_info_row.addWidget(HelpButton(
            "Invasiveness Analysis",
            "Measures invasiveness of the selected immune cell type(s) against "
            "one or more targets, over time and summarized per movie. Enable "
            "'invasiveness' during feature extraction to compute the underlying "
            "values.\n\n"
            "Feature extraction records, for every immune cell at every "
            "timepoint, the % of its surface touching a chosen organoid "
            "target (a '{target}_invasiveness_perc' value, 0-100) and a "
            "boolean invasive flag that is True once that contact reaches "
            "≥50%. This step turns those per-cell, per-timepoint values "
            "into three views:\n\n"
            "• Fraction over time: % of cells flagged invasive (≥50% "
            "contact) at each timepoint.\n"
            "• Mean / median % over time: average / typical contact % "
            "across ALL cells (including non-invasive 0% ones) at each "
            "timepoint.\n"
            "• Per-movie summary: the chosen over-time curve collapsed to "
            "one dot per movie (see the 'Per-movie summary stat' help).\n\n"
            "If organoid types are selected for fate cross-referencing, "
            "fate violins (one point per organoid, requires organoid filtered "
            "CSVs with death data) additionally show contact %/invasive-cell "
            "counts per organoid, split by whether that organoid died.",
        ))
        inv_lay.addLayout(inv_info_row)

        self.inv_immune_group = QGroupBox("Immune cell type(s)")
        self.inv_immune_lay = QVBoxLayout(self.inv_immune_group)
        self.inv_immune_lay.setSpacing(2)
        self._inv_immune_checks: dict[str, QCheckBox] = {}
        self._inv_immune_placeholder = QLabel("No immune cell types detected yet.")
        self._inv_immune_placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.inv_immune_lay.addWidget(self._inv_immune_placeholder)
        inv_lay.addWidget(self.inv_immune_group)

        self.inv_targets_group = QGroupBox("🎯 Targets to compare")
        self.inv_targets_lay = QVBoxLayout(self.inv_targets_group)
        self.inv_targets_lay.setSpacing(2)
        self._inv_target_checks: dict[str, QCheckBox] = {}
        self._inv_targets_placeholder = QLabel(
            "No invasiveness targets detected yet."
        )
        self._inv_targets_placeholder.setStyleSheet(
            "color: #888; font-style: italic;"
        )
        self.inv_targets_lay.addWidget(self._inv_targets_placeholder)
        inv_lay.addWidget(self.inv_targets_group)

        inv_stat_row = QHBoxLayout()
        inv_stat_row.addWidget(QLabel("Per-movie summary stat:"))
        self.inv_summary_combo = QComboBox()
        for _s in _INV_SUMMARY_STATS:
            self.inv_summary_combo.addItem(_s, userData=_s)
        self.inv_summary_combo.setMaximumWidth(160)
        inv_stat_row.addWidget(self.inv_summary_combo)
        inv_stat_row.addWidget(HelpButton(
            "Per-movie summary stat",
            "Collapses each movie's time series to one value (mean/median/max/AUC). "
            "Does not affect the separate mean/median-over-time curves.\n\n"
            "Defines the single number plotted as one dot per movie in the "
            "'Per-movie' plots (for both the boolean fraction and the "
            "percentage). For each movie (per immune type / target / line "
            "condition), the per-timepoint curve is first built by averaging "
            "across that movie's cells at each timepoint (mean %, or the "
            "per-timepoint median % when 'median' is selected), then "
            "collapsed over time using the selected statistic:\n\n"
            "• mean — average of the per-timepoint values.\n"
            "• median — median of the per-timepoint values.\n"
            "• max — highest per-timepoint value reached.\n"
            "• AUC — area under the per-timepoint curve, normalized by the "
            "time span (trapezoidal integral / duration), so it stays on "
            "the same 0-1 / 0-100 scale as the other statistics.\n\n"
            "This choice only affects the per-movie summary plots/CSV, not "
            "the separate fraction/mean/median over-time curves.",
        ))
        inv_stat_row.addStretch()
        inv_lay.addLayout(inv_stat_row)

        self.check_inv_group_by_lc = QCheckBox(
            "Separate by immune line condition"
        )
        inv_lc_row = QHBoxLayout()
        inv_lc_row.addWidget(self.check_inv_group_by_lc)
        inv_lc_row.addWidget(HelpButton(
            "Separate by immune line condition",
            "When enabled, results are additionally split by the immune "
            "cell's line condition (im_{type}_line_condition): a distinct "
            "box/line per line condition in the per-movie summary and "
            "over-time plots, instead of pooling all conditions together.",
        ))
        inv_lc_row.addStretch()
        inv_lay.addLayout(inv_lc_row)

        # Invasiveness temporal range: radio + Start T / End T
        from qtpy.QtWidgets import QRadioButton, QButtonGroup
        self.inv_radio_period_all = QRadioButton("All timepoints")
        self.inv_radio_period_custom = QRadioButton("Custom time range")
        self.inv_radio_period_all.setChecked(True)
        self.inv_period_radio_group = QButtonGroup()
        self.inv_period_radio_group.addButton(self.inv_radio_period_all)
        self.inv_period_radio_group.addButton(self.inv_radio_period_custom)

        self.spin_inv_start_t = QSpinBox()
        self.spin_inv_start_t.setRange(0, 999999)
        self.spin_inv_start_t.setValue(0)
        self.spin_inv_start_t.setPrefix("Start T: ")
        self.spin_inv_start_t.setMaximumWidth(130)
        self.spin_inv_start_t.setEnabled(False)
        self.spin_inv_end_t = QSpinBox()
        self.spin_inv_end_t.setRange(0, 999999)
        self.spin_inv_end_t.setValue(100)
        self.spin_inv_end_t.setPrefix("End T: ")
        self.spin_inv_end_t.setMaximumWidth(130)
        self.spin_inv_end_t.setEnabled(False)

        def _on_inv_period_radio_toggled(checked):
            custom = self.inv_radio_period_custom.isChecked()
            self.spin_inv_start_t.setEnabled(custom)
            self.spin_inv_end_t.setEnabled(custom)
        self.inv_radio_period_custom.toggled.connect(_on_inv_period_radio_toggled)

        inv_period_radio_row = QHBoxLayout()
        inv_period_radio_row.addWidget(self.inv_radio_period_all)
        inv_period_radio_row.addWidget(self.inv_radio_period_custom)
        inv_period_radio_row.addStretch()
        inv_period_t_row = QHBoxLayout()
        inv_period_t_row.addWidget(self.spin_inv_start_t)
        inv_period_t_row.addWidget(QLabel("to"))
        inv_period_t_row.addWidget(self.spin_inv_end_t)
        inv_period_t_row.addStretch()

        inv_temporal_label = QLabel("Temporal Range (in timepoints):")
        inv_temporal_label.setStyleSheet("color: #aaa; font-size: 10px;")
        inv_lay.addWidget(inv_temporal_label)
        inv_lay.addLayout(inv_period_radio_row)
        inv_lay.addLayout(inv_period_t_row)

        # Queue + view buttons for the invasiveness step.
        self.btn_queue_inv = QPushButton("+🛒")
        self.btn_queue_inv.setFixedSize(36, 28)
        self.btn_queue_inv.setToolTip("Add to Processing Queue")
        self.btn_queue_inv.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; "
            "border: 1px solid #ffc107; border-radius: 4px; "
            "font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )
        self.btn_view_inv = QPushButton("👁")
        self.btn_view_inv.setFixedSize(36, 28)
        self.btn_view_inv.setToolTip(
            "View the invasiveness result PDF in napari"
        )
        self.btn_view_inv.setStyleSheet(
            "QPushButton { background: #1a2e1a; color: #6dd56d; "
            "border: 1px solid #6dd56d; border-radius: 4px; "
            "font-size: 12px; font-weight: bold; }"
            "QPushButton:hover { background: #6dd56d; color: #1a2e1a; }"
            "QPushButton:disabled { background: #1a1a1a; color: #555; "
            "border: 1px solid #333; }"
        )
        self.btn_view_inv.setEnabled(False)
        self.btn_view_inv.clicked.connect(
            lambda: self._on_view_clicked("invasiveness")
        )

        inv_btn_row = QHBoxLayout()
        self.btn_run_inv = QPushButton("▶ Run Invasiveness Analysis")
        self._style_primary(self.btn_run_inv)
        self.btn_run_inv.clicked.connect(self._on_run_invasiveness)
        inv_btn_row.addWidget(self.btn_run_inv, stretch=1)
        inv_btn_row.addWidget(self.btn_queue_inv)
        inv_btn_row.addWidget(self.btn_view_inv)
        inv_lay.addLayout(inv_btn_row)

        layout.addWidget(self.inv_group)

        # ── Run All Available + Log ───────────────────────────────
        # Wrapped in a container widget so it can be hidden as a unit while
        # a single step is focused (running "all" doesn't make sense when
        # only one step's settings are on screen).
        run_all_row = QHBoxLayout()
        self.btn_run_all = QPushButton("▶▶  Run All Available")
        self.btn_run_all.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_all.clicked.connect(self._on_run_all_clicked)
        run_all_row.addWidget(self.btn_run_all, stretch=1)
        run_all_row.addWidget(self.btn_queue_all)
        self.run_all_container = QWidget()
        self.run_all_container.setLayout(run_all_row)
        layout.addWidget(self.run_all_container)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(140)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        layout.addStretch()

        # Try a first build (in case metadata is already loaded)
        self._rebuild()

        # Always land on the Guided overview; a specific analysis's settings
        # are only reached via that analysis's Start button.
        self._stack.setCurrentIndex(0)

    # ── Guided overview / focused settings ──────────────────────────────────
    def _on_guided_start(self, analysis_id: str):
        """Start from a Guided card: show only that analysis's settings."""
        self._focus_step(analysis_id)
        self._stack.setCurrentIndex(1)
        self._form_scroll.verticalScrollBar().setValue(0)

    def _focus_step(self, analysis_id: str):
        """Show only the group(s) relevant to ``analysis_id``, hide the rest."""
        from behav3d.napari.analysis_guided_copy import (
            DEATH_DYNAMICS, INTERACTION, INVASIVENESS,
        )

        visible = {
            "death_dynamics": {self.targets_group, self.dd_group},
            "interaction": {self.targets_group, self.ia_group},
            "invasiveness": {self.inv_group},
        }.get(analysis_id, set())
        for group in (self.targets_group, self.dd_group, self.ia_group, self.inv_group):
            group.setVisible(group in visible)
        self.run_all_container.setVisible(False)

        title = {
            "death_dynamics": DEATH_DYNAMICS["title"],
            "interaction": INTERACTION["title"],
            "invasiveness": INVASIVENESS["title"],
        }.get(analysis_id, "")
        self._focus_title.setText(title)

    def _style_primary(self, btn: QPushButton):
        btn.setStyleSheet(
            "QPushButton { background: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px; font-size: 12px; } "
            "QPushButton:hover:!disabled { background: #218838; } "
            "QPushButton:disabled { background: #3a3a3a; color: #888; }"
        )

    # ── Log ─────────────────────────────────────────────────────────────
    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    # ── Metadata wiring ─────────────────────────────────────────────────
    def _on_metadata_updated(self, *_):
        self._rebuild()

    def _clear_layout(self, lay):
        while lay.count():
            item = lay.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)

    def _rebuild(self):
        """Re-detect cell types and rebuild the target/interaction selectors."""
        self._target_checks.clear()
        self._interaction_checks.clear()
        self._interaction_labels.clear()
        self._target_meta.clear()
        self._clear_layout(self.targets_lay)
        self._clear_layout(self.interaction_lay)
        self._clear_layout(self.dd_thresholds_form_layout)

        if (
            self.metadata_loader is None
            or self.metadata_loader.metadata is None
            or not self.metadata_loader.output_dir
        ):
            self.targets_lay.addWidget(self._targets_placeholder)
            self.interaction_lay.addWidget(self._interaction_placeholder)
            self.list_dd_group_cols.clear()
            self._refresh_button_states()
            return

        organoid_types, immune_types, other_types = _detect_cell_types(
            self.metadata_loader
        )
        out_dir = Path(self.metadata_loader.output_dir)
        md = self.metadata_loader.metadata
        has_dead = _has_dead_channel(md)

        # ---- Condition-column grouping options -----------------------
        prev_selected = {item.text() for item in self.list_dd_group_cols.selectedItems()}
        self.list_dd_group_cols.clear()
        group_col_options = [
            c for c in md.columns
            if c in ("exp_nr", "well") or c.endswith("_line_condition")
        ]
        for col in group_col_options:
            self.list_dd_group_cols.addItem(col)
            if col in prev_selected:
                self.list_dd_group_cols.item(self.list_dd_group_cols.count() - 1).setSelected(True)

        # ---- Target cells (organoid + other) -----------------------
        target_types = list(organoid_types) + list(other_types)
        if not target_types:
            self.targets_lay.addWidget(
                QLabel("⚠ No organoid or other cell types detected in metadata.")
            )
        else:
            for ct in target_types:
                csv = _filtered_csv(out_dir, ct)
                has_filtered = csv.exists()
                has_dead_col = _csv_has_column(csv, "dead") if has_filtered else False
                self._target_meta[ct] = {
                    "filtered": has_filtered,
                    "dead_col": has_dead_col,
                    "csv": csv,
                }
                cb = QCheckBox(ct)
                cb.setEnabled(has_filtered)
                cb.toggled.connect(lambda _checked, _ct=ct: self._on_target_toggled(_ct))
                badge_parts = []
                badge_parts.append(
                    "✅ filtered" if has_filtered else "⚠ no filtered CSV (run Filtering)"
                )
                if has_filtered:
                    badge_parts.append(
                        "✅ dead-col" if has_dead_col else "⚠ no dead column"
                    )
                badge = QLabel(" · ".join(badge_parts))
                badge.setStyleSheet("color: #aaa; font-size: 10px;")
                row = QHBoxLayout()
                row.addWidget(cb)
                row.addWidget(badge)
                row.addStretch()
                container = QWidget()
                container.setLayout(row)
                self.targets_lay.addWidget(container)
                if not has_filtered:
                    cb.setToolTip(
                        "This target has no filtered features CSV. Run "
                        "Filtering for this cell type first."
                    )
                self._target_checks[ct] = cb

        # ---- Per-target read-only thresholds (inline info block) --
        params = (
            self.metadata_loader.behav3d_parameters
            if hasattr(self.metadata_loader, "behav3d_parameters")
            else {}
        ) or {}
        features = (params.get("features") or {})
        if not target_types:
            empty_lbl = QLabel("(no targets detected)")
            empty_lbl.setStyleSheet("color: #888; font-style: italic; font-size: 11px;")
            self.dd_thresholds_form_layout.addRow(empty_lbl)
        else:
            for ct in target_types:
                thr = features.get(ct, {}).get("dead_mask_percentage_threshold")
                thr_str = (
                    f"{float(thr) * 100.0:.4g} %" if thr not in (None, 0, 0.0)
                    else "<not configured>"
                )
                lbl = QLabel(thr_str)
                lbl.setStyleSheet("color: #ccc; font-size: 11px;")
                self.dd_thresholds_form_layout.addRow(f"{ct}:", lbl)

        # ---- Interaction cells (immune + other) ----------------------
        interaction_types = list(immune_types) + list(other_types)
        # Avoid listing an organoid that is also being treated as a
        # target (e.g. a merged grouped output). Targets and interactions
        # are conceptually different, but we still surface "other" types
        # in both panels; the user can choose.
        if not interaction_types:
            self.interaction_lay.addWidget(
                QLabel("⚠ No immune or other cell types detected in metadata.")
            )
        else:
            for ct in interaction_types:
                cb = QCheckBox(ct)
                cb.toggled.connect(
                    lambda _checked, _ct=ct: self._on_interaction_toggled(_ct)
                )
                badge = QLabel("")
                badge.setStyleSheet("color: #aaa; font-size: 10px;")
                row = QHBoxLayout()
                row.addWidget(cb)
                row.addWidget(badge)
                row.addStretch()
                container = QWidget()
                container.setLayout(row)
                self.interaction_lay.addWidget(container)
                self._interaction_checks[ct] = cb
                self._interaction_labels[ct] = badge

        # ---- Dead-channel gating -----------------------------------
        self.dd_warning_label.setVisible(not has_dead)
        self.ia_warning_label.setVisible(not has_dead)

        # Show / hide Step 1 controls (warning replaces buttons when off)
        for w in (
            self.btn_run_dd_single,
            self.btn_queue_dd_single,
            self.btn_run_dd_combined,
            self.btn_queue_dd_combined,
            self.dd_thresholds_frame,
        ):
            w.setVisible(has_dead)
        # Disclaimer visibility is governed by selection state in
        # ``_refresh_button_states``; hide it whenever Step 1 is fully off.
        if not has_dead:
            self.dd_disclaimer_label.setVisible(False)

        self._refresh_interaction_availability()
        self._refresh_button_states()
        self._rebuild_invasiveness()

    # ── State helpers ──────────────────────────────────────────────────
    def _selected_targets(self) -> list[str]:
        return [
            ct for ct, cb in self._target_checks.items()
            if cb.isEnabled() and cb.isChecked()
        ]

    def _selected_interactions(self) -> list[str]:
        return [ct for ct, cb in self._interaction_checks.items() if cb.isChecked()]

    def _on_target_toggled(self, _ct: str):
        self._refresh_interaction_availability()
        self._refresh_button_states()

    def _on_interaction_toggled(self, _ct: str):
        self._refresh_button_states()

    def _refresh_interaction_availability(self):
        """Annotate interaction-cell rows with per-target contact-column status."""
        targets = self._selected_targets()
        for ct, cb in self._interaction_checks.items():
            label = self._interaction_labels.get(ct)
            if not targets:
                if label is not None:
                    label.setText("(select target(s) above)")
                cb.setEnabled(False)
                continue
            parts: list[str] = []
            any_available = False
            for tgt in targets:
                meta = self._target_meta.get(tgt, {})
                csv: Path | None = meta.get("csv") if meta else None
                if csv and csv.exists():
                    contacts = _contact_cols_in_csv(csv)
                else:
                    contacts = set()
                ok = ct in contacts
                any_available = any_available or ok
                parts.append(f"{tgt}: {'✅' if ok else '⚠'}")
            cb.setEnabled(any_available)
            if not any_available and cb.isChecked():
                cb.blockSignals(True)
                cb.setChecked(False)
                cb.blockSignals(False)
            if label is not None:
                label.setText(", ".join(parts))

    def _refresh_button_states(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        has_dead = _has_dead_channel(
            self.metadata_loader.metadata if self.metadata_loader else None
        )
        # Targets whose filtered CSV actually contains a ``dead`` column —
        # i.e. Feature Extraction has been run with a dead-mask % threshold
        # so death classification is available for Death Dynamics.
        targets_with_dead = [
            ct for ct in targets
            if self._target_meta.get(ct, {}).get("dead_col")
        ]
        missing_dead_targets = [
            ct for ct in targets if ct not in targets_with_dead
        ]

        # Step 1 — Death Dynamics requires both the dead channel in
        # metadata AND a ``dead`` column in the selected target's CSV.
        dd_ready_single = bool(targets_with_dead) and has_dead
        dd_ready_combined = (len(targets_with_dead) >= 2) and has_dead
        self.btn_run_dd_single.setEnabled(dd_ready_single)
        self.btn_queue_dd_single.setEnabled(dd_ready_single)
        self.btn_run_dd_combined.setEnabled(dd_ready_combined)
        self.btn_queue_dd_combined.setEnabled(dd_ready_combined)

        # Disclaimer for missing death classification (only relevant when
        # the metadata says a dead channel is configured — otherwise the
        # whole Step 1 is hidden by ``_rebuild``).
        if has_dead and targets and missing_dead_targets:
            missing = ", ".join(missing_dead_targets)
            self.dd_disclaimer_label.setText(
                f"⚠ No death classification found for: {missing}.\n"
                "Re-run Feature Extraction with a Dead-mask % threshold > 0 "
                "(Feature Extraction tab) to generate the 'dead' column."
            )
            self.dd_disclaimer_label.setVisible(True)
        else:
            self.dd_disclaimer_label.setVisible(False)

        dd_tip = ""
        if has_dead and targets and not targets_with_dead:
            dd_tip = (
                "Selected target(s) have no 'dead' column in their filtered "
                "CSV. Re-run Feature Extraction with a Dead-mask % "
                "threshold > 0."
            )
        elif has_dead and not targets:
            dd_tip = "Select at least one target above."
        self.btn_run_dd_single.setToolTip(dd_tip)
        self.btn_queue_dd_single.setToolTip(dd_tip)
        if has_dead and len(targets_with_dead) < 2 and len(targets) >= 2:
            dd_tip_combined = (
                "Combined Death Dynamics needs at least 2 targets that have "
                "the 'dead' column in their filtered CSV."
            )
        else:
            dd_tip_combined = dd_tip
        self.btn_run_dd_combined.setToolTip(dd_tip_combined)
        self.btn_queue_dd_combined.setToolTip(dd_tip_combined)

        # Step 2 — Interaction Analysis runs with or without dead data.
        # The underlying ``run_interaction_analysis`` and
        # ``run_multi_organoid_interaction_comparison`` functions detect the
        # absence of a ``dead`` column and silently skip the alive/dead
        # comparison plots, cumulative-to-death curves and active-killing
        # dashboard — see ``_process_death_classification`` and the
        # ``has_dead_data`` branches in
        # ``behav3d/analysis/interaction_analysis.py``.
        ia_ready_single = bool(targets) and bool(interactions)
        self.btn_run_ia_single.setEnabled(ia_ready_single)
        self.btn_queue_ia_single.setEnabled(ia_ready_single)
        self.check_group_by_lc.setEnabled(ia_ready_single)
        ia_ready_overview = ia_ready_single
        self.btn_run_ia_combined.setEnabled(ia_ready_overview)
        self.btn_queue_ia_combined.setEnabled(ia_ready_overview)
        for w in self._ia_overview_setting_widgets:
            w.setEnabled(ia_ready_overview)
        # Before-death window only meaningful when death data exists.
        overview_has_dead = bool(targets_with_dead) and has_dead
        self.spin_time_window.setEnabled(ia_ready_overview and overview_has_dead)
        # Temporal T spinboxes: respect custom-range toggle when overview runs.
        if ia_ready_overview:
            custom_period = self.radio_period_custom.isChecked()
            self.spin_period_start_t.setEnabled(custom_period)
            self.spin_period_end_t.setEnabled(custom_period)
        tip_overview = ""
        if not targets:
            tip_overview = "Select at least one target above."
        elif not interactions:
            tip_overview = "Select at least one interaction cell type."
        elif not has_dead:
            tip_overview = (
                "Will run without dead data — fate split, cumulative-to-death "
                "curves and active-killing dashboard will be skipped."
            )
        elif len(targets) == 1:
            tip_overview = (
                "Summary plots for the selected target: violin, "
                "before-death curves, active-killing dashboard."
            )
        else:
            tip_overview = (
                "Summary plots across selected targets on the same figure."
            )
        self.btn_run_ia_combined.setToolTip(tip_overview)
        self.btn_queue_ia_combined.setToolTip(tip_overview)
        if ia_ready_overview and not overview_has_dead:
            self.spin_time_window.setToolTip(
                "Requires death classification in the filtered CSV."
            )
        else:
            self.spin_time_window.setToolTip(
                "Lookback window (minutes before Time of Death) for the "
                "cumulative-to-death curve only."
            )

        # IA disclaimer — shown when targets are selected but their filtered
        # CSV(s) lack a ``dead`` column (Feature Extraction not yet run with
        # a dead-mask % > 0) AND the dead channel IS configured in metadata.
        # When there is no dead channel at all, ``ia_warning_label`` already
        # tells the user what gets skipped, so we don't duplicate that here.
        if has_dead and targets and missing_dead_targets:
            missing = ", ".join(missing_dead_targets)
            self.ia_disclaimer_label.setText(
                f"ⓘ No 'dead' column found in filtered CSV for: {missing}.\n"
                "Interaction Analysis will still run, but alive vs dead "
                "comparisons, cumulative-to-death curves and the active "
                "killing dashboard will be skipped. Re-run Feature "
                "Extraction with a Dead-mask % threshold > 0 to enable them."
            )
            self.ia_disclaimer_label.setVisible(True)
        else:
            self.ia_disclaimer_label.setVisible(False)

        # "Run All Available" is *strict*: it runs immediately so its
        # inputs (filtered CSV with 'dead' column, contact columns, …)
        # must already be on disk. Enabled only when at least one
        # analysis step is runnable right now.
        run_all_ready = bool(targets) and (
            bool(targets_with_dead) or bool(interactions)
        )
        self.btn_run_all.setEnabled(run_all_ready)
        if not targets:
            run_all_tip = "Select at least one target above."
        elif not (targets_with_dead or interactions):
            run_all_tip = (
                "Nothing to run for the current selection — either pick at "
                "least one interaction cell type, or select targets that "
                "have a 'dead' column in their filtered CSV."
            )
        else:
            run_all_tip = ""
        self.btn_run_all.setToolTip(run_all_tip)

        # "+🛒" next to Run All Available is *lenient*: the steps it
        # queues run later, by which time Filter / Feature Extraction
        # steps earlier in the queue may have produced the missing
        # inputs. So we only require static conditions (dead channel in
        # metadata OR at least one interaction selected) — per-step
        # runners skip individual targets gracefully if their filtered
        # CSV / dead column / contact column is still missing at
        # execution time.
        queue_all_ready = bool(targets) and (has_dead or bool(interactions))
        self.btn_queue_all.setEnabled(queue_all_ready)
        if not targets:
            queue_all_tip = "Select at least one target above."
        elif not (has_dead or interactions):
            queue_all_tip = (
                "Nothing to queue for the current selection — either "
                "configure a dead channel in metadata, or pick at least "
                "one interaction cell type."
            )
        else:
            queue_all_tip = (
                "Add every applicable Death Dynamics / Interaction step "
                "for the current selection to the Processing Queue.\n"
                "Individual steps will skip themselves at run time if "
                "their inputs (filtered CSV, 'dead' column, contact "
                "columns) are still missing."
            )
        self.btn_queue_all.setToolTip(queue_all_tip)

        # Per-step "View result" buttons: enabled only when the matching
        # PDF already exists on disk.
        self._update_view_buttons()

    # ── View-result buttons & Results-panel notifications ──────────────
    def _view_pdf_candidates(self, kind: str) -> list[tuple[str, Path]]:
        """Return ``(label, path)`` tuples for every PDF the View button
        for *kind* could open given the current selection.

        For single-target / single-interaction selections this yields one
        entry, matching the original behaviour. For multi-selections we
        return one entry per (target) or (target, interaction) pair so
        the click handler can offer a chooser menu.
        """
        if self.metadata_loader is None or not self.metadata_loader.output_dir:
            return []
        out_dir = Path(self.metadata_loader.output_dir).expanduser()
        targets = self._selected_targets()
        interactions = self._selected_interactions()

        if kind == "dd_single":
            return [
                (
                    ct,
                    out_dir / "analysis" / ct / "results"
                    / f"combined_general_{ct}_dynamics_analysis.pdf",
                )
                for ct in targets
            ]
        if kind == "dd_combined":
            return [
                (
                    "Combined Death Dynamics",
                    out_dir / "analysis" / "multi_organoid_comparison"
                    / "multi_organoid_death_dynamics_comparison.pdf",
                )
            ]
        if kind == "ia_single":
            # Per-target interaction PDFs land in
            # ``analysis/{ct}/interaction_analysis/`` (see
            # ``run_interaction_analysis`` in
            # ``behav3d/analysis/interaction_analysis.py``), NOT in the
            # ``results`` sub-folder used for Death Dynamics.
            return [
                (
                    f"{ct} ↔ {it}",
                    out_dir / "analysis" / ct / "interaction_analysis"
                    / f"interaction_analysis_{ct}_vs_{it}.pdf",
                )
                for ct in targets
                for it in interactions
            ]
        if kind == "ia_combined":
            comp_dir = out_dir / "analysis" / "multi_organoid_comparison"
            pdfs = sorted(comp_dir.glob("multi_organoid_interaction_comparison*.pdf"))
            if not pdfs:
                return [
                    (
                        "Interaction Overview",
                        comp_dir / "multi_organoid_interaction_comparison.pdf",
                    )
                ]
            return [
                (p.stem.replace("multi_organoid_interaction_comparison", "Overview").strip("_") or "Overview", p)
                for p in pdfs
            ]
        if kind == "invasiveness":
            immune_list = self._selected_invasiveness_immune()
            if not immune_list:
                return []
            # Look in each selected immune type's folder plus the combined
            # multi-immune folder, matching all suffixed output PDFs.
            search_dirs = [
                out_dir / "analysis" / im / "invasiveness_analysis"
                for im in immune_list
            ]
            search_dirs.append(out_dir / "analysis" / "invasiveness_analysis")
            pdfs = []
            for d in search_dirs:
                if d.exists():
                    pdfs.extend(sorted(d.glob("invasiveness_analysis_*.pdf")))
            # De-duplicate while preserving order.
            seen = set()
            uniq = []
            for p in pdfs:
                if p not in seen:
                    seen.add(p)
                    uniq.append(p)
            if not uniq:
                im0 = immune_list[0]
                return [(
                    f"Invasiveness — {im0}",
                    out_dir / "analysis" / im0 / "invasiveness_analysis"
                    / f"invasiveness_analysis_{im0}.pdf",
                )]
            return [(p.stem, p) for p in uniq]
        return []

    def _update_view_buttons(self):
        def _set(btn: QPushButton, kind: str, *, empty_tip: str = ""):
            candidates = self._view_pdf_candidates(kind)
            if not candidates:
                btn.setEnabled(False)
                btn.setToolTip(
                    empty_tip
                    or "Use the Results panel below to browse outputs."
                )
                return
            existing = [(lbl, p) for (lbl, p) in candidates if p.exists()]
            btn.setEnabled(bool(existing))
            if existing:
                if len(existing) == 1:
                    btn.setToolTip(f"Open in napari:\n{existing[0][1]}")
                else:
                    listed = "\n".join(
                        f"  • {lbl}" for lbl, _p in existing
                    )
                    btn.setToolTip(
                        "Click to choose which PDF to open in napari:\n"
                        f"{listed}"
                    )
            else:
                missing_paths = "\n".join(str(p) for _lbl, p in candidates)
                btn.setToolTip(
                    "Result PDF not found yet (will be enabled once "
                    f"produced):\n{missing_paths}"
                )

        _set(
            self.btn_view_dd_single, "dd_single",
            empty_tip=(
                "Select at least one target to view its combined "
                "death-dynamics PDF."
            ),
        )
        _set(self.btn_view_dd_combined, "dd_combined")
        _set(
            self.btn_view_ia_single, "ia_single",
            empty_tip=(
                "Select at least one target and one interaction cell type "
                "to view the matching interaction-analysis PDF(s)."
            ),
        )
        _set(self.btn_view_ia_combined, "ia_combined")
        _set(
            self.btn_view_inv, "invasiveness",
            empty_tip=(
                "Select an immune cell type and target(s) to view its "
                "invasiveness-analysis PDF."
            ),
        )

    def _open_pdf_in_napari(self, path: Path):
        """Open *path* in napari using the shared Results panel's settings.

        Centralised so all View buttons (including the chooser menu) share
        the same DPI / reuse-layer behaviour.
        """
        try:
            panel = self._results_panel()
            dpi = int(panel.spin_dpi.value()) if panel else 150
            chk_reuse = getattr(panel, "chk_reuse", None) if panel else None
            reuse = chk_reuse.isChecked() if chk_reuse is not None else True
            open_pdf_in_napari(path, dpi=dpi, reuse=reuse)
        except Exception as e:
            traceback.print_exc()
            self._log(f"Could not open PDF in napari: {e}")

    def _on_view_clicked(self, kind: str):
        candidates = self._view_pdf_candidates(kind)
        existing = [(lbl, p) for (lbl, p) in candidates if p.exists()]
        if not existing:
            return
        if len(existing) == 1:
            self._open_pdf_in_napari(existing[0][1])
            return

        # Multiple PDFs available — show a chooser menu anchored under the
        # clicked button so the user can pick which one to open.
        btn = {
            "dd_single": self.btn_view_dd_single,
            "dd_combined": self.btn_view_dd_combined,
            "ia_single": self.btn_view_ia_single,
            "ia_combined": self.btn_view_ia_combined,
            "invasiveness": self.btn_view_inv,
        }.get(kind)
        menu = QMenu(self)
        for lbl, path in existing:
            act = menu.addAction(f"👁  {lbl}")
            act.triggered.connect(
                lambda _checked=False, _p=path: self._open_pdf_in_napari(_p)
            )
        if btn is not None:
            menu.exec_(btn.mapToGlobal(btn.rect().bottomLeft()))
        else:
            menu.exec_()

    def _results_panel(self):
        """Climb the parent chain to find the shared :class:`ResultsPanel`."""
        node = self.parent()
        while node is not None and not hasattr(node, "results_panel"):
            node = node.parent()
        return getattr(node, "results_panel", None) if node else None

    def _notify_results_changed(self):
        """Ask the shared Results panel (if any) to re-scan disk."""
        panel = self._results_panel()
        if panel is not None:
            try:
                panel.refresh()
            except Exception:
                # Never let a UI refresh failure mask an analysis error.
                traceback.print_exc()
        self._update_view_buttons()

    # ── Param persistence ──────────────────────────────────────────────
    def _persist_advanced(self):
        if self.metadata_loader is None:
            return
        params = self.metadata_loader.behav3d_parameters
        if params is None:
            return
        cfg = params.setdefault("analysis_tab", {}).setdefault("interaction", {})
        cfg["time_window_min"] = int(self.spin_time_window.value())
        cfg["group_by"] = self.combo_group_by.currentData() or "organoid_type"
        cfg["annotate_line_condition"] = bool(self.check_annotate_lc.isChecked())
        cfg["group_by_line_condition"] = bool(self.check_group_by_lc.isChecked())
        cfg["analysis_period_t"] = _period_from_t_radios(
            self.period_radio_group, self.spin_period_start_t, self.spin_period_end_t,
        )
        inv_cfg = params.setdefault("analysis_tab", {}).setdefault("invasiveness", {})
        inv_cfg["group_by_line_condition"] = bool(self.check_inv_group_by_lc.isChecked())
        inv_cfg["analysis_period_t"] = _period_from_t_radios(
            self.inv_period_radio_group, self.spin_inv_start_t, self.spin_inv_end_t,
        )
        params_path = getattr(self.metadata_loader, "behav3d_parameters_path", None)
        if params_path is not None:
            try:
                with open(params_path, "w", encoding="utf-8") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self._log(f"Warning: could not persist analysis params: {e}")
        else:
            out_dir = self.metadata_loader.output_dir
            if out_dir:
                try:
                    p = Path(out_dir) / "behav3d_parameters.yml"
                    with open(p, "w", encoding="utf-8") as f:
                        yaml.safe_dump(params, f, sort_keys=False)
                except Exception as e:
                    self._log(f"Warning: could not persist analysis params: {e}")

    # ── Open-folder pop-up ─────────────────────────────────────────────
    def _offer_open_results_folder(self, folder: Path, what: str = "Death Dynamics"):
        """Pop a dialog offering to open the results folder in the OS file manager.

        Silently no-ops if ``folder`` does not exist (e.g. analysis failed).
        """
        try:
            if not folder or not Path(folder).exists():
                return
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Information)
            box.setWindowTitle(f"{what} complete")
            box.setText(f"{what} analysis is complete.")
            box.setInformativeText(
                f"Results have been saved to:\n{folder}\n\n"
                "Would you like to open the results folder?"
            )
            btn_open = box.addButton("Open folder", QMessageBox.AcceptRole)
            box.addButton("Close", QMessageBox.RejectRole)
            box.setDefaultButton(btn_open)
            box.exec_()
            if box.clickedButton() is btn_open:
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))
        except Exception as e:
            self._log(f"Could not show results-folder dialog: {e}")

    # ── Click handlers ─────────────────────────────────────────────────
    def _on_run_dd_single(self):
        targets = self._selected_targets()
        if not targets:
            return
        try:
            ok = self.run_death_dynamics_for(targets, interactive=True)
            if ok:
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                if len(targets) == 1:
                    folder = out_dir / "analysis" / targets[0] / "results"
                else:
                    folder = out_dir / "analysis"
                self._offer_open_results_folder(folder, what="Death Dynamics")
        finally:
            self._notify_results_changed()

    def _on_run_dd_combined(self):
        targets = self._selected_targets()
        if len(targets) < 2:
            return
        try:
            ok = self.run_multi_organoid_death_for(targets, interactive=True)
            if ok:
                folder = (
                    Path(self.metadata_loader.output_dir).expanduser()
                    / "analysis"
                    / "multi_organoid_comparison"
                )
                self._offer_open_results_folder(
                    folder, what="Combined Death Dynamics"
                )
        finally:
            self._notify_results_changed()

    def _on_run_ia_single(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        if not targets or not interactions:
            return
        self._persist_advanced()
        try:
            self.run_interaction_for(
                targets,
                interactions,
                group_by_line_condition=self.check_group_by_lc.isChecked(),
                interactive=True,
            )
        finally:
            self._notify_results_changed()

    def _on_run_ia_combined(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        if not targets or not interactions:
            return
        self._persist_advanced()
        try:
            self.run_multi_organoid_interaction_for(
                targets,
                interactions,
                time_window_min=int(self.spin_time_window.value()),
                group_by=self.combo_group_by.currentData() or "organoid_type",
                annotate_line_condition=self.check_annotate_lc.isChecked(),
                analysis_period_t=_period_from_t_radios(
                    self.period_radio_group, self.spin_period_start_t, self.spin_period_end_t,
                ),
                interactive=True,
            )
        finally:
            self._notify_results_changed()

    def _on_run_all_clicked(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        if not targets:
            self._log("Select at least one target.")
            return
        has_dead = _has_dead_channel(self.metadata_loader.metadata)
        dd_any_ok = False
        if has_dead:
            dd_any_ok |= bool(self.run_death_dynamics_for(targets, interactive=True))
            if len(targets) >= 2:
                dd_any_ok |= bool(
                    self.run_multi_organoid_death_for(targets, interactive=True)
                )
        if interactions:
            self._persist_advanced()
            self.run_interaction_for(
                targets,
                interactions,
                group_by_line_condition=self.check_group_by_lc.isChecked(),
                interactive=True,
            )
            self.run_multi_organoid_interaction_for(
                targets,
                interactions,
                time_window_min=int(self.spin_time_window.value()),
                group_by=self.combo_group_by.currentData() or "organoid_type",
                annotate_line_condition=self.check_annotate_lc.isChecked(),
                analysis_period_t=_period_from_t_radios(
                    self.period_radio_group, self.spin_period_start_t, self.spin_period_end_t,
                ),
                interactive=True,
            )
        if dd_any_ok:
            folder = Path(self.metadata_loader.output_dir).expanduser() / "analysis"
            self._offer_open_results_folder(folder, what="Death Dynamics")
        self._notify_results_changed()

    # ── Public run methods (called by queue) ───────────────────────────
    def run_death_dynamics_for(
        self,
        cell_types: Iterable[str],
        *,
        interactive: bool = True,
    ) -> bool:
        """Run per-target Death Dynamics. Returns True if at least one
        target completed successfully.

        Skips any target whose filtered CSV is missing or lacks the
        ``dead`` column with a friendly log message, so the queue can
        continue with the next step instead of erroring out when
        upstream Filter / Feature Extraction steps have not produced the
        expected outputs yet.
        """
        from behav3d.analysis.organoid_analysis import run_organoid_analysis

        cts = list(cell_types)
        if not cts:
            return False
        if not _has_dead_channel(self.metadata_loader.metadata):
            self._log("⚠ Cannot run Death Dynamics: no dead channel configured.")
            return False
        out_dir_path = Path(self.metadata_loader.output_dir).expanduser()
        out_dir = str(out_dir_path)
        selected_group_cols = [
            item.text() for item in self.list_dd_group_cols.selectedItems()
        ] or None
        any_ok = False
        for ct in cts:
            csv = _filtered_csv(out_dir_path, ct)
            if not csv.exists():
                self._log(
                    f"⏭ Skipping Death Dynamics for {ct}: filtered CSV "
                    f"not found ({csv.name}). Run Filtering for this "
                    f"cell type first."
                )
                continue
            if not _csv_has_column(csv, "dead"):
                self._log(
                    f"⏭ Skipping Death Dynamics for {ct}: 'dead' column "
                    f"missing in {csv.name}. Re-run Feature Extraction "
                    f"with a Dead-mask % threshold > 0."
                )
                continue
            try:
                self._log(f"Death Dynamics for {ct}…")
                print(f"\n▶ Death Dynamics: {ct}", file=sys.stderr)
                run_organoid_analysis(
                    output_dir=out_dir,
                    df_tracks_path=None,
                    org_type=ct,
                    metadata=self.metadata_loader.metadata,
                    group_cols=selected_group_cols,
                    show_in_notebook=False,
                )
                self._log(f"✅ Death Dynamics complete for {ct}.")
                any_ok = True
            except Exception as e:
                traceback.print_exc()
                self._log(f"❌ Death Dynamics failed for {ct}: {e}")
                if not interactive:
                    raise
        return any_ok

    def run_multi_organoid_death_for(
        self,
        cell_types: Iterable[str],
        *,
        interactive: bool = True,
    ) -> bool:
        """Run Combined Death Dynamics. Returns True on success.

        Pre-filters *cell_types* to those whose filtered CSV exists and
        contains a ``dead`` column. If fewer than 2 usable targets
        remain we log a skip message and return ``False`` instead of
        raising — that way the queue can continue with the next step
        when an earlier Filter / Feature Extraction step has not yet
        produced the expected outputs.
        """
        from behav3d.analysis.organoid_analysis import (
            plot_multi_organoid_death_dynamics,
        )

        cts = list(cell_types)
        if len(cts) < 2:
            self._log("⚠ Need at least 2 targets for combined death dynamics.")
            return False
        if not _has_dead_channel(self.metadata_loader.metadata):
            self._log(
                "⚠ Cannot run Combined Death Dynamics: no dead channel configured."
            )
            return False
        out_dir_path = Path(self.metadata_loader.output_dir).expanduser()
        out_dir = str(out_dir_path)

        # Pre-flight: only keep targets that have a usable filtered CSV
        # with a 'dead' column.
        usable: list[str] = []
        skipped: list[str] = []
        for ct in cts:
            csv = _filtered_csv(out_dir_path, ct)
            if csv.exists() and _csv_has_column(csv, "dead"):
                usable.append(ct)
            else:
                skipped.append(ct)
        if skipped:
            self._log(
                f"⏭ Combined Death Dynamics: skipping {', '.join(skipped)} "
                f"(missing filtered CSV or 'dead' column)."
            )
        if len(usable) < 2:
            self._log(
                f"⏭ Skipping Combined Death Dynamics: only "
                f"{len(usable)} target(s) with usable data — need at "
                f"least 2."
            )
            return False

        # Build threshold map from feature-extraction config
        params = self.metadata_loader.behav3d_parameters or {}
        feat = params.get("features", {}) or {}
        thr_map = {}
        for ct in usable:
            v = feat.get(ct, {}).get("dead_mask_percentage_threshold")
            if v is not None:
                try:
                    thr_map[ct] = float(v)
                except Exception:
                    pass
        try:
            self._log(f"Combined Death Dynamics for {', '.join(usable)}…")
            plot_multi_organoid_death_dynamics(
                output_dir=out_dir,
                organoid_types=usable,
                dead_perc_threshold_map=thr_map,
                show_in_notebook=False,
            )
            self._log("✅ Combined Death Dynamics complete.")
            return True
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Combined Death Dynamics failed: {e}")
            if not interactive:
                raise
            return False

    def run_interaction_for(
        self,
        cell_types: Iterable[str],
        interaction_cell_types: Iterable[str],
        *,
        group_by_line_condition: bool = False,
        interactive: bool = True,
    ):
        from behav3d.analysis.interaction_analysis import run_interaction_analysis

        cts = list(cell_types)
        ints = list(interaction_cell_types)
        if not cts or not ints:
            return
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        for ct in cts:
            csv = _filtered_csv(Path(out_dir), ct)
            if not csv.exists():
                self._log(f"⚠ Skipping {ct}: filtered CSV not found ({csv}).")
                continue
            available = _contact_cols_in_csv(csv)
            selected = [c for c in ints if c in available]
            if not selected:
                self._log(
                    f"⚠ Skipping {ct}: none of the selected interaction cell "
                    f"types have a contact column in {csv.name}."
                )
                continue
            try:
                self._log(
                    f"Interaction Analysis: {ct} ↔ {', '.join(selected)}…"
                )
                run_interaction_analysis(
                    output_dir=out_dir,
                    cell_type=ct,
                    interacting_cell_types=selected,
                    df_tracks_path=str(csv),
                    show_plots=interactive,
                    group_by_line_condition=group_by_line_condition,
                )
                self._log(f"✅ Interaction Analysis complete for {ct}.")
            except Exception as e:
                traceback.print_exc()
                self._log(f"❌ Interaction Analysis failed for {ct}: {e}")
                if not interactive:
                    raise

    def run_multi_organoid_interaction_for(
        self,
        cell_types: Iterable[str],
        interaction_cell_types: Iterable[str],
        *,
        time_window_min: float = 60.0,
        group_by: str = "organoid_type",
        annotate_line_condition: bool = False,
        analysis_period_t: tuple = None,
        interactive: bool = True,
    ):
        """Run Interaction Overview.

        Pre-filters *cell_types* to those whose filtered CSV exists and
        carries a contact column for at least one of the selected
        interaction cell types. If no usable targets remain we log a skip
        message and return without raising, so the queue can continue when
        an earlier Filter step has not yet produced the expected outputs.
        """
        from behav3d.analysis.interaction_analysis import (
            run_multi_organoid_interaction_comparison,
        )

        cts = list(cell_types)
        ints = list(interaction_cell_types)
        if not cts or not ints:
            return
        if not _has_dead_channel(self.metadata_loader.metadata):
            self._log(
                "ⓘ Running Interaction Overview without dead "
                "data — alive/dead split, cumulative-to-death curves and "
                "active-killing dashboard will be skipped."
            )
        out_dir_path = Path(self.metadata_loader.output_dir).expanduser()
        out_dir = str(out_dir_path)

        # Pre-flight: only keep targets whose filtered CSV exists AND
        # carries a contact column for at least one selected interaction
        # cell type.
        usable: list[str] = []
        skipped: list[str] = []
        for ct in cts:
            csv = _filtered_csv(out_dir_path, ct)
            if not csv.exists():
                skipped.append(f"{ct} (no filtered CSV)")
                continue
            contacts = _contact_cols_in_csv(csv)
            if not any(it in contacts for it in ints):
                skipped.append(f"{ct} (no matching contact columns)")
                continue
            usable.append(ct)
        if skipped:
            self._log(
                "⏭ Interaction Overview: skipping "
                f"{', '.join(skipped)}."
            )
        if len(usable) < 1:
            self._log(
                f"⏭ Skipping Interaction Overview: no targets with "
                f"usable data."
            )
            return

        try:
            self._log(
                f"Interaction Overview: {', '.join(usable)} "
                f"× {', '.join(ints)}…"
            )
            run_multi_organoid_interaction_comparison(
                output_dir=out_dir,
                organoid_types=usable,
                immune_types=ints,
                metadata=self.metadata_loader.metadata,
                group_by=group_by,
                time_window_min=float(time_window_min),
                analysis_period_t=analysis_period_t,
                annotate_line_condition=annotate_line_condition,
                show_plots=interactive,
            )
            self._log("✅ Interaction Overview complete.")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Interaction Overview failed: {e}")
            if not interactive:
                raise

    # ── Invasiveness Analysis (Step 3, immune-cell perspective) ─────────
    def _selected_invasiveness_immune(self) -> list[str]:
        return [
            im for im, cb in getattr(self, "_inv_immune_checks", {}).items()
            if cb.isChecked()
        ]

    def _selected_invasiveness_targets(self) -> list[str]:
        return [
            t for t, cb in getattr(self, "_inv_target_checks", {}).items()
            if cb.isChecked()
        ]

    def _rebuild_invasiveness(self):
        """(Re)populate the immune-cell checkboxes, then rebuild target
        checkboxes. Additive: called whenever metadata changes."""
        self._inv_immune_checks = {}
        self._clear_layout(self.inv_immune_lay)
        _, immune_types, other_types = _detect_cell_types(self.metadata_loader)
        candidates = list(immune_types) + list(other_types)
        if not candidates:
            self.inv_immune_lay.addWidget(self._inv_immune_placeholder)
        else:
            for i, ct in enumerate(candidates):
                cb = QCheckBox(ct)
                cb.setChecked(i == 0)
                cb.toggled.connect(
                    lambda _checked=False: self._rebuild_invasiveness_targets()
                )
                self.inv_immune_lay.addWidget(cb)
                self._inv_immune_checks[ct] = cb
        self._rebuild_invasiveness_targets()

    def _rebuild_invasiveness_targets(self):
        self._inv_target_checks = {}
        self._clear_layout(self.inv_targets_lay)

        immune_list = self._selected_invasiveness_immune()
        if (
            self.metadata_loader is None
            or self.metadata_loader.metadata is None
            or not self.metadata_loader.output_dir
            or not immune_list
        ):
            self.inv_targets_lay.addWidget(self._inv_targets_placeholder)
            self.btn_run_inv.setEnabled(False)
            self.btn_queue_inv.setEnabled(False)
            self._update_view_buttons()
            return

        # Union of targets across all selected immune types; note any whose
        # filtered CSV is missing.
        out_dir = Path(self.metadata_loader.output_dir)
        targets: list[str] = []
        missing: list[str] = []
        for im in immune_list:
            csv = _filtered_csv(out_dir, im)
            tg = _invasiveness_targets_in_csv(csv)
            if not tg:
                missing.append(im)
            for t in tg:
                if t not in targets:
                    targets.append(t)
        if not targets:
            msg = QLabel(
                f"⚠ No invasiveness columns for: {', '.join(missing or immune_list)}. "
                "Run Feature Extraction with 'invasiveness' enabled, then Filtering."
            )
            msg.setWordWrap(True)
            msg.setStyleSheet("color: #888; font-style: italic;")
            self.inv_targets_lay.addWidget(msg)
            self.btn_run_inv.setEnabled(False)
            self.btn_queue_inv.setEnabled(False)
            self._update_view_buttons()
            return

        for t in sorted(targets):
            cb = QCheckBox(t)
            cb.setChecked(True)
            cb.toggled.connect(
                lambda _checked=False: self._refresh_invasiveness_buttons()
            )
            self.inv_targets_lay.addWidget(cb)
            self._inv_target_checks[t] = cb
        self._refresh_invasiveness_buttons()

    def _refresh_invasiveness_buttons(self):
        immune_list = self._selected_invasiveness_immune()
        ready = bool(immune_list) and bool(self._selected_invasiveness_targets())
        self.btn_run_inv.setEnabled(ready)
        self.btn_queue_inv.setEnabled(ready)
        if not immune_list:
            tip = "Select at least one immune cell type."
        elif not self._selected_invasiveness_targets():
            tip = "Select at least one target to compare."
        else:
            tip = ""
        self.btn_run_inv.setToolTip(tip)
        self.btn_queue_inv.setToolTip(tip)
        self._update_view_buttons()

    def _on_run_invasiveness(self):
        immune_list = self._selected_invasiveness_immune()
        targets = self._selected_invasiveness_targets()
        if not immune_list or not targets:
            return
        try:
            ok = self.run_invasiveness_for(
                immune_list,
                targets,
                summary_stat=self.inv_summary_combo.currentData() or "mean",
                group_by_line_condition=self.check_inv_group_by_lc.isChecked(),
                analysis_period_t=_period_from_t_radios(
                    self.inv_period_radio_group, self.spin_inv_start_t, self.spin_inv_end_t,
                ),
                interactive=True,
            )
            if ok:
                base = Path(self.metadata_loader.output_dir).expanduser() / "analysis"
                folder = (
                    base / "invasiveness_analysis" if len(immune_list) > 1
                    else base / immune_list[0] / "invasiveness_analysis"
                )
                self._offer_open_results_folder(
                    folder, what="Invasiveness Analysis"
                )
        finally:
            self._notify_results_changed()

    def run_invasiveness_for(
        self,
        immune_cell_types,
        targets: Iterable[str],
        *,
        summary_stat: str = "mean",
        group_by_line_condition: bool = False,
        analysis_period_t: tuple = None,
        interactive: bool = True,
    ) -> bool:
        """Run invasiveness analysis for one or more immune cell types vs
        targets. ``immune_cell_types`` may be a single string or a list.

        Returns True on success. Skips gracefully (logs) when none of the
        immune cell types' filtered CSVs exist, so the queue can continue.
        """
        from behav3d.analysis.invasiveness_analysis import (
            run_invasiveness_analysis,
        )

        if isinstance(immune_cell_types, str):
            immune_list = [immune_cell_types]
        else:
            immune_list = [str(c) for c in immune_cell_types if c]
        tgts = list(targets)
        if not immune_list or not tgts:
            return False
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        present = [im for im in immune_list if _filtered_csv(Path(out_dir), im).exists()]
        if not present:
            self._log(
                f"⚠ Skipping invasiveness: no filtered CSV found for "
                f"{', '.join(immune_list)}."
            )
            return False
        try:
            self._log(
                f"Invasiveness Analysis: {', '.join(present)} → {', '.join(tgts)}…"
            )
            organoid_types, _, _ = _detect_cell_types(self.metadata_loader)
            run_invasiveness_analysis(
                output_dir=out_dir,
                immune_cell_types=present,
                targets=tgts,
                summary_stat=summary_stat,
                show_plots=interactive,
                group_by_line_condition=group_by_line_condition,
                analysis_period_t=analysis_period_t,
                organoid_types=list(organoid_types),
            )
            self._log(f"✅ Invasiveness Analysis complete for {', '.join(present)}.")
            return True
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Invasiveness Analysis failed for {', '.join(present)}: {e}")
            if not interactive:
                raise
            return False


# ═══════════════════════════════════════════════════════════════════════════
# Single Cell sub-tab — real implementation imported from _single_cell.py
# ═══════════════════════════════════════════════════════════════════════════
from behav3d.napari._single_cell import SingleCellTab  # noqa: E402  (re-exported)


# NOTE: ResultsPanel used to be defined here; it now lives in
# behav3d/napari/_results_panel.py so the Filtering and Feature
# Extraction tabs can embed an instance without importing the
# (much heavier) Analysis module.

# ═══════════════════════════════════════════════════════════════════════════
# Outer Analysis tab
# ═══════════════════════════════════════════════════════════════════════════
class AnalysisTab(QWidget):
    """Outer tab wrapping Death Dynamics + (placeholder) Single Cell.

    Hosts a shared :class:`ResultsPanel` under a vertical splitter so the
    panel is visible from both inner sub-tabs.
    """

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # The whole Analysis GUI is gated behind a "metadata loaded" check,
        # mirroring the Segmentation / Filtering / Tracking / Feature
        # Extraction tabs. Until metadata is loaded we show a single
        # centered placeholder; after that the inner tabs + Results panel
        # are revealed.
        self.stack = QStackedWidget(self)
        outer.addWidget(self.stack)

        # -- Page 0: Placeholder --
        self.placeholder_page = QWidget()
        place_lay = QVBoxLayout(self.placeholder_page)
        place_lay.setContentsMargins(20, 20, 20, 20)
        place_lay.setAlignment(Qt.AlignCenter)
        self.placeholder_label = QLabel(
            "Load metadata in the Data Preparation tab to see analysis options."
        )
        self.placeholder_label.setAlignment(Qt.AlignCenter)
        self.placeholder_label.setWordWrap(True)
        self.placeholder_label.setStyleSheet(
            "color: #888; font-style: italic; font-size: 14px;"
        )
        place_lay.addWidget(self.placeholder_label)
        self.stack.addWidget(self.placeholder_page)

        # -- Page 1: Main Content --
        self.main_content = QWidget()
        main_lay = QVBoxLayout(self.main_content)
        main_lay.setContentsMargins(0, 0, 0, 0)
        main_lay.setSpacing(0)

        # ── Cell type grouping ──────────────────────────────────────────
        # Sits above both subtabs: grouping merges already-*filtered*
        # populations (see behav3d.analysis.grouping), so it belongs after
        # Filtering and before Death Dynamics / Single Cell, not before
        # Feature Extraction — that would duplicate feature
        # extraction/filtering work for every group.
        grouping_row = QHBoxLayout()
        grouping_row.setContentsMargins(6, 6, 6, 4)
        self.btn_open_grouping = QPushButton("🧬  Group cells together for analysis…")
        self.btn_open_grouping.setToolTip(
            "Merge several already-filtered cell types into a single "
            "'{name}_merged' population for Death Dynamics / Single Cell.\n"
            "Metadata is not modified — only behav3d_parameters.yml and a "
            "merged filtered CSV are written."
        )
        self.btn_open_grouping.setStyleSheet(
            "QPushButton { background:#1f3a5f; color:#bbdefb; font-weight:bold; "
            "border:1px solid #64b5f6; border-radius:4px; padding:8px 14px; } "
            "QPushButton:hover { background:#2c4f7f; }"
            "QPushButton:disabled { color:#666; border-color:#444; }"
        )
        self.btn_open_grouping.clicked.connect(self._on_open_grouping_dialog)
        self.btn_open_grouping.setEnabled(False)
        grouping_row.addWidget(self.btn_open_grouping)
        grouping_row.addWidget(HelpButton(
            "Group cells together for analysis",
            "Combine several already-filtered cell types (e.g. CD4 T cells + "
            "CD8 T cells) into one merged population you can run Death "
            "Dynamics / Single Cell on as a single group.\n\n"
            "Each cell keeps its original type too — a new 'origin_cell_type' "
            "column records which type it came from, so you can still split "
            "or colour results by the original identity afterwards.\n\n"
            "Your metadata.csv is never modified; only a new merged results "
            "file is created.",
        ))
        grouping_row.addStretch()
        main_lay.addLayout(grouping_row)

        self.splitter = QSplitter(Qt.Vertical, self.main_content)
        self.splitter.setChildrenCollapsible(True)
        main_lay.addWidget(self.splitter)

        self.inner_tabs = QTabWidget()
        self.splitter.addWidget(self.inner_tabs)

        self.death_dynamics_tab = DeathDynamicsTab(
            viewer=viewer, metadata_loader=metadata_loader, parent=self
        )
        self.inner_tabs.addTab(self.death_dynamics_tab, "💀 Death Dynamics")

        self.single_cell_tab = SingleCellTab(
            viewer=viewer, metadata_loader=metadata_loader, parent=self
        )
        self.inner_tabs.addTab(self.single_cell_tab, "🧬 Single Cell")

        self.results_panel = ResultsPanel(
            viewer=viewer, metadata_loader=metadata_loader, parent=self
        )
        self.splitter.addWidget(self.results_panel)

        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, True)
        self.splitter.setSizes([600, 400])

        # Re-scan when the user switches sub-tab so new outputs show up
        # without an explicit refresh click.
        self.inner_tabs.currentChanged.connect(
            lambda _i: self.results_panel.refresh()
        )

        self.stack.addWidget(self.main_content)

        # Wire metadata signal and initialise the visible page.
        if metadata_loader is not None and hasattr(
            metadata_loader, "metadata_loaded"
        ):
            metadata_loader.metadata_loaded.connect(self._on_metadata_loaded)

        if (
            metadata_loader is not None
            and getattr(metadata_loader, "metadata", None) is not None
        ):
            self.stack.setCurrentIndex(1)
            self.btn_open_grouping.setEnabled(True)
        else:
            self.stack.setCurrentIndex(0)

    def _on_metadata_loaded(self, *_):
        """Reveal the analysis GUI once metadata is available."""
        self.stack.setCurrentIndex(1)
        self.btn_open_grouping.setEnabled(
            self.metadata_loader is not None
            and getattr(self.metadata_loader, "metadata", None) is not None
        )

    def _on_metadata_updated(self, *_):
        """Cascade metadata updates to inner tabs and results panel."""
        if hasattr(self, "death_dynamics_tab"):
            self.death_dynamics_tab._on_metadata_updated()
        if hasattr(self, "single_cell_tab"):
            self.single_cell_tab._on_metadata_updated()
        if hasattr(self, "results_panel"):
            self.results_panel.refresh()

    # ── Cell type grouping ────────────────────────────────────────────────
    def _on_open_grouping_dialog(self):
        """Open the cell-type grouping dialog and refresh cell-type
        dropdowns (Death Dynamics / Single Cell) on success."""
        from behav3d.napari._grouping_dialog import GroupBuilderDialog

        if self.metadata_loader is None or self.metadata_loader.metadata is None:
            return

        dlg = GroupBuilderDialog(self.metadata_loader, parent=self)
        dlg.group_created.connect(self._on_group_changed)
        dlg.group_removed.connect(self._on_group_changed)
        dlg.exec_()

    def _on_group_changed(self, group_id: str):
        """Slot: refresh Death Dynamics / Single Cell dropdowns after a
        group is created or removed."""
        self._on_metadata_updated()
