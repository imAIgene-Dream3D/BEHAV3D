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
    QScrollArea, QSpinBox, QComboBox, QFrame, QSizePolicy,
    QToolButton, QMessageBox,
)
from qtpy.QtCore import Qt, QUrl
from qtpy.QtGui import QDesktopServices


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════
def _detect_cell_types(metadata_loader):
    """Return (organoid_types, immune_types, other_types) lists.

    Includes merged/grouped outputs routed to their detected category.
    Returns three empty lists when metadata is missing.
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
        self._content.setVisible(expanded)
        outer.addWidget(self._content)

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
    - Run All Selected + log
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
        for btn in (
            self.btn_queue_dd_single,
            self.btn_queue_dd_combined,
            self.btn_queue_ia_single,
            self.btn_queue_ia_combined,
        ):
            btn.setFixedSize(36, 28)
            btn.setToolTip("Add to Processing Queue")
            btn.setStyleSheet(
                "QPushButton { background: #1a1a2e; color: #ffc107; "
                "border: 1px solid #ffc107; border-radius: 4px; "
                "font-size: 11px; font-weight: bold; }"
                "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
            )

        self._init_ui()

        if metadata_loader is not None and hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    # ── UI ──────────────────────────────────────────────────────────────
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

        dd_btn_row = QHBoxLayout()
        self.btn_run_dd_single = QPushButton("▶ Run Death Dynamics (per target)")
        self._style_primary(self.btn_run_dd_single)
        self.btn_run_dd_single.clicked.connect(self._on_run_dd_single)
        dd_btn_row.addWidget(self.btn_run_dd_single, stretch=1)
        dd_btn_row.addWidget(self.btn_queue_dd_single)
        dd_lay.addLayout(dd_btn_row)

        dd_combined_row = QHBoxLayout()
        self.btn_run_dd_combined = QPushButton(
            "▶ Run Combined Death Dynamics (≥2 targets)"
        )
        self._style_primary(self.btn_run_dd_combined)
        self.btn_run_dd_combined.clicked.connect(self._on_run_dd_combined)
        dd_combined_row.addWidget(self.btn_run_dd_combined, stretch=1)
        dd_combined_row.addWidget(self.btn_queue_dd_combined)
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
            "  • Cumulative-to-death curves (Combined run)\n"
            "  • Active Killing dashboard (Combined run)"
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

        ia_btn_row = QHBoxLayout()
        self.btn_run_ia_single = QPushButton("▶ Run Interaction Analysis (per target)")
        self._style_primary(self.btn_run_ia_single)
        self.btn_run_ia_single.clicked.connect(self._on_run_ia_single)
        ia_btn_row.addWidget(self.btn_run_ia_single, stretch=1)
        ia_btn_row.addWidget(self.btn_queue_ia_single)
        ia_lay.addLayout(ia_btn_row)

        ia_combined_row = QHBoxLayout()
        self.btn_run_ia_combined = QPushButton(
            "▶ Run Combined Interaction Comparison (≥2 targets)"
        )
        self._style_primary(self.btn_run_ia_combined)
        self.btn_run_ia_combined.clicked.connect(self._on_run_ia_combined)
        ia_combined_row.addWidget(self.btn_run_ia_combined, stretch=1)
        ia_combined_row.addWidget(self.btn_queue_ia_combined)
        ia_lay.addLayout(ia_combined_row)

        # ── Advanced settings (Interaction Analysis) ─────────────
        # Genuinely collapsible (arrow toggle, no checkbox).
        self.ia_adv = CollapsibleSection(
            "Advanced settings (Interaction Analysis)", expanded=False
        )
        ia_adv_form = QFormLayout()
        ia_adv_form.setContentsMargins(0, 0, 0, 0)
        ia_adv_form.setSpacing(4)

        self.spin_time_window = QSpinBox()
        self.spin_time_window.setRange(1, 99999)
        self.spin_time_window.setValue(60)
        self.spin_time_window.setSuffix(" min")
        self.spin_time_window.setMaximumWidth(100)
        ia_adv_form.addRow("Time window before TOD:", self.spin_time_window)

        self.combo_group_by = QComboBox()
        self.combo_group_by.addItem("By target (organoid) type", userData="organoid_type")
        self.combo_group_by.addItem("By treatment (immune cell)", userData="treatment")
        self.combo_group_by.setMaximumWidth(280)
        ia_adv_form.addRow("Combined group by:", self.combo_group_by)

        self.ia_adv.addLayout(ia_adv_form)
        ia_lay.addWidget(self.ia_adv)

        layout.addWidget(self.ia_group)

        # ── Run All Selected + Log ───────────────────────────────
        self.btn_run_all = QPushButton("▶▶  Run All Selected")
        self.btn_run_all.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_all.clicked.connect(self._on_run_all_clicked)
        layout.addWidget(self.btn_run_all)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(140)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        layout.addStretch()

        # Try a first build (in case metadata is already loaded)
        self._rebuild()

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
            self._refresh_button_states()
            return

        organoid_types, immune_types, other_types = _detect_cell_types(
            self.metadata_loader
        )
        out_dir = Path(self.metadata_loader.output_dir)
        md = self.metadata_loader.metadata
        has_dead = _has_dead_channel(md)

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
                    f"{float(thr):.4g} %" if thr not in (None, 0, 0.0)
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
        ia_ready_combined = ia_ready_single and (len(targets) >= 2)
        self.btn_run_ia_combined.setEnabled(ia_ready_combined)
        self.btn_queue_ia_combined.setEnabled(ia_ready_combined)
        tip_combined = ""
        if len(targets) < 2:
            tip_combined = "Select at least 2 targets to compare."
        elif not interactions:
            tip_combined = "Select at least one interaction cell type."
        elif not has_dead:
            tip_combined = (
                "Will run without dead data — alive/dead comparison panels "
                "and cumulative-to-death curves will be skipped."
            )
        self.btn_run_ia_combined.setToolTip(tip_combined)
        self.btn_queue_ia_combined.setToolTip(tip_combined)

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
        ok = self.run_death_dynamics_for(targets, interactive=True)
        if ok:
            out_dir = Path(self.metadata_loader.output_dir).expanduser()
            if len(targets) == 1:
                folder = out_dir / "analysis" / targets[0] / "results"
            else:
                folder = out_dir / "analysis"
            self._offer_open_results_folder(folder, what="Death Dynamics")

    def _on_run_dd_combined(self):
        targets = self._selected_targets()
        if len(targets) < 2:
            return
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

    def _on_run_ia_single(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        if not targets or not interactions:
            return
        self.run_interaction_for(targets, interactions, interactive=True)

    def _on_run_ia_combined(self):
        targets = self._selected_targets()
        interactions = self._selected_interactions()
        if len(targets) < 2 or not interactions:
            return
        self._persist_advanced()
        self.run_multi_organoid_interaction_for(
            targets,
            interactions,
            time_window_min=int(self.spin_time_window.value()),
            group_by=self.combo_group_by.currentData() or "organoid_type",
            interactive=True,
        )

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
            self.run_interaction_for(targets, interactions, interactive=True)
            if len(targets) >= 2:
                self._persist_advanced()
                self.run_multi_organoid_interaction_for(
                    targets,
                    interactions,
                    time_window_min=int(self.spin_time_window.value()),
                    group_by=self.combo_group_by.currentData() or "organoid_type",
                    interactive=True,
                )
        if dd_any_ok:
            folder = Path(self.metadata_loader.output_dir).expanduser() / "analysis"
            self._offer_open_results_folder(folder, what="Death Dynamics")

    # ── Public run methods (called by queue) ───────────────────────────
    def run_death_dynamics_for(
        self,
        cell_types: Iterable[str],
        *,
        interactive: bool = True,
    ) -> bool:
        """Run per-target Death Dynamics. Returns True if at least one
        target completed successfully."""
        from behav3d.analysis.organoid_analysis import run_organoid_analysis

        cts = list(cell_types)
        if not cts:
            return False
        if not _has_dead_channel(self.metadata_loader.metadata):
            self._log("⚠ Cannot run Death Dynamics: no dead channel configured.")
            return False
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        any_ok = False
        for ct in cts:
            try:
                self._log(f"Death Dynamics for {ct}…")
                print(f"\n▶ Death Dynamics: {ct}", file=sys.stderr)
                run_organoid_analysis(
                    output_dir=out_dir,
                    df_tracks_path=None,
                    org_type=ct,
                    metadata=self.metadata_loader.metadata,
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
        """Run Combined Death Dynamics. Returns True on success."""
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
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        # Build threshold map from feature-extraction config
        params = self.metadata_loader.behav3d_parameters or {}
        feat = params.get("features", {}) or {}
        thr_map = {}
        for ct in cts:
            v = feat.get(ct, {}).get("dead_mask_percentage_threshold")
            if v is not None:
                try:
                    thr_map[ct] = float(v)
                except Exception:
                    pass
        try:
            self._log(f"Combined Death Dynamics for {', '.join(cts)}…")
            plot_multi_organoid_death_dynamics(
                output_dir=out_dir,
                organoid_types=cts,
                dead_perc_threshold_map=thr_map,
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
        interactive: bool = True,
    ):
        from behav3d.analysis.interaction_analysis import (
            run_multi_organoid_interaction_comparison,
        )

        cts = list(cell_types)
        ints = list(interaction_cell_types)
        if len(cts) < 2 or not ints:
            return
        if not _has_dead_channel(self.metadata_loader.metadata):
            self._log(
                "ⓘ Running Combined Interaction Comparison without dead "
                "data — alive/dead split, cumulative-to-death curves and "
                "active-killing dashboard will be skipped."
            )
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        try:
            self._log(
                f"Combined Interaction Comparison: {', '.join(cts)} "
                f"× {', '.join(ints)}…"
            )
            run_multi_organoid_interaction_comparison(
                output_dir=out_dir,
                organoid_types=cts,
                immune_types=ints,
                metadata=self.metadata_loader.metadata,
                group_by=group_by,
                time_window_min=float(time_window_min),
                show_plots=interactive,
            )
            self._log("✅ Combined Interaction Comparison complete.")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Combined Interaction Comparison failed: {e}")
            if not interactive:
                raise


# ═══════════════════════════════════════════════════════════════════════════
# Single Cell sub-tab (placeholder)
# ═══════════════════════════════════════════════════════════════════════════
class SingleCellTab(QWidget):
    """Placeholder for the upcoming Single Cell analysis sub-tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setAlignment(Qt.AlignCenter)

        icon = QLabel("🧬")
        icon.setStyleSheet("font-size: 48px;")
        icon.setAlignment(Qt.AlignCenter)
        lay.addWidget(icon)

        title = QLabel("Single Cell")
        title.setStyleSheet("font-size: 18px; font-weight: bold; color: #aaa;")
        title.setAlignment(Qt.AlignCenter)
        lay.addWidget(title)

        sub = QLabel(
            "Coming soon — DTW/UMAP/clustering analyses for motile cell "
            "types will live here."
        )
        sub.setStyleSheet("color: #888; font-size: 11px;")
        sub.setAlignment(Qt.AlignCenter)
        sub.setWordWrap(True)
        lay.addWidget(sub)


# ═══════════════════════════════════════════════════════════════════════════
# Outer Analysis tab
# ═══════════════════════════════════════════════════════════════════════════
class AnalysisTab(QWidget):
    """Outer tab wrapping Death Dynamics + (placeholder) Single Cell."""

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self.inner_tabs = QTabWidget()
        outer.addWidget(self.inner_tabs)

        self.death_dynamics_tab = DeathDynamicsTab(
            viewer=viewer, metadata_loader=metadata_loader, parent=self
        )
        self.inner_tabs.addTab(self.death_dynamics_tab, "💀 Death Dynamics")

        self.single_cell_tab = SingleCellTab(parent=self)
        self.inner_tabs.addTab(self.single_cell_tab, "🧬 Single Cell")
