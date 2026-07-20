"""
BEHAV3D napari plugin – Rename Cluster Dialog.

Provides a modal QDialog for renaming behavioral state clusters produced
by state classification (intrinsic / full binary-group) and track
classification steps.  Three modes share a single implementation:

- ``"intrinsic"``  — rename primary dynamic state clusters
                     (obs column ``intrinsic_behavioral_cluster``)
- ``"full"``       — rename clusters produced by crossing the intrinsic
                     state with binary group values, e.g. organoid contact
                     (obs column ``full_behavioral_cluster``).
                     Adds checkbox selection + combine-to-one-label feature.
- ``"track"``      — rename DTW trajectory clusters
                     (obs column ``behavioral_trajectory_cluster`` or
                     ``dtaidistance_cluster``, whichever is present)

All modes:
  • Save   → applies ``relabel_cluster_ids()`` in-place on the adata,
             writes the adata to disk, emits ``clusters_renamed`` signal.
  • Discard → ``reject()`` — no side effects.
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Dict, Optional

from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _unique_sorted(values) -> list:
    """Return unique non-None string values in sorted order."""
    seen = set()
    out = []
    for v in values:
        s = str(v)
        if s not in seen:
            seen.add(s)
            out.append(s)
    return sorted(out)


def _track_cluster_col(adata) -> Optional[str]:
    """Return the first recognised track-cluster obs column in *adata*, or None."""
    if adata is None:
        return None
    uns_key = adata.uns.get("dtai_trajectory_clustering", {}).get("cluster_key")
    if uns_key and uns_key in adata.obs.columns:
        return uns_key
    for col in (
        "ClusterID",
        "behavioral_trajectory_cluster",
        "dtaidistance_cluster",
        "track_cluster",
        "trajectory_cluster",
    ):
        if col in adata.obs.columns:
            return col
    return None


# ═══════════════════════════════════════════════════════════════════════════
# Dialog
# ═══════════════════════════════════════════════════════════════════════════

class RenameClusterDialog(QDialog):
    """Modal dialog for renaming behavioral clusters.

    Parameters
    ----------
    mode:
        ``"intrinsic"``, ``"full"``, or ``"track"``.
    model_adata:
        AnnData object whose obs column will be relabelled.
    adata_path:
        Path used to save the adata after renaming.
    parent:
        Optional parent widget.
    """

    # Emitted on Save with the mapping dict {old_name: new_name}.
    clusters_renamed = Signal(dict)

    _TITLES = {
        "intrinsic": "Rename Primary Dynamic State Clusters",
        "full":      "Rename Full Behavioral Clusters (Binary Groups)",
        "track":     "Rename Track Trajectory Clusters",
    }
    _SUBTITLES = {
        "intrinsic": (
            "Assign meaningful names to each primary movement-pattern "
            "cluster (derived from speed, displacement, sphericity, etc.)."
        ),
        "full": (
            "Rename clusters that combine the primary state with binary "
            "context labels (e.g. organoid contact = 0 or 1).  "
            "Use <b>Combine selected →</b> to merge multiple clusters into "
            "one label before saving."
        ),
        "track": (
            "Assign meaningful names to each trajectory cluster "
            "(derived from DTW distance on state sequences)."
        ),
    }
    _OBS_COLS = {
        "intrinsic": "intrinsic_behavioral_cluster",
        "full":      "full_behavioral_cluster",
    }

    def __init__(
        self,
        mode: str,
        model_adata,
        adata_path: Optional[Path] = None,
        parent=None,
    ):
        super().__init__(parent)
        if mode not in ("intrinsic", "full", "track"):
            raise ValueError(f"Unknown mode: {mode!r}")
        self._mode = mode
        self._adata = model_adata
        self._adata_path = Path(adata_path) if adata_path else None

        self._name_edits: Dict[str, QLineEdit] = {}
        self._select_checks: Dict[str, QCheckBox] = {}  # full mode only
        self._color_buttons: Dict[str, QPushButton] = {}

        self.setWindowTitle(self._TITLES[mode])
        self.setMinimumWidth(580)
        self.setMinimumHeight(380)
        self.setSizeGripEnabled(True)

        self._init_ui()
        self._populate()

    # ── UI ──────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setSpacing(8)
        outer.setContentsMargins(12, 12, 12, 12)

        # Header
        title_lbl = QLabel(self._TITLES[self._mode])
        title_lbl.setStyleSheet(
            "font-size: 14px; font-weight: bold; color: #ddd;"
        )
        title_lbl.setWordWrap(True)
        outer.addWidget(title_lbl)

        sub_lbl = QLabel(self._SUBTITLES[self._mode])
        sub_lbl.setStyleSheet("font-size: 11px; color: #999;")
        sub_lbl.setWordWrap(True)
        outer.addWidget(sub_lbl)

        # Status
        self._status_lbl = QLabel()
        self._status_lbl.setStyleSheet(
            "font-size: 11px; color: #aaa; font-style: italic;"
        )
        outer.addWidget(self._status_lbl)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("color: #444;")
        outer.addWidget(sep)

        # Column header
        hdr = QHBoxLayout()
        hdr.setContentsMargins(0, 0, 0, 0)
        if self._mode == "full":
            sel_hdr = QLabel("Sel")
            sel_hdr.setFixedWidth(30)
            sel_hdr.setStyleSheet("color:#888; font-size:10px;")
            hdr.addWidget(sel_hdr)
        old_hdr = QLabel("Current name")
        old_hdr.setStyleSheet("color:#888; font-size:10px; font-weight:bold;")
        hdr.addWidget(old_hdr, stretch=2)
        new_hdr = QLabel("New name")
        new_hdr.setStyleSheet("color:#888; font-size:10px; font-weight:bold;")
        hdr.addWidget(new_hdr, stretch=3)
        color_hdr = QLabel("Color")
        color_hdr.setFixedWidth(40)
        color_hdr.setStyleSheet("color:#888; font-size:10px; font-weight:bold;")
        hdr.addWidget(color_hdr)
        outer.addLayout(hdr)

        # Scrollable rows
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setStyleSheet("QScrollArea { border: none; }")

        self._rows_widget = QWidget()
        self._rows_layout = QVBoxLayout(self._rows_widget)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(3)
        self._rows_layout.addStretch(1)

        scroll.setWidget(self._rows_widget)
        outer.addWidget(scroll, stretch=1)

        # Combine row (full mode only)
        if self._mode == "full":
            combine_frame = QFrame()
            combine_frame.setFrameShape(QFrame.StyledPanel)
            combine_frame.setStyleSheet(
                "QFrame { background: #2a2a3a; border: 1px solid #444; "
                "border-radius: 3px; }"
            )
            combine_lay = QHBoxLayout(combine_frame)
            combine_lay.setContentsMargins(8, 4, 8, 4)
            combine_lay.setSpacing(6)
            combine_lbl = QLabel("Combine selected →")
            combine_lbl.setStyleSheet("color: #aaa; font-size: 11px;")
            combine_lay.addWidget(combine_lbl)
            self._combine_edit = QLineEdit()
            self._combine_edit.setPlaceholderText("New combined name…")
            self._combine_edit.setStyleSheet(
                "background: #1a1a2e; color: #ddd; border: 1px solid #555; "
                "border-radius: 2px; padding: 2px 4px;"
            )
            combine_lay.addWidget(self._combine_edit, stretch=1)
            self._btn_combine = QPushButton("Combine")
            self._btn_combine.setFixedHeight(26)
            self._btn_combine.setStyleSheet(
                "QPushButton { background: #2d5a8e; color: #ddd; "
                "border-radius: 3px; padding: 2px 10px; font-size: 11px; }"
                "QPushButton:hover { background: #3a72b3; }"
            )
            self._btn_combine.clicked.connect(self._on_combine_clicked)
            combine_lay.addWidget(self._btn_combine)
            outer.addWidget(combine_frame)

        # Standard Save / Discard buttons
        btn_box = QDialogButtonBox()
        self._btn_save = QPushButton("💾  Save")
        self._btn_save.setStyleSheet(
            "QPushButton { background: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px 18px; }"
            "QPushButton:hover { background: #218838; }"
        )
        self._btn_discard = QPushButton("✖  Discard")
        self._btn_discard.setStyleSheet(
            "QPushButton { background: #6c757d; color: white; "
            "border-radius: 4px; padding: 6px 18px; }"
            "QPushButton:hover { background: #5a6268; }"
        )
        btn_box.addButton(self._btn_save, QDialogButtonBox.AcceptRole)
        btn_box.addButton(self._btn_discard, QDialogButtonBox.RejectRole)
        btn_box.accepted.connect(self._on_save)
        btn_box.rejected.connect(self.reject)
        outer.addWidget(btn_box)

    def _make_row(self, old_name: str) -> QHBoxLayout:
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)

        if self._mode == "full":
            chk = QCheckBox()
            chk.setFixedWidth(30)
            self._select_checks[old_name] = chk
            row.addWidget(chk)

        lbl = QLabel(old_name)
        lbl.setStyleSheet("color: #bbb; font-size: 11px;")
        lbl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        row.addWidget(lbl, stretch=2)

        edit = QLineEdit(old_name)
        edit.setStyleSheet(
            "background: #1a1a2e; color: #ddd; border: 1px solid #555; "
            "border-radius: 2px; padding: 2px 4px; font-size: 11px;"
        )
        self._name_edits[old_name] = edit
        row.addWidget(edit, stretch=3)

        color_btn = QPushButton()
        color_btn.setFixedSize(36, 22)
        color_btn.setToolTip("Click to pick a color for this cluster")
        self._set_button_color(color_btn, "#808080")
        color_btn.clicked.connect(lambda _checked, n=old_name: self._on_color_btn_clicked(n))
        self._color_buttons[old_name] = color_btn
        row.addWidget(color_btn)

        return row

    # ── Color helpers ───────────────────────────────────────────────────

    def _set_button_color(self, btn: QPushButton, hex_color: str):
        btn.setStyleSheet(
            f"QPushButton {{ background: {hex_color}; border: 1px solid #555; "
            f"border-radius: 2px; }}"
            f"QPushButton:hover {{ border: 1px solid #aaa; }}"
        )
        btn.setProperty("hex_color", hex_color)

    def _get_button_color(self, btn: QPushButton) -> str:
        return btn.property("hex_color") or "#808080"

    def _on_color_btn_clicked(self, old_name: str):
        btn = self._color_buttons.get(old_name)
        if btn is None:
            return
        current = QColor(self._get_button_color(btn))
        chosen = QColorDialog.getColor(current, self, f"Color for '{old_name}'")
        if chosen.isValid():
            self._set_button_color(btn, chosen.name())

    # ── Population ──────────────────────────────────────────────────────

    def _populate(self):
        """Read cluster names from the adata and build the rename rows."""
        self._name_edits.clear()
        self._select_checks.clear()
        self._color_buttons.clear()

        # Clear existing rows (keep the trailing stretch)
        while self._rows_layout.count() > 1:
            item = self._rows_layout.takeAt(0)
            if item and item.widget():
                item.widget().deleteLater()
            elif item and item.layout():
                # clean up sub-layout
                while item.layout().count():
                    sub = item.layout().takeAt(0)
                    if sub and sub.widget():
                        sub.widget().deleteLater()

        clusters = self._get_cluster_names()
        if not clusters:
            self._status_lbl.setText("No clusters detected (run analysis first).")
            self._btn_save.setEnabled(False)
            return

        self._status_lbl.setText(
            f"{len(clusters)} cluster(s) detected."
        )
        self._btn_save.setEnabled(True)

        # Insert before the trailing stretch
        stretch_item = self._rows_layout.takeAt(self._rows_layout.count() - 1)
        for name in clusters:
            container = QWidget()
            row_lay = self._make_row(name)
            container.setLayout(row_lay)
            self._rows_layout.addWidget(container)
        self._rows_layout.addItem(stretch_item)

        # Apply saved colors to buttons
        self._apply_saved_colors_to_buttons(clusters)

    def _state_col(self) -> Optional[str]:
        if self._mode == "track":
            return _track_cluster_col(self._adata)
        return self._OBS_COLS.get(self._mode)

    def _apply_saved_colors_to_buttons(self, clusters: list):
        from behav3d.analysis.behavior.state.utils import (
            _get_classification_state_colors,
            _normalize_label_color_map,
        )
        if self._adata is None:
            return
        state_col = self._state_col()
        if state_col is None:
            return
        saved = _get_classification_state_colors(self._adata, state_col)
        colors = _normalize_label_color_map(clusters, colors=saved)
        for name, btn in self._color_buttons.items():
            if name in colors:
                self._set_button_color(btn, colors[name])

    def _get_cluster_names(self) -> list:
        if self._adata is None:
            return []
        col = self._state_col()
        if col is None or col not in self._adata.obs.columns:
            return []
        return _unique_sorted(
            self._adata.obs[col].dropna().astype(str).unique()
        )

    # ── Combine (full mode) ─────────────────────────────────────────────

    def _on_combine_clicked(self):
        target = self._combine_edit.text().strip()
        if not target:
            QMessageBox.warning(
                self, "Combine", "Enter a combined name first."
            )
            return
        selected = [
            name for name, chk in self._select_checks.items()
            if chk.isChecked()
        ]
        if not selected:
            QMessageBox.warning(
                self, "Combine", "Select at least one cluster to combine."
            )
            return
        # Use the color of the first selected cluster for all combined rows
        first_color = None
        for name in selected:
            if name in self._color_buttons:
                first_color = self._get_button_color(self._color_buttons[name])
                break
        for name in selected:
            if name in self._name_edits:
                self._name_edits[name].setText(target)
            if first_color is not None and name in self._color_buttons:
                self._set_button_color(self._color_buttons[name], first_color)
        self._combine_edit.clear()

    # ── Save ────────────────────────────────────────────────────────────

    def _on_save(self):
        mapping = {}
        for old_name, edit in self._name_edits.items():
            new_name = edit.text().strip()
            mapping[old_name] = new_name if new_name else old_name

        # Validate: no empty names
        empties = [k for k, v in mapping.items() if not v]
        if empties:
            QMessageBox.warning(
                self,
                "Empty Name",
                f"These clusters have no name:\n{', '.join(empties)}\n\n"
                "Fill in a name or leave the original.",
            )
            return

        try:
            self._apply_mapping(mapping)
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(
                self, "Rename Failed", f"Could not apply rename:\n\n{exc}"
            )
            return

        self.clusters_renamed.emit(mapping)
        self.accept()

    def _apply_mapping(self, mapping: dict):
        from behav3d.analysis.behavior.general import relabel_cluster_ids
        from behav3d.analysis.behavior.state.utils import (
            _coerce_hex_color,
            _set_classification_state_colors,
        )

        col = self._state_col()
        if col is None or col not in self._adata.obs.columns:
            return  # nothing to do

        # Save colors first (remap: color follows the renamed label)
        new_colors = {}
        for old_name, btn in self._color_buttons.items():
            new_label = mapping.get(old_name, old_name)
            new_colors[new_label] = _coerce_hex_color(self._get_button_color(btn))
        if new_colors:
            _set_classification_state_colors(self._adata, col, new_colors)

        # Only apply renaming if there are actual label changes
        existing = set(self._adata.obs[col].astype(str).unique())
        changed = any(
            str(mapping.get(lbl, lbl)) != str(lbl) for lbl in existing
        )
        if changed:
            relabel_cluster_ids(
                adata=self._adata,
                mapping=mapping,
                cluster_key=col,
                overwrite_original=True,
                keep_unmapped=True,
            )

        # Rebuild full_behavioral_cluster so its prefix reflects the renamed intrinsic labels
        if self._mode == "intrinsic" and "full_behavioral_cluster" in self._adata.obs.columns:
            from behav3d.analysis.behavior.state.utils import (
                _rebuild_full_behavioral_cluster_from_intrinsic,
            )
            import copy
            clustering_meta = self._adata.uns.get("clustering", {})
            binary_cols = clustering_meta.get("binary_cols_to_merge", [])
            bgc = copy.deepcopy(clustering_meta.get("binary_group_constraints", None))
            enforce = isinstance(bgc, dict) and "forbidden_binary_combinations" in bgc
            _rebuild_full_behavioral_cluster_from_intrinsic(
                self._adata,
                binary_cols_to_merge=binary_cols,
                binary_group_constraints=bgc,
                enforce_binary_group_constraints=enforce,
            )

        # Persist to disk
        if self._adata_path is not None:
            self._adata_path.parent.mkdir(parents=True, exist_ok=True)
            self._adata.write_h5ad(str(self._adata_path), compression="lzf")
