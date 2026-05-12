"""Qt port of :class:`behav3d.widgets.grouping.GroupBuilder`.

Lets the user group several detected cell types under a single name and
persist the resulting columns to the metadata CSV. The dialog is opened
from :class:`FeatureExtractionTab` via a button placed above the
per-cell-type tab widget.
"""

from __future__ import annotations

from pathlib import Path
import traceback

from qtpy.QtCore import Qt, Signal
from qtpy.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

import pandas as pd

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
)


class GroupBuilderDialog(QDialog):
    """Modal dialog letting the user create a ``{name}_merged`` cell type.

    Emits :pyattr:`group_created` with the merged cell-type name when the
    user successfully creates a new group, and :pyattr:`group_removed` when
    an existing ``*_merged`` group is deleted. The host (typically the
    Feature Extraction tab) listens to these signals to rebuild its tabs.
    """

    group_created = Signal(str)
    group_removed = Signal(str)

    def __init__(self, metadata_loader, parent=None):
        super().__init__(parent)
        self.metadata_loader = metadata_loader
        self.setWindowTitle("Cell Type Grouping")
        self.setMinimumWidth(420)
        self._checkboxes: dict[str, QCheckBox] = {}
        self._category_for: dict[str, str] = {}
        self._build_ui()

    # ── UI construction ──────────────────────────────────────────────
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        header = QLabel("<h3>Group Builder</h3>")
        layout.addWidget(header)

        instructions = QLabel(
            "Select existing cell types to group, give the group a name, and "
            "click <b>Create Group</b> to add it to the metadata. Groups are "
            "stored as <code>{name}_merged</code> columns and become available "
            "alongside the original cell types in downstream steps."
        )
        instructions.setWordWrap(True)
        instructions.setStyleSheet("color:#bbb; font-size:11px;")
        layout.addWidget(instructions)

        md = getattr(self.metadata_loader, "metadata", None)
        if md is None or md.empty:
            err = QLabel(
                "<b style='color:#e57373;'>No metadata loaded.</b><br>"
                "Open the Data Preparation tab and load metadata first."
            )
            err.setWordWrap(True)
            layout.addWidget(err)
            buttons = QDialogButtonBox(QDialogButtonBox.Close)
            buttons.rejected.connect(self.reject)
            layout.addWidget(buttons)
            return

        organoid_types = detect_organoid_types_from_metadata(md)
        immune_types = detect_immune_cell_types_from_metadata(md)
        other_types = detect_other_cell_types_from_metadata(md)
        all_types = sorted(organoid_types + immune_types + other_types)

        if not all_types:
            warn = QLabel(
                "<b style='color:#ffb74d;'>No cell types found in metadata.</b>"
            )
            warn.setWordWrap(True)
            layout.addWidget(warn)
            buttons = QDialogButtonBox(QDialogButtonBox.Close)
            buttons.rejected.connect(self.reject)
            layout.addWidget(buttons)
            return

        # Existing merged groups summary
        existing = detect_merged_cell_types_from_metadata(md)
        if existing:
            ex_label = QLabel(
                f"<i>Existing groups:</i> {', '.join(existing)}"
            )
            ex_label.setWordWrap(True)
            ex_label.setStyleSheet("color:#888; font-size:10px;")
            layout.addWidget(ex_label)

        # Checkbox list (scrollable when many types)
        select_group = QGroupBox("Select cell types to group")
        sel_lay = QVBoxLayout(select_group)
        sel_lay.setSpacing(2)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        sel_content = QWidget()
        sel_content_lay = QVBoxLayout(sel_content)
        sel_content_lay.setContentsMargins(2, 2, 2, 2)
        sel_content_lay.setSpacing(2)

        self._category_for = {ct: "organoid" for ct in organoid_types}
        self._category_for.update({ct: "immune" for ct in immune_types})
        self._category_for.update({ct: "other" for ct in other_types})

        for ct in all_types:
            cat = self._category_for.get(ct, "?")
            cb = QCheckBox(f"{ct}  ({cat})")
            cb.stateChanged.connect(self._on_selection_changed)
            sel_content_lay.addWidget(cb)
            self._checkboxes[ct] = cb
        sel_content_lay.addStretch()
        scroll.setWidget(sel_content)
        scroll.setMinimumHeight(140)
        sel_lay.addWidget(scroll)
        layout.addWidget(select_group)

        # Group name
        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Group name*:"))
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("e.g. tcells, tumors")
        name_row.addWidget(self.name_edit, stretch=1)
        layout.addLayout(name_row)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Buttons
        btn_row = QHBoxLayout()
        self.btn_remove = QPushButton("Remove Group")
        self.btn_remove.setToolTip(
            "Remove the selected *_merged group and its metadata columns.\n"
            "Only available when exactly one _merged cell type is selected."
        )
        self.btn_remove.setStyleSheet(
            "background-color:#c62828; color:white; font-weight:bold; "
            "padding:6px 14px; border-radius:4px;"
        )
        self.btn_remove.setEnabled(False)
        self.btn_remove.clicked.connect(self._on_remove_clicked)
        btn_row.addWidget(self.btn_remove)
        btn_row.addStretch()
        self.btn_create = QPushButton("Create Group")
        self.btn_create.setStyleSheet(
            "background-color:#28a745; color:white; font-weight:bold; "
            "padding:6px 14px; border-radius:4px;"
        )
        self.btn_create.clicked.connect(self._on_create_clicked)
        btn_row.addWidget(self.btn_create)
        self.btn_close = QPushButton("Close")
        self.btn_close.clicked.connect(self.reject)
        btn_row.addWidget(self.btn_close)
        layout.addLayout(btn_row)

    # ── Validation ───────────────────────────────────────────────────
    def _validate_group_name(self, group_name: str) -> tuple[bool, str]:
        group_name = (group_name or "").strip()
        if not group_name:
            return False, "Group name cannot be empty."
        if not all(c.isalnum() or c in "_-" for c in group_name):
            return False, (
                "Group name can only contain alphanumeric characters, hyphens, "
                "and underscores."
            )
        md = self.metadata_loader.metadata
        existing = detect_merged_cell_types_from_metadata(md)
        if group_name in existing:
            return False, f"Group '{group_name}' already exists in metadata."
        all_types = (
            detect_organoid_types_from_metadata(md)
            + detect_immune_cell_types_from_metadata(md)
            + detect_other_cell_types_from_metadata(md)
        )
        if group_name in all_types:
            return False, (
                f"Group name '{group_name}' conflicts with an existing cell type."
            )
        return True, ""

    def _selected_cell_types(self) -> list[str]:
        return [ct for ct, cb in self._checkboxes.items() if cb.isChecked()]

    # ── Action ───────────────────────────────────────────────────────
    def _on_create_clicked(self):
        try:
            group_name = self.name_edit.text().strip()
            ok, err = self._validate_group_name(group_name)
            if not ok:
                self._show_error(err)
                return

            selected = self._selected_cell_types()
            if not selected:
                self._show_error("Select at least one cell type.")
                return

            md = self.metadata_loader.metadata
            if md is None or md.empty:
                self._show_error("No metadata loaded.")
                return

            organoid_types = detect_organoid_types_from_metadata(md)
            immune_types = detect_immune_cell_types_from_metadata(md)
            other_types = detect_other_cell_types_from_metadata(md)

            sel_organoid = [t for t in selected if t in organoid_types]
            sel_immune = [t for t in selected if t in immune_types]

            if sel_organoid:
                prefix = "or_"
            elif sel_immune:
                prefix = "im_"
            else:
                prefix = "ot_"

            merged_name = f"{group_name}_merged"
            col_prefix = f"{prefix}{merged_name}"
            suffixes = [
                "line_condition",
                "segments_image_path",
                "tracks_image_path",
                "tracks_csv_path",
            ]

            added_cols: list[str] = []
            for suffix in suffixes:
                col_name = f"{col_prefix}_{suffix}"
                if col_name not in md.columns:
                    md[col_name] = pd.NA
                    added_cols.append(col_name)

            # Persist when the loader knows the CSV path.
            csv_path = (
                getattr(self.metadata_loader, "metadata_csv_path", None)
                or getattr(self.metadata_loader, "metadata_csv", None)
            )
            saved_to: str | None = None
            if not csv_path:
                # Fall back to the loader's behav3d_parameters paths entry.
                params = getattr(self.metadata_loader, "behav3d_parameters", {}) or {}
                csv_path = params.get("paths", {}).get("metadata_csv")
            if csv_path:
                try:
                    md.to_csv(csv_path, index=False)
                    saved_to = str(csv_path)
                except Exception as exc:
                    self._show_error(f"Could not save metadata CSV: {exc}")
                    return

            summary = (
                f"<b style='color:#81c784;'>✓ Group created:</b> "
                f"<code>{merged_name}</code><br>"
                f"<i>Members:</i> {', '.join(selected)}<br>"
                f"<i>Columns added:</i> {len(added_cols)}"
            )
            if saved_to:
                summary += f"<br><span style='color:#888;'>Saved to: {saved_to}</span>"
            self.status_label.setText(summary)

            try:
                self.group_created.emit(merged_name)
            except Exception:
                pass

            # Auto-close after a short delay so the user sees the status.
            # We use accept() which keeps the dialog modal-friendly.
            QMessageBox.information(
                self,
                "Group Created",
                f"Group '{merged_name}' was created with {len(selected)} member(s)."
                + (f"\nMetadata CSV updated: {saved_to}" if saved_to else ""),
            )
            self.accept()

        except Exception as exc:
            traceback.print_exc()
            self._show_error(f"Unexpected error: {exc}")

    # ── Selection change ─────────────────────────────────────────────
    def _on_selection_changed(self):
        """Enable the Remove Group button only when exactly one *_merged type is checked."""
        if not hasattr(self, "btn_remove"):
            return
        selected = self._selected_cell_types()
        can_remove = (
            len(selected) == 1
            and selected[0].endswith("_merged")
        )
        self.btn_remove.setEnabled(can_remove)

    # ── Remove action ────────────────────────────────────────────────
    def _on_remove_clicked(self):
        """Remove the single selected *_merged group from metadata."""
        try:
            selected = self._selected_cell_types()
            if len(selected) != 1 or not selected[0].endswith("_merged"):
                self._show_error(
                    "Select exactly one <code>*_merged</code> cell type to remove."
                )
                return

            merged_name = selected[0]
            md = self.metadata_loader.metadata
            if md is None or md.empty:
                self._show_error("No metadata loaded.")
                return

            category = self._category_for.get(merged_name, "")
            prefix_map = {"organoid": "or_", "immune": "im_", "other": "ot_"}
            prefix = prefix_map.get(category, "")

            # Collect all columns that belong to this merged group.
            # Cover both prefixed (or_/im_/ot_) and non-prefixed forms.
            cols_to_drop = [
                col for col in md.columns
                if col == f"{prefix}{merged_name}"
                or col.startswith(f"{prefix}{merged_name}_")
                or col == merged_name
                or col.startswith(f"{merged_name}_")
            ]

            if not cols_to_drop:
                self._show_error(
                    f"No metadata columns found for group <code>{merged_name}</code>."
                )
                return

            confirm = QMessageBox.question(
                self,
                "Confirm Removal",
                f"Remove group '{merged_name}' and its {len(cols_to_drop)} "
                f"metadata column(s)?\n\n"
                + "\n".join(f"  • {c}" for c in cols_to_drop)
                + "\n\nThis cannot be undone.",
                QMessageBox.Yes | QMessageBox.Cancel,
                QMessageBox.Cancel,
            )
            if confirm != QMessageBox.Yes:
                return

            md.drop(columns=cols_to_drop, inplace=True)

            csv_path = (
                getattr(self.metadata_loader, "metadata_csv_path", None)
                or getattr(self.metadata_loader, "metadata_csv", None)
            )
            if not csv_path:
                params = getattr(self.metadata_loader, "behav3d_parameters", {}) or {}
                csv_path = params.get("paths", {}).get("metadata_csv")

            saved_to: str | None = None
            if csv_path:
                try:
                    md.to_csv(csv_path, index=False)
                    saved_to = str(csv_path)
                except Exception as exc:
                    self._show_error(f"Could not save metadata CSV: {exc}")
                    return

            summary = (
                f"<b style='color:#81c784;'>✓ Group removed:</b> "
                f"<code>{merged_name}</code><br>"
                f"<i>Columns removed:</i> {len(cols_to_drop)}"
            )
            if saved_to:
                summary += f"<br><span style='color:#888;'>Saved to: {saved_to}</span>"
            self.status_label.setText(summary)

            try:
                self.group_removed.emit(merged_name)
            except Exception:
                pass

            QMessageBox.information(
                self,
                "Group Removed",
                f"Group '{merged_name}' was removed ({len(cols_to_drop)} column(s) deleted)."
                + (f"\nMetadata CSV updated: {saved_to}" if saved_to else ""),
            )
            self.accept()

        except Exception as exc:
            traceback.print_exc()
            self._show_error(f"Unexpected error: {exc}")

    def _show_error(self, message: str):
        self.status_label.setText(
            f"<b style='color:#e57373;'>Error:</b> {message}"
        )
