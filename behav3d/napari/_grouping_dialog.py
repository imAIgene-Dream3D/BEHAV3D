"""Cell-type grouping dialog for the Analysis tab.

Lets the user merge several already-*filtered* cell-type populations under
a single name for use in Death Dynamics / Single Cell. Group membership is
persisted to ``behav3d_parameters.yml`` (``cell_type_groups``) — metadata.csv
is never touched, so multicolor merged-tracking outputs (which must still
be processed through Feature Extraction + Filtering as their own real cell
types) are unaffected. The merge itself concatenates each member's filtered
track-features CSV into ``analysis/{name}_merged/track_features/...``.
"""

from __future__ import annotations

from pathlib import Path
import traceback

import yaml
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

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.analysis.grouping import (
    create_cell_type_group,
    create_group_tracked_segments,
    group_cell_type_id,
    group_tracked_segments_status,
    list_cell_type_groups,
    remove_cell_type_group,
    validate_group_name,
)
from behav3d.core.qt_help import HelpButton
from behav3d.napari._background_runner import BackgroundOperation, ProgressBarRow


class GroupBuilderDialog(QDialog):
    """Modal dialog letting the user merge filtered cell-type populations.

    Emits :pyattr:`group_created` with the merged cell-type id when the
    user successfully creates a new group, and :pyattr:`group_removed` when
    an existing group is deleted. The host (the Analysis tab) listens to
    these signals to refresh the Death Dynamics / Single Cell cell-type
    dropdowns.
    """

    group_created = Signal(str)
    group_removed = Signal(str)

    def __init__(self, metadata_loader, parent=None):
        super().__init__(parent)
        self.metadata_loader = metadata_loader
        self.setWindowTitle("Cell Type Grouping")
        self.setMinimumWidth(440)
        self._checkboxes: dict[str, QCheckBox] = {}
        self._category_for: dict[str, str] = {}
        self._all_types: list[str] = []
        self._bg = BackgroundOperation(self)
        self._build_ui()

    # ── UI construction ──────────────────────────────────────────────
    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)

        header = QLabel("<h3>Group Builder</h3>")
        layout.addWidget(header)

        instructions = QLabel(
            "Select cell types that have already been <b>filtered</b>, give "
            "the group a name, and click <b>Create Group</b>. The group "
            "merges each member's filtered track-features CSV into a single "
            "population, selectable in Death Dynamics and Single Cell "
            "alongside the original cell types. Metadata is not modified."
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
        self._all_types = all_types

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

        # Existing groups summary (from behav3d_parameters.yml).
        existing_groups = list_cell_type_groups(self._params())
        if existing_groups:
            ex_label = QLabel(
                "<i>Existing groups:</i> "
                + ", ".join(
                    f"{gid} ({', '.join(members)})"
                    for gid, members in existing_groups.items()
                )
            )
            ex_label.setWordWrap(True)
            ex_label.setStyleSheet("color:#888; font-size:10px;")
            layout.addWidget(ex_label)

        # Checkbox list (scrollable when many types)
        select_group_header = QHBoxLayout()
        select_group_header.setContentsMargins(0, 0, 0, 0)
        select_group_header.addWidget(HelpButton(
            "Selecting Cell Types vs. Groups",
            "The checkboxes below serve two purposes depending on what you "
            "select:\n\n"
            "- Check one or more real cell types (grayed out if not yet "
            "filtered) and click 'Create Group' to merge their filtered "
            "track-features into a new named population.\n\n"
            "- Check exactly one existing '(group)' entry to enable "
            "'Track Length Distribution' (preview) and 'Remove Group' for "
            "that group instead.\n\n"
            "A cell type only becomes selectable once it has been run "
            "through Filtering — 'not yet filtered' entries are disabled "
            "until then."
        ))
        select_group_header.addStretch()

        select_group = QGroupBox("Select cell types to group (must be filtered already)")
        sel_lay = QVBoxLayout(select_group)
        sel_lay.setSpacing(2)
        sel_lay.addLayout(select_group_header)

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

        output_dir = getattr(self.metadata_loader, "output_dir", None)
        for ct in all_types:
            cat = self._category_for.get(ct, "?")
            has_filtered = (
                output_dir is not None
                and self._filtered_csv_exists(output_dir, ct)
            )
            label = f"{ct}  ({cat})" if has_filtered else f"{ct}  ({cat}, not yet filtered)"
            cb = QCheckBox(label)
            cb.setEnabled(has_filtered)
            cb.stateChanged.connect(self._on_selection_changed)
            sel_content_lay.addWidget(cb)
            self._checkboxes[ct] = cb

        # Existing group ids double as removable "members" of one.
        for gid in existing_groups:
            label = f"{gid}  (group)"
            if output_dir is not None:
                status = group_tracked_segments_status(output_dir, gid, md)
                if status["total"] > 0:
                    label += f"  [tracked: {status['built']}/{status['total']}]"
            cb = QCheckBox(label)
            cb.stateChanged.connect(self._on_selection_changed)
            sel_content_lay.addWidget(cb)
            self._checkboxes[gid] = cb

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

        tracked_row = QHBoxLayout()
        self.build_tracked_checkbox = QCheckBox(
            "Also build tracked segments (slower)"
        )
        self.build_tracked_checkbox.setToolTip(
            "Rebuilds a per-sample tracked label image + tracks CSV for the "
            "group by relabeling each member's existing tracked images. "
            "Needed to see the group in Visualization/backprojection. "
            "This rewrites a label image per sample and can take a while."
        )
        tracked_row.addWidget(self.build_tracked_checkbox)
        tracked_row.addStretch()
        layout.addLayout(tracked_row)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Progress bar (shown while tracked segments build in the background)
        self.progress_row = ProgressBarRow(self)
        layout.addWidget(self.progress_row)

        # Buttons
        btn_row = QHBoxLayout()
        self.btn_preview_lengths = QPushButton("📊 Track Length Distribution")
        self.btn_preview_lengths.setToolTip(
            "Preview the track-length distribution of the selected existing "
            "group (select exactly one '(group)' entry above)."
        )
        self.btn_preview_lengths.setEnabled(False)
        self.btn_preview_lengths.clicked.connect(self._on_preview_lengths_clicked)
        btn_row.addWidget(self.btn_preview_lengths)
        self.btn_build_tracked = QPushButton("🧩 Build Tracked Segments")
        self.btn_build_tracked.setToolTip(
            "Build (or rebuild) per-sample tracked label images + tracks CSVs "
            "for the selected existing group, needed for it to appear in "
            "Visualization, Backprojection, and Exemplar export. Select "
            "exactly one '(group)' entry above. Rewrites a label image per "
            "sample and can take a while."
        )
        self.btn_build_tracked.setEnabled(False)
        self.btn_build_tracked.clicked.connect(self._on_build_tracked_clicked)
        btn_row.addWidget(self.btn_build_tracked)
        self.btn_remove = QPushButton("Remove Group")
        self.btn_remove.setToolTip(
            "Remove the selected group from the config.\n"
            "Only available when exactly one existing group is selected.\n"
            "The merged CSV on disk is kept."
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

    # ── Helpers ──────────────────────────────────────────────────────
    def _params(self) -> dict:
        return getattr(self.metadata_loader, "behav3d_parameters", {}) or {}

    @staticmethod
    def _filtered_csv_exists(output_dir, cell_type: str) -> bool:
        from behav3d.analysis.grouping import filtered_track_features_csv
        return filtered_track_features_csv(output_dir, cell_type).exists()

    def _selected_cell_types(self) -> list[str]:
        return [ct for ct, cb in self._checkboxes.items() if cb.isChecked()]

    def _existing_group_ids(self) -> set[str]:
        return set(list_cell_type_groups(self._params()).keys())

    def _persist_params(self):
        params = self._params()
        out_dir = getattr(self.metadata_loader, "output_dir", None)
        if not out_dir:
            return None
        params_path = Path(out_dir) / "behav3d_parameters.yml"
        try:
            with open(params_path, "w") as f:
                yaml.safe_dump(params, f, sort_keys=False)
            return str(params_path)
        except Exception as exc:
            self._show_error(f"Could not save behav3d_parameters.yml: {exc}")
            return None

    # ── Action: create ───────────────────────────────────────────────
    def _on_create_clicked(self):
        try:
            group_name = self.name_edit.text().strip()
            existing_groups = list_cell_type_groups(self._params())
            ok, err = validate_group_name(group_name, self._all_types, existing_groups)
            if not ok:
                self._show_error(err)
                return

            selected = self._selected_cell_types()
            if not selected:
                self._show_error("Select at least one cell type.")
                return

            output_dir = getattr(self.metadata_loader, "output_dir", None)
            if not output_dir:
                self._show_error("No output directory set.")
                return

            group_id = group_cell_type_id(group_name)
            try:
                out_csv = create_cell_type_group(output_dir, group_name, selected)
            except (FileNotFoundError, ValueError) as exc:
                self._show_error(str(exc))
                return

            params = self._params()
            groups = params.setdefault("cell_type_groups", {})
            groups[group_id] = list(selected)
            saved_to = self._persist_params()

            summary = (
                f"<b style='color:#81c784;'>✓ Group created:</b> "
                f"<code>{group_id}</code><br>"
                f"<i>Members:</i> {', '.join(selected)}<br>"
                f"<i>Merged CSV:</i> {out_csv}"
            )
            if saved_to:
                summary += f"<br><span style='color:#888;'>Saved to: {saved_to}</span>"
            self.status_label.setText(summary)

            try:
                self.group_created.emit(group_id)
            except Exception:
                pass

            if self.build_tracked_checkbox.isChecked():
                # Run in the background so the dialog stays responsive and
                # shows live progress; the sanity plot / info dialog / close
                # are deferred to _on_tracked_build_done.
                self._start_tracked_build(
                    output_dir, group_id, selected, out_csv, summary_base=summary
                )
                return

            # No tracked-segments build requested: finish immediately as before.
            self._show_length_distribution(output_dir, group_id)
            QMessageBox.information(
                self,
                "Group Created",
                f"Group '{group_id}' was created with {len(selected)} member(s).\n"
                f"Merged CSV: {out_csv}",
            )
            self.accept()

        except Exception as exc:
            traceback.print_exc()
            self._show_error(f"Unexpected error: {exc}")

    # ── Background: tracked-segment build ─────────────────────────────
    def _start_tracked_build(self, output_dir, group_id, selected=None, out_csv=None, summary_base=""):
        """Run ``create_group_tracked_segments`` off the UI thread with a
        live progress bar, reused by both group-creation (checkbox) and the
        existing-group "Build Tracked Segments" retrofit button."""
        if self._bg.is_running():
            self._show_error("Another background operation is already running.")
            return

        self._bg.run(
            fn=create_group_tracked_segments,
            kwargs={
                "output_dir": output_dir,
                "group_id": group_id,
                "metadata": self.metadata_loader.metadata,
            },
            desc=f"Building tracked segments for '{group_id}'…",
            progress_row=self.progress_row,
            buttons=[self.btn_create, self.btn_remove, self.btn_build_tracked, self.btn_close],
            on_done=lambda written: self._on_tracked_build_done(
                group_id, written, selected, out_csv, summary_base
            ),
            on_failed=lambda err: self._on_tracked_build_failed(group_id, err, summary_base),
        )

    def _on_tracked_build_done(self, group_id, written, selected, out_csv, summary_base):
        status = group_tracked_segments_status(
            getattr(self.metadata_loader, "output_dir", None), group_id, self.metadata_loader.metadata
        )
        summary = summary_base + (
            f"<br><b style='color:#81c784;'>✓ Tracked segments built for "
            f"{len(written)} sample(s)</b> ({status['built']}/{status['total']} total)."
        )
        self.status_label.setText(summary)

        if selected is not None and out_csv is not None:
            # This was the group-creation flow: finish with the same sanity
            # plot / info dialog / close that the synchronous path used to.
            output_dir = getattr(self.metadata_loader, "output_dir", None)
            self._show_length_distribution(output_dir, group_id)
            QMessageBox.information(
                self,
                "Group Created",
                f"Group '{group_id}' was created with {len(selected)} member(s).\n"
                f"Merged CSV: {out_csv}\n"
                f"Tracked segments: {len(written)} sample(s) written.",
            )
            self.accept()

    def _on_tracked_build_failed(self, group_id, error, summary_base):
        summary = summary_base + (
            f"<br><b style='color:#e57373;'>⚠ Could not build tracked segments "
            f"for '{group_id}':</b> {error}"
        )
        self.status_label.setText(summary)

    # ── Selection change ─────────────────────────────────────────────
    def _on_selection_changed(self):
        selected = self._selected_cell_types()
        existing = self._existing_group_ids()
        selected_groups = [s for s in selected if s in existing]
        can_remove = len(selected) == 1 and selected[0] in existing
        if hasattr(self, "btn_remove"):
            self.btn_remove.setEnabled(can_remove)
        if hasattr(self, "btn_preview_lengths"):
            self.btn_preview_lengths.setEnabled(len(selected_groups) == 1)
        if hasattr(self, "btn_build_tracked"):
            self.btn_build_tracked.setEnabled(len(selected_groups) == 1)

    # ── Action: preview track length distribution ───────────────────
    def _on_preview_lengths_clicked(self):
        selected = [s for s in self._selected_cell_types() if s in self._existing_group_ids()]
        if len(selected) != 1:
            self._show_error("Select exactly one existing group to preview.")
            return
        output_dir = getattr(self.metadata_loader, "output_dir", None)
        if not output_dir:
            self._show_error("No output directory set.")
            return
        self._show_length_distribution(output_dir, selected[0])

    def _show_length_distribution(self, output_dir, group_id: str):
        from behav3d.analysis.grouping import merged_track_features_csv
        from behav3d.widgets.utils import build_track_length_distribution_figure
        from behav3d.napari._pdf_view import show_matplotlib_figure

        csv_path = merged_track_features_csv(output_dir, group_id)
        fig = build_track_length_distribution_figure(
            csv_path, title=f"Track length distribution — {group_id}"
        )
        if fig is None:
            self._show_error(
                f"Could not build a track-length distribution for '{group_id}' "
                f"(missing or malformed CSV at {csv_path})."
            )
            return
        show_matplotlib_figure(fig, title=f"Track lengths — {group_id}", parent=self)

    # ── Action: build tracked segments (existing group) ──────────────
    def _on_build_tracked_clicked(self):
        """Build/rebuild tracked segments for the single selected existing group."""
        try:
            selected = [s for s in self._selected_cell_types() if s in self._existing_group_ids()]
            if len(selected) != 1:
                self._show_error("Select exactly one existing group to build tracked segments for.")
                return

            group_id = selected[0]
            output_dir = getattr(self.metadata_loader, "output_dir", None)
            if not output_dir:
                self._show_error("No output directory set.")
                return

            self.status_label.setText(
                f"<i>Building tracked segments for <code>{group_id}</code>… "
                "this may take a while.</i>"
            )
            self._start_tracked_build(output_dir, group_id)

        except Exception as exc:
            traceback.print_exc()
            self._show_error(f"Unexpected error: {exc}")

    # ── Action: remove ────────────────────────────────────────────────
    def _on_remove_clicked(self):
        """Remove the single selected group from ``cell_type_groups``."""
        try:
            selected = [s for s in self._selected_cell_types() if s in self._existing_group_ids()]
            if len(selected) != 1:
                self._show_error("Select exactly one existing group to remove.")
                return

            group_id = selected[0]
            confirm = QMessageBox.question(
                self,
                "Confirm Removal",
                f"Remove group '{group_id}' from the configuration?\n\n"
                "The merged CSV on disk will be kept, but the group will no "
                "longer be offered as a population in Death Dynamics / "
                "Single Cell.\n\nThis cannot be undone.",
                QMessageBox.Yes | QMessageBox.Cancel,
                QMessageBox.Cancel,
            )
            if confirm != QMessageBox.Yes:
                return

            params = self._params()
            remove_cell_type_group(params, group_id)
            saved_to = self._persist_params()

            summary = (
                f"<b style='color:#81c784;'>✓ Group removed:</b> <code>{group_id}</code>"
            )
            if saved_to:
                summary += f"<br><span style='color:#888;'>Saved to: {saved_to}</span>"
            self.status_label.setText(summary)

            try:
                self.group_removed.emit(group_id)
            except Exception:
                pass

            QMessageBox.information(
                self, "Group Removed", f"Group '{group_id}' was removed from the config."
            )
            self.accept()

        except Exception as exc:
            traceback.print_exc()
            self._show_error(f"Unexpected error: {exc}")

    def _show_error(self, message: str):
        self.status_label.setText(
            f"<b style='color:#e57373;'>Error:</b> {message}"
        )

    # ── Close guard ────────────────────────────────────────────────────
    def reject(self):
        """Block closing (Escape / Close button / window X) while the
        background tracked-segments build is running — the QThread outlives
        Python references to this dialog otherwise."""
        if self._bg.is_running():
            QMessageBox.warning(
                self,
                "Operation in progress",
                "Building tracked segments is still running.\n"
                "Please wait for it to finish before closing this dialog.",
            )
            return
        super().reject()
