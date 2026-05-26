"""
BEHAV3D napari plugin – shared Results panel widget.

Hosts the tree view of every file produced under
``output_dir/analysis/``. Instantiated independently inside each
top-level BEHAV3D tab that wants to surface results (Filtering,
Feature Extraction, Analysis); every instance points at the same
output directory through its ``metadata_loader`` so all trees stay in
sync after a ``refresh()``.

Public surface
--------------
- :class:`ResultsPanel` — the QWidget to embed.
- :func:`notify_results_changed` — walk the parent chain to the
  closest object that owns a ``results_panel`` attribute and call its
  ``refresh()``. Used by per-cell-type panels that don't know which
  outer tab they belong to.
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Optional

from qtpy.QtCore import QPointF, Qt, QTimer
from qtpy.QtGui import QBrush, QColor, QPainter, QPolygonF
from qtpy.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from behav3d.napari._pdf_view import (
    open_externally,
    open_image_in_napari,
    open_pdf_in_napari,
    reveal_in_folder,
)
from behav3d.napari._results_catalog import (
    CATEGORY_LABELS,
    SUBCATEGORY_LABELS,
    ResultFile,
    describe,
    group_by_tree,
    scan_outputs,
)


# ---------------------------------------------------------------------------
# Cross-widget refresh helper
# ---------------------------------------------------------------------------
def notify_results_changed(widget) -> bool:
    """Walk up ``widget``'s parent chain and refresh the first
    :class:`ResultsPanel` we encounter.

    Returns ``True`` if a panel was found and refreshed. Used by
    nested per-cell-type panels (FilteringTab/FeatureExtractionTab
    have per-CT sub-panels that don't directly own the panel).
    """
    node = widget
    while node is not None and not hasattr(node, "results_panel"):
        node = node.parent() if hasattr(node, "parent") else None
    panel = getattr(node, "results_panel", None) if node else None
    if panel is None:
        return False
    try:
        panel.refresh()
    except Exception:
        traceback.print_exc()
        return False
    return True


# ---------------------------------------------------------------------------
# Tree widget with visible expand/collapse arrows
# ---------------------------------------------------------------------------
class _ResultsTreeWidget(QTreeWidget):
    """QTreeWidget that draws bright, oversized branch indicators.

    Qt's default branch arrows are tiny and dark — on napari's dark
    theme they're almost invisible. We over-paint a filled triangle on
    top of the default indicator so the user can actually see whether
    a node is expanded or collapsed.
    """

    _ARROW_COLOR = QColor("#e6e6e6")
    _ARROW_SIZE = 9  # px

    def drawBranches(self, painter, rect, index):
        # Let Qt paint guide lines / indicator rect first (also keeps
        # hit-testing for clicks correct).
        super().drawBranches(painter, rect, index)

        if not index.isValid():
            return
        model = self.model()
        if model is None or not model.hasChildren(index):
            return

        size = self._ARROW_SIZE
        cx = rect.right() - self.indentation() / 2.0
        cy = rect.center().y() + 0.5

        painter.save()
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(self._ARROW_COLOR))

        if self.isExpanded(index):
            # Down-pointing triangle.
            triangle = QPolygonF([
                QPointF(cx - size / 2.0, cy - size / 3.0),
                QPointF(cx + size / 2.0, cy - size / 3.0),
                QPointF(cx, cy + size / 1.8),
            ])
        else:
            # Right-pointing triangle.
            triangle = QPolygonF([
                QPointF(cx - size / 2.5, cy - size / 2.0),
                QPointF(cx + size / 2.0, cy),
                QPointF(cx - size / 2.5, cy + size / 2.0),
            ])
        painter.drawPolygon(triangle)
        painter.restore()


# ---------------------------------------------------------------------------
# Widget
# ---------------------------------------------------------------------------
class ResultsPanel(QWidget):
    """Tree view of every file produced under ``output_dir/analysis/``.

    Grouped by analysis stage (Filtering, Feature Extraction,
    Analysis → Death Dynamics / Single Cell) and then by cell type.
    Each leaf carries a tooltip from
    :data:`behav3d.napari._results_catalog.FILE_CATALOG`.

    Double-click a PDF/PNG to open it in a dedicated secondary napari
    viewer (rasterized via ``pypdfium2``). Right-click anything for
    *Open externally* / *Reveal in folder* actions.
    """

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._collapsed = False
        self._init_ui()

        if metadata_loader is not None and hasattr(
            metadata_loader, "metadata_loaded"
        ):
            metadata_loader.metadata_loaded.connect(self.refresh)

        self.refresh()

        # Start in collapsed mode. Deferred so the panel is already
        # parented into its QSplitter when we ask the splitter to give
        # us the minimum slice of the screen.
        QTimer.singleShot(0, lambda: self._set_collapsed(True))

    # ── UI ─────────────────────────────────────────────────────────────
    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(4)

        # ── Header / toolbar ───────────────────────────────────────
        header = QHBoxLayout()
        header.setSpacing(6)

        title = QLabel("📄 Results")
        title.setStyleSheet("font-weight: bold; font-size: 12px; color: #ddd;")
        header.addWidget(title)
        header.addSpacing(8)

        header.addWidget(QLabel("DPI:"))
        self.spin_dpi = QSpinBox()
        self.spin_dpi.setRange(72, 300)
        self.spin_dpi.setValue(150)
        self.spin_dpi.setSingleStep(25)
        self.spin_dpi.setToolTip("Render resolution for PDF preview.")
        self.spin_dpi.setMaximumWidth(70)
        header.addWidget(self.spin_dpi)

        self.chk_show_extras = QCheckBox("Show non-viewable")
        self.chk_show_extras.setChecked(True)
        self.chk_show_extras.setToolTip(
            "Include CSV / HTML / JSON / Zarr entries. When unchecked, "
            "only PDFs and images are listed."
        )
        self.chk_show_extras.toggled.connect(self.refresh)
        header.addWidget(self.chk_show_extras)

        header.addStretch(1)

        self.btn_toggle = QToolButton()
        self.btn_toggle.setText("▾ Hide")
        self.btn_toggle.setToolTip("Collapse the Results panel.")
        self.btn_toggle.clicked.connect(self._on_toggle_collapsed)
        header.addWidget(self.btn_toggle)

        outer.addLayout(header)

        # ── Tree ───────────────────────────────────────────────────
        self.tree = _ResultsTreeWidget()
        self.tree.setHeaderLabels(["File"])
        self.tree.setColumnCount(1)
        # Alternating row colors clash with napari's dark theme (the
        # "alternate" tint ends up nearly white and unreadable), so we
        # leave them off and rely on indentation/expand markers instead.
        self.tree.setAlternatingRowColors(False)
        # Make the indentation a hair wider so our custom arrows have
        # room to breathe.
        self.tree.setIndentation(18)
        self.tree.setUniformRowHeights(False)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._on_context_menu)
        self.tree.itemDoubleClicked.connect(self._on_item_double_clicked)
        self.tree.itemSelectionChanged.connect(self._on_selection_changed)
        header_view = self.tree.header()
        try:
            header_view.setSectionResizeMode(0, QHeaderView.Stretch)
        except Exception:
            pass
        outer.addWidget(self.tree, stretch=1)

        # ── Selected-item details echo ─────────────────────────────
        self.details_label = QLabel(
            "Select a file to see what it is. Double-click a PDF or "
            "image to open it inside napari."
        )
        self.details_label.setWordWrap(True)
        self.details_label.setStyleSheet(
            "color: #aaa; font-size: 11px; padding: 2px 4px;"
        )
        outer.addWidget(self.details_label)

        # ── Per-item action buttons ────────────────────────────────
        actions = QHBoxLayout()
        actions.setSpacing(6)
        self.btn_open_napari = QPushButton("👁 Open in napari")
        self.btn_open_napari.setEnabled(False)
        self.btn_open_napari.clicked.connect(self._on_open_napari)
        actions.addWidget(self.btn_open_napari)
        self.btn_open_ext = QPushButton("Open externally")
        self.btn_open_ext.setEnabled(False)
        self.btn_open_ext.clicked.connect(self._on_open_external)
        actions.addWidget(self.btn_open_ext)
        self.btn_reveal = QPushButton("Reveal in folder")
        self.btn_reveal.setEnabled(False)
        self.btn_reveal.clicked.connect(self._on_reveal)
        actions.addWidget(self.btn_reveal)
        actions.addStretch(1)
        outer.addLayout(actions)

        self._body_widgets = (
            self.tree, self.details_label,
            self.btn_open_napari, self.btn_open_ext, self.btn_reveal,
        )

    # ── Hide / show toggle ────────────────────────────────────────────
    def _on_toggle_collapsed(self):
        self._set_collapsed(not self._collapsed)

    def _set_collapsed(self, collapsed: bool):
        """Show / hide the body and ask the parent splitter to resize.

        When collapsed the panel is constrained to its header height so
        the parent ``QSplitter`` gives the bulk of the screen back to
        the tab content; expanding lifts the constraint and restores a
        sensible default split.
        """
        self._collapsed = collapsed
        for w in self._body_widgets:
            w.setVisible(not collapsed)
        self.btn_toggle.setText("▴ Show" if collapsed else "▾ Hide")

        if collapsed:
            # Force the layout to recompute now that the body widgets
            # are hidden, then pin the height to whatever the header
            # alone wants.
            self.adjustSize()
            target_h = max(self.sizeHint().height(),
                           self.minimumSizeHint().height())
            self.setMaximumHeight(target_h)
            self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        else:
            self.setMaximumHeight(16777215)  # Qt's QWIDGETSIZE_MAX
            self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)

        self._resize_parent_splitter(collapsed)

    def _resize_parent_splitter(self, collapsed: bool):
        """Reallocate space in the surrounding ``QSplitter`` (if any).

        When collapsed we give the panel just enough for its header
        and hand the rest to the sibling widgets. When expanded we
        restore a 60/40-ish split between tab content and results.
        """
        splitter = self._find_parent_splitter()
        if splitter is None:
            return
        idx = splitter.indexOf(self)
        if idx < 0:
            return

        sizes = list(splitter.sizes())
        total = sum(sizes)
        if total <= 0:
            return

        if collapsed:
            header_h = max(self.sizeHint().height(),
                           self.minimumSizeHint().height(),
                           1)
            sizes[idx] = header_h
            remaining = max(total - header_h, 1)
            others = [i for i in range(len(sizes)) if i != idx]
            other_total = sum(splitter.sizes()[i] for i in others) or 1
            for i in others:
                sizes[i] = max(
                    int(remaining * splitter.sizes()[i] / other_total),
                    1,
                )
        else:
            # Restore a comfortable expanded split — sibling tabs keep
            # the majority of the space, the panel takes ~40%.
            target = int(total * 0.4)
            sizes[idx] = max(target, 200)
            remaining = max(total - sizes[idx], 1)
            others = [i for i in range(len(sizes)) if i != idx]
            other_total = sum(splitter.sizes()[i] for i in others) or 1
            for i in others:
                sizes[i] = max(
                    int(remaining * splitter.sizes()[i] / other_total),
                    1,
                )
        splitter.setSizes(sizes)

    def _find_parent_splitter(self) -> Optional[QSplitter]:
        node = self.parent()
        while node is not None:
            if isinstance(node, QSplitter):
                return node
            node = node.parent() if hasattr(node, "parent") else None
        return None

    # ── Public API ─────────────────────────────────────────────────────
    def refresh(self):
        """Re-scan the output directory and rebuild the tree."""
        self.tree.clear()
        out_dir = self._output_dir()
        if out_dir is None:
            placeholder = QTreeWidgetItem(
                self.tree,
                ["⚠ No output directory configured"],
            )
            placeholder.setToolTip(
                0, "Load metadata in the Data Preparation tab."
            )
            placeholder.setFlags(placeholder.flags() & ~Qt.ItemIsSelectable)
            return

        files = scan_outputs(out_dir)
        if not self.chk_show_extras.isChecked():
            files = [f for f in files if f.is_viewable]
        if not files:
            empty = QTreeWidgetItem(
                self.tree,
                ["(No outputs found yet)"],
            )
            empty.setToolTip(0, f"Searched: {out_dir / 'analysis'}")
            empty.setFlags(empty.flags() & ~Qt.ItemIsSelectable)
            return

        # Top-level categories in deterministic order. ``Analysis`` is the
        # only one with subcategories.
        category_order = ["filtering", "feature_extraction", "analysis"]
        subcategory_order = ["death_dynamics", "single_cell", None]

        tree = group_by_tree(files)
        for cat in category_order:
            cat_files = tree.get(cat)
            if not cat_files:
                continue
            cat_item = QTreeWidgetItem(
                self.tree, [CATEGORY_LABELS.get(cat, cat)]
            )
            cat_item.setToolTip(0, self._category_tooltip(cat))
            cat_item.setExpanded(True)

            if cat == "analysis":
                for sub in subcategory_order:
                    sub_files = cat_files.get(sub)
                    if not sub_files:
                        continue
                    sub_label = (
                        SUBCATEGORY_LABELS.get(sub, sub or "Other")
                    )
                    sub_item = QTreeWidgetItem(cat_item, [sub_label])
                    sub_item.setToolTip(0, self._subcategory_tooltip(sub))
                    sub_item.setExpanded(True)
                    self._populate_cell_types(sub_item, sub_files, out_dir)
            else:
                # Filtering / Feature Extraction — no subcategory.
                merged: dict = {}
                for sub_files in cat_files.values():
                    for ct, lst in sub_files.items():
                        merged.setdefault(ct, []).extend(lst)
                self._populate_cell_types(cat_item, merged, out_dir)

        # Restore the action-button state to "no selection".
        self._set_action_state(None)

    def set_metadata_loader(self, loader):
        """Swap in a new metadata loader (used by ``_widget``)."""
        if self.metadata_loader is loader:
            return
        self.metadata_loader = loader
        if loader is not None and hasattr(loader, "metadata_loaded"):
            try:
                loader.metadata_loaded.connect(self.refresh)
            except Exception:
                pass
        self.refresh()

    # ── Helpers ────────────────────────────────────────────────────────
    def _output_dir(self) -> Optional[Path]:
        if self.metadata_loader is None:
            return None
        out_dir = getattr(self.metadata_loader, "output_dir", None)
        if not out_dir:
            return None
        try:
            return Path(out_dir).expanduser()
        except Exception:
            return None

    def _category_tooltip(self, cat: str) -> str:
        if cat == "filtering":
            return ("Per-target QC plots produced during track filtering "
                    "(touching distributions, dead-dye distributions, "
                    "filter counts).")
        if cat == "feature_extraction":
            return ("Combined per-track / per-timepoint feature tables "
                    "produced by feature extraction.")
        if cat == "analysis":
            return ("Downstream analyses: death dynamics, interaction "
                    "analysis, single-cell clustering.")
        return ""

    def _subcategory_tooltip(self, sub: Optional[str]) -> str:
        if sub == "death_dynamics":
            return ("Death-dynamics outputs: per-organoid analyses, "
                    "interaction analyses, multi-organoid comparisons, "
                    "morpho-dead clustering.")
        if sub == "single_cell":
            return ("Single-cell outputs: UMAP clustering, behavioral "
                    "states, state transitions, exemplar tracks.")
        return ""

    def _populate_cell_types(
        self,
        parent_item: QTreeWidgetItem,
        cell_buckets: dict,
        out_dir: Path,
    ):
        analysis_root = out_dir / "analysis"
        for cell_type in sorted(
            cell_buckets.keys(),
            key=lambda c: ("" if c is None else c.lower()),
        ):
            files: list[ResultFile] = cell_buckets[cell_type]
            if not files:
                continue
            ct_label = cell_type if cell_type else "(unknown)"
            ct_item = QTreeWidgetItem(parent_item, [ct_label])
            ct_item.setToolTip(0, f"Cell type / folder: {ct_label}")
            ct_item.setExpanded(True)

            # Group files within this cell type by relative directory so
            # the tree mirrors the on-disk layout (results/, per_sample/,
            # behavioral_states/, etc.).
            by_dir: dict[tuple[str, ...], list[ResultFile]] = {}
            for f in files:
                try:
                    rel = f.path.relative_to(
                        analysis_root / (cell_type or "")
                    )
                except ValueError:
                    rel = Path(f.path.name)
                parts = rel.parts[:-1]  # subdirectories only
                by_dir.setdefault(parts, []).append(f)

            # Sort directories: shallow first, then alphabetical.
            dir_keys = sorted(by_dir.keys(), key=lambda p: (len(p), p))
            for parts in dir_keys:
                if parts:
                    folder_item = ct_item
                    for segment in parts:
                        folder_item = self._ensure_folder_child(
                            folder_item, segment
                        )
                        # Collapse track_intermediates / per_sample by
                        # default to keep the tree readable.
                        if segment in {"track_intermediates", "per_sample",
                                       "intermediate_results"}:
                            folder_item.setExpanded(False)
                else:
                    folder_item = ct_item

                for f in sorted(
                    by_dir[parts], key=lambda r: r.label.lower()
                ):
                    leaf = QTreeWidgetItem(folder_item, [f.label])
                    leaf.setData(0, Qt.UserRole, f)
                    leaf.setToolTip(
                        0, f.description or f"{f.label} ({f.path})"
                    )
                    if not f.is_viewable:
                        try:
                            from qtpy.QtGui import QBrush, QColor
                            leaf.setForeground(0, QBrush(QColor("#bbb")))
                        except Exception:
                            pass

    @staticmethod
    def _ensure_folder_child(
        parent_item: QTreeWidgetItem, name: str
    ) -> QTreeWidgetItem:
        """Return an existing folder child with ``name`` or create one."""
        for i in range(parent_item.childCount()):
            child = parent_item.child(i)
            stored = child.data(0, Qt.UserRole)
            if stored is None and child.text(0) == name:
                return child
        node = QTreeWidgetItem(parent_item, [name])
        node.setToolTip(0, f"Folder: {name}")
        node.setExpanded(True)
        return node

    # ── Selection / action plumbing ────────────────────────────────────
    def _current_file(self) -> Optional[ResultFile]:
        item = self.tree.currentItem()
        if item is None:
            return None
        data = item.data(0, Qt.UserRole)
        return data if isinstance(data, ResultFile) else None

    def _on_selection_changed(self):
        f = self._current_file()
        self._set_action_state(f)
        if f is None:
            self.details_label.setText(
                "Select a file to see what it is. Double-click a PDF or "
                "image to open it inside napari."
            )
        else:
            desc = f.description or describe(f.path) or "(no description)"
            self.details_label.setText(f"<b>{f.label}</b><br>{desc}")

    def _set_action_state(self, f: Optional[ResultFile]):
        has_file = f is not None
        self.btn_open_napari.setEnabled(bool(has_file and f.is_viewable))
        self.btn_open_ext.setEnabled(has_file)
        self.btn_reveal.setEnabled(has_file)

    def _on_item_double_clicked(self, item: QTreeWidgetItem, _column: int):
        f = item.data(0, Qt.UserRole)
        if not isinstance(f, ResultFile):
            item.setExpanded(not item.isExpanded())
            return
        if f.is_viewable:
            self._open_in_napari(f)
        else:
            open_externally(f.path)

    def _on_open_napari(self):
        f = self._current_file()
        if f is not None and f.is_viewable:
            self._open_in_napari(f)

    def _on_open_external(self):
        f = self._current_file()
        if f is not None:
            open_externally(f.path)

    def _on_reveal(self):
        f = self._current_file()
        if f is not None:
            reveal_in_folder(f.path)

    def _on_context_menu(self, pos):
        item = self.tree.itemAt(pos)
        if item is None:
            return
        f = item.data(0, Qt.UserRole)
        menu = QMenu(self)
        if isinstance(f, ResultFile):
            if f.is_viewable:
                act_napari = menu.addAction("Open in napari")
                act_napari.triggered.connect(lambda: self._open_in_napari(f))
            act_ext = menu.addAction("Open externally")
            act_ext.triggered.connect(lambda: open_externally(f.path))
            act_rev = menu.addAction("Reveal in folder")
            act_rev.triggered.connect(lambda: reveal_in_folder(f.path))
        else:
            # Folder / category node — only "Expand all / Collapse all".
            act_expand = menu.addAction(
                "Collapse all" if item.isExpanded() else "Expand all"
            )
            act_expand.triggered.connect(
                lambda: self._toggle_subtree(item, not item.isExpanded())
            )
        menu.exec_(self.tree.viewport().mapToGlobal(pos))

    def _toggle_subtree(self, item: QTreeWidgetItem, expand: bool):
        item.setExpanded(expand)
        for i in range(item.childCount()):
            self._toggle_subtree(item.child(i), expand)

    def _open_in_napari(self, f: ResultFile):
        try:
            if f.kind == "pdf":
                open_pdf_in_napari(
                    f.path, dpi=int(self.spin_dpi.value())
                )
            elif f.kind == "image":
                open_image_in_napari(f.path)
        except Exception as e:
            traceback.print_exc()
            QMessageBox.warning(
                self,
                "Could not open file",
                f"Failed to open {f.label}:\n{e}",
            )


__all__ = ["ResultsPanel", "notify_results_changed"]
