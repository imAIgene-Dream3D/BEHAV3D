"""
BEHAV3D napari plugin – main dock widget.
Provides a QTabWidget with tabs for the full BEHAV3D pipeline.
"""
from qtpy.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QTabWidget, QPushButton
from qtpy.QtCore import Qt, QSize, QEvent
from qtpy.QtGui import QIcon
import os
import napari

from behav3d.napari._queue import ProcessingQueuePanel, StepType
from behav3d.napari._global_workers import GlobalWorkersController
from behav3d.core.qt_help import disable_spinbox_wheel_scroll
from pathlib import Path


class FloatingAssistantButton(QPushButton):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setCheckable(True)
        self.setChecked(True)
        self.setToolTip("Show or hide the BEHAV3D Assistant panel")
        self.setFixedSize(60, 60)
        self.setCursor(Qt.PointingHandCursor)
        self.setStyleSheet("""
            QPushButton {
                background-color: rgba(42, 45, 48, 200);
                border: 2px solid #555;
                border-radius: 30px;
            }
            QPushButton:hover {
                background-color: rgba(80, 85, 90, 220);
            }
            QPushButton:pressed {
                background-color: rgba(100, 105, 110, 240);
                border: 2px solid #888;
            }
            QPushButton:checked {
                background-color: rgba(60, 120, 180, 200);
                border: 2px solid #6ba;
            }
            QPushButton:checked:pressed {
                background-color: rgba(70, 130, 190, 220);
                border: 2px solid #7cb;
            }
        """)
        icon_path = str(Path(__file__).resolve().parents[1] / "resources" / "assistant_icon.png")
        if os.path.exists(icon_path):
            self.setIcon(QIcon(icon_path))
            self.setIconSize(QSize(50, 50))
        else:
            self.setText("💬")
        
        self.setCheckable(True)
        self.setChecked(True)
        
        self._drag_active = False
        self._drag_start_global_pos = None
        self._drag_start_window_pos = None
        self._anchor_corner = "br"
        self._anchor_offset = (20, 20)

        if parent is not None:
            parent.installEventFilter(self)
            # Use a timer to initialize position once layout is computed
            from qtpy.QtCore import QTimer
            QTimer.singleShot(100, self._update_position)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._drag_start_global_pos = event.globalPos()
            self._drag_start_window_pos = self.pos()
            self._drag_active = False
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.LeftButton and self._drag_start_global_pos is not None:
            # Add a small distance threshold to differentiate a click from a drag
            diff = event.globalPos() - self._drag_start_global_pos
            if not self._drag_active and diff.manhattanLength() < 5:
                super().mouseMoveEvent(event)
                return

            if self.parent() is not None:
                new_pos = self._drag_start_window_pos + diff
                
                # Constrain within parent widget bounds
                p_rect = self.parent().rect()
                x = max(0, min(new_pos.x(), p_rect.width() - self.width()))
                y = max(0, min(new_pos.y(), p_rect.height() - self.height()))
                self.move(x, y)
                self._drag_active = True

                # Since the button follows the cursor, Qt keeps thinking the
                # mouse is hovering a pressed button and re-applies the
                # ":pressed" style on every move. Force it back to the
                # "up" visual state so only :checked (active/hidden) drives
                # the color while dragging, and skip the base-class handler
                # so it doesn't immediately re-set "down" itself.
                self.setDown(False)
                return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self._drag_active:
            # We dragged the button. Call super to reset visual state (un-press),
            # but block signals so it doesn't trigger a 'click' action.
            self._update_anchor_from_position()
            self.blockSignals(True)
            super().mouseReleaseEvent(event)
            self.blockSignals(False)
        else:
            # Normal click
            super().mouseReleaseEvent(event)

        self._drag_start_global_pos = None
        self._drag_start_window_pos = None
        self._drag_active = False

    def eventFilter(self, obj, event):
        if obj == self.parent() and event.type() == QEvent.Resize:
            self._update_position()
        return super().eventFilter(obj, event)

    def _update_anchor_from_position(self):
        """Remember which corner the button is closest to, and its offset
        from that corner, so _update_position can keep it there as the
        parent resizes (e.g. when the assistant panel opens/closes)."""
        if self.parent() is None:
            return
        p_rect = self.parent().rect()
        x, y = self.x(), self.y()
        left = x
        right = p_rect.width() - self.width() - x
        top = y
        bottom = p_rect.height() - self.height() - y

        horizontal = "l" if left <= right else "r"
        vertical = "t" if top <= bottom else "b"
        self._anchor_corner = vertical + horizontal
        self._anchor_offset = (
            left if horizontal == "l" else right,
            top if vertical == "t" else bottom,
        )

    def _update_position(self):
        if self.parent() is None:
            return
        p_rect = self.parent().rect()
        dx, dy = self._anchor_offset
        vertical, horizontal = self._anchor_corner[0], self._anchor_corner[1]

        x = dx if horizontal == "l" else p_rect.width() - self.width() - dx
        y = dy if vertical == "t" else p_rect.height() - self.height() - dy

        x = max(0, min(x, p_rect.width() - self.width()))
        y = max(0, min(y, p_rect.height() - self.height()))

        self.move(x, y)
        self.raise_()


from behav3d.napari._analysis import _period_from_t_radios


class BEHAV3DWidget(QWidget):
    """Main BEHAV3D dock widget registered as a napari plugin contribution."""

    def __init__(self, napari_viewer: "napari.Viewer", parent=None):
        super().__init__(parent)
        self.viewer = napari_viewer
        self.dev_mode = os.environ.get("BEHAV3D_DEV_MODE") == "1"
        self.setMinimumWidth(300)

        # Stop the mouse wheel from silently changing spin box / combo box
        # values when the user is just scrolling a tab past them.
        disable_spinbox_wheel_scroll()

        # --- Global QGroupBox styling -----------------------------------------
        # Qt's default QGroupBox leaves no room between the title and the top
        # border, so content rendered near the top of the box clips into the
        # title text.  margin-top reserves space *above* the border for the
        # title; padding-top pushes the content away from the border beneath it.
        self.setStyleSheet("""
            QGroupBox {
                border: 1px solid palette(mid);
                border-radius: 4px;
                margin-top: 1ex;
                padding-top: 8px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 8px;
                padding: 0 4px;
            }
        """)

        # --- Outer layout -------------------------------------------------
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)
        layout = outer_layout

        # --- Assistant show/hide toggle (Floating) ------------------------
        try:
            # Parent to BEHAV3DWidget itself so it overlaps this plugin and 
            # closes properly when this plugin closes.
            self.btn_toggle_assistant = FloatingAssistantButton(parent=self)
            self.btn_toggle_assistant.clicked.connect(self._toggle_assistant)
            self.btn_toggle_assistant.show()
        except Exception:

            # Fallback
            self.btn_toggle_assistant = QPushButton("💬 Assistant")
            self.btn_toggle_assistant.setCheckable(True)
            self.btn_toggle_assistant.setChecked(True)
            self.btn_toggle_assistant.clicked.connect(self._toggle_assistant)
            layout.addWidget(self.btn_toggle_assistant)


        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, stretch=1)

        # --- Tab 1: Data Preparation (metadata + dim order + zarr) --------
        from behav3d.napari._data_preparation import DataPreparationTab
        self.data_prep_tab = DataPreparationTab(parent=self)
        self.tabs.addTab(self.data_prep_tab, "📋 Data Preparation")

        # --- Tab 2: Visualization (dask-backed viewer) --------------------
        from behav3d.napari._visualization import VisualizationTab
        self.visualization_tab = VisualizationTab(viewer=self.viewer, data_prep=self.data_prep_tab, parent=self)
        self.tabs.addTab(self.visualization_tab, "👁 Visualization")

        # --- Tab 3: Segmentation ------------------------------------------
        from behav3d.napari._segmentation import SegmentationTab
        self.segmentation_tab = SegmentationTab(viewer=self.viewer, metadata_loader=self.data_prep_tab)
        self.tabs.addTab(self.segmentation_tab, "🦠 Segmentation")

        # --- Tab 4: Tracking ----------------------------------------------
        from behav3d.napari._tracking import TrackingTab
        self.tracking_tab = TrackingTab(viewer=self.viewer, metadata_loader=self.data_prep_tab)
        self.tabs.addTab(self.tracking_tab, "📍 Tracking")
        
        # Wire tracking completion signal to refresh metadata in data prep
        self.tracking_tab.tracking_completed.connect(self.data_prep_tab._on_tracking_completed)

        # --- Tab 5: Feature Extraction ------------------------------------
        from behav3d.napari._feature_extraction import FeatureExtractionTab
        self.feature_extraction_tab = FeatureExtractionTab(
            viewer=self.viewer, metadata_loader=self.data_prep_tab, parent=self
        )
        self.tabs.addTab(self.feature_extraction_tab, "🧪 Feature Extraction")

        # --- Tab 6: Filtering ---------------------------------------------
        from behav3d.napari._filtering import FilteringTab
        self.filtering_tab = FilteringTab(
            viewer=self.viewer, metadata_loader=self.data_prep_tab, parent=self
        )
        self.tabs.addTab(self.filtering_tab, "🧹 Filtering")

        # --- Tab 7: Analysis ----------------------------------------------
        from behav3d.napari._analysis import AnalysisTab
        self.analysis_tab = AnalysisTab(
            viewer=self.viewer,
            metadata_loader=self.data_prep_tab,
            parent=self,
        )
        self.tabs.addTab(self.analysis_tab, "📊 Analysis")

        # Backprojection is now integrated into the Single Cell tab
        # (Step 4 in State Classification, Step 5 in Track Classification)

        # --- Processing Queue (collapsible bottom panel) ------------------
        self.queue_panel = ProcessingQueuePanel(
            segmentation_tab=self.segmentation_tab,
            tracking_tab=self.tracking_tab,
            metadata_loader=self.data_prep_tab,
            feature_extraction_tab=self.feature_extraction_tab,
            filtering_tab=self.filtering_tab,
            analysis_tab=self.analysis_tab,
        )
        layout.addWidget(self.queue_panel)
        self.feature_extraction_tab.set_queue_panel(self.queue_panel)

        # Wire up +🛒 buttons to queue
        pc = self.segmentation_tab.pixel_classifier_page
        pc.btn_queue_train.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.TRAIN)
        )
        pc.btn_queue_segment.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.SEGMENT)
        )
        self.tracking_tab.btn_queue_track.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.TRACK)
        )
        self.feature_extraction_tab.btn_queue_feature.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.FEATURE_EXTRACT)
        )
        self.filtering_tab.btn_queue_filter.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.FILTER)
        )

        # Analysis tab +🛒 buttons (Death Dynamics + Interaction)
        dd = self.analysis_tab.death_dynamics_tab
        dd.btn_queue_dd_single.clicked.connect(self._add_death_dynamics_to_queue)
        dd.btn_queue_dd_combined.clicked.connect(self._add_multi_org_death_to_queue)
        dd.btn_queue_ia_single.clicked.connect(self._add_interaction_to_queue)
        dd.btn_queue_ia_combined.clicked.connect(self._add_multi_org_interaction_to_queue)
        dd.btn_queue_inv.clicked.connect(self._add_invasiveness_to_queue)
        dd.btn_queue_all.clicked.connect(self._add_all_analysis_to_queue)

        # Single Cell Analysis +🛒 buttons
        sc = self.analysis_tab.single_cell_tab
        sc.state_tab.btn_queue_state_cluster.clicked.connect(
            lambda: self._add_sc_step_to_queue(StepType.SC_STATE_CLUSTER)
        )
        sc.state_tab.btn_queue_train_state.clicked.connect(
            lambda: self._add_sc_step_to_queue(StepType.SC_TRAIN_STATE)
        )
        sc.state_tab.btn_queue_apply_state.clicked.connect(
            lambda: self._add_sc_step_to_queue(StepType.SC_APPLY_STATE)
        )
        sc.track_tab.btn_queue_track_cluster.clicked.connect(
            lambda: self._add_sc_step_to_queue(StepType.SC_TRACK_CLUSTER)
        )

        # Cellpose +🛒 button
        cp = self.segmentation_tab.cellpose_page
        cp.btn_queue_cellpose.clicked.connect(self._add_cellpose_to_queue)
        cp.btn_queue_otsu.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.DEAD_MASK)
        )

        # Cellpose-SAM +🛒 button
        cpsam = self.segmentation_tab.cellpose_sam_page
        cpsam.btn_queue_cellpose_sam.clicked.connect(self._add_cellpose_sam_to_queue)

        # APOC +🛒 buttons
        apoc = self.segmentation_tab.apoc_page
        apoc.btn_queue_apoc_segment.clicked.connect(self._add_apoc_segment_to_queue)

        # ConvPaint +🛒 buttons
        convpaint = self.segmentation_tab.convpaint_page
        convpaint.btn_queue_segment.clicked.connect(self._add_convpaint_segment_to_queue)

        # --- Tab Switch Logic ---------------------------------------------
        self._last_tab_index = self.tabs.currentIndex()
        self.tabs.currentChanged.connect(self._on_tab_changed)

        # --- Global Workers Controller ------------------------------------
        # Created AFTER all tabs so that every spin_workers spinbox already
        # exists when we link them.  The controller reads the saved
        # top-level ``n_workers`` value from behav3d_parameters.yml (via
        # data_prep_tab.behav3d_parameters) and propagates it to all spinboxes.
        self.workers_ctrl = GlobalWorkersController(metadata_loader=self.data_prep_tab)
        # Segmentation tab — pixel classifier spin_workers
        seg_pc = self.segmentation_tab.pixel_classifier_page
        if hasattr(seg_pc, "spin_workers"):
            self.workers_ctrl.link(seg_pc.spin_workers)
        # Segmentation tab — APOC spin_workers
        seg_apoc = self.segmentation_tab.apoc_page
        if hasattr(seg_apoc, "spin_workers"):
            self.workers_ctrl.link(seg_apoc.spin_workers)
        # Segmentation tab — ConvPaint spin_workers
        seg_cp = self.segmentation_tab.convpaint_page
        if hasattr(seg_cp, "spin_workers"):
            self.workers_ctrl.link(seg_cp.spin_workers)
        # Data Preparation tab — zarr conversion spin_workers
        if hasattr(self.data_prep_tab, "zarr_workers_spin"):
            self.workers_ctrl.link(self.data_prep_tab.zarr_workers_spin)
        # Feature extraction tab
        if hasattr(self.feature_extraction_tab, "spin_workers"):
            self.workers_ctrl.link(self.feature_extraction_tab.spin_workers)
        # When a new parameters file is loaded, reload the global value.
        self.data_prep_tab.metadata_loaded.connect(
            lambda _: self.workers_ctrl.reload()
        )
        # --- Co-pilot Assistant dock (right side) -------------------------
        # Added as a *separate* dock so it stays visible across all tabs. Wrapped
        # in try/except so an assistant failure can never block the main UI.
        self.assistant = None
        self._assistant_dock = None
        try:
            from behav3d.napari._assistant import AssistantDock
            self.assistant = AssistantDock(main_widget=self)
            # self._assistant_dock = self.viewer.window.add_dock_widget(
            #     self.assistant, area="right", name="BEHAV3D Assistant"
            # )
            # # Remove the native 'x' button so napari doesn't destroy the widget
            # from qtpy.QtWidgets import QDockWidget
            # self._assistant_dock.setFeatures(
            #     self._assistant_dock.features() & ~QDockWidget.DockWidgetClosable
            # )

            from napari._qt.widgets.qt_viewer_dock_widget import QtViewerDockWidget

            self._assistant_dock = QtViewerDockWidget(
                self.viewer.window._qt_viewer,
                self.assistant,
                name="BEHAV3D Assistant",
                area="right",
                close_btn=False,   # <-- this is the one that actually works
            )
            self.viewer.window._add_viewer_dock_widget(self._assistant_dock, menu=self.viewer.window.window_menu)
            
            # Both this pipeline and the assistant dock to the "right" area, so Qt
            # stacks them vertically by default. Defer to the next event-loop tick
            # (by which point this widget's own dock exists, for either launch
            # path) and re-split them *horizontally* → pipeline | assistant.
            from qtpy.QtCore import QTimer
            QTimer.singleShot(0, self._place_assistant_beside)
            # Keep the toggle button in sync if the dock is shown/hidden by other
            # means (e.g. its native close button → hide, not destroy).
            try:
                self._assistant_dock.visibilityChanged.connect(
                    self._on_assistant_visibility_changed
                )
            except Exception:
                pass
            # Keep the assistant's context bar in sync with workflow state.
            self.tabs.currentChanged.connect(
                lambda *_: self.assistant.refresh_context_bar()
            )
            if hasattr(self.data_prep_tab, "metadata_loaded"):
                self.data_prep_tab.metadata_loaded.connect(
                    lambda *_: self.assistant.refresh_context_bar()
                )
        except Exception as e:  # pragma: no cover - defensive
            import warnings
            warnings.warn(f"BEHAV3D Assistant could not be initialised: {e}")
            # No dock to toggle — disable the button so it isn't a dead control.
            try:
                self.btn_toggle_assistant.setEnabled(False)
            except Exception:
                pass

    def _toggle_assistant(self):
        """Show or hide the assistant dock. Closing only hides the dock, so the
        conversation is preserved and the panel can be reopened instantly."""
        dock = self._assistant_dock
        if dock is None:
            return
        try:
            if hasattr(dock, "toggleViewAction"):
                # Use the native QAction trigger to cleanly toggle the dock. 
                # This works even if Napari's internal layout manager hid it via the 'X' button.
                dock.toggleViewAction().trigger()
            else:
                visible = dock.isVisible()
                dock.setVisible(not visible)
            
            if dock.isVisible():
                dock.raise_()
        except Exception:
            pass

    def _on_assistant_visibility_changed(self, visible: bool):
        """Mirror the dock's visibility onto the toggle button without re-emitting."""
        try:
            self.btn_toggle_assistant.blockSignals(True)
            self.btn_toggle_assistant.setChecked(bool(visible))
        finally:
            self.btn_toggle_assistant.blockSignals(False)

    def _place_assistant_beside(self):
        """Move the assistant dock to sit side-by-side (horizontally) with this
        pipeline dock instead of stacked above it. Best-effort; never fatal."""
        try:
            from qtpy.QtWidgets import QDockWidget, QMainWindow
            from qtpy.QtCore import Qt
            if self._assistant_dock is None:
                return
            # Find the QMainWindow and this pipeline's enclosing QDockWidget.
            qmain = self._assistant_dock.parent()
            while qmain is not None and not isinstance(qmain, QMainWindow):
                qmain = qmain.parent()
            if qmain is None:
                return
            my_dock = None
            for dock in qmain.findChildren(QDockWidget):
                w = dock.widget()
                # the pipeline dock wraps `self` (possibly via a layout wrapper)
                if w is self or (w is not None and self in w.findChildren(type(self))):
                    if dock is not self._assistant_dock:
                        my_dock = dock
                        break
            if my_dock is not None:
                qmain.splitDockWidget(my_dock, self._assistant_dock, Qt.Horizontal)
        except Exception:
            pass

    def _add_cellpose_to_queue(self):
        """Validate cellpose config and add a CELLPOSE_SEGMENT step to the queue."""
        cp = self.segmentation_tab.cellpose_page
        ok, msg = cp.validate_for_queue()
        if not ok:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Cannot Add to Queue", msg)
            return
        params = cp.get_queue_params()
        self.queue_panel.add_step(StepType.CELLPOSE_SEGMENT, params=params)

    def _add_cellpose_sam_to_queue(self):
        """Validate Cellpose-SAM config and add a CELLPOSE_SAM_SEGMENT step to the queue."""
        cpsam = self.segmentation_tab.cellpose_sam_page
        ok, msg = cpsam.validate_for_queue()
        if not ok:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Cannot Add to Queue", msg)
            return
        params = cpsam.get_queue_params()
        self.queue_panel.add_step(StepType.CELLPOSE_SAM_SEGMENT, params=params)

    def _add_apoc_segment_to_queue(self):
        """Validate APOC config and add an APOC_SEGMENT step to the queue."""
        apoc = self.segmentation_tab.apoc_page
        ok, msg = apoc.validate_for_queue()
        if not ok:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Cannot Add to Queue", msg)
            return
        params = apoc.get_queue_params()
        self.queue_panel.add_step(StepType.APOC_SEGMENT, params=params)

    def _add_convpaint_segment_to_queue(self):
        """Validate ConvPaint config and add a CONVPAINT_SEGMENT step to the queue."""
        convpaint = self.segmentation_tab.convpaint_page
        ok, msg = convpaint.validate_for_queue()
        if not ok:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Cannot Add to Queue", msg)
            return
        params = convpaint.get_queue_params()
        self.queue_panel.add_step(StepType.CONVPAINT_SEGMENT, params=params)

    # ── Analysis tab queue helpers ─────────────────────────────────────
    def _add_death_dynamics_to_queue(self):
        dd = self.analysis_tab.death_dynamics_tab
        targets = dd._selected_targets()
        if not targets:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Cannot Add to Queue", "Select at least one target cell type."
            )
            return
        self.queue_panel.add_step(
            StepType.DEATH_DYNAMICS, params={"cell_types": list(targets)}
        )

    def _add_multi_org_death_to_queue(self):
        dd = self.analysis_tab.death_dynamics_tab
        targets = dd._selected_targets()
        if len(targets) < 2:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Cannot Add to Queue",
                "Select at least two targets for combined death dynamics."
            )
            return
        self.queue_panel.add_step(
            StepType.MULTI_ORG_DEATH, params={"cell_types": list(targets)}
        )

    def _add_interaction_to_queue(self):
        dd = self.analysis_tab.death_dynamics_tab
        targets = dd._selected_targets()
        interactions = dd._selected_interactions()
        if not targets or not interactions:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Cannot Add to Queue",
                "Select at least one target and one interaction cell type."
            )
            return
        dd._persist_advanced()
        self.queue_panel.add_step(
            StepType.INTERACTION,
            params={
                "cell_types": list(targets),
                "interaction_cell_types": list(interactions),
                "group_by_line_condition": bool(dd.check_group_by_lc.isChecked()),
            },
        )

    def _add_multi_org_interaction_to_queue(self):
        dd = self.analysis_tab.death_dynamics_tab
        targets = dd._selected_targets()
        interactions = dd._selected_interactions()
        if not targets or not interactions:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Cannot Add to Queue",
                "Select at least one target and one interaction cell type."
            )
            return
        # Persist advanced settings before queueing
        dd._persist_advanced()
        self.queue_panel.add_step(
            StepType.MULTI_ORG_INTERACTION,
            params={
                "cell_types": list(targets),
                "interaction_cell_types": list(interactions),
                "time_window_min": int(dd.spin_time_window.value()),
                "group_by": dd.combo_group_by.currentData() or "organoid_type",
                "annotate_line_condition": bool(dd.check_annotate_lc.isChecked()),
                "analysis_period_t": _period_from_t_radios(
                    dd.period_radio_group, dd.spin_period_start_t, dd.spin_period_end_t,
                ),
            },
        )

    def _add_invasiveness_to_queue(self):
        dd = self.analysis_tab.death_dynamics_tab
        immune_list = dd._selected_invasiveness_immune()
        targets = dd._selected_invasiveness_targets()
        if not immune_list or not targets:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(
                self, "Cannot Add to Queue",
                "Select at least one immune cell type and one target."
            )
            return
        self.queue_panel.add_step(
            StepType.INVASIVENESS,
            params={
                "immune_cell_types": list(immune_list),
                "targets": list(targets),
                "summary_stat": dd.inv_summary_combo.currentData() or "mean",
                "group_by_line_condition": bool(dd.check_inv_group_by_lc.isChecked()),
                "analysis_period_t": _period_from_t_radios(
                    dd.inv_period_radio_group, dd.spin_inv_start_t, dd.spin_inv_end_t,
                ),
            },
        )

    def _add_all_analysis_to_queue(self):
        """Queue every applicable Death Dynamics / Interaction step for the
        current Analysis-tab selection — mirrors ``_on_run_all_clicked``.

        Queue-time checks are deliberately *optimistic*: we only inspect
        static conditions (dead channel configured in metadata + selection
        counts). Intermediate-file availability (filtered CSV, dead
        column, contact columns) is checked by the per-step runners just
        before they execute, so users can chain a fresh
        Filter → Death Dynamics → Interaction pipeline in one click even
        when no analysis outputs exist yet.
        """
        from qtpy.QtWidgets import QMessageBox

        dd = self.analysis_tab.death_dynamics_tab
        targets = dd._selected_targets()
        interactions = dd._selected_interactions()
        if not targets:
            QMessageBox.warning(
                self, "Cannot Add to Queue",
                "Select at least one target cell type."
            )
            return

        # Static metadata check — does NOT depend on filtered CSV presence.
        has_dead = bool(
            dd.metadata_loader is not None
            and dd.metadata_loader.metadata is not None
            and "dead_channel" in dd.metadata_loader.metadata.columns
            and dd.metadata_loader.metadata["dead_channel"].notna().any()
        )

        added_any = False
        if has_dead:
            self.queue_panel.add_step(
                StepType.DEATH_DYNAMICS,
                params={"cell_types": list(targets)},
            )
            added_any = True
            if len(targets) >= 2:
                self.queue_panel.add_step(
                    StepType.MULTI_ORG_DEATH,
                    params={"cell_types": list(targets)},
                )

        if interactions:
            dd._persist_advanced()
            self.queue_panel.add_step(
                StepType.INTERACTION,
                params={
                    "cell_types": list(targets),
                    "interaction_cell_types": list(interactions),
                    "group_by_line_condition": bool(dd.check_group_by_lc.isChecked()),
                },
            )
            added_any = True
            self.queue_panel.add_step(
                StepType.MULTI_ORG_INTERACTION,
                params={
                    "cell_types": list(targets),
                    "interaction_cell_types": list(interactions),
                    "time_window_min": int(dd.spin_time_window.value()),
                    "group_by": (
                        dd.combo_group_by.currentData() or "organoid_type"
                    ),
                    "annotate_line_condition": bool(
                        dd.check_annotate_lc.isChecked()
                    ),
                    "analysis_period_t": _period_from_t_radios(
                        dd.period_radio_group,
                        dd.spin_period_start_t,
                        dd.spin_period_end_t,
                    ),
                },
            )

        # Invasiveness (immune-cell perspective) uses its own selectors,
        # independent of the target/interaction pickers above.
        inv_immune = dd._selected_invasiveness_immune()
        inv_targets = dd._selected_invasiveness_targets()
        if inv_immune and inv_targets:
            self.queue_panel.add_step(
                StepType.INVASIVENESS,
                params={
                    "immune_cell_types": list(inv_immune),
                    "targets": list(inv_targets),
                    "summary_stat": dd.inv_summary_combo.currentData() or "mean",
                    "group_by_line_condition": bool(
                        dd.check_inv_group_by_lc.isChecked()
                    ),
                    "analysis_period_t": _period_from_t_radios(
                        dd.inv_period_radio_group,
                        dd.spin_inv_start_t,
                        dd.spin_inv_end_t,
                    ),
                },
            )
            added_any = True

        if not added_any:
            QMessageBox.warning(
                self, "Nothing to Queue",
                "No applicable steps for the current selection. Either "
                "configure a dead channel in metadata or select at least "
                "one interaction cell type."
            )

    # ── Single Cell queue helpers ───────────────────────────────────────
    def _add_sc_step_to_queue(self, step_type):
        """Add a single-cell analysis step to the queue, tagging the
        currently selected cell type in the SingleCellTab."""
        try:
            sc = self.analysis_tab.single_cell_tab
            ct = sc.cell_type_combo.currentText()
        except Exception:
            ct = ""
        params = {"cell_type": ct} if ct else {}
        self.queue_panel.add_step(step_type, params=params)

    def _on_tab_changed(self, index):
        """Intercept tab switches to handle exit warnings and missing output dir."""
        # 1. Handle missing output directory (indices >= 2: Segmentation, Tracking, etc.)
        if index >= 2:
            has_metadata = self.data_prep_tab.metadata is not None
            has_out_dir = bool(self.data_prep_tab.output_dir)
            
            if has_metadata and not has_out_dir:
                from qtpy.QtWidgets import QMessageBox
                QMessageBox.warning(
                    self, 
                    "Output Directory Required",
                    "No output directory has been established. Please go to the Data Preparation tab to define it."
                )
                # Switch back to last tab (likely Data Prep or Visualization)
                self.tabs.blockSignals(True)
                self.tabs.setCurrentIndex(self._last_tab_index)
                self.tabs.blockSignals(False)
                return

        # 2. Handle exit logic (e.g. Segmentation Tab training session)
        if self._last_tab_index == 2:
            if hasattr(self, 'segmentation_tab'):
                if not self.segmentation_tab.request_tab_exit():
                    # If user chose 'No', switch back to Segmentation Tab
                    self.tabs.blockSignals(True)
                    self.tabs.setCurrentIndex(2)
                    self.tabs.blockSignals(False)
                    return

        # 3. Handle exit from the Visualization tab while a manual-correction
        # editor has unsaved changes — prompt Save / Discard / Cancel.
        if self._last_tab_index == 1:
            if hasattr(self, 'visualization_tab') and hasattr(
                self.visualization_tab, 'request_tab_exit'
            ):
                if not self.visualization_tab.request_tab_exit():
                    self.tabs.blockSignals(True)
                    self.tabs.setCurrentIndex(1)
                    self.tabs.blockSignals(False)
                    return

        # 4. Generic background-operation guards for tabs that adopted
        # the BackgroundOperation pattern (Tracking, Feature Extraction,
        # Filtering).  Each tab exposes ``request_tab_exit`` returning
        # False when it has work in flight that should block tab switch.
        guarded_attrs = {
            3: "tracking_tab",
            4: "feature_extraction_tab",
            5: "filtering_tab",
            6: "analysis_tab",
        }
        guarded_attr = guarded_attrs.get(self._last_tab_index)
        if guarded_attr is not None and hasattr(self, guarded_attr):
            tab_obj = getattr(self, guarded_attr)
            if hasattr(tab_obj, "request_tab_exit") and not tab_obj.request_tab_exit():
                self.tabs.blockSignals(True)
                self.tabs.setCurrentIndex(self._last_tab_index)
                self.tabs.blockSignals(False)
                return

        self._last_tab_index = index

        # Auto-refresh Analysis tab when switched to
        if index == 6 and hasattr(self, 'analysis_tab'):
            if hasattr(self.analysis_tab, '_on_metadata_updated'):
                self.analysis_tab._on_metadata_updated()
            # Also notify inner single cell tab to reload (picks up Track tab auto-fill)
            if hasattr(self.analysis_tab, 'single_cell_tab') and hasattr(self.analysis_tab.single_cell_tab, '_on_metadata_updated'):
                self.analysis_tab.single_cell_tab._on_metadata_updated()

    def sizeHint(self):
        return QSize(440, 650)