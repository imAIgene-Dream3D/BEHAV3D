"""
BEHAV3D napari plugin – main dock widget.
Provides a QTabWidget with tabs for the full BEHAV3D pipeline.
"""
from qtpy.QtWidgets import QWidget, QVBoxLayout, QTabWidget
from qtpy.QtCore import Qt, QSize
import napari

from behav3d.napari._queue import ProcessingQueuePanel, StepType


class BEHAV3DWidget(QWidget):
    """Main BEHAV3D dock widget registered as a napari plugin contribution."""

    def __init__(self, napari_viewer: "napari.Viewer", parent=None):
        super().__init__(parent)
        self.viewer = napari_viewer
        self.setMinimumWidth(300)

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

        # --- Tab 7: Analysis (Stub) ---------------------------------------
        from behav3d.napari._stubs import AnalysisTab
        self.tabs.addTab(AnalysisTab(parent=self), "📊 Analysis")

        # --- Processing Queue (collapsible bottom panel) ------------------
        self.queue_panel = ProcessingQueuePanel(
            segmentation_tab=self.segmentation_tab,
            tracking_tab=self.tracking_tab,
            metadata_loader=self.data_prep_tab,
            feature_extraction_tab=self.feature_extraction_tab,
            filtering_tab=self.filtering_tab,
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

        # Cellpose +🛒 button
        cp = self.segmentation_tab.cellpose_page
        cp.btn_queue_cellpose.clicked.connect(self._add_cellpose_to_queue)
        cp.btn_queue_otsu.clicked.connect(
            lambda: self.queue_panel.add_step(StepType.DEAD_MASK)
        )

        # APOC +🛒 buttons
        apoc = self.segmentation_tab.apoc_page
        apoc.btn_queue_apoc_segment.clicked.connect(self._add_apoc_segment_to_queue)

        # ConvPaint +🛒 buttons
        convpaint = self.segmentation_tab.convpaint_page
        convpaint.btn_queue_segment.clicked.connect(self._add_convpaint_segment_to_queue)

        # --- Tab Switch Logic ---------------------------------------------
        self._last_tab_index = self.tabs.currentIndex()
        self.tabs.currentChanged.connect(self._on_tab_changed)

        # --- Co-pilot Assistant dock (right side) -------------------------
        # Added as a *separate* dock so it stays visible across all tabs. Wrapped
        # in try/except so an assistant failure can never block the main UI.
        self.assistant = None
        self._assistant_dock = None
        try:
            from behav3d.napari._assistant import AssistantDock
            self.assistant = AssistantDock(main_widget=self)
            self._assistant_dock = self.viewer.window.add_dock_widget(
                self.assistant, area="right", name="🤖 BEHAV3D Assistant"
            )
            # Both this pipeline and the assistant dock to the "right" area, so Qt
            # stacks them vertically by default. Defer to the next event-loop tick
            # (by which point this widget's own dock exists, for either launch
            # path) and re-split them *horizontally* → pipeline | assistant.
            from qtpy.QtCore import QTimer
            QTimer.singleShot(0, self._place_assistant_beside)
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

        self._last_tab_index = index

    def sizeHint(self):
        return QSize(440, 650)
