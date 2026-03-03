"""
BEHAV3D napari plugin – main dock widget.
Provides a QTabWidget with tabs for the full BEHAV3D pipeline.
"""
from qtpy.QtWidgets import QWidget, QVBoxLayout, QTabWidget, QScrollArea
from qtpy.QtCore import Qt
import napari


class BEHAV3DWidget(QWidget):
    """Main BEHAV3D dock widget registered as a napari plugin contribution."""

    def __init__(self, napari_viewer: "napari.Viewer", parent=None):
        super().__init__(parent)
        self.viewer = napari_viewer
        self.setMinimumWidth(400)

        # --- Outer layout with scroll area --------------------------------
        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer_layout.addWidget(scroll)

        # --- Inner content widget -----------------------------------------
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        scroll.setWidget(content)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

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

        # --- Tab 5–7: Additional Steps (Stubs) ---------------------------
        from behav3d.napari._stubs import (
            FeatureExtractionTab,
            FilteringTab,
            AnalysisTab
        )
        self.tabs.addTab(FeatureExtractionTab(parent=self), "🧪 Feature Extraction")
        self.tabs.addTab(FilteringTab(parent=self), "🧹 Filtering")
        self.tabs.addTab(AnalysisTab(parent=self), "📊 Analysis")

        # --- Processing Queue (collapsible bottom panel) ------------------
        from behav3d.napari._queue import ProcessingQueuePanel, StepType
        self.queue_panel = ProcessingQueuePanel(
            segmentation_tab=self.segmentation_tab,
            tracking_tab=self.tracking_tab,
            metadata_loader=self.data_prep_tab,
        )
        layout.addWidget(self.queue_panel)

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

        # --- Tab Switch Logic ---------------------------------------------
        self._last_tab_index = self.tabs.currentIndex()
        self.tabs.currentChanged.connect(self._on_tab_changed)

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

        self._last_tab_index = index
