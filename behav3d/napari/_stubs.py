"""
Placeholder tab widgets for BEHAV3D napari plugin.
These will be replaced with full implementations as each module is ported.
"""
from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel
from qtpy.QtCore import Qt


class _StubTab(QWidget):
    """Base placeholder for tabs not yet implemented."""

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignCenter)

        icon = QLabel("🚧")
        icon.setStyleSheet("font-size: 48px;")
        icon.setAlignment(Qt.AlignCenter)
        layout.addWidget(icon)

        label = QLabel(f"{title}")
        label.setStyleSheet("font-size: 18px; font-weight: bold; color: #555;")
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        sub = QLabel("Coming soon — this tab is not yet implemented.")
        sub.setStyleSheet("font-size: 12px; color: #999;")
        sub.setAlignment(Qt.AlignCenter)
        layout.addWidget(sub)


class SegmentationTab(_StubTab):
    def __init__(self, parent=None):
        super().__init__("Segmentation", parent)


class TrackingTab(_StubTab):
    def __init__(self, parent=None):
        super().__init__("Tracking", parent)


class SingleCellAnalysisTab(_StubTab):
    def __init__(self, parent=None):
        super().__init__("Single-Cell Analysis", parent)


class ApplicationTab(_StubTab):
    def __init__(self, parent=None):
        super().__init__("Application", parent)
