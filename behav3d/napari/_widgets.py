"""
Reusable widgets for the BEHAV3D napari plugin.
"""
from qtpy.QtWidgets import QPushButton, QMessageBox, QHBoxLayout, QWidget
from qtpy.QtCore import Qt


class HelpButton(QPushButton):
    """Small circular '?' button that shows a help popup on click."""

    def __init__(self, title: str, description: str, parent=None):
        super().__init__("?", parent)
        self._title = title
        self._description = description

        self.setFixedSize(20, 20)
        self.setCursor(Qt.WhatsThisCursor)
        self.setStyleSheet(
            "QPushButton {"
            "  background-color: #5a9bd5;"
            "  color: white;"
            "  border-radius: 10px;"
            "  font-weight: bold;"
            "  font-size: 12px;"
            "  padding: 0px;"
            "}"
            "QPushButton:hover {"
            "  background-color: #3a7bc8;"
            "}"
        )
        self.clicked.connect(self._show_help)

    def _show_help(self):
        QMessageBox.information(self, self._title, self._description)


def make_help_row(widget, title: str, description: str) -> QHBoxLayout:
    """Return an HBoxLayout containing *widget* and a HelpButton."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.addWidget(widget, stretch=1)
    row.addWidget(HelpButton(title, description))
    return row
