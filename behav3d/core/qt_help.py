"""Shared Qt helpers for parameter-clarification "?" icons.

This module provides a small, framework-free Qt helper (``HelpButton``) and a
companion layout helper (``make_help_row``) used by both:

- the napari plugin (``behav3d.napari._widgets`` re-exports them), and
- the standalone Qt training widgets that live in
  ``behav3d.preprocessing.segmentation`` (e.g. ``apoc_train``,
  ``convpaint_train``).

Keeping the implementation here (rather than under ``behav3d.napari``) lets
``preprocessing.segmentation`` use the same "?" button without depending on
napari, which is important because those widgets are also opened from the
notebooks where ``behav3d.napari`` should not necessarily be imported.
"""

from qtpy.QtWidgets import QPushButton, QMessageBox, QHBoxLayout, QWidget
from qtpy.QtCore import Qt


class HelpButton(QPushButton):
    """Small circular '?' button that shows a help popup on click."""

    def __init__(self, title: str, description: str, parent=None):
        super().__init__("?", parent)
        self._title = title
        self._description = description

        self.setFixedSize(15, 15)
        self.setCursor(Qt.WhatsThisCursor)
        self.setStyleSheet(
            "QPushButton {"
            "  background-color: #5a9bd5;"
            "  color: white;"
            "  border-radius: 7px;"
            "  font-weight: bold;"
            "  font-size: 10px;"
            "  padding: 0px;"
            "}"
            "QPushButton:hover {"
            "  background-color: #3a7bc8;"
            "}"
        )
        self.clicked.connect(self._show_help)

    def set_help(self, title: str, description: str) -> None:
        """Replace the title/description shown when the user clicks the button."""
        self._title = title
        self._description = description

    def _show_help(self):
        QMessageBox.information(self, self._title, self._description)


def make_help_row(widget: QWidget, title: str, description: str) -> QHBoxLayout:
    """Return an HBoxLayout containing ``widget`` and a HelpButton."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.addWidget(widget, stretch=1)
    row.addWidget(HelpButton(title, description))
    return row


__all__ = ["HelpButton", "make_help_row"]
