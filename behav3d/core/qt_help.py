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

It also provides :func:`disable_spinbox_wheel_scroll`, a fix for accidental
spin-box/combo-box value changes while scrolling past them.
"""

from qtpy.QtWidgets import (
    QPushButton,
    QMessageBox,
    QHBoxLayout,
    QWidget,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QScrollArea,
    QApplication,
)
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


_wheel_scroll_disabled = False


def disable_spinbox_wheel_scroll() -> None:
    """Make QSpinBox/QDoubleSpinBox/QComboBox ignore mouse-wheel events
    instead of changing their value, application-wide and idempotently.

    Without this, scrolling down a tab that contains many parameter spin
    boxes silently changes whichever box the cursor passes over, which is
    easy to miss and hard to notice after the fact.

    The event is not simply dropped: it is forwarded to the nearest
    ``QScrollArea`` ancestor's viewport so the page keeps scrolling right
    through the box, instead of the scroll gesture going dead over it.

    Two things this deliberately avoids, learned from two earlier attempts
    at this fix:

    - It does NOT touch ``QAbstractSlider``/``QScrollBar`` — ``QScrollBar``
      *is* a ``QAbstractSlider``, and ``QScrollArea`` scrolls internally by
      forwarding wheel events to its own ``QScrollBar``. An earlier attempt
      that also intercepted ``QAbstractSlider`` swallowed those internal
      scrollbar events too, breaking scrolling everywhere, not just over
      spin/combo boxes.
    - It does NOT do the forwarding from a global, application-installed
      event filter. An earlier attempt did, and because the filter is
      global it also intercepted the very event it re-sent to the scroll
      area's viewport (which forwards to the scrollbar internally), causing
      an infinite viewport <-> scrollbar bounce and a native stack overflow
      crash. Forwarding instead from each widget's own (patched)
      ``wheelEvent`` is safe: the forward target (a ``QScrollArea``
      viewport) is never one of the patched classes, so it can't re-enter
      this code.
    """
    global _wheel_scroll_disabled
    if _wheel_scroll_disabled:
        return
    def _ignore_wheel(self, event):
        event.ignore()
        parent = self.parentWidget()
        while parent is not None and not isinstance(parent, QScrollArea):
            parent = parent.parentWidget()
        if parent is not None:
            QApplication.sendEvent(parent.viewport(), event)
    for cls in (QSpinBox, QDoubleSpinBox, QComboBox):
        cls.wheelEvent = _ignore_wheel
    _wheel_scroll_disabled = True


__all__ = ["HelpButton", "make_help_row", "disable_spinbox_wheel_scroll"]
