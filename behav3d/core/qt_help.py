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

from pathlib import Path

from qtpy.QtWidgets import (
    QPushButton,
    QMessageBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
    QWidget,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QScrollArea,
    QApplication,
    QStackedWidget,
    QTabWidget,
)
from qtpy.QtGui import QPixmap
from qtpy.QtCore import Qt
from typing import NamedTuple, Optional, Sequence, Union

# Figures live alongside the tutorial diagrams, resolved by plain filesystem
# path relative to this file (same convention as
# ``behav3d.napari._tutorial_dialog``).
_RESOURCES_DIR = Path(__file__).resolve().parents[1] / "resources"


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


class HelpSection(NamedTuple):
    """One block of an :class:`IllustratedHelpButton` popup.

    ``image`` is a path relative to ``behav3d/resources/`` (e.g.
    ``"tracking/bounded_propagation.png"``) or ``None`` for a text-only
    section. A figure that is missing from disk degrades to text only, so
    the dialog still works in a checkout where the images were never added.
    """

    heading: str
    text: str
    image: Optional[str] = None


class HelpDialog(QDialog):
    """Scrollable help window: heading + text + optional figure per section.

    The text-only ``QMessageBox`` used by :class:`HelpButton` cannot show
    figures at a usable size, which is the whole point for topics (like the
    propagation tracking variants) where a diagram explains the mechanism
    far better than a paragraph does.
    """

    def __init__(
        self,
        title: str,
        sections: Sequence[HelpSection],
        parent=None,
        image_width: int = 760,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self._image_width = int(image_width)

        outer = QVBoxLayout(self)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        body_lay = QVBoxLayout(body)
        body_lay.setSpacing(6)

        for i, section in enumerate(sections):
            if i:
                rule = QFrame()
                rule.setFrameShape(QFrame.HLine)
                rule.setStyleSheet("color:#444;")
                body_lay.addSpacing(6)
                body_lay.addWidget(rule)

            heading = QLabel(section.heading)
            heading.setStyleSheet("font-weight:bold; font-size:13px;")
            heading.setWordWrap(True)
            body_lay.addWidget(heading)

            if section.text:
                text = QLabel(section.text)
                text.setWordWrap(True)
                body_lay.addWidget(text)

            pixmap = self._load_image(section.image)
            if pixmap is not None:
                figure = QLabel()
                figure.setAlignment(Qt.AlignCenter)
                figure.setPixmap(pixmap)
                body_lay.addWidget(figure)

        body_lay.addStretch(1)
        scroll.setWidget(body)
        outer.addWidget(scroll)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        outer.addLayout(buttons)

        self.resize(self._image_width + 80, 640)

    def _load_image(self, name: Optional[str]) -> Optional[QPixmap]:
        """Load a figure by resources-relative name, scaled down to fit.

        Returns ``None`` when there is no image or the file is absent /
        unreadable, so a caller never has to guard the call.
        """
        if not name:
            return None
        pixmap = QPixmap(str(_RESOURCES_DIR / name))
        if pixmap.isNull():
            return None
        if pixmap.width() > self._image_width:
            # Only ever downscale; upscaling a small diagram just blurs it.
            pixmap = pixmap.scaledToWidth(self._image_width, Qt.SmoothTransformation)
        return pixmap


class IllustratedHelpButton(HelpButton):
    """``HelpButton`` whose popup can carry figures as well as text."""

    def __init__(
        self,
        title: str,
        sections: Sequence[HelpSection],
        parent=None,
        image_width: int = 760,
    ):
        # Give the base class a plain-text rendering too, so anything that
        # reads ``_description`` (or a future text-only fallback) still works.
        super().__init__(
            title,
            "\n\n".join(
                f"{s.heading} - {s.text}" if s.text else s.heading for s in sections
            ),
            parent,
        )
        self._sections = list(sections)
        self._image_width = int(image_width)

    def _show_help(self):
        HelpDialog(
            self._title, self._sections, parent=self, image_width=self._image_width
        ).exec_()


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


def reset_scroll_on_page_change(container: Union[QStackedWidget, QTabWidget]) -> None:
    """Scroll any ``QScrollArea`` inside a page back to the top whenever
    ``container`` (a ``QStackedWidget`` or ``QTabWidget``) switches to it.

    Qt keeps each page's scroll position between visits since the page
    widgets are never recreated, only hidden/shown. Without this, a page
    left scrolled down (e.g. a settings form) can appear empty the next
    time it's shown, because the visible viewport lands mid-content
    instead of at the top.
    """
    def _reset(index: int) -> None:
        page = container.widget(index)
        if page is None:
            return
        scrolls = page.findChildren(QScrollArea)
        if isinstance(page, QScrollArea):
            scrolls = [page, *scrolls]
        for scroll in scrolls:
            scroll.verticalScrollBar().setValue(0)

    container.currentChanged.connect(_reset)


__all__ = [
    "HelpButton",
    "HelpSection",
    "HelpDialog",
    "IllustratedHelpButton",
    "make_help_row",
    "disable_spinbox_wheel_scroll",
    "reset_scroll_on_page_change",
]
