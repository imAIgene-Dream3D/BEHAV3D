"""Step-by-step image-based tutorial dialog.

Sibling to :mod:`behav3d.core.qt_help`'s text-only ``HelpButton`` — this
variant shows a small diagram + caption per step with Next/Back
navigation, for workflows (like manual-editing seed placement) that users
found hard to follow from tooltips alone.

Diagram images are loaded from ``behav3d/resources/tutorials/`` following
the same convention as ``behav3d/napari/_widget.py``'s assistant icon: a
plain filesystem path resolved relative to this file, with a graceful
text-only fallback if an image is missing.
"""
from __future__ import annotations

from pathlib import Path
from typing import NamedTuple, Sequence

from qtpy.QtCore import Qt
from qtpy.QtGui import QPixmap
from qtpy.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
)

_TUTORIALS_DIR = Path(__file__).resolve().parents[1] / "resources" / "tutorials"


class TutorialStep(NamedTuple):
    heading: str
    image_name: str
    caption: str


class TutorialDialog(QDialog):
    """Modal wizard showing one diagram + caption per step."""

    def __init__(self, title: str, steps: Sequence[TutorialStep], parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(440)
        self._steps = list(steps)
        self._index = 0

        layout = QVBoxLayout(self)

        self._step_label = QLabel()
        self._step_label.setStyleSheet("color:#888; font-size:10px;")
        layout.addWidget(self._step_label)

        self._heading_label = QLabel()
        self._heading_label.setStyleSheet("font-weight:bold; font-size:13px;")
        self._heading_label.setWordWrap(True)
        layout.addWidget(self._heading_label)

        self._image_label = QLabel()
        self._image_label.setAlignment(Qt.AlignCenter)
        self._image_label.setMinimumHeight(220)
        self._image_label.setStyleSheet("color:#999; font-style:italic;")
        layout.addWidget(self._image_label)

        self._caption_label = QLabel()
        self._caption_label.setWordWrap(True)
        self._caption_label.setStyleSheet("color:#555;")
        layout.addWidget(self._caption_label)

        nav = QHBoxLayout()
        self._btn_back = QPushButton("< Back")
        self._btn_back.clicked.connect(self._go_back)
        self._btn_next = QPushButton("Next >")
        self._btn_next.clicked.connect(self._go_next)
        self._btn_close = QPushButton("Close")
        self._btn_close.clicked.connect(self.close)
        nav.addWidget(self._btn_back)
        nav.addWidget(self._btn_next)
        nav.addStretch(1)
        nav.addWidget(self._btn_close)
        layout.addLayout(nav)

        self._refresh()

    def _refresh(self) -> None:
        n = len(self._steps)
        step = self._steps[self._index]
        self._step_label.setText(f"Step {self._index + 1} of {n}")
        self._heading_label.setText(step.heading)
        self._caption_label.setText(step.caption)

        pix = QPixmap(str(_TUTORIALS_DIR / step.image_name))
        if not pix.isNull():
            self._image_label.setPixmap(
                pix.scaledToWidth(400, Qt.SmoothTransformation)
            )
            self._image_label.setText("")
        else:
            self._image_label.setPixmap(QPixmap())
            self._image_label.setText("(diagram unavailable)")

        self._btn_back.setEnabled(self._index > 0)
        self._btn_next.setEnabled(self._index < n - 1)

    def _go_back(self) -> None:
        if self._index > 0:
            self._index -= 1
            self._refresh()

    def _go_next(self) -> None:
        if self._index < len(self._steps) - 1:
            self._index += 1
            self._refresh()


class TutorialButton(QPushButton):
    """Small button that opens a :class:`TutorialDialog` on click."""

    def __init__(
        self,
        label: str,
        title: str,
        steps: Sequence[TutorialStep],
        parent=None,
    ):
        super().__init__(label, parent)
        self._title = title
        self._steps = list(steps)
        self.setCursor(Qt.PointingHandCursor)
        self.setStyleSheet(
            "QPushButton {"
            "  color:#1565C0;"
            "  border:1px solid #90CAF9;"
            "  border-radius:3px;"
            "  padding:2px 6px;"
            "}"
            "QPushButton:hover { background:#E3F2FD; }"
        )
        self.clicked.connect(self._show)

    def _show(self) -> None:
        dlg = TutorialDialog(self._title, self._steps, parent=self)
        dlg.exec_()


__all__ = ["TutorialStep", "TutorialDialog", "TutorialButton"]
