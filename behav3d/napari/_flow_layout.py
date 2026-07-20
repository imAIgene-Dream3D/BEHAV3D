"""Responsive flow layout for Qt widgets.

A ``FlowLayout`` arranges child widgets in a left-to-right flow that wraps
to a new row when the available width is exhausted — similar to CSS
``flexbox`` with ``flex-wrap: wrap``.

Adapted from the Qt flow-layout example
(https://doc.qt.io/qt-5/qtwidgets-layouts-flowlayout-example.html).
"""
from __future__ import annotations

from qtpy.QtCore import QPoint, QRect, QSize, Qt
from qtpy.QtWidgets import QLayout, QSizePolicy, QWidget


class FlowLayout(QLayout):
    """Layout that wraps children to a new row when width is exhausted.

    Parameters
    ----------
    parent:
        Parent widget (or ``None`` for a layout that is added to another
        layout).
    margin:
        Content margin applied on all four sides (px).
    h_spacing:
        Horizontal spacing between items (px).  ``-1`` uses the style
        default.
    v_spacing:
        Vertical spacing between rows (px).  ``-1`` uses the style default.
    """

    def __init__(
        self,
        parent: QWidget | None = None,
        margin: int = 0,
        h_spacing: int = -1,
        v_spacing: int = -1,
    ) -> None:
        super().__init__(parent)
        self._h_space = h_spacing
        self._v_space = v_spacing
        self._items: list = []
        self.setContentsMargins(margin, margin, margin, margin)

    # ------------------------------------------------------------------
    # QLayout overrides
    # ------------------------------------------------------------------
    def addItem(self, item) -> None:  # type: ignore[override]
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int):
        if 0 <= index < len(self._items):
            return self._items[index]
        return None

    def takeAt(self, index: int):
        if 0 <= index < len(self._items):
            return self._items.pop(index)
        return None

    def expandingDirections(self):
        return Qt.Orientations()  # type: ignore[attr-defined]

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._do_layout(QRect(0, 0, width, 0), test_only=True)

    def setGeometry(self, rect: QRect) -> None:
        super().setGeometry(rect)
        self._do_layout(rect, test_only=False)

    def sizeHint(self) -> QSize:
        return self.minimumSize()

    def minimumSize(self) -> QSize:
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        size += QSize(m.left() + m.right(), m.top() + m.bottom())
        return size

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _h_spacing_resolved(self) -> int:
        if self._h_space >= 0:
            return self._h_space
        return self._smart_spacing(QStyle_PM_LayoutHorizontalSpacing)

    def _v_spacing_resolved(self) -> int:
        if self._v_space >= 0:
            return self._v_space
        return self._smart_spacing(QStyle_PM_LayoutVerticalSpacing)

    def _smart_spacing(self, pm) -> int:  # noqa: ANN001
        parent = self.parent()
        if parent is None:
            return -1
        if parent.isWidgetType():
            style = parent.style()  # type: ignore[attr-defined]
            return style.pixelMetric(pm, None, parent)
        return self.spacing()

    def _do_layout(self, rect: QRect, *, test_only: bool) -> int:
        m = self.contentsMargins()
        effective = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        x = effective.x()
        y = effective.y()
        line_height = 0

        for item in self._items:
            widget = item.widget()
            h_sp = self._h_spacing_resolved()
            if h_sp == -1:
                h_sp = widget.style().layoutSpacing(  # type: ignore[union-attr]
                    QSizePolicy.PushButton, QSizePolicy.PushButton, Qt.Horizontal
                )
            v_sp = self._v_spacing_resolved()
            if v_sp == -1:
                v_sp = widget.style().layoutSpacing(  # type: ignore[union-attr]
                    QSizePolicy.PushButton, QSizePolicy.PushButton, Qt.Vertical
                )
            item_w = item.sizeHint().width()
            item_h = item.sizeHint().height()
            next_x = x + item_w + h_sp
            if next_x - h_sp > effective.right() and line_height > 0:
                x = effective.x()
                y = y + line_height + v_sp
                next_x = x + item_w + h_sp
                line_height = 0
            if not test_only:
                item.setGeometry(QRect(QPoint(x, y), item.sizeHint()))
            x = next_x
            line_height = max(line_height, item_h)

        return y + line_height - rect.y() + m.bottom()


# ---------------------------------------------------------------------------
# Qt style pixel metrics — resolved lazily to avoid import-time Qt calls
# ---------------------------------------------------------------------------
try:
    from qtpy.QtWidgets import QStyle
    QStyle_PM_LayoutHorizontalSpacing = QStyle.PM_LayoutHorizontalSpacing
    QStyle_PM_LayoutVerticalSpacing = QStyle.PM_LayoutVerticalSpacing
except Exception:
    QStyle_PM_LayoutHorizontalSpacing = 38  # fallback int value
    QStyle_PM_LayoutVerticalSpacing = 39
