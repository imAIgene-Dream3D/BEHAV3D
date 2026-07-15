"""
Per-parameter-group physical (µm) ⇄ pixel unit toggle for the BEHAV3D
napari plugin.

The backend stores segmentation/tracking size parameters in their *native*
units (pixels / voxels for segmentation; µm for tracking, whose linking runs
on micron-scaled coordinates). This module lets a parameter group *display*
those values in physical units without changing what is persisted: the
spinboxes show converted values, and :meth:`UnitGroupManager.get_native`
always returns the canonical native value for saving/running.

Conversions use the metadata pixel sizes (``pixel_distance_xy`` /
``pixel_distance_z``); see :func:`behav3d.core.utils.unit_conversion_factor`.
"""
from __future__ import annotations

from qtpy.QtCore import (
    Property,
    QEasingCurve,
    QPropertyAnimation,
    Qt,
)
from qtpy.QtGui import QColor, QPainter
from qtpy.QtWidgets import (
    QAbstractButton,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSpinBox,
    QWidget,
)

from behav3d.core.utils import (
    native_to_physical,
    physical_to_native,
    resolution_from_metadata,
)
from behav3d.core.qt_help import HelpButton


# ----------------------------------------------------------------------
# Custom toggle switch (moved here so both the units toggle and the
# visualization tab can share it).
# ----------------------------------------------------------------------
class QSwitch(QAbstractButton):
    """Modern left↔right toggle switch widget."""

    def __init__(self, parent=None, track_radius=10, thumb_radius=8):
        super().__init__(parent)
        self.setCheckable(True)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        self._track_radius = track_radius
        self._thumb_radius = thumb_radius

        self._margin = 3
        self._base_width = 40
        self._base_height = 22

        self._offset = 0.0
        self._animation = QPropertyAnimation(self, b"offset", self)
        self._animation.setDuration(150)
        self._animation.setEasingCurve(QEasingCurve.InOutSine)

        self.setFixedSize(self._base_width, self._base_height)

    @Property(float)
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, pos):
        self._offset = pos
        self.update()

    def sizeHint(self):
        return self.size()

    def nextCheckState(self):
        super().nextCheckState()
        start = self._offset
        end = 1.0 if self.isChecked() else 0.0
        self._animation.stop()
        self._animation.setStartValue(start)
        self._animation.setEndValue(end)
        self._animation.start()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        p.setPen(Qt.NoPen)
        if self.isChecked():
            p.setBrush(QColor("#4CAF50"))
        else:
            p.setBrush(QColor("#CCCCCC"))
        p.drawRoundedRect(0, 0, self.width(), self.height(), self._track_radius, self._track_radius)
        p.setBrush(QColor("white"))
        thumb_x = self._margin + self._offset * (self.width() - 2 * self._thumb_radius - 2 * self._margin)
        p.drawEllipse(int(thumb_x), self._margin, 2 * self._thumb_radius, 2 * self._thumb_radius)

    def setChecked(self, checked: bool):
        super().setChecked(checked)
        self._offset = 1.0 if checked else 0.0
        self.update()


# Suffixes shown on the spinboxes for each parameter kind + unit.
_SUFFIX = {
    ("distance", True): " µm",
    ("distance", False): " px",
    ("volume", True): " µm³",
    ("volume", False): " vox",
}


class UnitGroupManager:
    """Manage the physical/pixel display of a group of spinboxes.

    Persistence stays native: :meth:`get_native` returns the canonical
    pixel/voxel value for every registered widget regardless of the current
    display unit. The per-widget native value is only recomputed when the
    user edits a spinbox, so repeated toggling never accumulates rounding
    drift.
    """

    def __init__(self, metadata=None, default_physical=True, xy_um=None, z_um=None):
        """``metadata`` is a DataFrame with ``pixel_distance_xy``/``_z``
        columns (the common case). Callers that only have a resolved
        ``{"xy_um": ..., "z_um": ...}`` dict (e.g. the APOC/ConvPaint
        training widgets, which receive pixel sizes directly rather than
        full metadata) can instead pass ``xy_um``/``z_um`` explicitly —
        these take precedence over ``metadata`` when both are given.
        """
        if xy_um is not None and z_um is not None and xy_um > 0 and z_um > 0:
            self.xy, self.z, self.valid = float(xy_um), float(z_um), True
        else:
            self.xy, self.z, self.valid = resolution_from_metadata(metadata)
        # Physical display only makes sense with a valid resolution.
        self.physical = bool(default_physical and self.valid)
        self._entries = []  # list of spinbox widgets (tagged with attrs)

        self.switch = QSwitch()
        self.switch.setChecked(self.physical)
        if not self.valid:
            self.switch.setEnabled(False)
        self.switch.toggled.connect(self._on_toggle)

    # -- widget construction -------------------------------------------
    def header_row(self, label="Units:"):
        """Return a small ``QWidget`` row: ``label  px [switch] µm``."""
        row = QWidget()
        lay = QHBoxLayout(row)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        if label:
            lay.addWidget(QLabel(label))
        lay.addWidget(QLabel("px"))
        lay.addWidget(self.switch)
        lay.addWidget(QLabel("µm"))
        if not self.valid:
            note = QLabel("(no resolution in metadata)")
            note.setStyleSheet("color:#999; font-style:italic; font-size:10px;")
            lay.addWidget(note)
        lay.addWidget(HelpButton(
            "Physical vs. Pixel Units",
            "Toggle how the value(s) above are displayed:\n\n"
            "- 'px' (off): the value is shown in native pixels/voxels — "
            "the same units used for segmentation parameters.\n"
            "- 'µm' (on): the value is shown converted to physical units "
            "(µm for distances, µm³ for volumes) using this sample's pixel "
            "size ('pixel_distance_xy' / 'pixel_distance_z' from metadata).\n\n"
            "This only changes how the number is displayed — the value that "
            "gets saved/used by the pipeline is always converted back to its "
            "native unit automatically. If the toggle is disabled, the "
            "metadata has no valid pixel size and only native units are "
            "available."
        ))
        lay.addStretch(1)
        return row

    # -- registration ---------------------------------------------------
    def register(self, widget, kind, native_value, native_unit="pixel"):
        """Register ``widget`` (a spinbox) holding a value of ``kind``.

        ``kind`` is ``"distance"`` or ``"volume"``. ``native_value`` is the
        value in the parameter's *native* unit, which is what gets persisted
        and passed to the backend. ``native_unit`` names that native unit:

        - ``"pixel"`` (default): the backend uses pixels/voxels (segmentation).
        - ``"physical"``: the backend already uses µm (tracking, whose linking
          runs on micron-scaled coordinates).
        """
        widget._unit_kind = kind
        widget._unit_native_unit = native_unit
        widget._unit_native = float(native_value)
        self._set_display(widget)
        self._update_suffix(widget)
        widget.valueChanged.connect(lambda _v, w=widget: self._on_edit(w))
        self._entries.append(widget)

    def get_native(self, widget):
        """Return the canonical native value for ``widget`` (its native unit)."""
        return float(getattr(widget, "_unit_native", widget.value()))

    def set_native(self, widget, native_value):
        """Set ``widget``'s canonical native value and refresh its display."""
        widget._unit_native = float(native_value)
        self._set_display(widget)

    # -- internal -------------------------------------------------------
    def _display_to_native(self, value, kind, native_unit):
        if not self.valid:
            return float(value)
        if native_unit == "physical":
            # native is µm: physical display == native; pixel display *f -> µm
            if self.physical:
                return float(value)
            return native_to_physical(value, kind, self.xy, self.z)
        # native is pixels/voxels
        if self.physical:
            return physical_to_native(value, kind, self.xy, self.z)
        return float(value)

    def _native_to_display(self, native, kind, native_unit):
        if not self.valid:
            return float(native)
        if native_unit == "physical":
            if self.physical:
                return float(native)
            return physical_to_native(native, kind, self.xy, self.z)
        if self.physical:
            return native_to_physical(native, kind, self.xy, self.z)
        return float(native)

    def _on_edit(self, widget):
        widget._unit_native = self._display_to_native(
            widget.value(), widget._unit_kind, widget._unit_native_unit
        )

    def _set_display(self, widget):
        disp = self._native_to_display(
            widget._unit_native, widget._unit_kind, widget._unit_native_unit
        )
        widget.blockSignals(True)
        if isinstance(widget, QSpinBox):
            widget.setValue(int(round(disp)))
        elif isinstance(widget, QDoubleSpinBox):
            widget.setValue(float(disp))
        widget.blockSignals(False)

    def _update_suffix(self, widget):
        try:
            widget.setSuffix(_SUFFIX.get((widget._unit_kind, self.physical), ""))
        except Exception:
            pass

    def _on_toggle(self, checked):
        self.physical = bool(checked) and self.valid
        for w in self._entries:
            self._set_display(w)
            self._update_suffix(w)
