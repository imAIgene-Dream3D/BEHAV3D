"""Reusable widgets for the BEHAV3D napari plugin.

The ``HelpButton`` and ``make_help_row`` helpers live in
``behav3d.core.qt_help`` so they can be reused by the standalone Qt training
widgets in ``behav3d.preprocessing.segmentation`` without depending on napari.
They are re-exported here for backwards compatibility with the rest of the
plugin code that already imports them from this module.
"""

from pathlib import Path

from qtpy.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QLabel,
    QMessageBox,
    QVBoxLayout,
)

from behav3d.core.qt_help import (
    HelpButton,
    HelpSection,
    IllustratedHelpButton,
    make_help_row,
)

__all__ = [
    "HelpButton",
    "HelpSection",
    "IllustratedHelpButton",
    "make_help_row",
    "LABEL_DIM_ORDER_OPTIONS",
    "browse_file_or_zarr",
    "prompt_axis_order",
    "resolve_external_path",
]


# Axis-order permutations for *label* images (segmentation/tracking/dead-mask
# masks) — unlike raw intensity images these never require a channel axis, so
# this list (mirroring behav3d/widgets/utils.py's DIM_ORDER_OPTIONS) includes
# both with- and without-C variants. Kept separate from
# behav3d.napari._data_preparation.DIM_ORDER_OPTIONS, which is 5D-only and
# assumes a channel axis for raw images.
LABEL_DIM_ORDER_OPTIONS = [
    # 5D
    "TCZYX", "TZCYX", "CTZYX", "CZTYX", "ZCTYX", "ZTCYX",
    # 4D (missing one of T / C / Z)
    "CZYX", "ZCYX", "TZYX", "ZTYX", "TCYX", "CTYX",
    # 3D / 2D
    "ZYX", "CYX", "TYX",
    "YX",
]


def browse_file_or_zarr(parent, title: str, filter_str: str = "All Files (*.*)", allow_zarr: bool = True) -> str:
    """Open a file/directory browser and return the selected path (or "").

    When ``allow_zarr`` is True, first asks the user whether they want to
    pick a single file or a Zarr directory/folder, then opens the matching
    ``QFileDialog``. Mirrors
    ``behav3d.napari._data_preparation.DataPreparationTab._browse_path``.
    """
    if allow_zarr:
        msg = QMessageBox(parent)
        msg.setWindowTitle("Select Input Type")
        msg.setText("Would you like to select a file or a Zarr directory/folder?")

        btn_file = msg.addButton("File", QMessageBox.AcceptRole)
        btn_dir = msg.addButton("Zarr / Folder", QMessageBox.AcceptRole)
        msg.addButton("Cancel", QMessageBox.RejectRole)

        msg.exec_()

        clicked = msg.clickedButton()
        if clicked == btn_file:
            path, _ = QFileDialog.getOpenFileName(parent, title, "", filter_str)
            return path
        elif clicked == btn_dir:
            return QFileDialog.getExistingDirectory(parent, title, "")
        return ""
    else:
        path, _ = QFileDialog.getOpenFileName(parent, title, "", filter_str)
        return path


class _AxisOrderDialog(QDialog):
    """Modal dialog: shows a probed shape/detected axis order and lets the
    user confirm or pick the correct one from a filtered dropdown."""

    def __init__(self, parent, shape, detected_order, options):
        super().__init__(parent)
        self.setWindowTitle("Confirm Axis Order")

        layout = QVBoxLayout(self)

        shape_str = " x ".join(str(s) for s in shape)
        det_str = detected_order or "not detected"
        layout.addWidget(QLabel(
            f"<b>Shape:</b> {shape_str} ({len(shape)}D)<br>"
            f"<b>Detected order:</b> {det_str}"
        ))
        layout.addWidget(QLabel(
            "Please confirm the axis order of this label image "
            "(e.g. TZYX = time, z, y, x):"
        ))

        self.combo = QComboBox()
        self.combo.addItems(options)
        if detected_order and detected_order in options:
            self.combo.setCurrentText(detected_order)
        layout.addWidget(self.combo)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def chosen_axis_order(self) -> str:
        return self.combo.currentText()


def prompt_axis_order(parent, path, *, log=print):
    """Probe ``path``'s shape/axis order and let the user confirm/pick it.

    Returns the chosen axis-order string, or ``None`` if the user cancels or
    probing the source file fails.
    """
    from behav3d.io.images import get_image_dimension_order, get_image_shape

    path = Path(path)
    try:
        shape = get_image_shape(path)
    except Exception as exc:
        log(f"  Could not read shape for {path.name}: {exc}")
        QMessageBox.warning(parent, "Cannot Read File", f"Could not read shape for {path.name}:\n{exc}")
        return None

    try:
        detected_order = get_image_dimension_order(path)
    except Exception:
        detected_order = None

    options = [o for o in LABEL_DIM_ORDER_OPTIONS if len(o) == len(shape)] or list(LABEL_DIM_ORDER_OPTIONS)

    dialog = _AxisOrderDialog(parent, shape, detected_order, options)
    if dialog.exec_() != QDialog.Accepted:
        return None
    return dialog.chosen_axis_order()


def resolve_external_path(path_str, metadata_csv_path):
    """Resolve a metadata path string, trying it as-is then relative to the
    metadata CSV's parent directory. Returns a ``Path`` (which may not
    exist — callers should check ``.exists()``), or ``None`` if
    ``path_str`` is empty."""
    if not path_str:
        return None
    p = Path(path_str)
    if p.exists():
        return p
    if metadata_csv_path:
        p_rel = Path(metadata_csv_path).parent / path_str
        if p_rel.exists():
            return p_rel
    return p
