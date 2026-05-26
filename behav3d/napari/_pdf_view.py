"""
BEHAV3D napari plugin – PDF / image preview helpers.

Renders BEHAV3D's result PDFs into RGBA page stacks via ``pypdfium2`` and
opens them as an Image layer in a dedicated secondary ``napari.Viewer``.

The secondary viewer is reused across calls (by default) so users can
quickly stack multiple result PDFs in a single window and compare them
via the layer list. Each ``Open in napari`` action adds a fresh layer
(name-collisions are auto-suffixed by napari); previously opened layers
are never silently removed.

Also provides a couple of cross-platform helpers used by the Results
panel:

- :func:`open_externally` — open a file with the OS default app.
- :func:`reveal_in_folder` — show the file in Explorer / Finder / file
  manager, ideally with the file selected.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from qtpy.QtCore import QUrl
from qtpy.QtGui import QDesktopServices
from qtpy.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QWidget,
)


# ---------------------------------------------------------------------------
# PDF rasterization
# ---------------------------------------------------------------------------
def is_pdf_viewable(path: Path) -> bool:
    """Return True when ``path`` is an existing ``.pdf`` file."""
    return path.suffix.lower() == ".pdf" and path.is_file()


def pdf_to_stack(
    path: Path,
    dpi: int = 150,
    max_pages: Optional[int] = None,
) -> np.ndarray:
    """Rasterize ``path`` to a ``(N, H, W, 4)`` ``uint8`` RGBA stack.

    Pages with different aspect ratios are padded to the largest page's
    bounding box (top-left aligned) so the output is a single rectangular
    array. Empty padding is opaque white.

    Parameters
    ----------
    path : Path
        Path to the PDF file.
    dpi : int, default 150
        Rendering resolution in dots per inch. PDF native is 72 dpi, so
        ``scale = dpi / 72``.
    max_pages : int, optional
        Hard cap on the number of pages rendered (useful for huge PDFs).
    """
    import pypdfium2 as pdfium  # local import keeps module load fast

    pdf = pdfium.PdfDocument(str(path))
    scale = dpi / 72.0
    pages: list[np.ndarray] = []
    try:
        for i in range(len(pdf)):
            if max_pages is not None and i >= max_pages:
                break
            page = pdf[i]
            try:
                img = page.render(scale=scale).to_pil().convert("RGBA")
                pages.append(np.asarray(img))
            finally:
                page.close()
    finally:
        pdf.close()

    if not pages:
        # Return a single-pixel white image so callers don't crash on
        # empty / corrupt PDFs.
        return np.full((1, 1, 1, 4), 255, dtype=np.uint8)

    h = max(p.shape[0] for p in pages)
    w = max(p.shape[1] for p in pages)
    out = np.full((len(pages), h, w, 4), 255, dtype=np.uint8)
    for i, p in enumerate(pages):
        out[i, : p.shape[0], : p.shape[1]] = p
    return out


# ---------------------------------------------------------------------------
# Secondary napari viewer (cached, reusable)
# ---------------------------------------------------------------------------
_RESULT_VIEWER = None  # type: ignore[assignment]


def _invalidate_result_viewer(*_args) -> None:
    """Clear the cached secondary viewer reference."""
    global _RESULT_VIEWER
    _RESULT_VIEWER = None


def _is_viewer_alive(viewer) -> bool:
    """Return True only if ``viewer``'s underlying Qt window is usable.

    Just accessing ``viewer.window`` is not enough: napari keeps the
    Python wrapper alive after the user closes the second window, so
    later ``viewer.add_image`` calls fail deep inside Qt with a
    ``RuntimeError`` whose ``str`` is just the layer's repr ("Image
    layer 'foo' at 0x…"). We probe an actual Qt method to detect the
    destroyed C++ object, and also treat hidden windows as dead.
    """
    try:
        qt_window = getattr(viewer.window, "_qt_window", None)
        if qt_window is None:
            return False
        if not qt_window.isVisible():
            return False
        return True
    except (RuntimeError, AttributeError):
        return False


def _result_viewer(reuse: bool = True):
    """Return the cached secondary viewer, creating one if needed.

    The viewer is lazily imported so this module doesn't pull napari on
    import (matches the rest of the plugin's lazy-import style).
    """
    global _RESULT_VIEWER
    import napari  # noqa: WPS433 – lazy import on purpose

    if reuse and _RESULT_VIEWER is not None:
        if _is_viewer_alive(_RESULT_VIEWER):
            return _RESULT_VIEWER
        # The viewer's Qt window was closed/destroyed. Best-effort
        # tear-down of the stale Python wrapper so it doesn't leak.
        try:
            _RESULT_VIEWER.close()
        except Exception:
            pass
        _RESULT_VIEWER = None

    _RESULT_VIEWER = napari.Viewer(title="BEHAV3D – Results", show=True)

    # Drop the cache as soon as the window is destroyed/closed so the
    # next call spawns a fresh viewer instead of trying to use a dead
    # Qt widget.
    try:
        qt_window = getattr(_RESULT_VIEWER.window, "_qt_window", None)
        if qt_window is not None:
            qt_window.destroyed.connect(_invalidate_result_viewer)
    except Exception:
        pass

    return _RESULT_VIEWER


def _install_result_dock(viewer, path: Path) -> None:
    """Attach a tiny dock widget with Open / Reveal buttons.

    Re-added every time a new file is shown so the buttons always point
    at the file currently visible.
    """
    # Remove the previous instance if any (we tag the widget so we can
    # find it without storing a reference on the napari viewer).
    try:
        for dock in list(viewer.window._dock_widgets.values()):  # type: ignore[attr-defined]
            inner = dock.widget()
            if getattr(inner, "_behav3d_results_dock", False):
                viewer.window.remove_dock_widget(dock)
    except Exception:
        # If napari changes its private API we just live with stacked
        # docks; not a fatal issue.
        pass

    bar = QWidget()
    bar._behav3d_results_dock = True  # type: ignore[attr-defined]
    lay = QHBoxLayout(bar)
    lay.setContentsMargins(6, 4, 6, 4)
    btn_open = QPushButton("Open externally")
    btn_reveal = QPushButton("Reveal in folder")
    btn_open.clicked.connect(lambda: open_externally(path))
    btn_reveal.clicked.connect(lambda: reveal_in_folder(path))
    lay.addWidget(btn_open)
    lay.addWidget(btn_reveal)
    lay.addStretch(1)
    viewer.window.add_dock_widget(bar, area="bottom", name=path.name)


def open_pdf_in_napari(
    path: Path,
    *,
    dpi: int = 150,
    reuse: bool = True,
    max_pages: Optional[int] = None,
):
    """Render ``path`` and add it as an Image layer in the result viewer.

    Returns the (possibly cached) secondary ``napari.Viewer`` so callers
    can keep a reference if needed.
    """
    path = Path(path)
    if not is_pdf_viewable(path):
        raise ValueError(f"Not a viewable PDF: {path}")

    stack = pdf_to_stack(path, dpi=dpi, max_pages=max_pages)
    viewer = _result_viewer(reuse=reuse)

    # Always add a fresh layer — keep previously opened BEHAV3D results
    # so the user can compare files side by side. napari handles name
    # collisions by appending " [1]", " [2]", etc.
    layer = viewer.add_image(
        stack,
        rgb=True,
        name=path.stem,
        contrast_limits=(0, 255),
        metadata={"behav3d_pdf": True, "source": str(path)},
    )
    # ``layer`` is RGB so napari treats the leading axis as a dim slider:
    # label it "page" for clarity.
    try:
        viewer.dims.axis_labels = ("page", "y", "x")
        viewer.dims.set_point(0, 0)  # jump back to first page
    except Exception:
        pass

    _install_result_dock(viewer, path)

    try:
        viewer.window.activate()
    except Exception:
        pass

    return viewer


def open_image_in_napari(path: Path, *, reuse: bool = True):
    """Open a static image file (PNG/JPEG/TIFF) in the result viewer."""
    import imageio.v3 as iio  # lazy

    path = Path(path)
    data = iio.imread(str(path))
    viewer = _result_viewer(reuse=reuse)
    is_rgb = data.ndim == 3 and data.shape[-1] in (3, 4)
    # Always add as a new layer — preserve anything the user already
    # opened so multiple results can be compared at once.
    viewer.add_image(
        data,
        rgb=is_rgb,
        name=path.stem,
        metadata={"behav3d_pdf": True, "source": str(path)},
    )
    _install_result_dock(viewer, path)
    try:
        viewer.window.activate()
    except Exception:
        pass
    return viewer


# ---------------------------------------------------------------------------
# Cross-platform "open externally" / "reveal in folder"
# ---------------------------------------------------------------------------
def open_externally(path: Path) -> bool:
    """Open ``path`` with the OS default application.

    Returns True if a request was issued (no guarantee it succeeded).
    """
    path = Path(path)
    if not path.exists():
        return False
    return bool(QDesktopServices.openUrl(QUrl.fromLocalFile(str(path))))


def reveal_in_folder(path: Path) -> bool:
    """Show ``path`` in the OS file manager, selected if possible."""
    path = Path(path)
    if not path.exists():
        return False
    target = str(path.resolve())
    try:
        if sys.platform.startswith("win"):
            # /select, requires the path as a single argument; spaces ok.
            subprocess.Popen(["explorer", f"/select,{target}"])
            return True
        if sys.platform == "darwin":
            subprocess.Popen(["open", "-R", target])
            return True
        # Linux / other – no portable "select" flag; open the parent dir.
        parent = str(Path(target).parent)
        opener = shutil.which("xdg-open") or shutil.which("gio")
        if opener:
            if opener.endswith("gio"):
                subprocess.Popen([opener, "open", parent])
            else:
                subprocess.Popen([opener, parent])
            return True
    except Exception:
        pass

    # Fallback: Qt's openUrl on the parent directory.
    return bool(
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path.parent)))
    )


# Public API explicitly re-exported for clarity in `from _pdf_view import *`.
__all__ = [
    "is_pdf_viewable",
    "pdf_to_stack",
    "open_pdf_in_napari",
    "open_image_in_napari",
    "open_externally",
    "reveal_in_folder",
]
