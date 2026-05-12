"""Reusable widgets for the BEHAV3D napari plugin.

The ``HelpButton`` and ``make_help_row`` helpers live in
``behav3d.core.qt_help`` so they can be reused by the standalone Qt training
widgets in ``behav3d.preprocessing.segmentation`` without depending on napari.
They are re-exported here for backwards compatibility with the rest of the
plugin code that already imports them from this module.
"""

from behav3d.core.qt_help import HelpButton, make_help_row

__all__ = ["HelpButton", "make_help_row"]
