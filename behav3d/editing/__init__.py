"""BEHAV3D manual correction backend.

Pure-numpy / pure-skimage primitives for editing 4D tracked-segment volumes
(``T x Z x Y x X``), plus a small in-memory ``EditBuffer`` that wraps a
zarr-backed dask array and exposes an atomic save that overwrites the 
original ``*_tracked.zarr`` and regenerates the matching ``*_tracks.csv``.

The primitives are agnostic of any GUI; the napari Qt editor and the
ipywidgets notebook panel both consume the same API.
"""
from behav3d.editing.tracked_segments import (
    split_label,
    create_label,
    merge_labels,
    fill_label,
    erode_label,
    dilate_label,
    delete_label,
    delete_labels,
    lifetime_of,
)
from behav3d.editing.edit_buffer import EditBuffer

__all__ = [
    "EditBuffer",
    "split_label",
    "create_label",
    "merge_labels",
    "fill_label",
    "erode_label",
    "dilate_label",
    "delete_label",
    "delete_labels",
    "lifetime_of",
]
