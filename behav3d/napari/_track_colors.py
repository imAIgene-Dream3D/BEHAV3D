"""Sync a napari Tracks layer's colors to a matching Labels layer.

BEHAV3D's "tracked segments" images are labeled *by track ID*, and the
tracks CSV's ``TrackID`` column becomes the Tracks layer's built-in
``track_id`` property. That means a "tracked segments" Labels layer and its
matching Tracks layer already share one ID space — we just need to point the
Tracks layer at the same per-ID colors the Labels layer is already showing,
so a cell's track matches its segment (e.g. a yellow cell gets a yellow
track) instead of napari's default continuous track_id→turbo gradient.
"""
from __future__ import annotations


def _install_track_colormap(tracks_layer, prop_key: str, colormap) -> bool:
    """Register ``colormap`` for coloring ``tracks_layer`` by ``prop_key``.

    Works around a napari 0.5.x bug where the ``colormaps_dict`` property
    setter is unreachable: the ``@colormaps_dict.setter``-decorated method in
    napari's own source is itself named ``colomaps_dict`` (typo, missing an
    "r"), so Python binds the writable property to that misspelled name and
    leaves the correctly-spelled ``colormaps_dict`` read-only.

    Tries the correct spelling first (self-heals if a future napari release
    fixes the typo), then the misspelled alias that actually works today,
    then a direct private-attribute write as a last resort. Never raises —
    callers get ``False`` back if none of the tiers work, and the tracks
    layer is simply left with napari's default coloring.
    """
    try:
        merged = {**tracks_layer.colormaps_dict, prop_key: colormap}
    except Exception:
        return False

    for attr in ("colormaps_dict", "colomaps_dict"):
        try:
            setattr(tracks_layer, attr, merged)
            return True
        except AttributeError:
            continue
        except Exception:
            return False

    try:
        tracks_layer._colormaps_dict = merged
        return True
    except Exception:
        return False


def sync_tracks_colors_to_labels(tracks_layer, labels_layer) -> bool:
    """Recolor ``tracks_layer`` (by track_id) using the same per-ID colors
    ``labels_layer`` already shows, so a track matches the cell it belongs to.

    Assumes ``tracks_layer``'s ``track_id`` values live in the same ID space
    as ``labels_layer``'s label values — true for BEHAV3D's tracked-segments
    images, which are labeled by track ID.

    Returns True if the sync was applied, False if it was skipped/failed
    (e.g. either layer is missing, or the tracks layer's colormap API has
    changed in a way this helper doesn't handle) — in which case the tracks
    layer is left with whatever coloring it already had.
    """
    if tracks_layer is None or labels_layer is None:
        return False
    try:
        colormap = labels_layer.colormap
    except Exception:
        return False
    if not _install_track_colormap(tracks_layer, "track_id", colormap):
        return False
    try:
        tracks_layer.color_by = "track_id"  # forces _recolor_tracks() to run
    except Exception:
        return False
    return True
