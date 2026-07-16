"""Runtime workarounds for napari bugs that affect BEHAV3D viewers.

Patches are applied to napari classes, so they take effect for every viewer in the
process, including layers that already exist. Installed from this package's
``__init__``.
"""
import functools
import logging

logger = logging.getLogger(__name__)

_GUARD_FLAG = "_behav3d_thumbnail_guard"


def _patch_labels_thumbnail():
    """Skip Labels thumbnail updates whose cached slice has a stale ``ndisplay``.

    napari#7918. ``Layer._slice_dims`` always updates ``_slice_input`` (which
    carries ``ndisplay``), but ``_refresh_sync`` returns early for invisible
    layers, so a hidden layer keeps the slice it was last computed at. In 3D,
    ``Labels._update_thumbnail`` trusts ``_slice_input.ndisplay`` and max-projects
    that stale slice down to 1-D, which ``scipy.ndimage.zoom`` then rejects with
    "sequence argument must have length equal to input rank". Any opacity write
    reaches it, and napari's layer controls set opacity whenever they are built,
    so selecting a hidden Labels layer in 3D is enough to raise.

    Skipping matches napari's own policy of not slicing hidden layers: ``refresh``
    recomputes the thumbnail as soon as the layer is shown again.

    napari's own fix (PR #8251, 0.6.5) only guards ``_slice.empty``, so it still
    raises for a layer that was sliced while visible and then hidden. Comparing
    ``ndisplay`` covers that case too, and stays correct once we move off 0.5.6.
    """
    from napari.layers import Labels

    original = Labels._update_thumbnail
    if getattr(original, _GUARD_FLAG, False):
        return

    @functools.wraps(original)
    def _update_thumbnail(self):
        cached = getattr(getattr(self, "_slice", None), "slice_input", None)
        current = getattr(self, "_slice_input", None)
        if cached is not None and current is not None:
            if cached.ndisplay != current.ndisplay:
                return
        original(self)

    setattr(_update_thumbnail, _GUARD_FLAG, True)
    Labels._update_thumbnail = _update_thumbnail


def install_napari_patches():
    """Apply all napari workarounds. Safe to call more than once."""
    try:
        _patch_labels_thumbnail()
    except Exception:
        # Never let a workaround break importing the plugin: if napari is absent
        # (headless behav3d use) or its internals moved, the bug just resurfaces.
        logger.debug("Could not install napari Labels thumbnail guard", exc_info=True)
