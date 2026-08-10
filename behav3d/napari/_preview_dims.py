"""Shared helpers for tabs that keep a live, viewer-time-slider-driven
preview layer on the single shared ``napari.Viewer``.

Several tabs (State Classification, Track Classification, Feature
Backprojection, Feature Extraction's Dead/Alive + organoid previews) connect
a callback to ``viewer.dims.events.current_step`` so their overlay layer
recomputes as the user scrubs the time slider. Because every tab shares the
*same* ``Viewer`` instance (see ``behav3d.napari._widget.BEHAV3DWidget``), a
listener left connected by one tab can fire while another tab is bulk-
mutating ``viewer.layers`` (e.g. ``viewer.layers.clear()``), reacting by
re-adding a layer reentrantly. That out-of-band mutation, happening mid-
iteration, desyncs napari's own vispy canvas bookkeeping and raises a
``KeyError`` deep inside ``_vispy/canvas.py::_reorder_layers``.

This module centralizes the "only known-live listeners should ever be
connected, and any of them can be dropped in one call" bookkeeping that
``_feature_extraction.py`` originally implemented for its own two panels
(see its former ``_ACTIVE_PREVIEW_DIMS``), generalized to cover every
listener owner across every tab.

Usage, per preview owner:

    from behav3d.napari._preview_dims import (
        register_preview_dims_listener,
        unregister_preview_dims_listener,
    )

    def _connect_dims_listener(self):
        ...
        self.viewer.dims.events.current_step.connect(_on_step)
        register_preview_dims_listener(self.viewer, self, _on_step)

    def _disconnect_dims_listener(self):
        ...
        self.viewer.dims.events.current_step.disconnect(self._dims_callback)
        unregister_preview_dims_listener(self.viewer, self)

Usage, before any bulk layer mutation (``viewer.layers.clear()`` or
equivalent):

    from behav3d.napari._preview_dims import disconnect_all_preview_dims_listeners

    disconnect_all_preview_dims_listeners(self.viewer)
    self.viewer.layers.clear()
"""
from __future__ import annotations

# Maps id(viewer) -> {owner_key: callback}. ``owner_key`` is typically the
# owning tab/panel instance (``self``), so a given owner can register at
# most one listener per viewer at a time, while distinct owners' listeners
# are tracked independently.
_ACTIVE_PREVIEW_DIMS: dict = {}


def register_preview_dims_listener(viewer, owner_key, callback) -> None:
    """Record ``callback`` as ``owner_key``'s live preview dims listener on
    ``viewer``. Does not connect it — the caller is expected to have already
    called ``viewer.dims.events.current_step.connect(callback)``."""
    if viewer is None:
        return
    _ACTIVE_PREVIEW_DIMS.setdefault(id(viewer), {})[owner_key] = callback


def unregister_preview_dims_listener(viewer, owner_key) -> None:
    """Drop the bookkeeping entry for ``owner_key`` on ``viewer``. Does not
    disconnect anything — the caller handles that."""
    if viewer is None:
        return
    owners = _ACTIVE_PREVIEW_DIMS.get(id(viewer))
    if owners is None:
        return
    owners.pop(owner_key, None)
    if not owners:
        _ACTIVE_PREVIEW_DIMS.pop(id(viewer), None)


def disconnect_all_preview_dims_listeners(viewer) -> None:
    """Disconnect and forget every preview dims listener currently
    registered on ``viewer``, regardless of which tab/panel owns it.

    Call this immediately before any ``viewer.layers.clear()`` (or other
    bulk layer replacement) so no stale listener can fire reentrantly while
    the layer list is mid-mutation.
    """
    if viewer is None:
        return
    owners = _ACTIVE_PREVIEW_DIMS.pop(id(viewer), None)
    if not owners:
        return
    for callback in owners.values():
        try:
            viewer.dims.events.current_step.disconnect(callback)
        except Exception:
            pass


def stop_dim_playback(viewer) -> None:
    """Best-effort stop of any running napari dims animation.

    Clearing/replacing layers while a 3-D movie is playing destroys the
    ``QtDimSliderWidget`` objects the animation thread still references,
    raising ``RuntimeError: wrapped C/C++ object … has been deleted``. We
    stop playback via the public ``QtDims.stop()`` and also quit the
    animation thread directly, since depending on the napari version either
    (or neither) attribute may be present while a movie is live.
    """
    try:
        qt_dims = viewer.window._qt_viewer.dims
    except Exception:
        return
    # Public stop first (handles the current napari playback implementation).
    try:
        qt_dims.stop()
    except Exception:
        pass
    # Fall back to tearing down the animation thread if it is still around.
    try:
        thread = getattr(qt_dims, "_animation_thread", None)
        if thread is not None:
            thread.quit()
            thread.wait()
    except Exception:
        pass
