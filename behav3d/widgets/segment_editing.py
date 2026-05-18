"""Notebook entry point for the manual tracked-segment correction widget.

Provides ``TrackedSegmentEditingPanel`` — an ipywidgets panel modelled on
``TrackingVisualizationPanel`` (see :mod:`behav3d.widgets.tracking`).  Picking
a sample + cell-type and clicking *Open editor in napari* spins up a fresh
``napari.Viewer`` that loads:

* raw image channels,
* the tracked-segments labels for the chosen cell type,
* the tracks layer for that cell type.

Untracked segment layers are deliberately not loaded — the editor only acts
on tracked segments.

The panel docks the same ``TrackedSegmentEditor`` Qt widget used by the
napari plugin's Visualization tab, so save / undo / discard / exit-prompt
behaviour is identical across both surfaces.
"""
from __future__ import annotations

import traceback
from pathlib import Path

import ipywidgets as widgets
import pandas as pd

from behav3d.napari._loaders import (
    cell_types_with_tracked_segments,
    load_raw_into_viewer,
    load_tracked_segments_into_viewer,
    load_tracks_into_viewer,
    tracked_segments_paths_for,
)


class TrackedSegmentEditingPanel:
    """Independent ipywidgets panel that launches the manual editor.

    Parameters
    ----------
    metadata_loader:
        Object exposing ``.metadata`` (a ``pd.DataFrame``) and an
        ``output_dir`` attribute.  Same contract as the other notebook
        panels (e.g. ``MetadataLoader``).
    """

    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self._viewer = None
        self._editor = None
        self._dock = None

        self._header = widgets.HTML(
            "<h4 style='margin:4px 0;color:#1976D2;'>"
            "✏️ Manual correction of tracked segments</h4>"
        )
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Re-read the loaded metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        self.sample_dropdown = widgets.Dropdown(
            options=[], value=None, description="Sample:",
            layout=widgets.Layout(width="350px"), disabled=True,
        )
        self.sample_dropdown.observe(self._on_sample_changed, names="value")
        self.celltype_dropdown = widgets.Dropdown(
            options=[], value=None, description="Cell type:",
            layout=widgets.Layout(width="350px"), disabled=True,
        )

        self.open_button = widgets.Button(
            description="Open editor in napari",
            button_style="primary",
            tooltip="Launch napari with the chosen sample loaded for editing",
            icon="edit",
            layout=widgets.Layout(width="260px"),
            disabled=True,
        )
        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active napari viewer (with a Save/Discard prompt if needed)",
            layout=widgets.Layout(width="180px", display="none"),
        )
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)

        self.msg = widgets.Output()
        self._panel = widgets.VBox([
            self._header,
            widgets.HBox([self._status, self._refresh_btn]),
            self.sample_dropdown,
            self.celltype_dropdown,
            widgets.HBox([self.open_button, self.close_button]),
            self.msg,
        ])
        self._panel.add_class("behav3d-segment-editing-panel")
        self._maybe_build_from_loader()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def display(self):
        widgets.display(self._panel)

    # ------------------------------------------------------------------
    # State helpers
    # ------------------------------------------------------------------
    def _get_metadata(self) -> pd.DataFrame | None:
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty:
            return df
        return None

    def _selected_row(self) -> pd.Series | None:
        df = self._get_metadata()
        if df is None or self.sample_dropdown.value is None:
            return None
        rows = df[df["sample_name"].astype(str) == str(self.sample_dropdown.value)]
        if rows.empty:
            return None
        return rows.iloc[0]

    def _maybe_build_from_loader(self):
        df = self._get_metadata()
        if df is None or "sample_name" not in df.columns:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            self.sample_dropdown.disabled = True
            self.celltype_dropdown.disabled = True
            self.open_button.disabled = True
            return
        sample_names = sorted(df["sample_name"].astype(str).unique().tolist())
        self.sample_dropdown.options = sample_names
        self.sample_dropdown.value = sample_names[0] if sample_names else None
        self.sample_dropdown.disabled = False
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        self._refresh_celltypes()

    def _refresh_celltypes(self):
        row = self._selected_row()
        opts: list[str] = []
        if row is not None:
            opts = cell_types_with_tracked_segments(row)
        self.celltype_dropdown.options = opts
        if opts:
            self.celltype_dropdown.value = opts[0]
            self.celltype_dropdown.disabled = False
            self.open_button.disabled = False
        else:
            self.celltype_dropdown.value = None
            self.celltype_dropdown.disabled = True
            self.open_button.disabled = True
            with self.msg:
                self.msg.clear_output()
                print(
                    "  ⚠️ No tracked-segments files found for this sample.  "
                    "Run tracking first."
                )

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------
    def _on_refresh_clicked(self, _):
        with self.msg:
            self.msg.clear_output()
            self._maybe_build_from_loader()

    def _on_sample_changed(self, _change):
        self._refresh_celltypes()

    def _on_open_clicked(self, _):
        with self.msg:
            self.msg.clear_output()
            try:
                self._open_editor()
            except Exception as exc:
                traceback.print_exc()
                print(f"❌ {exc}")

    def _on_close_clicked(self, _):
        with self.msg:
            try:
                self._close_editor(prompt=True)
            except Exception as exc:
                traceback.print_exc()
                print(f"❌ {exc}")

    # ------------------------------------------------------------------
    # Editor lifecycle
    # ------------------------------------------------------------------
    def _open_editor(self):
        import napari  # local import: keeps notebook startup fast

        row = self._selected_row()
        if row is None:
            print("⚠️ No sample selected.")
            return
        ct_name = self.celltype_dropdown.value
        if not ct_name:
            print("⚠️ No cell type with tracked segments available.")
            return
        paths = tracked_segments_paths_for(row, ct_name)
        if paths is None:
            print(f"⚠️ No tracked segments path for {ct_name}.")
            return

        sample_name = str(row["sample_name"])
        output_dir = str(getattr(self.metadata_loader, "output_dir", "") or "")

        if self._viewer is not None:
            print("Closing previous viewer first…")
            self._close_editor(prompt=True)
            if self._viewer is not None:
                # User cancelled the close prompt.
                return

        viewer = napari.Viewer(title=f"BEHAV3D – Edit tracked segments ({sample_name})")
        self._viewer = viewer

        def _log(msg: str):
            print(msg)

        load_raw_into_viewer(viewer, sample_name, row, output_dir=output_dir, log=_log)
        added = load_tracked_segments_into_viewer(
            viewer, sample_name, row, log=_log,
            only_cell_type=ct_name, as_numpy=False,
        )
        load_tracks_into_viewer(
            viewer, sample_name, row, visible=True, log=_log,
            only_cell_type=ct_name,
        )

        # Pick the layer matching the chosen cell type.
        target = next(
            (n for n in added if n.endswith(f" – {ct_name} tracked segments")),
            None,
        )
        if target is None:
            print(f"❌ Could not find tracked-segments layer for {ct_name}.")
            viewer.close()
            self._viewer = None
            return

        from behav3d.napari._segment_editor import TrackedSegmentEditor

        editor = TrackedSegmentEditor(
            viewer=viewer,
            metadata_loader=self.metadata_loader,
            layer_name=target,
        )
        # When the user clicks "Stop editing" inside the editor, also
        # close this notebook-spawned viewer.
        editor.stop_callback = lambda: self._close_editor(prompt=False)
        self._editor = editor
        self._dock = viewer.window.add_dock_widget(
            editor, name="Manual edition", area="right"
        )

        # Intercept the napari window close so we can run the unsaved-edit
        # prompt before the viewer disappears.
        self._close_filter = _NapariCloseFilter(self)
        try:
            qt_window = viewer.window.qt_viewer.window()
            if qt_window is not None:
                qt_window.installEventFilter(self._close_filter)
                self._close_filter.watched = qt_window
        except Exception:
            pass

        self.close_button.layout.display = "inline-block"
        print(f"✅ Editor ready for '{sample_name} – {ct_name}'.")

    def _close_editor(self, prompt: bool = True) -> bool:
        if self._editor is not None:
            if prompt and not self._editor.request_exit():
                return False
            try:
                if not prompt:
                    self._editor._cleanup()
            except Exception:
                pass
            self._editor = None
        if self._dock is not None and self._viewer is not None:
            try:
                self._viewer.window.remove_dock_widget(self._dock)
            except Exception:
                pass
            self._dock = None
        if self._viewer is not None:
            try:
                self._viewer.close()
            except Exception:
                pass
            self._viewer = None
        self.close_button.layout.display = "none"
        return True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
from qtpy.QtCore import QEvent, QObject


class _NapariCloseFilter(QObject):
    """Qt event filter that runs the editor's exit prompt when the user
    closes the napari window via its title-bar X.
    """

    def __init__(self, panel: TrackedSegmentEditingPanel):
        super().__init__()
        self.panel = panel
        self.watched = None

    def eventFilter(self, obj, event):  # type: ignore[override]
        if event.type() == QEvent.Close and self.panel._editor is not None:
            if not self.panel._editor.request_exit():
                event.ignore()
                return True
            # Editor accepted the close — let the viewer close cleanly.
            self.panel._editor = None
            self.panel._dock = None
            self.panel._viewer = None
            self.panel.close_button.layout.display = "none"
        return False
