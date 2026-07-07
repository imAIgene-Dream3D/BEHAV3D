"""
Visualization tab for the BEHAV3D napari plugin.

After metadata is loaded (from the Data Preparation tab), this tab lets the
user choose a dataset/sample and load it into the napari viewer:
  - Raw zarr image → each channel as a separate Image layer (dask-backed)
  - Segment images → Labels layers for each cell type
  - Tracks CSVs    → Tracks layers for each cell type
All loaded lazily via dask to save memory.
"""
from __future__ import annotations

import traceback
from pathlib import Path

import yaml
import dask.array as da
import numpy as np
import pandas as pd
from qtpy.QtCore import Qt, QPropertyAnimation, QEasingCurve, Property, QRect, QPoint, QSize, QTimer
from qtpy.QtGui import QPainter, QColor
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QLabel,
    QComboBox,
    QPushButton,
    QCheckBox,
    QTextEdit,
    QStackedWidget,
    QAbstractButton,
    QSizePolicy,
    QMessageBox,
    QScrollArea,
    QFrame,
    QLayout,
)

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    is_multicolor_celltype,
)

# Channel colormaps (cycled if there are many channels)
_CHANNEL_COLORS = ["cyan", "yellow", "green", "red", "blue", "magenta"]
# Label colormaps per cell-type category
_LABEL_CMAP = {
    "or": "viridis",
    "im": "inferno",
    "ot": "plasma",
}
# ----------------------------------------------------------------------
# Dask cache-buster
# ----------------------------------------------------------------------
def _bust_dask_cache(arr: "da.Array") -> "da.Array":
    """Return a dask array with a unique task-graph name.

    ``da.from_zarr`` generates task-graph keys from the store path, so a
    second call for the same path may yield keys that overlap with results
    still held in dask's or napari's internal slice cache.  Wrapping with a
    UUID-based name forces dask (and napari) to treat the result as a brand-
    new array and re-read every block from disk.
    """
    import uuid
    try:
        return arr.map_blocks(
            lambda b: b,
            dtype=arr.dtype,
            name=f"zarr_load_{uuid.uuid4().hex}",
        )
    except Exception:
        return arr


# ----------------------------------------------------------------------
# Flow layout — wraps children to a new row when horizontal space is tight
# ----------------------------------------------------------------------
class _FlowLayout(QLayout):
    """Geometry manager that flows items left-to-right and wraps to a new row.

    Drop-in replacement for ``QHBoxLayout`` when you want automatic wrapping:
    items are placed left-to-right until the row would overflow, at which point
    a new row is started.  The enclosing widget's height expands automatically
    because :meth:`hasHeightForWidth` / :meth:`heightForWidth` are implemented.
    """

    def __init__(self, parent=None, h_spacing: int = 8, v_spacing: int = 4):
        super().__init__(parent)
        self._h_spacing = h_spacing
        self._v_spacing = v_spacing
        self._items: list = []

    # ------------------------------------------------------------------
    # QLayout interface
    # ------------------------------------------------------------------
    def addItem(self, item) -> None:
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
        return Qt.Orientations(Qt.Orientation(0))

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        return self._do_layout(QRect(0, 0, width, 0), test=True)

    def setGeometry(self, rect: QRect) -> None:
        super().setGeometry(rect)
        self._do_layout(rect, test=False)

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
    # Internal geometry calculation
    # ------------------------------------------------------------------
    def _do_layout(self, rect: QRect, test: bool) -> int:
        m = self.contentsMargins()
        eff = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        x, y, line_h = eff.x(), eff.y(), 0
        for item in self._items:
            hint = item.sizeHint()
            next_x = x + hint.width() + self._h_spacing
            if next_x - self._h_spacing > eff.right() and line_h > 0:
                # Wrap to the next row.
                x = eff.x()
                y += line_h + self._v_spacing
                next_x = x + hint.width() + self._h_spacing
                line_h = 0
            if not test:
                item.setGeometry(QRect(QPoint(x, y), hint))
            x = next_x
            line_h = max(line_h, hint.height())
        return y + line_h - eff.y() + m.bottom()


# ----------------------------------------------------------------------
# Custom Toggle Switch (QSwitch)
# ----------------------------------------------------------------------
class QSwitch(QAbstractButton):
    """Modern toggle switch widget."""
    def __init__(self, parent=None, track_radius=10, thumb_radius=8):
        super().__init__(parent)
        self.setCheckable(True)
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        self._track_radius = track_radius
        self._thumb_radius = thumb_radius
        
        self._margin = 3
        self._base_width = 40
        self._base_height = 22
        
        # Position of the thumb (0 to 1 float)
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
        
        # Draw background track
        if self.isChecked():
            # Green for ON
            p.setBrush(QColor("#4CAF50"))
        else:
            # Grey for OFF
            p.setBrush(QColor("#CCCCCC"))
            
        p.drawRoundedRect(0, 0, self.width(), self.height(), self._track_radius, self._track_radius)
        
        # Draw thumb
        p.setBrush(QColor("white"))
        thumb_x = self._margin + self._offset * (self.width() - 2 * self._thumb_radius - 2 * self._margin)
        p.drawEllipse(int(thumb_x), self._margin, 2 * self._thumb_radius, 2 * self._thumb_radius)

    def setChecked(self, checked: bool):
        super().setChecked(checked)
        self._offset = 1.0 if checked else 0.0
        self.update()



class VisualizationTab(QWidget):
    """Napari-integrated dataset viewer backed by dask arrays."""

    def __init__(self, viewer, data_prep, parent=None):
        """
        Parameters
        ----------
        viewer : napari.Viewer
            The napari viewer instance.
        data_prep : DataPreparationTab
            Reference to the data-preparation tab so we can read its metadata.
        """
        super().__init__(parent)
        self.viewer = viewer
        self.data_prep = data_prep

        # Wire up to the metadata_loaded signal
        self.data_prep.metadata_loaded.connect(self._on_metadata_loaded)

        self._metadata: pd.DataFrame | None = None

        # Debounced save for viewer display settings (contrast / colormap)
        self._display_save_timer = QTimer(self)
        self._display_save_timer.setSingleShot(True)
        self._display_save_timer.setInterval(1000)
        self._display_save_timer.timeout.connect(self._persist_viewer_display)

        # ------- Build UI -------
        self.main_layout = QVBoxLayout(self)
        self.main_layout.setContentsMargins(0, 0, 0, 0)

        self.stack = QStackedWidget()
        self.main_layout.addWidget(self.stack)

        # -- Page 0: Placeholder --
        self.placeholder_page = QWidget()
        place_lay = QVBoxLayout(self.placeholder_page)
        self.placeholder_label = QLabel("Load metadata in the Data Preparation tab to see visualization options.")
        self.placeholder_label.setAlignment(Qt.AlignCenter)
        self.placeholder_label.setStyleSheet("color: #888; font-style: italic; font-size: 14px; padding: 20px;")
        place_lay.addWidget(self.placeholder_label)
        self.stack.addWidget(self.placeholder_page)

        # -- Page 1: Main Content --
        self.main_content = QWidget()
        layout = QVBoxLayout(self.main_content)
        layout.setAlignment(Qt.AlignTop)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # Dataset selector
        sel_grp = QGroupBox("Dataset")
        sel_lay = QVBoxLayout(sel_grp)

        row = QHBoxLayout()
        row.addWidget(QLabel("Sample:"))
        self.sample_combo = QComboBox()
        self.sample_combo.setMinimumWidth(200)
        row.addWidget(self.sample_combo, stretch=1)
        sel_lay.addLayout(row)

        # Options
        self.clear_layers_cb = QCheckBox("Clear existing layers before loading")
        self.clear_layers_cb.setChecked(True)
        # sel_lay.addWidget(self.clear_layers_cb) # Removed as it does not make sense right now

        self.info_label = QLabel("")
        self.info_label.setStyleSheet("color: #666; font-style: italic;")
        sel_lay.addWidget(self.info_label)

        btn_load = QPushButton("Load Dataset into Napari")
        btn_load.setStyleSheet(
            "background-color: #4CAF50; color: white; font-weight: bold; padding: 8px;"
        )
        btn_load.clicked.connect(self._on_load_dataset)
        sel_lay.addWidget(btn_load)

        layout.addWidget(sel_grp)

        # Visibility Toggles
        vis_grp = QGroupBox("Visibility Toggles")
        vis_lay = _FlowLayout(vis_grp, h_spacing=12, v_spacing=4)

        def make_switch(label, group_type):
            container_widget = QWidget()
            container = QHBoxLayout(container_widget)
            container.setContentsMargins(0, 0, 0, 0)
            container.setSpacing(4)
            container.addWidget(QLabel(label))
            sw = QSwitch()
            sw.setChecked(True)
            sw.clicked.connect(lambda: self._on_toggle_layer_group(sw.isChecked(), group_type))
            container.addWidget(sw)
            return container_widget, sw

        wgt_raw, self.toggle_raw = make_switch("Raw:", "raw")
        wgt_seg, self.toggle_seg = make_switch("Segments:", "segments")
        wgt_tseg, self.toggle_track_seg = make_switch("Tracked Segments:", "tracked_segments")
        wgt_tracks, self.toggle_tracks = make_switch("Tracks:", "tracks")

        # Hide tracks by default
        self.toggle_tracks.setChecked(False)

        vis_lay.addWidget(wgt_raw)
        vis_lay.addWidget(wgt_seg)
        vis_lay.addWidget(wgt_tseg)
        vis_lay.addWidget(wgt_tracks)

        layout.addWidget(vis_grp)

        # ------------------------------------------------------------------
        # Manual edition (visible only when a tracked-segments layer exists)
        # ------------------------------------------------------------------
        self._editor = None  # active TrackedSegmentEditor instance (or None)
        self._editor_dock = None  # napari dock holding the editor (notebook only)

        self.edition_grp = QGroupBox("Manual Edition")
        edit_lay = QVBoxLayout(self.edition_grp)
        edit_lay.setContentsMargins(6, 4, 6, 4)
        edit_lay.setSpacing(4)

        ct_row = QHBoxLayout()
        ct_row.addWidget(QLabel("Tracked segments:"))
        self.edit_combo = QComboBox()
        self.edit_combo.setMinimumWidth(220)
        ct_row.addWidget(self.edit_combo, stretch=1)
        edit_lay.addLayout(ct_row)

        self.edit_toggle_btn = QPushButton("✏️  Edit tracked segments")
        self.edit_toggle_btn.setCheckable(True)
        self._style_edit_toggle(checked=False)
        self.edit_toggle_btn.toggled.connect(self._on_edit_toggle)
        edit_lay.addWidget(self.edit_toggle_btn)

        # Slot reserved for the inline TrackedSegmentEditor; populated on
        # demand and torn down when editing stops.
        self.edit_container = QWidget()
        self.edit_container_layout = QVBoxLayout(self.edit_container)
        self.edit_container_layout.setContentsMargins(0, 0, 0, 0)
        self.edit_container_layout.setSpacing(0)
        edit_lay.addWidget(self.edit_container)

        self.edition_grp.setVisible(False)
        layout.addWidget(self.edition_grp)

        # Log
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log)

        # Wrap the main content in a scroll area so the manual-correction
        # editor (which can grow tall) does not push controls out of view.
        self.main_scroll = QScrollArea()
        self.main_scroll.setWidgetResizable(True)
        self.main_scroll.setFrameShape(QFrame.NoFrame)
        self.main_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.main_scroll.setWidget(self.main_content)
        self.stack.addWidget(self.main_scroll)
        self.stack.setCurrentIndex(0)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _log(self, msg: str):
        self.log.append(msg)
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    # ------------------------------------------------------------------
    # Slot: metadata loaded
    # ------------------------------------------------------------------
    def _on_metadata_loaded(self, metadata: pd.DataFrame):
        self._metadata = metadata
        samples = sorted(str(s) for s in metadata["sample_name"].unique())
        self.sample_combo.clear()
        self.sample_combo.addItems(samples)
        self.info_label.setText(f"{len(samples)} sample(s) available.  Select one and click Load.")
        self._log(f"Metadata received — {len(samples)} sample(s)")
        self.stack.setCurrentIndex(1)

    # ------------------------------------------------------------------
    # Load dataset into napari
    # ------------------------------------------------------------------
    def _on_load_dataset(self):
        if self._metadata is None:
            self.info_label.setText("⚠️  No metadata loaded.")
            return

        sample_name = self.sample_combo.currentText()
        if not sample_name:
            return

        row = self._metadata[self._metadata["sample_name"] == sample_name]
        if row.empty:
            self.info_label.setText(f"⚠️  Sample '{sample_name}' not found in metadata.")
            return
        row = row.iloc[0]

        # Reset toggles to default before loading
        self.toggle_raw.setChecked(True)
        self.toggle_seg.setChecked(True)
        self.toggle_track_seg.setChecked(True)
        self.toggle_tracks.setChecked(False) # Tracks are hidden by default

        # if self.clear_layers_cb.isChecked():
        # Stop any running napari dim animation before clearing layers
        # (prevents RuntimeError: wrapped C/C++ object has been deleted)
        try:
            qt_dims = self.viewer.window._qt_viewer.dims
            if hasattr(qt_dims, '_animation_thread') and qt_dims._animation_thread is not None:
                qt_dims._animation_thread.quit()
                qt_dims._animation_thread.wait()
        except Exception:
            pass
        self._display_save_timer.stop()
        self.viewer.layers.clear()

        output_dir = self.data_prep.output_dir or ""

        # ---- 1. Raw image (zarr preferred, dask-backed) ------------------
        self._load_raw_image(sample_name, row, output_dir)

        # ---- 2. Segments -------------------------------------------------
        self._load_segments(sample_name, row)

        # ---- 3. Tracked Segments -----------------------------------------
        self._load_tracked_segments(sample_name, row)

        # ---- 4. Tracks ---------------------------------------------------
        self._load_tracks(sample_name, row)
        self._log(f"✅ Dataset '{sample_name}' loaded into napari")

        # ---- 5. Refresh the Manual Edition selector ---------------------
        self._refresh_edition_panel()

    # ------------------------------------------------------------------
    # Raw image loader
    # ------------------------------------------------------------------
    def _load_raw_image(self, sample_name: str, row: pd.Series, output_dir: str):
        """Load the raw image as dask array, splitting channels into layers."""
        from behav3d.io.formats.zarr import load_zarr

        raw_path = str(row.get("raw_image_path", "")).strip()
        if not raw_path:
            self._log(f"  ⚠️ No raw_image_path for {sample_name}")
            return

        raw_p = Path(raw_path)

        # Try zarr first (either the raw path is .zarr, or a converted version exists)
        zarr_candidates = [
            raw_p if raw_p.suffix == ".zarr" else None,
            Path(output_dir, f"{sample_name}.zarr") if output_dir else None,
            raw_p.with_suffix(".zarr"),
        ]
        loaded = False
        for zp in zarr_candidates:
            if zp and zp.exists():
                try:
                    dask_img = load_zarr(zp)
                    self._log(f"  📂 Loaded raw zarr: {zp}  shape={dask_img.shape}")
                    self._add_channels(dask_img, sample_name, row)
                    loaded = True
                    break
                except Exception as e:
                    self._log(f"  ⚠️ Error loading {zp}: {e}")

        if not loaded:
            # Fall back to full image load (may be large!)
            if raw_p.exists():
                try:
                    from behav3d.io.images import load_image
                    img = load_image(raw_p)
                    if not isinstance(img, da.Array):
                        img = da.from_array(img, chunks=(1,) + img.shape[1:])
                    self._log(f"  📂 Loaded raw image: {raw_p}  shape={img.shape}")
                    self._add_channels(img, sample_name, row)
                except Exception as e:
                    self._log(f"  ❌ Could not load raw image: {e}")
            else:
                self._log(f"  ⚠️ Raw image not found: {raw_path}")

    def _add_channels(self, dask_img: da.Array, sample_name: str, row: pd.Series):
        """Split a TCZYX array along C and add each channel as a napari Image layer."""
        # Determine dim order
        dim_order = str(row.get("dimension_order", "TCZYX")).strip() or "TCZYX"

        # Normalise to TCZYX for consistent indexing
        if dim_order != "TCZYX" and len(dim_order) == 5:
            try:
                axes = [dim_order.index(d) for d in "TCZYX"]
                dask_img = da.transpose(dask_img, axes)
            except ValueError:
                pass  # fallback: assume already TCZYX

        n_channels = dask_img.shape[1] if dask_img.ndim >= 5 else 1

        saved_channels = (
            self.data_prep.behav3d_parameters
            .get("viewer_display", {})
            .get("channels", {})
        )

        for c in range(n_channels):
            if dask_img.ndim >= 5:
                ch_data = dask_img[:, c, :, :, :]  # T Z Y X
            else:
                ch_data = dask_img

            saved = saved_channels.get(c) or saved_channels.get(str(c))
            color = (saved["colormap"] if saved and "colormap" in saved
                     else _CHANNEL_COLORS[c % len(_CHANNEL_COLORS)])
            layer_name = f"{sample_name} – Ch{c}"

            add_kwargs = dict(
                name=layer_name,
                colormap=color,
                blending="additive",
                visible=True,
            )
            if saved and "contrast_limits" in saved:
                add_kwargs["contrast_limits"] = tuple(saved["contrast_limits"])

            layer = self.viewer.add_image(ch_data, **add_kwargs)

            layer.events.contrast_limits.connect(self._on_layer_display_changed)
            layer.events.colormap.connect(self._on_layer_display_changed)

            self._log(f"    + Image layer: {layer_name}  ({color})")

    # ------------------------------------------------------------------
    # Viewer display persistence
    # ------------------------------------------------------------------
    def _on_layer_display_changed(self, event=None):
        self._display_save_timer.start()

    def _persist_viewer_display(self):
        """Snapshot current channel display settings and save to YAML."""
        out_dir = self.data_prep.output_dir
        if not out_dir:
            return

        params = self.data_prep.behav3d_parameters
        channels_cfg = params.setdefault("viewer_display", {}).setdefault("channels", {})

        for layer in self.viewer.layers:
            if " – Ch" not in layer.name:
                continue
            try:
                ch_idx = int(layer.name.rsplit("Ch", 1)[1])
            except (ValueError, IndexError):
                continue
            channels_cfg[ch_idx] = {
                "colormap": str(layer.colormap.name),
                "contrast_limits": [float(layer.contrast_limits[0]),
                                    float(layer.contrast_limits[1])],
            }

        params_path = Path(out_dir) / "behav3d_parameters.yml"
        try:
            with open(params_path, "w") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        except Exception as e:
            self._log(f"Warning: Could not save display settings: {e}")

    # ------------------------------------------------------------------
    # Segment loader
    # ------------------------------------------------------------------
    def _load_segments(self, sample_name: str, row: pd.Series):
        """Load segment images for each detected cell type."""
        from behav3d.io.formats.zarr import load_zarr

        ct_map = self._detect_cell_type_columns(row)

        for ct_name, prefix in ct_map.items():
            seg_col = f"{prefix}_{ct_name}_segments_image_path" if prefix else f"{ct_name}_segments_image_path"
            seg_path_val = row.get(seg_col)
            
            # Robust check for missing values (NaN or empty)
            if pd.isna(seg_path_val) or not str(seg_path_val).strip():
                self._log(f"    - Skipping {ct_name} segments: Path not defined in metadata")
                continue
            
            seg_path = str(seg_path_val).strip()
            if not Path(seg_path).exists():
                # Also try .zarr version in output dir as a fallback
                if self.data_prep.output_dir:
                    zarr_guess = Path(self.data_prep.output_dir, f"{Path(seg_path).stem}.zarr")
                    if zarr_guess.exists():
                        seg_path = str(zarr_guess)
                    else:
                        self._log(f"    - Skipping {ct_name} segments: File not found ({seg_path})")
                        continue
                else:
                    self._log(f"    - Skipping {ct_name} segments: File not found ({seg_path})")
                    continue

            try:
                seg_p = Path(seg_path)
                if seg_p.suffix == ".zarr":
                    seg_data = _bust_dask_cache(load_zarr(seg_p))
                else:
                    from behav3d.io.images import load_image
                    seg_data = load_image(seg_p)
                    if not isinstance(seg_data, da.Array):
                        seg_data = da.from_array(seg_data, chunks=(1,) + seg_data.shape[1:])

                layer_name = f"{sample_name} – {ct_name} segments"
                self.viewer.add_labels(
                    seg_data,
                    name=layer_name,
                    visible=True,
                )
                self._log(f"    + Labels layer: {layer_name}")
            except Exception as e:
                self._log(f"    ⚠️ Could not load segments for {ct_name}: {e}")

        # Load death mask if present
        dead_col = "dead_mask_path"
        dead_path_val = row.get(dead_col)
        if not pd.isna(dead_path_val) and str(dead_path_val).strip():
            dead_path = str(dead_path_val).strip()
            if Path(dead_path).exists():
                try:
                    dead_p = Path(dead_path)
                    if dead_p.suffix == ".zarr":
                        dead_data = _bust_dask_cache(load_zarr(dead_p))
                    else:
                        from behav3d.io.images import load_image
                        dead_data = load_image(dead_p)
                        if not isinstance(dead_data, da.Array):
                            dead_data = da.from_array(dead_data, chunks=(1,) + dead_data.shape[1:])

                    layer_name = f"{sample_name} – dead mask"
                    self.viewer.add_labels(
                        dead_data,
                        name=layer_name,
                        visible=False, # Hidden by default
                    )
                    self._log(f"    + Dead mask layer: {layer_name}")
                except Exception as e:
                    self._log(f"    ⚠️ Could not load dead mask: {e}")
            else:
                self._log(f"    - Skipping dead mask: File not found ({dead_path})")

    # ------------------------------------------------------------------
    # Track loader
    # ------------------------------------------------------------------
    def _load_tracks(self, sample_name: str, row: pd.Series):
        """Load tracks CSVs and add them as napari Tracks layers."""
        ct_map = self._detect_cell_type_columns(row)

        for ct_name, prefix in ct_map.items():
            tracks_csv_col = f"{prefix}_{ct_name}_tracks_csv_path" if prefix else f"{ct_name}_tracks_csv_path"
            tracks_csv_val = row.get(tracks_csv_col)
            
            # Robust check for missing values (NaN or empty)
            if pd.isna(tracks_csv_val) or not str(tracks_csv_val).strip():
                self._log(f"    - Skipping {ct_name} tracks: Path not defined in metadata")
                continue

            csv_path_str = str(tracks_csv_val).strip()
            if not Path(csv_path_str).exists():
                self._log(f"    - Skipping {ct_name} tracks: File not found ({csv_path_str})")
                continue

            try:
                tracks_df = pd.read_csv(csv_path_str)

                # Expect columns: track_id, t, z, y, x (at minimum)
                required = {"TrackID", "position_t", "pixel_position_y", "pixel_position_x"}
                if not required.issubset(set(tracks_df.columns)):
                    self._log(f"    ⚠️ Tracks CSV for {ct_name} missing columns {required - set(tracks_df.columns)}")
                    continue

                if "pixel_position_z" in tracks_df.columns:
                    track_data = tracks_df[["TrackID", "position_t", "pixel_position_z", "pixel_position_y", "pixel_position_x"]].values
                else:
                    track_data = tracks_df[["TrackID", "position_t", "pixel_position_y", "pixel_position_x"]].values

                layer_name = f"{sample_name} – {ct_name} tracks"
                self.viewer.add_tracks(
                    track_data,
                    name=layer_name,
                    visible=False,
                )
                self._log(f"    + Tracks layer: {layer_name}")
            except Exception as e:
                self._log(f"    ⚠️ Could not load tracks for {ct_name}: {e}")

    # ------------------------------------------------------------------
    # Tracked Segment loader
    # ------------------------------------------------------------------
    def _load_tracked_segments(self, sample_name: str, row: pd.Series):
        """Load tracked segments (labeled images with track IDs)."""
        from behav3d.io.formats.zarr import load_zarr
        
        ct_map = self._detect_cell_type_columns(row)

        for ct_name, prefix in ct_map.items():
            track_seg_col = f"{prefix}_{ct_name}_tracks_image_path" if prefix else f"{ct_name}_tracks_image_path"
            track_seg_path_val = row.get(track_seg_col)
            
            if pd.isna(track_seg_path_val) or not str(track_seg_path_val).strip():
                continue
                
            track_seg_path = str(track_seg_path_val).strip()
            if not Path(track_seg_path).exists():
                continue

            try:
                track_seg_p = Path(track_seg_path)
                if track_seg_p.suffix == ".zarr":
                    seg_data = _bust_dask_cache(load_zarr(track_seg_p))
                else:
                    from behav3d.io.images import load_image
                    seg_data = load_image(track_seg_p)
                    if not isinstance(seg_data, da.Array):
                        seg_data = da.from_array(seg_data, chunks=(1,) + seg_data.shape[1:])

                layer_name = f"{sample_name} – {ct_name} tracked segments"
                self.viewer.add_labels(
                    seg_data,
                    name=layer_name,
                    visible=not is_multicolor_celltype(ct_name),
                )
                self._log(f"    + Tracked Labels layer: {layer_name}")
                
                # If we have tracked segments, hide the regular segments
                reg_seg_layer_name = f"{sample_name} – {ct_name} segments"
                if reg_seg_layer_name in self.viewer.layers:
                    self.viewer.layers[reg_seg_layer_name].visible = False
                    self._log(f"    - Hiding regular segments for {ct_name}")
                    # Update toggle state to reflect that some segments are hidden
                    # Note: This affects ALL segments if multiple cell types exist, 
                    # but usually it's consistent.
                    self.toggle_seg.setChecked(False)

            except Exception as e:
                self._log(f"    ⚠️ Could not load tracked segments for {ct_name}: {e}")

    def _on_toggle_layer_group(self, state, group_type):
        """Batch show/hide layers based on their name suffix/pattern."""
        visible = bool(state)
        
        for layer in self.viewer.layers:
            name = layer.name
            match = False
            
            if group_type == "raw":
                # Raw images usually: "sample_name – Ch0"
                if " – Ch" in name:
                    match = True
            elif group_type == "segments":
                # Segments usually: "sample_name – celltype segments"
                if name.endswith(" segments") and not name.endswith(" tracked segments"):
                    match = True
            elif group_type == "tracked_segments":
                # Tracked segments: "sample_name – celltype tracked segments"
                if name.endswith(" tracked segments"):
                    body = name[:-17].rstrip()
                    if " – " in body:
                        ct_name = body.rsplit(" – ", 1)[1].strip()
                        if not is_multicolor_celltype(ct_name):
                            match = True
            elif group_type == "tracks":
                # Tracks: "sample_name – celltype tracks"
                if name.endswith(" tracks"):
                    body = name[:-7].rstrip()
                    if " – " in body:
                        ct_name = body.rsplit(" – ", 1)[1].strip()
                        if not is_multicolor_celltype(ct_name):
                            match = True
                    
            if match:
                layer.visible = visible

    # ------------------------------------------------------------------
    # Manual Edition section
    # ------------------------------------------------------------------
    def _style_edit_toggle(self, checked: bool) -> None:
        if checked:
            self.edit_toggle_btn.setText("⏹  Stop editing")
            self.edit_toggle_btn.setStyleSheet(
                "background-color:#c62828; color:white; font-weight:bold; padding:6px;"
            )
        else:
            self.edit_toggle_btn.setText("✏️  Edit tracked segments")
            self.edit_toggle_btn.setStyleSheet(
                "background-color:#28a745; color:white; font-weight:bold; padding:6px;"
            )

    def _tracked_layer_names(self) -> list[str]:
        names = []
        for layer in self.viewer.layers:
            name = layer.name
            if isinstance(name, str) and name.endswith(" tracked segments"):
                body = name[:-17].rstrip()
                if " – " in body:
                    ct_name = body.rsplit(" – ", 1)[1].strip()
                    if not is_multicolor_celltype(ct_name):
                        names.append(name)
        return names

    def _refresh_edition_panel(self) -> None:
        """Show the Manual Edition group iff a tracked-segments layer exists."""
        names = self._tracked_layer_names()
        # If we changed sample while editing, force-close the editor.
        if self._editor is not None and self._editor.layer_name not in names:
            self.edit_toggle_btn.setChecked(False)
        prev = self.edit_combo.currentText()
        self.edit_combo.blockSignals(True)
        self.edit_combo.clear()
        self.edit_combo.addItems(names)
        if prev in names:
            self.edit_combo.setCurrentText(prev)
        self.edit_combo.blockSignals(False)
        self.edition_grp.setVisible(bool(names))

    def _on_edit_toggle(self, checked: bool) -> None:
        if checked:
            if not self._open_editor():
                # Failed to open — revert.
                self.edit_toggle_btn.blockSignals(True)
                self.edit_toggle_btn.setChecked(False)
                self.edit_toggle_btn.blockSignals(False)
                return
            self._style_edit_toggle(True)
        else:
            ok = self._close_editor(prompt=True)
            if not ok:
                # User cancelled the close prompt — keep the toggle on.
                self.edit_toggle_btn.blockSignals(True)
                self.edit_toggle_btn.setChecked(True)
                self.edit_toggle_btn.blockSignals(False)
                return
            self._style_edit_toggle(False)

    def _open_editor(self) -> bool:
        from behav3d.napari._segment_editor import (
            TrackedSegmentEditor,
            _cell_type_from_layer_name,
            _sample_from_layer_name,
        )
        from behav3d.napari._loaders import (
            load_raw_into_viewer,
            load_tracked_segments_into_viewer,
            load_tracks_into_viewer,
        )

        layer_name = self.edit_combo.currentText()
        if not layer_name:
            self._log("⚠️ No tracked-segments layer to edit.")
            return False
        sample_name = _sample_from_layer_name(layer_name)
        cell_type = _cell_type_from_layer_name(layer_name)
        if not sample_name or not cell_type or self._metadata is None:
            self._log("⚠️ Could not infer sample/cell-type from layer name.")
            return False
        rows = self._metadata[
            self._metadata["sample_name"].astype(str) == str(sample_name)
        ]
        if rows.empty:
            self._log(f"⚠️ Sample '{sample_name}' not in metadata.")
            return False
        row = rows.iloc[0]
        output_dir = self.data_prep.output_dir or ""

        # Hard reset the viewer: keep ONLY the raw images, the chosen
        # tracked-segments layer (loaded as numpy so it is editable in
        # place) and the matching tracks for that cell type.
        self._log(
            f"  Entering edit mode — clearing viewer and reloading "
            f"raw + '{cell_type} tracked segments' + '{cell_type} tracks' only."
        )
        try:
            qt_dims = self.viewer.window._qt_viewer.dims
            if hasattr(qt_dims, "_animation_thread") and qt_dims._animation_thread is not None:
                qt_dims._animation_thread.quit()
                qt_dims._animation_thread.wait()
        except Exception:
            pass
        self._display_save_timer.stop()
        self.viewer.layers.clear()
        saved_channels = (
            self.data_prep.behav3d_parameters
            .get("viewer_display", {})
            .get("channels", {})
        )
        load_raw_into_viewer(
            self.viewer, sample_name, row,
            output_dir=output_dir, log=self._log,
            display_settings=saved_channels,
        )
        # Load the tracked-segments layer dask-backed (as_numpy=False) so
        # napari can display it immediately while materialisation happens in
        # the background.  The editor will kick off a QThread worker that
        # converts it to a writable numpy array before enabling editing.
        added = load_tracked_segments_into_viewer(
            self.viewer, sample_name, row, log=self._log,
            only_cell_type=cell_type, as_numpy=False,
        )
        load_tracks_into_viewer(
            self.viewer, sample_name, row, visible=True, log=self._log,
            only_cell_type=cell_type,
        )
        # Refresh combo (layer names may have shifted, but the chosen one
        # is still present if the load succeeded).
        self._refresh_edition_panel()
        if not added:
            self._log("❌ Tracked segments failed to load — editor not opened.")
            return False
        # The just-added tracked-segments layer is the one to edit.
        target = added[0]

        try:
            editor = TrackedSegmentEditor(
                viewer=self.viewer,
                metadata_loader=self.data_prep,
                layer_name=target,
                parent=self,
            )
        except Exception as exc:
            traceback.print_exc()
            QMessageBox.critical(
                self, "Cannot open editor",
                f"Failed to start the manual-correction editor:\n{exc}",
            )
            return False
        # Hook the editor's "Stop editing" footer button to the toggle.
        editor.stop_callback = lambda: self.edit_toggle_btn.setChecked(False)
        self.edit_container_layout.addWidget(editor)
        self._editor = editor
        self._log(f"  Editing started on '{target}'.")
        # Kick off background materialisation now that the editor is shown.
        editor.start_materialisation()
        return True

    def _close_editor(self, prompt: bool = True) -> bool:
        if self._editor is None:
            return True
        if prompt:
            if not self._editor.request_exit():
                return False
        else:
            self._editor._cleanup()
        # Remove from layout + delete.
        self.edit_container_layout.removeWidget(self._editor)
        self._editor.setParent(None)
        try:
            self._editor.deleteLater()
        except Exception:
            pass
        self._editor = None
        self._log("  Editing stopped. Reloading visualizer...")
        self._on_load_dataset()
        return True

    def request_tab_exit(self) -> bool:
        """Return False to veto leaving the Visualization tab while edits
        are pending.  Called by ``BEHAV3DWidget._on_tab_changed``.
        """
        if self._editor is None:
            return True
        if self._editor.request_exit():
            # request_exit() already saved/discarded; tear down the editor.
            self._close_editor(prompt=False)
            self.edit_toggle_btn.blockSignals(True)
            self.edit_toggle_btn.setChecked(False)
            self._style_edit_toggle(False)
            self.edit_toggle_btn.blockSignals(False)
            return True
        return False

    # ------------------------------------------------------------------
    # Detect cell-type columns from metadata row
    # ------------------------------------------------------------------
    @staticmethod
    def _detect_cell_type_columns(row: pd.Series) -> dict[str, str]:
        """Return {cell_type_name: prefix} for all detected cell types."""
        import re

        ct_map: dict[str, str] = {}
        pattern = re.compile(r"^(or|im|ot)_(.+?)_line_condition$")
        for col in row.index:
            m = pattern.match(col)
            if m:
                prefix, ct_name = m.group(1), m.group(2)
                ct_map[ct_name] = prefix
                
        merged_pattern = re.compile(r"^(.+?)_(segments_image_path|tracks_image_path|tracks_csv_path)$")
        for col in row.index:
            if col.startswith(("or_", "im_", "ot_")):
                continue
            m = merged_pattern.match(col)
            if m:
                ct_name = m.group(1)
                if ct_name.endswith("_merged") or ct_name.endswith("_grouped"):
                    if ct_name not in ct_map:
                        ct_map[ct_name] = ""
                        
        return ct_map
