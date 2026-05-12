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

import dask.array as da
import numpy as np
import pandas as pd
from qtpy.QtCore import Qt, QPropertyAnimation, QEasingCurve, Property, QRect, QPoint
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
)

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
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
        vis_lay = QHBoxLayout(vis_grp)
        
        def make_switch(label, group_type):
            container = QHBoxLayout()
            container.addWidget(QLabel(label))
            sw = QSwitch()
            sw.setChecked(True)
            sw.clicked.connect(lambda: self._on_toggle_layer_group(sw.isChecked(), group_type))
            container.addWidget(sw)
            return container, sw

        lay_raw, self.toggle_raw = make_switch("Raw:", "raw")
        lay_seg, self.toggle_seg = make_switch("Segments:", "segments")
        lay_tseg, self.toggle_track_seg = make_switch("Tracked Segments:", "tracked_segments")
        lay_tracks, self.toggle_tracks = make_switch("Tracks:", "tracks")
        
        # Hide tracks by default
        self.toggle_tracks.setChecked(False)

        vis_lay.addLayout(lay_raw)
        vis_lay.addLayout(lay_seg)
        vis_lay.addLayout(lay_tseg)
        vis_lay.addLayout(lay_tracks)
        
        layout.addWidget(vis_grp)

        # Log
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log)

        self.stack.addWidget(self.main_content)
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

        for c in range(n_channels):
            if dask_img.ndim >= 5:
                ch_data = dask_img[:, c, :, :, :]  # T Z Y X
            else:
                ch_data = dask_img

            color = _CHANNEL_COLORS[c % len(_CHANNEL_COLORS)]
            layer_name = f"{sample_name} – Ch{c}"

            self.viewer.add_image(
                ch_data,
                name=layer_name,
                colormap=color,
                blending="additive",
                visible=True,
            )
            self._log(f"    + Image layer: {layer_name}  ({color})")

    # ------------------------------------------------------------------
    # Segment loader
    # ------------------------------------------------------------------
    def _load_segments(self, sample_name: str, row: pd.Series):
        """Load segment images for each detected cell type."""
        from behav3d.io.formats.zarr import load_zarr

        ct_map = self._detect_cell_type_columns(row)

        for ct_name, prefix in ct_map.items():
            seg_col = f"{prefix}_{ct_name}_segments_image_path"
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
                    seg_data = load_zarr(seg_p)
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

    # ------------------------------------------------------------------
    # Track loader
    # ------------------------------------------------------------------
    def _load_tracks(self, sample_name: str, row: pd.Series):
        """Load tracks CSVs and add them as napari Tracks layers."""
        ct_map = self._detect_cell_type_columns(row)

        for ct_name, prefix in ct_map.items():
            csv_col = f"{prefix}_{ct_name}_tracks_csv_path"
            csv_path_val = row.get(csv_col)
            
            # Robust check for missing values (NaN or empty)
            if pd.isna(csv_path_val) or not str(csv_path_val).strip():
                self._log(f"    - Skipping {ct_name} tracks: Path not defined in metadata")
                continue

            csv_path_str = str(csv_path_val).strip()
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
            track_seg_col = f"{prefix}_{ct_name}_tracks_image_path"
            track_seg_path_val = row.get(track_seg_col)
            
            if pd.isna(track_seg_path_val) or not str(track_seg_path_val).strip():
                continue
                
            track_seg_path = str(track_seg_path_val).strip()
            if not Path(track_seg_path).exists():
                continue

            try:
                track_seg_p = Path(track_seg_path)
                if track_seg_p.suffix == ".zarr":
                    seg_data = load_zarr(track_seg_p)
                else:
                    from behav3d.io.images import load_image
                    seg_data = load_image(track_seg_p)
                    if not isinstance(seg_data, da.Array):
                        seg_data = da.from_array(seg_data, chunks=(1,) + seg_data.shape[1:])

                layer_name = f"{sample_name} – {ct_name} tracked segments"
                self.viewer.add_labels(
                    seg_data,
                    name=layer_name,
                    visible=True,
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
                    match = True
            elif group_type == "tracks":
                # Tracks: "sample_name – celltype tracks"
                if name.endswith(" tracks"):
                    match = True
                    
            if match:
                layer.visible = visible

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
        return ct_map
