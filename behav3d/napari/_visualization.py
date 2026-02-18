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
from qtpy.QtCore import Qt
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
)

from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)

# Channel colormaps (cycled if there are many channels)
_CHANNEL_COLORS = ["cyan", "yellow", "magenta", "green", "red", "blue"]
# Label colormaps per cell-type category
_LABEL_CMAP = {
    "or": "viridis",
    "im": "inferno",
    "ot": "plasma",
}


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
        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignTop)

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
        sel_lay.addWidget(self.clear_layers_cb)

        btn_load = QPushButton("Load Dataset into Napari")
        btn_load.setStyleSheet(
            "background-color: #4CAF50; color: white; font-weight: bold; padding: 8px;"
        )
        btn_load.clicked.connect(self._on_load_dataset)
        sel_lay.addWidget(btn_load)

        layout.addWidget(sel_grp)

        # Info panel
        self.info_label = QLabel("Load metadata in the Data Preparation tab first.")
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        # Log
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(150)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log)

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

        if self.clear_layers_cb.isChecked():
            self.viewer.layers.clear()

        output_dir = self.data_prep.output_dir or ""

        # ---- 1. Raw image (zarr preferred, dask-backed) ------------------
        self._load_raw_image(sample_name, row, output_dir)

        # ---- 2. Segments -------------------------------------------------
        self._load_segments(sample_name, row)

        # ---- 3. Tracks ---------------------------------------------------
        self._load_tracks(sample_name, row)

        self.info_label.setText(f"✅  Loaded dataset for '{sample_name}'")
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
                visible=(c == 0),  # only first channel visible by default
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
            seg_path = str(row.get(seg_col, "")).strip()
            if not seg_path or not Path(seg_path).exists():
                # Also try .zarr version in output dir
                if self.data_prep.output_dir:
                    zarr_guess = Path(self.data_prep.output_dir, f"{Path(seg_path).stem}.zarr") if seg_path else None
                    if zarr_guess and zarr_guess.exists():
                        seg_path = str(zarr_guess)
                    else:
                        continue
                else:
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
                    visible=False,
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
            csv_path_str = str(row.get(csv_col, "")).strip()
            if not csv_path_str or not Path(csv_path_str).exists():
                continue

            try:
                tracks_df = pd.read_csv(csv_path_str)

                # Expect columns: track_id, t, z, y, x (at minimum)
                required = {"track_id", "t", "y", "x"}
                if not required.issubset(set(tracks_df.columns)):
                    self._log(f"    ⚠️ Tracks CSV for {ct_name} missing columns {required - set(tracks_df.columns)}")
                    continue

                if "z" in tracks_df.columns:
                    track_data = tracks_df[["track_id", "t", "z", "y", "x"]].values
                else:
                    track_data = tracks_df[["track_id", "t", "y", "x"]].values

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
