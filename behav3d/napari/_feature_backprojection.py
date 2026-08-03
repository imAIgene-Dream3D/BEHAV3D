"""
BEHAV3D napari plugin – Feature Backprojection preview tab.

Quick, current-timepoint-only preview: pick a cell type, sample, and
numeric feature (e.g. ``elongation``, ``nr_dead_mask_pixels``) from the
combined track-features CSV, and see it backprojected onto the
segmentation labels for whichever timepoint napari's own time slider is
currently on.

Unlike the full backprojection already available under the Single Cell
tab (state/track cluster backprojection — writes a full per-sample zarr,
can take several minutes) or the notebook's "Backproject results"
section, this recomputes only the single currently-viewed frame on every
``viewer.dims.events.current_step`` change, mirroring the Feature
Extraction tab's Dead/Alive threshold preview (``_feature_extraction.py``,
``_refresh_preview_dead_layers``).
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from qtpy.QtWidgets import (
    QComboBox,
    QFormLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from behav3d.core.qt_help import make_help_row
from behav3d.napari._analysis import _detect_cell_types
from behav3d.napari._single_cell import _bp_add_raw_channels

_FEATURE_BP_EXCLUDED_COLUMNS = {
    "TrackID", "original_TrackID", "sample_name", "position_t",
    "position_x", "position_y", "position_z", "ClusterID", "UMAP1", "UMAP2",
}

_FEATURE_LAYER_PREFIX = "[Feature BP]"


class FeatureBackprojectionTab(QWidget):
    """Pick a cell type + numeric feature, preview it backprojected for
    whichever timepoint the viewer is currently showing."""

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        self._df_features: Optional[pd.DataFrame] = None
        self._df_features_csv_path: Optional[Path] = None
        self._df_features_cell_type: Optional[str] = None

        self._tracked_path: Optional[Path] = None
        self._sample: Optional[str] = None
        self._cell_type: Optional[str] = None
        self._feature_col: Optional[str] = None
        self._feature_layer_name: Optional[str] = None
        self._added_layer_names: list[str] = []
        self._dims_callback = None

        self._init_ui()

        if metadata_loader is not None and hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)
        self._on_metadata_updated()

    # ── UI ──────────────────────────────────────────────────────────────
    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(6)

        outer.addWidget(QLabel(
            "<div style='font-size:14px;font-weight:700;'>🔬 Feature Backprojection</div>"
        ))
        desc = QLabel(
            "Pick a cell type, sample, and numeric feature, then preview it "
            "backprojected onto the segmentation for whichever timepoint the "
            "viewer's time slider is currently on. Only the current frame is "
            "recomputed on each scrub — fast, and nothing is written to disk."
        )
        desc.setWordWrap(True)
        desc.setStyleSheet("color:#888;font-size:11px;")
        outer.addWidget(desc)

        form = QFormLayout()
        form.setSpacing(4)

        self.combo_cell_type = QComboBox()
        self.combo_cell_type.currentTextChanged.connect(self._on_cell_type_changed)
        form.addRow("Cell type:", self.combo_cell_type)

        self.combo_sample = QComboBox()
        self.combo_sample.currentTextChanged.connect(self._on_sample_changed)
        form.addRow("Sample:", self.combo_sample)

        self.combo_feature = QComboBox()
        form.addRow("Feature:", make_help_row(
            self.combo_feature, "Feature",
            "Numeric column from the combined track-features CSV for the "
            "selected cell type (the filtered CSV is used if it exists, "
            "otherwise the unfiltered one). Examples: elongation, "
            "nr_dead_mask_pixels, percentage_dead_mask, speed.",
        ))

        outer.addLayout(form)

        self.btn_show = QPushButton("▶ Show Preview")
        self.btn_show.setEnabled(False)
        self.btn_show.clicked.connect(self._on_show_clicked)
        outer.addWidget(self.btn_show)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color:#c66;font-size:11px;")
        outer.addWidget(self.status_label)

        outer.addStretch()

    # ── metadata / options ─────────────────────────────────────────────
    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _feature_csv_path(self, cell_type: str) -> Optional[Path]:
        out_dir = self._out_dir()
        if out_dir is None or not cell_type:
            return None
        base = out_dir / "analysis" / cell_type / "track_features"
        filtered = base / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"
        if filtered.exists():
            return filtered
        return base / f"BEHAV3D_{cell_type}_combined_track_features.csv"

    def _resolve_raw_path(self, row, out_dir: Path, sample: str) -> Optional[Path]:
        """Raw image path for `sample`.

        Prefers the ``raw_image_path`` column in metadata — the only
        source of truth for *combined* experiments (see
        ``behav3d.core.combine_experiments``), where samples keep pointing
        at their original per-experiment locations and are never copied
        under the combined run's `output_dir/images/`. Falls back to the
        `output_dir/images/{sample}/...` convention for older metadata
        that doesn't carry this column.
        """
        if row is not None and "raw_image_path" in row and pd.notna(row["raw_image_path"]):
            p = Path(str(row["raw_image_path"]).strip())
            if p.exists():
                return p

        from behav3d.analysis.behavior.state.visualization.backprojection import (
            _resolve_raw_image_path,
        )
        fallback = _resolve_raw_image_path(out_dir, sample, verbose=False)
        return fallback if fallback is not None and Path(fallback).exists() else None

    def _resolve_tracked_path(self, row, out_dir: Path, sample: str, cell_type: str) -> Optional[Path]:
        """Tracked-labels image path for `sample`/`cell_type` — same
        metadata-column-first, output_dir-glob-fallback strategy as
        `_resolve_raw_path` (and matching `DeathThresholdPreview`'s
        `_resolve_label_path` in the notebook GUI)."""
        if row is not None:
            for prefix in ("or", "im", "ot"):
                col = f"{prefix}_{cell_type}_tracks_image_path"
                if col in row and pd.notna(row[col]) and str(row[col]).strip():
                    p = Path(str(row[col]).strip())
                    if p.exists():
                        return p

        from behav3d.analysis.behavior.state.visualization.backprojection import (
            _resolve_tracked_image_path,
        )
        fallback = _resolve_tracked_image_path(out_dir, sample, cell_type, verbose=False)
        return fallback if fallback is not None and Path(fallback).exists() else None

    def _on_metadata_updated(self, *_):
        organoid_types, immune_types, other_types = _detect_cell_types(self.metadata_loader)
        cell_types = organoid_types + immune_types + other_types

        prev_ct = self.combo_cell_type.currentText()
        self.combo_cell_type.blockSignals(True)
        self.combo_cell_type.clear()
        self.combo_cell_type.addItems(cell_types)
        if prev_ct in cell_types:
            self.combo_cell_type.setCurrentText(prev_ct)
        self.combo_cell_type.blockSignals(False)

        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        samples = []
        if md is not None and "sample_name" in md.columns:
            samples = sorted(md["sample_name"].astype(str).unique())

        prev_sample = self.combo_sample.currentText()
        self.combo_sample.blockSignals(True)
        self.combo_sample.clear()
        self.combo_sample.addItems(samples)
        if prev_sample in samples:
            self.combo_sample.setCurrentText(prev_sample)
        self.combo_sample.blockSignals(False)

        self._on_cell_type_changed()

    def _on_cell_type_changed(self, *_):
        self._teardown_preview()
        cell_type = self.combo_cell_type.currentText()
        self._load_features_for_cell_type(cell_type)
        self._populate_feature_combo()

    def _on_sample_changed(self, *_):
        self._teardown_preview()

    def _load_features_for_cell_type(self, cell_type: str):
        self._df_features = None
        self._df_features_csv_path = None
        self._df_features_cell_type = cell_type
        if not cell_type:
            return
        csv_path = self._feature_csv_path(cell_type)
        self._df_features_csv_path = csv_path
        if csv_path is None or not csv_path.exists():
            return
        try:
            self._df_features = pd.read_csv(csv_path, low_memory=False)
        except Exception:
            traceback.print_exc()
            self._df_features = None

    def _populate_feature_combo(self):
        prev = self.combo_feature.currentText()
        self.combo_feature.blockSignals(True)
        self.combo_feature.clear()
        numeric_cols = []
        if self._df_features is not None:
            numeric_cols = sorted(
                c for c in self._df_features.columns
                if c not in _FEATURE_BP_EXCLUDED_COLUMNS
                and pd.api.types.is_numeric_dtype(self._df_features[c])
            )
            self.combo_feature.addItems(numeric_cols)
            if prev in numeric_cols:
                self.combo_feature.setCurrentText(prev)
        self.combo_feature.blockSignals(False)

        has_features = len(numeric_cols) > 0
        self.combo_feature.setEnabled(has_features)
        self.btn_show.setEnabled(has_features and self.combo_sample.count() > 0)

        if has_features:
            self.status_label.setText("")
        elif not self.combo_cell_type.currentText():
            self.status_label.setText("Load metadata in the Data Preparation tab to see cell types.")
        elif self._df_features_csv_path is not None and not self._df_features_csv_path.exists():
            self.status_label.setText(
                f"No track-features CSV found for this cell type "
                f"(expected at {self._df_features_csv_path}). Run Feature Extraction first."
            )
        elif self._df_features is None:
            self.status_label.setText("Could not read the track-features CSV.")
        else:
            self.status_label.setText("No numeric feature columns found in the track-features CSV.")

    # ── preview lifecycle ────────────────────────────────────────────────
    def _current_viewer_frame(self) -> int:
        if self.viewer is None:
            return 0
        try:
            return int(self.viewer.dims.current_step[0])
        except Exception:
            return 0

    def _disconnect_dims_listener(self):
        if self.viewer is None or self._dims_callback is None:
            self._dims_callback = None
            return
        try:
            self.viewer.dims.events.current_step.disconnect(self._dims_callback)
        except Exception:
            pass
        self._dims_callback = None

    def _connect_dims_listener(self):
        self._disconnect_dims_listener()
        if self.viewer is None:
            return

        def _on_step(*_):
            self._refresh_feature_layer()

        try:
            self.viewer.dims.events.current_step.connect(_on_step)
            self._dims_callback = _on_step
        except Exception:
            self._dims_callback = None

    def _teardown_preview(self):
        self._disconnect_dims_listener()
        if self.viewer is not None:
            for layer in list(self.viewer.layers):
                if layer.name in self._added_layer_names:
                    try:
                        self.viewer.layers.remove(layer)
                    except Exception:
                        pass
        self._added_layer_names = []
        self._tracked_path = None

    def _on_show_clicked(self):
        self.status_label.setStyleSheet("color:#c66;font-size:11px;")
        self.status_label.setText("")

        if self.viewer is None:
            self.status_label.setText("No viewer available.")
            return

        cell_type = self.combo_cell_type.currentText()
        sample = self.combo_sample.currentText()
        feature_col = self.combo_feature.currentText()
        if not cell_type or not sample or not feature_col:
            self.status_label.setText("Select a cell type, sample, and feature first.")
            return

        out_dir = self._out_dir()
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        if out_dir is None or md is None or md.empty:
            self.status_label.setText("No metadata loaded.")
            return

        if self._df_features is None or self._df_features_cell_type != cell_type:
            self._load_features_for_cell_type(cell_type)
        if self._df_features is None:
            self.status_label.setText(
                f"No track-features CSV found for '{cell_type}'. Run Feature Extraction first."
            )
            return
        if feature_col not in self._df_features.columns:
            self.status_label.setText(
                f"Feature '{feature_col}' not found in {self._df_features_csv_path.name}."
            )
            return

        from behav3d.io.images import load_image

        sample_rows = md[md["sample_name"].astype(str) == str(sample)]
        row = sample_rows.iloc[0] if not sample_rows.empty else None

        try:
            raw_path = self._resolve_raw_path(row, out_dir, sample)
            if raw_path is None:
                self.status_label.setText(f"Raw image not found for sample '{sample}'.")
                return
            tracked_path = self._resolve_tracked_path(row, out_dir, sample, cell_type)
            if tracked_path is None:
                self.status_label.setText(
                    f"Tracked image not found for sample '{sample}', cell type '{cell_type}'."
                )
                return
        except Exception as exc:
            traceback.print_exc()
            self.status_label.setText(f"Could not resolve image paths: {exc}")
            return

        self._teardown_preview()

        try:
            raw_img = load_image(raw_path)
            try:
                dim_order = str(row.get("dimension_order", "TCZYX")).strip() or "TCZYX" if row is not None else "TCZYX"
            except Exception:
                dim_order = "TCZYX"
            ch_names = _bp_add_raw_channels(self.viewer, raw_img, sample, {}, dim_order)
            self._added_layer_names.extend(ch_names)

            tracked_dask = load_image(tracked_path)
            tracked_layer_name = f"{_FEATURE_LAYER_PREFIX} {cell_type} TrackID"
            self.viewer.add_labels(
                tracked_dask, name=tracked_layer_name, opacity=0.4, visible=False,
            )
            self._added_layer_names.append(tracked_layer_name)
        except Exception as exc:
            traceback.print_exc()
            self.status_label.setText(f"Could not load raw/tracked images: {exc}")
            return

        self._tracked_path = Path(tracked_path)
        self._sample = sample
        self._cell_type = cell_type
        self._feature_col = feature_col
        self._feature_layer_name = f"{_FEATURE_LAYER_PREFIX} {feature_col}"
        self._added_layer_names.append(self._feature_layer_name)

        self._refresh_feature_layer()
        self._connect_dims_listener()

    def _refresh_feature_layer(self):
        from behav3d.analysis.backprojection import backproject_feature_at_timepoint
        from behav3d.io.images import load_image_timepoint

        if self._tracked_path is None or self._df_features is None or self.viewer is None:
            return

        t = self._current_viewer_frame()
        try:
            labels_frame = np.asarray(load_image_timepoint(self._tracked_path, t))
        except Exception as exc:
            self.status_label.setStyleSheet("color:#c66;font-size:11px;")
            self.status_label.setText(f"Could not read frame {t}: {exc}")
            return
        if labels_frame.ndim == 4:
            labels_frame = labels_frame[0]

        df = self._df_features
        pos_mask = pd.to_numeric(df["position_t"], errors="coerce").fillna(-1).astype(int) == int(t)
        if "sample_name" in df.columns:
            sample_mask = df["sample_name"].astype(str) == str(self._sample)
        else:
            sample_mask = pd.Series(True, index=df.index)
        df_t = df[sample_mask & pos_mask]

        if df_t.empty:
            mapped = np.zeros_like(labels_frame)
            ids_with_value = np.array([], dtype=np.int64)
        else:
            try:
                mapped, ids_with_value = backproject_feature_at_timepoint(
                    labels_frame=labels_frame,
                    df_features=df_t,
                    feature_col=self._feature_col,
                    track_col="TrackID",
                    background_value=0,
                )
            except ValueError as exc:
                self.status_label.setStyleSheet("color:#c66;font-size:11px;")
                self.status_label.setText(f"Cannot backproject '{self._feature_col}': {exc}")
                return

        present_mask = (
            np.isin(labels_frame, ids_with_value) if ids_with_value.size > 0
            else np.zeros_like(labels_frame, dtype=bool)
        )
        vals_for_contrast = mapped[present_mask]
        if vals_for_contrast.size > 0:
            vmin, vmax = float(np.percentile(vals_for_contrast, 2)), float(np.percentile(vals_for_contrast, 98))
            if vmin == vmax:
                vmax = vmin + 1e-6
        else:
            vmin, vmax = 0.0, 1.0

        name = self._feature_layer_name
        try:
            # Layer already exists: update only the data. Deliberately leave
            # contrast_limits/opacity untouched so manual adjustments made in
            # napari's own layer controls survive timepoint scrubbing instead
            # of being reset by every recompute.
            layer = self.viewer.layers[name]
            layer.data = mapped
            layer.refresh()
        except (KeyError, ValueError):
            self.viewer.add_image(
                mapped, name=name, colormap="inferno",
                contrast_limits=(vmin, vmax), opacity=0.5, blending="translucent",
            )

        self.status_label.setStyleSheet("color:#6a6;font-size:11px;")
        self.status_label.setText(
            f"Showing '{self._feature_col}' for {self._cell_type} / {self._sample} at t={t}."
        )
