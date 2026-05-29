"""
BEHAV3D napari plugin – Feature Extraction Tab.

Provides per-cell-type sub-tabs with feature checkboxes (movement, intensity,
morphology, contact, death), thresholds, workers, and batch-run capability.

Dead threshold logic
--------------------
- **Organoids**: a single *global* threshold is shared across all organoid types.
  The threshold is set in the "Death Classification — Organoids" group at the top
  of this tab, and the same value is used for every organoid cell type.
- **Immune / other cell types**: each type has its *own* threshold spinner, set
  inside its individual sub-tab panel.
"""

import sys
import os
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QTabWidget, QTextEdit, QCheckBox,
    QDoubleSpinBox, QSpinBox, QGroupBox, QMessageBox, QScrollArea,
    QComboBox, QToolTip, QSplitter,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QCursor

from behav3d.core.qt_help import HelpButton, make_help_row
from behav3d.napari._results_panel import (
    ResultsPanel,
    notify_results_changed,
)
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    fire_extra_callback,
)

# Colormaps for raw channel layers (same order as the Visualization tab)
_CHANNEL_COLORS = ["cyan", "yellow", "green", "red", "blue", "magenta"]

# Prefix used for all temporary preview layers so they can be cleaned up easily
_PREVIEW_PREFIX = "[Preview]"


def _find_main_widget(start):
    """Walk up the parent chain to the top-level ``BEHAV3DWidget``."""
    w = start
    while w is not None and not hasattr(w, "tabs"):
        w = w.parent() if hasattr(w, "parent") else None
    return w


def _notify_post_extraction(source_widget):
    """Show a 'Run Filtering before Analysis' reminder after extraction.

    Honours a session-level opt-out flag stored on the top ``BEHAV3DWidget``
    (``_skip_filter_reminder``). The reminder is never raised during queued
    or non-interactive runs (those should call this only from the
    interactive paths).
    """
    main = _find_main_widget(source_widget)
    if main is None:
        return
    if getattr(main, "_skip_filter_reminder", False):
        return

    box = QMessageBox(source_widget)
    box.setWindowTitle("Feature Extraction Complete")
    box.setIcon(QMessageBox.Information)
    box.setText(
        "Feature extraction is complete.\n\n"
        "Remember to run Filtering before any Analysis step \u2014 "
        "analysis reads the filtered track-features CSV."
    )

    btn_run = box.addButton("Run Filtering with current setup", QMessageBox.AcceptRole)
    btn_goto = box.addButton("Go to Filtering tab", QMessageBox.ActionRole)
    btn_ok = box.addButton("OK", QMessageBox.RejectRole)
    box.setDefaultButton(btn_ok)

    optout = QCheckBox("Don't remind me this session")
    box.setCheckBox(optout)

    box.exec_()
    if optout.isChecked():
        try:
            main._skip_filter_reminder = True
        except Exception:
            pass

    clicked = box.clickedButton()
    if clicked is btn_run:
        filt = getattr(main, "filtering_tab", None)
        if filt is not None and hasattr(filt, "run_batch_filtering"):
            try:
                filt.run_batch_filtering(interactive=True)
            except Exception:
                traceback.print_exc()
    elif clicked is btn_goto:
        try:
            tabs = getattr(main, "tabs", None)
            filt = getattr(main, "filtering_tab", None)
            if tabs is not None and filt is not None:
                idx = tabs.indexOf(filt)
                if idx >= 0:
                    tabs.setCurrentIndex(idx)
        except Exception:
            traceback.print_exc()

# ═══════════════════════════════════════════════════════════════════════════
# Shared preview helpers (module-level so panels + tab can both use them)
# ═══════════════════════════════════════════════════════════════════════════

def _find_sample_with_segments(md: pd.DataFrame, cell_type: str):
    """Return (row, seg_path_str) for the first sample that has a valid
    segments file for *cell_type*, or (None, None) if none found."""
    for _, row in md.iterrows():
        for prefix in ("or", "im", "ot"):
            col = f"{prefix}_{cell_type}_segments_image_path"
            if col in row.index and pd.notna(row.get(col)):
                p = str(row[col]).strip()
                if p and Path(p).exists():
                    return row, p
    return None, None


def _resolve_dead_mask(sample_row: pd.Series, output_dir: Path, log_fn=None):
    """Try to find the dead mask array for a sample.

    Resolution order:
      1. ``dead_mask_path`` column in metadata (set by APOC segmentation).
      2. Several constructed paths under ``output_dir/images/{sample}/``.
      3. Raw dead channel -> Otsu threshold on first timepoint (last resort).

    Parameters
    ----------
    log_fn : callable | None
        Optional log function; each tried path is reported through it.

    Returns
    -------
    arr : dask.Array | None
    method : str | None   'zarr' | 'raw' | None
    tried : list[str]     All paths checked (for diagnostics).
    """
    import dask.array as da
    _log = log_fn or (lambda m: None)
    tried: list = []

    def _try_load(p):
        p = Path(p)
        tried.append(str(p))
        if not p.exists():
            return None
        try:
            from behav3d.io.images import load_image
            arr = load_image(p)
            if not isinstance(arr, da.Array):
                arr = da.from_array(np.asarray(arr))
            return arr
        except Exception as exc:
            _log(f"  Found {p.name} but failed to load: {exc}")
            return None

    # 1. PRIMARY: dead_mask_path column in metadata
    #    Written by ALL segmentation pipelines:
    #      - APOC          (apoc_segment.py)
    #      - Cellpose      (cellpose_prediction.py)
    #      - Pixel classifier (napari_pixelclassifier.py)
    dm_val = sample_row.get("dead_mask_path")
    if dm_val and pd.notna(dm_val):
        _log(f"  [dead mask] metadata dead_mask_path = {dm_val}")
        arr = _try_load(str(dm_val).strip())
        if arr is not None:
            return arr, "zarr", tried
        _log("  [dead mask] path from metadata does not exist on disk — trying fallbacks")
    else:
        _log(
            "  [dead mask] dead_mask_path column is empty or missing in metadata. "
            "Make sure segmentation has been run and the metadata CSV was saved afterwards."
        )

    # 2. FALLBACK: Constructed path using the canonical naming convention
    #    ({sample}_mask_dead.zarr is what all pipelines actually write)
    sample_name = str(sample_row.get("sample_name", ""))
    img_dir = output_dir / "images" / sample_name
    for name in (
        f"{sample_name}_mask_dead.zarr",   # APOC / pixelclassifier standard
        f"{sample_name}_dead_mask.zarr",   # Cellpose alternate naming
    ):
        arr = _try_load(img_dir / name)
        if arr is not None:
            _log(f"  [dead mask] found via fallback path: {img_dir / name}")
            return arr, "zarr", tried


    # 3. Raw dead channel -> Otsu (last resort; only first timepoint)
    raw_path = sample_row.get("raw_image_path")
    dead_ch_idx = sample_row.get("dead_channel")
    if raw_path and pd.notna(raw_path) and pd.notna(dead_ch_idx):
        label = f"raw_channel[{int(dead_ch_idx)}] via Otsu"
        tried.append(label)
        try:
            from behav3d.io.images import load_image
            raw = load_image(raw_path)
            ch = int(dead_ch_idx)
            dead_ch = raw[:, ch, ...] if raw.ndim >= 5 else raw[ch, ...]
            t0 = np.asarray(dead_ch[0] if dead_ch.ndim >= 4 else dead_ch)
            from skimage.filters import threshold_otsu
            try:
                otsu_thr = threshold_otsu(t0[t0 > 0]) if (t0 > 0).any() else float(t0.mean())
            except Exception:
                otsu_thr = float(t0.mean())
            dead_mask_t0 = (t0 > otsu_thr).astype(np.uint8)
            return da.from_array(dead_mask_t0[np.newaxis, ...]), "raw", tried
        except Exception:
            pass

    return None, None, tried


def _load_raw_dask(sample_row: pd.Series, output_dir: Path):
    """Load the raw image as a dask array (T, C, Z, Y, X).  Returns None on failure."""
    import dask.array as da

    sample_name = str(sample_row.get("sample_name", ""))
    raw_path_val = sample_row.get("raw_image_path", "")
    if not raw_path_val or pd.isna(raw_path_val):
        return None

    raw_p = Path(str(raw_path_val).strip())

    # Prefer zarr (lazy loading)
    zarr_candidates = [
        raw_p if raw_p.suffix == ".zarr" else None,
        output_dir / "images" / sample_name / f"{sample_name}.zarr",
        raw_p.with_suffix(".zarr"),
    ]
    for zp in zarr_candidates:
        if zp and zp.exists():
            try:
                from behav3d.io.formats.zarr import load_zarr
                return load_zarr(zp)
            except Exception:
                pass

    # Fallback: load full image
    if raw_p.exists():
        try:
            from behav3d.io.images import load_image
            img = load_image(raw_p)
            if not isinstance(img, da.Array):
                img = da.from_array(np.asarray(img), chunks=(1,) + np.asarray(img).shape[1:])
            return img
        except Exception:
            pass

    return None


def _load_all_segments_for_sample(
    sample_row: pd.Series,
    all_cell_types: list,
    output_dir: Path,
) -> "tuple[dict[str, np.ndarray], dict[str, str]]":
    """Load segment arrays for every cell type that has data for *sample_row*.

    Returns a dict mapping cell_type -> ndarray, typically:
      - 4-D (T, Z, Y, X) for timelapse labels, or
      - 3-D (Z, Y, X) for single-volume labels.
    Missing or unreadable types are silently omitted.

    Segment source priority is tracked-first:
      1) metadata tracked path
      2) constructed tracked path
      3) metadata untracked path
      4) constructed untracked path

    Returns
    -------
    segs : dict[str, np.ndarray]
        Loaded segment arrays by cell type.
    sources : dict[str, str]
        Source tag per loaded cell type: "tracked" or "untracked".
    """
    from behav3d.io.images import load_image as _li
    result = {}
    sources = {}
    sample_name = str(sample_row.get("sample_name", ""))

    def _normalize(arr):
        arr_np = np.asarray(arr)
        if arr_np.ndim == 5:
            # Defensive: collapse channel axis if present.
            arr_np = arr_np[:, 0, ...]
        return arr_np

    def _try_load(p: Path):
        if not p.exists():
            return None
        try:
            return _normalize(_li(str(p)))
        except Exception:
            return None

    for ct in all_cell_types:
        # 1. Metadata tracked path (preferred)
        tracked_cols = [
            f"{pfx}_{ct}_tracks_image_path" for pfx in ("or", "im", "ot")
        ] + [f"{ct}_tracks_image_path"]
        for col in tracked_cols:
            val = sample_row.get(col)
            if val and pd.notna(val):
                arr_np = _try_load(Path(str(val).strip()))
                if arr_np is not None:
                    result[ct] = arr_np
                    sources[ct] = "tracked"
                    break
        if ct in result:
            continue

        # 2. Constructed tracked fallback
        for suffix in (
            f"{sample_name}_{ct}_tracked.zarr",
            f"{ct}_tracked.zarr",
        ):
            p = output_dir / "images" / sample_name / suffix
            arr_np = _try_load(p)
            if arr_np is not None:
                result[ct] = arr_np
                sources[ct] = "tracked"
                break
        if ct in result:
            continue

        # 3. Metadata untracked path
        for pfx in ("or", "im", "ot"):
            col = f"{pfx}_{ct}_segments_image_path"
            val = sample_row.get(col)
            if val and pd.notna(val):
                arr_np = _try_load(Path(str(val).strip()))
                if arr_np is not None:
                    result[ct] = arr_np
                    sources[ct] = "untracked"
                    break
        if ct in result:
            continue

        # 4. Constructed untracked fallback
        for suffix in (
            f"{sample_name}_{ct}_segments.zarr",
            f"{ct}_segments.zarr",
        ):
            p = output_dir / "images" / sample_name / suffix
            arr_np = _try_load(p)
            if arr_np is not None:
                result[ct] = arr_np
                sources[ct] = "untracked"
                break
    return result, sources


def _add_channel_layers(viewer, dask_img, sample_name: str, prefix: str = _PREVIEW_PREFIX):
    """Split a (T, C, Z, Y, X) dask array along C and add each channel as a
    napari Image layer with additive blending."""
    if dask_img is None:
        return
    n_channels = dask_img.shape[1] if dask_img.ndim >= 5 else 1
    for c in range(n_channels):
        ch_data = dask_img[:, c, ...] if dask_img.ndim >= 5 else dask_img
        color = _CHANNEL_COLORS[c % len(_CHANNEL_COLORS)]
        layer = viewer.add_image(
            ch_data,
            name=f"{prefix} {sample_name} – Ch{c}",
            colormap=color,
            blending="additive",
            opacity=0.7,
            visible=True,
        )
        # Auto-contrast on non-zero pixels
        try:
            flat = np.asarray(ch_data[0]).ravel()
            flat = flat[flat > 0]
            if flat.size > 0:
                layer.contrast_limits = (0, float(np.percentile(flat, 99.8)))
        except Exception:
            pass


def _overlay_for_volume(seg_vol, dead_vol, threshold_pct, frame_label="", log_fn=None):
    """Single-volume primitive: compute Dead/Alive overlay + per-region stats.

    Returns
    -------
    overlay_vol : np.ndarray (uint8)
        Same shape as ``seg_vol``. 0 = bg, 1 = alive, 2 = dead.
    region_stats : list[dict]
        One entry per labelled segment with keys ``label``, ``centroid``,
        ``pct_dead``, ``extent_x``.
    """
    import datetime
    from skimage.measure import regionprops

    _log = log_fn or (lambda m: None)

    def _ts():
        return datetime.datetime.now().strftime("%H:%M:%S")

    overlay_vol = np.zeros_like(seg_vol, dtype=np.uint8)
    dead_binary = dead_vol > 0
    region_stats = []

    seg_int = seg_vol.astype(np.int32)
    ndim = seg_int.ndim
    regions = regionprops(seg_int)
    n_regions = len(regions)
    if frame_label:
        msg = f"[{_ts()}] [Preview] Computing overlay{frame_label}: {n_regions} segments…"
        print(msg)
        _log(f"  ⏳ Overlay{frame_label}: {n_regions} segments…")
    for region in regions:
        label_val = region.label
        mask = seg_int == label_val
        total_pixels = int(mask.sum())
        if total_pixels == 0:
            continue
        dead_pixels = int((mask & dead_binary).sum())
        pct = (dead_pixels / total_pixels) * 100.0
        overlay_vol[mask] = 2 if (pct >= threshold_pct and threshold_pct > 0) else 1

        bbox_min = region.bbox[:ndim]
        bbox_max = region.bbox[ndim:]
        extent_x = float(bbox_max[-1] - bbox_min[-1])
        region_stats.append(
            {
                "label": int(label_val),
                "centroid": tuple(float(v) for v in region.centroid),
                "pct_dead": float(pct),
                "extent_x": extent_x,
            }
        )
    return overlay_vol, region_stats


def _dead_alive_colormap():
    """Return a DirectLabelColormap for the Dead/Alive overlay, or ``None``."""
    from napari.utils.colormaps import DirectLabelColormap
    try:
        return DirectLabelColormap(
            color_dict={
                None: [0, 0, 0, 0],      # default for any unlisted label
                0:    [0, 0, 0, 0],      # background → transparent
                1:    [0, 0.8, 0, 0.6],  # alive → green
                2:    [0.9, 0, 0, 0.6],  # dead → red
            }
        )
    except Exception:
        return None


def _update_dead_alive_overlay(
    viewer,
    seg_data: np.ndarray,
    dead_data: np.ndarray,
    threshold_pct: float,
    *,
    frame_idx: int | None = None,
    overlay_arr: np.ndarray | None = None,
    layer_name: str = f"{_PREVIEW_PREFIX} Dead/Alive",
    log_fn=None,
):
    """Compute the Dead/Alive overlay and update (or create) the labels layer.

    On-demand single-frame mode
    ---------------------------
    For timelapse data (``seg_data.ndim >= 4``) only the timepoint indexed by
    ``frame_idx`` is computed. ``overlay_arr`` is an in/out full ``T,Z,Y,X``
    (uint8) array that is allocated lazily on the first call and re-used
    across calls. Only the slice ``overlay_arr[frame_idx]`` is written; other
    frames remain whatever the caller previously set them to (typically zero
    = transparent until that frame is visited).

    For non-timelapse data (3D), behavior is the historical "single shot":
    the whole volume is computed in one go.

    Labels
    ------
    0  transparent (background)
    1  green  (alive)
    2  red    (dead)

    Returns
    -------
    overlay_arr : np.ndarray
        The full overlay array passed in (timelapse) or freshly computed (3D).
    frame_stats : list[dict]
        For timelapse mode: region stats for the single computed frame.
        For 3D mode: region stats for the whole volume.
    """
    import datetime

    def _ts():
        return datetime.datetime.now().strftime("%H:%M:%S")

    seg_arr = np.asarray(seg_data)
    dead_arr = np.asarray(dead_data)

    is_timelapse = seg_arr.ndim >= 4

    if is_timelapse:
        if dead_arr.ndim == 3:
            dead_arr = np.repeat(dead_arr[np.newaxis, ...], seg_arr.shape[0], axis=0)
        t_common = min(seg_arr.shape[0], dead_arr.shape[0])

        if frame_idx is None:
            frame_idx = 0
        frame_idx = int(frame_idx)
        if frame_idx < 0:
            frame_idx = 0
        if frame_idx >= t_common:
            frame_idx = t_common - 1

        overlay_shape = tuple(seg_arr[:t_common].shape)
        if overlay_arr is None:
            overlay_arr = np.zeros(overlay_shape, dtype=np.uint8)
        elif overlay_arr.shape != overlay_shape:
            overlay_arr = np.zeros(overlay_shape, dtype=np.uint8)

        seg_vol = np.asarray(seg_arr[frame_idx])
        dead_vol = np.asarray(dead_arr[frame_idx])
        if seg_vol.ndim > 3:
            seg_vol = seg_vol[0]
        if dead_vol.ndim > 3:
            dead_vol = dead_vol[0]
        if dead_vol.ndim < seg_vol.ndim:
            dead_vol = np.broadcast_to(dead_vol, seg_vol.shape)

        print(
            f"[{_ts()}] [Preview] Dead/Alive overlay (on-demand): "
            f"frame={frame_idx}/{t_common - 1}, threshold={threshold_pct:.2f}%"
        )
        overlay_t, stats_t = _overlay_for_volume(
            seg_vol, dead_vol, threshold_pct,
            frame_label=f" t={frame_idx}", log_fn=log_fn,
        )
        overlay_arr[frame_idx] = overlay_t
        n_dead = sum(
            1 for r in stats_t
            if r["pct_dead"] >= threshold_pct and threshold_pct > 0
        )
        print(
            f"[{_ts()}] [Preview]   t={frame_idx} done — "
            f"{len(stats_t)} cells, {n_dead} dead"
        )
        overlay = overlay_arr
        frame_stats = stats_t
    else:
        seg_vol = np.asarray(seg_arr)
        dead_vol = np.asarray(dead_arr)
        if dead_vol.ndim > seg_vol.ndim:
            dead_vol = dead_vol.max(axis=0)
        if dead_vol.ndim < seg_vol.ndim:
            dead_vol = np.broadcast_to(dead_vol, seg_vol.shape)
        print(
            f"[{_ts()}] [Preview] Dead/Alive overlay: single volume, "
            f"threshold={threshold_pct:.2f}%"
        )
        overlay, stats_0 = _overlay_for_volume(
            seg_vol, dead_vol, threshold_pct, frame_label="", log_fn=log_fn,
        )
        n_dead = sum(
            1 for r in stats_0
            if r["pct_dead"] >= threshold_pct and threshold_pct > 0
        )
        print(f"[{_ts()}] [Preview]   Done — {len(stats_0)} cells, {n_dead} dead")
        frame_stats = stats_0

    cmap = _dead_alive_colormap()

    print(f"[{_ts()}] [Preview] Updating Dead/Alive layer in viewer…")
    try:
        layer = viewer.layers[layer_name]
        # In timelapse mode the layer already references ``overlay_arr`` (same
        # object), so the in-place slice write above is enough; just refresh.
        # In 3D mode we replace ``layer.data`` with the freshly computed array.
        if not is_timelapse or layer.data is not overlay:
            layer.data = overlay
        layer.opacity = 1.0
        layer.refresh()
    except (KeyError, ValueError):
        kw = dict(name=layer_name, opacity=1.0)
        if cmap is not None:
            kw["colormap"] = cmap
        viewer.add_labels(overlay, **kw)

    print(f"[{_ts()}] [Preview] Dead/Alive layer ready.")
    return overlay, frame_stats


def _remove_layer_if_exists(viewer, layer_name: str):
    try:
        layer = viewer.layers[layer_name]
        viewer.layers.remove(layer)
    except Exception:
        pass



def _build_dead_pct_map(region_stats):
    """Build a {label_id: pct_dead} map from a single frame's region stats."""
    frame_map: dict = {}
    for r in region_stats or []:
        label_id = int(r.get("label", 0))
        if label_id > 0:
            frame_map[label_id] = float(r.get("pct_dead", 0.0))
    return frame_map


def _merge_org_segments(
    segs_dict: dict, org_cell_types: list
) -> "tuple[np.ndarray | None, dict]":
    """Merge segment arrays for multiple organoid types into one array with
    non-overlapping labels so the Dead/Alive overlay can cover all types.

    Returns
    -------
    merged_arr : ndarray | None
    label_type_map : dict
        merged_label_id -> (cell_type, original_label_id)
    """
    available = [ct for ct in org_cell_types if ct in segs_dict]
    if not available:
        return None, {}
    if len(available) == 1:
        ct = available[0]
        arr = np.asarray(segs_dict[ct])
        label_map = {
            int(lbl): (ct, int(lbl))
            for lbl in np.unique(arr)
            if lbl > 0
        }
        return arr.copy(), label_map

    ref = np.asarray(segs_dict[available[0]])
    merged = np.zeros(ref.shape, dtype=np.int32)
    label_map: dict = {}
    offset = 0
    for ct in available:
        arr = np.asarray(segs_dict[ct]).astype(np.int32)
        if arr.shape != merged.shape:
            continue
        max_lbl = int(arr.max()) if arr.size > 0 else 0
        for lbl in np.unique(arr):
            if lbl <= 0:
                continue
            label_map[int(lbl + offset)] = (ct, int(lbl))
        shifted = np.where(arr > 0, arr + offset, 0)
        merged = np.where(shifted > 0, shifted, merged)
        offset += max_lbl + 1
    return merged, label_map


# ═══════════════════════════════════════════════════════════════════════════
# Per-cell-type panel
# ═══════════════════════════════════════════════════════════════════════════
class CellTypeFeaturePanel(QWidget):
    """Feature extraction controls for one cell type.

    Parameters
    ----------
    is_organoid : bool
        When *True* the dead threshold is controlled by the global spinner in
        ``FeatureExtractionTab``; no per-panel spinner is shown.
        When *False* a per-panel threshold spinner + preview button are shown.
    threshold_getter : callable | None
        Callable that returns the effective dead threshold value (float).
        Required when ``is_organoid=True``; otherwise the panel uses its own
        ``spin_dead_threshold`` widget.
    """

    def __init__(
        self,
        cell_type: str,
        category: str,
        metadata_loader,
        all_cell_types: list,
        category_types: list,
        log_callback=None,
        is_organoid: bool = False,
        threshold_getter=None,
        parent=None,
        tab_progress_row=None,
    ):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types
        self.log = log_callback or (lambda m: None)
        self.viewer = None
        self._is_organoid = is_organoid
        self._threshold_getter = threshold_getter  # kept for backward compat
        self._preview_connected = False
        # Background-execution infrastructure.
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        # List of organoid cell types (set by FeatureExtractionTab for org panels)
        self._org_cell_types: list = []

        # Shared cache dict (for organoid panels, set by FeatureExtractionTab)
        self._org_preview_cache: dict | None = None

        # Cached arrays for live overlay updates
        self._preview_seg_t = None
        self._preview_dead_t = None
        # merged_label_id -> (cell_type, original_label_id)  — organoid panels
        self._preview_label_type_map: dict = {}
        self._preview_hover_callback = None

        # On-demand per-timepoint caches (populated lazily as the user
        # navigates the napari time slider or changes the threshold).
        # ``_preview_overlay_arr`` is the full T,Z,Y,X (or 3D) uint8 array
        # that backs the Dead/Alive labels layer — slices are filled in
        # place as frames are visited.
        self._preview_overlay_arr = None
        self._preview_pct_maps_by_frame: dict = {}
        self._preview_stats_by_frame: dict = {}
        self._preview_computed_frames: set = set()
        self._preview_current_thr: float | None = None
        self._preview_dims_callback = None
        self._preview_is_timelapse: bool = False

        # Read saved config
        params = self.metadata_loader.behav3d_parameters
        fcfg = params.get("features", {}).get(self.cell_type, {}) or {}

        # Detect dead channel presence
        md = metadata_loader.metadata
        self._has_dead = False
        if md is not None:
            self._has_dead = (
                "dead_channel" in md.columns and md["dead_channel"].notna().any()
            )

        # Will be created inside ``_init_ui`` if a dead channel is configured.
        self.btn_rerun_death = None

        self._init_ui(fcfg)
        self._refresh_rerun_death_button()

    # ── UI ──────────────────────────────────────────────────────────────────
    def _init_ui(self, fcfg):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # ── Feature checkboxes ────────────────────────────────────────────
        feat_group = QGroupBox("Features to extract")
        feat_lay = QVBoxLayout(feat_group)
        feat_lay.setSpacing(2)

        all_feats = ["movement", "intensity", "morphology", "contact", "invasiveness", "death"]

        # ── Mandatory features (same logic as the notebook widget) ─────────
        # intensity + contact are mandatory for all types.
        # movement is mandatory for immune / other cells.
        # death is mandatory whenever a dead channel is present.
        _mandatory: set = {"intensity", "contact"}
        if self.category in {"immune", "other"}:
            _mandatory.add("movement")
        if self._has_dead:
            _mandatory.add("death")
        self._mandatory_features: set = _mandatory

        default_feats = fcfg.get("features_choice", all_feats)
        if not isinstance(default_feats, (list, tuple)):
            default_feats = all_feats
        # Ensure mandatory features are always present in the saved config
        default_feats = list(default_feats)
        for _mf in self._mandatory_features:
            if _mf not in default_feats:
                default_feats.append(_mf)

        self.feature_checks: dict[str, QCheckBox] = {}
        for f in all_feats:
            if f == "death" and not self._has_dead:
                continue
            if f == "invasiveness" and self.category != "immune":
                # Invasiveness only makes sense for immune cells (counts immune
                # voxels inside organoid masks). Mirrors the notebook widget.
                continue
            label = (
                "Organoid Invasiveness (Advanced)"
                if f == "invasiveness" else f.capitalize()
            )
            cb = QCheckBox(label)
            cb.setChecked(f in default_feats)
            self.feature_checks[f] = cb

        # Force-check and disable mandatory checkboxes (cannot be toggled off)
        for _mf in self._mandatory_features:
            if _mf in self.feature_checks:
                self.feature_checks[_mf].setChecked(True)
                self.feature_checks[_mf].setEnabled(False)

        # Explanatory note shown at the top of the feature group
        _mandatory_note = QLabel(
            "ℹ  Greyed-out features are always extracted for this cell type."
        )
        _mandatory_note.setWordWrap(True)
        _mandatory_note.setStyleSheet(
            "color: #90A4AE; font-style: italic; font-size: 10px;"
        )
        feat_lay.addWidget(_mandatory_note)

        # Contact threshold (enabled only when Contact is checked)
        self.contact_threshold = QDoubleSpinBox()
        self.contact_threshold.setRange(0.0, 500.0)
        self.contact_threshold.setSingleStep(0.5)
        self.contact_threshold.setDecimals(2)
        self.contact_threshold.setValue(float(fcfg.get("contact_threshold", 0.0)))
        self.contact_threshold.setMaximumWidth(90)
        self.contact_threshold.setSuffix(" µm")

        contact_row = QHBoxLayout()
        contact_row.setSpacing(6)
        contact_row.addWidget(self.feature_checks.get("contact", QCheckBox()))
        contact_help = make_help_row(
            self.contact_threshold,
            "Contact Threshold (µm)",
            "Distance threshold (in micrometers) for detecting intensity transfer\n"
            "between neighbouring cells (contact features).\n\n"
            "Segments closer than this will be checked for contact signal.",
        )
        contact_row.addLayout(contact_help)

        # Group the optional 'Organoid Invasiveness' checkbox inline with the
        # contact row (immune cells only, mirroring the notebook widget).
        if self.category == "immune" and "invasiveness" in self.feature_checks:
            inv_cb = self.feature_checks["invasiveness"]
            inv_cb.setToolTip(
                "Count how many immune voxels lie inside organoid masks at each "
                "timepoint. Requires both immune segments and organoid masks; "
                "automatically aggregated downstream."
            )
            contact_row.addWidget(inv_cb)
            contact_row.addWidget(HelpButton(
                "Organoid Invasiveness (Advanced)",
                "Per-timepoint count of voxels of this immune cell type that fall "
                "inside any organoid mask.\n\n"
                "Helps quantify how aggressively immune cells penetrate organoids. "
                "Only available when Contact is enabled."
            ))
        contact_row.addStretch()

        def _toggle_ct(state=None):
            contact_on = self.feature_checks.get("contact", QCheckBox()).isChecked()
            self.contact_threshold.setEnabled(contact_on)
            # Invasiveness depends on contact (notebook parity).
            if self.category == "immune" and "invasiveness" in self.feature_checks:
                self.feature_checks["invasiveness"].setVisible(contact_on)

        _toggle_ct()
        if "contact" in self.feature_checks:
            self.feature_checks["contact"].stateChanged.connect(_toggle_ct)

        for f in all_feats:
            if f == "contact":
                feat_lay.addLayout(contact_row)
            elif f == "invasiveness":
                # Already rendered inline within contact_row (or skipped entirely
                # for non-immune cell types above).
                continue
            elif f in self.feature_checks:
                feat_lay.addWidget(self.feature_checks[f])

        layout.addWidget(feat_group)

        # ── Death Threshold ───────────────────────────────────────────────

        dead_group = QGroupBox("Death Threshold")
        dead_lay = QVBoxLayout(dead_group)
        dead_lay.setSpacing(4)

        # -- Death threshold spinner (ALL panels own one; organoid panels sync)
        dead_desc = QLabel(
            f"Percentage of dead-mask pixels overlapping a segment\n"
            f"above which {self.cell_type} cells are classified as 'dead'."
        )
        dead_desc.setWordWrap(True)
        dead_desc.setStyleSheet("color: #888; font-size: 10px;")
        dead_lay.addWidget(dead_desc)

        if self._is_organoid:
            shared_note = QLabel(
                "\u26a0\ufe0f  One threshold applies to ALL organoid types equally.\n"
                "Adjusting it here updates every other organoid tab in real time."
            )
            shared_note.setWordWrap(True)
            shared_note.setStyleSheet(
                "color: #FFB74D; font-style: italic; font-size: 10px;"
            )
            dead_lay.addWidget(shared_note)

        dead_form = QFormLayout()
        dead_form.setContentsMargins(0, 0, 0, 0)
        dead_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.spin_dead_threshold = QDoubleSpinBox()
        self.spin_dead_threshold.setRange(0.0, 100.0)
        self.spin_dead_threshold.setSingleStep(0.1)
        self.spin_dead_threshold.setDecimals(3)
        saved_thr = fcfg.get("dead_mask_percentage_threshold", 0.1)
        self.spin_dead_threshold.setValue(float(saved_thr) if saved_thr else 0.1)
        self.spin_dead_threshold.setSuffix(" %")
        self.spin_dead_threshold.setMaximumWidth(100)

        # Track the threshold that is currently baked into the combined CSV
        # on disk. ``None`` means feature extraction has never run for this
        # cell type (or no death classification was performed). Updated by
        # ``_persist()`` and by ``_on_rerun_death_clicked``.
        self._last_persisted_threshold = (
            float(saved_thr) if saved_thr not in (None, 0, 0.0) else None
        )

        # "Re-run death features" button — always visible when a combined
        # feature CSV exists; dimmed when the spinner value matches the
        # persisted threshold. Placed to the right of the spinner+help row.
        self.btn_rerun_death = QPushButton("\U0001F504  Re-run death")
        self.btn_rerun_death.setMaximumWidth(150)
        self.btn_rerun_death.setStyleSheet(
            "QPushButton { background: #1976D2; color: white; padding: 4px 8px; "
            "border-radius: 3px; font-size: 11px; font-weight: bold; } "
            "QPushButton:disabled { background: #455A64; color: #B0BEC5; } "
            "QPushButton:hover:!disabled { background: #1565C0; }"
        )
        self.btn_rerun_death.clicked.connect(self._on_rerun_death_clicked)

        spinner_row = make_help_row(
            self.spin_dead_threshold,
            "Dead Mask Percentage Threshold",
            "Percentage of dead-mask pixels overlapping a segment's volume\n"
            "required to classify the cell as dead.\n\n"
            "Set to 0 to skip dead classification.\n"
            "Typical range: 0.05\u20130.5 %.",
        )
        spinner_row.addWidget(self.btn_rerun_death)
        dead_form.addRow("Dead mask % threshold:", spinner_row)
        dead_lay.addLayout(dead_form)

        # Keep button state in sync with the spinner value.
        self.spin_dead_threshold.valueChanged.connect(
            lambda _v: self._refresh_rerun_death_button()
        )

        # Per-panel sample selector
        prev_sample_row_layout = QHBoxLayout()
        prev_sample_row_layout.addWidget(QLabel("Preview sample:"))
        self.preview_sample_combo = QComboBox()
        self.preview_sample_combo.setMinimumWidth(180)
        self.preview_sample_combo.setToolTip(
            "Select which sample to load for the dead threshold preview."
        )
        prev_sample_row_layout.addWidget(self.preview_sample_combo, stretch=1)
        dead_lay.addLayout(prev_sample_row_layout)

        # Populate combo
        _md = self.metadata_loader.metadata if self.metadata_loader else None
        if _md is not None and not _md.empty:
            _samples = sorted(str(s) for s in _md["sample_name"].unique())
            self.preview_sample_combo.addItems(_samples)

        btn_preview = QPushButton("\U0001F441  Preview Dead Threshold in Viewer")
        btn_preview.setStyleSheet(
            "QPushButton { background: #37474F; color: white; padding: 5px 10px; "
            "border-radius: 3px; font-size: 11px; } "
            "QPushButton:hover { background: #546E7A; }"
        )
        btn_preview.setToolTip(
            f"Load raw channels + dead mask for {self.cell_type} and show\n"
            "a green (alive) / red (dead) per-cell overlay.\n"
            "Adjust the threshold above to see live changes."
        )
        btn_preview.clicked.connect(self._on_preview_dead_clicked)
        dead_lay.addWidget(btn_preview)

        # Organoid spinners notify parent tab -> all org panels stay in sync
        if self._is_organoid:
            def _org_sync(value, _ct=self.cell_type):
                pt = self.parent()
                while pt and not hasattr(pt, '_notify_organoid_threshold_changed'):
                    pt = pt.parent()
                if pt:
                    pt._notify_organoid_threshold_changed(_ct, value)
            self.spin_dead_threshold.valueChanged.connect(_org_sync)
        dead_group.setVisible(self._has_dead)
        layout.addWidget(dead_group)
        self._dead_group = dead_group

        def _toggle_dead_group(state=None):
            cb = self.feature_checks.get("death")
            self._dead_group.setVisible(
                self._has_dead and (cb is not None and cb.isChecked())
            )

        if "death" in self.feature_checks:
            self.feature_checks["death"].stateChanged.connect(_toggle_dead_group)
        _toggle_dead_group()

        # ── Workers ───────────────────────────────────────────────────────
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, max_allowed)
        default_workers = min(int(fcfg.get("n_workers", max(4, n_cores // 2))), max_allowed)
        self.spin_workers.setValue(default_workers)
        self.spin_workers.setMaximumWidth(60)
        self.spin_workers.valueChanged.connect(self._on_workers_changed)

        workers_form = QFormLayout()
        workers_form.setContentsMargins(0, 0, 0, 0)
        workers_form.addRow(
            "Workers:",
            make_help_row(
                self.spin_workers,
                "Number of Workers",
                f"Number of CPU cores to use for parallel processing.\n\n"
                f"Your machine has {n_cores} cores.\n"
                f"Recommendation: Use at most {max(1, n_cores - 1)} cores to keep the system responsive.",
            ),
        )
        layout.addLayout(workers_form)

        # ── Sync-to-others buttons ────────────────────────────────────────
        cat_label = (
            self.category.capitalize() + "s" if self.category != "other" else "Other types"
        )
        sync_row = QHBoxLayout()
        btn_apply_cat = QPushButton(f"Apply to all {cat_label}")
        btn_apply_cat.clicked.connect(lambda: self._apply_to_others(category_only=True))
        btn_apply_all = QPushButton("Apply to all")
        btn_apply_all.clicked.connect(lambda: self._apply_to_others(category_only=False))
        sync_row.addWidget(btn_apply_cat)
        sync_row.addWidget(btn_apply_all)
        layout.addLayout(sync_row)

        # ── Run button ────────────────────────────────────────────────────
        self.btn_run = QPushButton(
            f"Run {self.cell_type.capitalize()} Feature Extraction"
        )
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        # ``clicked`` emits a ``bool`` (checked state) that would override the
        # ``interactive=True`` default of ``_on_run_clicked``; wrap in a lambda
        # so the overwrite prompt is not silently skipped on button presses.
        self.btn_run.clicked.connect(lambda: self._on_run_clicked(interactive=True))
        layout.addWidget(self.btn_run)

        layout.addStretch()

    # ── Helpers ──────────────────────────────────────────────────────────────
    def _get_threshold(self) -> float:
        """Return the effective dead-mask percentage threshold for this panel."""
        if self.spin_dead_threshold is not None:
            return round(float(self.spin_dead_threshold.value()), 2)
        if self._threshold_getter is not None:
            return round(float(self._threshold_getter()), 2)
        return 0.0

    def _disconnect_preview_dead_hover(self):
        if self.viewer is None or self._preview_hover_callback is None:
            self._preview_hover_callback = None
            return
        try:
            self.viewer.mouse_move_callbacks.remove(self._preview_hover_callback)
        except Exception:
            pass
        self._preview_hover_callback = None

    def _disconnect_preview_dims(self):
        """Disconnect the napari time-slider listener if connected."""
        viewer = self.viewer
        if viewer is None or self._preview_dims_callback is None:
            self._preview_dims_callback = None
            return
        try:
            viewer.dims.events.current_step.disconnect(self._preview_dims_callback)
        except Exception:
            pass
        self._preview_dims_callback = None

    def _connect_preview_dims(self):
        """Connect a viewer.dims.current_step listener that triggers an
        on-demand recomputation for whatever frame the slider moved to.

        Safe to call multiple times — any previous listener is dropped first.
        """
        viewer = self.viewer
        if viewer is None:
            return
        self._disconnect_preview_dims()

        def _on_step(*_):
            # Only act if we still have valid preview state
            if (
                self._preview_seg_t is None
                or self._preview_dead_t is None
                or not self._preview_is_timelapse
            ):
                return
            self._refresh_preview_dead_layers(value=None)

        try:
            viewer.dims.events.current_step.connect(_on_step)
            self._preview_dims_callback = _on_step
        except Exception:
            self._preview_dims_callback = None

    def _attach_preview_dead_hover(self, layer_name: str = f"{_PREVIEW_PREFIX} Dead/Alive"):
        viewer = self.viewer
        self._disconnect_preview_dead_hover()
        if (
            viewer is None
            or self._preview_seg_t is None
            or not self._preview_pct_maps_by_frame
        ):
            return
        try:
            viewer.layers[layer_name]
        except Exception:
            return

        seg_data = self._preview_seg_t
        seg_layer_name = f"{_PREVIEW_PREFIX} {self.cell_type} segments"

        def _segment_label_at(position, frame_idx: int) -> int:
            try:
                seg_arr = np.asarray(seg_data)
                if seg_arr.ndim >= 4:
                    if frame_idx < 0 or frame_idx >= seg_arr.shape[0]:
                        return 0
                    seg_vol = np.asarray(seg_arr[frame_idx])
                else:
                    seg_vol = np.asarray(seg_arr)
                if seg_vol.ndim > 3:
                    seg_vol = seg_vol[0]

                try:
                    seg_layer = self.viewer.layers[seg_layer_name]
                except Exception:
                    return 0
                data_position = seg_layer.world_to_data(position)
                coords = tuple(int(round(float(c))) for c in data_position[-seg_vol.ndim:])
                for i, c in enumerate(coords):
                    if c < 0 or c >= seg_vol.shape[i]:
                        return 0
                return int(seg_vol[coords])
            except Exception:
                return 0

        def _show_tooltip(text: str):
            try:
                QToolTip.showText(QCursor.pos(), text, self.viewer.window._qt_window)
            except Exception:
                try:
                    QToolTip.showText(QCursor.pos(), text)
                except Exception:
                    pass

        def _hide_tooltip():
            try:
                QToolTip.hideText()
            except Exception:
                pass

        # Capture label-type map so the closure always uses the current map
        _label_type_map = getattr(self, "_preview_label_type_map", {})

        def _on_mouse_move(*args):
            """Works with both napari ≤0.4.18 (event,) and ≥0.4.19 (viewer,event)."""
            if self.viewer is None:
                return
            position = getattr(self.viewer.cursor, "position", None)
            if position is None:
                self.viewer.status = ""
                _hide_tooltip()
                return

            frame_idx = self._current_viewer_frame() if self._preview_is_timelapse else 0

            label_id = _segment_label_at(position, frame_idx)
            if label_id <= 0:
                self.viewer.status = ""
                _hide_tooltip()
                return

            # Look up the per-frame pct map (may be missing if the user is
            # hovering before the on-demand recompute for this frame has
            # finished — in that case fall back to the label only).
            frame_map = self._preview_pct_maps_by_frame.get(frame_idx)
            pct_dead = frame_map.get(label_id) if frame_map else None
            # Resolve cell type + original label from merged organoid map
            type_info = _label_type_map.get(label_id)
            if type_info is not None:
                ct_name, orig_lbl = type_info
                if pct_dead is None:
                    status_text = f"{ct_name} #{orig_lbl}"
                else:
                    status_text = f"{ct_name} #{orig_lbl} | {pct_dead:.2f}% dead"
            else:
                if pct_dead is None:
                    status_text = f"Cell #{label_id}"
                else:
                    status_text = f"Cell #{label_id} | {pct_dead:.2f}% dead"
            self.viewer.status = status_text
            _show_tooltip(status_text)

        self._preview_hover_callback = _on_mouse_move
        # viewer.mouse_move_callbacks is the napari-version-agnostic API
        viewer.mouse_move_callbacks.append(_on_mouse_move)

    def _current_viewer_frame(self) -> int:
        """Return the napari viewer's current time-axis index (0 if no time)."""
        viewer = self.viewer
        if viewer is None:
            return 0
        try:
            return int(viewer.dims.current_step[0])
        except Exception:
            return 0

    def _invalidate_preview_cache(self):
        """Drop all cached per-frame overlays/stats (e.g. on threshold change)."""
        self._preview_computed_frames.clear()
        self._preview_pct_maps_by_frame.clear()
        self._preview_stats_by_frame.clear()
        if self._preview_overlay_arr is not None:
            # Zero so previously-visited frames don't render with stale colors
            # until they are recomputed at the new threshold.
            self._preview_overlay_arr.fill(0)

    def _refresh_preview_dead_layers(
        self,
        value: float | None = None,
        frame_idx: int | None = None,
    ):
        """Compute (or refresh) the Dead/Alive overlay for a single timepoint.

        Parameters
        ----------
        value : float | None
            New threshold value. If ``None``, uses the panel's spinner value.
        frame_idx : int | None
            Timepoint to compute. If ``None``, uses the viewer's current step.
        """
        viewer = self.viewer
        if (
            viewer is None
            or self._preview_seg_t is None
            or self._preview_dead_t is None
        ):
            return

        thr = self._get_threshold() if value is None else round(float(value), 2)

        # Threshold changed -> invalidate every cached frame
        if self._preview_current_thr is None or thr != self._preview_current_thr:
            self._invalidate_preview_cache()
            self._preview_current_thr = thr

        if frame_idx is None:
            frame_idx = self._current_viewer_frame() if self._preview_is_timelapse else 0
        frame_idx = int(frame_idx)

        # Fast path: this frame is already computed at the current threshold
        if frame_idx in self._preview_computed_frames:
            try:
                layer = viewer.layers[f"{_PREVIEW_PREFIX} Dead/Alive"]
                layer.refresh()
            except Exception:
                pass
            self._attach_preview_dead_hover()
            if self._is_organoid and self._org_preview_cache is not None:
                self._org_preview_cache["stats_by_frame"] = self._preview_stats_by_frame
                self._org_preview_cache["pct_maps_by_frame"] = self._preview_pct_maps_by_frame
                self._org_preview_cache["computed_frames"] = self._preview_computed_frames
                self._org_preview_cache["current_thr"] = thr
            return

        import datetime
        print(
            f"[{datetime.datetime.now().strftime('%H:%M:%S')}] [Preview] "
            f"Refreshing Dead/Alive overlay — threshold={thr:.2f}% "
            f"(frame={frame_idx}, {self.cell_type})"
        )
        overlay_arr, frame_stats = _update_dead_alive_overlay(
            viewer, self._preview_seg_t, self._preview_dead_t, thr,
            frame_idx=frame_idx if self._preview_is_timelapse else None,
            overlay_arr=self._preview_overlay_arr if self._preview_is_timelapse else None,
            log_fn=self.log,
        )
        if self._preview_is_timelapse:
            self._preview_overlay_arr = overlay_arr
        self._preview_stats_by_frame[frame_idx] = frame_stats
        self._preview_pct_maps_by_frame[frame_idx] = _build_dead_pct_map(frame_stats)
        self._preview_computed_frames.add(frame_idx)
        self._attach_preview_dead_hover()

        if self._is_organoid and self._org_preview_cache is not None:
            self._org_preview_cache["overlay_arr"] = self._preview_overlay_arr
            self._org_preview_cache["stats_by_frame"] = self._preview_stats_by_frame
            self._org_preview_cache["pct_maps_by_frame"] = self._preview_pct_maps_by_frame
            self._org_preview_cache["computed_frames"] = self._preview_computed_frames
            self._org_preview_cache["current_thr"] = thr

    def _on_workers_changed(self, value):
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        if value > max_allowed:
            self.spin_workers.setValue(max_allowed)
            self.log(
                f"⚠️ Workers capped to {max_allowed} (system has {n_cores} cores). "
                "Using all cores can freeze the system."
            )

    def _selected_features(self) -> list:
        result = [f for f, cb in self.feature_checks.items() if cb.isChecked()]
        # Defensive: mandatory features are always included regardless of checkbox state
        for _mf in getattr(self, "_mandatory_features", set()):
            if _mf in self.feature_checks and _mf not in result:
                result.append(_mf)
        return result

    def _collect_params(self) -> dict:
        thr = self._get_threshold()
        return {
            "features_choice": self._selected_features(),
            "contact_threshold": float(self.contact_threshold.value()),
            "dead_mask_percentage_threshold": thr if thr > 0 else None,
            "n_workers": int(self.spin_workers.value()),
        }

    def _apply_to_others(self, category_only=False):
        parent_tab = self.parent()
        while parent_tab and not hasattr(parent_tab, "panels"):
            parent_tab = parent_tab.parent()
        if not parent_tab:
            return

        settings = self._collect_params()
        targets = self.category_types if category_only else self.all_cell_types
        count = 0
        for ct in targets:
            if ct == self.cell_type:
                continue
            if ct in parent_tab.panels:
                p = parent_tab.panels[ct]
                for f, cb in p.feature_checks.items():
                    cb.setChecked(f in settings["features_choice"])
                p.contact_threshold.setValue(settings["contact_threshold"])
                p.spin_workers.setValue(settings["n_workers"])
                # Only sync dead threshold for non-organoid panels that own a spinner
                if not p._is_organoid and p.spin_dead_threshold is not None:
                    thr_val = settings.get("dead_mask_percentage_threshold")
                    p.spin_dead_threshold.setValue(thr_val if thr_val is not None else 0.0)
                count += 1
        scope = "category" if category_only else "all"
        self.log(f"Applied feature settings to {count} other cell type(s) ({scope}).")

    def _persist(self):
        params = self.metadata_loader.behav3d_parameters
        features = params.setdefault("features", {})
        features[self.cell_type] = self._collect_params()

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

        # Remember the threshold value that is now baked into the on-disk
        # config (and, after run, into the combined CSV).
        thr = self._get_threshold()
        self._last_persisted_threshold = thr if thr > 0 else None
        self._refresh_rerun_death_button()

    def _combined_csv_path(self, cell_type: str | None = None) -> Path:
        """Path to the combined features CSV for a given cell type."""
        ct = cell_type or self.cell_type
        return (
            Path(self.metadata_loader.output_dir)
            / "analysis"
            / ct
            / "track_features"
            / f"BEHAV3D_{ct}_combined_track_features.csv"
        )

    def _has_combined_csv(self) -> bool:
        try:
            return self._combined_csv_path().exists()
        except Exception:
            return False

    def _threshold_changed(self) -> bool:
        """Return True when current spinner value differs from the persisted one."""
        if not self._has_combined_csv():
            return False
        try:
            current = round(self._get_threshold(), 2)
        except Exception:
            return False
        persisted = self._last_persisted_threshold
        if persisted is None:
            return current > 0
        return abs(current - round(float(persisted), 2)) > 1e-6

    def _refresh_rerun_death_button(self):
        """Show/hide and enable/disable the Re-run death button per current state."""
        if not hasattr(self, "btn_rerun_death") or self.btn_rerun_death is None:
            return
        if not self._has_dead:
            self.btn_rerun_death.setVisible(False)
            return
        has_csv = self._has_combined_csv()
        self.btn_rerun_death.setVisible(has_csv)
        if not has_csv:
            return
        changed = self._threshold_changed()
        self.btn_rerun_death.setEnabled(changed)
        if changed:
            self.btn_rerun_death.setToolTip(
                "Recompute the 'dead' column using the new threshold "
                "without re-extracting all features."
            )
        else:
            self.btn_rerun_death.setToolTip(
                "Current threshold matches the value used in the existing "
                "features. Nothing to recompute."
            )

    def _on_rerun_death_clicked(self):
        from behav3d.features.timepoint_features import rerun_death_classification

        if not self._has_combined_csv():
            self.log(
                f"\u26a0\ufe0f No combined features CSV for {self.cell_type} \u2014 "
                "run full feature extraction first."
            )
            return
        new_thr = self._get_threshold()
        if new_thr <= 0:
            self.log(
                "\u26a0\ufe0f Dead threshold is 0 \u2014 set a positive value before re-running."
            )
            return
        try:
            rerun_death_classification(
                output_dir=str(Path(self.metadata_loader.output_dir).expanduser()),
                cell_type=self.cell_type,
                new_threshold=float(new_thr),
            )
            self._persist()
            self.log(
                f"\u2705 Re-ran death classification for {self.cell_type} "
                f"with threshold={new_thr}."
            )
            _notify_post_extraction(self)
        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during death re-run: {e}")

    def _run_feature_extraction_for(self, cell_type: str, overwrite: bool = False,
                                    params: dict = None, progress_cb=None):
        """Run feature extraction for a single cell type.

        ``params`` is an optional Qt-thread-safe snapshot from
        :meth:`_collect_params`; when supplied the widget values are not
        re-read (required when called from a background worker).
        ``progress_cb`` is forwarded to ``run_feature_extraction``.
        """
        from behav3d.features.timepoint_features import run_feature_extraction

        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        if params is None:
            params = self._collect_params()

        run_feature_extraction(
            metadata=self.metadata_loader.metadata,
            output_dir=out_dir,
            cell_type=cell_type,
            features_choice=list(params["features_choice"]),
            contact_threshold=float(params["contact_threshold"]),
            dead_mask_percentage_threshold=params["dead_mask_percentage_threshold"],
            n_workers=int(params["n_workers"]),
            overwrite=overwrite,
            progress_cb=progress_cb,
        )

    def _run_death_only_for(self, cell_type: str, params: dict = None):
        """Re-apply only the death classification without recomputing features.

        ``params`` is an optional snapshot from :meth:`_collect_params`;
        when supplied the threshold widget is not re-read.
        """
        from behav3d.features.timepoint_features import rerun_death_classification

        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        if params is not None:
            new_thr = params["dead_mask_percentage_threshold"] or 0.0
        else:
            new_thr = self._get_threshold()
        if new_thr <= 0:
            self.log(
                f"\u26a0\ufe0f Skipping death-only re-run for {cell_type}: "
                "threshold is 0."
            )
            return
        rerun_death_classification(
            output_dir=out_dir,
            cell_type=cell_type,
            new_threshold=float(new_thr),
        )

    def _check_existing_features(self, cell_types: list) -> list:
        warnings = []
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in cell_types:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            combined = feat_dir / f"BEHAV3D_{ct}_combined_track_features.csv"
            if combined.exists():
                warnings.append(f"{ct} feature data ({combined.name})")
        return warnings

    # ── Click handler ────────────────────────────────────────────────────────
    def _on_run_clicked(self, interactive=True):
        """Background-execute feature extraction for this cell type.

        Sample-level progress is forwarded from ``run_feature_extraction``.
        The death-only re-run path uses indeterminate progress because
        ``rerun_death_classification`` has no top-level sample loop.
        """
        from behav3d.napari._overwrite_prompt import prompt_overwrite_single
        from qtpy.QtWidgets import QMessageBox as _QMB

        if self._bg.is_running():
            self.log("⚠️ A feature-extraction run is already in progress for this panel.")
            return

        self.log(f"Running feature extraction for: {self.cell_type}")

        overwrite = False
        death_only = False
        existing = self._check_existing_features([self.cell_type])
        threshold_changed = self._threshold_changed()

        if existing:
            if interactive:
                extra = None
                if threshold_changed:
                    extra = [
                        (
                            "Re-run death only",
                            "death_only",
                            _QMB.ActionRole,
                        ),
                    ]
                choice = prompt_overwrite_single(
                    self,
                    "Overwrite Existing Features?",
                    existing,
                    extra_buttons=extra,
                )
                if choice == "cancel" or choice == "skip":
                    self.log(f"Feature extraction for {self.cell_type} cancelled.")
                    return
                if choice == "death_only":
                    death_only = True
                else:
                    overwrite = True
            else:
                overwrite = True

        # Persist + snapshot widget state on the Qt thread.
        self._persist()
        params_snapshot = self._collect_params()
        cell_type = self.cell_type

        if death_only:
            def _do_death_only(progress_cb=None):
                return self._run_death_only_for(cell_type, params=params_snapshot)

            def _on_done(_r):
                self.log(f"\u2705 {cell_type} death re-run finished.")
                _notify_post_extraction(self)
                notify_results_changed(self)

            def _on_failed(err: str):
                self.log(f"Error during death re-run: {err}")
                notify_results_changed(self)

            self._bg.run(
                fn=_do_death_only,
                desc=f"Death re-run \u2014 {cell_type}\u2026",
                progress_row=self.tab_progress_row,
                buttons=[self.btn_run],
                viewer=self.viewer,
                on_done=_on_done,
                on_failed=_on_failed,
                inject_progress=False,
                indeterminate=True,
            )
            return

        def _do_extraction(progress_cb=None):
            return self._run_feature_extraction_for(
                cell_type, overwrite=overwrite,
                params=params_snapshot, progress_cb=progress_cb,
            )

        def _on_done(_r):
            self.log(f"✅ {cell_type} feature extraction finished.")
            _notify_post_extraction(self)
            notify_results_changed(self)

        def _on_failed(err: str):
            self.log(f"Error during feature extraction: {err}")
            notify_results_changed(self)

        self._bg.run(
            fn=_do_extraction,
            desc=f"Feature extraction \u2014 {cell_type}\u2026",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    # ── Non-organoid Dead Threshold Preview ──────────────────────────────────
    # -- Dead Threshold Preview (all panel types) -------------------------
    def _on_preview_dead_clicked(self):
        """Load raw channels, dead mask, and ALL cell-type segments for the
        selected sample, then show a green/red dead-alive overlay.

        For organoid panels the overlay uses the shared threshold; adjusting
        the spinner in any organoid tab updates all organoid panels live.
        """
        try:
            viewer = self.viewer
            if viewer is None:
                self.log("\u26a0\ufe0f No viewer available for dead threshold preview.")
                return

            md = self.metadata_loader.metadata
            if md is None or md.empty:
                self.log("\u26a0\ufe0f No metadata loaded.")
                return

            # Sample selector
            sample_name = (
                self.preview_sample_combo.currentText()
                if hasattr(self, "preview_sample_combo")
                else ""
            )
            if not sample_name:
                self.log("\u26a0\ufe0f No sample selected in the preview combo.")
                return

            sample_rows = md[md["sample_name"] == sample_name]
            if sample_rows.empty:
                self.log(f"\u26a0\ufe0f Sample '{sample_name}' not in metadata.")
                return
            sample_row = sample_rows.iloc[0]

            output_dir = Path(self.metadata_loader.output_dir)
            self.log(f"Loading preview \u2014 {self.cell_type} / {sample_name}\u2026")

            import datetime
            def _ts():
                return datetime.datetime.now().strftime("%H:%M:%S")

            print(f"\n[{_ts()}] [Preview] {'='*46}")
            print(f"[{_ts()}] [Preview] Dead threshold preview: {self.cell_type} / {sample_name}")
            print(f"[{_ts()}] [Preview] {'='*46}")

            # Resolve dead mask (with diagnostic logging)
            print(f"[{_ts()}] [Preview] Step 1/5 \u2014 Resolving dead mask...")
            dead_arr, dead_method, tried_paths = _resolve_dead_mask(
                sample_row, output_dir, log_fn=self.log
            )
            if dead_arr is None:
                print(f"[{_ts()}] [Preview] \u274c Dead mask not found.")
                self.log(
                    "\u26a0\ufe0f Dead mask not found. Paths tried:\n"
                    + "\n".join(f"    {p}" for p in tried_paths)
                    + "\nRun dead mask segmentation first (Segmentation tab)."
                )
                return
            print(f"[{_ts()}] [Preview]   Dead mask loaded (method={dead_method}).")
            if dead_method == "raw":
                self.log(
                    "\u2139\ufe0f Dead mask estimated from raw dead channel (Otsu). "
                    "Run dedicated dead mask segmentation for best results."
                )

            # Load raw image (all channels, lazy)
            print(f"[{_ts()}] [Preview] Step 2/5 \u2014 Loading raw image (lazy)...")
            raw_dask = _load_raw_dask(sample_row, output_dir)
            if raw_dask is not None:
                print(f"[{_ts()}] [Preview]   Raw image shape: {raw_dask.shape}")
            else:
                print(f"[{_ts()}] [Preview]   Raw image not found \u2014 skipping channel layers.")

            # Collect ALL available segments for this sample
            all_types = getattr(self, "all_cell_types", [self.cell_type])
            print(
                f"[{_ts()}] [Preview] Step 3/5 \u2014 Loading segments for "
                f"{len(all_types)} type(s): {all_types}..."
            )
            segs_dict, seg_sources = _load_all_segments_for_sample(
                sample_row, all_types, output_dir
            )
            for ct, arr in segs_dict.items():
                shape_str = str(np.asarray(arr).shape)
                print(f"[{_ts()}] [Preview]   {ct}: shape={shape_str}, source={seg_sources.get(ct,'?')}")
            if self.cell_type not in segs_dict:
                print(f"[{_ts()}] [Preview] \u274c No segments for '{self.cell_type}' \u2014 aborting.")
                self.log(
                    f"\u26a0\ufe0f No segments found for '{self.cell_type}' "
                    f"in '{sample_name}'. Run segmentation first."
                )
                return

            # Clear all viewer layers entirely
            print(f"[{_ts()}] [Preview] Clearing viewer layers...")
            self._disconnect_preview_dead_hover()
            self._disconnect_preview_dims()
            self.viewer.layers.clear()

            # Raw channels
            if raw_dask is not None:
                n_ch = raw_dask.shape[1] if raw_dask.ndim >= 5 else 1
                _add_channel_layers(viewer, raw_dask, sample_name)
                self.log(f"  Added {n_ch} raw channel(s).")
            else:
                self.log("  \u26a0\ufe0f Could not load raw image.")

            # Dead mask layer (timelapse-aware)
            dead_stack = np.asarray(dead_arr)
            if dead_stack.ndim == 2:
                dead_stack = dead_stack[np.newaxis, np.newaxis, ...]
            elif dead_stack.ndim == 3:
                dead_stack = dead_stack[np.newaxis, ...]
            print(f"[{_ts()}] [Preview]   Dead mask stack shape: {dead_stack.shape}")

            # Mirror the feature-extraction-time cleaning: when the previewed
            # cell type is NOT itself an immune type, zero dead-mask voxels
            # that fall inside any immune segment available for this sample.
            from behav3d.core.metadata import detect_immune_cell_types_from_metadata
            try:
                all_immune_types = detect_immune_cell_types_from_metadata(
                    self.metadata_loader.metadata
                ) or []
            except Exception:
                all_immune_types = []
            immune_types_present = [
                ct for ct in all_immune_types if ct in segs_dict and ct != self.cell_type
            ]
            if self.cell_type not in all_immune_types and immune_types_present:
                any_immune = np.zeros(dead_stack.shape, dtype=bool)
                for im_type in immune_types_present:
                    im_arr = np.asarray(segs_dict[im_type])
                    if im_arr.ndim == 2:
                        im_arr = im_arr[np.newaxis, np.newaxis, ...]
                    elif im_arr.ndim == 3:
                        im_arr = im_arr[np.newaxis, ...]
                    if im_arr.shape != dead_stack.shape:
                        print(
                            f"[{_ts()}] [Preview]   Skipping immune '{im_type}' for "
                            f"dead-mask cleaning (shape {im_arr.shape} != {dead_stack.shape})."
                        )
                        continue
                    any_immune |= (im_arr > 0)
                if any_immune.any():
                    dead_stack = np.where(any_immune, 0, dead_stack)
                    self.log(
                        "  [Preview] Cleaned dead mask: zeroed voxels under immune "
                        f"segments {immune_types_present}"
                    )
                    print(
                        f"[{_ts()}] [Preview]   Cleaned dead mask: zeroed voxels under "
                        f"immune segments {immune_types_present}"
                    )

            viewer.add_image(
                dead_stack.astype(float),
                name=f"{_PREVIEW_PREFIX} Dead Mask",
                colormap="magenta",
                blending="additive",
                opacity=0.4,
            )

            # All segment layers (this type full opacity, others lighter)
            primary_seg = segs_dict[self.cell_type]
            for ct, seg_arr in segs_dict.items():
                is_primary = (ct == self.cell_type)
                viewer.add_labels(
                    seg_arr,
                    name=f"{_PREVIEW_PREFIX} {ct} segments",
                    opacity=0.7 if is_primary else 0.3,
                    visible=False,
                )
            self.log(
                "  Segments: "
                + ", ".join(
                    (
                        f"{ct}[primary,{seg_sources.get(ct, 'untracked')}]"
                        if ct == self.cell_type
                        else f"{ct}[{seg_sources.get(ct, 'untracked')}]"
                    )
                    for ct in segs_dict
                )
            )

            # For organoid panels: merge ALL organoid types into a single
            # segmentation so the Dead/Alive overlay covers every organoid type.
            # The shared threshold applies equally to all.
            if self._is_organoid and self._org_cell_types:
                print(
                    f"[{_ts()}] [Preview] Step 4/5 \u2014 Merging "
                    f"{len(self._org_cell_types)} organoid type(s)..."
                )
                overlay_seg, label_type_map = _merge_org_segments(
                    segs_dict, self._org_cell_types
                )
                if overlay_seg is None:
                    overlay_seg = primary_seg
                    label_type_map = {}
                self._preview_label_type_map = label_type_map
                org_names = [ct for ct in self._org_cell_types if ct in segs_dict]
                n_labels = len(label_type_map)
                print(f"[{_ts()}] [Preview]   Merged: {len(org_names)} type(s), {n_labels} total labels")
                self.log(
                    f"  Dead/Alive overlay merges {len(org_names)} organoid type(s): "
                    + ", ".join(org_names)
                )
            else:
                print(f"[{_ts()}] [Preview] Step 4/5 \u2014 Single type ({self.cell_type}), no merge needed.")
                overlay_seg = primary_seg
                self._preview_label_type_map = {}

            # Compute dead/alive overlay (on-demand: current frame only)
            thr = self._get_threshold()
            print(
                f"[{_ts()}] [Preview] Step 5/5 \u2014 Computing Dead/Alive overlay "
                f"on-demand for current frame (threshold={thr:.2f}%)..."
            )

            # Cache the source arrays + reset per-frame caches
            self._preview_seg_t = overlay_seg
            self._preview_dead_t = dead_stack
            self._preview_overlay_arr = None
            self._preview_pct_maps_by_frame = {}
            self._preview_stats_by_frame = {}
            self._preview_computed_frames = set()
            self._preview_current_thr = None

            overlay_seg_arr = np.asarray(overlay_seg)
            self._preview_is_timelapse = overlay_seg_arr.ndim >= 4

            # Pre-allocate the full overlay so the labels layer aligns with
            # the time slider; on-demand frames are written in place.
            if self._preview_is_timelapse:
                t_common = min(overlay_seg_arr.shape[0], dead_stack.shape[0])
                overlay_shape = tuple(overlay_seg_arr[:t_common].shape)
                self._preview_overlay_arr = np.zeros(overlay_shape, dtype=np.uint8)

            # Tear down any previous time-slider listener and create the
            # layer/compute the current frame.
            self._disconnect_preview_dims()
            current_frame = self._current_viewer_frame() if self._preview_is_timelapse else 0
            self._refresh_preview_dead_layers(thr, frame_idx=current_frame)

            # Organoid panels fill the shared tab cache too
            if self._is_organoid and self._org_preview_cache is not None:
                self._org_preview_cache["seg_t"] = overlay_seg
                self._org_preview_cache["dead_t"] = dead_stack
                self._org_preview_cache["overlay_arr"] = self._preview_overlay_arr
                self._org_preview_cache["pct_maps_by_frame"] = self._preview_pct_maps_by_frame
                self._org_preview_cache["stats_by_frame"] = self._preview_stats_by_frame
                self._org_preview_cache["computed_frames"] = self._preview_computed_frames
                self._org_preview_cache["current_thr"] = self._preview_current_thr
                self._org_preview_cache["label_type_map"] = self._preview_label_type_map
                self._org_preview_cache["is_timelapse"] = self._preview_is_timelapse

            # Connect time-slider listener for on-demand frame recompute.
            # Organoid panels delegate the listener to FeatureExtractionTab
            # (a single shared listener covers all organoid panels).
            if self._preview_is_timelapse and not self._is_organoid:
                self._connect_preview_dims()
            elif self._is_organoid:
                pt = self.parent()
                while pt and not hasattr(pt, "_connect_org_preview_dims"):
                    pt = pt.parent()
                if pt is not None:
                    if self._preview_is_timelapse:
                        pt._connect_org_preview_dims()
                    # Mirror this panel's freshly-loaded preview state onto
                    # every other organoid panel so they share the same
                    # per-frame caches and hover tooltip data.
                    pt._propagate_org_preview_to_panels()

            # Wire non-organoid spinner -> live overlay (once)
            if not self._preview_connected and self.spin_dead_threshold is not None:
                if not self._is_organoid:
                    self.spin_dead_threshold.valueChanged.connect(
                        self._on_threshold_spin_changed
                    )
                self._preview_connected = True

            self.log(
                f"\u2705 Preview loaded (on-demand per-frame). "
                f"Adjust threshold (currently {thr:.2f}%) and use the viewer time slider; "
                f"each frame is computed when you visit it."
            )

        except Exception as exc:
            import traceback as _tb
            _tb.print_exc()
            self.log(f"\u274c Error loading dead threshold preview: {exc}")


    def _on_threshold_spin_changed(self, value):
        """Live-update the overlay (current frame only) for non-organoid panels."""
        viewer = self.viewer
        if (
            viewer is None
            or self._preview_seg_t is None
            or self._preview_dead_t is None
        ):
            return
        # Threshold change invalidates every cached frame; recompute the
        # one the user is currently viewing.
        self._refresh_preview_dead_layers(value)



# ═══════════════════════════════════════════════════════════════════════════
# ActiveKillingPanel — Extended analysis for immune cells
# ═══════════════════════════════════════════════════════════════════════════
class ActiveKillingPanel(QWidget):
    """
    Extended analysis: Active Killing Analysis for immune cell types only.

    Runs AFTER baseline feature extraction (requires the combined_track_features
    CSV to be present for the chosen immune type).  Detects functional killing
    events by comparing the death-signal change in target (organoid) cells after
    immune-cell contact to a per-sample background death rate.

    Output
    ------
    CSVs  -> analysis/<immune>/active_killing/
    Plots -> analysis/<immune>/active_killing/plots/
    Viewer-> Top-N active killers shown as coloured Points layers (ALL killing
              timepoints per killer, one layer per cell).
    """

    _LAYER_PREFIX = "[Active Killing]"

    def __init__(
        self,
        immune_types: list,
        metadata_loader,
        viewer=None,
        log_callback=None,
        queue_callback=None,
        parent=None,
        tab_progress_row=None,
    ):
        super().__init__(parent)
        self.immune_types = list(immune_types)
        self.metadata_loader = metadata_loader
        self.viewer = viewer
        self.log = log_callback or (lambda m: None)
        self._queue_callback = queue_callback
        # Background-execution infrastructure.
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)
        self._init_ui()

    # ── UI ──────────────────────────────────────────────────────────────────
    def _init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        desc = QLabel(
            "Detects functional killing events by analysing how much the death "
            "signal in target (organoid) cells increases after immune cell contact, "
            "relative to the per-sample background death rate.\n"
            "\u26a0\ufe0f  Run baseline feature extraction for the immune type FIRST."
        )
        desc.setWordWrap(True)
        desc.setStyleSheet("color: #90A4AE; font-size: 10px; padding: 2px 0;")
        layout.addWidget(desc)

        # ── Immune cell type selector ──────────────────────────────────────
        imm_row = QHBoxLayout()
        imm_row.addWidget(QLabel("Immune cell type:"))
        self.immune_combo = QComboBox()
        if self.immune_types:
            self.immune_combo.addItems(self.immune_types)
        else:
            self.immune_combo.addItem("(no immune types detected)")
            self.immune_combo.setEnabled(False)
        self.immune_combo.currentTextChanged.connect(self._validate)
        imm_row.addWidget(self.immune_combo, stretch=1)
        layout.addLayout(imm_row)

        self.validation_label = QLabel("")
        self.validation_label.setWordWrap(True)
        self.validation_label.setStyleSheet("font-size: 10px;")
        layout.addWidget(self.validation_label)

        # ── Parameters form ────────────────────────────────────────────────
        params_group = QGroupBox("Parameters")
        params_form = QFormLayout(params_group)
        params_form.setContentsMargins(6, 6, 6, 6)
        params_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.spin_obs_window = QSpinBox()
        self.spin_obs_window.setRange(1, 100)
        self.spin_obs_window.setValue(5)
        self.spin_obs_window.setMaximumWidth(70)
        params_form.addRow(
            "Observation window (tp):",
            make_help_row(
                self.spin_obs_window,
                "Observation Window",
                "Timepoints after each contact timepoint to measure death-signal change.\n\n"
                "A contact is 'active killing' if the death-signal increase exceeds:\n"
                "  background_rate x observation_window x threshold_multiplier",
            ),
        )

        self.death_signal_combo = QComboBox()
        self.death_signal_combo.addItems(
            ["percentage_dead_mask", "mean_dead_dye", "nr_dead_mask_pixels"]
        )
        self.death_signal_combo.setCurrentText("percentage_dead_mask")
        self.death_signal_combo.setMaximumWidth(220)
        params_form.addRow("Death signal column:", self.death_signal_combo)

        self.spin_threshold_mult = QDoubleSpinBox()
        self.spin_threshold_mult.setRange(0.1, 20.0)
        self.spin_threshold_mult.setSingleStep(0.1)
        self.spin_threshold_mult.setDecimals(2)
        self.spin_threshold_mult.setValue(1.5)
        self.spin_threshold_mult.setMaximumWidth(90)
        params_form.addRow(
            "Killing threshold multiplier:",
            make_help_row(
                self.spin_threshold_mult,
                "Killing Threshold Multiplier",
                "Multiplier applied to the per-sample background death rate.\n\n"
                "A contact is classified as 'active killing' only when:\n"
                "  death_increase > background_rate x window x multiplier\n\n"
                "Default: 1.5  (signal must exceed 150% of expected background).",
            ),
        )

        self.check_abs_threshold = QCheckBox("Use absolute threshold instead of multiplier")
        self.check_abs_threshold.stateChanged.connect(self._on_abs_toggle)
        params_form.addRow("", self.check_abs_threshold)

        self.spin_abs_threshold = QDoubleSpinBox()
        self.spin_abs_threshold.setRange(0.0, 100.0)
        self.spin_abs_threshold.setSingleStep(0.01)
        self.spin_abs_threshold.setDecimals(4)
        self.spin_abs_threshold.setValue(0.0)
        self.spin_abs_threshold.setMaximumWidth(100)
        self.spin_abs_threshold.setEnabled(False)
        params_form.addRow(
            "Absolute threshold:",
            make_help_row(
                self.spin_abs_threshold,
                "Absolute Killing Threshold",
                "Fixed minimum death-signal increase required to classify a contact\n"
                "as active killing.  Bypasses the multiplier-based threshold.\n\n"
                "Useful when the background death rate is near zero.",
            ),
        )

        self.spin_min_contact = QSpinBox()
        self.spin_min_contact.setRange(1, 50)
        self.spin_min_contact.setValue(1)
        self.spin_min_contact.setMaximumWidth(70)
        params_form.addRow(
            "Min contact duration (tp):",
            make_help_row(
                self.spin_min_contact,
                "Minimum Contact Duration",
                "Minimum number of consecutive timepoints an immune cell must remain\n"
                "in contact with a target cell for the event to qualify.\n\n"
                "Default: 1 (all contacts are included).",
            ),
        )
        layout.addWidget(params_group)

        # ── Viewer preview ─────────────────────────────────────────────────
        viewer_group = QGroupBox("Viewer Preview — Top Active Killers")
        viewer_form = QFormLayout(viewer_group)
        viewer_form.setContentsMargins(6, 6, 6, 6)

        self.spin_top_n = QSpinBox()
        self.spin_top_n.setRange(1, 50)
        self.spin_top_n.setValue(5)
        self.spin_top_n.setMaximumWidth(70)
        viewer_form.addRow(
            "Top-N killers to display:",
            make_help_row(
                self.spin_top_n,
                "Top-N Active Killers",
                "Number of top-ranking immune cells (by total active-killing timepoints)\n"
                "to visualise in the napari viewer.  Each cell gets a colour-coded Points\n"
                "layer containing ALL its active-killing timepoints.",
            ),
        )

        self.btn_load_viewer = QPushButton("\U0001f441  Load Top Killers in Viewer")
        self.btn_load_viewer.setStyleSheet(
            "QPushButton { background: #37474F; color: white; padding: 5px 10px; "
            "border-radius: 3px; font-size: 11px; } "
            "QPushButton:hover { background: #546E7A; }"
        )
        self.btn_load_viewer.setToolTip(
            "Loads the top-N most active immune cells as separate colour-coded\n"
            "Points layers in the napari viewer.\n\n"
            "Each layer shows ALL timepoints where that cell was actively killing.\n"
            "Requires Active Killing Analysis to have been run at least once."
        )
        self.btn_load_viewer.clicked.connect(self._on_load_viewer_clicked)
        viewer_form.addRow("", self.btn_load_viewer)
        layout.addWidget(viewer_group)

        # ── Run button ─────────────────────────────────────────────────────
        action_row = QHBoxLayout()
        action_row.setSpacing(6)

        self.btn_run = QPushButton("\u25b6  Run Active Killing Analysis")
        self.btn_run.setStyleSheet(
            "background-color: #c0392b; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 8px; font-size: 13px;"
        )
        self.btn_run.clicked.connect(self._on_run_clicked)
        action_row.addWidget(self.btn_run, stretch=1)

        self.btn_queue = QPushButton("+\U0001f6d2")
        self.btn_queue.setFixedSize(40, 34)
        self.btn_queue.setToolTip("Add Active Killing Analysis to Processing Queue")
        self.btn_queue.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 12px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
            "QPushButton:disabled { background: #2a2a2a; color: #666; border: 1px solid #555; }"
        )
        self.btn_queue.clicked.connect(self._on_queue_clicked)
        action_row.addWidget(self.btn_queue)

        layout.addLayout(action_row)
        layout.addStretch()

        self._validate()

    # ── Helpers ──────────────────────────────────────────────────────────────
    def _on_abs_toggle(self, state):
        enabled = (state == Qt.Checked)
        self.spin_abs_threshold.setEnabled(enabled)
        self.spin_threshold_mult.setEnabled(not enabled)

    def _get_immune_type(self) -> str:
        return self.immune_combo.currentText()

    def _feature_csv_path(self, immune_type: str) -> Path:
        out = Path(self.metadata_loader.output_dir)
        feat_dir = out / "analysis" / immune_type / "track_features"
        filtered = feat_dir / f"BEHAV3D_{immune_type}_combined_track_features_filtered.csv"
        return filtered if filtered.exists() else (
            feat_dir / f"BEHAV3D_{immune_type}_combined_track_features.csv"
        )

    def _active_killing_dir(self, immune_type: str) -> Path:
        return Path(self.metadata_loader.output_dir) / "analysis" / immune_type / "active_killing"

    def set_queue_callback(self, callback):
        """Set the callback used by the cart button to enqueue this analysis."""
        self._queue_callback = callback
        if hasattr(self, "btn_queue"):
            self.btn_queue.setEnabled(callback is not None and self.btn_run.isEnabled())

    def get_queue_params(self) -> dict:
        """Snapshot the current Active Killing analysis settings for the queue."""
        return {
            "immune_type": self._get_immune_type(),
            **self._collect_params(),
        }

    def _on_queue_clicked(self):
        if self._queue_callback is None:
            self.log("⚠️ Active Killing queue is not connected.")
            return
        self._validate()
        if not self.btn_run.isEnabled():
            return
        self._queue_callback()

    def _validate(self):
        immune = self._get_immune_type()
        if not immune or "(no immune" in immune:
            self.validation_label.setText("\u26a0\ufe0f No immune cell types detected in metadata.")
            self.validation_label.setStyleSheet("color: #E57373; font-size: 10px;")
            self.btn_run.setEnabled(False)
            if hasattr(self, "btn_queue"):
                self.btn_queue.setEnabled(False)
            return
        csv = self._feature_csv_path(immune)
        if not csv.exists():
            self.validation_label.setText(
                f"\u26a0\ufe0f {immune} feature CSV not found \u2014 run feature extraction first.\n"
                f"Expected: .../analysis/{immune}/track_features/"
                f"BEHAV3D_{immune}_combined_track_features.csv"
            )
            self.validation_label.setStyleSheet("color: #E57373; font-size: 10px;")
            self.btn_run.setEnabled(False)
            if hasattr(self, "btn_queue"):
                self.btn_queue.setEnabled(self._queue_callback is not None and False)
        else:
            self.validation_label.setText(f"\u2713 Ready  \u2014  using: {csv.name}")
            self.validation_label.setStyleSheet("color: #66BB6A; font-size: 10px;")
            self.btn_run.setEnabled(True)
            if hasattr(self, "btn_queue"):
                self.btn_queue.setEnabled(self._queue_callback is not None)

    def _collect_params(self) -> dict:
        return {
            "observation_window": int(self.spin_obs_window.value()),
            "death_signal_column": self.death_signal_combo.currentText(),
            "killing_threshold_multiplier": float(self.spin_threshold_mult.value()),
            "absolute_killing_threshold": (
                float(self.spin_abs_threshold.value())
                if self.check_abs_threshold.isChecked() else None
            ),
            "min_contact_duration": int(self.spin_min_contact.value()),
        }

    def refresh_immune_types(self, immune_types: list):
        """Update the dropdown when metadata is reloaded."""
        self.immune_types = list(immune_types)
        current = self.immune_combo.currentText()
        self.immune_combo.blockSignals(True)
        self.immune_combo.clear()
        if immune_types:
            self.immune_combo.addItems(immune_types)
            if current in immune_types:
                self.immune_combo.setCurrentText(current)
            self.immune_combo.setEnabled(True)
        else:
            self.immune_combo.addItem("(no immune types detected)")
            self.immune_combo.setEnabled(False)
        self.immune_combo.blockSignals(False)
        self._validate()

    # ── Run ──────────────────────────────────────────────────────────────────
    def _on_run_clicked(self):
        """Run Active Killing Analysis in the background (indeterminate)."""
        self._validate()
        if not self.btn_run.isEnabled():
            return
        if self._bg.is_running():
            self.log("\u26a0\ufe0f An active killing run is already in progress.")
            return

        from behav3d.features.advanced_timepoint_features import run_active_killing_analysis
        from behav3d.core.metadata import detect_organoid_types_from_metadata

        immune = self._get_immune_type()
        params = self._collect_params()
        md = self.metadata_loader.metadata
        target_types = detect_organoid_types_from_metadata(md)
        if not target_types:
            self.log("\u26a0\ufe0f No organoid target types detected \u2014 cannot run active killing analysis.")
            return

        self.btn_run.setText("\u23f3 Running\u2026")
        self.log(
            f"\u25b6 Active Killing Analysis: {immune} vs {target_types}  "
            f"(window={params['observation_window']}, "
            f"signal={params['death_signal_column']}, "
            f"multiplier={params['killing_threshold_multiplier']})\u2026"
        )

        output_dir_str = str(self.metadata_loader.output_dir)

        def _do_active_killing(progress_cb=None):
            df_killing, df_summary, stats = run_active_killing_analysis(
                metadata=md,
                output_dir=output_dir_str,
                immune_cell_type=immune,
                target_cell_types=target_types,
                observation_window=params["observation_window"],
                death_signal_column=params["death_signal_column"],
                killing_threshold_multiplier=params["killing_threshold_multiplier"],
                absolute_killing_threshold=params["absolute_killing_threshold"],
                min_contact_duration=params["min_contact_duration"],
                save_results=True,
            )
            return (df_killing, df_summary, stats)

        def _on_done(result):
            df_killing, _df_summary, stats = result
            results_dir = self._active_killing_dir(immune)
            self._save_plots(df_killing, immune, results_dir)
            n_active = int(stats.get("total_active_killing_timepoints", 0))
            rate = stats.get("overall_killing_rate", 0.0)
            self.log(
                f"\u2705 Active Killing Analysis complete \u2014 "
                f"{n_active} active killing timepoints ({rate:.1%} of contact timepoints)."
            )
            self.btn_run.setText("\u25b6  Run Active Killing Analysis")
            self._offer_open_folder(results_dir)

        def _on_failed(err: str):
            self.log(f"\u274c Active Killing Analysis error: {err}")
            self.btn_run.setText("\u25b6  Run Active Killing Analysis")

        self._bg.run(
            fn=_do_active_killing,
            desc=f"Active Killing \u2014 {immune}\u2026",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
            inject_progress=False,
            indeterminate=True,
        )

    # ── Plot saving ───────────────────────────────────────────────────────────
    def _save_plots(self, df_killing, immune: str, results_dir: Path):
        """Save per-sample 4-panel kinetics plots + combined efficiency plot.
        Output layout mirrors the notebook ActiveKillingPanel exactly."""
        import matplotlib.pyplot as plt

        if df_killing is None or df_killing.empty:
            self.log("  \u2139\ufe0f No killing events \u2014 plots skipped.")
            return

        # Normalise column names
        if "TrackID" in df_killing.columns and "immune_track_id" not in df_killing.columns:
            df_killing = df_killing.copy()
            df_killing["immune_track_id"] = df_killing["TrackID"]

        try:
            import seaborn as sns
            _has_sns = True
        except ImportError:
            _has_sns = False

        df_active_all = df_killing[df_killing["is_active_killing"]].copy()
        if not df_active_all.empty and "killing_efficiency" in df_active_all.columns:
            fig, ax = plt.subplots(figsize=(10, 6))
            if _has_sns:
                sns.histplot(df_active_all["killing_efficiency"], kde=True, ax=ax, color="purple")
            else:
                ax.hist(df_active_all["killing_efficiency"].dropna(), bins=20, color="purple", alpha=0.8)
            ax.set_title("Combined Killing Efficiency Distribution", fontsize=16, fontweight="bold")
            ax.set_xlabel("Efficiency Score (signal increase / expected background)")
            ax.set_ylabel("Active Killing Events")
            combined_dir = results_dir / "plots"
            combined_dir.mkdir(parents=True, exist_ok=True)
            combined_path = combined_dir / "combined_killing_efficiency_distribution.png"
            plt.savefig(combined_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            try:
                self.log(f"  📊 Saved: {combined_path.relative_to(Path(self.metadata_loader.output_dir))}")
            except ValueError:
                self.log(f"  📊 Saved: {combined_path}")

    def run_analysis(self, interactive: bool = True, extra_callbacks=None):
        """Synchronous Active Killing analysis used by the processing queue.

        Kept blocking on purpose: the queue iterates steps serially and
        relies on each step completing before moving on.  GUI button
        clicks go through :meth:`_on_run_clicked` instead, which routes
        through ``BackgroundOperation`` for napari responsiveness.

        ``extra_callbacks`` is the queue's chaining hook fired once the
        synchronous analysis finishes.
        """
        self._validate()
        if not self.btn_run.isEnabled():
            fire_extra_callback(extra_callbacks, "on_failed", "run button disabled")
            return

        from behav3d.features.advanced_timepoint_features import run_active_killing_analysis
        from behav3d.core.metadata import detect_organoid_types_from_metadata

        immune = self._get_immune_type()
        params = self._collect_params()
        md = self.metadata_loader.metadata
        target_types = detect_organoid_types_from_metadata(md)
        if not target_types:
            self.log("⚠️ No organoid target types detected — cannot run active killing analysis.")
            fire_extra_callback(extra_callbacks, "on_failed", "no organoid target types")
            return

        self.btn_run.setEnabled(False)
        self.btn_run.setText("⏳ Running…")
        if hasattr(self, "btn_queue"):
            self.btn_queue.setEnabled(False)
        try:
            self.log(
                f"▶ Active Killing Analysis: {immune} vs {target_types}  "
                f"(window={params['observation_window']}, "
                f"signal={params['death_signal_column']}, "
                f"multiplier={params['killing_threshold_multiplier']})…"
            )
            df_killing, df_summary, stats = run_active_killing_analysis(
                metadata=md,
                output_dir=str(self.metadata_loader.output_dir),
                immune_cell_type=immune,
                target_cell_types=target_types,
                observation_window=params["observation_window"],
                death_signal_column=params["death_signal_column"],
                killing_threshold_multiplier=params["killing_threshold_multiplier"],
                absolute_killing_threshold=params["absolute_killing_threshold"],
                min_contact_duration=params["min_contact_duration"],
                save_results=True,
            )
            results_dir = self._active_killing_dir(immune)
            self._save_plots(df_killing, immune, results_dir)
            n_active = int(stats.get("total_active_killing_timepoints", 0))
            rate = stats.get("overall_killing_rate", 0.0)
            self.log(
                f"✅ Active Killing Analysis complete — "
                f"{n_active} active killing timepoints ({rate:.1%} of contact timepoints)."
            )
            if interactive:
                self._offer_open_folder(results_dir)
            fire_extra_callback(extra_callbacks, "on_done", (df_killing, df_summary, stats))
        except Exception as e:
            import traceback as _tb
            _tb.print_exc()
            self.log(f"❌ Active Killing Analysis error: {e}")
            fire_extra_callback(extra_callbacks, "on_failed", str(e))
        finally:
            self.btn_run.setEnabled(True)
            if hasattr(self, "btn_queue"):
                self.btn_queue.setEnabled(self._queue_callback is not None and self.btn_run.isEnabled())
            self.btn_run.setText("▶  Run Active Killing Analysis")
            notify_results_changed(self)

    # ── Folder open popup ──────────────────────────────────────────────────────
    def _offer_open_folder(self, results_dir: Path):
        """Show a modal dialog offering to open the output folder in the OS file manager."""
        box = QMessageBox(self)
        box.setWindowTitle("Active Killing Analysis Complete")
        box.setText(
            "\u2705  Active Killing Analysis finished!\n\n"
            f"Outputs saved to:\n{results_dir}\n\n"
            "Open output folder in file manager?"
        )
        btn_open = box.addButton("Open Folder", QMessageBox.AcceptRole)
        box.addButton("Close", QMessageBox.RejectRole)
        box.exec_()
        if box.clickedButton() == btn_open:
            import subprocess
            import sys as _sys
            try:
                if _sys.platform == "win32":
                    subprocess.Popen(["explorer", str(results_dir)])
                elif _sys.platform == "darwin":
                    subprocess.Popen(["open", str(results_dir)])
                else:
                    subprocess.Popen(["xdg-open", str(results_dir)])
            except Exception as e:
                self.log(f"Could not open folder: {e}")

    # ── Viewer loading ────────────────────────────────────────────────────────
    def _on_load_viewer_clicked(self):
        """Load top-N active killers as colour-coded Points layers in the napari viewer.

        Each killer gets its own layer containing ALL timepoints where it was
        actively killing (is_active_killing == True), coloured by rank order.
        """
        if self.viewer is None:
            self.log("\u26a0\ufe0f No viewer available.")
            return

        immune = self._get_immune_type()
        if not immune or "(no immune" in immune:
            self.log("\u26a0\ufe0f No immune type selected.")
            return

        advanced_path = (
            self._active_killing_dir(immune)
            / f"BEHAV3D_{immune}_advanced_track_features.csv"
        )
        if not advanced_path.exists():
            self.log(
                f"\u26a0\ufe0f Advanced features CSV not found at:\n  {advanced_path}\n"
                "Run Active Killing Analysis first."
            )
            return

        try:
            df = pd.read_csv(advanced_path)

            # Normalise TrackID -> immune_track_id
            if "TrackID" in df.columns and "immune_track_id" not in df.columns:
                df["immune_track_id"] = df["TrackID"]

            # Normalise centroid-* -> position_*
            for old, new in [
                ("centroid-0", "position_z"),
                ("centroid-1", "position_y"),
                ("centroid-2", "position_x"),
            ]:
                if old in df.columns and new not in df.columns:
                    df[new] = df[old]

            df_active = df[df["is_active_killing"] == True].copy()
            if df_active.empty:
                self.log("\u2139\ufe0f No active killing events found \u2014 nothing to display.")
                return

            # Remove previous active killing layers
            to_remove = [
                l for l in list(self.viewer.layers)
                if l.name.startswith(self._LAYER_PREFIX)
            ]
            for l in to_remove:
                try:
                    self.viewer.layers.remove(l)
                except Exception:
                    pass

            # Rank by total active-killing timepoints across all samples
            n_top = int(self.spin_top_n.value())
            top_killers = (
                df_active.groupby(["sample_name", "immune_track_id"])
                .size()
                .sort_values(ascending=False)
                .head(n_top)
            )

            # Distinct palette for up to 8 top killers (cycles after that)
            _palette = [
                [1.00, 0.18, 0.18, 0.90],  # red
                [1.00, 0.55, 0.00, 0.90],  # orange
                [0.93, 0.83, 0.00, 0.90],  # yellow
                [0.13, 0.70, 0.27, 0.90],  # green
                [0.13, 0.47, 0.90, 0.90],  # blue
                [0.60, 0.13, 0.90, 0.90],  # purple
                [0.90, 0.13, 0.54, 0.90],  # pink
                [0.00, 0.84, 0.84, 0.90],  # cyan
            ]

            n_loaded = 0
            for i, ((sample_name, track_id), count) in enumerate(top_killers.items()):
                mask = (
                    (df_active["sample_name"] == sample_name)
                    & (df_active["immune_track_id"] == track_id)
                )
                rows = df_active[mask].sort_values("position_t")
                coords = [
                    [
                        float(r.get("position_t", 0)),
                        float(r.get("position_z", 0)),
                        float(r.get("position_y", 0)),
                        float(r.get("position_x", 0)),
                    ]
                    for _, r in rows.iterrows()
                ]
                if not coords:
                    continue
                color = _palette[i % len(_palette)]
                layer_name = (
                    f"{self._LAYER_PREFIX} {immune} #{track_id} "
                    f"({sample_name}, {count} events)"
                )
                self.viewer.add_points(
                    np.array(coords),
                    name=layer_name,
                    face_color=[color],
                    size=8,
                    symbol="disc",
                    opacity=0.85,
                    out_of_slice_display=True,
                )
                n_loaded += 1

            self.log(
                f"\u2705 Loaded {n_loaded} top active killer(s) as Points layers. "
                "Each layer = ALL killing timepoints for that immune cell."
            )

        except Exception as e:
            import traceback as _tb
            _tb.print_exc()
            self.log(f"\u274c Error loading top killers in viewer: {e}")


# ═══════════════════════════════════════════════════════════════════════════
# FeatureExtractionTab - main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class FeatureExtractionTab(QWidget):
    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeFeaturePanel] = {}
        self._queue_panel = None

        # Background-execution infrastructure for batch + active-killing.
        self._bg = BackgroundOperation(self)

        # Track which cell types are organoids (needed for sync)
        self._org_types: list = []

        # Shared preview cache for organoid panels (all org panels share same ref).
        # Per-frame caches are populated lazily as the user navigates timepoints.
        self._org_preview_cache: dict = {
            "seg_t": None,
            "dead_t": None,
            "overlay_arr": None,
            "stats_by_frame": {},
            "pct_maps_by_frame": {},
            "computed_frames": set(),
            "current_thr": None,
            "label_type_map": None,
            "is_timelapse": False,
        }

        # Single shared napari time-slider listener for ALL organoid panels.
        # Attached the first time any organoid panel loads a preview.
        self._org_preview_dims_callback = None

        self._init_ui()

        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    # ── UI ──────────────────────────────────────────────────────────────────
    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Vertical splitter: tab content on top, shared Results panel
        # underneath (re-scans the same output directory used by the
        # other BEHAV3D tabs).
        self.splitter = QSplitter(Qt.Vertical, self)
        self.splitter.setChildrenCollapsible(True)
        outer.addWidget(self.splitter)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.splitter.addWidget(scroll)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)
        scroll.setWidget(content)

        self.results_panel = ResultsPanel(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            parent=self,
        )
        self.splitter.addWidget(self.results_panel)
        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, True)
        self.splitter.setSizes([600, 400])

        # ── Global Death Classification group (Organoids only) ─────────────

        # ── Cell type grouping ─────────────────────────────────────────────
        grouping_row = QHBoxLayout()
        grouping_row.setContentsMargins(0, 0, 0, 4)
        self.btn_open_grouping = QPushButton("🧬  Cell type grouping…")
        self.btn_open_grouping.setToolTip(
            "Group several detected cell types under a single '{name}_merged' "
            "label.\nThe new group is added as metadata columns and processed "
            "downstream alongside the original cell types."
        )
        self.btn_open_grouping.setStyleSheet(
            "QPushButton { background:#1f3a5f; color:#bbdefb; font-weight:bold; "
            "border:1px solid #64b5f6; border-radius:4px; padding:6px 12px; } "
            "QPushButton:hover { background:#2c4f7f; }"
            "QPushButton:disabled { color:#666; border-color:#444; }"
        )
        self.btn_open_grouping.clicked.connect(self._on_open_grouping_dialog)
        self.btn_open_grouping.setEnabled(False)
        grouping_row.addWidget(self.btn_open_grouping)
        grouping_row.addStretch()
        layout.addLayout(grouping_row)

        # ── Per-cell-type sub-tabs ─────────────────────────────────────────
        self.cell_tabs = QTabWidget()
        self.cell_tabs.setTabPosition(QTabWidget.West)
        layout.addWidget(self.cell_tabs)

        # ── Global Run + Queue ─────────────────────────────────────────────
        self.btn_run_batch = QPushButton(
            "Run Batch Feature Extraction (All Cell Types)"
        )
        self.btn_run_batch.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_batch.clicked.connect(self._on_run_batch_clicked)

        self.btn_queue_feature = QPushButton("+🛒")
        self.btn_queue_feature.setFixedSize(36, 32)
        self.btn_queue_feature.setToolTip("Add Feature Extraction to Processing Queue")
        self.btn_queue_feature.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_row = QHBoxLayout()
        batch_row.setSpacing(4)
        batch_row.addWidget(self.btn_run_batch, stretch=1)
        batch_row.addWidget(self.btn_queue_feature)
        layout.addLayout(batch_row)

        self.btn_run_batch.setVisible(False)
        self.btn_queue_feature.setVisible(False)

        # Shared progress row fed by every Run button on this tab.
        self.progress_row = ProgressBarRow()
        layout.addWidget(self.progress_row)

        # -- Active Killing Extended Analysis (collapsible) --------------------------
        # Panel body is hidden by default; user clicks the toggle to expand.
        self._ak_toggle_label_base = (
            "\U0001f9ec  Extended Analysis \u2014 Active Killing (Immune Cells)"
        )
        self._ak_toggle_btn = QPushButton(
            f"\u25b6  {self._ak_toggle_label_base}"
        )
        self._ak_toggle_btn.setCheckable(True)
        self._ak_toggle_btn.setChecked(False)
        self._ak_toggle_btn.setStyleSheet(
            "QPushButton { background: #2c1810; color: #e74c3c; font-weight: bold; "
            "border: 1px solid #c0392b; border-radius: 4px; padding: 6px 10px; "
            "text-align: left; } "
            "QPushButton:checked { background: #3d1f15; } "
            "QPushButton:hover { background: #3d1f15; }"
        )
        self._ak_toggle_btn.clicked.connect(self._on_ak_toggle)
        self._ak_toggle_btn.setVisible(False)
        layout.addWidget(self._ak_toggle_btn)

        self._ak_body = QWidget()
        self._ak_body.setVisible(False)
        self._ak_inner_layout = QVBoxLayout(self._ak_body)
        self._ak_inner_layout.setContentsMargins(4, 4, 4, 4)
        self.active_killing_panel = None  # created in _rebuild_tabs when immune types are known
        layout.addWidget(self._ak_body)

        # Alias kept for _rebuild_tabs compat
        self._active_killing_group = self._ak_body


        # ── Log ────────────────────────────────────────────────────────────
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(120)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        # Placeholder shown before metadata is loaded
        self._placeholder = QLabel(
            "Load metadata in the Data Preparation tab to see feature extraction options."
        )
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.cell_tabs.addTab(self._placeholder, "—")

    # ── Helpers ──────────────────────────────────────────────────────────────
    def _log(self, msg):
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing feature extraction tabs…")
        self._rebuild_tabs()

    def request_tab_exit(self) -> bool:
        """Block tab switching while a background feature run is in flight."""
        from qtpy.QtWidgets import QMessageBox

        if self._bg.is_running():
            QMessageBox.information(
                self,
                "Operation in progress",
                "A feature extraction run is still in progress. Please "
                "wait for it to finish before switching tabs.",
            )
            return False
        for panel in self.panels.values():
            bg = getattr(panel, "_bg", None)
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "A feature extraction run is still in progress. "
                    "Please wait for it to finish before switching tabs.",
                )
                return False
        ak = getattr(self, "active_killing_panel", None)
        if ak is not None:
            bg = getattr(ak, "_bg", None)
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "An active-killing analysis is still in progress. "
                    "Please wait for it to finish before switching tabs.",
                )
                return False
        return True

    def _detect_cell_types(self):
        """Return (organoid, immune, other) suitable for feature extraction.

        Mirrors :func:`behav3d.widgets.analysis._detect_downstream_cell_types`:
        per-channel multicolor inputs (``*_N_multicolor``) are stripped because
        feature extraction operates on the aggregated outputs (``*_merged`` or
        ``*_grouped``) emitted by the multicolor merging step. Merged /
        grouped names returned by the metadata helpers are preserved as-is.
        """
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            filter_multicolor_inputs,
        )
        md = self.metadata_loader.metadata
        if md is None:
            return [], [], []
        return (
            filter_multicolor_inputs(detect_organoid_types_from_metadata(md)),
            filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)),
            filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)),
        )

    def _rebuild_tabs(self):
        self.cell_tabs.clear()
        self.panels.clear()

        org, imm, oth = self._detect_cell_types()
        all_types = org + imm + oth
        self._org_types = list(org)

        # Dead channel detection is handled per-panel (CellTypeFeaturePanel).
        # Organoid threshold sync handled via _notify_organoid_threshold_changed.

        if not all_types:
            self.cell_tabs.addTab(self._placeholder, "—")
            self.btn_run_batch.setVisible(False)
            self.btn_queue_feature.setVisible(False)
            if hasattr(self, "btn_open_grouping"):
                self.btn_open_grouping.setEnabled(False)
            self._ak_toggle_btn.setVisible(False)
            self._ak_toggle_btn.blockSignals(True)
            self._ak_toggle_btn.setChecked(False)
            self._ak_toggle_btn.blockSignals(False)
            self._set_ak_collapsible_state(False)
            return

        self.btn_run_batch.setVisible(True)
        self.btn_queue_feature.setVisible(True)
        if hasattr(self, "btn_open_grouping"):
            self.btn_open_grouping.setEnabled(True)

        color_map = {"organoid": "🟣", "immune": "🔵", "other": "🟡"}
        for ct in org:
            self._add_panel(ct, "organoid", all_types, org, color_map, is_organoid=True)
        for ct in imm:
            self._add_panel(ct, "immune", all_types, imm, color_map, is_organoid=False)
        for ct in oth:
            self._add_panel(ct, "other", all_types, oth, color_map, is_organoid=False)

        # -- Refresh / create Active Killing panel ----------------------------
        if imm:
            if self.active_killing_panel is None:
                # First load: create and embed the widget
                self.active_killing_panel = ActiveKillingPanel(
                    immune_types=list(imm),
                    metadata_loader=self.metadata_loader,
                    viewer=self.viewer,
                    log_callback=self._log,
                    tab_progress_row=self.progress_row,
                )
                self._ak_inner_layout.addWidget(self.active_killing_panel)
            else:
                # Subsequent metadata reloads: refresh state in-place
                self.active_killing_panel.metadata_loader = self.metadata_loader
                self.active_killing_panel.viewer = self.viewer
                self.active_killing_panel.refresh_immune_types(list(imm))
            self.active_killing_panel.set_queue_callback(
                self._queue_active_killing if self._queue_panel is not None else None
            )
            self._ak_toggle_btn.setVisible(True)
            self._set_ak_collapsible_state(self._ak_toggle_btn.isChecked())
        else:
            self._ak_toggle_btn.setVisible(False)
            self._ak_toggle_btn.blockSignals(True)
            self._ak_toggle_btn.setChecked(False)
            self._ak_toggle_btn.blockSignals(False)
            self._set_ak_collapsible_state(False)

    def _add_panel(self, ct, category, all_types, cat_types, color_map, is_organoid: bool):
        panel = CellTypeFeaturePanel(
            cell_type=ct,
            category=category,
            metadata_loader=self.metadata_loader,
            all_cell_types=all_types,
            category_types=cat_types,
            log_callback=self._log,
            is_organoid=is_organoid,
            threshold_getter=None,
            tab_progress_row=self.progress_row,
        )
        panel.viewer = self.viewer
        if is_organoid:
            panel._org_preview_cache = self._org_preview_cache
            panel._org_cell_types = list(self._org_types)
        self.panels[ct] = panel
        icon = color_map.get(category, "")
        self.cell_tabs.addTab(panel, f"{icon} {ct}")

    def _on_open_grouping_dialog(self):
        """Open the cell-type grouping dialog and rebuild tabs on success."""
        from behav3d.napari._grouping_dialog import GroupBuilderDialog

        if self.metadata_loader is None or self.metadata_loader.metadata is None:
            self._log("\u26a0\ufe0f Load metadata before creating groups.")
            return

        dlg = GroupBuilderDialog(self.metadata_loader, parent=self)
        dlg.group_created.connect(self._on_group_created)
        dlg.group_removed.connect(self._on_group_removed)
        dlg.exec_()

    def _on_group_created(self, merged_name: str):
        """Slot: rebuild feature extraction tabs after a new group is added."""
        self._log(
            f"\U0001f9ec New cell-type group '{merged_name}' added \u2014 "
            "rebuilding feature extraction tabs."
        )
        # Notify other tabs (segmentation / tracking / filtering) when wired.
        loader = self.metadata_loader
        if loader is not None and hasattr(loader, "metadata_loaded"):
            try:
                loader.metadata_loaded.emit(loader.metadata)
            except Exception:
                pass
        else:
            self._rebuild_tabs()

    def _on_group_removed(self, merged_name: str):
        """Slot: rebuild feature extraction tabs after a *_merged group is removed."""
        self._log(
            f"\U0001f9f9 Cell-type group '{merged_name}' removed \u2014 "
            "rebuilding feature extraction tabs."
        )
        loader = self.metadata_loader
        if loader is not None and hasattr(loader, "metadata_loaded"):
            try:
                loader.metadata_loaded.emit(loader.metadata)
            except Exception:
                pass
        else:
            self._rebuild_tabs()

    def set_queue_panel(self, queue_panel):
        """Attach the processing queue so Active Killing can enqueue itself."""
        self._queue_panel = queue_panel
        if self.active_killing_panel is not None:
            self.active_killing_panel.set_queue_callback(self._queue_active_killing)

    def _queue_active_killing(self):
        """Enqueue the currently selected Active Killing analysis."""
        if self._queue_panel is None or self.active_killing_panel is None:
            return
        from behav3d.napari._queue import StepType
        params = self.active_killing_panel.get_queue_params()
        self._queue_panel.add_step(StepType.ACTIVE_KILLING, params=params)

    def _set_ak_collapsible_state(self, checked: bool):
        """Keep Active Killing toggle text and body visibility in sync."""
        self._ak_body.setVisible(bool(checked))
        arrow = "\u25bc" if checked else "\u25b6"
        self._ak_toggle_btn.setText(f"{arrow}  {self._ak_toggle_label_base}")

    # ── Global (Organoid) Dead Threshold Preview ─────────────────────────────
    def _on_ak_toggle(self, checked: bool):
        """Toggle body visibility of the Active Killing collapsible section."""
        self._set_ak_collapsible_state(checked)

    # ── Shared organoid preview helpers ─────────────────────────────────────
    def _current_viewer_frame(self) -> int:
        viewer = self.viewer
        if viewer is None:
            return 0
        try:
            return int(viewer.dims.current_step[0])
        except Exception:
            return 0

    def _disconnect_org_preview_dims(self):
        viewer = self.viewer
        if viewer is None or self._org_preview_dims_callback is None:
            self._org_preview_dims_callback = None
            return
        try:
            viewer.dims.events.current_step.disconnect(self._org_preview_dims_callback)
        except Exception:
            pass
        self._org_preview_dims_callback = None

    def _connect_org_preview_dims(self):
        """Attach a single shared dims listener that triggers an on-demand
        recompute of the Dead/Alive overlay (current frame only) for the
        organoid preview path. Safe to call multiple times — any previous
        listener is dropped first."""
        viewer = self.viewer
        if viewer is None:
            return
        self._disconnect_org_preview_dims()

        def _on_step(*_):
            cache = self._org_preview_cache
            if (
                cache.get("seg_t") is None
                or cache.get("dead_t") is None
                or not cache.get("is_timelapse")
            ):
                return
            self._refresh_org_preview_for_current_frame()

        try:
            viewer.dims.events.current_step.connect(_on_step)
            self._org_preview_dims_callback = _on_step
        except Exception:
            self._org_preview_dims_callback = None

    def _invalidate_org_preview_cache(self):
        """Drop per-frame caches and zero the shared overlay (e.g. on
        threshold change), so the layer doesn't show stale colours for
        previously-visited frames at the old threshold."""
        cache = self._org_preview_cache
        cache["computed_frames"] = set()
        cache["pct_maps_by_frame"] = {}
        cache["stats_by_frame"] = {}
        overlay_arr = cache.get("overlay_arr")
        if overlay_arr is not None:
            overlay_arr.fill(0)
        # Mirror invalidation onto every organoid panel
        for ct, panel in self.panels.items():
            if ct in self._org_types:
                panel._preview_computed_frames = cache["computed_frames"]
                panel._preview_pct_maps_by_frame = cache["pct_maps_by_frame"]
                panel._preview_stats_by_frame = cache["stats_by_frame"]

    def _refresh_org_preview_for_current_frame(self, thr: float | None = None):
        """Recompute the shared organoid Dead/Alive overlay for the
        currently displayed timepoint (only). Propagates the resulting
        per-frame caches to every organoid panel so hover tooltips stay
        in sync."""
        viewer = self.viewer
        cache = self._org_preview_cache
        if (
            viewer is None
            or cache.get("seg_t") is None
            or cache.get("dead_t") is None
        ):
            return

        if thr is None:
            # Use any organoid panel's current spinner value
            thr_val = None
            for ct, panel in self.panels.items():
                if ct in self._org_types and panel.spin_dead_threshold is not None:
                    thr_val = round(float(panel.spin_dead_threshold.value()), 2)
                    break
            if thr_val is None:
                return
            thr = thr_val
        else:
            thr = round(float(thr), 2)

        # Threshold changed -> invalidate all cached frames
        cached_thr = cache.get("current_thr")
        if cached_thr is None or thr != cached_thr:
            self._invalidate_org_preview_cache()
            cache["current_thr"] = thr

        is_timelapse = bool(cache.get("is_timelapse"))
        frame_idx = self._current_viewer_frame() if is_timelapse else 0

        # Fast path: already computed for this frame at the current threshold
        if frame_idx in cache["computed_frames"]:
            try:
                viewer.layers[f"{_PREVIEW_PREFIX} Dead/Alive"].refresh()
            except Exception:
                pass
            self._propagate_org_preview_to_panels()
            return

        overlay_arr, frame_stats = _update_dead_alive_overlay(
            viewer, cache["seg_t"], cache["dead_t"], thr,
            frame_idx=frame_idx if is_timelapse else None,
            overlay_arr=cache.get("overlay_arr") if is_timelapse else None,
        )
        if is_timelapse:
            cache["overlay_arr"] = overlay_arr

        cache["stats_by_frame"][frame_idx] = frame_stats
        cache["pct_maps_by_frame"][frame_idx] = _build_dead_pct_map(frame_stats)
        cache["computed_frames"].add(frame_idx)

        self._propagate_org_preview_to_panels()

    def _propagate_org_preview_to_panels(self):
        """Mirror the shared cache state onto every organoid panel and
        re-attach the hover tooltip so they stay in sync."""
        cache = self._org_preview_cache
        label_type_map = cache.get("label_type_map") or {}
        for ct, panel in self.panels.items():
            if ct in self._org_types and panel._preview_seg_t is not None:
                panel._preview_overlay_arr = cache.get("overlay_arr")
                panel._preview_stats_by_frame = cache["stats_by_frame"]
                panel._preview_pct_maps_by_frame = cache["pct_maps_by_frame"]
                panel._preview_computed_frames = cache["computed_frames"]
                panel._preview_current_thr = cache.get("current_thr")
                panel._preview_label_type_map = label_type_map
                panel._preview_is_timelapse = bool(cache.get("is_timelapse"))
                panel._attach_preview_dead_hover()

    def _notify_organoid_threshold_changed(self, source_ct: str, value: float):
        """Called by an organoid panel when its threshold spinner changes.

        1. Syncs all other organoid panels' spinners.
        2. Triggers an on-demand single-frame Dead/Alive overlay refresh.
        """
        for ct, panel in self.panels.items():
            if ct != source_ct and ct in self._org_types:
                if panel.spin_dead_threshold is not None:
                    panel.spin_dead_threshold.blockSignals(True)
                    panel.spin_dead_threshold.setValue(value)
                    panel.spin_dead_threshold.blockSignals(False)

        viewer = self.viewer
        cache = self._org_preview_cache
        if (
            viewer is not None
            and cache.get("seg_t") is not None
            and cache.get("dead_t") is not None
        ):
            self._refresh_org_preview_for_current_frame(float(value))

    def _on_run_batch_clicked(self):
        """User-triggered batch run — asynchronous."""
        self.run_batch_feature_extraction(interactive=True, block=False)

    def run_batch_feature_extraction(self, interactive=True, skip_existing=False,
                                     block=True, extra_callbacks=None):
        """Run feature extraction for all cell types sequentially.

        ``block=True`` (default, queue) runs synchronously.  ``block=False``
        (GUI) moves the per-cell-type loop to a background worker and
        shows a determinate progress bar scaled at cell-type granularity
        with sample-level updates rendered as labels.

        ``extra_callbacks`` is the queue's chaining hook
        (``{"on_done": cb, "on_failed": cb}``).
        """
        from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
        from qtpy.QtWidgets import QMessageBox as _QMB

        if not self.panels:
            self._log("No cell type panels available.")
            fire_extra_callback(extra_callbacks, "on_failed", "no cell type panels")
            return

        if not block and self._bg.is_running():
            self._log("⚠️ A batch feature-extraction run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        # Persist global organoid threshold before any panel runs
        self._sync_global_threshold_to_params()

        if self.active_killing_panel is not None:
            self.active_killing_panel.set_queue_callback(
                self._queue_active_killing if self._queue_panel is not None else None
            )

        total = len(self.panels)
        self._log(f"Starting batch feature extraction for {total} cell type(s)…")

        all_cts = list(self.panels.keys())
        existing = []
        existing_cts = set()
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in all_cts:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            combined = feat_dir / f"BEHAV3D_{ct}_combined_track_features.csv"
            if combined.exists():
                existing.append(f"{ct} feature data ({combined.name})")
                existing_cts.add(ct)

        changed_cts = [
            ct for ct in all_cts
            if ct in existing_cts and self.panels[ct]._threshold_changed()
        ]

        skip_existing_flag = skip_existing
        overwrite = not skip_existing
        death_only_cts: set[str] = set()
        if existing:
            if interactive:
                extra = None
                if changed_cts:
                    extra = [
                        (
                            "Re-run death only (changed thresholds)",
                            "death_only",
                            _QMB.ActionRole,
                        ),
                    ]
                choice = prompt_overwrite_batch(
                    self,
                    "Overwrite Existing Features?",
                    existing,
                    extra_buttons=extra,
                )
                if choice == "cancel":
                    self._log("Batch feature extraction cancelled.")
                    fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                    return
                if choice == "death_only":
                    death_only_cts = set(changed_cts)
                    skip_existing_flag = False
                    overwrite = False
                else:
                    skip_existing_flag = choice == "skip"
                    overwrite = not skip_existing_flag
            else:
                overwrite = True

        # Persist + snapshot every panel's Qt widget state on the Qt
        # thread — the worker must not read widgets.
        panel_params: dict[str, dict] = {}
        for ct, panel in self.panels.items():
            panel._persist()
            panel_params[ct] = panel._collect_params()

        SCALE = 100  # sub-units per cell-type step (for ETA granularity)

        def _do_batch(progress_cb=None):
            for i, (ct, panel) in enumerate(self.panels.items()):
                step_label_prefix = f"[{i + 1}/{total}] {ct}"

                def _make_step_cb(step_idx, label_prefix):
                    def _cb(curr, sub_total, label):
                        if progress_cb is None:
                            return
                        if sub_total and sub_total > 0:
                            frac = min(max(curr / sub_total, 0.0), 1.0)
                        else:
                            frac = 0.0
                        try:
                            progress_cb(
                                step_idx * SCALE + int(frac * SCALE),
                                total * SCALE,
                                f"{label_prefix}: {label}" if label else label_prefix,
                            )
                        except Exception:
                            pass
                    return _cb

                if death_only_cts:
                    if ct in death_only_cts:
                        self._log(f"--- [{i + 1}/{total}] Death-only re-run: {ct} ---")
                        panel._run_death_only_for(ct, params=panel_params[ct])
                        self._log(f"Done (death only): {ct}")
                    else:
                        self._log(
                            f"--- [{i + 1}/{total}] Skipping {ct} "
                            "(threshold unchanged) ---"
                        )
                    if progress_cb is not None:
                        try:
                            progress_cb((i + 1) * SCALE, total * SCALE, step_label_prefix)
                        except Exception:
                            pass
                    continue
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- [{i + 1}/{total}] Skipping {ct} (existing data) ---")
                    if progress_cb is not None:
                        try:
                            progress_cb((i + 1) * SCALE, total * SCALE, step_label_prefix)
                        except Exception:
                            pass
                    continue
                self._log(f"--- [{i + 1}/{total}] Feature extraction: {ct} ---")
                panel._run_feature_extraction_for(
                    ct,
                    overwrite=overwrite,
                    params=panel_params[ct],
                    progress_cb=_make_step_cb(i, step_label_prefix),
                )
                self._log(f"Done: {ct}")

            self._log("✅ Batch feature extraction finished.")

        if block:
            try:
                _do_batch(progress_cb=None)
                if interactive:
                    _notify_post_extraction(self)
                fire_extra_callback(extra_callbacks, "on_done", None)
            except Exception as e:
                traceback.print_exc()
                self._log(f"❌ Batch feature extraction error: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
            finally:
                try:
                    self.results_panel.refresh()
                except Exception:
                    traceback.print_exc()
            return

        def _on_done(result):
            if interactive:
                _notify_post_extraction(self)
            try:
                self.results_panel.refresh()
            except Exception:
                traceback.print_exc()
            fire_extra_callback(extra_callbacks, "on_done", result)

        def _on_failed(err: str):
            self._log(f"❌ Batch feature extraction error: {err}")
            try:
                self.results_panel.refresh()
            except Exception:
                traceback.print_exc()
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_batch,
            desc=f"Batch feature extraction ({total} cell types)\u2026",
            progress_row=self.progress_row,
            buttons=[self.btn_run_batch],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _sync_global_threshold_to_params(self):
        """Write the organoid dead threshold into behav3d_parameters
        for every organoid cell type, then save to YAML."""
        # Read from first available organoid panel spinner
        thr = 0.0
        for ct in self._org_types:
            if ct in self.panels and self.panels[ct].spin_dead_threshold is not None:
                thr = float(self.panels[ct].spin_dead_threshold.value())
                break
        thr_val = thr if thr > 0 else None
        params = self.metadata_loader.behav3d_parameters
        features = params.setdefault("features", {})
        for ct in self._org_types:
            ct_cfg = features.setdefault(ct, {})
            ct_cfg["dead_mask_percentage_threshold"] = thr_val

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self._log(f"Warning: Could not save parameters: {e}")

    def get_queue_params(self) -> dict:
        """Collect feature extraction params for all panels (used by queue)."""
        # Make sure organoid threshold is synced before snapshotting
        self._sync_global_threshold_to_params()
        return {ct: panel._collect_params() for ct, panel in self.panels.items()}
