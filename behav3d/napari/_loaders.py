"""Shared napari-layer loaders used by both the Visualization tab and the
notebook entry point of the manual-correction editor.

These helpers were extracted from
``behav3d.napari._visualization.VisualizationTab`` so that the
``TrackedSegmentEditingPanel`` (notebook) can spin up a fresh napari window
with the same exact layers — raw channels, tracked segments and tracks only,
**without** the regular untracked segments — without duplicating dataset-
loading logic.

All functions accept a napari ``viewer`` and a metadata ``row`` (a
``pd.Series``) and return the list of layer names they added.  Errors are
swallowed and reported through the optional ``log`` callback so that a
missing optional file does not abort the whole load.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Callable, Dict, List, Optional

import dask.array as da
import pandas as pd

from behav3d.core.metadata import is_multicolor_celltype
from behav3d.napari._track_colors import sync_tracks_colors_to_labels


_CHANNEL_COLORS = ["cyan", "yellow", "green", "red", "blue", "magenta"]


def _noop(_msg: str) -> None:
    return None


def detect_cell_type_columns(row: pd.Series) -> Dict[str, str]:
    """Return ``{cell_type_name: prefix}`` for all detected cell types in a row.

    Mirrors ``VisualizationTab._detect_cell_type_columns`` so that the
    notebook helper sees the same set of cell types as the napari plugin.
    """
    pattern = re.compile(r"^(or|im|ot)_(.+?)_line_condition$")
    out: Dict[str, str] = {}
    for col in row.index:
        m = pattern.match(col)
        if m:
            prefix, ct_name = m.group(1), m.group(2)
            out[ct_name] = prefix
            
    merged_pattern = re.compile(r"^(.+?)_(segments_image_path|tracks_image_path|tracks_csv_path)$")
    for col in row.index:
        if col.startswith(("or_", "im_", "ot_")):
            continue
        m = merged_pattern.match(col)
        if m:
            ct_name = m.group(1)
            if ct_name.endswith("_merged") or ct_name.endswith("_grouped"):
                if ct_name not in out:
                    out[ct_name] = ""
                    
    return out


def cell_types_with_tracked_segments(row: pd.Series) -> List[str]:
    """Return cell-type names that have a non-empty ``*_tracks_image_path``.

    Per-channel multicolor types (``*_N_multicolor``) are excluded because
    only the merged/grouped output is editable.  The individual channel
    segments should not be opened in the track editor.
    """
    out: List[str] = []
    for ct_name, prefix in detect_cell_type_columns(row).items():
        if is_multicolor_celltype(ct_name):
            continue  # skip individual channels — only merged/grouped is editable
        col = f"{prefix}_{ct_name}_tracks_image_path" if prefix else f"{ct_name}_tracks_image_path"
        v = row.get(col)
        if pd.notna(v) and str(v).strip():
            p = Path(str(v).strip())
            if p.exists():
                out.append(ct_name)
    return out


def tracked_segments_paths_for(row: pd.Series, cell_type: str) -> Optional[Dict[str, Path]]:
    """Return ``{"zarr": ..., "csv": ..., "prefix": ...}`` for the given
    cell type, or ``None`` if no tracked image is configured."""
    for prefix in ("or", "im", "ot"):
        zcol = f"{prefix}_{cell_type}_tracks_image_path"
        ccol = f"{prefix}_{cell_type}_tracks_csv_path"
        v = row.get(zcol)
        if pd.notna(v) and str(v).strip():
            return {
                "zarr": Path(str(v).strip()),
                "csv": Path(str(row.get(ccol)).strip())
                if pd.notna(row.get(ccol)) and str(row.get(ccol)).strip()
                else None,
                "prefix": prefix,
            }
    return None


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------
def load_raw_into_viewer(
    viewer,
    sample_name: str,
    row: pd.Series,
    output_dir: str = "",
    log: Callable[[str], None] = _noop,
    display_settings: Optional[Dict] = None,
) -> List[str]:
    """Load the raw image as one Image layer per channel.  Returns layer names."""
    from behav3d.io.formats.zarr import load_zarr

    raw_path = str(row.get("raw_image_path", "")).strip()
    if not raw_path:
        log(f"  ⚠️ No raw_image_path for {sample_name}")
        return []

    raw_p = Path(raw_path)
    candidates = [
        raw_p if raw_p.suffix == ".zarr" else None,
        Path(output_dir, f"{sample_name}.zarr") if output_dir else None,
        raw_p.with_suffix(".zarr"),
    ]
    dask_img = None
    for zp in candidates:
        if zp and zp.exists():
            try:
                dask_img = load_zarr(zp)
                log(f"  📂 Loaded raw zarr: {zp}  shape={dask_img.shape}")
                break
            except Exception as exc:
                log(f"  ⚠️ Error loading {zp}: {exc}")
    if dask_img is None and raw_p.exists():
        try:
            from behav3d.io.images import load_image

            img = load_image(raw_p)
            if not isinstance(img, da.Array):
                img = da.from_array(img, chunks=(1,) + img.shape[1:])
            dask_img = img
            log(f"  📂 Loaded raw image: {raw_p}  shape={dask_img.shape}")
        except Exception as exc:
            log(f"  ❌ Could not load raw image: {exc}")
            return []
    if dask_img is None:
        log(f"  ⚠️ Raw image not found: {raw_path}")
        return []

    dim_order = str(row.get("dimension_order", "TCZYX")).strip() or "TCZYX"
    if dim_order != "TCZYX" and len(dim_order) == 5:
        try:
            axes = [dim_order.index(d) for d in "TCZYX"]
            dask_img = da.transpose(dask_img, axes)
        except ValueError:
            pass

    n_channels = dask_img.shape[1] if dask_img.ndim >= 5 else 1
    saved_channels = display_settings or {}
    added: List[str] = []
    for c in range(n_channels):
        ch_data = dask_img[:, c, :, :, :] if dask_img.ndim >= 5 else dask_img
        saved = saved_channels.get(c) or saved_channels.get(str(c))
        color = (saved["colormap"] if saved and "colormap" in saved
                 else _CHANNEL_COLORS[c % len(_CHANNEL_COLORS)])
        name = f"{sample_name} – Ch{c}"
        add_kwargs = dict(
            name=name, colormap=color, blending="additive", visible=True,
        )
        if saved and "contrast_limits" in saved:
            add_kwargs["contrast_limits"] = tuple(saved["contrast_limits"])
        viewer.add_image(ch_data, **add_kwargs)
        log(f"    + Image layer: {name}  ({color})")
        added.append(name)
    return added


def load_tracked_segments_into_viewer(
    viewer,
    sample_name: str,
    row: pd.Series,
    log: Callable[[str], None] = _noop,
    only_cell_type: Optional[str] = None,
    as_numpy: bool = False,
) -> List[str]:
    """Load only the *tracked* segments (one Labels layer per cell type).

    Parameters
    ----------
    only_cell_type:
        If given, restrict the load to the matching cell type — used by the
        manual-correction editor when it wants exactly one layer.
    as_numpy:
        If True, eagerly materialise the labels array as ``np.ndarray``.
        The editor needs this so per-frame mutations work in-place.
    """
    import numpy as np
    from behav3d.io.formats.zarr import load_zarr

    added: List[str] = []
    for ct_name, prefix in detect_cell_type_columns(row).items():
        if only_cell_type is not None and ct_name != only_cell_type:
            continue
        col = f"{prefix}_{ct_name}_tracks_image_path" if prefix else f"{ct_name}_tracks_image_path"
        val = row.get(col)
        if pd.isna(val) or not str(val).strip():
            continue
        p = Path(str(val).strip())
        if not p.exists():
            log(f"    - Skipping {ct_name} tracked segments: file not found ({p})")
            continue
        try:
            if p.suffix == ".zarr":
                seg = load_zarr(p)
            else:
                from behav3d.io.images import load_image

                seg = load_image(p)
                if not isinstance(seg, da.Array):
                    seg = da.from_array(seg, chunks=(1,) + seg.shape[1:])
            if as_numpy:
                log(
                    f"    Materialising tracked segments for editing "
                    f"(shape={tuple(seg.shape)}, dtype={seg.dtype})…"
                )
                seg = np.asarray(seg)
            name = f"{sample_name} – {ct_name} tracked segments"
            viewer.add_labels(seg, name=name, visible=not is_multicolor_celltype(ct_name))
            log(f"    + Tracked Labels layer: {name}")
            added.append(name)
        except Exception as exc:
            log(f"    ⚠️ Could not load tracked segments for {ct_name}: {exc}")
    return added


def load_tracks_into_viewer(
    viewer,
    sample_name: str,
    row: pd.Series,
    visible: bool = False,
    log: Callable[[str], None] = _noop,
    only_cell_type: Optional[str] = None,
) -> List[str]:
    """Load the per-cell-type Tracks layers (CSV-derived)."""
    added: List[str] = []
    for ct_name, prefix in detect_cell_type_columns(row).items():
        if only_cell_type is not None and ct_name != only_cell_type:
            continue
        col = f"{prefix}_{ct_name}_tracks_csv_path" if prefix else f"{ct_name}_tracks_csv_path"
        val = row.get(col)
        if pd.isna(val) or not str(val).strip():
            continue
        p = Path(str(val).strip())
        if not p.exists():
            log(f"    - Skipping {ct_name} tracks: file not found ({p})")
            continue
        try:
            df = pd.read_csv(p)
            required = {"TrackID", "position_t", "pixel_position_y", "pixel_position_x"}
            if not required.issubset(set(df.columns)):
                log(
                    f"    ⚠️ Tracks CSV for {ct_name} missing columns "
                    f"{required - set(df.columns)}"
                )
                continue
            if "pixel_position_z" in df.columns:
                arr = df[
                    ["TrackID", "position_t", "pixel_position_z",
                     "pixel_position_y", "pixel_position_x"]
                ].values
            else:
                arr = df[
                    ["TrackID", "position_t", "pixel_position_y", "pixel_position_x"]
                ].values
            name = f"{sample_name} – {ct_name} tracks"
            viewer.add_tracks(arr, name=name, visible=visible)
            log(f"    + Tracks layer: {name}")
            added.append(name)

            # Color each track the same as its cell's tracked-segments
            # label, instead of napari's default track_id gradient.
            seg_layer_name = f"{sample_name} – {ct_name} tracked segments"
            if seg_layer_name in viewer.layers:
                sync_tracks_colors_to_labels(
                    viewer.layers[name], viewer.layers[seg_layer_name]
                )
        except Exception as exc:
            log(f"    ⚠️ Could not load tracks for {ct_name}: {exc}")
    return added
