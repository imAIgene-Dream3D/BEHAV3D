"""Cell-type grouping: merge already-filtered populations for analysis.

Grouping lives entirely after Filtering, in the Analysis step: it merges
several cell types' *filtered* track-features CSVs into one pseudo
cell-type population (``{name}_merged``) for use in Death Dynamics / Single
Cell. Group membership is recorded in ``behav3d_parameters.yml`` under
``cell_type_groups`` — metadata.csv is never touched, so multicolor
merged-tracking outputs (which must still be processed through Feature
Extraction + Filtering as their own real cell types) are unaffected.

Optionally, :func:`create_group_tracked_segments` builds per-sample tracked
label images + tracks CSVs for a group (for Visualization/backprojection)
by relabeling each member's existing tracked images to the group's merged
TrackIDs. These live under deterministic paths (``group_tracked_image_path``
/ ``group_tracked_csv_path``); callers that want a group to appear like a
real cell type in the viewer resolve those paths on demand rather than
writing them into metadata.csv.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

GROUP_SUFFIX = "_merged"
_PREFIXES = ("or", "im", "ot")


def group_cell_type_id(group_name: str) -> str:
    """Return the pseudo cell-type id used for a group's outputs.

    e.g. ``"tcells"`` -> ``"tcells_merged"``. Idempotent if already suffixed.
    """
    name = str(group_name).strip()
    return name if name.endswith(GROUP_SUFFIX) else f"{name}{GROUP_SUFFIX}"


def list_cell_type_groups(behav3d_parameters: dict) -> dict:
    """Return ``{group_id: [member_cell_types]}`` from config."""
    return dict((behav3d_parameters or {}).get("cell_type_groups", {}) or {})


def validate_group_name(group_name: str, existing_types, existing_groups: dict):
    """Return ``(ok, error_message)`` for a candidate group name."""
    name = (group_name or "").strip()
    if not name:
        return False, "Group name cannot be empty."
    if not all(c.isalnum() or c in "_-" for c in name):
        return False, (
            "Group name can only contain alphanumeric characters, hyphens, "
            "and underscores."
        )
    group_id = group_cell_type_id(name)
    if group_id in (existing_groups or {}):
        return False, f"Group '{group_id}' already exists."
    if group_id in (existing_types or []) or name in (existing_types or []):
        return False, f"Group name '{name}' conflicts with an existing cell type."
    return True, ""


def filtered_track_features_csv(output_dir, cell_type: str) -> Path:
    return (
        Path(output_dir) / "analysis" / cell_type / "track_features"
        / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"
    )


def merged_track_features_csv(output_dir, group_id: str) -> Path:
    return (
        Path(output_dir) / "analysis" / group_id / "track_features"
        / f"BEHAV3D_{group_id}_combined_track_features_filtered.csv"
    )


def group_tracked_image_path(output_dir, sample_name: str, group_id: str) -> Path:
    return (
        Path(output_dir) / "images" / sample_name
        / f"{sample_name}_{group_id}_tracked.zarr"
    )


def group_tracked_csv_path(output_dir, sample_name: str, group_id: str) -> Path:
    return (
        Path(output_dir) / "trackdata" / sample_name / group_id
        / f"{sample_name}_{group_id}_tracks.csv"
    )


def _member_tracked_image_path(metadata_row, cell_type: str):
    """Return the tracked-image path for ``cell_type`` on a metadata row.

    Tries each category prefix (``or``/``im``/``ot``) since a group's
    members can span categories.
    """
    for prefix in _PREFIXES:
        col = f"{prefix}_{cell_type}_tracks_image_path"
        if col in metadata_row.index:
            val = metadata_row[col]
            if pd.notna(val) and str(val).strip():
                return str(val).strip()
    return None


def create_cell_type_group(output_dir, group_name: str, member_cell_types) -> Path:
    """Merge each member's filtered track-features CSV into one population.

    Adds ``origin_cell_type`` (which member each row came from) and
    ``origin_TrackID`` (the member's original TrackID, kept for
    traceability), then reassigns ``TrackID`` so members can never collide
    even if their original numbering overlapped for the same sample —
    downstream code groups strictly by ``(TrackID, sample_name)``.

    Raises ``FileNotFoundError`` naming every member missing a filtered CSV
    (i.e. not yet run through Filtering), and ``ValueError`` if no member
    cell types are given.
    """
    output_dir = Path(output_dir)
    group_id = group_cell_type_id(group_name)

    if not member_cell_types:
        raise ValueError("No member cell types selected.")

    missing = []
    frames = []
    for ct in member_cell_types:
        csv_path = filtered_track_features_csv(output_dir, ct)
        if not csv_path.exists():
            missing.append(ct)
            continue
        df = pd.read_csv(csv_path, low_memory=False)
        if "TrackID" not in df.columns:
            missing.append(f"{ct} (missing TrackID column)")
            continue
        df["origin_cell_type"] = ct
        df["origin_TrackID"] = df["TrackID"]
        frames.append(df)

    if missing:
        raise FileNotFoundError(
            "Missing filtered track-features CSV for: " + ", ".join(missing)
            + ". Run Filtering for these cell types first."
        )

    # Reassign TrackID with a globally unique running counter across all
    # members so two members that happen to reuse the same TrackID (very
    # common, since each cell type is tracked independently) never collide
    # once concatenated into a single (TrackID, sample_name)-keyed table.
    next_id = 1
    for df in frames:
        unique_ids = sorted(df["TrackID"].unique())
        remap = {old: new for new, old in enumerate(unique_ids, start=next_id)}
        df["TrackID"] = df["TrackID"].map(remap)
        next_id = max(remap.values(), default=next_id - 1) + 1

    merged = pd.concat(frames, ignore_index=True, sort=False)

    out_csv = merged_track_features_csv(output_dir, group_id)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)
    return out_csv


def create_group_tracked_segments(output_dir, group_id: str, metadata, log=None, progress_cb=None) -> dict:
    """Build per-sample tracked label images + tracks CSVs for a group.

    Reuses each member cell type's existing tracked zarr (label value ==
    that member's ``TrackID``) together with the group's merged/filtered
    CSV — which records, for every surviving track, its
    ``origin_cell_type``, ``origin_TrackID`` and new group ``TrackID`` — to
    relabel and combine members into one tracked zarr per sample, then
    derives a matching tracks CSV from it (same schema as
    ``convert_tracked_image_to_csv``). This lets the group be visualized
    and backprojected like a real cell type, without ever touching
    metadata.csv: callers resolve the output paths on demand via
    :func:`group_tracked_image_path` / :func:`group_tracked_csv_path`.

    ``progress_cb(current, total, label)``, if given, is called once before
    each sample starts (``total`` = sample count in the merged CSV) and
    once more at completion — matching the ``BackgroundOperation``
    convention in ``behav3d.napari._background_runner`` so this can be run
    off the UI thread with a live progress bar.

    Returns ``{sample_name: {"image": Path, "csv": Path}}`` for samples
    that produced output; samples missing every member's tracked image are
    skipped (logged via ``log`` if given).
    """
    from behav3d.io.images import load_image, append_to_zarr
    from behav3d.preprocessing.tracking import convert_tracked_image_to_csv

    def _log(msg):
        if log:
            log(msg)

    def _progress(current, total, label):
        if progress_cb:
            try:
                progress_cb(current, total, label)
            except Exception:
                pass

    output_dir = Path(output_dir)
    merged_csv = merged_track_features_csv(output_dir, group_id)
    if not merged_csv.exists():
        raise FileNotFoundError(
            f"Merged CSV not found for group '{group_id}': {merged_csv}. "
            "Create the group first."
        )

    df = pd.read_csv(merged_csv, low_memory=False)
    required = {"TrackID", "origin_cell_type", "origin_TrackID", "sample_name", "position_t"}
    missing_cols = required - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Merged CSV for '{group_id}' is missing columns: {sorted(missing_cols)}"
        )

    written = {}
    sample_groups = list(df.groupby("sample_name", sort=False))
    n_samples = len(sample_groups)
    for sample_idx, (sample_name, sdf) in enumerate(sample_groups):
        _progress(sample_idx, n_samples, f"{sample_name} ({sample_idx + 1}/{n_samples})")
        sample_rows = metadata.loc[metadata["sample_name"].astype(str) == str(sample_name)]
        if sample_rows.empty:
            _log(f"  - Skipping {sample_name}: not found in metadata")
            continue
        row = sample_rows.iloc[0]

        member_images = {}
        n_timepoints = None
        for cell_type in sorted(sdf["origin_cell_type"].astype(str).unique()):
            img_path = _member_tracked_image_path(row, cell_type)
            if not img_path or not Path(img_path).exists():
                _log(f"  - {sample_name}/{cell_type}: tracked image not found, skipping member")
                continue
            arr = load_image(img_path)
            member_images[cell_type] = arr
            n_timepoints = len(arr) if n_timepoints is None else min(n_timepoints, len(arr))

        if not member_images:
            _log(f"  - Skipping {sample_name}: no member tracked images available")
            continue

        img_out = group_tracked_image_path(output_dir, sample_name, group_id)
        csv_out = group_tracked_csv_path(output_dir, sample_name, group_id)
        img_out.parent.mkdir(parents=True, exist_ok=True)
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        if img_out.exists():
            shutil.rmtree(img_out) if img_out.is_dir() else img_out.unlink()

        # Per member cell type, per timepoint: pairs of (origin_TrackID in
        # that member's tracked image, new group TrackID to relabel it to).
        relabel_by_ct_t: dict[str, dict[int, np.ndarray]] = {}
        for cell_type in member_images:
            ct_df = sdf[sdf["origin_cell_type"].astype(str) == cell_type]
            for t, tdf in ct_df.groupby("position_t"):
                pairs = tdf[["origin_TrackID", "TrackID"]].dropna().to_numpy()
                relabel_by_ct_t.setdefault(cell_type, {})[int(t)] = pairs

        for t in tqdm(range(n_timepoints), desc=f"{sample_name} [{group_id}]"):
            combined = None
            for cell_type, arr in member_images.items():
                t_img = np.asarray(arr[t])
                if combined is None:
                    combined = np.zeros_like(t_img)
                for old_id, new_id in relabel_by_ct_t.get(cell_type, {}).get(t, ()):
                    combined[t_img == old_id] = new_id
            append_to_zarr(img=np.expand_dims(combined, axis=0), outpath=img_out)

        elsize_x = row.get("pixel_distance_xy", 1) or 1
        elsize_y = row.get("pixel_distance_xy", 1) or 1
        elsize_z = row.get("pixel_distance_z", 1) or 1
        convert_tracked_image_to_csv(
            img_path=img_out,
            outpath=csv_out,
            element_size_x=elsize_x,
            element_size_y=elsize_y,
            element_size_z=elsize_z,
        )

        written[str(sample_name)] = {"image": img_out, "csv": csv_out}
        _log(f"  + {sample_name}: {img_out.name}, {csv_out.name}")

    _progress(n_samples, n_samples, "Done")
    return written


def group_tracked_segments_status(output_dir, group_id: str, metadata) -> dict:
    """Return ``{"built": int, "total": int}`` sample counts for a group's
    tracked segments/files (see :func:`create_group_tracked_segments`).

    A sample counts as "built" when both its tracked image and tracks CSV
    already exist on disk. Used by the Group Builder UI (to show/retrofit
    status) and by the Single Cell tab (to gate backprojection features).
    """
    output_dir = Path(output_dir)
    if metadata is None or "sample_name" not in getattr(metadata, "columns", []):
        return {"built": 0, "total": 0}

    sample_names = metadata["sample_name"].astype(str).unique().tolist()
    built = sum(
        1
        for sample_name in sample_names
        if group_tracked_image_path(output_dir, sample_name, group_id).exists()
        and group_tracked_csv_path(output_dir, sample_name, group_id).exists()
    )
    return {"built": built, "total": len(sample_names)}


def remove_cell_type_group(behav3d_parameters: dict, group_id: str) -> None:
    """Drop ``group_id`` from ``cell_type_groups`` in-place.

    The merged CSV on disk is intentionally left untouched (non-destructive);
    only removed from the config so it stops being offered as a population.
    """
    groups = behav3d_parameters.setdefault("cell_type_groups", {})
    groups.pop(group_id, None)


def group_category(behav3d_parameters: dict, metadata, group_id: str) -> str:
    """Best-effort organoid/immune/other category for a group's dropdown icon.

    Uses the category of the group's first member cell type.
    """
    from behav3d.widgets.utils import detect_cell_type_category

    members = list_cell_type_groups(behav3d_parameters).get(group_id, [])
    for member in members:
        try:
            return detect_cell_type_category(member, metadata)
        except Exception:
            continue
    return "immune"
