from pathlib import Path

import numpy as np
import pandas as pd


def _resolve_output_dir(output_dir):
    if output_dir is None:
        raise ValueError("output_dir is required.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return output_dir_path


def resolve_exemplar_positions_csv_path(output_dir, cell_type):
    """Resolve canonical track-features CSV path used for exemplar coordinate enrichment."""
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")

    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    feature_outdir = analysis_outdir / "track_features"
    filtered_path = feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"
    combined_path = feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv"
    if filtered_path.exists():
        return filtered_path
    if combined_path.exists():
        return combined_path
    raise FileNotFoundError(
        "Could not find track-features CSV required for exemplar coordinate enrichment. "
        f"Expected one of: '{filtered_path}' or '{combined_path}'."
    )


def merge_coordinate_columns_into_obs(
    adata,
    positions_csv_path,
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
):
    """Merge missing coordinate columns from a positions CSV into adata.obs."""
    key_cols = [str(sample_col), str(track_col), str(time_col)]
    missing_keys_obs = [c for c in key_cols if c not in adata.obs.columns]
    if len(missing_keys_obs) > 0:
        raise ValueError(
            f"Cannot enrich coordinates: adata.obs missing merge key columns {missing_keys_obs}."
        )

    df_pos = pd.read_csv(Path(positions_csv_path), low_memory=False)
    missing_keys_pos = [c for c in key_cols if c not in df_pos.columns]
    if len(missing_keys_pos) > 0:
        raise ValueError(
            "Cannot enrich coordinates: positions CSV missing merge key columns "
            f"{missing_keys_pos}. csv='{positions_csv_path}'"
        )

    candidate_coord_cols = [
        "pixel_position_x",
        "pixel_position_y",
        "pixel_position_z",
        "position_x",
        "position_y",
        "position_z",
    ]
    coord_cols_present = [c for c in candidate_coord_cols if c in df_pos.columns]
    if len(coord_cols_present) == 0:
        raise ValueError(
            "Cannot enrich coordinates: positions CSV has none of the expected coordinate columns "
            "['pixel_position_x','pixel_position_y','pixel_position_z','position_x','position_y','position_z']. "
            f"csv='{positions_csv_path}'"
        )

    missing_coord_cols_in_obs = [c for c in coord_cols_present if c not in adata.obs.columns]
    if len(missing_coord_cols_in_obs) == 0:
        return {
            "csv_path": str(positions_csv_path),
            "added_columns": [],
        }

    merge_key_cols = ["__k_sample", "__k_track", "__k_time"]
    obs = adata.obs.copy()
    obs["__k_sample"] = obs[str(sample_col)].astype("string")
    obs["__k_track"] = obs[str(track_col)].astype("string")
    obs["__k_time"] = pd.to_numeric(obs[str(time_col)], errors="coerce")

    pos = df_pos[key_cols + coord_cols_present].copy()
    pos["__k_sample"] = pos[str(sample_col)].astype("string")
    pos["__k_track"] = pos[str(track_col)].astype("string")
    pos["__k_time"] = pd.to_numeric(pos[str(time_col)], errors="coerce")

    pos = pos.dropna(subset=merge_key_cols)
    pos = (
        pos[merge_key_cols + coord_cols_present]
        .groupby(merge_key_cols, as_index=False, observed=False)
        .first()
    )

    obs["__orig_order"] = np.arange(len(obs), dtype=int)
    merged = obs.merge(
        pos[merge_key_cols + missing_coord_cols_in_obs],
        on=merge_key_cols,
        how="left",
        sort=False,
    )
    merged = merged.sort_values("__orig_order")
    merged.index = adata.obs.index

    drop_cols = ["__k_sample", "__k_track", "__k_time", "__orig_order"]
    merged = merged.drop(columns=[c for c in drop_cols if c in merged.columns])
    adata.obs = merged

    return {
        "csv_path": str(positions_csv_path),
        "added_columns": [str(c) for c in missing_coord_cols_in_obs],
    }


def ensure_exemplar_coordinate_columns(
    adata,
    *,
    output_dir,
    cell_type,
    require_pixel_for_video=False,
):
    """
    Ensure exemplar rendering prerequisites:
      - PDF trajectory panel: either position_* or pixel_position_* triplet.
      - Video panel: pixel_position_* triplet.
    """
    pos_triplet = ["position_x", "position_y", "position_z"]
    pix_triplet = ["pixel_position_x", "pixel_position_y", "pixel_position_z"]

    has_pos = all(c in adata.obs.columns for c in pos_triplet)
    has_pix = all(c in adata.obs.columns for c in pix_triplet)

    needs_enrichment = (not has_pos and not has_pix) or (bool(require_pixel_for_video) and not has_pix)
    merge_info = {"csv_path": None, "added_columns": []}
    if needs_enrichment:
        csv_path = resolve_exemplar_positions_csv_path(output_dir=output_dir, cell_type=cell_type)
        merge_info = merge_coordinate_columns_into_obs(
            adata=adata,
            positions_csv_path=csv_path,
            sample_col="sample_name",
            track_col="TrackID",
            time_col="position_t",
        )

    has_pos = all(c in adata.obs.columns for c in pos_triplet)
    has_pix = all(c in adata.obs.columns for c in pix_triplet)

    missing_pos = [c for c in pos_triplet if c not in adata.obs.columns]
    missing_pix = [c for c in pix_triplet if c not in adata.obs.columns]

    if not (has_pos or has_pix):
        raise ValueError(
            "Exemplar PDF export requires either coordinate triplet "
            "['position_x','position_y','position_z'] or "
            "['pixel_position_x','pixel_position_y','pixel_position_z'] in adata.obs. "
            f"Missing position triplet columns: {missing_pos}; missing pixel triplet columns: {missing_pix}. "
            f"Enrichment csv={merge_info.get('csv_path')}."
        )
    if bool(require_pixel_for_video) and (not has_pix):
        missing_pixel = [c for c in pix_triplet if c not in adata.obs.columns]
        raise ValueError(
            "Exemplar backprojection video export requires pixel coordinates "
            "['pixel_position_x','pixel_position_y','pixel_position_z'] in adata.obs. "
            "PDF export can still run with either position_* or pixel_position_* coordinates. "
            f"Missing pixel columns: {missing_pixel}. Enrichment csv={merge_info.get('csv_path')}."
        )

    return {
        "enriched": bool(needs_enrichment),
        "csv_path": merge_info.get("csv_path"),
        "added_columns": list(merge_info.get("added_columns", [])),
        "has_position_triplet": bool(has_pos),
        "has_pixel_triplet": bool(has_pix),
    }
