import pandas as pd
from pathlib import Path
import re


def is_multicolor_celltype(cell_type_name: str) -> bool:
    """Return True if a cell type name matches the per-channel multicolor pattern.

    Pattern used by the metadata builder: <base>_<index>_multicolor (e.g., tcells_1_multicolor)
    """
    if not isinstance(cell_type_name, str):
        return False
    return re.search(r'_[0-9]+_multicolor$', cell_type_name) is not None


def multicolor_base_name(cell_type_name: str) -> str:
    """If `cell_type_name` is a multicolor per-channel name, return the base name.

    Example: tcells_2_multicolor -> tcells
    If not multicolor, return the original name.
    """
    m = re.search(r'^(.*)_([0-9]+)_multicolor$', str(cell_type_name))
    if m:
        return m.group(1)
    return cell_type_name


def is_merged_celltype(cell_type_name: str) -> bool:
    """Return True if the cell type name looks like a merged multicolor output.

    We use the suffix `_merged` for the merged output (e.g., tcells_merged).
    """
    if not isinstance(cell_type_name, str):
        return False
    return cell_type_name.endswith('_merged')


def is_grouped_celltype(cell_type_name: str) -> bool:
    """Return True if the cell type name looks like a grouped multicolor output.

    We also support the suffix `_grouped` (e.g., tcells_grouped).
    """
    if not isinstance(cell_type_name, str):
        return False
    return cell_type_name.endswith('_grouped')


def is_combined_multicolor_celltype(cell_type_name: str) -> bool:
    """Return True for derived multicolor outputs (`*_merged` or `*_grouped`)."""
    return is_merged_celltype(cell_type_name) or is_grouped_celltype(cell_type_name)


def filter_multicolor_inputs(celltype_list):
    """Filter out per-channel multicolor entries from a list of cell type names.

    Keeps merged names (ending with `_merged`) and non-multicolor names.
    """
    out = []
    for ct in celltype_list:
        if is_multicolor_celltype(ct):
            continue
        out.append(ct)
    return out


def detect_merged_cell_types_from_metadata(metadata):
    """Detect derived multicolor outputs from metadata columns.

    Looks for standalone columns like ``tcells_merged_tracks_image_path`` /
    ``tcells_grouped_tracks_image_path`` or matching ``*_tracks_csv_path`` columns
    that are not prefixed with ``or_``, ``im_``, or ``ot_``.
    """
    if metadata is None:
        return []

    merged_types = set()
    suffixes = ("_tracks_image_path", "_tracks_csv_path", "_segments_image_path")
    for col in metadata.columns:
        if not col.endswith(suffixes):
            continue
        if col.startswith(("or_", "im_", "ot_")):
            continue
        base = col
        for suffix in suffixes:
            if base.endswith(suffix):
                base = base[: -len(suffix)]
                break
        if is_combined_multicolor_celltype(base):
            merged_types.add(base)

    return sorted(list(merged_types))


def expand_multicolor_celltype_names(base_name, n_channels):
    """Expand a base cell type into per-channel multicolor names.

    Example: expand_multicolor_celltype_names("tcells", 3) ->
    ["tcells_1_multicolor", "tcells_2_multicolor", "tcells_3_multicolor"]
    """
    base_name = str(base_name).strip()
    n_channels = int(n_channels)
    if not base_name:
        raise ValueError("base_name must be a non-empty string")
    if n_channels < 1:
        raise ValueError("n_channels must be at least 1")
    return [f"{base_name}_{idx}_multicolor" for idx in range(1, n_channels + 1)]

def detect_organoid_types_from_metadata(metadata):
    """
    Detect organoid types from metadata column names.
    Looks for columns matching pattern: or_{organoidtype}_line_condition, or_{organoidtype}_segments_image_path, etc.
    Returns list of organoid type names (e.g., ['organoid1', 'organoid2'])
    """
    if metadata is None:
        return []
    
    organoid_types = set()
    # Pattern: columns starting with 'or_' prefix
    suffixes = ("_line_condition", "_segments_image_path", "_tracks_image_path", "_tracks_csv_path")
    for col in metadata.columns:
        if col.startswith('or_'):
            base = col[3:]
            for suffix in suffixes:
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            if base:
                organoid_types.add(base)
    
    return sorted(list(organoid_types))

def detect_immune_cell_types_from_metadata(metadata):
    """
    Detect immune cell types from metadata column names.
    Looks for columns matching pattern: im_{celltype}_line_condition, im_{celltype}_segments_image_path, etc.
    Returns list of immune cell type names (e.g., ['tcell', 'macro', 'nk'])
    """
    if metadata is None:
        return []
    
    immune_types = set()
    # Pattern: columns starting with 'im_' prefix
    suffixes = ("_line_condition", "_segments_image_path", "_tracks_image_path", "_tracks_csv_path")
    for col in metadata.columns:
        if col.startswith('im_'):
            base = col[3:]
            for suffix in suffixes:
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            if base:
                immune_types.add(base)
    
    return sorted(list(immune_types))

def detect_other_cell_types_from_metadata(metadata):
    """
    Detect 'other' cell types from metadata column names.
    Looks for columns matching pattern: ot_{celltype}_line_condition, ot_{celltype}_segments_image_path, etc.
    Returns list of other cell type names (e.g., ['tumor1', 'fibroblast'])
    """
    if metadata is None:
        return []
    
    other_types = set()
    # Pattern: columns starting with 'ot_' prefix
    suffixes = ("_line_condition", "_segments_image_path", "_tracks_image_path", "_tracks_csv_path")
    for col in metadata.columns:
        if col.startswith('ot_'):
            base = col[3:]
            for suffix in suffixes:
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            if base:
                other_types.add(base)
    
    return sorted(list(other_types))

def has_dead_channel(metadata):
    """
    Check if dead_channel column exists and has non-null values in metadata.
    Returns True if dead channel is present, False otherwise.
    """
    if metadata is None:
        return False
    
    if 'dead_channel' not in metadata.columns:
        return False
    
    # Check if any row has a non-null dead_channel value
    return metadata['dead_channel'].notna().any()

def has_dead_mask(metadata):
    """
    Check if dead_mask_path column exists and has non-null, non-empty values in metadata.
    Returns True if dead mask path is present, False otherwise.
    """
    if metadata is None:
        return False
    
    if 'dead_mask_path' not in metadata.columns:
        return False
    
    # Check if any row has a non-null AND non-empty dead_mask_path value
    # (empty strings '' are not NaN but should be treated as missing)
    mask_values = metadata['dead_mask_path'].fillna('')
    return (mask_values.str.strip() != '').any()

def load_behav3d_metadata(
    metadata_path
    ):
    """
    Load metadata CSV with dynamic cell type support.
    Uses prefixed columns: or_*, im_*, ot_* for organoid, immune, and other cell types.
    """
    # Basic dtype dict for common columns (cell-type specific columns are dynamic)
    dtype_dict = {
        "sample_name": str,
        "exp_nr": "Int64",
        "well": str,
        "number_of_channels": "Int64",
        "dead_channel": "Int64",
        "dead_mask_path": str,
        "pixel_distance_xy": float,
        "pixel_distance_z": float,
        "distance_unit": str,
        "time_interval": float,
        "time_unit": str,
        "dimension_order": str, 
        "raw_image_path": str,
        "signal_unmixing_image_path": str, #FUNC
        "original_raw_image_path": str, # FUNC
    }

    metadata = pd.read_csv(metadata_path)
    
    # Apply dtypes only to columns that exist
    for col, dtype in dtype_dict.items():
        if col in metadata.columns:
            try:
                metadata[col] = metadata[col].astype(dtype)
            except:
                pass  # Skip if conversion fails

    metadata = metadata.dropna(how="all").reset_index(drop=True)
    return metadata 

def check_behav3d_metadata(
    metadata,
    func=False
    ):
    """
    Validate metadata with dynamic cell type support.
    Uses prefixed columns: or_*, im_*, ot_* for organoid, immune, and other cell types.
    """
    # Detect cell types dynamically.
    # Organoid and other types keep the existing schema, while immune types may
    # appear either in legacy form (im_tcell_...) or as multicolor-expanded
    # entries (im_tcell_1_multicolor_...).
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)

    immune_columns = [col for col in metadata.columns if col.startswith("im_")]
    immune_validation_types = sorted(
        {
            col[3 : -len(suffix)]
            for col in immune_columns
            for suffix in ("_line_condition", "_segments_image_path", "_tracks_image_path", "_tracks_csv_path")
            if col.endswith(suffix)
        }
    )
    if immune_validation_types:
        immune_types = immune_validation_types
    
    # Basic required columns (always needed)
    required_columns = [
        "sample_name", 
        "exp_nr", 
        "well", 
        "pixel_distance_xy", 
        "pixel_distance_z", 
        "distance_unit", 
        "time_interval", 
        "time_unit",
        "raw_image_path",
    ]
    
    # Build dynamic required columns for each detected cell type
    path_columns = []  # Columns that can be empty (filled by pipeline)
    
    for org_type in organoid_types:
        if is_combined_multicolor_celltype(org_type):
            continue
        required_columns.append(f"or_{org_type}_line_condition")
        path_columns.extend([
            f"or_{org_type}_segments_image_path",
            f"or_{org_type}_tracks_image_path",
            f"or_{org_type}_tracks_csv_path"
        ])
    
    for immune_type in immune_types:
        if is_combined_multicolor_celltype(immune_type):
            continue
        required_columns.append(f"im_{immune_type}_line_condition")
        path_columns.extend([
            f"im_{immune_type}_segments_image_path",
            f"im_{immune_type}_tracks_image_path",
            f"im_{immune_type}_tracks_csv_path"
        ])
    
    for other_type in other_types:
        if is_combined_multicolor_celltype(other_type):
            continue
        required_columns.append(f"ot_{other_type}_line_condition")
        path_columns.extend([
            f"ot_{other_type}_segments_image_path",
            f"ot_{other_type}_tracks_image_path",
            f"ot_{other_type}_tracks_csv_path"
        ])
    
    # Optional func columns
    func_columns = []
    if func:
        if "signal_unmixing_image_path" in metadata.columns:
            func_columns.append("signal_unmixing_image_path")
        if "original_raw_image_path" in metadata.columns:
            func_columns.append("original_raw_image_path")
    
    # Check required columns exist
    missing_columns = set(required_columns) - set(metadata.columns)
    if missing_columns:
        missing_string = '\n'.join(missing_columns)
        raise AssertionError(f"Not all required columns are present in the metadata .csv file\n{missing_string}")
    
    # Check that non-path columns are filled - collect ALL missing columns
    columns_to_check = [col for col in required_columns if col not in path_columns]
    missing_fields = []
    
    for col in columns_to_check:
        # Check for NaN, empty string, or literal 'nan' string
        is_missing = (
            metadata[col].isna() | 
            (metadata[col].astype(str).str.strip() == '') |
            (metadata[col].astype(str).str.strip().str.lower() == 'nan')
        )
        if is_missing.any():
            # Get row indices where values are missing
            missing_indices = metadata[is_missing].index.tolist()
            # Try to get sample names, fall back to row numbers if sample_name is NaN/empty
            missing_identifiers = []
            for idx in missing_indices:
                sample_name = metadata.loc[idx, 'sample_name']
                # Check for NaN, empty string, or literal 'nan' string
                if pd.notna(sample_name) and str(sample_name).strip() and str(sample_name).strip().lower() != 'nan':
                    missing_identifiers.append(str(sample_name))
                else:
                    missing_identifiers.append(f"Row {idx + 1}")
            missing_fields.append((col, missing_identifiers))
    
    # Check channel columns: if number_of_channels is set, validate channel columns
    # 0-indexed format
    channel_cols = sorted([col for col in metadata.columns if re.match(r'^channel\d+$', col)],
                          key=lambda x: int(re.match(r'^channel(\d+)$', x).group(1)))
    
    # Collect all channel-related errors to show them all at once
    channel_errors = []
    
    # If number_of_channels exists, validate it
    if 'number_of_channels' in metadata.columns:
        # 1. Check that number_of_channels is not empty in any row
        n_channels_missing = (
            metadata['number_of_channels'].isna() | 
            (metadata['number_of_channels'].astype(str).str.strip() == '') |
            (metadata['number_of_channels'].astype(str).str.strip().str.lower() == 'nan')
        )
        if n_channels_missing.any():
            missing_indices = metadata[n_channels_missing].index.tolist()
            missing_identifiers = []
            for idx in missing_indices:
                sample_name = metadata.loc[idx, 'sample_name']
                if pd.notna(sample_name) and str(sample_name).strip() and str(sample_name).strip().lower() != 'nan':
                    missing_identifiers.append(str(sample_name))
                else:
                    missing_identifiers.append(f"Row {idx + 1}")
            channel_errors.append(
                f"• number_of_channels is missing or empty in: {', '.join(missing_identifiers)}"
            )
        
        # Only continue with other checks if number_of_channels is filled in all rows
        if not n_channels_missing.any():
            # 2. Check that number_of_channels is the same in all rows
            n_channels_values = metadata['number_of_channels'].astype(int).unique()
            if len(n_channels_values) > 1:
                channel_errors.append(
                    f"• Different rows have different number_of_channels values: {list(n_channels_values)}. All rows must have the same value."
                )
            
            n_channels = int(n_channels_values[0])
            
            # Build list of valid channel labels from detected cell types FIRST
            # to calculate expected number of channels
            valid_channel_labels = set()
            valid_channel_labels.update([ct for ct in organoid_types if not is_combined_multicolor_celltype(ct)])
            valid_channel_labels.update([ct for ct in immune_types if not is_combined_multicolor_celltype(ct)])
            valid_channel_labels.update([ct for ct in other_types if not is_combined_multicolor_celltype(ct)])
            
            # Check if dead_channel is used (add "dead" to valid labels if any row has dead_channel set)
            has_dead_in_metadata = (
                'dead_channel' in metadata.columns and
                metadata['dead_channel'].notna().any() and
                (metadata['dead_channel'].astype(str).str.strip() != '').any() and
                (metadata['dead_channel'].astype(str).str.strip().str.lower() != 'nan').any()
            )
            if has_dead_in_metadata:
                valid_channel_labels.add('dead')
            
            # Calculate expected number of channels based on cell types
            expected_n_channels = len(valid_channel_labels)
            expected_channel_cols = [f'channel{i}' for i in range(expected_n_channels)]
            
            # 3. Check that number_of_channels matches expected (based on cell types)
            if n_channels != expected_n_channels:
                cell_types_breakdown = []
                if organoid_types:
                    cell_types_breakdown.append(f"{len(organoid_types)} organoid(s): {sorted(organoid_types)}")
                if immune_types:
                    cell_types_breakdown.append(f"{len(immune_types)} immune(s): {sorted(immune_types)}")
                if other_types:
                    cell_types_breakdown.append(f"{len(other_types)} other(s): {sorted(other_types)}")
                if has_dead_in_metadata:
                    cell_types_breakdown.append("1 dead channel")
                channel_errors.append(
                    f"• number_of_channels is {n_channels} but should be {expected_n_channels} based on cell types defined. "
                    f"Cell types: {', '.join(cell_types_breakdown)}"
                )
            
            # 4. Check that number of channel columns matches expected
            if len(channel_cols) != expected_n_channels:
                channel_errors.append(
                    f"• Found {len(channel_cols)} channel columns but expected {expected_n_channels}. "
                    f"Expected columns: {expected_channel_cols}, Found: {channel_cols}"
                )
            
            # Check that channel columns are named correctly (channel0, channel1, ..., channel{n-1})
            if channel_cols != expected_channel_cols and len(channel_cols) == expected_n_channels:
                channel_errors.append(
                    f"• Channel columns do not match expected format. Expected: {expected_channel_cols}, Found: {channel_cols}"
                )
            
            # 4. PER-ROW VALIDATION: Check that labels in each row match the cell type names
            # Labels can be in different order per row, but must match the valid cell types
            for rowidx, sample_metadata in metadata.iterrows():
                sample_name = sample_metadata.get('sample_name', f'Row {rowidx + 1}')
                
                # Get all channel labels for this row (can be in any order)
                row_channel_labels = []
                missing_channel_labels = []
                for ch_col in channel_cols:
                    ch_label = sample_metadata.get(ch_col)
                    is_missing = (
                        pd.isna(ch_label) or
                        str(ch_label).strip() == '' or
                        str(ch_label).strip().lower() == 'nan'
                    )
                    if is_missing:
                        missing_channel_labels.append(ch_col)
                    else:
                        row_channel_labels.append(str(ch_label).strip())
                
                # Check that all channel labels are filled
                if missing_channel_labels:
                    channel_errors.append(
                        f"• Sample '{sample_name}' (Row {rowidx + 1}): Missing channel labels in {missing_channel_labels}"
                    )
                    continue  # Skip further checks for this row if labels are missing
                
                # Check that all labels match the valid cell type names
                invalid_labels = [ch_label for ch_label in row_channel_labels if ch_label not in valid_channel_labels]
                if invalid_labels:
                    channel_errors.append(
                        f"• Sample '{sample_name}' (Row {rowidx + 1}): Invalid channel labels {invalid_labels}. "
                        f"Valid options: {sorted(valid_channel_labels)}"
                    )
                
                # Check that all required cell types are present in the labels
                present_labels = set(row_channel_labels)
                missing_types = valid_channel_labels - present_labels
                if missing_types and not invalid_labels:  # Only report if no invalid labels (to avoid duplicate info)
                    channel_errors.append(
                        f"• Sample '{sample_name}' (Row {rowidx + 1}): Missing required cell types {sorted(missing_types)}. "
                        f"Present labels: {sorted(present_labels)}"
                    )
    
    # Raise all channel errors at once if any
    if channel_errors:
        # Check for dead channel presence for the error message
        has_dead_for_msg = (
            'dead_channel' in metadata.columns and
            metadata['dead_channel'].notna().any() and
            (metadata['dead_channel'].astype(str).str.strip() != '').any() and
            (metadata['dead_channel'].astype(str).str.strip().str.lower() != 'nan').any()
        )
        error_msg = "Channel configuration errors:\n" + "\n".join(channel_errors)
        error_msg += f"\n\nValid cell type labels based on metadata columns:\n"
        error_msg += f"  Organoid types (or_*): {sorted(organoid_types) if organoid_types else 'None'}\n"
        error_msg += f"  Immune types (im_*): {sorted(immune_types) if immune_types else 'None'}\n"
        error_msg += f"  Other types (ot_*): {sorted(other_types) if other_types else 'None'}\n"
        error_msg += f"  Dead channel: {'dead' if has_dead_for_msg else 'N/A'}"
        raise AssertionError(error_msg)
    
    # If there are missing fields, raise error with ALL of them (BLOCKING)
    if missing_fields:
        error_lines = ["⚠️ The following required fields have missing values:\n"]
        for col, identifiers in missing_fields:
            ids_str = ', '.join(identifiers[:3])
            if len(identifiers) > 3:
                ids_str += f" (+{len(identifiers) - 3} more)"
            error_lines.append(f"  • {col}: missing in [{ids_str}]")
        error_lines.append("\nPlease fill in all required fields for all samples.")
        raise AssertionError('\n'.join(error_lines))
    
    
    # Validate paths per row
    ok = True
    all_cell_types = [(f"or_{t}", t) for t in organoid_types if not is_combined_multicolor_celltype(t)] + \
                     [(f"im_{t}", t) for t in immune_types if not is_combined_multicolor_celltype(t)] + \
                     [(f"ot_{t}", t) for t in other_types if not is_combined_multicolor_celltype(t)]
    
    for rowidx, sample_metadata in metadata.iterrows():
        print(f"Row {rowidx+1}: {sample_metadata['sample_name']}")
        sample_name = sample_metadata['sample_name']
        
        # Check raw image exists
        raw_path = str(sample_metadata["raw_image_path"]).strip().strip('"').strip("'")
        assert Path(raw_path).exists(), \
            f"The raw_image_path supplied for 'row {rowidx+1}: {sample_name}' does not exist: '{raw_path}'"
        
        # Check dead_mask_path (required if dead_channel is set)
        has_dead_channel_val = (
            'dead_channel' in sample_metadata.index and 
            pd.notna(sample_metadata['dead_channel']) and
            str(sample_metadata['dead_channel']).strip() and
            str(sample_metadata['dead_channel']).strip().lower() != 'nan'
        )
        
        if 'dead_mask_path' in sample_metadata.index:
            dead_mask_val = sample_metadata['dead_mask_path']
            # Check if it's a valid non-empty value (not NaN, not empty string, not 'nan' string)
            is_valid_value = (
                pd.notna(dead_mask_val) and 
                str(dead_mask_val).strip() and 
                str(dead_mask_val).strip().lower() != 'nan'
            )
            if is_valid_value:
                # Add robustness: strip quotes
                dead_mask_val = str(dead_mask_val).strip().strip('"').strip("'")
                if not Path(dead_mask_val).exists():
                    print(f"!!! dead_mask_path '{dead_mask_val}' does not exist. Please run segmentation below.")
                    ok = False
            elif has_dead_channel_val:
                # dead_channel is set but dead_mask_path is missing
                print(f"!!! No dead mask image. Please run segmentation below.")
                ok = False
        
        # Check each cell type's paths
        for prefix_type, display_type in all_cell_types:
            segments_col = f"{prefix_type}_segments_image_path"
            tracks_img_col = f"{prefix_type}_tracks_image_path"
            tracks_csv_col = f"{prefix_type}_tracks_csv_path"
            
            # Check segments
            if not pd.isna(sample_metadata[segments_col]) and sample_metadata[segments_col]:
                p = str(sample_metadata[segments_col]).strip().strip('"').strip("'")
                if not Path(p).exists():
                    print(f"⚠️ {segments_col} path does not exist: '{p}'")
            else:
                print(f"!!! No segmented {display_type} image. Please run segmentation below.")
                ok = False
            
            # Check tracks image
            if not pd.isna(sample_metadata[tracks_img_col]) and sample_metadata[tracks_img_col]:
                p = str(sample_metadata[tracks_img_col]).strip().strip('"').strip("'")
                if not Path(p).exists():
                    print(f"⚠️ {tracks_img_col} path does not exist: '{p}'")
            else:
                print(f"!!! No tracked {display_type} image. Please run tracking below.")
                ok = False
            
            # Check tracks CSV
            if not pd.isna(sample_metadata[tracks_csv_col]) and sample_metadata[tracks_csv_col]:
                p = str(sample_metadata[tracks_csv_col]).strip().strip('"').strip("'")
                if not Path(p).exists():
                    print(f"⚠️ {tracks_csv_col} path does not exist: '{p}'")
            else:
                print(f"!!! No tracked {display_type} CSV. Please run tracking below.")
                ok = False
        
        # Check func columns if needed
        if func:
            if "signal_unmixing_image_path" in func_columns:
                if not pd.isna(sample_metadata["signal_unmixing_image_path"]):
                    assert Path(sample_metadata["signal_unmixing_image_path"]).exists(), \
                        f"signal_unmixing_image_path for Row {rowidx+1} '{sample_name}' does not exist"
                else:
                    print(f"!!! No signal unmixed image supplied.")
                    ok = False
        
        print("")
    
    if ok:
        print("✅ Metadata file is complete and correct")
    else:
        print("⚠️ Some processing steps need to be run (see warnings above)")


METADATA_TECHNICAL_PATTERNS = (
    "_path", "_dir", "pixel_distance", "time_interval", "unit_",
    "channel_", "dead_channel", "zarr", "dimension_order",
)


def list_grouping_candidate_columns(metadata):
    """Return metadata columns usable as grouping/condition columns, excluding
    technical columns (paths, dirs, pixel/time units, channel/zarr internals)."""
    if metadata is None or not hasattr(metadata, "columns"):
        return []
    exclude = METADATA_TECHNICAL_PATTERNS
    return [c for c in metadata.columns if not any(pat in c.lower() for pat in exclude)]


def merge_condition_columns_into_obs(adata, metadata, columns, key_col="sample_name"):
    """Left-join `columns` from `metadata` into `adata.obs` on `key_col`, for any
    column not already present in `adata.obs`. Returns adata (mutated in place)."""
    if adata is None or not hasattr(adata, "obs"):
        return adata
    if metadata is None or key_col not in getattr(metadata, "columns", []):
        return adata
    cols_to_add = [c for c in columns if c not in adata.obs.columns and c in metadata.columns]
    if not cols_to_add:
        return adata
    meta_subset = metadata[[key_col] + cols_to_add].drop_duplicates(subset=[key_col])
    orig_index = adata.obs.index
    merged = adata.obs.merge(meta_subset, on=key_col, how="left")
    merged.index = orig_index
    adata.obs = merged
    return adata
