import pandas as pd
from pathlib import Path

from behav3d.core.metadata import load_behav3d_metadata


def _resolve_source(source):
    """Resolve a source path to (metadata_df, experiment_root, label).

    `source` is either a BEHAV3D output_dir (containing metadata.csv) or a
    direct path to a metadata csv file (any filename).
    """
    path = Path(source)
    if path.is_dir():
        metadata_path = path / "metadata.csv"
        experiment_root = path
    else:
        metadata_path = path
        experiment_root = path.parent

    if not metadata_path.exists():
        raise FileNotFoundError(f"metadata csv not found: '{metadata_path}' (from source '{source}')")

    metadata_df = load_behav3d_metadata(metadata_path)
    label = experiment_root.name or str(experiment_root)
    return metadata_df, experiment_root, label


def _find_combined_track_features(experiment_root):
    """Return {cell_type: path} for every BEHAV3D_<cell_type>_combined_track_features.csv
    found under experiment_root/analysis/<cell_type>/track_features/."""
    found = {}
    analysis_dir = Path(experiment_root, "analysis")
    if not analysis_dir.is_dir():
        return found

    for ct_dir in analysis_dir.iterdir():
        if not ct_dir.is_dir():
            continue
        candidate = Path(ct_dir, "track_features", f"BEHAV3D_{ct_dir.name}_combined_track_features.csv")
        if candidate.exists():
            found[ct_dir.name] = candidate

    return found


def combine_behav3d_experiments(sources, output_dir, overwrite=False):
    """
    Combine multiple BEHAV3D experiments into a single merged dataset.

    Merges each source's metadata.csv into one metadata.csv under `output_dir`
    (sample_names must be unique across sources; image/track paths are not
    copied, they stay pointing at their original absolute locations). For any
    cell type where every source has already completed Feature Extraction,
    also concatenates their BEHAV3D_<cell_type>_combined_track_features.csv
    files into one under `output_dir`.

    Parameters
    ----------
    sources : list of str/Path
        Each entry is either a BEHAV3D output_dir (containing metadata.csv)
        or a direct path to a metadata csv file.
    output_dir : str/Path
        Directory to write the merged metadata.csv (and merged
        combined_track_features.csv files, where applicable) into.
    overwrite : bool
        If False (default), raise if output_dir/metadata.csv already exists.

    Returns
    -------
    dict with keys:
        metadata : pd.DataFrame            merged metadata
        metadata_path : Path               where the merged metadata was written
        combined_features : dict           {cell_type: Path} of merged feature csvs written
        skipped_cell_types : dict          {cell_type: [labels of sources missing it]}
        warnings : list of str             non-fatal issues found during the merge
    """
    output_dir = Path(output_dir)
    metadata_out_path = output_dir / "metadata.csv"
    if metadata_out_path.exists() and not overwrite:
        raise FileExistsError(
            f"'{metadata_out_path}' already exists. Pass overwrite=True to replace it."
        )

    resolved = [_resolve_source(source) for source in sources]

    # Hard stop on any sample_name collision across sources.
    sample_owners = {}
    for metadata_df, _experiment_root, label in resolved:
        for sample_name in metadata_df["sample_name"]:
            sample_owners.setdefault(sample_name, []).append(label)

    duplicates = {name: labels for name, labels in sample_owners.items() if len(labels) > 1}
    if duplicates:
        lines = [f"  '{name}' appears in: {', '.join(labels)}" for name, labels in duplicates.items()]
        raise ValueError(
            "Sample names must be unique across combined experiments. "
            "Found overlapping sample_name(s):\n" + "\n".join(lines)
        )

    warnings = []

    merged_metadata = pd.concat(
        [metadata_df for metadata_df, _, _ in resolved], ignore_index=True, sort=False
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_metadata.to_csv(metadata_out_path, index=False)

    # Discover per-cell-type combined_track_features.csv files on disk for each source.
    per_source_features = {
        label: _find_combined_track_features(experiment_root)
        for _metadata_df, experiment_root, label in resolved
    }
    all_labels = [label for _, _, label in resolved]
    all_cell_types = sorted({ct for feats in per_source_features.values() for ct in feats})

    combined_features = {}
    skipped_cell_types = {}

    for cell_type in all_cell_types:
        present_labels = [label for label in all_labels if cell_type in per_source_features[label]]
        missing_labels = [label for label in all_labels if label not in present_labels]

        if missing_labels:
            skipped_cell_types[cell_type] = missing_labels
            warnings.append(
                f"Skipped merging '{cell_type}' combined_track_features: "
                f"Feature Extraction has not been run for it in: {', '.join(missing_labels)}"
            )
            continue

        dfs = {}
        for label in present_labels:
            dfs[label] = pd.read_csv(per_source_features[label][cell_type], low_memory=False)

        column_sets = {label: set(df.columns) for label, df in dfs.items()}
        reference_columns = next(iter(column_sets.values()))
        if any(cols != reference_columns for cols in column_sets.values()):
            mismatch_lines = []
            all_columns = set().union(*column_sets.values())
            for label, cols in column_sets.items():
                missing = sorted(all_columns - cols)
                if missing:
                    mismatch_lines.append(f"  {label} is missing: {missing}")
            warnings.append(
                f"'{cell_type}' combined_track_features do not have identical columns "
                f"across experiments (merging anyway, missing values filled with NaN):\n"
                + "\n".join(mismatch_lines)
            )

        merged_features_df = pd.concat(list(dfs.values()), ignore_index=True, sort=False)
        feature_outdir = Path(output_dir, "analysis", cell_type, "track_features")
        feature_outdir.mkdir(parents=True, exist_ok=True)
        feature_out_path = feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv"
        merged_features_df.to_csv(feature_out_path, index=False)
        combined_features[cell_type] = feature_out_path

    for warning in warnings:
        print(f"WARNING: {warning}")

    return {
        "metadata": merged_metadata,
        "metadata_path": metadata_out_path,
        "combined_features": combined_features,
        "skipped_cell_types": skipped_cell_types,
        "warnings": warnings,
    }
