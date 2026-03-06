from pathlib import Path
import time
import pandas as pd
import numpy as np
from behav3d.core.utils import format_time

def smooth_value_over_time(
    df, 
    column, 
    rolling_meanspeed_window, 
    min_periods=1,
    groupby=["TrackID", "sample_name"]
    ):
    smoothed_column = df.groupby(groupby)[column].transform(
        lambda x: x.rolling(
            window=rolling_meanspeed_window, 
            min_periods=min_periods
        ).mean()
    )
    return smoothed_column


def summarize_track_features(
    config=None,
    output_dir=None,
    imaris=False,
    cell_type=None
    ):
    """
    This code calculates summarized features (e.g. mean speed of the whole track) 
    for each TrackID for every experiment specified in the provided metadata.csv
    
    Output:
    - A .csv file containing all tracks from all experiments with their track-summarized features
    """
    
    assert config is not None or output_dir is not None, "Either 'config' or 'output_dir' must be supplied"
    
    start_time = time.time()
            
    print(f"--------------- Summarizing track features ---------------")
    
    if output_dir is None:
        output_dir = config['output_dir']
        if "imaris" in config.keys():
            imaris = config["imaris"]
        else:
            imaris = False
            
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    qc_outdir = Path(analysis_outdir, "quality_control")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not qc_outdir.exists():
        qc_outdir.mkdir(parents=True)

    df_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    df_tracks = pd.read_csv(df_tracks_path)

    # Force touching_* columns to string so they are never treated as numeric
    for col in df_tracks.columns:
        if col.startswith("touching_"):
            df_tracks[col] = df_tracks[col].astype(str)

    # Calculate mean values of track features over the whole track
    grouped_df_tracks=df_tracks.groupby(['sample_name','TrackID'])
    df_summarized_tracks = grouped_df_tracks.size().reset_index(name="track_length")
    
    # Select all numeric + boolean columns, excluding group keys
    num_bool_cols = df_tracks.select_dtypes(include=[np.number, bool]).columns.drop(
        ['TrackID'], errors='ignore'
    )
    means = grouped_df_tracks[num_bool_cols].mean().reset_index()
    means = means.rename(columns={col: f"mean_{col}" for col in num_bool_cols})
    df_summarized_tracks = pd.merge(df_summarized_tracks, means, on=['sample_name','TrackID'], how='left')

    # Dies column (if dead column exists)
    if 'dead' in df_tracks.columns:
        df_summarized_tracks['dies'] = grouped_df_tracks['dead'].any().reset_index()["dead"]
    
    # For cumulative/displacement features: max for displacement_from_origin, last for cumulative_displacement
    if 'displacement_from_origin' in df_tracks.columns:
        df_summarized_tracks['displacement_from_origin'] = grouped_df_tracks['displacement_from_origin'].max().reset_index()['displacement_from_origin']
    if 'cumulative_displacement' in df_tracks.columns:
        df_summarized_tracks['cumulative_displacement'] = grouped_df_tracks['cumulative_displacement'].last().reset_index()['cumulative_displacement']
    
    if not imaris:
        # Calculate active contact percentage for all *_contact columns (generic for any cell type)
        # For each contact type, calculates what % of contact time was "active" contact
        # Example: mean contact could be 1% (0.01) while active contact is 100% (1.0) of that contact time
        contact_cols = [col for col in df_tracks.columns if col.endswith('_contact') and not col.startswith('active_')]
        for contact_col in contact_cols:
            active_col = f"active_{contact_col}"
            if active_col in df_tracks.columns:
                def calculate_active_contact_when_contact(group, contact_col=contact_col, active_col=active_col):
                    if group[contact_col].any():
                        return group[group[contact_col]][active_col].mean()
                    else:
                        return 0
                # Use include_groups=False for pandas >= 2.1.0 compatibility
                result = grouped_df_tracks.apply(calculate_active_contact_when_contact, include_groups=False)
                # Extract just the values (in case it returns a DataFrame)
                if isinstance(result, pd.DataFrame):
                    result = result.iloc[:, 0]
                df_summarized_tracks[active_col] = result.values

    # Dynamically detect metadata columns (well, exp_nr, and any *_line_condition columns)
    metadata_cols = ['TrackID', 'sample_name']
    for col in ['well', 'exp_nr']:
        if col in df_tracks.columns:
            metadata_cols.append(col)
    # Add all *_line_condition columns (e.g., or_organoid1_line_condition, im_tcell_line_condition, etc.)
    line_cols = [col for col in df_tracks.columns if col.endswith('_line_condition')]
    metadata_cols.extend(line_cols)
    
    df_trackinfo = df_tracks[metadata_cols].drop_duplicates()
    df_summarized_tracks = pd.merge(df_trackinfo, df_summarized_tracks, how="left")
    # Write the summarized features to a .csv
    summ_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_summarized.csv")
    print(f"- Writing summarized tracks to {summ_tracks_out_path}")
    df_summarized_tracks.to_csv(summ_tracks_out_path, sep=",", index=False)
    
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_summarized_tracks)