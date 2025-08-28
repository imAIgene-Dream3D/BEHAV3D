from pathlib import Path
import time
import pandas as pd
import numpy as np
from behav3d.utils import format_time

def summarize_track_features(
    config=None,
    output_dir=None,
    imaris=False,
    cell_type="tcell"
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

    # df_summarized_tracks['mean_dead_dye'] = grouped_df_tracks['mean_dead_dye'].mean().reset_index()["mean_dead_dye"]
    # df_summarized_tracks['mean_MSD'] =  grouped_df_tracks['mean_square_displacement'].mean().reset_index()['mean_square_displacement']
    # df_summarized_tracks['mean_speed'] =  grouped_df_tracks['speed'].mean().reset_index()['speed']
    # df_summarized_tracks['mean_organoid_contact'] =  grouped_df_tracks['organoid_contact'].mean().reset_index()['organoid_contact']
    # df_summarized_tracks['mean_tcell_contact'] =  grouped_df_tracks['tcell_contact'].mean().reset_index()['tcell_contact']
    # df_summarized_tracks['mean_displacement'] =  grouped_df_tracks['displacement'].mean().reset_index()['displacement']
    df_summarized_tracks['dies'] =  grouped_df_tracks['dead'].any().reset_index()["dead"]
    
    # For some values, take the maximum of the track such as "displacement_from_origin"
    df_summarized_tracks['displacement_from_origin'] =  grouped_df_tracks['displacement_from_origin'].last().reset_index()['displacement_from_origin']
    df_summarized_tracks['cumulative_displacement'] =  grouped_df_tracks['cumulative_displacement'].last().reset_index()['cumulative_displacement']
    
    if not imaris:
        # Calculate for the contact that occurs, what percentage has been active contact
        # As it only takes points of contact, this can mean the mean contact is 1% (0.01)
        # while the active contact can then still be 100% (1.0)
        def calculate_active_contact_when_contact(group):
            if group['tcell_contact'].any():
                return group[group['tcell_contact']]['active_tcell_contact'].mean()
            else:
                return 0
        df_summarized_tracks['active_tcell_contact'] = grouped_df_tracks.apply(calculate_active_contact_when_contact).reset_index(drop=True)

    df_trackinfo = df_tracks[['TrackID', 'sample_name','well', 'exp_nr', 'organoid_line', 'tcell_line']].drop_duplicates()
    df_summarized_tracks = pd.merge(df_trackinfo, df_summarized_tracks, how="left")
    # Write the summarized features to a .csv
    summ_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_summarized.csv")
    print(f"- Writing summarized tracks to {summ_tracks_out_path}")
    df_summarized_tracks.to_csv(summ_tracks_out_path, sep=",", index=False)
    
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_summarized_tracks)