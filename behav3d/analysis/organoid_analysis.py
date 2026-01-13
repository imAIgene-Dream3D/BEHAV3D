from pathlib import Path
import time

import pandas as pd
import numpy as np

from behav3d.core.utils import format_time
from behav3d.analysis import smooth_value_over_time
from behav3d.preprocessing.filtering import plot_filter_count, filter_by_full_duration, filter_minimal_track_length, trim_to_maximal_track_length
from behav3d.preprocessing import calc_z_projection
from behav3d.io.images import load_image

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colormaps


def run_organoid_analysis(
    dead_perc_threshold,
    config=None,
    output_dir=None,
    metadata=None,
    df_tracks_path=None,
    org_type="organoid",
    # df_tracks_summarized_path=None,
    ):
    print(f"--------------- Performing {org_type} Death Dynamics analysis ---------------")
    start_time = time.time()
    assert config is not None or all(
        [output_dir]
    ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
        
    
    analysis_outdir = Path(output_dir, "analysis", org_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    results_outdir = Path(analysis_outdir, "results")
    sample_outdir = Path(results_outdir, "per_sample")
    
    if not sample_outdir.exists():
        sample_outdir.mkdir(parents=True)
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)    
    if not results_outdir.exists():
        results_outdir.mkdir(parents=True)    
        
    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
    # if df_tracks_summarized_path is None:
    #     df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_organoid_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path)
    # df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    df_tracks["TrackID"] = df_tracks["TrackID"].astype(str)
    
    # Check if dead channel data exists in tracks
    dead_columns = ["nr_dead_mask_pixels", "percentage_dead_mask", "mean_dead_dye"]
    has_dead_data = all(col in df_tracks.columns for col in dead_columns)
    
    if not has_dead_data:
        print(f"\n⚠️  WARNING: Dead channel data not found in features.")
        print(f"   Missing columns: {[col for col in dead_columns if col not in df_tracks.columns]}")
        print(f"   To run death dynamics analysis, enable dead_channel in MetadataBuilder before running feature extraction.")
        print(f"   Skipping death dynamics analysis.\n")
        
        # Cannot run death dynamics analysis without dead channel data
        end_time = time.time()
        h,m,s = format_time(start_time, end_time)
        print(f"### DONE (no death data) - elapsed time: {h}:{m:02}:{s:02}\n")
        return
    
    # death analysis (only runs if has_dead_data)
    df_tracks["smoothed_nr_dead_mask_pixels"] = smooth_value_over_time(
            df_tracks, 
            column="nr_dead_mask_pixels", 
            rolling_meanspeed_window=20,
            min_periods=20,
            groupby=["TrackID", "sample_name"]
        )
    
    df_tracks["smoothed_percentage_dead_mask"] = smooth_value_over_time(
            df_tracks, 
            column="percentage_dead_mask", 
            rolling_meanspeed_window=20,
            min_periods=20,
            groupby=["TrackID", "sample_name"]
        )
    
    df_tracks["smoothed_mean_dead_dye"] = smooth_value_over_time(
            df_tracks, 
            column="mean_dead_dye", 
            rolling_meanspeed_window=20,
            min_periods=20,
            groupby=["TrackID", "sample_name"]
        )
    
    df_tracks["dead"] = (df_tracks["smoothed_percentage_dead_mask"] > dead_perc_threshold)
    df_tracks["dead"] = df_tracks.groupby(["sample_name","TrackID"])["dead"].transform(lambda x: x.cummax())

    # Summarize tracks to their time of death
    summ = (
        df_tracks
        .groupby(["sample_name","TrackID"])
        .agg(
            t_first=("position_t","min"),
            t_last =("position_t","max"),
            # first time it is ever dead (NaN if never dead)
            t_dead =("position_t", lambda s: s[df_tracks.loc[s.index, "dead"]].min()
                    if (df_tracks.loc[s.index,"dead"]).any() else np.nan)
        )
        .reset_index()
    )

    timeline = df_tracks[["sample_name", "position_t"]].drop_duplicates()
    
    grid = timeline.merge(summ, on="sample_name", how="left")
    t = grid["position_t"]

    dead_at_t = grid["t_dead"].notna() & (t >= grid["t_dead"])
    seen_at_t = (t >= grid["t_first"]) & (t <= grid["t_last"])
    never_dead = grid["t_dead"].isna()
    disappeared_by_t = (t > grid["t_last"]) & never_dead
    alive_at_t = seen_at_t & (~dead_at_t)
    
    counts = (
        grid.assign(alive=alive_at_t, dead=dead_at_t, disappeared=disappeared_by_t)
            .groupby(["sample_name","position_t"])
            .agg(
                nr_alive=("alive","sum"),
                nr_dead=("dead","sum"),
                nr_disappeared=("disappeared","sum"),
            )
            .reset_index()
    )

    df_t0 = (
        df_tracks[df_tracks["position_t"]==0]
        .groupby("sample_name")
        .agg(nr_organoids_t0=("TrackID","nunique"))
        .reset_index()
    )

    df_general = counts.merge(df_t0, on="sample_name", how="left")
    df_general["percentage_dead"]        = df_general["nr_dead"]        / df_general["nr_organoids_t0"]
    df_general["percentage_alive"]       = df_general["nr_alive"]       / df_general["nr_organoids_t0"]
    df_general["percentage_disappeared"] = df_general["nr_disappeared"] / df_general["nr_organoids_t0"]
    
    df_general_outpath = Path(results_outdir, f"combined_general_{org_type}_dynamics_analysis.csv")
    df_general.to_csv(
        df_general_outpath,
        sep=",",
        index=False
    )
        
    general_pdf_outpath = Path(results_outdir, f"combined_general_{org_type}_dynamics_analysis.pdf")
    
    plot_general_organoid_analysis(
            df_general=df_general,
            outpath=general_pdf_outpath,
            figsize=(8.27, 11.69)
            )
    
    ## TODO PLOT A STACKED BARPLOT OVER TIME WHERE THE TOTAL BAR IS 
    # THE DETECTED ORGANOIDS AT THAT TIMEPOINT. 
    # STACK ALIVE AND DEAD. 
    # THIS CAN SHOW IF ORGANOIDS DISSAPPEAR OR ACTUALLY DIE OVER TIME
    
    for sample_name in df_tracks["sample_name"].unique():
        """
        Create a analysis pdf for each sample separately
        """
        analysis_sample_outdir = Path(sample_outdir, sample_name, org_type)
        if not analysis_sample_outdir.exists():
            analysis_sample_outdir.mkdir(parents=True)
            
        df_tracks_sample = df_tracks[df_tracks["sample_name"] == sample_name]
        
        img_outdir = Path(output_dir, "images", sample_name)
        
        df_tracks_sample_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_track_analysis.csv")
        df_tracks_sample.to_csv(
            df_tracks_sample_outpath,
            sep=",",
            index=False
        )
        
        df_experiment_sample = df_tracks_sample.groupby(["sample_name", "position_t"]).agg(
            nr_organoids_t0=("TrackID", lambda x: x.nunique()),
            nr_dead=("TrackID", lambda x: x[df_tracks_sample.loc[x.index, "dead"]].nunique()),
        ).reset_index()
        df_experiment_sample["nr_alive"]=df_experiment_sample["nr_organoids_t0"] - df_experiment_sample["nr_dead"]
        df_experiment_sample["percentage_dead"]=df_experiment_sample["nr_dead"] / df_experiment_sample["nr_organoids_t0"]
        df_experiment_sample["percentage_alive"]= 1.0 - df_experiment_sample["percentage_dead"]
        
        df_experiment_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_general_analysis.csv")
        df_experiment_sample.to_csv(
            df_experiment_outpath,
            sep=",",
            index=False
        )
        
        sample_pdf_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_analysis.pdf")
    
        raw_img_path = Path(img_outdir, f"{sample_name}.zarr")
        organoid_seg_path = Path(img_outdir, f"{sample_name}_{org_type}_tracked.zarr")
        
        # Get dead_mask_path from metadata
        # only runs if has_dead_data=True, so dead_channel must exist
        dead_mask_path = None
        
        if metadata is not None:
            sample_metadata = metadata[metadata['sample_name'] == sample_name]
            if not sample_metadata.empty:
                # dead_mask_path directly from metadata
                if 'dead_mask_path' in sample_metadata.columns:
                    mask_path_val = sample_metadata['dead_mask_path'].values[0]
                    if pd.notna(mask_path_val) and str(mask_path_val).strip():
                        dead_mask_path = Path(mask_path_val)
                        if dead_mask_path.exists():
                            print(f"  Using dead_mask_path from metadata: {dead_mask_path}")
                        else:
                            print(f"  ⚠️  dead_mask_path in metadata does not exist: {dead_mask_path}")
                            dead_mask_path = None
        
        # Fallback: construct path if not found in metadata
        if dead_mask_path is None:
            fallback_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
            if fallback_path.exists():
                dead_mask_path = fallback_path
                print(f"  Using fallback dead_mask_path: {dead_mask_path}")
            else:
                print(f"  ⚠️  Fallback dead_mask_path does not exist: {fallback_path}")
        
        # Skip PDF plotting if dead_mask_path is still None
        if dead_mask_path is None:
            print(f"  ⚠️  Skipping PDF generation for {sample_name} - no dead_mask_path available")
            print(f"      (has_dead_channel={has_dead_channel}, metadata columns: {list(sample_metadata.columns) if metadata is not None and not sample_metadata.empty else 'N/A'})")
            continue

        ### Plot analysis results per sample
        print(f"- Writing analysis pdf to {sample_pdf_outpath}")
        plot_sample_organoid_analysis(
            df_experiment_sample=df_experiment_sample,
            df_sample_tracks=df_tracks_sample,
            outpath=sample_pdf_outpath,
            raw_img_path=raw_img_path,
            organoid_seg_path=organoid_seg_path,
            dead_mask_path=dead_mask_path,
            dead_perc_threshold=dead_perc_threshold,
            figsize=(8.27, 11.69)
            )
    
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return()

def filter_organoid_tracks(
    metadata,
    config=None,
    output_dir=None,
    exp_duration=None,
    min_track_length=None,
    max_track_length=None,
    min_size=None,
    time_type="frames",
    cell_type="organoid",
    df_input_path=None  # Optional: path to input CSV (e.g., advanced features CSV)
    ):
    """
    This code filters tracks based on supplied parameters in the config.yml
    
    Filtering is based on:
    - Maximum experiment length (tcell_exp_duration)    
    - Minimum track length (tcell_min_track_length)
    - Tracks starting at timepoint 1 with a dead dye mean over the dead_dye_threshold (dead_dye_threshold)
    
    Additonally, all tracks are cut down to:
    - Maximum track length (tcell_max_track_length)
    
    Parameters
    ----------
    df_input_path : Path or str, optional
        Path to input CSV file. If provided, reads from this file instead of the default
        combined_track_features.csv. Useful for reading advanced features CSV with active killing data.
    
    Output:
    - A .csv file containing filtered tracks from all experiments
    """
    
    # assert config is not None or all(
    #     [output_dir]
    # ), "Either 'config' or 'output_dir parameters must be supplied"
    
    start_time = time.time()
    
    print(f"--------------- Filtering tracks ---------------")
    
    # if not all([output_dir]):
    #     output_dir = config['output_dir']
    #     exp_duration = config[f'{cell_type}_exp_duration']
    #     min_track_length = config[f'{cell_type}_min_track_length']
    #     max_track_length = config[f'{cell_type}_max_track_length']
    
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    qc_outdir = Path(analysis_outdir, "quality_control")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not qc_outdir.exists():
        qc_outdir.mkdir(parents=True)
    
    # Use provided input path or default to combined_track_features.csv
    if df_input_path is not None:
        df_all_tracks_path = Path(df_input_path)
        print(f"  Using input file: {df_all_tracks_path.name}")
    else:
        df_all_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    
    if not df_all_tracks_path.exists():
        raise FileNotFoundError(
            f"Track features file not found: {df_all_tracks_path}\n"
            f"Please run Feature Extraction for {cell_type} first."
        )
    
    df_all_tracks = pd.read_csv(df_all_tracks_path)
    
    # Dynamically detect *_line_condition columns from metadata
    line_condition_cols = [c for c in metadata.columns if c.endswith('_line_condition')]
    group_cols = ['TrackID', 'sample_name'] + line_condition_cols + ['exp_nr', 'well']
    cols_present = [c for c in group_cols if c in metadata.columns]
    metadata_info = metadata.loc[:, cols_present].copy()
    
    df_all_tracks_filt = pd.merge(df_all_tracks, metadata_info, how="left", on="sample_name")
    
    # Function to count the number of unique tracks in the DataFrame
    def count_tracks(df_all_tracks, col_name="nr_tracks", df_track_counts=None):
        # Dynamically detect *_line_condition columns
        line_cols_in_tracks = [c for c in df_all_tracks.columns if c.endswith('_line_condition')]
        potential_group_cols = ['sample_name'] + line_cols_in_tracks + ['exp_nr', 'well']
        group_cols_for_count = [c for c in potential_group_cols if c in df_all_tracks.columns]
        
        nr_tracks=df_all_tracks.groupby(group_cols_for_count).agg(
            nr_tracks=pd.NamedAgg(column='TrackID', aggfunc='nunique')
        ).reset_index()
        nr_tracks=nr_tracks.rename(columns={"nr_tracks":col_name})
        if df_track_counts is None:
            return(nr_tracks)
        else:
            return(pd.merge(df_track_counts, nr_tracks, how="left"))
    
    # Counting the nr of tracks before filtering
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_before_filtering")
    
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    if time_type=="real_time":
        time_column="time"
    else:
        time_column="position_t"
        if exp_duration is not None:
            exp_duration = exp_duration-1
        if min_track_length is not None:
            min_track_length = min_track_length-1
        if max_track_length is not None:
            max_track_length = max_track_length-1
        
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    if exp_duration is not None:
        df_all_tracks_filt=filter_by_full_duration(
            df=df_all_tracks_filt,
            time_column=time_column,
            exp_duration=exp_duration
            )
        
        df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_exp_duration", df_track_counts=df_track_counts)
    else:
        df_track_counts["nr_tracks_exp_duration"] = df_track_counts["nr_tracks_before_filtering"]

    # Filtering out tracks under specific track length and cutting them down to specified max track length
    df_all_tracks_filt = filter_minimal_track_length(
        df=df_all_tracks_filt,
        min_track_length=min_track_length,
        time_column=time_column,
        group_cols=group_cols
    )
    
    df_all_tracks_filt = trim_to_maximal_track_length(
        df=df_all_tracks_filt,
        max_track_length=max_track_length,
        time_column=time_column,
        group_cols=group_cols
    )
    
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_min_track_length", df_track_counts=df_track_counts)

    if min_size is not None:
        segment_start_frames = df_all_tracks_filt[df_all_tracks_filt['relative_time'] == 1]
        size_filtered_track_ids = segment_start_frames[segment_start_frames['nr_pixels'] >= min_size]['TrackID']
        df_all_tracks_filt = df_all_tracks_filt[df_all_tracks_filt['TrackID'].isin(size_filtered_track_ids)]
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_min_size", df_track_counts=df_track_counts)
    
    plot_filter_count_outpath = Path(qc_outdir, f"BEHAV3D_filter_counts.pdf")
    print(f"- Plotting track counts after filtering steps to {plot_filter_count_outpath}")
    plot_filter_count(
        df_track_counts,
        outpath=plot_filter_count_outpath,
        nr_cols=3,
        rows_per_page = 3,
        filter_cols=["nr_tracks_before_filtering", "nr_tracks_exp_duration", "nr_tracks_min_track_length", "nr_tracks_min_size"]
    )
    
    # Plot touching distributions for all contact types (same as immune/other cells)
    from behav3d.analysis.tcell_analysis import plot_touching_nontouching_distribution
    contact_columns = [col for col in df_all_tracks_filt.columns if col.endswith('_contact') and not col.startswith('active_')]
    for contact_col in contact_columns:
        target_type = contact_col.replace('_contact', '')
        plot_touching_outpath = Path(qc_outdir, f"BEHAV3D_{target_type}_touching_distribution.pdf")
        print(f"- Plotting {target_type} touching distribution to {plot_touching_outpath}")
        plot_touching_nontouching_distribution(
            df_all_tracks_filt, 
            outpath=plot_touching_outpath,
            contact_column=contact_col,
            nr_cols=3,
            rows_per_page=3,
            )
    
    # Plot dead dye distribution if dead channel exists
    if "mean_dead_dye" in df_all_tracks_filt.columns:
        from behav3d.analysis.tcell_analysis import plot_dead_dye_distribution
        
        # Plot distribution at all timepoints
        plot_dead_dye_distr_outpath = Path(qc_outdir, f"BEHAV3D_dead_dye_distribution.pdf")
        print(f"- Plotting dead dye distribution at all timepoints to {plot_dead_dye_distr_outpath}")
        plot_dead_dye_distribution(
            df_all_tracks_filt,
            outpath=plot_dead_dye_distr_outpath,
            nr_cols=2,
            rows_per_page=2
        )
        
        # Plot distribution at timepoint 1
        plot_dead_dye_distr_t0_outpath = Path(qc_outdir, f"BEHAV3D_dead_dye_distribution_t0.pdf")
        print(f"- Plotting dead dye distribution at timepoint 1 to {plot_dead_dye_distr_t0_outpath}")
        plot_dead_dye_distribution(
            df_all_tracks_filt[df_all_tracks_filt["relative_time"]==1],
            outpath=plot_dead_dye_distr_t0_outpath,
            nr_cols=2,
            rows_per_page=2
        )
    
    filt_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    print(f"- Writing filtered tracks to {filt_tracks_out_path}")
    df_all_tracks_filt.to_csv(filt_tracks_out_path, sep=",", index=False)
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_all_tracks_filt)

def plot_sample_organoid_analysis(
    df_experiment_sample,
    df_sample_tracks,
    outpath,
    raw_img_path,
    organoid_seg_path,
    dead_mask_path,
    dead_perc_threshold=None,
    figsize=(8.27, 11.69)
    ):
    
    raw_img = load_image(raw_img_path)
    seg_img = load_image(organoid_seg_path)
    dead_mask = load_image(dead_mask_path)
    
    raw_img = calc_z_projection(raw_img, projection="max", z_axis=-3)
    seg_img = calc_z_projection(seg_img, projection="max", z_axis=-3)
    dead_mask = calc_z_projection(dead_mask, projection="max", z_axis=-3)
    
    # Use segmentation image time dimension as reference (since it was actually processed)
    seg_timepoints = seg_img.shape[0]
    
    # Truncate raw image to match segmentation timepoints
    raw_img = raw_img[:seg_timepoints]
    
    t_middle = raw_img.shape[0] // 2
    
    t_start_img = np.asarray(raw_img[0])
    t_mid_img = np.asarray(raw_img[t_middle])
    t_end_img = np.asarray(raw_img[-1])
    
    t_start_img_seg = np.asarray(seg_img[0])
    t_mid_img_seg = np.asarray(seg_img[t_middle])
    t_end_img_seg = np.asarray(seg_img[-1])

    t_start_img_dead = np.asarray(dead_mask[0])
    t_mid_img_dead = np.asarray(dead_mask[t_middle])
    t_end_img_dead = np.asarray(dead_mask[-1])
    
    colors = ['black'] + [colormaps.get_cmap('tab20')(i) for i in range(1, 20)]*10
    cmap = LinearSegmentedColormap.from_list('custom_tab20', colors, 200)
    
    with PdfPages(outpath) as pdf:
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(4, 3, figure=fig, wspace=0.1, hspace=0.3)
        
        ax = fig.add_subplot(gs[0, 0])
        ax.imshow(t_start_img[1], cmap="grey", vmin=0, vmax=np.percentile(t_start_img[1], 99))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel('Live Channel', fontsize=10)
        ax.set_title(f'T0')
        
        ax = fig.add_subplot(gs[0, 1])
        ax.imshow(t_mid_img[1], cmap="grey", vmin=0, vmax=np.percentile(t_mid_img[1], 99))
        ax.axis('off')
        ax.set_title(f'T{t_middle}')
        
        ax = fig.add_subplot(gs[0, 2])
        ax.imshow(t_end_img[1], cmap="gray", vmin=0, vmax=np.percentile(t_end_img[1], 99))
        ax.axis('off')
        ax.set_title(f'T{raw_img.shape[0]}')

        ax = fig.add_subplot(gs[1, 0])
        ax.imshow(t_start_img_dead, cmap="grey")
        ax.imshow(np.ma.masked_where(t_start_img_seg == 0, t_start_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel('Dead Mask / Segmentation', fontsize=10)
        # ax.set_title(f'T0')
        
        ax = fig.add_subplot(gs[1, 1])
        ax.imshow(t_mid_img_dead, cmap="grey")
        ax.imshow(np.ma.masked_where(t_mid_img_seg == 0, t_mid_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.axis('off')
        # ax.set_title(f'T{t_middle} - Dead channel/Segmentation')
        
        ax = fig.add_subplot(gs[1, 2])
        ax.imshow(t_end_img_dead, cmap="gray")
        ax.imshow(np.ma.masked_where(t_end_img_seg == 0, t_end_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.axis('off')
        # ax.set_title(f'T{raw_img.shape[0]} - Dead channel/Segmentation')
        
        max_alive = df_experiment_sample["nr_alive"].max()
        ymax = max_alive + 1.0
        
        ax = fig.add_subplot(gs[2:, :])
        sns.lineplot(
        data = df_experiment_sample,
        x = 'position_t', 
        y = 'nr_alive', 
        ax = ax
        )
        ax.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax.set_ylim(0, ymax)
        ax.set_xlim(0)
        ax.set_ylabel('')
        ax.set_title(f'Alive Organoids')
        
        # Secondary Y axis (percent scale)
        ax_percent = ax.twinx()
        ax_percent.set_ylim(0, 100 * ymax / max_alive)
        ax_percent.set_ylabel('% Alive')

        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        pdf.savefig(fig)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 2, figure=fig, wspace=0.5, hspace=0.3)
        
        ax = fig.add_subplot(gs[0, 0])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_nr_dead_mask_pixels',
            title = f'# Dead Mask Pixels (smoothed)',
            ax = ax,
        )
        
        ax = fig.add_subplot(gs[0, 1])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_percentage_dead_mask',
            title = f'# Percentage Dead Organoid (smoothed)',
            ax = ax,
        )
        if dead_perc_threshold is not None:
            ax.axhline(dead_perc_threshold, color='black', linestyle=':', linewidth=1)
            # Add small label on the right
            ax.text(
                x=ax.get_xlim()[1],  # far right of x-axis
                y=dead_perc_threshold,
                s=f'  Dead Threshold ({dead_perc_threshold})',
                va='center', ha='left',
                fontsize=5, color='black',
            )
            
        ax = fig.add_subplot(gs[1, 0])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_mean_dead_dye',
            title = f'# Mean Dead Dye (smoothed)',
            ax = ax,
        )
        
        ax = fig.add_subplot(gs[1, 1])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'volume',
            title = f'Volume',
            ax = ax,
        )
        
        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        
        # plt.show()
        pdf.savefig(fig)
        plt.close(fig)
    
def plot_feature_lineplot(
        df_tracks,
        feature,
        title,
        ax,
        t_colname = 'position_t',
        hue = "TrackID"
    ):
    sns.lineplot(
        data = df_tracks,
        x = t_colname, 
        y = feature, 
        hue=hue, 
        ax = ax
        )
    ax.set_ylabel('')
    ax.set_title(title)
    ax.legend(
        bbox_to_anchor=(1.05, 1),  # Position to the right
        loc='upper left',
        borderaxespad=0,
        prop={'size': 5},  # Smaller font size
        title=hue  # Optional: keep or remove the title
    ).set_title(hue, prop={'size': 6})
    return ax

def plot_general_organoid_analysis(
    df_general,
    outpath,
    figsize=(8.27, 11.69)
    ):
    
    with PdfPages(outpath) as pdf:
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(4, 1, figure=fig, wspace=0.1, hspace=0.3)
        
        ### Plot the percentage organoid death of all samples
        ax = fig.add_subplot(gs[:2, :])
        sns.lineplot(
            data = df_general,
            x = 'position_t', 
            y = 'percentage_dead', 
            hue = 'sample_name',
            ax = ax,
            )
        ax.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax.set_ylim(0, 1.1)
        ax.set_xlim(0)
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('')
        ax.set_title(f'Percentage Dead Organoids')
        ax.legend(
            bbox_to_anchor=(1.05, 1),  # Position to the right
            loc='upper left',
            borderaxespad=0,
            prop={'size': 5},  # Smaller font size
            title=''  # Optional: keep or remove the title
            ).set_title('Sample Name', prop={'size': 6})
        
        ### Plot the barplot organoid death of all samples
        ax = fig.add_subplot(gs[2, :])
        # df_end = df_general[df_general["position_t"] == df_general["position_t"].max()]
        # df_end = df_end[["sample_name", "percentage_dead", "percentage_alive"]]
        
        # df_end.set_index("sample_name").plot(
        #     kind='bar', 
        #     stacked=True, 
        #     ax=ax, 
        #     color=["#CC6666", "#6699CC"], 
        #     alpha=1.0,
        #     width=0.25,
        #     zorder=2
        #     )
        # ax.grid(True, linestyle=':', linewidth=1, alpha=0.5, zorder=1)
        # ax.set_ylim(0, 1.1)
        # ax.set_xlabel('')
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', size=6)
        # ax.set_ylabel('')
        # ax.set_title(f'Percentage Alive Organoids at end of experiment')
        # handles, labels = ax.get_legend_handles_labels()
        # new_order = [1, 0]  # Dead first (red, index 0), then Alive (blue, index 1)
        # ax.legend(
        #     [handles[idx] for idx in new_order], [labels[idx] for idx in new_order],
        #     bbox_to_anchor=(1.05, 1),  # Position to the right
        #     loc='upper left',
        #     borderaxespad=0,
        #     prop={'size': 5},  # Smaller font size
        #     title=''  # Optional: keep or remove the title
        #     )
        
        # Take each sample's own endpoint
        idx = df_general.groupby("sample_name")["position_t"].idxmax()
        df_end = (
            df_general.loc[idx, ["sample_name", "percentage_alive", "percentage_dead", "percentage_disappeared"]]
            .copy()
            .sort_values("sample_name")
        )

        # If percentage_dead is missing, compute it from alive
        if "percentage_dead" not in df_end.columns:
            df_end["percentage_dead"] = 1.0 - df_end["percentage_alive"]

        # Keep the plotting columns in the intended order: dead (red) then alive (blue)
        plot_cols = ["percentage_alive", "percentage_dead", "percentage_disappeared"]

        df_end.set_index("sample_name")[plot_cols].plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=["#6699CC", "#CC6666", "#898989"],  # Red for Dead, Blue for Alive, Grey for Disappeared
            alpha=1.0,
            width=0.25,
            zorder=2,
        )

        ax.grid(True, linestyle=":", linewidth=1, alpha=0.5, zorder=1)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", size=6)
        ax.set_ylabel("")
        ax.set_title("Percentage Alive Organoids at end of experiment")

        # Legend: Dead first, then Alive (matches bar order/colors)
        handles, labels = ax.get_legend_handles_labels()
        order = [2, 0, 1]  # 0=Dead (red), 1=Alive (blue), 2=Disappeared (grey)
        ax.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0,
            prop={"size": 5},
            title="",
        )

        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        
        plt.show()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    