from pathlib import Path
import time

import pandas as pd
import numpy as np

from behav3d.utils import format_time
from behav3d.utils.analysis import plot_filter_count, smooth_value_over_time
from behav3d.utils.preprocessing import calc_z_projection
from behav3d.utils.fileio import load_image

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
    df_tracks_path=None,
    # df_tracks_summarized_path=None,
    ):
    print(f"--------------- Performing T-cell behavioral analysis ---------------")
    start_time = time.time()
    assert config is not None or all(
        [output_dir]
    ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
        
    
    analysis_outdir = Path(output_dir, "analysis", "organoid")
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
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_organoid_combined_track_features_filtered.csv")
    # if df_tracks_summarized_path is None:
    #     df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_organoid_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path)
    # df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    df_tracks["TrackID"] = df_tracks["TrackID"].astype(str)
    
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
    
    df_tracks["dead"]=df_tracks["smoothed_percentage_dead_mask"] > dead_perc_threshold
    df_tracks["dead"] = df_tracks.groupby(["sample_name", "TrackID"])["dead"].transform(lambda x: x.cummax())
    df_tracks = df_tracks[[
            "TrackID",
            "sample_name",
            "organoid_line",
            "tcell_line",
            "exp_nr",
            "well",
            "relative_time",
            "position_t",
            "nr_pixels",
            "nr_dead_mask_pixels",
            "percentage_dead_mask",
            "mean_dead_dye",
            "smoothed_nr_dead_mask_pixels",
            "smoothed_percentage_dead_mask",
            "smoothed_mean_dead_dye",
            "volume",
            "dead"]
            ]
    
    df_t0 = df_tracks[df_tracks["position_t"] == 0].groupby("sample_name").agg(
        nr_organoids_t0=("TrackID", "nunique"),
    ).reset_index()

    df_dead = df_tracks.groupby(["sample_name", "position_t"]).agg(
        nr_organoids=("TrackID", lambda x: x.nunique()),
        nr_dead=("TrackID", lambda x: x[df_tracks.loc[x.index, "dead"]].nunique())
    ).reset_index()
    
    df_general = df_dead.merge(df_t0, on="sample_name", how="left")
    df_general["nr_dead"] = df_general["nr_dead"] + df_general["nr_organoids_t0"] - df_general["nr_organoids"]

    df_general["nr_alive"]=df_general["nr_organoids_t0"] - df_general["nr_dead"]
    df_general["percentage_dead"]=df_general["nr_dead"]  / df_general["nr_organoids_t0"]
    df_general["percentage_alive"]= 1.0 - df_general["percentage_dead"]
    
    df_general_outpath = Path(results_outdir, f"combined_general_organoid_dynamics_analysis.csv")
    df_general.to_csv(
        df_general_outpath,
        sep=",",
        index=False
    )
        
    general_pdf_outpath = Path(results_outdir, f"combined_general_organoid_dynamics_analysis.pdf")
    
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
        analysis_sample_outdir = Path(sample_outdir, sample_name, "organoid")
        if not analysis_sample_outdir.exists():
            analysis_sample_outdir.mkdir(parents=True)
            
        df_tracks_sample = df_tracks[df_tracks["sample_name"] == sample_name]
        
        img_outdir = Path(output_dir, "images", sample_name)
        
        df_tracks_sample_outpath = Path(analysis_sample_outdir, f"{sample_name}_organoid_track_analysis.csv")
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
        
        df_experiment_outpath = Path(analysis_sample_outdir, f"{sample_name}_organoid_general_analysis.csv")
        df_experiment_sample.to_csv(
            df_experiment_outpath,
            sep=",",
            index=False
        )
        
        sample_pdf_outpath = Path(analysis_sample_outdir, f"{sample_name}_organoid_analysis.pdf")
    
        raw_img_path = Path(img_outdir, f"{sample_name}.zarr")
        organoid_seg_path = Path(img_outdir, f"{sample_name}_organoid_tracked.zarr")
        dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")

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
    cell_type="organoid"
    ):
    """
    This code filters tracks based on supplied parameters in the config.yml
    
    Filtering is based on:
    - Maximum experiment length (tcell_exp_duration)    
    - Minimum track length (tcell_min_track_length)
    - Tracks starting at timepoint 1 with a dead dye mean over the dead_dye_threshold (dead_dye_threshold)
    
    Additonally, all tracks are cut down to:
    - Maximum track length (tcell_max_track_length)
    
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
    
    df_all_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_all_tracks = pd.read_csv(df_all_tracks_path)
    
    group_cols = ['TrackID', 'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']
    df_all_tracks_filt = pd.merge(df_all_tracks, metadata, how="left", on="sample_name")

    # Function to count the number of unique tracks in the DataFrame
    def count_tracks(df_all_tracks, col_name="nr_tracks", df_track_counts=None):
        nr_tracks=df_all_tracks.groupby([
            'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']
            ).agg(nr_tracks=pd.NamedAgg(column='TrackID', aggfunc='nunique')).reset_index()
        nr_tracks=nr_tracks.rename(columns={"nr_tracks":col_name})
        if df_track_counts is None:
            return(nr_tracks)
        else:
            return(pd.merge(df_track_counts, nr_tracks, how="left"))
    
    # Counting the nr of tracks before filtering
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_before_filtering")
    
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    if exp_duration is not None:
        df_all_tracks_filt=df_all_tracks_filt.drop(df_all_tracks_filt[df_all_tracks_filt["position_t"]>exp_duration].index)
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_exp_duration", df_track_counts=df_track_counts)

    # Filtering out tracks under specific track length and cutting them down to specified max track length
    if min_track_length is not None:
        df_all_tracks_filt=df_all_tracks_filt.groupby(group_cols).filter(lambda group: len(group) >= min_track_length).reset_index(drop=True)
    if max_track_length is not None:
        df_all_tracks_filt=df_all_tracks_filt.groupby(group_cols).apply(lambda group: group.iloc[:max_track_length]).reset_index(drop=True)
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
            y = 'percentage_alive', 
            hue = 'sample_name',
            ax = ax,
            )
        ax.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax.set_ylim(0, 1.1)
        ax.set_xlim(0)
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('')
        ax.set_title(f'Percentage Alive Organoids')
        ax.legend(
            bbox_to_anchor=(1.05, 1),  # Position to the right
            loc='upper left',
            borderaxespad=0,
            prop={'size': 5},  # Smaller font size
            title=''  # Optional: keep or remove the title
            ).set_title('Sample Name', prop={'size': 6})
        
        ### Plot the barplot organoid death of all samples
        ax = fig.add_subplot(gs[2, :])
        df_end = df_general[df_general["position_t"] == df_general["position_t"].max()]
        df_end = df_end[["sample_name", "percentage_dead", "percentage_alive"]]
        df_end.set_index("sample_name").plot(
            kind='bar', 
            stacked=True, 
            ax=ax, 
            color=["#CC6666", "#6699CC"], 
            alpha=1.0,
            width=0.25,
            zorder=2
            )
        ax.grid(True, linestyle=':', linewidth=1, alpha=0.5, zorder=1)
        ax.set_ylim(0, 1.1)
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', size=6)
        ax.set_ylabel('')
        ax.set_title(f'Percentage Alive Organoids at end of experiment')
        handles, labels = ax.get_legend_handles_labels()
        new_order = [1, 0]  # Dead first (red, index 0), then Alive (blue, index 1)
        ax.legend(
            [handles[idx] for idx in new_order], [labels[idx] for idx in new_order],
            bbox_to_anchor=(1.05, 1),  # Position to the right
            loc='upper left',
            borderaxespad=0,
            prop={'size': 5},  # Smaller font size
            title=''  # Optional: keep or remove the title
            )
        
        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        
        plt.show()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    