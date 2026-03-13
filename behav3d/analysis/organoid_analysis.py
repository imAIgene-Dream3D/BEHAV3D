from pathlib import Path
import time

import json
import pandas as pd
import numpy as np
import scanpy as sc

from behav3d.core.utils import format_time
from behav3d.analysis import smooth_value_over_time


from behav3d.preprocessing import calc_z_projection
from behav3d.io.images import load_image, load_zarr

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colormaps
from sklearn.preprocessing import MinMaxScaler


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
    dead_columns = ["nr_dead_mask_pixels", "percentage_dead_mask"]
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
            groupby=["TrackID", "sample_name"]
        )
    
    df_tracks["smoothed_percentage_dead_mask"] = smooth_value_over_time(
            df_tracks,
            column="percentage_dead_mask",
            groupby=["TrackID", "sample_name"]
        )
    
    if "mean_dead_dye" in df_tracks.columns:
        df_tracks["smoothed_mean_dead_dye"] = smooth_value_over_time(
                df_tracks,
                column="mean_dead_dye",
                groupby=["TrackID", "sample_name"]
            )
    else:
        df_tracks["smoothed_mean_dead_dye"] = np.nan

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


def plot_multi_organoid_death_dynamics(
    output_dir,
    organoid_types,
    figsize=(12, 10)
    ):
    """
    Generate comparison plots for death dynamics across multiple organoid types.
    
    Parameters
    ----------
    output_dir : str or Path
        Base output directory containing analysis results
    organoid_types : list
        List of organoid type names (e.g., ['organoid1', 'organoid2'])
    figsize : tuple
        Figure size for plots
    
    Returns
    -------
    Path or None
        Path to generated PDF, or None if insufficient data
    """
    print(f"--------------- Multi-Organoid Death Dynamics Comparison ---------------")
    
    output_dir = Path(output_dir)
    
    # Check which organoid types have death dynamics data available
    available_data = {}
    for org_type in organoid_types:
        csv_path = Path(output_dir, "analysis", org_type, "results", f"combined_general_{org_type}_dynamics_analysis.csv")
        if csv_path.exists():
            available_data[org_type] = csv_path
            print(f"Found data for {org_type}")
        else:
            print(f"Missing data for {org_type}: {csv_path}")
    
    if len(available_data) < 2:
        print(f"\nNeed at least 2 organoid types with death dynamics data.")
        print(f"   Available: {list(available_data.keys())}")
        print(f"   Run death dynamics analysis for each organoid type first.")
        return None
    
    # Load and combine data from all organoid types
    all_data = []
    for org_type, csv_path in available_data.items():
        df = pd.read_csv(csv_path)
        df["organoid_type"] = org_type
        all_data.append(df)
    
    df_combined = pd.concat(all_data, ignore_index=True)
    
    # Combined label for each sample + organoid_type combination
    df_combined["sample_organoid"] = df_combined["sample_name"] + "_" + df_combined["organoid_type"]
    
    # Output directory
    results_outdir = Path(output_dir, "analysis", "multi_organoid_comparison")
    results_outdir.mkdir(parents=True, exist_ok=True)
    
    # Save combined data
    combined_csv_path = results_outdir / "combined_multi_organoid_death_dynamics.csv"
    df_combined.to_csv(combined_csv_path, index=False)
    print(f"\nCombined data saved to: {combined_csv_path}")
    
    # Plots
    pdf_path = results_outdir / "multi_organoid_death_dynamics_comparison.pdf"
    
    # Unique combinations for color palette
    unique_combinations = df_combined["sample_organoid"].unique()
    n_combinations = len(unique_combinations)
    
    # Enough distinct colors
    if n_combinations <= 10:
        colors = sns.color_palette("tab10", n_combinations)
    elif n_combinations <= 20:
        colors = sns.color_palette("tab20", n_combinations)
    else:
        # Many combinations: continuous colormap
        colors = sns.color_palette("husl", n_combinations)
    
    color_map = dict(zip(unique_combinations, colors))
    
    with PdfPages(pdf_path) as pdf:
        # Plot 1: Line plot - Percentage Dead per sample + organoid_type
        fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
        
        sns.lineplot(
            data=df_combined,
            x='position_t',
            y='percentage_dead',
            hue='sample_organoid',
            ax=ax1,
            linewidth=2,
            palette=color_map
        )
        
        ax1.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax1.set_ylim(0, 1.05)
        ax1.set_xlim(0)
        ax1.set_xlabel('Timepoint')
        ax1.set_ylabel('')
        ax1.set_title('Percentage of Dead Multi-Organoids')
        
        # Format y-axis as percentage
       #ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
        
        ax1.legend(
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            borderaxespad=0,
            prop={'size': 7},
            title='Sample Name'
        )
        
        plt.tight_layout()
        plt.show()
        pdf.savefig(fig1, bbox_inches='tight')
        plt.close(fig1)
        
        # Plot 2: Bar plot - End state per sample + organoid_type
        fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
        
        # Final timepoint for each sample and organoid type
        idx = df_combined.groupby(["sample_organoid"])["position_t"].idxmax()
        df_end = df_combined.loc[idx, ["organoid_type", "sample_name", "sample_organoid", 
                                        "percentage_dead", "percentage_alive", "percentage_disappeared"]].copy()
        
        # Sort by sample_name then organoid_type: consistent ordering
        df_end = df_end.sort_values(["sample_name", "organoid_type"]).reset_index(drop=True)
        
        # x-axis labels: "organoid_type_sample_name"
        df_end["bar_label"] = df_end["organoid_type"] + "_" + df_end["sample_name"]
        
        x = np.arange(len(df_end))
        width = 0.6
        
        # Stacked bar plot: Alive (bottom), Dead (middle), Disappeared (top)
        bars_alive = ax2.bar(x, df_end["percentage_alive"], width, 
                             label="percentage_alive", color="#6699CC", zorder=2)
        bars_dead = ax2.bar(x, df_end["percentage_dead"], width, 
                            bottom=df_end["percentage_alive"], 
                            label="percentage_dead", color="#CC6666", zorder=2)
        bars_disappeared = ax2.bar(x, df_end["percentage_disappeared"], width, 
                                   bottom=df_end["percentage_alive"] + df_end["percentage_dead"], 
                                   label="percentage_disappeared", color="#898989", zorder=2)
        
        ax2.set_xticks(x)
        ax2.set_xticklabels(df_end["bar_label"], rotation=45, ha='right', fontsize=8)
        ax2.grid(True, linestyle=':', linewidth=1, alpha=0.5, zorder=1)
        ax2.set_ylim(0, 1.05)
        ax2.set_xlabel('Organoid Name and Sample Name')
        ax2.set_ylabel('')
        ax2.set_title('Percentage of Multiple Organoids States at End of Experiment')
        
        # Format y-axis as percentage
        #x2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
        
        # Legend with correct order: disappeared, alive, dead (matching single organoid plot)
        handles, labels = ax2.get_legend_handles_labels()
        order = [2, 0, 1]  # disappeared, alive, dead
        ax2.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            bbox_to_anchor=(1.02, 1), 
            loc='upper left', 
            borderaxespad=0, 
            prop={'size': 8}
        )
        
        plt.tight_layout()
        plt.show()
        pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig2)
    
    print(f"PDF saved to: {pdf_path}")
    print(f"### DONE\n")
    
    return pdf_path

def run_organoid_morphology_dead_analysis(
    output_dir,
    organoid_types,
    initial_window=5,
    umap_n_neighbors=15,
    umap_min_dist=0.1,
    umap_spread=1.0,
    distance_metric = "euclidean",
    leiden_resolution=1.0,
    metadata=None,
):
    """
    Initial Morphology Related to Death: PCA + Leiden + UMAP

    Analyses whether the initial morphology of organoids (averaged over a
    selected time window) is predictive of whether they eventually die.

    Parameters
    ----------
    output_dir : str or Path
        Base output directory containing track features.
    organoid_types : list of str
        List of organoid types to analyze (e.g. ["organoid1"]).
    initial_window : int or tuple(int, int), default=5
        Window for baseline morphology. If int, uses 0 <= t <= N.
        If tuple, uses t_start <= t <= t_end. Features are averaged over this window.
    umap_n_neighbors : int, default=15
        Number of local neighbors used for graph construction.
    umap_min_dist : float, default=0.1
        Minimum distance between points in UMAP embedding.
    umap_spread : float, default=1.0
        The effective scale of embedded points.
    distance_metric : str, default="euclidean"
        Distance metric for UMAP/neighbors. Valid options include: euclidean,
        manhattan, chebyshev, cosine, correlation, canberra, braycurtis, etc.
        See umap-learn documentation for full list.
    leiden_resolution : float, default=1.0
        Controls cluster granularity. Higher values lead to more clusters.
    metadata : pd.DataFrame or None, default=None
        Metadata table (needed for organoid snapshot generation).
        If None, snapshots are skipped.

    Notes
    -----
    - Features: 4 key morphology descriptors (volume, sphericity,
      solidity, surface-to-volume ratio).
    - Normalization: Global Z-score across all features and samples.
    - Pipeline: PCA -> Neighborhood Graph (on top PCs) -> Leiden Clustering -> UMAP.
    - Data: All organoids from all types and samples are mixed for a shared morphological space.

    Outputs (in {output_dir}/analysis/morpho-dead/)
    ------------------------------------------------
    - morpho_dead_umap.pdf          : UMAP cluster + feature grid (T-cell style)
    - morpho_dead_heatmap.pdf       : Cluster × feature heatmap + violin plots
    - morpho_dead_dynamics.pdf      : Death dynamics per organoid, colored by cluster
    - morpho_dead_snapshots.pdf     : Organoid image crops per cluster at window start/end
    - morpho_dead_organoid_level.csv: Per-organoid results
    - morpho_dead_cluster_summary.csv: Cluster-level summary
    """
    # Hardcoded parameters for consistency
    random_state = 42
    #distance_metric = "euclidean"
    print("--------------- Initial Morphology Related to Death ---------------")

    output_dir = Path(output_dir)
    results_outdir = Path(output_dir, "analysis", "morpho-dead")
    results_outdir.mkdir(parents=True, exist_ok=True)

    # morphology features
    base_feature_cols = [
        'volume',
        'sphericity',
        'solidity',
        'surfrace_to_volume_ratio',
    ]

    required_cols = ["TrackID", "sample_name", "relative_time", "dead"] + base_feature_cols
    if not organoid_types:
        print("WARNING: No organoid types provided.")
        return None

    # Load data
    dfs = []
    missing = []
    for org_type in organoid_types:
        csv_path = Path(
            output_dir,
            "analysis",
            org_type,
            "track_features",
            f"BEHAV3D_{org_type}_combined_track_features_filtered.csv",
        )
        if not csv_path.exists():
            missing.append(org_type)
            continue
        df = pd.read_csv(csv_path, low_memory=False)
        df["organoid_type"] = org_type
        dfs.append(df)

    if missing:
        print(f"WARNING: Missing track features for: {', '.join(missing)}")
    if not dfs:
        print("WARNING: No track features found. Run track filtering/feature extraction first.")
        return None

    df_all = pd.concat(dfs, ignore_index=True)

    missing_cols = [c for c in required_cols if c not in df_all.columns]
    if missing_cols:
        print(f"WARNING: Missing required columns: {missing_cols}")
        return None

    if "SegmentID" not in df_all.columns:
        df_all["SegmentID"] = 0
    df_all["SegmentID"] = df_all["SegmentID"].fillna(0)

    # dead column boolean
    if df_all["dead"].dtype != bool:
        df_all["dead"] = (
            df_all["dead"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "t", "yes"])
        )

    df_all["relative_time"] = pd.to_numeric(df_all["relative_time"], errors="coerce")
    df_all = df_all[df_all["relative_time"].notna()].copy()

    df_all["TrackID"] = df_all["TrackID"].astype(str)
    df_all["SegmentID"] = df_all["SegmentID"].astype(str)
    df_all["organoid_id"] = (
        df_all["sample_name"].astype(str)
        + "_"
        + df_all["organoid_type"].astype(str)
        + "_"
        + df_all["TrackID"]
        + "_"
        + df_all["SegmentID"]
    )
    df_all = df_all.sort_values(by=["organoid_id", "relative_time"])
    df_all["dead"] = df_all.groupby("organoid_id")["dead"].transform(lambda x: x.cummax())

    # Select initial morphology window
    if isinstance(initial_window, (list, tuple)):
        w_start, w_end = initial_window
        print(f"\n  === Selecting Window Range ({w_start} <= t <= {w_end}) ===")
        df_init = df_all[(df_all["relative_time"] >= w_start) & (df_all["relative_time"] <= w_end)].copy()
        df_after = df_all[df_all["relative_time"] > w_end].copy()
    else:
        w_start = 0
        w_end = int(initial_window)
        print(f"\n  === Selecting Window Range ({w_start} <= t <= {w_end}) ===")
        df_init = df_all[df_all["relative_time"] <= w_end].copy()
        df_after = df_all[df_all["relative_time"] > w_end].copy()
    
    if df_init.empty:
        print("WARNING: No data within the initial morphology window.")
        return None

    # Organoids dead within initial window are excluded
    dead_in_init = df_init.groupby("organoid_id")["dead"].any()

    # Death after initial window
    if df_after.empty:
        print(f"WARNING: No data found after the initial window. "
              "All organoids will be labelled as 'alive'. Consider reducing the initial window.")
    dead_after = df_after.groupby("organoid_id")["dead"].any()

    # Feature extraction (mean morphology over window)
    print("\nCalculating baseline morphology (mean over window)")
    agg_dict = {
        "sample_name": ("sample_name", "first"),
        "organoid_type": ("organoid_type", "first"),
        "TrackID": ("TrackID", "first"),
        "SegmentID": ("SegmentID", "first"),
        "n_timepoints": ("relative_time", "nunique"),
    }
    
    # static features (mean)
    for col in base_feature_cols:
        agg_dict[f"{col}_mean"] = (col, "mean")
    
    df_feats = (
        df_init.groupby("organoid_id")
        .agg(**agg_dict)
        .reset_index()
    )

    feature_cols = [f"{col}_mean" for col in base_feature_cols]

    df_feats["dead_in_initial"] = df_feats["organoid_id"].map(dead_in_init).fillna(False)
    df_feats["dead_after"] = df_feats["organoid_id"].map(dead_after).fillna(False)
    
    total_orgs = len(df_feats)
    dead_init_count = int(df_feats["dead_in_initial"].sum())
    print(f"\n  Organoids in initial window: {total_orgs}")
    print(f"  Excluded (dead in initial window): {dead_init_count}")
    df_feats = df_feats[~df_feats["dead_in_initial"]].copy()

    # Drop NaN/inf feature rows
    df_feats[feature_cols] = df_feats[feature_cols].replace([np.inf, -np.inf], np.nan)
    before = len(df_feats)
    df_feats = df_feats.dropna(subset=feature_cols)
    if df_feats.empty:
        print("WARNING: No organoids left after filtering initial-dead and NaNs.")
        return None
    if before != len(df_feats):
        print(f"  Dropped {before - len(df_feats)} organoids due to NaN/inf features.")
    print(f"  Organoids used for clustering: {len(df_feats)}")

    n_obs = len(df_feats)
    if n_obs < 3:
        print("WARNING: Need at least 3 organoids for PCA/UMAP/Leiden clustering.")
        return None

    # Sample distribution
    print("\nSample distribution")
    sample_dist = df_feats.groupby(["organoid_type", "sample_name"]).size().reset_index(name="n_organoids")
    samples_per_type = df_feats.groupby("organoid_type")["sample_name"].nunique()
    
    for org_type in df_feats["organoid_type"].unique():
        df_type_dist = sample_dist[sample_dist["organoid_type"] == org_type]
        n_samples = samples_per_type[org_type]
        total_orgs = df_type_dist["n_organoids"].sum()
        print(f"  {org_type}: {n_samples} sample(s), {total_orgs} organoids")
        for _, row in df_type_dist.iterrows():
            print(f"    - {row['sample_name']}: {row['n_organoids']} organoids")

    #log transformation of any feature????
    '''print("\nLog transformation")
    df_feats[feature_cols] = np.log1p(df_feats[feature_cols])

    # Verify log transformation
    print("\nLog transformation verification")
    for col in feature_cols:
        mean_log = df_feats[col].mean()
        std_log = df_feats[col].std(ddof=0)
        status = "✓" if abs(mean_log) < 0.01 and abs(std_log - 1.0) < 0.1 else "⚠"
        print(f"    {status} {col}_log: μ={mean_log:10.7f}, σ={std_log:10.7f}")
    '''
    # Z-score standardization
    print("\nZ-score standardization")
    
    def _zscore(x):
        """Z-score normalization with safe handling of zero std."""
        std = x.std(ddof=0)
        if std == 0 or np.isnan(std):
            return x * 0.0
        return (x - x.mean()) / std

    df_scaled = df_feats.copy()
    
    for col in feature_cols:
        df_scaled[col + "_z"] = _zscore(df_scaled[col])

    # Verify standardization
    print("\nStandardization verification")
    for col in feature_cols:
        mean_z = df_scaled[col + "_z"].mean()
        std_z = df_scaled[col + "_z"].std(ddof=0)
        status = "✓" if abs(mean_z) < 0.01 and abs(std_z - 1.0) < 0.1 else "⚠"
        print(f"    {status} {col}_z: μ={mean_z:10.7f}, σ={std_z:10.7f}")

    # PCA > Neighbors > Leiden > UMAP
    z_cols = [c + "_z" for c in feature_cols]
    X = df_scaled[z_cols].to_numpy(dtype=float)

    obs_cols = ["organoid_id", "sample_name", "organoid_type", "TrackID", "SegmentID", "n_timepoints", "dead_after"]
    adata = sc.AnnData(X=X, obs=df_scaled[obs_cols].copy())
    adata.var_names = z_cols
    # Store raw (unscaled) features for reference
    adata.layers["raw"] = df_scaled[feature_cols].to_numpy(dtype=float)

    n_neighbors = int(umap_n_neighbors)
    if n_neighbors >= n_obs:
        n_neighbors = max(1, n_obs - 1)
        print(f"\n  Adjusted neighbors to {n_neighbors} due to small dataset.")

    # Step 1: PCA, compute all components for diagnostics
    print("\nRunning PCA -> Neighbors -> Leiden -> UMAP")
    n_comps = min(X.shape[0] - 1, X.shape[1] - 1, 50)
    sc.tl.pca(adata, n_comps=n_comps)

    # Determine n_pcs: enough PCs to explain >= 90% cumulative variance
    var_ratio = adata.uns["pca"]["variance_ratio"]
    cumvar = np.cumsum(var_ratio)
    n_pcs = int(np.searchsorted(cumvar, 0.90) + 1)
    n_pcs = max(2, min(n_pcs, n_comps))  # at least 2, at most all
    print(f"  PCA: {n_comps} components computed")
    print(f"  PCA Variance Ratio (top 5): {', '.join([f'{r:.1%}' for r in var_ratio[:5]])}")
    print(f"  Cumulative Variance: {cumvar[min(n_pcs-1, len(cumvar)-1)]:.1%} in {n_pcs} PCs (threshold: 90%)")

    # Step 2: Neighborhood graph, uses only top n_pcs components
    sc.pp.neighbors(
        adata,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        metric=str(distance_metric),
        random_state=int(random_state),
    )

    # Step 3: Leiden clustering, on the neighborhood graph
    sc.tl.leiden(
        adata,
        flavor="igraph",
        resolution=float(leiden_resolution),
        key_added="cluster",
        random_state=int(random_state),
    )

    # Step 4: UMAP, on the neighborhood graph (for visualization)
    sc.tl.umap(
        adata,
        min_dist=float(umap_min_dist),
        spread=float(umap_spread),
        random_state=int(random_state),
    )

    umap_coords = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["UMAP1"] = umap_coords[:, 0]
    obs["UMAP2"] = umap_coords[:, 1]
    # Leiden cluster labels
    obs["ClusterID"] = (obs["cluster"].astype(int) + 1).astype(str)

    df_out = df_scaled.merge(obs[["organoid_id", "ClusterID", "UMAP1", "UMAP2"]], on="organoid_id", how="left")
    df_out["ClusterID"] = df_out["ClusterID"].astype("category")
    df_out["status"] = df_out["dead_after"].map({True: "Dead", False: "Alive"})


    # Save results
    per_org_csv = results_outdir / "morpho_dead_organoid_level.csv"
    df_out.to_csv(per_org_csv, index=False)

    # Cluster summary with all features
    summary_agg = {
        "n": ("ClusterID", "size"),
        "pct_dead": ("dead_after", "mean"),
    }
    
    # Add mean and std for base morphology features
    for col in base_feature_cols:
        col_name = f"{col}_mean"  # actual column name in df_out
        summary_agg[f"{col}_mean"] = (col_name, "mean")
        summary_agg[f"{col}_std"] = (col_name, "std")
    
    summary = (
        df_out.groupby("ClusterID", observed=True)
        .agg(**summary_agg)
        .reset_index()
    )
    summary["pct_dead"] = summary["pct_dead"] * 100.0
    summary.to_csv(results_outdir / "morpho_dead_cluster_summary.csv", index=False)

    # Print summary table
    print("\n  Cluster Summary")
    for _, row in summary.iterrows():
        feat_str = ", ".join([f"{feat}={row[f'{feat}_mean']:.3f}±{row[f'{feat}_std']:.3f}" 
                              for feat in base_feature_cols])
        print(f"  Cluster {row['ClusterID']}: n={int(row['n'])}, "
              f"dead={row['pct_dead']:.1f}%, {feat_str}")
    print()

    # Composition tables
    comp_dead = df_out.groupby(["dead_after", "ClusterID"], observed=True).size().unstack(fill_value=0)
    comp_dead_pct = (comp_dead.div(comp_dead.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_dead_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_dead_alive.csv", index=False)

    comp_type = df_out.groupby(["ClusterID", "organoid_type"], observed=True).size().unstack(fill_value=0)
    comp_type_pct = (comp_type.div(comp_type.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_type_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_organoid_type.csv", index=False)

    comp_sample = df_out.groupby(["ClusterID", "sample_name"], observed=True).size().unstack(fill_value=0)
    comp_sample_pct = (comp_sample.div(comp_sample.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_sample_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_sample.csv", index=False)

    # Plots
    clusters = sorted(
        df_out["ClusterID"].dropna().unique().tolist(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    def _get_color_palette(n):
        if n <= 10: return sns.color_palette("Set1", n)
        elif n <= 20: return sns.color_palette("tab20", n)
        return sns.color_palette("husl", n)
    cluster_colors = _get_color_palette(len(clusters))
    cluster_color_map = dict(zip(clusters, cluster_colors))

    # Feature columns (continuous) vs metadata columns (categorical)
    umap_feature_cols = [f"{c}_mean" for c in base_feature_cols if f"{c}_mean" in df_out.columns]
    umap_meta_cols = [c for c in ["organoid_type", "sample_name", "status"] if c in df_out.columns]
    info_cols = umap_feature_cols + umap_meta_cols
    sample_cols_set = set(umap_meta_cols)

    print("\nGenerating Plots")

    # PDF 1: UMAP, cluster plot + feature/metadata subplots
    umap_pdf_path = results_outdir / "morpho_dead_umap.pdf"
    print("  Generating UMAP plots...")
    _rows_per_page = 4
    _nr_cols = 2
    _rows_first = 2
    n_sub_plots = len(info_cols)
    _n_rows_total = (n_sub_plots + _nr_cols - 1) // _nr_cols + _rows_first
    _nr_pages = max(1, (_n_rows_total + _rows_per_page - 1) // _rows_per_page)

    with PdfPages(umap_pdf_path) as pdf:
        plot_idx = 0
        for page in range(_nr_pages):
            fig = plt.figure(figsize=(8.27, 11.69))
            gs = GridSpec(_rows_per_page, _nr_cols, figure=fig, wspace=0.3)

            # Page 0: large cluster UMAP spanning top rows
            if page == 0:
                ax = fig.add_subplot(gs[:_rows_first, :])
                sns.scatterplot(
                    data=df_out, x="UMAP1", y="UMAP2", hue="ClusterID",
                    s=20, alpha=0.5, palette="Set1", ax=ax,
                )
                ax.legend(
                    loc="upper left", prop={"size": 8}, bbox_to_anchor=(1, 1),
                    borderpad=0.3, labelspacing=0.4, columnspacing=0.1, frameon=False,
                )
                ax.set_title("ClusterID", fontsize=10, loc="center")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_xlabel(""); ax.set_ylabel("")

            # Smaller subplots on remaining rows
            start_row = _rows_first if page == 0 else 0
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(start_row, _rows_per_page)
                for c in range(_nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_sub_plots:
                    ax.remove()
                    continue
                colorcol = info_cols[plot_idx]

                if colorcol in sample_cols_set or df_out.dtypes.get(colorcol) == bool:
                    # Categorical: colored points + legend
                    sns.scatterplot(
                        data=df_out, x="UMAP1", y="UMAP2", hue=colorcol,
                        s=20, alpha=0.5, palette="Set2", ax=ax,
                    )
                    ax.legend(
                        loc="upper left", prop={"size": 6}, bbox_to_anchor=(1, 1),
                        borderpad=0.3, labelspacing=0.3, columnspacing=0.1, frameon=False,
                    )
                else:
                    # Continuous: scatter + colorbar (clipped to 1st–99th percentile)
                    vals = pd.to_numeric(df_out[colorcol], errors="coerce")
                    vmin = np.nanpercentile(vals, 1)
                    vmax = np.nanpercentile(vals, 99)
                    sc_plot = ax.scatter(
                        df_out["UMAP1"], df_out["UMAP2"],
                        c=vals, cmap="viridis", s=20, alpha=0.5,
                        vmin=vmin, vmax=vmax,
                    )
                    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
                    cbar.ax.tick_params(labelsize=6)

                ax.set_title(colorcol, fontsize=10, loc="center")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_xlabel(""); ax.set_ylabel("")
                plot_idx += 1

            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    print(f"  UMAP PDF saved to: {umap_pdf_path}")
    print("  ClusterID")
    print(
        f"    PCs={n_pcs}, Leiden={leiden_resolution}, "
        f"UMAP(neighbors={umap_n_neighbors}, min_dist={umap_min_dist}, spread={umap_spread})"
    )

    # PDF 2: Heatmap + violin plots
    heatmap_pdf_path = results_outdir / "morpho_dead_heatmap.pdf"
    print("  Generating heatmap + violins...")

    # Cluster means (numeric features only)
    value_cols = [c for c in info_cols if c not in sample_cols_set]
    df_for_means = df_out[value_cols + ["ClusterID"]].copy()
    cluster_means = (
        df_for_means.groupby("ClusterID", observed=True)
        .mean(numeric_only=True).reset_index()
    )
    # Min-max scaling for heatmap
    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]
    X_hm = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X_hm = X_hm.replace([np.inf, -np.inf], np.nan)
    all_nan_cols = X_hm.columns[X_hm.isna().all()].tolist()
    if all_nan_cols:
        X_hm = X_hm.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]
    if len(scale_columns) > 0:
        X_filled = X_hm.fillna(X_hm.median(numeric_only=True))
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])
        df_hm_long = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_hm_long.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    # Violin data
    df_violin = df_out[["ClusterID"] + value_cols].copy()
    for c in value_cols:
        df_violin[c] = pd.to_numeric(df_violin[c], errors="coerce")
    df_violin.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_violin.melt(id_vars="ClusterID", var_name="var", value_name="value")
    cluster_order = sorted(df_violin["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_violin = len(feat_names)
    _vrows = 7
    _vcols = 2
    _v_rows_needed = (n_violin + _vcols - 1) // _vcols
    _v_pages = max(1, (_v_rows_needed + _vrows - 1) // _vrows)

    with PdfPages(heatmap_pdf_path) as pdf:
        # Page 1: full-page heatmap
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        if not overall_heatmap_data.empty:
            try:
                hm = sns.heatmap(
                    overall_heatmap_data, ax=ax, cmap="viridis",
                    cbar=True, yticklabels=True,
                )
                ax.set_title("Min–Max Scaled Cluster Means", fontsize=14, pad=14)
                ax.set_xlabel("ClusterID", fontsize=10)
                ax.set_ylabel("", fontsize=10)
                ax.tick_params(axis="y", labelsize=6)
                ax.tick_params(axis="x", labelsize=8, rotation=0)
                cbar_hm = hm.collections[0].colorbar
                cbar_hm.ax.tick_params(labelsize=8)
                fig.tight_layout(pad=2.0)
            except Exception:
                ax.text(0.5, 0.5, "Heatmap unavailable", ha="center", va="center")
                ax.axis("off")
        else:
            ax.text(0.5, 0.5, "No features for heatmap", ha="center", va="center")
            ax.axis("off")
        pdf.savefig(fig, dpi=600)
        plt.show()
        plt.close(fig)

        # Other pages: violin plots
        v_idx = 0
        for vpage in range(_v_pages):
            fig = plt.figure(figsize=(8.27, 11.69))
            gs = GridSpec(_vrows, _vcols, figure=fig, hspace=1.5, wspace=0.3)
            vaxes = [
                fig.add_subplot(gs[r, c])
                for r in range(_vrows) for c in range(_vcols)
            ]
            for ax in vaxes:
                if v_idx >= n_violin:
                    ax.remove(); continue
                feat = feat_names[v_idx]
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].dropna()
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no data)", ha="center", va="center")
                    ax.axis("off"); v_idx += 1; continue
                try:
                    sns.violinplot(
                        data=sub, x="ClusterID", y="value",
                        order=cluster_order, inner=None, ax=ax, cut=0,
                        color="steelblue",
                    )
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(unavailable)", ha="center", va="center")
                    ax.axis("off"); v_idx += 1; continue
                # orange = mean
                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)), means.values,
                    s=60, color="orange", edgecolor="black", linewidths=0.8, zorder=3,
                )
                ax.set_title(feat, fontsize=9)
                ax.set_xlabel("ClusterID", fontsize=8)
                ax.set_ylabel("Value", fontsize=8)
                ax.tick_params(axis="x", rotation=0, labelsize=7)
                ax.tick_params(axis="y", labelsize=7)
                v_idx += 1

            fig.subplots_adjust(left=0.20, right=0.98, top=0.95, bottom=0.08)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    print(f"  Heatmap PDF saved to: {heatmap_pdf_path}")

    # PDF 3: Death ---------------
    dynamics_pdf_path = results_outdir / "morpho_dead_dynamics.pdf"

    has_dead_feature = ("percentage_dead_mask" in df_all.columns
                        or "smoothed_percentage_dead_mask" in df_all.columns)
    has_dead_col = "dead" in df_all.columns

    if has_dead_feature or has_dead_col:
        print("  Generating death dynamics plots...")
        clustered_ids = set(df_out["organoid_id"].unique())
        df_dynamics = df_all[df_all["organoid_id"].isin(clustered_ids)].copy()
        id_to_cluster = df_out.set_index("organoid_id")["ClusterID"].to_dict()
        df_dynamics["ClusterID"] = df_dynamics["organoid_id"].map(id_to_cluster)

        # Individual curves: local smoothing with min_periods=1 so curves start at t=0
        if has_dead_feature and "percentage_dead_mask" in df_dynamics.columns:
            df_dynamics["viable_fraction"] = (
                df_dynamics
                .groupby("organoid_id")["percentage_dead_mask"]
                .transform(lambda s: s.rolling(window=20, min_periods=1, center=False).mean())
            )
            df_dynamics["viable_fraction"] = (1.0 - df_dynamics["viable_fraction"].clip(0, 1))

        # Cluster mean: binary alive from dead column
        if has_dead_col:
            df_dynamics["alive"] = (~df_dynamics["dead"]).astype(int)

        with PdfPages(dynamics_pdf_path) as pdf:
            # --- Individual organoid curves (viable fraction per cluster) ---
            if "viable_fraction" in df_dynamics.columns:
                for cl in clusters:
                    cl_data = df_dynamics[df_dynamics["ClusterID"] == cl]
                    n_orgs = cl_data["organoid_id"].nunique()
                    fig, ax = plt.subplots(figsize=(12, 6))
                    for oid in cl_data["organoid_id"].unique():
                        org_data = cl_data[cl_data["organoid_id"] == oid].dropna(subset=["viable_fraction"])
                        org_data = org_data.sort_values("relative_time")
                        ax.plot(
                            org_data["relative_time"], org_data["viable_fraction"],
                            color=cluster_color_map[cl], alpha=0.5, linewidth=0.8,
                        )
                    ax.axvspan(w_start, w_end, alpha=0.1, color="gray", label="Initial Window")
                    ax.set_xlabel("Relative Time", fontsize=11)
                    ax.set_ylabel("Viable Fraction (1 \u2212 dead mask coverage)", fontsize=11)
                    ax.set_title(f"Death Dynamics \u2014 Cluster {cl} ({n_orgs} organoids)", fontsize=13)
                    ax.set_ylim(-0.05, 1.05)
                    ax.set_xlim(left=0)
                    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
                    ax.legend(loc="upper right", prop={"size": 8})
                    plt.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

            # --- Cluster mean survival curve (binary alive, fixed cohort) ---
            if has_dead_col:
                fig, ax = plt.subplots(figsize=(12, 7))
                for cl in clusters:
                    cl_members = sorted(df_out.loc[df_out["ClusterID"] == cl, "organoid_id"].unique().tolist())
                    if not cl_members:
                        continue
                    cl_data = df_dynamics[df_dynamics["organoid_id"].isin(cl_members)].copy()
                    cl_data = cl_data[["organoid_id", "relative_time", "dead"]].dropna(subset=["relative_time"])
                    if cl_data.empty:
                        continue
                    # First timepoint where each organoid is dead (if never dead -> +inf)
                    first_dead_t = (
                        cl_data[cl_data["dead"] == True]
                        .groupby("organoid_id")["relative_time"]
                        .min()
                    )
                    first_dead_t = first_dead_t.reindex(cl_members).fillna(np.inf)
                    rel_times = np.sort(cl_data["relative_time"].unique())
                    frac_alive = [((first_dead_t > t).sum() / len(cl_members)) for t in rel_times]
                    ax.plot(
                        rel_times, frac_alive,
                        color=cluster_color_map[cl], linewidth=2.5, label=f"Cluster {cl}",
                    )

                ax.axvspan(w_start, w_end, alpha=0.1, color="gray", label="Initial Window")
                ax.set_xlabel("Relative Time", fontsize=11)
                ax.set_ylabel("% Alive", fontsize=11)
                ax.set_title("Death Dynamics per Morphology Cluster", fontsize=13)
                ax.set_ylim(-0.05, 1.05)
                ax.set_xlim(left=0)
                ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
                ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 8})
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"  Death dynamics PDF saved to: {dynamics_pdf_path}")
    else:
        print("  WARNING: No death feature found — skipping dynamics plot.")

    # ===== PDF 4: Organoid Snapshots — separate channels per organoid =====
    snapshot_pdf_path = results_outdir / "morpho_dead_snapshots.pdf"
    n_examples_per_cluster = 3
    t_col = "position_t" if "position_t" in df_init.columns else "relative_time"

    org_t_bounds = (
        df_init.groupby("organoid_id")[t_col]
        .agg(["min", "max"])
        .rename(columns={"min": "t_img_start", "max": "t_img_end"})
    )

    # Channel display colormaps (one per channel)
    from matplotlib.colors import LinearSegmentedColormap as _LSC
    _ch_cmaps = [
        _LSC.from_list("_g", ["black", "#00FF00"]),   # Ch 0: green
        _LSC.from_list("_m", ["black", "#FF00FF"]),   # Ch 1: magenta
        _LSC.from_list("_c", ["black", "#00FFFF"]),   # Ch 2: cyan
        _LSC.from_list("_y", ["black", "#FFFF00"]),   # Ch 3: yellow
        _LSC.from_list("_o", ["black", "#FF8000"]),   # Ch 4: orange
        _LSC.from_list("_b", ["black", "#8080FF"]),   # Ch 5: light blue
    ]
    _ch_color_names = ["green", "magenta", "cyan", "yellow", "orange", "light blue"]

    _zarr_cache = {}

    def _load_zarr_cached(zpath, cache_key):
        if cache_key in _zarr_cache:
            return _zarr_cache[cache_key]
        if not Path(zpath).exists():
            _zarr_cache[cache_key] = None
            return None
        try:
            z = load_zarr(zpath)
            _zarr_cache[cache_key] = z
            return z
        except Exception as e:
            print(f"    WARNING: Could not load {zpath}: {e}")
            _zarr_cache[cache_key] = None
            return None

    def _crop_per_channel(raw_zarr_data, seg_zarr_data, t_idx, track_id, margin=30):
        """Return (crops_list, contour_mask) for one organoid at one timepoint."""
        from scipy.ndimage import binary_dilation, binary_erosion
        t_idx = int(t_idx)
        if seg_zarr_data is None or raw_zarr_data is None:
            return None
        if t_idx < 0 or t_idx >= seg_zarr_data.shape[0]:
            return None

        seg_vol = np.asarray(seg_zarr_data[t_idx])
        seg_2d = seg_vol.max(axis=0) if seg_vol.ndim == 3 else seg_vol
        mask = (seg_2d == track_id)
        if not mask.any():
            return None
        ys, xs = np.where(mask)
        y0 = max(0, ys.min() - margin)
        y1 = min(seg_2d.shape[0], ys.max() + margin + 1)
        x0 = max(0, xs.min() - margin)
        x1 = min(seg_2d.shape[1], xs.max() + margin + 1)

        # Contour: dilated − eroded boundary of the mask
        mask_crop = mask[y0:y1, x0:x1]
        contour = binary_dilation(mask_crop, iterations=1) ^ binary_erosion(mask_crop, iterations=1)

        has_ch = (raw_zarr_data.ndim == 5)
        n_ch = raw_zarr_data.shape[1] if has_ch else 1
        crops = []
        for ci in range(n_ch):
            if has_ch:
                ch_vol = np.asarray(raw_zarr_data[t_idx, ci])
            else:
                ch_vol = np.asarray(raw_zarr_data[t_idx])
            ch_2d = ch_vol.max(axis=0) if ch_vol.ndim == 3 else ch_vol
            crops.append(ch_2d[y0:y1, x0:x1])
        return crops, contour

    # Build organoid_type → channel index map from metadata
    _type_to_channel = {}
    _channel_cols = []
    if metadata is not None:
        _channel_cols = sorted([c for c in metadata.columns if c.startswith("channel") and c[7:].isdigit()])
        if _channel_cols:
            for _, mrow in metadata.iterrows():
                for col in _channel_cols:
                    ch_idx = int(col.replace("channel", ""))
                    label = str(mrow[col]).strip() if pd.notna(mrow[col]) else ""
                    if label and label not in _type_to_channel:
                        _type_to_channel[label] = ch_idx

    def _get_organoid_channel(sample_name, organoid_type):
        """Resolve organoid channel per sample, fallback to global mapping."""
        if metadata is None or not _channel_cols:
            return _type_to_channel.get(organoid_type, -1)
        row = metadata[metadata["sample_name"] == sample_name]
        if row.empty:
            return _type_to_channel.get(organoid_type, -1)
        row0 = row.iloc[0]
        for col in _channel_cols:
            label = str(row0[col]).strip() if pd.notna(row0[col]) else ""
            if label == organoid_type:
                return int(col.replace("channel", ""))
        return _type_to_channel.get(organoid_type, -1)

    # Collect snapshot data per cluster — diverse sampling across types & samples
    snapshot_data = {}
    n_channels_found = 1
    try:
        print("  Generating organoid snapshots (separate channels)...")
        for cl in clusters:
            cl_orgs = df_out[df_out["ClusterID"] == cl]
            # Diverse sampling: one per (organoid_type, sample_name), then fill randomly
            selected = (
                cl_orgs
                .groupby(["organoid_type", "sample_name"], group_keys=False)
                .sample(n=1, random_state=random_state)
                .reset_index(drop=True)
            )
            selected = selected.drop_duplicates(subset=["organoid_id"]).copy()
            remaining = cl_orgs[~cl_orgs["organoid_id"].isin(selected["organoid_id"])]
            n_extra = max(0, n_examples_per_cluster - len(selected))
            if n_extra > 0 and len(remaining) > 0:
                extra = remaining.sample(min(n_extra, len(remaining)), random_state=random_state)
                selected = pd.concat([selected, extra], ignore_index=True)
            selected = selected.drop_duplicates(subset=["organoid_id"]).copy()
            selected = selected.head(n_examples_per_cluster)

            snapshot_data[cl] = []
            for _, org in selected.iterrows():
                oid = org["organoid_id"]
                if oid not in org_t_bounds.index:
                    continue
                ts = int(org_t_bounds.loc[oid, "t_img_start"])
                te = int(org_t_bounds.loc[oid, "t_img_end"])
                try:
                    track_id = int(float(org["TrackID"]))
                except (ValueError, TypeError):
                    continue
                sname = org["sample_name"]
                otype = org["organoid_type"]
                seg_path = output_dir / "images" / sname / f"{sname}_{otype}_tracked.zarr"
                raw_path = output_dir / "images" / sname / f"{sname}.zarr"
                seg_z = _load_zarr_cached(seg_path, f"seg_{sname}_{otype}")
                raw_z = _load_zarr_cached(raw_path, f"raw_{sname}")
                if seg_z is None or raw_z is None:
                    continue
                result_s = _crop_per_channel(raw_z, seg_z, ts, track_id)
                result_e = _crop_per_channel(raw_z, seg_z, te, track_id)
                crops_s, contour_s = result_s if result_s else (None, None)
                crops_e, contour_e = result_e if result_e else (None, None)
                if crops_s is not None:
                    n_channels_found = len(crops_s)
                elif crops_e is not None:
                    n_channels_found = len(crops_e)
                org_ch_idx = _get_organoid_channel(sname, otype)
                if crops_s is not None or crops_e is not None:
                    snapshot_data[cl].append({
                        "crops_s": crops_s, "crops_e": crops_e,
                        "contour_s": contour_s, "contour_e": contour_e,
                        "oid": oid, "otype": otype, "sname": sname,
                        "ts": ts, "te": te, "org_ch": org_ch_idx,
                    })
    except Exception as e:
        print(f"  WARNING: Snapshot generation failed: {e}")
        import traceback; traceback.print_exc()
        snapshot_data = {}
    finally:
        _zarr_cache.clear()

    clusters_with_data = [cl for cl in clusters if snapshot_data.get(cl)]
    if clusters_with_data:
        n_ch = n_channels_found

        def _show_crop(ax, img, cmap, contour=None, highlight=False):
            """Display a single channel crop. Red border if highlight=True."""
            if img is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="white")
                ax.set_facecolor("black")
                return
            img_f = img.astype(float)
            pos = img_f[img_f > 0]
            vmin = float(np.percentile(pos, 1)) if pos.size > 0 else 0
            vmax = float(np.percentile(pos, 99)) if pos.size > 0 else max(1, img_f.max())
            ax.imshow(img_f, cmap=cmap, vmin=vmin, vmax=vmax)
            if highlight and contour is not None:
                from matplotlib.colors import ListedColormap
                red_cmap = ListedColormap(["red"])
                contour_overlay = np.ma.masked_where(~contour, contour.astype(float))
                ax.imshow(contour_overlay, cmap=red_cmap, vmin=0, vmax=1, alpha=0.9)

        # Build channel labels from metadata for column headers
        _ch_labels = {}
        if metadata is not None:
            ch_cols = sorted([c for c in metadata.columns if c.startswith("channel") and c[7:].isdigit()])
            if ch_cols:
                first_row = metadata.iloc[0]
                for col in ch_cols:
                    ci = int(col.replace("channel", ""))
                    label = str(first_row[col]).strip() if pd.notna(first_row[col]) else ""
                    if label:
                        _ch_labels[ci] = label

        with PdfPages(snapshot_pdf_path) as pdf:
            for cl in clusters_with_data:
                examples = snapshot_data[cl]
                n_ex = len(examples)
                n_cols = n_ch * 2
                fig, axes = plt.subplots(
                    n_ex, n_cols,
                    figsize=(2.2 * n_cols + 1.5, 2.5 * n_ex + 1.5),
                    squeeze=False,
                )
                fig.suptitle(
                    f"Cluster {cl} \u2014 Organoid Snapshots",
                    fontsize=13, y=0.99, fontweight="bold",
                    color=cluster_color_map.get(cl, "black"),
                )
                fig.subplots_adjust(left=0.14, top=0.88, bottom=0.04, hspace=0.3, wspace=0.08)

                # Column headers (black text) — channel index + label if available
                for ci in range(n_ch):
                    ch_label = _ch_labels.get(ci, "")
                    header = f"Ch {ci}" + (f" ({ch_label})" if ch_label else "")
                    axes[0, ci].set_title(header, fontsize=8, color="black")
                    axes[0, n_ch + ci].set_title(header, fontsize=8, color="black")

                # Timepoint group headers
                fig.text(0.14 + (0.86 - 0.14) * 0.25, 0.92,
                         f"Window Start (t={w_start})", ha="center", fontsize=10, style="italic")
                fig.text(0.14 + (0.86 - 0.14) * 0.75, 0.92,
                         f"Window End (t={w_end})", ha="center", fontsize=10, style="italic")

                for ex_i, snap in enumerate(examples):
                    crops_s = snap["crops_s"]
                    crops_e = snap["crops_e"]
                    contour_s = snap["contour_s"]
                    contour_e = snap["contour_e"]
                    otype = snap["otype"]
                    sname = snap["sname"]
                    org_ch = snap["org_ch"]

                    for ci in range(n_ch):
                        cmap = _ch_cmaps[ci % len(_ch_cmaps)]
                        is_org_ch = (ci == org_ch)
                        img_s = crops_s[ci] if (crops_s and ci < len(crops_s)) else None
                        img_e = crops_e[ci] if (crops_e and ci < len(crops_e)) else None
                        _show_crop(axes[ex_i, ci], img_s, cmap,
                                   contour=contour_s, highlight=is_org_ch)
                        _show_crop(axes[ex_i, n_ch + ci], img_e, cmap,
                                   contour=contour_e, highlight=is_org_ch)
                        axes[ex_i, ci].axis("off")
                        axes[ex_i, n_ch + ci].axis("off")

                        # Red border around the organoid's channel
                        if is_org_ch:
                            from matplotlib.patches import Rectangle
                            for ax, img in [
                                (axes[ex_i, ci], img_s),
                                (axes[ex_i, n_ch + ci], img_e),
                            ]:
                                if img is not None:
                                    h, w = img.shape[:2]
                                    ax.add_patch(
                                        Rectangle(
                                            (-1.0, -1.0), w + 2.0, h + 2.0,
                                            transform=ax.transData,
                                            fill=False, edgecolor="red",
                                            linewidth=2.5, zorder=30, clip_on=False
                                        )
                                    )

                    # Row label
                    row_y = axes[ex_i, 0].get_position().y0 + axes[ex_i, 0].get_position().height / 2
                    fig.text(
                        0.01, row_y, f"{otype}\n{sname}",
                        fontsize=7, va="center", ha="left", style="italic",
                    )

                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"  Snapshots PDF saved to: {snapshot_pdf_path}")
    else:
        print("  WARNING: Could not generate snapshots — no images found.")

    # ===== PDF 5: Composition bar plots =====
    comp_pdf_path = results_outdir / "morpho_dead_composition.pdf"
    print("  Generating composition bar plots...")

    with PdfPages(comp_pdf_path) as pdf:
        # Bar 1: Dead vs Alive composition by cluster
        comp_dead_plot = comp_dead.div(comp_dead.sum(axis=1), axis=0) * 100.0
        comp_dead_plot = comp_dead_plot.reindex([True, False]).fillna(0.0)
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x_labels = ["Dead", "Alive"]
        x = np.arange(len(x_labels))
        bottom = np.zeros(len(x_labels))
        for cl in clusters:
            vals = comp_dead_plot.get(cl, pd.Series([0.0, 0.0])).values
            ax.bar(x, vals, bottom=bottom, label=f"Cluster {cl}", color=cluster_color_map.get(cl))
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels(x_labels)
        ax.set_ylabel("Percent of group"); ax.set_ylim(0, 100)
        ax.set_title("Cluster Composition: Dead vs Alive")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 7})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Bar 2: Organoid type composition per cluster
        comp_type_plot = comp_type.div(comp_type.sum(axis=1), axis=0) * 100.0
        comp_type_plot = comp_type_plot.reindex(clusters).fillna(0.0)
        types = comp_type_plot.columns.tolist()
        type_colors = _get_color_palette(len(types))
        type_color_map = dict(zip(types, type_colors))
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x = np.arange(len(comp_type_plot.index))
        bottom = np.zeros(len(comp_type_plot.index))
        for t in types:
            vals = comp_type_plot[t].values
            ax.bar(x, vals, bottom=bottom, label=t, color=type_color_map[t])
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels([str(c) for c in comp_type_plot.index])
        ax.set_ylabel("Percent in cluster"); ax.set_ylim(0, 100)
        ax.set_title("Organoid Type Composition per Cluster")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 7})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Bar 3: Sample composition per cluster
        comp_sample_plot = comp_sample.div(comp_sample.sum(axis=1), axis=0) * 100.0
        comp_sample_plot = comp_sample_plot.reindex(clusters).fillna(0.0)
        samples = comp_sample_plot.columns.tolist()
        sample_colors = _get_color_palette(len(samples))
        sample_color_map = dict(zip(samples, sample_colors))
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        x = np.arange(len(comp_sample_plot.index))
        bottom = np.zeros(len(comp_sample_plot.index))
        for s in samples:
            vals = comp_sample_plot[s].values
            ax.bar(x, vals, bottom=bottom, label=s, color=sample_color_map[s])
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels([str(c) for c in comp_sample_plot.index])
        ax.set_ylabel("Percent in cluster"); ax.set_ylim(0, 100)
        ax.set_title("Sample Composition per Cluster")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 6})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Box 4: Time of death per cluster (boxplot + individual organoid points)
        dead_orgs = df_out[df_out["dead_after"] == True].copy()
        if not dead_orgs.empty:
            dead_ids = set(dead_orgs["organoid_id"].unique())
            df_death_times = (
                df_all[df_all["organoid_id"].isin(dead_ids) & (df_all["dead"] == True)]
                .groupby("organoid_id")["relative_time"]
                .min()
                .reset_index()
                .rename(columns={"relative_time": "death_time"})
            )
            df_death_times = df_death_times.merge(
                df_out[["organoid_id", "ClusterID"]].drop_duplicates(),
                on="organoid_id", how="left",
            )
            df_death_times = df_death_times.dropna(subset=["ClusterID", "death_time"])

            if not df_death_times.empty:
                fig, ax = plt.subplots(1, 1, figsize=(10, 7))
                sns.boxplot(
                    data=df_death_times, x="ClusterID", y="death_time",
                    order=clusters, color="white", fliersize=0,
                    boxprops=dict(edgecolor="black"),
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(color="black"),
                    capprops=dict(color="black"),
                    ax=ax,
                )
                sns.stripplot(
                    data=df_death_times, x="ClusterID", y="death_time",
                    order=clusters, ax=ax, size=5, alpha=0.6, jitter=0.2,
                    palette=cluster_color_map, hue="ClusterID", hue_order=clusters,
                    legend=False,
                )
                ax.set_xlabel("Morphology Cluster", fontsize=11)
                ax.set_ylabel("Death Timepoint (relative time)", fontsize=11)
                ax.set_title("Organoid Death Timepoint per Morphology Cluster", fontsize=13)
                ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.show()
                plt.close(fig)

    print(f"  Composition PDF saved to: {comp_pdf_path}")

    print(f"\n  Results saved to: {results_outdir}")
    print("\n### ANALYSIS COMPLETE ###\n")
    return results_outdir
