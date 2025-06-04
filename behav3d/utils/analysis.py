
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# def smooth_value_over_time(
#     df, 
#     column, 
#     rolling_meanspeed_window, 
#     min_periods=1,
#     groupby=["TrackID", "sample_name"]
#     ):
#     smoothed_column = df.groupby(groupby)[column].apply(
#         lambda x: x.rolling(
#             window=rolling_meanspeed_window, 
#             min_periods=min_periods
#             ).mean()
#         ).reset_index(drop=True)
#     return(smoothed_column)
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
   
def plot_filter_count(
    df_track_counts,
    outpath,
    nr_cols=3,
    rows_per_page = 3,
    figsize=(8.27, 11.69),
    filter_cols=["nr_tracks_before_filtering", "nr_tracks_exp_duration", "nr_tracks_min_track_length", "nr_tracks_dead_t1"]
    ):
    
    sample_names = df_track_counts["sample_name"].unique()
    n_plots = len(sample_names)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0)
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)
    
    with PdfPages(outpath) as pdf:
        plot_idx = 0
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, wspace=0.5, hspace=0.3)
            remaining_axes = [
                fig.add_subplot(gs[i, j])
                for i in range(rows_per_page)
                for j in range(nr_cols)
            ]
            
            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()
                    continue
                
                sample = sample_names[plot_idx]
                df_subset = df_track_counts[df_track_counts["sample_name"] == sample]
                
                sns.barplot(
                    x=filter_cols, 
                    y=df_subset[filter_cols].values.flatten(),
                    ax=ax
                )
                
                ax.set_title(sample, fontsize=10, loc='center')
                ax.set_xlabel("")
                ax.set_ylabel("Count")
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                
                plot_idx += 1
            
            fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
            # plt.show()
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)