
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def filter_by_full_duration(df: pd.DataFrame,
                            exp_duration=None,
                            time_column="position_t"
                            ) -> pd.DataFrame:
    """
    Cut rows where time_column > exp_duration.
    """
    if exp_duration is not None:
        df=df.drop(df[df[time_column]>exp_duration].index)
    return df

def filter_minimal_track_length(
    df,
    min_track_length=None,
    time_column="position_t",
    group_cols=["TrackID", "sample_name"]
    ):

    first_t = (
            df.groupby(group_cols)[time_column]
            .transform("first")
        )

    # --- 2) Min track length on the chosen axis (keep groups with sufficient span) ---
    if min_track_length is not None:
        # span = last - first on the chosen axis
        span = (
            df.groupby(group_cols)[time_column]
            .transform("last") - first_t
        )
        df = df[span >= float(min_track_length)].reset_index(drop=True)
    return df
        
def trim_to_maximal_track_length(
    df,
    max_track_length=None,
    time_column="position_t",
    group_cols=["TrackID", "sample_name"]
    ):

    if max_track_length is not None:
        last_t = df.groupby(group_cols)[time_column].transform("last")
        within_last_window = (last_t - df[time_column]) <= float(max_track_length)
        df = df[within_last_window].reset_index(drop=True)
    return df

def plot_filter_count(
    df_track_counts,
    outpath,
    nr_cols=3,
    rows_per_page = 3,
    figsize=(8.27, 11.69),
    filter_cols=["nr_tracks_before_filtering", "nr_tracks_exp_duration", "nr_tracks_min_track_length", "nr_tracks_dead_t1"],
    plot_results=False
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
                
                # Force y-axis to use integer formatting (track counts should always be whole numbers)
                ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
                
                # Use plt.setp for better compatibility with matplotlib
                plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
                
                plot_idx += 1
            
            fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
            # plt.show()
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)