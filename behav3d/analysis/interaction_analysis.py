"""
Interaction Analysis for BEHAV3D

Analyzes interactions from the organoid's point of view:
- Cumulative interactions over time with selected cell types
- Comparison between organoids that survive vs die
- Per-sample breakdowns

Outputs:
- PDF plots to analysis/{organoid_type}/interaction_analysis/
- CSV statistics per sample
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_interaction_analysis(
    output_dir: str,
    cell_type: str,
    interacting_cell_types: list,
    dead_threshold: float = 0.02,
    df_tracks_path: str = None,
    show_plots: bool = True,
):
    """
    Run interaction analysis for an organoid cell type.
    
    Parameters
    ----------
    output_dir : str
        Base output directory for BEHAV3D results
    cell_type : str
        Name of the organoid cell type (e.g., "organoid")
    interacting_cell_types : list
        List of cell types to analyze interactions with (e.g., ["tcell", "macrophage"])
    dead_threshold : float
        Threshold for percentage_dead_mask to classify organoid as dead (default: 0.02)
    df_tracks_path : str, optional
        Path to filtered track features CSV. If None, uses default location.
    show_plots : bool
        Whether to display plots inline (default: True)
        
    Returns
    -------
    dict
        Dictionary with results per interacting cell type:
        {
            "cell_type": {
                "stats_df": pd.DataFrame,
                "pdf_path": Path,
                "figs": dict of matplotlib figures (if show_plots=True)
            }
        }
    """
    from behav3d.utils.analysis import smooth_value_over_time
    
    print(f"--------------- Performing {cell_type} Interaction Analysis ---------------")
    start_time = time.time()
    
    # Setup paths
    output_dir = Path(output_dir)
    analysis_outdir = output_dir / "analysis" / cell_type
    feature_outdir = analysis_outdir / "track_features"
    results_dir = analysis_outdir / "interaction_analysis"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    if df_tracks_path is None:
        df_tracks_path = feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"
    else:
        df_tracks_path = Path(df_tracks_path)
    
    if not df_tracks_path.exists():
        raise FileNotFoundError(
            f"Filtered track features not found at: {df_tracks_path}\n"
            "Run Track Filtering (Step 1) first."
        )
    
    # Load data
    print("   Loading track features...")
    df = pd.read_csv(df_tracks_path)
    df = df.sort_values(by=["sample_name", "TrackID", "position_t"])
    df["TrackID"] = df["TrackID"].astype(str)
    
    n_organoids = df.groupby("sample_name")["TrackID"].nunique().sum()
    n_samples = df["sample_name"].nunique()
    print(f"   Loaded {len(df)} timepoints from {n_organoids} {cell_type}s across {n_samples} samples")
    
    # Process death classification
    df, has_dead_column = _process_death_classification(df, dead_threshold)
    
    if has_dead_column:
        n_alive = (df.groupby(["sample_name", "TrackID"])["survives"].first()).sum()
        n_dead = n_organoids - n_alive
        print(f"   Death status: {n_alive} survive, {n_dead} die (threshold={dead_threshold})")
    else:
        print("   ⚠️ No death data available - skipping alive vs dead comparisons")
    
    # Process each interacting cell type
    results = {}
    
    for ct in interacting_cell_types:
        contact_col = f"{ct}_contact"
        if contact_col not in df.columns:
            print(f"   ⚠️ Skipping {ct} - no contact column '{contact_col}' found")
            continue
        
        print(f"\n   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        print(f"   Processing: {cell_type} ↔ {ct}")
        print(f"   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        
        # Calculate cumulative contact
        df[contact_col] = df[contact_col].astype(int)
        cumulative_col = f"cumulative_{ct}_contact"
        df[cumulative_col] = df.groupby(["sample_name", "TrackID"])[contact_col].cumsum()
        
        # Calculate statistics
        stats_df = calculate_interaction_stats(df, cell_type, ct, cumulative_col, has_dead_column)
        
        # Save statistics CSV
        stats_csv_path = results_dir / f"interaction_stats_{cell_type}_vs_{ct}.csv"
        stats_df.to_csv(stats_csv_path, index=False)
        print(f"   📊 Statistics saved: {stats_csv_path.name}")
        
        # Generate plots and PDF
        figs = generate_interaction_plots(
            df=df,
            cell_type=cell_type,
            interacting_type=ct,
            cumulative_col=cumulative_col,
            has_dead_column=has_dead_column,
            n_organoids=n_organoids,
            n_samples=n_samples,
        )
        
        # Save PDF
        pdf_path = results_dir / f"interaction_analysis_{cell_type}_vs_{ct}.pdf"
        save_plots_to_pdf(figs, pdf_path)
        print(f"   📈 PDF saved: {pdf_path.name}")
        
        # Display overall plots if requested
        if show_plots:
            print(f"\n   📈 Displaying overall plots for {ct}:")
            from IPython.display import display
            
            if "cumulative_overall" in figs:
                display(figs["cumulative_overall"])
            if "alive_vs_dead_overall" in figs:
                display(figs["alive_vs_dead_overall"])
        
        # Close all figures
        for fig in figs.values():
            plt.close(fig)
        
        results[ct] = {
            "stats_df": stats_df,
            "pdf_path": pdf_path,
            "stats_csv_path": stats_csv_path,
        }
    
    elapsed = time.time() - start_time
    print(f"\n   ✅ Interaction Analysis complete! ({elapsed:.1f}s)")
    print(f"   Results saved to: {results_dir}")
    
    return results


def _process_death_classification(df: pd.DataFrame, dead_threshold: float):
    """
    Process death classification for the dataframe.
    
    Returns
    -------
    tuple
        (df with death columns, has_dead_column: bool)
    """
    from behav3d.utils.analysis import smooth_value_over_time
    
    has_dead_column = "dead" in df.columns
    has_percentage_dead = "percentage_dead_mask" in df.columns
    
    if not has_dead_column and has_percentage_dead:
        # Create dead classification based on threshold
        print(f"   Calculating death classification (threshold={dead_threshold})...")
        
        # Smooth the percentage_dead_mask
        df["smoothed_percentage_dead_mask"] = smooth_value_over_time(
            df, 
            column="percentage_dead_mask", 
            rolling_meanspeed_window=20,
            min_periods=20,
            groupby=["TrackID", "sample_name"]
        )
        
        df["dead"] = (df["smoothed_percentage_dead_mask"] > dead_threshold)
        df["dead"] = df.groupby(["sample_name", "TrackID"])["dead"].transform(lambda x: x.cummax())
        has_dead_column = True
    
    # Determine final alive/dead status for each organoid
    if has_dead_column:
        df_final_status = df.groupby(["sample_name", "TrackID"]).agg(
            final_dead=("dead", "last"),
            max_time=("position_t", "max")
        ).reset_index()
        df_final_status["survives"] = ~df_final_status["final_dead"]
        
        # Merge back
        df = df.merge(
            df_final_status[["sample_name", "TrackID", "survives"]], 
            on=["sample_name", "TrackID"], 
            how="left"
        )
    
    return df, has_dead_column


def calculate_interaction_stats(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
    has_dead_column: bool,
) -> pd.DataFrame:
    """
    Calculate summary statistics for interactions with a specific cell type.
    
    Parameters
    ----------
    df : pd.DataFrame
        Track features dataframe with contact columns
    cell_type : str
        Organoid cell type name
    interacting_type : str
        Interacting cell type name
    cumulative_col : str
        Name of cumulative contact column
    has_dead_column : bool
        Whether death data is available
        
    Returns
    -------
    pd.DataFrame
        Per-sample summary statistics
    """
    contact_col = f"{interacting_type}_contact"
    
    # Per-track summary
    track_stats = df.groupby(["sample_name", "TrackID"]).agg(
        total_contact_timepoints=(contact_col, "sum"),
        total_timepoints=("position_t", "count"),
        max_cumulative_contacts=(cumulative_col, "max"),
    ).reset_index()
    
    track_stats["contact_percentage"] = (
        track_stats["total_contact_timepoints"] / track_stats["total_timepoints"] * 100
    )
    
    # Add survival status if available
    if has_dead_column and "survives" in df.columns:
        survival_info = df.groupby(["sample_name", "TrackID"])["survives"].first().reset_index()
        track_stats = track_stats.merge(survival_info, on=["sample_name", "TrackID"], how="left")
    
    # Per-sample summary
    sample_stats = []
    for sample in df["sample_name"].unique():
        track_sample = track_stats[track_stats["sample_name"] == sample]
        
        stats = {
            "sample_name": sample,
            "organoid_type": cell_type,
            "interacting_cell_type": interacting_type,
            "n_organoids": track_sample["TrackID"].nunique(),
            "mean_contact_percentage": track_sample["contact_percentage"].mean(),
            "std_contact_percentage": track_sample["contact_percentage"].std(),
            "mean_total_contacts": track_sample["max_cumulative_contacts"].mean(),
            "std_total_contacts": track_sample["max_cumulative_contacts"].std(),
        }
        
        if has_dead_column and "survives" in track_sample.columns:
            stats["n_survives"] = track_sample["survives"].sum()
            stats["n_dies"] = (~track_sample["survives"]).sum()
            
            # Mean contacts for alive vs dead
            alive_contacts = track_sample[track_sample["survives"]]["max_cumulative_contacts"]
            dead_contacts = track_sample[~track_sample["survives"]]["max_cumulative_contacts"]
            
            stats["mean_contacts_survives"] = alive_contacts.mean() if len(alive_contacts) > 0 else np.nan
            stats["mean_contacts_dies"] = dead_contacts.mean() if len(dead_contacts) > 0 else np.nan
        
        sample_stats.append(stats)
    
    return pd.DataFrame(sample_stats)


def generate_interaction_plots(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
    has_dead_column: bool,
    n_organoids: int,
    n_samples: int,
) -> dict:
    """
    Generate all interaction plots.
    
    Returns
    -------
    dict
        Dictionary of matplotlib figures:
        - "cumulative_overall": Cumulative contacts over time (all samples)
        - "cumulative_per_sample": Cumulative contacts per sample
        - "alive_vs_dead_overall": Comparison by survival status (all samples)
        - "alive_vs_dead_per_sample": Comparison by survival status per sample
    """
    figs = {}
    
    # Plot 1: Cumulative overall
    figs["cumulative_overall"] = plot_cumulative_overall(
        df, cell_type, interacting_type, cumulative_col, n_organoids, n_samples
    )
    
    # Plot 2: Alive vs Dead overall
    if has_dead_column and "survives" in df.columns:
        fig = plot_alive_vs_dead_overall(
            df, cell_type, interacting_type, cumulative_col, n_organoids
        )
        if fig is not None:
            figs["alive_vs_dead_overall"] = fig
    
    # Plot 3: Cumulative per sample
    figs["cumulative_per_sample"] = plot_cumulative_per_sample(
        df, cell_type, interacting_type, cumulative_col
    )
    
    # Plot 4: Alive vs Dead per sample
    if has_dead_column and "survives" in df.columns:
        fig = plot_alive_vs_dead_per_sample(
            df, cell_type, interacting_type, cumulative_col
        )
        if fig is not None:
            figs["alive_vs_dead_per_sample"] = fig
    
    return figs


def save_plots_to_pdf(figs: dict, pdf_path: Path):
    """Save all figures to a PDF file."""
    with PdfPages(pdf_path) as pdf:
        # Order: overall plots first, then per-sample
        for key in ["cumulative_overall", "alive_vs_dead_overall", 
                    "cumulative_per_sample", "alive_vs_dead_per_sample"]:
            if key in figs:
                pdf.savefig(figs[key])


def plot_cumulative_overall(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
    n_organoids: int,
    n_samples: int,
) -> plt.Figure:
    """Plot cumulative interactions over time (all samples combined)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate mean and SEM per timepoint
    stats = df.groupby("position_t")[cumulative_col].agg(["mean", "std", "count"]).reset_index()
    stats["sem"] = stats["std"] / np.sqrt(stats["count"])
    
    ax.plot(stats["position_t"], stats["mean"], linewidth=2, color="steelblue")
    ax.fill_between(
        stats["position_t"],
        stats["mean"] - stats["sem"],
        stats["mean"] + stats["sem"],
        alpha=0.3,
        color="steelblue"
    )
    
    ax.set_xlabel("Time (frames)", fontsize=12)
    ax.set_ylabel(f"Cumulative {interacting_type} contacts", fontsize=12)
    ax.set_title(
        f"Average Cumulative {interacting_type} Interactions with {cell_type}s Over Time\n"
        f"(n = {n_organoids} {cell_type}s, {n_samples} samples)",
        fontsize=14
    )
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_cumulative_per_sample(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
) -> plt.Figure:
    """Plot cumulative interactions over time per sample."""
    samples = df["sample_name"].unique()
    n_cols = min(3, len(samples))
    n_rows = (len(samples) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows * n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, sample in enumerate(samples):
        ax = axes[idx]
        df_sample = df[df["sample_name"] == sample]
        n_org_sample = df_sample["TrackID"].nunique()
        
        stats_sample = df_sample.groupby("position_t")[cumulative_col].agg(
            ["mean", "std", "count"]
        ).reset_index()
        stats_sample["sem"] = stats_sample["std"] / np.sqrt(stats_sample["count"])
        
        ax.plot(stats_sample["position_t"], stats_sample["mean"], linewidth=2, color="steelblue")
        ax.fill_between(
            stats_sample["position_t"],
            stats_sample["mean"] - stats_sample["sem"],
            stats_sample["mean"] + stats_sample["sem"],
            alpha=0.3,
            color="steelblue"
        )
        
        ax.set_xlabel("Time", fontsize=10)
        ax.set_ylabel(f"Cumulative {interacting_type} contacts", fontsize=10)
        ax.set_title(f"{sample}\n(n = {n_org_sample})", fontsize=11)
        ax.grid(True, alpha=0.3)
    
    # Hide unused axes
    for idx in range(len(samples), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(
        f"Cumulative {interacting_type} Interactions with {cell_type}s Per Sample",
        fontsize=14, y=1.02
    )
    plt.tight_layout()
    return fig


def plot_alive_vs_dead_overall(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
    n_organoids: int,
) -> plt.Figure:
    """Plot alive vs dead comparison (all samples combined)."""
    if "survives" not in df.columns:
        return None
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for survives, label, color in [(True, "Survives", "forestgreen"), (False, "Dies", "crimson")]:
        df_subset = df[df["survives"] == survives]
        n_subset = df_subset["TrackID"].nunique()
        
        if n_subset == 0:
            continue
        
        stats = df_subset.groupby("position_t")[cumulative_col].agg(
            ["mean", "std", "count"]
        ).reset_index()
        stats["sem"] = stats["std"] / np.sqrt(stats["count"])
        
        ax.plot(
            stats["position_t"], stats["mean"], 
            linewidth=2, color=color, label=f"{label} (n={n_subset})"
        )
        ax.fill_between(
            stats["position_t"],
            stats["mean"] - stats["sem"],
            stats["mean"] + stats["sem"],
            alpha=0.2,
            color=color
        )
    
    ax.set_xlabel("Time (frames)", fontsize=12)
    ax.set_ylabel(f"Cumulative {interacting_type} contacts", fontsize=12)
    ax.set_title(
        f"Cumulative {interacting_type} Interactions: Surviving vs Dying {cell_type}s\n"
        f"(Total n = {n_organoids})",
        fontsize=14
    )
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def plot_alive_vs_dead_per_sample(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
) -> plt.Figure:
    """Plot alive vs dead comparison per sample."""
    if "survives" not in df.columns:
        return None
    
    samples = df["sample_name"].unique()
    n_cols = min(3, len(samples))
    n_rows = (len(samples) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows * n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for idx, sample in enumerate(samples):
        ax = axes[idx]
        df_sample = df[df["sample_name"] == sample]
        
        for survives, label, color in [(True, "Survives", "forestgreen"), (False, "Dies", "crimson")]:
            df_subset = df_sample[df_sample["survives"] == survives]
            n_subset = df_subset["TrackID"].nunique()
            
            if n_subset == 0:
                continue
            
            stats = df_subset.groupby("position_t")[cumulative_col].agg(
                ["mean", "std", "count"]
            ).reset_index()
            stats["sem"] = stats["std"] / np.sqrt(stats["count"])
            
            ax.plot(
                stats["position_t"], stats["mean"],
                linewidth=2, color=color, label=f"{label} (n={n_subset})"
            )
            ax.fill_between(
                stats["position_t"],
                stats["mean"] - stats["sem"],
                stats["mean"] + stats["sem"],
                alpha=0.2,
                color=color
            )
        
        ax.set_xlabel("Time", fontsize=10)
        ax.set_ylabel(f"Cumulative {interacting_type} contacts", fontsize=10)
        ax.set_title(f"{sample}", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    # Hide unused axes
    for idx in range(len(samples), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(
        f"Cumulative {interacting_type} Interactions: Surviving vs Dying {cell_type}s (Per Sample)",
        fontsize=14, y=1.02
    )
    plt.tight_layout()
    return fig
