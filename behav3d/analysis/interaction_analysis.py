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
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


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
        Retained for backward compatibility. Interaction analysis now reads the
        final death classification directly from the ``dead`` column in the CSV.
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
    
    # Process death classification from the final feature-extraction output.
    df, has_dead_column = _process_death_classification(df)
    
    if has_dead_column:
        n_alive = (df.groupby(["sample_name", "TrackID"])["survives"].first()).sum()
        n_dead = n_organoids - n_alive
        print(f"   Death status: {n_alive} survive, {n_dead} die")
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


def _process_death_classification(df: pd.DataFrame):
    """
    Process death classification for the dataframe.
    
    Returns
    -------
    tuple
        (df with death columns, has_dead_column: bool)
    """
    
    
    has_dead_column = "dead" in df.columns
    if not has_dead_column:
        return df, False

    if pd.api.types.is_bool_dtype(df["dead"]):
        df["dead"] = df["dead"].fillna(False)
    else:
        df["dead"] = (
            df["dead"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "1.0", "yes"})
        )
    
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
            "n_organoids": track_sample.groupby(["sample_name", "TrackID"]).ngroups,
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
        n_org_sample = df_sample.groupby(["sample_name", "TrackID"]).ngroups
        
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
        n_subset = df_subset.groupby(["sample_name", "TrackID"]).ngroups
        
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
            n_subset = df_subset.groupby(["sample_name", "TrackID"]).ngroups
                
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


# ---------------------------------------------------------------------------
# Multi-Organoid Interaction Comparison
# ---------------------------------------------------------------------------

def run_multi_organoid_interaction_comparison(
    output_dir: str,
    organoid_types: list,
    immune_types: list,
    metadata: pd.DataFrame,
    group_by: str = "organoid_type",
    time_window_min: float = 60,
    show_plots: bool = True,
):
    """
    Compare interactions across multiple organoid types.

    Parameters
    ----------
    output_dir : str
        Base output directory for BEHAV3D results.
    organoid_types : list
        Organoid type names (e.g. ["healthy", "tumor"]).
    immune_types : list
        Immune / other cell type names to include (e.g. ["cd4", "cd8"]).
    metadata : pd.DataFrame
        Full metadata table (used to resolve time_interval).
    group_by : str
        "organoid_type" or "treatment".
    time_window_min : float
        Time window in minutes before TOD for the cumulative curve.
    show_plots : bool
        Whether to display figures inline.

    Returns
    -------
    dict  with keys "violin_fig", "curve_fig", "pdf_path", "summary_csv_path".
    """
    print("--------------- Multi-Organoid Interaction Comparison ---------------")
    start = time.time()

    output_dir = Path(output_dir)
    results_dir = output_dir / "analysis" / "multi_organoid_comparison"
    results_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load filtered tracks for every organoid type
    # ------------------------------------------------------------------
    frames = []
    for org_type in organoid_types:
        csv_path = (
            output_dir / "analysis" / org_type / "track_features"
            / f"BEHAV3D_{org_type}_combined_track_features_filtered.csv"
        )
        if not csv_path.exists():
            print(f"   Skipping {org_type}: filtered tracks not found")
            continue

        df = pd.read_csv(csv_path)
        df = df.sort_values(by=["sample_name", "TrackID", "position_t"])
        df["TrackID"] = df["TrackID"].astype(str)
        df["organoid_type"] = org_type

        # Keep only immune types whose contact column exists
        available_im = [ct for ct in immune_types if f"{ct}_contact" in df.columns]
        if not available_im:
            print(f"   Skipping {org_type}: no contact columns for {immune_types}")
            continue

        for ct in available_im:
            col = f"{ct}_contact"
            df[col] = df[col].astype(int)
            df[f"cumulative_{ct}_contact"] = (
                df.groupby(["sample_name", "TrackID"])[col].cumsum()
            )

        df, has_dead = _process_death_classification(df)

        # Compute t_dead per track (first position_t where dead)
        if has_dead:
            t_dead_map = (
                df[df["dead"]]
                .groupby(["sample_name", "TrackID"])["position_t"]
                .min()
                .rename("t_dead")
            )
            df = df.merge(
                t_dead_map.reset_index(),
                on=["sample_name", "TrackID"],
                how="left",
            )
        else:
            df["t_dead"] = np.nan
            df["survives"] = True

        frames.append(df)

    if not frames:
        print("   No data available for any organoid type.")
        return None

    df_all = pd.concat(frames, ignore_index=True)

    n_organoids = df_all.groupby(["organoid_type", "sample_name", "TrackID"]).ngroups
    print(f"   Loaded {n_organoids} organoids across {len(frames)} organoid type(s)")

    # Resolve which immune types actually have data
    active_immune = [ct for ct in immune_types if f"cumulative_{ct}_contact" in df_all.columns]
    if not active_immune:
        print("   No contact columns found across loaded data.")
        return None

    # ------------------------------------------------------------------
    # 2. Build per-track summary for violin
    # ------------------------------------------------------------------
    track_keys = ["organoid_type", "sample_name", "TrackID"]
    track_summary = df_all.groupby(track_keys).agg(
        survives=("survives", "first") if "survives" in df_all.columns else ("position_t", lambda _: True),
        t_dead=("t_dead", "first"),
        t_last=("position_t", "max"),
    ).reset_index()

    for ct in active_immune:
        cum_col = f"cumulative_{ct}_contact"
        max_vals = (
            df_all.groupby(track_keys)[cum_col]
            .max()
            .rename(f"max_{ct}_contact")
            .reset_index()
        )
        track_summary = track_summary.merge(max_vals, on=track_keys, how="left")

    # Total contacts across all selected immune types (for "by organoid_type" grouping)
    im_max_cols = [f"max_{ct}_contact" for ct in active_immune]
    track_summary["total_contacts"] = track_summary[im_max_cols].sum(axis=1)
    track_summary["fate"] = track_summary["survives"].map({True: "Live", False: "Dying"})

    # Save summary
    summary_csv = results_dir / "multi_organoid_interaction_summary.csv"
    track_summary.to_csv(summary_csv, index=False)
    print(f"   Summary saved: {summary_csv.name}")

    # ------------------------------------------------------------------
    # 3. Generate plots
    # ------------------------------------------------------------------
    figs = {}

    has_dead_data = "survives" in track_summary.columns and not track_summary["survives"].all()

    fig_violin = plot_interaction_violin_comparison(
        track_summary, group_by, active_immune, has_dead_data,
    )
    figs["violin"] = fig_violin

    fig_curve = plot_cumulative_to_death_curves(
        df_all, group_by, active_immune, time_window_min, track_summary,
    )
    if fig_curve is not None:
        figs["cumulative_to_death"] = fig_curve

    # Save PDF
    pdf_path = results_dir / "multi_organoid_interaction_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        for key in ["violin", "cumulative_to_death"]:
            if key in figs:
                pdf.savefig(figs[key])
    print(f"   PDF saved: {pdf_path.name}")

    if show_plots:
        from IPython.display import display
        for fig in figs.values():
            display(fig)

    for fig in figs.values():
        plt.close(fig)

    elapsed = time.time() - start
    print(f"\n   Multi-Organoid Interaction Comparison complete! ({elapsed:.1f}s)")
    return {
        "summary_csv_path": summary_csv,
        "pdf_path": pdf_path,
    }


def plot_interaction_violin_comparison(
    track_summary: pd.DataFrame,
    group_by: str,
    immune_types: list,
    has_dead_data: bool,
) -> plt.Figure:
    if group_by == "treatment":
        rows = []
        for ct in immune_types:
            col = f"max_{ct}_contact"
            if col not in track_summary.columns:
                continue
            sub = track_summary[["organoid_type", "sample_name", "TrackID", "fate", col]].copy()
            sub = sub.rename(columns={col: "cumulative_interactions"})
            sub["condition"] = ct
            rows.append(sub)
        df_plot = pd.concat(rows, ignore_index=True)
        x_label = "Treatment (immune cell type)"
    else:
        df_plot = track_summary.copy()
        df_plot["cumulative_interactions"] = df_plot["total_contacts"]
        df_plot["condition"] = df_plot["organoid_type"]
        x_label = "Organoid type"

    conditions = df_plot["condition"].unique()
    n_conditions = len(conditions)

    fig, ax = plt.subplots(figsize=(max(6, 3 * n_conditions), 7))

    order = sorted(conditions)

    if has_dead_data:
        hue_col = "fate"
        hue_order = ["Dying", "Live"]
        palette = {"Dying": "#e74c3c", "Live": "#3498db"}
        # Darker versions for boxes
        dark_palette = {"Dying": "#922b21", "Live": "#1a5276"}
    else:
        hue_col = None
        hue_order = None
        palette = None
        dark_palette = {}

    common_kw = dict(
        data=df_plot, x="condition", y="cumulative_interactions",
        hue=hue_col, hue_order=hue_order, order=order, palette=palette,
        ax=ax, dodge=True, width=0.8,
    )

    sns.violinplot(
        **common_kw,
        inner=None, linewidth=1, alpha=0.35,
        density_norm="count", cut=0,
    )

    # Use a copy of common_kw with the darker palette + gap for narrow centered box
    box_kw = common_kw.copy()
    box_kw["palette"] = dark_palette
    sns.boxplot(
        **box_kw,
        gap=0.6,          # Large gap makes the box narrow
        fliersize=0, legend=False,
        boxprops=dict(alpha=0.8),
        whiskerprops=dict(alpha=0.8),
        capprops=dict(alpha=0.8),
        showcaps=True,
        medianprops=dict(color="white", linewidth=2),
    )

    # Stripplot doesn't support width; remove it to avoid crash
    strip_kw = common_kw.copy()
    strip_kw.pop("width", None)
    sns.stripplot(
        **strip_kw,
        size=3, alpha=0.5, jitter=True, legend=False,
    )

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        n_unique = len(hue_order) if hue_order else 1
        ax.legend(handles[:n_unique], labels[:n_unique], fontsize=11, loc="best")

    # Annotate zero counts
    if has_dead_data:
        y_range = ax.get_ylim()[1] - max(ax.get_ylim()[0], 0)
        annot_y = y_range * 0.02
        for cond in order:
            for fate_val in ["Dying", "Live"]:
                sub = df_plot[(df_plot["condition"] == cond) & (df_plot["fate"] == fate_val)]
                n_zeros = (sub["cumulative_interactions"] == 0).sum()
                if n_zeros > 0:
                    x_idx = list(order).index(cond)
                    offset = -0.2 if fate_val == "Dying" else 0.2
                    ax.annotate(
                        f"{n_zeros} zeros",
                        xy=(x_idx + offset, annot_y),
                        fontsize=9, ha="center", va="bottom",
                        color="#000000", # Pure black for maximum visibility
                        fontweight="normal",
                        fontstyle="italic",
                    )

    # X-axis
    tick_labels = []
    for c in order:
        cond_data = df_plot[df_plot["condition"] == c]
        n_total = len(cond_data)
        if has_dead_data:
            n_dying = (cond_data["fate"] == "Dying").sum()
            pct = (n_dying / n_total * 100) if n_total > 0 else 0
            tick_labels.append(f"{c}\n({n_dying}/{n_total} dying, {pct:.0f}%)")
        else:
            tick_labels.append(f"{c}\n(n={n_total})")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(tick_labels)

    ax.set_ylim(bottom=0)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel("Cumulative interactions before death / end", fontsize=12)
    title_detail = " + ".join(immune_types) if group_by == "organoid_type" else "all organoid types"
    ax.set_title(
        f"Interaction Distribution ({title_detail})",
        fontsize=14,
    )
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    return fig


def plot_cumulative_to_death_curves(
    df_all: pd.DataFrame,
    group_by: str,
    immune_types: list,
    time_window_min: float,
    track_summary: pd.DataFrame = None,
) -> plt.Figure:
    """
    Mean +/- SEM cumulative interaction curves aligned to time of death.

    Only dying organoids are included.  X-axis runs from ``-time_window_min``
    to 0 (TOD) with three tick marks.  Per-track curves are interpolated onto
    a common time grid before averaging to produce smooth output.

    Parameters
    ----------
    df_all : pd.DataFrame
        Combined timepoint-level data with ``organoid_type``, ``t_dead``,
        ``survives``, ``time`` (hours), and ``cumulative_{ct}_contact`` columns.
    group_by : str
        "organoid_type" or "treatment".
    immune_types : list
        Immune cell type names with data.
    time_window_min : float
        Time window in minutes before death.
    track_summary : pd.DataFrame, optional
        Per-track summary (one row per organoid) used to compute total counts
        for the legend labels.
    """
    if "survives" not in df_all.columns:
        print("   No death data -- skipping cumulative-to-death curves.")
        return None

    df_dying = df_all[df_all["survives"] == False].copy()
    if df_dying.empty:
        print("   No dying organoids -- skipping cumulative-to-death curves.")
        return None

    if "time" not in df_dying.columns:
        print("   No 'time' column (hours) -- skipping cumulative-to-death curves.")
        return None

    track_keys = ["organoid_type", "sample_name", "TrackID"]

    # Compute the real time at death for each track
    t_dead_time = (
        df_dying[df_dying["dead"]]
        .groupby(track_keys)["time"]
        .min()
        .rename("time_at_death")
    )
    df_dying = df_dying.merge(t_dead_time.reset_index(), on=track_keys, how="left")

    # Filter to window before death
    df_dying["relative_time_min"] = (df_dying["time"] - df_dying["time_at_death"]) * 60.0
    df_dying = df_dying[
        (df_dying["relative_time_min"] >= -time_window_min)
        & (df_dying["relative_time_min"] <= 0)
    ]

    if df_dying.empty:
        print("   No timepoints within the specified window -- skipping curve.")
        return None

    # Recompute cumulative contacts within the window only
    for ct in immune_types:
        contact_col = f"{ct}_contact"
        if contact_col not in df_dying.columns:
            continue
        df_dying[f"window_cum_{ct}"] = (
            df_dying
            .sort_values(track_keys + ["position_t"])
            .groupby(track_keys)[contact_col]
            .cumsum()
        )

    # Build long-form data for plotting
    if group_by == "treatment":
        curve_frames = []
        for ct in immune_types:
            cum_col = f"window_cum_{ct}"
            if cum_col not in df_dying.columns:
                continue
            sub = df_dying[track_keys + ["relative_time_min", cum_col]].copy()
            sub = sub.rename(columns={cum_col: "cumulative_interactions"})
            sub["condition"] = ct
            curve_frames.append(sub)
        if not curve_frames:
            return None
        df_curve = pd.concat(curve_frames, ignore_index=True)
    else:
        cum_cols = [f"window_cum_{ct}" for ct in immune_types if f"window_cum_{ct}" in df_dying.columns]
        if not cum_cols:
            return None
        df_dying["cumulative_interactions"] = df_dying[cum_cols].sum(axis=1)
        df_curve = df_dying[track_keys + ["relative_time_min", "cumulative_interactions"]].copy()
        df_curve["condition"] = df_curve["organoid_type"]

    # ---- Interpolate per-track curves onto a common time grid ----
    n_grid = max(int(time_window_min), 30)
    time_grid = np.linspace(-time_window_min, 0, n_grid)

    interp_rows = []
    for (cond, sn, tid), grp in df_curve.groupby(["condition", "sample_name", "TrackID"]):
        grp = grp.sort_values("relative_time_min")
        t_vals = grp["relative_time_min"].values
        y_vals = grp["cumulative_interactions"].values
        if len(t_vals) < 2:
            continue
        y_interp = np.interp(time_grid, t_vals, y_vals, left=0.0, right=y_vals[-1])
        for t, y in zip(time_grid, y_interp):
            interp_rows.append({"condition": cond, "time_grid": t, "cumulative_interactions": y})

    if not interp_rows:
        print("   Not enough data for interpolation -- skipping curve.")
        return None

    df_interp = pd.DataFrame(interp_rows)

    # ---- Compute totals for legend ----
    total_per_cond = {}
    if track_summary is not None:
        if group_by == "treatment":
            for ct in immune_types:
                total_per_cond[ct] = len(track_summary)
        else:
            for org in track_summary["organoid_type"].unique():
                total_per_cond[org] = (track_summary["organoid_type"] == org).sum()

    # ---- Plot ----
    conditions = sorted(df_interp["condition"].unique())
    n_conds = len(conditions)
    palette = sns.color_palette("tab10", n_conds)
    color_map = dict(zip(conditions, palette))

    fig, ax = plt.subplots(figsize=(10, 6))

    for cond in conditions:
        sub = df_interp[df_interp["condition"] == cond]
        stats = (
            sub.groupby("time_grid")["cumulative_interactions"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        stats["sem"] = stats["std"] / np.sqrt(stats["count"])
        stats = stats.sort_values("time_grid")

        n_dying = int(stats["count"].iloc[0]) if len(stats) else 0
        n_total = total_per_cond.get(cond, n_dying)
        pct = (n_dying / n_total * 100) if n_total > 0 else 0
        label = f"{cond} ({n_dying}/{n_total} dying, {pct:.0f}%)"

        color = color_map[cond]
        ax.plot(
            stats["time_grid"], stats["mean"],
            linewidth=2, color=color, label=label,
        )
        ax.fill_between(
            stats["time_grid"],
            np.maximum(stats["mean"] - stats["sem"], 0),
            stats["mean"] + stats["sem"],
            alpha=0.2, color=color,
        )

    ax.set_ylim(bottom=0)
    ax.set_xticks([-time_window_min, -time_window_min / 2, 0])
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: "TOD" if v == 0 else f"{v:.0f}"
    ))
    ax.set_xlabel("Time relative to death (min)", fontsize=12)
    ax.set_ylabel("Cumulative interactions", fontsize=12)
    ax.set_title(
        f"Cumulative Interactions Before Death (window = {time_window_min:.0f} min)",
        fontsize=14,
    )
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig
