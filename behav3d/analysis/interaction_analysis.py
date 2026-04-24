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
    # 2b. Load active killing data (if available)
    # ------------------------------------------------------------------
    df_contact_events = _load_active_killing_data(
        output_dir, active_immune, organoid_types,
        track_summary=track_summary,
    )
    has_killing_data = not df_contact_events.empty
    if has_killing_data:
        _n_ev = df_contact_events["contact_event_id"].nunique()
        _n_targeted = int(df_contact_events["has_active_killing"].sum())
        _n_kill_events = (
            df_contact_events.loc[
                df_contact_events["has_active_killing"],
                "contact_event_id",
            ].nunique()
        )
        print(f"   Active killing data loaded: {_n_ev} contact events "
              f"({_n_kill_events} with active killing, "
              f"{_n_targeted} targeted organoid hits)")
    else:
        print("   \u2139\ufe0f  No active killing data found \u2014 "
              "panels C & D require running Active Killing analysis first")

    # ------------------------------------------------------------------
    # 3. Generate plots
    # ------------------------------------------------------------------
    figs = {}

    has_dead_data = "survives" in track_summary.columns and not track_summary["survives"].all()

    fig_violin = plot_interaction_violin_comparison(
        track_summary, group_by, active_immune, has_dead_data,
    )
    figs["violin"] = fig_violin

    curve_data_container = []
    fig_curve = plot_cumulative_to_death_curves(
        df_all,
        group_by,
        active_immune,
        time_window_min,
        track_summary,
        curve_data_out=curve_data_container,
    )
    if fig_curve is not None:
        figs["cumulative_to_death"] = fig_curve

    if curve_data_container:
        curve_csv = results_dir / "multi_organoid_cumulative_to_death_curves_min.csv"
        curve_df = curve_data_container[0]
        curve_df.to_csv(curve_csv, index=False)
        print(f"   Curve data saved (time in minutes): {curve_csv.name}")

    # Active killing dashboard (only generated when killing data exists)
    fig_dashboard = plot_interaction_overview_dashboard(
        track_summary, df_contact_events, active_immune, has_dead_data,
    )
    if fig_dashboard is not None:
        figs["dashboard"] = fig_dashboard

    # Save PDF. Order matches the notebook display order (violin and
    # cumulative-to-death first, active-killing dashboard last).
    pdf_path = results_dir / "multi_organoid_interaction_comparison.pdf"
    with PdfPages(pdf_path) as pdf:
        for key in ["violin", "cumulative_to_death", "dashboard"]:
            if key in figs:
                pdf.savefig(
                    figs[key],
                    bbox_inches="tight",
                    pad_inches=0.25,
                )
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

    fig, ax = plt.subplots(figsize=(max(6, 3.5 * n_conditions), 7), dpi=120)
    fig.patch.set_facecolor("#fafafa")
    ax.set_facecolor("#fafafa")

    order = sorted(conditions)

    if has_dead_data:
        hue_col = "fate"
        # Internal fate labels ("Dying" / "Live") stay as-is to match the
        # track_summary column; display labels are set via the legend below.
        # Order is Dead first, Live second (matches dashboard panels C/D).
        hue_order = ["Dying", "Live"]
        palette = {"Dying": "#e8a0a0", "Live": "#a0c4e8"}
        # Darker, richer versions for boxes
        dark_palette = {"Dying": "#b03a3a", "Live": "#2b6a99"}
        strip_palette = {"Dying": "#c0392b", "Live": "#2471a3"}
    else:
        hue_col = None
        hue_order = None
        palette = None
        dark_palette = {}
        strip_palette = None

    common_kw = dict(
        data=df_plot, x="condition", y="cumulative_interactions",
        hue=hue_col, hue_order=hue_order, order=order, palette=palette,
        ax=ax, dodge=True, width=0.8,
    )

    sns.violinplot(
        **common_kw,
        inner=None, linewidth=0.8, alpha=0.30,
        density_norm="count", cut=0,
    )

    # Use a copy of common_kw with the darker palette + gap for narrow centered box
    box_kw = common_kw.copy()
    box_kw["palette"] = dark_palette
    sns.boxplot(
        **box_kw,
        gap=0.6,
        fliersize=0, legend=False,
        boxprops=dict(alpha=0.85, linewidth=0.8),
        whiskerprops=dict(alpha=0.8, linewidth=0.8),
        capprops=dict(alpha=0.8, linewidth=0.8),
        showcaps=True,
        medianprops=dict(color="white", linewidth=2),
    )

    # Stripplot — slightly more muted dots
    strip_kw = common_kw.copy()
    strip_kw.pop("width", None)
    if strip_palette:
        strip_kw["palette"] = strip_palette
    sns.stripplot(
        **strip_kw,
        size=3, alpha=0.45, jitter=True, legend=False,
        edgecolor="white", linewidth=0.3,
    )

    # Legend -- framed, semi-transparent. Re-label "Dying" -> "Dead" at
    # display time so internal track_summary values stay unchanged.
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        n_unique = len(hue_order) if hue_order else 1
        display_labels = [
            _FATE_DISPLAY.get(l, l) for l in labels[:n_unique]
        ] if hue_order else labels[:n_unique]
        leg = ax.legend(
            handles[:n_unique], display_labels,
            fontsize=10, loc="best",
            frameon=True, framealpha=0.85, edgecolor="#cccccc",
            fancybox=True, shadow=False,
            title="Organoid fate" if hue_order else None,
            title_fontsize=9,
        )
        leg.get_frame().set_linewidth(0.6)

    # Annotate organoids with zero contacts with a subtle badge.
    if has_dead_data:
        y_range = ax.get_ylim()[1] - max(ax.get_ylim()[0], 0)
        annot_y = y_range * 0.015
        for cond in order:
            for fate_val in ["Dying", "Live"]:
                sub = df_plot[
                    (df_plot["condition"] == cond)
                    & (df_plot["fate"] == fate_val)
                ]
                n_zeros = (sub["cumulative_interactions"] == 0).sum()
                if n_zeros > 0:
                    x_idx = list(order).index(cond)
                    offset = -0.2 if fate_val == "Dying" else 0.2
                    badge_color = dark_palette.get(fate_val, "#555555")
                    ax.annotate(
                        f" {n_zeros} no contacts ",
                        xy=(x_idx + offset, annot_y),
                        fontsize=7.5, ha="center", va="bottom",
                        color="white",
                        fontweight="bold",
                        bbox=dict(
                            boxstyle="round,pad=0.25",
                            facecolor=badge_color,
                            edgecolor="none",
                            alpha=0.75,
                        ),
                    )

    # X-axis tick labels: show dead / total per condition.
    tick_labels = []
    for c in order:
        cond_data = df_plot[df_plot["condition"] == c]
        n_total = len(cond_data)
        if has_dead_data:
            n_dead = (cond_data["fate"] == "Dying").sum()
            pct = (n_dead / n_total * 100) if n_total > 0 else 0
            tick_labels.append(
                f"{c}\n({n_dead}/{n_total} dead, {pct:.0f}%)",
            )
        else:
            tick_labels.append(f"{c}\n(n={n_total})")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(tick_labels, fontsize=10)

    ax.set_ylim(bottom=0)
    ax.set_xlabel(x_label, fontsize=12, labelpad=8)
    # y-label: cumulative contact timepoints per organoid (integrated over
    # the full observation window; for dead organoids the cumulative stops
    # at time of death, matching how track features are computed).
    ax.set_ylabel(
        "Cumulative contact timepoints per organoid",
        fontsize=12, labelpad=8,
    )
    title_detail = (
        " + ".join(immune_types) if group_by == "organoid_type"
        else "all organoid types"
    )
    # Title makes the unit of observation explicit: one point = one
    # organoid, not one immune cell.
    ax.set_title(
        f"Cumulative interactions per organoid -- contacts with "
        f"{title_detail}",
        fontsize=14, fontweight="semibold", pad=12,
    )

    # Clean up spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["left"].set_color("#888888")
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["bottom"].set_color("#888888")
    ax.tick_params(axis="both", which="both", length=4, width=0.6, colors="#555555")

    ax.grid(True, axis="y", alpha=0.25, linewidth=0.5, color="#cccccc")
    ax.set_axisbelow(True)

    plt.tight_layout()
    return fig


def _time_unit_to_minutes_factor(df: pd.DataFrame) -> float:
    """
    Return the multiplier that converts values in the ``time`` column to minutes.

    Reads the ``time_unit`` column if present (expected "s", "m" or "h").
    Falls back to ``"h"`` (the legacy default of ``generalize_units_of_track_features``)
    if the column is missing, so upstream CSVs without an explicit unit keep
    working.
    """
    unit = "h"
    if "time_unit" in df.columns:
        vals = df["time_unit"].dropna().astype(str).str.lower().unique()
        if len(vals) == 1:
            unit = vals[0]
        elif len(vals) > 1:
            # Mixed units across tracks: refuse to guess silently.
            raise ValueError(
                f"Mixed time_unit values in data: {sorted(vals)}. "
                "Re-run feature extraction with a consistent time unit."
            )
    factors = {"s": 1.0 / 60.0, "m": 1.0, "h": 60.0}
    if unit not in factors:
        raise ValueError(
            f"Unsupported time_unit '{unit}'. Expected one of {sorted(factors)}."
        )
    return factors[unit]


def plot_cumulative_to_death_curves(
    df_all: pd.DataFrame,
    group_by: str,
    immune_types: list,
    time_window_min: float,
    track_summary: pd.DataFrame = None,
    curve_data_out: list = None,
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
        ``survives``, ``time`` and ``cumulative_{ct}_contact`` columns. The
        ``time`` column may be in seconds, minutes or hours; its unit is
        detected from the ``time_unit`` column (defaults to ``"h"`` for
        backward compatibility) and converted to minutes internally.
    group_by : str
        "organoid_type" or "treatment".
    immune_types : list
        Immune cell type names with data.
    time_window_min : float
        Time window in minutes before death.
    track_summary : pd.DataFrame, optional
        Per-track summary (one row per organoid) used to compute total counts
        for the legend labels.
    curve_data_out : list, optional
        If provided, the per-condition aggregated curve (time in minutes,
        mean, sem, n) is appended as a single DataFrame to this list so the
        caller can persist it to disk.
    """
    if "survives" not in df_all.columns:
        print("   No death data -- skipping cumulative-to-death curves.")
        return None

    df_dying = df_all[df_all["survives"] == False].copy()
    if df_dying.empty:
        print("   No dying organoids -- skipping cumulative-to-death curves.")
        return None

    if "time" not in df_dying.columns:
        print("   No 'time' column -- skipping cumulative-to-death curves.")
        return None

    track_keys = ["organoid_type", "sample_name", "TrackID"]

    to_min = _time_unit_to_minutes_factor(df_dying)

    # Compute the real time at death for each track (in minutes)
    t_dead_time = (
        df_dying[df_dying["dead"]]
        .groupby(track_keys)["time"]
        .min()
        .rename("time_at_death")
    )
    df_dying = df_dying.merge(t_dead_time.reset_index(), on=track_keys, how="left")

    # Time relative to each track's time-of-death, in minutes, regardless of
    # the upstream ``time_unit``. ``to_min`` maps s/m/h onto a minutes scale.
    df_dying["relative_time_min"] = (df_dying["time"] - df_dying["time_at_death"]) * to_min
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
    skipped_tracks = 0
    for (cond, sn, tid), grp in df_curve.groupby(["condition", "sample_name", "TrackID"]):
        grp = grp.sort_values("relative_time_min")
        t_vals = grp["relative_time_min"].values
        y_vals = grp["cumulative_interactions"].values
        if len(t_vals) == 0:
            # No observations in window -- still contribute a flat zero curve so
            # this dying organoid is counted in the legend, consistent with the
            # violin plot.
            y_interp = np.zeros_like(time_grid)
        elif len(t_vals) == 1:
            # Single observation inside the window: flat step at y_vals[0] from
            # that timepoint onward; 0 before it.
            y_interp = np.where(time_grid < t_vals[0], 0.0, float(y_vals[0]))
        else:
            y_interp = np.interp(time_grid, t_vals, y_vals, left=0.0, right=y_vals[-1])
        for t, y in zip(time_grid, y_interp):
            interp_rows.append({"condition": cond, "time_grid": t, "cumulative_interactions": y})
    if skipped_tracks:
        print(f"   Skipped {skipped_tracks} track(s) with no timepoints in window.")

    if not interp_rows:
        print("   Not enough data for interpolation -- skipping curve.")
        return None

    df_interp = pd.DataFrame(interp_rows)

    # ---- Compute totals for legend ----
    total_per_cond = {}
    dying_per_cond = {}
    if track_summary is not None:
        dying_mask = track_summary["fate"] == "Dying" if "fate" in track_summary.columns else ~track_summary["survives"].fillna(True)
        if group_by == "treatment":
            for ct in immune_types:
                total_per_cond[ct] = len(track_summary)
                dying_per_cond[ct] = int(dying_mask.sum())
        else:
            for org in track_summary["organoid_type"].unique():
                org_mask = track_summary["organoid_type"] == org
                total_per_cond[org] = int(org_mask.sum())
                dying_per_cond[org] = int((org_mask & dying_mask).sum())

    # ---- Plot ----
    conditions = sorted(df_interp["condition"].unique())
    n_conds = len(conditions)
    palette = sns.color_palette("tab10", n_conds)
    color_map = dict(zip(conditions, palette))

    fig, ax = plt.subplots(figsize=(10, 6))

    aggregated_frames = []

    for cond in conditions:
        sub = df_interp[df_interp["condition"] == cond]
        stats = (
            sub.groupby("time_grid")["cumulative_interactions"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        stats["sem"] = stats["std"] / np.sqrt(stats["count"])
        stats = stats.sort_values("time_grid")

        curve_count = int(stats["count"].iloc[0]) if len(stats) else 0
        n_dying = dying_per_cond.get(cond, curve_count)
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

        out = stats.rename(
            columns={"time_grid": "relative_time_min", "count": "n_tracks"}
        ).copy()
        out.insert(0, "condition", cond)
        out["n_dying_total"] = n_dying
        out["n_organoids_total"] = n_total
        aggregated_frames.append(out)

    if curve_data_out is not None and aggregated_frames:
        curve_data_out.append(pd.concat(aggregated_frames, ignore_index=True))

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


# ---------------------------------------------------------------------------
# Interaction Overview Dashboard
# ---------------------------------------------------------------------------

def _load_active_killing_data(
    output_dir: Path,
    immune_types: list,
    organoid_types: list,
    track_summary: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Load contact events from active-killing CSVs and link each event to
    its target organoid type (and fate, when ``track_summary`` is given).

    For each immune type whose active-killing output exists, reads
    ``contact_events_{im}.csv`` and ``active_killing_per_timepoint_{im}.csv``,
    determines which events contain active killing, and explodes the
    ``target_track_ids`` field so that each row represents one
    (contact_event \u00d7 target_organoid) pair.

    Parameters
    ----------
    track_summary : pd.DataFrame, optional
        Per-organoid summary with columns ``sample_name``, ``TrackID`` and
        ``fate`` (``"Live"`` / ``"Dying"``). When provided, the returned
        DataFrame gains a ``fate`` column per (event, target-organoid) row
        used by the dashboard panels to split by fate. When missing, ``fate``
        defaults to ``"Live"`` so downstream plotting still works.

    Returns
    -------
    pd.DataFrame
        Columns: contact_event_id, sample_name, immune_track_id,
        immune_type, target_track_id, organoid_type, contact_duration,
        has_active_killing, n_killing_tp, fate.
        Empty DataFrame if no active-killing data is found.
    """
    # Build (sample_name, TrackID-str) -> organoid_type lookup
    track_type_map: dict = {}
    for org_type in organoid_types:
        csv_path = (
            output_dir / "analysis" / org_type / "track_features"
            / f"BEHAV3D_{org_type}_combined_track_features_filtered.csv"
        )
        if not csv_path.exists():
            continue
        df_org = pd.read_csv(csv_path, usecols=["sample_name", "TrackID"])
        for _, row in (
            df_org[["sample_name", "TrackID"]].drop_duplicates().iterrows()
        ):
            key = (str(row["sample_name"]), str(row["TrackID"]))
            track_type_map.setdefault(key, org_type)

    all_events: list = []

    for im_type in immune_types:
        ak_dir = output_dir / "analysis" / im_type / "active_killing"
        events_path = ak_dir / f"contact_events_{im_type}.csv"
        killing_path = ak_dir / f"active_killing_per_timepoint_{im_type}.csv"

        if not events_path.exists():
            continue

        df_ev = pd.read_csv(events_path)
        if df_ev.empty:
            continue
        df_ev["immune_type"] = im_type

        # Per-target active-killing attribution from the per-timepoint CSV.
        #
        # ``advanced_timepoint_features.analyze_active_killing_per_timepoint``
        # writes ``targeted_track_id`` only on timepoints where
        # ``is_active_killing=True`` (it is the specific organoid whose
        # death signal increase triggered the classification; see
        # behav3d/features/advanced_timepoint_features.py line 437). So
        # each killing event is attributable to exactly one target.
        per_target_kill = None
        if killing_path.exists():
            df_tp = pd.read_csv(killing_path)
            required_cols = {
                "is_active_killing", "contact_event_id", "sample_name",
                "targeted_track_id",
            }
            if not df_tp.empty and required_cols.issubset(df_tp.columns):
                kill_mask = df_tp["is_active_killing"].astype(bool)
                df_tp_kill = df_tp.loc[
                    kill_mask
                    & df_tp["targeted_track_id"].notna(),
                    ["contact_event_id", "sample_name",
                     "targeted_track_id"],
                ].copy()
                if not df_tp_kill.empty:
                    df_tp_kill["targeted_track_id"] = (
                        df_tp_kill["targeted_track_id"]
                        .astype(float).astype("Int64").astype(str)
                    )
                    per_target_kill = (
                        df_tp_kill.groupby(
                            ["sample_name", "contact_event_id",
                             "targeted_track_id"],
                        ).size().rename("n_killing_tp").reset_index()
                    )

        # Explode target_track_ids -> one row per (event, target organoid)
        if "target_track_ids" not in df_ev.columns:
            continue
        df_ev["_targets"] = (
            df_ev["target_track_ids"].astype(str).str.split(",")
        )
        df_ex = df_ev.explode("_targets")
        df_ex["target_track_id"] = df_ex["_targets"].str.strip()
        df_ex = df_ex[
            df_ex["target_track_id"].ne("")
            & df_ex["target_track_id"].ne("nan")
        ]

        # Attach per-target killing (True only for the actual targeted
        # organoid, not for co-contacted organoids during the same event).
        df_ex["sample_name"] = df_ex["sample_name"].astype(str)
        df_ex["target_track_id"] = df_ex["target_track_id"].astype(str)
        if per_target_kill is not None and not per_target_kill.empty:
            df_ex = df_ex.merge(
                per_target_kill.rename(
                    columns={"targeted_track_id": "target_track_id"},
                ),
                on=["sample_name", "contact_event_id", "target_track_id"],
                how="left",
            )
        else:
            df_ex["n_killing_tp"] = np.nan
        df_ex["n_killing_tp"] = (
            df_ex["n_killing_tp"].fillna(0).astype(int)
        )
        df_ex["has_active_killing"] = df_ex["n_killing_tp"] > 0

        # Map to organoid type
        df_ex["organoid_type"] = [
            track_type_map.get(
                (str(sn), str(tid)), "unknown"
            )
            for sn, tid in zip(df_ex["sample_name"], df_ex["target_track_id"])
        ]
        df_ex = df_ex[df_ex["organoid_type"] != "unknown"]

        all_events.append(df_ex)

    if not all_events:
        return pd.DataFrame()

    df_all = pd.concat(all_events, ignore_index=True)

    # Attach fate of the target organoid from track_summary when available.
    # Default to "Live" for any organoid missing from track_summary (e.g. when
    # death data is not present) so downstream grouping still works.
    if (
        track_summary is not None
        and not track_summary.empty
        and "fate" in track_summary.columns
    ):
        ts = track_summary[["sample_name", "TrackID", "fate"]].copy()
        ts["sample_name"] = ts["sample_name"].astype(str)
        ts["TrackID"] = ts["TrackID"].astype(str)
        df_all["sample_name"] = df_all["sample_name"].astype(str)
        df_all["target_track_id"] = df_all["target_track_id"].astype(str)
        df_all = df_all.merge(
            ts.rename(columns={"TrackID": "target_track_id"}),
            on=["sample_name", "target_track_id"],
            how="left",
        )
        df_all["fate"] = df_all["fate"].fillna("Live")
    else:
        df_all["fate"] = "Live"

    keep = [
        "contact_event_id", "sample_name", "immune_track_id", "immune_type",
        "target_track_id", "organoid_type", "fate", "contact_duration",
        "has_active_killing", "n_killing_tp",
    ]
    return df_all[[c for c in keep if c in df_all.columns]]


def _style_dashboard_ax(ax: plt.Axes):
    """Apply consistent styling to a dashboard subplot."""
    ax.set_facecolor("#fafafa")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.6)
        ax.spines[spine].set_color("#888888")
    ax.tick_params(
        axis="both", which="both", length=4, width=0.6, colors="#555555",
    )
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.5, color="#cccccc")
    ax.set_axisbelow(True)


# -- Dashboard group-axis helpers (dynamic 2 x N_org_types layout) ----------

# Internal fate labels match the values written into ``track_summary["fate"]``
# by run_multi_organoid_interaction_comparison (``"Live"`` / ``"Dying"``).
# Order is Dead-then-Live so every group appears Dead -> Live -> Dead -> Live.
_FATE_ORDER = ("Dying", "Live")

# Display labels used on axes, legends and tick text. We show "Dead" rather
# than "Dying" since these organoids have crossed the death threshold by the
# end of the experiment.
_FATE_DISPLAY = {"Dying": "Dead", "Live": "Live"}

_GROUP_SEP = " | "
_FATE_COLORS = {"Live": "#66BB6A", "Dying": "#EF5350"}


def _build_org_fate_groups(
    org_types: list, fates: tuple = _FATE_ORDER,
) -> list:
    """Return the canonical ordered list of (organoid_type, fate) groups.

    Used by Panels C and D so axis positions are consistent regardless of
    which (org, fate) combinations actually contain data. Scales to
    ``2 * len(org_types)`` groups (N organoid types -> 2N slots). Fate
    order within each organoid type is Dead then Live.
    """
    return [(org, fate) for org in org_types for fate in fates]


def _group_label(org: str, fate: str) -> str:
    """Compose a unique categorical label for an (org, fate) group."""
    return f"{org}{_GROUP_SEP}{fate}"


def _apply_org_fate_xticks(
    ax: plt.Axes,
    group_order: list,
    *,
    show_separators: bool = True,
):
    """
    Set two-line x-tick labels (``organoid_type`` on top, display fate
    below) and draw faint vertical separators between successive organoid
    types. Fate is re-labeled via ``_FATE_DISPLAY`` (Dying -> Dead).
    """
    xticks = list(range(len(group_order)))
    xlabels = [
        f"{org}\n{_FATE_DISPLAY.get(fate, fate)}"
        for (org, fate) in group_order
    ]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=10)

    if show_separators and len(group_order) >= 2:
        for i in range(1, len(group_order)):
            prev_org = group_order[i - 1][0]
            this_org = group_order[i][0]
            if this_org != prev_org:
                ax.axvline(
                    i - 0.5, color="#bbbbbb", linewidth=0.7,
                    linestyle=(0, (2, 2)), alpha=0.8, zorder=0,
                )


def plot_interaction_overview_dashboard(
    track_summary: pd.DataFrame,
    df_contact_events: pd.DataFrame,
    immune_types: list,
    has_dead_data: bool,
):
    """
    Active-killing dashboard: two panels derived from the per-event / per-
    timepoint active-killing CSVs written by ``advanced_timepoint_features``
    (``contact_events_{im}.csv`` and ``active_killing_per_timepoint_{im}.csv``).

    Panels
    ------
    Left  -- Contact Event Duration per (organoid_type x fate); one point
             per immune-cell <-> organoid contact event, colored by whether
             that event actively killed the targeted organoid.
    Right -- Active Killing Efficiency per (organoid_type x fate); for each
             (organoid_type, fate) group, the % of its contact events that
             actively killed the targeted organoid, aggregated across all
             selected immune types.

    Returns ``None`` when no active-killing data is available -- the
    dashboard has nothing to show without it.
    """
    if df_contact_events.empty:
        return None

    org_types = sorted(df_contact_events["organoid_type"].unique())
    if not org_types:
        return None

    # Dynamic group count drives figure width: 2 x N_present_org_types.
    n_groups = max(2 * len(org_types), 4)
    figwidth = max(14.0, 3.8 * n_groups + 4.0)

    fig, (axC, axD) = plt.subplots(1, 2, figsize=(figwidth, 7), dpi=120)
    fig.patch.set_facecolor("#fafafa")

    _panel_contact_duration(axC, df_contact_events, org_types)
    _panel_killing_proportion(
        axD, df_contact_events, org_types, immune_types,
    )

    fig.suptitle(
        "Active Killing Dashboard",
        fontsize=15, fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# -- Dashboard sub-panels ---------------------------------------------------

def _panel_contact_duration(
    ax: plt.Axes,
    df_contact_events: pd.DataFrame,
    org_types: list,
):
    """
    Contact Event Duration panel.

    Each point = one immune-cell <-> organoid contact event (as written to
    ``contact_events_{im}.csv`` by the active-killing feature step). Points
    are grouped on the x-axis by (organoid_type, fate of the contacted
    organoid) and colored by whether that event actively killed the
    targeted organoid. x-axis scales dynamically to ``2 * len(org_types)``
    slots; fate order is Dead then Live per organoid type.
    """
    _style_dashboard_ax(ax)

    panel_title = "Contact event duration per organoid fate"

    if (
        df_contact_events.empty
        or "contact_duration" not in df_contact_events.columns
    ):
        ax.text(
            0.5, 0.5, "No contact event data",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color="#999999",
        )
        ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
        return

    df = df_contact_events.copy()
    if "fate" not in df.columns:
        df["fate"] = "Live"
    df = df[df["fate"].isin(_FATE_ORDER)].copy()
    if df.empty:
        ax.text(
            0.5, 0.5, "No matched organoid data",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color="#999999",
        )
        ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
        return

    df["killing_label"] = df["has_active_killing"].map(
        {True: "Active killing event", False: "Non-killing event"},
    )
    df["group"] = [
        _group_label(org, fate)
        for org, fate in zip(df["organoid_type"], df["fate"])
    ]

    group_pairs = _build_org_fate_groups(org_types)
    group_order = [_group_label(o, f) for (o, f) in group_pairs]
    group_to_idx = {g: i for i, g in enumerate(group_order)}

    # Per-group violin background color follows fate (soft-red for Dead,
    # slate for Live) so the fate axis is readable at a glance.
    fate_violin_color = {"Live": "#90A4AE", "Dying": "#E57373"}
    fate_violin_dark = {"Live": "#546E7A", "Dying": "#C62828"}
    palette = {
        _group_label(o, f): fate_violin_color[f]
        for (o, f) in group_pairs
    }
    dark_palette = {
        _group_label(o, f): fate_violin_dark[f]
        for (o, f) in group_pairs
    }

    common = dict(
        data=df, x="group", y="contact_duration",
        order=group_order, ax=ax, width=0.75,
    )

    sns.violinplot(
        **common, hue="group", palette=palette, legend=False,
        inner=None, linewidth=0.8, alpha=0.30,
        density_norm="count", cut=0,
    )
    sns.boxplot(
        **common, hue="group", palette=dark_palette, legend=False,
        fliersize=0,
        boxprops=dict(alpha=0.6, linewidth=0.8),
        whiskerprops=dict(alpha=0.7, linewidth=0.8),
        capprops=dict(alpha=0.7, linewidth=0.8),
        showcaps=True,
        medianprops=dict(color="white", linewidth=2),
    )

    # Manual scatter: color encodes whether this contact event actively
    # killed its targeted organoid (per-target attribution from
    # active_killing_per_timepoint_{im}.csv::targeted_track_id).
    kill_color = {
        "Active killing event": "#C62828",
        "Non-killing event": "#1565C0",
    }
    rng = np.random.default_rng(1)
    for g_label, sub in df.groupby("group"):
        if g_label not in group_to_idx or len(sub) == 0:
            continue
        xi = group_to_idx[g_label]
        jitter = rng.uniform(-0.18, 0.18, size=len(sub))
        colors = [kill_color[lab] for lab in sub["killing_label"]]
        ax.scatter(
            xi + jitter, sub["contact_duration"],
            s=10, c=colors, alpha=0.55,
            edgecolors="white", linewidths=0.25, zorder=3,
        )

    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=kill_color["Active killing event"],
            markeredgecolor="white", markeredgewidth=0.5,
            label="Active killing event",
        ),
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=kill_color["Non-killing event"],
            markeredgecolor="white", markeredgewidth=0.5,
            label="Non-killing event",
        ),
    ]
    ax.legend(
        handles=handles, fontsize=9,
        frameon=True, framealpha=0.85, edgecolor="#cccccc",
        title="Each point = one contact event",
        title_fontsize=9, loc="upper right",
    )

    _apply_org_fate_xticks(ax, group_pairs)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Organoid type / fate", fontsize=11, labelpad=8)
    ax.set_ylabel(
        "Contact duration (timepoints)",
        fontsize=11, labelpad=8,
    )
    ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)


def _panel_killing_proportion(
    ax: plt.Axes,
    df_contact_events: pd.DataFrame,
    org_types: list,
    immune_types: list,
):
    """
    Active Killing Efficiency panel.

    One bar per (organoid_type, fate): the percentage of that group's
    contact events that actively killed the targeted organoid. Numerator
    and denominator count contact events (not organoids, not immune cells);
    killing attribution uses the per-target ``targeted_track_id`` column of
    ``active_killing_per_timepoint_{im}.csv``. Aggregates across all
    selected immune types (per-immune breakdown already lives in
    ``active_killing_summary_{im}.csv``).
    """
    _style_dashboard_ax(ax)

    panel_title = "Active killing efficiency per organoid fate"

    if df_contact_events.empty:
        ax.text(
            0.5, 0.5, "No active killing data",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color="#999999",
        )
        ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
        return

    df = df_contact_events.copy()
    if "fate" not in df.columns:
        df["fate"] = "Live"
    df = df[df["fate"].isin(_FATE_ORDER)].copy()
    if df.empty:
        ax.text(
            0.5, 0.5, "No matched organoid data",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color="#999999",
        )
        ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
        return

    # Aggregate across all selected immune types. ``has_active_killing`` is
    # per (event, target) thanks to the per-target attribution done in
    # _load_active_killing_data, so summing gives the number of contact
    # events that killed this specific targeted organoid without double
    # counting co-contacted organoids.
    summary = (
        df.groupby(["organoid_type", "fate"])
          .agg(
              n_events=("contact_event_id", "count"),
              n_killing=("has_active_killing", "sum"),
          )
          .reset_index()
    )
    summary["killing_pct"] = np.where(
        summary["n_events"] > 0,
        summary["n_killing"] / summary["n_events"] * 100,
        0.0,
    )

    group_pairs = _build_org_fate_groups(org_types)

    x = np.arange(len(group_pairs))
    vals = []
    n_events_list = []
    n_kill_list = []
    colors = []
    for (org, fate) in group_pairs:
        row = summary[
            (summary["organoid_type"] == org) & (summary["fate"] == fate)
        ]
        if len(row):
            vals.append(float(row["killing_pct"].iloc[0]))
            n_events_list.append(int(row["n_events"].iloc[0]))
            n_kill_list.append(int(row["n_killing"].iloc[0]))
        else:
            vals.append(0.0)
            n_events_list.append(0)
            n_kill_list.append(0)
        colors.append(_FATE_COLORS[fate])

    bar_w = 0.6
    ax.bar(
        x, vals, bar_w, color=colors, edgecolor="white",
        linewidth=0.8, alpha=0.88,
    )

    # Annotate each bar with killing_events / total_events.
    for xi, v, ne, nk in zip(x, vals, n_events_list, n_kill_list):
        txt = f"{nk}/{ne} events" if ne > 0 else "0/0 events"
        ax.text(
            xi, v + 1.2, txt,
            ha="center", va="bottom", fontsize=8.5,
            color="#444444",
        )

    # Legend for fate colors (Dead first, Live second).
    handles = [
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=_FATE_COLORS["Dying"],
            edgecolor="white", label="Dead",
        ),
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=_FATE_COLORS["Live"],
            edgecolor="white", label="Live",
        ),
    ]
    ax.legend(
        handles=handles, fontsize=9,
        frameon=True, framealpha=0.85, edgecolor="#cccccc",
        title="Organoid fate", title_fontsize=9, loc="upper right",
    )

    _apply_org_fate_xticks(ax, group_pairs)
    top = max(max(vals) * 1.25 if vals else 10, 10)
    ax.set_ylim(bottom=0, top=min(top, 105))
    ax.set_xlabel("Organoid type / fate", fontsize=11, labelpad=8)
    ax.set_ylabel(
        "Contact events that actively killed (%)",
        fontsize=11, labelpad=8,
    )
    ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
