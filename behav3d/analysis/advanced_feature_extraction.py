"""
Advanced Feature Extraction for BEHAV3D

This module provides advanced analysis features that go beyond basic track feature extraction.
It focuses on interaction-based analyses such as active killing detection.

-------------------------------------
--------------- FEATURES ------------
-------------------------------------

### Active Killing Analysis
Detects when an immune cell (e.g., T cell) makes contact with a target cell (e.g., organoid)
and subsequently causes an increase in death signal (dead dye/dead mask) above background levels.

Key concepts:
- Contact event: When immune cell touches target cell
- Observation window: N timepoints after contact to observe death signal changes
- Background death rate: Average death signal increase across all timepoints (regardless of contact)
- Active killing: Death signal increase after contact that exceeds background rate

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Per-contact interaction features
- contact_event_id: Unique identifier for each contact event
- immune_track_id: TrackID of the immune cell
- target_track_id: TrackID of the target cell  
- contact_start_t: Timepoint when contact began
- death_signal_before: Death signal at contact start
- death_signal_after: Death signal after observation window
- death_signal_increase: Change in death signal
- is_active_killing: Whether this exceeds background threshold

# Summary statistics
- total_contact_events: Number of contact events
- active_killing_events: Number classified as active killing
- active_killing_rate: Proportion of contacts resulting in killing
- background_death_rate: Average death increase per timepoint
"""

import numpy as np
import pandas as pd
from pathlib import Path
import time
from typing import Optional, List, Tuple, Dict, Union

from behav3d.utils import get_current_time, format_time


def calculate_background_death_rate(
    df_target_tracks: pd.DataFrame,
    death_signal_column: str = "mean_dead_dye",
    groupby: List[str] = ["sample_name", "TrackID"]
) -> Tuple[float, pd.DataFrame]:
    """
    Calculate the average background death rate across all target cells.
    
    This represents the baseline increase in death signal per timepoint,
    regardless of immune cell contact. Active killing should exceed this rate.
    
    Parameters
    ----------
    df_target_tracks : pd.DataFrame
        DataFrame containing target cell tracks with death signal information
    death_signal_column : str
        Column name containing the death signal (e.g., 'mean_dead_dye', 'percentage_dead_mask')
    groupby : List[str]
        Columns to group by for per-track calculations
        
    Returns
    -------
    background_rate : float
        Average death signal increase per timepoint across all tracks
    df_death_rates : pd.DataFrame
        Per-track death rate statistics
    """
    df = df_target_tracks.copy()
    df = df.sort_values(by=groupby + ["position_t"])
    
    # Calculate per-timepoint change in death signal for each track
    df["death_signal_change"] = df.groupby(groupby)[death_signal_column].diff()
    
    # Calculate statistics per track
    df_death_rates = df.groupby(groupby).agg(
        mean_death_increase=(death_signal_column, lambda x: x.iloc[-1] - x.iloc[0] if len(x) > 1 else 0),
        track_duration=("position_t", lambda x: x.max() - x.min() + 1),
        mean_death_rate_per_t=("death_signal_change", "mean"),
        total_death_signal_change=("death_signal_change", "sum"),
        initial_death_signal=(death_signal_column, "first"),
        final_death_signal=(death_signal_column, "last"),
    ).reset_index()
    
    # Calculate rate per timepoint
    df_death_rates["death_rate_per_timepoint"] = (
        df_death_rates["mean_death_increase"] / df_death_rates["track_duration"].clip(lower=1)
    )
    
    # Overall background rate (mean across all tracks)
    background_rate = df_death_rates["death_rate_per_timepoint"].mean()
    
    return background_rate, df_death_rates


def identify_contact_events(
    df_immune_tracks: pd.DataFrame,
    contact_column: str = "organoid_contact",
    touching_column: str = "touching_organoids",
    min_contact_duration: int = 1
) -> pd.DataFrame:
    """
    Identify distinct contact events between immune cells and target cells.
    
    A contact event is defined as a continuous period where an immune cell
    is in contact with a specific target cell.
    
    Parameters
    ----------
    df_immune_tracks : pd.DataFrame
        DataFrame containing immune cell tracks with contact information
    contact_column : str
        Boolean column indicating contact (e.g., 'organoid_contact')
    touching_column : str
        Column containing IDs of touched cells (e.g., 'touching_organoids')
    min_contact_duration : int
        Minimum number of timepoints for a valid contact event
        
    Returns
    -------
    df_contact_events : pd.DataFrame
        DataFrame with one row per contact event
    """
    df = df_immune_tracks.copy()
    df = df.sort_values(by=["sample_name", "TrackID", "position_t"])
    
    contact_events = []
    event_id = 0
    
    for (sample_name, track_id), group in df.groupby(["sample_name", "TrackID"]):
        group = group.reset_index(drop=True)
        
        # Skip if no contact column or touching column
        if contact_column not in group.columns or touching_column not in group.columns:
            continue
            
        in_contact = False
        current_targets = set()
        contact_start_t = None
        contact_start_idx = None
        
        for idx, row in group.iterrows():
            is_contacting = row.get(contact_column, False)
            touching_str = row.get(touching_column, "")
            
            # Parse touching targets
            if pd.notna(touching_str) and str(touching_str).strip():
                current_touching = set(str(touching_str).split(","))
            else:
                current_touching = set()
            
            if is_contacting and current_touching:
                if not in_contact:
                    # Start of new contact event
                    in_contact = True
                    contact_start_t = row["position_t"]
                    contact_start_idx = idx
                    current_targets = current_touching
                else:
                    # Continue contact - add any new targets
                    current_targets.update(current_touching)
            else:
                if in_contact:
                    # End of contact event
                    contact_end_t = group.loc[idx - 1, "position_t"] if idx > 0 else contact_start_t
                    contact_duration = int(contact_end_t - contact_start_t + 1)
                    
                    if contact_duration >= min_contact_duration:
                        for target_id in current_targets:
                            event_id += 1
                            contact_events.append({
                                "contact_event_id": event_id,
                                "sample_name": sample_name,
                                "immune_track_id": track_id,
                                "target_track_id": target_id.strip() if isinstance(target_id, str) else target_id,
                                "contact_start_t": contact_start_t,
                                "contact_end_t": contact_end_t,
                                "contact_duration": contact_duration,
                            })
                    
                    in_contact = False
                    current_targets = set()
                    contact_start_t = None
        
        # Handle contact that extends to end of track
        if in_contact and current_targets:
            contact_end_t = group.iloc[-1]["position_t"]
            contact_duration = int(contact_end_t - contact_start_t + 1)
            
            if contact_duration >= min_contact_duration:
                for target_id in current_targets:
                    event_id += 1
                    contact_events.append({
                        "contact_event_id": event_id,
                        "sample_name": sample_name,
                        "immune_track_id": track_id,
                        "target_track_id": target_id.strip() if isinstance(target_id, str) else target_id,
                        "contact_start_t": contact_start_t,
                        "contact_end_t": contact_end_t,
                        "contact_duration": contact_duration,
                    })
    
    df_contact_events = pd.DataFrame(contact_events)
    return df_contact_events


def analyze_active_killing(
    df_immune_tracks: pd.DataFrame,
    df_target_tracks: pd.DataFrame,
    observation_window: int = 5,
    death_signal_column: str = "mean_dead_dye",
    contact_column: str = "organoid_contact",
    touching_column: str = "touching_organoids",
    min_contact_duration: int = 1,
    killing_threshold_multiplier: float = 1.5,
    absolute_killing_threshold: Optional[float] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Analyze active killing events between immune cells and target cells.
    
    Active killing is defined as: immune cell contacts target cell, and within
    the observation window, the target's death signal increases above the
    background death rate.
    
    Parameters
    ----------
    df_immune_tracks : pd.DataFrame
        DataFrame containing immune cell tracks (e.g., T cells)
    df_target_tracks : pd.DataFrame
        DataFrame containing target cell tracks (e.g., organoids)
    observation_window : int
        Number of timepoints after contact start to observe death signal changes
    death_signal_column : str
        Column in target tracks containing death signal
    contact_column : str
        Boolean column in immune tracks indicating contact
    touching_column : str
        Column in immune tracks containing IDs of touched targets
    min_contact_duration : int
        Minimum timepoints for valid contact event
    killing_threshold_multiplier : float
        Multiplier for background rate to determine killing threshold.
        A contact is "active killing" if death increase > background_rate * multiplier * window
    absolute_killing_threshold : float, optional
        If provided, use this absolute threshold instead of relative to background
        
    Returns
    -------
    df_killing_events : pd.DataFrame
        DataFrame with killing analysis for each contact event
    df_summary : pd.DataFrame
        Summary statistics per sample
    stats : dict
        Overall statistics dictionary
    """
    print(f"{get_current_time()} - Calculating background death rate...")
    
    # Calculate background death rate from target cells
    background_rate, df_death_rates = calculate_background_death_rate(
        df_target_tracks,
        death_signal_column=death_signal_column
    )
    
    print(f"{get_current_time()} - Background death rate: {background_rate:.4f} per timepoint")
    
    # Identify contact events
    print(f"{get_current_time()} - Identifying contact events...")
    df_contact_events = identify_contact_events(
        df_immune_tracks,
        contact_column=contact_column,
        touching_column=touching_column,
        min_contact_duration=min_contact_duration
    )
    
    if df_contact_events.empty:
        print(f"{get_current_time()} - No contact events found")
        return pd.DataFrame(), pd.DataFrame(), {"total_contacts": 0}
    
    print(f"{get_current_time()} - Found {len(df_contact_events)} contact events")
    
    # Analyze death signal changes for each contact event
    print(f"{get_current_time()} - Analyzing death signal changes...")
    
    killing_analysis = []
    
    for _, event in df_contact_events.iterrows():
        sample_name = event["sample_name"]
        target_id = event["target_track_id"]
        contact_start = event["contact_start_t"]
        
        # Try to convert target_id to match the format in target tracks
        try:
            target_id_int = int(float(target_id))
        except (ValueError, TypeError):
            target_id_int = target_id
        
        # Get target track data
        target_mask = (
            (df_target_tracks["sample_name"] == sample_name) & 
            (df_target_tracks["TrackID"].astype(str) == str(target_id_int))
        )
        df_target = df_target_tracks[target_mask].sort_values("position_t")
        
        if df_target.empty:
            # Target track not found - skip this event
            continue
        
        # Get death signal at contact start
        at_contact = df_target[df_target["position_t"] == contact_start]
        if at_contact.empty:
            # Find closest timepoint before or at contact
            before_contact = df_target[df_target["position_t"] <= contact_start]
            if before_contact.empty:
                continue
            at_contact = before_contact.iloc[[-1]]
        
        death_signal_at_contact = at_contact[death_signal_column].values[0]
        
        # Get death signal after observation window
        observation_end = contact_start + observation_window
        after_window = df_target[df_target["position_t"] >= observation_end]
        
        if after_window.empty:
            # Use last available timepoint
            after_window = df_target.iloc[[-1]]
        else:
            after_window = after_window.iloc[[0]]
        
        death_signal_after = after_window[death_signal_column].values[0]
        actual_window = after_window["position_t"].values[0] - contact_start
        
        # Calculate death signal increase
        death_signal_increase = death_signal_after - death_signal_at_contact
        
        # Determine killing threshold
        if absolute_killing_threshold is not None:
            killing_threshold = absolute_killing_threshold
        else:
            # Expected increase based on background rate over the observation window
            expected_increase = background_rate * actual_window
            killing_threshold = expected_increase * killing_threshold_multiplier
        
        is_active_killing = death_signal_increase > killing_threshold
        
        killing_analysis.append({
            **event.to_dict(),
            "target_track_id_matched": target_id_int,
            "death_signal_at_contact": death_signal_at_contact,
            "death_signal_after_window": death_signal_after,
            "death_signal_increase": death_signal_increase,
            "observation_window_actual": actual_window,
            "expected_background_increase": background_rate * actual_window,
            "killing_threshold": killing_threshold,
            "is_active_killing": is_active_killing,
            "killing_efficiency": death_signal_increase / (background_rate * actual_window + 1e-10)
        })
    
    df_killing_events = pd.DataFrame(killing_analysis)
    
    if df_killing_events.empty:
        print(f"{get_current_time()} - No matched contact events found")
        return pd.DataFrame(), pd.DataFrame(), {"total_contacts": 0}
    
    # Calculate summary statistics
    print(f"{get_current_time()} - Calculating summary statistics...")
    
    df_summary = df_killing_events.groupby("sample_name").agg(
        total_contact_events=("contact_event_id", "count"),
        active_killing_events=("is_active_killing", "sum"),
        mean_death_signal_increase=("death_signal_increase", "mean"),
        median_death_signal_increase=("death_signal_increase", "median"),
        mean_contact_duration=("contact_duration", "mean"),
        mean_killing_efficiency=("killing_efficiency", "mean"),
    ).reset_index()
    
    df_summary["active_killing_rate"] = (
        df_summary["active_killing_events"] / df_summary["total_contact_events"]
    )
    
    # Overall statistics
    stats = {
        "total_contacts": len(df_killing_events),
        "total_active_killing": df_killing_events["is_active_killing"].sum(),
        "overall_killing_rate": df_killing_events["is_active_killing"].mean(),
        "background_death_rate": background_rate,
        "mean_death_signal_increase": df_killing_events["death_signal_increase"].mean(),
        "observation_window": observation_window,
        "killing_threshold_multiplier": killing_threshold_multiplier,
    }
    
    print(f"{get_current_time()} - Active killing analysis complete:")
    print(f"    Total contact events: {stats['total_contacts']}")
    print(f"    Active killing events: {stats['total_active_killing']}")
    print(f"    Overall killing rate: {stats['overall_killing_rate']:.2%}")
    
    return df_killing_events, df_summary, stats


def run_active_killing_analysis(
    metadata: pd.DataFrame,
    output_dir: Union[str, Path],
    immune_cell_type: str = "tcell",
    target_cell_type: str = "organoid",
    observation_window: int = 5,
    death_signal_column: str = "mean_dead_dye",
    contact_column: Optional[str] = None,
    touching_column: Optional[str] = None,
    min_contact_duration: int = 1,
    killing_threshold_multiplier: float = 1.5,
    absolute_killing_threshold: Optional[float] = None,
    save_results: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Run active killing analysis on BEHAV3D processed data.
    
    This is the main entry point for active killing analysis. It loads the
    appropriate track feature files and runs the analysis.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        BEHAV3D metadata DataFrame
    output_dir : str or Path
        BEHAV3D output directory containing analysis results
    immune_cell_type : str
        Type of immune cell (e.g., 'tcell', 'macro', 'nk')
    target_cell_type : str
        Type of target cell (e.g., 'organoid', 'tumorcell')
    observation_window : int
        Timepoints after contact to observe death signal
    death_signal_column : str
        Column containing death signal in target tracks
    contact_column : str, optional
        Column indicating contact. Defaults to '{target_cell_type}_contact'
    touching_column : str, optional
        Column with touched cell IDs. Defaults to 'touching_{target_cell_type}s'
    min_contact_duration : int
        Minimum contact duration in timepoints
    killing_threshold_multiplier : float
        Multiplier for background rate threshold
    absolute_killing_threshold : float, optional
        Absolute threshold (overrides multiplier if set)
    save_results : bool
        Whether to save results to CSV files
        
    Returns
    -------
    df_killing_events : pd.DataFrame
        Per-contact killing analysis
    df_summary : pd.DataFrame
        Per-sample summary statistics
    stats : dict
        Overall statistics
    """
    print(f"--------------- Running Active Killing Analysis ---------------")
    print(f"Immune cell type: {immune_cell_type}")
    print(f"Target cell type: {target_cell_type}")
    print(f"Observation window: {observation_window} timepoints")
    start_time = time.time()
    
    output_dir = Path(output_dir)
    
    # Set default column names based on target cell type
    if contact_column is None:
        contact_column = f"{target_cell_type}_contact"
    if touching_column is None:
        touching_column = f"touching_{target_cell_type}s"
    
    # Load immune cell tracks
    immune_feature_dir = output_dir / "analysis" / immune_cell_type / "track_features"
    immune_tracks_path = immune_feature_dir / f"BEHAV3D_{immune_cell_type}_combined_track_features_filtered.csv"
    
    if not immune_tracks_path.exists():
        # Try non-filtered version
        immune_tracks_path = immune_feature_dir / f"BEHAV3D_{immune_cell_type}_combined_track_features.csv"
    
    if not immune_tracks_path.exists():
        raise FileNotFoundError(f"Could not find immune cell tracks at {immune_tracks_path}")
    
    print(f"{get_current_time()} - Loading immune cell tracks from {immune_tracks_path}")
    df_immune_tracks = pd.read_csv(immune_tracks_path)
    
    # Load target cell tracks
    target_feature_dir = output_dir / "analysis" / target_cell_type / "track_features"
    target_tracks_path = target_feature_dir / f"BEHAV3D_{target_cell_type}_combined_track_features_filtered.csv"
    
    if not target_tracks_path.exists():
        target_tracks_path = target_feature_dir / f"BEHAV3D_{target_cell_type}_combined_track_features.csv"
    
    if not target_tracks_path.exists():
        raise FileNotFoundError(f"Could not find target cell tracks at {target_tracks_path}")
    
    print(f"{get_current_time()} - Loading target cell tracks from {target_tracks_path}")
    df_target_tracks = pd.read_csv(target_tracks_path)
    
    # Run the analysis
    df_killing_events, df_summary, stats = analyze_active_killing(
        df_immune_tracks=df_immune_tracks,
        df_target_tracks=df_target_tracks,
        observation_window=observation_window,
        death_signal_column=death_signal_column,
        contact_column=contact_column,
        touching_column=touching_column,
        min_contact_duration=min_contact_duration,
        killing_threshold_multiplier=killing_threshold_multiplier,
        absolute_killing_threshold=absolute_killing_threshold,
    )
    
    # Save results
    if save_results and not df_killing_events.empty:
        results_dir = output_dir / "analysis" / "active_killing"
        results_dir.mkdir(parents=True, exist_ok=True)
        
        killing_events_path = results_dir / f"active_killing_events_{immune_cell_type}_vs_{target_cell_type}.csv"
        summary_path = results_dir / f"active_killing_summary_{immune_cell_type}_vs_{target_cell_type}.csv"
        
        df_killing_events.to_csv(killing_events_path, index=False)
        df_summary.to_csv(summary_path, index=False)
        
        print(f"{get_current_time()} - Results saved to {results_dir}")
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    
    return df_killing_events, df_summary, stats


def calculate_killing_dynamics_over_time(
    df_killing_events: pd.DataFrame,
    df_immune_tracks: pd.DataFrame,
    time_bins: Optional[List[int]] = None,
    time_column: str = "contact_start_t"
) -> pd.DataFrame:
    """
    Calculate how active killing rate changes over the course of the experiment.
    
    Parameters
    ----------
    df_killing_events : pd.DataFrame
        Output from analyze_active_killing
    df_immune_tracks : pd.DataFrame
        Immune cell tracks for time reference
    time_bins : List[int], optional
        Custom time bin edges. If None, uses quartiles of experiment duration.
    time_column : str
        Column to use for time binning
        
    Returns
    -------
    df_dynamics : pd.DataFrame
        Killing statistics per time bin
    """
    if df_killing_events.empty:
        return pd.DataFrame()
    
    if time_bins is None:
        max_t = df_immune_tracks["position_t"].max()
        time_bins = [0, max_t // 4, max_t // 2, 3 * max_t // 4, max_t + 1]
    
    df = df_killing_events.copy()
    df["time_bin"] = pd.cut(
        df[time_column], 
        bins=time_bins, 
        labels=[f"{time_bins[i]}-{time_bins[i+1]}" for i in range(len(time_bins)-1)]
    )
    
    df_dynamics = df.groupby(["sample_name", "time_bin"]).agg(
        n_contacts=("contact_event_id", "count"),
        n_active_killing=("is_active_killing", "sum"),
        mean_death_increase=("death_signal_increase", "mean"),
    ).reset_index()
    
    df_dynamics["killing_rate"] = df_dynamics["n_active_killing"] / df_dynamics["n_contacts"]
    
    return df_dynamics
