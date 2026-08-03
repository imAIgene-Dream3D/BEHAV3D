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
- Contact event: Continuous period where immune cell touches ANY target cell
- Total contact duration: Full length of continuous contact (must exceed min_contact_duration)
- Per-organoid, per-timepoint calculation: Each touched organoid is evaluated independently
  at every timepoint of the contact (sliding window, not just once at contact start)
- Observation window: N timepoints after EACH contact timepoint to measure death signal change
- Active killing: Death signal increase (t -> t + window) exceeds either a multiplier of the
  organoid's own signal at t, or a fixed absolute increase, depending on the selected mode

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Per-timepoint features (added to track features CSV)
- is_active_killing: Boolean, True if this timepoint's contact caused killing
- death_signal_increase_{N}tp: Death signal change over observation window
- killing_efficiency: Ratio of actual death signal increase vs the organoid's own threshold

# Where no contact occurs:
- is_active_killing = False
- death_signal_increase = 0.0
- killing_efficiency = 0.0

# Summary statistics per sample
- total_contact_events: Number of qualifying contact events (â‰¥ min_duration)
- active_killing_events: Number of timepoints classified as active killing
- active_killing_rate: Proportion of contact timepoints with active killing
"""

import pandas as pd
from pathlib import Path
import time
from typing import Optional, List, Tuple, Dict, Union
import numpy as np

from behav3d.core.metadata import detect_organoid_types_from_metadata
from behav3d.core.utils import get_current_time, format_time
from behav3d.io.images import load_image


# If an organoid's death signal is exactly 0 at a given window's start timepoint, a
# multiplier threshold would be trivially 0 (any nonzero signal would count as active
# killing). Substitute this small epsilon as the effective starting signal in that case only.
ZERO_BASELINE_EPSILON = 0.1


def identify_contact_events_global(
    df_immune_tracks: pd.DataFrame,
    target_cell_types: List[str],
    min_contact_duration: int = 1
) -> pd.DataFrame:
    """
    Identify distinct contact events between immune cells and ANY target cell.
    
    A contact event is defined as a continuous period where an immune cell
    is in contact with ANY target cell type. The total duration is measured
    across the entire continuous contact period.
    
    Parameters
    ----------
    df_immune_tracks : pd.DataFrame
        DataFrame containing immune cell tracks with contact information
    target_cell_types : List[str]
        List of target cell types to check for contact (e.g., ['organoid', 'organoid2'])
    min_contact_duration : int
        Minimum number of timepoints for a valid contact event
        
    Returns
    -------
    df_contact_events : pd.DataFrame
        DataFrame with one row per contact event (contact with any target)
    """
    df = df_immune_tracks.copy()
    df = df.sort_values(by=["sample_name", "TrackID", "position_t"])
    
    # Create a global contact column (True if contacting ANY target type)
    contact_columns = [f"{target_type}_contact" for target_type in target_cell_types]
    existing_contact_cols = [c for c in contact_columns if c in df.columns]
    
    if not existing_contact_cols:
        print(f"  Warning: No contact columns found for target types: {target_cell_types}")
        return pd.DataFrame()
    
    # Combine all contact columns
    df["_any_contact"] = df[existing_contact_cols].any(axis=1)
    
    # Combine all touching columns
    touching_columns = [f"touching_{target_type}s" for target_type in target_cell_types]
    existing_touching_cols = [c for c in touching_columns if c in df.columns]
    
    def combine_touching(row):
        """Combine all touching targets into one string"""
        all_targets = []
        for col in existing_touching_cols:
            val = row.get(col, "")
            if pd.notna(val) and str(val).strip():
                all_targets.extend([t.strip() for t in str(val).split(",") if t.strip()])
        return ",".join(all_targets) if all_targets else ""
    
    df["_all_touching"] = df.apply(combine_touching, axis=1)
    
    contact_events = []
    event_id = 0
    
    for (sample_name, track_id), group in df.groupby(["sample_name", "TrackID"]):
        group = group.reset_index(drop=True)
        
        in_contact = False
        contact_start_t = None
        contact_timepoints = []  # List of timepoints in this contact
        all_targets = set()
        
        for idx, row in group.iterrows():
            is_contacting = row["_any_contact"]
            touching_str = row["_all_touching"]
            
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
                    contact_timepoints = [row["position_t"]]
                    all_targets = current_touching
                else:
                    # Continue contact
                    contact_timepoints.append(row["position_t"])
                    all_targets.update(current_touching)
            else:
                if in_contact:
                    # End of contact event
                    contact_end_t = contact_timepoints[-1]
                    contact_duration = len(contact_timepoints)
                    
                    if contact_duration >= min_contact_duration:
                        event_id += 1
                        contact_events.append({
                            "contact_event_id": event_id,
                            "sample_name": sample_name,
                            "immune_track_id": track_id,
                            "target_track_ids": ",".join(sorted([str(t) for t in all_targets])),
                            "contact_start_t": contact_start_t,
                            "contact_end_t": contact_end_t,
                            "contact_duration": contact_duration,
                            "contact_timepoints": contact_timepoints.copy(),
                        })
                    
                    in_contact = False
                    contact_start_t = None
                    contact_timepoints = []
                    all_targets = set()
        
        # Handle contact that extends to end of track
        if in_contact and contact_timepoints:
            contact_end_t = contact_timepoints[-1]
            contact_duration = len(contact_timepoints)
            
            if contact_duration >= min_contact_duration:
                event_id += 1
                contact_events.append({
                    "contact_event_id": event_id,
                    "sample_name": sample_name,
                    "immune_track_id": track_id,
                    "target_track_ids": ",".join(sorted([str(t) for t in all_targets])),
                    "contact_start_t": contact_start_t,
                    "contact_end_t": contact_end_t,
                    "contact_duration": contact_duration,
                    "contact_timepoints": contact_timepoints.copy(),
                })
    
    df_contact_events = pd.DataFrame(contact_events)
    return df_contact_events


def analyze_active_killing_per_timepoint(
    df_immune_tracks: pd.DataFrame,
    df_target_tracks: pd.DataFrame,
    df_contact_events: pd.DataFrame,
    observation_window: int = 5,
    death_signal_column: str = "mean_dead_dye",
    killing_threshold_multiplier: float = 1.5,
    absolute_killing_threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Calculate active killing status for each contact timepoint using a
    sliding observation window anchored at that timepoint.

    Each touched target organoid is evaluated independently at every
    timepoint of the contact event (not just once at contact start):
    1. Get the target's death signal at timepoint t
    2. Get the target's death signal at t + observation_window
       (or the last available timepoint on the target's track, if the
       track ends before then -- this can reach past the contact event's
       own end)
    3. Compare the increase to a threshold:
       - multiplier mode: threshold = death_at_t x (killing_threshold_multiplier - 1)
         (i.e. active if death_at_t_plus_window >= death_at_t x killing_threshold_multiplier).
         If death_at_t is exactly 0, ZERO_BASELINE_EPSILON is substituted so the
         threshold isn't trivially 0.
       - absolute mode: threshold = absolute_killing_threshold (flat, in death_signal_column units)
    4. Among all touched targets, the one with the highest efficiency (increase relative to
       its own threshold) AT THAT TIMEPOINT is reported as the timepoint's representative
       target/values -- so the representative target can change from one contact timepoint
       to the next. This is reported whether or not it is active; targeted_track_id is only
       set when active.
    5. Each contact timepoint gets its own independently computed row: a long contact can
       show killing at some timepoints and not others instead of one value broadcast
       across the whole event.

    Parameters
    ----------
    df_immune_tracks : pd.DataFrame
        DataFrame containing immune cell tracks
    df_target_tracks : pd.DataFrame
        DataFrame containing target cell tracks with death signal
    df_contact_events : pd.DataFrame
        Contact events from identify_contact_events_global
    observation_window : int
        Number of timepoints after each contact timepoint to measure death signal change
    death_signal_column : str
        Column in target tracks containing death signal
    killing_threshold_multiplier : float
        Multiplier applied to a target's own signal at each timepoint (used when
        absolute_killing_threshold is None)
    absolute_killing_threshold : float, optional
        If provided, use this flat value as the killing threshold instead of the
        multiplier-based threshold. Useful when organoids start with near-zero signal.

    Returns
    -------
    df_killing_per_timepoint : pd.DataFrame
        DataFrame with one row per (contact_event, timepoint) with killing status
    """
    if df_contact_events.empty:
        return pd.DataFrame()

    # Get max timepoint per sample (to know when the observation window runs past track end)
    max_timepoints = df_target_tracks.groupby("sample_name")["position_t"].max().to_dict()

    killing_results = []

    for _, event in df_contact_events.iterrows():
        sample_name = event["sample_name"]
        immune_track_id = event["immune_track_id"]
        target_ids = [t.strip() for t in str(event["target_track_ids"]).split(",")]
        contact_timepoints = event["contact_timepoints"]

        max_t = max_timepoints.get(sample_name, contact_timepoints[-1])

        # Get target tracks for this sample
        target_mask = df_target_tracks["sample_name"] == sample_name
        df_targets_sample = df_target_tracks[target_mask]

        # Precompute a sorted (position_t, death_signal) lookup per touched
        # target once per event, reused across every contact timepoint below
        # (avoids re-filtering/sorting the target dataframe per timepoint).
        target_lookups = {}
        for target_id in target_ids:
            try:
                target_id_int = int(float(target_id))
            except (ValueError, TypeError):
                target_id_int = target_id

            target_rows = df_targets_sample[
                df_targets_sample["TrackID"].astype(str) == str(target_id_int)
            ].sort_values("position_t")

            if target_rows.empty:
                continue

            target_lookups[target_id_int] = (
                target_rows["position_t"].to_numpy(),
                target_rows[death_signal_column].to_numpy(),
            )

        for t in contact_timepoints:
            observation_end_t = t + observation_window
            can_observe = observation_end_t <= max_t

            target_evals = []

            for target_id_int, (pt_arr, death_arr) in target_lookups.items():
                # Death signal at or before t (mirrors the old exact-match,
                # else-nearest-before lookup for contact_start_t).
                start_idx = np.searchsorted(pt_arr, t, side="right") - 1
                if start_idx < 0:
                    continue
                death_at_start = death_arr[start_idx]

                # Death signal at the first frame >= observation_end_t,
                # falling back to the target's last available frame.
                end_idx = np.searchsorted(pt_arr, observation_end_t, side="left")
                if end_idx >= len(pt_arr):
                    end_idx = len(pt_arr) - 1
                death_at_end = death_arr[end_idx]

                death_increase = death_at_end - death_at_start

                if absolute_killing_threshold is not None:
                    threshold_increase = absolute_killing_threshold
                else:
                    effective_start = death_at_start if death_at_start != 0 else ZERO_BASELINE_EPSILON
                    threshold_increase = effective_start * (killing_threshold_multiplier - 1)

                is_active = death_increase > threshold_increase
                efficiency = death_increase / (threshold_increase + 1e-10) if threshold_increase > 0 else (
                    1.0 if death_increase > 0 else 0.0
                )

                target_evals.append({
                    "target_id": target_id_int,
                    "death_increase": death_increase,
                    "threshold_increase": threshold_increase,
                    "is_active": is_active,
                    "efficiency": efficiency,
                })

            # Representative target for this timepoint = highest efficiency
            # among all targets evaluated at this timepoint.
            best = max(target_evals, key=lambda d: d["efficiency"]) if target_evals else None

            if best is not None:
                is_active_killing = best["is_active"]
                targeted_track_id = best["target_id"] if is_active_killing else None
                death_signal_increase = best["death_increase"]
                killing_efficiency = best["efficiency"]
                killing_threshold_used = best["threshold_increase"]
            else:
                is_active_killing = False
                targeted_track_id = None
                death_signal_increase = 0.0
                killing_efficiency = 0.0
                killing_threshold_used = 0.0

            killing_results.append({
                "contact_event_id": event["contact_event_id"],
                "sample_name": sample_name,
                "immune_track_id": immune_track_id,
                "position_t": t,
                "death_signal_increase": death_signal_increase,
                "is_active_killing": is_active_killing,
                "killing_efficiency": killing_efficiency,
                "targeted_track_id": targeted_track_id,
                "killing_threshold_used": killing_threshold_used,
                "observation_complete": can_observe,
            })

    return pd.DataFrame(killing_results)


def find_advanced_features_csv(output_dir: Union[str, Path], cell_type: str) -> Optional[Path]:
    """
    Locate the active-killing advanced-features CSV for a cell type.

    run_active_killing_analysis() writes into a per-target subfolder (the
    target cell type name, or "combined" for multi-target runs), so this
    searches analysis/<cell_type>/active_killing/*/ rather than assuming a
    flat layout. Falls back to the legacy flat path for older runs. Returns
    the most recently modified match, or None if none exist.
    """
    active_killing_dir = Path(output_dir) / "analysis" / cell_type / "active_killing"
    if not active_killing_dir.exists():
        return None
    filename = f"BEHAV3D_{cell_type}_advanced_track_features.csv"
    candidates = list(active_killing_dir.glob(f"*/{filename}"))
    legacy = active_killing_dir / filename
    if legacy.exists():
        candidates.append(legacy)
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def run_active_killing_analysis(
    metadata: pd.DataFrame,
    output_dir: Union[str, Path],
    immune_cell_type: str = "tcell",
    target_cell_types: Optional[List[str]] = None,
    observation_window: int = 5,
    death_signal_column: str = "mean_dead_dye",
    min_contact_duration: int = 1,
    killing_threshold_multiplier: float = 1.5,
    absolute_killing_threshold: Optional[float] = None,
    save_results: bool = True,
    output_subfolder: str = ""
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """
    Run active killing analysis on BEHAV3D processed data.
    
    This function analyzes GLOBAL active killing (across all target types combined).
    It calculates active killing PER TIMEPOINT during qualifying contact events.
    
    Algorithm:
    1. Identify continuous contact events with ANY target cell (across all target types)
    2. Filter contacts by total duration (must be >= min_contact_duration)
    3. For each contact event, evaluate every touched target organoid independently at
       EACH contact timepoint t (sliding window):
       - Measure that organoid's own death signal increase from t to t + observation_window
       - Compare to either a multiplier of its signal at t, or a flat absolute increase
       - Classify as active killing if exceeds threshold; the highest-efficiency touched
         organoid at that timepoint is reported as the timepoint's target (this can differ
         from one contact timepoint to the next)
    
    Parameters
    ----------
    metadata : pd.DataFrame
        BEHAV3D metadata DataFrame
    output_dir : str or Path
        BEHAV3D output directory containing analysis results
    immune_cell_type : str
        Type of immune cell (e.g., 'tcell', 'macro', 'nk')
    target_cell_types : List[str], optional
        List of target cell types (e.g., ['organoid', 'organoid2']).
        If None, auto-detects all organoid types from metadata.
    observation_window : int
        Timepoints after each contact timepoint to measure death signal
    death_signal_column : str
        Column containing death signal in target tracks
    min_contact_duration : int
        Minimum TOTAL contact duration (continuous) in timepoints
    killing_threshold_multiplier : float
        Multiplier applied to a target organoid's own signal at each timepoint.
        Used when absolute_killing_threshold is None.
    absolute_killing_threshold : float, optional
        If provided, use this flat value as the killing threshold instead of the
        multiplier-based threshold. Useful when organoids start with near-zero signal.
        Default is None (uses killing_threshold_multiplier instead).
    save_results : bool
        Whether to save results to CSV files
        
    Returns
    -------
    df_killing_per_tp : pd.DataFrame
        Per-timepoint killing analysis for all qualifying contacts
    df_summary : pd.DataFrame
        Per-sample summary statistics
    stats : dict
        Overall statistics
    """
    print(f"--------------- Running Active Killing Analysis ---------------")
    print(f"Immune cell type: {immune_cell_type}")
    print(f"Algorithm: Global contacts (any target), per-organoid sliding-window evaluation")
    start_time = time.time()
    
    output_dir = Path(output_dir)
    
    # Auto-detect target cell types if not provided
    if target_cell_types is None:
        target_cell_types = detect_organoid_types_from_metadata(metadata)
        if not target_cell_types:
            raise ValueError("No organoid types detected in metadata. Please specify target_cell_types.")
    
    print(f"Target cell types: {target_cell_types}")
    print(f"Observation window: {observation_window} timepoints")
    print(f"Min contact duration: {min_contact_duration} timepoints")
    if absolute_killing_threshold is not None:
        print(f"Killing threshold mode: ABSOLUTE ({absolute_killing_threshold})")
    else:
        print(f"Killing threshold mode: MULTIPLIER ({killing_threshold_multiplier}x organoid's own signal at each timepoint)")
    
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
    
    # Load and combine ALL target cell tracks
    all_target_tracks = []
    
    for target_type in target_cell_types:
        target_feature_dir = output_dir / "analysis" / target_type / "track_features"
        target_tracks_path = target_feature_dir / f"BEHAV3D_{target_type}_combined_track_features_filtered.csv"
        
        if not target_tracks_path.exists():
            target_tracks_path = target_feature_dir / f"BEHAV3D_{target_type}_combined_track_features.csv"
        
        if target_tracks_path.exists():
            print(f"{get_current_time()} - Loading {target_type} tracks from {target_tracks_path}")
            df_target = pd.read_csv(target_tracks_path)
            df_target["_target_cell_type"] = target_type
            all_target_tracks.append(df_target)
        else:
            print(f"{get_current_time()} - Warning: {target_type} tracks not found at {target_tracks_path}")
    
    if not all_target_tracks:
        raise FileNotFoundError(f"Could not find any target cell tracks for types: {target_cell_types}")
    
    # Combine all target tracks
    df_target_tracks = pd.concat(all_target_tracks, ignore_index=True)
    print(f"{get_current_time()} - Combined {len(df_target_tracks)} target cell track rows")

    # Identify GLOBAL contact events (contact with ANY target type)
    print(f"{get_current_time()} - Identifying global contact events...")
    df_contact_events = identify_contact_events_global(
        df_immune_tracks=df_immune_tracks,
        target_cell_types=target_cell_types,
        min_contact_duration=min_contact_duration
    )
    
    if df_contact_events.empty:
        print(f"{get_current_time()} - No qualifying contact events found (duration >= {min_contact_duration})")
        # Still create advanced features with all zeros
        df_advanced = create_advanced_features_csv(
            df_immune_tracks=df_immune_tracks,
            df_killing_per_timepoint=pd.DataFrame(),
            observation_window=observation_window,
        )
        
        if save_results:
            if output_subfolder:
                results_dir = output_dir / "analysis" / immune_cell_type / "active_killing" / output_subfolder
            else:
                results_dir = output_dir / "analysis" / immune_cell_type / "active_killing"
            results_dir.mkdir(parents=True, exist_ok=True)
            advanced_features_path = results_dir / f"BEHAV3D_{immune_cell_type}_advanced_track_features.csv"
            df_advanced.to_csv(advanced_features_path, index=False)
            print(f"{get_current_time()} - Advanced features saved to {advanced_features_path}")
        
        return pd.DataFrame(), pd.DataFrame(), {"total_contacts": 0, "total_active_killing": 0}
    
    n_events = len(df_contact_events)
    n_timepoints = sum(len(e) for e in df_contact_events["contact_timepoints"])
    print(f"{get_current_time()} - Found {n_events} qualifying contact events ({n_timepoints} total contact timepoints)")
    
    # Analyze active killing PER TIMEPOINT
    print(f"{get_current_time()} - Analyzing active killing per timepoint...")
    df_killing_per_tp = analyze_active_killing_per_timepoint(
        df_immune_tracks=df_immune_tracks,
        df_target_tracks=df_target_tracks,
        df_contact_events=df_contact_events,
        observation_window=observation_window,
        death_signal_column=death_signal_column,
        killing_threshold_multiplier=killing_threshold_multiplier,
        absolute_killing_threshold=absolute_killing_threshold,
    )
    
    # Calculate summary statistics
    print(f"{get_current_time()} - Calculating summary statistics...")
    
    if not df_killing_per_tp.empty:
        df_summary = df_killing_per_tp.groupby("sample_name").agg(
            total_contact_timepoints=("position_t", "count"),
            active_killing_timepoints=("is_active_killing", "sum"),
            mean_death_signal_increase=("death_signal_increase", "mean"),
            median_death_signal_increase=("death_signal_increase", "median"),
            mean_killing_efficiency=("killing_efficiency", "mean"),
            n_contact_events=("contact_event_id", "nunique"),
        ).reset_index()
        
        df_summary["active_killing_rate"] = (
            df_summary["active_killing_timepoints"] / df_summary["total_contact_timepoints"]
        )
    else:
        df_summary = pd.DataFrame()
    
    # Overall statistics
    stats = {
        "total_contact_events": n_events,
        "total_contact_timepoints": len(df_killing_per_tp) if not df_killing_per_tp.empty else 0,
        "total_active_killing_timepoints": int(df_killing_per_tp["is_active_killing"].sum()) if not df_killing_per_tp.empty else 0,
        "overall_killing_rate": df_killing_per_tp["is_active_killing"].mean() if not df_killing_per_tp.empty else 0.0,
        "observation_window": observation_window,
        "min_contact_duration": min_contact_duration,
        "killing_threshold_multiplier": killing_threshold_multiplier,
        "absolute_killing_threshold": absolute_killing_threshold,
        "threshold_mode": "absolute" if absolute_killing_threshold is not None else "multiplier",
        "targeted_organoids_tracked": True
    }
    
    print(f"{get_current_time()} - Active killing analysis complete:")
    print(f"    Total qualifying contact events: {stats['total_contact_events']}")
    print(f"    Total contact timepoints analyzed: {stats['total_contact_timepoints']}")
    print(f"    Active killing timepoints: {stats['total_active_killing_timepoints']}")
    print(f"    Overall killing rate: {stats['overall_killing_rate']:.2%}")
    
    # Save results
    if save_results:
        if output_subfolder:
            results_dir = output_dir / "analysis" / immune_cell_type / "active_killing" / output_subfolder
        else:
            results_dir = output_dir / "analysis" / immune_cell_type / "active_killing"
        results_dir.mkdir(parents=True, exist_ok=True)
        
        if not df_killing_per_tp.empty:
            killing_events_path = results_dir / f"active_killing_per_timepoint_{immune_cell_type}.csv"
            summary_path = results_dir / f"active_killing_summary_{immune_cell_type}.csv"
            contact_events_path = results_dir / f"contact_events_{immune_cell_type}.csv"
            
            df_killing_per_tp.to_csv(killing_events_path, index=False)
            df_summary.to_csv(summary_path, index=False)
            df_contact_events.drop(columns=["contact_timepoints"]).to_csv(contact_events_path, index=False)
            
            print(f"{get_current_time()} - Results saved to {results_dir}")
        
        # Create and save advanced features CSV
        print(f"{get_current_time()} - Creating advanced features CSV...")
        df_advanced = create_advanced_features_csv(
            df_immune_tracks=df_immune_tracks,
            df_killing_per_timepoint=df_killing_per_tp,
            observation_window=observation_window,
        )
        
        advanced_features_path = results_dir / f"BEHAV3D_{immune_cell_type}_advanced_track_features.csv"
        df_advanced.to_csv(advanced_features_path, index=False)
        print(f"{get_current_time()} - Advanced features saved to {advanced_features_path}")
        
        # Print summary of new columns
        n_killing_contacts = df_advanced["is_active_killing"].sum()
        print(f"    Rows with active killing: {n_killing_contacts}")
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    
    return df_killing_per_tp, df_summary, stats


def create_advanced_features_csv(
    df_immune_tracks: pd.DataFrame,
    df_killing_per_timepoint: pd.DataFrame,
    observation_window: int = 5,
) -> pd.DataFrame:
    """
    Create an advanced features CSV that enriches the original immune cell track 
    features with GLOBAL per-timepoint active killing information.
    
    This function adds killing information for each timepoint:
    - is_active_killing: Boolean (True if active killing at this timepoint)
    - death_signal_increase_{N}tp: Death signal change over observation window
    - killing_efficiency: Ratio of actual death signal increase vs the targeted
      organoid's own threshold (its signal at that timepoint, scaled by the multiplier,
      or the flat absolute threshold)
    
    NO NAs: Where no contact occurs, is_active_killing=False and numeric values=0.0
    
    Parameters
    ----------
    df_immune_tracks : pd.DataFrame
        Original immune cell track features (per-timepoint data)
    df_killing_per_timepoint : pd.DataFrame
        Output from analyze_active_killing_per_timepoint
    observation_window : int
        Observation window size (for column naming)
        
    Returns
    -------
    df_advanced : pd.DataFrame
        Enhanced track features with killing information (NO NAs)
    """
    df = df_immune_tracks.copy()
    
    # Create suffix with observation window size
    window_suffix = f"_{observation_window}tp"
    
    # Initialize columns with default values (NO NAs)
    df["is_active_killing"] = False
    df["killing_efficiency"] = 0.0
    df[f"death_signal_increase{window_suffix}"] = 0.0
    df["targeted_track_id"] = -1
    df["contact_event_id"] = -1
    
    if df_killing_per_timepoint.empty:
        print("No killing events to merge - returning original features with zero killing columns")
        return df
    
    # Create a lookup key for merging (vectorized approach)
    df["_merge_key"] = (
        df["sample_name"].astype(str) + "_" + 
        df["TrackID"].astype(str) + "_" + 
        df["position_t"].astype(str)
    )
    
    df_killing_per_timepoint = df_killing_per_timepoint.copy()
    df_killing_per_timepoint["_merge_key"] = (
        df_killing_per_timepoint["sample_name"].astype(str) + "_" + 
        df_killing_per_timepoint["immune_track_id"].astype(str) + "_" + 
        df_killing_per_timepoint["position_t"].astype(str)
    )
    
    # Handle potential duplicates (same immune cell, same timepoint, multiple contact events)
    # Keep last to match original loop behavior which would overwrite with later values
    df_killing_deduped = df_killing_per_timepoint.drop_duplicates(
        subset=["_merge_key"], keep="last"
    )
    
    # Create indexed Series for vectorized lookup
    killing_indexed = df_killing_deduped.set_index("_merge_key")[
        ["is_active_killing", "killing_efficiency", "death_signal_increase", "targeted_track_id", "contact_event_id"]
    ]
    
    # Vectorized assignment using reindex - matches df's merge keys to killing data
    matched_values = killing_indexed.reindex(df["_merge_key"])
    
    # Update columns - use where/fillna to preserve defaults for non-matches
    df["is_active_killing"] = matched_values["is_active_killing"].fillna(False).astype(bool).values
    df["killing_efficiency"] = matched_values["killing_efficiency"].fillna(0.0).values
    df[f"death_signal_increase{window_suffix}"] = matched_values["death_signal_increase"].fillna(0.0).values
    df["targeted_track_id"] = matched_values["targeted_track_id"].fillna(-1).values
    df["contact_event_id"] = matched_values["contact_event_id"].fillna(-1).astype(int).values
    
    # Clean up merge key
    df.drop(columns=["_merge_key"], inplace=True)
    
    return df


def calculate_killing_dynamics_over_time(
    df_killing_per_timepoint: pd.DataFrame,
    df_immune_tracks: pd.DataFrame,
    time_bins: Optional[List[int]] = None,
    time_column: str = "position_t"
) -> pd.DataFrame:
    """
    Calculate how active killing rate changes over the course of the experiment.
    
    Parameters
    ----------
    df_killing_per_timepoint : pd.DataFrame
        Output from analyze_active_killing_per_timepoint
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
    if df_killing_per_timepoint.empty:
        return pd.DataFrame()
    
    if time_bins is None:
        max_t = df_immune_tracks["position_t"].max()
        time_bins = [0, max_t // 4, max_t // 2, 3 * max_t // 4, max_t + 1]
    
    df = df_killing_per_timepoint.copy()
    df["time_bin"] = pd.cut(
        df[time_column], 
        bins=time_bins, 
        labels=[f"{time_bins[i]}-{time_bins[i+1]}" for i in range(len(time_bins)-1)]
    )
    
    df_dynamics = df.groupby(["sample_name", "time_bin"]).agg(
        n_contact_timepoints=("position_t", "count"),
        n_active_killing=("is_active_killing", "sum"),
        mean_death_increase=("death_signal_increase", "mean"),
    ).reset_index()
    
    df_dynamics["killing_rate"] = df_dynamics["n_active_killing"] / df_dynamics["n_contact_timepoints"]
    
    return df_dynamics




def calculate_invasiveness_single_timepoint(args):
    """
    Calculate organoid invasiveness features for immune cells at a single timepoint.
    
    Invasiveness measures what percentage of an immune cell's surface is in contact
    with organoid surfaces. A cell is considered invasive if >50% of its surface
    is in contact.
    
    Parameters
    ----------
    args : tuple
        Contains:
        - t: timepoint index
        - current_cell_segments_path: path to immune cell segments
        - organoid_segments_paths: dict of organoid segment paths
        - element_size_x, element_size_y, element_size_z: voxel spacing (um)
        - contact_threshold: distance threshold (um) for determining contact
        - calculate_from: cell type being analyzed (should be immune type)
    
    Returns
    -------
    pd.DataFrame
        Columns: TrackID, position_t, <type>_invasiveness, <type>_invasiveness_perc, any_org_invasiveness_perc
    """
    from scipy.ndimage import distance_transform_edt
    import math
    
    (
        t,
        current_cell_segments_path,
        organoid_segments_paths,
        element_size_x,
        element_size_y,
        element_size_z,
        contact_threshold,
        calculate_from
    ) = args

    # Load current cell type's segments for this timepoint
    current_segments = np.asarray(load_image(current_cell_segments_path)[t])
    
    # Load all organoid types' segments
    organoid_segments_dict = {}
    for org_type, org_path in organoid_segments_paths.items():
        organoid_segments_dict[org_type] = np.asarray(load_image(org_path)[t])
    
    df_invasiveness = []
    segment_ids = np.unique(current_segments)
    
    # If there are no foreground segments at this timepoint, return an empty DataFrame
    foreground_ids = [sid for sid in segment_ids if sid != 0]
    if not foreground_ids:
        return pd.DataFrame(columns=["TrackID", "position_t"])
    
    for segment_id in foreground_ids:
        stack_max_z, stack_max_y, stack_max_x = current_segments.shape
        seg_locs = np.argwhere(current_segments == segment_id)
        min_z, min_y, min_x = seg_locs.min(axis=0)
        max_z, max_y, max_x = seg_locs.max(axis=0)
        
        z_ext = 2 * math.ceil(contact_threshold / element_size_z)
        y_ext = 2 * math.ceil(contact_threshold / element_size_y)
        x_ext = 2 * math.ceil(contact_threshold / element_size_x)
        
        slicer = (
            slice(max(0, min_z - z_ext), min(stack_max_z, max_z + z_ext + 1)),
            slice(max(0, min_y - y_ext), min(stack_max_y, max_y + y_ext + 1)),
            slice(max(0, min_x - x_ext), min(stack_max_x, max_x + x_ext + 1))
        )
        
        seg_cutout = current_segments[slicer]
        
        # Calculate distance transform from current cell boundary
        real_distances = distance_transform_edt(
            seg_cutout != segment_id,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        
        # Define "surface" as pixels within 2 um of cell boundary
        surface_threshold = 2.0  # um
        surface_mask = real_distances <= surface_threshold
        total_surface_pixels = np.sum(surface_mask)
        
        if total_surface_pixels == 0:
            # Cell has no surface (single pixel?), skip
            continue
        
        invasiveness_data = {
            'TrackID': segment_id,
            'position_t': t,
        }
        
        any_invasiveness_list = []
        invasiveness_perc_list = []
        
        # Calculate invasiveness for each organoid type
        for org_type, org_segments in organoid_segments_dict.items():
            org_cutout = org_segments[slicer]
            
            # Count surface pixels in contact with this organoid type
            org_contact_mask = (real_distances <= contact_threshold) & (org_cutout != 0)
            contacted_surface_pixels = np.sum(org_contact_mask)
            
            # Calculate percentage
            invasiveness_perc = (contacted_surface_pixels / total_surface_pixels) * 100.0
            
            # Boolean: invasive if >= 50% of surface is in contact
            invasiveness_bool = invasiveness_perc >= 50.0
            
            invasiveness_data[f'{org_type}_invasiveness'] = invasiveness_bool
            invasiveness_data[f'{org_type}_invasiveness_perc'] = invasiveness_perc
            
            any_invasiveness_list.append(invasiveness_bool)
            invasiveness_perc_list.append(invasiveness_perc)
        
        # Add aggregate: True if invasive against ANY organoid type
        invasiveness_data['any_org_invasiveness'] = any(any_invasiveness_list)
        invasiveness_data['any_org_invasiveness_perc'] = max(invasiveness_perc_list) if invasiveness_perc_list else 0.0
        
        df_invasiveness.append(pd.DataFrame([invasiveness_data]))
    
    if df_invasiveness:
        return pd.concat(df_invasiveness, ignore_index=True)
    else:
        return pd.DataFrame(columns=["TrackID", "position_t"])

