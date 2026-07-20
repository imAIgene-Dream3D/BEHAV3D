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
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


# ---------------------------------------------------------------------------
# Shared helpers: line-condition grouping and absolute analysis-time window.
#
# All of these are additive and OFF by default: callers that do not pass a
# line-condition flag / analysis period get exactly the previous behaviour.
# They are defined here (and imported by invasiveness_analysis) so the two
# analyses share one implementation instead of repeating the logic.
# ---------------------------------------------------------------------------

def resolve_line_condition_map(df: pd.DataFrame, immune_types: list):
    """
    Map each ``sample_name`` to its immune *line condition* value.

    Uses the metadata convention ``im_{immune_type}_line_condition``. These
    per-sample columns are carried into the filtered track-feature CSVs, so a
    lookup on the loaded DataFrame is enough (no separate metadata read). The
    first immune type in ``immune_types`` with a usable column wins; a couple
    of fallbacks keep it robust to naming drift.

    Returns
    -------
    (mapping, column_name) : (dict, str | None)
        ``mapping`` is ``{sample_name: condition_value}`` (values as str).
        ``column_name`` is the resolved source column (``None`` if none found).
    """
    candidates = [f"im_{ct}_line_condition" for ct in (immune_types or [])]
    # Fallbacks: any immune line_condition column, then any *_line_condition.
    candidates += [
        c for c in df.columns
        if c.startswith("im_") and c.endswith("_line_condition")
    ]
    candidates += [c for c in df.columns if c.endswith("_line_condition")]

    seen = set()
    for col in candidates:
        if col in seen or col not in df.columns:
            continue
        seen.add(col)
        sub = df.dropna(subset=[col])
        if sub.empty:
            continue
        mapping = (
            sub.drop_duplicates("sample_name")
            .set_index("sample_name")[col]
            .astype(str)
            .to_dict()
        )
        if mapping:
            return mapping, col
    return {}, None


def add_line_condition_column(
    df: pd.DataFrame, immune_types: list, target_col: str = "line_condition",
):
    """
    Attach a per-sample ``line_condition`` column to ``df`` in place.

    Missing values (samples with no condition) become ``"(no line condition)"``
    so grouping never drops rows silently. Returns ``(df, source_column)``.
    """
    mapping, col = resolve_line_condition_map(df, immune_types)
    if mapping:
        df[target_col] = (
            df["sample_name"].map(mapping).fillna("(no line condition)")
        )
    else:
        df[target_col] = "(no line condition)"
    return df, col


def add_minutes_from_start(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add ``analysis_time_min`` = minutes elapsed from the start of each movie.

    Per sample, ``t0`` is that movie's minimum ``time``; the offset is converted
    to minutes with :func:`_time_unit_to_minutes_factor` so movies with
    different frame intervals are directly comparable. Falls back to
    ``position_t`` (treated as minutes) when no ``time`` column exists so the
    filter still works on legacy CSVs. Operates in place and returns ``df``.
    """
    if "time" in df.columns:
        to_min = _time_unit_to_minutes_factor(df)
        t0 = df.groupby("sample_name")["time"].transform("min")
        df["analysis_time_min"] = (df["time"] - t0) * to_min
    else:
        t0 = df.groupby("sample_name")["position_t"].transform("min")
        df["analysis_time_min"] = (df["position_t"] - t0).astype(float)
    return df


def _normalize_period(analysis_period_min):
    """Coerce a period spec to ``(start, end)`` floats/None; tolerant of junk."""
    if analysis_period_min is None:
        return None, None
    try:
        start, end = analysis_period_min
    except (TypeError, ValueError):
        return None, None
    start = None if start is None or start == "" else float(start)
    end = None if end is None or end == "" else float(end)
    return start, end


def apply_analysis_period(df: pd.DataFrame, analysis_period_min) -> pd.DataFrame:
    """
    Keep only rows inside the absolute analysis window (minutes from movie
    start). ``analysis_period_min`` is ``(start, end)`` in minutes; either bound
    may be ``None`` (open-ended). Returns a filtered copy; a no-op (returns the
    same object) when the period is ``None`` / fully open.
    """
    start, end = _normalize_period(analysis_period_min)
    if start is None and end is None:
        return df
    if "analysis_time_min" not in df.columns:
        df = add_minutes_from_start(df.copy())
    mask = pd.Series(True, index=df.index)
    if start is not None:
        mask &= df["analysis_time_min"] >= start
    if end is not None:
        mask &= df["analysis_time_min"] <= end
    return df[mask].copy()


def period_bounds_by_sample(df: pd.DataFrame, analysis_period_min):
    """
    Translate the absolute minutes window into per-sample ``position_t`` bounds.

    Needed to filter the active-killing event CSVs, which only carry
    ``position_t`` (no real ``time``). Returns ``{sample_name: (t_start, t_end)}``
    inclusive, omitting samples with no timepoints in the window. Returns
    ``None`` when the period is open (i.e. no filtering required).
    """
    start, end = _normalize_period(analysis_period_min)
    if start is None and end is None:
        return None
    if "analysis_time_min" not in df.columns:
        df = add_minutes_from_start(df.copy())
    bounds = {}
    for sample, grp in df.groupby("sample_name"):
        sub = grp
        if start is not None:
            sub = sub[sub["analysis_time_min"] >= start]
        if end is not None:
            sub = sub[sub["analysis_time_min"] <= end]
        if sub.empty:
            continue
        bounds[str(sample)] = (
            int(sub["position_t"].min()), int(sub["position_t"].max()),
        )
    return bounds


def _in_bounds(period_bounds: dict, sample_name, position_t) -> bool:
    """True if ``position_t`` falls within ``period_bounds`` for that sample.

    Samples absent from ``period_bounds`` have no timepoints in the window and
    return False. Non-finite timepoints return False.
    """
    b = period_bounds.get(str(sample_name))
    if b is None:
        return False
    try:
        t = float(position_t)
    except (TypeError, ValueError):
        return False
    if np.isnan(t):
        return False
    return b[0] <= t <= b[1]


def _sample_period_bounds(period_bounds: dict, sample_name):
    """Return inclusive ``(lo, hi)`` for *sample_name*, or ``None``."""
    if period_bounds is None:
        return None
    return period_bounds.get(str(sample_name))


def _event_end_t(row) -> float:
    """Best-effort inclusive end timepoint for a contact-event row."""
    if "contact_end_t" in row.index and pd.notna(row["contact_end_t"]):
        return float(row["contact_end_t"])
    if "contact_duration" in row.index and pd.notna(row["contact_duration"]):
        try:
            dur = int(row["contact_duration"])
        except (TypeError, ValueError):
            dur = 0
        if dur > 0 and pd.notna(row.get("contact_start_t")):
            return float(row["contact_start_t"]) + dur - 1
    if pd.notna(row.get("contact_start_t")):
        return float(row["contact_start_t"])
    return float("nan")


def _interval_overlaps_window(start_t, end_t, lo, hi) -> bool:
    """True when inclusive ``[start_t, end_t]`` overlaps ``[lo, hi]``."""
    try:
        s, e = float(start_t), float(end_t)
    except (TypeError, ValueError):
        return False
    if np.isnan(s) or np.isnan(e):
        return False
    if s > e:
        s, e = e, s
    return s <= hi and e >= lo


def _overlap_duration(start_t, end_t, lo, hi) -> int:
    """Inclusive timepoint count in the overlap of event and window."""
    try:
        s, e = float(start_t), float(end_t)
    except (TypeError, ValueError):
        return 0
    if np.isnan(s) or np.isnan(e):
        return 0
    if s > e:
        s, e = e, s
    overlap_lo = max(s, lo)
    overlap_hi = min(e, hi)
    if overlap_lo > overlap_hi:
        return 0
    return int(overlap_hi - overlap_lo + 1)


def _clip_contact_events_to_period(df_ev: pd.DataFrame, period_bounds: dict) -> pd.DataFrame:
    """Keep events overlapping the per-sample window; clip ``contact_duration``.

    Events that start before the window but extend into it are included (with
  duration counted only inside the window). Events with zero overlap are dropped.
    """
    if period_bounds is None or df_ev.empty or "contact_start_t" not in df_ev.columns:
        return df_ev

    rows = []
    for _, row in df_ev.iterrows():
        bounds = _sample_period_bounds(period_bounds, row["sample_name"])
        if bounds is None:
            continue
        lo, hi = bounds
        start_t = row["contact_start_t"]
        end_t = _event_end_t(row)
        if not _interval_overlaps_window(start_t, end_t, lo, hi):
            continue
        dur = _overlap_duration(start_t, end_t, lo, hi)
        if dur <= 0:
            continue
        r = row.copy()
        r["contact_duration"] = dur
        rows.append(r)
    if not rows:
        return df_ev.iloc[0:0].copy()
    return pd.DataFrame(rows)


def format_period_label(analysis_period_min) -> str:
    """Human-readable label for the analysis window, e.g. '10-60 min'."""
    start, end = _normalize_period(analysis_period_min)
    if start is None and end is None:
        return "full movie"
    lo = "start" if start is None else f"{start:g}"
    hi = "end" if end is None else f"{end:g}"
    return f"{lo}-{hi} min"


# ---------------------------------------------------------------------------
# Timepoint-based temporal range (position_t) — new API
# ---------------------------------------------------------------------------

def _normalize_period_t(analysis_period_t):
    """Coerce a timepoint period spec to ``(start, end)`` ints/None."""
    if analysis_period_t is None:
        return None, None
    try:
        start, end = analysis_period_t
    except (TypeError, ValueError):
        return None, None
    start = None if start is None or start == "" else int(start)
    end = None if end is None or end == "" else int(end)
    return start, end


def apply_analysis_period_timepoints(
    df: pd.DataFrame, analysis_period_t,
) -> pd.DataFrame:
    """Keep only rows where ``position_t`` is inside ``[start_t, end_t]``.

    ``analysis_period_t`` is ``(start_t, end_t)`` in raw timepoints; either
    bound may be ``None`` (open-ended).  Returns a filtered copy; a no-op
    when the period is ``None`` / fully open.
    """
    start, end = _normalize_period_t(analysis_period_t)
    if start is None and end is None:
        return df
    mask = pd.Series(True, index=df.index)
    if start is not None:
        mask &= df["position_t"] >= start
    if end is not None:
        mask &= df["position_t"] <= end
    return df[mask].copy()


def format_period_label_timepoints(analysis_period_t) -> str:
    """Human-readable label for a timepoint window, e.g. 'T 10–50'."""
    start, end = _normalize_period_t(analysis_period_t)
    if start is None and end is None:
        return "all timepoints"
    lo = "start" if start is None else str(start)
    hi = "end" if end is None else str(end)
    return f"T {lo}–{hi}"


def period_bounds_by_sample_t(df: pd.DataFrame, analysis_period_t):
    """Per-sample ``position_t`` bounds from a timepoint-based window.

    Returns ``{sample_name: (t_start, t_end)}`` inclusive, or ``None``
    when no filtering is required.
    """
    start, end = _normalize_period_t(analysis_period_t)
    if start is None and end is None:
        return None
    bounds = {}
    for sample, grp in df.groupby("sample_name"):
        sub = grp
        if start is not None:
            sub = sub[sub["position_t"] >= start]
        if end is not None:
            sub = sub[sub["position_t"] <= end]
        if sub.empty:
            continue
        bounds[str(sample)] = (
            int(sub["position_t"].min()), int(sub["position_t"].max()),
        )
    return bounds


def _build_output_suffix(
    group_by: str,
    analysis_period_t,
    annotate_line_condition: bool = False,
) -> str:
    """Build a filename suffix encoding key settings so runs don't overwrite."""
    parts = []
    gb_short = {"organoid_type": "orgtype", "treatment": "treatment"}
    parts.append(gb_short.get(group_by, group_by))
    start, end = _normalize_period_t(analysis_period_t)
    if start is None and end is None:
        parts.append("Tall")
    else:
        lo = "s" if start is None else str(start)
        hi = "e" if end is None else str(end)
        parts.append(f"T{lo}-{hi}")
    if annotate_line_condition:
        parts.append("annoLC")
    return "_" + "_".join(parts)


def run_interaction_analysis(
    output_dir: str,
    cell_type: str,
    interacting_cell_types: list,
    dead_threshold: float = 0.02,
    df_tracks_path: str = None,
    show_plots: bool = True,
    group_by_line_condition: bool = False,
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
    group_by_line_condition : bool
        When True, the two *overall* plots are split by the interacting immune
        cell's line condition (``im_{type}_line_condition``): the cumulative
        curve draws one line per condition, and the alive-vs-dead plot colours
        by condition (solid = survives, dashed = dies). Default False keeps the
        original single-curve behaviour. Per-sample plots are unaffected.

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

        # Optional per-immune-type line-condition split column for the overall
        # plots. Attached per-ct so each interacting type uses its own
        # im_{ct}_line_condition values.
        split_col = None
        if group_by_line_condition:
            df, _lc_used = add_line_condition_column(
                df, [ct], target_col=f"__lc_{ct}",
            )
            split_col = f"__lc_{ct}"
            if _lc_used is None:
                print(f"   ⚠️ No line_condition column found for {ct}; "
                      "falling back to combined curves.")
                split_col = None

        # Generate plots and PDF
        figs = generate_interaction_plots(
            df=df,
            cell_type=cell_type,
            interacting_type=ct,
            cumulative_col=cumulative_col,
            has_dead_column=has_dead_column,
            n_organoids=n_organoids,
            n_samples=n_samples,
            split_col=split_col,
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
    split_col: str = None,
) -> dict:
    """
    Generate all interaction plots.

    Parameters
    ----------
    split_col : str, optional
        Column name to split the overall plots by (used for line-condition
        grouping). When None, the overall plots are single combined curves.

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
        df, cell_type, interacting_type, cumulative_col, n_organoids, n_samples,
        split_col=split_col,
    )
    
    # Plot 2: Alive vs Dead overall
    if has_dead_column and "survives" in df.columns:
        fig = plot_alive_vs_dead_overall(
            df, cell_type, interacting_type, cumulative_col, n_organoids,
            split_col=split_col,
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
    split_col: str = None,
) -> plt.Figure:
    """Plot cumulative interactions over time (all samples combined).

    When ``split_col`` is given (line-condition grouping), one mean +/- SEM
    curve is drawn per value of that column instead of a single combined curve.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    if split_col and split_col in df.columns and df[split_col].notna().any():
        groups = sorted(df[split_col].dropna().unique())
        palette = sns.color_palette("tab10", max(len(groups), 1))
        for gi, gval in enumerate(groups):
            sub = df[df[split_col] == gval]
            stats = (
                sub.groupby("position_t")[cumulative_col]
                .agg(["mean", "std", "count"]).reset_index()
            )
            stats["sem"] = stats["std"] / np.sqrt(stats["count"])
            n_g = sub.groupby(["sample_name", "TrackID"]).ngroups
            ax.plot(
                stats["position_t"], stats["mean"],
                linewidth=2, color=palette[gi], label=f"{gval} (n={n_g})",
            )
            ax.fill_between(
                stats["position_t"],
                stats["mean"] - stats["sem"],
                stats["mean"] + stats["sem"],
                alpha=0.25, color=palette[gi],
            )
        ax.legend(title="Line condition", fontsize=9)
    else:
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
        fontsize=14, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def plot_alive_vs_dead_overall(
    df: pd.DataFrame,
    cell_type: str,
    interacting_type: str,
    cumulative_col: str,
    n_organoids: int,
    split_col: str = None,
) -> plt.Figure:
    """Plot alive vs dead comparison (all samples combined).

    When ``split_col`` is given (line-condition grouping), curves are coloured
    by that column's value with survival encoded as line style
    (solid = survives, dashed = dies).
    """
    if "survives" not in df.columns:
        return None

    if split_col and split_col in df.columns and df[split_col].notna().any():
        fig, ax = plt.subplots(figsize=(10, 6))
        groups = sorted(df[split_col].dropna().unique())
        palette = sns.color_palette("tab10", max(len(groups), 1))
        color_map = dict(zip(groups, palette))
        for gval in groups:
            for survives, style, fate_lbl in [
                (True, "-", "Survives"), (False, "--", "Dies"),
            ]:
                sub = df[(df[split_col] == gval) & (df["survives"] == survives)]
                n_sub = sub.groupby(["sample_name", "TrackID"]).ngroups
                if n_sub == 0:
                    continue
                stats = (
                    sub.groupby("position_t")[cumulative_col]
                    .agg(["mean", "std", "count"]).reset_index()
                )
                stats["sem"] = stats["std"] / np.sqrt(stats["count"])
                ax.plot(
                    stats["position_t"], stats["mean"],
                    linewidth=2, color=color_map[gval], linestyle=style,
                    label=f"{gval} — {fate_lbl} (n={n_sub})",
                )
                ax.fill_between(
                    stats["position_t"],
                    stats["mean"] - stats["sem"],
                    stats["mean"] + stats["sem"],
                    alpha=0.15, color=color_map[gval],
                )
        ax.set_xlabel("Time (frames)", fontsize=12)
        ax.set_ylabel(f"Cumulative {interacting_type} contacts", fontsize=12)
        ax.set_title(
            f"Cumulative {interacting_type} Interactions by line condition: "
            f"Surviving vs Dying {cell_type}s\n(Total n = {n_organoids})",
            fontsize=14,
        )
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        return fig

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
        fontsize=14, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


# ---------------------------------------------------------------------------
# Interaction Overview (multi-organoid / cross-target)
# ---------------------------------------------------------------------------

def run_multi_organoid_interaction_comparison(
    output_dir: str,
    organoid_types: list,
    immune_types: list,
    metadata: pd.DataFrame,
    group_by: str = "organoid_type",
    time_window_min: float = 60,
    show_plots: bool = True,
    analysis_period_t: tuple = None,
    annotate_line_condition: bool = False,
    # Backward compat: old callers may pass analysis_period_min (minutes).
    # When set AND analysis_period_t is None, it is used instead.
    analysis_period_min: tuple = None,
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
        Time window in minutes before TOD for the *cumulative-to-death* curve
        only. This is independent of ``analysis_period_t``.
    show_plots : bool
        Whether to display figures inline.
    analysis_period_t : tuple, optional
        Temporal range ``(start_t, end_t)`` in raw timepoints (``position_t``).
        Restricts data for the violin ("cumulative interactions per organoid")
        and the active-killing dashboard. Either bound may be ``None`` for
        open-ended. The cumulative-to-death curve keeps its own
        ``time_window_min`` and is NOT affected. Default ``None`` (all
        timepoints).
    annotate_line_condition : bool
        When True, the violin and dashboard add the immune line condition as
        an annotation (hue / line style) within the selected ``group_by``
        grouping. Default False.
    analysis_period_min : tuple, optional
        **Deprecated** — kept for backward compatibility. When
        ``analysis_period_t`` is None and this is set, the old minutes-based
        window is used instead.

    Returns
    -------
    dict  with keys "violin_fig", "curve_fig", "pdf_path", "summary_csv_path".
    """
    print("--------------- Interaction Overview ---------------")
    start = time.time()

    # Resolve analysis_period_t vs legacy analysis_period_min.
    if analysis_period_t is None and analysis_period_min is not None:
        # Backward compat — convert old minutes-based to same variable
        analysis_period_t = analysis_period_min
        _use_minutes_fallback = True
    else:
        _use_minutes_fallback = False

    output_dir = Path(output_dir)
    results_dir = output_dir / "analysis" / "multi_organoid_comparison"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Build output suffix for non-overwriting filenames.
    out_suffix = _build_output_suffix(
        group_by, analysis_period_t, annotate_line_condition,
    )

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
    # 1b. Line-condition annotation + temporal range (both optional)
    # ------------------------------------------------------------------
    # Resolve line condition for annotation (NOT as primary grouping).
    lc_map = {}
    lc_col = None
    lc_maps = {}  # per-immune-type: {immune_type: {sample_name: line_condition}}
    if annotate_line_condition:
        lc_map, lc_col = resolve_line_condition_map(df_all, active_immune)
        if lc_map:
            df_all["line_condition"] = (
                df_all["sample_name"].map(lc_map).fillna("(no line condition)")
            )
            print(f"   Annotating by line condition (from '{lc_col}')")
        else:
            print("   ⚠️ No line_condition column found — "
                  "skipping line condition annotation.")
        # Per-immune-type line condition maps so each immune type's own
        # ``im_{type}_line_condition`` column drives its own violins/lines.
        for ct in active_immune:
            col = f"im_{ct}_line_condition"
            if col in df_all.columns:
                m = (
                    df_all.dropna(subset=[col])
                    .drop_duplicates("sample_name")
                    .set_index("sample_name")[col]
                    .astype(str)
                    .to_dict()
                )
                if m:
                    lc_maps[ct] = m

    # Determine the temporal range (timepoint-based or legacy minutes).
    if _use_minutes_fallback:
        # Legacy path: convert minutes to internal filtering
        start_min, end_min = _normalize_period(analysis_period_min)
        has_period = not (start_min is None and end_min is None)
        if has_period:
            df_all = add_minutes_from_start(df_all)
            print(f"   Analysis window applied to violin + dashboard: "
                  f"{format_period_label(analysis_period_min)}")
    else:
        start_t, end_t = _normalize_period_t(analysis_period_t)
        has_period = not (start_t is None and end_t is None)
        if has_period:
            print(f"   Temporal range applied to violin + dashboard: "
                  f"{format_period_label_timepoints(analysis_period_t)}")

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

    # Attach line condition (per sample) to the track summary for annotation.
    if annotate_line_condition and lc_map:
        track_summary["line_condition"] = (
            track_summary["sample_name"].map(lc_map).fillna("(no line condition)")
        )

    # Total contacts across all selected immune types (for "by organoid_type" grouping)
    im_max_cols = [f"max_{ct}_contact" for ct in active_immune]
    track_summary["total_contacts"] = track_summary[im_max_cols].sum(axis=1)
    track_summary["fate"] = track_summary["survives"].map({True: "Live", False: "Dying"})

    # Save summary
    summary_csv = results_dir / f"multi_organoid_interaction_summary{out_suffix}.csv"
    track_summary.to_csv(summary_csv, index=False)
    print(f"   Summary saved: {summary_csv.name}")

    # ------------------------------------------------------------------
    # 2a. Period-restricted per-track summary for the violin (optional)
    # ------------------------------------------------------------------
    track_summary_violin = track_summary
    if has_period:
        if _use_minutes_fallback:
            df_period = apply_analysis_period(df_all, analysis_period_min)
        else:
            df_period = apply_analysis_period_timepoints(df_all, analysis_period_t)
        ts_v = track_summary.copy()
        for ct in active_immune:
            contact_col = f"{ct}_contact"
            mcol = f"max_{ct}_contact"
            if contact_col in df_period.columns:
                windowed = (
                    df_period.groupby(track_keys)[contact_col]
                    .sum()
                    .rename(mcol)
                    .reset_index()
                )
                ts_v = ts_v.drop(columns=[mcol], errors="ignore").merge(
                    windowed, on=track_keys, how="left",
                )
                ts_v[mcol] = ts_v[mcol].fillna(0)
        ts_v["total_contacts"] = ts_v[im_max_cols].sum(axis=1)
        track_summary_violin = ts_v

    # Per-sample bounds for filtering event CSVs to the window.
    if has_period:
        if _use_minutes_fallback:
            dash_period_bounds = period_bounds_by_sample(df_all, analysis_period_min)
        else:
            dash_period_bounds = period_bounds_by_sample_t(df_all, analysis_period_t)
    else:
        dash_period_bounds = None

    # ------------------------------------------------------------------
    # 2b. Load active killing data (if available)
    # ------------------------------------------------------------------
    df_contact_events = _load_active_killing_data(
        output_dir, active_immune, organoid_types,
        track_summary=track_summary,
        period_bounds=dash_period_bounds,
        line_condition_map=(lc_map if annotate_line_condition and lc_map else None),
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

    # Build period label for plot annotation.
    if has_period:
        if _use_minutes_fallback:
            period_label = format_period_label(analysis_period_min)
        else:
            period_label = format_period_label_timepoints(analysis_period_t)
    else:
        period_label = None

    fig_violin = plot_interaction_violin_comparison(
        track_summary_violin, group_by, active_immune, has_dead_data,
        period_label=period_label,
        annotate_line_condition=annotate_line_condition and bool(lc_map),
        lc_maps=lc_maps,
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
        annotate_line_condition=annotate_line_condition and bool(lc_map),
        lc_maps=lc_maps,
    )
    if fig_curve is not None:
        figs["cumulative_to_death"] = fig_curve

    if curve_data_container:
        curve_csv = results_dir / f"multi_organoid_cumulative_to_death_curves_min{out_suffix}.csv"
        curve_df = curve_data_container[0]
        curve_df.to_csv(curve_csv, index=False)
        print(f"   Curve data saved (time in minutes): {curve_csv.name}")

    fig_dashboard = plot_interaction_overview_dashboard(
        track_summary, df_contact_events, active_immune, has_dead_data,
        group_by=group_by,
        period_label=period_label,
        annotate_line_condition=annotate_line_condition and bool(lc_map),
        lc_maps=lc_maps,
    )
    if fig_dashboard is not None:
        figs["dashboard"] = fig_dashboard

    # Save PDF with non-overwriting filename.
    pdf_path = results_dir / f"multi_organoid_interaction_comparison{out_suffix}.pdf"
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
    print(f"\n   Interaction Overview complete! ({elapsed:.1f}s)")
    return {
        "summary_csv_path": summary_csv,
        "pdf_path": pdf_path,
    }


def plot_interaction_violin_comparison(
    track_summary: pd.DataFrame,
    group_by: str,
    immune_types: list,
    has_dead_data: bool,
    period_label: str = None,
    annotate_line_condition: bool = False,
    lc_maps: dict = None,
) -> plt.Figure:
    # ``lc_maps`` = {immune_type: {sample_name: line_condition}}. When
    # annotation is requested we DISAGGREGATE: instead of pooling immune types
    # / line conditions into one violin (and marking them with symbols), we
    # draw a SEPARATE violin for each organoid_type x immune_type x
    # line_condition combination, side by side in the same figure.
    lc_maps = lc_maps or {}

    def _lc_for(ct):
        m = {str(k): str(v) for k, v in lc_maps.get(ct, {}).items()}
        return m

    # Decide whether disaggregation is meaningful: more than one immune type,
    # or at least one immune type spanning more than one line condition.
    disaggregate = False
    if annotate_line_condition:
        if len([c for c in immune_types if f"max_{c}_contact" in track_summary.columns]) > 1:
            disaggregate = True
        else:
            for ct in immune_types:
                if len(set(_lc_for(ct).values())) > 1:
                    disaggregate = True
                    break

    per_immune = False
    if disaggregate:
        # One point per (organoid, immune_type); condition label carries the
        # immune type and line condition so seaborn renders distinct violins.
        per_immune = True
        rows = []
        for ct in immune_types:
            col = f"max_{ct}_contact"
            if col not in track_summary.columns:
                continue
            sub = track_summary[["organoid_type", "sample_name", "TrackID", "fate", col]].copy()
            sub = sub.rename(columns={col: "cumulative_interactions"})
            lc_map = _lc_for(ct)
            sub["line_condition"] = (
                sub["sample_name"].astype(str).map(lc_map).fillna("(no LC)")
            )
            sub["immune_type"] = ct
            if group_by == "treatment":
                sub["condition"] = ct + "\n" + sub["line_condition"].astype(str)
            else:
                sub["condition"] = (
                    sub["organoid_type"].astype(str)
                    + "\n" + ct
                    + "\n" + sub["line_condition"].astype(str)
                )
            rows.append(sub)
        df_plot = pd.concat(rows, ignore_index=True)
        x_label = (
            "Immune type / line condition" if group_by == "treatment"
            else "Organoid type / immune / line condition"
        )
    elif group_by == "treatment":
        per_immune = True
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
        # Light fills for violins; full-strength shared red/green for box
        # edges and strip dots so the colour vocabulary matches the dashboard.
        palette = {k: _FATE_COLORS_LIGHT[k] for k in hue_order}
        dark_palette = {k: _FATE_COLORS[k] for k in hue_order}
        strip_palette = {k: _FATE_COLORS[k] for k in hue_order}
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

    # Boxplot: filled by seaborn with `dark_palette`, then we hollow each box
    # so the strip points underneath stay visible while the per-fate colour
    # is preserved on the edges, whiskers and median tick.
    # Use a copy of common_kw with the darker palette + gap for narrow centered box.
    # Only override the palette when we actually have a hue mapping; otherwise
    # seaborn falls back to looking up `x` values in the palette dict and an
    # empty dict raises "palette is missing keys" (and triggers a FutureWarning
    # about passing palette without hue).
    box_kw = common_kw.copy()
    if has_dead_data and dark_palette:
        box_kw["palette"] = dark_palette
    sns.boxplot(
        **box_kw,
        gap=0.6,
        fliersize=0, legend=False,
        boxprops=dict(linewidth=1.4),
        whiskerprops=dict(linewidth=1.0),
        capprops=dict(linewidth=1.0),
        showcaps=True,
        medianprops=dict(linewidth=2.2),
    )

    # Hollow out the boxes so the underlying stripplot is visible. Each
    # PathPatch added by sns.boxplot keeps its palette colour as the
    # edgecolor; only the face becomes transparent.
    for patch in ax.patches:
        if isinstance(patch, mpatches.PathPatch):
            face = patch.get_facecolor()
            patch.set_edgecolor(face)
            patch.set_facecolor("none")

    # Stripplot — slightly larger / punchier dots now that the box is open.
    strip_kw = common_kw.copy()
    strip_kw.pop("width", None)
    if strip_palette:
        strip_kw["palette"] = strip_palette
    sns.stripplot(
        **strip_kw,
        size=3.5, alpha=0.6, jitter=True, legend=False,
        edgecolor="white", linewidth=0.3,
    )

    # Legend: organoid fate (colour). Immune type / line condition are now
    # disaggregated onto the x-axis (separate violins), so no marker legend.
    legend_handles = []
    if hue_order:
        legend_handles.extend([
            mpatches.Patch(
                facecolor=_FATE_COLORS[f],
                edgecolor="white",
                label=_FATE_DISPLAY.get(f, f),
            )
            for f in hue_order
        ])
    if legend_handles:
        existing = ax.get_legend()
        if existing is not None:
            existing.remove()
        ax.legend(
            handles=legend_handles,
            fontsize=9, loc="best",
            frameon=True, framealpha=0.85, edgecolor="#cccccc",
            fancybox=True, shadow=False,
            title="Organoid fate",
            title_fontsize=9,
        )

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
    # at time of death, matching how track features are computed). When
    # several immune types are selected they are summed -- noted on the
    # label so the unit is unambiguous.
    ax.set_ylabel(
        "Cumulative contact timepoints per organoid\n"
        + ("(per immune type shown separately)" if per_immune
           else "(sum over selected immune types)"),
        fontsize=12, labelpad=8,
    )
    if group_by == "organoid_type":
        title_detail = " + ".join(immune_types)
        title = (
            f"Cumulative interactions per organoid — contacts with {title_detail}"
        )
    else:
        title = (
            f"Cumulative interactions per organoid — by immune type "
            f"({' + '.join(immune_types)}; all organoid types pooled)"
        )
    if period_label:
        title += f"\n(analysis window: {period_label})"
    ax.set_title(
        title,
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
    annotate_line_condition: bool = False,
    lc_maps: dict = None,
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

    # Per-immune line-condition maps for disaggregation.
    lc_maps = lc_maps or {}

    def _lc_series(ct):
        m = {str(k): str(v) for k, v in lc_maps.get(ct, {}).items()}
        return df_dying["sample_name"].astype(str).map(m).fillna("(no LC)")

    # Decide whether to disaggregate lines by immune type x line condition.
    avail_immune = [ct for ct in immune_types if f"window_cum_{ct}" in df_dying.columns]
    disaggregate = False
    if annotate_line_condition:
        if len(avail_immune) > 1:
            disaggregate = True
        else:
            for ct in avail_immune:
                if len(set(str(v) for v in lc_maps.get(ct, {}).values())) > 1:
                    disaggregate = True
                    break

    # Build long-form data for plotting
    if disaggregate:
        # Separate line per (organoid_type ×) immune_type × line_condition.
        curve_frames = []
        for ct in avail_immune:
            cum_col = f"window_cum_{ct}"
            sub = df_dying[track_keys + ["relative_time_min", cum_col]].copy()
            sub = sub.rename(columns={cum_col: "cumulative_interactions"})
            lc = _lc_series(ct).values
            if group_by == "treatment":
                sub["condition"] = [f"{ct} | {v}" for v in lc]
            else:
                sub["condition"] = [
                    f"{ot} | {ct} | {v}"
                    for ot, v in zip(sub["organoid_type"].astype(str), lc)
                ]
            curve_frames.append(sub)
        if not curve_frames:
            return None
        df_curve = pd.concat(curve_frames, ignore_index=True)
    elif group_by == "treatment":
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
        cond_col = "organoid_type"
        df_dying["cumulative_interactions"] = df_dying[cum_cols].sum(axis=1)
        keep_cols = track_keys + ["relative_time_min", "cumulative_interactions"]
        if cond_col not in keep_cols:
            keep_cols = keep_cols + [cond_col]
        df_curve = df_dying[keep_cols].copy()
        df_curve["condition"] = df_curve[cond_col]

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
    if track_summary is not None and not track_summary.empty:
        ts = track_summary.copy()
        if "fate" in ts.columns:
            dying_mask = ts["fate"] == "Dying"
        else:
            dying_mask = ~ts["survives"].fillna(True)

        conditions_for_legend = sorted(df_interp["condition"].unique())

        if disaggregate:
            # Keys match composite labels like "13T | blue | teg_cd4".
            for cond in conditions_for_legend:
                parts = [p.strip() for p in str(cond).split(" | ")]
                if group_by == "treatment" and len(parts) >= 2:
                    ct, lc = parts[0], parts[-1]
                    lc_map = {str(k): str(v) for k, v in lc_maps.get(ct, {}).items()}
                    samples = {s for s, v in lc_map.items() if v == lc}
                    mask = ts["sample_name"].astype(str).isin(samples)
                elif len(parts) >= 3:
                    ot, ct, lc = parts[0], parts[1], parts[-1]
                    lc_map = {str(k): str(v) for k, v in lc_maps.get(ct, {}).items()}
                    samples = {s for s, v in lc_map.items() if v == lc}
                    mask = (
                        (ts["organoid_type"].astype(str) == ot)
                        & ts["sample_name"].astype(str).isin(samples)
                    )
                else:
                    mask = pd.Series(False, index=ts.index)
                total_per_cond[cond] = int(mask.sum())
                dying_per_cond[cond] = int((mask & dying_mask).sum())
        elif group_by == "treatment":
            for ct in immune_types:
                total_per_cond[ct] = len(ts)
                dying_per_cond[ct] = int(dying_mask.sum())
        else:
            legend_col = (
                "line_condition"
                if group_by == "line_condition" and "line_condition" in ts.columns
                else "organoid_type"
            )
            for cond_val in ts[legend_col].unique():
                cond_mask = ts[legend_col] == cond_val
                total_per_cond[cond_val] = int(cond_mask.sum())
                dying_per_cond[cond_val] = int((cond_mask & dying_mask).sum())

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

def _resolve_active_killing_csv_paths(
    ak_dir: Path,
    im_type: str,
    organoid_types: list = None,
) -> tuple:
    """
    Locate active-killing CSVs under ``active_killing/`` or a target subfolder.

    Active killing writes to ``active_killing/<target>/`` (or ``combined/`` when
    multiple targets are selected). Legacy runs used a flat ``active_killing/``
    layout. Returns ``(events_path, killing_path)``; either may be ``None``.
    """
    events_name = f"contact_events_{im_type}.csv"
    killing_name = f"active_killing_per_timepoint_{im_type}.csv"

    flat_events = ak_dir / events_name
    if flat_events.exists():
        flat_killing = ak_dir / killing_name
        return flat_events, flat_killing if flat_killing.exists() else None

    combined_events = ak_dir / "combined" / events_name
    if combined_events.exists():
        combined_killing = ak_dir / "combined" / killing_name
        return combined_events, (
            combined_killing if combined_killing.exists() else None
        )

    for org_type in organoid_types or []:
        cand_events = ak_dir / org_type / events_name
        if cand_events.exists():
            cand_killing = ak_dir / org_type / killing_name
            return cand_events, (
                cand_killing if cand_killing.exists() else None
            )

    event_candidates = sorted(
        ak_dir.glob(f"*/{events_name}"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if len(event_candidates) == 1:
        ev = event_candidates[0]
        k = ev.parent / killing_name
        return ev, k if k.exists() else None
    if len(event_candidates) > 1:
        ev = event_candidates[0]
        k = ev.parent / killing_name
        return ev, k if k.exists() else None

    return None, None


def _load_active_killing_data(
    output_dir: Path,
    immune_types: list,
    organoid_types: list,
    track_summary: pd.DataFrame = None,
    period_bounds: dict = None,
    line_condition_map: dict = None,
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
    period_bounds : dict, optional
        ``{sample_name: (t_start, t_end)}`` inclusive ``position_t`` bounds for
        the absolute analysis window. When given, contact events that *overlap*
        their sample's bounds are kept and ``contact_duration`` is recomputed as
        the number of timepoints inside the window only (not the full event
        length). Per-timepoint killing is counted only within those bounds.
        Samples not present in the dict have no timepoints in the window and
        are dropped. ``None`` (default) means no window restriction.
    line_condition_map : dict, optional
        ``{sample_name: line_condition}``. When given, the returned DataFrame
        gains a ``line_condition`` column (per sample) so the dashboard can
        group by line condition instead of organoid type.

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
        events_path, killing_path = _resolve_active_killing_csv_paths(
            ak_dir, im_type, organoid_types=organoid_types,
        )

        if events_path is None or not events_path.exists():
            continue

        df_ev = pd.read_csv(events_path)
        if df_ev.empty:
            continue
        df_ev["immune_type"] = im_type

        # Absolute analysis window: keep overlapping events and clip duration
        # to timepoints inside the per-sample bounds.
        if period_bounds is not None:
            df_ev = _clip_contact_events_to_period(df_ev, period_bounds)
            if df_ev.empty:
                continue

        # Per-target active-killing attribution from the per-timepoint CSV.
        #
        # ``advanced_timepoint_features.analyze_active_killing_per_timepoint``
        # writes ``targeted_track_id`` only on timepoints where
        # ``is_active_killing=True`` (it is the specific organoid whose
        # death signal increase triggered the classification; see
        # behav3d/features/advanced_timepoint_features.py line 437). So
        # each killing event is attributable to exactly one target.
        per_target_kill = None
        if killing_path is not None and killing_path.exists():
            df_tp = pd.read_csv(killing_path)
            required_cols = {
                "is_active_killing", "contact_event_id", "sample_name",
                "targeted_track_id",
            }
            if not df_tp.empty and required_cols.issubset(df_tp.columns):
                # Restrict killing timepoints to the analysis window so the
                # killing-efficiency numerator only counts kills inside it.
                if period_bounds is not None and "position_t" in df_tp.columns:
                    df_tp = df_tp[
                        df_tp.apply(
                            lambda r: _in_bounds(
                                period_bounds, r["sample_name"], r["position_t"],
                            ),
                            axis=1,
                        )
                    ].copy()
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

    # Attach immune line condition (per sample) when requested, so the
    # dashboard can group by line_condition instead of organoid type.
    if line_condition_map:
        df_all["line_condition"] = (
            df_all["sample_name"].astype(str)
            .map({str(k): v for k, v in line_condition_map.items()})
            .fillna("(no line condition)")
        )

    keep = [
        "contact_event_id", "sample_name", "immune_track_id", "immune_type",
        "target_track_id", "organoid_type", "line_condition", "fate",
        "contact_duration", "has_active_killing", "n_killing_tp",
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

# Shared palettes used by all multi-organoid interaction plots so the same
# semantic colour appears everywhere:
#   - fate: dead = red, live = green (full-strength for box edges, strips, bars)
#   - fate (light): for violin fills behind the per-fate dot clouds
#   - killing status: purple = active killing event, blue = non-killing event
_FATE_COLORS = {"Dying": "#E53935", "Live": "#43A047"}
_FATE_COLORS_LIGHT = {"Dying": "#EF9A9A", "Live": "#A5D6A7"}
_KILLING_COLORS = {
    "Active killing event": "#7E57C2",
    "Non-killing event": "#1976D2",
}


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
    group_by: str = "organoid_type",
    period_label: str = None,
    annotate_line_condition: bool = False,
    lc_maps: dict = None,
):
    """
    Active-killing dashboard: two panels derived from the per-event / per-
    timepoint active-killing CSVs written by ``advanced_timepoint_features``
    (``contact_events_{im}.csv`` and ``active_killing_per_timepoint_{im}.csv``).

    Panels
    ------
    Left  -- Contact Event Duration per (primary x fate); one point
             per immune-cell <-> organoid contact event, colored by whether
             that event actively killed the targeted organoid.
    Right -- Active Killing Efficiency per (primary x fate); for each
             (primary, fate) group, the % of its contact events that
             actively killed the targeted organoid.

    ``group_by`` mirrors the violin / before-death curve:
    ``"organoid_type"`` keeps organoid types separate on the x-axis;
    ``"treatment"`` pools organoid types and groups by immune type (and,
    when annotating, by line condition). ``period_label`` (when set) is
    added to the figure title to note the absolute analysis window.

    Returns ``None`` when no active-killing data is available -- the
    dashboard has nothing to show without it.
    """
    if df_contact_events.empty:
        return None

    lc_maps = lc_maps or {}
    df = df_contact_events.copy()

    # Per-row line condition from the event's own immune type (multi-immune safe).
    if annotate_line_condition and "immune_type" in df.columns and lc_maps:
        def _row_lc(row):
            m = {str(k): str(v) for k, v in lc_maps.get(row["immune_type"], {}).items()}
            return m.get(str(row["sample_name"]), None)
        mapped = df.apply(_row_lc, axis=1)
        if "line_condition" in df.columns:
            df["line_condition"] = mapped.fillna(df["line_condition"])
        else:
            df["line_condition"] = mapped.fillna("(no LC)")

    # Same disaggregation trigger as the violin plot.
    disaggregate = False
    if annotate_line_condition:
        n_immune = df["immune_type"].nunique() if "immune_type" in df.columns else 1
        n_lc = df["line_condition"].nunique() if "line_condition" in df.columns else 1
        if n_immune > 1 or n_lc > 1:
            disaggregate = True

    # Fold the x-axis grouping into ``organoid_type`` (the column the panels read).
    if disaggregate:
        if group_by == "treatment":
            parts = df["immune_type"].astype(str)
            if "line_condition" in df.columns:
                parts = parts + "\n" + df["line_condition"].astype(str)
            primary_label = "Immune / line condition"
        else:
            parts = df["organoid_type"].astype(str)
            if "immune_type" in df.columns and df["immune_type"].nunique() > 1:
                parts = parts + "\n" + df["immune_type"].astype(str)
            if "line_condition" in df.columns:
                parts = parts + "\n" + df["line_condition"].astype(str)
            primary_label = "Organoid / immune / line condition"
        df["organoid_type"] = parts
    elif group_by == "treatment":
        df["organoid_type"] = df["immune_type"].astype(str)
        primary_label = "Immune type"
    else:
        primary_label = "Organoid"

    org_types = sorted(df["organoid_type"].unique())
    if not org_types:
        return None

    # Dynamic group count drives figure width: 2 x N_present_org_types.
    # Panels are now stacked vertically (contact duration on top as the
    # main view, active-killing efficiency as a smaller sub-panel below),
    # so the figure no longer needs to stretch horizontally to fit two
    # side-by-side panels.
    n_groups = max(2 * len(org_types), 4)
    figwidth = max(11.0, 3.2 * n_groups + 3.0)

    fig = plt.figure(figsize=(figwidth, 12), dpi=120)
    fig.patch.set_facecolor("#fafafa")
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.35)
    axC = fig.add_subplot(gs[0, 0])
    axD = fig.add_subplot(gs[1, 0])

    _panel_contact_duration(
        axC, df, org_types, primary_label=primary_label,
        # Disaggregation already places immune x line condition on the x-axis,
        # so no marker-based annotation is needed here.
        annotate_line_condition=annotate_line_condition and not disaggregate,
        period_label=period_label,
    )
    _panel_killing_proportion(
        axD, df, org_types, immune_types,
        primary_label=primary_label,
    )

    suptitle = "Interaction -- Active Killing"
    if period_label:
        suptitle += f"  (analysis window: {period_label})"
    fig.suptitle(
        suptitle,
        fontsize=15, fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


# -- Dashboard sub-panels ---------------------------------------------------

def _panel_contact_duration(
    ax: plt.Axes,
    df_contact_events: pd.DataFrame,
    org_types: list,
    primary_label: str = "organoid",
    annotate_line_condition: bool = False,
    period_label: str = None,
):
    """
    Contact Event Duration panel.

    Each point = one immune-cell <-> organoid contact event (as written to
    ``contact_events_{im}.csv`` by the active-killing feature step). Points
    are grouped on the x-axis by (primary, fate of the contacted
    organoid) and colored by whether that event actively killed the
    targeted organoid. x-axis scales dynamically to ``2 * len(org_types)``
    slots; fate order is Dead then Live per group.
    """
    _style_dashboard_ax(ax)

    panel_title = f"Contact event duration per {primary_label} fate"

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

    # Per-group violin background colour follows fate (light red for Dead,
    # light green for Live -- pulled from the shared module-level palette
    # so the fate vocabulary stays consistent across plots). Box edges use
    # the full-strength fate colour from `_FATE_COLORS`.
    palette = {
        _group_label(o, f): _FATE_COLORS_LIGHT[f]
        for (o, f) in group_pairs
    }
    dark_palette = {
        _group_label(o, f): _FATE_COLORS[f]
        for (o, f) in group_pairs
    }

    common = dict(
        data=df, x="group", y="contact_duration",
        order=group_order, ax=ax,
    )

    sns.violinplot(
        **common, hue="group", palette=palette, legend=False,
        width=0.75,
        inner=None, linewidth=0.8, alpha=0.30,
        density_norm="count", cut=0,
    )
    sns.boxplot(
        **common, hue="group", palette=dark_palette, legend=False,
        width=0.28,
        fliersize=0,
        boxprops=dict(linewidth=1.2),
        whiskerprops=dict(linewidth=1.0),
        capprops=dict(linewidth=1.0),
        showcaps=True,
        medianprops=dict(linewidth=2.0),
    )

    # Hollow the centered box so neither point cloud is hidden behind a
    # solid rectangle; per-fate colour is preserved on the edge / whiskers
    # / median tick because we re-use the box's own facecolor as edgecolor.
    for patch in ax.patches:
        if isinstance(patch, mpatches.PathPatch):
            face = patch.get_facecolor()
            patch.set_edgecolor(face)
            patch.set_facecolor("none")

    # Manual scatter: active-killing points sit on the LEFT half of each
    # violin (purple), non-killing points on the RIGHT half (blue). One-
    # sided jitter ranges in (-0.22, -0.03) / (+0.03, +0.22) reserve the
    # centre for the box without overlapping it.
    kill_color = _KILLING_COLORS
    rng = np.random.default_rng(1)
    use_lc_markers = (
        annotate_line_condition
        and "line_condition" in df.columns
        and df["line_condition"].nunique() > 1
    )
    lc_marker_map = {}
    if use_lc_markers:
        marker_cycle = ["o", "s", "^", "D", "v", "P", "*", "X"]
        for i, lc in enumerate(sorted(df["line_condition"].dropna().unique())):
            lc_marker_map[lc] = marker_cycle[i % len(marker_cycle)]

    for g_label, sub in df.groupby("group"):
        if g_label not in group_to_idx or len(sub) == 0:
            continue
        xi = group_to_idx[g_label]
        active_mask = sub["killing_label"] == "Active killing event"
        active = sub[active_mask]
        non_active = sub[~active_mask]
        if len(active) > 0:
            x_left = xi + rng.uniform(-0.22, -0.03, size=len(active))
            if use_lc_markers:
                active_plot = active.copy()
                active_plot["_x"] = x_left
                for lc, chunk in active_plot.groupby("line_condition"):
                    mk = lc_marker_map.get(lc, "o")
                    ax.scatter(
                        chunk["_x"], chunk["contact_duration"],
                        s=10, c=kill_color["Active killing event"],
                        marker=mk,
                        alpha=0.55, edgecolors="white", linewidths=0.25, zorder=3,
                    )
            else:
                ax.scatter(
                    x_left, active["contact_duration"],
                    s=10, c=kill_color["Active killing event"],
                    alpha=0.55, edgecolors="white", linewidths=0.25, zorder=3,
                )
        if len(non_active) > 0:
            x_right = xi + rng.uniform(0.03, 0.22, size=len(non_active))
            if use_lc_markers:
                non_plot = non_active.copy()
                non_plot["_x"] = x_right
                for lc, chunk in non_plot.groupby("line_condition"):
                    mk = lc_marker_map.get(lc, "o")
                    ax.scatter(
                        chunk["_x"], chunk["contact_duration"],
                        s=10, c=kill_color["Non-killing event"],
                        marker=mk,
                        alpha=0.55, edgecolors="white", linewidths=0.25, zorder=3,
                    )
            else:
                ax.scatter(
                    x_right, non_active["contact_duration"],
                    s=10, c=kill_color["Non-killing event"],
                    alpha=0.55, edgecolors="white", linewidths=0.25, zorder=3,
                )

    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=kill_color["Active killing event"],
            markeredgecolor="white", markeredgewidth=0.5,
            label="Active killing event (left)",
        ),
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=kill_color["Non-killing event"],
            markeredgecolor="white", markeredgewidth=0.5,
            label="Non-killing event (right)",
        ),
    ]
    if use_lc_markers:
        for lc, mk in lc_marker_map.items():
            handles.append(
                plt.Line2D(
                    [0], [0], marker=mk, linestyle="", markersize=6,
                    markerfacecolor="#666666", markeredgecolor="white",
                    markeredgewidth=0.5, label=str(lc),
                )
            )
    ax.legend(
        handles=handles, fontsize=9,
        frameon=True, framealpha=0.85, edgecolor="#cccccc",
        title="Each point = one contact event",
        title_fontsize=9, loc="upper right",
    )

    _apply_org_fate_xticks(ax, group_pairs)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(
        f"{primary_label.capitalize()} / fate", fontsize=11, labelpad=8,
    )
    ax.set_ylabel(
        "Contact duration in window (timepoints)"
        if period_label
        else "Contact duration (timepoints)",
        fontsize=11, labelpad=8,
    )
    ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)


def _panel_killing_proportion(
    ax: plt.Axes,
    df_contact_events: pd.DataFrame,
    org_types: list,
    immune_types: list,
    primary_label: str = "organoid",
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

    panel_title = f"Active killing efficiency per {primary_label} fate"

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
    ax.set_xlabel(
        f"{primary_label.capitalize()} / fate", fontsize=11, labelpad=8,
    )
    ax.set_ylabel(
        "Contact events that actively killed (%)",
        fontsize=11, labelpad=8,
    )
    ax.set_title(panel_title, fontsize=13, fontweight="semibold", pad=10)
