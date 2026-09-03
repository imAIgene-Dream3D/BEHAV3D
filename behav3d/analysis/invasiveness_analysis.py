"""
Invasiveness Analysis for BEHAV3D

Analyzes invasiveness from the immune cell's point of view:
- Fraction of invasive immune cells over time (boolean invasiveness, >= 50% surface contact)
- Mean / median invasiveness percentage over time
- Per-movie summary collapsing each sample into one statistic (mean/median/max/AUC)

The analysis reads the *immune* cell type's filtered track features and compares
invasiveness against one or more organoid *targets* (e.g. "27T" vs "MDO"), where
each target has ``{target}_invasiveness`` (bool) and ``{target}_invasiveness_perc``
(float, 0-100) columns produced during feature extraction.

Outputs:
- PDF plots to analysis/{immune_cell_type}/invasiveness_analysis/
- CSV statistics per sample / per target

This module is fully self-contained and does not modify any existing pipeline
code. Boolean coercion is kept local (see ``coerce_boolean_columns``) so we do
not repeat ``== "True"`` logic elsewhere.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from behav3d.analysis.interaction_analysis import (
    add_line_condition_column,
    apply_analysis_period,
    apply_analysis_period_timepoints,
    format_period_label,
    format_period_label_timepoints,
    _normalize_period_t,
)


# Values (case-insensitive, stripped) that count as boolean True.
_TRUE_STRINGS = {"true", "1", "1.0", "yes", "t"}

# Summary statistics available for the per-movie collapse.
SUMMARY_STATS = ("mean", "median", "max", "auc")


def coerce_boolean_series(s: pd.Series) -> pd.Series:
    """Robustly coerce a Series to boolean."""
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0) != 0
    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .isin(_TRUE_STRINGS)
    )


def coerce_boolean_columns(df: pd.DataFrame, cols) -> pd.DataFrame:
    """Coerce the given columns of ``df`` to boolean in place and return df."""
    for col in cols:
        if col in df.columns:
            df[col] = coerce_boolean_series(df[col])
    return df


def detect_invasiveness_targets(df: pd.DataFrame) -> list:
    """
    Return the list of available invasiveness targets in ``df``.

    A target ``T`` is available if a ``{T}_invasiveness_perc`` column exists.
    The ``any_org`` pseudo-target is excluded — only real organoid targets are
    returned, selectable individually or in combination.
    """
    targets = []
    for col in df.columns:
        if col.endswith("_invasiveness_perc"):
            target = col[: -len("_invasiveness_perc")]
            if target and target != "any_org":
                targets.append(target)
    return sorted(set(targets))


def _bool_col(target: str) -> str:
    return f"{target}_invasiveness"


def _perc_col(target: str) -> str:
    return f"{target}_invasiveness_perc"


def _time_group_keys(df: pd.DataFrame) -> list:
    """Group keys for over-time aggregation, carrying immune_type / line
    condition when present so multiple immune types (and line conditions)
    stay separated instead of being pooled."""
    keys = []
    if "immune_type" in df.columns:
        keys.append("immune_type")
    keys.append("sample_name")
    if "line_condition" in df.columns:
        keys.append("line_condition")
    keys.append("position_t")
    return keys


def aggregate_fraction_over_time(df: pd.DataFrame, targets: list) -> pd.DataFrame:
    """Fraction of invasive immune cells per (immune_type, sample_name,
    line_condition, position_t) per target. immune_type / line_condition are
    only added to the grouping when those columns exist."""
    keys = _time_group_keys(df)
    frames = []
    for target in targets:
        bcol = _bool_col(target)
        if bcol not in df.columns:
            continue
        grouped = (
            df.groupby(keys)[bcol]
            .agg(fraction="mean", n_cells="count")
            .reset_index()
        )
        grouped["target"] = target
        frames.append(grouped)
    if not frames:
        return pd.DataFrame(
            columns=keys + ["target", "fraction", "n_cells"]
        )
    return pd.concat(frames, ignore_index=True)


def aggregate_perc_over_time(df: pd.DataFrame, targets: list) -> pd.DataFrame:
    """Mean & median invasiveness percentage per (immune_type, sample_name,
    line_condition, position_t) per target."""
    keys = _time_group_keys(df)
    frames = []
    for target in targets:
        pcol = _perc_col(target)
        if pcol not in df.columns:
            continue
        grouped = (
            df.groupby(keys)[pcol]
            .agg(mean_perc="mean", median_perc="median", n_cells="count")
            .reset_index()
        )
        grouped["target"] = target
        frames.append(grouped)
    if not frames:
        return pd.DataFrame(
            columns=keys + ["target", "mean_perc", "median_perc", "n_cells"]
        )
    return pd.concat(frames, ignore_index=True)


def _collapse_series(times: np.ndarray, values: np.ndarray, summary_stat: str) -> float:
    """Collapse a per-timepoint series into a single summary statistic."""
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.all(np.isnan(values)):
        return np.nan
    if summary_stat == "mean":
        return float(np.nanmean(values))
    if summary_stat == "median":
        return float(np.nanmedian(values))
    if summary_stat == "max":
        return float(np.nanmax(values))
    if summary_stat == "auc":
        order = np.argsort(times)
        t = np.asarray(times, dtype=float)[order]
        y = values[order]
        mask = ~np.isnan(y)
        t, y = t[mask], y[mask]
        if t.size < 2:
            return float(y[0]) if y.size else np.nan
        span = t[-1] - t[0]
        if span <= 0:
            return float(np.nanmean(y))
        return float(np.trapezoid(y, t) / span)
    raise ValueError(
        f"Unknown summary_stat '{summary_stat}'. Expected one of {SUMMARY_STATS}."
    )


def _perc_time_col_for_collapse(summary_stat: str) -> str:
    """Column from ``aggregate_perc_over_time`` to collapse for per-movie %.

    ``median`` uses the per-timepoint median across cells; all other stats use
    the per-timepoint mean across cells (then collapse over time).
    """
    return "median_perc" if summary_stat == "median" else "mean_perc"


def summarize_per_movie(
    df: pd.DataFrame,
    targets: list,
    summary_stat: str = "mean",
) -> pd.DataFrame:
    """Collapse each (immune_type, sample, line_condition, target) time series
    into a single summary value. Each (movie x immune type x target) stays a
    separate row so nothing is pooled across immune types or line conditions."""
    frac = aggregate_fraction_over_time(df, targets)
    perc = aggregate_perc_over_time(df, targets)
    perc_time_col = _perc_time_col_for_collapse(summary_stat)

    # Enumerate distinct (immune_type, sample, line_condition) combinations
    # present in the data so each becomes its own per-movie point.
    ctx_cols = [c for c in ("immune_type", "sample_name", "line_condition") if c in df.columns]
    ctx = df[ctx_cols].drop_duplicates()

    rows = []
    for target in targets:
        f_t = frac[frac["target"] == target]
        p_t = perc[perc["target"] == target]
        for _, ctx_row in ctx.iterrows():
            mask_f = pd.Series(True, index=f_t.index)
            mask_p = pd.Series(True, index=p_t.index)
            for c in ctx_cols:
                mask_f &= (f_t[c] == ctx_row[c])
                mask_p &= (p_t[c] == ctx_row[c])
            f_s = f_t[mask_f].sort_values("position_t")
            p_s = p_t[mask_p].sort_values("position_t")
            if f_s.empty and p_s.empty:
                continue
            fraction_summary = _collapse_series(
                f_s["position_t"].values, f_s["fraction"].values, summary_stat
            )
            perc_summary = _collapse_series(
                p_s["position_t"].values, p_s[perc_time_col].values, summary_stat
            )
            n_cells = int(f_s["n_cells"].max()) if not f_s.empty else int(p_s["n_cells"].max())
            row = {
                "sample_name": ctx_row["sample_name"],
                "target": target,
                "summary_stat": summary_stat,
                "fraction_summary": fraction_summary,
                "perc_summary": perc_summary,
                "n_cells": n_cells,
            }
            if "immune_type" in ctx_cols:
                row["immune_type"] = ctx_row["immune_type"]
            if "line_condition" in ctx_cols:
                row["line_condition"] = ctx_row["line_condition"]
            rows.append(row)
    return pd.DataFrame(rows)


def _samples_grid(n_samples: int):
    n_cols = min(3, max(1, n_samples))
    n_rows = (n_samples + n_cols - 1) // n_cols
    return n_rows, n_cols


def _series_setup(long_df: pd.DataFrame):
    """Build a per-row ``series`` label and colour palette.

    When several immune types and/or line conditions are present, each
    combination gets its own legend entry (e.g. ``blue · teg_cd4 · 13T``).
    Single-immune, single-condition plots keep the plain target label.
    """
    long_df = long_df.copy()
    multi_immune = (
        "immune_type" in long_df.columns and long_df["immune_type"].nunique() > 1
    )
    has_lc = (
        "line_condition" in long_df.columns
        and long_df["line_condition"].notna().any()
    )
    if has_lc:
        lc_vals = long_df["line_condition"].astype(str).str.strip()
        has_lc = not lc_vals.isin(["(no line condition)", "", "nan"]).all()

    def _series_label(row):
        parts = []
        if multi_immune and "immune_type" in row.index:
            parts.append(str(row["immune_type"]))
        if has_lc and "line_condition" in row.index:
            lc = str(row["line_condition"]).strip()
            if lc not in ("(no line condition)", "nan"):
                parts.append(lc)
        parts.append(str(row["target"]))
        return " \u00b7 ".join(parts) if len(parts) > 1 else parts[0]

    long_df["series"] = long_df.apply(_series_label, axis=1)
    series_order = sorted(long_df["series"].unique())
    palette = dict(zip(series_order, sns.color_palette("tab10", max(len(series_order), 1))))
    return long_df, series_order, palette, multi_immune


def _facet_title(df_s: pd.DataFrame, sample: str, multi_immune: bool = False) -> str:
    """Sample title. Line conditions are shown in the legend when multi-immune."""
    if multi_immune:
        return str(sample)
    if "line_condition" in df_s.columns:
        lcs = sorted(
            str(v) for v in df_s["line_condition"].dropna().unique()
            if str(v).strip() not in ("(no line condition)", "nan", "")
        )
        if len(lcs) == 1:
            return f"{sample}\n[{lcs[0]}]"
    return str(sample)


def _plot_over_time(
    long_df: pd.DataFrame,
    value_col: str,
    ylabel: str,
    title: str,
    period_label: str = None,
    ylim=None,
    legend_title: str = "Series",
    subtitle: str = None,
) -> plt.Figure:
    """Shared over-time renderer: one coloured line per series, one subplot
    per sample. Immune type and line condition appear in the legend label."""
    long_df, series_order, palette, multi_immune = _series_setup(long_df)
    samples = sorted(long_df["sample_name"].unique())
    n_rows, n_cols = _samples_grid(len(samples))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4.2 * n_rows), squeeze=False)
    axes = axes.flatten()

    for idx, sample in enumerate(samples):
        ax = axes[idx]
        df_s = long_df[long_df["sample_name"] == sample]
        for series in series_order:
            df_t = df_s[df_s["series"] == series].sort_values("position_t")
            if df_t.empty:
                continue
            ax.plot(
                df_t["position_t"], df_t[value_col],
                linewidth=2, color=palette[series], label=series,
            )
        ax.set_xlabel("Time point", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.set_title(_facet_title(df_s, sample, multi_immune=multi_immune), fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(title=legend_title, fontsize=8)

    for idx in range(len(samples), len(axes)):
        axes[idx].set_visible(False)

    full_title = title
    if subtitle:
        full_title += f"\n{subtitle}"
    if period_label and period_label != "full movie":
        full_title += f"\n(window: {period_label})"
    fig.suptitle(full_title, fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


def plot_fraction_over_time(
    frac_long: pd.DataFrame,
    immune_cell_type: str,
    period_label: str = None,
) -> plt.Figure:
    """Fraction of invasive cells over time; separate line per immune x target."""
    return _plot_over_time(
        frac_long, "fraction", "Fraction invasive cells",
        f"Fraction of invasive {immune_cell_type} cells over time",
        period_label=period_label, ylim=(0, 1),
        legend_title="Immune \u00b7 line condition \u00b7 target",
        subtitle="% of cells with \u226550% surface contact (boolean invasive)",
    )


def plot_perc_over_time(
    perc_long: pd.DataFrame,
    immune_cell_type: str,
    statistic: str = "mean",
    period_label: str = None,
) -> plt.Figure:
    """Invasiveness percentage over time; separate line per immune x target."""
    value_col = "mean_perc" if statistic == "mean" else "median_perc"
    if statistic == "mean":
        subtitle = (
            "Average contact % across all cells at each timepoint "
            "(including non-invasive zeros)"
        )
    else:
        subtitle = (
            "Typical cell contact % at each timepoint; stays at 0 when "
            "most cells are non-invasive"
        )
    return _plot_over_time(
        perc_long, value_col, f"{statistic.capitalize()} invasiveness (%)",
        f"{statistic.capitalize()} {immune_cell_type} invasiveness (%) over time",
        period_label=period_label,
        legend_title="Immune \u00b7 line condition \u00b7 target",
        subtitle=subtitle,
    )


def plot_per_movie_summary(
    summary_df: pd.DataFrame,
    immune_cell_type: str,
    summary_stat: str,
    value: str = "fraction",
    hue_col: str = None,
    split_line_condition: bool = False,
) -> plt.Figure:
    """Per-movie summary: one point per movie, separated by target and (when
    present) by immune type and line condition into distinct dodged boxes."""
    value_col = "fraction_summary" if value == "fraction" else "perc_summary"
    ylabel = (
        "Fraction invasive cells" if value == "fraction"
        else "Invasiveness (%)"
    )

    targets = sorted(summary_df["target"].unique())
    plot_df = summary_df.dropna(subset=[value_col]).copy()

    # Build the sub-group (hue) that separates immune type x line condition.
    hue_parts = []
    if "immune_type" in plot_df.columns and plot_df["immune_type"].nunique() > 1:
        hue_parts.append("immune_type")
    if (
        split_line_condition
        and "line_condition" in plot_df.columns
        and plot_df["line_condition"].nunique() > 1
    ):
        hue_parts.append("line_condition")
    # Backward-compat: honour an explicit hue_col request.
    if not hue_parts and hue_col and hue_col in plot_df.columns:
        hue_parts.append(hue_col)

    use_hue = bool(hue_parts)
    legend_title = " \u00b7 ".join(
        {"immune_type": "immune", "line_condition": "line condition"}.get(p, p)
        for p in hue_parts
    ) if use_hue else None
    if use_hue:
        plot_df["_grp"] = plot_df[hue_parts].astype(str).agg(" \u00b7 ".join, axis=1)
        n_hue = max(plot_df["_grp"].nunique(), 1)
        fig_w = max(7, 2.6 * len(targets) * min(n_hue, 4))
    else:
        fig_w = max(6, 2.2 * len(targets))
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    if not plot_df.empty:
        box_hue = "_grp" if use_hue else "target"
        hue_order = (
            sorted(plot_df["_grp"].dropna().unique()) if use_hue else None
        )
        plot_width = 0.8 if use_hue else 0.5
        common_kw = dict(
            data=plot_df, x="target", y=value_col, order=targets,
            hue=box_hue, hue_order=hue_order, palette="tab10",
            ax=ax, dodge=use_hue,
        )
        if use_hue:
            common_kw["width"] = plot_width
        sns.boxplot(
            **common_kw,
            legend=use_hue,
            gap=0.6 if use_hue else 0,
            fliersize=0,
        )
        for patch in ax.patches:
            if hasattr(patch, "set_facecolor"):
                fc = patch.get_facecolor()
                if len(fc) == 4:
                    r, g, b, _ = fc
                    patch.set_facecolor((r, g, b, 0.35))
        strip_kw = common_kw.copy()
        strip_kw.pop("width", None)
        sns.stripplot(
            **strip_kw,
            legend=False,
            size=7,
            # jitter must be 0 when dodging so strip x-offsets match boxplot
            # (seaborn 0.13); otherwise points sit left/right of their boxes.
            jitter=0 if use_hue else True,
            edgecolor="white", linewidth=0.5,
        )
        if use_hue:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(
                    handles[:len(hue_order)], labels[:len(hue_order)],
                    title=legend_title, fontsize=9,
                )

    stat_label = {
        "mean": "Mean", "median": "Median", "max": "Max",
        "auc": "Time-normalized AUC",
    }.get(summary_stat, summary_stat)

    ax.set_xlabel("Target", fontsize=12)
    ax.set_ylabel(f"{stat_label} {ylabel} (per movie)", fontsize=12)
    if value == "fraction":
        ax.set_ylim(0, 1)
    else:
        ax.set_ylim(bottom=0)

    n_movies = plot_df["sample_name"].nunique() if not plot_df.empty else 0
    n_points = len(plot_df)
    if use_hue and n_points != n_movies:
        n_line = (
            f"{n_points} points ({n_movies} movie"
            f"{'s' if n_movies != 1 else ''}; one per movie \u00d7 group)"
        )
    else:
        n_line = (
            f"{n_points} point{'s' if n_points != 1 else ''} "
            f"({n_movies} movie{'s' if n_movies != 1 else ''})"
        )

    perc_note = ""
    if value == "perc" and summary_stat == "median":
        perc_note = " — collapsed from median-% over time"
    elif value == "perc":
        perc_note = " — collapsed from mean-% over time"

    ax.set_title(
        f"Per-movie {stat_label.lower()} invasiveness of {immune_cell_type} "
        f"({'boolean fraction' if value == 'fraction' else 'percentage'}"
        f"{perc_note})\n"
        f"(one point = one movie \u00d7 group, {n_line})",
        fontsize=12,
    )
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# New violin-by-fate plots: cross-reference immune invasiveness with organoid
# death outcome.  Faceted per sample, user selects target(s).
# ---------------------------------------------------------------------------

def _load_organoid_fate(
    output_dir: Path, organoid_types: list,
) -> pd.DataFrame:
    """Load organoid fate (survives / dies) from each organoid type's filtered CSV.

    Returns a DataFrame with columns:
        sample_name, TrackID, organoid_type, organoid_died  (bool)
    """
    frames = []
    for ot in organoid_types:
        csv_path = (
            output_dir / "analysis" / ot / "track_features"
            / f"BEHAV3D_{ot}_combined_track_features_filtered.csv"
        )
        if not csv_path.exists():
            continue
        dfo = pd.read_csv(csv_path)
        if "dead" not in dfo.columns:
            continue
        dfo["dead"] = coerce_boolean_series(dfo["dead"])
        # A track dies if it ever has dead == True
        fate = (
            dfo.groupby(["sample_name", "TrackID"])["dead"]
            .any()
            .reset_index()
            .rename(columns={"dead": "organoid_died"})
        )
        fate["organoid_type"] = ot
        fate["TrackID"] = fate["TrackID"].map(_normalize_track_id)
        frames.append(fate)
    if not frames:
        return pd.DataFrame(
            columns=["sample_name", "TrackID", "organoid_type", "organoid_died"]
        )
    return pd.concat(frames, ignore_index=True)


def _find_touching_column(df: pd.DataFrame, target: str):
    """Return the column that lists the organoid TrackID(s) each immune cell
    touches for ``target``.

    Feature extraction writes this as ``touching_{target}s`` (plural of the
    organoid-type name, e.g. ``touching_13Ts``). Older/other layouts may use
    ``touching_{target}``. Returns the first match or ``None``.
    """
    for cand in (f"touching_{target}s", f"touching_{target}", f"touching_{target}es"):
        if cand in df.columns:
            return cand
    # Last resort: any column starting with touching_ and containing the target.
    for c in df.columns:
        if c.startswith("touching_") and target.lower() in c.lower():
            return c
    return None


def _normalize_track_id(val) -> str:
    """Coerce a touched-track id (possibly '28.0') to a clean string id."""
    s = str(val).strip()
    if s in ("", "nan", "None", "0", "0.0", "False", "false"):
        return ""
    try:
        return str(int(float(s)))
    except (TypeError, ValueError):
        return s


def aggregate_invasiveness_by_organoid_fate(
    df_immune: pd.DataFrame,
    df_organoid_fate: pd.DataFrame,
    target: str,
) -> pd.DataFrame:
    """Cross-reference immune invasiveness with organoid fate.

    For each organoid in ``df_organoid_fate``, find immune cells contacting
    it and compute:
      - **mean_invasiveness_perc**: mean of ``{target}_invasiveness_perc``
        across all immune-cell timepoints in contact.
      - **n_invasive_cells**: count of unique immune TrackIDs with
        ``{target}_invasiveness == True`` at any timepoint.

    Each row is ONE organoid (never pooled across organoids), so the violin
    shows one point per organoid — matching the per-organoid logic requested.

    Returns one row per organoid with columns:
        sample_name, organoid_TrackID, organoid_type, organoid_died,
        mean_invasiveness_perc, n_invasive_cells
    """
    perc_c = _perc_col(target)
    bool_c = _bool_col(target)
    contact_col = _find_touching_column(df_immune, target)

    # Need contact column to know which organoid each immune cell touches.
    if contact_col is None:
        return pd.DataFrame()  # Gracefully return empty if not linkable.

    # Filter immune data to rows where the cell is touching a specific
    # organoid TrackID of the given type. The touching column may hold a
    # single id or a comma-separated list, so we explode it.
    df_c = df_immune.copy()
    df_c["_touch_raw"] = df_c[contact_col].astype(str)
    df_c = df_c[
        df_c["_touch_raw"].notna()
        & ~df_c["_touch_raw"].str.strip().str.lower().isin(
            ["", "nan", "none", "0", "0.0", "false"]
        )
    ].copy()
    if df_c.empty:
        return pd.DataFrame()

    df_c["_touch_list"] = df_c["_touch_raw"].str.split(",")
    df_contact = df_c.explode("_touch_list")
    df_contact["organoid_TrackID"] = df_contact["_touch_list"].map(_normalize_track_id)
    df_contact = df_contact[df_contact["organoid_TrackID"] != ""].copy()
    if df_contact.empty:
        return pd.DataFrame()

    # Ensure numeric perc column
    df_contact[perc_c] = pd.to_numeric(df_contact[perc_c], errors="coerce").fillna(0)
    if bool_c in df_contact.columns:
        df_contact[bool_c] = coerce_boolean_series(df_contact[bool_c])
    else:
        df_contact[bool_c] = df_contact[perc_c] >= 50.0

    # Aggregate per organoid (and per immune type when several are present, so
    # e.g. cd4 vs cd8 invasion of the same organoid stay separate points).
    org_keys = ["sample_name", "organoid_TrackID"]
    if "immune_type" in df_contact.columns:
        org_keys = ["immune_type"] + org_keys
    grouped = df_contact.groupby(org_keys)
    agg = grouped.agg(
        mean_invasiveness_perc=(perc_c, "mean"),
    ).reset_index()

    # Count of unique invasive immune cells (invasiveness == True at any t)
    invasive_only = df_contact[df_contact[bool_c]]
    if not invasive_only.empty:
        inv_counts = (
            invasive_only.groupby(org_keys)["TrackID"]
            .nunique()
            .rename("n_invasive_cells")
            .reset_index()
        )
        agg = agg.merge(inv_counts, on=org_keys, how="left")
    else:
        agg["n_invasive_cells"] = 0
    agg["n_invasive_cells"] = agg["n_invasive_cells"].fillna(0).astype(int)

    # Merge with organoid fate. Restrict to the target organoid type first:
    # TrackID numbering can repeat across organoid types (e.g. 13T id 10 and
    # green id 10 are different organoids), so merging without this filter
    # would duplicate rows and assign conflicting fates.
    fate = df_organoid_fate.copy()
    if "organoid_type" in fate.columns:
        fate_t = fate[fate["organoid_type"] == target]
        if not fate_t.empty:
            fate = fate_t
    # Normalise the merge key on both sides to the same string form so the
    # join never fails on dtype mismatch (int vs object), regardless of how
    # the caller built the fate table.
    fate = fate.copy()
    fate["TrackID"] = fate["TrackID"].map(_normalize_track_id)
    agg["organoid_TrackID"] = agg["organoid_TrackID"].map(_normalize_track_id)
    fate = fate.drop_duplicates(subset=["sample_name", "TrackID"])
    agg = agg.merge(
        fate.rename(columns={"TrackID": "organoid_TrackID"}),
        on=["sample_name", "organoid_TrackID"],
        how="left",
    )
    agg["organoid_died"] = agg["organoid_died"].fillna(False)
    agg["target"] = target
    return agg


def _plot_fate_violin(
    df_by_fate: pd.DataFrame,
    y_col: str,
    ylabel: str,
    title: str,
    period_label: str = None,
    ymin: float = None,
):
    """Shared fate violin renderer, faceted per sample. When several immune
    types are present, separate violins are dodged by immune type within each
    died/survived group; otherwise fate colours (pink/teal) are used."""
    if df_by_fate.empty:
        return None
    samples = sorted(df_by_fate["sample_name"].unique())
    n = len(samples)
    multi_immune = (
        "immune_type" in df_by_fate.columns
        and df_by_fate["immune_type"].nunique() > 1
    )
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5), squeeze=False)
    axes = axes.flatten()
    fate_palette = {False: "#f8a19a", True: "#7ecfc0"}  # pink=survived, teal=died

    for i, sample in enumerate(samples):
        ax = axes[i]
        ds = df_by_fate[df_by_fate["sample_name"] == sample]
        if ds.empty:
            ax.set_title(sample)
            continue
        if multi_immune:
            im_order = sorted(ds["immune_type"].unique())
            common = dict(
                data=ds, x="organoid_died", y=y_col, order=[False, True],
                hue="immune_type", hue_order=im_order, palette="tab10",
            )
            sns.violinplot(**common, inner=None, alpha=0.3, ax=ax,
                           legend=(i == 0), dodge=True)
            sns.stripplot(**common, dodge=True, jitter=True, size=4,
                          alpha=0.8, ax=ax, legend=False)
            if i == 0:
                ax.legend(title="Immune", fontsize=8)
        else:
            common = dict(
                data=ds, x="organoid_died", y=y_col, order=[False, True],
                hue="organoid_died", palette=fate_palette,
            )
            sns.violinplot(**common, inner=None, alpha=0.3, ax=ax, legend=False)
            sns.stripplot(**common, dodge=False, jitter=True, size=5,
                          alpha=0.8, ax=ax, legend=False)
        ax.set_title(sample, fontsize=10)
        ax.set_xlabel("Organoid died")
        ax.set_ylabel(ylabel if i == 0 else "")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["FALSE", "TRUE"])
        if ymin is not None:
            ax.set_ylim(bottom=ymin)

    if period_label:
        title += f"\n({period_label})"
    fig.suptitle(title, fontsize=13, y=1.02)
    fig.tight_layout()
    return fig


def plot_invasiveness_by_fate_perc(
    df_by_fate: pd.DataFrame,
    target: str,
    immune_cell_type: str,
    period_label: str = None,
):
    """Violin+strip: mean invasiveness (%) by organoid fate, faceted per sample."""
    return _plot_fate_violin(
        df_by_fate, "mean_invasiveness_perc",
        "Mean invasiveness (%) in window",
        f"Mean invasiveness (%) of {immune_cell_type} by {target} fate",
        period_label=period_label,
        ymin=0,
    )


def plot_invasiveness_by_fate_count(
    df_by_fate: pd.DataFrame,
    target: str,
    immune_cell_type: str,
    period_label: str = None,
):
    """Violin+strip: count of invasive cells by organoid fate, faceted per sample."""
    return _plot_fate_violin(
        df_by_fate, "n_invasive_cells",
        f"Invasive {immune_cell_type} count",
        f"Invasive {immune_cell_type} count per {target} by fate",
        period_label=period_label,
        ymin=0,
    )


def _build_invasiveness_output_suffix(
    targets: list,
    summary_stat: str,
    analysis_period_t,
    group_by_line_condition: bool = False,
    immune_types: list = None,
) -> str:
    """Filename suffix so different user selections do not overwrite."""
    tgt = "-".join(sorted(targets)) if targets else "all"
    parts = []
    if immune_types and len(immune_types) > 1:
        parts.append("im-" + "-".join(sorted(immune_types)))
    parts += [tgt, summary_stat]
    start, end = _normalize_period_t(analysis_period_t)
    if start is None and end is None:
        parts.append("Tall")
    else:
        lo = "s" if start is None else str(start)
        hi = "e" if end is None else str(end)
        parts.append(f"T{lo}-{hi}")
    if group_by_line_condition:
        parts.append("annoLC")
    return "_" + "_".join(parts)


def run_invasiveness_analysis(
    output_dir: str,
    immune_cell_type=None,
    targets: list = None,
    df_tracks_path: str = None,
    summary_stat: str = "mean",
    sample_names: list = None,
    show_plots: bool = True,
    group_by_line_condition: bool = False,
    analysis_period_t: tuple = None,
    organoid_types: list = None,
    immune_cell_types: list = None,
    # Backward compat:
    analysis_period_min: tuple = None,
):
    """
    Run invasiveness analysis for one or more immune cell types against one or
    more organoid targets. Any combination of immune types, targets and line
    conditions is shown in the SAME figures but kept separated (distinct
    lines / boxes / violins) — nothing is pooled across those dimensions.

    Parameters
    ----------
    immune_cell_type : str, optional
        Single immune cell type (backward-compatible). Ignored when
        ``immune_cell_types`` is given.
    immune_cell_types : list, optional
        Multiple immune cell types to combine into one analysis. Each type's
        own filtered CSV is loaded and tagged with an ``immune_type`` column so
        plots can separate by it.
    group_by_line_condition : bool
        When True, results are additionally separated by immune line condition
        (``im_{immune}_line_condition``) — a distinct box per line condition in
        the per-movie summary and an annotation in the over-time facets.
    analysis_period_t : tuple, optional
        Temporal range ``(start_t, end_t)`` in raw timepoints. Replaces the
        old ``analysis_period_min`` (which is kept for backward compat).
    organoid_types : list, optional
        When provided, the organoid filtered CSVs are loaded to cross-reference
        immune invasiveness with organoid fate, producing the violin-by-fate
        plots.
    analysis_period_min : tuple, optional
        **Deprecated** — kept for backward compatibility.
    """
    # Resolve the immune type list (multi takes precedence over single).
    if immune_cell_types:
        immune_list = [str(c) for c in immune_cell_types if c]
    elif immune_cell_type:
        immune_list = [str(immune_cell_type)]
    else:
        raise ValueError("Provide immune_cell_type or immune_cell_types.")
    immune_list = list(dict.fromkeys(immune_list))  # de-dup, keep order
    multi_immune = len(immune_list) > 1
    immune_label = "-".join(immune_list)

    print(f"--------------- Performing {immune_label} Invasiveness Analysis ---------------")
    start_time = time.time()

    # Resolve analysis_period_t vs legacy analysis_period_min.
    if analysis_period_t is None and analysis_period_min is not None:
        analysis_period_t = analysis_period_min
        _use_minutes_fallback = True
    else:
        _use_minutes_fallback = False

    if summary_stat not in SUMMARY_STATS:
        raise ValueError(
            f"summary_stat '{summary_stat}' invalid. Expected one of {SUMMARY_STATS}."
        )

    output_dir = Path(output_dir)
    # Combined output dir when several immune types are analysed together; keep
    # the per-immune location for single-immune runs (backward compatible).
    if multi_immune:
        results_dir = output_dir / "analysis" / "invasiveness_analysis"
    else:
        results_dir = output_dir / "analysis" / immune_list[0] / "invasiveness_analysis"
    results_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load + combine each immune type's filtered track features.
    # ------------------------------------------------------------------
    print("   Loading track features...")
    frames = []
    for im in immune_list:
        if df_tracks_path is not None and not multi_immune:
            csv_path = Path(df_tracks_path)
        else:
            csv_path = (
                output_dir / "analysis" / im / "track_features"
                / f"BEHAV3D_{im}_combined_track_features_filtered.csv"
            )
        if not csv_path.exists():
            raise FileNotFoundError(
                f"Filtered track features not found at: {csv_path}\n"
                "Run Feature Extraction (with invasiveness) and Track Filtering first."
            )
        df_im = pd.read_csv(csv_path)
        df_im["TrackID"] = df_im["TrackID"].astype(str)
        df_im["immune_type"] = im
        # Per-immune line condition (its own im_{immune}_line_condition column).
        if group_by_line_condition:
            lc_src = f"im_{im}_line_condition"
            if lc_src in df_im.columns:
                df_im["line_condition"] = df_im[lc_src].astype(str)
            else:
                df_im, _src = add_line_condition_column(df_im, [im])
        frames.append(df_im)

    df = pd.concat(frames, ignore_index=True)
    df = df.sort_values(by=["immune_type", "sample_name", "TrackID", "position_t"])

    if group_by_line_condition and "line_condition" not in df.columns:
        print("   ⚠️ No line_condition column found — results will not split "
              "by line condition.")
    elif group_by_line_condition:
        df["line_condition"] = df["line_condition"].fillna("(no line condition)")

    if sample_names:
        df = df[df["sample_name"].isin(list(sample_names))].copy()
        if df.empty:
            raise ValueError(
                f"No rows left after filtering to samples {list(sample_names)}."
            )

    available = detect_invasiveness_targets(df)
    if not available:
        raise ValueError(
            "No invasiveness columns (e.g. '{target}_invasiveness_perc') found for "
            f"{immune_label}. Enable 'invasiveness' during feature extraction."
        )

    if targets is None:
        targets = available
    else:
        targets = list(targets)
        missing = [t for t in targets if _perc_col(t) not in df.columns]
        if missing:
            print(f"   ⚠️ Skipping targets without invasiveness columns: {missing}")
        targets = [t for t in targets if _perc_col(t) in df.columns]
        if not targets:
            raise ValueError(
                f"None of the requested targets have invasiveness columns. "
                f"Available: {available}"
            )

    print(f"   Immune types: {immune_list}")
    print(f"   Targets: {targets}")

    out_suffix = _build_invasiveness_output_suffix(
        targets, summary_stat, analysis_period_t, group_by_line_condition,
        immune_types=immune_list,
    )

    coerce_boolean_columns(df, [_bool_col(t) for t in targets])

    n_cells = df.groupby(["immune_type", "sample_name"])["TrackID"].nunique().sum()
    n_samples = df["sample_name"].nunique()
    print(f"   Loaded {len(df)} timepoints from {n_cells} {immune_label} cells across {n_samples} samples")

    # Apply temporal range.
    period_label = None
    if _use_minutes_fallback and analysis_period_min is not None:
        df = apply_analysis_period(df, analysis_period_min)
        period_label = format_period_label(analysis_period_min)
        if df.empty:
            raise ValueError(
                f"No timepoints left after applying analysis window ({period_label})."
            )
        print(f"   Analysis window: {period_label} ({len(df)} timepoints retained)")
    elif analysis_period_t is not None:
        start_t, end_t = _normalize_period_t(analysis_period_t)
        if not (start_t is None and end_t is None):
            df = apply_analysis_period_timepoints(df, analysis_period_t)
            period_label = format_period_label_timepoints(analysis_period_t)
            if df.empty:
                raise ValueError(
                    f"No timepoints left after applying temporal range ({period_label})."
                )
            print(f"   Temporal range: {period_label} ({len(df)} timepoints retained)")

    if group_by_line_condition and "line_condition" in df.columns:
        n_lc = df["line_condition"].nunique()
        print(f"   Separating by immune line condition ({n_lc} condition(s))")

    frac_long = aggregate_fraction_over_time(df, targets)
    perc_long = aggregate_perc_over_time(df, targets)
    summary_df = summarize_per_movie(df, targets, summary_stat=summary_stat)

    frac_csv = results_dir / f"invasiveness_fraction_over_time_{immune_label}{out_suffix}.csv"
    perc_csv = results_dir / f"invasiveness_perc_over_time_{immune_label}{out_suffix}.csv"
    summary_csv = results_dir / f"invasiveness_per_movie_summary_{immune_label}{out_suffix}.csv"
    frac_long.to_csv(frac_csv, index=False)
    perc_long.to_csv(perc_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)
    print(f"   Statistics saved: {frac_csv.name}, {perc_csv.name}, {summary_csv.name}")

    figs = {}
    if not frac_long.empty:
        figs["fraction_over_time"] = plot_fraction_over_time(
            frac_long, immune_label, period_label=period_label,
        )
    if not perc_long.empty:
        figs["perc_over_time_mean"] = plot_perc_over_time(
            perc_long, immune_label, statistic="mean", period_label=period_label,
        )
        figs["perc_over_time_median"] = plot_perc_over_time(
            perc_long, immune_label, statistic="median", period_label=period_label,
        )
    if not summary_df.empty:
        figs["per_movie_fraction"] = plot_per_movie_summary(
            summary_df, immune_label, summary_stat,
            value="fraction", split_line_condition=group_by_line_condition,
        )
        figs["per_movie_perc"] = plot_per_movie_summary(
            summary_df, immune_label, summary_stat,
            value="perc", split_line_condition=group_by_line_condition,
        )

    # ------------------------------------------------------------------
    # NEW: Violin-by-fate plots (cross-referencing with organoid death)
    # ------------------------------------------------------------------
    fate_plot_keys = []
    if organoid_types:
        print("   Loading organoid fate data for violin-by-fate plots...")
        df_organoid_fate = _load_organoid_fate(output_dir, organoid_types)
        if not df_organoid_fate.empty:
            for target in targets:
                df_by_fate = aggregate_invasiveness_by_organoid_fate(
                    df, df_organoid_fate, target,
                )
                if df_by_fate.empty:
                    print(f"   ⚠️ No contact linkage for target '{target}' — "
                          "skipping violin-by-fate plots (need touching_{target} column).")
                    continue
                key_perc = f"fate_perc_{target}"
                key_count = f"fate_count_{target}"
                fig_perc = plot_invasiveness_by_fate_perc(
                    df_by_fate, target, immune_label,
                    period_label=period_label,
                )
                fig_count = plot_invasiveness_by_fate_count(
                    df_by_fate, target, immune_label,
                    period_label=period_label,
                )
                if fig_perc is not None:
                    figs[key_perc] = fig_perc
                    fate_plot_keys.append(key_perc)
                if fig_count is not None:
                    figs[key_count] = fig_count
                    fate_plot_keys.append(key_count)
                # Save fate CSV
                fate_csv = results_dir / f"invasiveness_by_fate_{target}_{immune_label}{out_suffix}.csv"
                df_by_fate.to_csv(fate_csv, index=False)
                print(f"   Fate data saved: {fate_csv.name}")
        else:
            print("   ⚠️ No organoid fate data available — skipping violin-by-fate plots.")

    # Build ordered list of plot keys for the PDF.
    pdf_keys = [
        "fraction_over_time", "perc_over_time_mean", "perc_over_time_median",
        "per_movie_fraction", "per_movie_perc",
    ] + fate_plot_keys

    pdf_path = results_dir / f"invasiveness_analysis_{immune_label}{out_suffix}.pdf"
    with PdfPages(pdf_path) as pdf:
        for key in pdf_keys:
            if key in figs:
                pdf.savefig(figs[key], bbox_inches="tight")
    print(f"   📈 PDF saved: {pdf_path.name}")

    if show_plots:
        from IPython.display import display
        for fig in figs.values():
            display(fig)

    for fig in figs.values():
        plt.close(fig)

    elapsed = time.time() - start_time
    print(f"\n   ✅ Invasiveness Analysis complete! ({elapsed:.1f}s)")
    print(f"   Results saved to: {results_dir}")

    return {
        "fraction_over_time": frac_long,
        "perc_over_time": perc_long,
        "per_movie_summary": summary_df,
        "targets": targets,
        "pdf_path": pdf_path,
        "fraction_csv_path": frac_csv,
        "perc_csv_path": perc_csv,
        "summary_csv_path": summary_csv,
    }
