from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    _build_state_color_map,
)
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    compute_class_by_stack_proportions,
    draw_stacked_proportion_barv,
    draw_diff_barh,
    _welch_diff_rows,
)
from behav3d.analysis.behavior.track.contact_state_shift import (
    compute_track_state_shift_features,
    summarize_state_shift_track_fractions,
    BEFORE_LABEL,
    AFTER_LABEL,
)

_CONTACT_GROUPS = ("contact", "no_contact")
_CONTACT_GROUP_TITLES = {"contact": "Contact tracks", "no_contact": "No-contact tracks (null)"}


def _plot_contact_state_shift_page(
    state_timepoints_df,
    track_fraction_df,
    *,
    state_col,
    state_order,
    colors,
    params_text,
):
    """Build the combined 2x2 report figure: rows = {diff-bars, stacked composition}, columns =
    {contact tracks, no-contact tracks (null)}.

    ``track_fraction_df`` is a flat (non-indexed) DataFrame with columns ``["contact_group",
    "period", *state_order]`` — one row per (track, period).

    Returns ``(fig, long_rows)`` where ``long_rows`` is a list of flat dicts (one per state x
    group x panel-type) for the combined CSV export.
    """
    fig = plt.figure(figsize=(11.0, 9.5))
    outer = GridSpec(nrows=3, ncols=2, height_ratios=[0.35, 1.0, 1.0], hspace=0.55, wspace=0.35)

    ax_header = fig.add_subplot(outer[0, :])
    ax_header.axis("off")
    ax_header.text(0.0, 0.5, params_text, ha="left", va="center", fontsize=8, family="monospace")

    long_rows = []

    for c, group in enumerate(_CONTACT_GROUPS):
        ax = fig.add_subplot(outer[1, c])
        group_fractions = track_fraction_df[track_fraction_df["contact_group"] == group]
        before_df = group_fractions[group_fractions["period"] == BEFORE_LABEL][state_order]
        after_df = group_fractions[group_fractions["period"] == AFTER_LABEL][state_order]

        if len(before_df) > 0 and len(after_df) > 0:
            diff_rows = _welch_diff_rows(state_order, before_df, after_df)
            diff_df = pd.DataFrame(diff_rows)
            draw_diff_barh(ax, state_order, diff_df, colors)
            for row in diff_rows:
                long_rows.append({
                    "panel": "diff_bar", "contact_group": group, "state": row["class"],
                    "mean_before": row["mean_a"], "mean_after": row["mean_b"],
                    "diff": row["diff"], "t_stat": row["t_stat"], "p_value": row["p_value"],
                    "stars": row["stars"], "n_before": row["n_a"], "n_after": row["n_b"],
                })
        else:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
        ax.set_title(f"{_CONTACT_GROUP_TITLES[group]}\nBefore → After state change", fontsize=10)

        ax2 = fig.add_subplot(outer[2, c])
        state_tp = state_timepoints_df[state_timepoints_df["contact_group"] == group]
        if len(state_tp) > 0:
            props_df = compute_class_by_stack_proportions(
                state_tp, class_col="period", stack_col=state_col,
                class_order=[BEFORE_LABEL, AFTER_LABEL], stack_order=state_order,
            )
            draw_stacked_proportion_barv(ax2, props_df, [BEFORE_LABEL, AFTER_LABEL], state_order, colors)
            for period in (BEFORE_LABEL, AFTER_LABEL):
                if period in props_df.index:
                    for state in state_order:
                        long_rows.append({
                            "panel": "stacked_composition", "contact_group": group,
                            "period": period, "state": state,
                            "proportion": float(props_df.loc[period, state]),
                        })
        else:
            ax2.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax2.transAxes)
            ax2.axis("off")
        ax2.set_title(f"{_CONTACT_GROUP_TITLES[group]}\nState composition", fontsize=10)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[s]) for s in state_order]
    fig.legend(handles, state_order, loc="lower center", ncol=min(len(state_order), 6), frameon=False, fontsize=8)
    fig.suptitle("Contact-triggered behavioral-state shift (before vs. after)", fontsize=13, fontweight="bold")
    return fig, long_rows


def save_track_contact_state_shift_report(
    adata_tracks,
    df_timepoints,
    adata_states,
    out_dir,
    *,
    contact_col,
    min_bout_length,
    state_col=FULL_STATE_COL,
    window_mode="fixed",
    fixed_window_length=10,
    min_window_timepoints=3,
    sample_col="sample_name",
    track_col="TrackID",
    state_order=None,
    state_colors=None,
    null_seed=0,
    verbose=False,
):
    """Compare each track's behavioral-state composition before vs. after its first sufficiently
    long contact bout (contact tracks), against a timing-matched null before/after split for
    no-contact tracks. Writes a single combined PDF (2x2 grid: diff-bars / stacked composition x
    contact / no-contact) plus CSVs, into the same ``{out_dir}/contact_analysis/{contact_col}/``
    folder used by ``save_track_contact_group_analysis`` (as sibling artifacts, not appended
    pages, so this analysis can be re-run independently with its own parameters).

    Returns a dict of artifact paths plus ``n_contact_tracks``/``n_no_contact_tracks``/
    ``n_excluded_tracks``.
    """
    groupby_cols = [str(sample_col), str(track_col)]

    if state_order is None or state_colors is None:
        resolved_order, resolved_colors = _build_state_color_map(adata_states, state_col)
        state_order = state_order if state_order is not None else resolved_order
        state_colors = state_colors if state_colors is not None else resolved_colors

    features = compute_track_state_shift_features(
        df_timepoints,
        adata_tracks,
        adata_states,
        contact_col=contact_col,
        min_bout_length=min_bout_length,
        state_col=state_col,
        window_mode=window_mode,
        fixed_window_length=fixed_window_length,
        min_window_timepoints=min_window_timepoints,
        groupby_cols=groupby_cols,
        null_seed=null_seed,
        verbose=verbose,
    )
    track_windows = features["track_windows"]
    state_timepoints = features["state_timepoints"]

    fraction_df = summarize_state_shift_track_fractions(
        state_timepoints, groupby_cols=groupby_cols, state_col=state_col, state_order=state_order,
    ).reset_index()
    # Attach contact_group (constant per track) as a plain column so the plotting step can filter
    # by it directly, without juggling a mixed-depth MultiIndex.
    contact_group_by_track = state_timepoints.drop_duplicates(subset=groupby_cols)[groupby_cols + ["contact_group"]]
    fraction_df = fraction_df.merge(contact_group_by_track, on=groupby_cols, how="left")

    out_dir = Path(out_dir) / "contact_analysis" / str(contact_col)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = out_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "contact_state_shift.pdf"

    params_text = (
        f"contact_col={contact_col}  min_bout_length={min_bout_length}  "
        f"state_col={state_col}  window_mode={window_mode}"
        + (f"  fixed_window_length={fixed_window_length}" if window_mode == "fixed" else "")
        + f"  min_window_timepoints={min_window_timepoints}  null_seed={null_seed}\n"
        f"n_contact_tracks={features['n_contact_tracks']}  "
        f"n_no_contact_tracks={features['n_no_contact_tracks']}  "
        f"n_excluded_tracks={features['n_excluded_tracks']}"
    )

    fig, long_rows = _plot_contact_state_shift_page(
        state_timepoints, fraction_df,
        state_col=state_col, state_order=state_order, colors=state_colors, params_text=params_text,
    )
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    windows_csv = csv_dir / "state_shift_track_windows.csv"
    track_windows.reset_index().to_csv(windows_csv, index=False)

    long_df = pd.DataFrame(long_rows)
    diff_csv = csv_dir / "state_shift_diff_bars.csv"
    stacked_csv = csv_dir / "state_shift_stacked_composition.csv"
    if len(long_df) > 0:
        long_df[long_df["panel"] == "diff_bar"].drop(columns=["panel"]).to_csv(diff_csv, index=False)
        long_df[long_df["panel"] == "stacked_composition"].drop(columns=["panel"]).to_csv(stacked_csv, index=False)
    else:
        pd.DataFrame().to_csv(diff_csv, index=False)
        pd.DataFrame().to_csv(stacked_csv, index=False)

    if verbose:
        print(f"Saved contact state-shift report: {pdf_path}")

    return {
        "contact_col": str(contact_col),
        "min_bout_length": int(min_bout_length),
        "window_mode": str(window_mode),
        "pdf_path": str(pdf_path),
        "csv_dir": str(csv_dir),
        "track_windows_csv": str(windows_csv),
        "diff_bars_csv": str(diff_csv),
        "stacked_composition_csv": str(stacked_csv),
        "n_contact_tracks": features["n_contact_tracks"],
        "n_no_contact_tracks": features["n_no_contact_tracks"],
        "n_excluded_tracks": features["n_excluded_tracks"],
    }
