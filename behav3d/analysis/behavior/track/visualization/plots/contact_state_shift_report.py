from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from behav3d.features.state_descriptive_features import rle_encode
from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.contact_grouping import (
    compute_track_contact_features,
    _contact_group_col_name,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    _build_state_color_map,
    _extract_track_window_obs,
    _compute_state_bar_segments,
    _plot_statebar_segments_on_ax,
)
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    compute_class_by_stack_proportions,
    draw_stacked_proportion_barv,
    draw_diff_barh,
    SIGNIFICANCE_LEGEND_TEXT,
    _welch_diff_rows,
    _chunk_list,
)
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    compute_relative_state_time_matrix,
)
from behav3d.analysis.behavior.track.contact_state_shift import (
    compute_track_state_shift_features,
    summarize_state_shift_track_fractions,
    BEFORE_LABEL,
    AFTER_LABEL,
)

_CONTACT_GROUPS = ("contact", "no_contact")
_CONTACT_GROUP_TITLES = {"contact": "Contact tracks", "no_contact": "No-contact tracks (null)"}


def _plot_state_composition_over_relative_time(ax, panel_df, *, state_col, state_order, colors, fixed_window_length):
    """Continuous stacked composition vs. time-relative-to-contact (`relative_t`), x-axis capped
    to `[-fixed_window_length, fixed_window_length]` — the "state composition over time" view,
    windowed around the contact bout instead of spanning the whole track."""
    if len(panel_df) == 0:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return
    mat = compute_relative_state_time_matrix(
        panel_df, time_col="relative_t", state_col=state_col, state_order=state_order,
    )
    x = mat.index.to_numpy(dtype=float)
    bottom = np.zeros(len(x))
    for state in reversed(list(state_order)):
        v = mat[state].to_numpy(dtype=float)
        ax.bar(x, v, bottom=bottom, width=1.0, align="edge", linewidth=0, color=colors.get(state), label=str(state))
        bottom += v
    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_xlim(-float(fixed_window_length), float(fixed_window_length))
    ax.set_ylim(0.0, 1.0)
    ax.margins(x=0)
    ax.set_xlabel("Time relative to contact (timepoints)", fontsize=8)
    ax.set_ylabel("Proportion", fontsize=8)


def _plot_contact_state_shift_page(
    state_timepoints_df,
    track_fraction_df,
    *,
    state_col,
    state_order,
    colors,
    params_text,
    fixed_window_length,
):
    """Build the combined report figure: rows = {composition over time (capped), diff-bars,
    stacked before/after composition}, columns = {contact tracks, no-contact tracks (null)}.

    ``track_fraction_df`` is a flat (non-indexed) DataFrame with columns ``["contact_group",
    "period", *state_order]`` — one row per (track, period).

    Returns ``(fig, long_rows)`` where ``long_rows`` is a list of flat dicts (one per state x
    group x panel-type) for the combined CSV export.
    """
    fig = plt.figure(figsize=(11.0, 12.5))
    outer = GridSpec(nrows=4, ncols=2, height_ratios=[0.35, 1.0, 1.0, 1.0], hspace=0.6, wspace=0.35)

    ax_header = fig.add_subplot(outer[0, :])
    ax_header.axis("off")
    ax_header.text(0.0, 0.5, params_text, ha="left", va="center", fontsize=8, family="monospace")
    ax_header.text(1.0, 0.5, SIGNIFICANCE_LEGEND_TEXT, ha="right", va="center", fontsize=8)

    long_rows = []

    for c, group in enumerate(_CONTACT_GROUPS):
        ax_time = fig.add_subplot(outer[1, c])
        group_state_tp = state_timepoints_df[state_timepoints_df["contact_group"] == group]
        _plot_state_composition_over_relative_time(
            ax_time, group_state_tp, state_col=state_col, state_order=state_order, colors=colors,
            fixed_window_length=fixed_window_length,
        )
        ax_time.set_title(
            f"{_CONTACT_GROUP_TITLES[group]}\nState composition over time (window={fixed_window_length})",
            fontsize=10,
        )

        ax = fig.add_subplot(outer[2, c])
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

        ax2 = fig.add_subplot(outer[3, c])
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

    legend_state_order = list(reversed(state_order))
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[s]) for s in legend_state_order]
    fig.legend(handles, legend_state_order, loc="lower center", ncol=min(len(state_order), 6), frameon=False, fontsize=8)
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
        fixed_window_length=fixed_window_length,
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


_TRACK_OVERVIEW_WINDOW_COL = "trajectory_window_id"
_CONTACT_BAR_COLOR = "#2ca02c"
_NO_CONTACT_BAR_COLOR = "#bbbbbb"


def _compute_contact_bar_segments(times, is_contact, min_bout_length):
    """RLE-encode a track's per-timepoint contact boolean series into colored bar segments,
    time-aligned with `_compute_state_bar_segments` (same width convention). A run is green only
    if it's truthy AND long enough to meet `min_bout_length` -- the same threshold
    `contact_grouping.compute_track_contact_features` uses for "contact" grouping -- so a
    too-short contact blip renders grey, same as no contact at all.
    """
    times = np.asarray(times, dtype=float)
    is_contact = np.asarray(is_contact, dtype=bool)
    if len(times) == 0:
        raise ValueError("Cannot render contact bar: no timepoints found.")

    dx = np.diff(times)
    default_w = np.median(dx[dx > 0]) if np.any(dx > 0) else 1.0
    widths = np.r_[dx, default_w]

    min_bout_length = int(min_bout_length)
    segments = []
    pos = 0
    for value, length in rle_encode(is_contact.tolist()):
        start = float(times[pos])
        end = float(times[pos + length - 1] + widths[pos + length - 1])
        color = _CONTACT_BAR_COLOR if (bool(value) and length >= min_bout_length) else _NO_CONTACT_BAR_COLOR
        segments.append((start, max(0.0, end - start), color))
        pos += length
    return segments


def _plot_track_contact_overview_page(page_rows, *, sample_name, contact_col, min_bout_length):
    """One PDF page: each track gets its state-over-time bar (top) with the matching contact bar
    (bottom, grey/green) directly underneath. All rows share one page-wide time axis (rather than
    each track being independently stretched to the same width), so a track's bar width in the
    page reflects its actual duration relative to the other tracks on that page.
    """
    n = len(page_rows)
    page_xlim = (
        min(r["xlim"][0] for r in page_rows),
        max(r["xlim"][1] for r in page_rows),
    )
    fig = plt.figure(figsize=(11.0, max(2.4, 1.7 * n)))
    outer = fig.add_gridspec(nrows=n, ncols=1, hspace=0.75, top=0.90, bottom=0.10)

    for i, row_data in enumerate(page_rows):
        inner = outer[i].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.1)
        ax_state = fig.add_subplot(inner[0])
        ax_contact = fig.add_subplot(inner[1], sharex=ax_state)

        _plot_statebar_segments_on_ax(
            ax_state,
            segments=row_data["state_segments"],
            xlim=page_xlim,
            title=f"{row_data['sample_name']} | Track {row_data['track_id']}",
        )
        ax_state.set_xlabel("")
        ax_state.tick_params(labelbottom=False)

        for start, width, color in row_data["contact_segments"]:
            ax_contact.broken_barh([(start, width)], (0.0, 1.0), facecolors=color)
        ax_contact.set_xlim(*page_xlim)
        ax_contact.set_ylim(0.0, 1.0)
        ax_contact.set_yticks([])
        ax_contact.set_xlabel("position_t", fontsize=8)

    legend_handles = [
        Patch(facecolor=_CONTACT_BAR_COLOR, edgecolor="k", label=f"contact (bout >= {min_bout_length})"),
        Patch(facecolor=_NO_CONTACT_BAR_COLOR, edgecolor="k", label="no contact"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle(
        f"Track contact overview | sample: {sample_name} | contact_col={contact_col}",
        fontsize=11, fontweight="bold",
    )
    return fig


def save_track_contact_overview_report(
    adata_tracks,
    df_timepoints,
    adata_states,
    out_dir,
    *,
    contact_col,
    min_bout_length,
    state_col=FULL_STATE_COL,
    sample_col="sample_name",
    track_col="TrackID",
    rows_per_page=6,
    state_order=None,
    state_colors=None,
    plot_dpi=200,
    verbose=False,
):
    """For every track whose contact meets `min_bout_length` (the same threshold used by
    `contact_grouping.compute_track_contact_features` for "contact" vs. "no_contact" grouping),
    plot its full (untrimmed, classified-window) behavioral-state trajectory as a colored bar,
    with a grey/green bar directly beneath it marking every `contact_col` bout of at least
    `min_bout_length` timepoints. Pages are grouped by sample -- a sample's tracks are never
    split across a page shared with the next sample's, even if that leaves the page under-full.

    Writes into the same `{out_dir}/contact_analysis/{contact_col}/` folder used by the other two
    contact reports, as a sibling PDF (`track_contact_overview.pdf`).

    Returns a dict with `pdf_path`, `n_tracks`, `n_samples`.
    """
    sample_col = str(sample_col)
    track_col = str(track_col)
    groupby_cols = [sample_col, track_col]
    time_col = "position_t"

    if state_order is None or state_colors is None:
        resolved_order, resolved_colors = _build_state_color_map(adata_states, state_col)
        state_order = state_order if state_order is not None else resolved_order
        state_colors = state_colors if state_colors is not None else resolved_colors

    out_dir = Path(out_dir) / "contact_analysis" / str(contact_col)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "track_contact_overview.pdf"

    contact_features = compute_track_contact_features(
        df_timepoints, adata_tracks,
        contact_col=contact_col, min_bout_length=min_bout_length,
        groupby_cols=groupby_cols, verbose=verbose,
    )
    group_col = _contact_group_col_name(contact_col)
    contact_tracks = contact_features[contact_features[group_col] == "contact"].reset_index()

    has_window_col = (
        _TRACK_OVERVIEW_WINDOW_COL in adata_tracks.obs.columns
        and _TRACK_OVERVIEW_WINDOW_COL in contact_tracks.columns
    )
    key_cols = groupby_cols + ([_TRACK_OVERVIEW_WINDOW_COL] if has_window_col else [])

    windows = (
        adata_tracks.obs[key_cols + ["position_t_min", "position_t_max"]]
        .drop_duplicates(subset=key_cols)
        .copy()
    )
    for col in key_cols:
        contact_tracks[col] = contact_tracks[col].astype(str)
        windows[col] = windows[col].astype(str)
    tracks_df = contact_tracks.merge(windows, on=key_cols, how="left")

    if len(tracks_df) == 0:
        fig, ax = plt.subplots(figsize=(10, 2.2))
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            f"No tracks met the contact threshold (min_bout_length={min_bout_length}) for '{contact_col}'.",
            ha="center", va="center", fontsize=10,
        )
        with PdfPages(pdf_path) as pdf:
            pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)
        return {"pdf_path": str(pdf_path), "n_tracks": 0, "n_samples": 0}

    prepared_rows = []
    sample_order = []
    for _, row in tracks_df.iterrows():
        sample_name = row[sample_col]
        track_id = row[track_col]
        tmin = row["position_t_min"]
        tmax = row["position_t_max"]
        if sample_name not in sample_order:
            sample_order.append(sample_name)

        state_track_df = _extract_track_window_obs(
            adata_states,
            sample_name=sample_name,
            track_id=track_id,
            tmin=tmin,
            tmax=tmax,
            sample_key=sample_col,
            track_key=track_col,
            time_key=time_col,
            extra_cols=[state_col],
        )
        state_segments, xlim = _compute_state_bar_segments(
            state_track_df, state_key=state_col, time_key=time_col, state_color_map=state_colors,
        )

        contact_track_df = df_timepoints[
            (df_timepoints[sample_col].astype(str) == str(sample_name))
            & (df_timepoints[track_col].astype(str) == str(track_id))
        ].copy()
        contact_track_df[time_col] = pd.to_numeric(contact_track_df[time_col], errors="coerce")
        contact_track_df = contact_track_df[
            (contact_track_df[time_col] >= float(tmin)) & (contact_track_df[time_col] <= float(tmax))
        ].sort_values(time_col)
        contact_segments = _compute_contact_bar_segments(
            contact_track_df[time_col].to_numpy(dtype=float),
            pd.to_numeric(contact_track_df[contact_col], errors="coerce").fillna(0).astype(bool).to_numpy(),
            min_bout_length,
        )

        prepared_rows.append(
            {
                "sample_name": sample_name,
                "track_id": track_id,
                "state_segments": state_segments,
                "xlim": xlim,
                "contact_segments": contact_segments,
            }
        )

    rows_per_page = max(1, int(rows_per_page))
    with PdfPages(pdf_path) as pdf:
        for sample_name in sample_order:
            sample_rows = [r for r in prepared_rows if r["sample_name"] == sample_name]
            for page_rows in _chunk_list(sample_rows, rows_per_page):
                fig = _plot_track_contact_overview_page(
                    page_rows, sample_name=sample_name, contact_col=contact_col,
                    min_bout_length=min_bout_length,
                )
                pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
                plt.close(fig)

    if verbose:
        print(f"Saved track contact overview report: {pdf_path}")

    return {
        "contact_col": str(contact_col),
        "min_bout_length": int(min_bout_length),
        "pdf_path": str(pdf_path),
        "n_tracks": int(len(prepared_rows)),
        "n_samples": int(len(sample_order)),
    }
