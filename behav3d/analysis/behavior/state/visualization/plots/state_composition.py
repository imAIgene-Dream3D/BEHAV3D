import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.analysis.behavior.utils import _natural_sort_key
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    draw_thin_stacked_proportion_barh,
    compute_condition_diff_stats_pairwise,
    legend_layout,
    plot_condition_diff_grid,
    plot_page_stacked_proportion_barh_grid,
    stacked_proportion_barh_rows_per_page,
    _resolve_effective_group_cols,
    _make_group_label,
)


A4_PORTRAIT = (8.27, 11.69)
_A4_CONTENT_H = 9.0    # usable height (in) after title + legend + margins
_A4_CONTENT_W = 7.2    # usable width (in) after side margins
_MIN_PANEL_H = 2.2     # minimum comfortable panel height → max 4 rows on A4
_MIN_PANEL_W = 2.0     # minimum comfortable panel width  → max 3 cols on A4


def _panels_per_a4_page(grid_ncols):
    """Return (panels_per_page, ncols, max_rows) fitting comfortably on A4."""
    ncols = max(1, min(int(grid_ncols), int(_A4_CONTENT_W / _MIN_PANEL_W)))
    max_rows = max(1, int(_A4_CONTENT_H / _MIN_PANEL_H))
    return ncols * max_rows, ncols, max_rows


def _chunk_list(lst, n):
    """Split list into chunks of at most n items."""
    n = max(1, int(n))
    return [lst[i : i + n] for i in range(0, max(1, len(lst)), n)]


def _pub_style_ax(ax):
    """Apply publication-ready style to a data axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=7, length=3, width=0.8)


def _wrap_title(text, max_chars=28):
    """Wrap a long label at underscore boundaries, capped at 2 lines."""
    if len(text) <= max_chars:
        return text
    parts = text.split("_")
    lines, line = [], ""
    for p in parts:
        candidate = (line + "_" + p) if line else p
        if line and len(candidate) > max_chars:
            lines.append(line)
            line = p
            if len(lines) >= 2:
                break
        else:
            line = candidate
    if line and len(lines) < 2:
        lines.append(line)
    return "\n".join(lines)


def plot_state_composition_over_time(
    adata,
    *,
    time_col: str = "position_t",
    state_col: str = "ClusterID",
    sample_col: str = "sample_name",
    relative: bool = True,
    group_by_sample: bool = False,
    state_order=None,
    figsize=(8, 4),
    show: bool = True,
    ax=None,
    legend: bool = True,
    state_colors=None,
    sample_title_fontsize: float = 8,
    sample_title_pad: float = 2,
):
    """
    Compute (and optionally plot) continuous stacked state composition over time.

    Returns
    -------
    If show == False:
        dict of DataFrames:
            {sample_name -> DataFrame(time x state)}
            or {"all": DataFrame} if not grouped
    If show == True:
        (data_dict, fig, ax/axes)
    """

    obs = adata.obs

    # --- validate ---
    for c in [time_col, state_col]:
        if c not in obs.columns:
            raise KeyError(f"{c} not found in adata.obs")
    if group_by_sample and sample_col not in obs.columns:
        raise KeyError(f"{sample_col} not found in adata.obs")

    # --- clean ---
    cols = [time_col, state_col] + ([sample_col] if group_by_sample else [])
    df = obs[cols].dropna().copy()
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    df = df.dropna(subset=[time_col])
    df[time_col] = df[time_col].astype(int)

    if state_order is None:
        state_order = (
            df[state_col]
            .astype(str)
            .value_counts()
            .index
            .tolist()
        )
    state_color_map = _normalize_label_color_map(state_order, colors=state_colors)

    # --- compute per-sample matrices ---
    def _compute(panel_df):
        panel_df = panel_df.copy()
        panel_df[state_col] = panel_df[state_col].astype(str)

        mat = (
            panel_df
            .groupby([time_col, state_col], observed=True)
            .size()
            .unstack(fill_value=0)
            .sort_index()
        )

        for s in state_order:
            if s not in mat.columns:
                mat[s] = 0
        mat = mat[state_order]

        if relative:
            mat = mat.div(
                mat.sum(axis=1).replace(0, np.nan),
                axis=0
            ).fillna(0.0)

        return mat

    if group_by_sample:
        data = {
            str(s): _compute(df[df[sample_col].astype(str) == str(s)])
            for s in df[sample_col].astype(str).unique()
        }
    else:
        data = {"all": _compute(df)}

    # --- early exit if no plotting ---
    if not show:
        return data

    # --- plotting ---
    def _plot_one(mat, ax, title=None, title_fontsize=None, title_pad=None):
        x = mat.index.to_numpy()
        bottom = np.zeros(len(x))

        for s in reversed(list(mat.columns)):
            v = mat[s].to_numpy()
            ax.bar(
                x,
                v,
                bottom=bottom,
                width=1.0,
                align="edge",
                linewidth=0,
                color=state_color_map.get(str(s), None),
                label=str(s),
            )
            bottom += v

        ax.set_xlim(x.min(), x.max() + 1)
        ax.margins(x=0)
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion" if relative else "Count")
        if title:
            ax.set_title(title, fontsize=title_fontsize, pad=title_pad)

    if group_by_sample:
        n = len(data)
        fig, axes = plt.subplots(
            nrows=n,
            ncols=1,
            sharex=True,
            figsize=(figsize[0], figsize[1] * n),
        )
        if n == 1:
            axes = [axes]

        for ax_i, (samp, mat) in zip(axes, data.items()):
            _plot_one(
                mat,
                ax_i,
                title=f"{sample_col}: {samp}",
                title_fontsize=sample_title_fontsize,
                title_pad=sample_title_pad,
            )

        if legend:
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower center",
                ncol=min(len(labels), 6),
                frameon=False,
            )
            fig.tight_layout(rect=(0, 0.06, 1, 1))
        else:
            fig.tight_layout()

        return data, fig, axes

    # single axis
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    _plot_one(data["all"], ax)

    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[::-1],
            labels[::-1],
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            frameon=False,
        )
        fig.tight_layout()

    return data, fig, ax


def _prepare_state_composition_df(
    adata,
    *,
    time_col="position_t",
    state_col="ClusterID",
    sample_col="sample_name",
    state_order=None,
):
    """Validate/clean obs and return normalized DataFrame + state/sample ordering."""
    obs = adata.obs
    required = [time_col, state_col, sample_col]
    missing = [c for c in required if c not in obs.columns]
    if len(missing) > 0:
        raise KeyError(f"Missing required columns in adata.obs: {missing}")

    df = obs[required].dropna().copy()
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    df = df.dropna(subset=[time_col]).copy()
    if len(df) == 0:
        raise ValueError("No valid rows remain after filtering NaNs in required columns.")
    df[time_col] = df[time_col].astype(int)
    df[state_col] = df[state_col].astype(str)
    df[sample_col] = df[sample_col].astype(str)

    observed_states = sorted(df[state_col].unique().tolist(), key=_natural_sort_key)
    if state_order is None:
        state_order = [str(s) for s in observed_states]
    else:
        state_order = [str(s) for s in list(state_order)]
        extras = [s for s in observed_states if str(s) not in state_order]
        state_order.extend([str(s) for s in extras])

    sample_order = df[sample_col].drop_duplicates().astype(str).tolist()
    return df, state_order, sample_order


def compute_relative_state_time_matrix(panel_df, *, time_col, state_col, state_order):
    """Build relative state-composition matrix (time x state) for a panel."""
    mat = (
        panel_df.groupby([time_col, state_col], observed=True)
        .size()
        .unstack(fill_value=0)
        .sort_index()
    )
    for state in state_order:
        if state not in mat.columns:
            mat[state] = 0
    mat = mat[state_order]
    mat = mat.div(mat.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return mat


def _compute_relative_matrices_by_sample(
    df,
    *,
    time_col,
    state_col,
    sample_col,
    state_order,
    sample_order,
):
    """Compute per-sample and pooled relative composition matrices."""
    relative_by_sample = {}
    for sample in sample_order:
        panel = df[df[sample_col] == sample]
        relative_by_sample[str(sample)] = compute_relative_state_time_matrix(
            panel,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )

    relative_pooled = compute_relative_state_time_matrix(
        df,
        time_col=time_col,
        state_col=state_col,
        state_order=state_order,
    )
    return relative_by_sample, relative_pooled


def _compute_overall_relative_composition_by_sample(
    df,
    *,
    state_col,
    sample_col,
    state_order,
    sample_order,
):
    """Compute pooled relative composition per sample across the full experiment."""
    overall_by_sample = {}
    for sample in sample_order:
        panel = df[df[sample_col] == sample]
        counts = panel[state_col].value_counts().reindex(state_order, fill_value=0).astype(float)
        total = float(counts.sum())
        if total > 0:
            counts = counts / total
        overall_by_sample[str(sample)] = counts
    return overall_by_sample


def _compute_count_matrix(panel_df, *, time_col, state_col, state_order):
    """Build raw count state-composition matrix (time x state) for a panel."""
    mat = (
        panel_df.groupby([time_col, state_col], observed=True)
        .size()
        .unstack(fill_value=0)
        .sort_index()
    )
    for state in state_order:
        if state not in mat.columns:
            mat[state] = 0
    mat = mat[state_order]
    return mat


def _compute_count_matrices_by_sample(
    df,
    *,
    time_col,
    state_col,
    sample_col,
    state_order,
    sample_order,
):
    """Compute per-sample and pooled raw count matrices."""
    count_by_sample = {}
    for sample in sample_order:
        panel = df[df[sample_col] == sample]
        count_by_sample[str(sample)] = _compute_count_matrix(
            panel,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )
    count_pooled = _compute_count_matrix(
        df,
        time_col=time_col,
        state_col=state_col,
        state_order=state_order,
    )
    return count_by_sample, count_pooled


def _compute_overall_count_by_sample(
    df,
    *,
    state_col,
    sample_col,
    state_order,
    sample_order,
):
    """Compute raw count per state per sample (pooled across time)."""
    overall_by_sample = {}
    for sample in sample_order:
        panel = df[df[sample_col] == sample]
        counts = panel[state_col].value_counts().reindex(state_order, fill_value=0).astype(int)
        overall_by_sample[str(sample)] = counts
    return overall_by_sample


def _compute_relative_auc_table(
    relative_by_sample,
    *,
    state_order,
    include_pooled_summary=True,
    relative_pooled=None,
    pooled_sample_name="__all__",
):
    """Compute trapezoidal AUC per sample and state from relative composition curves."""

    def _rows_for_matrix(sample_name, mat):
        rows = []
        x = mat.index.to_numpy(dtype=float)
        n_timepoints = int(len(x))
        time_min = float(x.min()) if n_timepoints > 0 else np.nan
        time_max = float(x.max()) if n_timepoints > 0 else np.nan
        for state in state_order:
            y = mat[state].to_numpy(dtype=float) if state in mat.columns else np.zeros(n_timepoints, dtype=float)
            auc = float(np.trapezoid(y=y, x=x)) if n_timepoints >= 2 else 0.0
            rows.append(
                {
                    "sample_name": str(sample_name),
                    "state_id": str(state),
                    "auc_relative": auc,
                    "n_timepoints": n_timepoints,
                    "time_min": time_min,
                    "time_max": time_max,
                }
            )
        return rows

    rows = []
    for sample_name, mat in relative_by_sample.items():
        rows.extend(_rows_for_matrix(sample_name, mat))
    if include_pooled_summary and relative_pooled is not None:
        rows.extend(_rows_for_matrix(pooled_sample_name, relative_pooled))

    return pd.DataFrame(rows)


def _build_relative_plot_data_table(relative_by_sample, plot_view=None, label_col_name="sample_name"):
    """Build long-form CSV table for all relative curves used in the report pages.

    ``label_col_name`` lets this same builder be reused for the grouped report
    rows (``relative_by_group``, keyed by group label rather than sample name).
    """
    rows = []
    for sample_name, mat in relative_by_sample.items():
        if mat is None or len(mat) == 0:
            continue
        for state_id in mat.columns:
            y = mat[state_id].to_numpy(dtype=float)
            x = mat.index.to_numpy(dtype=float)
            for t, v in zip(x.tolist(), y.tolist()):
                rows.append(
                    {
                        label_col_name: str(sample_name),
                        "time": float(t),
                        "state_id": str(state_id),
                        "relative_proportion": float(v),
                    }
                )
    out = pd.DataFrame(rows)
    if plot_view is not None:
        out["plot_view"] = str(plot_view)
    return out


def _build_overall_summary_plot_data_table(overall_by_sample, plot_view=None, label_col_name="sample_name"):
    """Build long-form CSV table for pooled per-sample summary bars.

    ``label_col_name`` lets this same builder be reused for the grouped report
    rows (``overall_by_group``, keyed by group label rather than sample name).
    """
    rows = []
    for sample_name, ser in overall_by_sample.items():
        if ser is None or len(ser) == 0:
            continue
        for state_id, value in ser.items():
            rows.append(
                {
                    label_col_name: str(sample_name),
                    "time": np.nan,
                    "state_id": str(state_id),
                    "relative_proportion": float(value),
                }
            )
    out = pd.DataFrame(rows)
    if plot_view is not None:
        out["plot_view"] = str(plot_view)
    return out


def _paginate_samples(sample_names, samples_per_page=4):
    """Split sample names into page-sized chunks while preserving order."""
    sample_names = [str(x) for x in list(sample_names)]
    samples_per_page = max(int(samples_per_page), 1)
    return [
        sample_names[i : i + samples_per_page]
        for i in range(0, len(sample_names), samples_per_page)
    ]


def _plot_overall_summary_bar(ax, overall, *, state_order, state_colors, show_ylabel=False):
    """Draw a narrow vertical stacked bar for pooled sample composition."""
    bottom = 0.0
    for state in reversed(list(state_order)):
        value = float(overall.get(state, 0.0))
        ax.bar(
            [0.0],
            [value],
            bottom=bottom,
            width=0.8,
            color=state_colors[state],
            linewidth=0,
        )
        bottom += value
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([])
    if show_ylabel:
        ax.set_ylabel("Proportion", fontsize=7)
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Overall", fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)


def _plot_overall_count_bar(ax, overall_counts, *, state_order, state_colors):
    """Draw a narrow vertical stacked bar showing raw cell counts, scaled to its own total."""
    total = float(sum(int(overall_counts.get(s, 0)) for s in state_order))
    ylim_top = max(total * 1.1, 1.0)
    bottom = 0.0
    for state in reversed(list(state_order)):
        value = float(int(overall_counts.get(state, 0)))
        ax.bar(
            [0.0],
            [value],
            bottom=bottom,
            width=0.8,
            color=state_colors[state],
            linewidth=0,
        )
        bottom += value
    ax.set_ylim(0.0, ylim_top)
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_xlabel("Overall", fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)


def _plot_page_relative_stacked_grid(
    relative_by_sample,
    *,
    overall_by_sample,
    state_order,
    time_col,
    sample_col,
    grid_ncols,
    figsize_per_panel,
    state_colors,
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """One A4 page: relative stacked composition per sample (grid)."""
    samples = list(relative_by_sample.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(samples)) / float(ncols))))
    fig = plt.figure(figsize=A4_PORTRAIT)
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    first_handles = []
    first_labels = []

    for i, sample in enumerate(samples):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1], sharey=ax)
        mat = relative_by_sample[sample]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(
                _wrap_title(sample) + "\n(empty)",
                fontsize=sample_title_fontsize,
                fontweight="bold",
                pad=sample_title_pad,
            )
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
        colors = [state_colors[state] for state in reversed(list(state_order))]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            _wrap_title(sample),
            fontsize=sample_title_fontsize,
            fontweight="bold",
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col, fontsize=8)
        ax.set_ylabel("Proportion", fontsize=8)
        _pub_style_ax(ax)
        _plot_overall_summary_bar(
            ax_bar,
            overall_by_sample[sample],
            state_order=state_order,
            state_colors=state_colors,
        )
        if len(first_labels) == 0:
            first_handles, first_labels = ax.get_legend_handles_labels()

    for i in range(len(samples), nrows * ncols):
        ax_empty = fig.add_subplot(outer[i])
        ax_empty.axis("off")

    if len(first_labels) > 0:
        legend_ncol, _, legend_margin_in = legend_layout(len(first_labels), base_margin_in=0.82)
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=legend_ncol,
            frameon=False,
            fontsize=7,
        )
        fig.tight_layout(rect=(0.03, legend_margin_in / A4_PORTRAIT[1], 1, 0.93))
    else:
        fig.tight_layout(rect=(0.03, 0, 1, 0.93))
    fig.suptitle("Relative State Composition (Stacked) by Sample", y=0.97, fontsize=12, fontweight="bold")
    return fig


def _plot_page_count_stacked_grid(
    count_by_sample,
    *,
    overall_count_by_sample,
    state_order,
    time_col,
    sample_col,
    grid_ncols,
    figsize_per_panel,
    state_colors,
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """One A4 page: absolute cell count stacked composition per sample (grid)."""
    samples = list(count_by_sample.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(samples)) / float(ncols))))

    all_maxes = []
    for mat in count_by_sample.values():
        if mat is not None and len(mat) > 0:
            row_totals = mat.sum(axis=1)
            if len(row_totals) > 0:
                all_maxes.append(float(row_totals.max()))
    global_ymax = max(all_maxes) * 1.1 if all_maxes else 1.0

    fig = plt.figure(figsize=A4_PORTRAIT)
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    first_handles = []
    first_labels = []

    for i, sample in enumerate(samples):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1])
        mat = count_by_sample[sample]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(
                _wrap_title(sample) + "\n(empty)",
                fontsize=sample_title_fontsize,
                fontweight="bold",
                pad=sample_title_pad,
            )
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
        colors = [state_colors[state] for state in reversed(list(state_order))]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)
        ax.set_ylim(0.0, global_ymax)
        ax.set_title(
            _wrap_title(sample),
            fontsize=sample_title_fontsize,
            fontweight="bold",
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col, fontsize=8)
        ax.set_ylabel("Cell count", fontsize=8)
        _pub_style_ax(ax)
        _plot_overall_count_bar(
            ax_bar,
            overall_count_by_sample[sample],
            state_order=state_order,
            state_colors=state_colors,
        )
        if len(first_labels) == 0:
            first_handles, first_labels = ax.get_legend_handles_labels()

    for i in range(len(samples), nrows * ncols):
        ax_empty = fig.add_subplot(outer[i])
        ax_empty.axis("off")

    if len(first_labels) > 0:
        legend_ncol, _, legend_margin_in = legend_layout(len(first_labels), base_margin_in=0.82)
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=legend_ncol,
            frameon=False,
            fontsize=7,
        )
        fig.tight_layout(rect=(0.03, legend_margin_in / A4_PORTRAIT[1], 1, 0.93))
    else:
        fig.tight_layout(rect=(0.03, 0, 1, 0.93))
    fig.suptitle("Absolute Cell Count (Stacked) by Sample", y=0.97, fontsize=12, fontweight="bold")
    return fig


def save_state_condition_comparison_report(
    adata,
    output_pdf_path,
    output_csv_path=None,
    *,
    state_col="ClusterID",
    sample_col="sample_name",
    condition_col,
    group_cols=None,
    group_x=None,
    group_x_levels_map=None,
    condition_groups=None,
    state_order=None,
    state_colors=None,
    verbose=True,
):
    """Per-cluster overall (pooled, non-per-timepoint) proportion difference between every
    pairwise combination of `condition_col`'s levels, with Welch's two-sided unpaired t-test
    significance stars — one row per pairwise comparison, one column per group, optionally
    split into multiple groups (columns) via `group_x`/`group_cols`.

    diff = mean of the second level minus mean of the first level in each pair.

    condition_groups : dict[str, str], optional
        Maps raw condition_col levels to merged group labels (see
        compute_condition_diff_stats_pairwise) - when given, compares the merged groups
        instead of every raw level pairwise.

    group_x_levels_map : dict[str, str], optional
        Maps raw `group_x` levels to merged group labels - when given, `group_x` is pooled
        into the merged labels (columns) instead of showing one column per raw level.
    """
    obs = adata.obs
    effective_group_cols, _ = _resolve_effective_group_cols(group_cols, group_x, None)
    # Grouping by the same column being compared is degenerate (every group would
    # contain a single condition level) and, worse, would duplicate condition_col
    # inside metadata_cols below, turning df[condition_col] into a 2-column
    # DataFrame instead of a Series and breaking downstream .unique() calls.
    if condition_col in effective_group_cols and bool(verbose):
        print(
            f"  Note: '{condition_col}' is both the comparison condition and a "
            "requested group column — ignoring it as a group."
        )
    valid_group_cols = [c for c in effective_group_cols if c in obs.columns and c != condition_col]
    required = [state_col, sample_col, condition_col] + valid_group_cols
    missing = [c for c in required if c not in obs.columns]
    if len(missing) > 0:
        raise KeyError(f"Missing required columns in adata.obs: {missing}")

    df = obs[required].dropna(subset=[state_col, sample_col, condition_col]).copy()
    df[state_col] = df[state_col].astype(str)
    df[sample_col] = df[sample_col].astype(str)
    df[condition_col] = df[condition_col].astype(str)
    for gc in valid_group_cols:
        df[gc] = df[gc].astype(str)
        if gc == group_x and group_x_levels_map:
            df[gc] = df[gc].map(group_x_levels_map)
    if group_x in valid_group_cols and group_x_levels_map:
        # A level_map applies only to whichever raw levels are still present - map() leaves
        # anything outside the mapping as NaN, so drop those rows rather than let them
        # silently poison the grouped comparison.
        df = df.dropna(subset=[group_x])
    if len(df) == 0:
        raise ValueError("No valid rows remain after filtering NaNs in required columns.")

    # The unit compared below is (sample, condition_col level), not the sample alone.
    # condition_col is usually constant within a sample (e.g. a per-well treatment),
    # in which case this is a no-op relabeling - each sample still maps to exactly
    # one unit. But for a merged cell-type group, condition_col can instead be a
    # per-cell tag (e.g. origin_cell_type) that varies *within* a sample. Collapsing
    # straight to one row per sample would arbitrarily keep whichever level happened
    # to appear first and silently leave a single degenerate condition level with
    # zero pairs left to compare (empty report). Using this composite unit lets a
    # sample that contains multiple condition levels contribute one row per level.
    unit_col = "__comparison_unit__"
    df[unit_col] = df[sample_col] + " | " + df[condition_col]

    observed_states = sorted(df[state_col].unique().tolist(), key=_natural_sort_key)
    if state_order is None:
        resolved_state_order = _apply_state_order(
            [str(s) for s in observed_states], _get_classification_state_order(adata, state_col)
        )
    else:
        resolved_state_order = [str(s) for s in list(state_order)]
        extras = [s for s in observed_states if str(s) not in resolved_state_order]
        resolved_state_order.extend(extras)

    unit_order = df[unit_col].drop_duplicates().tolist()

    overall_by_unit = _compute_overall_relative_composition_by_sample(
        df,
        state_col=state_col,
        sample_col=unit_col,
        state_order=resolved_state_order,
        sample_order=unit_order,
    )
    per_unit_df = pd.DataFrame(overall_by_unit).T.reindex(columns=resolved_state_order, fill_value=0.0)

    metadata_cols = [unit_col, condition_col] + valid_group_cols
    unit_metadata = (
        df[metadata_cols].drop_duplicates(subset=[unit_col]).set_index(unit_col)
    )

    if state_colors is None:
        state_colors = _get_classification_state_colors(adata, state_col)
    resolved_colors = _normalize_label_color_map(resolved_state_order, colors=state_colors, cmap_name="tab20")

    diff_stats_by_group = compute_condition_diff_stats_pairwise(
        per_unit_df,
        unit_metadata,
        class_order=resolved_state_order,
        condition_col=condition_col,
        group_cols=valid_group_cols,
        condition_groups=condition_groups,
    )

    output_pdf_path = Path(output_pdf_path)
    output_pdf_path.parent.mkdir(parents=True, exist_ok=True)
    title = f"{state_col} — {condition_col} pairwise comparison"
    result = plot_condition_diff_grid(
        diff_stats_by_group,
        class_order=resolved_state_order,
        colors=resolved_colors,
        title=title,
        out_pdf=output_pdf_path,
        out_csv=output_csv_path,
    )
    if bool(verbose):
        print(f"Saved condition comparison report: {result['pdf_path']}")
    return result


def _compute_grouped_relative_matrices(
    df,
    *,
    group_cols,
    time_col,
    state_col,
    state_order,
):
    """Compute relative composition matrices per unique group label."""
    group_labels = _make_group_label(df, group_cols)
    tmp = df.copy()
    tmp["_group_label"] = group_labels.values
    unique_groups = tmp["_group_label"].dropna().unique().tolist()

    relative_by_group = {}
    overall_by_group = {}
    for grp in unique_groups:
        panel = tmp[tmp["_group_label"] == grp]
        relative_by_group[str(grp)] = compute_relative_state_time_matrix(
            panel,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )
        counts = panel[state_col].value_counts().reindex(state_order, fill_value=0).astype(float)
        total = float(counts.sum())
        if total > 0:
            counts = counts / total
        overall_by_group[str(grp)] = counts

    return relative_by_group, overall_by_group


def _compute_grouped_count_matrices(
    df,
    *,
    group_cols,
    time_col,
    state_col,
    state_order,
):
    """Compute raw count matrices per unique group label."""
    group_labels = _make_group_label(df, group_cols)
    tmp = df.copy()
    tmp["_group_label"] = group_labels.values
    unique_groups = tmp["_group_label"].dropna().unique().tolist()

    count_by_group = {}
    overall_count_by_group = {}
    for grp in unique_groups:
        panel = tmp[tmp["_group_label"] == grp]
        count_by_group[str(grp)] = _compute_count_matrix(
            panel,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )
        counts = panel[state_col].value_counts().reindex(state_order, fill_value=0).astype(int)
        overall_count_by_group[str(grp)] = counts

    return count_by_group, overall_count_by_group


def _plot_page_grouped_stacked_grid(
    relative_by_group,
    *,
    overall_by_group,
    state_order,
    time_col,
    group_label_title,
    grid_ncols,
    figsize_per_panel,
    state_colors,
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """One A4 page: relative stacked composition per group (flat grid, for 3+ group_cols)."""
    groups = list(relative_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))
    fig = plt.figure(figsize=A4_PORTRAIT)
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    first_handles = []
    first_labels = []

    for i, grp in enumerate(groups):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1], sharey=ax)
        mat = relative_by_group[grp]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(f"{grp} (empty)", fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
        colors = [state_colors[state] for state in reversed(list(state_order))]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(grp, fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
        ax.set_xlabel(time_col, fontsize=8)
        ax.set_ylabel("Proportion", fontsize=8)
        _pub_style_ax(ax)
        _plot_overall_summary_bar(
            ax_bar,
            overall_by_group[grp],
            state_order=state_order,
            state_colors=state_colors,
        )
        if len(first_labels) == 0:
            first_handles, first_labels = ax.get_legend_handles_labels()

    for i in range(len(groups), nrows * ncols):
        ax_empty = fig.add_subplot(outer[i])
        ax_empty.axis("off")

    if len(first_labels) > 0:
        legend_ncol, _, legend_margin_in = legend_layout(len(first_labels), base_margin_in=0.82)
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=legend_ncol,
            frameon=False,
            fontsize=7,
        )
        fig.tight_layout(rect=(0.03, legend_margin_in / A4_PORTRAIT[1], 1, 0.93))
    else:
        fig.tight_layout(rect=(0.03, 0, 1, 0.93))
    fig.suptitle(
        f"Grouped Relative State Composition — {group_label_title}",
        y=0.97,
        fontsize=12,
        fontweight="bold",
    )
    return fig


def _plot_page_grouped_count_stacked_grid(
    count_by_group,
    *,
    overall_count_by_group,
    state_order,
    time_col,
    group_label_title,
    grid_ncols,
    figsize_per_panel,
    state_colors,
    global_ymax=None,
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """One A4 page: absolute cell count stacked composition per group (flat grid, for 3+ group_cols)."""
    groups = list(count_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))

    if global_ymax is None:
        all_maxes = []
        for mat in count_by_group.values():
            if mat is not None and len(mat) > 0:
                row_totals = mat.sum(axis=1)
                if len(row_totals) > 0:
                    all_maxes.append(float(row_totals.max()))
        global_ymax = max(all_maxes) * 1.1 if all_maxes else 1.0

    fig = plt.figure(figsize=A4_PORTRAIT)
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    first_handles = []
    first_labels = []

    for i, grp in enumerate(groups):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1])
        mat = count_by_group[grp]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(f"{grp} (empty)", fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
        colors = [state_colors[state] for state in reversed(list(state_order))]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)
        ax.set_ylim(0.0, global_ymax)
        ax.set_title(grp, fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
        ax.set_xlabel(time_col, fontsize=8)
        ax.set_ylabel("Cell count", fontsize=8)
        _pub_style_ax(ax)
        _plot_overall_count_bar(
            ax_bar,
            overall_count_by_group[grp],
            state_order=state_order,
            state_colors=state_colors,
        )
        if len(first_labels) == 0:
            first_handles, first_labels = ax.get_legend_handles_labels()

    for i in range(len(groups), nrows * ncols):
        ax_empty = fig.add_subplot(outer[i])
        ax_empty.axis("off")

    if len(first_labels) > 0:
        legend_ncol, _, legend_margin_in = legend_layout(len(first_labels), base_margin_in=0.82)
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=legend_ncol,
            frameon=False,
            fontsize=7,
        )
        fig.tight_layout(rect=(0.03, legend_margin_in / A4_PORTRAIT[1], 1, 0.93))
    else:
        fig.tight_layout(rect=(0.03, 0, 1, 0.93))
    fig.suptitle(
        f"Grouped Absolute Cell Count — {group_label_title}",
        y=0.97,
        fontsize=12,
        fontweight="bold",
    )
    return fig


def _plot_page_grouped_2d_grid(
    data_by_group,
    *,
    overall_by_group,
    group_cols,
    unique_vals_per_col,
    state_order,
    time_col,
    group_label_title,
    state_colors,
    mode="relative",
    global_ymax=None,
    row_slice=None,
    axis_cols=None,
    sample_title_fontsize=8,
    sample_title_pad=2,
    page_size=A4_PORTRAIT,
):
    """
    One A4 page: grouped state composition arranged in a 2D grid.

    For 1 group_col: single column of panels, one row per unique value.
    For 2 group_cols: by default the column with more unique values goes to
    rows (Y) and the column with fewer unique values goes to columns (X);
    pass ``axis_cols=(col_y, col_x)`` to pick the axes explicitly instead.
    Header labels are drawn on each axis. row_slice=(start, end) selects a
    subset of Y rows for pagination. ``mode="count"`` plots raw cell counts
    (using ``global_ymax`` for the shared y-axis) instead of proportions.
    """
    if len(group_cols) == 2:
        if axis_cols is not None:
            col_y, col_x = axis_cols
        else:
            # Assign high-cardinality column to rows for A4 portrait fit
            col_y, col_x = sorted(group_cols, key=lambda c: -len(unique_vals_per_col[c]))
        col_x_vals = unique_vals_per_col[col_x]
        col_y_vals_all = unique_vals_per_col[col_y]
        if row_slice is not None:
            page_col_y_vals = col_y_vals_all[row_slice[0] : row_slice[1]]
        else:
            page_col_y_vals = col_y_vals_all

        nrows_data = max(1, len(page_col_y_vals))
        ncols_data = max(1, len(col_x_vals))

        header_h = 0.06
        header_w = 0.06
        fig = plt.figure(figsize=page_size)
        outer = GridSpec(
            nrows=nrows_data + 1,
            ncols=ncols_data + 1,
            figure=fig,
            height_ratios=[header_h] + [1.0] * nrows_data,
            width_ratios=[header_w] + [1.0] * ncols_data,
            hspace=0.45,
            wspace=0.35,
        )

        # Corner cell
        ax_corner = fig.add_subplot(outer[0, 0])
        ax_corner.axis("off")

        # Column header labels (col_x values)
        for c, col_x_val in enumerate(col_x_vals):
            ax_ch = fig.add_subplot(outer[0, c + 1])
            ax_ch.axis("off")
            ax_ch.text(
                0.5, 0.5, str(col_x_val),
                ha="center", va="center",
                fontsize=9, fontweight="bold",
                transform=ax_ch.transAxes,
            )

        # Row header labels (col_y values) + data cells
        first_handles = []
        first_labels = []

        for r, col_y_val in enumerate(page_col_y_vals):
            ax_rh = fig.add_subplot(outer[r + 1, 0])
            ax_rh.axis("off")
            ax_rh.text(
                0.5, 0.5, str(col_y_val),
                ha="center", va="center",
                fontsize=9, fontweight="bold",
                rotation=90,
                transform=ax_rh.transAxes,
            )

            for c, col_x_val in enumerate(col_x_vals):
                # Reconstruct key in group_cols order
                vals_map = {col_y: col_y_val, col_x: col_x_val}
                key = " | ".join(str(vals_map[gc]) for gc in group_cols)

                inner = outer[r + 1, c + 1].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.08)
                ax = fig.add_subplot(inner[0, 0])
                ax_bar = fig.add_subplot(inner[0, 1])

                if key not in data_by_group:
                    ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                    ax.axis("off")
                    ax_bar.axis("off")
                    continue

                mat = data_by_group[key]
                x = mat.index.to_numpy(dtype=float)
                if len(x) == 0:
                    ax.text(0.5, 0.5, "empty", ha="center", va="center", fontsize=8, transform=ax.transAxes)
                    ax.axis("off")
                    ax_bar.axis("off")
                    continue

                y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
                colors = [state_colors[state] for state in reversed(list(state_order))]
                ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)

                if mode == "relative":
                    ax.set_ylim(0.0, 1.0)
                    ax.set_ylabel("Proportion", fontsize=7)
                else:
                    ax.set_ylim(0.0, global_ymax if global_ymax is not None else 1.0)
                    ax.set_ylabel("Cell count", fontsize=7)

                ax.set_xlabel(time_col, fontsize=7)
                _pub_style_ax(ax)
                ax.tick_params(labelsize=6)

                if mode == "relative":
                    _plot_overall_summary_bar(
                        ax_bar, overall_by_group[key],
                        state_order=state_order, state_colors=state_colors,
                    )
                else:
                    _plot_overall_count_bar(
                        ax_bar, overall_by_group[key],
                        state_order=state_order, state_colors=state_colors,
                    )

                if len(first_labels) == 0:
                    first_handles, first_labels = ax.get_legend_handles_labels()

        mode_label = "Relative State Composition" if mode == "relative" else "Absolute Cell Count"
        title_str = f"Grouped {mode_label} — {col_x} (columns) × {col_y} (rows)"
        if row_slice is not None and row_slice[0] > 0:
            total_rows = len(col_y_vals_all)
            title_str += f"\n(rows {row_slice[0]+1}–{row_slice[1]} of {total_rows})"

    else:
        # 1 group_col: single-column layout
        col_y = group_cols[0]
        col_y_vals_all = unique_vals_per_col[col_y]
        if row_slice is not None:
            page_col_y_vals = col_y_vals_all[row_slice[0] : row_slice[1]]
        else:
            page_col_y_vals = col_y_vals_all

        nrows_data = max(1, len(page_col_y_vals))
        fig = plt.figure(figsize=page_size)
        outer = GridSpec(nrows=nrows_data, ncols=1, figure=fig, hspace=0.45)
        first_handles = []
        first_labels = []

        for r, col_y_val in enumerate(page_col_y_vals):
            key = str(col_y_val)
            inner = outer[r].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.08)
            ax = fig.add_subplot(inner[0, 0])
            ax_bar = fig.add_subplot(inner[0, 1])

            if key not in data_by_group:
                ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                ax.axis("off")
                ax_bar.axis("off")
                continue

            mat = data_by_group[key]
            x = mat.index.to_numpy(dtype=float)
            if len(x) == 0:
                ax.set_title(f"{col_y_val} (empty)", fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
                ax.axis("off")
                ax_bar.axis("off")
                continue

            y_arrays = [mat[state].to_numpy(dtype=float) for state in reversed(list(state_order))]
            colors = [state_colors[state] for state in reversed(list(state_order))]
            ax.stackplot(x, *y_arrays, labels=[str(s) for s in reversed(list(state_order))], colors=colors, alpha=0.85)

            if mode == "relative":
                ax.set_ylim(0.0, 1.0)
                ax.set_ylabel("Proportion", fontsize=8)
            else:
                ax.set_ylim(0.0, global_ymax if global_ymax is not None else 1.0)
                ax.set_ylabel("Cell count", fontsize=8)

            ax.set_title(f"{col_y}: {col_y_val}", fontsize=sample_title_fontsize, fontweight="bold", pad=sample_title_pad)
            ax.set_xlabel(time_col, fontsize=8)
            _pub_style_ax(ax)

            if mode == "relative":
                _plot_overall_summary_bar(
                    ax_bar, overall_by_group[key],
                    state_order=state_order, state_colors=state_colors,
                )
            else:
                _plot_overall_count_bar(
                    ax_bar, overall_by_group[key],
                    state_order=state_order, state_colors=state_colors,
                )

            if len(first_labels) == 0:
                first_handles, first_labels = ax.get_legend_handles_labels()

        mode_label = "Relative State Composition" if mode == "relative" else "Absolute Cell Count"
        title_str = f"Grouped {mode_label} — {group_label_title}"
        if row_slice is not None and row_slice[0] > 0:
            total_rows = len(col_y_vals_all)
            title_str += f"\n(groups {row_slice[0]+1}–{row_slice[1]} of {total_rows})"

    if len(first_labels) > 0:
        legend_ncol, _, legend_margin_in = legend_layout(len(first_labels), base_margin_in=0.82)
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=legend_ncol,
            frameon=False,
            fontsize=7,
        )
        fig.tight_layout(rect=(0.03, legend_margin_in / page_size[1], 1, 0.92))
    else:
        fig.tight_layout(rect=(0.03, 0, 1, 0.92))
    fig.suptitle(title_str, y=0.97, fontsize=11, fontweight="bold")
    return fig


def save_state_composition_report(
    adata,
    output_pdf_path,
    output_csv_path=None,
    *,
    time_col="position_t",
    state_col="ClusterID",
    sample_col="sample_name",
    state_order=None,
    grid_ncols=3,
    figsize_per_panel=(4.0, 2.8),
    include_pooled_summary=True,
    dpi=300,
    verbose=True,
    state_colors=None,
    sample_title_fontsize: float = 8,
    sample_title_pad: float = 2,
    group_cols=None,
    group_x=None,
    group_y=None,
    group_x_levels_map=None,
    group_y_levels_map=None,
):
    """
    Save a combined multi-page PDF report + merged plot-data CSV for relative state composition.

    ``group_x``/``group_y`` explicitly pick the grouped-grid's axes; ``group_cols``
    is the "group per page" column list (unaffected in meaning).

    ``group_x_levels_map``/``group_y_levels_map`` : dict[str, str], optional
        Maps raw ``group_x``/``group_y`` levels to merged group labels (same semantics
        as ``condition_groups`` in ``compute_condition_diff_stats_pairwise``) - when
        given, that axis is pooled into the merged labels instead of showing one panel
        per raw level. Rows whose raw level isn't in the mapping are dropped.

    Outputs:
      1) one combined PDF with all report pages
      2) one merged plot-data CSV for all report views
      3) one AUC CSV
    """
    output_pdf_path = Path(output_pdf_path)
    stem = output_pdf_path.stem
    parent = output_pdf_path.parent
    merged_pdf_path = output_pdf_path
    merged_plot_csv_path = output_pdf_path.with_suffix(".csv")

    if output_csv_path is None:
        output_auc_csv_path = parent / f"{stem}_auc.csv"
    else:
        output_auc_csv_path = Path(output_csv_path)
    if output_auc_csv_path == merged_plot_csv_path:
        output_auc_csv_path = parent / f"{stem}_auc.csv"

    output_pdf_path.parent.mkdir(parents=True, exist_ok=True)
    output_auc_csv_path.parent.mkdir(parents=True, exist_ok=True)

    _state_order_explicit = state_order is not None
    df, state_order, sample_order = _prepare_state_composition_df(
        adata,
        time_col=time_col,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
    )
    if not _state_order_explicit:
        state_order = _apply_state_order(state_order, _get_classification_state_order(adata, state_col))

    relative_by_sample, relative_pooled = _compute_relative_matrices_by_sample(
        df,
        time_col=time_col,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
        sample_order=sample_order,
    )
    overall_by_sample = _compute_overall_relative_composition_by_sample(
        df,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
        sample_order=sample_order,
    )

    count_by_sample, _count_pooled = _compute_count_matrices_by_sample(
        df,
        time_col=time_col,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
        sample_order=sample_order,
    )
    overall_count_by_sample = _compute_overall_count_by_sample(
        df,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
        sample_order=sample_order,
    )

    effective_group_cols, requested_axis_cols = _resolve_effective_group_cols(
        group_cols, group_x, group_y,
    )

    valid_group_cols = []
    if effective_group_cols:
        obs_cols_available = list(adata.obs.columns)
        for gc in effective_group_cols:
            if gc in obs_cols_available:
                valid_group_cols.append(gc)
                df[gc] = adata.obs.loc[df.index, gc].astype(str).fillna("(unknown)").values
                if gc == group_x and group_x_levels_map:
                    df[gc] = df[gc].map(group_x_levels_map)
                elif gc == group_y and group_y_levels_map:
                    df[gc] = df[gc].map(group_y_levels_map)
            else:
                if verbose:
                    print(f"  Warning: group_col '{gc}' not found in adata.obs — skipping.")
    if valid_group_cols:
        # A level_map applies only to whichever raw levels are still present in this
        # dataset - map() leaves anything outside the mapping as NaN, so drop those rows
        # here rather than let them silently poison unique_vals_per_col/group matrices.
        df = df.dropna(subset=valid_group_cols)

    axis_cols = (
        requested_axis_cols
        if requested_axis_cols is not None
        and requested_axis_cols[0] in valid_group_cols
        and requested_axis_cols[1] in valid_group_cols
        else None
    )

    auc_table = _compute_relative_auc_table(
        relative_by_sample,
        state_order=state_order,
        include_pooled_summary=bool(include_pooled_summary),
        relative_pooled=relative_pooled,
    )
    auc_table.to_csv(output_auc_csv_path, index=False)

    if state_colors is None:
        state_colors = _get_classification_state_colors(adata, state_col)
    state_colors = _normalize_label_color_map(
        state_order,
        colors=state_colors,
        cmap_name="tab20",
    )
    # Pagination: compute how many panels fit on one A4 page
    panels_per_page, ncols_eff, max_rows = _panels_per_a4_page(grid_ncols)

    relative_by_group = {}
    overall_by_group = {}
    count_by_group = {}
    overall_count_by_group = {}
    group_label_title = ""
    unique_vals_per_col = {}
    if len(valid_group_cols) > 0:
        group_label_title = ", ".join(valid_group_cols)
        relative_by_group, overall_by_group = _compute_grouped_relative_matrices(
            df,
            group_cols=valid_group_cols,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )
        count_by_group, overall_count_by_group = _compute_grouped_count_matrices(
            df,
            group_cols=valid_group_cols,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )
        if len(valid_group_cols) in (1, 2):
            for col in valid_group_cols:
                unique_vals_per_col[col] = sorted(
                    df[col].astype(str).dropna().unique().tolist()
                )

    _grp_count_maxes = [
        float(mat.sum(axis=1).max())
        for mat in count_by_group.values()
        if mat is not None and len(mat) > 0 and len(mat.sum(axis=1)) > 0
    ]
    global_group_ymax = max(_grp_count_maxes) * 1.1 if _grp_count_maxes else 1.0

    with PdfPages(merged_pdf_path) as pdf:
        # --- Grouped pages (when group_cols provided) ---
        if len(valid_group_cols) > 0:
            if len(valid_group_cols) in (1, 2):
                # 2D grid layout with pagination on the Y (row) axis
                if len(valid_group_cols) == 2:
                    row_col = axis_cols[0] if axis_cols is not None else max(
                        valid_group_cols, key=lambda c: len(unique_vals_per_col[c])
                    )
                    col_y_vals = unique_vals_per_col[row_col]
                else:
                    col_y_vals = unique_vals_per_col[valid_group_cols[0]]

                for row_start in range(0, max(1, len(col_y_vals)), max_rows):
                    row_end = min(row_start + max_rows, len(col_y_vals))
                    _rs = (row_start, row_end)

                    fig_g = _plot_page_grouped_2d_grid(
                        relative_by_group,
                        overall_by_group=overall_by_group,
                        group_cols=valid_group_cols,
                        unique_vals_per_col=unique_vals_per_col,
                        state_order=state_order,
                        time_col=time_col,
                        group_label_title=group_label_title,
                        state_colors=state_colors,
                        row_slice=_rs,
                        axis_cols=axis_cols,
                        sample_title_fontsize=sample_title_fontsize,
                        sample_title_pad=sample_title_pad,
                    )
                    pdf.savefig(fig_g, dpi=dpi)
                    plt.close(fig_g)

                    fig_g_cnt = _plot_page_grouped_2d_grid(
                        count_by_group,
                        overall_by_group=overall_count_by_group,
                        group_cols=valid_group_cols,
                        unique_vals_per_col=unique_vals_per_col,
                        state_order=state_order,
                        time_col=time_col,
                        group_label_title=group_label_title,
                        state_colors=state_colors,
                        mode="count",
                        global_ymax=global_group_ymax,
                        row_slice=_rs,
                        axis_cols=axis_cols,
                        sample_title_fontsize=sample_title_fontsize,
                        sample_title_pad=sample_title_pad,
                    )
                    pdf.savefig(fig_g_cnt, dpi=dpi)
                    plt.close(fig_g_cnt)

            else:
                # 3+ group_cols: flat grid, paginated
                group_keys = list(relative_by_group.keys())
                for page_groups in _chunk_list(group_keys, panels_per_page):
                    page_rel = {k: relative_by_group[k] for k in page_groups}
                    page_overall_g = {k: overall_by_group[k] for k in page_groups}
                    fig_flat = _plot_page_grouped_stacked_grid(
                        page_rel,
                        overall_by_group=page_overall_g,
                        state_order=state_order,
                        time_col=time_col,
                        group_label_title=group_label_title,
                        grid_ncols=ncols_eff,
                        figsize_per_panel=figsize_per_panel,
                        state_colors=state_colors,
                        sample_title_fontsize=sample_title_fontsize,
                        sample_title_pad=sample_title_pad,
                    )
                    pdf.savefig(fig_flat, dpi=dpi)
                    plt.close(fig_flat)

                for page_groups in _chunk_list(group_keys, panels_per_page):
                    page_cnt = {k: count_by_group[k] for k in page_groups}
                    page_overall_cnt = {k: overall_count_by_group[k] for k in page_groups}
                    fig_flat_cnt = _plot_page_grouped_count_stacked_grid(
                        page_cnt,
                        overall_count_by_group=page_overall_cnt,
                        state_order=state_order,
                        time_col=time_col,
                        group_label_title=group_label_title,
                        grid_ncols=ncols_eff,
                        figsize_per_panel=figsize_per_panel,
                        state_colors=state_colors,
                        global_ymax=global_group_ymax,
                        sample_title_fontsize=sample_title_fontsize,
                        sample_title_pad=sample_title_pad,
                    )
                    pdf.savefig(fig_flat_cnt, dpi=dpi)
                    plt.close(fig_flat_cnt)

            # --- Grouped overall relative composition — horizontal bars ---
            group_row_order = list(overall_by_group.keys())
            rows_per_page = stacked_proportion_barh_rows_per_page()
            for page_groups in _chunk_list(group_row_order, rows_per_page):
                page_overall_g = {g: overall_by_group[g] for g in page_groups}
                fig_group_overall_h = plot_page_stacked_proportion_barh_grid(
                    page_overall_g,
                    row_order=page_groups,
                    class_order=state_order,
                    colors=state_colors,
                    title=f"Grouped Overall Relative State Composition — {group_label_title}",
                    row_label_fontsize=sample_title_fontsize,
                )
                pdf.savefig(fig_group_overall_h, dpi=dpi)
                plt.close(fig_group_overall_h)

        # --- Per-sample relative stacked grid (paginated) ---
        for page_samples in _paginate_samples(sample_order, samples_per_page=panels_per_page):
            page_rel = {s: relative_by_sample[s] for s in page_samples if s in relative_by_sample}
            page_overall = {s: overall_by_sample[s] for s in page_samples if s in overall_by_sample}
            fig1 = _plot_page_relative_stacked_grid(
                page_rel,
                overall_by_sample=page_overall,
                state_order=state_order,
                time_col=time_col,
                sample_col=sample_col,
                grid_ncols=ncols_eff,
                figsize_per_panel=figsize_per_panel,
                state_colors=state_colors,
                sample_title_fontsize=sample_title_fontsize,
                sample_title_pad=sample_title_pad,
            )
            pdf.savefig(fig1, dpi=dpi)
            plt.close(fig1)

        # --- Per-sample absolute count stacked grid (paginated) ---
        for page_samples in _paginate_samples(sample_order, samples_per_page=panels_per_page):
            page_cnt = {s: count_by_sample[s] for s in page_samples if s in count_by_sample}
            page_overall_cnt = {s: overall_count_by_sample[s] for s in page_samples if s in overall_count_by_sample}
            fig_counts = _plot_page_count_stacked_grid(
                page_cnt,
                overall_count_by_sample=page_overall_cnt,
                state_order=state_order,
                time_col=time_col,
                sample_col=sample_col,
                grid_ncols=ncols_eff,
                figsize_per_panel=figsize_per_panel,
                state_colors=state_colors,
                sample_title_fontsize=sample_title_fontsize,
                sample_title_pad=sample_title_pad,
            )
            pdf.savefig(fig_counts, dpi=dpi)
            plt.close(fig_counts)

        # --- Overall relative composition — horizontal bars, all samples (paginated) ---
        rows_per_page = stacked_proportion_barh_rows_per_page()
        for page_samples in _chunk_list(sample_order, rows_per_page):
            page_overall = {s: overall_by_sample[s] for s in page_samples if s in overall_by_sample}
            fig_overall_h = plot_page_stacked_proportion_barh_grid(
                page_overall,
                row_order=page_samples,
                class_order=state_order,
                colors=state_colors,
                title="Overall Relative State Composition — All Samples",
                row_label_fontsize=sample_title_fontsize,
            )
            pdf.savefig(fig_overall_h, dpi=dpi)
            plt.close(fig_overall_h)

    tables_to_concat = [
        _build_relative_plot_data_table(relative_by_sample, plot_view="stacked_by_sample").assign(
            plot_component="timecourse"
        ),
        _build_overall_summary_plot_data_table(overall_by_sample, plot_view="stacked_by_sample").assign(
            plot_component="overall_summary"
        ),
    ]
    # The PDF's grouped pages come from relative_by_group/overall_by_group, but
    # the CSV previously only ever included the plain per-sample tables above -
    # a group_cols selection changed what the PDF showed with no matching
    # breakdown in the CSV at all. Add the same grouped data here, with each
    # requested group column broken back out (via the group label each row was
    # built from) so it can be filtered/pivoted directly.
    if len(valid_group_cols) > 0:
        group_label_lookup = (
            df.assign(_group_label=_make_group_label(df, valid_group_cols).values)
            .drop_duplicates(subset=["_group_label"])
            .set_index("_group_label")[valid_group_cols]
        )
        group_timecourse = _build_relative_plot_data_table(
            relative_by_group, plot_view="stacked_by_group", label_col_name="group_label"
        ).assign(plot_component="group_timecourse")
        group_overall = _build_overall_summary_plot_data_table(
            overall_by_group, plot_view="stacked_by_group", label_col_name="group_label"
        ).assign(plot_component="group_overall_summary")
        for gtable in (group_timecourse, group_overall):
            if len(gtable) > 0:
                for col in valid_group_cols:
                    gtable[col] = gtable["group_label"].map(group_label_lookup[col])
        tables_to_concat.extend([group_timecourse, group_overall])

    merged_plot_data = pd.concat(tables_to_concat, ignore_index=True)
    merged_plot_data.to_csv(merged_plot_csv_path, index=False)

    if verbose:
        print(f"Saved state composition report PDF: {merged_pdf_path}")
        print(f"Saved state composition report CSV: {merged_plot_csv_path}")
        print(f"Saved state composition AUC CSV: {output_auc_csv_path}")

    return {
        "pdf_path": str(merged_pdf_path),
        "csv_path": str(output_auc_csv_path),
        "plot_data_csv_path": str(merged_plot_csv_path),
        "pdf_paths": {
            "combined": str(merged_pdf_path),
        },
        "plot_data_csv_paths": {
            "combined": str(merged_plot_csv_path),
        },
        "auc_table": auc_table,
        "plot_data_table": merged_plot_data,
        "relative_by_sample": relative_by_sample,
        "relative_pooled": relative_pooled if include_pooled_summary else None,
    }


# df_fig, fig, ax = plot_state_composition_over_time(
#     adata_full, 
#     time_col="position_t", 
#     state_col="ClusterID", 
#     relative=False
#     )

# df_fig, fig, axes= plot_state_composition_over_time(
#     adata_full, 
#     time_col="position_t", 
#     state_col="ClusterID", 
#     relative=True,
#     group_by_sample=True
#     )

# fig, axes = plot_state_composition_over_time(adata, group_by_sample=True, relative=True)
# plt.show()
