import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

from behav3d.analysis.behavior.state.utils import _normalize_label_color_map


A4_PORTRAIT = (8.27, 11.69)


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

        for s in mat.columns:
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
        ax.legend(
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

    observed_states = df[state_col].value_counts().index.tolist()
    if state_order is None:
        state_order = [str(s) for s in observed_states]
    else:
        state_order = [str(s) for s in list(state_order)]
        extras = [s for s in observed_states if str(s) not in state_order]
        state_order.extend([str(s) for s in extras])

    sample_order = df[sample_col].drop_duplicates().astype(str).tolist()
    return df, state_order, sample_order


def _compute_relative_matrix(panel_df, *, time_col, state_col, state_order):
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
        relative_by_sample[str(sample)] = _compute_relative_matrix(
            panel,
            time_col=time_col,
            state_col=state_col,
            state_order=state_order,
        )

    relative_pooled = _compute_relative_matrix(
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


def _build_color_map(labels, cmap_name="tab20"):
    """Create deterministic label->color mapping."""
    labels = [str(x) for x in list(labels)]
    if len(labels) == 0:
        return {}
    cmap = plt.get_cmap(cmap_name)
    if len(labels) <= getattr(cmap, "N", 256):
        color_values = [cmap(i / max(len(labels) - 1, 1)) for i in range(len(labels))]
    else:
        hsv = plt.get_cmap("hsv")
        color_values = [hsv(i / len(labels)) for i in range(len(labels))]
    return {lab: color_values[i] for i, lab in enumerate(labels)}


def _make_panel_grid(n_panels, ncols, figsize_per_panel):
    """Create a dynamic subplot grid and return flattened axes."""
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(float(n_panels) / float(ncols)))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=False,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
    )
    axes = np.atleast_1d(axes).ravel()
    for i in range(n_panels, len(axes)):
        axes[i].axis("off")
    return fig, axes


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


def _plot_overall_count_bar(ax, overall_counts, *, state_order, state_colors):
    """Draw a narrow vertical stacked bar showing raw cell counts."""
    total = float(sum(int(overall_counts.get(s, 0)) for s in state_order))
    ylim_top = max(total * 1.1, 1.0)
    bottom = 0.0
    for state in state_order:
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
    ax.set_xlabel("Overall")


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
    """New page: absolute cell count stacked composition per sample (grid)."""
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

    fig = plt.figure(figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
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
                f"{sample_col}: {sample} (empty)",
                fontsize=sample_title_fontsize,
                pad=sample_title_pad,
            )
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in state_order]
        colors = [state_colors[state] for state in state_order]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in state_order], colors=colors, alpha=0.9)
        ax.set_ylim(0.0, global_ymax)
        ax.set_title(
            f"{sample_col}: {sample}",
            fontsize=sample_title_fontsize,
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col)
        ax.set_ylabel("Cell count")
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

    fig.subplots_adjust(hspace=0.35)
    if len(first_labels) > 0:
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=min(len(first_labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle("Absolute Cell Count (Stacked) by Sample", y=0.998, fontsize=12)
    return fig


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
            auc = float(np.trapz(y=y, x=x)) if n_timepoints >= 2 else 0.0
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


def _build_relative_plot_data_table(relative_by_sample, plot_view=None):
    """Build long-form CSV table for all relative curves used in the report pages."""
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
                        "sample_name": str(sample_name),
                        "time": float(t),
                        "state_id": str(state_id),
                        "relative_proportion": float(v),
                    }
                )
    out = pd.DataFrame(rows)
    if plot_view is not None:
        out["plot_view"] = str(plot_view)
    return out


def _build_overall_summary_plot_data_table(overall_by_sample, plot_view=None):
    """Build long-form CSV table for pooled per-sample summary bars."""
    rows = []
    for sample_name, ser in overall_by_sample.items():
        if ser is None or len(ser) == 0:
            continue
        for state_id, value in ser.items():
            rows.append(
                {
                    "sample_name": str(sample_name),
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


def _compute_state_panel_ylims(
    relative_by_sample,
    *,
    state_order,
    headroom=0.10,
    min_positive_ylim=0.05,
):
    """Compute per-state y-limits from the largest observed relative proportion."""
    ylims = {}
    for state in state_order:
        max_value = 0.0
        for mat in relative_by_sample.values():
            if mat is None or len(mat) == 0 or state not in mat.columns:
                continue
            state_max = float(np.nanmax(mat[state].to_numpy(dtype=float))) if len(mat.index) > 0 else 0.0
            if np.isfinite(state_max):
                max_value = max(max_value, state_max)
        if max_value <= 0.0:
            ylims[str(state)] = (0.0, 1.0)
            continue
        y_max = min(1.0, max(max_value * (1.0 + float(headroom)), float(min_positive_ylim)))
        ylims[str(state)] = (0.0, float(y_max))
    return ylims


def _plot_overall_summary_bar(ax, overall, *, state_order, state_colors, show_ylabel=False):
    """Draw a narrow vertical stacked bar for pooled sample composition."""
    bottom = 0.0
    for state in state_order:
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
        ax.set_ylabel("Proportion")
    else:
        ax.set_yticklabels([])
    ax.set_xlabel("Overall")


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
    """Page 1: relative stacked composition per sample (grid)."""
    samples = list(relative_by_sample.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(samples)) / float(ncols))))
    fig = plt.figure(figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    plot_axes = []
    first_handles = []
    first_labels = []

    for i, sample in enumerate(samples):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1], sharey=ax)
        plot_axes.append(ax)
        mat = relative_by_sample[sample]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(
                f"{sample_col}: {sample} (empty)",
                fontsize=sample_title_fontsize,
                pad=sample_title_pad,
            )
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in state_order]
        colors = [state_colors[state] for state in state_order]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in state_order], colors=colors, alpha=0.9)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            f"{sample_col}: {sample}",
            fontsize=sample_title_fontsize,
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion")
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

    fig.subplots_adjust(hspace=0.35)
    if len(first_labels) > 0:
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=min(len(first_labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle("Relative State Composition (Stacked) by Sample", y=0.998, fontsize=12)
    return fig


def _plot_page_nonstacked_state_lines_grid(
    relative_by_sample,
    *,
    state_order,
    time_col,
    sample_col,
    grid_ncols,
    figsize_per_panel,
    state_colors,
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """Page 2: non-stacked line plots by sample (one line per state)."""
    samples = list(relative_by_sample.keys())
    fig, axes = _make_panel_grid(len(samples), grid_ncols, figsize_per_panel)
    for i, sample in enumerate(samples):
        ax = axes[i]
        mat = relative_by_sample[sample]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(
                f"{sample_col}: {sample} (empty)",
                fontsize=sample_title_fontsize,
                pad=sample_title_pad,
            )
            ax.axis("off")
            continue
        for state in state_order:
            ax.plot(
                x,
                mat[state].to_numpy(dtype=float),
                color=state_colors[state],
                linewidth=1.4,
                alpha=0.9,
                label=str(state),
            )
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            f"{sample_col}: {sample}",
            fontsize=sample_title_fontsize,
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion")

    fig.subplots_adjust(hspace=0.35)
    handles, labels = axes[0].get_legend_handles_labels() if len(samples) > 0 else ([], [])
    if len(labels) > 0:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=min(len(labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle("Relative State Composition (Non-Stacked Lines) by Sample", y=0.998, fontsize=12)
    return fig


def _plot_nonstacked_state_lines_page(
    relative_by_sample,
    *,
    overall_by_sample,
    state_order,
    time_col,
    sample_col,
    state_colors,
    sample_title_fontsize=8,
    sample_title_pad=2,
    page_size=A4_PORTRAIT,
):
    """Render one portrait PDF page with one sample panel per row."""
    samples = list(relative_by_sample.keys())
    nrows = max(len(samples), 1)
    fig = plt.figure(figsize=page_size)
    outer = GridSpec(nrows=nrows, ncols=1, figure=fig)
    plot_axes = []
    first_handles = []
    first_labels = []

    for i, sample in enumerate(samples):
        inner = outer[i].subgridspec(1, 2, width_ratios=[8, 1], wspace=0.10)
        ax = fig.add_subplot(inner[0, 0])
        ax_bar = fig.add_subplot(inner[0, 1], sharey=ax)
        plot_axes.append(ax)
        mat = relative_by_sample[sample]
        x = mat.index.to_numpy(dtype=float)
        if len(x) == 0:
            ax.set_title(
                f"{sample_col}: {sample} (empty)",
                fontsize=sample_title_fontsize,
                pad=sample_title_pad,
            )
            ax.axis("off")
            ax_bar.axis("off")
            continue
        for state in state_order:
            ax.plot(
                x,
                mat[state].to_numpy(dtype=float),
                color=state_colors[state],
                linewidth=1.4,
                alpha=0.9,
                label=str(state),
            )
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            f"{sample_col}: {sample}",
            fontsize=sample_title_fontsize,
            pad=sample_title_pad,
        )
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion")
        _plot_overall_summary_bar(
            ax_bar,
            overall_by_sample[sample],
            state_order=state_order,
            state_colors=state_colors,
        )
        if len(first_labels) == 0:
            first_handles, first_labels = ax.get_legend_handles_labels()

    if len(first_labels) > 0:
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=min(len(first_labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 0.98))
    else:
        fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.suptitle("Relative State Composition (Non-Stacked Lines) by Sample", y=0.995, fontsize=12)
    return fig


def _plot_page_per_state_sample_lines_grid(
    relative_by_sample,
    *,
    state_order,
    time_col,
    sample_col,
    grid_ncols,
    figsize_per_panel,
    sample_colors,
    state_ylims=None,
):
    """Page 3: one panel per state, with sample-colored lines."""
    fig, axes = _make_panel_grid(len(state_order), grid_ncols, figsize_per_panel)
    samples = list(relative_by_sample.keys())
    for i, state in enumerate(state_order):
        ax = axes[i]
        for sample in samples:
            mat = relative_by_sample[sample]
            x = mat.index.to_numpy(dtype=float)
            if len(x) == 0:
                continue
            ax.plot(
                x,
                mat[state].to_numpy(dtype=float),
                color=sample_colors[sample],
                linewidth=1.4,
                alpha=0.9,
                label=str(sample),
            )
        if isinstance(state_ylims, dict) and str(state) in state_ylims:
            ax.set_ylim(*state_ylims[str(state)])
        else:
            ax.set_ylim(0.0, 1.0)
        ax.set_title(f"{state}")
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion")

    handles, labels = axes[0].get_legend_handles_labels() if len(state_order) > 0 else ([], [])
    if len(labels) > 0:
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=min(len(labels), 8),
            frameon=False,
            title=sample_col,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle("Per-State Relative Composition with Sample-Colored Lines", y=0.998, fontsize=12)
    return fig


def _make_group_label(df, group_cols):
    """Return a Series of concatenated group-column values (index aligned with df)."""
    parts = []
    for col in group_cols:
        if col in df.columns:
            parts.append(df[col].fillna("(unknown)").astype(str))
        else:
            parts.append(pd.Series(["(unknown)"] * len(df), index=df.index))
    if len(parts) == 0:
        return pd.Series(["(all)"] * len(df), index=df.index)
    result = parts[0]
    for p in parts[1:]:
        result = result + " | " + p
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
        relative_by_group[str(grp)] = _compute_relative_matrix(
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
    """One page: relative stacked composition per group, same layout as per-sample page."""
    groups = list(relative_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))
    fig = plt.figure(figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
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
            ax.set_title(f"{grp} (empty)", fontsize=sample_title_fontsize, pad=sample_title_pad)
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in state_order]
        colors = [state_colors[state] for state in state_order]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in state_order], colors=colors, alpha=0.9)
        ax.set_ylim(0.0, 1.0)
        ax.set_title(grp, fontsize=sample_title_fontsize, pad=sample_title_pad)
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion")
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

    fig.subplots_adjust(hspace=0.35)
    if len(first_labels) > 0:
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=min(len(first_labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle(
        f"Grouped Relative State Composition — grouped by: {group_label_title}",
        y=0.998,
        fontsize=12,
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
    sample_title_fontsize=8,
    sample_title_pad=2,
):
    """One page: absolute cell count stacked composition per group."""
    groups = list(count_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))

    all_maxes = []
    for mat in count_by_group.values():
        if mat is not None and len(mat) > 0:
            row_totals = mat.sum(axis=1)
            if len(row_totals) > 0:
                all_maxes.append(float(row_totals.max()))
    global_ymax = max(all_maxes) * 1.1 if all_maxes else 1.0

    fig = plt.figure(figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
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
            ax.set_title(f"{grp} (empty)", fontsize=sample_title_fontsize, pad=sample_title_pad)
            ax.axis("off")
            ax_bar.axis("off")
            continue
        y_arrays = [mat[state].to_numpy(dtype=float) for state in state_order]
        colors = [state_colors[state] for state in state_order]
        ax.stackplot(x, *y_arrays, labels=[str(s) for s in state_order], colors=colors, alpha=0.9)
        ax.set_ylim(0.0, global_ymax)
        ax.set_title(grp, fontsize=sample_title_fontsize, pad=sample_title_pad)
        ax.set_xlabel(time_col)
        ax.set_ylabel("Cell count")
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

    fig.subplots_adjust(hspace=0.35)
    if len(first_labels) > 0:
        fig.legend(
            first_handles,
            first_labels,
            loc="lower center",
            ncol=min(len(first_labels), 8),
            frameon=False,
        )
        fig.tight_layout(rect=(0, 0.06, 1, 1))
    else:
        fig.tight_layout()
    fig.suptitle(
        f"Grouped Absolute Cell Count — grouped by: {group_label_title}",
        y=0.998,
        fontsize=12,
    )
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
):
    """
    Save a combined multi-page PDF report + merged plot-data CSV for relative state composition.

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

    df, state_order, sample_order = _prepare_state_composition_df(
        adata,
        time_col=time_col,
        state_col=state_col,
        sample_col=sample_col,
        state_order=state_order,
    )

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

    valid_group_cols = []
    if group_cols:
        obs_cols_available = list(adata.obs.columns)
        for gc in group_cols:
            if gc in obs_cols_available:
                valid_group_cols.append(gc)
                df[gc] = adata.obs.loc[df.index, gc].astype(str).fillna("(unknown)").values
            else:
                if verbose:
                    print(f"  Warning: group_col '{gc}' not found in adata.obs — skipping.")

    auc_table = _compute_relative_auc_table(
        relative_by_sample,
        state_order=state_order,
        include_pooled_summary=bool(include_pooled_summary),
        relative_pooled=relative_pooled,
    )
    auc_table.to_csv(output_auc_csv_path, index=False)

    state_colors = _normalize_label_color_map(
        state_order,
        colors=state_colors,
        cmap_name="tab20",
    )
    sample_colors = _build_color_map(sample_order, cmap_name="tab20")
    nonstacked_samples_per_page = 4
    per_state_ylim_headroom = 0.10

    nonstacked_pages = _paginate_samples(sample_order, samples_per_page=nonstacked_samples_per_page)

    relative_by_group = {}
    overall_by_group = {}
    count_by_group = {}
    overall_count_by_group = {}
    group_label_title = ""
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

    with PdfPages(merged_pdf_path) as pdf:
        # Grouped pages first (when group_cols provided)
        if len(valid_group_cols) > 0:
            fig_grouped = _plot_page_grouped_stacked_grid(
                relative_by_group,
                overall_by_group=overall_by_group,
                state_order=state_order,
                time_col=time_col,
                group_label_title=group_label_title,
                grid_ncols=grid_ncols,
                figsize_per_panel=figsize_per_panel,
                state_colors=state_colors,
                sample_title_fontsize=sample_title_fontsize,
                sample_title_pad=sample_title_pad,
            )
            pdf.savefig(fig_grouped, dpi=dpi, bbox_inches="tight")
            plt.close(fig_grouped)

            fig_grouped_counts = _plot_page_grouped_count_stacked_grid(
                count_by_group,
                overall_count_by_group=overall_count_by_group,
                state_order=state_order,
                time_col=time_col,
                group_label_title=group_label_title,
                grid_ncols=grid_ncols,
                figsize_per_panel=figsize_per_panel,
                state_colors=state_colors,
                sample_title_fontsize=sample_title_fontsize,
                sample_title_pad=sample_title_pad,
            )
            pdf.savefig(fig_grouped_counts, dpi=dpi, bbox_inches="tight")
            plt.close(fig_grouped_counts)

        # Per-sample relative stacked grid
        fig1 = _plot_page_relative_stacked_grid(
            relative_by_sample,
            overall_by_sample=overall_by_sample,
            state_order=state_order,
            time_col=time_col,
            sample_col=sample_col,
            grid_ncols=grid_ncols,
            figsize_per_panel=figsize_per_panel,
            state_colors=state_colors,
            sample_title_fontsize=sample_title_fontsize,
            sample_title_pad=sample_title_pad,
        )
        pdf.savefig(fig1, dpi=dpi, bbox_inches="tight")
        plt.close(fig1)

        # Per-sample absolute count stacked grid
        fig_counts = _plot_page_count_stacked_grid(
            count_by_sample,
            overall_count_by_sample=overall_count_by_sample,
            state_order=state_order,
            time_col=time_col,
            sample_col=sample_col,
            grid_ncols=grid_ncols,
            figsize_per_panel=figsize_per_panel,
            state_colors=state_colors,
            sample_title_fontsize=sample_title_fontsize,
            sample_title_pad=sample_title_pad,
        )
        pdf.savefig(fig_counts, dpi=dpi, bbox_inches="tight")
        plt.close(fig_counts)

        # Non-stacked per-sample line plots
        for page_samples in nonstacked_pages:
            page_relative = {
                sample: relative_by_sample[sample]
                for sample in page_samples
                if sample in relative_by_sample
            }
            page_overall = {
                sample: overall_by_sample[sample]
                for sample in page_samples
                if sample in overall_by_sample
            }
            fig2 = _plot_nonstacked_state_lines_page(
                page_relative,
                overall_by_sample=page_overall,
                state_order=state_order,
                time_col=time_col,
                sample_col=sample_col,
                state_colors=state_colors,
                sample_title_fontsize=sample_title_fontsize,
                sample_title_pad=sample_title_pad,
            )
            pdf.savefig(fig2, dpi=dpi, bbox_inches="tight")
            plt.close(fig2)

        # Last page: per-state lines colored by sample
        state_ylims = _compute_state_panel_ylims(
            relative_by_sample,
            state_order=state_order,
            headroom=per_state_ylim_headroom,
        )
        fig3 = _plot_page_per_state_sample_lines_grid(
            relative_by_sample,
            state_order=state_order,
            time_col=time_col,
            sample_col=sample_col,
            grid_ncols=grid_ncols,
            figsize_per_panel=figsize_per_panel,
            sample_colors=sample_colors,
            state_ylims=state_ylims,
        )
        pdf.savefig(fig3, dpi=dpi, bbox_inches="tight")
        plt.close(fig3)

    merged_plot_data = pd.concat(
        [
            _build_relative_plot_data_table(relative_by_sample, plot_view="stacked_by_sample").assign(
                plot_component="timecourse"
            ),
            _build_overall_summary_plot_data_table(overall_by_sample, plot_view="stacked_by_sample").assign(
                plot_component="overall_summary"
            ),
            _build_relative_plot_data_table(relative_by_sample, plot_view="nonstacked_by_sample").assign(
                plot_component="timecourse"
            ),
            _build_overall_summary_plot_data_table(overall_by_sample, plot_view="nonstacked_by_sample").assign(
                plot_component="overall_summary"
            ),
            _build_relative_plot_data_table(relative_by_sample, plot_view="per_state_by_sample").assign(
                plot_component="timecourse"
            ),
        ],
        ignore_index=True,
    )
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
