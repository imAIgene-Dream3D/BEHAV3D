import hashlib
from contextlib import nullcontext
from itertools import combinations

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy import stats

A4_PORTRAIT = (8.27, 11.69)


def _resolve_effective_group_cols(group_cols, group_x=None, group_y=None):
    """Combine explicit X/Y axis choices with a "group per page" column list.

    Returns ``(effective_group_cols, axis_cols)``. ``axis_cols`` is a
    ``(col_y, col_x)`` tuple that forces explicit 2D-grid axis assignment when
    the grid is driven purely by ``group_x``/``group_y`` (no additional
    page-grouping columns selected); otherwise it is ``None`` and the caller
    falls back to automatic cardinality-based axis assignment.
    """
    xy = [c for c in (group_x, group_y) if c]
    page_cols = list(group_cols) if group_cols else []
    effective_group_cols = list(dict.fromkeys(xy + page_cols))
    axis_cols = (group_y, group_x) if (group_x and group_y and not page_cols) else None
    return effective_group_cols, axis_cols


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


def hash_stable_label_color(label, cmap_name="tab20"):
    """Deterministic per-label color independent of sort order or subset filtering."""
    cmap = plt.get_cmap(cmap_name)
    n = max(int(getattr(cmap, "N", 20)), 1)
    idx = int(hashlib.md5(str(label).encode("utf-8")).hexdigest(), 16) % n
    return str(mcolors.to_hex(cmap(idx / max(n - 1, 1))))


def hash_stable_label_color_map(labels, colors=None, cmap_name="tab20"):
    """Like _normalize_label_color_map, but the default fill is hash-stable, not position-based."""
    labels = [str(x) for x in list(labels or [])]
    saved = {str(k): v for k, v in dict(colors or {}).items()} if colors else {}
    return {label: saved.get(label, hash_stable_label_color(label, cmap_name=cmap_name)) for label in labels}


def draw_thin_stacked_proportion_barh(
    ax,
    values,
    class_order,
    colors,
    *,
    xmax=1.0,
    bar_height_frac=0.22,
):
    """Draw one thin horizontal stacked bar of class proportions on `ax`, always summing to `xmax`.

    Caller is expected to set a title above the axes (e.g. ax.set_title(..., pad=1)) for the
    row/group label.
    """
    half = float(bar_height_frac) / 2.0
    left = 0.0
    for class_name in class_order:
        val = float(values.get(class_name, 0.0))
        if val <= 0.0:
            continue
        ax.barh(
            [0.0], [val], left=left, color=colors[class_name],
            height=bar_height_frac, edgecolor="none", linewidth=0.0,
        )
        left += val
    ax.set_xlim(0.0, xmax if xmax and xmax > 0 else 1.0)
    ax.set_ylim(-half - 0.05, half + 0.05)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_stacked_proportion_barv(
    ax,
    props_df,
    class_order,
    stack_order,
    colors,
    *,
    ymax=1.0,
    bar_width=0.6,
    xtick_fontsize=7,
):
    """Draw a grid-cell panel of vertical stacked bars: one bar per `class_order` entry on the
    x-axis, each stacked bottom-up by `stack_order`, height summing to `ymax`.

    `props_df` is indexed by class (rows) with `stack_order` columns holding proportions (e.g.
    from `compute_class_by_stack_proportions`); missing (class, stack) combos default to 0.
    """
    x = np.arange(len(class_order))
    bottom = np.zeros(len(class_order))
    for stack_name in stack_order:
        if stack_name in props_df.columns:
            vals = props_df[stack_name].reindex(class_order).fillna(0.0).to_numpy(dtype=float)
        else:
            vals = np.zeros(len(class_order))
        ax.bar(x, vals, bottom=bottom, width=bar_width, color=colors[stack_name], edgecolor="none", linewidth=0.0)
        bottom += vals
    ax.set_xlim(-0.5, len(class_order) - 0.5)
    ax.set_ylim(0.0, ymax if ymax and ymax > 0 else 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([str(c) for c in class_order], rotation=45, ha="right", fontsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=xtick_fontsize)
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def compute_class_by_stack_proportions(df, *, class_col, stack_col, class_order, stack_order):
    """Crosstab of `class_col` x `stack_col`, normalized so each class's row sums to 1."""
    props = pd.crosstab(df[class_col], df[stack_col], normalize="index")
    return props.reindex(index=list(class_order), columns=list(stack_order), fill_value=0.0)


def _chunk_list(lst, n):
    """Split list into chunks of at most n items."""
    n = max(1, int(n))
    lst = list(lst)
    return [lst[i : i + n] for i in range(0, max(1, len(lst)), n)]


def _wrap_row_label(text, max_chars=60):
    """Wrap a long label at underscore boundaries, capped at 2 lines."""
    text = str(text)
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


def stacked_proportion_barh_rows_per_page(page_size=A4_PORTRAIT, row_height_in=0.20, base_height_in=1.1):
    """Recommended max rows per page before `plot_page_stacked_proportion_barh_grid`'s height caps out."""
    return max(1, int((page_size[1] - base_height_in) / row_height_in))


def plot_page_stacked_proportion_barh_grid(
    data_by_key,
    *,
    row_order,
    class_order,
    colors,
    title,
    page_size=A4_PORTRAIT,
    bar_height_frac=0.22,
    row_height_in=0.20,
    base_height_in=1.1,
    row_label_fontsize=8,
):
    """One page of thin stacked-proportion bars, one per row (sample/group).

    Figure height grows with the number of rows but is capped at ``page_size[1]`` so a
    dense page stays compact instead of growing indefinitely; page ``row_order`` with
    `_chunk_list` (sized via `stacked_proportion_barh_rows_per_page`) when there are more
    rows than comfortably fit on one page.
    """
    rows = [k for k in row_order if k in data_by_key]
    n_rows = max(1, len(rows))
    fig_h = min(page_size[1], row_height_in * n_rows + base_height_in)
    fig, axes = plt.subplots(nrows=n_rows, ncols=1, figsize=(page_size[0], fig_h), squeeze=False)
    axes = axes.flatten()

    for i, key in enumerate(rows):
        values = data_by_key.get(key)
        if values is None:
            values = {}
        draw_thin_stacked_proportion_barh(
            axes[i], values, class_order, colors, xmax=1.0, bar_height_frac=bar_height_frac,
        )
        axes[i].set_title(_wrap_row_label(key), fontsize=row_label_fontsize, pad=1, loc="left")

    handles = [Patch(facecolor=colors[c], label=str(c)) for c in class_order]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 8), frameon=False, fontsize=7)
    fig.tight_layout(rect=(0.03, 0.05, 1, 0.94), h_pad=0.3)
    fig.suptitle(title, y=0.97, fontsize=12, fontweight="bold")
    return fig


def welch_ttest_stars(p_value):
    """'' if p is NaN/>=0.05, else '*' p<0.05, '**' p<0.01, '***' p<0.001, '****' p<0.0001."""
    if p_value is None or not np.isfinite(p_value):
        return ""
    if p_value < 1e-4:
        return "****"
    if p_value < 1e-3:
        return "***"
    if p_value < 1e-2:
        return "**"
    if p_value < 5e-2:
        return "*"
    return ""


def _welch_diff_rows(class_order, vals_a_panel, vals_b_panel):
    """Per-class Welch's two-sided unpaired t-test rows comparing two panels of samples."""
    rows = []
    for class_name in class_order:
        vals_a = pd.to_numeric(vals_a_panel[class_name], errors="coerce").dropna().to_numpy()
        vals_b = pd.to_numeric(vals_b_panel[class_name], errors="coerce").dropna().to_numpy()
        mean_a = float(np.mean(vals_a)) if len(vals_a) > 0 else float("nan")
        mean_b = float(np.mean(vals_b)) if len(vals_b) > 0 else float("nan")
        if len(vals_a) >= 2 and len(vals_b) >= 2:
            t_stat, p_value = stats.ttest_ind(vals_b, vals_a, equal_var=False, nan_policy="omit")
            t_stat = float(t_stat)
            p_value = float(p_value)
        else:
            t_stat, p_value = float("nan"), float("nan")
        rows.append({
            "class": class_name,
            "mean_a": mean_a,
            "mean_b": mean_b,
            "diff": (mean_b - mean_a) if np.isfinite(mean_a) and np.isfinite(mean_b) else float("nan"),
            "t_stat": t_stat,
            "p_value": p_value,
            "stars": welch_ttest_stars(p_value),
            "n_a": int(len(vals_a)),
            "n_b": int(len(vals_b)),
        })
    return rows


def compute_condition_diff_stats_pairwise(
    per_sample_class_props,
    sample_metadata,
    *,
    class_order,
    condition_col,
    group_cols=None,
):
    """Per group (or a single '(all)' bucket if group_cols is empty), per pair of
    condition_col levels present anywhere in the dataset — every combination, ordered by
    first appearance in the data, e.g. levels [None, M21, M23] -> (None, M21), (None,
    M23), (M21, M23) — per class: Welch's two-sided unpaired t-test. The same set of
    pairs is used for every group so rows stay aligned in the resulting grid; a pair
    with too few samples on one side within a given group simply gets NaN stats for it.

    diff = mean of the second level minus mean of the first level in each pair.

    Parameters
    ----------
    per_sample_class_props : pd.DataFrame
        index=sample, columns=class_order, values=proportion (0-1).
    sample_metadata : pd.DataFrame
        index=sample, columns include condition_col (+ group_cols if given).

    Returns
    -------
    dict[str, dict[tuple[str, str], pd.DataFrame]]
        Keyed by group label (or "(all)"), then by (level_a, level_b); each DataFrame
        has columns: class, mean_a, mean_b, diff, t_stat, p_value, stars, n_a, n_b.
    """
    class_order = [str(c) for c in list(class_order)]
    joined = per_sample_class_props.join(sample_metadata, how="inner")

    # Excluding condition_col here too (not just at the caller) guards against
    # any other caller passing an overlapping group_cols/condition_col combo -
    # joined[condition_col] would otherwise be ambiguous if condition_col were
    # duplicated among the selected columns, since DataFrame.unique() doesn't exist.
    valid_group_cols = [c for c in (group_cols or []) if c in joined.columns and c != condition_col]
    if valid_group_cols:
        group_labels = _make_group_label(joined, valid_group_cols).astype(str)
        groups = {label: joined[group_labels == label] for label in group_labels.unique().tolist()}
    else:
        groups = {"(all)": joined}

    all_levels = sample_metadata[condition_col].astype(str).drop_duplicates().tolist()

    out = {}
    for group_label, panel in groups.items():
        cond = panel[condition_col].astype(str)
        pair_stats = {}
        for level_a, level_b in combinations(all_levels, 2):
            vals_a_panel = panel[cond == level_a]
            vals_b_panel = panel[cond == level_b]
            rows = _welch_diff_rows(class_order, vals_a_panel, vals_b_panel)
            pair_stats[(level_a, level_b)] = pd.DataFrame(rows)
        out[str(group_label)] = pair_stats
    return out


def draw_diff_barh(
    ax,
    class_order,
    diff_df,
    colors,
    *,
    xlim=None,
    bar_height_frac=0.6,
    star_fontsize=8,
    label_fontsize=8,
):
    """One horizontal bar per class showing signed proportion difference (diff_df['diff'] in
    fractional units, plotted as percent), with a significance star at the bar tip. class_order[0]
    is drawn at the top (matches typical reference-figure convention: first cluster on top)."""
    class_order = [str(c) for c in list(class_order)]
    by_class = diff_df.set_index("class")
    rows = list(reversed(class_order))
    y_positions = list(range(len(rows)))

    max_abs = 0.0
    for class_name in rows:
        diff_pct = float(by_class.loc[class_name, "diff"]) * 100.0 if class_name in by_class.index else 0.0
        max_abs = max(max_abs, abs(diff_pct))
    xmax = xlim if xlim is not None else max(max_abs * 1.35, 5.0)

    for y, class_name in zip(y_positions, rows):
        if class_name not in by_class.index:
            continue
        row = by_class.loc[class_name]
        diff_pct = float(row["diff"]) * 100.0 if np.isfinite(row["diff"]) else 0.0
        color = colors.get(class_name, "#808080")
        ax.barh([y], [diff_pct], height=bar_height_frac, color=color, edgecolor="none", linewidth=0.0)
        stars = str(row.get("stars", ""))
        if stars:
            tip = diff_pct + (xmax * 0.03 if diff_pct >= 0 else -xmax * 0.03)
            ha = "left" if diff_pct >= 0 else "right"
            ax.text(tip, y, stars, ha=ha, va="center", fontsize=star_fontsize)

    ax.axvline(0.0, color="black", linewidth=0.8)
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(rows, fontsize=label_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Cluster size difference (%)", fontsize=8)


def _diff_df_has_data(diff_df):
    """True if at least one class has a finite diff (i.e. both condition sides had samples)."""
    return diff_df is not None and bool(np.isfinite(diff_df["diff"].to_numpy(dtype=float)).any())


def plot_condition_diff_grid(
    diff_stats_by_group,
    *,
    class_order,
    colors,
    title,
    out_pdf,
    out_csv=None,
    figsize_per_panel=(3.2, None),
    max_cols_per_page=4,
    label_fontsize=9,
    pdf_pages=None,
):
    """Multi-page PDF: one row per pairwise level comparison (the same set of pairs is shared
    across every group), one column per group (or a single column if diff_stats_by_group has
    only the '(all)' key). Saves a PDF and (optionally) a merged long-form CSV.

    If ``pdf_pages`` (an open ``matplotlib.backends.backend_pdf.PdfPages``) is given, pages are
    appended to it instead of a standalone ``out_pdf`` file being created; ``out_pdf`` is still
    used as the nominal path reported in the return value.
    """
    class_order = [str(c) for c in list(class_order)]
    group_labels = list(diff_stats_by_group.keys())
    pair_keys = list(next(iter(diff_stats_by_group.values()), {}).keys())
    n_pairs = max(len(pair_keys), 1)
    panel_w = float(figsize_per_panel[0])
    panel_h = figsize_per_panel[1] if figsize_per_panel[1] is not None else max(1.2 + 0.45 * len(class_order), 3.0)
    header_h = 0.35

    global_xmax = 0.0
    for pairs in diff_stats_by_group.values():
        for df in pairs.values():
            if len(df) > 0:
                global_xmax = max(global_xmax, float(np.nanmax(np.abs(df["diff"].to_numpy(dtype=float))) * 100.0))
    global_xmax = max(global_xmax * 1.35, 5.0)

    out_pdf = str(out_pdf)
    long_rows = []

    with (nullcontext(pdf_pages) if pdf_pages is not None else PdfPages(out_pdf)) as pdf:
        for page_groups in _chunk_list(group_labels, max(1, int(max_cols_per_page))):
            n_cols = len(page_groups)
            fig = plt.figure(figsize=(panel_w * n_cols, header_h + panel_h * n_pairs))
            outer = GridSpec(
                nrows=1 + n_pairs, ncols=n_cols, figure=fig,
                width_ratios=[panel_w] * n_cols,
                height_ratios=[header_h] + [1.0] * n_pairs,
                hspace=0.55, wspace=0.25,
            )

            for c, group_label in enumerate(page_groups):
                ax_ch = fig.add_subplot(outer[0, c])
                ax_ch.axis("off")
                if group_label != "(all)":
                    ax_ch.text(
                        0.5, 0.5, str(group_label),
                        ha="center", va="center", fontsize=label_fontsize,
                        fontweight="bold", transform=ax_ch.transAxes,
                    )

            for r, (level_a, level_b) in enumerate(pair_keys):
                for c, group_label in enumerate(page_groups):
                    ax = fig.add_subplot(outer[r + 1, c])
                    ax.set_title(
                        _wrap_row_label(f"{level_a}  vs  {level_b}"),
                        fontsize=max(label_fontsize - 2, 6), fontweight="normal", pad=4,
                    )
                    diff_df = diff_stats_by_group[group_label].get((level_a, level_b))
                    if not _diff_df_has_data(diff_df):
                        ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                        ax.axis("off")
                        continue
                    draw_diff_barh(ax, class_order, diff_df, colors, xlim=global_xmax)
                    if c > 0:
                        ax.set_yticklabels([])

                    df_long = diff_df.copy()
                    df_long["group"] = group_label
                    df_long["level_a"] = level_a
                    df_long["level_b"] = level_b
                    long_rows.append(df_long)

            fig.suptitle(title, fontsize=11, fontweight="bold", y=0.99)
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    csv_path = None
    if out_csv is not None:
        merged = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
        csv_path = str(out_csv)
        merged.to_csv(csv_path, index=False)

    return {"pdf_path": out_pdf, "csv_path": csv_path}


def plot_condition_diff_grid_2d(
    per_sample_class_props,
    sample_metadata,
    *,
    class_order,
    colors,
    condition_col,
    grid_axis_cols,
    extra_group_cols=None,
    title,
    out_pdf,
    out_csv=None,
    figsize_per_panel=(3.2, None),
    label_fontsize=9,
    pdf_pages=None,
):
    """Multi-page PDF: a true 2D grid of diff-bar panels — columns = ``grid_axis_cols[0]``,
    rows = ``grid_axis_cols[1]`` — for a binary ``condition_col`` (exactly one pairwise
    comparison, e.g. "contact" vs "no_contact").

    Unlike ``plot_condition_diff_grid`` (which always dedicates the row axis to the pairwise
    condition-level comparisons and pools every group column into a single column axis), this
    variant only makes sense when there's exactly one comparison to show per cell — freeing both
    ``grid_axis_cols`` to be true grid axes instead. ``extra_group_cols``, if given, paginate into
    one full grid page per unique combination rather than adding a third grid dimension.
    """
    class_order = [str(c) for c in list(class_order)]
    extra_group_cols = [c for c in (extra_group_cols or []) if c in sample_metadata.columns]
    grid_axis_cols = [c for c in grid_axis_cols if c in sample_metadata.columns]
    if len(grid_axis_cols) != 2:
        raise ValueError("plot_condition_diff_grid_2d requires exactly 2 grid_axis_cols (group_x and group_y).")
    col_x, col_y = grid_axis_cols[0], grid_axis_cols[1]

    panel_w = float(figsize_per_panel[0])
    panel_h = figsize_per_panel[1] if figsize_per_panel[1] is not None else max(1.2 + 0.45 * len(class_order), 3.0)
    header_w = 1.6
    header_h = 0.35

    out_pdf = str(out_pdf)
    long_rows = []

    if extra_group_cols:
        page_labels = _make_group_label(sample_metadata, extra_group_cols).astype(str)
        page_keys = sorted(page_labels.dropna().unique().tolist())
    else:
        page_labels = pd.Series(["(all)"] * len(sample_metadata), index=sample_metadata.index)
        page_keys = ["(all)"]

    with (nullcontext(pdf_pages) if pdf_pages is not None else PdfPages(out_pdf)) as pdf:
        for page_key in page_keys:
            page_metadata = sample_metadata[page_labels == page_key]
            if len(page_metadata) == 0:
                continue

            diff_stats_by_facet = compute_condition_diff_stats_pairwise(
                per_sample_class_props,
                page_metadata,
                class_order=class_order,
                condition_col=condition_col,
                group_cols=grid_axis_cols,
            )
            pair_keys = list(next(iter(diff_stats_by_facet.values()), {}).keys())
            if len(pair_keys) > 1:
                raise ValueError(
                    f"plot_condition_diff_grid_2d expects a binary condition_col='{condition_col}' "
                    f"(exactly one pairwise comparison); found {len(pair_keys)}."
                )
            level_a, level_b = pair_keys[0] if pair_keys else (None, None)

            col_x_vals = sorted(page_metadata[col_x].astype(str).dropna().unique().tolist())
            col_y_vals = sorted(page_metadata[col_y].astype(str).dropna().unique().tolist())
            nrows_data = max(1, len(col_y_vals))
            ncols_data = max(1, len(col_x_vals))

            global_xmax = 0.0
            if level_a is not None:
                for pairs in diff_stats_by_facet.values():
                    d = pairs.get((level_a, level_b))
                    if d is not None and len(d) > 0:
                        global_xmax = max(global_xmax, float(np.nanmax(np.abs(d["diff"].to_numpy(dtype=float))) * 100.0))
            global_xmax = max(global_xmax * 1.35, 5.0)

            fig = plt.figure(figsize=(header_w + panel_w * ncols_data, header_h + panel_h * nrows_data))
            outer = GridSpec(
                nrows=1 + nrows_data, ncols=1 + ncols_data, figure=fig,
                width_ratios=[header_w] + [panel_w] * ncols_data,
                height_ratios=[header_h] + [panel_h] * nrows_data,
                hspace=0.55, wspace=0.25,
            )
            fig.add_subplot(outer[0, 0]).axis("off")

            for c, col_x_val in enumerate(col_x_vals):
                ax_ch = fig.add_subplot(outer[0, c + 1])
                ax_ch.axis("off")
                ax_ch.text(
                    0.5, 0.5, str(col_x_val),
                    ha="center", va="center", fontsize=label_fontsize, fontweight="bold",
                    transform=ax_ch.transAxes,
                )

            for r, col_y_val in enumerate(col_y_vals):
                ax_rh = fig.add_subplot(outer[r + 1, 0])
                ax_rh.axis("off")
                ax_rh.text(
                    0.93, 0.5, _wrap_row_label(str(col_y_val)),
                    ha="right", va="center", ma="right", fontsize=label_fontsize, fontweight="bold",
                    transform=ax_rh.transAxes,
                )

                for c, col_x_val in enumerate(col_x_vals):
                    facet_label = f"{col_x_val} | {col_y_val}"
                    ax = fig.add_subplot(outer[r + 1, c + 1])
                    diff_df = None
                    if level_a is not None:
                        ax.set_title(
                            f"{level_a}  vs  {level_b}",
                            fontsize=max(label_fontsize - 2, 6), fontweight="normal", pad=4,
                        )
                        diff_df = diff_stats_by_facet.get(facet_label, {}).get((level_a, level_b))
                    if not _diff_df_has_data(diff_df):
                        ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                        ax.axis("off")
                        continue
                    draw_diff_barh(ax, class_order, diff_df, colors, xlim=global_xmax)
                    if c > 0:
                        ax.set_yticklabels([])

                    df_long = diff_df.copy()
                    df_long[col_x] = col_x_val
                    df_long[col_y] = col_y_val
                    if extra_group_cols:
                        for gc, gv in zip(extra_group_cols, str(page_key).split(" | ")):
                            df_long[gc] = gv
                    df_long["level_a"] = level_a
                    df_long["level_b"] = level_b
                    long_rows.append(df_long)

            page_title = title if level_a is not None else f"{title} — no {condition_col} comparison available"
            if page_key != "(all)":
                page_title += f" ({', '.join(extra_group_cols)}: {page_key})"
            fig.suptitle(page_title, fontsize=11, fontweight="bold", y=0.99)
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    csv_path = None
    if out_csv is not None:
        merged = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
        csv_path = str(out_csv)
        merged.to_csv(csv_path, index=False)

    return {"pdf_path": out_pdf, "csv_path": csv_path}
