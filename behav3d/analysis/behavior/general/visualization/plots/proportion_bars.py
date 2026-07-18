import hashlib

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from scipy import stats

A4_PORTRAIT = (8.27, 11.69)


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
        draw_thin_stacked_proportion_barh(
            axes[i], data_by_key.get(key) or {}, class_order, colors, xmax=1.0, bar_height_frac=bar_height_frac,
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


def compute_condition_diff_stats(
    per_sample_class_props,
    sample_metadata,
    *,
    class_order,
    condition_col,
    level_a,
    level_b,
    facet_col=None,
):
    """Per facet value (or a single '(all)' bucket if facet_col is None), per class: Welch's
    two-sided unpaired t-test comparing per-sample proportions between level_a and level_b.

    diff = mean_b - mean_a (level_b minus level_a).

    Parameters
    ----------
    per_sample_class_props : pd.DataFrame
        index=sample, columns=class_order, values=proportion (0-1).
    sample_metadata : pd.DataFrame
        index=sample, columns include condition_col (+ facet_col if given).

    Returns
    -------
    dict[str, pd.DataFrame]
        Keyed by facet value (or "(all)"); each DataFrame has columns:
        class, mean_a, mean_b, diff, t_stat, p_value, stars, n_a, n_b.
    """
    class_order = [str(c) for c in list(class_order)]
    joined = per_sample_class_props.join(sample_metadata, how="inner")

    if facet_col is not None and facet_col in joined.columns:
        facet_values = sorted(joined[facet_col].dropna().astype(str).unique().tolist())
        facet_groups = {fv: joined[joined[facet_col].astype(str) == fv] for fv in facet_values}
    else:
        facet_groups = {"(all)": joined}

    out = {}
    for facet_value, panel in facet_groups.items():
        rows = []
        cond = panel[condition_col].astype(str)
        vals_a_panel = panel[cond == str(level_a)]
        vals_b_panel = panel[cond == str(level_b)]
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
        out[str(facet_value)] = pd.DataFrame(rows)
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


def plot_condition_diff_composite(
    diff_stats_by_facet,
    *,
    class_order,
    colors,
    level_a_label,
    level_b_label,
    title,
    out_pdf,
    out_csv=None,
    figsize_per_panel=(3.2, None),
):
    """One composite figure: one side-by-side column panel per facet value (or a single panel if
    diff_stats_by_facet has only the '(all)' key), one draw_diff_barh row per class. Saves a PDF
    and (optionally) a merged long-form CSV."""
    class_order = [str(c) for c in list(class_order)]
    facet_values = list(diff_stats_by_facet.keys())
    n_panels = max(1, len(facet_values))
    panel_w = float(figsize_per_panel[0])
    panel_h = figsize_per_panel[1] if figsize_per_panel[1] is not None else max(1.2 + 0.45 * len(class_order), 3.0)
    fig, axes = plt.subplots(
        nrows=1, ncols=n_panels, figsize=(panel_w * n_panels, panel_h), squeeze=False,
        sharex=False,
    )
    axes = axes.flatten()

    global_xmax = 0.0
    for df in diff_stats_by_facet.values():
        if len(df) > 0:
            global_xmax = max(global_xmax, float(np.nanmax(np.abs(df["diff"].to_numpy(dtype=float))) * 100.0))
    global_xmax = max(global_xmax * 1.35, 5.0)

    for i, facet_value in enumerate(facet_values):
        ax = axes[i]
        diff_df = diff_stats_by_facet[facet_value]
        draw_diff_barh(ax, class_order, diff_df, colors, xlim=global_xmax)
        if i > 0:
            ax.set_yticklabels([])
        panel_title = facet_value if facet_value != "(all)" else ""
        header = f"{panel_title}\n" if panel_title else ""
        ax.set_title(f"{header}{level_a_label}          {level_b_label}", fontsize=8, pad=6)

    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.99)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))

    out_pdf = str(out_pdf)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    csv_path = None
    if out_csv is not None:
        long_rows = []
        for facet_value, df in diff_stats_by_facet.items():
            df = df.copy()
            df["facet_value"] = facet_value
            long_rows.append(df)
        merged = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
        csv_path = str(out_csv)
        merged.to_csv(csv_path, index=False)

    return {"pdf_path": out_pdf, "csv_path": csv_path}
