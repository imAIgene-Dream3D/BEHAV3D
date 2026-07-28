"""Sankey report of window-to-window trajectory-cluster transitions.

When a track is longer than the configured trajectory size, trajectory
clustering can split it into fixed-length, non-overlapping windows while
keeping the original ``TrackID`` (see ``trajectory_window_id`` /
``split_long_tracks`` in :mod:`behav3d.analysis.behavior.track.utils`), and
DTW-clusters each window independently. This module reconnects those split
windows back through their shared parent (``sample_name``, ``TrackID``) and
shows how the assigned cluster changes from one window to the next - unlike
:mod:`behav3d.analysis.behavior.state.visualization.plots.state_transitions`,
which pools per-timepoint state transitions across the whole population, this
keeps each sample's per-track window sequence intact.
"""
import warnings

import pandas as pd
import plotly.graph_objects as go

from behav3d.analysis.behavior.utils import _natural_sort_key, _sanitize_filename_token
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
    _export_plotly_static,
    _is_main_thread,
    _write_sankey_html,
)
from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths, _winfo


def compute_window_transition_links(
    adata_tracks,
    cluster_key,
    *,
    window_col="trajectory_window_id",
    id_cols=("sample_name", "TrackID"),
    sample_col="sample_name",
):
    """Per-track window-to-next-window cluster transition counts.

    ``adata_tracks.obs`` has one row per (sample, TrackID, window) - so a
    ``groupby(id_cols)`` + ``shift(-1)`` links each track's window directly to
    its own next window, the same way `compute_cluster_transition_matrix`
    links consecutive timepoints, but keyed on window order instead of raw
    time. Only genuinely consecutive windows (i, i+1) are linked - a track
    missing a window in the middle does not get its neighbours silently
    connected.

    Returns
    -------
    dict with keys "per_sample" and "pooled", each a DataFrame with columns
    [("sample_name" for per_sample only), window_from, window_to,
    cluster_from, cluster_to, count].
    """
    id_cols = list(id_cols)
    obs_cols = adata_tracks.obs.columns
    if str(window_col) not in obs_cols:
        raise KeyError(
            f"'{window_col}' not found in adata_tracks.obs. This report needs tracks split "
            "into windows - re-run trajectory clustering with 'Divide long tracks' enabled."
        )
    if str(cluster_key) not in obs_cols:
        raise KeyError(f"Missing cluster column '{cluster_key}' in adata_tracks.obs.")

    df = adata_tracks.obs.loc[:, id_cols + [str(window_col), str(cluster_key)]].copy()
    df[str(window_col)] = pd.to_numeric(df[str(window_col)], errors="coerce")
    df = df.dropna(subset=[str(window_col)])
    df[str(window_col)] = df[str(window_col)].astype(int)
    df[str(cluster_key)] = df[str(cluster_key)].astype(str).str.strip()
    df = df.sort_values(id_cols + [str(window_col)])

    df["_next_cluster"] = df.groupby(id_cols, observed=True)[str(cluster_key)].shift(-1)
    df["_next_window"] = df.groupby(id_cols, observed=True)[str(window_col)].shift(-1)
    transitions = df.dropna(subset=["_next_cluster", "_next_window"]).copy()
    transitions["_next_window"] = transitions["_next_window"].astype(int)
    transitions = transitions[transitions["_next_window"] == transitions[str(window_col)] + 1]

    pooled = (
        transitions.groupby(
            [str(window_col), str(cluster_key), "_next_cluster"], observed=True
        )
        .size()
        .reset_index(name="count")
        .rename(columns={str(window_col): "window_from", str(cluster_key): "cluster_from", "_next_cluster": "cluster_to"})
    )
    pooled["window_to"] = pooled["window_from"] + 1
    pooled = pooled[["window_from", "window_to", "cluster_from", "cluster_to", "count"]]

    per_sample = (
        transitions.groupby(
            [str(sample_col), str(window_col), str(cluster_key), "_next_cluster"], observed=True
        )
        .size()
        .reset_index(name="count")
        .rename(columns={
            str(sample_col): "sample_name", str(window_col): "window_from",
            str(cluster_key): "cluster_from", "_next_cluster": "cluster_to",
        })
    )
    per_sample["window_to"] = per_sample["window_from"] + 1
    per_sample = per_sample[["sample_name", "window_from", "window_to", "cluster_from", "cluster_to", "count"]]

    return {"per_sample": per_sample, "pooled": pooled}


def plot_window_transition_sankey(
    df_links,
    *,
    state_colors=None,
    state_order=None,
    link_alpha=0.45,
    show_node_labels=True,
    transparent_bg=True,
    font_size=14,
    title=None,
):
    """Sankey of window_i -> window_i+1 cluster transitions.

    ``df_links`` columns: window_from, window_to, cluster_from, cluster_to, count -
    e.g. one scope (one sample, or pooled) from `compute_window_transition_links`.
    Node columns are the real window indices in natural left-to-right order (no
    path-length padding/right-alignment, unlike `plot_sankey_diagram_between_states` -
    every track contributes exactly one real node per window it has).
    """
    if df_links is None or len(df_links) == 0:
        raise ValueError("No window transitions to plot (empty link table).")

    df = df_links.copy()
    df["cluster_from"] = df["cluster_from"].astype(str)
    df["cluster_to"] = df["cluster_to"].astype(str)
    df["window_from"] = df["window_from"].astype(int)
    df["window_to"] = df["window_to"].astype(int)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0.0)
    df = df[df["count"] > 0]
    if len(df) == 0:
        raise ValueError("No window transitions to plot after filtering.")

    windows = sorted(set(df["window_from"]).union(df["window_to"]))
    n_cols = len(windows)
    col_index = {w: i for i, w in enumerate(windows)}

    category_order = _apply_state_order(
        sorted(set(df["cluster_from"]) | set(df["cluster_to"]), key=_natural_sort_key),
        state_order,
    )
    color_map = _normalize_label_color_map(category_order, colors=state_colors)

    def to_rgba(hex_color, a):
        hx = str(hex_color).lstrip("#")
        r, g, b = int(hx[0:2], 16), int(hx[2:4], 16), int(hx[4:6], 16)
        return f"rgba({r},{g},{b},{a})"

    node_id = {}
    node_colors = []
    node_x = []

    def add_node(window, cluster):
        key = (window, cluster)
        if key in node_id:
            return node_id[key]
        nid = len(node_colors)
        node_id[key] = nid
        node_colors.append(to_rgba(color_map[cluster], 1.0))
        node_x.append(col_index[window] / (n_cols - 1) if n_cols > 1 else 0.0)
        return nid

    link_rows = []
    for _, r in df.iterrows():
        s_id = add_node(r["window_from"], r["cluster_from"])
        t_id = add_node(r["window_to"], r["cluster_to"])
        link_rows.append((s_id, t_id, r["cluster_from"], float(r["count"])))

    # Node label shows the flow actually explained by the drawn links: incoming
    # total for any node with predecessors, else its outgoing total (window 0).
    # Incoming/outgoing sums for the same interior node can legitimately differ
    # (some tracks stop after reaching it - not every track has the same number
    # of windows), so this is a documented choice, not an approximation.
    incoming_total, outgoing_total = {}, {}
    for s_id, t_id, _cf, v in link_rows:
        outgoing_total[s_id] = outgoing_total.get(s_id, 0.0) + v
        incoming_total[t_id] = incoming_total.get(t_id, 0.0) + v

    node_labels = [""] * len(node_colors)
    if show_node_labels:
        for (window, cluster), nid in node_id.items():
            total = incoming_total.get(nid, outgoing_total.get(nid, 0.0))
            node_labels[nid] = f"{cluster}  (n={int(round(total))})"

    sources = [s for s, _, _, _ in link_rows]
    targets = [t for _, t, _, _ in link_rows]
    values = [v for _, _, _, v in link_rows]
    colors = [to_rgba(color_map[cf], link_alpha) for _, _, cf, _ in link_rows]

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=15,
                    thickness=18,
                    line=dict(color="rgba(0,0,0,0)", width=0),
                    label=node_labels,
                    color=node_colors,
                    x=node_x,
                ),
                link=dict(source=sources, target=targets, value=values, color=colors),
            )
        ]
    )

    for st in category_order:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=to_rgba(color_map[st], 1.0)),
                name=st,
                hoverinfo="skip",
                showlegend=True,
            )
        )

    bg = "rgba(0,0,0,0)" if transparent_bg else "white"
    legend_bg = "rgba(0,0,0,0)" if transparent_bg else "rgba(255,255,255,0.8)"
    fig.update_layout(
        title_text=title or "Window-to-window trajectory transitions",
        font_size=font_size,
        margin=dict(l=20, r=200, t=50, b=20),
        legend=dict(
            x=1.02,
            y=0.5,
            xanchor="left",
            yanchor="middle",
            bgcolor=legend_bg,
            bordercolor="rgba(0,0,0,0)" if transparent_bg else "rgba(0,0,0,0.15)",
            borderwidth=0 if transparent_bg else 1,
        ),
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
        paper_bgcolor=bg,
        plot_bgcolor=bg,
    )
    return fig


def save_window_transition_report(
    adata_tracks,
    output_dir,
    cell_type,
    *,
    cluster_key=None,
    window_col="trajectory_window_id",
    id_cols=("sample_name", "TrackID"),
    sample_col="sample_name",
    include_pooled=True,
    min_count=1,
    state_colors=None,
    state_order=None,
    verbose=True,
):
    """Save one Sankey PDF per sample (+ optional pooled page) of window-to-window
    trajectory-cluster transitions, plus the underlying link-count CSVs.

    Unlike `save_state_transition_report` (per-timepoint HMM state, pooled across
    every track), this respects each track's own window sequence within its own
    sample - answering "did *this* track end up in a different trajectory cluster
    by its next window?" rather than a population-wide flicker pattern.
    """
    meta = adata_tracks.uns.get("dtai_trajectory_clustering", {})
    meta = meta if isinstance(meta, dict) else {}
    if cluster_key is None:
        cluster_key = str(meta.get("cluster_key", "ClusterID"))
    if state_order is None:
        state_order = _get_classification_state_order(adata_tracks, cluster_key)
    if state_colors is None:
        state_colors = _get_classification_state_colors(adata_tracks, cluster_key)

    paths = _resolve_dtaidistance_paths(output_dir, cell_type)
    out_dir = paths["outfolder"] / "window_transitions"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_pages_dir = out_dir / "sankey_pdf_pages"
    html_dir = out_dir / "sankey_html"

    links = compute_window_transition_links(
        adata_tracks, cluster_key, window_col=window_col, id_cols=id_cols, sample_col=sample_col,
    )
    per_sample_csv = out_dir / "window_transition_links_by_sample.csv"
    pooled_csv = out_dir / "window_transition_links_pooled.csv"
    links["per_sample"].to_csv(per_sample_csv, index=False)
    links["pooled"].to_csv(pooled_csv, index=False)

    samples = sorted(
        links["per_sample"]["sample_name"].astype(str).unique().tolist(), key=_natural_sort_key
    )

    use_spawn = not _is_main_thread()
    if use_spawn and verbose:
        _winfo(
            "trajectory-dtai",
            "Sankey static export delegated to a spawned child process (background thread).",
        )

    index_rows = []
    pdf_page_paths = []
    kaleido_available = True
    kaleido_warned = False

    def _make_and_export(df_links_sub, label, page_stem):
        nonlocal kaleido_available, kaleido_warned
        n_transitions = int(pd.to_numeric(df_links_sub["count"], errors="coerce").fillna(0).sum())
        row = {
            "scope": label, "n_transitions": n_transitions,
            "status": "no_transitions", "output_path": None, "error": None,
        }
        work = df_links_sub[df_links_sub["count"] >= int(min_count)] if min_count > 0 else df_links_sub
        if len(work) == 0:
            index_rows.append(row)
            return
        try:
            fig = plot_window_transition_sankey(
                work, state_colors=state_colors, state_order=state_order,
                title=f"{cell_type} — {label}: window-to-window trajectory transitions",
            )
        except Exception as exc:
            row["status"] = "error"
            row["error"] = str(exc)
            index_rows.append(row)
            return

        if kaleido_available:
            try:
                pdf_pages_dir.mkdir(parents=True, exist_ok=True)
                page_path = pdf_pages_dir / f"sankey_{_sanitize_filename_token(page_stem)}.pdf"
                result = _export_plotly_static(
                    fig, export_format="pdf", output_path=page_path,
                    width=1200, height=700, scale=2, use_spawn=use_spawn,
                )
                if not bool(result.get("ok", False)):
                    raise RuntimeError(str(result.get("error", "unknown Plotly export error")))
                row["status"] = "pdf"
                row["output_path"] = str(page_path)
                pdf_page_paths.append(page_path)
            except Exception as exc:
                kaleido_available = False
                if not kaleido_warned:
                    warnings.warn(
                        "Static Plotly export unavailable for Sankey PDF; falling back to HTML files. "
                        f"Details: {exc}",
                        RuntimeWarning,
                    )
                    kaleido_warned = True

        if not kaleido_available:
            try:
                html_path = _write_sankey_html(fig, html_dir=html_dir, start_state=page_stem, end_state="windows")
                row["status"] = "html"
                row["output_path"] = str(html_path)
            except Exception as exc:
                row["status"] = "error"
                row["error"] = str(exc)

        index_rows.append(row)

    for sample in samples:
        sub = links["per_sample"][links["per_sample"]["sample_name"].astype(str) == sample]
        _make_and_export(sub, sample, sample)

    if include_pooled:
        _make_and_export(links["pooled"], "All samples (pooled)", "pooled")

    merged_pdf_path = out_dir / "window_transitions_all_samples.pdf"
    merged = False
    if len(pdf_page_paths) > 0:
        try:
            # pypdf >=5 removed PdfMerger in favour of PdfWriter.append(); try the
            # current API first and fall back to PdfMerger for older installs.
            try:
                from pypdf import PdfWriter
                writer = PdfWriter()
                for page_path in pdf_page_paths:
                    writer.append(str(page_path))
                with open(merged_pdf_path, "wb") as f:
                    writer.write(f)
                writer.close()
            except ImportError:
                from pypdf import PdfMerger
                merger = PdfMerger()
                for page_path in pdf_page_paths:
                    merger.append(str(page_path))
                with open(merged_pdf_path, "wb") as f:
                    merger.write(f)
                merger.close()
            merged = True
            for row in index_rows:
                if row["status"] == "pdf":
                    row["merged_pdf"] = str(merged_pdf_path)
        except Exception as exc:
            if verbose:
                _winfo(
                    "trajectory-dtai",
                    f"could not merge Sankey pages into one PDF ({exc}); keeping per-page PDFs.",
                )

    index_df = pd.DataFrame(index_rows)
    index_csv = out_dir / "window_transitions_index.csv"
    index_df.to_csv(index_csv, index=False)

    if verbose:
        _winfo("trajectory-dtai", f"saved window-transition links CSV (per sample): {per_sample_csv}")
        _winfo("trajectory-dtai", f"saved window-transition links CSV (pooled): {pooled_csv}")
        if merged:
            _winfo("trajectory-dtai", f"saved window-transitions Sankey PDF: {merged_pdf_path}")
        elif pdf_pages_dir.exists():
            _winfo("trajectory-dtai", f"saved window-transitions Sankey PDF pages: {pdf_pages_dir}")
        elif html_dir.exists():
            _winfo("trajectory-dtai", f"saved window-transitions Sankey HTML fallback: {html_dir}")
        _winfo("trajectory-dtai", f"saved window-transitions index CSV: {index_csv}")

    return {
        "output_dir": str(out_dir),
        "cluster_key": str(cluster_key),
        "samples": samples,
        "per_sample_links_csv": str(per_sample_csv),
        "pooled_links_csv": str(pooled_csv),
        "index_csv": str(index_csv),
        "sankey_pdf": str(merged_pdf_path) if merged else None,
        "sankey_pdf_pages_dir": str(pdf_pages_dir) if pdf_pages_dir.exists() else None,
        "sankey_html_dir": str(html_dir) if html_dir.exists() else None,
    }
