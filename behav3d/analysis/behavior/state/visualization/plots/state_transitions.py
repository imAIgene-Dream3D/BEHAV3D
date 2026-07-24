import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import multiprocessing
import re
import threading
import warnings
from io import BytesIO
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from behav3d.analysis.behavior.utils import _natural_sort_key, _sanitize_filename_token
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)

# -----------------------------
# Plot timepoint>timepoint state transition matrix
# -----------------------------


def _is_main_thread():
    """Return True when running on the interpreter's main thread."""
    return threading.current_thread() is threading.main_thread()


def _plotly_static_export_worker(payload, conn):
    """Render a Plotly figure inside a spawned child process."""
    try:
        fig = go.Figure(payload["figure"])
        export_format = str(payload["export_format"])
        width = int(payload.get("width", 1400))
        height = int(payload.get("height", 700))
        scale = int(payload.get("scale", 2))
        output_path = payload.get("output_path", None)

        if export_format == "pdf":
            if output_path is None:
                raise ValueError("output_path is required for PDF static export.")
            fig.write_image(str(output_path), format="pdf", width=width, height=height, scale=scale)
            conn.send({"ok": True, "output_path": str(output_path)})
        elif export_format == "png":
            png_bytes = fig.to_image(format="png", width=width, height=height, scale=scale)
            conn.send({"ok": True, "png_bytes": png_bytes})
        else:
            raise ValueError(f"Unsupported Plotly static export format: {export_format}")
    except Exception as exc:
        conn.send({"ok": False, "error": str(exc)})
    finally:
        conn.close()


def _export_plotly_static_with_spawn(
    fig,
    *,
    export_format,
    output_path=None,
    width=1400,
    height=700,
    scale=2,
):
    """Run Plotly static export in a dedicated spawned child process."""
    ctx = multiprocessing.get_context("spawn")
    parent_conn, child_conn = ctx.Pipe(duplex=False)
    payload = {
        "figure": fig.to_dict(),
        "export_format": str(export_format),
        "output_path": (str(output_path) if output_path is not None else None),
        "width": int(width),
        "height": int(height),
        "scale": int(scale),
    }
    proc = ctx.Process(target=_plotly_static_export_worker, args=(payload, child_conn))
    try:
        proc.start()
        child_conn.close()
        try:
            result = parent_conn.recv()
        except EOFError:
            result = {"ok": False, "error": "spawned Plotly export worker exited without returning a result"}
        proc.join()
        if (proc.exitcode not in (0, None)) and bool(result.get("ok", False)):
            result = {
                "ok": False,
                "error": f"spawned Plotly export worker exited with code {proc.exitcode}",
            }
        return result
    finally:
        try:
            parent_conn.close()
        except Exception:
            pass
        if proc.is_alive():
            proc.terminate()
            proc.join()


def _export_plotly_static(
    fig,
    *,
    export_format,
    output_path=None,
    width=1400,
    height=700,
    scale=2,
    use_spawn=False,
):
    """Export Plotly figure either directly or through a spawned child process."""
    if use_spawn:
        return _export_plotly_static_with_spawn(
            fig,
            export_format=export_format,
            output_path=output_path,
            width=width,
            height=height,
            scale=scale,
        )

    if str(export_format) == "pdf":
        if output_path is None:
            raise ValueError("output_path is required for PDF static export.")
        fig.write_image(str(output_path), format="pdf", width=width, height=height, scale=scale)
        return {"ok": True, "output_path": str(output_path)}
    if str(export_format) == "png":
        png_bytes = fig.to_image(format="png", width=width, height=height, scale=scale)
        return {"ok": True, "png_bytes": png_bytes}
    raise ValueError(f"Unsupported Plotly static export format: {export_format}")


def _write_sankey_html(fig, *, html_dir, start_state, end_state):
    """Persist Sankey plot as HTML and return the written path."""
    html_dir = Path(html_dir)
    html_dir.mkdir(parents=True, exist_ok=True)
    html_name = (
        "sankey_"
        + _sanitize_filename_token(start_state, "start")
        + "_to_"
        + _sanitize_filename_token(end_state, "end")
        + ".html"
    )
    html_path = html_dir / html_name
    fig.write_html(str(html_path), include_plotlyjs="cdn", full_html=True)
    return html_path

def compute_cluster_transition_matrix(
    adata,
    cluster_key: str,
    id_cols=("sample_name", "TrackID"),
    time_key="position_t",
    normalize: bool = True,
    plot: bool = False,
    ax: plt.Axes = None,
    only_transitions: bool = False,
    state_order=None,
):
    """
    Compute a transition matrix between cluster states from tracked objects over time.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing tracking and clustering info in .obs.
    cluster_key : str
        Column in adata.obs with cluster labels (e.g. 'ClusterID', 'leiden').
    id_cols : sequence of str, default ("sample_name", "TrackID")
        Columns in adata.obs that together identify each track/object.
    time_key : str, default "position_t"
        Column in adata.obs giving time or frame index (must be sortable).
    normalize : bool, default True
        If True, return row-normalized probabilities (HMM-style).
        If False, returns raw transition counts.
    plot : bool, default False
        If True, plot the transition matrix as a heatmap.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None and plot=True, a new figure/axis is created.
    only_transitions : bool, default False
        If True, remove self-transitions (diagonal) from the *returned* matrices
        by setting diagonal counts to 0 and re-normalizing across off-diagonal
        entries (so rows sum to 1 when there is at least one off-diagonal transition).
        Also makes the diagonal appear empty in the plot.
    state_order : list of str, optional
        Saved display order for the states (e.g. from `_get_classification_state_order`). States
        not present in this list are appended afterwards. Defaults to alphabetical.

    Returns
    -------
    transition_counts : pandas.DataFrame
        Matrix of transition counts, shape (n_states, n_states).
        Rows = current state, columns = next state.
    transition_probs : pandas.DataFrame
        Row-normalized transition matrix.
        If only_transitions=True, this is P(next | current, next != current).
    """
    if state_order is None:
        state_order = _get_classification_state_order(adata, cluster_key)

    df = adata.obs[list(id_cols) + [cluster_key, time_key]].copy()
    df = df.sort_values(list(id_cols) + [time_key])

    # Next state within each track
    df["next_state"] = df.groupby(list(id_cols))[cluster_key].shift(-1)

    # Drop last timepoint of each track (no "next" state)
    transitions = df.dropna(subset=["next_state"]).copy()

    # Build transition counts
    transition_counts = pd.crosstab(
        transitions[cluster_key],
        transitions["next_state"]
    )

    # Ensure all states present on both axes, in the saved display order (falling back to
    # alphabetical for any state absent from `state_order`). Order matching is done via string
    # keys (matching `_get_classification_state_order`'s stored labels) but the original
    # `cluster_key` dtype (e.g. int cluster ids) is preserved for the reindex/tick labels below.
    states_raw = sorted(df[cluster_key].unique())
    ordered_labels = _apply_state_order([str(s) for s in states_raw], state_order)
    label_to_raw = {str(s): s for s in states_raw}
    states = [label_to_raw[label] for label in ordered_labels]
    transition_counts = transition_counts.reindex(
        index=states, columns=states, fill_value=0
    )

    # Optionally remove self-transitions in the returned matrices
    if only_transitions:
        counts_arr = transition_counts.to_numpy(copy=True)
        np.fill_diagonal(counts_arr, 0)
        transition_counts = pd.DataFrame(
            counts_arr,
            index=transition_counts.index,
            columns=transition_counts.columns,
        )

    # Row-normalize to probabilities
    row_sums = transition_counts.sum(axis=1)
    transition_probs = transition_counts.div(row_sums.replace(0, np.nan), axis=0)

    if plot:
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        data_to_plot = transition_probs if normalize else transition_counts
        plot_arr = data_to_plot.to_numpy(dtype=float, copy=True)

        # Make diagonal appear empty ONLY when only_transitions=True
        if only_transitions:
            diag_mask = np.eye(plot_arr.shape[0], dtype=bool)
            plot_arr = np.ma.array(plot_arr, mask=diag_mask)

        imshow_kwargs = {}
        if normalize:  # probability heatmap
            imshow_kwargs.update(dict(vmin=0.0, vmax=1.0))
            
        im = ax.imshow(plot_arr, aspect="auto", **imshow_kwargs)

        # Make masked values transparent (so diagonal looks empty)
        if only_transitions:
            cmap = im.get_cmap()
            cmap.set_bad(alpha=0)

        ax.set_xticks(np.arange(len(states)))
        ax.set_yticks(np.arange(len(states)))
        ax.set_xticklabels(states, rotation=90)
        ax.set_yticklabels(states)

        plt.colorbar(im, ax=ax, label="P(next | current)" if normalize else "Count")
        ax.set_xlabel("Next state")
        ax.set_ylabel("Current state")
        ax.set_title("Cluster transition matrix")
        plt.tight_layout()

    return transition_counts, transition_probs
 
##########################################################
# N-gram calculations to describe state behavioral transition
##########################################################

# Collapse state runs into bouts (111222333 -> 123)
def collapse_runs(seq):
    """Collapse consecutive identical states into bouts."""
    out = []
    prev = None
    first = True
    for x in seq:
        if first or x != prev:
            out.append(x)
            prev = x
            first = False
    return out

# Calculate all n-grams across all tracks (e.g. 1>2, 1>2>3, 2>3>2, 5>4>5>2 etc)
def all_ngrams(
    adata_full,
    n_values=(2, 3),
    state_col="ClusterID",
    group_cols=("sample_name", "TrackID"),
    time_col="position_t",
    collapse_bouts=True,
):
    """
    Compute pooled n-gram counts across all tracks for arbitrary n.

    Returns a DataFrame with columns:
      n, ngram (tuple), ngram_str, end_state, count
    """
    obs = adata_full.obs

    needed = [state_col, *group_cols, time_col]
    missing = [c for c in needed if c not in obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    df = obs.loc[:, [*group_cols, time_col, state_col]].copy()
    df = df.sort_values([*group_cols, time_col], kind="stable")
    df["_state"] = df[state_col].astype(str).str.strip()

    # counters[n] is a plain dict: {ngram_tuple: count}
    counters = {n: {} for n in n_values}

    for _, g in df.groupby(list(group_cols), sort=False, observed=False):
        seq = g["_state"].tolist()
        if collapse_bouts:
            seq = collapse_runs(seq)

        L = len(seq)
        for n in n_values:
            if L < n:
                continue
            cnt_map = counters[n]
            for t in range(n - 1, L):
                gram = tuple(seq[t - n + 1 : t + 1])
                cnt_map[gram] = cnt_map.get(gram, 0) + 1

    rows = []
    for n in n_values:
        for gram, cnt in counters[n].items():
            rows.append(
                {
                    "n": n,
                    "ngram": gram,
                    "ngram_str": " > ".join(gram),
                    "end_state": gram[-1],
                    "count": cnt,
                }
            )

    if len(rows) == 0:
        return pd.DataFrame(
            columns=["n", "ngram", "ngram_str", "end_state", "count"]
        )

    out = pd.DataFrame(rows)
    out = out.sort_values(["n", "count"], ascending=[True, False], kind="stable").reset_index(drop=True)
    return out

# Get top-N n-grams (of order n) that lead to end_state
def top_ngrams_per_end_state(df_ngrams, n, top_n=10):
    """Get top-N n-grams (of order n) per end-state."""
    sub = df_ngrams[df_ngrams["n"] == n]
    out = {}
    for end_state, g in sub.groupby("end_state", sort=True):
        out[end_state] = (
            g.sort_values("count", ascending=False, kind="stable")
            .head(top_n)
            .reset_index(drop=True)
        )
    return out

# Plot top-N n-grams (of order n) pooled across end-states
def plot_top_ngrams(df_ngrams, n, top_n=30, title=None):
    sub = df_ngrams[df_ngrams["n"] == n].head(top_n)
    if sub.empty:
        print(f"No {n}-grams to plot.")
        return

    plt.figure(figsize=(12, max(4, 0.35 * len(sub))))
    plt.barh(sub["ngram_str"][::-1], sub["count"][::-1])
    plt.xlabel("Count (pooled across tracks)")
    plt.ylabel(f"{n}-gram")
    plt.title(title or f"Top {top_n} {n}-grams")
    plt.tight_layout()
    plt.show()

def plot_top_ngrams_per_end_state(
    df_ngrams,
    n,
    top_n=10,
    min_total_end_occurrences=1,
):
    sub = df_ngrams[df_ngrams["n"] == n]
    totals = sub.groupby("end_state")["count"].sum().sort_values(ascending=False)

    for end_state in totals.index:
        if totals[end_state] < min_total_end_occurrences:
            continue

        g = (
            sub[sub["end_state"] == end_state]
            .sort_values("count", ascending=False, kind="stable")
            .head(top_n)
        )
        if g.empty:
            continue

        plt.figure(figsize=(12, max(4, 0.35 * len(g))))
        plt.barh(g["ngram_str"][::-1], g["count"][::-1])
        plt.xlabel("Count (pooled across tracks)")
        plt.ylabel(f"{n}-gram")
        plt.title(
            f"Top {top_n} {n}-grams ending in {end_state} "
            f"(total={int(totals[end_state])})"
        )
        plt.tight_layout()
        plt.show()

# -----------------------------
# Calculate all paths from a beginning state to the first encountered selected end state
# -----------------------------
def paths_between_states(
    adata_full,
    start_state,
    end_state,
    state_col="ClusterID",
    group_cols=("sample_name", "TrackID"),
    time_col="position_t",
    collapse_bouts=True,
    mode="next_end",
    allow_revisit_end=False,
    max_len=None,
):
    """
    Extract observed paths from start_state to end_state across all tracks.

    mode:
      - "next_end": for each start occurrence, take the path to the *next* end occurrence after it.
      - "any_end":  for each start occurrence, take paths to *all* later end occurrences.
    """
    obs = adata_full.obs
    needed = [state_col, *group_cols, time_col]
    missing = [c for c in needed if c not in obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    df = obs.loc[:, [*group_cols, time_col, state_col]].copy()
    df = df.sort_values([*group_cols, time_col], kind="stable")
    df["_state"] = df[state_col].astype(str).str.strip()

    s = str(start_state).strip()
    e = str(end_state).strip()

    # counter is a plain dict: {path_tuple: count}
    counter = {}

    for _, g in df.groupby(list(group_cols), sort=False, observed=False):
        seq = g["_state"].tolist()
        if collapse_bouts:
            seq = collapse_runs(seq)

        start_idx = [i for i, x in enumerate(seq) if x == s]
        end_idx = [i for i, x in enumerate(seq) if x == e]
        if not start_idx or not end_idx:
            continue

        for i in start_idx:
            ends_after = [j for j in end_idx if j > i]
            if not ends_after:
                continue

            if mode == "next_end":
                ends_after = [ends_after[0]]

            for j in ends_after:
                path = seq[i : j + 1]

                if not allow_revisit_end and len(path) > 1:
                    # if end appears earlier in path (besides last), truncate at first
                    for k in range(1, len(path) - 1):
                        if path[k] == e:
                            path = path[: k + 1]
                            break

                if max_len is not None and len(path) > max_len:
                    continue

                tpath = tuple(path)
                counter[tpath] = counter.get(tpath, 0) + 1

    out = pd.DataFrame(
        [{"path": p, "path_str": " > ".join(p), "length": len(p), "count": c}
         for p, c in counter.items()]
    )

    if out.empty:
        return out.assign(
            path_str=pd.Series(dtype=str),
            length=pd.Series(dtype=int),
            count=pd.Series(dtype=int),
        )

    out = out.sort_values(
        ["length", "count"], ascending=[True, False], kind="stable"
    ).reset_index(drop=True)
    return out

# Plot a ranking of paths from beginning to end state by occurrence count
def plot_paths_by_count(
    df_paths,
    top_n=30,
    min_count=1,
    title=None,
):
    if df_paths.empty:
        print("No paths to plot.")
        return

    sub = (
        df_paths[df_paths["count"] >= min_count]
        .sort_values("count", ascending=False, kind="stable")
        .head(top_n)
        .copy()
    )

    if sub.empty:
        print("No paths pass min_count filter.")
        return

    plt.figure(figsize=(12, max(4, 0.35 * len(sub))))
    plt.barh(sub["path_str"][::-1], sub["count"][::-1])
    plt.xlabel("Count (number of occurrences)")
    plt.ylabel("Path")
    plt.title(title or "Paths ranked by occurrence")
    plt.tight_layout()
    plt.show()


def _filter_paths_for_sankey(
    df,
    count_col="count",
    min_count=0,
    relative_count=0.0,
):
    """Filter paths by absolute and/or relative occurrence thresholds."""
    if count_col not in df.columns:
        raise KeyError(f"Missing '{count_col}' column in path table.")

    work = df.copy()
    work[count_col] = pd.to_numeric(work[count_col], errors="coerce").fillna(0.0)
    work = work[work[count_col] > 0].copy()

    min_count = float(min_count)
    relative_count = float(relative_count)
    if min_count < 0:
        raise ValueError("min_count must be >= 0.")
    if not (0.0 <= relative_count <= 1.0):
        raise ValueError("relative_count must be in [0, 1].")

    total_count = float(work[count_col].sum())
    if min_count > 0:
        work = work[work[count_col] >= min_count].copy()
    if relative_count > 0 and total_count > 0:
        work = work[(work[count_col] / total_count) >= relative_count].copy()

    return work, total_count


# Plot a Sankey diagram of paths between states, proportions of most common behavioral routes
def plot_sankey_diagram_between_states(
    df,
    path_col="path_str",
    count_col="count",
    sep=" > ",
    min_count=0,
    arrangement="snap",
    link_alpha=0.45,
    state_colors=None,
    show_node_labels=False,
    transparent_bg=True,
    font_size=25,
    relative_count=0.0,
):
    df, _total_count = _filter_paths_for_sankey(
        df,
        count_col=count_col,
        min_count=min_count,
        relative_count=relative_count,
    )
    if df.empty:
        raise ValueError("No rows left after filtering by min_count/relative_count")

    def to_rgba(color, a):
        if color is None:
            return None
        c = str(color).strip()
        m = re.match(r"rgba\((\d+),(\d+),(\d+),([0-9.]+)\)", c.replace(" ", ""), re.I)
        if m:
            r, g, b, _ = m.groups()
            return f"rgba({r},{g},{b},{a})"
        m = re.match(r"rgb\((\d+),(\d+),(\d+)\)", c.replace(" ", ""), re.I)
        if m:
            r, g, b = m.groups()
            return f"rgba({r},{g},{b},{a})"
        m = re.match(r"^#([0-9a-fA-F]{6})$", c)
        if m:
            hx = m.group(1)
            r, g, b = int(hx[0:2], 16), int(hx[2:4], 16), int(hx[4:6], 16)
            return f"rgba({r},{g},{b},{a})"
        return c

    def parse_states(s):
        return [x.strip() for x in str(s).split(sep) if str(x).strip()]

    raw = []
    max_len = 0
    seen = set()
    category_order = []
    for _, r in df.iterrows():
        states = parse_states(r[path_col])
        if len(states) < 2:
            continue
        w = float(r[count_col])
        raw.append((states, w))
        if len(states) > max_len:
            max_len = len(states)
        for st in states:
            if st not in seen:
                seen.add(st)
                category_order.append(st)

    if not raw:
        raise ValueError("No valid paths (need at least 2 states per row)")

    n_spots = max_len

    # Build color map, filling in any state missing from `state_colors` with the same
    # hash-stable default used by every other state report (instead of a local hash scheme).
    color_map = {}
    if isinstance(state_colors, dict):
        color_map = {str(k): str(v) for k, v in state_colors.items()}
    elif isinstance(state_colors, (list, tuple)):
        cols = list(state_colors)
        if cols:
            for i, st in enumerate(category_order):
                color_map[st] = cols[i % len(cols)]
    elif state_colors is not None:
        raise TypeError("state_colors must be a dict, list/tuple, or None")
    color_map = _normalize_label_color_map(category_order, colors=color_map)

    def base_color(state):
        return to_rgba(color_map[str(state)], 1.0)

    def link_color(key):
        if key is None:
            return f"rgba(150,150,150,{link_alpha})"
        return to_rgba(base_color(str(key)), link_alpha)

    def blank_name(c):
        return f"__BLANK_{c}__"

    def is_blank(s):
        return str(s).startswith("__BLANK_")

    def expand_to_cols(states):
        # right-align path in n_spots columns, keep first at col 0 and last at col n_spots-1
        L = len(states)
        cols = [0] + [None] * (L - 2) + [n_spots - 1]
        for i in range(1, L - 1):
            cols[i] = n_spots - (L - i)
        expanded = [None] * n_spots
        for st, c in zip(states, cols):
            expanded[c] = st
        for c in range(n_spots):
            if expanded[c] is None:
                expanded[c] = blank_name(c)
        return expanded

    node_id = {}  # (col, state) -> node index
    node_labels = []
    node_colors = []
    node_x = []

    def add_node(col, state):
        key = (col, state)
        if key in node_id:
            return node_id[key]
        nid = len(node_labels)
        node_id[key] = nid

        if is_blank(state):
            node_labels.append("")
            node_colors.append("rgba(0,0,0,0)")
        else:
            node_labels.append(str(state) if show_node_labels else "")
            node_colors.append(base_color(str(state)))

        node_x.append(col / (n_spots - 1) if n_spots > 1 else 0)
        return nid

    # link_value: (s, t, key_color_state) -> value
    link_value = {}

    # blank_weight: blank_node_id -> {key_state: weight}
    blank_weight = {}

    for states, w in raw:
        expanded = expand_to_cols(states)
        last_real = None
        for c in range(n_spots - 1):
            cur_state = expanded[c]
            nxt_state = expanded[c + 1]

            if not is_blank(cur_state):
                last_real = cur_state

            if last_real is None:
                # fallback: first non-blank in expanded
                for s0 in expanded:
                    if not is_blank(s0):
                        last_real = s0
                        break

            key_state = last_real

            s_id = add_node(c, cur_state)
            t_id = add_node(c + 1, nxt_state)

            lk = (s_id, t_id, key_state)
            link_value[lk] = link_value.get(lk, 0.0) + w

            if is_blank(cur_state):
                by_key = blank_weight.get(s_id)
                if by_key is None:
                    by_key = {}
                    blank_weight[s_id] = by_key
                by_key[key_state] = by_key.get(key_state, 0.0) + w

            if is_blank(nxt_state):
                by_key = blank_weight.get(t_id)
                if by_key is None:
                    by_key = {}
                    blank_weight[t_id] = by_key
                by_key[key_state] = by_key.get(key_state, 0.0) + w

    # Color blank nodes using the most common key_state passing through them
    for blank_nid, by_key in blank_weight.items():
        best_key = max(by_key.items(), key=lambda kv: kv[1])[0]
        node_colors[blank_nid] = link_color(best_key)

    sources, targets, values, colors = [], [], [], []
    for (s_id, t_id, key_state), v in link_value.items():
        sources.append(s_id)
        targets.append(t_id)
        values.append(v)
        colors.append(link_color(key_state))

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement=arrangement,
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="rgba(0,0,0,0)", width=0),
                    label=node_labels,
                    color=node_colors,
                    x=node_x,
                ),
                link=dict(source=sources, target=targets, value=values, color=colors),
            )
        ]
    )

    # legend entries
    for st in category_order:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(size=10, color=to_rgba(base_color(st), 1.0)),
                name=str(st),
                hoverinfo="skip",
                showlegend=True,
            )
        )

    bg = "rgba(0,0,0,0)" if transparent_bg else "white"
    legend_bg = "rgba(0,0,0,0)" if transparent_bg else "rgba(255,255,255,0.8)"

    threshold_text = []
    if float(min_count) > 0:
        threshold_text.append(f"count >= {min_count}")
    if float(relative_count) > 0:
        threshold_text.append(f"share >= {100.0 * float(relative_count):.1f}%")
    if len(threshold_text) == 0:
        threshold_suffix = "no path filter"
    else:
        threshold_suffix = ", ".join(threshold_text)

    fig.update_layout(
        title_text=f"Right-aligned Sankey ({threshold_suffix})",
        font_size=font_size,
        margin=dict(l=20, r=220, t=50, b=20),
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


def _plot_transition_heatmap_page(df_matrix, title, colorbar_label, normalized, mask_diagonal=False):
    fig, ax = plt.subplots(figsize=(8.27, 11.69), constrained_layout=True)
    arr = df_matrix.to_numpy(dtype=float, copy=True)
    if mask_diagonal:
        arr = np.ma.array(arr, mask=np.eye(arr.shape[0], dtype=bool))
    if normalized:
        im = ax.imshow(arr, aspect="equal", vmin=0.0, vmax=1.0)
    else:
        im = ax.imshow(arr, aspect="equal")
    if mask_diagonal:
        cmap = plt.colormaps[im.cmap.name].copy()
        cmap.set_bad(color="white")
        im.set_cmap(cmap)

    n_rows, n_cols = arr.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))
    ax.set_xticklabels([str(c) for c in df_matrix.columns], rotation=90)
    ax.set_yticklabels([str(i) for i in df_matrix.index])
    ax.set_xlabel("Next state")
    ax.set_ylabel("Current state")
    ax.set_title(str(title))
    plt.colorbar(im, ax=ax, label=str(colorbar_label))
    return fig


def _plot_heatmap_on_ax(ax, df_matrix, title, colorbar_label, normalized, mask_diagonal=False):
    """Plot a square transition heatmap onto an existing axes object."""
    arr = df_matrix.to_numpy(dtype=float, copy=True)
    if mask_diagonal:
        arr = np.ma.array(arr, mask=np.eye(arr.shape[0], dtype=bool))
    if normalized:
        im = ax.imshow(arr, aspect="equal", vmin=0.0, vmax=1.0)
    else:
        im = ax.imshow(arr, aspect="equal")
    if mask_diagonal:
        cmap = plt.colormaps[im.cmap.name].copy()
        cmap.set_bad(color="white")
        im.set_cmap(cmap)

    n_rows, n_cols = arr.shape
    ax.set_xticks(np.arange(n_cols))
    ax.set_yticks(np.arange(n_rows))
    ax.set_xticklabels([str(c) for c in df_matrix.columns], rotation=90, fontsize=8)
    ax.set_yticklabels([str(i) for i in df_matrix.index], fontsize=8)
    ax.set_xlabel("Next state", fontsize=9)
    ax.set_ylabel("Current state", fontsize=9)
    ax.set_title(str(title), fontsize=10, pad=6)
    ax.get_figure().colorbar(im, ax=ax, label=str(colorbar_label), fraction=0.046, pad=0.04)


def _plot_transition_heatmap_quad_page(
    df_probs, df_probs_no_self,
    df_counts, df_counts_no_self,
    state_col="",
):
    """A4 portrait figure with all four transition heatmaps in a 2×2 square grid.

    Row 0: transition probabilities (full | no-self with hollow diagonal)
    Row 1: transition counts       (full | no-self with hollow diagonal)
    """
    fig, axes = plt.subplots(2, 2, figsize=(8.27, 11.69), constrained_layout=True)
    _plot_heatmap_on_ax(
        axes[0, 0], df_probs,
        "Transition Probabilities\n(all transitions incl. self)",
        "P(next | current)", normalized=True, mask_diagonal=False,
    )
    _plot_heatmap_on_ax(
        axes[0, 1], df_probs_no_self,
        "Transition Probabilities\n(self-transitions excluded)",
        "P(next | current, next ≠ current)", normalized=True, mask_diagonal=True,
    )
    _plot_heatmap_on_ax(
        axes[1, 0], df_counts,
        "Transition Counts\n(all transitions incl. self)",
        "Count", normalized=False, mask_diagonal=False,
    )
    _plot_heatmap_on_ax(
        axes[1, 1], df_counts_no_self,
        "Transition Counts\n(self-transitions excluded)",
        "Count  (diagonal = no data)", normalized=False, mask_diagonal=True,
    )
    if state_col:
        fig.suptitle(f"Transition Matrices  ·  {state_col}", fontsize=12, fontweight="bold")
    return fig


def _plot_ngram_batch_page(items, y_col="ngram_str", x_col="count"):
    """A4 portrait figure with 1 or 2 n-gram bar charts stacked vertically.

    items: list of (df_ranking, title) tuples — 1 or 2 entries.
    Returns fig or None if all items are empty.
    """
    non_empty = [(df, t) for df, t in items if df is not None and len(df) > 0]
    if not non_empty:
        return None
    n_rows = len(non_empty)
    fig, axes = plt.subplots(n_rows, 1, figsize=(8.27, 11.69), squeeze=False)
    for i, (df, title) in enumerate(non_empty):
        ax = axes[i, 0]
        work = df.copy()
        work[x_col] = pd.to_numeric(work[x_col], errors="coerce").fillna(0.0)
        work = work.sort_values(x_col, ascending=False, kind="stable")
        labels = work[y_col].astype(str).tolist()[::-1]
        values = work[x_col].astype(float).tolist()[::-1]
        ax.barh(labels, values)
        ax.set_xlabel("Count (pooled across tracks)", fontsize=9)
        ax.set_ylabel("N-gram", fontsize=9)
        ax.set_title(str(title), fontsize=10)
        ax.tick_params(axis="y", labelsize=8)
    fig.subplots_adjust(left=0.28, right=0.95, top=0.93, bottom=0.05, hspace=0.50)
    return fig


def _plot_ngram_ranking_page(df_ranking, title, y_col="ngram_str", x_col="count"):
    work = df_ranking.copy()
    if len(work) == 0:
        return None
    work[x_col] = pd.to_numeric(work[x_col], errors="coerce").fillna(0.0)
    work = work.sort_values(x_col, ascending=False, kind="stable")
    if len(work) == 0:
        return None

    fig_h = max(4.0, min(18.0, 0.35 * float(len(work)) + 1.5))
    fig, ax = plt.subplots(figsize=(11.0, fig_h))
    labels = work[y_col].astype(str).tolist()[::-1]
    values = work[x_col].astype(float).tolist()[::-1]
    ax.barh(labels, values)
    ax.set_xlabel("Count (pooled across tracks)")
    ax.set_ylabel("N-gram")
    ax.set_title(str(title))
    fig.tight_layout()
    return fig


def _coerce_ngram_orders(ngram_orders):
    if isinstance(ngram_orders, (int, np.integer)):
        vals = [int(ngram_orders)]
    else:
        vals = [int(v) for v in list(ngram_orders)]
    vals = sorted({int(v) for v in vals if int(v) >= 2})
    if len(vals) == 0:
        raise ValueError("ngram_orders must contain at least one integer >= 2.")
    return tuple(vals)


def _build_ngram_tables(
    adata,
    state_col,
    id_cols,
    time_col,
    ngram_orders,
    collapse_bouts,
    ngram_min_count,
    ngram_pooled_top_n,
    ngram_per_state_top_n,
):
    df_ngrams = all_ngrams(
        adata_full=adata,
        n_values=tuple(ngram_orders),
        state_col=str(state_col),
        group_cols=tuple(id_cols),
        time_col=str(time_col),
        collapse_bouts=bool(collapse_bouts),
    ).copy()

    if len(df_ngrams) == 0:
        df_ngrams = pd.DataFrame(
            columns=["n", "ngram", "ngram_str", "end_state", "count", "start_state"]
        )
        return df_ngrams, df_ngrams.copy(), df_ngrams.copy(), df_ngrams.copy()

    df_ngrams["count"] = pd.to_numeric(df_ngrams["count"], errors="coerce").fillna(0).astype(int)
    df_ngrams = df_ngrams[df_ngrams["count"] > 0].copy()
    df_ngrams["start_state"] = df_ngrams["ngram"].apply(
        lambda g: str(g[0]) if isinstance(g, tuple) and len(g) > 0 else ""
    )
    df_ngrams["end_state"] = df_ngrams["end_state"].astype(str)
    df_ngrams["ngram_str"] = df_ngrams["ngram_str"].astype(str)

    work = df_ngrams[df_ngrams["count"] >= int(ngram_min_count)].copy()

    pooled_rows = []
    for n in tuple(ngram_orders):
        sub = work[work["n"] == int(n)].sort_values("count", ascending=False, kind="stable").head(int(ngram_pooled_top_n))
        if len(sub) > 0:
            sub = sub.copy()
            sub["rank_scope"] = "pooled"
            pooled_rows.append(sub)
    pooled_df = pd.concat(pooled_rows, ignore_index=True) if len(pooled_rows) > 0 else work.iloc[0:0].copy()

    per_end_rows = []
    for n in tuple(ngram_orders):
        sub_n = work[work["n"] == int(n)]
        for end_state, g in sub_n.groupby("end_state", sort=False, observed=False):
            g = g.sort_values("count", ascending=False, kind="stable").head(int(ngram_per_state_top_n))
            if len(g) == 0:
                continue
            g = g.copy()
            g["rank_scope"] = "per_end_state"
            g["group_state"] = str(end_state)
            per_end_rows.append(g)
    per_end_df = pd.concat(per_end_rows, ignore_index=True) if len(per_end_rows) > 0 else work.iloc[0:0].copy()

    per_start_rows = []
    for n in tuple(ngram_orders):
        sub_n = work[work["n"] == int(n)]
        for start_state, g in sub_n.groupby("start_state", sort=False, observed=False):
            g = g.sort_values("count", ascending=False, kind="stable").head(int(ngram_per_state_top_n))
            if len(g) == 0:
                continue
            g = g.copy()
            g["rank_scope"] = "per_start_state"
            g["group_state"] = str(start_state)
            per_start_rows.append(g)
    per_start_df = pd.concat(per_start_rows, ignore_index=True) if len(per_start_rows) > 0 else work.iloc[0:0].copy()

    return df_ngrams, pooled_df, per_end_df, per_start_df


def save_state_transition_report(
    adata,
    output_dir,
    *,
    state_col="full_behavioral_cluster",
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    include_self_pairs=True,
    rows_per_page=3,
    sankey_min_count=100,
    sankey_relative_count=0.0,
    sankey_vector_export=True,
    sankey_one_plot_per_page=True,
    sankey_merge_backend="pypdf_optional",
    sankey_mode="next_end",
    collapse_bouts=True,
    allow_revisit_end=False,
    max_path_len=None,
    include_no_self_matrices=True,
    include_ngram_rankings=True,
    ngram_orders=(2, 3),
    ngram_pooled_top_n=30,
    ngram_per_state_top_n=10,
    ngram_include_per_end_state=True,
    ngram_include_per_start_state=True,
    ngram_min_count=1,
    state_colors=None,
    state_order=None,
    verbose=True,
):
    """Save transition matrix + all-pairs Sankey report for one state column."""
    if state_order is None:
        state_order = _get_classification_state_order(adata, state_col)
    if state_colors is None:
        state_colors = _get_classification_state_colors(adata, state_col)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir = output_dir / "transition_matrix"
    matrix_data_dir = matrix_dir / "data"
    matrix_data_dir.mkdir(parents=True, exist_ok=True)
    html_dir = output_dir / "sankey_html"

    include_no_self_matrices = bool(include_no_self_matrices)
    include_ngram_rankings = bool(include_ngram_rankings)
    ngram_include_per_end_state = bool(ngram_include_per_end_state)
    ngram_include_per_start_state = bool(ngram_include_per_start_state)
    ngram_orders = _coerce_ngram_orders(ngram_orders)
    ngram_pooled_top_n = max(1, int(ngram_pooled_top_n))
    ngram_per_state_top_n = max(1, int(ngram_per_state_top_n))
    ngram_min_count = max(1, int(ngram_min_count))

    counts, probs = compute_cluster_transition_matrix(
        adata,
        cluster_key=str(state_col),
        id_cols=tuple(id_cols),
        time_key=str(time_col),
        normalize=True,
        plot=False,
        only_transitions=False,
        state_order=state_order,
    )

    matrix_counts_csv = matrix_data_dir / "transition_matrix_counts.csv"
    matrix_probs_csv = matrix_data_dir / "transition_matrix_probs.csv"
    matrix_counts_no_self_csv = matrix_data_dir / "transition_matrix_counts_no_self.csv"
    matrix_probs_no_self_csv = matrix_data_dir / "transition_matrix_probs_no_self.csv"
    matrix_heatmap_pdf = matrix_dir / "transition_matrix_heatmap.pdf"
    counts.to_csv(matrix_counts_csv)
    probs.to_csv(matrix_probs_csv)

    counts_no_self = None
    probs_no_self = None
    if include_no_self_matrices:
        counts_no_self, probs_no_self = compute_cluster_transition_matrix(
            adata,
            cluster_key=str(state_col),
            id_cols=tuple(id_cols),
            time_key=str(time_col),
            normalize=True,
            plot=False,
            only_transitions=True,
            state_order=state_order,
        )
        counts_no_self.to_csv(matrix_counts_no_self_csv)
        probs_no_self.to_csv(matrix_probs_no_self_csv)

    if str(state_col) not in adata.obs.columns:
        raise KeyError(f"Missing state column '{state_col}' in adata.obs.")
    states = (
        pd.Series(adata.obs[str(state_col)])
        .astype("string")
        .dropna()
        .astype(str)
        .str.strip()
    )
    states = _apply_state_order(sorted(states.unique().tolist(), key=_natural_sort_key), state_order)

    ngrams_all_csv = matrix_data_dir / "transition_ngrams_all.csv"
    ngrams_pooled_csv = matrix_data_dir / "transition_ngrams_pooled.csv"
    ngrams_per_end_csv = matrix_data_dir / "transition_ngrams_per_end.csv"
    ngrams_per_start_csv = matrix_data_dir / "transition_ngrams_per_start.csv"

    df_ngrams = None
    df_ngrams_pooled = None
    df_ngrams_per_end = None
    df_ngrams_per_start = None
    if include_ngram_rankings:
        df_ngrams, df_ngrams_pooled, df_ngrams_per_end, df_ngrams_per_start = _build_ngram_tables(
            adata=adata,
            state_col=state_col,
            id_cols=id_cols,
            time_col=time_col,
            ngram_orders=ngram_orders,
            collapse_bouts=collapse_bouts,
            ngram_min_count=ngram_min_count,
            ngram_pooled_top_n=ngram_pooled_top_n,
            ngram_per_state_top_n=ngram_per_state_top_n,
        )
        df_ngrams.to_csv(ngrams_all_csv, index=False)
        df_ngrams_pooled.to_csv(ngrams_pooled_csv, index=False)
        df_ngrams_per_end.to_csv(ngrams_per_end_csv, index=False)
        df_ngrams_per_start.to_csv(ngrams_per_start_csv, index=False)

    with PdfPages(matrix_heatmap_pdf) as pdf:
        # Page 1 — all four heatmaps in a 2×2 square grid (when no-self matrices available)
        if include_no_self_matrices and probs_no_self is not None and counts_no_self is not None:
            fig_quad = _plot_transition_heatmap_quad_page(
                probs, probs_no_self, counts, counts_no_self, state_col=str(state_col),
            )
            pdf.savefig(fig_quad, bbox_inches="tight")
            plt.close(fig_quad)
        else:
            # Fallback: one A4 page per metric when no-self matrices are not computed
            fig_probs = _plot_transition_heatmap_page(
                probs,
                title=f"Transition Probabilities — all transitions  ({state_col})",
                colorbar_label="P(next | current)",
                normalized=True,
            )
            pdf.savefig(fig_probs, bbox_inches="tight")
            plt.close(fig_probs)

            fig_counts = _plot_transition_heatmap_page(
                counts,
                title=f"Transition Counts — all transitions  ({state_col})",
                colorbar_label="Count",
                normalized=False,
            )
            pdf.savefig(fig_counts, bbox_inches="tight")
            plt.close(fig_counts)

        # Pages 3+ — N-gram rankings (A4, batched 2 per page)
        if include_ngram_rankings and df_ngrams is not None:
            for n in ngram_orders:
                pooled_n = df_ngrams_pooled[df_ngrams_pooled["n"] == int(n)].copy()
                fig = _plot_ngram_batch_page(
                    [(pooled_n, f"Top {ngram_pooled_top_n} {n}-grams (pooled across all tracks)")],
                )
                if fig is not None:
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

            if ngram_include_per_end_state:
                end_states = sorted(
                    df_ngrams_per_end.get("group_state", pd.Series([], dtype=str)).astype(str).unique().tolist(),
                    key=_natural_sort_key,
                )
                per_end_items = []
                for end_state in end_states:
                    for n in ngram_orders:
                        sub = df_ngrams_per_end[
                            (df_ngrams_per_end["n"] == int(n))
                            & (df_ngrams_per_end["group_state"].astype(str) == str(end_state))
                        ].copy()
                        per_end_items.append((sub, f"Top {ngram_per_state_top_n} {n}-grams ending in state {end_state}"))
                for batch_start in range(0, len(per_end_items), 2):
                    batch = per_end_items[batch_start: batch_start + 2]
                    fig = _plot_ngram_batch_page(batch)
                    if fig is not None:
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)

            if ngram_include_per_start_state:
                start_states = sorted(
                    df_ngrams_per_start.get("group_state", pd.Series([], dtype=str)).astype(str).unique().tolist(),
                    key=_natural_sort_key,
                )
                per_start_items = []
                for start_state in start_states:
                    for n in ngram_orders:
                        sub = df_ngrams_per_start[
                            (df_ngrams_per_start["n"] == int(n))
                            & (df_ngrams_per_start["group_state"].astype(str) == str(start_state))
                        ].copy()
                        per_start_items.append((sub, f"Top {ngram_per_state_top_n} {n}-grams starting from state {start_state}"))
                for batch_start in range(0, len(per_start_items), 2):
                    batch = per_start_items[batch_start: batch_start + 2]
                    fig = _plot_ngram_batch_page(batch)
                    if fig is not None:
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)

    rows_per_page = max(1, int(rows_per_page))
    include_self_pairs = bool(include_self_pairs)
    sankey_min_count = int(sankey_min_count)
    sankey_relative_count = float(sankey_relative_count)
    sankey_vector_export = bool(sankey_vector_export)
    sankey_one_plot_per_page = bool(sankey_one_plot_per_page)
    sankey_merge_backend = str(sankey_merge_backend).strip().lower()
    if sankey_merge_backend not in {"pypdf_optional", "pypdf_required", "none"}:
        raise ValueError("sankey_merge_backend must be one of: 'pypdf_optional', 'pypdf_required', 'none'.")
    if not (0.0 <= sankey_relative_count <= 1.0):
        raise ValueError("sankey_relative_count must be in [0, 1].")
    if sankey_vector_export and (not sankey_one_plot_per_page) and verbose:
        print("Sankey vector export uses one plot per page; ignoring sankey_one_plot_per_page=False.")

    pair_rows = []
    sankey_pdf_pages_dir = output_dir / "sankey_pdf_pages"
    sankey_pdf_page_rows = []
    sankey_png_rows = []
    sankey_pdf_path = output_dir / "sankey_all_pairs.pdf"
    kaleido_available = True
    kaleido_warned = False
    use_spawn_for_static_export = not _is_main_thread()
    if use_spawn_for_static_export and verbose:
        print(
            "Sankey static export delegated to a spawned child process because the current thread is not the main thread."
        )

    for start_state in states:
        for end_state in states:
            if (not include_self_pairs) and (str(start_state) == str(end_state)):
                continue

            row = {
                "start_state": str(start_state),
                "end_state": str(end_state),
                "n_paths": 0,
                "count_sum": 0,
                "status": "no_paths",
                "sankey_output_path": None,
                "sankey_pdf_page": None,
                "error": None,
            }

            df_paths = paths_between_states(
                adata,
                start_state=str(start_state),
                end_state=str(end_state),
                state_col=str(state_col),
                group_cols=tuple(id_cols),
                time_col=str(time_col),
                collapse_bouts=bool(collapse_bouts),
                mode=str(sankey_mode),
                allow_revisit_end=bool(allow_revisit_end),
                max_len=max_path_len,
            )

            if not df_paths.empty:
                row["n_paths"] = int(len(df_paths))
                row["count_sum"] = int(pd.to_numeric(df_paths["count"], errors="coerce").fillna(0).sum())

            if len(df_paths) == 0:
                pair_rows.append(row)
                continue

            df_paths_filtered, _total_count = _filter_paths_for_sankey(
                df_paths,
                count_col="count",
                min_count=sankey_min_count,
                relative_count=sankey_relative_count,
            )
            row["n_paths_filtered"] = int(len(df_paths_filtered))
            row["count_sum_filtered"] = int(
                pd.to_numeric(df_paths_filtered.get("count", pd.Series([], dtype=float)), errors="coerce")
                .fillna(0)
                .sum()
            )
            if len(df_paths_filtered) == 0:
                row["status"] = "no_paths_after_filter"
                pair_rows.append(row)
                continue

            try:
                fig = plot_sankey_diagram_between_states(
                    df_paths,
                    state_colors=state_colors,
                    min_count=sankey_min_count,
                    relative_count=sankey_relative_count,
                )
            except Exception as exc:
                row["status"] = "error"
                row["error"] = str(exc)
                pair_rows.append(row)
                continue

            if kaleido_available and sankey_vector_export:
                try:
                    sankey_pdf_pages_dir.mkdir(parents=True, exist_ok=True)
                    page_filename = (
                        "sankey_"
                        + _sanitize_filename_token(start_state, "start")
                        + "_to_"
                        + _sanitize_filename_token(end_state, "end")
                        + ".pdf"
                    )
                    page_path = sankey_pdf_pages_dir / page_filename
                    export_result = _export_plotly_static(
                        fig,
                        export_format="pdf",
                        output_path=page_path,
                        width=1400,
                        height=700,
                        scale=2,
                        use_spawn=use_spawn_for_static_export,
                    )
                    if not bool(export_result.get("ok", False)):
                        raise RuntimeError(str(export_result.get("error", "unknown Plotly export error")))
                    row["status"] = "pdf_unmerged"
                    row["sankey_output_path"] = str(page_path)
                    sankey_pdf_page_rows.append({"row": row, "pdf_path": page_path})
                except Exception as exc:
                    kaleido_available = False
                    if not kaleido_warned:
                        warnings.warn(
                            "Static Plotly export unavailable for Sankey vector PDF; falling back to HTML files. "
                            f"Details: {exc}",
                            RuntimeWarning,
                        )
                        kaleido_warned = True

            if kaleido_available and (not sankey_vector_export):
                try:
                    export_result = _export_plotly_static(
                        fig,
                        export_format="png",
                        width=1400,
                        height=700,
                        scale=2,
                        use_spawn=use_spawn_for_static_export,
                    )
                    if not bool(export_result.get("ok", False)):
                        raise RuntimeError(str(export_result.get("error", "unknown Plotly export error")))
                    png_bytes = export_result["png_bytes"]
                    row["status"] = "pdf_pending"
                    sankey_png_rows.append({"row": row, "png_bytes": png_bytes})
                except Exception as exc:
                    kaleido_available = False
                    if not kaleido_warned:
                        warnings.warn(
                            "Kaleido unavailable for Sankey PDF export; falling back to HTML files. "
                            f"Details: {exc}",
                            RuntimeWarning,
                        )
                        kaleido_warned = True

            if not kaleido_available:
                try:
                    html_path = _write_sankey_html(
                        fig,
                        html_dir=html_dir,
                        start_state=start_state,
                        end_state=end_state,
                    )
                    row["status"] = "html"
                    row["sankey_output_path"] = str(html_path)
                except Exception as exc:
                    row["status"] = "error"
                    row["error"] = str(exc)

            pair_rows.append(row)

    sankey_pdf_written = False
    if sankey_vector_export and len(sankey_pdf_page_rows) > 0:
        if sankey_merge_backend == "none":
            sankey_pdf_written = False
        else:
            merger_cls = None
            try:
                from pypdf import PdfMerger as _PdfMerger
                merger_cls = _PdfMerger
            except Exception as exc:
                if sankey_merge_backend == "pypdf_required":
                    raise RuntimeError(
                        "sankey_merge_backend='pypdf_required' but pypdf is not available."
                    ) from exc
                if not kaleido_warned:
                    warnings.warn(
                        "pypdf unavailable for merging Sankey vector pages; keeping per-pair PDF files only.",
                        RuntimeWarning,
                    )
                    kaleido_warned = True
                merger_cls = None

            if merger_cls is not None:
                merger = merger_cls()
                for i, item in enumerate(sankey_pdf_page_rows):
                    row = item["row"]
                    page_path = item["pdf_path"]
                    merger.append(str(page_path))
                    row["status"] = "pdf"
                    row["sankey_pdf_page"] = int(i) + 1
                    row["sankey_output_path"] = str(sankey_pdf_path)
                with open(sankey_pdf_path, "wb") as f:
                    merger.write(f)
                merger.close()
                sankey_pdf_written = True

    if (not sankey_vector_export) and len(sankey_png_rows) > 0:
        with PdfPages(sankey_pdf_path) as pdf:
            for i in range(0, len(sankey_png_rows), rows_per_page):
                chunk = sankey_png_rows[i : i + rows_per_page]
                n_chunk = len(chunk)
                fig, axes = plt.subplots(n_chunk, 1, figsize=(11.69, 3.6 * n_chunk))
                if n_chunk == 1:
                    axes = [axes]

                for j, item in enumerate(chunk):
                    row = item["row"]
                    page_number = int((i + j) // rows_per_page) + 1
                    row["status"] = "pdf"
                    row["sankey_pdf_page"] = page_number
                    row["sankey_output_path"] = str(sankey_pdf_path)
                    img = plt.imread(BytesIO(item["png_bytes"]), format="png")
                    ax = axes[j]
                    ax.imshow(img)
                    ax.axis("off")
                    ax.set_title(
                        f"{row['start_state']} -> {row['end_state']} "
                        f"(n_paths={row['n_paths']}, count_sum={row['count_sum']})",
                        fontsize=10,
                    )

                fig.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
        sankey_pdf_written = True

    pair_index_df = pd.DataFrame(pair_rows)
    pair_index_csv = output_dir / "sankey_pairs_index.csv"
    pair_index_df.to_csv(pair_index_csv, index=False)

    if verbose:
        print(f"Saved transition matrix counts CSV: {matrix_counts_csv}")
        print(f"Saved transition matrix probabilities CSV: {matrix_probs_csv}")
        if include_no_self_matrices:
            print(f"Saved no-self transition matrix counts CSV: {matrix_counts_no_self_csv}")
            print(f"Saved no-self transition matrix probabilities CSV: {matrix_probs_no_self_csv}")
        print(f"Saved transition matrix heatmap PDF: {matrix_heatmap_pdf}")
        if include_ngram_rankings:
            print(f"Saved transition n-grams CSV (all): {ngrams_all_csv}")
            print(f"Saved transition n-grams CSV (pooled): {ngrams_pooled_csv}")
            print(f"Saved transition n-grams CSV (per-end): {ngrams_per_end_csv}")
            print(f"Saved transition n-grams CSV (per-start): {ngrams_per_start_csv}")
        if sankey_pdf_written:
            print(f"Saved Sankey all-pairs PDF: {sankey_pdf_path}")
        elif html_dir.exists():
            print(f"Saved Sankey HTML fallback directory: {html_dir}")
        print(f"Saved Sankey pair index CSV: {pair_index_csv}")

    return {
        "output_dir": str(output_dir),
        "state_col": str(state_col),
        "states": list(states),
        "sankey_min_count": int(sankey_min_count),
        "sankey_relative_count": float(sankey_relative_count),
        "transition_matrix_counts_csv": str(matrix_counts_csv),
        "transition_matrix_probs_csv": str(matrix_probs_csv),
        "transition_matrix_counts_no_self_csv": (
            str(matrix_counts_no_self_csv) if include_no_self_matrices else None
        ),
        "transition_matrix_probs_no_self_csv": (
            str(matrix_probs_no_self_csv) if include_no_self_matrices else None
        ),
        "transition_matrix_heatmap_pdf": str(matrix_heatmap_pdf),
        "transition_ngrams_all_csv": (str(ngrams_all_csv) if include_ngram_rankings else None),
        "transition_ngrams_pooled_csv": (str(ngrams_pooled_csv) if include_ngram_rankings else None),
        "transition_ngrams_per_end_csv": (str(ngrams_per_end_csv) if include_ngram_rankings else None),
        "transition_ngrams_per_start_csv": (str(ngrams_per_start_csv) if include_ngram_rankings else None),
        "sankey_pairs_index_csv": str(pair_index_csv),
        "sankey_pdf_path": (str(sankey_pdf_path) if sankey_pdf_written else None),
        "sankey_pdf_pages_dir": (str(sankey_pdf_pages_dir) if sankey_pdf_pages_dir.exists() else None),
        "sankey_html_dir": (str(html_dir) if html_dir.exists() else None),
    }
