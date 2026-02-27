import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import hashlib
import re

# -----------------------------
# Plot timepoint>timepoint state transition matrix
# -----------------------------

def compute_cluster_transition_matrix(
    adata,
    cluster_key: str,
    id_cols=("sample_name", "TrackID"),
    time_key="position_t",
    normalize: bool = True,
    plot: bool = False,
    ax: plt.Axes = None,
    only_transitions: bool = False,
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

    Returns
    -------
    transition_counts : pandas.DataFrame
        Matrix of transition counts, shape (n_states, n_states).
        Rows = current state, columns = next state.
    transition_probs : pandas.DataFrame
        Row-normalized transition matrix.
        If only_transitions=True, this is P(next | current, next != current).
    """
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

    # Ensure all states present on both axes
    states = sorted(df[cluster_key].unique())
    transition_counts = transition_counts.reindex(
        index=states, columns=states, fill_value=0
    )

    # Optionally remove self-transitions in the returned matrices
    if only_transitions:
        np.fill_diagonal(transition_counts.values, 0)

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

    for _, g in df.groupby(list(group_cols), sort=False):
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

    out = (
        pd.DataFrame(rows)
        .sort_values(["n", "count"], ascending=[True, False], kind="stable")
        .reset_index(drop=True)
    )
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

    for _, g in df.groupby(list(group_cols), sort=False):
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
):
    df = df[df[count_col] > min_count].copy()
    if df.empty:
        raise ValueError("No rows left after filtering by min_count")

    def stable_rgba(label, a=1.0):
        h = hashlib.md5(str(label).encode("utf-8")).hexdigest()
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        r = int(0.55 * r + 0.45 * 255)
        g = int(0.55 * g + 0.45 * 255)
        b = int(0.55 * b + 0.45 * 255)
        return f"rgba({r},{g},{b},{a})"

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

    # Build color map (optional)
    color_map = {}
    if isinstance(state_colors, dict):
        color_map = {str(k): str(v) for k, v in state_colors.items()}
    elif isinstance(state_colors, (list, tuple)):
        cols = list(state_colors)
        if cols:
            for i, st in enumerate(category_order):
                color_map[st] = cols[i % len(cols)]
    elif state_colors is None:
        color_map = {}
    else:
        raise TypeError("state_colors must be a dict, list/tuple, or None")

    def base_color(state):
        st = str(state)
        if st in color_map:
            return to_rgba(color_map[st], 1.0)
        return stable_rgba(st, 1.0)

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

    fig.update_layout(
        title_text=f"Right-aligned Sankey (count > {min_count})",
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


test():
    # -----------------------------
    # Usage (same behavior as your script)
    # -----------------------------
    # compute 2-, 3-, and 4-grams
    cluster_column = "full_behavioral_cluster"  # matches your snippet; could also be "ClusterID"
    df_ngrams = all_ngrams(
        adata_full,
        n_values=(2, 3, 4),
        state_col=cluster_column,
        group_cols=("sample_name", "TrackID"),
        time_col="position_t",
        collapse_bouts=True,  # strongly recommended
    )

    plot_top_ngrams(df_ngrams, n=3, top_n=30)
    plot_top_ngrams_per_end_state(df_ngrams, n=4, top_n=10)

    df_paths_1_5 = paths_between_states(
        adata_full,
        start_state="static",
        end_state="organoid_contact",
        state_col=cluster_column,
        collapse_bouts=True,
        mode="next_end",
    )

    plot_paths_by_count(
        df_paths_1_5,
        top_n=25,
        min_count=1,
        title="Most common paths from state 1 to state 5",
    )

    colors = [
        "#d62728",
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    fig = plot_sankey_diagram_between_states(
        df_paths_1_5,
        state_colors=colors,
        min_count=100,
    )

    # export + show
    # (requires kaleido installed for write_image / to_image)
    fig.write_image(
        "/Users/s.deblank-3/Downloads/newplot.pdf",
        width=1400,
        height=700,
        scale=2,
    )
    fig.show()
