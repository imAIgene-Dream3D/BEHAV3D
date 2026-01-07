import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def _collapse_runs(seq):
    out, prev = [], object()
    for x in seq:
        if x != prev:
            out.append(x)
            prev = x
    return out

def all_trigram_counts(
    adata_full,
    state_col="ClusterID",
    group_cols=("sample_name", "TrackID"),
    time_col="position_t",
    collapse_bouts=True,
):
    """
    Pooled trigram counts across ALL tracks.
    Returns df with columns: trigram (tuple), trigram_str, end_state, count
    """
    obs = adata_full.obs
    needed = [state_col, *group_cols, time_col]
    missing = [c for c in needed if c not in obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    df = obs.loc[:, [*group_cols, time_col, state_col]].copy()
    df = df.sort_values([*group_cols, time_col], kind="stable")
    df["_state"] = df[state_col].astype(str).str.strip()

    trigram_counter = Counter()

    for _, g in df.groupby(list(group_cols), sort=False):
        seq = g["_state"].tolist()
        if collapse_bouts:
            seq = _collapse_runs(seq)
        if len(seq) < 3:
            continue
        for t in range(2, len(seq)):
            tri = (seq[t-2], seq[t-1], seq[t])
            trigram_counter[tri] += 1

    out = pd.DataFrame(
        [{"trigram": tri,
          "trigram_str": " > ".join(tri),
          "end_state": tri[2],
          "count": cnt}
         for tri, cnt in trigram_counter.items()]
    ).sort_values("count", ascending=False, kind="stable").reset_index(drop=True)

    return out

def top_trigrams_per_end_state(df_trigrams, top_n=15):
    """
    Returns a dict: end_state -> top-N DataFrame (sorted by count desc)
    """
    out = {}
    for end_state, sub in df_trigrams.groupby("end_state", sort=True):
        out[end_state] = sub.sort_values("count", ascending=False, kind="stable").head(top_n).reset_index(drop=True)
    return out

def plot_top_trigrams_per_end_state(df_trigrams, top_n=15, min_total_end_occurrences=1):
    """
    Makes one ranked plot per end_state (separate figures).
    min_total_end_occurrences filters out end_states with too few total counts.
    """
    # total number of trigrams ending in each state (i.e., how often that state is an end-state)
    totals = df_trigrams.groupby("end_state")["count"].sum().sort_values(ascending=False)

    for end_state in totals.index:
        if totals[end_state] < min_total_end_occurrences:
            continue

        sub = (
            df_trigrams[df_trigrams["end_state"] == end_state]
            .sort_values("count", ascending=False, kind="stable")
            .head(top_n)
            .copy()
        )
        if sub.empty:
            continue

        plt.figure(figsize=(12, max(4, 0.35 * len(sub))))
        plt.barh(sub["trigram_str"][::-1], sub["count"][::-1])
        plt.xlabel("Count (pooled across tracks)")
        plt.ylabel("Trigram")
        plt.title(f"Top {top_n} trigrams ending in state {end_state} (total={int(totals[end_state])})")
        plt.tight_layout()
        plt.show()

# -----------------------
# RUN (assumes adata_full already exists)
# -----------------------
df_trigrams = all_trigram_counts(
    adata_full,
    state_col="ClusterID",
    group_cols=("sample_name", "TrackID"),
    time_col="position_t",
    collapse_bouts=True,   # bout-based; set False for per-timepoint
)

# (A) Table: top-N per end-state (as a dict of DataFrames)
top_dict = top_trigrams_per_end_state(df_trigrams, top_n=15)

# Example: inspect one end-state
# display(top_dict["5"].head(15))

# (B) Plots: one figure per end-state
plot_top_trigrams_per_end_state(df_trigrams, top_n=15, min_total_end_occurrences=1)