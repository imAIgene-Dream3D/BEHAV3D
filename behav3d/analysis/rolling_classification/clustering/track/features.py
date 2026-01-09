import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Literal, List, Iterable
from pandas.api.types import is_numeric_dtype

from tqdm import tqdm

import numpy as np

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
from sklearn.feature_selection import VarianceThreshold


from behav3d.utils.rolling_classification import *

"""
Full track features based on states
"""
def rle_encode(states):
    """Return list of (state, run_length)."""
    if len(states) == 0:
        return []
    runs = []
    cur = states[0]
    length = 1
    for s in states[1:]:
        if s == cur:
            length += 1
        else:
            runs.append((cur, length))
            cur = s
            length = 1
    runs.append((cur, length))
    return runs

def fractions_from_runs(runs, states):
    total = sum(l for _, l in runs) if runs else 0
    cols = [f"overall_fraction_{s}" for s in states]
    out = {c: 0.0 for c in cols}
    if total == 0:
        return out, cols
    for s, l in runs:
        out[f"overall_fraction_{s}"] += l / total
    return out, cols

def bout_stats_from_runs(runs, states):
    lengths = {s: [] for s in states}
    for s, l in runs:
        lengths[s].append(l)

    total_bouts = len(runs)
    track_length = sum(l for _, l in runs)

    cols = []
    for s in states:
        cols.extend([f"bouts_nr_{s}", f"bouts_mean_length_{s}", f"bouts_max_length_{s}"])

    if total_bouts == 0 or track_length == 0:
        return {c: 0.0 for c in cols}, cols

    out = {}
    for s in states:
        arr = np.asarray(lengths[s], dtype=float)
        out[f"bouts_nr_{s}"] = float(arr.size) / float(total_bouts)
        out[f"bouts_mean_length_{s}"] = float(arr.mean()) / float(track_length) if arr.size else 0.0
        out[f"bouts_max_length_{s}"] = float(arr.max()) / float(track_length) if arr.size else 0.0

    return out, cols

def transition_probs_from_runs(runs, state_to_idx, states):
    K = len(states)
    M = np.zeros((K, K), dtype=float)
    labels = [s for s, _ in runs]
    if len(labels) < 2:
        return M  # all zeros
    for a, b in zip(labels[:-1], labels[1:]):
        if a in state_to_idx and b in state_to_idx:
            M[state_to_idx[a], state_to_idx[b]] += 1.0
    row_sums = M.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return M / row_sums

def ngram_counts_from_runs(runs, n=3, weight="count"):
    labels = [s for s, _ in runs]
    lens = [l for _, l in runs]
    out = {}
    if len(labels) < n:
        return out
    for i in range(len(labels) - n + 1):
        g = tuple(labels[i : i + n])
        w = float(lens[i]) if weight == "duration" else 1.0
        out[g] = out.get(g, 0.0) + w
    return out

def extract_descibing_track_state_features(
    adata,
    group_col=("sample_name", "TrackID"),
    time_col="position_t",
    state_col="ClusterID",
    use_fractions=True,
    use_bout_stats=True,
    use_transitions=True,
    use_bigrams=True,
    use_trigrams=True,
    ngram_weight="count",
):
    """
    Track-level feature extraction -> returns (track_adata, blocks)

    track_adata.X        = features (n_tracks x n_features)
    track_adata.obs      = per-track metadata (sample_name, TrackID, position_t_min/max, n_timepoints)
    track_adata.var_names= feature names
    """

    group_col = list(group_col) if isinstance(group_col, (list, tuple)) else [group_col]

    obs = adata.obs[group_col + [time_col, state_col]].copy()
    obs = obs.sort_values(group_col + [time_col])

    # Stable state universe
    states = pd.Index(obs[state_col].astype("category").cat.categories).tolist()
    state_to_idx = {s: i for i, s in enumerate(states)}

    # Use observed=True to avoid unobserved categorical groups exploding results
    gb = obs.groupby(group_col, sort=False, observed=True)

    # Collect runs + n-gram vocab
    runs_by_track = []
    bigram_set, trigram_set = set(), set()

    for tid, df_t in gb:
        seq = df_t[state_col].tolist()
        runs = rle_encode(seq)  # <-- your RLE function

        runs_by_track.append((tid, runs))

        if use_bigrams:
            bigram_set |= set(ngram_counts_from_runs(runs, n=2, weight=ngram_weight).keys())
        if use_trigrams:
            trigram_set |= set(ngram_counts_from_runs(runs, n=3, weight=ngram_weight).keys())

    bigram_list = sorted(bigram_set)
    trigram_list = sorted(trigram_set)

    def bg_col(g): return f"bigram_{g[0]}>{g[1]}"
    def tg_col(g): return f"trigram_{g[0]}>{g[1]}>{g[2]}"

    # Track block columns explicitly while building features
    fraction_cols = []
    bout_cols = []
    transition_cols = []
    bigram_cols = [bg_col(g) for g in bigram_list] if use_bigrams else []
    trigram_cols = [tg_col(g) for g in trigram_list] if use_trigrams else []

    rows = []
    id_rows = []

    for tid, runs in runs_by_track:
        # tid is scalar if len(group_col)==1, else tuple
        if len(group_col) == 1:
            id_rows.append([tid])
        else:
            id_rows.append(list(tid))

        feats = {}

        if use_fractions:
            f, fcols = fractions_from_runs(runs, states)
            feats.update(f)
            if not fraction_cols:
                fraction_cols = fcols

        if use_bout_stats:
            b, bcols = bout_stats_from_runs(runs, states)
            feats.update(b)
            if not bout_cols:
                bout_cols = bcols

        if use_transitions:
            T = transition_probs_from_runs(runs, state_to_idx, states)
            if not transition_cols:
                transition_cols = [f"transitions_{a}>{b}" for a in states for b in states]
            for a in states:
                ia = state_to_idx[a]
                for b in states:
                    ib = state_to_idx[b]
                    feats[f"transitions_{a}>{b}"] = float(T[ia, ib])

        if use_bigrams:
            bg = ngram_counts_from_runs(runs, n=2, weight=ngram_weight)
            total = sum(bg.values()) or 1.0
            for g in bigram_list:
                feats[bg_col(g)] = bg.get(g, 0.0) / total

        if use_trigrams:
            tg = ngram_counts_from_runs(runs, n=3, weight=ngram_weight)
            total = sum(tg.values()) or 1.0
            for g in trigram_list:
                feats[tg_col(g)] = tg.get(g, 0.0) / total

        rows.append(feats)

    # Features
    df_feat = pd.DataFrame(rows).fillna(0.0)

    # IDs in the exact same order as df_feat rows
    df_ids = pd.DataFrame(id_rows, columns=group_col)

    # Time summaries (also observed=True to avoid unobserved category combos)
    df_time = (
        obs.groupby(group_col, sort=False, observed=True)[time_col]
           .agg(position_t_min="min", position_t_max="max", n_timepoints="size")
           .reset_index()
    )

    # Build obs table in the same row order as features
    df_meta = df_ids.merge(df_time, on=group_col, how="left")

    # Combine for conversion using your helper
    df_out = pd.concat([df_meta.reset_index(drop=True), df_feat.reset_index(drop=True)], axis=1)

    blocks = [fraction_cols, bout_cols, transition_cols, bigram_cols, trigram_cols]

    feature_cols = df_feat.columns.tolist()
    obs_cols = df_meta.columns.tolist()

    track_adata = df_to_adata(df_out, feature_cols=feature_cols, obs_cols=obs_cols)

    # Optional bookkeeping
    track_adata.uns["feature_blocks"] = {
        name: cols for name, cols in zip(
            ["fractions", "bout_stats", "transitions", "bigrams", "trigrams"], blocks
        ) if cols
    }
    track_adata.uns["feature_params"] = {
        "group_col": group_col,
        "time_col": time_col,
        "state_col": state_col,
        "ngram_weight": ngram_weight,
        "use_fractions": use_fractions,
        "use_bout_stats": use_bout_stats,
        "use_transitions": use_transitions,
        "use_bigrams": use_bigrams,
        "use_trigrams": use_trigrams,
    }

    return track_adata, blocks


"""
Accessory feautre selection or scaling methods
"""

def scale_feature_blocks(
    adata,
    blocks,
    layer=None,
    mode="sqrt"  # "sqrt" or "linear"
):
    for block in blocks:
        d = len(block)
        if d == 0:
            continue

        scale = np.sqrt(d) if mode == "sqrt" else float(d)

        idx = adata.var_names.get_indexer(block)

        if np.any(idx < 0):
            raise ValueError("Some block features not found in adata.var_names")

        if layer is not None:
            adata.layers[layer][:, idx] /= scale
        else:
            adata.X[:, idx] /= scale

    return adata

def l2_normalize_block(adata, cols, layer=None):
    """Row-wise L2 normalize selected columns. Leaves all-zero rows unchanged."""
    adata_norm = adata.copy()
    if layer is not None:
        X = adata_norm[:, cols].layers[layer]
    else:
        X = adata_norm[:, cols].X
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata_norm[:, cols] = X / norms
    return adata_norm

def l2_normalize_all(adata):
    """Row-wise L2 normalize all columns. Leaves all-zero rows unchanged."""
   
    X = adata.X.copy()
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata.X = X / norms
    return(adata)

def l2_normalize_features_blocks(
    adata,
    blocks,
    layer=None
    ):
    for block in blocks:
        adata = l2_normalize_block(adata, block, layer=layer)
    return adata
 
def drop_highly_correlated_features(
    data,
    feature_cols: list,
    threshold: float = 0.98,
    prefer_keep: list | None = None,
    also_drop_constant: bool = True,
    redundancy_tolerance_abs: float = 0.02,       # apply prefer_keep only if |Δredundancy| <= this
    redundancy_tolerance_rel: float = 0.05,
) -> tuple[object, list, pd.DataFrame]:
    """
    Remove one feature from any pair with |corr| > threshold.

    Accepts either:
      - pandas DataFrame (features are columns)
      - AnnData (features are var_names; uses .X via .to_df())

    Returns:
      reduced_data : same type as input (DataFrame or AnnData), with dropped features removed
      dropped      : list of dropped feature names
      report       : DataFrame logging each decision: kept, dropped, correlation, reason
    """
    is_adata = hasattr(data, "var_names") and hasattr(data, "obs_names")

    # Build working feature frame X (pandas) and restrict to existing features
    if is_adata:
        existing = [f for f in feature_cols if f in data.var_names]
        X = data[:, existing].to_df() if existing else pd.DataFrame(index=data.obs_names)
        feature_cols_eff = existing
    else:
        existing = [f for f in feature_cols if f in data.columns]
        X = data[existing].copy() if existing else pd.DataFrame(index=data.index)
        feature_cols_eff = existing

    # Your current behavior: ignore passed prefer_keep and auto-prefer median/mean features
    prefer_keep = [f for f in feature_cols_eff if "median" in f]
    prefer_keep += [f for f in feature_cols_eff if "mean" in f]
    prefer_keep = set(prefer_keep)

    # Numeric-only
    cols_ok = [c for c in X.columns if is_numeric_dtype(X[c])]
    X = X[cols_ok]

    # Replace infs and optionally drop constants
    X = X.replace([np.inf, -np.inf], np.nan)
    dropped: list[str] = []
    decisions: list[dict] = []

    if also_drop_constant and X.shape[1] > 0:
        nunique = X.nunique(dropna=True)
        const_cols = nunique[nunique <= 1].index.tolist()
        for c in const_cols:
            decisions.append({
                "kept_feature": None,
                "dropped_feature": c,
                "abs_correlation": np.nan,
                "reason": "constant_or_single_value"
            })
        if const_cols:
            X = X.drop(columns=const_cols)
            dropped.extend(const_cols)

    if X.shape[1] <= 1:
        report = pd.DataFrame(decisions, columns=["kept_feature", "dropped_feature", "abs_correlation", "reason"])
        if is_adata:
            keep_vars = [v for v in data.var_names if v not in set(dropped)]
            return data[:, keep_vars].copy(), dropped, report
        else:
            reduced_df = data.drop(columns=[c for c in dropped if c in data.columns])
            return reduced_df, dropped, report

    # Impute NaNs so corr works
    X_impute = X.fillna(X.mean(numeric_only=True))

    # Abs correlation and upper triangle
    corr = X_impute.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

    to_drop = set()

    for col in upper.columns:
        partners = upper.index[upper[col] > threshold].tolist()
        for partner in partners:
            if col in to_drop or partner in to_drop:
                continue

            corr_val = float(upper.loc[partner, col])

            col_redund = corr[col].drop(index=[col]).mean()
            part_redund = corr[partner].drop(index=[partner]).mean()
            diff = abs(col_redund - part_redund)

            rel_ok = diff <= redundancy_tolerance_rel * max(col_redund, part_redund, 1e-12)
            abs_ok = diff <= redundancy_tolerance_abs

            if abs_ok or rel_ok:
                if (col in prefer_keep) ^ (partner in prefer_keep):
                    kept = col if col in prefer_keep else partner
                    loser = partner if col in prefer_keep else col
                    reason = "prefer_keep_tie"
                else:
                    if col_redund > part_redund:
                        kept, loser = partner, col
                    elif part_redund > col_redund:
                        kept, loser = col, partner
                    else:
                        kept, loser = sorted([col, partner])[0], sorted([col, partner])[1]
                    reason = "more_redundant"
            else:
                if col_redund >= part_redund:
                    kept, loser = partner, col
                else:
                    kept, loser = col, partner
                reason = "more_redundant"

            to_drop.add(loser)
            decisions.append({
                "kept_feature": kept,
                "dropped_feature": loser,
                "abs_correlation": corr_val,
                "reason": reason,
                "col_redundancy": float(col_redund),
                "partner_redundancy": float(part_redund),
                "redundancy_diff": float(diff),
            })

    dropped.extend(sorted(to_drop))

    report = pd.DataFrame(decisions, columns=["kept_feature", "dropped_feature", "abs_correlation", "reason"])
    report = report.sort_values(
        by=["reason", "abs_correlation"],
        ascending=[True, False],
        na_position="last"
    ).reset_index(drop=True)

    if is_adata:
        keep_vars = [v for v in data.var_names if v not in set(dropped)]
        return data[:, keep_vars].copy(), dropped, report
    else:
        reduced_df = data.drop(columns=[c for c in dropped if c in data.columns])
        return reduced_df, dropped, report


def drop_low_variance_features(
    data,
    feature_cols: list[str],
    low_var_threshold= 1e-4,
) -> tuple[object, list[str], list[str]]:
    """
    Drop low-variance features using sklearn VarianceThreshold.

    Accepts either:
      - pandas DataFrame (features are columns)
      - AnnData (features are var_names; uses .X via .to_df())

    Returns:
      reduced_data      : same type as input (DataFrame or AnnData)
      kept_features     : list of kept feature names
      dropped_low_var   : list of dropped feature names
    """
    is_adata = hasattr(data, "var_names") and hasattr(data, "obs_names")

    if is_adata:
        existing = [f for f in feature_cols if f in data.var_names]
        X = data[:, existing].to_df() if existing else pd.DataFrame(index=data.obs_names)
    else:
        existing = [f for f in feature_cols if f in data.columns]
        X = data[existing].copy() if existing else pd.DataFrame(index=data.index)

    if not existing or X.shape[1] == 0:
        return data, [], []

    selector = VarianceThreshold(threshold=low_var_threshold)
    selector.fit(X)

    keep_mask = selector.get_support()
    kept_features = X.columns[keep_mask].tolist()
    dropped_low_var = X.columns[~keep_mask].tolist()

    if is_adata:
        keep_vars = [v for v in data.var_names if v not in set(dropped_low_var)]
        return data[:, keep_vars].copy(), kept_features, dropped_low_var
    else:
        reduced_df = data.drop(columns=dropped_low_var)
        return reduced_df, kept_features, dropped_low_var
    
# def drop_highly_correlated_features(
#     df: pd.DataFrame,
#     feature_cols: list,
#     threshold: float = 0.98,
#     prefer_keep: list | None = None,
#     also_drop_constant: bool = True,
#     redundancy_tolerance_abs: float = 0.02,       # apply prefer_keep only if |Δredundancy| <= this
#     redundancy_tolerance_rel: float = 0.05,
# ) -> tuple[pd.DataFrame, list, pd.DataFrame]:
#     """
#     Remove one feature from any pair with |corr| > threshold.
#     Returns:
#       X_reduced : DataFrame of remaining features
#       dropped   : list of dropped feature names
#       report    : DataFrame logging each decision: kept, dropped, correlation, reason
#     """
#     # prefer_keep = set(prefer_keep or [])
#     prefer_keep = [f for f in feature_cols if "median" in f]
#     prefer_keep += [f for f in feature_cols if "mean" in f]
#     # Work on a numeric-only copy (keep order)
#     X = df[feature_cols].copy()
#     # Cast non-numeric to numeric where possible, otherwise drop
#     cols_ok = [c for c in X.columns if is_numeric_dtype(X[c])]
#     X = X[cols_ok]

#     # Replace infs and optionally drop constants
#     X = X.replace([np.inf, -np.inf], np.nan)
#     dropped: list[str] = []
#     decisions: list[dict] = []

#     if also_drop_constant:
#         nunique = X.nunique(dropna=True)
#         const_cols = nunique[nunique <= 1].index.tolist()
#         for c in const_cols:
#             decisions.append({
#                 "kept_feature": None,
#                 "dropped_feature": c,
#                 "abs_correlation": np.nan,
#                 "reason": "constant_or_single_value"
#             })
#         if const_cols:
#             X = X.drop(columns=const_cols)
#             dropped.extend(const_cols)

#     if X.shape[1] <= 1:
#         report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
#         return X, dropped, report

#     # Impute NaNs (mean) so corr works; PCA later will use its own scaling
#     X_impute = X.fillna(X.mean(numeric_only=True))

#     # Absolute correlation matrix and its upper triangle
#     corr = X_impute.corr().abs()
#     upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

#     to_drop = set()

#     # Greedy pass over pairs above threshold
#     for col in upper.columns:
#         partners = upper.index[upper[col] > threshold].tolist()
#         for partner in partners:
#             if col in to_drop or partner in to_drop:
#                 continue

#             corr_val = float(upper.loc[partner, col])

#             # compute redundancies (mean abs corr to others)
#             col_redund    = corr[col].drop(index=[col]).mean()
#             part_redund   = corr[partner].drop(index=[partner]).mean()
#             diff = abs(col_redund - part_redund)
#             rel_ok = diff <= redundancy_tolerance_rel * max(col_redund, part_redund, 1e-12)
#             abs_ok = diff <= redundancy_tolerance_abs

#             if abs_ok or rel_ok:
#                 # tie → apply prefer_keep if it singles out one; else fall back to redundancy anyway
#                 if (col in prefer_keep) ^ (partner in prefer_keep):
#                     kept  = col if col in prefer_keep else partner
#                     loser = partner if col in prefer_keep else col
#                     reason = "prefer_keep_tie"
#                 else:
#                     # tie but no clear preference → arbitrarily drop the slightly more redundant (or lexical tiebreak)
#                     if col_redund >= part_redund:
#                         kept, loser = partner, col
#                     else:
#                         kept, loser = col, partner
#                     reason = "more_redundant"
#             else:
#                 # not a tie → drop the more redundant one (ignore prefer_keep)
#                 if col_redund >= part_redund:
#                     kept, loser = partner, col
#                 else:
#                     kept, loser = col, partner
#                 reason = "more_redundant"

#             to_drop.add(loser)
#             decisions.append({
#                 "kept_feature": kept,
#                 "dropped_feature": loser,
#                 "abs_correlation": corr_val,
#                 "reason": reason,
#                 "col_redundancy": float(col_redund),
#                 "partner_redundancy": float(part_redund),
#                 "redundancy_diff": float(diff),
#             })

#     if to_drop:
#         df = df.drop(columns=list(to_drop))
#         dropped.extend(list(to_drop))

#     report = pd.DataFrame(decisions, columns=["kept_feature","dropped_feature","abs_correlation","reason"])
#     # Sort report: constants first, then by correlation desc
#     report = report.sort_values(
#         by=["reason", "abs_correlation"],
#         ascending=[True, False],
#         na_position="last"
#     ).reset_index(drop=True)

#     return df, dropped, report

# def drop_low_variance_features(
#     df: pd.DataFrame,
#     feature_cols: list[str],
#     low_var_threshold: float,
# ) -> tuple[pd.DataFrame, list[str], list[str]]:
#     """
#     Drop low-variance features using sklearn VarianceThreshold.

#     Returns:
#       df_reduced        : df with low-variance features removed
#       kept_features     : list of kept feature names
#       dropped_low_var   : list of dropped feature names
#     """
#     selector = VarianceThreshold(threshold=low_var_threshold)
#     selector.fit(df[feature_cols])

#     keep_mask = selector.get_support()
#     kept_features = df[feature_cols].columns[keep_mask].tolist()
#     dropped_low_var = df[feature_cols].columns[~keep_mask].tolist()

#     df_reduced = df.drop(columns=dropped_low_var)

#     return df_reduced, kept_features, dropped_low_var

