import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering

try:
    from dtaidistance import dtw_ndim
except Exception:
    dtw_ndim = None

try:
    import umap
except Exception:
    umap = None

from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.utils import _winfo


def _normalize_state_cols(state_cols):
    if isinstance(state_cols, str):
        state_cols = [state_cols]
    out = [str(c).strip() for c in list(state_cols) if str(c).strip() != ""]
    if len(out) == 0:
        raise ValueError("At least one categorical state column is required for distance clustering.")
    return out


def _format_state_token(row_values, missing_token="missing"):
    tokens = []
    for value in row_values:
        if pd.isna(value):
            tokens.append(str(missing_token))
        else:
            tokens.append(str(value))
    return "|".join(tokens)


def _local_dissimilarity_name(value):
    if callable(value):
        return getattr(value, "__name__", value.__class__.__name__)
    return str(value)


def extract_categorical_track_sequences(
    adata,
    *,
    state_cols=(FULL_STATE_COL,),
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    missing_policy="keep",
    missing_token="missing",
    extra_meta_cols=(),
):
    """Return per-track categorical state sequences and track-level metadata.

    extra_meta_cols are carried through as track-level metadata (first value per
    track) without affecting how tracks are grouped, e.g. a per-track-constant
    population tag such as ``origin_cell_type``.
    """
    state_cols = _normalize_state_cols(state_cols)
    groupby_cols = [str(c) for c in list(groupby_cols)]
    extra_meta_cols = [
        str(c) for c in list(extra_meta_cols) if str(c) not in groupby_cols and str(c) in adata.obs.columns
    ]
    required_cols = list(groupby_cols) + [str(time_col)] + list(state_cols) + extra_meta_cols
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required columns for distance clustering: {missing}")

    obs = adata.obs[required_cols].copy()
    obs["_orig_idx"] = np.arange(len(obs))
    obs = obs.sort_values(groupby_cols + [str(time_col), "_orig_idx"])

    missing_policy = str(missing_policy).strip().lower()
    if missing_policy not in {"keep", "drop"}:
        raise ValueError("missing_policy must be 'keep' or 'drop'.")

    sequences = []
    meta_rows = []
    grouped = obs.groupby(groupby_cols, sort=False, observed=True)
    for track_id, df_track in grouped:
        if missing_policy == "drop":
            df_seq = df_track.dropna(subset=state_cols).copy()
        else:
            df_seq = df_track.copy()
        if len(df_seq) == 0:
            continue

        if len(state_cols) == 1:
            series = df_seq[state_cols[0]]
            seq = [str(missing_token) if pd.isna(value) else str(value) for value in series.tolist()]
        else:
            seq = [
                _format_state_token(row, missing_token=missing_token)
                for row in df_seq[state_cols].itertuples(index=False, name=None)
            ]

        if len(groupby_cols) == 1:
            track_values = [track_id]
        else:
            track_values = list(track_id)
        meta = {col: track_values[i] for i, col in enumerate(groupby_cols)}
        for col in extra_meta_cols:
            vals = df_seq[col].dropna().unique()
            meta[col] = vals[0] if len(vals) > 0 else pd.NA
        meta["position_t_min"] = df_seq[str(time_col)].min()
        meta["position_t_max"] = df_seq[str(time_col)].max()
        meta["n_timepoints"] = int(len(df_seq))
        meta["distance_sequence"] = " ".join(seq)
        meta_rows.append(meta)
        sequences.append(seq)

    if len(sequences) == 0:
        raise ValueError("No valid track sequences were available for clustering.")
    return sequences, pd.DataFrame(meta_rows)


def _one_hot_encode_sequences(sequences, *, dtype=np.double):
    n_tracks = len(sequences)
    if n_tracks == 0:
        raise ValueError("No sequences available for distance computation.")

    categories = sorted({str(value) for seq in sequences for value in seq})
    if len(categories) == 0:
        raise ValueError("No categorical labels available for distance computation.")

    category_to_index = {label: index for index, label in enumerate(categories)}
    encoded = []
    for seq in sequences:
        seq_values = [str(value) for value in seq]
        arr = np.zeros((len(seq_values), len(categories)), dtype=dtype)
        for row_index, label in enumerate(seq_values):
            arr[row_index, category_to_index[label]] = 1.0
        encoded.append(arr)
    return encoded, categories


def _validate_distance_matrix(distances, n_expected):
    distances = np.asarray(distances, dtype=float)
    if distances.shape != (int(n_expected), int(n_expected)):
        raise ValueError(
            f"Distance matrix has shape {distances.shape}, expected {(int(n_expected), int(n_expected))}."
        )
    if not np.isfinite(distances).all():
        raise ValueError("Distance matrix contains non-finite values.")
    distances = 0.5 * (distances + distances.T)
    np.fill_diagonal(distances, 0.0)
    return distances


def compute_dtaidistance_onehot_distance_matrix(
    sequences,
    *,
    window=None,
    max_dist=None,
    max_length_diff=None,
    penalty=None,
    psi=None,
    parallel=True,
    inner_dist="squared euclidean",
    verbose=True,
):
    """Compute a DTW distance matrix using dtaidistance over one-hot vectors."""
    if dtw_ndim is None:
        raise ImportError(
            "dtaidistance is required for behavioral trajectory classification. "
            "Install the BEHAV3D requirements in the active notebook kernel."
        )
    encoded_sequences, categories = _one_hot_encode_sequences(sequences)
    if bool(verbose):
        _winfo(
            "trajectory-dtai",
            "dtaidistance distance matrix | "
            f"tracks={len(encoded_sequences)} | states={len(categories)} | "
            f"inner_dist={inner_dist} | window={window} | parallel={bool(parallel)}",
        )
    distances = dtw_ndim.distance_matrix_fast(
        encoded_sequences,
        compact=False,
        parallel=bool(parallel),
        inner_dist=str(inner_dist),
        window=None if window is None else int(window),
        max_dist=max_dist,
        max_length_diff=max_length_diff,
        penalty=penalty,
        psi=psi,
    )
    return _validate_distance_matrix(distances, len(encoded_sequences)), categories


def compute_transition_profile_distance_matrix(sequences, *, use_bigrams=True, use_trigrams=False, verbose=True):
    """Distance matrix based on which state transitions a track ever exhibits.

    Unlike the one-hot DTW distance above (which aligns sequences timepoint by
    timepoint and is dominated by how much time is spent in each state), this
    represents each track as the *set* of state-to-state transitions (bigrams
    and/or trigrams) it ever exhibits, then measures Jaccard distance between
    those sets. Two tracks that both contain a brief A->B->A blip end up at
    distance 0 from each other regardless of when the blip occurs or how long
    the surrounding runs of A are, and both are pulled away from a track that
    never transitions at all — the timepoint-alignment approach instead scores
    a single-timepoint blip as nearly identical to "no event" and often fails
    to separate them.
    """
    from behav3d.features.state_descriptive_features import rle_encode, ngram_counts_from_runs
    from scipy.spatial.distance import pdist, squareform

    profiles = []
    vocab = set()
    for seq in sequences:
        runs = rle_encode([str(value) for value in seq])
        observed = set()
        if use_bigrams:
            observed.update(ngram_counts_from_runs(runs, n=2, weight="count").keys())
        if use_trigrams:
            observed.update(ngram_counts_from_runs(runs, n=3, weight="count").keys())
        profiles.append(observed)
        vocab.update(observed)

    n = len(sequences)
    vocab = sorted(vocab, key=str)
    if bool(verbose):
        _winfo(
            "trajectory-dtai",
            "transition-profile distance matrix | "
            f"tracks={n} | vocabulary_size={len(vocab)} | "
            f"use_bigrams={use_bigrams} | use_trigrams={use_trigrams}",
        )
    if not vocab:
        return np.zeros((n, n)), []

    vocab_index = {g: i for i, g in enumerate(vocab)}
    presence = np.zeros((n, len(vocab)), dtype=float)
    for i, observed in enumerate(profiles):
        for g in observed:
            presence[i, vocab_index[g]] = 1.0

    distances = squareform(pdist(presence, metric="jaccard"))
    distances = np.where(np.isnan(distances), 0.0, distances)
    np.fill_diagonal(distances, 0.0)
    return _validate_distance_matrix(distances, n), vocab


def _cluster_precomputed_distances(distances, *, n_clusters=6, linkage="average"):
    n_clusters = int(n_clusters)
    if n_clusters < 2:
        raise ValueError("n_clusters must be at least 2.")
    if distances.shape[0] < n_clusters:
        raise ValueError(
            f"n_clusters={n_clusters} exceeds the number of available tracks ({distances.shape[0]})."
        )
    if str(linkage).strip().lower() == "ward":
        raise ValueError("Ward linkage is not valid for precomputed DTAI/DTW distances. Use average, complete, or single.")
    try:
        model = AgglomerativeClustering(
            n_clusters=n_clusters,
            metric="precomputed",
            linkage=str(linkage),
        )
    except TypeError:
        model = AgglomerativeClustering(
            n_clusters=n_clusters,
            affinity="precomputed",
            linkage=str(linkage),
        )
    raw_labels = model.fit_predict(distances)
    return raw_labels, model


def _cluster_precomputed_distances_leiden(
    distances,
    *,
    n_neighbors=15,
    resolution=1.0,
    random_state=123,
):
    """Cluster a precomputed DTW distance matrix with Leiden graph clustering.

    Unlike agglomerative clustering with a fixed n_clusters, Leiden finds
    community structure on the same k-NN graph that the QC UMAP embedding
    (`_ensure_dtaidistance_umap`) is built from, so partitions tend to track
    the visual blobs in that plot rather than cutting through them.
    """
    from behav3d.analysis.behavior.general.leiden import run_leiden_clustering

    distances = np.asarray(distances, dtype=float)
    n_obs = int(distances.shape[0])
    if n_obs < 3:
        raise ValueError("At least three tracks are required for Leiden clustering.")
    resolved_n_neighbors = max(2, min(int(n_neighbors), n_obs - 1))
    raw_labels = np.asarray(
        run_leiden_clustering(
            distances,
            n_neighbors=resolved_n_neighbors,
            metric="precomputed",
            method="umap",
            use_rep="X",
            resolution=resolution,
            random_state=int(random_state),
            key_added="leiden_cluster",
        )
    )
    if len(set(raw_labels.tolist())) < 2:
        raise ValueError(
            f"Leiden clustering collapsed to a single cluster at resolution={resolution}; "
            "try a higher resolution."
        )
    return raw_labels, None


def _ensure_dtaidistance_umap(adata_tracks, distances, *, random_state=123, n_neighbors=15, min_dist=0.1):
    """Compute (or reuse) the 2D UMAP embedding for a dtaidistance track adata.

    The cached embedding is only reused when the requested UMAP parameters
    match those it was computed with — otherwise a stale embedding from a
    previous n_neighbors/min_dist would silently be returned even after the
    user changed those controls.
    """
    requested_params = {
        "n_neighbors": int(n_neighbors),
        "min_dist": float(min_dist),
        "random_state": int(random_state),
    }
    cached_params = adata_tracks.uns.get("_umap_params")
    if "X_umap" in adata_tracks.obsm and cached_params == requested_params:
        return np.asarray(adata_tracks.obsm["X_umap"], dtype=float)
    if umap is None:
        raise ImportError("umap-learn is required to create DTAI UMAP quality-control plots.")
    n_obs = int(adata_tracks.n_obs)
    if n_obs < 2:
        raise ValueError("At least two tracks are required for UMAP.")
    reducer = umap.UMAP(
        n_components=2,
        metric="precomputed",
        n_neighbors=max(2, min(requested_params["n_neighbors"], n_obs - 1)),
        min_dist=requested_params["min_dist"],
        random_state=requested_params["random_state"],
    )
    embedding = reducer.fit_transform(np.asarray(distances, dtype=float))
    adata_tracks.obsm["X_umap"] = np.asarray(embedding, dtype=float)
    adata_tracks.uns["_umap_params"] = requested_params
    return adata_tracks.obsm["X_umap"]

def _relabel_by_cluster_size(raw_labels):
    labels = pd.Series(raw_labels).astype(str)
    ranked = labels.value_counts().sort_values(ascending=False)
    mapping = {old: str(i + 1) for i, old in enumerate(ranked.index.tolist())}
    return labels.map(mapping).astype(str).tolist(), mapping


def _add_cluster_medoids(adata_tracks, distances, cluster_key="ClusterID"):
    labels = pd.Series(adata_tracks.obs[cluster_key], index=adata_tracks.obs.index).astype(str)
    medoid_flags = pd.Series(False, index=adata_tracks.obs.index)
    medoid_rank = pd.Series(pd.NA, index=adata_tracks.obs.index, dtype="Int64")
    for cluster in sorted(labels.unique(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x))):
        idx = np.flatnonzero(labels.to_numpy() == cluster)
        if len(idx) == 0:
            continue
        within = distances[np.ix_(idx, idx)]
        order = np.argsort(within.mean(axis=1))
        for rank, local_idx in enumerate(order, start=1):
            obs_idx = adata_tracks.obs.index[idx[local_idx]]
            medoid_rank.loc[obs_idx] = int(rank)
        medoid_flags.loc[adata_tracks.obs.index[idx[order[0]]]] = True
    adata_tracks.obs[f"{cluster_key}_medoid"] = medoid_flags
    adata_tracks.obs[f"{cluster_key}_medoid_rank"] = medoid_rank
    return adata_tracks
