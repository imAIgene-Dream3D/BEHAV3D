"""One-hot categorical track clustering preview utilities using dtaidistance.

This module mirrors the DTW preview widget, but encodes categorical behavioral
states as one-hot vectors and computes pairwise distances with dtaidistance.
"""

from pathlib import Path
import re
import time
import traceback

import anndata as ad
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

try:
    from dtaidistance import dtw_ndim
except Exception:
    dtw_ndim = None

try:
    import umap
except Exception:
    umap = None

from behav3d.analysis.clustering.state.classification import FULL_STATE_COL
from behav3d.analysis.clustering.track.classification import (
    export_track_cluster_backprojection,
    rename_track_clusters,
    show_track_cluster_backprojection,
)
from behav3d.analysis.tcell_analysis import (
    plot_cluster_percentage_bars,
    plot_clustering_feature_heatmap,
    plot_feature_umap,
    run_tcell_analysis,
)
from behav3d.analysis.clustering.track.visualization.plots.exemplar_coordinate_utils import (
    ensure_exemplar_coordinate_columns,
)
from behav3d.analysis.clustering.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    select_exemplar_tracks_by_cluster,
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.widgets.utils import spinning_loader


def _winfo(prefix, message):
    print(f"[{prefix}] INFO {message}")


def _sanitize_filename_token(value, fallback="track"):
    token = str(value).strip()
    if token == "":
        token = str(fallback)
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", token).strip("._-")
    return token if token else str(fallback)


def _ordered_unique(values):
    out = []
    seen = set()
    for value in list(values or []):
        text = str(value).strip()
        if text == "" or text in seen:
            continue
        out.append(text)
        seen.add(text)
    return out


def _resolve_dtaidistance_paths(output_dir, cell_type, output_subdir_name="behavorial_trajectories"):
    root = Path(output_dir).expanduser()
    analysis_outdir = root / "analysis" / str(cell_type)
    state_outdir = analysis_outdir / "behavioral_states"
    outfolder = analysis_outdir / str(output_subdir_name)
    clustering_outfolder = outfolder / "clustering"
    quality_control_outfolder = outfolder / "quality_control"
    outfolder.mkdir(parents=True, exist_ok=True)
    clustering_outfolder.mkdir(parents=True, exist_ok=True)
    quality_control_outfolder.mkdir(parents=True, exist_ok=True)
    return {
        "root": root,
        "analysis_outdir": analysis_outdir,
        "state_outdir": state_outdir,
        "outfolder": outfolder,
        "clustering_outfolder": clustering_outfolder,
        "quality_control_outfolder": quality_control_outfolder,
    }


def get_dtaidistance_track_trajectories_filename(cell_type):
    cell_token = _sanitize_filename_token(cell_type, fallback="cell")
    return f"BEHAV3D_{cell_token}_behavioral_trajectories_dtaidistance.h5ad"


def _default_behavioral_states_path(output_dir, cell_type):
    return (
        Path(output_dir).expanduser()
        / "analysis"
        / str(cell_type)
        / "behavioral_states"
        / f"BEHAV3D_{cell_type}_behavioral_states.h5ad"
    )


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


def _resolve_optional_int(value):
    text = str(value).strip()
    return None if text == "" else int(text)


def _resolve_optional_float(value):
    text = str(value).strip()
    return None if text == "" else float(text)


def _filter_tracks_for_dtaidistance(
    adata,
    *,
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    trajectory_size=None,
    min_length=None,
    trim_mode="last",
):
    groupby_cols = [str(c) for c in list(groupby_cols)]
    missing = [col for col in groupby_cols if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing groupby_cols in adata.obs: {missing}")
    if str(time_col) not in adata.obs.columns:
        raise KeyError(f"'{time_col}' not found in adata.obs.")

    obs = adata.obs.copy()
    obs["_orig_idx"] = np.arange(len(obs))
    obs = obs.sort_values(groupby_cols + [str(time_col), "_orig_idx"])

    if min_length is not None:
        group_sizes = obs.groupby(groupby_cols, observed=True).size()
        keep_groups = group_sizes[group_sizes >= int(min_length)].index
        obs["_keep"] = obs.set_index(groupby_cols).index.isin(keep_groups)
        obs = obs.loc[obs["_keep"]].copy()

    if trajectory_size is not None:
        trajectory_size = int(trajectory_size)
        if trajectory_size <= 0:
            raise ValueError("trajectory_size must be positive when supplied.")
        trim_mode = str(trim_mode).strip().lower()
        if trim_mode not in {"first", "last"}:
            raise ValueError("trim_mode must be 'first' or 'last'.")
        obs["_rank"] = obs.groupby(groupby_cols, observed=True).cumcount()
        sizes = obs.groupby(groupby_cols, observed=True)["_rank"].transform("max") + 1
        if trim_mode == "first":
            obs = obs.loc[obs["_rank"] < trajectory_size]
        else:
            obs = obs.loc[obs["_rank"] >= (sizes - trajectory_size)]

    idx = obs.index
    adata_out = adata[idx].copy()
    for col in ["_orig_idx", "_keep", "_rank"]:
        if col in adata_out.obs.columns:
            adata_out.obs.drop(columns=[col], inplace=True)
    return adata_out


def extract_categorical_track_sequences(
    adata,
    *,
    state_cols=(FULL_STATE_COL,),
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    missing_policy="keep",
    missing_token="missing",
):
    """Return per-track categorical state sequences and track-level metadata."""
    state_cols = _normalize_state_cols(state_cols)
    groupby_cols = [str(c) for c in list(groupby_cols)]
    required_cols = list(groupby_cols) + [str(time_col)] + list(state_cols)
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


def _ensure_dtaidistance_umap(adata_tracks, distances, *, random_state=123, n_neighbors=15, min_dist=0.1):
    if "X_umap" in adata_tracks.obsm:
        return np.asarray(adata_tracks.obsm["X_umap"], dtype=float)
    if umap is None:
        raise ImportError("umap-learn is required to create DTAI UMAP quality-control plots.")
    n_obs = int(adata_tracks.n_obs)
    if n_obs < 2:
        raise ValueError("At least two tracks are required for UMAP.")
    reducer = umap.UMAP(
        n_components=2,
        metric="precomputed",
        n_neighbors=max(2, min(int(n_neighbors), n_obs - 1)),
        min_dist=float(min_dist),
        random_state=int(random_state),
    )
    embedding = reducer.fit_transform(np.asarray(distances, dtype=float))
    adata_tracks.obsm["X_umap"] = np.asarray(embedding, dtype=float)
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


def _save_diagnostics(
    adata_tracks,
    distances,
    outfolder,
    *,
    cluster_key="ClusterID",
    max_heatmap_tracks=200,
    random_state=123,
):
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)
    pdf_path = outfolder / "dtaidistance_clustering_diagnostics.pdf"
    counts_csv = outfolder / "dtaidistance_cluster_counts.csv"
    medoids_csv = outfolder / "dtaidistance_cluster_medoids.csv"
    umap_csv = outfolder / "dtaidistance_umap_clusters.csv"

    labels = pd.Series(adata_tracks.obs[cluster_key]).astype(str)
    counts = labels.value_counts().rename_axis(cluster_key).reset_index(name="n_tracks")
    counts.to_csv(counts_csv, index=False)

    medoid_cols = [
        c
        for c in ["sample_name", "TrackID", cluster_key, f"{cluster_key}_medoid", f"{cluster_key}_medoid_rank"]
        if c in adata_tracks.obs.columns
    ]
    adata_tracks.obs.loc[adata_tracks.obs[f"{cluster_key}_medoid"].astype(bool), medoid_cols].to_csv(
        medoids_csv,
        index=False,
    )
    umap_embedding = None
    umap_error = None
    try:
        umap_embedding = _ensure_dtaidistance_umap(
            adata_tracks,
            distances,
            random_state=int(random_state),
        )
        umap_df = adata_tracks.obs.copy()
        umap_df["UMAP1"] = umap_embedding[:, 0]
        umap_df["UMAP2"] = umap_embedding[:, 1]
        umap_df.to_csv(umap_csv, index=False)
    except Exception as exc:
        umap_error = str(exc)

    with PdfPages(pdf_path) as pdf:
        if umap_embedding is not None:
            fig, ax = plt.subplots(figsize=(7, 6))
            cluster_order = sorted(labels.unique().tolist(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
            cmap = plt.get_cmap("tab20")
            color_map = {cluster: cmap(i % cmap.N) for i, cluster in enumerate(cluster_order)}
            for cluster in cluster_order:
                mask = labels.to_numpy() == str(cluster)
                ax.scatter(
                    umap_embedding[mask, 0],
                    umap_embedding[mask, 1],
                    s=18,
                    alpha=0.8,
                    color=color_map[cluster],
                    label=str(cluster),
                    linewidths=0,
                )
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            ax.set_title("One-hot dtaidistance UMAP by cluster")
            ax.legend(title=str(cluster_key), loc="best", frameon=False, markerscale=1.4)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        elif umap_error is not None:
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.axis("off")
            ax.text(0.02, 0.6, f"UMAP plot skipped: {umap_error}", ha="left", va="center", wrap=True)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(counts) + 1.5)))
        ax.barh(counts[cluster_key].astype(str), counts["n_tracks"], color="#2E6FBA")
        ax.invert_yaxis()
        ax.set_xlabel("N tracks")
        ax.set_ylabel(str(cluster_key))
        ax.set_title("One-hot dtaidistance cluster counts")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        n = distances.shape[0]
        if n > 0:
            rng = np.random.default_rng(int(random_state))
            if n > int(max_heatmap_tracks):
                keep = np.sort(rng.choice(n, size=int(max_heatmap_tracks), replace=False))
                heat = distances[np.ix_(keep, keep)]
                title_suffix = f"random {int(max_heatmap_tracks)}/{n} tracks"
            else:
                heat = distances
                title_suffix = f"{n} tracks"
            fig, ax = plt.subplots(figsize=(8, 7))
            im = ax.imshow(heat, aspect="auto", interpolation="nearest", cmap="viridis")
            ax.set_title(f"One-hot dtaidistance matrix ({title_suffix})")
            ax.set_xlabel("Track")
            ax.set_ylabel("Track")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    return {
        "diagnostics_pdf": str(pdf_path),
        "cluster_counts_csv": str(counts_csv),
        "cluster_medoids_csv": str(medoids_csv),
        "umap_csv": str(umap_csv) if umap_embedding is not None else None,
        "umap_error": umap_error,
    }


def _dtai_meta(adata_tracks):
    meta = adata_tracks.uns.get("dtai_trajectory_clustering", {})
    return meta if isinstance(meta, dict) else {}


def _resolve_cluster_key(adata_tracks, cluster_key=None):
    if cluster_key is not None:
        return str(cluster_key)
    return str(_dtai_meta(adata_tracks).get("cluster_key", "ClusterID"))


def _write_model_if_requested(adata_tracks, output_dir, cell_type, *, save_outputs=True):
    output_path = (
        _resolve_dtaidistance_paths(output_dir, cell_type)["outfolder"]
        / get_dtaidistance_track_trajectories_filename(cell_type)
    )
    if bool(save_outputs):
        adata_tracks.write(output_path, compression="gzip")
    return output_path


def _build_identity_cluster_mapping_from_obs(obs, cluster_col="ClusterID"):
    if cluster_col not in obs.columns:
        return {}
    values = pd.Series(obs[cluster_col]).dropna().astype(str).unique().tolist()
    values = sorted(values, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
    return {str(value): str(value) for value in values}


def _feature_dtw_outdir(output_dir, cell_type):
    return Path(output_dir).expanduser() / "analysis" / str(cell_type) / "timepoint_feature_dtw"


def _feature_dtw_clustered_csv_path(output_dir, cell_type):
    return _feature_dtw_outdir(output_dir, cell_type) / f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv"


def _feature_dtw_umap_csv_path(output_dir, cell_type):
    return _feature_dtw_outdir(output_dir, cell_type) / f"BEHAV3D_{cell_type}_UMAP_clusters.csv"


def _feature_dtw_output_csv_paths(output_dir, cell_type):
    outdir = _feature_dtw_outdir(output_dir, cell_type)
    return [
        outdir / f"BEHAV3D_{cell_type}_UMAP_clusters.csv",
        outdir / f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv",
        outdir / f"BEHAV3D_{cell_type}_UMAP_cluster_percentages.csv",
    ]


def _feature_dtw_rename_mapping_path(output_dir, cell_type):
    return _feature_dtw_outdir(output_dir, cell_type) / f"feature_dtw_cluster_names_{cell_type}.yml"


def _load_feature_dtw_name_mapping(output_dir, cell_type):
    mapping_path = _feature_dtw_rename_mapping_path(output_dir, cell_type)
    if not mapping_path.exists():
        return {}
    try:
        with mapping_path.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
        mapping = data.get("cluster_names", data)
        if not isinstance(mapping, dict):
            return {}
        return {str(k): str(v) for k, v in mapping.items()}
    except Exception:
        return {}


def _feature_dtw_plot_info_cols(df):
    preferred = ["sample_name", "TrackID", "ClusterID", "ClusterName", "UMAP1", "UMAP2"]
    return [col for col in preferred if col in df.columns]


def _feature_dtw_sample_cols(df):
    sample_cols = [col for col in ["sample_name"] if col in df.columns]
    return sample_cols if len(sample_cols) > 0 else None


def _feature_dtw_cluster_percentage_groups(df):
    if "ClusterID" not in df.columns:
        return None
    group_cols = [col for col in ["sample_name", "condition", "treatment", "well"] if col in df.columns]
    if len(group_cols) == 0:
        counts = df["ClusterID"].astype(str).value_counts().rename_axis("ClusterID").reset_index(name="count")
        total = counts["count"].sum()
        counts["percentage"] = counts["count"] / total if total else 0.0
        return counts
    counts = df.groupby(group_cols + ["ClusterID"], observed=True).size().reset_index(name="count")
    totals = counts.groupby(group_cols, observed=True)["count"].transform("sum")
    counts["percentage"] = np.where(totals > 0, counts["count"] / totals, 0.0)
    return counts


def _save_feature_dtw_quality_control(output_dir, cell_type, *, outfolder=None):
    outdir = _feature_dtw_outdir(output_dir, cell_type)
    qc_outdir = outdir / "quality_control" if outfolder is None else Path(outfolder)
    qc_outdir.mkdir(parents=True, exist_ok=True)
    umap_csv = _feature_dtw_umap_csv_path(output_dir, cell_type)
    if not umap_csv.exists():
        raise FileNotFoundError(f"Original BEHAV3D UMAP CSV not found: {umap_csv}")

    df_umap = pd.read_csv(umap_csv, low_memory=False)
    df_plot = df_umap.copy()
    if "ClusterName" in df_plot.columns:
        df_plot["ClusterID"] = pd.Categorical(df_plot["ClusterName"].astype(str))
    elif "ClusterID" in df_plot.columns:
        df_plot["ClusterID"] = pd.Categorical(df_plot["ClusterID"].astype(str))
    else:
        raise ValueError("Original BEHAV3D UMAP CSV is missing ClusterID.")

    sample_cols = _feature_dtw_sample_cols(df_plot)
    info_cols = _feature_dtw_plot_info_cols(df_plot)
    cluster_umap_path = qc_outdir / f"BEHAV3D_{cell_type}_UMAP_clusters.pdf"
    heatmap_path = qc_outdir / f"BEHAV3D_{cell_type}_UMAP_cluster_feature_heatmap.pdf"
    perc_prefix = qc_outdir / f"BEHAV3D_{cell_type}_UMAP_cluster_percentages"

    plot_feature_umap(
        df_umap=df_plot,
        info_cols=info_cols,
        sample_cols=sample_cols,
        outpath=cluster_umap_path,
        rows_per_page=4,
        nr_cols=2,
        rows_first_img=2,
        figsize=(8.27, 11.69),
        plot_results=True,
    )
    plot_clustering_feature_heatmap(
        df_plot,
        info_cols,
        sample_cols,
        heatmap_path,
        rows_per_page=7,
        nr_cols=2,
        figsize=(8.27, 11.69),
        plot_results=True,
    )
    df_clust_perc = _feature_dtw_cluster_percentage_groups(df_plot)
    if df_clust_perc is not None:
        df_clust_perc.to_csv(perc_prefix.with_suffix(".csv"), index=False)
        plot_cluster_percentage_bars(df_clust_perc, perc_prefix, group_by_columns=None, plot_results=True)

    return {
        "umap_pdf": str(cluster_umap_path),
        "heatmap_pdf": str(heatmap_path),
        "percentages_pdf": str(perc_prefix.with_suffix(".pdf")),
    }


def save_dtaidistance_diagnostics(
    adata_tracks,
    output_dir,
    cell_type,
    *,
    cluster_key=None,
    outfolder=None,
    max_heatmap_tracks=200,
    random_state=123,
    save_outputs=True,
    verbose=True,
):
    """Write diagnostics for a saved one-hot dtaidistance clustering model."""
    paths = _resolve_dtaidistance_paths(output_dir, cell_type)
    resolved_cluster_key = _resolve_cluster_key(adata_tracks, cluster_key=cluster_key)
    distances = _validate_distance_matrix(adata_tracks.X, adata_tracks.n_obs)
    plot_paths = _save_diagnostics(
        adata_tracks,
        distances,
        paths["quality_control_outfolder"] if outfolder is None else outfolder,
        cluster_key=resolved_cluster_key,
        max_heatmap_tracks=int(max_heatmap_tracks),
        random_state=int(random_state),
    )
    adata_tracks.uns.setdefault("visualization", {})
    adata_tracks.uns["visualization"].update(plot_paths)
    output_path = _write_model_if_requested(
        adata_tracks,
        output_dir,
        cell_type,
        save_outputs=save_outputs,
    )
    if bool(verbose):
        _winfo("trajectory-dtai", f"saved diagnostics: {plot_paths.get('diagnostics_pdf')}")
        if bool(save_outputs):
            _winfo("trajectory-dtai", f"updated one-hot dtaidistance model: {output_path}")
    return plot_paths


def _load_filtered_state_adata_for_model(
    adata_tracks,
    output_dir,
    cell_type,
    *,
    adata_full_path=None,
    verbose=True,
):
    meta = _dtai_meta(adata_tracks)
    if adata_full_path is None:
        adata_full_path = meta.get("source_adata_full_path", None)
    if adata_full_path is None:
        adata_full_path = _default_behavioral_states_path(output_dir, cell_type)
    adata_full_path = Path(adata_full_path).expanduser()
    if not adata_full_path.exists():
        raise FileNotFoundError(f"Could not find behavioral-state h5ad: {adata_full_path}")

    groupby_cols = list(meta.get("groupby_cols", ["sample_name", "TrackID"]))
    time_col = str(meta.get("time_col", "position_t"))
    behavioral_trajectory_size = meta.get("behavioral_trajectory_size", 100)
    if behavioral_trajectory_size is not None:
        behavioral_trajectory_size = int(behavioral_trajectory_size)
    min_track_length = meta.get("min_track_length", behavioral_trajectory_size)
    if min_track_length is not None:
        min_track_length = int(min_track_length)
    trim_mode = str(meta.get("trajectory_trim_mode", "last"))
    if bool(verbose):
        _winfo("trajectory-dtai", f"loading behavioral states for exemplar plots: {adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)
    return _filter_tracks_for_dtaidistance(
        adata_full,
        groupby_cols=groupby_cols,
        time_col=time_col,
        trajectory_size=behavioral_trajectory_size,
        min_length=min_track_length,
        trim_mode=trim_mode,
    )


def save_dtaidistance_exemplar_plots(
    adata_tracks,
    output_dir,
    cell_type,
    *,
    adata_full_path=None,
    cluster_key=None,
    n_per_cluster=None,
    exemplar_pdf_rows_per_page=6,
    exemplar_layout_mode="both",
    exemplar_num_example_ranks=5,
    make_overview_statebars=True,
    make_backprojection_pdf=False,
    make_backprojection_mp4=False,
    random_state=None,
    save_outputs=True,
    verbose=True,
):
    """Write exemplar overview and per-cluster statebar PDFs for a DTAI model."""
    paths = _resolve_dtaidistance_paths(output_dir, cell_type)
    meta = _dtai_meta(adata_tracks)
    resolved_cluster_key = _resolve_cluster_key(adata_tracks, cluster_key=cluster_key)
    state_col = str(meta.get("state_col", FULL_STATE_COL))
    time_col = str(meta.get("time_col", "position_t"))
    if n_per_cluster is None:
        n_per_cluster = int(meta.get("n_per_cluster", 10))
    if random_state is None:
        random_state = int(meta.get("random_state", 123))

    adata_filt = _load_filtered_state_adata_for_model(
        adata_tracks,
        output_dir,
        cell_type,
        adata_full_path=adata_full_path,
        verbose=verbose,
    )
    exemplar_root = paths["clustering_outfolder"] / "example_tracks"
    exemplar_root.mkdir(parents=True, exist_ok=True)

    coord_enrichment = ensure_exemplar_coordinate_columns(
        adata_filt,
        output_dir=output_dir,
        cell_type=cell_type,
        require_pixel_for_video=False,
    )
    chosen_exemplars, _ = select_exemplar_tracks_by_cluster(
        adata_tracks=adata_tracks,
        n_per_cluster=int(n_per_cluster),
        sample_key="sample_name",
        track_key="TrackID",
        cluster_key=resolved_cluster_key,
        tmin_key="position_t_min",
        tmax_key="position_t_max",
        seed=int(random_state),
    )
    exemplar_selection_csv = exemplar_root / (
        f"example_track_selection_cluster_{_sanitize_filename_token(resolved_cluster_key, fallback='cluster')}_"
        f"state_{_sanitize_filename_token(state_col, fallback='state')}_dtai.csv"
    )
    chosen_exemplars.to_csv(exemplar_selection_csv, index=False)

    overview_pdf = None
    per_cluster_pdf_out = {}
    backprojection_pdf_out = {}
    backprojection_mp4_out = {}
    if bool(make_overview_statebars):
        fig_exemplar, _, _ = plot_exemplar_tracks_by_cluster(
            adata_filt,
            adata_tracks,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            time_key=time_col,
            state_key=state_col,
            cluster_key=resolved_cluster_key,
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            seed=int(random_state),
        )
        overview_pdf = exemplar_root / "example_tracks_overview.pdf"
        with PdfPages(overview_pdf) as pdf:
            pdf.savefig(fig_exemplar, bbox_inches="tight", dpi=300)
        plt.close(fig_exemplar)

        per_cluster_pdf_out = save_exemplar_statebar_track_pdf_per_cluster(
            adata_full=adata_filt,
            out_dir=exemplar_root,
            chosen_df=chosen_exemplars,
            adata_tracks=None,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            time_key=time_col,
            state_key=state_col,
            cluster_key=resolved_cluster_key,
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            rows_per_page=int(exemplar_pdf_rows_per_page),
            plot_dpi=300,
            seed=int(random_state),
            cmap_name="tab20",
            layout_mode=str(exemplar_layout_mode),
            num_example_ranks=int(exemplar_num_example_ranks),
        )

    if bool(make_backprojection_pdf):
        backprojection_pdf_out = save_exemplar_statebar_backprojection_pdf(
            adata_full=adata_filt,
            output_dir=output_dir,
            cell_type=cell_type,
            out_dir=exemplar_root / "backprojection",
            chosen_df=chosen_exemplars,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            time_key=time_col,
            state_key=state_col,
            cluster_key=resolved_cluster_key,
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            plot_dpi=220,
            seed=int(random_state),
            examples_per_cluster=int(n_per_cluster),
            num_example_ranks=int(exemplar_num_example_ranks),
            verbose=verbose,
        )

    if bool(make_backprojection_mp4):
        backprojection_mp4_out = save_exemplar_statebar_backprojection_video_per_cluster(
            adata_full=adata_filt,
            output_dir=output_dir,
            cell_type=cell_type,
            out_dir=exemplar_root / "backprojection",
            chosen_df=chosen_exemplars,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            time_key=time_col,
            state_key=state_col,
            cluster_key=resolved_cluster_key,
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            dpi=180,
            seed=int(random_state),
            examples_per_cluster=int(n_per_cluster),
            num_example_ranks=int(exemplar_num_example_ranks),
            verbose=verbose,
        )

    plot_paths = {
        "exemplar_selection_csv": str(exemplar_selection_csv),
        "exemplar_tracks_overview_pdf": str(overview_pdf) if overview_pdf is not None else None,
        "exemplar_statebar_track_pdf_by_cluster": dict(per_cluster_pdf_out.get("pdf_paths_by_cluster", {})),
        "exemplar_statebar_track_pdf_by_example_rank": dict(
            per_cluster_pdf_out.get("pdf_paths_by_example_rank", {})
        ),
        "exemplar_backprojection_pdf_by_cluster": dict(backprojection_pdf_out.get("pdf_paths_by_cluster", {})),
        "exemplar_backprojection_pdf_by_example_rank": dict(
            backprojection_pdf_out.get("pdf_paths_by_example_rank", {})
        ),
        "exemplar_backprojection_mp4_by_cluster": dict(backprojection_mp4_out.get("video_paths_by_cluster", {})),
        "exemplar_backprojection_mp4_by_example_rank": dict(
            backprojection_mp4_out.get("video_paths_by_example_rank", {})
        ),
        "exemplar_render_config": {
            "stage": "after_clustering",
            "cluster_key": resolved_cluster_key,
            "state_key": state_col,
            "n_per_cluster": int(n_per_cluster),
            "coordinate_enrichment": dict(coord_enrichment),
            "overview_statebars_enabled": bool(make_overview_statebars),
            "backprojection_pdfs_enabled": bool(make_backprojection_pdf),
            "backprojection_videos_enabled": bool(make_backprojection_mp4),
        },
    }
    adata_tracks.uns.setdefault("visualization", {})
    adata_tracks.uns["visualization"].update(plot_paths)
    output_path = _write_model_if_requested(
        adata_tracks,
        output_dir,
        cell_type,
        save_outputs=save_outputs,
    )
    if bool(verbose):
        _winfo("trajectory-dtai", f"saved exemplar track PDFs: {exemplar_root}")
        if bool(save_outputs):
            _winfo("trajectory-dtai", f"updated one-hot dtaidistance model: {output_path}")
    return plot_paths


def run_categorical_dtaidistance_trajectory_clustering(
    output_dir,
    cell_type="tcell",
    *,
    adata_full_path=None,
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    behavioral_trajectory_size=100,
    min_track_length=None,
    trajectory_trim_mode="last",
    max_tracks=None,
    n_clusters=6,
    window=None,
    max_dist=None,
    max_length_diff=None,
    penalty=None,
    psi=None,
    parallel=True,
    linkage="average",
    missing_policy="keep",
    cluster_key="ClusterID",
    save_outputs=True,
    save_distance_matrix=False,
    plot_results=False,
    max_heatmap_tracks=200,
    plot_exemplars=False,
    n_per_cluster=10,
    exemplar_pdf_rows_per_page=6,
    exemplar_layout_mode="both",
    exemplar_num_example_ranks=5,
    output_subdir_name="behavorial_trajectories",
    random_state=123,
    verbose=True,
):
    """Cluster categorical state trajectories with dtaidistance over one-hot encodings."""
    started = time.perf_counter()
    paths = _resolve_dtaidistance_paths(output_dir, cell_type, output_subdir_name=output_subdir_name)
    if adata_full_path is None:
        adata_full_path = _default_behavioral_states_path(output_dir, cell_type)
    adata_full_path = Path(adata_full_path).expanduser()
    if not adata_full_path.exists():
        raise FileNotFoundError(f"Could not find behavioral-state h5ad: {adata_full_path}")

    state_cols = [str(FULL_STATE_COL)]
    trajectory_size = None if behavioral_trajectory_size is None else int(behavioral_trajectory_size)
    min_length = trajectory_size if min_track_length is None else int(min_track_length)
    trim_mode = str(trajectory_trim_mode).strip().lower()

    if bool(verbose):
        _winfo("trajectory-dtai", f"loading behavioral states: {adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)

    if bool(verbose):
        length_text = "all timepoints" if trajectory_size is None else f"{trim_mode} {trajectory_size} timepoints"
        _winfo(
            "trajectory-dtai",
            f"filtering tracks with min_length={min_length}, keeping {length_text} | "
            f"state_col={FULL_STATE_COL}",
        )
    adata_filt = _filter_tracks_for_dtaidistance(
        adata_full,
        groupby_cols=list(groupby_cols),
        time_col=str(time_col),
        trajectory_size=trajectory_size,
        min_length=min_length,
        trim_mode=trim_mode,
    )

    sequences, track_obs = extract_categorical_track_sequences(
        adata_filt,
        state_cols=state_cols,
        groupby_cols=groupby_cols,
        time_col=time_col,
        missing_policy=missing_policy,
    )

    if max_tracks is not None and int(max_tracks) > 0 and len(sequences) > int(max_tracks):
        rng = np.random.default_rng(int(random_state))
        keep = np.sort(rng.choice(len(sequences), size=int(max_tracks), replace=False))
        sequences = [sequences[i] for i in keep]
        track_obs = track_obs.iloc[keep].reset_index(drop=True)
        if bool(verbose):
            _winfo("trajectory-dtai", f"sampled {len(sequences)} tracks for distance clustering")

    if bool(verbose):
        _winfo(
            "trajectory-dtai",
            f"computing pairwise one-hot dtaidistance distances for {len(sequences)} tracks",
        )
    distances, categories = compute_dtaidistance_onehot_distance_matrix(
        sequences,
        window=window,
        max_dist=max_dist,
        max_length_diff=max_length_diff,
        penalty=penalty,
        psi=psi,
        parallel=parallel,
        verbose=verbose,
    )
    dtw_backend = "dtaidistance"

    if bool(verbose):
        _winfo("trajectory-dtai", f"clustering precomputed distances with n_clusters={int(n_clusters)}")
    raw_labels, _ = _cluster_precomputed_distances(
        distances,
        n_clusters=int(n_clusters),
        linkage=str(linkage),
    )
    labels, size_mapping = _relabel_by_cluster_size(raw_labels)
    track_obs[cluster_key] = pd.Categorical(labels)

    var_names = [f"distance_to_track_{i}" for i in range(distances.shape[0])]
    adata_tracks = ad.AnnData(
        X=distances,
        obs=track_obs.copy(),
        var=pd.DataFrame(index=pd.Index(var_names, name="distance_target")),
    )
    adata_tracks = _add_cluster_medoids(adata_tracks, distances, cluster_key=cluster_key)

    silhouette = None
    try:
        if len(set(labels)) > 1 and len(set(labels)) < len(labels):
            silhouette = float(silhouette_score(distances, labels, metric="precomputed"))
    except Exception:
        silhouette = None

    adata_tracks.uns["dtai_trajectory_clustering"] = {
        "method": "categorical_onehot_dtaidistance_agglomerative",
        "dtw_backend": str(dtw_backend),
        "local_encoding": "one_hot",
        "inner_dist": "squared euclidean",
        "one_hot_categories": [str(c) for c in categories],
        "groupby_cols": [str(c) for c in list(groupby_cols)],
        "time_col": str(time_col),
        "state_col": str(FULL_STATE_COL),
        "state_cols": list(state_cols),
        "behavioral_trajectory_size": None if trajectory_size is None else int(trajectory_size),
        "min_track_length": None if min_length is None else int(min_length),
        "trajectory_trim_mode": trim_mode,
        "max_tracks": None if max_tracks is None else int(max_tracks),
        "n_clusters": int(n_clusters),
        "window": None if window is None else int(window),
        "max_dist": None if max_dist is None else float(max_dist),
        "max_length_diff": None if max_length_diff is None else int(max_length_diff),
        "penalty": None if penalty is None else float(penalty),
        "psi": None if psi is None else int(psi),
        "parallel": bool(parallel),
        "linkage": str(linkage),
        "missing_policy": str(missing_policy),
        "cluster_key": str(cluster_key),
        "raw_label_size_mapping": dict(size_mapping),
        "silhouette_score_precomputed": silhouette,
        "n_per_cluster": int(n_per_cluster),
        "random_state": int(random_state),
        "source_adata_full_path": str(adata_full_path),
    }

    if bool(plot_results):
        plot_paths = _save_diagnostics(
            adata_tracks,
            distances,
            paths["quality_control_outfolder"],
            cluster_key=cluster_key,
            max_heatmap_tracks=int(max_heatmap_tracks),
            random_state=int(random_state),
        )
        adata_tracks.uns.setdefault("visualization", {})
        adata_tracks.uns["visualization"].update(plot_paths)

    if bool(plot_exemplars):
        exemplar_root = paths["clustering_outfolder"] / "example_tracks"
        exemplar_root.mkdir(parents=True, exist_ok=True)
        try:
            ensure_exemplar_coordinate_columns(
                adata_filt,
                output_dir=output_dir,
                cell_type=cell_type,
                require_pixel_for_video=False,
            )
            chosen_exemplars, _ = select_exemplar_tracks_by_cluster(
                adata_tracks=adata_tracks,
                n_per_cluster=int(n_per_cluster),
                sample_key="sample_name",
                track_key="TrackID",
                cluster_key=str(cluster_key),
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                seed=int(random_state),
            )
            exemplar_selection_csv = exemplar_root / (
                f"example_track_selection_cluster_{_sanitize_filename_token(cluster_key, fallback='cluster')}_"
                f"state_{_sanitize_filename_token(state_cols[0], fallback='state')}_dtai.csv"
            )
            chosen_exemplars.to_csv(exemplar_selection_csv, index=False)

            fig_exemplar, _, _ = plot_exemplar_tracks_by_cluster(
                adata_filt,
                adata_tracks,
                n_per_cluster=int(n_per_cluster),
                sample_key="sample_name",
                track_key="TrackID",
                time_key=str(time_col),
                state_key=str(state_cols[0]),
                cluster_key=str(cluster_key),
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                seed=int(random_state),
            )
            overview_pdf = exemplar_root / "example_tracks_overview.pdf"
            with PdfPages(overview_pdf) as pdf:
                pdf.savefig(fig_exemplar, bbox_inches="tight", dpi=300)
            plt.close(fig_exemplar)

            try:
                per_cluster_pdf_out = save_exemplar_statebar_track_pdf_per_cluster(
                    adata_full=adata_filt,
                    out_dir=exemplar_root,
                    chosen_df=chosen_exemplars,
                    adata_tracks=None,
                    n_per_cluster=int(n_per_cluster),
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key=str(time_col),
                    state_key=str(state_cols[0]),
                    cluster_key=str(cluster_key),
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    rows_per_page=int(exemplar_pdf_rows_per_page),
                    plot_dpi=300,
                    seed=int(random_state),
                    cmap_name="tab20",
                    layout_mode=str(exemplar_layout_mode),
                    num_example_ranks=int(exemplar_num_example_ranks),
                )
                pdf_paths_by_cluster = dict(per_cluster_pdf_out.get("pdf_paths_by_cluster", {}))
                pdf_paths_by_example_rank = dict(per_cluster_pdf_out.get("pdf_paths_by_example_rank", {}))
                exemplar_warning = None
            except Exception as exc:
                raise RuntimeError(f"Failed to save exemplar PDFs: {exc}") from exc
            adata_tracks.uns.setdefault("visualization", {})
            adata_tracks.uns["visualization"].update(
                {
                    "exemplar_selection_csv": str(exemplar_selection_csv),
                    "exemplar_tracks_overview_pdf": str(overview_pdf),
                    "exemplar_statebar_track_pdf_by_cluster": dict(pdf_paths_by_cluster),
                    "exemplar_statebar_track_pdf_by_example_rank": dict(pdf_paths_by_example_rank),
                    "exemplar_statebar_warning": exemplar_warning,
                }
            )
            if bool(verbose):
                _winfo("trajectory-dtai", f"saved exemplar track PDFs: {exemplar_root}")
        except Exception as exc:
            adata_tracks.uns.setdefault("visualization", {})
            adata_tracks.uns["visualization"]["exemplar_error"] = str(exc)
            if bool(verbose):
                _winfo("trajectory-dtai", f"skipping exemplar PDFs due to error: {exc}")

    if bool(save_distance_matrix):
        distance_csv = paths["clustering_outfolder"] / "categorical_dtai_distance_matrix.csv"
        pd.DataFrame(distances, index=adata_tracks.obs.index, columns=adata_tracks.obs.index).to_csv(distance_csv)
        adata_tracks.uns.setdefault("dtai_trajectory_clustering", {})
        adata_tracks.uns["dtai_trajectory_clustering"]["distance_matrix_csv"] = str(distance_csv)

    output_path = paths["outfolder"] / get_dtaidistance_track_trajectories_filename(cell_type)
    if bool(save_outputs):
        adata_tracks.write(output_path, compression="gzip")
        if bool(verbose):
            _winfo("trajectory-dtai", f"saved one-hot dtaidistance model: {output_path}")

    elapsed = time.perf_counter() - started
    if bool(verbose):
        _winfo("trajectory-dtai", f"finished in {elapsed:.2f}s")
    return adata_tracks


class TrajectoryDTAIClassificationPanel:
    """Notebook widget for one-hot categorical track clustering."""

    def __init__(self, metadata_loader, cell_type=None):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(getattr(metadata_loader, "output_dir", "")).expanduser())
        self.model_adata = None

        cell_types = self._detect_cell_types()
        if len(cell_types) == 0:
            cell_types = ["tcell"]
        initial_cell_type = str(cell_type) if cell_type is not None else cell_types[0]
        if initial_cell_type not in cell_types:
            cell_types.append(initial_cell_type)

        self.cell_type_dd = widgets.Dropdown(
            description="Cell type",
            options=list(cell_types),
            value=initial_cell_type,
            layout=widgets.Layout(width="260px"),
            style={"description_width": "90px"},
        )
        self.refresh_btn = widgets.Button(description="Refresh", icon="refresh", layout=widgets.Layout(width="110px"))
        self.refresh_spinner = widgets.HTML(value=spinning_loader)
        self.refresh_spinner.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.output_dir_html = widgets.HTML("")
        self.use_original_behav3d = widgets.Checkbox(
            description="Use original feature-based BEHAV3D DTW clustering",
            value=False,
            indent=False,
            layout=widgets.Layout(width="420px"),
        )
        self.use_original_behav3d.observe(self._on_use_original_behav3d_changed, names="value")

        self.state_col_html = widgets.HTML("")
        self.behavioral_trajectory_size = widgets.Text(
            description="Trajectory size",
            value="100",
            placeholder="blank = all timepoints",
            layout=widgets.Layout(width="230px"),
            style={"description_width": "110px"},
        )
        self.n_clusters = widgets.IntText(
            description="Clusters",
            value=6,
            layout=widgets.Layout(width="170px"),
            style={"description_width": "80px"},
        )
        self.n_per_cluster = widgets.IntText(
            description="Exemplars/cluster",
            value=10,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "125px"},
        )
        self.random_state = widgets.IntText(
            description="Seed",
            value=123,
            layout=widgets.Layout(width="150px"),
            style={"description_width": "60px"},
        )
        self.advanced = widgets.Checkbox(description="Advanced", value=False, indent=False)
        self.advanced.observe(self._on_advanced_changed, names="value")

        self.window = widgets.Text(
            description="Window",
            value="",
            placeholder="blank = unconstrained",
            layout=widgets.Layout(width="220px"),
            style={"description_width": "80px"},
        )
        self.max_dist = widgets.Text(
            description="Max dist",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="180px"),
            style={"description_width": "80px"},
        )
        self.max_length_diff = widgets.Text(
            description="Max len diff",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="200px"),
            style={"description_width": "95px"},
        )
        self.penalty = widgets.Text(
            description="Penalty",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="170px"),
            style={"description_width": "70px"},
        )
        self.psi = widgets.Text(
            description="Psi",
            value="",
            placeholder="blank = off",
            layout=widgets.Layout(width="150px"),
            style={"description_width": "45px"},
        )
        self.parallel = widgets.Checkbox(description="Parallel", value=True, indent=False)
        self.linkage = widgets.Dropdown(
            description="Linkage",
            options=["average", "complete", "single"],
            value="average",
            layout=widgets.Layout(width="190px"),
            style={"description_width": "80px"},
        )
        self.missing_policy = widgets.Dropdown(
            description="Missing",
            options=[("Keep as category", "keep"), ("Drop missing timepoints", "drop")],
            value="keep",
            layout=widgets.Layout(width="260px"),
            style={"description_width": "80px"},
        )
        self.save_distance_matrix = widgets.Checkbox(description="Save distance matrix CSV", value=False, indent=False)
        self.trajectory_trim_mode = widgets.Dropdown(
            description="Trim",
            options=[("Last timepoints", "last"), ("First timepoints", "first")],
            value="last",
            layout=widgets.Layout(width="210px"),
            style={"description_width": "60px"},
        )

        self.btn_run = widgets.Button(
            description="Run one-hot dtaidistance clustering",
            button_style="success",
            layout=widgets.Layout(width="330px"),
        )
        self.run_spinner = widgets.HTML(value=spinning_loader)
        self.run_spinner.layout.display = "none"
        self.out_run = widgets.Output()
        self.plot_status_html = widgets.HTML("<i>Run clustering first, then create plots as needed.</i>")
        self.btn_diagnostics = widgets.Button(
            description="Create diagnostics plots",
            button_style="info",
            layout=widgets.Layout(width="230px"),
        )
        self.diagnostics_spinner = widgets.HTML(value=spinning_loader)
        self.diagnostics_spinner.layout.display = "none"
        self.btn_exemplars = widgets.Button(
            description="Create exemplar PDFs",
            button_style="info",
            layout=widgets.Layout(width="210px"),
        )
        self.exemplar_spinner = widgets.HTML(value=spinning_loader)
        self.exemplar_spinner.layout.display = "none"
        self.make_overview_statebars = widgets.Checkbox(
            description="Overview statebars",
            value=True,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.make_backprojection_pdf = widgets.Checkbox(
            description="Backprojection PDFs",
            value=False,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.make_backprojection_mp4 = widgets.Checkbox(
            description="Backprojection MP4",
            value=False,
            indent=False,
            layout=widgets.Layout(width="190px"),
        )
        self.state_outputs_warning = widgets.HTML("")
        self.out_plots = widgets.Output()
        self.rename_status = widgets.HTML("<i>Run clustering first to rename clusters.</i>")
        self.rename_rows = widgets.VBox([])
        self.btn_refresh_rename = widgets.Button(
            description="Refresh clusters",
            icon="refresh",
            layout=widgets.Layout(width="170px"),
        )
        self.btn_rename = widgets.Button(
            description="Apply cluster names",
            button_style="warning",
            layout=widgets.Layout(width="210px"),
        )
        self.rename_spinner = widgets.HTML(value=spinning_loader)
        self.rename_spinner.layout.display = "none"
        self.out_rename = widgets.Output()
        self._name_boxes = {}
        self.backprojection_status = widgets.HTML("<i>No samples detected yet.</i>")
        self.backproj_sample_dd = widgets.Dropdown(
            description="Sample",
            options=[],
            value=None,
            layout=widgets.Layout(width="360px"),
            style={"description_width": "90px"},
        )
        self.btn_backproject = widgets.Button(
            description="Open backprojection",
            button_style="success",
            layout=widgets.Layout(width="220px"),
        )
        self.backprojection_spinner = widgets.HTML(value=spinning_loader)
        self.backprojection_spinner.layout.display = "none"
        self.out_backprojection = widgets.Output()

        self.advanced_row_1 = widgets.HBox(
            [self.window, self.max_dist, self.max_length_diff, self.penalty, self.psi],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.advanced_row_2 = widgets.HBox(
            [
                self.trajectory_trim_mode,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
                self.use_original_behav3d,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.backend_summary_html = widgets.HTML(
            "<b>dtaidistance backend:</b> one-hot vectors with <code>inner_dist='squared euclidean'</code>; "
            "the DTW window is the main speed constraint."
        )

        self.run_section = widgets.VBox(
            [
                self.state_col_html,
                widgets.HBox(
                    [
                        self.behavioral_trajectory_size,
                        self.n_clusters,
                        self.random_state,
                        self.advanced,
                    ],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                self.advanced_row_1,
                self.advanced_row_2,
                self.backend_summary_html,
                widgets.HBox([self.btn_run, self.run_spinner]),
                self.out_run,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.plot_section = widgets.VBox(
            [
                self.plot_status_html,
                widgets.HTML(
                    "<b>Diagnostics</b><br>"
                    "<span style='color:#555;'>Writes cluster counts, medoids, and a distance-matrix heatmap.</span>"
                ),
                widgets.HBox([self.btn_diagnostics, self.diagnostics_spinner], layout=widgets.Layout(gap="8px")),
                widgets.HTML(
                    "<b>Exemplar tracks</b><br>"
                    "<span style='color:#555;'>Writes exemplar overview and per-cluster statebar PDFs.</span>"
                ),
                widgets.HBox(
                    [
                        self.n_per_cluster,
                        self.make_overview_statebars,
                        self.make_backprojection_pdf,
                        self.make_backprojection_mp4,
                        self.btn_exemplars,
                        self.exemplar_spinner,
                    ],
                    layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
                ),
                self.state_outputs_warning,
                self.out_plots,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.rename_section = widgets.VBox(
            [
                self.rename_status,
                self.rename_rows,
                widgets.HBox(
                    [self.btn_refresh_rename, self.btn_rename, self.rename_spinner],
                    layout=widgets.Layout(gap="8px"),
                ),
                self.out_rename,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.backprojection_section = widgets.VBox(
            [
                self.backprojection_status,
                self.backproj_sample_dd,
                widgets.HBox([self.btn_backproject, self.backprojection_spinner], layout=widgets.Layout(gap="8px")),
                self.out_backprojection,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.steps = widgets.Accordion(
            children=[self.run_section, self.rename_section, self.plot_section, self.backprojection_section],
            selected_index=0,
        )
        self.steps.set_title(0, "Run clustering")
        self.steps.set_title(1, "Rename clusters")
        self.steps.set_title(2, "Create plots")
        self.steps.set_title(3, "Backprojection")

        self.original_trajectory_size = widgets.IntText(
            description="Trajectory size",
            value=100,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "110px"},
        )
        self.original_n_clusters = widgets.IntText(
            description="Clusters",
            value=5,
            layout=widgets.Layout(width="170px"),
            style={"description_width": "80px"},
        )
        self.original_umap_n_neighbors = widgets.IntText(
            description="UMAP n_neighbors",
            value=15,
            layout=widgets.Layout(width="230px"),
            style={"description_width": "130px"},
        )
        self.original_umap_min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=0.1,
            layout=widgets.Layout(width="220px"),
            style={"description_width": "110px"},
        )
        self.btn_run_original = widgets.Button(
            description="Run original BEHAV3D",
            button_style="success",
            layout=widgets.Layout(width="230px"),
        )
        self.original_spinner = widgets.HTML(value=spinning_loader)
        self.original_spinner.layout.display = "none"
        self.out_original = widgets.Output()
        self.original_mode_switch_row = widgets.HBox(
            [self.use_original_behav3d],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.original_description_html = widgets.HTML(
            "<span style='color:#555;'>Original BEHAV3D clusters tracks from continuous timepoint features. "
            "For each track, selected features are normalized, DTW distances are calculated across feature "
            "trajectories, UMAP projects the distance structure, and K-means assigns the requested clusters.</span>"
        )
        self.original_controls_row = widgets.HBox(
            [
                self.original_trajectory_size,
                self.original_n_clusters,
                self.original_umap_n_neighbors,
                self.original_umap_min_dist,
            ],
            layout=widgets.Layout(flex_flow="row wrap", gap="8px"),
        )
        self.original_run_row = widgets.HBox(
            [self.btn_run_original, self.original_spinner],
            layout=widgets.Layout(gap="8px"),
        )
        self.original_run_section = widgets.VBox(
            [
                self.original_description_html,
                self.original_controls_row,
                self.original_run_row,
                self.out_original,
            ],
            layout=widgets.Layout(gap="8px"),
        )
        self.mode_body = widgets.VBox([self.steps], layout=widgets.Layout(gap="8px"))

        self.ui = widgets.VBox(
            [
                widgets.HTML('<div style="font-size:20px;font-weight:700;">Behavioral trajectory classification</div>'),
                widgets.HTML(
                    "<span style='color:#555;'>Classifies track super-behaviors using dynamic time warping "
                    "based on state transitions and proportions</span>"
                ),
                widgets.HBox([self.cell_type_dd, self.refresh_btn, self.refresh_spinner], layout=widgets.Layout(gap="8px")),
                self.status_html,
                self.mode_body,
            ],
            layout=widgets.Layout(gap="8px"),
        )

        self.cell_type_dd.observe(self._on_cell_type_changed, names="value")
        self.refresh_btn.on_click(self._on_refresh_clicked)
        self.btn_run.on_click(self._on_run_clicked)
        self.btn_run_original.on_click(self._on_run_original_clicked)
        self.btn_refresh_rename.on_click(self._on_refresh_rename_clicked)
        self.btn_rename.on_click(self._on_rename_clicked)
        self.btn_diagnostics.on_click(self._on_diagnostics_clicked)
        self.btn_exemplars.on_click(self._on_exemplars_clicked)
        self.make_overview_statebars.observe(self._on_exemplar_output_changed, names="value")
        self.make_backprojection_pdf.observe(self._on_exemplar_output_changed, names="value")
        self.make_backprojection_mp4.observe(self._on_exemplar_output_changed, names="value")
        self.btn_backproject.on_click(self._on_backproject_clicked)
        self._refresh_context()
        self._sync_advanced_visibility()
        self._sync_mode()

    def _detect_cell_types(self):
        md = getattr(self.metadata_loader, "metadata", None)
        cell_types = []
        if md is not None:
            try:
                cell_types.extend(detect_immune_cell_types_from_metadata(md))
                cell_types.extend(detect_organoid_types_from_metadata(md))
                cell_types.extend(detect_other_cell_types_from_metadata(md))
            except Exception:
                pass
        analysis_dir = Path(self.output_dir) / "analysis" if self.output_dir else None
        if analysis_dir is not None and analysis_dir.exists():
            for path in analysis_dir.iterdir():
                if path.is_dir():
                    cell_types.append(path.name)
        return _ordered_unique(cell_types)

    def _current_cell_type(self):
        return str(self.cell_type_dd.value).strip()

    def _model_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return (
            Path(self.output_dir)
            / "analysis"
            / ct
            / "behavorial_trajectories"
            / get_dtaidistance_track_trajectories_filename(ct)
        )

    def _original_track_features_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        feature_outdir = Path(self.output_dir) / "analysis" / ct / "track_features"
        filtered = feature_outdir / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        if filtered.exists():
            return filtered
        return feature_outdir / f"BEHAV3D_{ct}_combined_track_features.csv"

    def _state_adata_path(self, cell_type=None):
        ct = self._current_cell_type() if cell_type is None else str(cell_type)
        return _default_behavioral_states_path(self.output_dir, ct)

    def _has_behavioral_states(self):
        return Path(self._state_adata_path()).exists()

    def _any_exemplar_outputs_selected(self):
        return any(
            [
                bool(self.make_overview_statebars.value),
                bool(self.make_backprojection_pdf.value),
                bool(self.make_backprojection_mp4.value),
            ]
        )

    def _sync_exemplar_state_controls(self):
        has_states = self._has_behavioral_states()
        for widget in [
            self.make_overview_statebars,
            self.make_backprojection_pdf,
            self.make_backprojection_mp4,
        ]:
            widget.disabled = not has_states

        if not has_states:
            self.state_outputs_warning.value = (
                "<b style='color:#a66;'>No states have been defined.</b> "
                "Run behavioral state classification before creating overview statebars, "
                "backprojection PDFs, or backprojection MP4s."
            )
            self.btn_exemplars.disabled = True
            return

        self.state_outputs_warning.value = ""
        if bool(self.use_original_behav3d.value):
            has_clusters = _feature_dtw_clustered_csv_path(self.output_dir, self._current_cell_type()).exists()
        else:
            has_clusters = self.model_adata is not None or self._model_adata_path().exists()
        self.btn_exemplars.disabled = not (has_clusters and self._any_exemplar_outputs_selected())

    def _on_exemplar_output_changed(self, _):
        self._sync_exemplar_state_controls()

    def _sync_advanced_visibility(self):
        display = None if bool(self.advanced.value) else "none"
        self.advanced_row_1.layout.display = "none"
        self.advanced_row_2.layout.display = display
        self.backend_summary_html.layout.display = "none"

    def _on_advanced_changed(self, _):
        self._sync_advanced_visibility()

    def _sync_mode(self):
        if bool(self.use_original_behav3d.value):
            self.advanced_row_2.children = [
                self.trajectory_trim_mode,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
            ]
            self.original_run_section.children = [
                self.original_mode_switch_row,
                self.original_description_html,
                self.original_controls_row,
                self.original_run_row,
                self.out_original,
            ]
            self.steps.children = [
                self.original_run_section,
                self.rename_section,
                self.plot_section,
                self.backprojection_section,
            ]
            self.steps.set_title(0, "Run original BEHAV3D")
            self.steps.set_title(1, "Rename clusters")
            self.steps.set_title(2, "Create plots")
            self.steps.set_title(3, "Backprojection")
            self.status_html.value = (
                "<b>Original BEHAV3D mode:</b> running feature-based DTW on continuous track features."
            )
            self.btn_diagnostics.disabled = not _feature_dtw_umap_csv_path(
                self.output_dir,
                self._current_cell_type(),
            ).exists()
            self.btn_exemplars.disabled = True
            self.plot_status_html.value = "<i>Run original BEHAV3D clustering first, then create QC or rename clusters.</i>"
            self._rebuild_rename_rows()
            self._sync_exemplar_state_controls()
        else:
            self.advanced_row_2.children = [
                self.trajectory_trim_mode,
                self.linkage,
                self.parallel,
                self.save_distance_matrix,
                self.use_original_behav3d,
            ]
            self.original_run_section.children = [
                self.original_description_html,
                self.original_controls_row,
                self.original_run_row,
                self.out_original,
            ]
            self.steps.children = [self.run_section, self.rename_section, self.plot_section, self.backprojection_section]
            self.steps.set_title(0, "Run clustering")
            self.steps.set_title(1, "Rename clusters")
            self.steps.set_title(2, "Create plots")
            self.steps.set_title(3, "Backprojection")
            self._refresh_context()

    def _on_use_original_behav3d_changed(self, _):
        self._sync_mode()

    def _set_busy(self, button, spinner, busy=True):
        button.disabled = bool(busy)
        spinner.layout.display = None if bool(busy) else "none"

    def _refresh_context(self):
        self.output_dir = str(Path(getattr(self.metadata_loader, "output_dir", "")).expanduser())
        self.output_dir_html.value = ""
        current = self._current_cell_type()
        cell_types = self._detect_cell_types()
        if len(cell_types) == 0:
            cell_types = [current or "tcell"]
        if current and current not in cell_types:
            cell_types.append(current)
        self.cell_type_dd.options = _ordered_unique(cell_types)
        if self.cell_type_dd.value not in self.cell_type_dd.options:
            self.cell_type_dd.value = self.cell_type_dd.options[0]
        model_path = self._model_adata_path()
        if model_path.exists():
            if not bool(self.use_original_behav3d.value):
                self.status_html.value = "<b style='color:#080;'>Ready for plots:</b> an existing model is available."
            self.plot_status_html.value = "<b>Ready for plots:</b> an existing DTAI model is available."
            self.btn_diagnostics.disabled = False
            self.btn_exemplars.disabled = False
        else:
            if not bool(self.use_original_behav3d.value):
                self.status_html.value = (
                    "<i>No one-hot dtaidistance model detected yet. Press Run to one-hot encode the "
                    f"<code>{FULL_STATE_COL}</code> trajectories and create it.</i>"
                )
            self.plot_status_html.value = "<i>Run clustering first, then create plots as needed.</i>"
            self.btn_diagnostics.disabled = True
            self.btn_exemplars.disabled = True
        self._refresh_backprojection_samples()
        self._rebuild_rename_rows()
        if bool(self.use_original_behav3d.value):
            self.btn_diagnostics.disabled = not _feature_dtw_umap_csv_path(
                self.output_dir,
                self._current_cell_type(),
            ).exists()
            self._sync_exemplar_state_controls()
        else:
            self._sync_exemplar_state_controls()

    def _on_refresh_clicked(self, _):
        self._set_busy(self.refresh_btn, self.refresh_spinner, busy=True)
        try:
            self._refresh_context()
        finally:
            self._set_busy(self.refresh_btn, self.refresh_spinner, busy=False)

    def _on_cell_type_changed(self, _):
        self._refresh_context()

    def _load_model_adata(self):
        if self.model_adata is not None:
            return self.model_adata
        model_path = self._model_adata_path()
        if not model_path.exists():
            raise FileNotFoundError(f"No one-hot dtaidistance model found at: {model_path}")
        self.model_adata = sc.read_h5ad(model_path)
        return self.model_adata

    def _detect_backprojection_sample_names(self):
        names = []
        adata_tracks = self.model_adata
        model_path = self._model_adata_path()
        if adata_tracks is None and model_path.exists():
            try:
                adata_tracks = sc.read_h5ad(model_path)
            except Exception:
                adata_tracks = None
        if adata_tracks is not None and "sample_name" in adata_tracks.obs.columns:
            names.extend(
                adata_tracks.obs["sample_name"].astype("string").dropna().str.strip().tolist()
            )
        return sorted({str(name).strip() for name in names if str(name).strip() != ""})

    def _refresh_backprojection_samples(self):
        sample_names = self._detect_backprojection_sample_names()
        self.backproj_sample_dd.options = sample_names
        if len(sample_names) == 0:
            self.backproj_sample_dd.value = None
            self.backprojection_status.value = "<i>No samples detected for backprojection.</i>"
            self.btn_backproject.disabled = True
            return
        if self.backproj_sample_dd.value not in sample_names:
            self.backproj_sample_dd.value = sample_names[0]
        self.backprojection_status.value = f"<b>Available samples:</b> {len(sample_names)}"
        self.btn_backproject.disabled = False

    def _rebuild_rename_rows(self):
        self._name_boxes = {}
        if bool(self.use_original_behav3d.value):
            clustered_path = _feature_dtw_clustered_csv_path(self.output_dir, self._current_cell_type())
            if not clustered_path.exists():
                self.rename_rows.children = []
                self.rename_status.value = "<i>Run original BEHAV3D clustering first to rename clusters.</i>"
                self.btn_rename.disabled = True
                return
            try:
                df = pd.read_csv(clustered_path, usecols=["ClusterID"], low_memory=False)
            except Exception as exc:
                self.rename_rows.children = []
                self.rename_status.value = f"<i>Could not load original BEHAV3D clusters: {exc}</i>"
                self.btn_rename.disabled = True
                return
            mapping = _build_identity_cluster_mapping_from_obs(df, cluster_col="ClusterID")
            saved_mapping = _load_feature_dtw_name_mapping(self.output_dir, self._current_cell_type())
            for key in list(mapping.keys()):
                mapping[key] = saved_mapping.get(key, key)
            label = "Original BEHAV3D clusters"
        else:
            try:
                adata_tracks = self.model_adata
                if adata_tracks is None and self._model_adata_path().exists():
                    adata_tracks = sc.read_h5ad(self._model_adata_path())
            except Exception:
                adata_tracks = None
            if adata_tracks is None or "ClusterID" not in adata_tracks.obs.columns:
                self.rename_rows.children = []
                self.rename_status.value = "<i>Run DTAI clustering first to rename clusters.</i>"
                self.btn_rename.disabled = True
                return
            mapping = _build_identity_cluster_mapping_from_obs(adata_tracks.obs, cluster_col="ClusterID")
            label = "DTAI clusters"

        rows = []
        for old_name, current_name in mapping.items():
            txt = widgets.Text(value=str(current_name), layout=widgets.Layout(width="300px"))
            self._name_boxes[str(old_name)] = txt
            rows.append(
                widgets.HBox(
                    [widgets.Label(str(old_name), layout=widgets.Layout(width="120px")), txt],
                    layout=widgets.Layout(align_items="center", gap="8px"),
                )
            )
        self.rename_rows.children = rows
        self.rename_status.value = f"<b>{label}:</b> {len(rows)}"
        self.btn_rename.disabled = len(rows) == 0

    def _on_refresh_rename_clicked(self, _):
        self._rebuild_rename_rows()

    def _on_rename_clicked(self, _):
        self._set_busy(self.btn_rename, self.rename_spinner, busy=True)
        self.out_rename.clear_output()
        with self.out_rename:
            try:
                mapping = {}
                for old_name, txt in self._name_boxes.items():
                    new_name = str(txt.value).strip()
                    mapping[str(old_name)] = new_name if new_name != "" else str(old_name)
                if bool(self.use_original_behav3d.value):
                    self._rename_original_behav3d_clusters(mapping)
                else:
                    self._rename_dtaidistance_clusters(mapping)
                self._rebuild_rename_rows()
                self.steps.selected_index = 2
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_rename, self.rename_spinner, busy=False)

    def _rename_dtaidistance_clusters(self, mapping):
        adata_tracks = self._load_model_adata()
        rename_track_clusters(
            adata=adata_tracks,
            mapping=mapping,
            cluster_col="ClusterID",
            keep_unmapped=True,
        )
        adata_tracks.write(self._model_adata_path(), compression="lzf")
        qc_dir = _resolve_dtaidistance_paths(self.output_dir, self._current_cell_type())[
            "quality_control_outfolder"
        ] / "after_renaming"
        plot_paths = save_dtaidistance_diagnostics(
            adata_tracks,
            output_dir=self.output_dir,
            cell_type=self._current_cell_type(),
            cluster_key="ClusterID",
            outfolder=qc_dir,
            random_state=int(self.random_state.value),
            save_outputs=True,
            verbose=True,
        )
        self.plot_status_html.value = "<b>Renamed QC ready:</b> diagnostics were written after renaming."
        _winfo("trajectory-dtai-widget", f"Renamed clusters and wrote QC: {plot_paths.get('diagnostics_pdf')}")

    def _rename_original_behav3d_clusters(self, mapping):
        ct = self._current_cell_type()
        outdir = _feature_dtw_outdir(self.output_dir, ct)
        mapping_path = _feature_dtw_rename_mapping_path(self.output_dir, ct)
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        with mapping_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(
                {
                    "cell_type": ct,
                    "cluster_id_column": "ClusterID",
                    "cluster_name_column": "ClusterName",
                    "cluster_names": dict(mapping),
                },
                f,
                sort_keys=False,
            )

        updated = []
        for csv_path in _feature_dtw_output_csv_paths(self.output_dir, ct):
            if not csv_path.exists():
                continue
            df = pd.read_csv(csv_path, low_memory=False)
            if "ClusterID" not in df.columns:
                continue
            df["ClusterName"] = df["ClusterID"].astype(str).map(mapping).fillna(df["ClusterID"].astype(str))
            df.to_csv(csv_path, index=False)
            updated.append(csv_path)

        plot_paths = _save_feature_dtw_quality_control(
            self.output_dir,
            ct,
            outfolder=outdir / "quality_control" / "after_renaming",
        )
        self.plot_status_html.value = "<b>Renamed QC ready:</b> original BEHAV3D diagnostics were written after renaming."
        _winfo("trajectory-dtai-widget", f"Saved original BEHAV3D cluster names: {mapping_path}")
        for path in updated:
            _winfo("trajectory-dtai-widget", f"Updated renamed CSV: {path}")
        for path in plot_paths.values():
            _winfo("trajectory-dtai-widget", f"Created renamed QC: {path}")

    def _load_original_exemplar_data(self):
        ct = self._current_cell_type()
        clustered_path = _feature_dtw_clustered_csv_path(self.output_dir, ct)
        if not clustered_path.exists():
            raise FileNotFoundError(f"Original BEHAV3D clustered CSV not found: {clustered_path}")
        state_path = self._state_adata_path(ct)
        if not Path(state_path).exists():
            raise FileNotFoundError(f"Behavioral-state h5ad not found: {state_path}")

        adata_full = sc.read_h5ad(state_path)
        adata_filt = filter_and_truncate_tracks_anndata(
            adata_full,
            groupby_cols=["sample_name", "TrackID"],
            time_col="position_t",
            min_length=int(self.original_trajectory_size.value),
            max_length=int(self.original_trajectory_size.value),
        )
        adata_filt = adata_filt.copy()
        adata_filt.obs["sample_name"] = adata_filt.obs["sample_name"].astype(str)
        adata_filt.obs["TrackID"] = adata_filt.obs["TrackID"].astype(str)
        clusters = pd.read_csv(clustered_path, low_memory=False)
        required = ["sample_name", "TrackID", "ClusterID"]
        missing = [col for col in required if col not in clusters.columns]
        if missing:
            raise ValueError(f"Original BEHAV3D clustered CSV missing required columns: {missing}")
        clusters = (
            clusters[required]
            .dropna(subset=["sample_name", "TrackID", "ClusterID"])
            .assign(
                sample_name=lambda df: df["sample_name"].astype(str),
                TrackID=lambda df: df["TrackID"].astype(str),
                ClusterID=lambda df: df["ClusterID"].astype(str),
            )
            .drop_duplicates(subset=["sample_name", "TrackID"])
        )
        track_obs = (
            adata_filt.obs[["sample_name", "TrackID", "position_t"]]
            .copy()
            .assign(
                sample_name=lambda df: df["sample_name"].astype(str),
                TrackID=lambda df: df["TrackID"].astype(str),
            )
            .groupby(["sample_name", "TrackID"], observed=True, as_index=False)
            .agg(position_t_min=("position_t", "min"), position_t_max=("position_t", "max"))
            .merge(clusters, on=["sample_name", "TrackID"], how="inner")
        )
        if len(track_obs) == 0:
            raise ValueError("No original BEHAV3D clustered tracks matched the behavioral-state h5ad.")
        adata_tracks = ad.AnnData(
            X=np.zeros((len(track_obs), 1), dtype=float),
            obs=track_obs,
            var=pd.DataFrame(index=["feature_dtw_cluster"]),
        )
        return adata_filt, adata_tracks

    def _save_original_exemplar_plots(self):
        adata_filt, adata_tracks = self._load_original_exemplar_data()
        ct = self._current_cell_type()
        outdir = _feature_dtw_outdir(self.output_dir, ct) / "clustering" / "example_tracks"
        outdir.mkdir(parents=True, exist_ok=True)
        n_per_cluster = int(self.n_per_cluster.value)
        num_example_ranks = 5
        chosen_exemplars, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=n_per_cluster,
            sample_key="sample_name",
            track_key="TrackID",
            cluster_key="ClusterID",
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            seed=int(self.random_state.value),
        )
        selection_csv = outdir / "example_track_selection_cluster_ClusterID_state_behavioral_state_original_behav3d.csv"
        chosen_exemplars.to_csv(selection_csv, index=False)

        plot_paths = {"exemplar_selection_csv": str(selection_csv)}
        if bool(self.make_overview_statebars.value):
            fig_exemplar, _, _ = plot_exemplar_tracks_by_cluster(
                adata_filt,
                adata_tracks,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                seed=int(self.random_state.value),
            )
            overview_pdf = outdir / "example_tracks_overview.pdf"
            with PdfPages(overview_pdf) as pdf:
                pdf.savefig(fig_exemplar, bbox_inches="tight", dpi=300)
            plt.close(fig_exemplar)
            statebar_out = save_exemplar_statebar_track_pdf_per_cluster(
                adata_full=adata_filt,
                out_dir=outdir,
                chosen_df=chosen_exemplars,
                adata_tracks=None,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                rows_per_page=6,
                plot_dpi=300,
                seed=int(self.random_state.value),
                cmap_name="tab20",
                layout_mode="both",
                num_example_ranks=num_example_ranks,
            )
            plot_paths["exemplar_tracks_overview_pdf"] = str(overview_pdf)
            plot_paths["exemplar_statebar_track_pdf_by_cluster"] = dict(
                statebar_out.get("pdf_paths_by_cluster", {})
            )

        if bool(self.make_backprojection_pdf.value):
            pdf_out = save_exemplar_statebar_backprojection_pdf(
                adata_full=adata_filt,
                output_dir=self.output_dir,
                cell_type=ct,
                out_dir=outdir / "backprojection",
                chosen_df=chosen_exemplars,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                plot_dpi=220,
                seed=int(self.random_state.value),
                examples_per_cluster=n_per_cluster,
                num_example_ranks=num_example_ranks,
                verbose=True,
            )
            plot_paths["exemplar_backprojection_pdf_by_cluster"] = dict(pdf_out.get("pdf_paths_by_cluster", {}))

        if bool(self.make_backprojection_mp4.value):
            mp4_out = save_exemplar_statebar_backprojection_video_per_cluster(
                adata_full=adata_filt,
                output_dir=self.output_dir,
                cell_type=ct,
                out_dir=outdir / "backprojection",
                chosen_df=chosen_exemplars,
                n_per_cluster=n_per_cluster,
                sample_key="sample_name",
                track_key="TrackID",
                time_key="position_t",
                state_key=FULL_STATE_COL,
                cluster_key="ClusterID",
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                dpi=180,
                seed=int(self.random_state.value),
                examples_per_cluster=n_per_cluster,
                num_example_ranks=num_example_ranks,
                verbose=True,
            )
            plot_paths["exemplar_backprojection_mp4_by_cluster"] = dict(mp4_out.get("video_paths_by_cluster", {}))

        return plot_paths

    def _original_behav3d_features(self, csv_path):
        df_head = pd.read_csv(csv_path, nrows=5)
        preferred = [
            "mean_square_displacement",
            "speed",
            "mean_dead_dye",
            f"{self._current_cell_type()}_contact",
            "organoid_contact",
        ]
        return [col for col in preferred if col in df_head.columns]

    def _on_run_clicked(self, _):
        self._set_busy(self.btn_run, self.run_spinner, busy=True)
        self.out_run.clear_output()
        with self.out_run:
            try:
                trajectory_size = _resolve_optional_int(self.behavioral_trajectory_size.value)
                _winfo(
                    "trajectory-dtai-widget",
                    f"Running one-hot dtaidistance on fixed state column='{FULL_STATE_COL}'",
                )
                self.model_adata = run_categorical_dtaidistance_trajectory_clustering(
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    behavioral_trajectory_size=trajectory_size,
                    min_track_length=trajectory_size,
                    trajectory_trim_mode=str(self.trajectory_trim_mode.value),
                    max_tracks=None,
                    n_clusters=int(self.n_clusters.value),
                    window=None,
                    max_dist=None,
                    max_length_diff=None,
                    penalty=None,
                    psi=None,
                    parallel=bool(self.parallel.value),
                    linkage=str(self.linkage.value),
                    missing_policy="keep",
                    save_distance_matrix=bool(self.save_distance_matrix.value),
                    plot_results=True,
                    plot_exemplars=False,
                    n_per_cluster=int(self.n_per_cluster.value),
                    random_state=int(self.random_state.value),
                    verbose=True,
                )
                self.status_html.value = (
                    f"<b style='color:#080;'>Saved one-hot dtaidistance model:</b> {self._model_adata_path().name} "
                    f"({self.model_adata.n_obs} tracks)"
                )
                self.plot_status_html.value = "<b>Ready for plots:</b> clustering finished."
                self.btn_diagnostics.disabled = False
                self.btn_exemplars.disabled = False
                self._rebuild_rename_rows()
                self._refresh_backprojection_samples()
                self._sync_exemplar_state_controls()
                self.steps.selected_index = 1
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_run, self.run_spinner, busy=False)

    def _on_run_original_clicked(self, _):
        self._set_busy(self.btn_run_original, self.original_spinner, busy=True)
        self.out_original.clear_output()
        with self.out_original:
            try:
                ct = self._current_cell_type()
                csv_path = self._original_track_features_path(ct)
                if not csv_path.exists():
                    raise FileNotFoundError(f"Original BEHAV3D track-features CSV not found: {csv_path}")
                features = self._original_behav3d_features(csv_path)
                if len(features) == 0:
                    raise ValueError(
                        "None of the default original BEHAV3D features were found in the track-features CSV."
                    )
                normalize = [
                    col
                    for col in ["mean_square_displacement", "speed", "mean_dead_dye"]
                    if col in features
                ]
                _winfo(
                    "trajectory-dtai-widget",
                    f"Running original BEHAV3D feature DTW with features={features}",
                )
                run_tcell_analysis(
                    cell_type=ct,
                    output_dir=self.output_dir,
                    df_tracks_path=str(csv_path),
                    columns_to_use=features,
                    columns_to_normalize=normalize,
                    umap_minimal_distance=float(self.original_umap_min_dist.value),
                    umap_n_neighbors=int(self.original_umap_n_neighbors.value),
                    nr_of_clusters=int(self.original_n_clusters.value),
                    plot_results=True,
                    seed=int(self.random_state.value),
                    output_subdir_name="timepoint_feature_dtw",
                    feature_scaling_preset="original_behav3d",
                    min_track_length=int(self.original_trajectory_size.value),
                    max_track_length=int(self.original_trajectory_size.value),
                )
                plot_paths = _save_feature_dtw_quality_control(self.output_dir, ct)
                _winfo("trajectory-dtai-widget", "Original BEHAV3D feature-DTW analysis complete.")
                for path in plot_paths.values():
                    _winfo("trajectory-dtai-widget", f"Created original BEHAV3D QC: {path}")
                self.btn_diagnostics.disabled = False
                self.plot_status_html.value = "<b>Ready for renaming:</b> raw original BEHAV3D QC was written."
                self._rebuild_rename_rows()
                self._sync_exemplar_state_controls()
                self.steps.selected_index = 1
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_run_original, self.original_spinner, busy=False)

    def _on_diagnostics_clicked(self, _):
        self._set_busy(self.btn_diagnostics, self.diagnostics_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if bool(self.use_original_behav3d.value):
                    plot_paths = _save_feature_dtw_quality_control(
                        self.output_dir,
                        self._current_cell_type(),
                    )
                    self.plot_status_html.value = "<b>Diagnostics ready:</b> original BEHAV3D QC was written."
                    for path in plot_paths.values():
                        _winfo("trajectory-dtai-widget", f"Created original BEHAV3D QC: {path}")
                else:
                    adata_tracks = self._load_model_adata()
                    plot_paths = save_dtaidistance_diagnostics(
                        adata_tracks,
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        random_state=int(self.random_state.value),
                        save_outputs=True,
                        verbose=True,
                    )
                    self.plot_status_html.value = "<b>Diagnostics ready:</b> cluster diagnostics were written."
                    _winfo(
                        "trajectory-dtai-widget",
                        f"Created diagnostics plots: {plot_paths.get('diagnostics_pdf')}",
                    )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_diagnostics, self.diagnostics_spinner, busy=False)

    def _on_exemplars_clicked(self, _):
        self._set_busy(self.btn_exemplars, self.exemplar_spinner, busy=True)
        self.out_plots.clear_output()
        with self.out_plots:
            try:
                if not self._any_exemplar_outputs_selected():
                    raise ValueError("Select at least one exemplar output.")
                if not self._has_behavioral_states():
                    raise ValueError(
                        "No states have been defined. Run behavioral state classification before creating exemplar outputs."
                    )
                if bool(self.use_original_behav3d.value):
                    plot_paths = self._save_original_exemplar_plots()
                else:
                    adata_tracks = self._load_model_adata()
                    plot_paths = save_dtaidistance_exemplar_plots(
                        adata_tracks,
                        output_dir=self.output_dir,
                        cell_type=self._current_cell_type(),
                        n_per_cluster=int(self.n_per_cluster.value),
                        make_overview_statebars=bool(self.make_overview_statebars.value),
                        make_backprojection_pdf=bool(self.make_backprojection_pdf.value),
                        make_backprojection_mp4=bool(self.make_backprojection_mp4.value),
                        random_state=int(self.random_state.value),
                        save_outputs=True,
                        verbose=True,
                    )
                self.plot_status_html.value = "<b>Exemplar PDFs ready:</b> exemplar track outputs were written."
                _winfo(
                    "trajectory-dtai-widget",
                    f"Created exemplar overview: {plot_paths.get('exemplar_tracks_overview_pdf')}",
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_exemplars, self.exemplar_spinner, busy=False)
                self._sync_exemplar_state_controls()

    def _on_backproject_clicked(self, _):
        self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=True)
        self.out_backprojection.clear_output()
        with self.out_backprojection:
            try:
                adata_tracks = self._load_model_adata()
                sample_name = self.backproj_sample_dd.value
                if sample_name is None or len(str(sample_name).strip()) == 0:
                    raise ValueError("Please select a sample.")
                adata_full_path = self._state_adata_path()
                if not Path(adata_full_path).exists():
                    raise FileNotFoundError(f"Behavioral-state h5ad not found: {adata_full_path}")

                adata_full = sc.read_h5ad(adata_full_path)
                output_col = "dtai_track_behavioral_cluster"
                manifest = export_track_cluster_backprojection(
                    adata_full=adata_full,
                    adata_tracks=adata_tracks,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    cluster_col="ClusterID",
                    output_col=output_col,
                    sample_name=str(sample_name),
                    verbose=True,
                )
                sample_key = str(sample_name).strip()
                state_img_path = manifest.get("output_paths", {}).get(sample_key, None)
                if state_img_path is None:
                    raise RuntimeError(
                        "Backprojection export finished but no state image was written for sample "
                        f"'{sample_key}'. manifest={manifest}"
                    )
                _winfo(
                    "trajectory-dtai-widget",
                    f"Opening backprojection for sample '{sample_key}' | state_image={state_img_path}",
                )
                show_track_cluster_backprojection(
                    sample_name=sample_key,
                    output_dir=self.output_dir,
                    cell_type=self._current_cell_type(),
                    adata_full=adata_full,
                    adata_tracks=adata_tracks,
                    cluster_col="ClusterID",
                    state_col=FULL_STATE_COL,
                    state_img_path=state_img_path,
                    output_col=output_col,
                    run=True,
                    verbose=True,
                )
            except Exception:
                traceback.print_exc()
            finally:
                self._set_busy(self.btn_backproject, self.backprojection_spinner, busy=False)
                self._refresh_backprojection_samples()
