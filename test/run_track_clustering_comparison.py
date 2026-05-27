"""Run track-clustering methods side by side and compare exemplar state tracks.

Edit the parameters below, then run this file in the BEHAV3D environment.
The script assumes HMM behavioral states have already been assigned.
"""

from pathlib import Path
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages

from behav3d.analysis.clustering.state.classification import FULL_STATE_COL
from behav3d.analysis.clustering.track.classification import run_state_based_analysis
from behav3d.analysis.clustering.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.analysis.tcell_analysis import run_tcell_analysis

TEST_DIR = Path(__file__).resolve().parent
if str(TEST_DIR) not in sys.path:
    sys.path.insert(0, str(TEST_DIR))

from trajectory_dtai_classification import run_categorical_dtaidistance_trajectory_clustering


# ---------------------------------------------------------------------------
# Edit these parameters before running.
# ---------------------------------------------------------------------------
OUTPUT_DIR = Path(".")
CELL_TYPE = "tcell"
STATE_COL = FULL_STATE_COL

TRAJECTORY_SIZE = 100
N_CLUSTERS = 6
N_PER_CLUSTER = 5
RANDOM_STATE = 123
MAX_TRACKS = None  # e.g. 20 for a fast smoke test; None runs all available tracks.

UMAP_MIN_DIST = 0.1
UMAP_N_NEIGHBORS = 15

TCELL_DTW_FEATURES = [
    "mean_square_displacement",
    "speed",
    "mean_dead_dye",
    "tcell_contact",
    "organoid_contact",
]
TCELL_DTW_NORMALIZE_FEATURES = [
    "mean_square_displacement",
    "speed",
    "mean_dead_dye",
]

RUN_TCELL_FEATURE_DTW = True
RUN_ONEHOT_DTAI = True
RUN_STATE_FEATURE_CLUSTERING = True


def _method_subdir(name):
    return f"track_clustering_comparison/{name}"


def _analysis_dir():
    return Path(OUTPUT_DIR, "analysis", CELL_TYPE)


def _comparison_root():
    path = _analysis_dir() / "track_clustering_comparison"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _state_adata_path():
    return _analysis_dir() / "behavioral_states" / f"BEHAV3D_{CELL_TYPE}_behavioral_states.h5ad"


def _track_features_path():
    feature_dir = _analysis_dir() / "track_features"
    filtered = feature_dir / f"BEHAV3D_{CELL_TYPE}_combined_track_features_filtered.csv"
    if filtered.exists():
        return filtered
    return feature_dir / f"BEHAV3D_{CELL_TYPE}_combined_track_features.csv"


def _track_features_summarized_path():
    return (
        _analysis_dir()
        / "track_features"
        / f"BEHAV3D_{CELL_TYPE}_combined_track_features_summarized.csv"
    )


def _key_frame(df, sample_col="sample_name", track_col="TrackID"):
    out = df[[sample_col, track_col]].copy()
    out[sample_col] = out[sample_col].astype(str)
    out[track_col] = out[track_col].astype(str)
    return out.drop_duplicates()


def _normalize_obs_keys(adata):
    adata = adata.copy()
    for col in ("sample_name", "TrackID"):
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype(str)
    return adata


def _sample_common_tracks(adata_full):
    if MAX_TRACKS is None or int(MAX_TRACKS) <= 0:
        return None

    hmm_keys = _key_frame(adata_full.obs)
    feature_path = _track_features_path()
    if feature_path.exists():
        feature_keys = _key_frame(pd.read_csv(feature_path, usecols=["sample_name", "TrackID"]))
        common = hmm_keys.merge(feature_keys, on=["sample_name", "TrackID"], how="inner")
    else:
        common = hmm_keys

    if len(common) == 0:
        raise ValueError("No common tracks found for sampling.")

    n_keep = min(int(MAX_TRACKS), len(common))
    sampled = common.sample(n=n_keep, random_state=int(RANDOM_STATE)).reset_index(drop=True)
    print(f"Sampled {len(sampled)} common tracks for comparison.")
    return sampled


def _filter_df_to_keys(df, keys):
    if keys is None:
        return df
    work = df.copy()
    work["sample_name"] = work["sample_name"].astype(str)
    work["TrackID"] = work["TrackID"].astype(str)
    return work.merge(keys, on=["sample_name", "TrackID"], how="inner")


def _prepare_inputs(adata_full):
    """Return paths to full/sampled state h5ad and feature CSVs."""
    root = _comparison_root()
    keys = _sample_common_tracks(adata_full)
    if keys is None:
        return {
            "adata_path": _state_adata_path(),
            "features_path": _track_features_path(),
            "features_summarized_path": _track_features_summarized_path(),
            "keys": None,
        }

    key_tuples = set(map(tuple, keys[["sample_name", "TrackID"]].to_numpy()))
    mask = [
        (str(sample), str(track)) in key_tuples
        for sample, track in adata_full.obs[["sample_name", "TrackID"]].itertuples(index=False, name=None)
    ]
    sampled_adata = adata_full[np.asarray(mask, dtype=bool)].copy()
    sampled_adata_path = root / "sampled_behavioral_states.h5ad"
    sampled_adata.write(sampled_adata_path, compression="gzip")

    sampled_features_path = root / "sampled_track_features.csv"
    sampled_summary_path = root / "sampled_track_features_summarized.csv"
    _filter_df_to_keys(pd.read_csv(_track_features_path(), low_memory=False), keys).to_csv(
        sampled_features_path,
        index=False,
    )
    _filter_df_to_keys(pd.read_csv(_track_features_summarized_path(), low_memory=False), keys).to_csv(
        sampled_summary_path,
        index=False,
    )

    return {
        "adata_path": sampled_adata_path,
        "features_path": sampled_features_path,
        "features_summarized_path": sampled_summary_path,
        "keys": keys,
    }


def _filtered_state_adata(adata_full):
    adata = filter_and_truncate_tracks_anndata(
        adata_full,
        groupby_cols=["sample_name", "TrackID"],
        time_col="position_t",
        min_length=int(TRAJECTORY_SIZE),
        max_length=int(TRAJECTORY_SIZE),
    )
    return _normalize_obs_keys(adata)


def _track_windows_from_states(adata_states):
    obs = adata_states.obs[["sample_name", "TrackID", "position_t"]].copy()
    obs["sample_name"] = obs["sample_name"].astype(str)
    obs["TrackID"] = obs["TrackID"].astype(str)
    obs["position_t"] = pd.to_numeric(obs["position_t"], errors="coerce")
    windows = (
        obs.dropna(subset=["position_t"])
        .groupby(["sample_name", "TrackID"], observed=True, as_index=False)
        .agg(
            position_t_min=("position_t", "min"),
            position_t_max=("position_t", "max"),
            n_timepoints=("position_t", "size"),
        )
    )
    return windows


def _adata_from_tcell_clusters(cluster_csv, adata_states):
    df = pd.read_csv(cluster_csv, low_memory=False)
    required = {"sample_name", "TrackID", "ClusterID"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"tcell_analysis cluster CSV missing columns: {missing}")

    df = df[["sample_name", "TrackID", "ClusterID"]].copy()
    df["sample_name"] = df["sample_name"].astype(str)
    df["TrackID"] = df["TrackID"].astype(str)
    df["ClusterID"] = df["ClusterID"].astype(str)
    df = df.drop_duplicates(subset=["sample_name", "TrackID"])
    df = df.merge(_track_windows_from_states(adata_states), on=["sample_name", "TrackID"], how="inner")
    if len(df) == 0:
        raise ValueError("No tcell_analysis clustered tracks overlapped the HMM state windows.")
    df["ClusterID"] = pd.Categorical(df["ClusterID"])
    df.index = [f"{s}__{t}" for s, t in df[["sample_name", "TrackID"]].itertuples(index=False, name=None)]
    return ad.AnnData(
        X=np.zeros((len(df), 1), dtype=float),
        obs=df,
        var=pd.DataFrame(index=["tcell_feature_dtw_cluster"]),
    )


def _run_tcell_feature_dtw(inputs):
    out_subdir = _method_subdir("timepoint_feature_dtw")
    run_tcell_analysis(
        output_dir=OUTPUT_DIR,
        df_tracks_path=inputs["features_path"],
        df_tracks_summarized_path=inputs["features_summarized_path"],
        cell_type=CELL_TYPE,
        columns_to_use=TCELL_DTW_FEATURES,
        columns_to_normalize=TCELL_DTW_NORMALIZE_FEATURES,
        umap_minimal_distance=float(UMAP_MIN_DIST),
        umap_n_neighbors=int(UMAP_N_NEIGHBORS),
        nr_of_clusters=int(N_CLUSTERS),
        plot_results=True,
        seed=int(RANDOM_STATE),
        output_subdir_name=out_subdir,
    )
    cluster_csv = _analysis_dir() / out_subdir / f"BEHAV3D_{CELL_TYPE}_UMAP_clusters.csv"
    return {
        "name": "timepoint_feature_dtw",
        "label": "Timepoint feature DTW (HMM state trajectories shown)",
        "cluster_csv": cluster_csv,
        "outdir": _analysis_dir() / out_subdir,
    }


def _run_onehot_dtaidistance(inputs):
    out_subdir = _method_subdir("onehot_dtaidistance")
    adata_tracks = run_categorical_dtaidistance_trajectory_clustering(
        output_dir=OUTPUT_DIR,
        cell_type=CELL_TYPE,
        adata_full_path=inputs["adata_path"],
        behavioral_trajectory_size=int(TRAJECTORY_SIZE),
        min_track_length=int(TRAJECTORY_SIZE),
        max_tracks=None,
        n_clusters=int(N_CLUSTERS),
        parallel=True,
        linkage="average",
        missing_policy="keep",
        save_outputs=True,
        save_distance_matrix=False,
        plot_results=True,
        plot_exemplars=False,
        n_per_cluster=int(N_PER_CLUSTER),
        output_subdir_name=out_subdir,
        random_state=int(RANDOM_STATE),
        verbose=True,
    )
    return {
        "name": "onehot_dtaidistance",
        "label": "One-hot DTAIDistance",
        "adata": _normalize_obs_keys(adata_tracks),
        "outdir": _analysis_dir() / out_subdir,
    }


def _run_state_feature_clustering(inputs):
    out_subdir = _method_subdir("state_feature_clustering")
    adata_tracks = run_state_based_analysis(
        output_dir=OUTPUT_DIR,
        cell_type=CELL_TYPE,
        adata_full_path=inputs["adata_path"],
        state_col=STATE_COL,
        behavioral_trajectory_size=int(TRAJECTORY_SIZE),
        n_neighbors=30,
        leiden_resolution=0.2,
        cluster_key="ClusterID",
        umap_min_dist=0.1,
        plot_results=True,
        plot_exemplars=False,
        plot_exemplar_backprojection_videos=False,
        plot_exemplar_backprojection_pdfs=False,
        n_per_cluster=int(N_PER_CLUSTER),
        save_outputs=True,
        output_subdir_name=out_subdir,
        random_state=int(RANDOM_STATE),
        verbose=True,
    )
    return {
        "name": "state_feature_clustering",
        "label": "State features: fractions, bouts, transitions, n-grams",
        "adata": _normalize_obs_keys(adata_tracks),
        "outdir": _analysis_dir() / out_subdir,
    }


def _write_title_page(pdf, title, lines=None):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.text(0.03, 0.92, title, fontsize=20, fontweight="bold", va="top")
    y = 0.82
    for line in lines or []:
        ax.text(0.04, y, str(line), fontsize=11, va="top")
        y -= 0.06
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _add_method_to_pdf(pdf, method, adata_states):
    if "adata" in method:
        adata_tracks = method["adata"]
    else:
        adata_tracks = _adata_from_tcell_clusters(method["cluster_csv"], adata_states)

    _write_title_page(
        pdf,
        method["label"],
        [
            f"Output folder: {method['outdir']}",
            f"Tracks: {adata_tracks.n_obs}",
            f"Displayed trajectory state column: {STATE_COL}",
        ],
    )
    fig, _, _ = plot_exemplar_tracks_by_cluster(
        adata_full=adata_states,
        adata_tracks=adata_tracks,
        n_per_cluster=int(N_PER_CLUSTER),
        sample_key="sample_name",
        track_key="TrackID",
        time_key="position_t",
        state_key=STATE_COL,
        cluster_key="ClusterID",
        tmin_key="position_t_min",
        tmax_key="position_t_max",
        seed=int(RANDOM_STATE),
    )
    pdf.savefig(fig, bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    state_path = _state_adata_path()
    if not state_path.exists():
        raise FileNotFoundError(f"Missing HMM behavioral-state h5ad: {state_path}")
    if not _track_features_path().exists() and RUN_TCELL_FEATURE_DTW:
        raise FileNotFoundError(f"Missing track feature CSV: {_track_features_path()}")
    if not _track_features_summarized_path().exists() and RUN_TCELL_FEATURE_DTW:
        raise FileNotFoundError(f"Missing summarized track feature CSV: {_track_features_summarized_path()}")

    print(f"Loading HMM behavioral states: {state_path}")
    adata_full = _normalize_obs_keys(sc.read_h5ad(state_path))
    if STATE_COL not in adata_full.obs.columns:
        raise ValueError(f"HMM behavioral-state h5ad missing '{STATE_COL}'.")

    inputs = _prepare_inputs(adata_full)
    adata_for_plots = _filtered_state_adata(sc.read_h5ad(inputs["adata_path"]))

    methods = []
    if RUN_TCELL_FEATURE_DTW:
        methods.append(_run_tcell_feature_dtw(inputs))
    if RUN_ONEHOT_DTAI:
        methods.append(_run_onehot_dtaidistance(inputs))
    if RUN_STATE_FEATURE_CLUSTERING:
        methods.append(_run_state_feature_clustering(inputs))

    combined_pdf = _comparison_root() / "track_clustering_example_clusters.pdf"
    with PdfPages(combined_pdf) as pdf:
        _write_title_page(
            pdf,
            "Track Clustering Comparison",
            [
                f"Output dir: {Path(OUTPUT_DIR).resolve()}",
                f"Cell type: {CELL_TYPE}",
                f"Trajectory size: {TRAJECTORY_SIZE}",
                f"N clusters: {N_CLUSTERS}",
                f"Examples per cluster: {N_PER_CLUSTER}",
                f"Max tracks: {MAX_TRACKS}",
                f"Feature DTW selected features: {', '.join(TCELL_DTW_FEATURES)}",
            ],
        )
        for method in methods:
            _add_method_to_pdf(pdf, method, adata_for_plots)

    print("\n### Track clustering comparison complete")
    print(f"Combined exemplar PDF: {combined_pdf}")
    for method in methods:
        print(f"- {method['label']}: {method['outdir']}")


if __name__ == "__main__":
    main()
