from pathlib import Path
import shutil
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import silhouette_score

from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.dtw import (
    _add_cluster_medoids,
    _cluster_precomputed_distances,
    _ensure_dtaidistance_umap,
    _relabel_by_cluster_size,
    _validate_distance_matrix,
    compute_dtaidistance_onehot_distance_matrix,
    extract_categorical_track_sequences,
)
from behav3d.analysis.behavior.track.utils import (
    _default_behavioral_states_path,
    _filter_tracks_for_dtaidistance,
    _resolve_dtaidistance_paths,
    _winfo,
    get_dtaidistance_track_trajectories_filename,
)
from behav3d.analysis.behavior.track.bouts import (
    _build_track_feature_adata,
    _resolve_track_classifier_path,
    train_track_classifier,
)
from behav3d.analysis.behavior.utils import _sanitize_filename_token
from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
    ensure_exemplar_coordinate_columns as _ensure_exemplar_coordinate_columns,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    select_exemplar_tracks_by_cluster,
)


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

    coord_enrichment = _ensure_exemplar_coordinate_columns(
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
    clear_outputs=True,
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
    if bool(clear_outputs):
        outfolder = Path(output_dir).expanduser() / "analysis" / str(cell_type) / str(output_subdir_name)
        if outfolder.exists():
            shutil.rmtree(outfolder, ignore_errors=True)
            if bool(verbose):
                _winfo("trajectory-dtai", f"cleared previous outputs: {outfolder}")
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
            _ensure_exemplar_coordinate_columns(
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


def train_dtaidistance_trajectory_classifier(
    output_dir,
    cell_type,
    model_adata,
    *,
    cluster_col="ClusterID",
    classifier_n_estimators=300,
    classifier_min_samples_leaf=2,
    classifier_min_samples_split=2,
    classifier_max_features="sqrt",
    classifier_max_depth=None,
    classifier_n_jobs=-1,
    validation_test_size=0.2,
    validation_stratify=True,
    random_state=123,
    save_classifier=True,
    verbose=True,
):
    """Train a RandomForest classifier from a dtaidistance trajectory model.

    Builds track-level state features (fractions, bigrams, bout stats) from the
    source behavioral-states file recorded in model_adata.uns, transfers the
    ClusterID labels, then delegates to train_track_classifier.
    """
    meta = model_adata.uns.get("dtai_trajectory_clustering", {})
    adata_full_path = meta.get("source_adata_full_path")
    if not adata_full_path:
        raise ValueError(
            "model_adata.uns['dtai_trajectory_clustering']['source_adata_full_path'] is missing. "
            "Re-run dtaidistance clustering to populate it."
        )
    adata_full_path = Path(adata_full_path).expanduser()
    if not adata_full_path.exists():
        raise FileNotFoundError(f"Source behavioral-states file not found: {adata_full_path}")

    trajectory_size = meta.get("behavioral_trajectory_size") or 100

    if bool(verbose):
        _winfo("trajectory-dtai-clf", f"loading behavioral states: {adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)

    if bool(verbose):
        _winfo("trajectory-dtai-clf", f"building track features (trajectory_size={trajectory_size})")
    feature_adata = _build_track_feature_adata(
        adata_full,
        state_col=str(FULL_STATE_COL),
        behavioral_trajectory_size=int(trajectory_size),
    )

    # Transfer ClusterID from dtaidistance model to feature adata by (sample_name, TrackID)
    dtai_obs = model_adata.obs[[cluster_col, "sample_name", "TrackID"]].copy()
    dtai_obs = dtai_obs.reset_index(drop=True)
    dtai_obs = dtai_obs.set_index(["sample_name", "TrackID"])[cluster_col]

    feat_obs = feature_adata.obs.copy()
    feat_obs["_cluster"] = (
        pd.MultiIndex.from_frame(feat_obs[["sample_name", "TrackID"]])
        .map(dtai_obs.to_dict())
    )
    has_label = feat_obs["_cluster"].notna()
    n_dropped = (~has_label).sum()
    if n_dropped > 0 and bool(verbose):
        _winfo("trajectory-dtai-clf", f"dropped {n_dropped} tracks with no matching ClusterID")
    feature_adata = feature_adata[has_label.values].copy()
    feature_adata.obs[cluster_col] = feat_obs.loc[has_label, "_cluster"].astype(str).values

    classifier_path = _resolve_track_classifier_path(
        output_dir, cell_type, output_subdir_name="behavorial_trajectories"
    )

    result = train_track_classifier(
        output_dir=output_dir,
        cell_type=cell_type,
        model_adata=feature_adata,
        cluster_col=cluster_col,
        classifier_n_estimators=int(classifier_n_estimators),
        classifier_min_samples_leaf=int(classifier_min_samples_leaf),
        classifier_min_samples_split=int(classifier_min_samples_split),
        classifier_max_features=classifier_max_features,
        classifier_max_depth=None if classifier_max_depth is None else int(classifier_max_depth),
        classifier_n_jobs=int(classifier_n_jobs),
        validation_test_size=float(validation_test_size),
        validation_stratify=bool(validation_stratify),
        random_state=int(random_state),
        save_classifier=bool(save_classifier),
        classifier_path=classifier_path,
        output_subdir_name="behavorial_trajectories",
        verbose=bool(verbose),
    )
    return result
