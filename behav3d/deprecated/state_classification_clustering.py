import warnings
import math
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from behav3d.analysis.behavior.state.legacy_clustering import *
from behav3d.core.utils import rmtree_ignore_missing

def rename_intrinsic_behavioral_clusters(
    adata,
    mapping,
    binary_cols_to_merge=None,
    keep_unmapped=True,
):
    """
    Relabel intrinsic_behavioral_cluster and rebuild full_behavioral_cluster from binary groups.
    """
    if "intrinsic_behavioral_cluster" not in adata.obs.columns:
        raise ValueError("Missing 'intrinsic_behavioral_cluster' in adata.obs.")

    relabel_cluster_ids(
        adata=adata,
        mapping=mapping,
        cluster_key="intrinsic_behavioral_cluster",
        new_key="intrinsic_behavioral_cluster",
        keep_unmapped=keep_unmapped,
        overwrite_original=True,
    )

    if binary_cols_to_merge is None:
        binary_cols_to_merge = []
        if isinstance(adata.uns, dict):
            binary_cols_to_merge = list(
                adata.uns.get("clustering", {}).get("binary_cols_to_merge", [])
            )
    binary_group_constraints = None
    enforce_binary_group_constraints = False
    if isinstance(adata.uns, dict):
        clustering_meta = adata.uns.get("clustering", {})
        if isinstance(clustering_meta, dict):
            binary_group_constraints = clustering_meta.get("binary_group_constraints", None)
            enforce_binary_group_constraints = (
                isinstance(binary_group_constraints, dict)
                and ("forbidden_binary_combinations" in binary_group_constraints)
            )

    return _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=enforce_binary_group_constraints,
    )


def subsample_with_temporal_spacing(
    data,
    id_cols=None,
    time_col="position_t",
    min_spacing=5,
    min_time_from_track_start=None,
    max_samples=None,
    random_state=123,
):
    """Subsample rows with temporal spacing per track for DataFrame or AnnData input."""
    if id_cols is None:
        id_cols = ["sample_name", "TrackID"]

    rng = np.random.RandomState(random_state)

    def _select_track_row_indices(track_df):
        n_rows = len(track_df)
        if n_rows == 0:
            return []

        eligible_positions = np.arange(n_rows, dtype=int)
        if min_time_from_track_start is not None:
            time_values = pd.to_numeric(track_df[time_col], errors="coerce").to_numpy(dtype=float, copy=False)
            finite_mask = np.isfinite(time_values)
            if not np.any(finite_mask):
                return []
            track_start_t = float(np.nanmin(time_values[finite_mask]))
            min_eligible_t = track_start_t + float(min_time_from_track_start)
            eligible_positions = np.flatnonzero(finite_mask & (time_values >= min_eligible_t))

        n_eligible = int(len(eligible_positions))
        if n_eligible == 0:
            return []

        start_offset = int(rng.randint(0, min(int(min_spacing), n_eligible)))
        return eligible_positions[start_offset:: int(min_spacing)].tolist()

    if hasattr(data, "obs") and hasattr(data, "n_obs"):
        obs = data.obs
        missing_cols = [c for c in list(id_cols) + [time_col] if c not in obs.columns]
        if len(missing_cols) > 0:
            raise ValueError(
                "AnnData input is missing required columns for temporal subsampling: "
                f"{missing_cols}"
            )

        obs_work = obs[list(id_cols) + [time_col]].copy()
        obs_work["_obs_idx"] = np.arange(data.n_obs, dtype=int)
        selected_obs_idx = []

        for _, track_obs in obs_work.groupby(id_cols, observed=False):
            track_obs = track_obs.sort_values(time_col).reset_index(drop=True)
            local_selected = _select_track_row_indices(track_obs)
            if len(local_selected) == 0:
                continue
            selected_obs_idx.extend(track_obs.iloc[local_selected]["_obs_idx"].tolist())

        if len(selected_obs_idx) == 0:
            subset_empty = data[:0].copy()
            subset_empty.uns = copy.deepcopy(dict(getattr(data, "uns", {})))
            return subset_empty

        selected_obs_idx = np.asarray(selected_obs_idx, dtype=int)
        if max_samples is not None and len(selected_obs_idx) > int(max_samples):
            selected_obs_idx = rng.choice(selected_obs_idx, size=int(max_samples), replace=False)
        subset = data[selected_obs_idx].copy()
        subset.uns = copy.deepcopy(dict(getattr(data, "uns", {})))
        return subset

    subsampled_rows = []
    for _, track_df in data.groupby(id_cols):
        track_df = track_df.sort_values(time_col).reset_index(drop=True)
        selected_indices = _select_track_row_indices(track_df)
        if len(selected_indices) == 0:
            continue
        subsampled_rows.append(track_df.iloc[selected_indices])

    if len(subsampled_rows) == 0:
        return data.iloc[:0].copy()

    df_subsampled = pd.concat(subsampled_rows, ignore_index=True)
    if max_samples is not None and len(df_subsampled) > int(max_samples):
        df_subsampled = df_subsampled.sample(n=int(max_samples), random_state=random_state).reset_index(drop=True)
    return df_subsampled


def _resolve_processing_outdir(base_outdir):
    """Return/create processing output root under a stage root folder."""
    if base_outdir is None:
        raise ValueError("base_outdir is required to resolve processing output.")
    processing_outdir = Path(base_outdir) / "processing"
    processing_outdir.mkdir(parents=True, exist_ok=True)
    return processing_outdir


def _resolve_processing_subdir(base_outdir, *parts):
    """Return/create a subfolder under processing output root."""
    processing_outdir = _resolve_processing_outdir(base_outdir)
    subdir = processing_outdir.joinpath(*parts)
    subdir.mkdir(parents=True, exist_ok=True)
    return subdir


def _resolve_state_clustering_outdir(base_outdir):
    """Return/create clustering diagnostics/proportion folder under processing root."""
    if base_outdir is None:
        raise ValueError("base_outdir is required to resolve state clustering output.")
    state_clustering_outdir = Path(base_outdir) / "state_clustering"
    state_clustering_outdir.mkdir(parents=True, exist_ok=True)
    return state_clustering_outdir


def _rmtree_ignore_missing(path):
    """Remove a directory tree while tolerating concurrent missing-file races."""
    path = Path(path)
    if not path.exists():
        return

    rmtree_ignore_missing(path)


def _select_exemplar_windows_by_cluster(
    adata_windows,
    *,
    cluster_key="intrinsic_behavioral_cluster",
    n_per_cluster=5,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    tmin_key="window_start_position_t",
    tmax_key="window_end_position_t",
    seed=0,
):
    required_cols = [sample_key, track_key, cluster_key, time_key, tmin_key, tmax_key]
    missing = [c for c in required_cols if c not in adata_windows.obs.columns]
    if len(missing) > 0:
        raise ValueError(f"adata_windows.obs missing required columns for exemplar selection: {missing}")

    rng = np.random.default_rng(seed)
    window_df = adata_windows.obs[required_cols].copy()
    window_df = window_df.dropna(subset=[cluster_key, tmin_key, tmax_key]).copy()
    if len(window_df) == 0:
        raise ValueError("No windows available after dropping rows missing cluster or window bounds.")

    cluster_total_counts = window_df.groupby(cluster_key, observed=True).size().to_dict()
    chosen_parts = []

    for _, df_cl in window_df.groupby(cluster_key, sort=False, observed=True):
        if len(df_cl) == 0:
            continue

        df_cl = df_cl.copy()
        df_cl["_draw"] = rng.random(len(df_cl))
        df_cl = df_cl.sort_values(
            ["_draw", str(time_key), str(tmin_key), str(tmax_key)],
            kind="stable",
        )
        df_cl = df_cl.drop_duplicates(
            subset=[sample_key, track_key],
            keep="first",
        ).copy()

        k = min(int(n_per_cluster), len(df_cl))
        if k <= 0:
            continue
        idx = rng.choice(len(df_cl), size=k, replace=False)
        chosen_parts.append(df_cl.iloc[idx].copy())

    if len(chosen_parts) == 0:
        raise ValueError("No exemplar windows were selected.")

    chosen = pd.concat(chosen_parts, axis=0, ignore_index=True)
    chosen = chosen.sort_values(
        [str(cluster_key), str(sample_key), str(track_key), str(tmin_key), str(tmax_key), str(time_key)],
        kind="stable",
    ).reset_index(drop=True)
    chosen["example_rank"] = (
        chosen.groupby(cluster_key, observed=True).cumcount() + 1
    ).astype(int)
    return chosen, cluster_total_counts


def _render_state_cluster_exemplar_videos(
    model_adata,
    *,
    state_paths,
    cell_type,
    state_clustering_outdir,
    cluster_key="intrinsic_behavioral_cluster",
    n_per_cluster=5,
    random_state=123,
    video_fps=6,
    video_dpi=180,
    video_margin=20,
    video_pmin=0.0,
    video_pmax=99.99,
    video_track_color="#63ff33",
    video_show_segments=True,
    video_segment_style="outline",
    video_segment_color="#ffffff",
    verbose=True,
):
    result = {
        "enabled": True,
        "cluster_key": str(cluster_key),
        "selection_csv": None,
        "video_paths_by_cluster": {},
        "global_xy_shape": None,
        "cluster_total_counts": {},
        "segment_outline_errors": {},
        "error": None,
        "render_config": {
            "stage": "after_clustering",
            "cluster_key": str(cluster_key),
            "state_key": str(cluster_key),
            "n_per_cluster": int(n_per_cluster),
            "selection_random_state": int(random_state),
            "layout_mode": "per_cluster",
            "video_fps": int(video_fps),
            "video_dpi": int(video_dpi),
            "video_margin": int(video_margin),
            "video_pmin": float(video_pmin),
            "video_pmax": float(video_pmax),
            "video_track_color": str(video_track_color),
            "video_show_segments": bool(video_show_segments),
            "video_segment_style": str(video_segment_style),
            "video_segment_color": str(video_segment_color),
            "show_state_legend": False,
            "positions_csv_path": None,
        },
    }

    if str(cluster_key) not in model_adata.obs.columns:
        result["error"] = f"Missing clustering exemplar column '{cluster_key}' in model_adata.obs."
        return result

    exemplar_started = _vstart(
        verbose,
        "state-clustering",
        f"render exemplar videos | cluster_key={cluster_key}",
    )
    exemplar_root = (
        Path(state_clustering_outdir)
        / "example_windows"
        / _sanitize_filename_token(cluster_key, fallback="cluster")
    )
    _rmtree_ignore_missing(exemplar_root)
    exemplar_root.mkdir(parents=True, exist_ok=True)

    try:
        positions_csv_path = _resolve_positions_csv_path(state_paths=state_paths)
        result["render_config"]["positions_csv_path"] = str(positions_csv_path)
        df_positions = pd.read_csv(positions_csv_path, low_memory=False)
        df_positions = df_positions.sort_values(["sample_name", "TrackID", "position_t"])
        required_pos_cols = [
            "sample_name",
            "TrackID",
            "position_t",
            "pixel_position_x",
            "pixel_position_y",
            "pixel_position_z",
        ]
        missing_pos_cols = [c for c in required_pos_cols if c not in df_positions.columns]
        if len(missing_pos_cols) > 0:
            raise ValueError(
                "State exemplar videos require pixel-position columns in the positions CSV. "
                f"Missing columns: {missing_pos_cols}. csv='{positions_csv_path}'"
            )

        chosen_exemplars, cluster_total_counts = _select_exemplar_windows_by_cluster(
            adata_windows=model_adata,
            n_per_cluster=int(n_per_cluster),
            cluster_key=str(cluster_key),
            seed=int(random_state),
        )
        result["cluster_total_counts"] = {
            str(k): int(v) for k, v in dict(cluster_total_counts).items()
        }
        chosen_exemplars["cluster_total_windows"] = chosen_exemplars[str(cluster_key)].astype(str).map(
            result["cluster_total_counts"]
        )

        selection_csv = exemplar_root / (
            f"example_window_selection_cluster_{_sanitize_filename_token(cluster_key, fallback='cluster')}.csv"
        )
        chosen_exemplars.to_csv(selection_csv, index=False)
        result["selection_csv"] = str(selection_csv)
        _vsave(verbose, "state-clustering", "example window selection CSV", selection_csv)

        backprojection_video_out = save_selected_fulltrack_cluster_videos(
            chosen_df=chosen_exemplars,
            df_positions=df_positions,
            output_folder=state_paths.output_dir,
            out_dir=exemplar_root,
            cluster_key=str(cluster_key),
            sample_key="sample_name",
            track_key="TrackID",
            tmin_key="window_start_position_t",
            tmax_key="window_end_position_t",
            fps=int(video_fps),
            dpi=int(video_dpi),
            margin=int(video_margin),
            pmin=float(video_pmin),
            pmax=float(video_pmax),
            track_color=str(video_track_color),
            show_segment_outlines=bool(video_show_segments),
            segment_style=str(video_segment_style),
            segment_color=str(video_segment_color),
            cell_type=cell_type,
            verbose=bool(verbose),
        )
        result["video_paths_by_cluster"] = dict(
            backprojection_video_out.get("video_paths_by_cluster", {})
        )
        result["segment_outline_errors"] = dict(
            backprojection_video_out.get("segment_outline_errors", {})
        )
    except Exception as exc:
        result["error"] = str(exc)
        _vinfo(verbose, "state-clustering", f"State exemplar videos skipped: {exc}")

    _vdone(verbose, "state-clustering", "render exemplar videos", exemplar_started)
    return result


def plot_clustering_diagnostics_pdf(
    adata,
    cluster_col,
    feature_cols,
    pdf_path,
    title="Clustering diagnostics",
):

    if adata is None or adata.n_obs == 0:
        return
    if cluster_col not in adata.obs.columns:
        return

    ad = _ensure_umap(adata, n_neighbors=30, random_state=123)

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(pdf_path) as pdf:
        fig, (ax_umap, ax_corr) = plt.subplots(
            1, 2, figsize=A4_LANDSCAPE, gridspec_kw={"width_ratios": [1.25, 0.75]}
        )
        sc.pl.umap(
            ad,
            color=cluster_col,
            ax=ax_umap,
            show=False,
            title="UMAP",
            legend_loc="on data",
            legend_fontsize=7,
        )
        ax_umap.set_aspect("equal", adjustable="box")

        # Correlation matrix on same page using explicit axis.  
        sc.pl.correlation_matrix(
            ad,
            groupby=cluster_col,
            show_correlation_numbers=True,
            ax=ax_corr,
            show=False,
        )
        ax_corr.set_title("Cluster correlations", fontsize=10)

        fig.suptitle(f"{title} | UMAP + Correlation", fontsize=12, fontweight="bold", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        umap_bbox = ax_umap.get_position().frozen()

        # Keep UMAP unchanged; only move/shrink correlation panel.
        bbox = ax_corr.get_position()
        extra_gap = 0.03
        new_x = min(0.98 - bbox.width, bbox.x0 + extra_gap)
        height_shrink = 0.5
        new_h = bbox.height * height_shrink
        ax_corr.set_position([new_x, bbox.y0 + (bbox.height - new_h) / 2, bbox.width, new_h])
        ax_umap.set_position(umap_bbox)

        _save_pdf_page_a4(pdf, fig, orientation="landscape")
        plt.close(fig)

        use_dendrogram = ad.obs[cluster_col].astype(str).nunique() > 1
        dendrogram_key = f"dendrogram_{cluster_col}"
        dendrogram_arg = False
        if use_dendrogram:
            try:
                if dendrogram_key not in ad.uns:
                    sc.tl.dendrogram(ad, groupby=cluster_col, key_added=dendrogram_key)
                if dendrogram_key in ad.uns:
                    dendrogram_arg = dendrogram_key
            except Exception:
                dendrogram_arg = False

        sc.pl.heatmap(
            ad,
            var_names=list(feature_cols),
            groupby=cluster_col,
            standard_scale="var",
            figsize=A4_LANDSCAPE,
            swap_axes=True,
            dendrogram=dendrogram_arg,
            show_gene_labels=True,
            show=False,
        )
  
        fig_hm = plt.gcf()
        fig_hm.suptitle(f"{title} | Heatmap", y=0.995)
        _save_pdf_page_a4(pdf, fig_hm, orientation="landscape")
        plt.close(fig_hm)


def plot_feature_distribution_pdf(
    adata,
    feature_cols,
    pdf_path,
    preprocessing_params=None,
    title="Feature distributions",
    features_per_page=6,
):
    """Write per-feature raw and z-scaled violin distributions to a PDF."""
    if adata is None or adata.n_obs == 0:
        return

    valid_features = [str(c) for c in list(feature_cols or []) if str(c) in adata.var_names]
    if len(valid_features) == 0:
        return

    X_scaled = _to_numpy_2d(adata[:, valid_features].X).astype(float, copy=False)
    X_raw = X_scaled
    raw_title = "Raw"
    preprocessing_meta = preprocessing_params or getattr(adata, "uns", {}).get("preprocessing", {})
    scaler_meta = _coerce_scaler_params(preprocessing_meta)
    if scaler_meta is not None:
        mean = np.asarray(scaler_meta["mean"], dtype=float)
        scale = np.asarray(scaler_meta["scale"], dtype=float)
        if mean.shape[0] == len(valid_features) and scale.shape[0] == len(valid_features):
            safe_scale = scale.copy()
            safe_scale[safe_scale == 0] = 1.0
            X_raw = (X_scaled * safe_scale[None, :]) + mean[None, :]
            X_raw = _invert_log_scaling_in_continuous_matrix(
                X_raw,
                preprocessing_params=preprocessing_meta,
                feature_cols=valid_features,
                inplace=True,
                strict=False,
            )
        else:
            raw_title = "Raw (scaler mismatch; showing scaled)"
    else:
        log_meta = _coerce_log_scaling_params(preprocessing_meta)
        matched_log_indices = []
        if log_meta is not None:
            matched_log_indices = _log_scaling_feature_indices(
                valid_features,
                log_meta.get("resolved_feature_cols", []),
                strict=False,
            )
        if len(matched_log_indices) > 0:
            X_raw = _invert_log_scaling_in_continuous_matrix(
                X_raw,
                preprocessing_params=preprocessing_meta,
                feature_cols=valid_features,
                inplace=True,
                strict=False,
            )
        else:
            raw_title = "Raw (scaler unavailable; showing current values)"

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    features_per_page = max(int(features_per_page), 1)
    n_features = len(valid_features)
    with PdfPages(pdf_path) as pdf:
        for start in range(0, n_features, features_per_page):
            page_features = valid_features[start : start + features_per_page]
            nrows = len(page_features)
            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=2,
                figsize=(A4_LANDSCAPE[0], max(2.2 * nrows, 4.5)),
                squeeze=False,
            )
            for row_idx, feature_name in enumerate(page_features):
                idx = start + row_idx
                raw_vals = pd.to_numeric(pd.Series(X_raw[:, idx]), errors="coerce").dropna().to_numpy(dtype=float)
                scaled_vals = pd.to_numeric(pd.Series(X_scaled[:, idx]), errors="coerce").dropna().to_numpy(dtype=float)
                feature_title = str(feature_name)

                for ax, vals, col_title, color in [
                    (axes[row_idx, 0], raw_vals, raw_title, "#4C7E9E"),
                    (axes[row_idx, 1], scaled_vals, "Z-scaled", "#D18951"),
                ]:
                    if vals.size > 0:
                        vp = ax.violinplot(vals, positions=[1], showmeans=True, showextrema=False, widths=0.8)
                        for body in vp["bodies"]:
                            body.set_facecolor(color)
                            body.set_edgecolor(color)
                            body.set_alpha(0.55)
                        if "cmeans" in vp:
                            vp["cmeans"].set_color("#222222")
                            vp["cmeans"].set_linewidth(1.0)
                        x_jitter = np.full(vals.size, 1.0)
                        if vals.size > 1:
                            x_jitter = 1.0 + np.linspace(-0.06, 0.06, vals.size)
                        ax.scatter(x_jitter, vals, s=6, alpha=0.18, color="#222222", linewidths=0)
                    else:
                        ax.text(0.5, 0.5, "No finite values", ha="center", va="center", transform=ax.transAxes)
                    ax.set_xlim(0.5, 1.5)
                    ax.set_xticks([1])
                    ax.set_xticklabels([feature_title], rotation=18, ha="right", fontsize=8)
                    ax.set_title(col_title, fontsize=9)
                    ax.grid(axis="y", alpha=0.2)

            fig.suptitle(
                f"{title} | features {start + 1}-{start + len(page_features)} of {n_features}",
                fontsize=12,
                fontweight="bold",
                y=0.995,
            )
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            _save_pdf_page_a4(pdf, fig, orientation="landscape")
            plt.close(fig)


def _export_clustering_diagnostics_csvs(
    adata,
    cluster_col,
    feature_cols,
    outdir,
    filename_prefix="behavioral_states_diagnostics",
):
    """Export CSV companions for clustering diagnostics panels."""
    if adata is None or adata.n_obs == 0:
        return {}
    if cluster_col not in adata.obs.columns:
        return {}

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ad = _ensure_umap(adata, n_neighbors=30, random_state=123)

    diagnostics_paths = {}

    # UMAP panel source data.
    umap = np.asarray(ad.obsm.get("X_umap", np.empty((ad.n_obs, 2))), dtype=float)
    if umap.ndim != 2 or umap.shape[0] != ad.n_obs:
        umap = np.zeros((ad.n_obs, 2), dtype=float)
    if umap.shape[1] < 2:
        umap = np.pad(umap, ((0, 0), (0, 2 - umap.shape[1])), mode="constant")
    umap_df = pd.DataFrame(
        {
            "obs_index": ad.obs.index.astype(str),
            "umap_1": umap[:, 0],
            "umap_2": umap[:, 1],
            "cluster_col": ad.obs[cluster_col].astype(str).to_numpy(),
        }
    )
    for c in ["sample_name", "TrackID", "position_t"]:
        if c in ad.obs.columns:
            umap_df[c] = ad.obs[c].to_numpy()
    umap_csv = outdir / f"{filename_prefix}_umap.csv"
    umap_df.to_csv(umap_csv, index=False)
    diagnostics_paths["umap_csv"] = str(umap_csv)

    valid_features = [str(c) for c in list(feature_cols or []) if str(c) in ad.var_names]
    cluster_order = sorted(ad.obs[cluster_col].astype(str).unique().tolist(), key=_mixed_label_sort_key)
    if len(valid_features) == 0:
        cluster_means = pd.DataFrame(index=cluster_order)
    else:
        X = _to_numpy_2d(ad[:, valid_features].X).astype(float, copy=False)
        feature_df = pd.DataFrame(X, columns=valid_features, index=ad.obs.index)
        feature_df["_cluster"] = ad.obs[cluster_col].astype(str).to_numpy()
        cluster_means = feature_df.groupby("_cluster", observed=True)[valid_features].mean()
        cluster_means = cluster_means.reindex(cluster_order)

    # Correlation panel source data (cluster-by-cluster correlation).
    if cluster_means.shape[0] > 0 and cluster_means.shape[1] > 0:
        corr_df = cluster_means.T.corr()
    else:
        corr_df = pd.DataFrame(index=cluster_means.index, columns=cluster_means.index, dtype=float)
    corr_csv = outdir / f"{filename_prefix}_cluster_correlation.csv"
    corr_df.to_csv(corr_csv)
    diagnostics_paths["cluster_correlation_csv"] = str(corr_csv)

    # Heatmap panel source data (cluster x feature, standard_scale='var' equivalent).
    heatmap_df = cluster_means.copy()
    if heatmap_df.shape[1] > 0:
        for col in heatmap_df.columns:
            series = pd.to_numeric(heatmap_df[col], errors="coerce")
            col_min = series.min()
            col_max = series.max()
            if pd.isna(col_min) or pd.isna(col_max) or float(col_max) == float(col_min):
                heatmap_df[col] = 0.0
            else:
                heatmap_df[col] = (series - col_min) / (col_max - col_min)
    heatmap_csv = outdir / f"{filename_prefix}_heatmap_matrix.csv"
    heatmap_df.to_csv(heatmap_csv)
    diagnostics_paths["heatmap_matrix_csv"] = str(heatmap_csv)

    return diagnostics_paths


def _ensure_umap(adata, n_neighbors=30, random_state=123):
    """Ensure AnnData has a UMAP embedding for plotting."""
    if "X_umap" in adata.obsm:
        return adata

    requested_n_neighbors = int(n_neighbors)
    neighbors_meta = adata.uns.get("neighbors", {}) if isinstance(adata.uns, dict) else {}
    params = neighbors_meta.get("params", {}) if isinstance(neighbors_meta, dict) else {}
    existing_n_neighbors = params.get("n_neighbors", None)
    conn_key = neighbors_meta.get("connectivities_key", "connectivities")
    dist_key = neighbors_meta.get("distances_key", "distances")
    has_graph = (conn_key in adata.obsp) and (dist_key in adata.obsp)

    can_reuse_neighbors = False
    try:
        can_reuse_neighbors = has_graph and (int(existing_n_neighbors) == requested_n_neighbors)
    except (TypeError, ValueError):
        can_reuse_neighbors = False

    if not can_reuse_neighbors:
        if "X_pca" in adata.obsm:
            sc.pp.neighbors(
                adata,
                n_neighbors=requested_n_neighbors,
                use_rep="X_pca",
                random_state=int(random_state),
            )
        else:
            sc.pp.neighbors(adata, n_neighbors=requested_n_neighbors, random_state=int(random_state))
    sc.tl.umap(adata, random_state=int(random_state))
    return adata


def check_model_clusterid_consistency(
    model_adata,
    adata_full,
    model_cluster_col="intrinsic_behavioral_cluster",
    full_cluster_col="intrinsic_behavioral_cluster",
    key_cols=("sample_name", "TrackID", "position_t"),
    max_examples=10,
    verbose=True,
):
    """
    Check whether model rows keep the same cluster label after label transfer into full data.

    Matching priority:
    1) key-based join on key_cols
    2) obs_names intersection (fallback only if keys are unavailable)
    """
    if model_cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{model_cluster_col}' in model_adata.obs.")
    if full_cluster_col not in adata_full.obs.columns:
        raise ValueError(f"Missing '{full_cluster_col}' in adata_full.obs.")

    model_obs = model_adata.obs.copy()
    full_obs = adata_full.obs.copy()
    model_obs["_cluster_model"] = model_obs[model_cluster_col].astype(str)
    full_obs["_cluster_full"] = full_obs[full_cluster_col].astype(str)

    missing_model_keys = [c for c in key_cols if c not in model_obs.columns]
    missing_full_keys = [c for c in key_cols if c not in full_obs.columns]
    keys_available = len(missing_model_keys) == 0 and len(missing_full_keys) == 0
    used_match_mode = "key_cols" if keys_available else "obs_names"

    if keys_available:
        keys = list(key_cols)
        model_cmp = model_obs[keys + ["_cluster_model"]].copy()
        full_cmp = full_obs[keys + ["_cluster_full"]].copy()
        for c in keys:
            model_cmp[c] = model_cmp[c].astype(str)
            full_cmp[c] = full_cmp[c].astype(str)

        model_dups = int(model_cmp.duplicated(subset=keys, keep=False).sum())
        full_dups = int(full_cmp.duplicated(subset=keys, keep=False).sum())
        if model_dups > 0 or full_dups > 0:
            raise ValueError(
                f"Cannot do one-to-one comparison: duplicate key rows found "
                f"(model duplicates={model_dups}, full duplicates={full_dups}) for key_cols={keys}."
            )

        merged = model_cmp.merge(full_cmp, on=keys, how="left", validate="one_to_one")
    else:
        overlap = model_obs.index.intersection(full_obs.index)
        if len(overlap) == 0:
            raise ValueError(
                "Cannot compare labels: no shared obs_names and key columns are missing for key-based matching. "
                f"missing in model_adata: {missing_model_keys}; missing in adata_full: {missing_full_keys}"
            )
        merged = pd.DataFrame(
            {
                "_cluster_model": model_obs.loc[overlap, "_cluster_model"].to_numpy(),
                "_cluster_full": full_obs.loc[overlap, "_cluster_full"].to_numpy(),
            },
            index=overlap,
        ).reset_index(drop=False).rename(columns={"index": "_row_id"})

    merged["_match"] = merged["_cluster_model"] == merged["_cluster_full"]
    missing_in_full = int(merged["_cluster_full"].isna().sum())
    n_total = int(len(merged))
    n_match = int(merged["_match"].sum())
    n_mismatch = int((~merged["_match"] & merged["_cluster_full"].notna()).sum())

    summary = {
        "match_mode": used_match_mode,
        "n_model_rows_checked": n_total,
        "n_matching_clusterid": n_match,
        "n_mismatched_clusterid": n_mismatch,
        "n_missing_rows_in_adata_full": missing_in_full,
        "match_rate": float(n_match / n_total) if n_total > 0 else np.nan,
    }

    _vinfo(
        verbose,
        "state-apply",
        "clusterid consistency check "
        f"(mode={summary['match_mode']}): "
        f"checked={summary['n_model_rows_checked']}, "
        f"matches={summary['n_matching_clusterid']}, "
        f"mismatches={summary['n_mismatched_clusterid']}, "
        f"missing_in_full={summary['n_missing_rows_in_adata_full']}",
    )
    if n_mismatch > 0:
        _vinfo(verbose, "state-apply", f"mismatch preview suppressed (n_mismatch={n_mismatch})")

    return summary, merged


def _matches_model_sampling_cache(
    adata,
    *,
    kept_features,
    min_spacing_used,
    start_cutoff_mode=None,
    start_cutoff_window_size=None,
    max_samples,
    random_state,
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    nan_policy="drop_any_nan_in_kept_features",
):
    """Check whether cached model adata matches sampling/cleaning stage settings."""
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["model cache uns missing/invalid"]
    pre_meta = adata.uns.get("preprocessing", {})
    if not isinstance(pre_meta, dict):
        return False, ["model cache preprocessing metadata missing/invalid"]
    model_cache = pre_meta.get("model_cache", {})
    if not isinstance(model_cache, dict):
        return False, ["model cache metadata missing"]

    features_meta = model_cache.get("features", {})
    if isinstance(features_meta, dict):
        kept_cached = features_meta.get("kept_features", [])
    else:
        kept_cached = features_meta
    if kept_cached is None:
        kept_cached = []
    kept_cached = [str(x) for x in list(kept_cached)]

    if kept_features is None:
        kept_features = []
    kept_expected = [str(x) for x in list(kept_features)]
    if kept_cached != kept_expected:
        reasons.append("model_cache.features.kept_features mismatch")
    var_names = [str(x) for x in list(getattr(adata, "var_names", []))]
    if var_names != kept_expected:
        reasons.append("model cache var_names mismatch")

    sampling_meta = model_cache.get("sampling", {})
    sampling_meta = sampling_meta if isinstance(sampling_meta, dict) else {}
    try:
        got_spacing = int(sampling_meta.get("min_spacing_used", -1))
    except Exception:
        got_spacing = -1
    if got_spacing != int(min_spacing_used):
        reasons.append("model_cache.sampling.min_spacing_used mismatch")
    got_max_samples = sampling_meta.get("max_samples", None)
    got_max_samples = None if got_max_samples is None else int(got_max_samples)
    exp_max_samples = None if max_samples is None else int(max_samples)
    if got_max_samples != exp_max_samples:
        reasons.append("model_cache.sampling.max_samples mismatch")
    try:
        got_random_state = int(sampling_meta.get("random_state", -1))
    except Exception:
        got_random_state = -1
    if got_random_state != int(random_state):
        reasons.append("model_cache.sampling.random_state mismatch")
    if [str(x) for x in list(sampling_meta.get("id_cols", []))] != [str(x) for x in list(id_cols)]:
        reasons.append("model_cache.sampling.id_cols mismatch")
    if str(sampling_meta.get("time_col", "")) != str(time_col):
        reasons.append("model_cache.sampling.time_col mismatch")
    got_start_cutoff_mode = sampling_meta.get("start_cutoff_mode", None)
    got_start_cutoff_mode = None if got_start_cutoff_mode is None else str(got_start_cutoff_mode)
    exp_start_cutoff_mode = None if start_cutoff_mode is None else str(start_cutoff_mode)
    if got_start_cutoff_mode != exp_start_cutoff_mode:
        reasons.append("model_cache.sampling.start_cutoff_mode mismatch")
    got_start_cutoff_window = sampling_meta.get("start_cutoff_window_size", None)
    got_start_cutoff_window = None if got_start_cutoff_window is None else int(got_start_cutoff_window)
    exp_start_cutoff_window = (
        None if start_cutoff_window_size is None else int(start_cutoff_window_size)
    )
    if got_start_cutoff_window != exp_start_cutoff_window:
        reasons.append("model_cache.sampling.start_cutoff_window_size mismatch")

    cleaning_meta = model_cache.get("cleaning", {})
    cleaning_meta = cleaning_meta if isinstance(cleaning_meta, dict) else {}
    if str(cleaning_meta.get("nan_policy", "")) != str(nan_policy):
        reasons.append("model_cache.cleaning.nan_policy mismatch")

    return len(reasons) == 0, reasons


def _representation_label(use_pca):
    return "X_pca" if bool(use_pca) else "X"


def _scanpy_use_rep(use_pca):
    return "X_pca" if bool(use_pca) else None


def _matches_cached_pca_stage(
    adata,
    *,
    use_pca,
    pca_var_selection,
    svd_solver,
    random_state,
    ncomps_requested,
):
    """Check whether cached PCA output matches current PCA settings."""
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    pca_meta = model_cache.get("pca", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(pca_meta, dict):
        return False, ["model_cache.pca missing"]
    got_enabled = pca_meta.get("enabled", None)
    if got_enabled is None:
        return False, ["model_cache.pca.enabled missing"]
    if bool(got_enabled) != bool(use_pca):
        reasons.append("model_cache.pca.enabled mismatch")

    if not bool(use_pca):
        if "X_pca" in getattr(adata, "obsm", {}):
            reasons.append("X_pca present when PCA disabled")
        if "pca" in adata.uns:
            reasons.append("pca metadata present when PCA disabled")
        return len(reasons) == 0, reasons

    if "X_pca" not in adata.obsm:
        return False, reasons + ["X_pca missing"]

    got_pca_var = pca_meta.get("pca_var_selection", np.nan)
    if not np.isclose(float(got_pca_var), float(pca_var_selection), rtol=1e-8, atol=1e-12):
        reasons.append("model_cache.pca.pca_var_selection mismatch")
    if str(pca_meta.get("svd_solver", "")) != str(svd_solver):
        reasons.append("model_cache.pca.svd_solver mismatch")
    try:
        got_random_state = int(pca_meta.get("random_state", -1))
    except Exception:
        got_random_state = -1
    if got_random_state != int(random_state):
        reasons.append("model_cache.pca.random_state mismatch")
    try:
        got_ncomps_requested = int(pca_meta.get("ncomps_requested", -1))
    except Exception:
        got_ncomps_requested = -1
    if got_ncomps_requested != int(ncomps_requested):
        reasons.append("model_cache.pca.ncomps_requested mismatch")
    x_pca = adata.obsm.get("X_pca", None)
    ncomps_effective = -1 if x_pca is None else int(x_pca.shape[1])
    if int(pca_meta.get("ncomps_effective", -1)) != ncomps_effective:
        reasons.append("model_cache.pca.ncomps_effective mismatch")

    return len(reasons) == 0, reasons


def _matches_cached_neighbors_stage(
    adata,
    *,
    n_neighbors,
    random_state,
    use_rep="X_pca",
):
    """Check whether cached neighbors graph matches current neighbor settings."""
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    neighbors_uns = adata.uns.get("neighbors", {})
    neighbors_uns = neighbors_uns if isinstance(neighbors_uns, dict) else {}
    params = neighbors_uns.get("params", {})
    params = params if isinstance(params, dict) else {}
    conn_key = neighbors_uns.get("connectivities_key", "connectivities")
    dist_key = neighbors_uns.get("distances_key", "distances")
    has_graph = (conn_key in adata.obsp) and (dist_key in adata.obsp)
    if not has_graph:
        reasons.append("neighbors graph missing")
    try:
        if int(params.get("n_neighbors", -1)) != int(n_neighbors):
            reasons.append("neighbors params.n_neighbors mismatch")
    except Exception:
        reasons.append("neighbors params.n_neighbors invalid")

    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    neighbors_meta = model_cache.get("neighbors", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(neighbors_meta, dict):
        return False, reasons + ["model_cache.neighbors missing"]
    if int(neighbors_meta.get("n_neighbors", -1)) != int(n_neighbors):
        reasons.append("model_cache.neighbors.n_neighbors mismatch")
    if int(neighbors_meta.get("random_state", -1)) != int(random_state):
        reasons.append("model_cache.neighbors.random_state mismatch")
    if str(neighbors_meta.get("use_rep", "")) != str(use_rep):
        reasons.append("model_cache.neighbors.use_rep mismatch")

    return len(reasons) == 0, reasons


def _matches_cached_umap_stage(
    adata,
    *,
    min_dist,
    random_state,
    use_rep="X_pca",
):
    """Check whether cached UMAP embedding matches current UMAP settings."""
    reasons = []
    if "X_umap" not in adata.obsm:
        return False, ["X_umap missing"]
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    umap_meta = model_cache.get("umap", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(umap_meta, dict):
        return False, ["model_cache.umap missing"]
    if float(umap_meta.get("min_dist", np.nan)) != float(min_dist):
        reasons.append("model_cache.umap.min_dist mismatch")
    if int(umap_meta.get("random_state", -1)) != int(random_state):
        reasons.append("model_cache.umap.random_state mismatch")
    if str(umap_meta.get("use_rep", "")) != str(use_rep):
        reasons.append("model_cache.umap.use_rep mismatch")
    return len(reasons) == 0, reasons


def train_random_forest_classifier(
    X_train,
    y_train,
    random_state=123,
    n_estimators=300,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
):
    """Train a RandomForestClassifier and return (classifier, fit_info)."""
    X_train = np.asarray(X_train, dtype=float)
    y_train = np.asarray(y_train).astype(str)
    if X_train.ndim != 2:
        raise ValueError(f"X_train must be 2D, got shape {X_train.shape}.")
    if len(y_train) == 0:
        raise ValueError("No rows available to train RandomForest classifier.")
    if X_train.shape[0] != len(y_train):
        raise ValueError(
            f"X/y row mismatch for RandomForest training: X rows={X_train.shape[0]}, y rows={len(y_train)}."
        )

    params = {
        "n_estimators": int(n_estimators),
        "min_samples_leaf": int(min_samples_leaf),
        "n_jobs": int(n_jobs),
        "random_state": int(random_state),
        "max_depth": None if max_depth is None else int(max_depth),
        "min_samples_split": int(min_samples_split),
        "max_features": max_features,
        "class_weight": class_weight,
    }

    clf = RandomForestClassifier(**params, oob_score=True)
    clf.fit(X_train, y_train)
    train_acc = float((clf.predict(X_train) == y_train).mean())
    fit_info = {
        "backend": "random_forest",
        "n_train_rows": int(X_train.shape[0]),
        "n_features": int(X_train.shape[1]),
        "classes": [str(c) for c in getattr(clf, "classes_", [])],
        "train_accuracy": train_acc,
        "params": {k: repr(v) for k, v in clf.get_params().items()},
        "oob_score": float(getattr(clf, "oob_score_", np.nan)),
    }
    return clf, fit_info


def train_clusterid_classifier_from_model_adata(
    model_adata,
    cluster_col="intrinsic_behavioral_cluster",
    classifier_backend="random_forest",
    random_state=123,
    n_estimators=300,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    verbose=True,
):
    """Train a supervised classifier on model_adata features -> cluster labels."""
    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = list(model_adata.var_names)
    X_train = _to_numpy_2d(model_adata[:, feature_cols].X)
    y_train = model_adata.obs[cluster_col].astype(str).to_numpy()

    backend = str(classifier_backend).lower()
    if backend not in {"random_forest", "rf", "randomforest", "randomforestclassifier"}:
        raise ValueError("Only RandomForestClassifier is currently supported (classifier_backend='random_forest').")

    clf, fit_info = train_random_forest_classifier(
        X_train=X_train,
        y_train=y_train,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
    )
    fit_info["cluster_col"] = cluster_col
    _vinfo(
        verbose,
        "state-training",
        (
            "Trained intrinsic classifier "
            f"(backend={fit_info['backend']}, rows={fit_info['n_train_rows']}, "
            f"features={fit_info['n_features']}, train_acc={fit_info['train_accuracy']:.4f})"
        ),
    )
    return clf, fit_info


def _resolve_validation_mode(validation_test_size):
    """Resolve validation mode from user-specified test size."""
    if validation_test_size is None or float(validation_test_size) == 0.0:
        return "oob_only"
    test_size = float(validation_test_size)
    if not (0.0 < test_size < 1.0):
        raise ValueError(
            "validation_test_size must be in (0, 1), or set to 0/None for OOB-only mode."
        )
    return "holdout"


def _validate_stratified_split_feasibility(y, test_size, target_name):
    """Raise ValueError if stratified holdout split cannot be constructed."""
    y = np.asarray(y).astype(str)
    class_counts = pd.Series(y).value_counts(dropna=False)
    n_samples = int(len(y))
    n_classes = int(class_counts.shape[0])

    if n_classes < 2:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': need at least 2 classes, got {n_classes}."
        )
    if int(class_counts.min()) < 2:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': each class needs >=2 rows; "
            f"min class count={int(class_counts.min())}. Class counts={class_counts.to_dict()}"
        )

    test_size = float(test_size)
    n_test = int(np.ceil(n_samples * test_size))
    n_train = n_samples - n_test
    if n_test < n_classes:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': test rows ({n_test}) "
            f"must be >= number of classes ({n_classes}). "
            f"n_samples={n_samples}, test_size={test_size}, class_counts={class_counts.to_dict()}"
        )
    if n_train < n_classes:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': train rows ({n_train}) "
            f"must be >= number of classes ({n_classes}). "
            f"n_samples={n_samples}, test_size={test_size}, class_counts={class_counts.to_dict()}"
        )


def _build_compact_metric_summary(y_true, y_pred):
    """Build compact scalar classification metrics."""
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    label_set = pd.Index(np.concatenate([y_true, y_pred])).unique()
    return {
        "n_rows": int(len(y_true)),
        "n_classes": int(len(label_set)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
    }


def _evaluate_holdout_predictions(
    *,
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="classifier_holdout_qc",
    row_normalized_decimals=5,
):
    """Evaluate holdout predictions and return compact metrics + saved artifacts."""
    summary = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=outfolder,
        filename_prefix=filename_prefix,
        row_normalized_decimals=row_normalized_decimals,
    )
    compact = _build_compact_metric_summary(y_true, y_pred)
    return {
        "metrics": compact,
        "artifacts": dict(summary.get("artifacts", {})),
        "summary": summary,
    }


def evaluate_classifier_predictions(
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="classifier_qc",
    row_normalized_decimals=5,
):
    """Evaluate predicted labels and optionally save confusion/report artifacts."""
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    label_set = pd.Index(np.concatenate([y_true, y_pred])).unique().tolist()
    labels = sorted([str(x) for x in label_set], key=_mixed_label_sort_key)

    cm_counts = confusion_matrix(y_true, y_pred, labels=labels)
    cm_true_norm = confusion_matrix(y_true, y_pred, labels=labels, normalize="true")
    report_dict = classification_report(y_true, y_pred, labels=labels, output_dict=True, zero_division=0)

    summary = {
        "n_rows": int(len(y_true)),
        "n_classes": int(len(labels)),
        "classes": [str(c) for c in labels],
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "artifacts": {},
    }

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)

        cm_counts_df = pd.DataFrame(cm_counts, index=labels, columns=labels)
        cm_norm_df = pd.DataFrame(cm_true_norm, index=labels, columns=labels)
        report_df = pd.DataFrame(report_dict).T

        cm_counts_csv = outfolder / f"{filename_prefix}_confusion_counts.csv"
        cm_norm_csv = outfolder / f"{filename_prefix}_confusion_true_normalized.csv"
        report_csv = outfolder / f"{filename_prefix}_classification_report.csv"
        cm_pdf = outfolder / f"{filename_prefix}_confusion_matrices.pdf"

        cm_counts_df.to_csv(cm_counts_csv)
        cm_norm_df.to_csv(cm_norm_csv)
        report_df.to_csv(report_csv)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

        im0 = axes[0].imshow(cm_counts, cmap="Blues", aspect="auto")
        axes[0].set_title("Confusion Matrix (Counts)")
        axes[0].set_xlabel("Predicted label")
        axes[0].set_ylabel("True label")
        axes[0].set_xticks(np.arange(len(labels)))
        axes[0].set_yticks(np.arange(len(labels)))
        axes[0].set_xticklabels(labels, rotation=90)
        axes[0].set_yticklabels(labels)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

        im1 = axes[1].imshow(cm_true_norm, cmap="Blues", aspect="auto", vmin=0.0, vmax=1.0)
        axes[1].set_title("Confusion Matrix (Row-normalized)")
        axes[1].set_xlabel("Predicted label")
        axes[1].set_ylabel("True label")
        axes[1].set_xticks(np.arange(len(labels)))
        axes[1].set_yticks(np.arange(len(labels)))
        axes[1].set_xticklabels(labels, rotation=90)
        axes[1].set_yticklabels(labels)
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

        if len(labels) <= 20:
            for i in range(len(labels)):
                for j in range(len(labels)):
                    axes[0].text(j, i, str(int(cm_counts[i, j])), ha="center", va="center", fontsize=7, color="black")
                    axes[1].text(
                        j,
                        i,
                        f"{cm_true_norm[i, j]:.{int(row_normalized_decimals)}f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="black",
                    )

        fig.suptitle("Classifier QC", fontsize=12)
        fig.savefig(cm_pdf, dpi=300, bbox_inches="tight")
        plt.close(fig)

        summary["artifacts"] = {
            "confusion_counts_csv": str(cm_counts_csv),
            "confusion_true_normalized_csv": str(cm_norm_csv),
            "classification_report_csv": str(report_csv),
            "confusion_matrices_pdf": str(cm_pdf),
        }

    return summary


def evaluate_clusterid_classifier_on_model_adata(
    model_adata,
    classifier,
    cluster_col="intrinsic_behavioral_cluster",
    outfolder=None,
    filename_prefix="behavioral_states_classifier_qc",
    row_normalized_decimals=5,
    verbose=True,
):
    """Evaluate classifier quality on model_adata and optionally save QC outputs."""
    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = list(model_adata.var_names)
    X_eval = _to_numpy_2d(model_adata[:, feature_cols].X)
    y_true = model_adata.obs[cluster_col].astype(str).to_numpy()
    y_pred = np.asarray(classifier.predict(X_eval)).astype(str)

    summary = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=outfolder,
        filename_prefix=filename_prefix,
        row_normalized_decimals=row_normalized_decimals,
    )

    _vinfo(
        verbose,
        "state-training",
        (
            "Intrinsic model-data QC: "
            f"n={summary['n_rows']}, classes={summary['n_classes']}, "
            f"accuracy={summary['accuracy']:.4f}, balanced_accuracy={summary['balanced_accuracy']:.4f}"
        ),
    )

    return summary


def _resolve_binary_classifier_feature_cols(adata, binary_cols_to_merge):
    """Resolve binary classifier features strictly from supplied binary columns."""
    cols = []
    for col in list(binary_cols_to_merge):
        if col in adata.obs.columns:
            cols.append(col)
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(cols))


def _build_classifier_matrix_from_adata(
    adata,
    continuous_feature_cols,
    binary_feature_cols,
    binary_expanded_feature_cols=None,
    return_binary_feature_names=False,
    row_indices=None,
):
    """Build a numeric feature matrix from adata.X (continuous) and obs (binary)."""
    cont_cols = list(continuous_feature_cols)
    bin_cols = list(binary_feature_cols)

    missing_cont = [c for c in cont_cols if c not in adata.var_names]
    if len(missing_cont) > 0:
        raise ValueError(f"Missing continuous feature columns in adata.var_names: {missing_cont[:10]}")

    if row_indices is None:
        adata_cont = adata[:, cont_cols]
        obs_df = adata.obs
    else:
        row_indices = np.asarray(row_indices)
        adata_cont = adata[row_indices, cont_cols]
        if np.issubdtype(row_indices.dtype, np.integer):
            obs_df = adata.obs.iloc[row_indices]
        else:
            obs_df = adata.obs.loc[row_indices]

    X_cont = _to_numpy_2d(adata_cont.X).astype(float, copy=False)
    if len(bin_cols) == 0:
        if return_binary_feature_names:
            return X_cont, []
        return X_cont

    bin_df = _build_binary_feature_dataframe(
        obs_df=obs_df,
        binary_feature_cols=bin_cols,
        expected_expanded_cols=binary_expanded_feature_cols,
    )
    X_bin = bin_df.to_numpy(dtype=float)
    X = np.hstack([X_cont, X_bin])
    if return_binary_feature_names:
        return X, list(bin_df.columns)
    return X


def train_full_clusterid_classifier_from_adata(
    adata,
    target_col="full_behavioral_cluster",
    binary_feature_cols=None,
    preprocessing_artifact=None,
    random_state=123,
    n_estimators=500,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    verbose=True,
):
    """
    Train final random-forest classifier using continuous + binary features.

    Returns classifier artifact and fit metadata.
    """
    if target_col not in adata.obs.columns:
        raise ValueError(f"Missing target column '{target_col}' in adata.obs.")

    continuous_feature_cols = list(adata.var_names)
    binary_cols = [] if binary_feature_cols is None else list(binary_feature_cols)
    X_train, expanded_binary_cols = _build_classifier_matrix_from_adata(
        adata=adata,
        continuous_feature_cols=continuous_feature_cols,
        binary_feature_cols=binary_cols,
        return_binary_feature_names=True,
    )
    y_train = adata.obs[target_col].astype(str).to_numpy()
    if len(y_train) == 0:
        raise ValueError("No rows available to train full ClusterID classifier.")

    clf, fit_info = train_random_forest_classifier(
        X_train=X_train,
        y_train=y_train,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
    )

    artifact = {
        "classifier": clf,
        "continuous_feature_cols": list(continuous_feature_cols),
        "binary_feature_cols": list(binary_cols),
        "binary_expanded_feature_cols": list(expanded_binary_cols),
        "target_col": str(target_col),
        "preprocessing": preprocessing_artifact,
    }
    fit_info["target_col"] = str(target_col)
    fit_info["continuous_feature_cols"] = list(continuous_feature_cols)
    fit_info["binary_feature_cols"] = list(binary_cols)
    fit_info["n_features_total"] = int(X_train.shape[1])
    fit_info["n_features_continuous"] = int(len(continuous_feature_cols))
    fit_info["n_features_binary"] = int(len(binary_cols))
    fit_info["n_features_binary_expanded"] = int(len(expanded_binary_cols))
    _vinfo(
        verbose,
        "state-training",
        (
            "Trained full classifier "
            f"(rows={fit_info['n_train_rows']}, total_features={fit_info['n_features_total']}, "
            f"train_acc={fit_info['train_accuracy']:.4f})"
        ),
    )
    return artifact, fit_info


def train_full_classifier_on_labeled_adata(
    adata,
    binary_features_to_group=None,
    preprocessing_artifact=None,
    random_state=123,
    n_estimators=500,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    output_col="full_behavioral_cluster",
    confidence_col="full_behavioral_cluster_confidence",
    outfolder=None,
    qc_filename_prefix="behavioral_states_full_classifier_qc_random_forest",
    windowing_artifact=None,
    verbose=True,
):
    """
    Train/apply/evaluate the full classifier (continuous + binary) on a labeled AnnData.

    This lets you run full-classifier training separately from the continuous-only stage.
    """
    preprocessing_params = preprocessing_artifact
    if preprocessing_params is None:
        preprocessing_params = _build_preprocessing_params_from_uns(
            adata.uns.get("preprocessing", {}) if isinstance(adata.uns, dict) else None,
            feature_cols=list(adata.var_names),
        )
    windowing_params = windowing_artifact
    if binary_features_to_group is None:
        binary_features_to_group = []
        if isinstance(adata.uns, dict):
            binary_features_to_group = list(
                adata.uns.get("clustering", {}).get("binary_cols_to_merge", [])
            )

    full_label_classifier_artifact, full_classifier_fit_info = train_full_clusterid_classifier_from_adata(
        adata=adata,
        target_col="full_behavioral_cluster",
        binary_feature_cols=binary_features_to_group,
        preprocessing_artifact=preprocessing_params,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
        verbose=verbose,
    )
    if windowing_params is not None:
        full_label_classifier_artifact["windowing"] = dict(windowing_params)

    predict_clusterids_with_full_classifier(
        adata=adata,
        classifier_artifact=full_label_classifier_artifact,
        output_col=output_col,
        confidence_col=confidence_col,
        inplace=True,
        apply_preprocessing=False,
    )
    y_true = adata.obs["full_behavioral_cluster"].astype(str).to_numpy()
    y_pred = adata.obs[output_col].astype(str).to_numpy()
    full_classifier_qc_model_data = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=(Path(outfolder) if outfolder is not None else None),
        filename_prefix=qc_filename_prefix,
    )
    full_cm_pdf = full_classifier_qc_model_data.get("artifacts", {}).get("confusion_matrices_pdf", None)
    if full_cm_pdf is not None:
        _vsave(verbose, "state-training", "full-classifier confusion matrices", full_cm_pdf)
    _vinfo(
        verbose,
        "state-training",
        (
            "Full model-data QC: "
            f"n={full_classifier_qc_model_data['n_rows']}, "
            f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
            f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
        ),
    )
    if "classification" not in adata.uns:
        adata.uns["classification"] = {}
    adata.uns["classification"]["full_classifier_fit_info"] = full_classifier_fit_info
    adata.uns["classification"]["full_classifier_qc_model_data"] = full_classifier_qc_model_data
    adata.uns["classification"]["full_classifier_prediction_col"] = str(output_col)
    adata.uns["classification"]["full_classifier_confidence_col"] = (
        None if confidence_col is None else str(confidence_col)
    )
    return full_label_classifier_artifact, full_classifier_fit_info, full_classifier_qc_model_data


def plot_binary_group_behavioral_cluster_grid(
    adata,
    binary_group_col="binary_group",
    behavioral_cluster_col="behavioral_clusterid",
    ncols=1,
    figsize_per_plot=(10.0, 1.8),
    csv_path=None,
    pdf_path=None,
    return_csv=True,
    verbose=True,
):
    """
    Plot one horizontal stacked bar per binary group, showing proportions of behavioral clusters.
    """
    obs = adata.obs.copy()
    if binary_group_col not in obs.columns:
        raise ValueError(f"Missing '{binary_group_col}' in adata.obs.")
    if behavioral_cluster_col not in obs.columns:
        raise ValueError(f"Missing '{behavioral_cluster_col}' in adata.obs.")

    plot_df = obs[[binary_group_col, behavioral_cluster_col]].copy()
    plot_df[binary_group_col] = plot_df[binary_group_col].astype("string").fillna("unassigned")
    plot_df[behavioral_cluster_col] = plot_df[behavioral_cluster_col].astype("string").fillna("unassigned")

    # Keep all behavioral clusters present in the dataset (even if concentrated in
    # one binary group), so no cluster silently disappears in the grid.
    cluster_counts = plot_df[behavioral_cluster_col].value_counts(dropna=False)
    all_clusters = cluster_counts.index.tolist()

    ctab = pd.crosstab(
        plot_df[binary_group_col],
        plot_df[behavioral_cluster_col],
        normalize="index",
        dropna=False,
    ).reindex(columns=all_clusters, fill_value=0.0)

    def _sort_key(x):
        s = str(x)
        return (0, int(s)) if s.isdigit() else (1, s)

    group_names = sorted(ctab.index.tolist(), key=_sort_key)
    cluster_names = sorted(ctab.columns.tolist(), key=_sort_key)

    n_groups = len(group_names)
    if n_groups == 0:
        raise ValueError("No groups available for plotting.")

    nrows = int(np.ceil(n_groups / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    # Build a larger distinct palette by combining qualitative colormaps.
    palette = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(cmap_name)
        palette.extend([cmap(i) for i in range(cmap.N)])
    if len(cluster_names) > len(palette):
        cmap = plt.get_cmap("hsv")
        palette.extend([cmap(i / max(len(cluster_names), 1)) for i in range(len(cluster_names) - len(palette))])
    colors = {cl: palette[i] for i, cl in enumerate(cluster_names)}

    for i, grp in enumerate(group_names):
        ax = axes[i]
        left = 0.0
        for cl in cluster_names:
            val = float(ctab.loc[grp, cl])
            if val <= 0:
                continue
            ax.barh([0], [val], left=left, color=colors[cl], height=0.92, edgecolor="none", linewidth=0.0)
            left += val
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(False)
        ax.text(
            -0.01,
            0.5,
            str(grp),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
        )
        grp_n = int((plot_df[binary_group_col] == grp).sum())
        ax.set_title(f"n={grp_n}", fontsize=9, loc="right")

    for j in range(n_groups, len(axes)):
        axes[j].axis("off")

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="none", color=colors[cl], label=str(cl), markersize=7)
        for cl in cluster_names
    ]
    fig.legend(
        handles=handles,
        labels=[str(c) for c in cluster_names],
        title="behavioral_clusterid",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
    )
    fig.suptitle("Behavioral Cluster Composition per Binary Group", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.82, 0.97])

    if verbose and "5" in ctab.columns:
        by_group = (ctab["5"] * 100.0).round(3)
        nonzero = int((by_group > 0).sum())
        _vinfo(
            verbose,
            "state-apply",
            f"Cluster 5 composition summary: total_n={int(cluster_counts.get('5', 0))}, nonzero_groups={nonzero}",
        )

    if csv_path is None and pdf_path is not None:
        csv_path = Path(pdf_path).with_suffix(".csv")
    if csv_path is not None:
        ctab.to_csv(csv_path)
    if pdf_path is not None:
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")

    if return_csv:
        return fig, ctab
    return fig


def plot_behavioral_cluster_binary_group_grid(
    adata,
    binary_group_col="binary_group",
    behavioral_cluster_col="behavioral_clusterid",
    ncols=1,
    figsize_per_plot=(10.0, 1.8),
    csv_path=None,
    pdf_path=None,
    return_csv=True,
):
    """
    Plot one horizontal stacked bar per behavioral cluster, showing proportions of binary groups.
    """
    obs = adata.obs.copy()
    if binary_group_col not in obs.columns:
        raise ValueError(f"Missing '{binary_group_col}' in adata.obs.")
    if behavioral_cluster_col not in obs.columns:
        raise ValueError(f"Missing '{behavioral_cluster_col}' in adata.obs.")

    plot_df = obs[[binary_group_col, behavioral_cluster_col]].copy()
    plot_df[binary_group_col] = plot_df[binary_group_col].astype("string").fillna("unassigned")
    plot_df[behavioral_cluster_col] = plot_df[behavioral_cluster_col].astype("string").fillna("unassigned")

    # Keep all binary groups present in the dataset so no group disappears.
    group_counts = plot_df[binary_group_col].value_counts(dropna=False)
    all_groups = group_counts.index.tolist()

    ctab = pd.crosstab(
        plot_df[behavioral_cluster_col],
        plot_df[binary_group_col],
        normalize="index",
        dropna=False,
    ).reindex(columns=all_groups, fill_value=0.0)

    def _sort_key(x):
        s = str(x)
        return (0, int(s)) if s.isdigit() else (1, s)

    cluster_names = sorted(ctab.index.tolist(), key=_sort_key)
    group_names = sorted(ctab.columns.tolist(), key=_sort_key)

    n_clusters = len(cluster_names)
    if n_clusters == 0:
        raise ValueError("No behavioral clusters available for plotting.")

    nrows = int(np.ceil(n_clusters / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    palette = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(cmap_name)
        palette.extend([cmap(i) for i in range(cmap.N)])
    if len(group_names) > len(palette):
        cmap = plt.get_cmap("hsv")
        palette.extend([cmap(i / max(len(group_names), 1)) for i in range(len(group_names) - len(palette))])
    colors = {grp: palette[i] for i, grp in enumerate(group_names)}

    for i, cl in enumerate(cluster_names):
        ax = axes[i]
        left = 0.0
        for grp in group_names:
            val = float(ctab.loc[cl, grp])
            if val <= 0:
                continue
            ax.barh([0], [val], left=left, color=colors[grp], height=0.92, edgecolor="none", linewidth=0.0)
            left += val
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(False)
        ax.text(
            -0.01,
            0.5,
            str(cl),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
        )
        cl_n = int((plot_df[behavioral_cluster_col] == cl).sum())
        ax.set_title(f"n={cl_n}", fontsize=9, loc="right")

    for j in range(n_clusters, len(axes)):
        axes[j].axis("off")

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="none", color=colors[grp], label=str(grp), markersize=7)
        for grp in group_names
    ]
    fig.legend(
        handles=handles,
        labels=[str(g) for g in group_names],
        title="binary_group",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
    )
    fig.suptitle("Binary Group Composition per Behavioral Cluster", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.82, 0.97])

    if csv_path is None and pdf_path is not None:
        csv_path = Path(pdf_path).with_suffix(".csv")
    if csv_path is not None:
        ctab.to_csv(csv_path)
    if pdf_path is not None:
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")

    if return_csv:
        return fig, ctab
    return fig


def run_state_clustering(
    features,
    binary_features_to_group,
    output_dir,
    cell_type,
    window_size=5,
    min_spacing=None,
    max_samples=None,
    n_neighbors=60,
    min_dist=0.1,
    resolution=0.2,
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    pca_var_selection=0.95,
    use_pca=True,
    clustering_method="leiden",
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    log_scale_features=None,
    incomplete_window_policy="partial",
    random_state=123,
    reuse_prepared_dataset=True,
    plot_exemplar_videos=True,
    exemplar_n_per_cluster=5,
    exemplar_video_fps=6,
    exemplar_video_dpi=180,
    exemplar_video_margin=20,
    exemplar_video_pmin=0.0,
    exemplar_video_pmax=99.99,
    exemplar_video_track_color="#63ff33",
    exemplar_video_show_segments=True,
    exemplar_video_segment_style="outline",
    exemplar_video_segment_color="#ffffff",
    # save_prepared_dataset=True,
    df_positions=None,
    state_paths=None,
    verbose=True,
):
    """Stage 1: load/prepare positions, cluster, return model_adata.

    Stage-1 metadata is stored in:
    - model_adata.uns["preprocessing"]
    - model_adata.uns["clustering"]
    - model_adata.uns["preprocessing"]["prepared_adata_full_*"] cache metadata
    """
    run_started = _vstart(verbose, "state-clustering", "run state clustering")
    _vinfo(
        verbose,
        "state-clustering",
        f"requested descriptive_features={list(descriptive_features)} | use_pca={bool(use_pca)}",
    )
    use_pca = bool(use_pca)
    representation_label = _representation_label(use_pca)
    scanpy_use_rep = _scanpy_use_rep(use_pca)
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    state_outdir = state_paths.state_outdir
    state_clustering_outdir = _resolve_state_clustering_outdir(state_paths.processing_outdir)
    prepared_adata_full_path = state_paths.prepared_full_adata_path
    legacy_prepared_adata_full_path = _resolve_legacy_prepared_full_adata_path(
        state_paths=state_paths,
        cell_type=cell_type,
    )
    prepared_adata_read_path = None
    if prepared_adata_full_path is not None and prepared_adata_full_path.exists():
        prepared_adata_read_path = prepared_adata_full_path
    elif legacy_prepared_adata_full_path.exists():
        prepared_adata_read_path = legacy_prepared_adata_full_path
    prepared_adata_model_path = state_paths.model_adata_path
    df_tracks_path = None
    stage_prepare_started = time.perf_counter()
    adata_prepared = None
    prepared_cache_used = False
    positions_loaded_from_csv = False

    def _load_positions_if_needed():
        nonlocal df_positions, df_tracks_path, positions_loaded_from_csv
        if df_positions is None:
            df_tracks_path = _resolve_positions_csv_path(state_paths=state_paths)
            df_positions = pd.read_csv(df_tracks_path)
            positions_loaded_from_csv = True
        _require_columns(
            df_positions,
            ["sample_name", "TrackID", "position_t"] + list(features),
            context="state clustering input",
        )
        return df_positions

    if bool(reuse_prepared_dataset) and prepared_adata_read_path is not None:
        try:
            cached_prepared = sc.read_h5ad(prepared_adata_read_path)
            pre_match_cached, pre_reasons_cached = _matches_requested_preprocessing_in_adata(
                cached_prepared,
                features=features,
                binary_features_to_group=binary_features_to_group,
                window_size=window_size,
                min_spacing=min_spacing,
                descriptive_features=descriptive_features,
                incomplete_window_policy=incomplete_window_policy,
                lower_quantile_cap=lower_quantile_cap,
                upper_quantile_cap=upper_quantile_cap,
                log_scale_features=log_scale_features,
                scale_features=True,
            )
            if pre_match_cached:
                adata_prepared = cached_prepared
                prepared_cache_used = True
                _vinfo(verbose, "state-clustering", f"reusing prepared full adata cache: {prepared_adata_read_path}")
                if prepared_adata_read_path != prepared_adata_full_path:
                    _vinfo(
                        verbose,
                        "state-clustering",
                        f"using legacy prepared cache for read; future writes use {prepared_adata_full_path}",
                    )
            else:
                _vinfo(
                    verbose,
                    "state-clustering",
                    "prepared full cache mismatch; rebuilding from raw positions | "
                    f"{_short_reasons(pre_reasons_cached)}",
                )
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"could not load prepared full cache; rebuilding ({exc})")

    if adata_prepared is None:
        _load_positions_if_needed()
        _vinfo(
            verbose,
            "state-clustering",
            "preparing dataset (windowed feature extraction, quantile capping, optional log scaling, scaling)",
        )
        if positions_loaded_from_csv:
            _vinfo(verbose, "state-clustering", f"input track-features CSV: {df_tracks_path}")
        else:
            _vinfo(verbose, "state-clustering", "input source: provided df_positions DataFrame")

        prepare_kwargs = {
            "df_positions": df_positions,
            "features": features,
            "binary_features_to_group": binary_features_to_group,
            "window_size": window_size,
            "min_spacing": min_spacing,
            "descriptive_features": list(descriptive_features),
            "lower_quantile_cap": lower_quantile_cap,
            "upper_quantile_cap": upper_quantile_cap,
            "log_scale_features": log_scale_features,
            "outfolder": state_outdir,
            "scale_features": True,
            "incomplete_window_policy": incomplete_window_policy,
            "reuse_prepared_dataset": reuse_prepared_dataset,
            "prepared_dataset_path": prepared_adata_full_path,
            "verbose": verbose,
        }
        adata_prepared = prepare_state_classification_dataset(**prepare_kwargs)

        pre_match, pre_mismatch_reasons = _matches_requested_preprocessing_in_adata(
            adata_prepared,
            features=features,
            binary_features_to_group=binary_features_to_group,
            window_size=window_size,
            min_spacing=min_spacing,
            descriptive_features=descriptive_features,
            incomplete_window_policy=incomplete_window_policy,
            lower_quantile_cap=lower_quantile_cap,
            upper_quantile_cap=upper_quantile_cap,
            log_scale_features=log_scale_features,
            scale_features=True,
        )
        if not pre_match:
            _vinfo(
                verbose,
                "state-clustering",
                "prepared full cache mismatch; recomputing full preprocessing | "
                f"{_short_reasons(pre_mismatch_reasons)}",
            )
            prepare_kwargs["reuse_prepared_dataset"] = False
            adata_prepared = prepare_state_classification_dataset(**prepare_kwargs)
        del df_positions
    else:
        _vinfo(
            verbose,
            "state-clustering",
            "preparing dataset (windowed feature extraction, quantile capping, optional log scaling, scaling)",
        )
    _vinfo(
        verbose,
        "state-clustering",
        (
            f"prepare stage source={'cache' if prepared_cache_used else 'recomputed'} | "
            f"descriptive_features={list(descriptive_features)}"
        ),
    )
    _vdone(verbose, "state-clustering", "prepare/reuse full dataset", stage_prepare_started)
    
    pre_meta = (
        adata_prepared.uns.get("preprocessing", {})
        if isinstance(getattr(adata_prepared, "uns", {}), dict)
        else {}
    )
    kept_features = pre_meta.get("kept_features", list(adata_prepared.var_names))
    if not isinstance(kept_features, list):
        kept_features = list(adata_prepared.var_names)
    kept_features = [str(c) for c in kept_features if str(c) in adata_prepared.var_names]
    if len(kept_features) == 0:
        kept_features = [str(c) for c in list(adata_prepared.var_names)]

    default_non_feature_cols = [
        "sample_name",
        "TrackID",
        "sub_TrackID",
        "position_t",
        "window_start_position_t",
        "window_end_position_t",
        "window_length_frames",
    ]
    non_feature_cols = pre_meta.get("non_feature_cols", default_non_feature_cols)
    if not isinstance(non_feature_cols, list):
        non_feature_cols = list(default_non_feature_cols)

    binary_cols_to_merge = pre_meta.get("binary_cols_to_merge", None)
    if not isinstance(binary_cols_to_merge, list):
        binary_cols_to_merge = [c for c in list(binary_features_to_group) if c in adata_prepared.obs.columns]
    prepared_preprocessing_meta = (
        dict(adata_prepared.uns.get("preprocessing", {}))
        if isinstance(getattr(adata_prepared, "uns", {}), dict)
        else {}
    )

    if adata_prepared.n_obs < 50:
        raise ValueError("Insufficient rows in full descriptive dataset for clustering.")
    _vinfo(
        verbose,
        "state-clustering",
        "prepared dataset summary | "
        f"rows={adata_prepared.n_obs}, continuous_features={len(kept_features)}, "
        f"binary_group_features={len(binary_cols_to_merge)}",
    )

    if min_spacing is None:
        spacing_to_use = int(window_size)
    else:
        spacing_to_use = int(min_spacing)
    sampling_id_cols = ["sample_name", "TrackID"]
    sampling_time_col = "position_t"
    sampling_start_cutoff_mode = "time_from_track_start"
    sampling_start_cutoff_window_size = None if window_size is None else int(window_size)
    sampling_min_time_from_track_start = (
        None
        if sampling_start_cutoff_window_size is None
        else max(int(sampling_start_cutoff_window_size) - 1, 0)
    )
    nan_policy = "drop_any_nan_in_kept_features"

    stage_sampling_started = _vstart(
        verbose,
        "state-clustering",
        f"sampling + NaN cleanup | min_spacing={spacing_to_use} | max_samples={max_samples}",
    )
    model_adata = None
    model_cache_loaded = False
    if (
        bool(reuse_prepared_dataset)
        and prepared_adata_model_path is not None
        and prepared_adata_model_path.exists()
    ):
        try:
            cached_model_adata = sc.read_h5ad(prepared_adata_model_path)
            pre_match_model, pre_reasons_model = _matches_requested_preprocessing_in_adata(
                cached_model_adata,
                features=features,
                binary_features_to_group=binary_features_to_group,
                window_size=window_size,
                min_spacing=min_spacing,
                descriptive_features=descriptive_features,
                incomplete_window_policy=incomplete_window_policy,
                lower_quantile_cap=lower_quantile_cap,
                upper_quantile_cap=upper_quantile_cap,
                log_scale_features=log_scale_features,
                scale_features=True,
            )
            sampling_match, sampling_reasons = _matches_model_sampling_cache(
                cached_model_adata,
                kept_features=kept_features,
                min_spacing_used=spacing_to_use,
                start_cutoff_mode=sampling_start_cutoff_mode,
                start_cutoff_window_size=sampling_start_cutoff_window_size,
                max_samples=max_samples,
                random_state=random_state,
                id_cols=sampling_id_cols,
                time_col=sampling_time_col,
                nan_policy=nan_policy,
            )
            if pre_match_model and sampling_match:
                model_adata = cached_model_adata
                model_cache_loaded = True
                _vinfo(
                    verbose,
                    "state-clustering",
                    f"reusing cached model adata: {prepared_adata_model_path}",
                )
            else:
                mismatch_reasons = []
                if not pre_match_model:
                    mismatch_reasons.extend([f"preprocessing: {r}" for r in pre_reasons_model])
                if not sampling_match:
                    mismatch_reasons.extend([f"sampling: {r}" for r in sampling_reasons])
                _vinfo(
                    verbose,
                    "state-clustering",
                    "cached model adata mismatch; rebuilding sampling/clean stage | "
                    f"{_short_reasons(mismatch_reasons)}",
                )
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"could not load cached model adata; rebuilding ({exc})")

    subsampled_rows_for_log = None
    rows_after_nan_for_pca = None
    if model_adata is None:
        adata_train = subsample_with_temporal_spacing(
            adata_prepared,
            id_cols=sampling_id_cols,
            time_col=sampling_time_col,
            min_spacing=spacing_to_use,
            min_time_from_track_start=sampling_min_time_from_track_start,
            max_samples=max_samples,
            random_state=random_state,
        )
        
        if adata_train.n_obs < 50:
            raise ValueError("Insufficient rows after global subsampling for clustering.")
        subsampled_rows_for_log = int(adata_train.n_obs)

        X_train = _to_numpy_2d(adata_train[:, kept_features].X).astype(float, copy=False)
        valid_mask = ~np.isnan(X_train).any(axis=1)
        n_valid = int(valid_mask.sum())
        if n_valid < 50:
            raise ValueError("Insufficient rows after dropping NaNs for PCA stage.")
        rows_after_nan_for_pca = int(n_valid)
        if n_valid != adata_train.n_obs:
            _vinfo(verbose, "state-clustering", f"dropped NaN rows before PCA: kept {n_valid}/{adata_train.n_obs}")
        model_adata = adata_train[valid_mask, kept_features].copy()
        del X_train, valid_mask, adata_train
    else:
        if model_adata.n_obs < 50:
            raise ValueError("Cached model dataset has insufficient rows for clustering.")
        subsampled_rows_for_log = int(model_adata.n_obs)
        rows_after_nan_for_pca = int(model_adata.n_obs)

    _vinfo(
        verbose,
        "state-clustering",
        f"model stage source={'cache' if model_cache_loaded else 'recomputed'}",
    )

    ncomps_requested = 0
    if use_pca:
        ncomps_requested = min(len(kept_features), model_adata.n_obs - 1)
        if ncomps_requested < 2:
            raise ValueError("Insufficient rows/features to run PCA stage.")

    _vdone(verbose, "state-clustering", "sampling + NaN cleanup", stage_sampling_started)
    _vinfo(verbose, "state-clustering", f"subsampled rows={subsampled_rows_for_log} (spacing={spacing_to_use})")

    pca_input_shape = (int(model_adata.n_obs), int(len(kept_features)))
    pca_recomputed = False
    neighbors_recomputed = False
    if use_pca:
        stage_pca_started = _vstart(
            verbose,
            "state-clustering",
            (
                "PCA stage | "
                f"matrix_shape={pca_input_shape}, kept_features={len(kept_features)}, "
                f"rows_after_nan={rows_after_nan_for_pca}, min_var_selection={pca_var_selection}"
            ),
        )
        if model_cache_loaded:
            pca_match, pca_reasons = _matches_cached_pca_stage(
                model_adata,
                use_pca=use_pca,
                pca_var_selection=pca_var_selection,
                svd_solver="full",
                random_state=random_state,
                ncomps_requested=ncomps_requested,
            )
            if pca_match:
                _vinfo(verbose, "state-clustering", "reusing cached PCA stage")
            else:
                _vinfo(
                    verbose,
                    "state-clustering",
                    f"PCA recompute required | reasons={_short_reasons(pca_reasons)}",
                )
                model_adata = run_pca(
                    model_adata,
                    pca_var_selection=pca_var_selection,
                    ncomps=ncomps_requested,
                    svd_solver="full",
                    random_state=random_state,
                )
                pca_recomputed = True
        else:
            model_adata = run_pca(
                model_adata,
                pca_var_selection=pca_var_selection,
                ncomps=ncomps_requested,
                svd_solver="full",
                random_state=random_state,
            )
            pca_recomputed = True
        _vdone(verbose, "state-clustering", "PCA stage", stage_pca_started)
    else:
        stage_pca_started = _vstart(
            verbose,
            "state-clustering",
            (
                "PCA stage skipped | "
                f"matrix_shape={pca_input_shape}, kept_features={len(kept_features)}, "
                f"rows_after_nan={rows_after_nan_for_pca}"
            ),
        )
        if model_cache_loaded:
            pca_match, pca_reasons = _matches_cached_pca_stage(
                model_adata,
                use_pca=use_pca,
                pca_var_selection=pca_var_selection,
                svd_solver="full",
                random_state=random_state,
                ncomps_requested=ncomps_requested,
            )
            if pca_match:
                _vinfo(verbose, "state-clustering", "reusing cached no-PCA representation")
            else:
                _vinfo(
                    verbose,
                    "state-clustering",
                    f"no-PCA mode requested; clearing cached PCA artifacts | reasons={_short_reasons(pca_reasons)}",
                )
                if "X_pca" in model_adata.obsm:
                    del model_adata.obsm["X_pca"]
                if "pca" in model_adata.uns:
                    del model_adata.uns["pca"]
                if "PCs" in model_adata.varm:
                    del model_adata.varm["PCs"]
                pca_recomputed = True
        else:
            _vinfo(verbose, "state-clustering", "PCA disabled; using scaled feature matrix directly")
        _vdone(verbose, "state-clustering", "PCA stage (skipped)", stage_pca_started)

    stage_neighbors_started = _vstart(
        verbose,
        "state-clustering",
        f"neighbors stage | n_neighbors={n_neighbors} | use_rep={representation_label}",
    )
    if model_cache_loaded:
        neighbors_match, neighbors_reasons = _matches_cached_neighbors_stage(
            model_adata,
            n_neighbors=n_neighbors,
            random_state=random_state,
            use_rep=representation_label,
        )
        if pca_recomputed:
            neighbors_match = False
            neighbors_reasons = list(neighbors_reasons) + ["upstream PCA recomputed"]
        if neighbors_match:
            _vinfo(
                verbose,
                "state-clustering",
                f"reusing cached neighbors graph (n_neighbors={n_neighbors}, use_rep={representation_label})",
            )
        else:
            _vinfo(
                verbose,
                "state-clustering",
                f"neighbors recompute required | reasons={_short_reasons(neighbors_reasons)}",
            )
            sc.pp.neighbors(model_adata, n_neighbors=n_neighbors, use_rep=scanpy_use_rep, random_state=random_state)
            neighbors_recomputed = True
    else:
        sc.pp.neighbors(model_adata, n_neighbors=n_neighbors, use_rep=scanpy_use_rep, random_state=random_state)
        neighbors_recomputed = True
    _vdone(verbose, "state-clustering", "neighbors stage", stage_neighbors_started)

    stage_umap_started = _vstart(
        verbose,
        "state-clustering",
        f"UMAP stage | min_dist={min_dist}, random_state={random_state}",
    )
    if model_cache_loaded:
        umap_match, umap_reasons = _matches_cached_umap_stage(
            model_adata,
            min_dist=min_dist,
            random_state=random_state,
            use_rep=representation_label,
        )
        if pca_recomputed or neighbors_recomputed:
            umap_match = False
            umap_reasons = list(umap_reasons) + ["upstream representation/graph recomputed"]
        if umap_match:
            _vinfo(verbose, "state-clustering", f"reusing cached UMAP embedding (min_dist={min_dist})")
        else:
            _vinfo(
                verbose,
                "state-clustering",
                f"UMAP recompute required | reasons={_short_reasons(umap_reasons)}",
            )
            sc.tl.umap(
                model_adata, 
                min_dist=min_dist, 
                random_state=random_state,
                n_components=2
                )
    else:
        sc.tl.umap(
            model_adata, 
            min_dist=min_dist, 
            random_state=random_state,
            n_components=2
            )
    _vdone(verbose, "state-clustering", "UMAP stage", stage_umap_started)

    stage_clustering_started = _vstart(
        verbose,
        "state-clustering",
        (
            f"{clustering_method} clustering | "
            f"n_neighbors={n_neighbors}, min_dist={min_dist}, resolution={resolution}"
        ),
    )
    if isinstance(resolution, (list, tuple, np.ndarray, pd.Series)):
        stability_resolutions = [float(r) for r in list(resolution)]
    else:
        stability_resolutions = [float(resolution)]
    method = str(clustering_method).lower()
    if method == "leiden":
        if ("neighbors" not in model_adata.uns) or ("connectivities" not in model_adata.obsp):
            raise ValueError(
                "Leiden clustering requires a precomputed neighbors graph when n_neighbors=None. "
                "Expected adata.uns['neighbors'] and adata.obsp['connectivities']."
            )
        model_adata = run_leiden_clustering(
            model_adata,
            n_neighbors=None,
            resolution=resolution,
            use_rep=scanpy_use_rep,
            metric="euclidean",
            method="umap",
            key_added="intrinsic_behavioral_cluster",
            random_state=random_state,
            stability_resolutions=stability_resolutions,
        )
        model_adata.obs["intrinsic_behavioral_cluster"] = model_adata.obs["intrinsic_behavioral_cluster"].astype("category")
    elif method == "kmeans":
        k = int(resolution) if isinstance(resolution, (int, float)) and float(resolution) >= 2 else 5
        kmeans_input = model_adata.obsm["X_pca"] if use_pca else _to_numpy_2d(model_adata.X).astype(float, copy=False)
        labels = KMeans(n_clusters=k, random_state=random_state, n_init="auto").fit_predict(kmeans_input)
        model_adata.obs["intrinsic_behavioral_cluster"] = pd.Categorical(
            pd.Series(labels.astype(str), index=model_adata.obs.index)
        )
    else:
        raise ValueError("clustering_method must be one of: 'leiden', 'kmeans'.")

    binary_group_constraints = _infer_binary_group_constraints(
        model_adata.obs,
        binary_cols_to_merge,
    )
    _vinfo(
        verbose,
        "state-clustering",
        (
            "binary-group constraints inferred | "
            f"allowed_groups={len(binary_group_constraints.get('allowed_binary_groups', []))}, "
            f"forbidden_combinations={len(binary_group_constraints.get('forbidden_binary_combinations', []))}"
        ),
    )
    _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=True,
    )
    _vdone(verbose, "state-clustering", "clustering stage", stage_clustering_started)

    stage_prepared_write_started = time.perf_counter()
    prepared_adata_full = adata_prepared[:, kept_features].copy()
    prepared_adata_full.uns["preprocessing"] = dict(prepared_preprocessing_meta)
    prepared_adata_full.uns["preprocessing"].pop("prepared_adata_full_signature", None)
    prepared_adata_full.uns["preprocessing"]["prepared_adata_full_path"] = str(prepared_adata_full_path)
    prepared_adata_full.write(prepared_adata_full_path, compression="gzip")
    del prepared_adata_full, adata_prepared
    _vdone(verbose, "state-clustering", "write prepared full adata cache", stage_prepared_write_started)

    n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
    n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique())
    _vinfo(
        verbose,
        "state-clustering",
        "stage-1 clustering summary | "
        f"intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}, model_rows={model_adata.n_obs}",
    )
    _vinfo(
        verbose,
        "state-clustering",
        (
            "final feature summary | "
            f"descriptive_features={list(descriptive_features)}, "
            f"continuous_features={len(kept_features)}"
        ),
    )

    diagnostics_csvs = {}
    diagnostics_pdf = None
    feature_distribution_pdf = None
    if state_clustering_outdir is not None:
        diagnostics_pdf = state_clustering_outdir / "behavioral_clustering_diagnostics.pdf"
        plot_clustering_diagnostics_pdf(
            adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            feature_cols=kept_features,
            pdf_path=diagnostics_pdf,
            title=f"all_data | {clustering_method} (resolution={resolution})",
        )
        diagnostics_csvs = _export_clustering_diagnostics_csvs(
            adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            feature_cols=kept_features,
            outdir=state_clustering_outdir,
            filename_prefix="behavioral_clustering_diagnostics",
        )
        feature_distribution_pdf = state_clustering_outdir / "behavioral_clustering_feature_distributions.pdf"
        plot_feature_distribution_pdf(
            adata=model_adata,
            feature_cols=kept_features,
            pdf_path=feature_distribution_pdf,
            preprocessing_params=prepared_preprocessing_meta,
            title=f"all_data | {clustering_method} (resolution={resolution})",
        )
        _vsave(
            verbose,
            "state-clustering",
            "feature distribution PDF",
            feature_distribution_pdf,
        )
        if ("binary_group" in model_adata.obs.columns) and ("behavioral_clusterid" in model_adata.obs.columns):
            pdf1 = state_clustering_outdir / "behavioral_clustering_binary_group_cluster_proportions.pdf"
            fig1, _ = plot_binary_group_behavioral_cluster_grid(
                model_adata,
                pdf_path=pdf1,
                return_csv=True,
                verbose=verbose,
            )
            plt.close(fig1)
            pdf2 = state_clustering_outdir / "behavioral_clustering_cluster_binary_group_proportions.pdf"
            fig2, _ = plot_behavioral_cluster_binary_group_grid(
                model_adata,
                pdf_path=pdf2,
                return_csv=True,
            )
            plt.close(fig2)

    ncomps_effective = (
        int(model_adata.obsm["X_pca"].shape[1])
        if ("X_pca" in model_adata.obsm and model_adata.obsm["X_pca"] is not None)
        else 0
    )
    model_cache_meta = {
        "sampling": {
            "min_spacing_used": int(spacing_to_use),
            "max_samples": None if max_samples is None else int(max_samples),
            "random_state": int(random_state),
            "time_col": str(sampling_time_col),
            "id_cols": [str(c) for c in sampling_id_cols],
            "start_cutoff_mode": str(sampling_start_cutoff_mode),
            "start_cutoff_window_size": (
                None
                if sampling_start_cutoff_window_size is None
                else int(sampling_start_cutoff_window_size)
            ),
        },
        "cleaning": {
            "nan_policy": str(nan_policy),
        },
        "features": {
            "kept_features": [str(c) for c in kept_features],
        },
        "pca": {
            "enabled": bool(use_pca),
            "pca_var_selection": (float(pca_var_selection) if use_pca else None),
            "svd_solver": ("full" if use_pca else None),
            "random_state": int(random_state),
            "ncomps_requested": int(ncomps_requested),
            "ncomps_effective": int(ncomps_effective),
        },
        "neighbors": {
            "n_neighbors": int(n_neighbors),
            "use_rep": str(representation_label),
            "random_state": int(random_state),
        },
        "umap": {
            "min_dist": float(min_dist),
            "random_state": int(random_state),
            "use_rep": str(representation_label),
        },
        "clustering": {
            "use_rep": str(representation_label),
        },
    }
    model_preprocessing_meta = dict(prepared_preprocessing_meta)
    model_preprocessing_meta["model_cache"] = dict(model_cache_meta)
    model_preprocessing_meta.pop("prepared_adata_full_signature", None)
    model_preprocessing_meta["prepared_adata_full_path"] = (
        None if prepared_adata_full_path is None else str(prepared_adata_full_path)
    )
    model_adata.uns["preprocessing"] = model_preprocessing_meta

    model_adata.uns["clustering"] = {
        "clustering_method": clustering_method,
        "use_pca": bool(use_pca),
        "use_rep": str(representation_label),
        "resolution": resolution,
        "n_neighbors": int(n_neighbors),
        "random_state": int(random_state),
        "non_feature_cols": list(non_feature_cols),
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "binary_group_constraints": copy.deepcopy(dict(binary_group_constraints)),
        "diagnostics_pdf": None if diagnostics_pdf is None else str(diagnostics_pdf),
        "diagnostics_csvs": dict(diagnostics_csvs),
        "feature_distribution_pdf": (
            None if feature_distribution_pdf is None else str(feature_distribution_pdf)
        ),
    }

    exemplar_render_result = {
        "enabled": bool(plot_exemplar_videos),
        "cluster_key": "intrinsic_behavioral_cluster",
        "selection_csv": None,
        "video_paths_by_cluster": {},
        "global_xy_shape": None,
        "cluster_total_counts": {},
        "error": None,
        "render_config": {
            "stage": "after_clustering",
            "cluster_key": "intrinsic_behavioral_cluster",
            "state_key": "intrinsic_behavioral_cluster",
            "n_per_cluster": int(exemplar_n_per_cluster),
            "selection_random_state": int(random_state),
            "layout_mode": "per_cluster",
            "video_fps": int(exemplar_video_fps),
            "video_dpi": int(exemplar_video_dpi),
            "video_margin": int(exemplar_video_margin),
            "video_pmin": float(exemplar_video_pmin),
            "video_pmax": float(exemplar_video_pmax),
            "video_track_color": str(exemplar_video_track_color),
            "video_show_segments": bool(exemplar_video_show_segments),
            "video_segment_style": str(exemplar_video_segment_style),
            "video_segment_color": str(exemplar_video_segment_color),
            "show_state_legend": False,
            "positions_csv_path": None,
        },
    }
    if bool(plot_exemplar_videos):
        exemplar_render_result = _render_state_cluster_exemplar_videos(
            model_adata=model_adata,
            state_paths=state_paths,
            cell_type=cell_type,
            state_clustering_outdir=state_clustering_outdir,
            cluster_key="intrinsic_behavioral_cluster",
            n_per_cluster=int(exemplar_n_per_cluster),
            random_state=int(random_state),
            video_fps=int(exemplar_video_fps),
            video_dpi=int(exemplar_video_dpi),
            video_margin=int(exemplar_video_margin),
            video_pmin=float(exemplar_video_pmin),
            video_pmax=float(exemplar_video_pmax),
            video_track_color=str(exemplar_video_track_color),
            video_show_segments=bool(exemplar_video_show_segments),
            video_segment_style=str(exemplar_video_segment_style),
            video_segment_color=str(exemplar_video_segment_color),
            verbose=verbose,
        )

    model_adata.uns.setdefault("visualization", {})
    model_adata.uns["visualization"].update(
        {
            "exemplar_render_stage": "after_clustering",
            "exemplar_cluster_key": str(exemplar_render_result.get("cluster_key", "intrinsic_behavioral_cluster")),
            "exemplar_selection_csv": exemplar_render_result.get("selection_csv", None),
            "exemplar_video_paths_by_cluster": dict(
                exemplar_render_result.get("video_paths_by_cluster", {})
            ),
            "exemplar_backprojection_video_paths_by_cluster": dict(
                exemplar_render_result.get("video_paths_by_cluster", {})
            ),
            "exemplar_video_canvas_shape": exemplar_render_result.get(
                "global_xy_shape", None
            ),
            "exemplar_cluster_total_counts": dict(
                exemplar_render_result.get("cluster_total_counts", {})
            ),
            "exemplar_render_config": dict(
                exemplar_render_result.get("render_config", {})
            ),
            "exemplar_segment_outline_errors": dict(
                exemplar_render_result.get("segment_outline_errors", {})
            ),
            "exemplar_render_error": exemplar_render_result.get("error", None),
        }
    )

    stage_model_write_started = time.perf_counter()
    model_adata.write(prepared_adata_model_path, compression="gzip")
    _vdone(verbose, "state-clustering", "write model adata cache", stage_model_write_started)
    _vsave(verbose, "state-clustering", "model-stage cache", prepared_adata_model_path)
    _vsave(verbose, "state-clustering", "prepared full adata cache", prepared_adata_full_path)
    _vdone(verbose, "state-clustering", "run state clustering", run_started)
    return model_adata


def train_state_classifiers(
    output_dir,
    cell_type,
    model_adata=None,
    label_transfer_method="classifier",
    classifier_backend="random_forest",
    classifier_n_estimators=300,
    classifier_min_samples_leaf=2,
    classifier_n_jobs=-1,
    classifier_max_depth=None,
    classifier_min_samples_split=2,
    classifier_max_features="sqrt",
    classifier_class_weight=None,
    classifier_confidence_col="intrinsic_behavioral_cluster_confidence",
    save_label_classifier=True,
    train_continuous_classifier=True,
    train_full_classifier=True,
    full_classifier_confidence_col="full_behavioral_cluster_confidence",
    save_full_label_classifier=True,
    df_prepared=None,
    validation_test_size=0.05,
    validation_random_state=None,
    validation_stratify=True,
    state_paths=None,
    verbose=True,
):
    """
    Stage 2: train intrinsic/full classifiers using model_adata (auto-loaded from canonical path if needed).

    Mutates model_adata in-place by adding/updating model_adata.uns["classification"].

    Returns:
        dict: {
            "full_classifier_path": str | None,
            "partial_classifier_path": str | None,
        }
    """
    train_started = _vstart(verbose, "state-training", "train state classifiers")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    model_adata_path = state_paths.model_adata_path
    if model_adata is None:
        load_started = _vstart(verbose, "state-training", "load model adata")
        if not model_adata_path.exists():
            raise FileNotFoundError(
                "Could not load model_adata for classifier training. "
                f"Expected '{model_adata_path}'. "
                "Run run_state_clustering(...) first or pass model_adata explicitly."
            )
        model_adata = sc.read_h5ad(model_adata_path)
        _vdone(verbose, "state-training", "load model adata", load_started)
        _vinfo(verbose, "state-training", f"Using model adata with rows={model_adata.n_obs}, features={len(model_adata.var_names)}")
    else:
        _vinfo(verbose, "state-training", f"Using provided model adata with rows={model_adata.n_obs}, features={len(model_adata.var_names)}")

    if isinstance(model_adata, dict):
        raise ValueError(
            "train_state_classifiers now expects model_adata (AnnData), not a stage-1 artifact dict."
        )
    if not (hasattr(model_adata, "uns") and hasattr(model_adata, "var_names")):
        raise ValueError("model_adata must be an AnnData-like object with .uns and .var_names.")
    pre_meta = model_adata.uns.get("preprocessing", None)
    clust_meta = model_adata.uns.get("clustering", None)
    if not isinstance(pre_meta, dict):
        raise ValueError(
            "model_adata.uns['preprocessing'] is missing or invalid. "
            "Run run_state_clustering first."
        )
    if not isinstance(clust_meta, dict):
        raise ValueError(
            "model_adata.uns['clustering'] is missing or invalid. "
            "Run run_state_clustering first."
        )

    windowing_params = pre_meta.get("windowing", None)
    if not isinstance(windowing_params, dict):
        raise ValueError("model_adata.uns['preprocessing']['windowing'] must be a dict.")

    transfer_mode = str(label_transfer_method).lower()
    if transfer_mode != "classifier":
        raise ValueError(
            "train_state_classifiers is training-only and supports only "
            "label_transfer_method='classifier'. Use apply_state_classifiers_to_full_dataset for full-data inference."
        )

    classifier_backend_norm = str(classifier_backend).lower()
    if classifier_backend_norm not in {"random_forest", "rf", "randomforest", "randomforestclassifier"}:
        raise ValueError("Only RandomForestClassifier is currently supported (classifier_backend='random_forest').")
    classifier_backend_name = "random_forest"

    cont_cols = list(model_adata.var_names)
    preprocessing_params = _build_preprocessing_params_from_uns(
        pre_meta,
        feature_cols=cont_cols,
    )
    if preprocessing_params is not None and not isinstance(preprocessing_params, dict):
        raise ValueError("model_adata.uns['preprocessing'] could not be converted to preprocessing params.")

    required_clustering_keys = ["non_feature_cols", "binary_cols_to_merge", "resolution", "n_neighbors", "random_state"]
    missing_clustering_keys = [k for k in required_clustering_keys if k not in clust_meta]
    if len(missing_clustering_keys) > 0:
        raise ValueError(
            "model_adata.uns['clustering'] is missing required keys: "
            f"{missing_clustering_keys}"
        )

    non_feature_cols = clust_meta.get("non_feature_cols", None)
    binary_cols_to_merge = clust_meta.get("binary_cols_to_merge", None)
    if not isinstance(non_feature_cols, list):
        raise ValueError("model_adata.uns['clustering']['non_feature_cols'] must be a list.")
    if not isinstance(binary_cols_to_merge, list):
        raise ValueError("model_adata.uns['clustering']['binary_cols_to_merge'] must be a list.")
    non_feature_cols = [str(c) for c in non_feature_cols]
    binary_cols_to_merge = [str(c) for c in binary_cols_to_merge]
    clustering_method = clust_meta.get("clustering_method", None)
    clustering_resolution = clust_meta.get("resolution", None)
    try:
        clustering_n_neighbors = int(clust_meta.get("n_neighbors", 0))
    except Exception as exc:
        raise ValueError(
            "model_adata.uns['clustering']['n_neighbors'] must be int-castable."
        ) from exc
    try:
        stage1_random_state = int(clust_meta.get("random_state", 123))
    except Exception as exc:
        raise ValueError(
            "model_adata.uns['clustering']['random_state'] must be int-castable."
        ) from exc

    validation_mode = _resolve_validation_mode(validation_test_size)
    if validation_mode == "holdout":
        validation_test_size = float(validation_test_size)
    else:
        validation_test_size = 0.0
    validation_random_state = (
        stage1_random_state if validation_random_state is None else int(validation_random_state)
    )
    validation_stratify = bool(validation_stratify)

    required_intrinsic_cols = ["intrinsic_behavioral_cluster"]
    missing_intrinsic_cols = [c for c in required_intrinsic_cols if c not in model_adata.obs.columns]
    if len(missing_intrinsic_cols) > 0:
        raise ValueError(
            "model_adata is missing required intrinsic label columns for training: "
            f"{missing_intrinsic_cols}. Run run_state_clustering before train_state_classifiers."
        )
    if train_full_classifier:
        required_full_cols = ["full_behavioral_cluster"]
        missing_full_cols = [c for c in required_full_cols if c not in model_adata.obs.columns]
        if len(missing_full_cols) > 0:
            raise ValueError(
                "model_adata is missing required full label columns for training: "
                f"{missing_full_cols}. Run run_state_clustering (and optionally "
                "rename_intrinsic_behavioral_clusters) before train_state_classifiers."
            )

    missing_binary_cols = [c for c in binary_cols_to_merge if c not in model_adata.obs.columns]
    if train_full_classifier and len(missing_binary_cols) > 0:
        raise ValueError(
            "model_adata is missing binary columns required for full-classifier training: "
            f"{missing_binary_cols}. Re-run run_state_clustering so binary columns are in model_adata.obs."
        )

    label_classifier = None  # canonical train-split intrinsic model
    label_classifier_artifact_train_split = None
    classifier_fit_info = None  # canonical train-split fit info
    classifier_fit_info_train_split = None
    classifier_qc_model_data = None
    classifier_model_path = None
    classifier_model_paths = {}
    intrinsic_validation = {
        "target_col": "intrinsic_behavioral_cluster",
        "n_train": 0,
        "n_test": 0,
        "train_metrics": None,
        "test_metrics": None,
        "oob_score_train_split": None,
        "oob_score_full_data": None,
        "holdout_artifacts": {},
    }

    outfolder_path = state_paths.processing_outdir
    outfolder_path.mkdir(parents=True, exist_ok=True)
    intrinsic_outfolder = state_paths.intrinsic_outdir
    full_outfolder = state_paths.full_outdir
    intrinsic_qc_outfolder = state_paths.intrinsic_qc_outdir
    full_qc_outfolder = state_paths.full_qc_outdir
    intrinsic_outfolder.mkdir(parents=True, exist_ok=True)
    full_outfolder.mkdir(parents=True, exist_ok=True)
    intrinsic_qc_outfolder.mkdir(parents=True, exist_ok=True)
    full_qc_outfolder.mkdir(parents=True, exist_ok=True)

    if train_continuous_classifier:
        intrinsic_started = _vstart(verbose, "state-training", "train intrinsic classifier")
        intrinsic_label_mask = pd.Series(
            model_adata.obs["intrinsic_behavioral_cluster"],
            index=model_adata.obs.index,
            dtype="string",
        ).notna()
        intrinsic_train_adata = model_adata[intrinsic_label_mask.to_numpy(dtype=bool)].copy()
        X_intrinsic_all = _to_numpy_2d(intrinsic_train_adata[:, cont_cols].X).astype(float, copy=False)
        y_intrinsic_all = intrinsic_train_adata.obs["intrinsic_behavioral_cluster"].astype(str).to_numpy()
        X_intrinsic_train = X_intrinsic_all
        y_intrinsic_train = y_intrinsic_all
        if validation_mode == "holdout":
            if validation_stratify:
                _validate_stratified_split_feasibility(
                    y_intrinsic_all, validation_test_size, target_name="intrinsic_behavioral_cluster"
                )
            idx_all = np.arange(intrinsic_train_adata.n_obs, dtype=int)
            idx_train, idx_test = train_test_split(
                idx_all,
                test_size=validation_test_size,
                random_state=validation_random_state,
                stratify=(y_intrinsic_all if validation_stratify else None),
            )
            X_intrinsic_train = X_intrinsic_all[idx_train]
            y_intrinsic_train = y_intrinsic_all[idx_train]
            X_intrinsic_test = X_intrinsic_all[idx_test]
            y_intrinsic_test = y_intrinsic_all[idx_test]
            intrinsic_validation["n_train"] = int(len(idx_train))
            intrinsic_validation["n_test"] = int(len(idx_test))
        else:
            intrinsic_validation["n_train"] = int(intrinsic_train_adata.n_obs)
            intrinsic_validation["n_test"] = 0

        label_classifier, classifier_fit_info = train_random_forest_classifier(
            X_train=X_intrinsic_train,
            y_train=y_intrinsic_train,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
        )
        classifier_fit_info["cluster_col"] = "intrinsic_behavioral_cluster"
        classifier_fit_info_train_split = classifier_fit_info

        label_classifier_artifact_train_split = {
            "classifier": label_classifier,
            "backend": classifier_backend_name,
            "feature_cols": cont_cols,
            "cluster_col": "intrinsic_behavioral_cluster",
            "preprocessing": preprocessing_params,
            "windowing": windowing_params,
        }

        y_pred_intrinsic_train = np.asarray(label_classifier.predict(X_intrinsic_train)).astype(str)
        intrinsic_validation["train_metrics"] = _build_compact_metric_summary(
            y_intrinsic_train, y_pred_intrinsic_train
        )
        intrinsic_validation["oob_score_train_split"] = float(
            classifier_fit_info.get("oob_score", np.nan)
        )

        if validation_mode == "holdout":
            y_pred_intrinsic_test = np.asarray(label_classifier.predict(X_intrinsic_test)).astype(str)
            holdout_intrinsic_eval = _evaluate_holdout_predictions(
                y_true=y_intrinsic_test,
                y_pred=y_pred_intrinsic_test,
                outfolder=intrinsic_qc_outfolder,
                filename_prefix=f"state_classification_intrinsic_classifier_holdout_qc_{classifier_backend_name}",
            )
            intrinsic_validation["test_metrics"] = holdout_intrinsic_eval["metrics"]
            intrinsic_validation["holdout_artifacts"] = holdout_intrinsic_eval["artifacts"]

        classifier_qc_model_data = evaluate_clusterid_classifier_on_model_adata(
            model_adata=intrinsic_train_adata,
            classifier=label_classifier,
            cluster_col="intrinsic_behavioral_cluster",
            outfolder=intrinsic_qc_outfolder,
            filename_prefix=f"state_classification_classifier_qc_{classifier_backend_name}",
            verbose=verbose,
        )
        intrinsic_validation["oob_score_full_data"] = (
            intrinsic_validation["oob_score_train_split"] if validation_mode != "holdout" else None
        )

        if validation_mode == "holdout":
            tm = intrinsic_validation["train_metrics"] or {}
            vm = intrinsic_validation["test_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic validation_set scores: "
                    f"n_val={intrinsic_validation['n_test']}, "
                    f"val_acc={vm.get('accuracy', np.nan):.4f}, "
                    f"val_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                    f"val_macro_f1={vm.get('macro_f1', np.nan):.4f}"
                ),
            )
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic training/OOB scores: "
                    f"n_train={intrinsic_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={intrinsic_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        else:
            tm = intrinsic_validation["train_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic training/OOB scores: "
                    f"n_train={intrinsic_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={intrinsic_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        _vdone(verbose, "state-training", "train intrinsic classifier", intrinsic_started)
        del X_intrinsic_all, y_intrinsic_all, X_intrinsic_train, y_intrinsic_train
        if validation_mode == "holdout":
            del X_intrinsic_test, y_intrinsic_test

    full_label_classifier_artifact = None  # canonical train-split full-label artifact
    full_label_classifier_artifact_train_split = None
    full_classifier_fit_info = None  # canonical train-split full fit info
    full_classifier_fit_info_train_split = None
    full_classifier_qc_model_data = None
    full_classifier_model_path = None
    full_classifier_model_paths = {}
    full_validation = {
        "target_col": "full_behavioral_cluster",
        "n_train": 0,
        "n_test": 0,
        "train_metrics": None,
        "test_metrics": None,
        "oob_score_train_split": None,
        "oob_score_full_data": None,
        "holdout_artifacts": {},
    }

    if train_full_classifier:
        full_started = _vstart(verbose, "state-training", "train full classifier")
        full_label_mask = pd.Series(
            model_adata.obs["full_behavioral_cluster"],
            index=model_adata.obs.index,
            dtype="string",
        ).notna()
        full_train_adata = model_adata[full_label_mask.to_numpy(dtype=bool)].copy()
        y_full_all = full_train_adata.obs["full_behavioral_cluster"].astype(str).to_numpy()

        if validation_mode == "holdout":
            if validation_stratify:
                _validate_stratified_split_feasibility(
                    y_full_all, validation_test_size, target_name="full_behavioral_cluster"
                )
            idx_all = np.arange(full_train_adata.n_obs, dtype=int)
            idx_train, idx_test = train_test_split(
                idx_all,
                test_size=validation_test_size,
                random_state=validation_random_state,
                stratify=(y_full_all if validation_stratify else None),
            )

            X_full_train, expanded_binary_train_cols = _build_classifier_matrix_from_adata(
                adata=full_train_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                return_binary_feature_names=True,
                row_indices=idx_train,
            )
            X_full_test = _build_classifier_matrix_from_adata(
                adata=full_train_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                binary_expanded_feature_cols=expanded_binary_train_cols,
                return_binary_feature_names=False,
                row_indices=idx_test,
            )
            y_full_train = y_full_all[idx_train]
            y_full_test = y_full_all[idx_test]
            full_validation["n_train"] = int(len(idx_train))
            full_validation["n_test"] = int(len(idx_test))
        else:
            X_full_train, expanded_binary_train_cols = _build_classifier_matrix_from_adata(
                adata=full_train_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                return_binary_feature_names=True,
            )
            y_full_train = y_full_all
            full_validation["n_train"] = int(full_train_adata.n_obs)
            full_validation["n_test"] = 0

        clf_full_train_split, fit_info_full_train_split = train_random_forest_classifier(
            X_train=X_full_train,
            y_train=y_full_train,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
        )
        fit_info_full_train_split["target_col"] = "full_behavioral_cluster"
        fit_info_full_train_split["continuous_feature_cols"] = list(cont_cols)
        fit_info_full_train_split["binary_feature_cols"] = list(binary_cols_to_merge)
        fit_info_full_train_split["n_features_total"] = int(X_full_train.shape[1])
        fit_info_full_train_split["n_features_continuous"] = int(len(cont_cols))
        fit_info_full_train_split["n_features_binary"] = int(len(binary_cols_to_merge))
        fit_info_full_train_split["n_features_binary_expanded"] = int(len(expanded_binary_train_cols))
        full_classifier_fit_info_train_split = fit_info_full_train_split
        full_classifier_fit_info = fit_info_full_train_split

        full_label_classifier_artifact_train_split = {
            "classifier": clf_full_train_split,
            "continuous_feature_cols": list(cont_cols),
            "binary_feature_cols": list(binary_cols_to_merge),
            "binary_expanded_feature_cols": list(expanded_binary_train_cols),
            "target_col": "full_behavioral_cluster",
            "preprocessing": preprocessing_params,
            "windowing": dict(windowing_params),
        }
        full_label_classifier_artifact = full_label_classifier_artifact_train_split

        y_pred_full_train = np.asarray(clf_full_train_split.predict(X_full_train)).astype(str)
        full_validation["train_metrics"] = _build_compact_metric_summary(y_full_train, y_pred_full_train)
        full_validation["oob_score_train_split"] = float(
            fit_info_full_train_split.get("oob_score", np.nan)
        )

        if validation_mode == "holdout":
            y_pred_full_test = np.asarray(clf_full_train_split.predict(X_full_test)).astype(str)
            holdout_full_eval = _evaluate_holdout_predictions(
                y_true=y_full_test,
                y_pred=y_pred_full_test,
                outfolder=full_qc_outfolder,
                filename_prefix=f"state_classification_full_classifier_holdout_qc_{classifier_backend_name}",
            )
            full_validation["test_metrics"] = holdout_full_eval["metrics"]
            full_validation["holdout_artifacts"] = holdout_full_eval["artifacts"]

        if validation_mode == "holdout":
            X_full_all = _build_classifier_matrix_from_adata(
                adata=full_train_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                binary_expanded_feature_cols=full_label_classifier_artifact_train_split.get(
                    "binary_expanded_feature_cols", []
                ),
                return_binary_feature_names=False,
            )
        else:
            X_full_all = X_full_train
        y_pred_full_all = np.asarray(
            full_label_classifier_artifact_train_split["classifier"].predict(X_full_all)
        ).astype(str)
        full_classifier_qc_model_data = evaluate_classifier_predictions(
            y_true=y_full_all,
            y_pred=y_pred_full_all,
            outfolder=full_qc_outfolder,
            filename_prefix="state_classification_full_classifier_qc_random_forest",
        )
        full_validation["oob_score_full_data"] = (
            full_validation["oob_score_train_split"] if validation_mode != "holdout" else None
        )
        _vinfo(
            verbose,
            "state-training",
            (
                "Full model-data QC: "
                f"n={full_classifier_qc_model_data['n_rows']}, "
                f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
                f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
            ),
        )
        if validation_mode == "holdout":
            tm = full_validation["train_metrics"] or {}
            vm = full_validation["test_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full validation_set scores: "
                    f"n_val={full_validation.get('n_test', 0)}, "
                    f"val_acc={vm.get('accuracy', np.nan):.4f}, "
                    f"val_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                    f"val_macro_f1={vm.get('macro_f1', np.nan):.4f}"
                ),
            )
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full training/OOB scores: "
                    f"n_train={full_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={full_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        else:
            tm = full_validation["train_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full training/OOB scores: "
                    f"n_train={full_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={full_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        _vdone(verbose, "state-training", "train full classifier", full_started)
        del X_full_train, y_full_train, X_full_all, y_pred_full_train, y_pred_full_all
        if validation_mode == "holdout":
            del X_full_test, y_full_test, y_pred_full_test

    intrinsic_metrics_csv = None
    full_metrics_csv = None
    if intrinsic_qc_outfolder is not None and train_continuous_classifier:
        intrinsic_tm = intrinsic_validation.get("train_metrics") or {}
        intrinsic_vm = intrinsic_validation.get("test_metrics") or {}
        intrinsic_metrics_row = {
            "mode": validation_mode,
            "n_train": intrinsic_validation.get("n_train", 0),
            "n_test": intrinsic_validation.get("n_test", 0),
            "oob_score_train_split": intrinsic_validation.get("oob_score_train_split", np.nan),
            "train_accuracy": intrinsic_tm.get("accuracy", np.nan),
            "train_balanced_accuracy": intrinsic_tm.get("balanced_accuracy", np.nan),
            "train_macro_f1": intrinsic_tm.get("macro_f1", np.nan),
            "test_accuracy": intrinsic_vm.get("accuracy", np.nan),
            "test_balanced_accuracy": intrinsic_vm.get("balanced_accuracy", np.nan),
            "test_macro_f1": intrinsic_vm.get("macro_f1", np.nan),
            "model_data_accuracy": (
                None if classifier_qc_model_data is None else classifier_qc_model_data.get("accuracy", np.nan)
            ),
            "model_data_balanced_accuracy": (
                None if classifier_qc_model_data is None else classifier_qc_model_data.get("balanced_accuracy", np.nan)
            ),
        }
        intrinsic_metrics_csv = (
            intrinsic_qc_outfolder / f"state_classification_intrinsic_classifier_metrics_{classifier_backend_name}.csv"
        )
        pd.DataFrame([intrinsic_metrics_row]).to_csv(intrinsic_metrics_csv, index=False)
        _vsave(verbose, "state-training", "intrinsic metrics csv", intrinsic_metrics_csv)

    if full_qc_outfolder is not None and train_full_classifier:
        full_tm = full_validation.get("train_metrics") or {}
        full_vm = full_validation.get("test_metrics") or {}
        full_metrics_row = {
            "mode": validation_mode,
            "n_train": full_validation.get("n_train", 0),
            "n_test": full_validation.get("n_test", 0),
            "oob_score_train_split": full_validation.get("oob_score_train_split", np.nan),
            "train_accuracy": full_tm.get("accuracy", np.nan),
            "train_balanced_accuracy": full_tm.get("balanced_accuracy", np.nan),
            "train_macro_f1": full_tm.get("macro_f1", np.nan),
            "test_accuracy": full_vm.get("accuracy", np.nan),
            "test_balanced_accuracy": full_vm.get("balanced_accuracy", np.nan),
            "test_macro_f1": full_vm.get("macro_f1", np.nan),
            "model_data_accuracy": (
                None if full_classifier_qc_model_data is None else full_classifier_qc_model_data.get("accuracy", np.nan)
            ),
            "model_data_balanced_accuracy": (
                None if full_classifier_qc_model_data is None else full_classifier_qc_model_data.get("balanced_accuracy", np.nan)
            ),
        }
        full_metrics_csv = (
            full_qc_outfolder / f"state_classification_full_classifier_metrics_{classifier_backend_name}.csv"
        )
        pd.DataFrame([full_metrics_row]).to_csv(full_metrics_csv, index=False)
        _vsave(verbose, "state-training", "full metrics csv", full_metrics_csv)

    model_adata.uns["classification"] = {
        "label_transfer_method": transfer_mode,
        "classifier_backend": classifier_backend_name,
        "classifier_fit_info": classifier_fit_info,
        "classifier_fit_info_train_split": classifier_fit_info_train_split,
        "classifier_qc_model_data": classifier_qc_model_data,
        "classifier_confidence_col": classifier_confidence_col,
        "classifier_model_path": None,
        "classifier_model_paths": {},
        "full_classifier_fit_info": None,
        "full_classifier_fit_info_train_split": None,
        "full_classifier_qc_model_data": None,
        "full_classifier_prediction_col": "full_behavioral_cluster",
        "full_classifier_confidence_col": full_classifier_confidence_col,
        "full_classifier_model_path": None,
        "full_classifier_model_paths": {},
        "model_vs_full_clusterid_check": None,
        "intrinsic_artifacts_dir": None if intrinsic_outfolder is None else str(intrinsic_outfolder),
        "full_artifacts_dir": None if full_outfolder is None else str(full_outfolder),
        "intrinsic_qc_dir": None if intrinsic_qc_outfolder is None else str(intrinsic_qc_outfolder),
        "full_qc_dir": None if full_qc_outfolder is None else str(full_qc_outfolder),
        "intrinsic_metrics_csv": None if intrinsic_metrics_csv is None else str(intrinsic_metrics_csv),
        "full_metrics_csv": None if full_metrics_csv is None else str(full_metrics_csv),
        "validation_config": {
            "mode": validation_mode,
            "test_size": None if validation_mode == "oob_only" else float(validation_test_size),
            "random_state": int(validation_random_state),
            "stratify": bool(validation_stratify),
            "split_unit": "row",
        },
        "intrinsic_validation": intrinsic_validation,
        "full_validation": full_validation,
    }

    if train_full_classifier:
        model_adata.uns["classification"]["full_classifier_fit_info"] = full_classifier_fit_info
        model_adata.uns["classification"]["full_classifier_fit_info_train_split"] = (
            full_classifier_fit_info_train_split
        )
        model_adata.uns["classification"]["full_classifier_qc_model_data"] = full_classifier_qc_model_data

    def _build_pipeline_metadata_payload():
        return {
            "preprocessing": copy.deepcopy(dict(model_adata.uns.get("preprocessing", {}))),
            "clustering": copy.deepcopy(dict(model_adata.uns.get("clustering", {}))),
            "classification": copy.deepcopy(dict(model_adata.uns.get("classification", {}))),
        }

    if label_classifier_artifact_train_split is not None:
        label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
    if full_label_classifier_artifact_train_split is not None:
        full_label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()

    if save_label_classifier and label_classifier_artifact_train_split is not None:
        classifier_model_path = state_paths.intrinsic_classifier_default_path
        classifier_model_path.parent.mkdir(parents=True, exist_ok=True)
        model_adata.uns["classification"]["classifier_model_path"] = str(classifier_model_path)
        classifier_model_paths["train_split"] = str(classifier_model_path)
        model_adata.uns["classification"]["classifier_model_paths"] = dict(classifier_model_paths)
        label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
        with open(classifier_model_path, "wb") as f:
            pickle.dump(label_classifier_artifact_train_split, f)
        _vsave(verbose, "state-training", "intrinsic classifier artifact", classifier_model_path)

    if train_full_classifier and save_full_label_classifier and full_label_classifier_artifact_train_split is not None:
        full_classifier_model_path = state_paths.full_classifier_default_path
        full_classifier_model_path.parent.mkdir(parents=True, exist_ok=True)
        model_adata.uns["classification"]["full_classifier_model_path"] = str(full_classifier_model_path)
        full_classifier_model_paths["train_split"] = str(full_classifier_model_path)
        model_adata.uns["classification"]["full_classifier_model_paths"] = dict(full_classifier_model_paths)
        full_label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
        with open(full_classifier_model_path, "wb") as f:
            pickle.dump(full_label_classifier_artifact_train_split, f)
        _vsave(verbose, "state-training", "full classifier artifact", full_classifier_model_path)

    write_started = _vstart(verbose, "state-training", "write model adata")
    model_adata.write(state_paths.model_adata_path, compression="gzip")
    _vdone(verbose, "state-training", "write model adata", write_started)
    _vsave(verbose, "state-training", "model adata", state_paths.model_adata_path)
    n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
    n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique()) if "full_behavioral_cluster" in model_adata.obs.columns else 0
    _vinfo(
        verbose,
        "state-training",
        f"Training summary: intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}",
    )
    _vdone(verbose, "state-training", "train state classifiers", train_started)
    return {
        "full_classifier_path": None if full_classifier_model_path is None else str(full_classifier_model_path),
        "partial_classifier_path": None if classifier_model_path is None else str(classifier_model_path),
    }


def test_pipeline():
    model_adata = run_state_clustering(
        features=features,
        binary_features_to_group=binary_features_to_group,
        output_dir=output_dir,
        cell_type=cell_type,
        window_size=window_size,
        max_samples=max_samples,
        min_spacing=min_spacing,
        n_neighbors=n_neighbors,
        resolution=resolution,
        descriptive_features=descriptive_features,
        pca_var_selection=pca_var_selection,
        clustering_method=clustering_method,
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        incomplete_window_policy=incomplete_window_policy,
        random_state=random_state,
        reuse_prepared_dataset=reuse_prepared_dataset,
        verbose=verbose,
    )
    
    model_adata.obs["intrinsic_behavioral_cluster"]
    model_adata.obs["full_behavioral_cluster"]
    
    mapping = {
        1: "dead",
        2: "static",
        3: "plastic_scanner",
        4: "static",
        5: "round_scanner",
    }
    rename_intrinsic_behavioral_clusters(
        adata=model_adata,
        mapping=mapping,
        binary_cols_to_merge=binary_features_to_group,
    )
    
    build_identity_cluster_mapping(
        model_adata,
        cluster_col="full_behavioral_cluster",
    )
    
    full_mapping = {
        'no_contact_dead': 'dead',
        'no_contact_plastic_scanner': 'plastic_scanner',
        'no_contact_round_scanner': 'round_scanner',
        'no_contact_static': 'static',
        'organoid_contact_pixels_and_tcell_contact_pixels_dead': 'dead',
        'organoid_contact_pixels_and_tcell_contact_pixels_plastic_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_and_tcell_contact_pixels_round_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_and_tcell_contact_pixels_static': 'organoid_contact_pixels',
        'organoid_contact_pixels_dead': 'dead',
        'organoid_contact_pixels_plastic_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_round_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_static': 'organoid_contact_pixels',
        'tcell_contact_pixels_dead': 'dead',
        'tcell_contact_pixels_plastic_scanner': 'tcell_contact_pixels',
        'tcell_contact_pixels_round_scanner': 'tcell_contact_pixels',
        'tcell_contact_pixels_static': 'tcell_contact_pixels'
        }
    
    model_adata = relabel_cluster_ids(
        adata=model_adata,
        mapping=full_mapping,
        cluster_key="full_behavioral_cluster",
        # Set overwrite_original=True to make reruns idempotent.
        overwrite_original=True,
    )
    
    classifier_paths = train_state_classifiers(
        output_dir=output_dir,
        cell_type=cell_type,
        model_adata=model_adata,
        label_transfer_method="classifier",
        classifier_backend=classifier_backend,
        classifier_n_estimators=classifier_n_estimators,
        classifier_min_samples_leaf=classifier_min_samples_leaf,
        classifier_n_jobs=classifier_n_jobs,
        classifier_max_depth=classifier_max_depth,
        classifier_min_samples_split=classifier_min_samples_split,
        classifier_max_features=classifier_max_features,
        classifier_class_weight=classifier_class_weight,
        classifier_confidence_col=classifier_confidence_col,
        save_label_classifier=save_label_classifier,
        train_continuous_classifier=True,
        train_full_classifier=True,
        )
    
    adata_full = apply_state_classifiers_to_full_dataset(
        output_dir=output_dir,
        cell_type=cell_type,
        # label_classifier_artifact=classifier_paths["partial_classifier_path"],
        # full_label_classifier_artifact=classifier_paths["full_classifier_path"],
        continuous_output_col="intrinsic_behavioral_cluster",
        full_output_col="full_behavioral_cluster",
        combine_binary_with_continuous=True,
        verbose=verbose,
    )

    if "sample_name" not in adata_full.obs.columns:
        raise ValueError(
            "test_pipeline requires 'sample_name' in adata_full.obs to visualize backprojection."
        )
    sample_names = (
        pd.Series(adata_full.obs["sample_name"])
        .dropna()
        .astype("string")
        .str.strip()
    )
    sample_names = sample_names[sample_names != ""]
    unique_samples = sorted(sample_names.unique().tolist())
    if len(unique_samples) == 0:
        raise ValueError(
            "test_pipeline could not resolve a sample_name from adata_full.obs for backprojection visualization."
        )
    selected_sample_name = "ROCHE_JM1_Exp042-8_Img03_10T_HER2II"
    _vinfo(verbose, "state-apply", f"Opening backprojection viewer for sample '{selected_sample_name}'")
    _viewer = show_behavioral_state_backprojection(
        sample_name=selected_sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        run=True,
        verbose=verbose,
    )
