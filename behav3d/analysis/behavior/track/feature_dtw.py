import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from dtaidistance import dtw_ndim
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from behav3d.analysis.behavior.general.umap import fit_umap
from behav3d.analysis.behavior.track.visualization.plots.feature_dtw import (
    plot_cluster_percentage_bars,
    plot_clustering_feature_heatmap,
    plot_feature_umap,
)
from behav3d.analysis.behavior.utils import _handle_nan_in_distance_matrix
from behav3d.core.utils import expand_column_patterns, format_time


# ---------------------------------------------------------------------------
# Scaling helpers
# ---------------------------------------------------------------------------

def _minmax_scale(series):
    series = pd.to_numeric(series, errors="coerce").astype("float64")
    min_value = series.min(skipna=True)
    max_value = series.max(skipna=True)
    if pd.isna(min_value) or pd.isna(max_value) or max_value == min_value:
        return pd.Series(0.0, index=series.index)
    return (series - min_value) / (max_value - min_value)


def _original_behav3d_quantile_scale(values):
    values = pd.to_numeric(values, errors="coerce").astype("float64")
    std = values.std()
    if pd.isna(std) or std == 0:
        z_values = pd.Series(0.0, index=values.index)
    else:
        z_values = (values - values.mean()) / std

    q75 = z_values.quantile(0.75)
    min_z = z_values.min(skipna=True)
    if pd.isna(q75) or pd.isna(min_z):
        q_values = pd.Series(0.0, index=values.index)
    else:
        q_values = z_values.where(z_values > q75, min_z)

    q_values = _minmax_scale(q_values)
    max_quantile = q_values.quantile(0.9999999)
    if pd.notna(max_quantile) and max_quantile != 0:
        q_values = q_values / max_quantile

    return q_values


def apply_original_behav3d_feature_scaling(df_tracks, cell_type="tcell"):
    """Add the original BEHAV3D/R-style scaled feature columns used for feature DTW."""
    contact_col = f"{cell_type}_contact"
    required_cols = [
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
        "organoid_contact",
        contact_col,
        "exp_nr",
    ]
    missing_cols = [col for col in required_cols if col not in df_tracks.columns]
    if missing_cols:
        raise ValueError(
            "feature_scaling_preset='original_behav3d' requires missing columns: "
            + ", ".join(missing_cols)
        )

    print("- Applying original BEHAV3D feature scaling preset")
    df_scaled = df_tracks.copy()
    quantile_features = {
        "mean_square_displacement": "q_mean_square_displacement",
        "speed": "q_speed",
        "mean_dead_dye": "q_mean_dead_dye",
    }

    group_key = "exp_nr"
    for source_col, out_col in quantile_features.items():
        df_scaled[out_col] = (
            df_scaled
            .groupby(group_key, group_keys=False)[source_col]
            .transform(_original_behav3d_quantile_scale)
        )

    contact_features = {
        "organoid_contact": "s_organoid_contact",
        contact_col: f"s_{contact_col}",
    }
    for source_col, out_col in contact_features.items():
        df_scaled[out_col] = (
            df_scaled
            .groupby(group_key, group_keys=False)[source_col]
            .transform(_minmax_scale)
        )

    dtw_features = [
        "q_mean_square_displacement",
        "q_speed",
        "q_mean_dead_dye",
        "s_organoid_contact",
        f"s_{contact_col}",
    ]
    return df_scaled, dtw_features


def normalize_track_features(
    df_tracks,
    columns=["mean_square_displacement", "speed", "mean_dead_dye"],
):
    for column_name in columns:
        if column_name not in df_tracks.columns:
            continue
        if df_tracks[column_name].dtype == "object":
            try:
                df_tracks[column_name] = pd.to_numeric(df_tracks[column_name], errors="coerce")
            except Exception:
                print(f"Skipping normalization of non-numeric column: {column_name}")
                continue
        if df_tracks[column_name].std() > 0:
            df_tracks.loc[:, f"z_{column_name}"] = df_tracks[column_name].transform(
                lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
            )
        else:
            print(f"Skipping normalization of zero-variance column: {column_name}")
            df_tracks[f"z_{column_name}"] = 0
    return df_tracks


# ---------------------------------------------------------------------------
# Track filtering & summarization
# ---------------------------------------------------------------------------

def filter_feature_dtw_track_lengths(
    df_tracks,
    min_track_length=None,
    max_track_length=None,
    group_cols=("sample_name", "TrackID"),
):
    """Filter tracks by number of timepoints and trim each track from its start."""
    group_cols = [col for col in group_cols if col in df_tracks.columns]
    if len(group_cols) == 0:
        raise ValueError("Track length filtering requires sample_name and/or TrackID columns.")

    time_cols = [col for col in ["relative_time", "position_t", "time"] if col in df_tracks.columns]
    sort_cols = group_cols + time_cols
    df_filtered = df_tracks.sort_values(sort_cols).reset_index(drop=True).copy()
    before_tracks = df_filtered[group_cols].drop_duplicates().shape[0]

    if min_track_length is not None:
        min_track_length = int(min_track_length)
        track_lengths = df_filtered.groupby(group_cols, observed=True)[group_cols[0]].transform("size")
        df_filtered = df_filtered.loc[track_lengths >= min_track_length].reset_index(drop=True)

    if max_track_length is not None:
        max_track_length = int(max_track_length)
        if max_track_length < 1:
            raise ValueError("max_track_length must be at least 1.")
        df_filtered = (
            df_filtered
            .groupby(group_cols, group_keys=False, observed=True)
            .head(max_track_length)
            .reset_index(drop=True)
        )

    after_tracks = df_filtered[group_cols].drop_duplicates().shape[0]
    print(
        "- Applied Feature DTW track length filtering: "
        f"{before_tracks} → {after_tracks} tracks"
    )
    if len(df_filtered) == 0:
        raise ValueError("No tracks remain after Feature DTW track length filtering.")
    return df_filtered


def summarize_feature_dtw_tracks(df_tracks):
    """Summarize the in-memory tracks used for Feature DTW plotting/clustering."""
    group_cols = ["sample_name", "TrackID"]
    missing = [col for col in group_cols if col not in df_tracks.columns]
    if missing:
        raise ValueError("Cannot summarize Feature DTW tracks; missing columns: " + ", ".join(missing))

    df_tracks = df_tracks.copy()
    for col in df_tracks.columns:
        if col.startswith("touching_"):
            df_tracks[col] = df_tracks[col].astype(str)

    grouped = df_tracks.groupby(group_cols, observed=True)
    df_summarized = grouped.size().reset_index(name="track_length")

    num_bool_cols = df_tracks.select_dtypes(include=[np.number, bool]).columns.drop(
        ["TrackID"], errors="ignore"
    )
    if len(num_bool_cols) > 0:
        means = grouped[num_bool_cols].mean().reset_index()
        means = means.rename(columns={col: f"mean_{col}" for col in num_bool_cols})
        df_summarized = pd.merge(df_summarized, means, on=group_cols, how="left")

    if "dead" in df_tracks.columns:
        df_summarized["dies"] = grouped["dead"].any().reset_index()["dead"]
    if "displacement_from_origin" in df_tracks.columns:
        df_summarized["displacement_from_origin"] = (
            grouped["displacement_from_origin"].max().reset_index()["displacement_from_origin"]
        )
    if "cumulative_displacement" in df_tracks.columns:
        df_summarized["cumulative_displacement"] = (
            grouped["cumulative_displacement"].last().reset_index()["cumulative_displacement"]
        )

    contact_cols = [
        col for col in df_tracks.columns
        if col.endswith("_contact") and not col.startswith("active_")
    ]
    for contact_col in contact_cols:
        active_col = f"active_{contact_col}"
        if active_col in df_tracks.columns:
            def _active_contact_mean(group, contact_col=contact_col, active_col=active_col):
                mask = group[contact_col].fillna(False).astype(bool)
                return group.loc[mask, active_col].mean() if mask.any() else 0

            result = grouped.apply(_active_contact_mean, include_groups=False)
            if isinstance(result, pd.DataFrame):
                result = result.iloc[:, 0]
            df_summarized[active_col] = result.values

    metadata_cols = ["TrackID", "sample_name"]
    for col in ["well", "exp_nr"]:
        if col in df_tracks.columns:
            metadata_cols.append(col)
    metadata_cols.extend([col for col in df_tracks.columns if col.endswith("_line_condition")])
    metadata_cols = list(dict.fromkeys(metadata_cols))

    df_trackinfo = df_tracks[metadata_cols].drop_duplicates()
    return pd.merge(df_trackinfo, df_summarized, how="left")


def _resolve_selected_plot_feature_columns(
    df_tracks_summarized,
    df_tracks,
    selected_features,
    summary_match_mode="all",
):
    """Resolve user-selected feature names to columns available for result PDFs."""
    if selected_features is None:
        selected_features = []
    elif isinstance(selected_features, str):
        selected_features = [selected_features]
    else:
        selected_features = list(selected_features)
    summary_cols = [str(c) for c in list(df_tracks_summarized.columns)]
    summary_set = set(summary_cols)
    track_cols = set(str(c) for c in list(df_tracks.columns))
    out = []
    seen = set()

    def _is_allowed_plot_column(col):
        return "_on_distance" not in str(col)

    def _add(col):
        col = str(col)
        if col in summary_set and col not in seen and _is_allowed_plot_column(col):
            out.append(col)
            seen.add(col)

    def _summary_matches_for_feature(raw_feature):
        token = str(raw_feature).strip()
        if token == "":
            return []
        escaped = re.escape(token)
        pattern = re.compile(rf"(^|_){escaped}($|_)")
        matched = [
            col for col in summary_cols
            if _is_allowed_plot_column(col) and bool(pattern.search(str(col)))
        ]
        if summary_match_mode == "all":
            return matched
        return matched[:1]

    for feature in selected_features:
        feature = str(feature)
        expanded = expand_column_patterns(feature, df_tracks_summarized.columns)
        exact_summary_matches = [
            str(col) for col in expanded
            if str(col) in summary_set and _is_allowed_plot_column(col)
        ]
        for col in exact_summary_matches:
            _add(col)

        base_feature = feature[2:] if feature.startswith("z_") else feature
        base_feature = str(base_feature).strip()
        if base_feature in summary_set:
            _add(base_feature)
            continue

        mean_feature = f"mean_{base_feature}"
        if mean_feature in summary_set:
            _add(mean_feature)
            continue

        if len(exact_summary_matches) == 0 and base_feature not in summary_set:
            for col in _summary_matches_for_feature(base_feature):
                _add(col)
            if base_feature in track_cols:
                for col in _summary_matches_for_feature(base_feature):
                    _add(col)
    return out


# ---------------------------------------------------------------------------
# DTW computation
# ---------------------------------------------------------------------------

def calculate_dtw(
    df_tracks,
    features=[
        "z_mean_square_displacement",
        "z_speed",
        "z_mean_dead_dye",
        "tcell_contact",
        "organoid_contact",
    ],
):
    print("- Calculating the dynamic time warping distance matrix")
    df_tracks = df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])

    nan_counts = df_tracks[features].isna().sum()
    features_with_nan = nan_counts[nan_counts > 0]
    if len(features_with_nan) > 0:
        print("  ⚠️ Warning: Found NaN values in selected features:")
        for feat, count in features_with_nan.items():
            print(f"     - {feat}: {count} NaN values")
        print("  → Filling NaN with 0 for DTW calculation")

    df_tracks_filled = df_tracks.copy()
    df_tracks_filled[features] = df_tracks_filled[features].fillna(0)

    dtw_input_tracks = []
    dtw_rownames = []
    for (TrackID, sample_name), group in df_tracks_filled.groupby(["TrackID", "sample_name"]):
        track_features = group[features].to_numpy().astype(np.double)
        dtw_rownames.append(f"{TrackID}--{sample_name}")
        dtw_input_tracks.append(track_features)

    dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input_tracks)
    dtw_distance_matrix = pd.DataFrame(dtw_distance_matrix, index=dtw_rownames, columns=dtw_rownames)
    dtw_distance_matrix = _handle_nan_in_distance_matrix(
        dtw_distance_matrix, context="DTW distance matrix"
    )
    return dtw_distance_matrix


# ---------------------------------------------------------------------------
# UMAP clustering
# ---------------------------------------------------------------------------

def cluster_umap(
    umap_embedding,
    config=None,
    nr_of_clusters=None,
    df_tracks=None,
    df_tracks_summarized=None,
    random_state=None,
    output_dir=None,
    cell_type="tcell",
    cluster_percentage_group_by=None,
    plot_results=True,
    output_subdir_name="results",
    plot_feature_cols=None,
):
    assert config is not None or all(
        [output_dir, nr_of_clusters]
    ), "Either 'config' or 'output_dir, nr_of_clusters' parameters must be supplied"

    print("- Performing clustering on the UMAP data")
    if all([output_dir, nr_of_clusters]) is None:
        output_dir = Path(config["output_dir"])
        nr_of_clusters = config["nr_of_clusters"]

    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    results_outdir = Path(analysis_outdir, str(output_subdir_name).strip() or "results")
    analysis_outdir.mkdir(parents=True, exist_ok=True)
    feature_outdir.mkdir(parents=True, exist_ok=True)
    results_outdir.mkdir(parents=True, exist_ok=True)

    if df_tracks is None:
        df_tracks = pd.read_csv(
            Path(feature_outdir, "BEHAV3D_combined_track_features_filtered.csv")
        )
    if df_tracks_summarized is None:
        df_tracks_summarized = pd.read_csv(
            Path(feature_outdir, "BEHAV3D_combined_track_features_summarized.csv")
        )

    df_umap = pd.merge(df_tracks_summarized, umap_embedding, how="left")

    kmeans = KMeans(n_clusters=nr_of_clusters, n_init=100, random_state=random_state)
    df_umap["ClusterID"] = kmeans.fit_predict(
        df_umap[umap_embedding.columns.difference(df_tracks_summarized.columns)]
    )
    df_umap["ClusterID"] = df_umap["ClusterID"] + 1
    df_umap["ClusterID"] = df_umap["ClusterID"].astype("category")

    df_umap_out_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv")
    print(f"- Writing clustered tracks to {df_umap_out_path}")
    df_umap.to_csv(df_umap_out_path, sep=",", index=False)

    line_cols = [c for c in df_tracks_summarized.columns if c.endswith("_line_condition")]
    sample_cols = line_cols + ["sample_name"]
    if "well" in df_tracks_summarized.columns:
        sample_cols.append("well")
    if "exp_nr" in df_tracks_summarized.columns:
        sample_cols.append("exp_nr")
    if cluster_percentage_group_by is None:
        cluster_percentage_context_cols = []
    elif isinstance(cluster_percentage_group_by, str):
        cluster_percentage_context_cols = [cluster_percentage_group_by]
    else:
        cluster_percentage_context_cols = list(cluster_percentage_group_by)
    metadata_cols = [
        c for c in sample_cols + cluster_percentage_context_cols
        if c in df_umap.columns
    ]
    if plot_feature_cols is None:
        info_cols = list(df_umap.drop(columns=["TrackID", "UMAP1", "UMAP2", "ClusterID"]).columns)
    else:
        selected_plot_cols = [c for c in list(plot_feature_cols) if c in df_umap.columns]
        if len(selected_plot_cols) == 0:
            print(
                "  ⚠️ Warning: None of the selected DTW features could be resolved to summarized columns for plotting. "
                "Output PDFs will show metadata/context columns only."
            )
        info_cols = list(dict.fromkeys(metadata_cols + selected_plot_cols))

    cluster_UMAP_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.pdf")
    print(f"- Plotting clustered UMAP plots with displayed Track features to {cluster_UMAP_path}")
    plot_feature_umap(
        df_umap=df_umap,
        info_cols=info_cols,
        sample_cols=sample_cols,
        outpath=cluster_UMAP_path,
        rows_per_page=4,
        nr_cols=2,
        rows_first_img=2,
        figsize=(8.27, 11.69),
        plot_results=plot_results,
    )

    cluster_features_heatmap_path = Path(
        results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_feature_heatmap.pdf"
    )
    print(f"- Plotting heatmaps with summarized cluster features to {cluster_features_heatmap_path}")
    plot_clustering_feature_heatmap(
        df_umap,
        info_cols,
        sample_cols,
        cluster_features_heatmap_path,
        rows_per_page=7,
        nr_cols=2,
        figsize=(8.27, 11.69),
        plot_results=plot_results,
    )

    all_line_cols = [c for c in df_umap.columns if c.endswith("_line_condition")]
    this_cell_line_cols = [c for c in all_line_cols if f"_{cell_type}_" in c]
    other_line_cols = [c for c in all_line_cols if c not in this_cell_line_cols]
    line_cols = this_cell_line_cols + other_line_cols

    if cluster_percentage_group_by is None:
        cluster_percentage_groups = []
    elif isinstance(cluster_percentage_group_by, str):
        cluster_percentage_groups = [cluster_percentage_group_by]
    else:
        cluster_percentage_groups = list(cluster_percentage_group_by)
    extra_group_cols = [
        col for col in cluster_percentage_groups
        if col in df_umap.columns and col not in line_cols
    ]

    group_cols_sample = line_cols + extra_group_cols + ["sample_name", "ClusterID"]
    group_cols_sample_only = line_cols + extra_group_cols + ["sample_name"]

    df_clust_perc = (
        df_umap
        .groupby(group_cols_sample, observed=True)
        .size()
        .reset_index(name="count")
    )
    sample_totals = (
        df_clust_perc
        .groupby(group_cols_sample_only, observed=True)["count"]
        .sum()
        .reset_index(name="sample_total")
    )
    combo_totals = (
        df_clust_perc
        .groupby(line_cols, observed=True)["count"]
        .sum()
        .reset_index(name="combo_total")
    ) if line_cols else None

    df_clust_perc = df_clust_perc.merge(sample_totals, how="left", on=group_cols_sample_only)
    if combo_totals is not None:
        df_clust_perc = df_clust_perc.merge(combo_totals, how="left", on=line_cols)
    df_clust_perc["percentage"] = df_clust_perc["count"] / df_clust_perc["sample_total"]

    cluster_percentage_plot_prefix = Path(
        results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_percentages"
    )
    print(f"- Plotting percentage plots of each cluster per sample to {cluster_percentage_plot_prefix}_*.pdf")
    plot_cluster_percentage_bars(
        df_clust_perc,
        cluster_percentage_plot_prefix,
        group_by_columns=cluster_percentage_group_by,
        plot_results=plot_results,
    )

    df_clust_perc = df_clust_perc.reset_index(drop=True)
    df_clust_perc_out_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_percentages.csv")
    print(f"- Writing cluster percentages to {df_clust_perc_out_path}")
    df_clust_perc.to_csv(df_clust_perc_out_path, sep=",", index=False)

    df_clust_tracks_out_path = Path(
        results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv"
    )
    print(f"- Writing clustered track info to {df_clust_tracks_out_path}")
    df_tracks = pd.merge(df_tracks, df_umap[["TrackID", "ClusterID"]], on="TrackID", how="left")
    df_tracks.to_csv(df_clust_tracks_out_path, sep=",", index=False)
    return()


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_tcell_analysis(
    config=None,
    output_dir=None,
    df_tracks_path=None,
    df_tracks_summarized_path=None,
    cell_type="tcell",
    columns_to_use=[
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
        "tcell_contact",
        "organoid_contact",
    ],
    columns_to_normalize=[
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
    ],
    umap_minimal_distance=None,
    umap_n_neighbors=None,
    nr_of_clusters=None,
    cluster_percentage_group_by=None,
    plot_results=True,
    seed=42,
    output_subdir_name="results",
    feature_scaling_preset=None,
    min_track_length=None,
    max_track_length=None,
):
    print(f"--------------- Performing {cell_type} behavioral analysis ---------------")
    start_time = time.time()

    if output_dir is None:
        output_dir = config["output_dir"]
    if feature_scaling_preset is None and config is not None:
        feature_scaling_preset = config.get("feature_scaling_preset", None)

    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    analysis_outdir.mkdir(parents=True, exist_ok=True)
    feature_outdir.mkdir(parents=True, exist_ok=True)

    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    if df_tracks_summarized_path is None:
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_summarized.csv")

    df_tracks = pd.read_csv(df_tracks_path, low_memory=False)
    if Path(df_tracks_summarized_path).exists():
        df_tracks_summarized = pd.read_csv(df_tracks_summarized_path, low_memory=False)
    else:
        df_tracks_summarized = summarize_feature_dtw_tracks(df_tracks)
    df_tracks = df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])

    if min_track_length is not None or max_track_length is not None:
        df_tracks = filter_feature_dtw_track_lengths(
            df_tracks,
            min_track_length=min_track_length,
            max_track_length=max_track_length,
        )
        df_tracks_summarized = summarize_feature_dtw_tracks(df_tracks)

    if df_tracks_summarized["track_length"].nunique() != 1:
        print("Warning: The track lengths are not cut to similar length, this might influence dynamic time warping")
        print("Set 'min_track_length' and 'max_track_length' to the same value to create equal tracks")

    if feature_scaling_preset is None:
        non_wildcard_cols_to_normalize = []
        for col in columns_to_normalize:
            non_wildcard_cols_to_normalize += expand_column_patterns(col, df_tracks.columns)
        df_tracks = normalize_track_features(df_tracks, columns=non_wildcard_cols_to_normalize)
        dtw_features = [x if x not in columns_to_normalize else f"z_{x}" for x in columns_to_use]
        plot_feature_cols = _resolve_selected_plot_feature_columns(
            df_tracks_summarized=df_tracks_summarized,
            df_tracks=df_tracks,
            selected_features=columns_to_use,
        )
    elif feature_scaling_preset == "original_behav3d":
        df_tracks, dtw_features = apply_original_behav3d_feature_scaling(
            df_tracks,
            cell_type=cell_type,
        )
        plot_feature_cols = _resolve_selected_plot_feature_columns(
            df_tracks_summarized=df_tracks_summarized,
            df_tracks=df_tracks,
            selected_features=[
                "mean_square_displacement",
                "speed",
                "mean_dead_dye",
                f"{cell_type}_contact",
                "organoid_contact",
            ],
        )
    else:
        raise ValueError(
            "Invalid feature_scaling_preset. Supported values are None and 'original_behav3d'."
        )

    dtw_distance_matrix = calculate_dtw(df_tracks, features=dtw_features)

    umap_embedding = fit_umap(
        dtw_distance_matrix=dtw_distance_matrix,
        umap_n_neighbors=umap_n_neighbors,
        umap_minimal_distance=umap_minimal_distance,
        random_state=seed,
    )

    df_clusters = cluster_umap(
        umap_embedding=umap_embedding,
        output_dir=output_dir,
        nr_of_clusters=nr_of_clusters,
        df_tracks=df_tracks,
        df_tracks_summarized=df_tracks_summarized,
        cell_type=cell_type,
        cluster_percentage_group_by=cluster_percentage_group_by,
        plot_results=plot_results,
        random_state=seed,
        output_subdir_name=output_subdir_name,
        plot_feature_cols=plot_feature_cols,
    )
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return df_clusters


# ---------------------------------------------------------------------------
# Path helpers and QC wrapper (unchanged from before)
# ---------------------------------------------------------------------------

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
        plot_results=False,
    )
    plot_clustering_feature_heatmap(
        df_plot,
        info_cols,
        sample_cols,
        heatmap_path,
        rows_per_page=7,
        nr_cols=2,
        figsize=(8.27, 11.69),
        plot_results=False,
    )
    df_clust_perc = _feature_dtw_cluster_percentage_groups(df_plot)
    if df_clust_perc is not None:
        df_clust_perc.to_csv(perc_prefix.with_suffix(".csv"), index=False)
        plot_cluster_percentage_bars(df_clust_perc, perc_prefix, group_by_columns=None, plot_results=False)

    return {
        "umap_pdf": str(cluster_umap_path),
        "heatmap_pdf": str(heatmap_path),
        "percentages_pdf": str(perc_prefix.with_suffix(".pdf")),
    }
