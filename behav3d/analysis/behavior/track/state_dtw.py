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
from scipy.stats import chi2_contingency, ttest_rel, wilcoxon

from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.features.state_descriptive_features import extract_descibing_track_state_features
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
from behav3d.analysis.behavior.utils import (
    _mixed_label_sort_key,
    _sanitize_filename_token,
    _save_adata_obs_csv,
)
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


_DIAGNOSTICS_GROUP_COL_CANDIDATES = ("origin_cell_type",)


def _resolve_diagnostics_group_col(adata_tracks, group_col=None):
    """Pick a per-track grouping column (e.g. population tag) to break diagnostics down by."""
    if group_col is not None:
        candidates = [str(group_col)]
    else:
        candidates = list(_DIAGNOSTICS_GROUP_COL_CANDIDATES)
    for col in candidates:
        if col not in adata_tracks.obs.columns:
            continue
        n_unique = adata_tracks.obs[col].dropna().nunique()
        if n_unique >= 2:
            return col
    return None


_DIAGNOSTICS_REPLICATE_COL_CANDIDATES = ("well", "sample_name")


def _resolve_diagnostics_replicate_col(adata_tracks, replicate_col=None, exclude=()):
    """Pick a per-track replicate column (e.g. well) to assess between-replicate variation."""
    exclude = set(exclude)
    if replicate_col is not None:
        candidates = [str(replicate_col)]
    else:
        candidates = [c for c in _DIAGNOSTICS_REPLICATE_COL_CANDIDATES if c not in exclude]
    for col in candidates:
        if col not in adata_tracks.obs.columns:
            continue
        n_unique = adata_tracks.obs[col].dropna().nunique()
        if n_unique >= 2:
            return col
    return None


def _chi2_group_cluster_test(df, *, cluster_key, group_col):
    """Chi-square test of independence between cluster membership and group_col."""
    table = pd.crosstab(df[cluster_key], df[group_col])
    if table.shape[0] < 2 or table.shape[1] < 2 or table.to_numpy().sum() == 0:
        return {"n_tracks": int(table.to_numpy().sum()), "chi2": np.nan, "p_value": np.nan, "dof": np.nan, "cramers_v": np.nan}
    chi2, p, dof, _ = chi2_contingency(table)
    n = table.to_numpy().sum()
    min_dim = min(table.shape) - 1
    cramers_v = float(np.sqrt((chi2 / n) / min_dim)) if min_dim > 0 and n > 0 else np.nan
    return {"n_tracks": int(n), "chi2": float(chi2), "p_value": float(p), "dof": int(dof), "cramers_v": cramers_v}


def _full_replicate_cluster_group_counts(df, *, replicate_col, cluster_key, group_col):
    """Cross-tab of (replicate, cluster, group) counts, with zero-count combinations filled in."""
    full_index = pd.MultiIndex.from_product(
        [
            sorted(df[replicate_col].unique().tolist()),
            sorted(df[cluster_key].unique().tolist(), key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x))),
            sorted(df[group_col].unique().tolist()),
        ],
        names=[replicate_col, cluster_key, group_col],
    )
    counts = (
        df.groupby([replicate_col, cluster_key, group_col], observed=True)
        .size()
        .reindex(full_index, fill_value=0)
        .rename("n_tracks")
        .reset_index()
    )
    denom = counts.groupby([replicate_col, group_col], observed=True)["n_tracks"].transform("sum")
    counts["frac_within_replicate_group"] = np.where(denom > 0, counts["n_tracks"] / denom, 0.0)
    return counts


def _paired_group_stats_by_replicate(replicate_group_counts, *, cluster_key, group_col, replicate_col):
    """For exactly two groups, compare per-cluster group fractions across replicates (paired)."""
    group_values = sorted(replicate_group_counts[group_col].unique().tolist())
    if len(group_values) != 2:
        return None
    grp_a, grp_b = group_values
    cluster_order = sorted(
        replicate_group_counts[cluster_key].unique().tolist(),
        key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)),
    )

    rows = []
    for cluster in cluster_order:
        pivot = (
            replicate_group_counts[replicate_group_counts[cluster_key] == cluster]
            .pivot(index=replicate_col, columns=group_col, values="frac_within_replicate_group")
            .reindex(columns=group_values, fill_value=0.0)
        )
        vals_a = pivot[grp_a].to_numpy(dtype=float)
        vals_b = pivot[grp_b].to_numpy(dtype=float)
        n_rep = len(pivot)
        row = {
            cluster_key: cluster,
            "n_replicates": n_rep,
            f"mean_frac_{grp_a}": float(np.mean(vals_a)) if n_rep > 0 else np.nan,
            f"mean_frac_{grp_b}": float(np.mean(vals_b)) if n_rep > 0 else np.nan,
            "mean_diff": float(np.mean(vals_a - vals_b)) if n_rep > 0 else np.nan,
            "std_diff": float(np.std(vals_a - vals_b, ddof=1)) if n_rep > 1 else np.nan,
            "paired_ttest_stat": np.nan,
            "paired_ttest_p": np.nan,
            "wilcoxon_stat": np.nan,
            "wilcoxon_p": np.nan,
        }
        if n_rep >= 2 and not np.allclose(vals_a, vals_b):
            try:
                t_stat, t_p = ttest_rel(vals_a, vals_b)
                row["paired_ttest_stat"] = float(t_stat)
                row["paired_ttest_p"] = float(t_p)
            except Exception:
                pass
            try:
                w_stat, w_p = wilcoxon(vals_a, vals_b)
                row["wilcoxon_stat"] = float(w_stat)
                row["wilcoxon_p"] = float(w_p)
            except Exception:
                pass
        rows.append(row)
    result = pd.DataFrame(rows)
    result.attrs["group_a"] = grp_a
    result.attrs["group_b"] = grp_b
    return result


def _save_diagnostics(
    adata_tracks,
    distances,
    outfolder,
    *,
    cluster_key="ClusterID",
    max_heatmap_tracks=200,
    random_state=123,
    group_col=None,
):
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)
    pdf_path = outfolder / "dtaidistance_clustering_diagnostics.pdf"
    counts_csv = outfolder / "dtaidistance_cluster_counts.csv"
    medoids_csv = outfolder / "dtaidistance_cluster_medoids.csv"
    umap_csv = outfolder / "dtaidistance_umap_clusters.csv"

    labels = pd.Series(adata_tracks.obs[cluster_key]).astype(str)
    resolved_cluster_order = _apply_state_order(
        sorted(labels.unique().tolist(), key=_mixed_label_sort_key),
        _get_classification_state_order(adata_tracks, cluster_key),
    )
    color_map = _normalize_label_color_map(
        resolved_cluster_order, colors=_get_classification_state_colors(adata_tracks, cluster_key),
    )
    counts = labels.value_counts().rename_axis(cluster_key).reset_index(name="n_tracks")
    counts.to_csv(counts_csv, index=False)

    resolved_group_col = _resolve_diagnostics_group_col(adata_tracks, group_col=group_col)
    group_counts = None
    group_counts_csv = None
    if resolved_group_col is not None:
        groups = pd.Series(adata_tracks.obs[resolved_group_col]).astype(str)
        group_df = pd.DataFrame({cluster_key: labels.to_numpy(), resolved_group_col: groups.to_numpy()})
        group_counts = (
            group_df.groupby([cluster_key, resolved_group_col], observed=True)
            .size()
            .rename("n_tracks")
            .reset_index()
        )
        group_counts["frac_of_cluster"] = group_counts["n_tracks"] / group_counts.groupby(cluster_key)[
            "n_tracks"
        ].transform("sum")
        group_counts["frac_of_group"] = group_counts["n_tracks"] / group_counts.groupby(resolved_group_col)[
            "n_tracks"
        ].transform("sum")
        group_counts_csv = outfolder / f"dtaidistance_cluster_counts_by_{resolved_group_col}.csv"
        group_counts.to_csv(group_counts_csv, index=False)

    resolved_replicate_col = None
    replicate_group_counts = None
    replicate_group_counts_csv = None
    chi2_stats = None
    chi2_stats_csv = None
    paired_stats = None
    paired_stats_csv = None
    if resolved_group_col is not None:
        resolved_replicate_col = _resolve_diagnostics_replicate_col(adata_tracks, exclude={resolved_group_col})
    if resolved_group_col is not None and resolved_replicate_col is not None:
        reps = pd.Series(adata_tracks.obs[resolved_replicate_col]).astype(str)
        rep_df = pd.DataFrame(
            {
                cluster_key: labels.to_numpy(),
                resolved_group_col: groups.to_numpy(),
                resolved_replicate_col: reps.to_numpy(),
            }
        )
        replicate_group_counts = _full_replicate_cluster_group_counts(
            rep_df,
            replicate_col=resolved_replicate_col,
            cluster_key=cluster_key,
            group_col=resolved_group_col,
        )
        replicate_group_counts_csv = (
            outfolder / f"dtaidistance_cluster_counts_by_{resolved_group_col}_by_{resolved_replicate_col}.csv"
        )
        replicate_group_counts.to_csv(replicate_group_counts_csv, index=False)

        overall_chi2 = _chi2_group_cluster_test(rep_df, cluster_key=cluster_key, group_col=resolved_group_col)
        overall_chi2["scope"] = "overall"
        chi2_rows = [overall_chi2]
        for rep_val, sub in rep_df.groupby(resolved_replicate_col, observed=True):
            row = _chi2_group_cluster_test(sub, cluster_key=cluster_key, group_col=resolved_group_col)
            row["scope"] = str(rep_val)
            chi2_rows.append(row)
        chi2_stats = pd.DataFrame(chi2_rows)[
            ["scope", "n_tracks", "chi2", "dof", "p_value", "cramers_v"]
        ]
        chi2_stats_csv = (
            outfolder / f"dtaidistance_{resolved_group_col}_vs_cluster_chi2_by_{resolved_replicate_col}.csv"
        )
        chi2_stats.to_csv(chi2_stats_csv, index=False)

        paired_stats = _paired_group_stats_by_replicate(
            replicate_group_counts,
            cluster_key=cluster_key,
            group_col=resolved_group_col,
            replicate_col=resolved_replicate_col,
        )
        if paired_stats is not None:
            paired_stats_csv = (
                outfolder / f"dtaidistance_{resolved_group_col}_paired_by_{resolved_replicate_col}.csv"
            )
            paired_stats.to_csv(paired_stats_csv, index=False)

    medoid_cols = [
        c
        for c in [
            "sample_name",
            "TrackID",
            "trajectory_window_id",
            "position_t_min",
            "position_t_max",
            cluster_key,
            f"{cluster_key}_medoid",
            f"{cluster_key}_medoid_rank",
        ]
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
            for cluster in resolved_cluster_order:
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

            if resolved_group_col is not None:
                groups = pd.Series(adata_tracks.obs[resolved_group_col]).astype(str)
                group_order = sorted(groups.unique().tolist())
                fig, ax = plt.subplots(figsize=(7, 6))
                gcmap = plt.get_cmap("Set2")
                gcolor_map = {grp: gcmap(i % gcmap.N) for i, grp in enumerate(group_order)}
                for grp in group_order:
                    mask = groups.to_numpy() == grp
                    ax.scatter(
                        umap_embedding[mask, 0],
                        umap_embedding[mask, 1],
                        s=18,
                        alpha=0.8,
                        color=gcolor_map[grp],
                        label=grp,
                        linewidths=0,
                    )
                ax.set_xlabel("UMAP1")
                ax.set_ylabel("UMAP2")
                ax.set_title(f"One-hot dtaidistance UMAP by {resolved_group_col}")
                ax.legend(title=str(resolved_group_col), loc="best", frameon=False, markerscale=1.4)
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

        counts_by_cluster = counts.set_index(counts[cluster_key].astype(str))["n_tracks"].reindex(
            resolved_cluster_order, fill_value=0
        )
        fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(counts_by_cluster) + 1.5)))
        ax.barh(
            resolved_cluster_order,
            counts_by_cluster.to_numpy(),
            color=[color_map[c] for c in resolved_cluster_order],
        )
        ax.invert_yaxis()
        ax.set_xlabel("N tracks")
        ax.set_ylabel(str(cluster_key))
        ax.set_title("One-hot dtaidistance cluster counts")
        ax.grid(axis="x", alpha=0.2)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        if group_counts is not None:
            cluster_order = resolved_cluster_order
            group_order = sorted(group_counts[resolved_group_col].unique().tolist())
            gcmap = plt.get_cmap("Set2")
            gcolor_map = {grp: gcmap(i % gcmap.N) for i, grp in enumerate(group_order)}
            pivot_counts = (
                group_counts.pivot(index=cluster_key, columns=resolved_group_col, values="n_tracks")
                .reindex(index=cluster_order, columns=group_order)
                .fillna(0)
            )

            y_pos = np.arange(len(cluster_order))
            fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(cluster_order) + 1.5)))
            left = np.zeros(len(cluster_order))
            for grp in group_order:
                vals = pivot_counts[grp].to_numpy()
                ax.barh(y_pos, vals, left=left, height=0.7, color=gcolor_map[grp], label=grp)
                left = left + vals
            ax.set_yticks(y_pos)
            ax.set_yticklabels(cluster_order)
            ax.invert_yaxis()
            ax.set_xlabel("N tracks")
            ax.set_ylabel(str(cluster_key))
            ax.set_title(f"One-hot dtaidistance cluster counts by {resolved_group_col} (stacked)")
            ax.legend(title=str(resolved_group_col), loc="best", frameon=False)
            ax.grid(axis="x", alpha=0.2)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

            pivot_frac = pivot_counts.div(pivot_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
            fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(cluster_order) + 1.5)))
            left = np.zeros(len(cluster_order))
            for grp in group_order:
                vals = pivot_frac[grp].to_numpy()
                ax.barh(y_pos, vals, left=left, height=0.7, color=gcolor_map[grp], label=grp)
                left = left + vals
            ax.set_yticks(y_pos)
            ax.set_yticklabels(cluster_order)
            ax.invert_yaxis()
            ax.set_xlim(0, 1)
            ax.set_xlabel(f"Fraction of cluster ({resolved_group_col})")
            ax.set_ylabel(str(cluster_key))
            ax.set_title(f"One-hot dtaidistance cluster composition by {resolved_group_col}")
            ax.legend(title=str(resolved_group_col), loc="lower right", frameon=False)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

            if replicate_group_counts is not None:
                replicate_order = sorted(replicate_group_counts[resolved_replicate_col].unique().tolist())
                group_order_rep = sorted(replicate_group_counts[resolved_group_col].unique().tolist())
                cluster_color_map = color_map

                fig, axes = plt.subplots(
                    1, len(replicate_order), figsize=(3.2 * len(replicate_order) + 1, 5), sharey=True
                )
                axes = np.atleast_1d(axes)
                for ax, rep_val in zip(axes, replicate_order):
                    sub = replicate_group_counts[replicate_group_counts[resolved_replicate_col] == rep_val]
                    pivot_rep = (
                        sub.pivot(index=resolved_group_col, columns=cluster_key, values="frac_within_replicate_group")
                        .reindex(index=group_order_rep, columns=cluster_order)
                        .fillna(0.0)
                    )
                    bottom = np.zeros(len(group_order_rep))
                    x_pos = np.arange(len(group_order_rep))
                    for c in cluster_order:
                        vals = pivot_rep[c].to_numpy()
                        ax.bar(x_pos, vals, bottom=bottom, width=0.6, color=cluster_color_map[c], label=c)
                        bottom = bottom + vals
                    ax.set_xticks(x_pos)
                    ax.set_xticklabels(group_order_rep, rotation=30, ha="right", fontsize=8)
                    ax.set_title(str(rep_val), fontsize=10, fontweight="bold")
                    ax.set_ylim(0, 1)
                axes[0].set_ylabel(f"Fraction of tracks (by {cluster_key})")
                handles, labels_leg = axes[0].get_legend_handles_labels()
                fig.legend(
                    handles,
                    labels_leg,
                    title=str(cluster_key),
                    loc="lower center",
                    ncol=min(len(cluster_order), 8),
                    frameon=False,
                    fontsize=7,
                )
                fig.suptitle(
                    f"Cluster composition by {resolved_group_col}, split by {resolved_replicate_col}",
                    fontsize=12,
                    fontweight="bold",
                )
                fig.tight_layout(rect=(0.02, 0.1, 1, 0.92))
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

            if chi2_stats is not None:
                fig, ax = plt.subplots(figsize=(8.27, 11.69))
                ax.axis("off")
                lines = [
                    f"{resolved_group_col} vs {cluster_key} — chi-square test of independence",
                    "",
                ]
                for _, row in chi2_stats.iterrows():
                    lines.append(
                        f"{str(row['scope']):>16s}:  n={int(row['n_tracks'])}, chi2={row['chi2']:.3g}, "
                        f"dof={row['dof']}, p={row['p_value']:.3g}, Cramer's V={row['cramers_v']:.3g}"
                    )
                if paired_stats is not None and len(paired_stats) > 0:
                    grp_a = paired_stats.attrs.get("group_a", "group_a")
                    grp_b = paired_stats.attrs.get("group_b", "group_b")
                    max_n_rep = int(paired_stats["n_replicates"].max())
                    lines += [
                        "",
                        f"Paired comparison across {resolved_replicate_col} ({grp_a} vs {grp_b}, "
                        "fraction of that well's population per cluster):",
                        "",
                    ]
                    for _, row in paired_stats.iterrows():
                        lines.append(
                            f"{str(cluster_key)} {row[cluster_key]}: "
                            f"mean {grp_a}={row[f'mean_frac_{grp_a}']:.3f}, "
                            f"mean {grp_b}={row[f'mean_frac_{grp_b}']:.3f}, "
                            f"diff={row['mean_diff']:.3f}, "
                            f"paired-t p={row['paired_ttest_p']:.3g}, wilcoxon p={row['wilcoxon_p']:.3g}"
                        )
                    lines += [
                        "",
                        f"Caveat: only {max_n_rep} replicate(s) ({resolved_replicate_col}) available — "
                        "treat these p-values as indicative, not confirmatory.",
                    ]
                ax.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", fontsize=9, family="monospace")
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
        "cluster_counts_by_group_csv": str(group_counts_csv) if group_counts_csv is not None else None,
        "diagnostics_group_col": resolved_group_col,
        "diagnostics_replicate_col": resolved_replicate_col,
        "cluster_counts_by_group_by_replicate_csv": (
            str(replicate_group_counts_csv) if replicate_group_counts_csv is not None else None
        ),
        "group_vs_cluster_chi2_csv": str(chi2_stats_csv) if chi2_stats_csv is not None else None,
        "group_paired_by_replicate_csv": str(paired_stats_csv) if paired_stats_csv is not None else None,
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
        _save_adata_obs_csv(adata_tracks, output_path)
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
    group_col=None,
    save_outputs=True,
    verbose=True,
):
    """Write diagnostics for a saved one-hot dtaidistance clustering model.

    group_col optionally breaks the cluster-count/UMAP diagnostics down by a
    per-track grouping column (default: auto-detects "origin_cell_type" if present).
    """
    paths = _resolve_dtaidistance_paths(output_dir, cell_type)
    resolved_cluster_key = _resolve_cluster_key(adata_tracks, cluster_key=cluster_key)
    if _dtai_meta(adata_tracks).get("method") == "original_behav3d_feature_dtw":
        raise ValueError(
            "Diagnostics are not available for the 'Original BEHAV3D' feature-DTW model — "
            "no pairwise distance matrix is stored for this clustering method. Use "
            "_save_feature_dtw_quality_control() instead."
        )
    distances = _validate_distance_matrix(adata_tracks.X, adata_tracks.n_obs)
    plot_paths = _save_diagnostics(
        adata_tracks,
        distances,
        paths["quality_control_outfolder"] if outfolder is None else outfolder,
        cluster_key=resolved_cluster_key,
        max_heatmap_tracks=int(max_heatmap_tracks),
        random_state=int(random_state),
        group_col=group_col,
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


def save_dtaidistance_exemplar_overview(
    adata_tracks,
    output_dir,
    cell_type,
    *,
    outfolder=None,
    n_per_cluster=10,
    random_state=None,
    verbose=True,
):
    """Save only the exemplar overview grid PDF to outfolder (default: quality_control/)."""
    paths = _resolve_dtaidistance_paths(output_dir, cell_type)
    meta = _dtai_meta(adata_tracks)
    resolved_cluster_key = _resolve_cluster_key(adata_tracks)
    state_col = str(meta.get("state_col", FULL_STATE_COL))
    time_col = str(meta.get("time_col", "position_t"))
    if random_state is None:
        random_state = int(meta.get("random_state", 123))
    adata_filt = _load_filtered_state_adata_for_model(
        adata_tracks, output_dir, cell_type, verbose=verbose,
    )
    _ensure_exemplar_coordinate_columns(
        adata_filt, output_dir=output_dir, cell_type=cell_type, require_pixel_for_video=False,
    )
    dest = Path(outfolder) if outfolder is not None else paths["quality_control_outfolder"]
    dest.mkdir(parents=True, exist_ok=True)
    fig, _, _ = plot_exemplar_tracks_by_cluster(
        adata_filt, adata_tracks,
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
    overview_pdf = dest / "example_tracks_overview.pdf"
    with PdfPages(overview_pdf) as pdf:
        pdf.savefig(fig, bbox_inches="tight", dpi=300)
    plt.close(fig)
    if bool(verbose):
        _winfo("trajectory-dtai", f"saved exemplar overview: {overview_pdf}")
    return str(overview_pdf)


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
    split_long_tracks = bool(meta.get("split_long_tracks", False))
    window_col = str(meta.get("trajectory_window_col", "trajectory_window_id"))
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
        split_long_tracks=split_long_tracks,
        window_col=window_col,
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
    exemplar_root = paths["outfolder"] / "example_tracks"
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
    split_long_tracks=False,
    trajectory_window_col="trajectory_window_id",
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
    split_long_tracks = bool(split_long_tracks)
    trajectory_window_col = str(trajectory_window_col)
    sequence_groupby_cols = list(groupby_cols)
    if split_long_tracks and trajectory_size is not None:
        sequence_groupby_cols = list(groupby_cols) + [trajectory_window_col]

    if bool(verbose):
        _winfo("trajectory-dtai", f"loading behavioral states: {adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)

    if bool(verbose):
        length_text = "all timepoints" if trajectory_size is None else f"{trim_mode} {trajectory_size} timepoints"
        action_text = "dividing" if split_long_tracks and trajectory_size is not None else "keeping"
        _winfo(
            "trajectory-dtai",
            f"filtering tracks with min_length={min_length}, {action_text} {length_text} | "
            f"state_col={FULL_STATE_COL}",
        )
    adata_filt = _filter_tracks_for_dtaidistance(
        adata_full,
        groupby_cols=list(groupby_cols),
        time_col=str(time_col),
        trajectory_size=trajectory_size,
        min_length=min_length,
        trim_mode=trim_mode,
        split_long_tracks=split_long_tracks,
        window_col=trajectory_window_col,
    )

    line_condition_cols = [c for c in adata_filt.obs.columns if c.endswith("_line_condition")]
    sequences, track_obs = extract_categorical_track_sequences(
        adata_filt,
        state_cols=state_cols,
        groupby_cols=sequence_groupby_cols,
        time_col=time_col,
        missing_policy=missing_policy,
        extra_meta_cols=("origin_cell_type", "well", *line_condition_cols),
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
        "sequence_groupby_cols": [str(c) for c in list(sequence_groupby_cols)],
        "time_col": str(time_col),
        "state_col": str(FULL_STATE_COL),
        "state_cols": list(state_cols),
        "behavioral_trajectory_size": None if trajectory_size is None else int(trajectory_size),
        "min_track_length": None if min_length is None else int(min_length),
        "trajectory_trim_mode": trim_mode,
        "split_long_tracks": bool(split_long_tracks),
        "trajectory_window_col": str(trajectory_window_col),
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
        raw_qc_dir = paths["quality_control_outfolder"] / "raw"
        raw_qc_dir.mkdir(parents=True, exist_ok=True)
        plot_paths = _save_diagnostics(
            adata_tracks,
            distances,
            raw_qc_dir,
            cluster_key=cluster_key,
            max_heatmap_tracks=int(max_heatmap_tracks),
            random_state=int(random_state),
        )
        adata_tracks.uns.setdefault("visualization", {})
        adata_tracks.uns["visualization"].update(plot_paths)

    if bool(plot_exemplars):
        exemplar_root = paths["outfolder"] / "example_tracks"
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
        paths["clustering_outfolder"].mkdir(parents=True, exist_ok=True)
        distance_csv = paths["clustering_outfolder"] / "categorical_dtai_distance_matrix.csv"
        pd.DataFrame(distances, index=adata_tracks.obs.index, columns=adata_tracks.obs.index).to_csv(distance_csv)
        adata_tracks.uns.setdefault("dtai_trajectory_clustering", {})
        adata_tracks.uns["dtai_trajectory_clustering"]["distance_matrix_csv"] = str(distance_csv)

    output_path = paths["outfolder"] / get_dtaidistance_track_trajectories_filename(cell_type)
    if bool(save_outputs):
        adata_tracks.write(output_path, compression="gzip")
        _save_adata_obs_csv(adata_tracks, output_path)
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
    groupby_cols = list(meta.get("groupby_cols", ["sample_name", "TrackID"]))
    sequence_groupby_cols = list(meta.get("sequence_groupby_cols", groupby_cols))
    time_col = str(meta.get("time_col", "position_t"))
    min_track_length = meta.get("min_track_length", trajectory_size)
    trim_mode = str(meta.get("trajectory_trim_mode", "last"))
    split_long_tracks = bool(meta.get("split_long_tracks", False))
    window_col = str(meta.get("trajectory_window_col", "trajectory_window_id"))

    if bool(verbose):
        _winfo("trajectory-dtai-clf", f"loading behavioral states: {adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)

    if bool(verbose):
        _winfo("trajectory-dtai-clf", f"building track features (trajectory_size={trajectory_size})")
    if split_long_tracks:
        adata_filt = _filter_tracks_for_dtaidistance(
            adata_full,
            groupby_cols=groupby_cols,
            time_col=time_col,
            trajectory_size=int(trajectory_size),
            min_length=None if min_track_length is None else int(min_track_length),
            trim_mode=trim_mode,
            split_long_tracks=True,
            window_col=window_col,
        )
        feature_adata, _ = extract_descibing_track_state_features(
            adata_filt,
            state_col=str(FULL_STATE_COL),
            group_col=sequence_groupby_cols,
        )
        feature_adata.uns["track_filtering"] = {
            "groupby_cols": [str(c) for c in list(groupby_cols)],
            "sequence_groupby_cols": [str(c) for c in list(sequence_groupby_cols)],
            "time_col": str(time_col),
            "state_col": str(FULL_STATE_COL),
            "behavioral_trajectory_size": int(trajectory_size),
            "min_track_length": None if min_track_length is None else int(min_track_length),
            "trajectory_trim_mode": trim_mode,
            "split_long_tracks": True,
            "trajectory_window_col": str(window_col),
        }
    else:
        feature_adata = _build_track_feature_adata(
            adata_full,
            state_col=str(FULL_STATE_COL),
            behavioral_trajectory_size=int(trajectory_size),
        )

    # Transfer ClusterID from dtaidistance model to feature adata by the
    # analysis trajectory key. Split windows keep original TrackID and add
    # trajectory_window_id, so TrackID alone is not unique in that mode.
    transfer_cols = [str(c) for c in sequence_groupby_cols if str(c) in model_adata.obs.columns]
    missing_transfer = [str(c) for c in transfer_cols if str(c) not in feature_adata.obs.columns]
    if missing_transfer:
        raise ValueError(f"Feature adata is missing dtaidistance transfer columns: {missing_transfer}")
    dtai_obs = model_adata.obs[[cluster_col] + transfer_cols].copy()
    dtai_obs = dtai_obs.reset_index(drop=True)
    dtai_obs = dtai_obs.set_index(transfer_cols)[cluster_col]

    feat_obs = feature_adata.obs.copy()
    feat_obs["_cluster"] = (
        pd.MultiIndex.from_frame(feat_obs[transfer_cols])
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
