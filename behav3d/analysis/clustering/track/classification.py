import time
import random
import re
import pickle
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_hex
import seaborn as sns
import umap
import scanpy as sc

from sklearn.cluster import KMeans, HDBSCAN
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    confusion_matrix,
    classification_report,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold

from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata
from behav3d.analysis.clustering.general import relabel_cluster_ids
from behav3d.features.state_descriptive_features import (
    extract_descibing_track_state_features, 
    scale_feature_blocks, 
    l2_normalize_features_blocks,
    drop_highly_correlated_features,
    drop_low_variance_features
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.analysis.clustering.general.leiden import (
    run_pca, 
    run_leiden_clustering
)
from behav3d.analysis.clustering.state.visualization.backprojection import (
    export_behavioral_state_backprojection_zarrs,
)
from behav3d.analysis.clustering.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    select_exemplar_tracks_by_cluster,
    _prepare_exemplar_backprojection_data,
)
from behav3d.analysis.clustering.track.visualization.plots.exemplar_coordinate_utils import (
    ensure_exemplar_coordinate_columns as _shared_ensure_exemplar_coordinate_columns,
    merge_coordinate_columns_into_obs as _shared_merge_coordinate_columns_into_obs,
    resolve_exemplar_positions_csv_path as _shared_resolve_exemplar_positions_csv_path,
)
from behav3d.analysis.clustering.general.visualization.plots import plot_top_ranking_features
# %matplotlib inline


def _vinfo(verbose, prefix, message):
    if not bool(verbose):
        return
    print(f"[{prefix}] INFO {str(message)}")


def _vstart(verbose, prefix, step_name):
    if bool(verbose):
        print(f"[{prefix}] START {step_name}")
    return time.perf_counter()


def _vdone(verbose, prefix, step_name, t_start):
    if bool(verbose):
        dt = time.perf_counter() - float(t_start)
        print(f"[{prefix}] DONE {step_name} | took {dt:.2f}s")


def _vsave(verbose, prefix, label, path):
    if not bool(verbose):
        return
    print(f"[{prefix}] SAVED {str(label)}: {path}")


# Backward-compatible shims (temporary internal use).
def _vprint(verbose, message):
    if bool(verbose):
        print(str(message))


def _vstep(verbose, prefix, step_name, t_start):
    _vdone(verbose, prefix, step_name, t_start)


def _resolve_output_dir(output_dir):
    """Resolve/create the global BEHAV3D output root directory."""
    if output_dir is None:
        raise ValueError("output_dir is required.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return output_dir_path


def _resolve_track_paths(output_dir, cell_type, output_subdir_name="behavioral_state_trajectories"):
    """Resolve canonical track-classification paths under analysis/<cell_type>/."""
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")

    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    analysis_outdir.mkdir(parents=True, exist_ok=True)

    state_outdir = analysis_outdir / "behavioral_states"
    state_outdir.mkdir(parents=True, exist_ok=True)

    outfolder = analysis_outdir / str(output_subdir_name)
    outfolder.mkdir(parents=True, exist_ok=True)

    return {
        "output_dir": root,
        "analysis_outdir": analysis_outdir,
        "state_outdir": state_outdir,
        "outfolder": outfolder,
    }


def _resolve_track_positions_csv_path(output_dir, cell_type):
    return _shared_resolve_exemplar_positions_csv_path(
        output_dir=output_dir,
        cell_type=cell_type,
    )


def _merge_coordinate_columns_into_obs(
    adata,
    positions_csv_path,
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
):
    return _shared_merge_coordinate_columns_into_obs(
        adata=adata,
        positions_csv_path=positions_csv_path,
        sample_col=sample_col,
        track_col=track_col,
        time_col=time_col,
    )


def _ensure_exemplar_coordinate_columns(
    adata,
    *,
    output_dir,
    cell_type,
    require_pixel_for_video=False,
):
    return _shared_ensure_exemplar_coordinate_columns(
        adata=adata,
        output_dir=output_dir,
        cell_type=cell_type,
        require_pixel_for_video=require_pixel_for_video,
    )


def _mixed_label_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def _sanitize_filename_token(value, fallback="plot"):
    token = str(value).strip()
    if token == "":
        token = str(fallback)
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", token)
    token = token.strip("._-")
    return token if token != "" else str(fallback)


def get_track_trajectories_filename(cell_type):
    cell_token = _sanitize_filename_token(cell_type, fallback="cell")
    return f"BEHAV3D_{cell_token}_behavioral_trajectories.h5ad"


def get_track_classifier_filename(cell_type):
    cell_token = _sanitize_filename_token(cell_type, fallback="cell")
    return f"track_classification_random_forest_{cell_token}.pkl"


def _resolve_track_classifier_path(output_dir, cell_type, output_subdir_name="behavioral_state_trajectories"):
    paths = _resolve_track_paths(
        output_dir=output_dir,
        cell_type=cell_type,
        output_subdir_name=output_subdir_name,
    )
    classifier_dir = paths["outfolder"] / "classification"
    classifier_dir.mkdir(parents=True, exist_ok=True)
    return classifier_dir / get_track_classifier_filename(cell_type)


def _to_numpy_2d(X):
    if hasattr(X, "toarray"):
        X = X.toarray()
    arr = np.asarray(X)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape={arr.shape}.")
    return arr


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


def evaluate_classifier_predictions(
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="track_classifier_qc",
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


def _evaluate_holdout_predictions(
    *,
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="track_classifier_holdout_qc",
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


def load_track_classifier_artifact(path):
    """Load a saved track-classifier artifact pickle."""
    with Path(path).open("rb") as f:
        artifact = pickle.load(f)
    if not isinstance(artifact, dict) or "classifier" not in artifact:
        raise ValueError("Invalid track classifier artifact: expected dict with key 'classifier'.")
    return artifact


def _coerce_track_classifier_artifact(classifier_artifact_or_path):
    if isinstance(classifier_artifact_or_path, (str, Path)):
        return load_track_classifier_artifact(classifier_artifact_or_path)
    if isinstance(classifier_artifact_or_path, dict) and "classifier" in classifier_artifact_or_path:
        return classifier_artifact_or_path
    raise ValueError(
        "classifier_artifact_or_path must be a classifier artifact dict or path to a pickled artifact."
    )


def _resolve_track_apply_input_adata(
    *,
    adata_full,
    adata_full_path,
    state_folder,
    cell_type,
):
    """
    Resolve apply input in priority order:
    1) in-memory adata_full
    2) explicit adata_full_path
    3) default behavioral_states path for cell_type
    """
    if adata_full is not None:
        return adata_full, "in_memory", None

    if adata_full_path is None:
        path = Path(state_folder, f"BEHAV3D_{cell_type}_behavioral_states.h5ad")
        source_kind = "default_behavioral_states_path"
    else:
        path = Path(adata_full_path)
        source_kind = "path"

    if not path.exists():
        raise FileNotFoundError(f"Full dataset h5ad not found: {path}")
    return sc.read_h5ad(path), source_kind, str(path)


def _build_track_preprocessing_spec(*, state_col, filter_cfg, feat_cfg):
    groupby_cols = filter_cfg.get("groupby_cols", ("sample_name", "TrackID"))
    if not isinstance(groupby_cols, (list, tuple)) or len(groupby_cols) == 0:
        groupby_cols = ("sample_name", "TrackID")
    return {
        "state_col": str(state_col),
        "groupby_cols": [str(c) for c in list(groupby_cols)],
        "time_col": str(filter_cfg.get("time_col", "position_t")),
        "behavioral_trajectory_size": int(filter_cfg.get("behavioral_trajectory_size", 100)),
        "use_fractions": bool(feat_cfg.get("use_fractions", True)),
        "use_bout_stats": bool(feat_cfg.get("use_bout_stats", True)),
        "include_bouts_mean_length": bool(feat_cfg.get("include_bouts_mean_length", True)),
        "include_bouts_nr": bool(feat_cfg.get("include_bouts_nr", True)),
        "include_bouts_max_length": bool(feat_cfg.get("include_bouts_max_length", True)),
        "use_transitions": bool(feat_cfg.get("use_transitions", True)),
        "use_bigrams": bool(feat_cfg.get("use_bigrams", True)),
        "use_trigrams": bool(feat_cfg.get("use_trigrams", True)),
        "ngram_weight": str(feat_cfg.get("ngram_weight", "count")),
    }


def _resolve_apply_state_col(artifact, requested_state_col):
    """
    Resolve state column for apply.
    If artifact provides training state column, enforce exact match.
    """
    requested = str(requested_state_col)
    artifact_state_col = str(artifact.get("state_col_used_for_training", "")).strip()
    if artifact_state_col == "":
        return requested, None
    if requested != artifact_state_col:
        mismatch = (
            "Requested state column does not match classifier training state column: "
            f"requested='{requested}', trained='{artifact_state_col}'."
        )
        raise ValueError(mismatch)
    return artifact_state_col, None


def _drop_disabled_bout_features(
    adata_state_features,
    include_bouts_mean_length=True,
    include_bouts_nr=True,
    include_bouts_max_length=True,
):
    drop_cols = []
    if not bool(include_bouts_mean_length):
        drop_cols.extend(
            [c for c in adata_state_features.var_names if str(c).startswith("bouts_mean_length_")]
        )
    if not bool(include_bouts_nr):
        drop_cols.extend(
            [c for c in adata_state_features.var_names if str(c).startswith("bouts_nr_")]
        )
    if not bool(include_bouts_max_length):
        drop_cols.extend(
            [c for c in adata_state_features.var_names if str(c).startswith("bouts_max_length_")]
        )
    if len(drop_cols) == 0:
        return adata_state_features

    keep_cols = [c for c in adata_state_features.var_names if c not in set(drop_cols)]
    adata_state_features = adata_state_features[:, keep_cols].copy()
    feature_blocks = dict(adata_state_features.uns.get("feature_blocks", {}))
    if "bout_stats" in feature_blocks:
        feature_blocks["bout_stats"] = [c for c in feature_blocks["bout_stats"] if c in keep_cols]
    adata_state_features.uns["feature_blocks"] = feature_blocks
    return adata_state_features


def _build_track_feature_adata(
    adata_full,
    *,
    state_col="full_behavioral_cluster",
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    behavioral_trajectory_size=100,
    use_fractions=True,
    use_bout_stats=True,
    include_bouts_mean_length=True,
    include_bouts_nr=True,
    include_bouts_max_length=True,
    use_transitions=True,
    use_bigrams=True,
    use_trigrams=True,
    ngram_weight="count",
):
    trajectory_size = int(behavioral_trajectory_size)
    adata_filt = filter_and_truncate_tracks_anndata(
        adata_full,
        groupby_cols=[str(c) for c in list(groupby_cols)],
        time_col=str(time_col),
        min_length=trajectory_size,
        max_length=trajectory_size,
    )
    adata_state_features, _ = extract_descibing_track_state_features(
        adata_filt,
        state_col=str(state_col),
        use_fractions=bool(use_fractions),
        use_bout_stats=bool(use_bout_stats),
        use_transitions=bool(use_transitions),
        use_bigrams=bool(use_bigrams),
        use_trigrams=bool(use_trigrams),
        ngram_weight=str(ngram_weight),
    )

    if bool(use_bout_stats):
        adata_state_features = _drop_disabled_bout_features(
            adata_state_features,
            include_bouts_mean_length=bool(include_bouts_mean_length),
            include_bouts_nr=bool(include_bouts_nr),
            include_bouts_max_length=bool(include_bouts_max_length),
        )

    adata_state_features.uns["track_feature_selection"] = {
        "state_col": str(state_col),
        "use_fractions": bool(use_fractions),
        "use_bout_stats": bool(use_bout_stats),
        "include_bouts_mean_length": bool(include_bouts_mean_length),
        "include_bouts_nr": bool(include_bouts_nr),
        "include_bouts_max_length": bool(include_bouts_max_length),
        "use_transitions": bool(use_transitions),
        "use_bigrams": bool(use_bigrams),
        "use_trigrams": bool(use_trigrams),
        "ngram_weight": str(ngram_weight),
    }
    adata_state_features.uns["track_filtering"] = {
        "groupby_cols": [str(c) for c in list(groupby_cols)],
        "time_col": str(time_col),
        "state_col": str(state_col),
        "behavioral_trajectory_size": trajectory_size,
    }
    return adata_state_features


def _align_feature_matrix_to_artifact(adata_tracks, feature_cols):
    expected_cols = [str(c) for c in list(feature_cols)]
    if len(expected_cols) == 0:
        raise ValueError("Classifier artifact has no feature columns.")

    present_cols = [str(c) for c in list(adata_tracks.var_names)]
    X_df = pd.DataFrame(
        _to_numpy_2d(adata_tracks.X).astype(float, copy=False),
        index=adata_tracks.obs.index,
        columns=present_cols,
    )
    missing = [c for c in expected_cols if c not in X_df.columns]
    for c in missing:
        X_df[c] = 0.0
    X_query = X_df[expected_cols].to_numpy(dtype=float, copy=False)
    extra = [c for c in present_cols if c not in set(expected_cols)]
    return X_query, missing, extra


def predict_track_clusterids_with_classifier(
    adata_tracks,
    classifier_artifact_or_path,
    output_col=None,
    confidence_col=None,
    inplace=True,
):
    """Predict track clusters for an AnnData track-feature matrix."""
    artifact = _coerce_track_classifier_artifact(classifier_artifact_or_path)
    target = adata_tracks if bool(inplace) else adata_tracks.copy()
    cluster_col = str(artifact.get("cluster_col", "ClusterID"))
    out_col = cluster_col if output_col is None else str(output_col)
    conf_col = f"{out_col}_confidence" if confidence_col is None else str(confidence_col)
    feature_cols = artifact.get("feature_cols", None)
    if feature_cols is None:
        feature_cols = list(target.var_names)
    clf = artifact["classifier"]

    X_query, missing_features, extra_features = _align_feature_matrix_to_artifact(
        target,
        feature_cols=feature_cols,
    )
    pred = np.asarray(clf.predict(X_query)).astype(str)
    target.obs[out_col] = pd.Categorical(pd.Series(pred, index=target.obs.index).astype(str))

    if hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(X_query)
        target.obs[conf_col] = pd.Series(proba.max(axis=1), index=target.obs.index, dtype=float)

    target.uns["classification_apply"] = {
        "classifier_backend": str(artifact.get("backend", "random_forest")),
        "classifier_cluster_col": cluster_col,
        "prediction_col": out_col,
        "confidence_col": conf_col,
        "n_rows": int(target.n_obs),
        "n_features_expected": int(len(feature_cols)),
        "n_features_missing_filled_zero": int(len(missing_features)),
        "n_features_extra_ignored": int(len(extra_features)),
        "missing_features": [str(c) for c in missing_features],
    }
    return target


def train_track_classifier(
    output_dir,
    cell_type="tcell",
    model_adata=None,
    cluster_col="ClusterID",
    classifier_n_estimators=300,
    classifier_min_samples_leaf=2,
    classifier_n_jobs=-1,
    classifier_max_depth=None,
    classifier_min_samples_split=2,
    classifier_max_features="sqrt",
    classifier_class_weight=None,
    output_subdir_name="behavioral_state_trajectories",
    save_classifier=True,
    classifier_path=None,
    random_state=123,
    validation_test_size=0.05,
    validation_random_state=None,
    validation_stratify=True,
    verbose=True,
):
    """Train a random-forest classifier on track features in model_adata."""
    train_started = _vstart(verbose, "trajectory-training", "train track classifier")
    resolved_paths = _resolve_track_paths(
        output_dir=output_dir,
        cell_type=cell_type,
        output_subdir_name=output_subdir_name,
    )
    outfolder = resolved_paths["outfolder"]
    model_adata_path = outfolder / get_track_trajectories_filename(cell_type)
    if model_adata is None:
        if not model_adata_path.exists():
            raise FileNotFoundError(
                f"Could not load model adata for classifier training: {model_adata_path}. "
                "Run run_state_based_analysis first or pass model_adata explicitly."
            )
        model_adata = sc.read_h5ad(model_adata_path)
        _vinfo(verbose, "trajectory-training", f"using model adata cache: {model_adata_path}")

    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = [str(c) for c in list(model_adata.var_names)]
    X_all = _to_numpy_2d(model_adata[:, feature_cols].X).astype(float, copy=False)
    y_all = model_adata.obs[cluster_col].astype(str).to_numpy()
    if len(y_all) == 0:
        raise ValueError("No rows available to train classifier.")

    stage1_random_state = int(random_state)
    try:
        clustering_meta = model_adata.uns.get("clustering", {})
        if isinstance(clustering_meta, dict) and ("random_state" in clustering_meta):
            stage1_random_state = int(clustering_meta.get("random_state", stage1_random_state))
    except Exception:
        stage1_random_state = int(random_state)

    classifier_backend_name = "random_forest"
    validation_mode = _resolve_validation_mode(validation_test_size)
    if validation_mode == "holdout":
        validation_test_size = float(validation_test_size)
    else:
        validation_test_size = 0.0
    validation_random_state = (
        int(stage1_random_state) if validation_random_state is None else int(validation_random_state)
    )
    validation_stratify = bool(validation_stratify)

    classifier_qc_outfolder = outfolder / "classification" / "qc"
    classifier_qc_outfolder.mkdir(parents=True, exist_ok=True)

    classifier_validation = {
        "mode": str(validation_mode),
        "test_size": None if validation_mode == "oob_only" else float(validation_test_size),
        "random_state": int(validation_random_state),
        "stratify": bool(validation_stratify),
        "target_col": str(cluster_col),
        "n_train": 0,
        "n_test": 0,
        "train_metrics": None,
        "test_metrics": None,
        "oob_score_train_split": None,
        "oob_score_full_data": None,
        "holdout_artifacts": {},
    }

    X_train = X_all
    y_train = y_all
    X_test = None
    y_test = None
    if validation_mode == "holdout":
        if validation_stratify:
            _validate_stratified_split_feasibility(y_all, validation_test_size, target_name=str(cluster_col))
        idx_all = np.arange(len(y_all), dtype=int)
        idx_train, idx_test = train_test_split(
            idx_all,
            test_size=validation_test_size,
            random_state=validation_random_state,
            stratify=(y_all if validation_stratify else None),
        )
        X_train = X_all[idx_train]
        y_train = y_all[idx_train]
        X_test = X_all[idx_test]
        y_test = y_all[idx_test]
        classifier_validation["n_train"] = int(len(idx_train))
        classifier_validation["n_test"] = int(len(idx_test))
    else:
        classifier_validation["n_train"] = int(len(y_all))
        classifier_validation["n_test"] = 0

    clf, fit_info = train_random_forest_classifier(
        X_train=X_train,
        y_train=y_train,
        random_state=int(stage1_random_state),
        n_estimators=int(classifier_n_estimators),
        min_samples_leaf=int(classifier_min_samples_leaf),
        n_jobs=int(classifier_n_jobs),
        max_depth=classifier_max_depth,
        min_samples_split=int(classifier_min_samples_split),
        max_features=classifier_max_features,
        class_weight=classifier_class_weight,
    )
    fit_info["cluster_col"] = str(cluster_col)

    y_pred_train = np.asarray(clf.predict(X_train)).astype(str)
    classifier_validation["train_metrics"] = _build_compact_metric_summary(y_train, y_pred_train)
    classifier_validation["oob_score_train_split"] = float(fit_info.get("oob_score", np.nan))

    if validation_mode == "holdout" and X_test is not None:
        y_pred_test = np.asarray(clf.predict(X_test)).astype(str)
        holdout_eval = _evaluate_holdout_predictions(
            y_true=y_test,
            y_pred=y_pred_test,
            outfolder=classifier_qc_outfolder,
            filename_prefix=f"track_classification_classifier_holdout_qc_{classifier_backend_name}",
        )
        classifier_validation["test_metrics"] = holdout_eval.get("metrics", None)
        classifier_validation["holdout_artifacts"] = holdout_eval.get("artifacts", {})

    X_eval = X_all if validation_mode == "holdout" else X_train
    y_eval = y_all if validation_mode == "holdout" else y_train
    y_pred_eval = np.asarray(clf.predict(X_eval)).astype(str)
    classifier_qc_model_data = evaluate_classifier_predictions(
        y_true=y_eval,
        y_pred=y_pred_eval,
        outfolder=classifier_qc_outfolder,
        filename_prefix=f"track_classification_classifier_qc_{classifier_backend_name}",
    )
    classifier_validation["oob_score_full_data"] = (
        classifier_validation["oob_score_train_split"] if validation_mode != "holdout" else None
    )

    _vinfo(
        verbose,
        "trajectory-training",
        "QC on model labels | "
        f"n={classifier_qc_model_data['n_rows']}, "
        f"accuracy={classifier_qc_model_data['accuracy']:.4f}, "
        f"balanced_accuracy={classifier_qc_model_data['balanced_accuracy']:.4f}",
    )
    if bool(verbose):
        tm = classifier_validation["train_metrics"] or {}
        _vinfo(
            verbose,
            "trajectory-training",
            "training scores | "
            f"n_train={classifier_validation['n_train']}, "
            f"train_acc={tm.get('accuracy', np.nan):.4f}, "
            f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
            f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
            f"oob_train_split={classifier_validation.get('oob_score_train_split', np.nan):.4f}",
        )
        if validation_mode == "holdout":
            vm = classifier_validation["test_metrics"] or {}
            _vinfo(
                verbose,
                "trajectory-training",
                "validation_set scores | "
                f"n_val={classifier_validation['n_test']}, "
                f"val_acc={vm.get('accuracy', np.nan):.4f}, "
                f"val_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                f"val_macro_f1={vm.get('macro_f1', np.nan):.4f}",
            )
        else:
            _vinfo(
                verbose,
                "trajectory-training",
                "validation_set scores: not available (validation_test_size=0, OOB-only mode).",
            )

    tm = classifier_validation.get("train_metrics") or {}
    vm = classifier_validation.get("test_metrics") or {}
    metrics_row = {
        "mode": validation_mode,
        "n_train": classifier_validation.get("n_train", 0),
        "n_test": classifier_validation.get("n_test", 0),
        "oob_score_train_split": classifier_validation.get("oob_score_train_split", np.nan),
        "train_accuracy": tm.get("accuracy", np.nan),
        "train_balanced_accuracy": tm.get("balanced_accuracy", np.nan),
        "train_macro_f1": tm.get("macro_f1", np.nan),
        "test_accuracy": vm.get("accuracy", np.nan),
        "test_balanced_accuracy": vm.get("balanced_accuracy", np.nan),
        "test_macro_f1": vm.get("macro_f1", np.nan),
        "model_data_accuracy": classifier_qc_model_data.get("accuracy", np.nan),
        "model_data_balanced_accuracy": classifier_qc_model_data.get("balanced_accuracy", np.nan),
    }
    metrics_csv = (
        classifier_qc_outfolder / f"track_classification_classifier_metrics_{classifier_backend_name}.csv"
    )
    pd.DataFrame([metrics_row]).to_csv(metrics_csv, index=False)

    filter_cfg = dict(model_adata.uns.get("track_filtering", {}))
    feat_cfg = dict(model_adata.uns.get("track_feature_selection", {}))
    training_state_col = str(
        feat_cfg.get(
            "state_col",
            filter_cfg.get("state_col", "full_behavioral_cluster"),
        )
    )
    preprocessing_spec = _build_track_preprocessing_spec(
        state_col=training_state_col,
        filter_cfg=filter_cfg,
        feat_cfg=feat_cfg,
    )

    artifact = {
        "classifier": clf,
        "backend": classifier_backend_name,
        "feature_cols": feature_cols,
        "cluster_col": str(cluster_col),
        "cell_type": str(cell_type),
        "fit_info": dict(fit_info),
        "validation": dict(classifier_validation),
        "track_filtering": dict(filter_cfg),
        "track_feature_selection": dict(feat_cfg),
        "state_col_used_for_training": str(training_state_col),
        "preprocessing_spec_version": 1,
        "preprocessing_spec": dict(preprocessing_spec),
        "classifier_qc_model_data": classifier_qc_model_data,
        "classifier_metrics_csv": str(metrics_csv),
    }

    save_path = None
    classifier_model_paths = {}
    if bool(save_classifier):
        save_path = (
            Path(classifier_path)
            if classifier_path is not None
            else _resolve_track_classifier_path(
                output_dir=output_dir,
                cell_type=cell_type,
                output_subdir_name=output_subdir_name,
            )
        )
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with save_path.open("wb") as f:
            pickle.dump(artifact, f)
        classifier_model_paths["train_split"] = str(save_path)
        _vsave(verbose, "trajectory-training", "track classifier artifact", save_path)

    model_adata.uns.setdefault("classification", {})
    model_adata.uns["classification"].update(
        {
            "classifier_backend": classifier_backend_name,
            "classifier_cluster_col": str(cluster_col),
            "classifier_fit_info": dict(fit_info),
            "classifier_fit_info_train_split": dict(fit_info),
            "classifier_qc_model_data": classifier_qc_model_data,
            "classifier_qc_dir": str(classifier_qc_outfolder),
            "classifier_metrics_csv": str(metrics_csv),
            "validation_config": {
                "mode": validation_mode,
                "test_size": None if validation_mode == "oob_only" else float(validation_test_size),
                "random_state": int(validation_random_state),
                "stratify": bool(validation_stratify),
                "split_unit": "row",
            },
            "classifier_validation": dict(classifier_validation),
            "classifier_model_path": None if save_path is None else str(save_path),
            "classifier_model_paths": dict(classifier_model_paths),
            "state_col_used_for_training": str(training_state_col),
            "preprocessing_spec_version": 1,
            "preprocessing_spec": dict(preprocessing_spec),
        }
    )
    _vdone(verbose, "trajectory-training", "train track classifier", train_started)
    return {
        "classifier_path": None if save_path is None else str(save_path),
        "fit_info": dict(fit_info),
        "classifier_validation": dict(classifier_validation),
        "classifier_qc_model_data": classifier_qc_model_data,
        "classifier_metrics_csv": str(metrics_csv),
        "classifier_artifact": artifact,
        "model_adata": model_adata,
    }


def apply_track_classifier_to_subtracks(
    output_dir,
    cell_type="tcell",
    classifier_artifact_or_path=None,
    adata_full=None,
    adata_full_path=None,
    state_col="full_behavioral_cluster",
    output_col="ClusterID",
    confidence_col=None,
    output_subdir_name="behavioral_state_trajectories",
    plot_exemplars=True,
    plot_exemplar_backprojection_videos=True,
    plot_exemplar_backprojection_pdfs=False,
    n_per_cluster=10,
    plot_dpi=300,
    plot_track_class_proportions=True,
    track_class_proportion_plot_dpi=300,
    exemplar_video_fps=6,
    exemplar_video_dpi=180,
    exemplar_video_margin=20,
    exemplar_video_pmin=0.0,
    exemplar_video_pmax=99.99,
    exemplar_video_track_color="#63ff33",
    exemplar_backprojection_show_segments=True,
    exemplar_backprojection_segment_style="outline",
    exemplar_backprojection_segment_color="#ffffff",
    exemplar_backprojection_layout_mode="both",
    exemplar_backprojection_examples_per_cluster=5,
    exemplar_backprojection_num_example_ranks=5,
    exemplar_backprojection_pdf_rows_per_page=6,
    exemplar_backprojection_pdf_show_raw_image=True,
    exemplar_backprojection_show_state_legend=True,
    save_outputs=True,
    save_as_model=True,
    random_state=123,
    verbose=True,
):
    """
    Build (sub)track feature rows from the full dataset and apply a trained classifier.
    """
    apply_started = _vstart(
        verbose,
        "trajectory-apply",
        f"apply classifier | cell_type={cell_type} | output_col={output_col}",
    )
    resolved_paths = _resolve_track_paths(
        output_dir=output_dir,
        cell_type=cell_type,
        output_subdir_name=output_subdir_name,
    )
    outfolder = resolved_paths["outfolder"]
    state_folder = resolved_paths["state_outdir"]

    if classifier_artifact_or_path is None:
        classifier_artifact_or_path = _resolve_track_classifier_path(
            output_dir=output_dir,
            cell_type=cell_type,
            output_subdir_name=output_subdir_name,
        )
    artifact = _coerce_track_classifier_artifact(classifier_artifact_or_path)
    _vinfo(
        verbose,
        "trajectory-apply",
        "loaded classifier artifact | "
        f"backend={artifact.get('backend', 'unknown')} | cluster_col={artifact.get('cluster_col', 'ClusterID')}",
    )

    filter_cfg = dict(artifact.get("track_filtering", {}))
    feat_cfg = dict(artifact.get("track_feature_selection", {}))
    artifact_preprocessing_spec_version = int(artifact.get("preprocessing_spec_version", 1))
    state_col_mismatch_error = None
    groupby_cols = filter_cfg.get("groupby_cols", ("sample_name", "TrackID"))
    if not isinstance(groupby_cols, (list, tuple)) or len(groupby_cols) == 0:
        groupby_cols = ("sample_name", "TrackID")
    time_col = str(filter_cfg.get("time_col", "position_t"))
    resolved_state_col, _ = _resolve_apply_state_col(
        artifact=artifact,
        requested_state_col=state_col,
    )
    if str(artifact.get("state_col_used_for_training", "")).strip() == "":
        _vinfo(
            verbose,
            "trajectory-apply",
            "legacy classifier artifact without 'state_col_used_for_training'; using requested state_col.",
        )
    if "behavioral_trajectory_size" not in filter_cfg:
        raise ValueError(
            "Classifier artifact is missing required 'track_filtering.behavioral_trajectory_size'. "
            "This appears to be a legacy model. Regenerate the behavioral state trajectory model."
        )
    behavioral_trajectory_size = int(filter_cfg["behavioral_trajectory_size"])

    apply_preprocessing_spec = dict(
        artifact.get(
            "preprocessing_spec",
            _build_track_preprocessing_spec(
                state_col=resolved_state_col,
                filter_cfg=filter_cfg,
                feat_cfg=feat_cfg,
            ),
        )
    )

    load_started = time.perf_counter()
    adata_full, apply_input_source_kind, apply_input_source_path = _resolve_track_apply_input_adata(
        adata_full=adata_full,
        adata_full_path=adata_full_path,
        state_folder=state_folder,
        cell_type=cell_type,
    )
    required_cols = set([str(time_col), str(resolved_state_col)] + [str(c) for c in list(groupby_cols)])
    missing_required_cols = sorted([c for c in required_cols if c not in adata_full.obs.columns])
    if len(missing_required_cols) > 0:
        raise ValueError(
            "Behavioral-state AnnData is missing required columns for classifier apply: "
            f"{missing_required_cols}."
        )
    _vinfo(
        verbose,
        "trajectory-apply",
        "loaded behavioral-state dataset | "
        f"source={apply_input_source_kind} | path={apply_input_source_path}",
    )
    _vdone(verbose, "trajectory-apply", "load behavioral-state dataset", load_started)

    use_bout_stats = bool(feat_cfg.get("use_bout_stats", True))
    predict_started = time.perf_counter()
    _vinfo(
        verbose,
        "trajectory-apply",
        "building trajectory feature matrix | "
        f"state_col={resolved_state_col} | behavioral_trajectory_size={behavioral_trajectory_size}",
    )
    adata_tracks = _build_track_feature_adata(
        adata_full,
        state_col=str(resolved_state_col),
        groupby_cols=groupby_cols,
        time_col=time_col,
        behavioral_trajectory_size=behavioral_trajectory_size,
        use_fractions=bool(feat_cfg.get("use_fractions", True)),
        use_bout_stats=bool(use_bout_stats),
        include_bouts_mean_length=bool(feat_cfg.get("include_bouts_mean_length", True)),
        include_bouts_nr=bool(feat_cfg.get("include_bouts_nr", True)),
        include_bouts_max_length=bool(feat_cfg.get("include_bouts_max_length", True)),
        use_transitions=bool(feat_cfg.get("use_transitions", True)),
        use_bigrams=bool(feat_cfg.get("use_bigrams", True)),
        use_trigrams=bool(feat_cfg.get("use_trigrams", True)),
        ngram_weight=str(feat_cfg.get("ngram_weight", "count")),
    )
    adata_tracks = predict_track_clusterids_with_classifier(
        adata_tracks=adata_tracks,
        classifier_artifact_or_path=artifact,
        output_col=str(output_col),
        confidence_col=confidence_col,
        inplace=True,
    )
    _vinfo(
        verbose,
        "trajectory-apply",
        "predicted clusters | "
        f"rows={adata_tracks.n_obs} | unique_labels={adata_tracks.obs[str(output_col)].astype(str).nunique()}",
    )
    _vdone(verbose, "trajectory-apply", "build features + predict", predict_started)
    track_class_proportions_plot_pdf = None
    track_class_proportions_table_csv = None
    track_class_proportions_sample_order = []
    track_class_proportions_class_order = []
    track_class_proportions_colors = {}
    track_class_proportions_error = None
    track_class_proportions_config = {
        "enabled": bool(plot_track_class_proportions),
        "sample_col": "sample_name",
        "class_col": str(output_col),
        "dpi": int(track_class_proportion_plot_dpi),
        "cmap_name": "tab20",
    }
    if bool(plot_track_class_proportions):
        proportions_started = time.perf_counter()
        try:
            class_outdir = outfolder / "classification"
            class_outdir.mkdir(parents=True, exist_ok=True)
            proportions_out = save_track_class_proportions_by_sample_plot(
                adata_tracks=adata_tracks,
                out_dir=class_outdir,
                sample_col="sample_name",
                class_col=str(output_col),
                dpi=int(track_class_proportion_plot_dpi),
                cmap_name="tab20",
            )
            track_class_proportions_plot_pdf = proportions_out.get("pdf_path")
            track_class_proportions_table_csv = proportions_out.get("csv_path")
            track_class_proportions_sample_order = list(
                proportions_out.get("sample_order", [])
            )
            track_class_proportions_class_order = list(
                proportions_out.get("class_order", [])
            )
            track_class_proportions_colors = dict(proportions_out.get("colors", {}))
            _vsave(verbose, "trajectory-apply", "track class proportions PDF", track_class_proportions_plot_pdf)
            _vsave(verbose, "trajectory-apply", "track class proportions CSV", track_class_proportions_table_csv)
            _vdone(verbose, "trajectory-apply", "write track class proportions outputs", proportions_started)
        except Exception as exc:
            track_class_proportions_error = str(exc)
            _vinfo(verbose, "trajectory-apply", f"skipping track class proportions export due to error: {exc}")

    exemplar_statebar_track_pdf_by_cluster = {}
    exemplar_statebar_track_pdf_by_example_rank = {}
    exemplar_backprojection_pdf_paths_by_cluster = {}
    exemplar_backprojection_pdf_paths_by_example_rank = {}
    exemplar_backprojection_video_paths_by_cluster = {}
    exemplar_backprojection_video_paths_by_example_rank = {}
    exemplar_backprojection_global_canvas_xy_shape = None
    exemplar_selection_csv = None
    exemplar_render_config = {
        "stage": "after_classification",
        "cluster_key": str(output_col),
        "state_key": str(resolved_state_col),
        "n_per_cluster": int(n_per_cluster),
        "plot_exemplars": bool(plot_exemplars),
        "plot_exemplar_backprojection_videos": bool(plot_exemplar_backprojection_videos),
        "plot_exemplar_backprojection_pdfs": bool(plot_exemplar_backprojection_pdfs),
        "plot_dpi": int(plot_dpi),
        "video_fps": int(exemplar_video_fps),
        "video_dpi": int(exemplar_video_dpi),
        "video_margin": int(exemplar_video_margin),
        "video_pmin": float(exemplar_video_pmin),
        "video_pmax": float(exemplar_video_pmax),
        "video_track_color": str(exemplar_video_track_color),
        "backprojection_show_segments": bool(exemplar_backprojection_show_segments),
        "backprojection_segment_style": str(exemplar_backprojection_segment_style),
        "backprojection_segment_color": str(exemplar_backprojection_segment_color),
        "backprojection_layout_mode": str(exemplar_backprojection_layout_mode),
        "backprojection_examples_per_cluster": int(exemplar_backprojection_examples_per_cluster),
        "backprojection_num_example_ranks": int(exemplar_backprojection_num_example_ranks),
        "backprojection_pdf_rows_per_page": int(exemplar_backprojection_pdf_rows_per_page),
        "backprojection_pdf_show_raw_image": bool(exemplar_backprojection_pdf_show_raw_image),
        "backprojection_show_state_legend": bool(exemplar_backprojection_show_state_legend),
        "coordinate_enrichment": None,
    }

    if bool(plot_exemplars):
        exemplar_started = time.perf_counter()
        _vinfo(
            verbose,
            "trajectory-apply",
            "rendering exemplar visualizations | "
            f"backprojection_videos={bool(plot_exemplar_backprojection_videos)} | "
            f"backprojection_pdfs={bool(plot_exemplar_backprojection_pdfs)}",
        )
        exemplar_root = outfolder / "classification" / "example_tracks"
        exemplar_root.mkdir(parents=True, exist_ok=True)

        adata_plot = filter_and_truncate_tracks_anndata(
            adata_full,
            groupby_cols=[str(c) for c in list(groupby_cols)],
            time_col=str(time_col),
            min_length=int(behavioral_trajectory_size),
            max_length=int(behavioral_trajectory_size),
        )
        coord_enrichment_pdf = _ensure_exemplar_coordinate_columns(
            adata=adata_plot,
            output_dir=output_dir,
            cell_type=cell_type,
            require_pixel_for_video=False,
        )
        exemplar_render_config["coordinate_enrichment"] = dict(coord_enrichment_pdf)

        chosen_exemplars, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            cluster_key=str(output_col),
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            seed=int(random_state),
            query=None,
        )

        exemplar_selection_csv = exemplar_root / (
            f"example_track_selection_cluster_{_sanitize_filename_token(output_col, fallback='cluster')}_"
            f"state_{_sanitize_filename_token(resolved_state_col, fallback='state')}_classified.csv"
        )
        chosen_exemplars.to_csv(exemplar_selection_csv, index=False)
        _vsave(verbose, "trajectory-apply", "example track selection CSV", exemplar_selection_csv)

        per_cluster_pdf_out = save_exemplar_statebar_track_pdf_per_cluster(
            adata_full=adata_plot,
            out_dir=exemplar_root,
            chosen_df=chosen_exemplars,
            adata_tracks=None,
            n_per_cluster=int(n_per_cluster),
            sample_key="sample_name",
            track_key="TrackID",
            time_key=str(time_col),
            state_key=str(resolved_state_col),
            cluster_key=str(output_col),
            tmin_key="position_t_min",
            tmax_key="position_t_max",
            rows_per_page=6,
            plot_dpi=int(plot_dpi),
            seed=int(random_state),
            cmap_name="tab20",
            layout_mode=str(exemplar_backprojection_layout_mode),
            num_example_ranks=int(exemplar_backprojection_num_example_ranks),
        )
        exemplar_statebar_track_pdf_by_cluster = dict(
            per_cluster_pdf_out.get("pdf_paths_by_cluster", {})
        )
        exemplar_statebar_track_pdf_by_example_rank = dict(
            per_cluster_pdf_out.get("pdf_paths_by_example_rank", {})
        )

        if bool(plot_exemplar_backprojection_videos) or bool(plot_exemplar_backprojection_pdfs):
            coord_enrichment_video = _ensure_exemplar_coordinate_columns(
                adata=adata_plot,
                output_dir=output_dir,
                cell_type=cell_type,
                require_pixel_for_video=True,
            )
            exemplar_render_config["coordinate_enrichment_video"] = dict(coord_enrichment_video)
            exemplar_backprojection_outdir = exemplar_root
            exemplar_backprojection_outdir.mkdir(parents=True, exist_ok=True)
            prep_t0 = time.perf_counter()
            exemplar_backprojection_prepared = _prepare_exemplar_backprojection_data(
                adata_full=adata_plot,
                output_dir=output_dir,
                cell_type=cell_type,
                chosen_df=chosen_exemplars.copy().reset_index(drop=True),
                sample_key="sample_name",
                track_key="TrackID",
                time_key=str(time_col),
                state_key=str(resolved_state_col),
                cluster_key=str(output_col),
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                margin=int(exemplar_video_margin),
                pmin=float(exemplar_video_pmin),
                pmax=float(exemplar_video_pmax),
                cmap_name="tab20",
                coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                show_segment_outlines=bool(exemplar_backprojection_show_segments),
                segment_style=str(exemplar_backprojection_segment_style),
                segment_color=str(exemplar_backprojection_segment_color),
            )
            exemplar_render_config["segment_outline_errors"] = dict(
                exemplar_backprojection_prepared.get("segment_outline_errors", {})
            )
            _vdone(verbose, "trajectory-apply", "prepare backprojection data", prep_t0)

            if bool(plot_exemplar_backprojection_videos):
                backprojection_video_out = save_exemplar_statebar_backprojection_video_per_cluster(
                    adata_full=adata_plot,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    out_dir=exemplar_backprojection_outdir,
                    chosen_df=chosen_exemplars,
                    adata_tracks=None,
                    n_per_cluster=int(n_per_cluster),
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key=str(time_col),
                    state_key=str(resolved_state_col),
                    cluster_key=str(output_col),
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    fps=int(exemplar_video_fps),
                    dpi=int(exemplar_video_dpi),
                    margin=int(exemplar_video_margin),
                    pmin=float(exemplar_video_pmin),
                    pmax=float(exemplar_video_pmax),
                    track_color=str(exemplar_video_track_color),
                    show_segment_outlines=bool(exemplar_backprojection_show_segments),
                    segment_style=str(exemplar_backprojection_segment_style),
                    segment_color=str(exemplar_backprojection_segment_color),
                    coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                    seed=int(random_state),
                    cmap_name="tab20",
                    layout_mode=str(exemplar_backprojection_layout_mode),
                    examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                    num_example_ranks=int(exemplar_backprojection_num_example_ranks),
                    show_state_legend=bool(exemplar_backprojection_show_state_legend),
                    ffmpeg_preset="veryfast",
                    ffmpeg_crf=21,
                    prepared_data=exemplar_backprojection_prepared,
                    verbose=bool(verbose),
                )
                exemplar_backprojection_video_paths_by_cluster = dict(
                    backprojection_video_out.get("video_paths_by_cluster", {})
                )
                exemplar_backprojection_video_paths_by_example_rank = dict(
                    backprojection_video_out.get("video_paths_by_example_rank", {})
                )
                exemplar_backprojection_global_canvas_xy_shape = backprojection_video_out.get(
                    "global_xy_shape", exemplar_backprojection_global_canvas_xy_shape
                )

            if bool(plot_exemplar_backprojection_pdfs):
                backprojection_pdf_out = save_exemplar_statebar_backprojection_pdf(
                    adata_full=adata_plot,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    out_dir=exemplar_backprojection_outdir,
                    chosen_df=chosen_exemplars,
                    adata_tracks=None,
                    n_per_cluster=int(n_per_cluster),
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key=str(time_col),
                    state_key=str(resolved_state_col),
                    cluster_key=str(output_col),
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    rows_per_page=int(exemplar_backprojection_pdf_rows_per_page),
                    plot_dpi=max(120, int(exemplar_video_dpi)),
                    margin=int(exemplar_video_margin),
                    pmin=float(exemplar_video_pmin),
                    pmax=float(exemplar_video_pmax),
                    track_color=str(exemplar_video_track_color),
                    show_segment_outlines=bool(exemplar_backprojection_show_segments),
                    segment_style=str(exemplar_backprojection_segment_style),
                    segment_color=str(exemplar_backprojection_segment_color),
                    coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                    seed=int(random_state),
                    cmap_name="tab20",
                    layout_mode=str(exemplar_backprojection_layout_mode),
                    examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                    num_example_ranks=int(exemplar_backprojection_num_example_ranks),
                    show_raw_image=bool(exemplar_backprojection_pdf_show_raw_image),
                    show_state_legend=bool(exemplar_backprojection_show_state_legend),
                    prepared_data=exemplar_backprojection_prepared,
                    verbose=bool(verbose),
                )
                exemplar_backprojection_pdf_paths_by_cluster = dict(
                    backprojection_pdf_out.get("pdf_paths_by_cluster", {})
                )
                exemplar_backprojection_pdf_paths_by_example_rank = dict(
                    backprojection_pdf_out.get("pdf_paths_by_example_rank", {})
                )
                exemplar_backprojection_global_canvas_xy_shape = backprojection_pdf_out.get(
                    "global_xy_shape", exemplar_backprojection_global_canvas_xy_shape
                )
                exemplar_render_config["segment_outline_errors"] = dict(
                    backprojection_pdf_out.get(
                        "segment_outline_errors",
                        exemplar_render_config.get("segment_outline_errors", {}),
                    )
                )
        _vdone(verbose, "trajectory-apply", "render exemplar outputs", exemplar_started)

    output_path = None
    if bool(save_outputs):
        if bool(save_as_model):
            output_path = outfolder / get_track_trajectories_filename(cell_type)
        else:
            output_path = outfolder / f"BEHAV3D_{_sanitize_filename_token(cell_type, fallback='cell')}_behavioral_trajectories_predicted.h5ad"
        output_path.parent.mkdir(parents=True, exist_ok=True)

    adata_tracks.uns.setdefault("classification", {})
    adata_tracks.uns["classification"].update(
        {
            "apply_random_state": int(random_state),
            "classifier_backend": str(artifact.get("backend", "random_forest")),
            "classifier_cluster_col": str(artifact.get("cluster_col", "ClusterID")),
            "classifier_model_path": (
                str(classifier_artifact_or_path)
                if isinstance(classifier_artifact_or_path, (str, Path))
                else None
            ),
            "output_col": str(output_col),
            "confidence_col": (
                f"{str(output_col)}_confidence" if confidence_col is None else str(confidence_col)
            ),
            "applied_on_adata_full_path": None if apply_input_source_path is None else str(apply_input_source_path),
            "apply_input_source_kind": str(apply_input_source_kind),
            "apply_input_source_path": None if apply_input_source_path is None else str(apply_input_source_path),
            "apply_preprocessing_spec": dict(apply_preprocessing_spec),
            "apply_preprocessing_spec_version": int(artifact_preprocessing_spec_version),
            "apply_state_col_resolved": str(resolved_state_col),
            "apply_state_col_mismatch_error": state_col_mismatch_error,
            "saved_output_path": None if output_path is None else str(output_path),
            "track_class_proportions_plot_pdf": track_class_proportions_plot_pdf,
            "track_class_proportions_table_csv": track_class_proportions_table_csv,
            "track_class_proportions_sample_order": list(track_class_proportions_sample_order),
            "track_class_proportions_class_order": list(track_class_proportions_class_order),
            "track_class_proportions_colors": dict(track_class_proportions_colors),
            "track_class_proportions_error": track_class_proportions_error,
            "track_class_proportions_config": dict(track_class_proportions_config),
            "exemplar_render_stage": "after_classification",
            "exemplar_statebar_track_pdf_by_cluster": dict(exemplar_statebar_track_pdf_by_cluster),
            "exemplar_statebar_track_pdf_by_example_rank": dict(
                exemplar_statebar_track_pdf_by_example_rank
            ),
            "exemplar_statebar_backprojection_video_by_cluster": dict(
                exemplar_backprojection_video_paths_by_cluster
            ),
            "exemplar_backprojection_pdf_paths_by_cluster": dict(
                exemplar_backprojection_pdf_paths_by_cluster
            ),
            "exemplar_backprojection_pdf_paths_by_example_rank": dict(
                exemplar_backprojection_pdf_paths_by_example_rank
            ),
            "exemplar_backprojection_video_paths_by_cluster": dict(
                exemplar_backprojection_video_paths_by_cluster
            ),
            "exemplar_backprojection_video_paths_by_example_rank": dict(
                exemplar_backprojection_video_paths_by_example_rank
            ),
            "exemplar_backprojection_global_canvas_xy_shape": exemplar_backprojection_global_canvas_xy_shape,
            "exemplar_selection_csv": (
                None if exemplar_selection_csv is None else str(exemplar_selection_csv)
            ),
            "exemplar_segment_outline_errors": dict(
                exemplar_render_config.get("segment_outline_errors", {})
            ),
            "exemplar_render_config": dict(exemplar_render_config),
        }
    )
    if bool(save_outputs) and output_path is not None:
        save_started = time.perf_counter()
        adata_tracks.write(output_path, compression="gzip")
        _vsave(verbose, "trajectory-apply", "classifier-applied trajectory model", output_path)
        _vdone(verbose, "trajectory-apply", "save classifier outputs", save_started)
    _vdone(verbose, "trajectory-apply", "apply classifier", apply_started)
    return {
        "adata_tracks": adata_tracks,
        "output_path": None if output_path is None else str(output_path),
        "track_class_proportions_plot_pdf": track_class_proportions_plot_pdf,
        "track_class_proportions_table_csv": track_class_proportions_table_csv,
        "apply_input_source_kind": str(apply_input_source_kind),
        "apply_input_source_path": None if apply_input_source_path is None else str(apply_input_source_path),
        "apply_preprocessing_spec": dict(apply_preprocessing_spec),
        "apply_state_col_resolved": str(resolved_state_col),
    }


def archive_track_clustering_pdfs(outfolder, archive_dir_name="clustering_originals"):
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)

    pdf_paths = sorted([p for p in outfolder.glob("*.pdf") if p.is_file()])
    if len(pdf_paths) == 0:
        return {"archive_dir": None, "archived_paths": []}

    ts = time.strftime("%Y%m%d_%H%M%S")
    archive_dir = outfolder / str(archive_dir_name) / ts
    archive_dir.mkdir(parents=True, exist_ok=True)

    archived_paths = []
    for src in pdf_paths:
        dst = archive_dir / src.name
        suffix_idx = 1
        while dst.exists():
            dst = archive_dir / f"{src.stem}_{suffix_idx}{src.suffix}"
            suffix_idx += 1
        src.replace(dst)
        archived_paths.append(dst)

    return {"archive_dir": archive_dir, "archived_paths": archived_paths}


def _apply_best_pdf_orientation(fig, default_orientation="landscape"):
    """
    Keep figure orientation matched to content shape before PDF export.
    """
    if fig is None:
        return None
    default_orientation = str(default_orientation).strip().lower()
    width, height = fig.get_size_inches()
    if abs(width - height) < 0.05:
        if default_orientation == "portrait":
            fig.set_size_inches(height, max(width, height * 1.25), forward=True)
        else:
            fig.set_size_inches(max(width, height * 1.25), height, forward=True)
    width, height = fig.get_size_inches()
    return "landscape" if float(width) >= float(height) else "portrait"


def _rank_cluster_counts(adata_tracks, cluster_key):
    if cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{cluster_key}' in adata_tracks.obs.")

    labels = (
        pd.Series(adata_tracks.obs[cluster_key])
        .astype("string")
        .dropna()
        .astype(str)
    )
    counts = labels.value_counts(dropna=False).to_dict()
    ranked = sorted(
        [(str(label), int(n)) for label, n in counts.items()],
        key=lambda x: (-int(x[1]), _mixed_label_sort_key(x[0])),
    )
    return ranked


def _plot_cluster_rankings(adata_tracks, cluster_key):
    ranked = _rank_cluster_counts(adata_tracks=adata_tracks, cluster_key=cluster_key)
    cluster_labels = [x[0] for x in ranked]
    cluster_counts = [x[1] for x in ranked]

    fig_h = max(4.5, 0.45 * len(cluster_labels) + 2.0)
    fig, ax = plt.subplots(figsize=(10.5, fig_h))
    y = np.arange(len(cluster_labels))
    bars = ax.barh(y, cluster_counts, color="#2E6FBA")
    ax.set_yticks(y)
    ax.set_yticklabels(cluster_labels)
    ax.invert_yaxis()
    ax.set_xlabel("N tracks")
    ax.set_ylabel(str(cluster_key))
    ax.set_title(f"Cluster occurrence ranking ({cluster_key})")
    ax.grid(axis="x", alpha=0.2)
    for bar, n in zip(bars, cluster_counts):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2.0,
            f" {int(n)}",
            va="center",
            ha="left",
            fontsize=9,
        )
    fig.tight_layout()
    return fig


def _build_label_color_map(labels, cmap_name="tab20"):
    """Create deterministic label->color mapping for stacked plots."""
    labels = [str(x) for x in list(labels)]
    if len(labels) == 0:
        return {}
    cmap = plt.get_cmap(cmap_name)
    if len(labels) <= getattr(cmap, "N", 256):
        denom = max(len(labels) - 1, 1)
        color_values = [cmap(i / denom) for i in range(len(labels))]
    else:
        hsv = plt.get_cmap("hsv")
        color_values = [hsv(i / max(len(labels), 1)) for i in range(len(labels))]
    return {lab: color_values[i] for i, lab in enumerate(labels)}


def save_track_class_proportions_by_sample_plot(
    adata_tracks,
    out_dir,
    *,
    sample_col="sample_name",
    class_col="ClusterID",
    dpi=300,
    cmap_name="tab20",
):
    """
    Save one horizontal stacked bar per sample showing track-class proportions.
    """
    if sample_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{sample_col}' in adata_tracks.obs.")
    if class_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{class_col}' in adata_tracks.obs.")

    plot_df = adata_tracks.obs[[sample_col, class_col]].copy()
    plot_df[sample_col] = (
        plot_df[sample_col]
        .astype("string")
        .fillna("unassigned")
        .astype(str)
        .str.strip()
        .replace("", "unassigned")
    )
    plot_df[class_col] = (
        plot_df[class_col]
        .astype("string")
        .fillna("unassigned")
        .astype(str)
        .str.strip()
        .replace("", "unassigned")
    )
    if len(plot_df) == 0:
        raise ValueError("No rows available to plot track class proportions.")

    sample_order = (
        pd.Series(plot_df[sample_col], dtype="string")
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    class_order = sorted(
        pd.Series(plot_df[class_col], dtype="string")
        .dropna()
        .astype(str)
        .unique()
        .tolist(),
        key=_mixed_label_sort_key,
    )
    if len(sample_order) == 0:
        raise ValueError("No samples available to plot track class proportions.")
    if len(class_order) == 0:
        raise ValueError("No classes available to plot track class proportions.")

    proportion_table = pd.crosstab(
        plot_df[sample_col],
        plot_df[class_col],
        normalize="index",
        dropna=False,
    ).reindex(index=sample_order, columns=class_order, fill_value=0.0)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    class_token = _sanitize_filename_token(class_col, fallback="ClusterID")
    pdf_path = out_dir / f"track_class_proportions_by_sample_{class_token}.pdf"
    csv_path = pdf_path.with_suffix(".csv")
    proportion_table.to_csv(csv_path)

    colors = _build_label_color_map(class_order, cmap_name=cmap_name)
    fig_h = max(2.0, 1.2 * len(sample_order))
    fig, axes = plt.subplots(
        nrows=len(sample_order),
        ncols=1,
        figsize=(10.0, fig_h),
        squeeze=False,
    )
    axes = axes.flatten()
    for i, sample_name in enumerate(sample_order):
        ax = axes[i]
        left = 0.0
        for class_name in class_order:
            val = float(proportion_table.loc[sample_name, class_name])
            if val <= 0.0:
                continue
            ax.barh(
                [0.0],
                [val],
                left=left,
                color=colors[class_name],
                height=0.92,
                edgecolor="none",
                linewidth=0.0,
            )
            left += val
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(str(sample_name), loc="left", fontsize=10, pad=2.0)

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            color=colors[class_name],
            label=str(class_name),
            markersize=7,
        )
        for class_name in class_order
    ]
    fig.legend(
        handles=legend_handles,
        labels=[str(class_name) for class_name in class_order],
        title=str(class_col),
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0, 0.82, 1.0])
    fig.savefig(pdf_path, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)

    color_hex = {str(k): str(to_hex(v)) for k, v in colors.items()}
    return {
        "pdf_path": str(pdf_path),
        "csv_path": str(csv_path),
        "sample_order": [str(x) for x in sample_order],
        "class_order": [str(x) for x in class_order],
        "colors": color_hex,
    }


def generate_track_clustering_report_pdfs(
    adata_tracks,
    outfolder,
    cluster_key="ClusterID",
    heatmap_figsize=(25, 20),
    matrixplot_figsize=(20, 40),
    umap_size=1,
    umap_alpha=0.5,
    plot_dpi=300,
    filename_suffix="",
    verbose=False,
):
    """
    Write track clustering diagnostics as one multi-page PDF:
    - page 1: UMAP + correlation matrix (side-by-side)
    - page 2: heatmap
    - page 3: matrixplot
    - page 4: cluster occurrence ranking
    """
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)

    if cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{cluster_key}' in adata_tracks.obs.")

    if "X_umap" not in adata_tracks.obsm:
        sc.tl.umap(adata_tracks, random_state=0)

    sc.tl.dendrogram(adata_tracks, groupby=cluster_key)

    suffix = str(filename_suffix).strip()
    diagnostics_pdf = Path(
        outfolder,
        f"clustering_diagnostics_{_sanitize_filename_token(cluster_key)}{suffix}.pdf",
    )
    diagnostics_started = _vstart(verbose, "trajectory-clustering", "write clustering diagnostics PDF")

    with PdfPages(diagnostics_pdf) as pdf:
        fig, (ax_umap, ax_corr) = plt.subplots(
            1,
            2,
            figsize=(14.5, 8.5),
            gridspec_kw={"width_ratios": [1.25, 0.75]},
        )
        sc.pl.umap(
            adata_tracks,
            color=cluster_key,
            size=umap_size,
            alpha=umap_alpha,
            ax=ax_umap,
            show=False,
            title="UMAP",
        )
        ax_umap.set_aspect("equal", adjustable="box")

        sc.pl.correlation_matrix(
            adata_tracks,
            groupby=cluster_key,
            show_correlation_numbers=True,
            ax=ax_corr,
            show=False,
        )
        ax_corr.set_title("Cluster correlations", fontsize=10)
        fig.suptitle("Track clustering diagnostics | UMAP + Correlation", fontsize=12, fontweight="bold", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        umap_bbox = ax_umap.get_position().frozen()
        bbox = ax_corr.get_position()
        extra_gap = 0.03
        new_x = min(0.98 - bbox.width, bbox.x0 + extra_gap)
        height_shrink = 0.5
        new_h = bbox.height * height_shrink
        ax_corr.set_position([new_x, bbox.y0 + (bbox.height - new_h) / 2, bbox.width, new_h])
        ax_umap.set_position(umap_bbox)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        sc.pl.heatmap(
            adata_tracks,
            var_names=adata_tracks.var_names,
            groupby=cluster_key,
            standard_scale="var",
            figsize=heatmap_figsize,
            swap_axes=True,
            dendrogram=True,
            show_gene_labels=True,
            show=False,
        )
        fig = plt.gcf()
        default_orientation = "landscape" if float(heatmap_figsize[0]) >= float(heatmap_figsize[1]) else "portrait"
        _apply_best_pdf_orientation(fig, default_orientation=default_orientation)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        sc.pl.matrixplot(
            adata_tracks,
            var_names=adata_tracks.var_names,
            groupby=cluster_key,
            standard_scale="var",
            figsize=matrixplot_figsize,
            swap_axes=True,
            dendrogram=True,
            show=False,
        )
        fig = plt.gcf()
        default_orientation = "landscape" if float(matrixplot_figsize[0]) >= float(matrixplot_figsize[1]) else "portrait"
        _apply_best_pdf_orientation(fig, default_orientation=default_orientation)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        fig = _plot_cluster_rankings(adata_tracks=adata_tracks, cluster_key=cluster_key)
        _apply_best_pdf_orientation(fig, default_orientation="landscape")
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

    _vsave(verbose, "trajectory-clustering", "clustering diagnostics PDF", diagnostics_pdf)
    _vdone(verbose, "trajectory-clustering", "write clustering diagnostics PDF", diagnostics_started)
    return {
        "diagnostics_pdf": diagnostics_pdf,
    }


def build_identity_cluster_mapping(adata, cluster_col="ClusterID"):
    """Build identity mapping from unique labels in `cluster_col`."""
    if not hasattr(adata, "obs"):
        raise ValueError("adata must have an .obs attribute.")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in adata.obs.")

    labels = (
        pd.Series(adata.obs[cluster_col])
        .astype("string")
        .dropna()
        .unique()
        .tolist()
    )
    labels = sorted([str(x) for x in labels], key=_mixed_label_sort_key)
    return {label: label for label in labels}


def rename_track_clusters(
    adata,
    mapping,
    cluster_col="ClusterID",
    keep_unmapped=True,
):
    """
    Relabel track clusters in-place using the same behavior as state classification.
    """
    relabel_cluster_ids(
        adata=adata,
        mapping=mapping,
        cluster_key=cluster_col,
        new_key=cluster_col,
        keep_unmapped=keep_unmapped,
        overwrite_original=True,
    )
    return adata


def assign_track_clusters_to_full_dataset(
    adata_full,
    adata_tracks,
    cluster_col="ClusterID",
    output_col="track_behavioral_cluster",
    id_cols=("sample_name", "TrackID"),
    unassigned_label="unassigned",
    inplace=False,
):
    """
    Broadcast track-level cluster labels to all rows in the full per-timepoint dataset.
    """
    if adata_full is None or not hasattr(adata_full, "obs"):
        raise ValueError("adata_full with .obs is required.")
    if adata_tracks is None or not hasattr(adata_tracks, "obs"):
        raise ValueError("adata_tracks with .obs is required.")

    required_full = [str(c) for c in id_cols]
    missing_full = [c for c in required_full if c not in adata_full.obs.columns]
    if len(missing_full) > 0:
        raise ValueError(f"adata_full.obs missing required id columns: {missing_full}")

    required_tracks = [str(c) for c in id_cols] + [str(cluster_col)]
    missing_tracks = [c for c in required_tracks if c not in adata_tracks.obs.columns]
    if len(missing_tracks) > 0:
        raise ValueError(
            "adata_tracks.obs missing required columns for track-label assignment: "
            f"{missing_tracks}"
        )

    adata_out = adata_full if bool(inplace) else adata_full.copy()
    id_cols = [str(c) for c in id_cols]

    track_labels = (
        adata_tracks.obs[id_cols + [cluster_col]]
        .copy()
        .dropna(subset=[cluster_col])
        .groupby(id_cols, observed=True, as_index=False)[cluster_col]
        .first()
    )

    full_obs = adata_out.obs[id_cols].copy()
    merged = full_obs.merge(track_labels, on=id_cols, how="left")
    labels = pd.Series(merged[cluster_col], index=adata_out.obs.index).astype("string")
    labels = labels.fillna(str(unassigned_label)).astype(str)
    adata_out.obs[output_col] = pd.Categorical(labels)
    return adata_out


def export_track_cluster_backprojection(
    adata_full,
    adata_tracks,
    output_dir,
    cell_type,
    cluster_col="ClusterID",
    output_col="track_behavioral_cluster",
    sample_name=None,
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    verbose=True,
):
    """
    Export track-cluster backprojections by assigning each predicted track label to the
    selected per-track window (typically the last N timepoints used for feature extraction),
    then matching by (sample_name, TrackID, position_t) during voxel mapping.
    """
    required_full = [str(sample_col), str(track_col), str(time_col)]
    missing_full = [c for c in required_full if c not in adata_full.obs.columns]
    if len(missing_full) > 0:
        raise ValueError(
            "Cannot export track-cluster backprojection: adata_full.obs missing required columns "
            f"{missing_full}"
        )

    required_tracks = [str(sample_col), str(track_col), str(cluster_col)]
    missing_tracks = [c for c in required_tracks if c not in adata_tracks.obs.columns]
    if len(missing_tracks) > 0:
        raise ValueError(
            "Cannot export track-cluster backprojection: adata_tracks.obs missing required columns "
            f"{missing_tracks}"
        )

    tracks_obs = adata_tracks.obs.copy()
    full_obs = adata_full.obs[[str(sample_col), str(track_col), str(time_col)]].copy()

    tracks_obs["__sample"] = tracks_obs[str(sample_col)].astype("string").str.strip()
    tracks_obs["__track"] = pd.to_numeric(tracks_obs[str(track_col)], errors="coerce")
    tracks_obs["__cluster"] = tracks_obs[str(cluster_col)].astype("string").str.strip()
    tracks_obs = tracks_obs.dropna(subset=["__sample", "__track", "__cluster"]).copy()
    tracks_obs = tracks_obs[(tracks_obs["__sample"] != "") & (tracks_obs["__cluster"] != "")].copy()

    if sample_name is not None and len(str(sample_name).strip()) > 0:
        selected = str(sample_name).strip()
        tracks_obs = tracks_obs[tracks_obs["__sample"] == selected].copy()
        if len(tracks_obs) == 0:
            raise ValueError(
                f"Could not export track-cluster backprojection: sample '{selected}' not found in adata_tracks."
            )

    full_obs["__sample"] = full_obs[str(sample_col)].astype("string").str.strip()
    full_obs["__track"] = pd.to_numeric(full_obs[str(track_col)], errors="coerce")
    full_obs["__time"] = pd.to_numeric(full_obs[str(time_col)], errors="coerce")
    full_obs = full_obs.dropna(subset=["__sample", "__track", "__time"]).copy()
    full_obs = full_obs[full_obs["__sample"] != ""].copy()
    full_obs["__track"] = full_obs["__track"].astype(np.int64)
    full_obs["__time"] = full_obs["__time"].astype(np.int64)
    full_obs = full_obs.sort_values(["__sample", "__track", "__time"], kind="mergesort")

    times_by_track = {
        (str(sample), int(track)): np.asarray(g["__time"].to_numpy(dtype=np.int64), dtype=np.int64)
        for (sample, track), g in full_obs.groupby(["__sample", "__track"], observed=True, sort=False)
    }

    has_tmin = "position_t_min" in tracks_obs.columns
    has_tmax = "position_t_max" in tracks_obs.columns
    has_n_timepoints = "n_timepoints" in tracks_obs.columns
    expanded_parts = []
    missing_windows = 0

    for _, row in tracks_obs.iterrows():
        s = str(row["__sample"])
        tr = int(row["__track"])
        cluster_label = str(row["__cluster"])
        times = times_by_track.get((s, tr), None)
        if times is None or times.size == 0:
            missing_windows += 1
            continue

        chosen = times
        if has_tmin and has_tmax:
            tmin = pd.to_numeric(pd.Series([row.get("position_t_min", np.nan)]), errors="coerce").iloc[0]
            tmax = pd.to_numeric(pd.Series([row.get("position_t_max", np.nan)]), errors="coerce").iloc[0]
            if pd.notna(tmin) and pd.notna(tmax):
                lo = int(min(tmin, tmax))
                hi = int(max(tmin, tmax))
                chosen = chosen[(chosen >= lo) & (chosen <= hi)]

        if has_n_timepoints:
            n_pos = pd.to_numeric(pd.Series([row.get("n_timepoints", np.nan)]), errors="coerce").iloc[0]
            if pd.notna(n_pos):
                n_pos = max(0, int(n_pos))
                if n_pos > 0 and chosen.size > n_pos:
                    chosen = chosen[-n_pos:]

        if chosen.size == 0:
            missing_windows += 1
            continue

        expanded_parts.append(
            pd.DataFrame(
                {
                    str(sample_col): np.repeat(s, chosen.size),
                    str(track_col): np.repeat(tr, chosen.size),
                    str(time_col): chosen.astype(np.int64, copy=False),
                    str(output_col): np.repeat(cluster_label, chosen.size),
                }
            )
        )

    if len(expanded_parts) == 0:
        raise ValueError(
            "Could not export track-cluster backprojection: no matching (sample, TrackID, position_t) "
            "rows were resolved from adata_tracks windows."
        )

    backproj_obs = pd.concat(expanded_parts, ignore_index=True)
    backproj_obs = backproj_obs.sort_values(
        [str(sample_col), str(track_col), str(time_col)],
        kind="mergesort",
    ).drop_duplicates(
        subset=[str(sample_col), str(track_col), str(time_col)],
        keep="last",
    )
    _vinfo(
        verbose,
        "trajectory-apply",
        "expanded track windows to frame assignments | "
        f"tracks={int(len(tracks_obs))} | rows={int(len(backproj_obs))} | "
        f"missing_tracks={int(missing_windows)}",
    )

    adata_export = SimpleNamespace(
        obs=backproj_obs,
        uns={},
        filename=getattr(adata_full, "filename", None),
    )

    manifest = export_behavioral_state_backprojection_zarrs(
        adata=adata_export,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=output_col,
        sample_col=sample_col,
        track_col=track_col,
        time_col=time_col,
        enforce_time_coverage=False,
        verbose=verbose,
    )
    manifest["window_backprojection_rows"] = int(len(backproj_obs))
    manifest["window_backprojection_missing_tracks"] = int(missing_windows)
    return manifest


def run_state_based_analysis(
    output_dir,
    cell_type="tcell",
    
    # Input
    adata_full_path=None,  # if None, will look under output_dir/analysis/<cell_type>/behavioral_states/BEHAV3D_<cell_type>_behavioral_states.h5ad
    state_col="full_behavioral_cluster",
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",

    # Track filtering / truncation
    behavioral_trajectory_size=100,

    # Feature extraction / normalization
    use_fractions=True,
    use_bout_stats=True,
    include_bouts_mean_length=True,
    include_bouts_nr=True,
    include_bouts_max_length=True,
    use_transitions=True,
    use_bigrams=True,
    use_trigrams=True,
    ngram_weight="count",
    do_block_scaling=False,
    do_l2_normalization=False,

    # Feature selection
    drop_highly_correlated=False,
    corr_threshold=0.95,
    drop_low_variance=False,
    low_var_threshold=1e-4,

    # PCA / Leiden / UMAP
    do_pca=True,
    pca_var_selection=0.95,
    n_neighbors=30,
    leiden_resolution=0.2,
    leiden_metric="euclidean",  # or cosine
    leiden_use_rep="X",   # "X" matches your snippet; could also be "X_pca"
    cluster_key="ClusterID",
    umap_min_dist=0.1,

    # Plotting
    plot_results=True,
    autosave_plots=True,
    plot_file_ext="pdf",
    plot_dpi=300,
    heatmap_figsize=(25, 20),
    matrixplot_figsize=(20, 40),
    umap_size=1,
    umap_alpha=0.5,

    # Exemplar track plotting
    plot_exemplars=True,
    n_per_cluster=10,
    exemplar_state_keys=("full_behavioral_cluster",),
    plot_exemplar_backprojection_videos=True,
    plot_exemplar_backprojection_pdfs=False,
    exemplar_video_fps=6,
    exemplar_video_dpi=180,
    exemplar_video_margin=20,
    exemplar_video_pmin=0.0,
    exemplar_video_pmax=99.99,
    exemplar_video_track_color="#63ff33",
    exemplar_backprojection_show_segments=True,
    exemplar_backprojection_segment_style="outline",
    exemplar_backprojection_segment_color="#ffffff",
    exemplar_backprojection_layout_mode="both",
    exemplar_backprojection_examples_per_cluster=5,
    exemplar_backprojection_num_example_ranks=5,
    exemplar_backprojection_pdf_rows_per_page=6,
    exemplar_backprojection_pdf_show_raw_image=True,
    exemplar_backprojection_show_state_legend=True,

    # Saving
    save_outputs=True,
    output_subdir_name="behavioral_state_trajectories",
    relabel_mapping=None,
    relabel_keep_unmapped=True,

    random_state=123,
    verbose=True,
):
    run_started = _vstart(verbose, "trajectory-clustering", "run trajectory clustering")
    start_time = time.time()

    resolved_paths = _resolve_track_paths(
        output_dir=output_dir,
        cell_type=cell_type,
        output_subdir_name=output_subdir_name,
    )
    outfolder = resolved_paths["outfolder"]
    clustering_outfolder = outfolder / "clustering"
    exemplar_root = outfolder / "clustering" / "example_tracks"
    state_folder = resolved_paths["state_outdir"]
    
    if adata_full_path is None:
        adata_full_path = Path(state_folder, f"BEHAV3D_{cell_type}_behavioral_states.h5ad")

    adata_full_path = Path(adata_full_path)
    if not adata_full_path.exists():
        raise FileNotFoundError(f"Could not find BEHAV3D_{cell_type}_behavioral_states.h5ad at: {adata_full_path}")
    load_started = _vstart(verbose, "trajectory-clustering", f"load behavioral-state adata | path={adata_full_path}")
    adata_full = sc.read_h5ad(adata_full_path)
    _vinfo(verbose, "trajectory-clustering", f"loaded dataset | rows={adata_full.n_obs} | features={adata_full.n_vars}")
    _vdone(verbose, "trajectory-clustering", "load dataset", load_started)

    # --------- Filter + truncate tracks ----------
    filter_started = _vstart(
        verbose,
        "trajectory-clustering",
        f"filter + truncate trajectories | behavioral_trajectory_size={int(behavioral_trajectory_size)}",
    )
    adata_filt = filter_and_truncate_tracks_anndata(
        adata_full,
        groupby_cols=list(groupby_cols),
        time_col=time_col,
        min_length=int(behavioral_trajectory_size),
        max_length=int(behavioral_trajectory_size),
    )
    _vdone(verbose, "trajectory-clustering", "filter + truncate trajectories", filter_started)

    # --------- Extract state-describing features ----------
    feature_started = _vstart(
        verbose,
        "trajectory-clustering",
        f"extract trajectory features | state_col={state_col}",
    )
    adata_state_features, blocks = extract_descibing_track_state_features(
        adata_filt,
        state_col=state_col,
        use_fractions=use_fractions,
        use_bout_stats=use_bout_stats,
        use_transitions=use_transitions,
        use_bigrams=use_bigrams,
        use_trigrams=use_trigrams,
        ngram_weight=ngram_weight,
    )
    _vdone(verbose, "trajectory-clustering", "extract trajectory features", feature_started)

    if bool(use_bout_stats):
        drop_cols = []
        if not bool(include_bouts_mean_length):
            drop_cols.extend(
                [c for c in adata_state_features.var_names if str(c).startswith("bouts_mean_length_")]
            )
        if not bool(include_bouts_nr):
            drop_cols.extend(
                [c for c in adata_state_features.var_names if str(c).startswith("bouts_nr_")]
            )
        if not bool(include_bouts_max_length):
            drop_cols.extend(
                [c for c in adata_state_features.var_names if str(c).startswith("bouts_max_length_")]
            )
        if len(drop_cols) > 0:
            keep_cols = [c for c in adata_state_features.var_names if c not in set(drop_cols)]
            adata_state_features = adata_state_features[:, keep_cols].copy()
            feature_blocks = dict(adata_state_features.uns.get("feature_blocks", {}))
            if "bout_stats" in feature_blocks:
                feature_blocks["bout_stats"] = [c for c in feature_blocks["bout_stats"] if c in keep_cols]
            adata_state_features.uns["feature_blocks"] = feature_blocks

    feature_blocks = dict(adata_state_features.uns.get("feature_blocks", {}))
    blocks = [
        [c for c in feature_blocks.get("fractions", []) if c in adata_state_features.var_names],
        [c for c in feature_blocks.get("bout_stats", []) if c in adata_state_features.var_names],
        [c for c in feature_blocks.get("transitions", []) if c in adata_state_features.var_names],
        [c for c in feature_blocks.get("bigrams", []) if c in adata_state_features.var_names],
        [c for c in feature_blocks.get("trigrams", []) if c in adata_state_features.var_names],
    ]

    adata_state_features.uns["track_feature_selection"] = {
        "state_col": str(state_col),
        "use_fractions": bool(use_fractions),
        "use_bout_stats": bool(use_bout_stats),
        "include_bouts_mean_length": bool(include_bouts_mean_length),
        "include_bouts_nr": bool(include_bouts_nr),
        "include_bouts_max_length": bool(include_bouts_max_length),
        "use_transitions": bool(use_transitions),
        "use_bigrams": bool(use_bigrams),
        "use_trigrams": bool(use_trigrams),
        "ngram_weight": str(ngram_weight),
    }
    adata_state_features.uns["track_filtering"] = {
        "groupby_cols": [str(c) for c in list(groupby_cols)],
        "time_col": str(time_col),
        "state_col": str(state_col),
        "behavioral_trajectory_size": int(behavioral_trajectory_size),
    }

    # --------- (Optional) Block scaling + L2 normalize ----------
    if do_block_scaling:
        adata_state_features = scale_feature_blocks(adata_state_features, blocks=blocks)

    if do_l2_normalization:
        adata_state_features = l2_normalize_features_blocks(adata_state_features, blocks=blocks)

    # --------- Feature selection ----------
    dropped_high_corr = []
    dropped_low_var = []

    if drop_highly_correlated:
        adata_state_features, dropped_high_corr, report = drop_highly_correlated_features(
            data=adata_state_features,
            feature_cols=adata_state_features.var_names,
            threshold=corr_threshold,
        )

    kept_features = list(adata_state_features.var_names)

    if drop_low_variance:
        adata_state_features, kept_features, dropped_low_var = drop_low_variance_features(
            data=adata_state_features,
            feature_cols=kept_features,
            low_var_threshold=low_var_threshold,
        )

    resolved_leiden_use_rep = str(leiden_use_rep)
    # --------- PCA ----------
    if do_pca:
        pca_started = _vstart(
            verbose,
            "trajectory-clustering",
            f"PCA | pca_var_selection={float(pca_var_selection)}",
        )
        adata_state_features = run_pca(
            adata_state_features,
            ncomps=len(adata_state_features.var_names),
            pca_var_selection=pca_var_selection,  # matches your snippet signature
            random_state=random_state,
        )
        if resolved_leiden_use_rep == "X":
            resolved_leiden_use_rep = "X_pca"
        _vdone(verbose, "trajectory-clustering", "PCA", pca_started)

    # --------- Leiden ----------
    # Keep embedding source configurable (X vs X_pca).
    leiden_started = _vstart(
        verbose,
        "trajectory-clustering",
        "Leiden clustering | "
        f"neighbors={int(n_neighbors)} | resolution={float(leiden_resolution)} | metric={leiden_metric}",
    )
    adata_state_features = run_leiden_clustering(
            adata_state_features, 
            n_neighbors=n_neighbors,
            resolution=leiden_resolution,
            method="umap",
            use_rep=resolved_leiden_use_rep,
            key_added=cluster_key,
            random_state=random_state,
            metric=leiden_metric
        )
    _vdone(verbose, "trajectory-clustering", "Leiden clustering", leiden_started)

    if isinstance(relabel_mapping, dict) and len(relabel_mapping) > 0:
        _vinfo(verbose, "trajectory-clustering", f"apply relabel mapping ({len(relabel_mapping)} entries)")
        rename_track_clusters(
            adata=adata_state_features,
            mapping=relabel_mapping,
            cluster_col=cluster_key,
            keep_unmapped=bool(relabel_keep_unmapped),
        )

    # --------- UMAP ----------
    umap_started = _vstart(
        verbose,
        "trajectory-clustering",
        f"UMAP embedding | min_dist={float(umap_min_dist)}",
    )
    sc.tl.umap(
        adata_state_features,
        min_dist=umap_min_dist,
        random_state=random_state,
    )
    _vdone(verbose, "trajectory-clustering", "UMAP embedding", umap_started)

    if plot_results:
        if bool(autosave_plots):
            diagnostics_started = time.perf_counter()
            clustering_outfolder.mkdir(parents=True, exist_ok=True)
            report_paths = generate_track_clustering_report_pdfs(
                adata_tracks=adata_state_features,
                outfolder=clustering_outfolder,
                cluster_key=cluster_key,
                heatmap_figsize=heatmap_figsize,
                matrixplot_figsize=matrixplot_figsize,
                umap_size=umap_size,
                umap_alpha=umap_alpha,
                plot_dpi=plot_dpi,
                verbose=verbose,
            )
            _vsave(verbose, "trajectory-clustering", "clustering diagnostics PDF", report_paths["diagnostics_pdf"])
            _vdone(verbose, "trajectory-clustering", "write clustering diagnostics", diagnostics_started)

    exemplar_statebar_track_pdf_by_cluster = {}
    exemplar_statebar_track_pdf_by_example_rank = {}
    exemplar_backprojection_pdf_paths_by_cluster = {}
    exemplar_backprojection_pdf_paths_by_example_rank = {}
    exemplar_backprojection_video_paths_by_cluster = {}
    exemplar_backprojection_video_paths_by_example_rank = {}
    exemplar_backprojection_global_canvas_xy_shape = None
    exemplar_selection_csv = None
    exemplar_render_config = {
        "stage": "after_clustering",
        "cluster_key": str(cluster_key),
        "state_key": str(state_col),
        "n_per_cluster": int(n_per_cluster),
        "autosave_plots": bool(autosave_plots),
        "plot_exemplars": bool(plot_exemplars),
        "plot_exemplar_backprojection_videos": bool(plot_exemplar_backprojection_videos),
        "plot_exemplar_backprojection_pdfs": bool(plot_exemplar_backprojection_pdfs),
        "plot_dpi": int(plot_dpi),
        "video_fps": int(exemplar_video_fps),
        "video_dpi": int(exemplar_video_dpi),
        "video_margin": int(exemplar_video_margin),
        "video_pmin": float(exemplar_video_pmin),
        "video_pmax": float(exemplar_video_pmax),
        "video_track_color": str(exemplar_video_track_color),
        "backprojection_show_segments": bool(exemplar_backprojection_show_segments),
        "backprojection_segment_style": str(exemplar_backprojection_segment_style),
        "backprojection_segment_color": str(exemplar_backprojection_segment_color),
        "backprojection_layout_mode": str(exemplar_backprojection_layout_mode),
        "backprojection_examples_per_cluster": int(exemplar_backprojection_examples_per_cluster),
        "backprojection_num_example_ranks": int(exemplar_backprojection_num_example_ranks),
        "backprojection_pdf_rows_per_page": int(exemplar_backprojection_pdf_rows_per_page),
        "backprojection_pdf_show_raw_image": bool(exemplar_backprojection_pdf_show_raw_image),
        "backprojection_show_state_legend": bool(exemplar_backprojection_show_state_legend),
        "coordinate_enrichment": None,
    }

    # --------- Exemplar tracks by cluster ----------
    if plot_exemplars:
        exemplar_started = time.perf_counter()
        # This assumes plot_exemplar_tracks_by_cluster signature: (adata_tracks, adata_clusters, n_per_cluster, state_key)
        # where `adata_clusters` contains the clustering in .obs[cluster_key].
        fig_exemplar, _, chosen_exemplars = plot_exemplar_tracks_by_cluster(
            adata_filt,
            adata_state_features,
            n_per_cluster=n_per_cluster,
            state_key=state_col,
            cluster_key=cluster_key,
        )
        if bool(autosave_plots) and fig_exemplar is not None:
            exemplar_root.mkdir(parents=True, exist_ok=True)
            exemplar_path = Path(
                exemplar_root,
                "example_tracks_overview.pdf",
            )
            exemplar_path.parent.mkdir(parents=True, exist_ok=True)
            with PdfPages(exemplar_path) as pdf:
                _apply_best_pdf_orientation(fig_exemplar, default_orientation="landscape")
                pdf.savefig(fig_exemplar, dpi=int(plot_dpi), bbox_inches="tight")

            coord_enrichment_pdf = _ensure_exemplar_coordinate_columns(
                adata=adata_filt,
                output_dir=output_dir,
                cell_type=cell_type,
                require_pixel_for_video=False,
            )
            exemplar_render_config["coordinate_enrichment"] = dict(coord_enrichment_pdf)

            exemplar_selection_csv = exemplar_root / (
                f"example_track_selection_cluster_{_sanitize_filename_token(cluster_key, fallback='cluster')}_"
                f"state_{_sanitize_filename_token(state_col, fallback='state')}_clustering.csv"
            )
            chosen_exemplars.to_csv(exemplar_selection_csv, index=False)

            per_cluster_pdf_out = save_exemplar_statebar_track_pdf_per_cluster(
                adata_full=adata_filt,
                out_dir=exemplar_root,
                chosen_df=chosen_exemplars,
                adata_tracks=None,
                n_per_cluster=int(n_per_cluster),
                sample_key="sample_name",
                track_key="TrackID",
                time_key=str(time_col),
                state_key=str(state_col),
                cluster_key=str(cluster_key),
                tmin_key="position_t_min",
                tmax_key="position_t_max",
                rows_per_page=6,
                plot_dpi=int(plot_dpi),
                seed=int(random_state),
                cmap_name="tab20",
                layout_mode=str(exemplar_backprojection_layout_mode),
                num_example_ranks=int(exemplar_backprojection_num_example_ranks),
            )
            exemplar_statebar_track_pdf_by_cluster = dict(
                per_cluster_pdf_out.get("pdf_paths_by_cluster", {})
            )
            exemplar_statebar_track_pdf_by_example_rank = dict(
                per_cluster_pdf_out.get("pdf_paths_by_example_rank", {})
            )

            if bool(plot_exemplar_backprojection_videos) or bool(plot_exemplar_backprojection_pdfs):
                coord_enrichment_video = _ensure_exemplar_coordinate_columns(
                    adata=adata_filt,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    require_pixel_for_video=True,
                )
                exemplar_render_config["coordinate_enrichment_video"] = dict(coord_enrichment_video)
                exemplar_backprojection_outdir = exemplar_root
                exemplar_backprojection_outdir.mkdir(parents=True, exist_ok=True)
                prep_t0 = time.perf_counter()
                exemplar_backprojection_prepared = _prepare_exemplar_backprojection_data(
                    adata_full=adata_filt,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    chosen_df=chosen_exemplars.copy().reset_index(drop=True),
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key=str(time_col),
                    state_key=str(state_col),
                    cluster_key=str(cluster_key),
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    margin=int(exemplar_video_margin),
                    pmin=float(exemplar_video_pmin),
                    pmax=float(exemplar_video_pmax),
                    cmap_name="tab20",
                    coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                    examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                    show_segment_outlines=bool(exemplar_backprojection_show_segments),
                    segment_style=str(exemplar_backprojection_segment_style),
                    segment_color=str(exemplar_backprojection_segment_color),
                )
                exemplar_render_config["segment_outline_errors"] = dict(
                    exemplar_backprojection_prepared.get("segment_outline_errors", {})
                )
                _vdone(verbose, "trajectory-clustering", "prepare backprojection data", prep_t0)

                if bool(plot_exemplar_backprojection_videos):
                    backprojection_video_out = save_exemplar_statebar_backprojection_video_per_cluster(
                        adata_full=adata_filt,
                        output_dir=output_dir,
                        cell_type=cell_type,
                        out_dir=exemplar_backprojection_outdir,
                        chosen_df=chosen_exemplars,
                        adata_tracks=None,
                        n_per_cluster=int(n_per_cluster),
                        sample_key="sample_name",
                        track_key="TrackID",
                        time_key=str(time_col),
                        state_key=str(state_col),
                        cluster_key=str(cluster_key),
                        tmin_key="position_t_min",
                        tmax_key="position_t_max",
                        fps=int(exemplar_video_fps),
                        dpi=int(exemplar_video_dpi),
                        margin=int(exemplar_video_margin),
                        pmin=float(exemplar_video_pmin),
                        pmax=float(exemplar_video_pmax),
                        track_color=str(exemplar_video_track_color),
                        show_segment_outlines=bool(exemplar_backprojection_show_segments),
                        segment_style=str(exemplar_backprojection_segment_style),
                        segment_color=str(exemplar_backprojection_segment_color),
                        coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                        seed=int(random_state),
                        cmap_name="tab20",
                        layout_mode=str(exemplar_backprojection_layout_mode),
                        examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                        num_example_ranks=int(exemplar_backprojection_num_example_ranks),
                        show_state_legend=bool(exemplar_backprojection_show_state_legend),
                        ffmpeg_preset="veryfast",
                        ffmpeg_crf=21,
                        prepared_data=exemplar_backprojection_prepared,
                        verbose=bool(verbose),
                    )
                    exemplar_backprojection_video_paths_by_cluster = dict(
                        backprojection_video_out.get("video_paths_by_cluster", {})
                    )
                    exemplar_backprojection_video_paths_by_example_rank = dict(
                        backprojection_video_out.get("video_paths_by_example_rank", {})
                    )
                    exemplar_backprojection_global_canvas_xy_shape = backprojection_video_out.get(
                        "global_xy_shape", exemplar_backprojection_global_canvas_xy_shape
                    )

                if bool(plot_exemplar_backprojection_pdfs):
                    backprojection_pdf_out = save_exemplar_statebar_backprojection_pdf(
                        adata_full=adata_filt,
                        output_dir=output_dir,
                        cell_type=cell_type,
                        out_dir=exemplar_backprojection_outdir,
                        chosen_df=chosen_exemplars,
                        adata_tracks=None,
                        n_per_cluster=int(n_per_cluster),
                        sample_key="sample_name",
                        track_key="TrackID",
                        time_key=str(time_col),
                        state_key=str(state_col),
                        cluster_key=str(cluster_key),
                        tmin_key="position_t_min",
                        tmax_key="position_t_max",
                        rows_per_page=int(exemplar_backprojection_pdf_rows_per_page),
                        plot_dpi=max(120, int(exemplar_video_dpi)),
                        margin=int(exemplar_video_margin),
                        pmin=float(exemplar_video_pmin),
                        pmax=float(exemplar_video_pmax),
                        track_color=str(exemplar_video_track_color),
                        show_segment_outlines=bool(exemplar_backprojection_show_segments),
                        segment_style=str(exemplar_backprojection_segment_style),
                        segment_color=str(exemplar_backprojection_segment_color),
                        coordinate_source_hint=coord_enrichment_video.get("csv_path"),
                        seed=int(random_state),
                        cmap_name="tab20",
                        layout_mode=str(exemplar_backprojection_layout_mode),
                        examples_per_cluster=int(exemplar_backprojection_examples_per_cluster),
                        num_example_ranks=int(exemplar_backprojection_num_example_ranks),
                        show_raw_image=bool(exemplar_backprojection_pdf_show_raw_image),
                        show_state_legend=bool(exemplar_backprojection_show_state_legend),
                        prepared_data=exemplar_backprojection_prepared,
                        verbose=bool(verbose),
                    )
                    exemplar_backprojection_pdf_paths_by_cluster = dict(
                        backprojection_pdf_out.get("pdf_paths_by_cluster", {})
                    )
                    exemplar_backprojection_pdf_paths_by_example_rank = dict(
                        backprojection_pdf_out.get("pdf_paths_by_example_rank", {})
                    )
                    exemplar_backprojection_global_canvas_xy_shape = backprojection_pdf_out.get(
                        "global_xy_shape", exemplar_backprojection_global_canvas_xy_shape
                    )
                    exemplar_render_config["segment_outline_errors"] = dict(
                        backprojection_pdf_out.get(
                            "segment_outline_errors",
                            exemplar_render_config.get("segment_outline_errors", {}),
                        )
                    )
        _vdone(verbose, "trajectory-clustering", "render exemplar outputs", exemplar_started)
        plt.close(fig_exemplar)

    adata_state_features.uns.setdefault("visualization", {})
    adata_state_features.uns["visualization"].update(
        {
            "exemplar_render_stage": "after_clustering",
            "exemplar_statebar_track_pdf_by_cluster": dict(exemplar_statebar_track_pdf_by_cluster),
            "exemplar_statebar_track_pdf_by_example_rank": dict(
                exemplar_statebar_track_pdf_by_example_rank
            ),
            "exemplar_statebar_backprojection_video_by_cluster": dict(
                exemplar_backprojection_video_paths_by_cluster
            ),
            "exemplar_backprojection_pdf_paths_by_cluster": dict(
                exemplar_backprojection_pdf_paths_by_cluster
            ),
            "exemplar_backprojection_pdf_paths_by_example_rank": dict(
                exemplar_backprojection_pdf_paths_by_example_rank
            ),
            "exemplar_backprojection_video_paths_by_cluster": dict(
                exemplar_backprojection_video_paths_by_cluster
            ),
            "exemplar_backprojection_video_paths_by_example_rank": dict(
                exemplar_backprojection_video_paths_by_example_rank
            ),
            "exemplar_backprojection_global_canvas_xy_shape": exemplar_backprojection_global_canvas_xy_shape,
            "exemplar_selection_csv": (
                None if exemplar_selection_csv is None else str(exemplar_selection_csv)
            ),
            "exemplar_segment_outline_errors": dict(
                exemplar_render_config.get("segment_outline_errors", {})
            ),
            "exemplar_render_config": dict(exemplar_render_config),
        }
    )


    # --------- Save outputs ----------
    adata_feat_out = Path(outfolder, get_track_trajectories_filename(cell_type))
    if save_outputs:
        save_started = time.perf_counter()
        adata_state_features.write(adata_feat_out, compression="gzip")
        _vsave(verbose, "trajectory-clustering", "clustered trajectory model", adata_feat_out)
        _vdone(verbose, "trajectory-clustering", "save clustered model", save_started)

    elapsed = time.time() - start_time
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    _vinfo(verbose, "trajectory-clustering", f"elapsed={h}:{m:02}:{s:02}")
    _vdone(verbose, "trajectory-clustering", "run trajectory clustering", run_started)

    return adata_state_features


def test():
    """Manual helper for quick local test runs of track state-based analysis."""
    output_dir = Path(r"F:/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    cell_type = "tcell"
    adata_full_path = None
    state_col = "full_behavioral_cluster"
    groupby_cols = ("sample_name", "TrackID")
    time_col = "position_t"
    behavioral_trajectory_size = 100
    use_fractions = True
    use_bout_stats = True
    include_bouts_mean_length = True
    include_bouts_nr = True
    include_bouts_max_length = True
    use_transitions = True
    use_bigrams = True
    use_trigrams = True
    ngram_weight = "count"
    do_block_scaling = False
    do_l2_normalization = False
    drop_highly_correlated = False
    corr_threshold = 0.95
    drop_low_variance = False
    low_var_threshold = 1e-4
    do_pca = True
    pca_var_selection = 0.95
    n_neighbors = 30
    leiden_resolution = 0.2
    leiden_metric = "euclidean"
    leiden_use_rep = "X"
    umap_min_dist = 0.1
    plot_results = True
    autosave_plots = True
    plot_file_ext = "pdf"
    plot_dpi = 300
    heatmap_figsize = (25, 20)
    matrixplot_figsize = (20, 40)
    umap_size = 1
    umap_alpha = 0.5
    plot_exemplars = True
    n_per_cluster = 10
    exemplar_state_keys = ("full_behavioral_cluster",)
    save_outputs = True
    output_subdir_name = "behavioral_state_trajectories"
    random_state = 123

    run_state_based_analysis(
        output_dir=output_dir,
        cell_type=cell_type,
        adata_full_path=adata_full_path,
        state_col=state_col,
        groupby_cols=groupby_cols,
        time_col=time_col,
        behavioral_trajectory_size=behavioral_trajectory_size,
        use_fractions=use_fractions,
        use_bout_stats=use_bout_stats,
        include_bouts_mean_length=include_bouts_mean_length,
        include_bouts_nr=include_bouts_nr,
        include_bouts_max_length=include_bouts_max_length,
        use_transitions=use_transitions,
        use_bigrams=use_bigrams,
        use_trigrams=use_trigrams,
        ngram_weight=ngram_weight,
        do_block_scaling=do_block_scaling,
        do_l2_normalization=do_l2_normalization,
        drop_highly_correlated=drop_highly_correlated,
        corr_threshold=corr_threshold,
        drop_low_variance=drop_low_variance,
        low_var_threshold=low_var_threshold,
        do_pca=do_pca,
        pca_var_selection=pca_var_selection,
        n_neighbors=n_neighbors,
        leiden_resolution=leiden_resolution,
        leiden_metric=leiden_metric,
        leiden_use_rep=leiden_use_rep,
        umap_min_dist=umap_min_dist,
        plot_results=plot_results,
        autosave_plots=autosave_plots,
        plot_file_ext=plot_file_ext,
        plot_dpi=plot_dpi,
        heatmap_figsize=heatmap_figsize,
        matrixplot_figsize=matrixplot_figsize,
        umap_size=umap_size,
        umap_alpha=umap_alpha,
        plot_exemplars=plot_exemplars,
        n_per_cluster=n_per_cluster,
        exemplar_state_keys=exemplar_state_keys,
        save_outputs=save_outputs,
        output_subdir_name=output_subdir_name,
        random_state=random_state,
    )
