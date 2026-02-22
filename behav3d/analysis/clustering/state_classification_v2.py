import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype

from pathlib import Path
import pickle

import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler

from behav3d.core.anndata import df_to_adata
from behav3d.analysis.clustering.state_classification import (
    _prepare_state_classification_dataset,
    cluster_group,
    subsample_with_temporal_spacing,
)
from behav3d.analysis.clustering.general import relabel_cluster_ids


def _binary_col_to_group_name(col: str) -> str:
    name = str(col)
    name = name.replace("_pixels", "")
    name = name.replace("_value", "")
    return name


def _assign_binary_group_labels(df: pd.DataFrame, binary_cols: list[str]) -> pd.Series:
    """Build a single categorical group label from binary indicator columns."""
    if len(binary_cols) == 0:
        return pd.Series(["no_contact"] * len(df), index=df.index, dtype="object")

    labels = []
    for _, row in df[binary_cols].iterrows():
        active = []
        for col in binary_cols:
            val = row[col]
            if pd.notna(val) and float(val) == 1.0:
                active.append(_binary_col_to_group_name(col))
        if len(active) == 0:
            labels.append("no_contact")
        elif len(active) == 1:
            labels.append(active[0])
        else:
            labels.append("_and_".join(active))
    return pd.Series(labels, index=df.index, dtype="object")


def _add_clean_binary_annotation_columns(df: pd.DataFrame, binary_cols: list[str]) -> pd.DataFrame:
    """Add user-facing binary annotation columns (e.g. organoid_contact) to obs."""
    out = df.copy()
    for col in binary_cols:
        if col not in out.columns:
            continue
        clean_name = _binary_col_to_group_name(col)
        if clean_name == col or clean_name in out.columns:
            continue
        vals = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
        out[clean_name] = (vals > 0).astype(np.int8)
    return out


def check_model_clusterid_consistency(
    model_adata,
    adata_full,
    model_cluster_col="ClusterID",
    full_cluster_col="ClusterID",
    key_cols=("sample_name", "TrackID", "position_t"),
    max_examples=10,
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

    print(
        "ClusterID consistency check "
        f"(mode={summary['match_mode']}): "
        f"checked={summary['n_model_rows_checked']}, "
        f"matches={summary['n_matching_clusterid']}, "
        f"mismatches={summary['n_mismatched_clusterid']}, "
        f"missing_in_full={summary['n_missing_rows_in_adata_full']}"
    )
    if n_mismatch > 0:
        cols = [c for c in merged.columns if c not in {"_match"}]
        print("First mismatches:")
        print(merged.loc[(~merged["_match"]) & merged["_cluster_full"].notna(), cols].head(max_examples))

    return summary, merged


def _to_numpy_2d(X):
    """Convert dense/sparse matrix-like data to a 2D numpy array."""
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def train_clusterid_classifier_from_model_adata(
    model_adata,
    cluster_col="ClusterID",
    classifier_backend="random_forest",
    classifier_kwargs=None,
    random_state=123,
):
    """Train a supervised classifier on model_adata features -> cluster labels."""
    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = list(model_adata.var_names)
    X_train = _to_numpy_2d(model_adata[:, feature_cols].X)
    y_train = model_adata.obs[cluster_col].astype(str).to_numpy()
    if len(y_train) == 0:
        raise ValueError("No rows in model_adata to train classifier.")

    backend = str(classifier_backend).lower()
    params = dict(classifier_kwargs) if classifier_kwargs is not None else {}

    if backend in {"random_forest", "rf", "randomforest", "randomforestclassifier"}:
        default_params = {
            "n_estimators": 300,
            "min_samples_leaf": 2,
            "n_jobs": -1,
            "random_state": int(random_state),
        }
        default_params.update(params)
        clf = RandomForestClassifier(**default_params)
        backend_name = "random_forest"
    elif backend in {"logistic_regression", "logreg", "lr"}:
        default_params = {
            "max_iter": 1000,
            "multi_class": "auto",
            "random_state": int(random_state),
        }
        default_params.update(params)
        clf = LogisticRegression(**default_params)
        backend_name = "logistic_regression"
    else:
        raise ValueError(
            "Unknown classifier_backend. Use one of: "
            "'random_forest', 'logistic_regression'."
        )

    clf.fit(X_train, y_train)
    train_acc = float((clf.predict(X_train) == y_train).mean())
    fit_info = {
        "backend": backend_name,
        "cluster_col": cluster_col,
        "n_train_rows": int(X_train.shape[0]),
        "n_features": int(X_train.shape[1]),
        "classes": [str(c) for c in getattr(clf, "classes_", [])],
        "train_accuracy": train_acc,
        "params": {k: repr(v) for k, v in clf.get_params().items()},
    }
    print(
        "Trained cluster classifier "
        f"(backend={backend_name}, rows={fit_info['n_train_rows']}, "
        f"features={fit_info['n_features']}, train_acc={train_acc:.4f})"
    )
    return clf, fit_info


def predict_clusterids_with_classifier(
    adata,
    classifier,
    feature_cols=None,
    output_col="ClusterID",
    confidence_col=None,
    inplace=True,
):
    """Predict cluster labels for an AnnData object using a trained classifier."""
    target = adata if inplace else adata.copy()
    cols = list(target.var_names) if feature_cols is None else list(feature_cols)

    missing = [c for c in cols if c not in target.var_names]
    if len(missing) > 0:
        raise ValueError(f"adata is missing classifier feature columns: {missing[:10]}")

    X_query = _to_numpy_2d(target[:, cols].X)
    pred = classifier.predict(X_query)
    target.obs[output_col] = pd.Categorical(pd.Series(pred, index=target.obs.index).astype(str))

    if confidence_col is not None and hasattr(classifier, "predict_proba"):
        proba = classifier.predict_proba(X_query)
        target.obs[confidence_col] = pd.Series(proba.max(axis=1), index=target.obs.index, dtype=float)

    return target


def _mixed_label_sort_key(x):
    s = str(x)
    return (0, int(s)) if s.isdigit() else (1, s)


def evaluate_classifier_predictions(
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="classifier_qc",
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
                    axes[1].text(j, i, f"{cm_true_norm[i, j]:.2f}", ha="center", va="center", fontsize=7, color="black")

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
    cluster_col="ClusterID",
    outfolder=None,
    filename_prefix="state_classification_v2_classifier_qc",
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
    )

    print(
        "Classifier QC on model data: "
        f"n={summary['n_rows']}, classes={summary['n_classes']}, "
        f"accuracy={summary['accuracy']:.4f}, balanced_accuracy={summary['balanced_accuracy']:.4f}"
    )

    return summary


def _resolve_binary_classifier_feature_cols(adata, binary_cols_to_merge):
    """Resolve which binary obs columns should be used in full-label classifier."""
    cols = []
    for col in list(binary_cols_to_merge):
        clean_col = _binary_col_to_group_name(col)
        if clean_col in adata.obs.columns:
            cols.append(clean_col)
        elif col in adata.obs.columns:
            cols.append(col)
    if "binary_group" in adata.obs.columns:
        cols.append("binary_group")
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(cols))


def _build_binary_feature_dataframe(
    obs_df,
    binary_feature_cols,
    expected_expanded_cols=None,
):
    """Build numeric binary/group feature frame with one-hot encoding for categoricals."""
    pieces = []
    for col in list(binary_feature_cols):
        if col not in obs_df.columns:
            raise ValueError(f"Missing binary/group feature column in adata.obs: '{col}'")
        s = obs_df[col]
        if is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s):
            piece = pd.DataFrame(
                {col: pd.to_numeric(s, errors="coerce").fillna(0.0).astype(float)},
                index=obs_df.index,
            )
        else:
            piece = pd.get_dummies(
                s.astype("string").fillna("missing"),
                prefix=str(col),
                dtype=float,
            )
            piece.index = obs_df.index
        pieces.append(piece)

    if len(pieces) == 0:
        out = pd.DataFrame(index=obs_df.index)
    else:
        out = pd.concat(pieces, axis=1)

    if expected_expanded_cols is not None:
        expected = list(expected_expanded_cols)
        out = out.reindex(columns=expected, fill_value=0.0)
    return out


def _build_classifier_matrix_from_adata(
    adata,
    continuous_feature_cols,
    binary_feature_cols,
    binary_expanded_feature_cols=None,
    return_binary_feature_names=False,
):
    """Build a numeric feature matrix from adata.X (continuous) and obs (binary)."""
    cont_cols = list(continuous_feature_cols)
    bin_cols = list(binary_feature_cols)

    missing_cont = [c for c in cont_cols if c not in adata.var_names]
    if len(missing_cont) > 0:
        raise ValueError(f"Missing continuous feature columns in adata.var_names: {missing_cont[:10]}")

    X_cont = _to_numpy_2d(adata[:, cont_cols].X).astype(float, copy=False)
    if len(bin_cols) == 0:
        if return_binary_feature_names:
            return X_cont, []
        return X_cont

    bin_df = _build_binary_feature_dataframe(
        obs_df=adata.obs,
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
    target_col="ClusterID",
    binary_feature_cols=None,
    classifier_kwargs=None,
    random_state=123,
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

    params = {
        "n_estimators": 500,
        "min_samples_leaf": 2,
        "n_jobs": -1,
        "random_state": int(random_state),
    }
    if classifier_kwargs is not None:
        params.update(dict(classifier_kwargs))

    clf = RandomForestClassifier(**params)
    clf.fit(X_train, y_train)
    y_pred_train = clf.predict(X_train)
    train_acc = float((y_pred_train == y_train).mean())

    artifact = {
        "classifier": clf,
        "continuous_feature_cols": list(continuous_feature_cols),
        "binary_feature_cols": list(binary_cols),
        "binary_expanded_feature_cols": list(expanded_binary_cols),
        "target_col": str(target_col),
    }
    fit_info = {
        "backend": "random_forest",
        "target_col": str(target_col),
        "continuous_feature_cols": list(continuous_feature_cols),
        "binary_feature_cols": list(binary_cols),
        "n_train_rows": int(X_train.shape[0]),
        "n_features_total": int(X_train.shape[1]),
        "n_features_continuous": int(len(continuous_feature_cols)),
        "n_features_binary": int(len(binary_cols)),
        "n_features_binary_expanded": int(len(expanded_binary_cols)),
        "classes": [str(c) for c in getattr(clf, "classes_", [])],
        "train_accuracy": train_acc,
        "params": {k: repr(v) for k, v in clf.get_params().items()},
    }
    print(
        "Trained full ClusterID classifier "
        f"(rows={fit_info['n_train_rows']}, total_features={fit_info['n_features_total']}, "
        f"train_acc={train_acc:.4f})"
    )
    return artifact, fit_info


def predict_clusterids_with_full_classifier(
    adata,
    classifier_artifact,
    output_col="ClusterID_pred",
    confidence_col=None,
    inplace=True,
):
    """Apply a trained full ClusterID classifier artifact to an AnnData object."""
    target = adata if inplace else adata.copy()
    clf = classifier_artifact["classifier"]
    cont_cols = classifier_artifact["continuous_feature_cols"]
    bin_cols = classifier_artifact.get("binary_feature_cols", [])
    bin_expanded_cols = classifier_artifact.get("binary_expanded_feature_cols", None)

    X_query = _build_classifier_matrix_from_adata(
        adata=target,
        continuous_feature_cols=cont_cols,
        binary_feature_cols=bin_cols,
        binary_expanded_feature_cols=bin_expanded_cols,
    )
    pred = np.asarray(clf.predict(X_query)).astype(str)
    target.obs[output_col] = pd.Categorical(pd.Series(pred, index=target.obs.index).astype(str))

    if confidence_col is not None and hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(X_query)
        target.obs[confidence_col] = pd.Series(proba.max(axis=1), index=target.obs.index, dtype=float)

    return target


def plot_binary_group_behavioral_cluster_grid(
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

    # Diagnostic print to quickly verify cluster presence in plotting table.
    if "5" in ctab.columns:
        by_group = (ctab["5"] * 100.0).round(3)
        print(f"Cluster 5 total n={int(cluster_counts.get('5', 0))}; per-group proportions (%):")
        print(by_group[by_group > 0].sort_values(ascending=False))

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


def run_state_classification_v2(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    max_samples=None,
    min_spacing=None,
    n_neighbors=60,
    resolution="auto",
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    pca_var_selection=0.95,
    clustering_method="leiden",
    resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    leiden_seed_tries=10,
    leiden_subsample_tries=10,
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    incomplete_window_policy="drop",
    outfolder=None,
    random_state=123,
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
    label_transfer_method="classifier",
    classifier_backend="random_forest",
    classifier_kwargs=None,
    classifier_confidence_col="ClusterID_confidence",
    save_label_classifier=True,
    train_full_label_classifier=True,
    full_classifier_kwargs=None,
    full_classifier_confidence_col="ClusterID_full_classifier_confidence",
    save_full_label_classifier=True,
    return_artifacts=False,
):
    """
    V2 state classification:
    - cluster all continuous descriptive data jointly (single global model)
    - annotate binary contact info and combined ClusterID afterwards
    - transfer model labels to full dataset via `ingest` or a supervised classifier
    - optionally train final RF on continuous + binary features to predict final labels
    """
    prepared = _prepare_state_classification_dataset(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=window_size,
        descriptive_features=list(descriptive_features),
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        outfolder=outfolder,
        scale_features=False,
        incomplete_window_policy=incomplete_window_policy,
        reuse_prepared_dataset=reuse_prepared_dataset,
        save_prepared_dataset=save_prepared_dataset,
    )
    df_analysis = prepared["df_analysis"]
    kept_features = prepared["kept_features"]
    non_feature_cols = prepared["non_feature_cols"]
    binary_cols_to_merge = prepared["binary_cols_to_merge"]

    scaler = None
    df_scaled = df_analysis.copy()
    if len(kept_features) > 0:
        scaler = StandardScaler().fit(df_scaled[kept_features])
        df_scaled[kept_features] = scaler.transform(df_scaled[kept_features])


    if len(df_scaled) < 50:
        raise ValueError("Insufficient rows in full descriptive dataset for clustering.")

    if min_spacing is None:
        target_samples = max(500, len(df_scaled) // 20)
        spacing_to_use = max(1, len(df_scaled) // target_samples // 10)
    else:
        spacing_to_use = int(min_spacing)

    df_train = subsample_with_temporal_spacing(
        df_scaled,
        id_cols=["sample_name", "TrackID"],
        time_col="position_t",
        min_spacing=spacing_to_use,
        max_samples=max_samples,
        random_state=random_state,
    )
    
    print(f"Subsampled {len(df_train)} rows for clustering (spacing={spacing_to_use}).")
    if len(df_train) < 50:
        raise ValueError("Insufficient rows after global subsampling for clustering.")

    model_adata = cluster_group(
        df_group=df_train,
        feature_cols=kept_features,
        non_feature_cols=non_feature_cols,
        n_neighbors=n_neighbors,
        resolution=resolution,
        pca_var_selection=pca_var_selection,
        outfolder=Path(outfolder) if outfolder else None,
        group_name="all_data",
        clustering_method=clustering_method,
        resolutions=resolutions,
        leiden_seed_tries=leiden_seed_tries,
        leiden_subsample_tries=leiden_subsample_tries,
        random_state=random_state,
    )
    if model_adata is None:
        raise RuntimeError("Global clustering model could not be fit.")

    if "preprocessing" not in model_adata.uns:
        model_adata.uns["preprocessing"] = {}
    model_adata.uns["preprocessing"]["kept_features"] = list(kept_features)
    if scaler is not None:
        model_adata.uns["preprocessing"]["scaler"] = {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }

    mapping =  {
            "1": "dead",
            "2": "plastic_static",
            "3": "round_static",
            "4": "scanning",
            "5": "round_static", 
        }
    model_adata = relabel_cluster_ids(
        model_adata,
        mapping,
        cluster_key="ClusterID",
    )
    
    
    # Inference to full dataset
    # df_full = df_analysis.copy()
    # if scaler is not None:
    #     x = df_full[kept_features].to_numpy(dtype=float)
    #     df_full.loc[:, kept_features] = (x - scaler.mean_) / scaler.scale_

    obs_cols = [c for c in non_feature_cols if c in df_scaled.columns] + [c for c in binary_cols_to_merge if c in df_scaled.columns]
    adata_full = df_to_adata(df_scaled, kept_features, obs_cols=obs_cols)
    adata_full = adata_full[:, model_adata.var_names].copy()
    transfer_mode = str(label_transfer_method).lower()
    label_classifier = None
    classifier_fit_info = None
    classifier_qc_model_data = None
    classifier_model_path = None
    if transfer_mode == "ingest":
        sc.tl.ingest(adata_full, model_adata, obs="ClusterID")
    elif transfer_mode == "classifier":
        label_classifier, classifier_fit_info = train_clusterid_classifier_from_model_adata(
            model_adata=model_adata,
            cluster_col="ClusterID",
            classifier_backend=classifier_backend,
            classifier_kwargs=classifier_kwargs,
            random_state=random_state,
        )
        predict_clusterids_with_classifier(
            adata=adata_full,
            classifier=label_classifier,
            feature_cols=list(model_adata.var_names),
            output_col="ClusterID",
            confidence_col=classifier_confidence_col,
            inplace=True,
        )
        classifier_qc_model_data = evaluate_clusterid_classifier_on_model_adata(
            model_adata=model_adata,
            classifier=label_classifier,
            cluster_col="ClusterID",
            outfolder=(Path(outfolder) if outfolder is not None else None),
            filename_prefix=f"state_classification_v2_classifier_qc_{classifier_backend}",
        )
    else:
        raise ValueError("label_transfer_method must be one of: 'ingest', 'classifier'.")

    cluster_check_summary, _ = check_model_clusterid_consistency(
        model_adata=model_adata,
        adata_full=adata_full,
        model_cluster_col="ClusterID",
        full_cluster_col="ClusterID",
        key_cols=("sample_name", "TrackID", "position_t"),
    )
    
    adata_full.obs = _add_clean_binary_annotation_columns(adata_full.obs, binary_cols_to_merge)
    adata_full.obs["behavioral_clusterid"] = adata_full.obs["ClusterID"].astype(str)
    adata_full.obs["binary_group"] = _assign_binary_group_labels(adata_full.obs, binary_cols_to_merge).astype("category")
    adata_full.obs["ClusterID"] = (
        adata_full.obs["binary_group"].astype(str) + "_" + adata_full.obs["behavioral_clusterid"].astype(str)
    ).astype("category")
    adata_full.obs["behavioral_clusterid"] = adata_full.obs["behavioral_clusterid"].astype("category")

    adata_full.uns["state_classification_v2"] = {
        "window_size": int(window_size),
        "features": list(features),
        "descriptive_features": list(descriptive_features),
        "binary_features_to_group": list(binary_cols_to_merge),
        "clustering_method": clustering_method,
        "resolution": resolution,
        "n_neighbors": int(n_neighbors),
        "label_transfer_method": transfer_mode,
        "classifier_backend": str(classifier_backend) if transfer_mode == "classifier" else None,
        "classifier_fit_info": classifier_fit_info,
        "classifier_qc_model_data": classifier_qc_model_data,
        "classifier_confidence_col": classifier_confidence_col if transfer_mode == "classifier" else None,
        "classifier_model_path": None,
        "full_classifier_fit_info": None,
        "full_classifier_qc_model_data": None,
        "full_classifier_prediction_col": "ClusterID_full_classifier" if train_full_label_classifier else None,
        "full_classifier_confidence_col": full_classifier_confidence_col if train_full_label_classifier else None,
        "full_classifier_model_path": None,
        "model_vs_full_clusterid_check": cluster_check_summary,
    }

    mapping = {
    "no_contact_dead": "dead",
    "no_contact_plastic_static": "plastic_static",
    "no_contact_round_static": "round_static",
    "no_contact_scanning": "scanning",
    "organoid_contact_and_tcell_contact_dead": "dead",
    "organoid_contact_and_tcell_contact_plastic_static": "organoid_and_tcell_contact",
    "organoid_contact_and_tcell_contact_round_static": "organoid_and_tcell_contact",
    "organoid_contact_and_tcell_contact_scanning": "organoid_and_tcell_contact",
    "organoid_contact_dead": "dead",
    "organoid_contact_plastic_static": "organoid_contact",
    "organoid_contact_round_static": "organoid_contact",
    "organoid_contact_scanning": "organoid_contact",
    "tcell_contact_dead": "dead",
    "tcell_contact_plastic_static": "tcell_contact",
    "tcell_contact_round_static": "tcell_contact",
    "tcell_contact_scanning": "tcell_contact"
    }
    
    adata_full = relabel_cluster_ids(
        adata_full,
        mapping,
        cluster_key="ClusterID",
    )

    full_label_classifier_artifact = None
    full_classifier_fit_info = None
    full_classifier_qc_model_data = None
    full_classifier_model_path = None
    if train_full_label_classifier:
        full_binary_feature_cols = _resolve_binary_classifier_feature_cols(
            adata=adata_full,
            binary_cols_to_merge=binary_cols_to_merge,
        )
        full_label_classifier_artifact, full_classifier_fit_info = train_full_clusterid_classifier_from_adata(
            adata=adata_full,
            target_col="ClusterID",
            binary_feature_cols=full_binary_feature_cols,
            classifier_kwargs=full_classifier_kwargs,
            random_state=random_state,
        )

        predict_clusterids_with_full_classifier(
            adata=adata_full,
            classifier_artifact=full_label_classifier_artifact,
            output_col="ClusterID_full_classifier",
            confidence_col=full_classifier_confidence_col,
            inplace=True,
        )
        y_full_true = adata_full.obs["ClusterID"].astype(str).to_numpy()
        y_full_pred = adata_full.obs["ClusterID_full_classifier"].astype(str).to_numpy()
        full_classifier_qc_model_data = evaluate_classifier_predictions(
            y_true=y_full_true,
            y_pred=y_full_pred,
            outfolder=(Path(outfolder) if outfolder is not None else None),
            filename_prefix="state_classification_v2_full_classifier_qc_random_forest",
        )
        print(
            "Full classifier QC on adata_full labels: "
            f"n={full_classifier_qc_model_data['n_rows']}, "
            f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
            f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
        )
        adata_full.uns["state_classification_v2"]["full_classifier_fit_info"] = full_classifier_fit_info
        adata_full.uns["state_classification_v2"]["full_classifier_qc_model_data"] = full_classifier_qc_model_data

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=False, exist_ok=True)

        if transfer_mode == "classifier" and save_label_classifier and label_classifier is not None:
            classifier_model_path = outfolder / f"state_classification_v2_label_classifier_{classifier_backend}.pkl"
            with open(classifier_model_path, "wb") as f:
                pickle.dump(label_classifier, f)
            adata_full.uns["state_classification_v2"]["classifier_model_path"] = str(classifier_model_path)
        if train_full_label_classifier and save_full_label_classifier and full_label_classifier_artifact is not None:
            full_classifier_model_path = outfolder / "state_classification_v2_full_label_classifier_random_forest.pkl"
            with open(full_classifier_model_path, "wb") as f:
                pickle.dump(full_label_classifier_artifact, f)
            adata_full.uns["state_classification_v2"]["full_classifier_model_path"] = str(full_classifier_model_path)

        adata_full.write(outfolder / "adata_state_classification_v2_full.h5ad", compression="gzip")
        model_adata.write(outfolder / "adata_state_classification_v2_model.h5ad", compression="gzip")

        pdf1 = outfolder / "state_classification_v2_binary_group_cluster_proportions.pdf"
        fig, _ = plot_binary_group_behavioral_cluster_grid(
            adata_full,
            pdf_path=pdf1,
            return_csv=True,
        )
        plt.close(fig)
        pdf2 = outfolder / "state_classification_v2_cluster_binary_group_proportions.pdf"
        fig2, _ = plot_behavioral_cluster_binary_group_grid(
            adata_full,
            pdf_path=pdf2,
            return_csv=True,
        )
        plt.close(fig2)

    if return_artifacts:
        return {
            "adata_full": adata_full,
            "model_adata": model_adata,
            "label_classifier": label_classifier,
            "classifier_fit_info": classifier_fit_info,
            "classifier_qc_model_data": classifier_qc_model_data,
            "full_label_classifier_artifact": full_label_classifier_artifact,
            "full_classifier_fit_info": full_classifier_fit_info,
            "full_classifier_qc_model_data": full_classifier_qc_model_data,
            "cluster_check_summary": cluster_check_summary,
        }

    return adata_full
