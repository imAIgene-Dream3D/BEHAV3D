import copy
import numbers
import pickle
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from sklearn.preprocessing import StandardScaler

from behav3d.core.anndata import df_to_adata
from behav3d.features.rolling_window_features import create_descriptive_track_dataset
from behav3d.analysis.behavior.state.hmm import (
    run_hmm_state_classification,
    run_sticky_hmm_state_classification,
)
from behav3d.analysis.behavior.state.utils import (
    A4_LANDSCAPE,
    _apply_log1p_to_feature_matrix,
    _get_classification_state_colors,
    _infer_binary_group_constraints,
    _invert_log_scaling_in_continuous_matrix,
    _mixed_label_sort_key,
    _normalize_label_color_map,
    _normalize_log_scale_feature_selectors,
    _rebuild_full_behavioral_cluster_from_intrinsic,
    _require_columns,
    _resolve_log_scale_feature_cols,
    _resolve_positions_csv_path,
    _resolve_state_paths,
    _save_pdf_page_a4,
    _sanitize_filename_token,
    _set_classification_state_colors,
    _to_numpy_2d,
    _vdone,
    _vinfo,
    _vsave,
    _vstart,
    cap_values_to_quantile,
)
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    save_state_composition_report,
)
from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
    save_state_transition_report,
)


HMM_OPTIONAL_WINDOW_FEATURES = (
    "net_displacement",
    "straightness",
    "mean_square_displacement",
)
INTRINSIC_STATE_COL = "hmm_intrinsic_behavioral_state"
HMM_INTRINSIC_RAW_STATE_COL = "hmm_intrinsic_behavioral_state_raw"
FULL_STATE_COL = "behavioral_state"
BINARY_GROUP_COL = "binary_group"


def _obs_label_values(adata, col):
    if adata is None or not hasattr(adata, "obs") or col not in adata.obs.columns:
        return []
    labels = pd.Series(adata.obs[col]).dropna().astype("string").str.strip()
    labels = labels[labels != ""]
    return sorted([str(v) for v in labels.unique().tolist()], key=_mixed_label_sort_key)


def _extract_state_color_mapping(colors, state_col):
    if not isinstance(colors, dict):
        return {}
    nested = colors.get(str(state_col), None)
    if isinstance(nested, dict):
        return {str(k): str(v) for k, v in nested.items()}
    if all(not isinstance(v, dict) for v in colors.values()):
        return {str(k): str(v) for k, v in colors.items()}
    return {}


def _current_hmm_full_state_colors(adata, state_col=FULL_STATE_COL):
    labels = _obs_label_values(adata, state_col)
    colors = _normalize_label_color_map(
        labels,
        colors=_get_classification_state_colors(adata, state_col),
    )
    if len(colors) > 0:
        _set_classification_state_colors(adata, state_col, colors)
    return colors


def _normalize_selection(selection):
    if selection is None:
        return []
    return [str(v) for v in list(selection)]


def _dedupe_preserve_order(values):
    out = []
    seen = set()
    if values is None:
        values = []
    for value in list(values):
        value = str(value)
        if value in seen:
            continue
        out.append(value)
        seen.add(value)
    return out


def _coerce_hmm_raw_state_series(raw_series, *, label):
    series = pd.Series(raw_series)
    stringified = series.astype("string").str.strip()
    empty_mask = stringified.isna() | (stringified == "")
    numeric = pd.to_numeric(series, errors="coerce")
    invalid_mask = (~empty_mask) & numeric.isna()
    if bool(invalid_mask.any()):
        bad_vals = sorted(set(stringified[invalid_mask].tolist()))
        raise ValueError(
            f"{label} contains non-numeric HMM raw state values: {bad_vals[:10]}"
        )
    numeric_valid = numeric[~numeric.isna()]
    if not numeric_valid.empty:
        fractional = np.mod(numeric_valid.astype(float), 1.0)
        non_integer_mask = ~np.isclose(fractional, 0.0)
        if bool(np.any(non_integer_mask)):
            bad_vals = sorted(set(numeric_valid[non_integer_mask].tolist()))
            raise ValueError(
                f"{label} contains non-integer HMM raw state values: {bad_vals[:10]}"
            )
    coerced = numeric.round(0).astype("Int64")
    coerced[empty_mask] = pd.NA
    coerced.index = series.index
    return coerced


def _format_hmm_raw_state_series_for_key(raw_series, *, label):
    raw_int = _coerce_hmm_raw_state_series(raw_series, label=label)
    return raw_int.astype("string").str.strip()


def _format_hmm_raw_state_value_for_key(raw_state, *, label):
    raw_int = _coerce_hmm_raw_state_series(pd.Series([raw_state]), label=label)
    value = raw_int.iloc[0]
    if pd.isna(value):
        return ""
    return str(int(value))


def _validate_hmm_intrinsic_label_mapping_keys(intrinsic_label_mapping, artifact_name):
    bad_keys = [
        k
        for k in intrinsic_label_mapping.keys()
        if isinstance(k, bool) or not isinstance(k, numbers.Integral)
    ]
    if bad_keys:
        raise ValueError(
            f"{artifact_name} intrinsic_label_mapping keys must be integers; "
            f"rebuild the artifact. Invalid keys: {bad_keys[:10]}"
        )
    return intrinsic_label_mapping


def _coerce_hmm_log_scaling_params(preprocessing_meta):
    if not isinstance(preprocessing_meta, dict):
        return None
    log_meta = preprocessing_meta.get("log_scaling", None)
    if not isinstance(log_meta, dict):
        return None

    requested = _normalize_log_scale_feature_selectors(log_meta.get("requested_features", []))
    resolved_raw = log_meta.get("resolved_feature_cols", None)
    resolved = [] if resolved_raw is None else [str(c) for c in list(resolved_raw)]
    transform = log_meta.get("transform", None)
    if transform is None and len(resolved) > 0:
        transform = "log1p"
    if transform is None and len(requested) == 0 and len(resolved) == 0:
        return None
    return {
        "transform": None if transform is None else str(transform),
        "requested_features": list(requested),
        "resolved_feature_cols": list(resolved),
    }


def _coerce_hmm_scaler_params(preprocessing_meta):
    if not isinstance(preprocessing_meta, dict):
        return None
    scaler_meta = preprocessing_meta.get("scaler", None)
    if not isinstance(scaler_meta, dict):
        return None
    mean_raw = scaler_meta.get("mean", None)
    scale_raw = scaler_meta.get("scale", None)
    if mean_raw is None or scale_raw is None:
        return None
    return {
        "mean": np.asarray(mean_raw, dtype=float),
        "scale": np.asarray(scale_raw, dtype=float),
    }


def _apply_hmm_saved_preprocessing_to_matrix(
    X,
    *,
    preprocessing_meta,
    feature_cols,
):
    out = np.asarray(X, dtype=float).copy()
    feature_cols = [str(c) for c in list(feature_cols)]

    qmeta = preprocessing_meta.get("quantile_capping", {}) if isinstance(preprocessing_meta, dict) else {}
    feature_limits = qmeta.get("feature_limits", {}) if isinstance(qmeta, dict) else {}
    if isinstance(feature_limits, dict) and len(feature_limits) > 0:
        for idx, feature_name in enumerate(feature_cols):
            limits = feature_limits.get(feature_name, None)
            if not isinstance(limits, dict):
                continue
            lower = limits.get("lower", None)
            upper = limits.get("upper", None)
            out[:, idx] = np.clip(out[:, idx], a_min=lower, a_max=upper)

    log_meta = _coerce_hmm_log_scaling_params(preprocessing_meta)
    if isinstance(log_meta, dict) and len(log_meta.get("resolved_feature_cols", [])) > 0:
        out = _apply_log1p_to_feature_matrix(
            out,
            feature_cols=feature_cols,
            resolved_feature_cols=log_meta["resolved_feature_cols"],
            inplace=True,
        )

    scaler_meta = _coerce_hmm_scaler_params(preprocessing_meta)
    if isinstance(scaler_meta, dict):
        mean = np.asarray(scaler_meta["mean"], dtype=float)
        scale = np.asarray(scaler_meta["scale"], dtype=float)
        if out.shape[1] != len(mean) or out.shape[1] != len(scale):
            raise ValueError(
                "HMM scaler parameter size mismatch: "
                f"n_features={out.shape[1]}, len(mean)={len(mean)}, len(scale)={len(scale)}"
            )
        safe_scale = scale.copy()
        safe_scale[safe_scale == 0] = 1.0
        out -= mean
        out /= safe_scale
    return out


def _resolve_hmm_deployment_artifact_path(output_dir=None, cell_type=None, state_paths=None):
    if state_paths is None:
        if cell_type is None or len(str(cell_type).strip()) == 0:
            raise ValueError("cell_type is required when state_paths is not provided.")
        state_paths = _resolve_state_paths(output_dir, cell_type)
    resolved_cell_type = (
        str(cell_type).strip()
        if cell_type is not None and len(str(cell_type).strip()) > 0
        else str(Path(state_paths.analysis_outdir).name)
    )
    return (
        Path(state_paths.processing_outdir)
        / "hmm_behavioral_classification"
        / f"hmm_state_deployment_{resolved_cell_type}.pkl"
    )


def _resolve_hmm_quality_control_outdir(output_dir=None, cell_type=None, state_paths=None):
    if state_paths is None:
        if cell_type is None or len(str(cell_type).strip()) == 0:
            raise ValueError("cell_type is required when state_paths is not provided.")
        state_paths = _resolve_state_paths(output_dir, cell_type)
    outdir = Path(state_paths.processing_outdir) / "hmm_behavioral_classification" / "quality_control"
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def _build_hmm_intrinsic_label_mapping(
    model_adata,
    *,
    hmm_model=None,
    intrinsic_state_col=HMM_INTRINSIC_RAW_STATE_COL,
    label_col=INTRINSIC_STATE_COL,
):
    if model_adata is None or not hasattr(model_adata, "obs"):
        raise ValueError("model_adata with .obs is required to build an HMM label mapping.")
    if label_col not in model_adata.obs.columns:
        raise ValueError(f"model_adata is missing '{label_col}'.")
    if intrinsic_state_col == label_col:
        raise ValueError(
            "HMM intrinsic label mapping requires distinct raw and curated columns; "
            f"got intrinsic_state_col='{intrinsic_state_col}' and label_col='{label_col}'."
        )

    obs = model_adata.obs.copy()
    mapping = {}
    n_components = int(getattr(hmm_model, "n_components", 0) or 0)
    if n_components > 0:
        mapping = {int(idx): str(idx) for idx in range(1, n_components + 1)}

    if intrinsic_state_col not in obs.columns:
        raise ValueError(
            f"model_adata is missing '{intrinsic_state_col}', so curated intrinsic labels cannot be linked back "
            "to raw HMM state ids. Re-run HMM clustering with the updated pipeline."
        )

    work = obs[[intrinsic_state_col, label_col]].copy()
    work[label_col] = work[label_col].astype("string").str.strip().fillna("")
    work[intrinsic_state_col] = _coerce_hmm_raw_state_series(
        work[intrinsic_state_col],
        label=intrinsic_state_col,
    )
    work = work[work[intrinsic_state_col].notna() & (work[label_col] != "")]

    per_raw = (
        work.groupby(intrinsic_state_col, observed=False)[label_col]
        .agg(lambda s: sorted(set(str(v) for v in s if str(v) != ""), key=_mixed_label_sort_key))
    )
    for raw_state_id, labels in per_raw.items():
        if len(labels) != 1:
            raise ValueError(
                "Curated intrinsic label mapping is inconsistent for raw HMM state "
                f"'{raw_state_id}': labels={labels}"
            )
        mapping[int(raw_state_id)] = str(labels[0])
    return {int(k): str(v) for k, v in mapping.items()}


def _build_hmm_full_label_mapping(
    model_adata,
    *,
    full_label_col=FULL_STATE_COL,
    binary_group_col=BINARY_GROUP_COL,
    raw_state_col=HMM_INTRINSIC_RAW_STATE_COL,
):
    if model_adata is None or not hasattr(model_adata, "obs"):
        raise ValueError("model_adata with .obs is required to build a full-label mapping.")
    obs = model_adata.obs.copy()
    if full_label_col not in obs.columns:
        raise ValueError(f"model_adata is missing '{full_label_col}'.")
    if binary_group_col not in obs.columns:
        raise ValueError(f"model_adata is missing '{binary_group_col}'.")
    if raw_state_col not in obs.columns:
        raise ValueError(f"model_adata is missing '{raw_state_col}'.")

    work = obs[[binary_group_col, raw_state_col, full_label_col]].copy()
    work[binary_group_col] = work[binary_group_col].astype("string").str.strip().fillna("")
    work[raw_state_col] = _coerce_hmm_raw_state_series(
        work[raw_state_col],
        label=raw_state_col,
    )
    work[full_label_col] = work[full_label_col].astype("string").str.strip().fillna("")
    missing_key_mask = (
        (work[full_label_col] != "")
        & ((work[binary_group_col] == "") | work[raw_state_col].isna())
    )
    if bool(missing_key_mask.any()):
        raise ValueError(
            "Curated full labels cannot be anchored to original HMM labels because some rows are missing "
            f"'{binary_group_col}' or '{raw_state_col}'."
        )
    work = work[
        (work[binary_group_col] != "")
        & work[raw_state_col].notna()
        & (work[full_label_col] != "")
    ].copy()
    raw_state_key = _format_hmm_raw_state_series_for_key(work[raw_state_col], label=raw_state_col)
    work["_raw_full_cluster_key"] = work[binary_group_col].astype(str) + "_" + raw_state_key.astype(str)

    per_generated = (
        work.groupby("_raw_full_cluster_key", observed=False)[full_label_col]
        .agg(lambda s: sorted(set(str(v) for v in s if str(v) != ""), key=_mixed_label_sort_key))
    )
    mapping = {}
    for raw_full_key, labels in per_generated.items():
        if len(labels) != 1:
            raise ValueError(
                "Curated full label mapping is inconsistent for canonical raw full cluster "
                f"'{raw_full_key}': labels={labels}"
            )
        mapping[str(raw_full_key)] = str(labels[0])
    return mapping


def _build_hmm_raw_full_cluster_keys(
    obs,
    *,
    binary_group_col=BINARY_GROUP_COL,
    raw_state_col=HMM_INTRINSIC_RAW_STATE_COL,
):
    if binary_group_col not in obs.columns:
        raise ValueError(f"obs is missing '{binary_group_col}'.")
    if raw_state_col not in obs.columns:
        raise ValueError(f"obs is missing '{raw_state_col}'.")

    binary_group = pd.Series(obs[binary_group_col], index=obs.index, dtype="string").str.strip()
    raw_state = _format_hmm_raw_state_series_for_key(
        pd.Series(obs[raw_state_col], index=obs.index),
        label=raw_state_col,
    )
    raw_full_keys = pd.Series(pd.NA, index=obs.index, dtype="string")
    valid = binary_group.notna() & raw_state.notna() & (binary_group != "") & (raw_state != "")
    raw_full_keys.loc[valid] = (
        binary_group.loc[valid].astype(str) + "_" + raw_state.loc[valid].astype(str)
    )
    return raw_full_keys


def _expand_hmm_full_label_mapping_with_raw_state_keys(
    full_label_mapping,
    intrinsic_label_mapping,
    binary_groups,
):
    mapping = {str(k): str(v) for k, v in dict(full_label_mapping).items()}
    for binary_group in binary_groups:
        binary_group = str(binary_group).strip()
        if binary_group == "":
            continue
        for raw_state, intrinsic_label in dict(intrinsic_label_mapping).items():
            intrinsic_label = str(intrinsic_label).strip()
            raw_state_key = _format_hmm_raw_state_value_for_key(
                raw_state,
                label="intrinsic_label_mapping",
            )
            if raw_state_key == "" or intrinsic_label == "":
                continue
            raw_key = f"{binary_group}_{raw_state_key}"
            intrinsic_key = f"{binary_group}_{intrinsic_label}"
            if raw_key not in mapping and intrinsic_key in mapping:
                mapping[raw_key] = mapping[intrinsic_key]
    return mapping


def _rebuild_hmm_full_behavioral_state_from_intrinsic(
    adata,
    *,
    binary_cols_to_merge,
    intrinsic_col=INTRINSIC_STATE_COL,
    full_state_col=FULL_STATE_COL,
    binary_group_constraints=None,
    enforce_binary_group_constraints=False,
):
    _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col=intrinsic_col,
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=enforce_binary_group_constraints,
    )
    if "full_behavioral_cluster" not in adata.obs.columns:
        raise ValueError("Could not rebuild HMM full behavioral state labels.")
    adata.obs[full_state_col] = pd.Categorical(
        pd.Series(adata.obs["full_behavioral_cluster"], index=adata.obs.index, dtype="string")
    )
    if intrinsic_col != "intrinsic_behavioral_cluster":
        adata.obs["intrinsic_behavioral_cluster"] = adata.obs[intrinsic_col].copy()
    return adata


def _build_hmm_deployment_artifact(
    *,
    model_adata,
    hmm_model,
    state_paths=None,
    output_dir=None,
    cell_type=None,
    source_model_adata_path=None,
):
    if hmm_model is None:
        raise ValueError("hmm_model is required to build an HMM deployment artifact.")
    if model_adata is None or not hasattr(model_adata, "obs"):
        raise ValueError("model_adata with .obs is required to build an HMM deployment artifact.")

    if state_paths is None:
        if cell_type is None or len(str(cell_type).strip()) == 0:
            raise ValueError("cell_type is required when state_paths is not provided.")
        state_paths = _resolve_state_paths(output_dir, cell_type)

    pipeline_metadata = {
        "preprocessing": copy.deepcopy(dict(getattr(model_adata, "uns", {}).get("preprocessing", {}))),
        "clustering": copy.deepcopy(dict(getattr(model_adata, "uns", {}).get("clustering", {}))),
        "classification": copy.deepcopy(dict(getattr(model_adata, "uns", {}).get("classification", {}))),
    }
    clustering_meta = pipeline_metadata["clustering"]
    binary_cols_to_merge = [str(c) for c in list(clustering_meta.get("binary_cols_to_merge", []))]
    binary_group_constraints = copy.deepcopy(clustering_meta.get("binary_group_constraints", None))
    observed_binary_groups = []
    if BINARY_GROUP_COL in model_adata.obs.columns:
        observed_binary_groups = sorted(
            [
                str(v)
                for v in (
                    pd.Series(model_adata.obs[BINARY_GROUP_COL])
                    .dropna()
                    .astype("string")
                    .str.strip()
                    .unique()
                    .tolist()
                )
                if v != ""
            ],
            key=_mixed_label_sort_key,
        )

    source_path = (
        str(Path(source_model_adata_path))
        if source_model_adata_path is not None
        else str(Path(state_paths.model_adata_path))
    )
    preprocessing_meta = pipeline_metadata["preprocessing"]
    raw_feature_cols = list(preprocessing_meta.get("raw_feature_cols", []))
    observation_feature_cols = list(preprocessing_meta.get("continuous_feature_cols", []))
    if len(observation_feature_cols) == 0:
        observation_feature_cols = list(getattr(model_adata, "var_names", []))
    if len(observation_feature_cols) == 0:
        raise ValueError("Could not determine HMM observation feature columns for deployment artifact.")

    full_state_colors = _current_hmm_full_state_colors(model_adata, state_col=FULL_STATE_COL)
    classification_meta = pipeline_metadata["classification"]
    existing_state_colors = classification_meta.get("state_colors", {})
    if not isinstance(existing_state_colors, dict):
        existing_state_colors = {}
    existing_state_colors[str(FULL_STATE_COL)] = dict(full_state_colors)
    classification_meta["state_colors"] = existing_state_colors

    _coerce_hmm_raw_state_series(
        model_adata.obs[HMM_INTRINSIC_RAW_STATE_COL],
        label=HMM_INTRINSIC_RAW_STATE_COL,
    )

    return {
        "artifact_kind": "hmm_state_deployment",
        "model": hmm_model,
        "pipeline_metadata": pipeline_metadata,
        "observation_feature_cols": [str(c) for c in observation_feature_cols],
        "raw_feature_cols": [str(c) for c in raw_feature_cols],
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "binary_group_constraints": binary_group_constraints,
        "intrinsic_label_mapping": _build_hmm_intrinsic_label_mapping(
            model_adata,
            hmm_model=hmm_model,
            intrinsic_state_col=HMM_INTRINSIC_RAW_STATE_COL,
            label_col=INTRINSIC_STATE_COL,
        ),
        "full_label_mapping": _build_hmm_full_label_mapping(
            model_adata,
            full_label_col=FULL_STATE_COL,
            binary_group_col=BINARY_GROUP_COL,
            raw_state_col=HMM_INTRINSIC_RAW_STATE_COL,
        ),
        "full_label_strategy": "binary_group_plus_intrinsic",
        "observed_binary_groups": list(observed_binary_groups),
        "state_colors": {str(FULL_STATE_COL): dict(full_state_colors)},
        "source_model_adata_path": source_path,
    }


def _validate_hmm_deployment_artifact(artifact, artifact_name="hmm_deployment_artifact"):
    if not isinstance(artifact, dict):
        raise ValueError(f"{artifact_name} must be a dict, got {type(artifact)}.")
    if str(artifact.get("artifact_kind", "")) != "hmm_state_deployment":
        raise ValueError(
            f"{artifact_name} has invalid artifact_kind='{artifact.get('artifact_kind', None)}'; "
            "expected 'hmm_state_deployment'."
        )
    model = artifact.get("model", None)
    if model is None or (not hasattr(model, "predict")) or (not hasattr(model, "predict_proba")):
        raise ValueError(f"{artifact_name} is missing a valid fitted HMM model.")
    if int(getattr(model, "n_components", 0) or 0) < 1:
        raise ValueError(f"{artifact_name} model must expose n_components >= 1.")

    pipeline_metadata = artifact.get("pipeline_metadata", None)
    if not isinstance(pipeline_metadata, dict):
        raise ValueError(f"{artifact_name} is missing required pipeline_metadata.")
    for key in ("preprocessing", "clustering", "classification"):
        if not isinstance(pipeline_metadata.get(key, None), dict):
            raise ValueError(f"{artifact_name} pipeline_metadata.{key} is missing/invalid.")

    observation_feature_cols = artifact.get("observation_feature_cols", None)
    if not isinstance(observation_feature_cols, list) or len(observation_feature_cols) == 0:
        raise ValueError(f"{artifact_name} observation_feature_cols must be a non-empty list.")
    raw_feature_cols = artifact.get("raw_feature_cols", None)
    if not isinstance(raw_feature_cols, list):
        raise ValueError(f"{artifact_name} raw_feature_cols must be a list.")
    binary_cols_to_merge = artifact.get("binary_cols_to_merge", None)
    if not isinstance(binary_cols_to_merge, list):
        raise ValueError(f"{artifact_name} binary_cols_to_merge must be a list.")
    intrinsic_label_mapping = artifact.get("intrinsic_label_mapping", None)
    if not isinstance(intrinsic_label_mapping, dict) or len(intrinsic_label_mapping) == 0:
        raise ValueError(f"{artifact_name} intrinsic_label_mapping must be a non-empty dict.")
    _validate_hmm_intrinsic_label_mapping_keys(intrinsic_label_mapping, artifact_name)
    full_label_mapping = artifact.get("full_label_mapping", None)
    if not isinstance(full_label_mapping, dict) or len(full_label_mapping) == 0:
        raise ValueError(f"{artifact_name} full_label_mapping must be a non-empty dict.")
    if str(artifact.get("full_label_strategy", "")) != "binary_group_plus_intrinsic":
        raise ValueError(
            f"{artifact_name} full_label_strategy must be 'binary_group_plus_intrinsic', "
            f"got '{artifact.get('full_label_strategy', None)}'."
        )
    observed_binary_groups = artifact.get("observed_binary_groups", None)
    if not isinstance(observed_binary_groups, list):
        raise ValueError(f"{artifact_name} observed_binary_groups must be a list.")
    state_colors = artifact.get("state_colors", {})
    if state_colors is not None and not isinstance(state_colors, dict):
        raise ValueError(f"{artifact_name} state_colors must be a dict when present.")
    return artifact


def load_hmm_deployment_artifact(path):
    with open(path, "rb") as f:
        artifact = pickle.load(f)
    _validate_hmm_deployment_artifact(artifact, artifact_name=f"HMM deployment artifact '{path}'")
    return artifact


def save_hmm_deployment_artifact(
    *,
    output_path,
    model_adata,
    hmm_model,
    state_paths=None,
    output_dir=None,
    cell_type=None,
    source_model_adata_path=None,
    verbose=True,
):
    artifact = _build_hmm_deployment_artifact(
        model_adata=model_adata,
        hmm_model=hmm_model,
        state_paths=state_paths,
        output_dir=output_dir,
        cell_type=cell_type,
        source_model_adata_path=source_model_adata_path,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as f:
        pickle.dump(artifact, f)
    _vsave(verbose, "state-hmm", "HMM deployment artifact", output_path)
    return artifact


def _coerce_hmm_deployment_artifact_input(hmm_deployment_artifact, artifact_name):
    if isinstance(hmm_deployment_artifact, (str, Path)):
        resolved_path = Path(hmm_deployment_artifact)
        return load_hmm_deployment_artifact(resolved_path), "explicit_path", str(resolved_path)
    _validate_hmm_deployment_artifact(hmm_deployment_artifact, artifact_name=artifact_name)
    return hmm_deployment_artifact, "explicit_artifact", None


def _summarize_hmm_apply_row_drops(
    df_base,
    valid_mask,
    *,
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    max_preview=5,
):
    total_rows = int(len(df_base))
    kept_rows = int(np.count_nonzero(np.asarray(valid_mask, dtype=bool)))
    dropped_rows = int(total_rows - kept_rows)
    dropped_pct = float((100.0 * dropped_rows / total_rows) if total_rows > 0 else 0.0)

    key_cols = [str(sample_col), str(track_col)]
    dropped_rows_df = df_base.loc[~pd.Series(valid_mask, index=df_base.index), key_cols + [str(time_col)]].copy()
    input_tracks = df_base[key_cols].drop_duplicates()
    kept_tracks = df_base.loc[pd.Series(valid_mask, index=df_base.index), key_cols].drop_duplicates()
    dropped_tracks = dropped_rows_df[key_cols].drop_duplicates()

    preview = [
        f"{str(row[sample_col])}/{str(row[track_col])}"
        for _, row in dropped_tracks.head(int(max_preview)).iterrows()
    ]
    omitted = max(0, int(len(dropped_tracks)) - len(preview))

    return {
        "n_timepoints_input": total_rows,
        "n_timepoints_kept": kept_rows,
        "n_timepoints_dropped_nan": dropped_rows,
        "dropped_pct": dropped_pct,
        "n_tracks_input": int(len(input_tracks)),
        "n_tracks_kept": int(len(kept_tracks)),
        "n_tracks_with_dropped_timepoints": int(len(dropped_tracks)),
        "dropped_track_preview": preview,
        "dropped_track_preview_omitted": int(omitted),
    }


def _warn_hmm_apply_row_drops(drop_summary):
    if int(drop_summary.get("n_timepoints_dropped_nan", 0)) <= 0:
        return

    preview = list(drop_summary.get("dropped_track_preview", []))
    preview_text = ", ".join(preview) if len(preview) > 0 else "none"
    omitted = int(drop_summary.get("dropped_track_preview_omitted", 0))
    if omitted > 0:
        preview_text = f"{preview_text} (+{omitted} more)"

    warnings.warn(
        "HMM apply removed "
        f"{int(drop_summary.get('n_timepoints_dropped_nan', 0))}/"
        f"{int(drop_summary.get('n_timepoints_input', 0))} rows "
        f"({float(drop_summary.get('dropped_pct', 0.0)):.1f}%) before state assignment. "
        "Rows were removed because one or more required HMM observation features became NaN "
        "after smoothing, optional window-feature construction, or numeric coercion. "
        f"Kept rows={int(drop_summary.get('n_timepoints_kept', 0))}. "
        f"Affected tracks={int(drop_summary.get('n_tracks_with_dropped_timepoints', 0))}/"
        f"{int(drop_summary.get('n_tracks_input', 0))}. "
        f"Preview: {preview_text}.",
        RuntimeWarning,
        stacklevel=2,
    )


def _validate_n_states(n_states):
    if n_states is None:
        raise ValueError("n_states is required for HMM state clustering. Pass an integer or 'auto'.")
    if isinstance(n_states, str):
        if n_states.lower() != "auto":
            raise ValueError("n_states must be an integer or 'auto'.")
        return "auto"
    n_states = int(n_states)
    if n_states < 1:
        raise ValueError("n_states must be >= 1.")
    return n_states


def _smooth_timepoint_features(
    df,
    feature_cols,
    *,
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    window=1,
    min_periods=1,
):
    window = max(1, int(window))
    min_periods = max(1, int(min_periods))
    if window <= 1 or len(feature_cols) == 0:
        return df.copy()

    sort_cols = list(id_cols) + [time_col]
    out = df.sort_values(sort_cols, kind="mergesort").copy()
    numeric = out[list(feature_cols)].apply(pd.to_numeric, errors="coerce")
    grouped = numeric.groupby([out[c] for c in id_cols], sort=False, observed=False)
    out.loc[:, list(feature_cols)] = grouped.transform(
        lambda g: g.rolling(window=window, min_periods=min_periods).mean()
    )
    return out


def _compute_optional_window_feature_frame(
    df,
    *,
    window_features,
    time_col="position_t",
    id_cols=("sample_name", "TrackID"),
    window_size=1,
    verbose=True,
):
    window_features = _dedupe_preserve_order(window_features)
    if len(window_features) == 0:
        keep_cols = list(id_cols) + [time_col]
        return pd.DataFrame(columns=keep_cols), 1

    unknown = [feat for feat in window_features if feat not in HMM_OPTIONAL_WINDOW_FEATURES]
    if len(unknown) > 0:
        raise ValueError(
            "Unsupported HMM additional_window_features: "
            f"{unknown}. Supported values: {list(HMM_OPTIONAL_WINDOW_FEATURES)}"
        )

    required_cols = list(id_cols) + [time_col, "position_x", "position_y", "position_z"]
    _require_columns(df, required_cols, "HMM window feature computation")
    effective_window = max(2, int(window_size))

    started = _vstart(
        verbose,
        "state-hmm",
        (
            "compute optional window features "
            f"| features={window_features}, window={effective_window}"
        ),
    )
    df_window = create_descriptive_track_dataset(
        df_tracks=df[required_cols].copy(),
        columns_to_summarize=[],
        window_size=effective_window,
        step_size=1,
        time_col=time_col,
        id_cols=id_cols,
        features_to_compute=window_features,
        only_nonbinary=True,
        incomplete_window_policy="partial",
    )
    keep_cols = list(id_cols) + [time_col] + window_features
    if len(df_window) == 0:
        df_window = pd.DataFrame(columns=keep_cols)
    else:
        df_window = df_window[keep_cols].copy()
    _vdone(verbose, "state-hmm", "compute optional window features", started)
    return df_window, effective_window


def _compute_hmm_sequence_lengths(df, *, id_cols=("sample_name", "TrackID")):
    if df is None or len(df) == 0:
        return []
    grouped = df.groupby(list(id_cols), sort=False, observed=False)
    return [int(len(group)) for _, group in grouped]


def _prepare_hmm_apply_adata_from_df_positions(
    df_positions,
    *,
    preprocessing_meta,
    binary_cols_to_merge,
    feature_cols,
    start_offset=None,
    start_offset_fill_mode=None,
    verbose=True,
):
    if not isinstance(preprocessing_meta, dict):
        raise ValueError("preprocessing_meta must be a dict for HMM apply.")

    windowing = preprocessing_meta.get("windowing", {})
    if not isinstance(windowing, dict):
        raise ValueError("HMM apply preprocessing.windowing must be a dict.")

    raw_feature_cols = [str(c) for c in list(windowing.get("features", feature_cols))]
    additional_window_features = _dedupe_preserve_order(windowing.get("additional_window_features", []))
    resolved_start_offset = max(0, int(windowing.get("start_offset", 0) if start_offset is None else start_offset))
    resolved_fill_mode = str(
        windowing.get("start_offset_fill_mode", "backfill")
        if start_offset_fill_mode is None
        else start_offset_fill_mode
    ).strip().lower()
    if resolved_fill_mode not in {"backfill", "leave_unassigned"}:
        raise ValueError(
            "start_offset_fill_mode must be 'backfill' or 'leave_unassigned', "
            f"got '{resolved_fill_mode}'."
        )
    required_cols = ["sample_name", "TrackID", "position_t"] + list(raw_feature_cols) + list(binary_cols_to_merge)
    if len(additional_window_features) > 0:
        required_cols.extend(["position_x", "position_y", "position_z"])
    required_cols = _dedupe_preserve_order(required_cols)
    _require_columns(df_positions, required_cols, context="HMM apply input")

    df_base = df_positions.sort_values(["sample_name", "TrackID", "position_t"], kind="mergesort").copy()
    df_base["_start_offset_track_idx"] = (
        df_base.groupby(["sample_name", "TrackID"], sort=False, observed=False).cumcount().astype(int)
    )
    if len(additional_window_features) > 0:
        df_window, _ = _compute_optional_window_feature_frame(
            df_base,
            window_features=additional_window_features,
            time_col="position_t",
            id_cols=("sample_name", "TrackID"),
            window_size=int(windowing.get("window_feature_window_size", windowing.get("feature_smoothing_window", 1))),
            verbose=verbose,
        )
        df_base = df_base.merge(
            df_window,
            on=["sample_name", "TrackID", "position_t"],
            how="left",
            validate="one_to_one",
        )

    df_base = _smooth_timepoint_features(
        df_base,
        raw_feature_cols,
        id_cols=("sample_name", "TrackID"),
        time_col="position_t",
        window=int(windowing.get("feature_smoothing_window", 1)),
        min_periods=int(windowing.get("smoothing_min_periods", 1)),
    )

    cont_cols = [str(c) for c in list(feature_cols)]
    numeric_frame = df_base[cont_cols].apply(pd.to_numeric, errors="coerce")
    valid_mask = numeric_frame.notna().all(axis=1)
    drop_summary = _summarize_hmm_apply_row_drops(
        df_base,
        valid_mask,
        sample_col="sample_name",
        track_col="TrackID",
        time_col="position_t",
    )
    _warn_hmm_apply_row_drops(drop_summary)
    if not bool(valid_mask.any()):
        raise ValueError("No valid timepoints remain after smoothing for HMM apply.")

    df_valid = df_base.loc[valid_mask].copy()
    df_valid.loc[:, cont_cols] = numeric_frame.loc[valid_mask, cont_cols].to_numpy(dtype=float)
    df_valid["_start_offset_scored"] = df_valid["_start_offset_track_idx"] >= int(resolved_start_offset)
    df_valid["_start_offset_skipped"] = ~df_valid["_start_offset_scored"]
    df_scored = df_valid[df_valid["_start_offset_scored"]].copy()

    input_track_sizes = df_base.groupby(["sample_name", "TrackID"], sort=False, observed=False).size()
    valid_track_sizes = df_valid.groupby(["sample_name", "TrackID"], sort=False, observed=False).size()
    start_offset_summary = {
        "start_offset": int(resolved_start_offset),
        "start_offset_fill_mode": str(resolved_fill_mode),
        "n_timepoints_skipped_start_offset": int(df_valid["_start_offset_skipped"].sum()),
        "n_timepoints_scored": int(df_valid["_start_offset_scored"].sum()),
        "n_tracks_too_short_after_offset": int((input_track_sizes <= int(resolved_start_offset)).sum()),
        "n_tracks_scored_after_offset": int(len(df_scored.groupby(["sample_name", "TrackID"], sort=False, observed=False))),
        "n_tracks_valid_for_offset": int(len(valid_track_sizes)),
    }

    _meta_obs_extra = [
        c for c in (["exp_nr", "well"] + [col for col in df_valid.columns if col.endswith("_line_condition")])
        if c in df_valid.columns and c not in set(binary_cols_to_merge)
    ]
    obs_cols = (
        ["sample_name", "TrackID", "position_t"]
        + list(binary_cols_to_merge)
        + _meta_obs_extra
        + ["_start_offset_track_idx", "_start_offset_scored", "_start_offset_skipped"]
    )
    adata_query = df_to_adata(df_valid, feature_cols=cont_cols, obs_cols=obs_cols)
    adata_query.uns["preprocessing"] = {
        **copy.deepcopy(preprocessing_meta),
        **dict(drop_summary),
        **dict(start_offset_summary),
    }
    return {
        "adata": adata_query,
        "df_valid": df_valid,
        "df_scored": df_scored,
        "drop_summary": dict(drop_summary),
        "start_offset_summary": dict(start_offset_summary),
        "resolved_start_offset": int(resolved_start_offset),
        "resolved_start_offset_fill_mode": str(resolved_fill_mode),
    }


def _compute_state_means(adata, cluster_col, feature_cols):
    feature_cols = [] if feature_cols is None else list(feature_cols)
    valid_features = [str(c) for c in feature_cols if str(c) in adata.var_names]
    if adata is None or adata.n_obs == 0 or cluster_col not in adata.obs.columns or len(valid_features) == 0:
        return pd.DataFrame()

    X = _to_numpy_2d(adata[:, valid_features].X).astype(float, copy=False)
    feature_df = pd.DataFrame(X, columns=valid_features, index=adata.obs.index)
    feature_df["_cluster"] = pd.Series(adata.obs[cluster_col], index=adata.obs.index, dtype="string")
    feature_df = feature_df[feature_df["_cluster"].notna()].copy()
    if feature_df.empty:
        return pd.DataFrame()
    cluster_order = sorted(feature_df["_cluster"].unique().tolist(), key=_mixed_label_sort_key)
    return feature_df.groupby("_cluster", observed=True)[valid_features].mean().reindex(cluster_order)


def _compute_cluster_correlation_df(adata, cluster_col, feature_cols):
    state_means = _compute_state_means(adata, cluster_col, feature_cols)
    if state_means.empty or state_means.shape[0] < 2 or state_means.shape[1] < 2:
        return pd.DataFrame()
    corr_df = state_means.transpose().corr()
    corr_df.index.name = cluster_col
    corr_df.columns.name = cluster_col
    return corr_df


def _compute_hmm_transition_matrix_df(model, *, raw_to_label=None, raw_weights=None):
    transmat = getattr(model, "transmat_", None)
    if transmat is None:
        return pd.DataFrame()
    transmat = np.asarray(transmat, dtype=float)
    if transmat.ndim != 2 or transmat.shape[0] == 0 or transmat.shape[0] != transmat.shape[1]:
        return pd.DataFrame()
    state_labels = [str(i + 1) for i in range(int(transmat.shape[0]))]
    if not isinstance(raw_to_label, dict) or len(raw_to_label) == 0:
        return pd.DataFrame(transmat, index=state_labels, columns=state_labels)

    label_map = {str(k): str(v) for k, v in raw_to_label.items()}
    labels = sorted({label_map.get(label, label) for label in state_labels}, key=_mixed_label_sort_key)
    weights = {str(k): float(v) for k, v in dict(raw_weights or {}).items()}
    counts = pd.DataFrame(0.0, index=labels, columns=labels)
    for i, raw_i in enumerate(state_labels):
        row_label = label_map.get(raw_i, raw_i)
        row_weight = weights.get(raw_i, 1.0)
        for j, raw_j in enumerate(state_labels):
            col_label = label_map.get(raw_j, raw_j)
            counts.loc[row_label, col_label] += float(transmat[i, j]) * float(row_weight)
    row_sums = counts.sum(axis=1).replace(0.0, np.nan)
    return counts.div(row_sums, axis=0).fillna(0.0)


def _build_hmm_transition_label_mapping(adata, *, cluster_col):
    if (
        adata is None
        or not hasattr(adata, "obs")
        or HMM_INTRINSIC_RAW_STATE_COL not in adata.obs.columns
        or cluster_col not in adata.obs.columns
    ):
        return {}, {}

    work = adata.obs[[HMM_INTRINSIC_RAW_STATE_COL, cluster_col]].copy()
    work[HMM_INTRINSIC_RAW_STATE_COL] = _format_hmm_raw_state_series_for_key(
        work[HMM_INTRINSIC_RAW_STATE_COL],
        label=HMM_INTRINSIC_RAW_STATE_COL,
    )
    work[cluster_col] = work[cluster_col].astype("string").str.strip().fillna("")
    work = work[(work[HMM_INTRINSIC_RAW_STATE_COL] != "") & (work[cluster_col] != "")]
    if work.empty:
        return {}, {}

    grouped = work.groupby(HMM_INTRINSIC_RAW_STATE_COL, observed=False)[cluster_col].agg(
        lambda s: sorted({str(v) for v in s if str(v) != ""}, key=_mixed_label_sort_key)
    )
    mapping = {
        str(raw_state): str(labels[0])
        for raw_state, labels in grouped.items()
        if len(labels) == 1
    }
    weights = work[HMM_INTRINSIC_RAW_STATE_COL].value_counts().to_dict()
    return mapping, {str(k): float(v) for k, v in weights.items()}


def save_hmm_quality_control_outputs(
    adata,
    *,
    feature_cols,
    output_dir,
    model=None,
    selection_df=None,
    cluster_col=INTRINSIC_STATE_COL,
    scaler_mean=None,
    scaler_scale=None,
    title=None,
    preprocessing_params=None,
    verbose=True,
):
    """Write HMM QC PDFs/CSVs for the requested state labels."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_cols = [str(c) for c in list(feature_cols or [])]
    preprocessing_meta = (
        dict(preprocessing_params)
        if isinstance(preprocessing_params, dict)
        else getattr(adata, "uns", {}).get("preprocessing", {})
    )
    if scaler_mean is None and isinstance(preprocessing_meta, dict):
        scaler_mean = preprocessing_meta.get("scaler", {}).get("mean", None)
    if scaler_scale is None and isinstance(preprocessing_meta, dict):
        scaler_scale = preprocessing_meta.get("scaler", {}).get("scale", None)
    if title is None:
        n_states = getattr(model, "n_components", None)
        state_text = f"states={int(n_states)}" if n_states is not None else f"cluster_col={cluster_col}"
        title = f"all_data | hmm ({state_text})"
    raw_to_label, raw_weights = _build_hmm_transition_label_mapping(adata, cluster_col=cluster_col)
    transition_df = _compute_hmm_transition_matrix_df(
        model,
        raw_to_label=raw_to_label,
        raw_weights=raw_weights,
    )

    diagnostics_pdf = output_dir / "behavioral_clustering_diagnostics.pdf"
    diagnostics_started = _vstart(verbose, "state-hmm", "write HMM diagnostics PDF")
    _plot_hmm_state_diagnostics_pdf(
        adata,
        cluster_col=cluster_col,
        feature_cols=feature_cols,
        pdf_path=diagnostics_pdf,
        model=model,
        transition_df=transition_df,
        selection_df=selection_df,
        title=title,
        verbose=verbose,
    )
    _vdone(verbose, "state-hmm", "write HMM diagnostics PDF", diagnostics_started)
    _vsave(verbose, "state-hmm", "diagnostics pdf", diagnostics_pdf)

    diagnostics_csvs_started = _vstart(verbose, "state-hmm", "export HMM diagnostics CSVs")
    diagnostics_csvs = {}
    state_means = _compute_state_means(adata, cluster_col, feature_cols)
    if not state_means.empty:
        state_means_csv = output_dir / "behavioral_clustering_hmm_state_means.csv"
        state_means.to_csv(state_means_csv)
        diagnostics_csvs["state_means_csv"] = str(state_means_csv)
        _vsave(verbose, "state-hmm", "state means csv", state_means_csv)

    cluster_corr_df = _compute_cluster_correlation_df(adata, cluster_col, feature_cols)
    if not cluster_corr_df.empty:
        cluster_corr_csv = output_dir / "behavioral_clustering_hmm_cluster_correlation.csv"
        cluster_corr_df.to_csv(cluster_corr_csv)
        diagnostics_csvs["cluster_correlation_csv"] = str(cluster_corr_csv)
        _vsave(verbose, "state-hmm", "cluster correlation csv", cluster_corr_csv)

    transition_csv = None
    if not transition_df.empty:
        transition_csv = output_dir / "behavioral_clustering_hmm_transition_matrix.csv"
        transition_df.to_csv(transition_csv)
        diagnostics_csvs["transition_matrix_csv"] = str(transition_csv)
        _vsave(verbose, "state-hmm", "transition matrix csv", transition_csv)
    _vdone(verbose, "state-hmm", "export HMM diagnostics CSVs", diagnostics_csvs_started)

    state_counts_csv = output_dir / "behavioral_clustering_hmm_state_counts.csv"
    (
        pd.Series(adata.obs[cluster_col], dtype="string")
        .dropna()
        .value_counts()
        .sort_index(key=lambda idx: [_mixed_label_sort_key(v) for v in idx])
        .rename_axis(cluster_col)
        .reset_index(name="count")
        .to_csv(state_counts_csv, index=False)
    )
    _vsave(verbose, "state-hmm", "state counts csv", state_counts_csv)

    selection_csv = None
    if isinstance(selection_df, pd.DataFrame) and (not selection_df.empty):
        selection_csv = output_dir / "behavioral_clustering_hmm_model_selection.csv"
        selection_df.to_csv(selection_csv, index=False)
        diagnostics_csvs["model_selection_csv"] = str(selection_csv)
        _vsave(verbose, "state-hmm", "HMM model selection csv", selection_csv)

    feature_distribution_pdf = output_dir / "behavioral_clustering_feature_distributions.pdf"
    _plot_hmm_feature_distribution_pdf(
        adata,
        feature_cols=feature_cols,
        pdf_path=feature_distribution_pdf,
        scaler_mean=scaler_mean,
        scaler_scale=scaler_scale,
        cluster_col=cluster_col,
        title=title,
        preprocessing_params=preprocessing_meta,
        verbose=verbose,
    )
    _vsave(verbose, "state-hmm", "feature distribution pdf", feature_distribution_pdf)

    return {
        "diagnostics_pdf": str(diagnostics_pdf),
        "diagnostics_csvs": dict(diagnostics_csvs),
        "feature_distribution_pdf": str(feature_distribution_pdf),
        "quality_control_dir": str(output_dir),
        "state_counts_csv": str(state_counts_csv),
        "transition_matrix_csv": None if transition_csv is None else str(transition_csv),
        "model_selection_csv": None if selection_csv is None else str(selection_csv),
    }

def _plot_hmm_feature_distribution_pdf(
    adata,
    *,
    feature_cols,
    pdf_path,
    scaler_mean,
    scaler_scale,
    cluster_col=INTRINSIC_STATE_COL,
    title="Feature distributions",
    features_per_page=4,
    preprocessing_params=None,
    verbose=True,
):
    if adata is None or adata.n_obs == 0:
        return

    feature_cols = [] if feature_cols is None else list(feature_cols)
    valid_features = [str(c) for c in feature_cols if str(c) in adata.var_names]
    if len(valid_features) == 0:
        return

    X_scaled = _to_numpy_2d(adata[:, valid_features].X).astype(float, copy=False)
    X_raw = X_scaled
    raw_title = "Raw"
    preprocessing_meta = (
        dict(preprocessing_params)
        if isinstance(preprocessing_params, dict)
        else getattr(adata, "uns", {}).get("preprocessing", {})
    )
    mean = np.asarray(scaler_mean, dtype=float) if scaler_mean is not None else np.asarray([], dtype=float)
    scale = np.asarray(scaler_scale, dtype=float) if scaler_scale is not None else np.asarray([], dtype=float)
    if mean.shape[0] == len(valid_features) and scale.shape[0] == len(valid_features):
        safe_scale = scale.copy()
        safe_scale[safe_scale == 0] = 1.0
        X_raw = (X_scaled * safe_scale[None, :]) + mean[None, :]
    else:
        raw_title = "Raw (scaler mismatch; showing scaled)"
    X_raw = _invert_log_scaling_in_continuous_matrix(
        X_raw,
        preprocessing_params=preprocessing_meta,
        feature_cols=valid_features,
        inplace=True,
        strict=False,
    )

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    n_features = len(valid_features)
    if cluster_col in adata.obs.columns:
        cluster_series = pd.Series(adata.obs[cluster_col], index=adata.obs.index, dtype="string")
        cluster_order = sorted(cluster_series.dropna().unique().tolist(), key=_mixed_label_sort_key)
    else:
        cluster_series = pd.Series(["all"] * adata.n_obs, index=adata.obs.index)
        cluster_order = []
    page_count = n_features

    overall_started = _vstart(
        verbose,
        "state-hmm",
        (
            "feature distribution pdf sections "
            f"| features={n_features}, pages={page_count}, clusters={len(cluster_order)}"
        ),
    )
    with PdfPages(pdf_path) as pdf:
        pages_started = _vstart(
            verbose,
            "state-hmm",
            f"feature distribution feature pages | pages={page_count}, panels_per_page={1 + len(cluster_order)}",
        )
        cluster_values = cluster_series.to_numpy()
        for idx, feature_name in enumerate(valid_features):
            page_specs = [("overall", np.ones(adata.n_obs, dtype=bool), "#C94C4C")]
            page_specs.extend(
                (f"cluster {cluster_name}", cluster_values == str(cluster_name), "#4C7E9E")
                for cluster_name in cluster_order
            )
            page_data = []
            finite_vals = []
            for panel_label, page_mask, color in page_specs:
                raw_vals = pd.to_numeric(
                    pd.Series(X_raw[page_mask, idx]),
                    errors="coerce",
                ).dropna().to_numpy(dtype=float)
                page_data.append((panel_label, color, raw_vals))
                if raw_vals.size > 0:
                    finite_vals.append(raw_vals)
            y_limits = None
            if len(finite_vals) > 0:
                all_vals = np.concatenate(finite_vals)
                y_min = float(np.nanmin(all_vals))
                y_max = float(np.nanmax(all_vals))
                if np.isfinite(y_min) and np.isfinite(y_max):
                    if np.isclose(y_min, y_max):
                        pad = 0.5 if np.isclose(y_min, 0.0) else max(abs(y_min) * 0.1, 1e-6)
                    else:
                        pad = max((y_max - y_min) * 0.05, 1e-6)
                    y_limits = (y_min - pad, y_max + pad)
            fig = Figure(figsize=(max(6.4, 1.15 * len(page_data) + 3.6), 7.2))
            ax = fig.add_subplot(1, 1, 1)
            positions = np.arange(1, len(page_data) + 1, dtype=float)
            tick_labels = []
            for pos, (panel_label, color, raw_vals) in zip(positions, page_data):
                n_vals = int(raw_vals.size)
                tick_labels.append(f"{panel_label}\n(n={n_vals})")
                if raw_vals.size > 0:
                    vp = ax.violinplot(raw_vals, positions=[pos], showmeans=True, showextrema=False, widths=0.72)
                    for body in vp["bodies"]:
                        body.set_facecolor(color)
                        body.set_edgecolor(color)
                        body.set_alpha(0.62)
                    if "cmeans" in vp:
                        vp["cmeans"].set_color("#222222")
                        vp["cmeans"].set_linewidth(1.0)
                else:
                    ax.text(pos, 0.5, "No finite values", ha="center", va="center", transform=ax.get_xaxis_transform())
            ax.set_xlim(0.35, float(len(page_data)) + 0.65)
            ax.set_xticks(positions)
            ax.set_xticklabels(tick_labels, rotation=20, ha="right", fontsize=8)
            ax.set_ylabel(str(feature_name))
            ax.set_xlabel("overall and intrinsic cluster number")
            ax.grid(axis="y", alpha=0.2)
            if y_limits is not None:
                ax.set_ylim(*y_limits)
            ax.set_title(f"{raw_title} distributions", fontsize=10)

            fig.suptitle(
                f"{title} | {raw_title} | feature {idx + 1} of {n_features}: {feature_name}",
                fontsize=12,
                fontweight="bold",
                y=0.995,
            )
            fig.tight_layout(rect=[0, 0, 1, 0.98])
            _save_pdf_page_a4(pdf, fig, orientation="portrait")
            plt.close(fig)
        _vdone(verbose, "state-hmm", "feature distribution feature pages", pages_started)
    _vdone(verbose, "state-hmm", "feature distribution pdf sections", overall_started)


def _plot_hmm_state_diagnostics_pdf(
    adata,
    *,
    cluster_col,
    feature_cols,
    pdf_path,
    model=None,
    transition_df=None,
    sample_col="sample_name",
    selection_df=None,
    title="HMM state diagnostics",
    verbose=True,
):
    if adata is None or adata.n_obs == 0 or cluster_col not in adata.obs.columns:
        return

    ad = adata.copy()
    if transition_df is None:
        transition_df = _compute_hmm_transition_matrix_df(model)
    cluster_corr_df = _compute_cluster_correlation_df(ad, cluster_col, feature_cols)

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    feature_cols = [] if feature_cols is None else list(feature_cols)
    valid_features = [str(c) for c in feature_cols if str(c) in ad.var_names]
    overall_started = _vstart(
        verbose,
        "state-hmm",
        (
            "diagnostics pdf sections "
            f"| features={len(valid_features)}, clusters={int(ad.obs[cluster_col].astype(str).nunique())}"
        ),
    )
    with PdfPages(pdf_path) as pdf:
        if len(valid_features) > 0:
            heatmap_started = _vstart(
                verbose,
                "state-hmm",
                f"diagnostics Scanpy heatmap | features={len(valid_features)}",
            )
            ad_heat = ad[:, valid_features].copy()
            cluster_order = sorted(
                pd.Series(ad_heat.obs[cluster_col], index=ad_heat.obs.index, dtype="string").dropna().unique().tolist(),
                key=_mixed_label_sort_key,
            )
            ad_heat.obs[cluster_col] = pd.Categorical(
                pd.Series(ad_heat.obs[cluster_col], index=ad_heat.obs.index, dtype="string"),
                categories=cluster_order,
                ordered=True,
            )
            dendrogram_arg = False
            if len(cluster_order) > 1:
                try:
                    dendrogram_key = f"dendrogram_{cluster_col}"
                    if dendrogram_key not in ad_heat.uns:
                        sc.tl.dendrogram(ad_heat, groupby=cluster_col, key_added=dendrogram_key)
                    if dendrogram_key in ad_heat.uns:
                        dendrogram_arg = dendrogram_key
                except Exception:
                    dendrogram_arg = False

            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Starting a Matplotlib GUI outside of the main thread",
                )
                sc.pl.heatmap(
                    ad_heat,
                    var_names=valid_features,
                    groupby=cluster_col,
                    standard_scale="var",
                    figsize=A4_LANDSCAPE,
                    swap_axes=True,
                    dendrogram=dendrogram_arg,
                    show_gene_labels=True,
                    show=False,
                )
                fig_sc = plt.gcf()
            fig_sc.suptitle(f"{title} | Scanpy heatmap", y=0.995)
            _save_pdf_page_a4(pdf, fig_sc, orientation="landscape")
            plt.close(fig_sc)
            _vdone(verbose, "state-hmm", "diagnostics Scanpy heatmap", heatmap_started)

        if not transition_df.empty:
            transition_started = _vstart(verbose, "state-hmm", "diagnostics transition matrix heatmap")
            fig_dyn = Figure(figsize=A4_LANDSCAPE)
            ax = fig_dyn.add_subplot(1, 1, 1)
            sns.heatmap(
                transition_df,
                cmap="mako",
                ax=ax,
                cbar_kws={"label": "Transition probability"},
            )
            ax.set_title("Learned transition matrix", fontsize=10)
            ax.set_xlabel(transition_df.columns.name or "State")
            ax.set_ylabel(transition_df.index.name or "State")
            ax.tick_params(axis="x", rotation=35, labelsize=8)
            ax.tick_params(axis="y", labelsize=8)
            fig_dyn.suptitle(f"{title} | HMM transition matrix", fontsize=12, fontweight="bold", y=0.995)
            fig_dyn.tight_layout(rect=[0, 0, 1, 0.97])
            _save_pdf_page_a4(pdf, fig_dyn, orientation="landscape")
            plt.close(fig_dyn)
            _vdone(verbose, "state-hmm", "diagnostics transition matrix heatmap", transition_started)

        if not cluster_corr_df.empty:
            cluster_corr_started = _vstart(verbose, "state-hmm", "diagnostics cluster correlation heatmap")
            fig_corr = Figure(figsize=A4_LANDSCAPE)
            ax = fig_corr.add_subplot(1, 1, 1)
            sns.heatmap(
                cluster_corr_df,
                cmap="vlag",
                vmin=-1.0,
                vmax=1.0,
                center=0.0,
                ax=ax,
                cbar_kws={"label": "Feature-profile correlation"},
            )
            ax.set_title("Cluster correlation by mean feature profile", fontsize=10)
            ax.set_xlabel(cluster_corr_df.columns.name or "State")
            ax.set_ylabel(cluster_corr_df.index.name or "State")
            ax.tick_params(axis="x", rotation=35, labelsize=8)
            ax.tick_params(axis="y", labelsize=8)
            fig_corr.suptitle(f"{title} | Cluster correlation", fontsize=12, fontweight="bold", y=0.995)
            fig_corr.tight_layout(rect=[0, 0, 1, 0.97])
            _save_pdf_page_a4(pdf, fig_corr, orientation="landscape")
            plt.close(fig_corr)
            _vdone(verbose, "state-hmm", "diagnostics cluster correlation heatmap", cluster_corr_started)

        if isinstance(selection_df, pd.DataFrame) and (not selection_df.empty):
            selection_started = _vstart(
                verbose,
                "state-hmm",
                f"diagnostics model selection plot | rows={len(selection_df)}",
            )
            fig_sel = Figure(figsize=A4_LANDSCAPE)
            axes = fig_sel.subplots(1, 2)
            axes[0].plot(selection_df["k"], selection_df["bic"], marker="o", color="#D18951")
            axes[0].set_title("BIC by number of states", fontsize=10)
            axes[0].set_xlabel("k")
            axes[0].set_ylabel("BIC")
            axes[0].grid(alpha=0.25)

            axes[1].plot(selection_df["k"], selection_df["logL"], marker="o", color="#4C7E9E")
            axes[1].set_title("Log-likelihood by number of states", fontsize=10)
            axes[1].set_xlabel("k")
            axes[1].set_ylabel("logL")
            axes[1].grid(alpha=0.25)

            fig_sel.suptitle(f"{title} | Model selection", fontsize=12, fontweight="bold", y=0.995)
            fig_sel.tight_layout(rect=[0, 0, 1, 0.97])
            _save_pdf_page_a4(pdf, fig_sel, orientation="landscape")
            plt.close(fig_sel)
            _vdone(verbose, "state-hmm", "diagnostics model selection plot", selection_started)
    _vdone(verbose, "state-hmm", "diagnostics pdf sections", overall_started)


def run_hmm_state_clustering(
    features,
    binary_features_to_group,
    output_dir,
    cell_type,
    *,
    n_states=None,
    k_min=1,
    k_max=8,
    covariance_type="full",
    n_iter=200,
    tol=1e-3,
    sticky=False,
    stickiness_kappa=8.0,
    transmat_alpha=1.0,
    min_covar=1e-3,
    feature_smoothing_window=1,
    smoothing_min_periods=1,
    start_offset=0,
    start_offset_fill_mode="backfill",
    additional_window_features=None,
    window_features_window=5,
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    log_scale_features=None,
    random_state=123,
    df_positions=None,
    state_paths=None,
    return_details=False,
    verbose=True,
):
    hmm_started = _vstart(verbose, "state-hmm", "run HMM state clustering")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    state_clustering_outdir = _resolve_hmm_quality_control_outdir(state_paths=state_paths) / "raw"
    state_clustering_outdir.mkdir(parents=True, exist_ok=True)

    raw_features = _normalize_selection(features)
    window_features = _dedupe_preserve_order(additional_window_features)
    kept_features = _dedupe_preserve_order(list(raw_features) + list(window_features))
    binary_cols_to_merge = _normalize_selection(binary_features_to_group)
    log_scale_selectors = _normalize_log_scale_feature_selectors(log_scale_features)
    n_states = _validate_n_states(n_states)
    start_offset = max(0, int(start_offset))
    start_offset_fill_mode = str(start_offset_fill_mode).strip().lower()
    if start_offset_fill_mode not in {"backfill", "leave_unassigned"}:
        raise ValueError(
            "start_offset_fill_mode must be 'backfill' or 'leave_unassigned', "
            f"got '{start_offset_fill_mode}'."
        )
    sort_cols = ["sample_name", "TrackID", "position_t"]
    required_cols = sort_cols + raw_features + binary_cols_to_merge
    if len(window_features) > 0:
        required_cols.extend(["position_x", "position_y", "position_z"])
    required_cols = _dedupe_preserve_order(required_cols)

    if df_positions is None:
        positions_csv_path = _resolve_positions_csv_path(state_paths=state_paths)
        load_started = _vstart(verbose, "state-hmm", "load track-features CSV")
        df_positions = pd.read_csv(positions_csv_path, low_memory=False)
        _vdone(verbose, "state-hmm", "load track-features CSV", load_started)
    else:
        positions_csv_path = None

    _require_columns(df_positions, required_cols, "HMM state clustering")
    df_base = df_positions.sort_values(sort_cols, kind="mergesort").copy()
    df_base["_start_offset_track_idx"] = (
        df_base.groupby(["sample_name", "TrackID"], sort=False, observed=False).cumcount().astype(int)
    )
    window_feature_window_size = 1
    if len(window_features) > 0:
        df_window, window_feature_window_size = _compute_optional_window_feature_frame(
            df_base,
            window_features=window_features,
            time_col="position_t",
            id_cols=("sample_name", "TrackID"),
            window_size=window_features_window,
            verbose=verbose,
        )
        df_base = df_base.merge(
            df_window,
            on=sort_cols,
            how="left",
            validate="one_to_one",
        )
    df_base = _smooth_timepoint_features(
        df_base,
        raw_features,
        id_cols=("sample_name", "TrackID"),
        time_col="position_t",
        window=feature_smoothing_window,
        min_periods=smoothing_min_periods,
    )

    cap_started = _vstart(verbose, "state-hmm", "cap observation features to quantiles")
    df_capped, cap_limits = cap_values_to_quantile(
        df_base,
        kept_features,
        lower_quantile=lower_quantile_cap,
        upper_quantile=upper_quantile_cap,
        return_limits=True,
    )
    _vdone(verbose, "state-hmm", "cap observation features to quantiles", cap_started)

    numeric_frame = df_capped[kept_features].apply(pd.to_numeric, errors="coerce")
    valid_mask = numeric_frame.notna().all(axis=1)
    n_dropped = int((~valid_mask).sum())
    if n_dropped > 0:
        _vinfo(
            verbose,
            "state-hmm",
            f"Dropping {n_dropped} timepoints with NaNs in the HMM observation matrix.",
        )
    if not bool(valid_mask.any()):
        raise ValueError("No valid timepoints remain after smoothing/quantile capping for HMM state clustering.")

    df_valid = df_capped.loc[valid_mask].copy()
    df_valid.loc[:, kept_features] = numeric_frame.loc[valid_mask, kept_features].to_numpy(dtype=float)
    df_valid["_start_offset_scored"] = df_valid["_start_offset_track_idx"] >= int(start_offset)
    df_valid["_start_offset_skipped"] = ~df_valid["_start_offset_scored"]
    df_scored = df_valid[df_valid["_start_offset_scored"]].copy()
    if len(df_scored) == 0:
        raise ValueError(
            "No valid HMM timepoints remain after applying start_offset="
            f"{int(start_offset)}."
        )

    log_scaling_params = _resolve_log_scale_feature_cols(
        kept_features,
        log_scale_selectors,
    )
    if len(log_scaling_params["unresolved_features"]) > 0:
        raise ValueError(
            "log_scale_features could not be resolved to HMM observation feature columns: "
            f"{log_scaling_params['unresolved_features'][:20]}"
        )
    if len(log_scaling_params["resolved_feature_cols"]) > 0:
        X_log = _apply_log1p_to_feature_matrix(
            df_valid[kept_features].to_numpy(dtype=float, copy=True),
            feature_cols=kept_features,
            resolved_feature_cols=log_scaling_params["resolved_feature_cols"],
            inplace=True,
        )
        df_valid.loc[:, kept_features] = X_log
        df_scored.loc[:, kept_features] = df_valid.loc[df_scored.index, kept_features].to_numpy(dtype=float)
        del X_log

    scale_started = _vstart(verbose, "state-hmm", "z-scale observation features")
    scaler = StandardScaler()
    df_valid.loc[:, kept_features] = scaler.fit_transform(df_valid[kept_features].to_numpy(dtype=float))
    df_scored.loc[:, kept_features] = df_valid.loc[df_scored.index, kept_features].to_numpy(dtype=float)
    _vdone(verbose, "state-hmm", "z-scale observation features", scale_started)

    hmm_helper = run_sticky_hmm_state_classification if bool(sticky) else run_hmm_state_classification
    fit_started = _vstart(
        verbose,
        "state-hmm",
        f"fit {'sticky ' if bool(sticky) else ''}HMM | n_states={n_states}",
    )
    df_states, model, selection_df = hmm_helper(
        df_features=df_scored,
        feature_cols=kept_features,
        id_cols=("sample_name", "TrackID"),
        time_col="position_t",
        n_states=n_states,
        k_min=int(k_min),
        k_max=int(k_max),
        covariance_type=str(covariance_type),
        n_iter=int(n_iter),
        tol=float(tol),
        random_state=int(random_state),
        verbose=bool(verbose),
        error_on_nans=True,
        out_col_name=HMM_INTRINSIC_RAW_STATE_COL,
        stickiness_kappa=float(stickiness_kappa),
        transmat_alpha=float(transmat_alpha),
        min_covar=float(min_covar),
    )
    _vdone(verbose, "state-hmm", "fit HMM", fit_started)

    df_model = df_valid.copy()
    df_model[HMM_INTRINSIC_RAW_STATE_COL] = pd.Series(pd.NA, index=df_model.index, dtype="Int64")
    df_model[INTRINSIC_STATE_COL] = pd.Series(pd.NA, index=df_model.index, dtype="string")
    scored_state_map = df_states.set_index(sort_cols)[HMM_INTRINSIC_RAW_STATE_COL]
    scored_index = pd.MultiIndex.from_frame(df_model[sort_cols])
    raw_states = pd.Series(
        pd.to_numeric(
            scored_state_map.reindex(scored_index).to_numpy(),
            errors="coerce",
        ),
        index=df_model.index,
        dtype="Int64",
    )
    df_model.loc[:, HMM_INTRINSIC_RAW_STATE_COL] = raw_states
    df_model.loc[:, INTRINSIC_STATE_COL] = raw_states.astype("string")
    if start_offset_fill_mode == "backfill" and start_offset > 0:
        for _, group_idx in df_model.groupby(["sample_name", "TrackID"], sort=False, observed=False).groups.items():
            group = df_model.loc[list(group_idx)]
            scored_group = group[group["_start_offset_scored"]]
            skipped_group = group[group["_start_offset_skipped"]]
            if len(scored_group) == 0 or len(skipped_group) == 0:
                continue
            first_raw = scored_group[HMM_INTRINSIC_RAW_STATE_COL].iloc[0]
            if pd.isna(first_raw):
                continue
            df_model.loc[skipped_group.index, HMM_INTRINSIC_RAW_STATE_COL] = first_raw
            first_intrinsic = _format_hmm_raw_state_value_for_key(
                first_raw,
                label=HMM_INTRINSIC_RAW_STATE_COL,
            )
            if first_intrinsic != "":
                df_model.loc[skipped_group.index, INTRINSIC_STATE_COL] = first_intrinsic

    obs_cols = sort_cols + binary_cols_to_merge + [HMM_INTRINSIC_RAW_STATE_COL, INTRINSIC_STATE_COL]
    model_adata = df_to_adata(df_model, feature_cols=kept_features, obs_cols=obs_cols)
    raw_labels = _coerce_hmm_raw_state_series(
        model_adata.obs[HMM_INTRINSIC_RAW_STATE_COL],
        label=HMM_INTRINSIC_RAW_STATE_COL,
    )
    intrinsic_labels = pd.Series(
        model_adata.obs[INTRINSIC_STATE_COL],
        index=model_adata.obs.index,
        dtype="string",
    )
    model_adata.obs[HMM_INTRINSIC_RAW_STATE_COL] = raw_labels
    model_adata.obs[INTRINSIC_STATE_COL] = pd.Categorical(intrinsic_labels)

    binary_group_constraints = _infer_binary_group_constraints(df_model, binary_cols_to_merge)
    _rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col=INTRINSIC_STATE_COL,
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=True,
    )
    if start_offset_fill_mode == "leave_unassigned":
        unassigned_mask = pd.Series(
            model_adata.obs[INTRINSIC_STATE_COL],
            index=model_adata.obs.index,
            dtype="string",
        ).isna()
        if bool(unassigned_mask.any()):
            for col in [
                INTRINSIC_STATE_COL,
                "behavioral_clusterid",
                "full_behavioral_cluster",
                FULL_STATE_COL,
            ]:
                if col in model_adata.obs.columns:
                    model_adata.obs[col] = pd.Categorical(
                        pd.Series(model_adata.obs[col], index=model_adata.obs.index, dtype="string").where(
                            ~unassigned_mask,
                            pd.NA,
                        )
                    )

    state_paths.model_adata_path.parent.mkdir(parents=True, exist_ok=True)
    early_write_started = _vstart(verbose, "state-hmm", "write preliminary model adata")
    model_adata.write(state_paths.model_adata_path, compression="gzip")
    _vdone(verbose, "state-hmm", "write preliminary model adata", early_write_started)
    _vsave(verbose, "state-hmm", "preliminary model adata", state_paths.model_adata_path)

    diagnostics_pdf = state_clustering_outdir / "behavioral_clustering_diagnostics.pdf"
    diagnostics_started = _vstart(verbose, "state-hmm", "write HMM diagnostics PDF")
    _plot_hmm_state_diagnostics_pdf(
        model_adata,
        cluster_col=INTRINSIC_STATE_COL,
        feature_cols=kept_features,
        pdf_path=diagnostics_pdf,
        model=model,
        selection_df=selection_df,
        title=f"all_data | hmm (states={int(model.n_components)})",
        verbose=verbose,
    )
    _vdone(verbose, "state-hmm", "write HMM diagnostics PDF", diagnostics_started)
    _vsave(verbose, "state-hmm", "diagnostics pdf", diagnostics_pdf)

    diagnostics_csvs_started = _vstart(verbose, "state-hmm", "export HMM diagnostics CSVs")
    diagnostics_csvs = {}
    state_means = _compute_state_means(model_adata, INTRINSIC_STATE_COL, kept_features)
    if not state_means.empty:
        state_means_started = _vstart(verbose, "state-hmm", "write HMM state means CSV")
        state_means_csv = state_clustering_outdir / "behavioral_clustering_hmm_state_means.csv"
        state_means.to_csv(state_means_csv)
        diagnostics_csvs["state_means_csv"] = str(state_means_csv)
        _vdone(verbose, "state-hmm", "write HMM state means CSV", state_means_started)
        _vsave(verbose, "state-hmm", "state means csv", state_means_csv)

    cluster_corr_df = _compute_cluster_correlation_df(
        model_adata,
        INTRINSIC_STATE_COL,
        kept_features,
    )
    if not cluster_corr_df.empty:
        cluster_corr_started = _vstart(verbose, "state-hmm", "write HMM cluster correlation CSV")
        cluster_corr_csv = state_clustering_outdir / "behavioral_clustering_hmm_cluster_correlation.csv"
        cluster_corr_df.to_csv(cluster_corr_csv)
        diagnostics_csvs["cluster_correlation_csv"] = str(cluster_corr_csv)
        _vdone(verbose, "state-hmm", "write HMM cluster correlation CSV", cluster_corr_started)
        _vsave(verbose, "state-hmm", "cluster correlation csv", cluster_corr_csv)

    transition_csv = None
    transition_df = _compute_hmm_transition_matrix_df(model)
    if not transition_df.empty:
        transition_csv_started = _vstart(verbose, "state-hmm", "write HMM transition matrix CSV")
        transition_csv = state_clustering_outdir / "behavioral_clustering_hmm_transition_matrix.csv"
        transition_df.to_csv(transition_csv)
        diagnostics_csvs["transition_matrix_csv"] = str(transition_csv)
        _vdone(verbose, "state-hmm", "write HMM transition matrix CSV", transition_csv_started)
        _vsave(verbose, "state-hmm", "transition matrix csv", transition_csv)
    _vdone(verbose, "state-hmm", "export HMM diagnostics CSVs", diagnostics_csvs_started)

    state_counts_started = _vstart(verbose, "state-hmm", "write HMM state counts CSV")
    state_counts_csv = state_clustering_outdir / "behavioral_clustering_hmm_state_counts.csv"
    (
        pd.Series(model_adata.obs[INTRINSIC_STATE_COL], dtype="string")
        .dropna()
        .value_counts()
        .sort_index(key=lambda idx: [_mixed_label_sort_key(v) for v in idx])
        .rename_axis(INTRINSIC_STATE_COL)
        .reset_index(name="count")
        .to_csv(state_counts_csv, index=False)
    )
    _vdone(verbose, "state-hmm", "write HMM state counts CSV", state_counts_started)
    _vsave(verbose, "state-hmm", "state counts csv", state_counts_csv)

    selection_csv = None
    if isinstance(selection_df, pd.DataFrame) and (not selection_df.empty):
        selection_csv_started = _vstart(verbose, "state-hmm", "write HMM model selection CSV")
        selection_csv = state_clustering_outdir / "behavioral_clustering_hmm_model_selection.csv"
        selection_df.to_csv(selection_csv, index=False)
        diagnostics_csvs["model_selection_csv"] = str(selection_csv)
        _vdone(verbose, "state-hmm", "write HMM model selection CSV", selection_csv_started)
        _vsave(verbose, "state-hmm", "HMM model selection csv", selection_csv)

    preprocessing_meta = {
        "observation_mode": "timepoint_hmm",
        "features": list(kept_features),
        "raw_feature_cols": list(raw_features),
        "kept_features": list(kept_features),
        "continuous_feature_cols": list(kept_features),
        "additional_window_features": list(window_features),
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "non_feature_cols": list(sort_cols + binary_cols_to_merge),
        "windowing": {
            "observation_mode": "timepoint_hmm",
            "features": list(raw_features),
            "binary_features_to_group": list(binary_cols_to_merge),
            "additional_window_features": list(window_features),
            "window_size": 1,
            "window_feature_window_size": int(window_feature_window_size),
            "min_spacing": None,
            "descriptive_features": ["timepoint_hmm"],
            "incomplete_window_policy": "timepoint",
            "lower_quantile_cap": None if lower_quantile_cap is None else float(lower_quantile_cap),
            "upper_quantile_cap": None if upper_quantile_cap is None else float(upper_quantile_cap),
            "feature_smoothing_window": int(feature_smoothing_window),
            "smoothing_min_periods": int(smoothing_min_periods),
            "start_offset": int(start_offset),
            "start_offset_fill_mode": str(start_offset_fill_mode),
        },
        "timepoint_smoothing": {
            "feature_smoothing_window": int(feature_smoothing_window),
            "smoothing_min_periods": int(smoothing_min_periods),
            "start_offset": int(start_offset),
            "start_offset_fill_mode": str(start_offset_fill_mode),
        },
        "quantile_capping": {
            "lower_quantile": None if lower_quantile_cap is None else float(lower_quantile_cap),
            "upper_quantile": None if upper_quantile_cap is None else float(upper_quantile_cap),
            "feature_limits": dict(cap_limits),
        },
        "log_scaling": {
            "transform": "log1p",
            "requested_features": list(log_scaling_params.get("requested_features", [])),
            "resolved_feature_cols": list(log_scaling_params.get("resolved_feature_cols", [])),
        },
        "scaler": {
            "mean": scaler.mean_.tolist(),
            "scale": scaler.scale_.tolist(),
        },
        "positions_csv_path": None if positions_csv_path is None else str(positions_csv_path),
        "n_timepoints_input": int(len(df_positions)),
        "n_timepoints_kept": int(model_adata.n_obs),
        "n_timepoints_dropped_nan": int(n_dropped),
        "start_offset": int(start_offset),
        "start_offset_fill_mode": str(start_offset_fill_mode),
        "n_timepoints_skipped_start_offset": int(df_valid["_start_offset_skipped"].sum()),
        "n_timepoints_scored": int(df_valid["_start_offset_scored"].sum()),
        "n_tracks_too_short_after_offset": int(
            (
                df_base.groupby(["sample_name", "TrackID"], sort=False, observed=False).size()
                <= int(start_offset)
            ).sum()
        ),
    }
    feature_distribution_pdf = state_clustering_outdir / "behavioral_clustering_feature_distributions.pdf"
    feature_distribution_started = _vstart(verbose, "state-hmm", "write HMM feature distribution PDF")
    _plot_hmm_feature_distribution_pdf(
        model_adata,
        feature_cols=kept_features,
        pdf_path=feature_distribution_pdf,
        scaler_mean=scaler.mean_.tolist(),
        scaler_scale=scaler.scale_.tolist(),
        title=f"all_data | hmm (states={int(model.n_components)})",
        preprocessing_params=preprocessing_meta,
        verbose=verbose,
    )
    _vdone(verbose, "state-hmm", "write HMM feature distribution PDF", feature_distribution_started)
    _vsave(verbose, "state-hmm", "feature distribution pdf", feature_distribution_pdf)
    clustering_meta = {
        "clustering_method": "hmm",
        "use_pca": False,
        "use_rep": "X",
        "random_state": int(random_state),
        "n_neighbors": int(min(30, max(2, model_adata.n_obs - 1))) if int(model_adata.n_obs) > 1 else 1,
        "resolution": None,
        "non_feature_cols": list(sort_cols + binary_cols_to_merge),
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "binary_group_constraints": copy.deepcopy(dict(binary_group_constraints)),
        "hmm": {
            "n_states_requested": n_states,
            "auto_selected": isinstance(n_states, str) and (str(n_states).lower() == "auto"),
            "selected_k": int(model.n_components),
            "covariance_type": str(covariance_type),
            "sticky": bool(sticky),
            "n_iter": int(n_iter),
            "tol": float(tol),
            "stickiness_kappa": float(stickiness_kappa),
            "transmat_alpha": float(transmat_alpha),
            "min_covar": float(min_covar),
            "k_min": int(k_min),
            "k_max": int(k_max),
            "start_offset": int(start_offset),
            "start_offset_fill_mode": str(start_offset_fill_mode),
        },
        "diagnostics_pdf": str(diagnostics_pdf),
        "diagnostics_csvs": dict(diagnostics_csvs),
        "feature_distribution_pdf": str(feature_distribution_pdf),
        "quality_control_dir": str(state_clustering_outdir.parent),
        "raw_quality_control_dir": str(state_clustering_outdir),
        "state_counts_csv": str(state_counts_csv),
        "transition_matrix_csv": None if transition_csv is None else str(transition_csv),
        "sample_state_occupancy_csv": None,
        "model_selection_csv": None if selection_csv is None else str(selection_csv),
        "state_composition_report_pdf": None,
        "state_composition_report_auc_csv": None,
        "state_composition_report_plot_csvs": [],
        "state_composition_report_error": None,
        "state_transition_report_dir": None,
        "state_transition_matrix_counts_csv": None,
        "state_transition_matrix_probs_csv": None,
        "state_transition_report_error": None,
        "state_reports_ready": False,
        "state_reports_reason": (
            "State reports have not been created yet. Use 'Create analysis plots' in the HMM widget after curated renaming."
        ),
        "exemplar_videos": {
            "enabled": False,
            "reason": "HMM stage-1 clustering uses timepoint observations and does not render window exemplars.",
        },
    }
    classification_meta = {
        "intrinsic_output_col": INTRINSIC_STATE_COL,
        "full_output_col": FULL_STATE_COL,
        "intrinsic_intrinsic_state_col": INTRINSIC_STATE_COL,
        "full_label_strategy": "binary_group_plus_intrinsic",
        "hmm_deployment_artifact_path": str(
            _resolve_hmm_deployment_artifact_path(state_paths=state_paths, cell_type=cell_type)
        ),
    }
    metadata_started = _vstart(verbose, "state-hmm", "attach HMM preprocessing/clustering metadata")
    model_adata.uns["preprocessing"] = preprocessing_meta
    model_adata.uns["clustering"] = clustering_meta
    model_adata.uns["classification"] = classification_meta
    _current_hmm_full_state_colors(model_adata, state_col=FULL_STATE_COL)
    _vdone(verbose, "state-hmm", "attach HMM preprocessing/clustering metadata", metadata_started)

    write_started = _vstart(verbose, "state-hmm", "write model adata")
    state_paths.model_adata_path.parent.mkdir(parents=True, exist_ok=True)
    model_adata.write(state_paths.model_adata_path, compression="gzip")
    _vdone(verbose, "state-hmm", "write model adata", write_started)
    _vsave(verbose, "state-hmm", "model adata", state_paths.model_adata_path)
    _vinfo(
        verbose,
        "state-hmm",
        (
            "HMM summary: "
            f"timepoints={model_adata.n_obs}, "
            f"features={model_adata.n_vars}, "
            f"states={int(model.n_components)}, "
            f"binary_groups={pd.Series(model_adata.obs['binary_group'], dtype='string').dropna().nunique()}"
        ),
    )
    _vdone(verbose, "state-hmm", "run HMM state clustering", hmm_started)
    if bool(return_details):
        return {
            "model_adata": model_adata,
            "hmm_model": model,
            "selection_df": selection_df,
            "state_paths": state_paths,
            "hmm_deployment_artifact": _build_hmm_deployment_artifact(
                model_adata=model_adata,
                hmm_model=model,
                state_paths=state_paths,
            ),
        }
    return model_adata


def _finalize_hmm_apply_outputs(
    adata_full,
    *,
    state_paths,
    preprocessing_meta,
    classification_meta,
    full_output_col,
    export_state_composition=False,
    export_state_transitions=False,
    state_transitions_rows_per_page=3,
    state_transitions_include_self_pairs=True,
    state_transitions_min_count=100,
    state_transitions_relative_count=0.0,
    verbose=True,
):
    report_adata = adata_full
    if full_output_col in adata_full.obs.columns:
        report_mask = pd.Series(adata_full.obs[full_output_col], index=adata_full.obs.index, dtype="string").notna()
        if not bool(report_mask.any()):
            report_adata = adata_full[:0].copy()
        elif not bool(report_mask.all()):
            report_adata = adata_full[report_mask.to_numpy(dtype=bool)].copy()

    adata_full.uns["preprocessing"] = dict(preprocessing_meta)
    adata_full.uns["classification"] = dict(classification_meta)
    report_state_colors = {}
    if full_output_col in report_adata.obs.columns:
        report_state_colors = _normalize_label_color_map(
            _obs_label_values(report_adata, full_output_col),
            colors=(
                _extract_state_color_mapping(classification_meta.get("state_colors", {}), full_output_col)
                or _extract_state_color_mapping(classification_meta.get("state_colors", {}), FULL_STATE_COL)
            ),
        )

    state_composition_report_pdf = None
    state_composition_report_pdfs = []
    state_composition_report_auc_csv = None
    state_composition_report_plot_csvs = []
    state_composition_report_error = None

    if bool(export_state_composition) and (full_output_col in report_adata.obs.columns) and int(report_adata.n_obs) > 0:
        report_dir = state_paths.state_composition_outdir
        report_dir.mkdir(parents=True, exist_ok=True)
        report_token = _sanitize_filename_token(
            full_output_col,
            fallback=FULL_STATE_COL,
        )
        report_pdf_path = report_dir / f"state_composition_report_{report_token}.pdf"
        report_auc_csv_path = report_dir / f"state_composition_report_{report_token}.csv"
        try:
            report_out = save_state_composition_report(
                adata=report_adata,
                output_pdf_path=report_pdf_path,
                output_csv_path=report_auc_csv_path,
                time_col="position_t",
                state_col=str(full_output_col),
                sample_col="sample_name",
                include_pooled_summary=True,
                state_colors=report_state_colors,
                verbose=verbose,
            )
            pdf_paths = report_out.get("pdf_paths", {})
            if isinstance(pdf_paths, dict) and len(pdf_paths) > 0:
                state_composition_report_pdfs = [str(v) for v in pdf_paths.values()]
                state_composition_report_pdf = str(
                    pdf_paths.get("stacked_by_sample", state_composition_report_pdfs[0])
                )
            else:
                state_composition_report_pdf = str(report_out.get("pdf_path", report_pdf_path))
                state_composition_report_pdfs = [state_composition_report_pdf]
            state_composition_report_auc_csv = str(report_out.get("csv_path", report_auc_csv_path))
            plot_csv_paths = report_out.get("plot_data_csv_paths", {})
            if isinstance(plot_csv_paths, dict) and len(plot_csv_paths) > 0:
                state_composition_report_plot_csvs = [str(v) for v in plot_csv_paths.values()]
            else:
                plot_csv_path = report_out.get("plot_data_csv_path", None)
                if plot_csv_path is not None:
                    state_composition_report_plot_csvs = [str(plot_csv_path)]
        except Exception as exc:
            state_composition_report_error = str(exc)
            warnings.warn(
                "Could not generate state composition report after HMM apply: "
                f"{exc}",
                RuntimeWarning,
            )
            _vinfo(verbose, "state-hmm-apply", f"State composition report skipped: {exc}")

    adata_full.uns["classification"]["state_composition_report_pdf"] = state_composition_report_pdf
    adata_full.uns["classification"]["state_composition_report_pdfs"] = list(state_composition_report_pdfs)
    adata_full.uns["classification"]["state_composition_report_auc_csv"] = state_composition_report_auc_csv
    adata_full.uns["classification"]["state_composition_report_plot_csvs"] = list(
        state_composition_report_plot_csvs
    )
    adata_full.uns["classification"]["state_composition_report_error"] = state_composition_report_error

    state_transition_report_dir = None
    state_transition_matrix_counts_csv = None
    state_transition_matrix_probs_csv = None
    state_transition_matrix_counts_no_self_csv = None
    state_transition_matrix_probs_no_self_csv = None
    state_transition_ngrams_all_csv = None
    state_transition_ngrams_pooled_csv = None
    state_transition_ngrams_per_end_csv = None
    state_transition_ngrams_per_start_csv = None
    state_transition_sankey_pdf = None
    state_transition_sankey_html_dir = None
    state_transition_report_error = None

    if bool(export_state_transitions) and (full_output_col in report_adata.obs.columns) and int(report_adata.n_obs) > 0:
        try:
            transitions_outdir = state_paths.state_transitions_outdir
            transitions_outdir.mkdir(parents=True, exist_ok=True)
            transition_out = save_state_transition_report(
                adata=report_adata,
                output_dir=transitions_outdir,
                state_col=str(full_output_col),
                id_cols=("sample_name", "TrackID"),
                time_col="position_t",
                include_self_pairs=bool(state_transitions_include_self_pairs),
                rows_per_page=max(1, int(state_transitions_rows_per_page)),
                sankey_min_count=int(state_transitions_min_count),
                sankey_relative_count=float(state_transitions_relative_count),
                state_colors=report_state_colors,
                verbose=verbose,
            )
            state_transition_report_dir = str(transition_out.get("output_dir", transitions_outdir))
            state_transition_matrix_counts_csv = transition_out.get("transition_matrix_counts_csv", None)
            state_transition_matrix_probs_csv = transition_out.get("transition_matrix_probs_csv", None)
            state_transition_matrix_counts_no_self_csv = transition_out.get(
                "transition_matrix_counts_no_self_csv", None
            )
            state_transition_matrix_probs_no_self_csv = transition_out.get(
                "transition_matrix_probs_no_self_csv", None
            )
            state_transition_ngrams_all_csv = transition_out.get("transition_ngrams_all_csv", None)
            state_transition_ngrams_pooled_csv = transition_out.get("transition_ngrams_pooled_csv", None)
            state_transition_ngrams_per_end_csv = transition_out.get("transition_ngrams_per_end_csv", None)
            state_transition_ngrams_per_start_csv = transition_out.get("transition_ngrams_per_start_csv", None)
            state_transition_sankey_pdf = transition_out.get("sankey_pdf_path", None)
            state_transition_sankey_html_dir = transition_out.get("sankey_html_dir", None)
        except Exception as exc:
            state_transition_report_error = str(exc)
            warnings.warn(
                "Could not generate state transition report after HMM apply: "
                f"{exc}",
                RuntimeWarning,
            )
            _vinfo(verbose, "state-hmm-apply", f"State transition report skipped: {exc}")

    adata_full.uns["classification"]["state_transition_report_dir"] = state_transition_report_dir
    adata_full.uns["classification"]["state_transition_matrix_counts_csv"] = state_transition_matrix_counts_csv
    adata_full.uns["classification"]["state_transition_matrix_probs_csv"] = state_transition_matrix_probs_csv
    adata_full.uns["classification"]["state_transition_matrix_counts_no_self_csv"] = (
        state_transition_matrix_counts_no_self_csv
    )
    adata_full.uns["classification"]["state_transition_matrix_probs_no_self_csv"] = (
        state_transition_matrix_probs_no_self_csv
    )
    adata_full.uns["classification"]["state_transition_ngrams_all_csv"] = state_transition_ngrams_all_csv
    adata_full.uns["classification"]["state_transition_ngrams_pooled_csv"] = state_transition_ngrams_pooled_csv
    adata_full.uns["classification"]["state_transition_ngrams_per_end_csv"] = state_transition_ngrams_per_end_csv
    adata_full.uns["classification"]["state_transition_ngrams_per_start_csv"] = (
        state_transition_ngrams_per_start_csv
    )
    adata_full.uns["classification"]["state_transition_sankey_pdf"] = state_transition_sankey_pdf
    adata_full.uns["classification"]["state_transition_sankey_html_dir"] = state_transition_sankey_html_dir
    adata_full.uns["classification"]["state_transition_report_error"] = state_transition_report_error

    adata_full.obs = adata_full.obs.drop(
        columns=["_start_offset_track_idx", "_start_offset_scored", "_start_offset_skipped"],
        errors="ignore",
    )

    write_started = _vstart(verbose, "state-hmm-apply", "write full output adata")
    adata_full.write(state_paths.full_output_adata_path, compression="gzip")
    _vdone(verbose, "state-hmm-apply", "write full output adata", write_started)
    _vsave(verbose, "state-hmm-apply", "full output adata", state_paths.full_output_adata_path)
    return adata_full


def apply_hmm_deployment_artifact_to_full_dataset(
    output_dir,
    cell_type,
    *,
    hmm_deployment_artifact,
    continuous_output_col=INTRINSIC_STATE_COL,
    full_output_col=FULL_STATE_COL,
    export_state_composition=False,
    export_state_transitions=False,
    state_transitions_rows_per_page=3,
    state_transitions_include_self_pairs=True,
    state_transitions_min_count=100,
    state_transitions_relative_count=0.0,
    start_offset=None,
    start_offset_fill_mode=None,
    state_paths=None,
    verbose=True,
):
    apply_started = _vstart(verbose, "state-hmm-apply", "apply HMM deployment artifact")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    artifact, artifact_source_type, artifact_source_path = _coerce_hmm_deployment_artifact_input(
        hmm_deployment_artifact,
        artifact_name="hmm_deployment_artifact",
    )

    pipeline_meta = copy.deepcopy(dict(artifact.get("pipeline_metadata", {})))
    pre_meta = dict(pipeline_meta.get("preprocessing", {}))
    class_meta = dict(pipeline_meta.get("classification", {}))
    binary_cols_to_merge = [str(c) for c in list(artifact.get("binary_cols_to_merge", []))]
    binary_group_constraints = copy.deepcopy(artifact.get("binary_group_constraints", None))
    enforce_binary_group_constraints = (
        isinstance(binary_group_constraints, dict)
        and ("forbidden_binary_combinations" in binary_group_constraints)
    )

    positions_csv_path = _resolve_positions_csv_path(state_paths=state_paths)
    load_started = _vstart(verbose, "state-hmm-apply", "load positions/features CSV")
    df_positions = pd.read_csv(positions_csv_path, low_memory=False)
    df_positions = df_positions.sort_values(by=["sample_name", "TrackID", "position_t"])
    _vdone(verbose, "state-hmm-apply", "load positions/features CSV", load_started)

    prepared = _prepare_hmm_apply_adata_from_df_positions(
        df_positions,
        preprocessing_meta=pre_meta,
        binary_cols_to_merge=binary_cols_to_merge,
        feature_cols=list(artifact["observation_feature_cols"]),
        start_offset=start_offset,
        start_offset_fill_mode=start_offset_fill_mode,
        verbose=verbose,
    )
    adata_full = prepared["adata"]
    df_scored = prepared["df_scored"]
    pre_meta = {
        **dict(pre_meta),
        **dict(prepared.get("drop_summary", {})),
        **dict(prepared.get("start_offset_summary", {})),
    }
    preprocess_started = _vstart(verbose, "state-hmm-apply", "rerun HMM preprocessing on observation matrix")
    X_query_all = _apply_hmm_saved_preprocessing_to_matrix(
        _to_numpy_2d(adata_full.X).astype(float, copy=True),
        preprocessing_meta=pre_meta,
        feature_cols=list(artifact["observation_feature_cols"]),
    )
    _vdone(verbose, "state-hmm-apply", "rerun HMM preprocessing on observation matrix", preprocess_started)
    score_mask = adata_full.obs["_start_offset_scored"].to_numpy(dtype=bool)
    X_query = X_query_all[score_mask, :]

    model = artifact["model"]
    adata_full.obs[HMM_INTRINSIC_RAW_STATE_COL] = pd.Series(
        pd.NA,
        index=adata_full.obs.index,
        dtype="Int64",
    )
    adata_full.obs[INTRINSIC_STATE_COL] = pd.Series(
        pd.NA,
        index=adata_full.obs.index,
        dtype="string",
    )
    adata_full.obs["intrinsic_behavioral_cluster_confidence"] = pd.Series(
        np.nan,
        index=adata_full.obs.index,
        dtype=float,
    )
    _validate_hmm_intrinsic_label_mapping_keys(
        artifact["intrinsic_label_mapping"],
        "hmm_deployment_artifact",
    )
    label_mapping = {int(k): str(v) for k, v in artifact["intrinsic_label_mapping"].items()}
    if bool(np.any(score_mask)):
        lengths = _compute_hmm_sequence_lengths(df_scored, id_cols=("sample_name", "TrackID"))
        raw_state_codes = np.asarray(model.predict(X_query, lengths=lengths), dtype=np.int64) + 1
        state_confidence = np.asarray(model.predict_proba(X_query, lengths=lengths), dtype=float).max(axis=1)
        raw_state_scored = pd.Series(
            raw_state_codes,
            index=adata_full.obs.index[score_mask],
            dtype="Int64",
        )
        intrinsic_scored = raw_state_scored.map(label_mapping)
        if intrinsic_scored.isna().any():
            missing_ids = sorted(set(raw_state_scored[intrinsic_scored.isna()].tolist()))
            raise ValueError(
                "HMM deployment artifact is missing curated intrinsic labels for predicted raw state id(s): "
                f"{missing_ids[:20]}"
            )
        adata_full.obs.loc[adata_full.obs.index[score_mask], HMM_INTRINSIC_RAW_STATE_COL] = raw_state_scored
        adata_full.obs.loc[adata_full.obs.index[score_mask], INTRINSIC_STATE_COL] = intrinsic_scored
        adata_full.obs.loc[adata_full.obs.index[score_mask], "intrinsic_behavioral_cluster_confidence"] = state_confidence

        if str(prepared["resolved_start_offset_fill_mode"]) == "backfill" and int(prepared["resolved_start_offset"]) > 0:
            for _, group_idx in adata_full.obs.groupby(["sample_name", "TrackID"], sort=False, observed=False).groups.items():
                group = adata_full.obs.loc[list(group_idx)]
                scored_group = group[group["_start_offset_scored"]]
                skipped_group = group[group["_start_offset_skipped"]]
                if len(scored_group) == 0 or len(skipped_group) == 0:
                    continue
                first_raw = scored_group[HMM_INTRINSIC_RAW_STATE_COL].iloc[0]
                first_intrinsic = scored_group[INTRINSIC_STATE_COL].iloc[0]
                first_conf = pd.to_numeric(pd.Series(scored_group["intrinsic_behavioral_cluster_confidence"]), errors="coerce").iloc[0]
                if not pd.isna(first_raw):
                    adata_full.obs.loc[skipped_group.index, HMM_INTRINSIC_RAW_STATE_COL] = first_raw
                if not pd.isna(first_intrinsic):
                    adata_full.obs.loc[skipped_group.index, INTRINSIC_STATE_COL] = first_intrinsic
                if np.isfinite(first_conf):
                    adata_full.obs.loc[skipped_group.index, "intrinsic_behavioral_cluster_confidence"] = first_conf

    adata_full.obs[HMM_INTRINSIC_RAW_STATE_COL] = _coerce_hmm_raw_state_series(
        adata_full.obs[HMM_INTRINSIC_RAW_STATE_COL],
        label=HMM_INTRINSIC_RAW_STATE_COL,
    )
    adata_full.obs[INTRINSIC_STATE_COL] = pd.Categorical(
        pd.Series(adata_full.obs[INTRINSIC_STATE_COL], index=adata_full.obs.index, dtype="string")
    )

    _rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=adata_full,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col=INTRINSIC_STATE_COL,
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=enforce_binary_group_constraints,
    )
    full_label_mapping = {str(k): str(v) for k, v in artifact["full_label_mapping"].items()}
    if str(prepared["resolved_start_offset_fill_mode"]) == "leave_unassigned":
        unassigned_mask = pd.Series(
            adata_full.obs[INTRINSIC_STATE_COL],
            index=adata_full.obs.index,
            dtype="string",
        ).isna()
        for col in ["behavioral_clusterid", "full_behavioral_cluster", FULL_STATE_COL]:
            if col in adata_full.obs.columns:
                adata_full.obs[col] = pd.Series(adata_full.obs[col], index=adata_full.obs.index, dtype="string").where(
                    ~unassigned_mask,
                    pd.NA,
                )
    raw_full_keys = _build_hmm_raw_full_cluster_keys(
        adata_full.obs,
        binary_group_col=BINARY_GROUP_COL,
        raw_state_col=HMM_INTRINSIC_RAW_STATE_COL,
    )
    observed_apply_binary_groups = sorted(
        [
            str(v)
            for v in (
                pd.Series(adata_full.obs[BINARY_GROUP_COL])
                .dropna()
                .astype("string")
                .str.strip()
                .unique()
                .tolist()
            )
            if v != ""
        ],
        key=_mixed_label_sort_key,
    )
    full_label_mapping = _expand_hmm_full_label_mapping_with_raw_state_keys(
        full_label_mapping,
        label_mapping,
        observed_apply_binary_groups,
    )
    curated_full_labels = pd.Series(pd.NA, index=adata_full.obs.index, dtype="string")
    nonmissing_full = raw_full_keys.notna()
    curated_full_labels.loc[nonmissing_full] = raw_full_keys.loc[nonmissing_full].map(full_label_mapping)
    if curated_full_labels.loc[nonmissing_full].isna().any():
        missing_labels = sorted(
            set(raw_full_keys.loc[nonmissing_full][curated_full_labels.loc[nonmissing_full].isna()].astype(str).tolist())
        )
        raise ValueError(
            "HMM deployment artifact is not applicable to this data because canonical original full-cluster "
            "labels are missing from the artifact mapping: "
            f"{missing_labels[:20]}"
        )
    adata_full.obs[FULL_STATE_COL] = pd.Categorical(curated_full_labels)
    adata_full.obs["full_behavioral_cluster"] = pd.Categorical(curated_full_labels)
    adata_full.obs["full_behavioral_cluster_confidence"] = adata_full.obs[
        "intrinsic_behavioral_cluster_confidence"
    ].astype(float)

    if continuous_output_col != INTRINSIC_STATE_COL:
        adata_full.obs[str(continuous_output_col)] = pd.Categorical(
            pd.Series(adata_full.obs[INTRINSIC_STATE_COL], index=adata_full.obs.index, dtype="string")
        )
        adata_full.obs[f"{str(continuous_output_col)}_confidence"] = adata_full.obs[
            "intrinsic_behavioral_cluster_confidence"
        ].astype(float)
    if full_output_col != FULL_STATE_COL:
        adata_full.obs[str(full_output_col)] = pd.Categorical(
            pd.Series(adata_full.obs[FULL_STATE_COL], index=adata_full.obs.index, dtype="string")
        )
        adata_full.obs[f"{str(full_output_col)}_confidence"] = adata_full.obs[
            "full_behavioral_cluster_confidence"
        ].astype(float)

    classification_meta = dict(class_meta)
    classification_meta["applied_continuous_output_col"] = str(continuous_output_col)
    classification_meta["applied_continuous_model_variant"] = "hmm_deployment"
    classification_meta["applied_full_output_col"] = str(full_output_col)
    classification_meta["applied_full_model_variant"] = "hmm_deployment"
    classification_meta["full_label_strategy"] = str(artifact.get("full_label_strategy", "binary_group_plus_intrinsic"))
    classification_meta["hmm_deployment_artifact_kind"] = str(artifact.get("artifact_kind", "hmm_state_deployment"))
    classification_meta["hmm_deployment_artifact_source_type"] = str(artifact_source_type)
    classification_meta["hmm_deployment_artifact_source_path"] = artifact_source_path
    classification_meta["hmm_deployment_source_model_adata_path"] = artifact.get("source_model_adata_path", None)
    classification_meta["intrinsic_label_mapping"] = {str(k): str(v) for k, v in label_mapping.items()}
    classification_meta["full_label_mapping"] = dict(full_label_mapping)
    classification_meta["observed_binary_groups_train"] = list(artifact.get("observed_binary_groups", []))
    artifact_state_colors = _extract_state_color_mapping(
        artifact.get("state_colors", {}),
        FULL_STATE_COL,
    )
    full_labels = _obs_label_values(adata_full, FULL_STATE_COL)
    full_state_colors = _normalize_label_color_map(full_labels, colors=artifact_state_colors)
    if len(full_state_colors) > 0:
        state_colors_meta = dict(classification_meta.get("state_colors", {}) or {})
        state_colors_meta[str(FULL_STATE_COL)] = dict(full_state_colors)
        if str(full_output_col) != FULL_STATE_COL:
            state_colors_meta[str(full_output_col)] = dict(full_state_colors)
        classification_meta["state_colors"] = state_colors_meta
        _set_classification_state_colors(adata_full, FULL_STATE_COL, full_state_colors)
        if str(full_output_col) != FULL_STATE_COL:
            _set_classification_state_colors(adata_full, str(full_output_col), full_state_colors)

    adata_full = _finalize_hmm_apply_outputs(
        adata_full,
        state_paths=state_paths,
        preprocessing_meta=pre_meta,
        classification_meta=classification_meta,
        full_output_col=str(full_output_col),
        export_state_composition=export_state_composition,
        export_state_transitions=export_state_transitions,
        state_transitions_rows_per_page=state_transitions_rows_per_page,
        state_transitions_include_self_pairs=state_transitions_include_self_pairs,
        state_transitions_min_count=state_transitions_min_count,
        state_transitions_relative_count=state_transitions_relative_count,
        verbose=verbose,
    )
    _vdone(verbose, "state-hmm-apply", "apply HMM deployment artifact", apply_started)
    return adata_full
