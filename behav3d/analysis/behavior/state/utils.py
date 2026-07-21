from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from behav3d.analysis.behavior.utils import (
    _mixed_label_sort_key,
    _resolve_output_dir,
    _save_adata_obs_csv,
    _sanitize_filename_token,
    _to_numpy_2d,
    _vdone,
    _vinfo,
    _vsave,
    _vstart,
)


A4_PORTRAIT = (8.27, 11.69)
A4_LANDSCAPE = (11.69, 8.27)

_GLOBAL_DESCRIPTOR_COLUMNS = {
    "summed_displacement": "summed_displacement",
    "net_displacement": "net_displacement",
    "straightness": "straightness",
    "directional_persistence": "directional_persistence",
    "median_turning_angle": "median_turning_angle",
    "fraction_reversed_movement": "fraction_reversed_movement",
    "mean_square_displacement": "mean_square_displacement",
}
_LOG_SCALING_NEGATIVE_TOL = 1e-12


def _coerce_hex_color(value, fallback="#808080"):
    try:
        if value is None:
            raise ValueError("missing color")
        return str(mcolors.to_hex(value, keep_alpha=False))
    except Exception:
        return str(fallback)


def _default_label_color_map(labels, cmap_name="tab20"):
    labels = [] if labels is None else [str(x) for x in list(labels)]
    if len(labels) == 0:
        return {}
    cmap = plt.get_cmap(cmap_name)
    if len(labels) <= getattr(cmap, "N", 256):
        color_values = [cmap(i / max(len(labels) - 1, 1)) for i in range(len(labels))]
    else:
        hsv = plt.get_cmap("hsv")
        color_values = [hsv(i / len(labels)) for i in range(len(labels))]
    return {label: _coerce_hex_color(color_values[i]) for i, label in enumerate(labels)}


def _normalize_label_color_map(labels, colors=None, cmap_name="tab20"):
    labels = [] if labels is None else [str(x) for x in list(labels)]
    defaults = _default_label_color_map(labels, cmap_name=cmap_name)
    raw = {}
    if isinstance(colors, dict):
        raw = {str(k): v for k, v in colors.items()}
    return {
        label: _coerce_hex_color(raw.get(label, defaults.get(label)), fallback=defaults.get(label, "#808080"))
        for label in labels
    }


def _get_classification_state_colors(adata, state_col):
    if adata is None or not hasattr(adata, "uns"):
        return {}
    classification_meta = adata.uns.get("classification", {})
    if not isinstance(classification_meta, dict):
        return {}
    state_colors = classification_meta.get("state_colors", {})
    if not isinstance(state_colors, dict):
        return {}
    nested = state_colors.get(str(state_col), None)
    if isinstance(nested, dict):
        return {str(k): _coerce_hex_color(v) for k, v in nested.items()}
    if all(not isinstance(v, dict) for v in state_colors.values()):
        return {str(k): _coerce_hex_color(v) for k, v in state_colors.items()}
    return {}


def _set_classification_state_colors(adata, state_col, colors):
    if adata is None or not hasattr(adata, "uns"):
        return {}
    classification_meta = adata.uns.get("classification", {})
    if not isinstance(classification_meta, dict):
        classification_meta = {}
    state_colors = classification_meta.get("state_colors", {})
    if not isinstance(state_colors, dict):
        state_colors = {}
    normalized = {
        str(k): _coerce_hex_color(v)
        for k, v in dict(colors or {}).items()
        if str(k).strip() != ""
    }
    state_colors[str(state_col)] = normalized
    classification_meta["state_colors"] = state_colors
    adata.uns["classification"] = classification_meta
    return normalized


def _get_classification_state_order(adata, state_col):
    """Return the user-defined display order for `state_col`'s labels, or [] if none is set."""
    if adata is None or not hasattr(adata, "uns"):
        return []
    classification_meta = adata.uns.get("classification", {})
    if not isinstance(classification_meta, dict):
        return []
    state_order = classification_meta.get("state_order", {})
    if not isinstance(state_order, dict):
        return []
    order = state_order.get(str(state_col), None)
    if isinstance(order, (list, tuple)):
        return [str(x) for x in order]
    return []


def _set_classification_state_order(adata, state_col, order):
    """Persist a user-defined display order for `state_col`'s labels."""
    if adata is None or not hasattr(adata, "uns"):
        return []
    classification_meta = adata.uns.get("classification", {})
    if not isinstance(classification_meta, dict):
        classification_meta = {}
    state_order_map = classification_meta.get("state_order", {})
    if not isinstance(state_order_map, dict):
        state_order_map = {}
    normalized = [str(x) for x in list(order or []) if str(x).strip() != ""]
    state_order_map[str(state_col)] = normalized
    classification_meta["state_order"] = state_order_map
    adata.uns["classification"] = classification_meta
    return normalized


def _apply_state_order(labels, order):
    """Sort `labels` by a stored display `order`.

    Labels found in `order` come first, in that order; any label not in `order` (e.g. added
    after the order was last saved) is appended afterwards, in its original relative order.
    """
    labels = [str(x) for x in list(labels or [])]
    order = [str(x) for x in list(order or [])]
    ranked = {label: i for i, label in enumerate(order)}
    known = sorted((label for label in labels if label in ranked), key=lambda label: ranked[label])
    unknown = [label for label in labels if label not in ranked]
    return known + unknown


def _require_columns(df, required_cols, context):
    missing = [c for c in required_cols if c not in df.columns]
    if len(missing) > 0:
        raise ValueError(
            f"Missing required columns for {context}: {missing}. "
            f"Available columns include: {list(df.columns[:20])}"
        )


def _save_pdf_page_a4(pdf, fig, orientation="portrait"):
    fig.set_size_inches(*(A4_PORTRAIT if orientation == "portrait" else A4_LANDSCAPE), forward=True)
    pdf.savefig(fig, bbox_inches="tight")


def _normalize_log_scale_feature_selectors(log_scale_features):
    if log_scale_features is None:
        return []
    if isinstance(log_scale_features, str):
        selectors = [log_scale_features]
    else:
        selectors = list(log_scale_features)

    normalized = []
    seen = set()
    for selector in selectors:
        s = str(selector).strip()
        if len(s) == 0:
            raise ValueError("log_scale_features cannot contain blank selector names.")
        if s in seen:
            continue
        normalized.append(s)
        seen.add(s)
    return normalized


def _resolve_log_scale_feature_cols(feature_cols, log_scale_features):
    if feature_cols is None:
        feature_cols = []
    feature_cols = [str(c) for c in list(feature_cols)]
    selectors = _normalize_log_scale_feature_selectors(log_scale_features)
    feature_set = set(feature_cols)
    resolved = []
    unresolved = []

    for selector in selectors:
        matched = []
        if selector in feature_set:
            matched.append(selector)
        if selector in _GLOBAL_DESCRIPTOR_COLUMNS:
            col = str(_GLOBAL_DESCRIPTOR_COLUMNS[selector])
            if col in feature_set and col not in matched:
                matched.append(col)
        prefix = f"{selector}_"
        for col in feature_cols:
            if col.startswith(prefix) and col not in matched:
                matched.append(col)

        if len(matched) == 0:
            unresolved.append(selector)
            continue
        for col in matched:
            if col not in resolved:
                resolved.append(col)

    return {
        "requested_features": list(selectors),
        "resolved_feature_cols": list(resolved),
        "unresolved_features": list(unresolved),
    }


def _coerce_log_scaling_params(preprocessing_params):
    if not isinstance(preprocessing_params, dict):
        return None
    log_meta = preprocessing_params.get("log_scaling", None)
    if not isinstance(log_meta, dict):
        return None

    requested = log_meta.get("requested_features", [])
    resolved = log_meta.get("resolved_feature_cols", [])
    transform = log_meta.get("transform", None)
    requested = _normalize_log_scale_feature_selectors(requested)
    resolved = [str(c) for c in list(resolved or [])]
    if transform is None and len(resolved) > 0:
        transform = "log1p"
    if len(resolved) == 0 and len(requested) == 0 and transform is None:
        return None
    return {
        "transform": None if transform is None else str(transform),
        "requested_features": list(requested),
        "resolved_feature_cols": list(resolved),
    }


def _log_scaling_feature_indices(feature_cols, resolved_feature_cols, *, strict=True):
    if feature_cols is None:
        feature_cols = []
    if resolved_feature_cols is None:
        resolved_feature_cols = []
    feature_cols = [str(c) for c in list(feature_cols)]
    resolved_feature_cols = [str(c) for c in list(resolved_feature_cols)]
    col_to_idx = {col: idx for idx, col in enumerate(feature_cols)}
    missing = [col for col in resolved_feature_cols if col not in col_to_idx]
    if strict and len(missing) > 0:
        raise ValueError(
            "Log-scaling metadata references feature columns that are missing from the current matrix: "
            f"{missing[:20]}"
        )
    return [col_to_idx[col] for col in resolved_feature_cols if col in col_to_idx]


def _apply_log1p_to_feature_matrix(
    X,
    feature_cols,
    resolved_feature_cols,
    *,
    inplace=False,
    strict=True,
):
    out = np.asarray(X, dtype=float)
    if not inplace:
        out = out.copy()

    feat_indices = _log_scaling_feature_indices(
        feature_cols,
        resolved_feature_cols,
        strict=strict,
    )
    for idx in feat_indices:
        vals = out[:, idx]
        finite = np.isfinite(vals)
        if not np.any(finite):
            continue
        finite_vals = vals[finite]
        min_val = float(np.min(finite_vals))
        if min_val < -_LOG_SCALING_NEGATIVE_TOL:
            raise ValueError(
                "Cannot apply log1p to feature values below zero. "
                f"feature='{feature_cols[idx]}', min_value={min_val}"
            )
        if min_val < 0:
            vals = vals.copy()
            vals[finite & (vals < 0)] = 0.0
            out[:, idx] = vals
        out[:, idx] = np.log1p(out[:, idx])
    return out


def _invert_log1p_in_feature_matrix(
    X,
    feature_cols,
    resolved_feature_cols,
    *,
    inplace=False,
    strict=True,
):
    out = np.asarray(X, dtype=float)
    if not inplace:
        out = out.copy()

    feat_indices = _log_scaling_feature_indices(
        feature_cols,
        resolved_feature_cols,
        strict=strict,
    )
    for idx in feat_indices:
        out[:, idx] = np.expm1(out[:, idx])
    return out


def _apply_log_scaling_to_continuous_matrix(
    X,
    preprocessing_params,
    feature_cols,
    *,
    inplace=False,
    strict=True,
):
    out = np.asarray(X, dtype=float)
    if not inplace:
        out = out.copy()

    log_meta = _coerce_log_scaling_params(preprocessing_params)
    if log_meta is None:
        return out
    transform = log_meta.get("transform", None)
    if transform != "log1p":
        raise ValueError(f"Unsupported log scaling transform: {transform}")
    return _apply_log1p_to_feature_matrix(
        out,
        feature_cols=feature_cols,
        resolved_feature_cols=log_meta.get("resolved_feature_cols", []),
        inplace=True,
        strict=strict,
    )


def _invert_log_scaling_in_continuous_matrix(
    X,
    preprocessing_params,
    feature_cols,
    *,
    inplace=False,
    strict=True,
):
    out = np.asarray(X, dtype=float)
    if not inplace:
        out = out.copy()

    log_meta = _coerce_log_scaling_params(preprocessing_params)
    if log_meta is None:
        return out
    transform = log_meta.get("transform", None)
    if transform != "log1p":
        raise ValueError(f"Unsupported log scaling transform: {transform}")
    return _invert_log1p_in_feature_matrix(
        out,
        feature_cols=feature_cols,
        resolved_feature_cols=log_meta.get("resolved_feature_cols", []),
        inplace=True,
        strict=strict,
    )


def cap_values_to_quantile(df, features, lower_quantile=None, upper_quantile=0.99, return_limits=False):
    if lower_quantile is None and upper_quantile is None:
        if return_limits:
            return df, {}
        return df

    df_capped = df.copy()
    cap_limits = {}

    for feature in features:
        if feature not in df_capped.columns:
            continue
        col = df_capped[feature]
        if pd.api.types.is_bool_dtype(col):
            if return_limits:
                cap_limits[feature] = {"lower": None, "upper": None}
            continue
        series = pd.to_numeric(col, errors="coerce").dropna()
        if len(series) == 0:
            if return_limits:
                cap_limits[feature] = {"lower": None, "upper": None}
            continue
        lower_limit = series.quantile(lower_quantile) if lower_quantile is not None else None
        upper_limit = series.quantile(upper_quantile) if upper_quantile is not None else None
        numeric_col = pd.to_numeric(df_capped[feature], errors="coerce")
        df_capped[feature] = numeric_col.clip(lower=lower_limit, upper=upper_limit)
        if return_limits:
            cap_limits[feature] = {
                "lower": None if lower_limit is None else float(lower_limit),
                "upper": None if upper_limit is None else float(upper_limit),
            }

    if return_limits:
        return df_capped, cap_limits
    return df_capped


def _binary_col_to_group_name(col: str) -> str:
    return str(col)


def _binary_value_is_active(val) -> bool:
    if pd.isna(val):
        return False
    try:
        return float(val) == 1.0
    except Exception:
        return False


def _active_binary_groups_from_row(row, binary_cols: list[str]) -> list[str]:
    active = []
    for col in binary_cols:
        if _binary_value_is_active(row[col]):
            active.append(_binary_col_to_group_name(col))
    return active


def _binary_groups_to_label(active_groups: list[str]) -> str:
    if len(active_groups) == 0:
        return "no_contact"
    if len(active_groups) == 1:
        return str(active_groups[0])
    return "_and_".join([str(x) for x in active_groups])


def _infer_binary_group_constraints(df: pd.DataFrame, binary_cols: list[str]) -> dict:
    cols = [str(c) for c in list(binary_cols or []) if str(c) in df.columns]
    clean_cols = [_binary_col_to_group_name(c) for c in cols]

    if len(cols) == 0:
        return {
            "inference_rule": "exact_zero_cooccurrence",
            "binary_cols": [],
            "binary_group_names": [],
            "allowed_binary_groups": ["no_contact"],
            "forbidden_binary_combinations": {},
            "support_counts": {},
        }

    support_counts = {name: 0 for name in clean_cols}
    cooccur_counts = {}
    allowed_groups = set()

    for _, row in df[cols].iterrows():
        active = _active_binary_groups_from_row(row, cols)
        allowed_groups.add(_binary_groups_to_label(active))
        for grp in active:
            support_counts[grp] = int(support_counts.get(grp, 0)) + 1
        for size in range(2, len(active) + 1):
            for combo in combinations(active, size):
                combo_key = tuple(sorted([str(x) for x in combo]))
                cooccur_counts[combo_key] = int(cooccur_counts.get(combo_key, 0)) + 1

    forbidden_combinations = []
    supported_groups = [str(name) for name, count in support_counts.items() if int(count) > 0]
    for size in range(2, len(supported_groups) + 1):
        for combo in combinations(supported_groups, size):
            combo_key = tuple(sorted([str(x) for x in combo]))
            if int(cooccur_counts.get(combo_key, 0)) == 0:
                forbidden_combinations.append(list(combo_key))

    return {
        "inference_rule": "exact_zero_cooccurrence",
        "binary_cols": list(cols),
        "binary_group_names": list(clean_cols),
        "allowed_binary_groups": sorted([str(x) for x in allowed_groups], key=_mixed_label_sort_key),
        "forbidden_binary_combinations": _group_forbidden_binary_combinations(forbidden_combinations),
        "support_counts": {str(k): int(v) for k, v in support_counts.items()},
    }


def _group_forbidden_binary_combinations(forbidden_combinations):
    grouped = {}
    for item in list(forbidden_combinations or []):
        if not isinstance(item, (list, tuple)) or len(item) < 2:
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] entries must be length-2+ lists/tuples."
            )
        combo = tuple(sorted([str(x) for x in item]))
        if len(set(combo)) != len(combo):
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] entries must not contain duplicate group names."
            )
        arity = str(len(combo))
        grouped.setdefault(arity, set()).add(combo)

    return {
        str(arity): [list(combo) for combo in sorted(combos, key=lambda c: tuple(c))]
        for arity, combos in sorted(grouped.items(), key=lambda kv: int(kv[0]))
    }


def _normalize_binary_group_constraints(binary_group_constraints):
    if binary_group_constraints is None:
        return None
    if not isinstance(binary_group_constraints, dict):
        raise ValueError("binary_group_constraints must be a dict when provided.")

    allowed_raw = binary_group_constraints.get("allowed_binary_groups", None)
    forbidden_raw = binary_group_constraints.get("forbidden_binary_combinations", None)
    if allowed_raw is None and forbidden_raw is None:
        raise ValueError("binary_group_constraints['forbidden_binary_combinations'] is required.")

    allowed_set = None
    if allowed_raw is not None:
        if not hasattr(allowed_raw, "__iter__") or isinstance(allowed_raw, (str, dict, set)):
            raise ValueError("binary_group_constraints['allowed_binary_groups'] must be a list or sequence.")
        allowed_set = set(str(x) for x in allowed_raw)

    if forbidden_raw is None:
        raise ValueError("binary_group_constraints['forbidden_binary_combinations'] is required.")
    if not isinstance(forbidden_raw, dict):
        raise ValueError("binary_group_constraints['forbidden_binary_combinations'] must be a dict grouped by arity.")
    forbidden_set = set()
    for arity, combos in forbidden_raw.items():
        try:
            arity_int = int(str(arity))
        except Exception as exc:
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] keys must be integer arities."
            ) from exc
        if arity_int < 2:
            raise ValueError("binary_group_constraints['forbidden_binary_combinations'] arity keys must be >= 2.")
        if not hasattr(combos, "__iter__") or isinstance(combos, (str, dict, set)):
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] values must be lists of combinations."
            )
        for item in combos:
            if not isinstance(item, (list, tuple, np.ndarray)) or len(item) != arity_int:
                raise ValueError(
                    "binary_group_constraints['forbidden_binary_combinations'] entries must match their arity key."
                )
            combo = tuple(sorted([str(x) for x in item]))
            if len(set(combo)) != len(combo):
                raise ValueError(
                    "binary_group_constraints['forbidden_binary_combinations'] entries must not contain duplicate group names."
                )
            forbidden_set.add(combo)

    return {
        "allowed_binary_groups": allowed_set,
        "forbidden_binary_combinations": forbidden_set,
    }


def _assign_binary_group_labels(
    df: pd.DataFrame,
    binary_cols: list[str],
    binary_group_constraints=None,
    enforce_binary_group_constraints: bool = False,
) -> pd.Series:
    if len(binary_cols) == 0:
        if enforce_binary_group_constraints and binary_group_constraints is not None:
            normalized = _normalize_binary_group_constraints(binary_group_constraints)
            if normalized is not None:
                allowed = normalized.get("allowed_binary_groups", None)
                if allowed is not None and "no_contact" not in allowed:
                    raise ValueError(
                        "Binary group constraints violated: 'no_contact' is required for empty binary_cols."
                    )
        return pd.Series(["no_contact"] * len(df), index=df.index, dtype="object")

    normalized_constraints = None
    if binary_group_constraints is not None:
        normalized_constraints = _normalize_binary_group_constraints(binary_group_constraints)

    labels = []
    forbidden_hits = {}
    unknown_hits = {}
    for _, row in df[binary_cols].iterrows():
        active = _active_binary_groups_from_row(row, binary_cols)
        label = _binary_groups_to_label(active)
        labels.append(label)

        if not enforce_binary_group_constraints or normalized_constraints is None:
            continue

        row_idx = str(row.name)
        forbidden_combinations = normalized_constraints.get("forbidden_binary_combinations", set())
        for size in range(2, len(active) + 1):
            for combo in combinations(active, size):
                combo_key = tuple(sorted([str(x) for x in combo]))
                if combo_key not in forbidden_combinations:
                    continue
                hit = forbidden_hits.setdefault(
                    combo_key,
                    {"count": 0, "first_index": row_idx},
                )
                hit["count"] = int(hit["count"]) + 1

        allowed_groups = normalized_constraints.get("allowed_binary_groups", None)
        if allowed_groups is not None and label not in allowed_groups:
            hit = unknown_hits.setdefault(
                label,
                {"count": 0, "first_index": row_idx},
            )
            hit["count"] = int(hit["count"]) + 1

    if enforce_binary_group_constraints and (len(forbidden_hits) > 0 or len(unknown_hits) > 0):
        parts = []
        if len(forbidden_hits) > 0:
            combo_msgs = []
            for combo, info in sorted(forbidden_hits.items(), key=lambda x: (len(x[0]), x[0])):
                combo_msgs.append(
                    f"{'+'.join(list(combo))} (n={int(info['count'])}, first_index={info['first_index']})"
                )
            parts.append("forbidden_combinations=[" + "; ".join(combo_msgs[:8]) + "]")
        if len(unknown_hits) > 0:
            group_msgs = []
            for group, info in sorted(unknown_hits.items(), key=lambda x: x[0]):
                group_msgs.append(
                    f"{group} (n={int(info['count'])}, first_index={info['first_index']})"
                )
            parts.append("unknown_binary_groups=[" + "; ".join(group_msgs[:8]) + "]")
        raise ValueError(
            "Binary group constraints violated while rebuilding full behavioral clusters. "
            + " | ".join(parts)
            + " | This indicates the supplied data contains binary-group combinations not present in training."
        )

    return pd.Series(labels, index=df.index, dtype="object")


def _add_clean_binary_annotation_columns(df: pd.DataFrame, binary_cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in binary_cols:
        if col not in out.columns:
            continue
        clean_name = _binary_col_to_group_name(col)
        if clean_name != col and clean_name in out.columns:
            continue
        vals = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
        out[clean_name] = (vals > 0).astype(np.int8)
    return out


def _rebuild_full_behavioral_cluster_from_intrinsic(
    adata,
    binary_cols_to_merge,
    intrinsic_col="intrinsic_behavioral_cluster",
    binary_group_constraints=None,
    enforce_binary_group_constraints=False,
):
    if intrinsic_col not in adata.obs.columns:
        raise ValueError(f"Missing '{intrinsic_col}' in adata.obs.")

    _bcols = binary_cols_to_merge if binary_cols_to_merge is not None else []
    binary_cols = [str(c) for c in list(_bcols)]
    adata.obs = _add_clean_binary_annotation_columns(adata.obs, binary_cols)
    adata.obs["behavioral_clusterid"] = adata.obs[intrinsic_col].astype(str)
    adata.obs["binary_group"] = _assign_binary_group_labels(
        adata.obs,
        binary_cols,
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=bool(enforce_binary_group_constraints),
    ).astype("category")
    adata.obs["full_behavioral_cluster"] = (
        adata.obs["binary_group"].astype(str) + "_" + adata.obs["behavioral_clusterid"].astype(str)
    ).astype("category")
    adata.obs["behavioral_clusterid"] = adata.obs["behavioral_clusterid"].astype("category")
    return adata


@dataclass(frozen=True)
class StatePaths:
    output_dir: Path
    analysis_outdir: Path
    feature_outdir: Path
    filtered_tracks_csv: Path
    combined_tracks_csv: Path
    state_outdir: Path
    processing_outdir: Path
    state_clustering_outdir: Path
    model_adata_path: Path
    prepared_full_adata_path: Path
    full_output_adata_path: Path
    intrinsic_outdir: Path
    full_outdir: Path
    intrinsic_qc_outdir: Path
    full_qc_outdir: Path
    intrinsic_classifier_default_path: Path
    full_classifier_default_path: Path
    state_composition_outdir: Path
    state_transitions_outdir: Path


def _resolve_state_outdir(output_dir, cell_type):
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    state_outdir = root / "analysis" / str(cell_type) / "behavioral_states"
    state_outdir.mkdir(parents=True, exist_ok=True)
    return state_outdir


def _resolve_state_paths(output_dir, cell_type):
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    feature_outdir = analysis_outdir / "track_features"

    state_outdir = root / "analysis" / str(cell_type) / "behavioral_states"
    state_outdir.mkdir(parents=True, exist_ok=True)
    processing_outdir = state_outdir / "processing"
    processing_outdir.mkdir(parents=True, exist_ok=True)

    intrinsic_outdir = processing_outdir / "intrinsic_behavioral_classification"
    full_outdir = processing_outdir / "full_behavioral_classification"
    intrinsic_qc_outdir = intrinsic_outdir / "quality_control"
    full_qc_outdir = full_outdir / "quality_control"

    return StatePaths(
        output_dir=root,
        analysis_outdir=analysis_outdir,
        feature_outdir=feature_outdir,
        filtered_tracks_csv=feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
        combined_tracks_csv=feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv",
        state_outdir=state_outdir,
        processing_outdir=processing_outdir,
        state_clustering_outdir=processing_outdir / "state_clustering",
        model_adata_path=processing_outdir / f"BEHAV3D_{cell_type}_behavioral_states_modeldata.h5ad",
        prepared_full_adata_path=processing_outdir / f"BEHAV3D_{cell_type}_behavioral_states_prepared_full.h5ad",
        full_output_adata_path=state_outdir / f"BEHAV3D_{cell_type}_behavioral_states.h5ad",
        intrinsic_outdir=intrinsic_outdir,
        full_outdir=full_outdir,
        intrinsic_qc_outdir=intrinsic_qc_outdir,
        full_qc_outdir=full_qc_outdir,
        intrinsic_classifier_default_path=intrinsic_outdir / f"intrinsic_state_random_forest_{cell_type}.pkl",
        full_classifier_default_path=full_outdir / f"state_random_forest_{cell_type}.pkl",
        state_composition_outdir=state_outdir / "state_composition",
        state_transitions_outdir=state_outdir / "state_transitions",
    )


def _resolve_analysis_paths(output_dir, cell_type):
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    feature_outdir = analysis_outdir / "track_features"
    return {
        "output_dir": root,
        "analysis_outdir": analysis_outdir,
        "feature_outdir": feature_outdir,
        "filtered_tracks_csv": feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
        "combined_tracks_csv": feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv",
    }


def _resolve_state_stage_paths(output_dir, cell_type):
    return asdict(_resolve_state_paths(output_dir, cell_type))


def _resolve_positions_csv_path(output_dir=None, cell_type=None, state_paths=None):
    if state_paths is None:
        analysis_paths = _resolve_analysis_paths(output_dir, cell_type)
        filtered_path = analysis_paths["filtered_tracks_csv"]
        combined_path = analysis_paths["combined_tracks_csv"]
    else:
        filtered_path = state_paths.filtered_tracks_csv
        combined_path = state_paths.combined_tracks_csv
    if filtered_path.exists():
        return filtered_path
    if combined_path.exists():
        return combined_path
    raise FileNotFoundError(
        "Could not find track-features CSV for state classification. "
        f"Expected one of: '{filtered_path}' or '{combined_path}'. "
        f"Run feature extraction for cell_type='{cell_type}' first, or place the CSV in "
        "output_dir/analysis/<cell_type>/track_features."
    )
