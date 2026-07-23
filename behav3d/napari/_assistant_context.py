"""
BEHAV3D assistant — workflow context serializer.

``build_context(main_widget)`` snapshots *where the user is* and *what they have
configured* into a plain JSON-serialisable dict that is sent to the model on
every turn. It is intentionally defensive: tabs may be half-initialised and the
metadata DataFrame may be empty or missing expected columns (see the historical
``KeyError: 'sample_name'`` bug), so every access is guarded.

The returned dict shape:

    {
      "current_step": "segmentation",
      "current_tab_label": "🦠 Segmentation",
      "tab_index": 2,
      "output_dir": "/path" | "",
      "output_dir_set": true,
      "metadata": {
          "loaded": true,
          "n_samples": 3,
          "sample_names": [...],
          "columns": [...],
          "cell_types": {"immune": [...], "organoid": [...], "other": [...], "merged": [...]},
      },
      "experiment_reference": {
          "notes": [{"source": "README_BEHAV3D_Exp010.md", "text": "..."}],
          "saved_configurations": [{"source": "...yml", "settings": {...}}],
      },
      "queue": [{"type": "track", "label": "📍 Batch Tracking", "status": "pending", "params": {...}}],
      "parameters": {                      # only values that differ from defaults
          "tracking.immune.method": {"current": "btrack", "default": "trackpy"},
          ...
      },
      "step_schema": [ <cards for current_step> ],
    }
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional

from behav3d.napari._assistant_controls import (
    CONTROL_CONTRACT_VERSION,
    active_cell_type,
    control_registry,
)

# Tab index -> workflow step key (mirrors the order tabs are added in _widget.py)
_TAB_INDEX_TO_STEP = {
    0: "data_preparation",
    1: "visualization",
    2: "segmentation",
    3: "tracking",
    4: "feature_extraction",
    5: "filtering",
    6: "analysis",
}


def _safe(fn, default=None):
    try:
        return fn()
    except Exception:
        return default


_DIMENSION_ORDERS = {"TCZYX", "TZCYX", "ZCTYX", "ZTCYX", "CZTYX", "CTZYX"}
_RAW_IMAGE_SUFFIXES = (".zarr", ".zarr.zip", ".czi", ".lif", ".liff",
                       ".tif", ".tiff", ".ims", ".h5")
_EXPERIMENT_NOTE_PATTERNS = (
    "README_BEHAV3D*.md",
    "EXPERIMENT_CONTEXT*.md",
    "BEHAV3D_CONTEXT*.md",
)
_EXPERIMENT_CONFIG_PATTERNS = (
    "behav3d_parameters.yml",
    "behav3d_parameters*.yml",
    "behav3d_parameters*.yaml",
)
_EXPERIMENT_NOTE_CHAR_LIMIT = 16_000


def _json_value(value):
    """Convert pandas/numpy/path values to JSON primitives without hiding data."""
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    try:
        import pandas as pd
        if bool(pd.isna(value)):
            return None
    except (ImportError, TypeError, ValueError):
        pass
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return _safe(value.item, str(value))
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _reference_roots(output_dir: str, parameters: dict) -> list[Path]:
    roots: list[Path] = []
    candidates = [output_dir]
    paths = parameters.get("paths") if isinstance(parameters, dict) else None
    metadata_path = paths.get("metadata_csv") if isinstance(paths, dict) else None
    if metadata_path:
        candidates.append(str(Path(str(metadata_path)).expanduser().parent))
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(str(candidate)).expanduser()
        if path.is_dir() and path not in roots:
            roots.append(path)
    return roots


def _compact_experiment_config(config: dict) -> dict:
    """Keep interpretation-relevant settings without paths or large model features."""
    if not isinstance(config, dict):
        return {}
    summary: dict[str, Any] = {}

    pixel = config.get("pixel_classifier")
    if isinstance(pixel, dict):
        summary["segmentation"] = {
            key: _json_value(pixel[key])
            for key in (
                "apoc_strategy", "convpaint_strategy", "examples_per_sample",
                "use_all_timepoints", "tp_start", "tp_end",
            )
            if key in pixel
        }

    apoc = config.get("apoc")
    if isinstance(apoc, dict):
        suffixes = (
            "prob_mask_threshold", "prob_seed_threshold", "edt_threshold",
            "segment_size_min", "opening_nr_pixels", "fill_holes",
            "max_depth", "num_ensembles",
        )
        by_cell_type: dict[str, dict] = {}
        for key, value in apoc.items():
            match = re.match(r"^apoc_(.+)_(%s)$" % "|".join(suffixes), str(key))
            if match:
                cell_type, setting = match.groups()
                by_cell_type.setdefault(cell_type, {})[setting] = _json_value(value)
        if by_cell_type:
            summary.setdefault("segmentation", {})["apoc_by_cell_type"] = by_cell_type

    tracking = config.get("tracking")
    if isinstance(tracking, dict):
        tracking_summary = {}
        if "track_organoids_together" in tracking:
            tracking_summary["track_organoids_together"] = _json_value(
                tracking["track_organoids_together"]
            )
        for cell_type, settings in tracking.items():
            if not isinstance(settings, dict) or "method" not in settings:
                continue
            item = {"method": _json_value(settings.get("method"))}
            btrack = settings.get("btrack")
            if isinstance(btrack, dict):
                item["btrack"] = {
                    key: _json_value(btrack[key])
                    for key in (
                        "max_search_radius", "use_visual_features", "use_optimize",
                        "dist_thresh", "time_thresh", "step_size",
                    )
                    if key in btrack
                }
            tracking_summary[str(cell_type)] = item
        if tracking_summary:
            summary["tracking"] = tracking_summary

    feature_summary = {}
    features = config.get("features")
    for cell_type, settings in (
        features.items() if isinstance(features, dict) else ()
    ):
        if not isinstance(settings, dict):
            continue
        item = {
            key: _json_value(settings[key])
            for key in (
                "features_choice", "contact_threshold",
                "dead_mask_percentage_threshold",
            )
            if key in settings
        }
        if item:
            feature_summary[str(cell_type)] = item
    if feature_summary:
        summary["features"] = feature_summary

    filtering_summary = {}
    filtering = config.get("track_filtering")
    for cell_type, settings in (
        filtering.items() if isinstance(filtering, dict) else ()
    ):
        if not isinstance(settings, dict):
            continue
        item = {
            key: _json_value(settings[key])
            for key in (
                "exp_duration_enabled", "exp_duration", "min_length_enabled",
                "min_track_length", "max_length_enabled", "max_track_length",
                "split_long_tracks", "filter_min_size_t1", "min_size_t1",
                "filter_t0_dead", "time_type",
            )
            if key in settings
        }
        if item:
            filtering_summary[str(cell_type)] = item
    if filtering_summary:
        summary["filtering"] = filtering_summary

    active_killing = config.get("active_killing")
    if isinstance(active_killing, dict):
        summary["active_killing"] = {
            key: _json_value(active_killing[key])
            for key in (
                "observation_window", "death_signal_column",
                "killing_threshold_multiplier", "use_absolute_threshold",
                "absolute_killing_threshold", "min_contact_duration",
            )
            if key in active_killing
        }
    death_dynamics = config.get("death_dynamics")
    if isinstance(death_dynamics, dict):
        summary["death_dynamics"] = _json_value(death_dynamics)

    module_specs = {
        "behavioral_state_classification": (
            "hmm_n_states_mode", "hmm_n_states",
            "hmm_feature_smoothing_window", "hmm_start_offset",
            "selected_features", "binary_features_to_group",
        ),
        "behavioral_track_classification": (
            "behavioral_trajectory_size", "n_clusters",
            "trajectory_trim_mode", "linkage",
        ),
    }
    for module, keys in module_specs.items():
        settings = config.get(module)
        if not isinstance(settings, dict):
            continue
        defaults = settings.get("defaults", settings)
        if isinstance(defaults, dict):
            summary[module] = {
                key: _json_value(defaults[key])
                for key in keys if key in defaults
            }
    for module in ("state_classification", "track_classification"):
        settings = config.get(module)
        if isinstance(settings, dict):
            summary[module] = {
                "cell_types_present": [str(key) for key in settings],
                "cell_types_with_nonempty_settings": [
                    str(key) for key, value in settings.items() if value
                ],
            }
    return summary


def _experiment_reference_context(output_dir: str, parameters: dict) -> dict | None:
    """Discover optional, dataset-specific notes and a compact saved configuration."""
    roots = _reference_roots(output_dir, parameters)
    if not roots:
        return None

    note_paths: list[Path] = []
    config_paths: list[Path] = []
    for root in roots:
        for pattern in _EXPERIMENT_NOTE_PATTERNS:
            note_paths.extend(sorted(root.glob(pattern)))
        for pattern in _EXPERIMENT_CONFIG_PATTERNS:
            config_paths.extend(sorted(root.glob(pattern)))
    note_paths = list(dict.fromkeys(path for path in note_paths if path.is_file()))
    config_paths = list(dict.fromkeys(path for path in config_paths if path.is_file()))

    remaining = _EXPERIMENT_NOTE_CHAR_LIMIT
    notes = []
    for path in note_paths:
        if remaining <= 0:
            break
        text = _safe(lambda p=path: p.read_text(encoding="utf-8"), "") or ""
        excerpt = text[:remaining]
        if excerpt:
            notes.append({
                "source": path.name,
                "text": excerpt,
                "truncated": len(excerpt) < len(text),
            })
            remaining -= len(excerpt)

    saved_configs = []
    for path in config_paths[:3]:
        try:
            import yaml
            loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        compact = _compact_experiment_config(loaded)
        if compact:
            saved_configs.append({"source": path.name, "settings": compact})

    if not notes and not saved_configs:
        return None
    return {
        "notes": notes,
        "saved_configurations": saved_configs,
        "provenance": (
            "User-provided, dataset-specific reference files discovered beside "
            "the selected output directory or metadata file."
        ),
        "configuration_caveat": (
            "Saved settings describe configured intent and are not proof that a "
            "module ran. Confirm execution from discovered result files."
        ),
    }


def validate_metadata_records(records: list[dict]) -> list[dict]:
    """Return factual metadata issues and review notes without biological guesses."""
    issues: list[dict] = []

    def add(row, field, severity, message):
        issues.append({"row": row, "field": field, "severity": severity,
                       "message": message})

    required = ("sample_name", "raw_image_path", "pixel_distance_xy",
                "pixel_distance_z", "time_interval")
    for index, record in enumerate(records):
        sample = record.get("sample_name") or f"row {index + 1}"
        for field in required:
            value = record.get(field)
            if value is None or (isinstance(value, str) and not value.strip()):
                add(index, field, "error", f"{sample}: {field.replace('_', ' ')} is missing.")

        raw = str(record.get("raw_image_path") or "").strip()
        if raw:
            path = Path(raw).expanduser()
            lower = raw.lower()
            if path.is_dir() and not lower.endswith(".zarr"):
                add(index, "raw_image_path", "error",
                    f"{sample}: the raw image path is a folder; select one multidimensional image file.")
            elif not any(lower.endswith(suffix) for suffix in _RAW_IMAGE_SUFFIXES):
                add(index, "raw_image_path", "error",
                    f"{sample}: the raw image format is not supported by the loader.")
            elif not path.exists():
                add(index, "raw_image_path", "warning",
                    f"{sample}: the raw image path does not exist on this machine.")

        order = str(record.get("dimension_order") or "").strip().upper()
        if order and order not in _DIMENSION_ORDERS:
            add(index, "dimension_order", "error",
                f"{sample}: dimension order must be one of {', '.join(sorted(_DIMENSION_ORDERS))}.")
        for field in ("pixel_distance_xy", "pixel_distance_z", "time_interval"):
            value = record.get(field)
            if value is None or value == "":
                continue
            try:
                if float(value) <= 0:
                    add(index, field, "error",
                        f"{sample}: {field.replace('_', ' ')} must be greater than zero.")
            except (TypeError, ValueError):
                add(index, field, "error",
                    f"{sample}: {field.replace('_', ' ')} must be numeric.")

    names = [str(r.get("sample_name") or "").strip() for r in records]
    duplicates = sorted({name for name in names if name and names.count(name) > 1})
    for name in duplicates:
        add(None, "sample_name", "error", f"Sample name '{name}' is used more than once.")

    for field in ("pixel_distance_xy", "pixel_distance_z", "time_interval",
                  "time_unit", "dimension_order"):
        values = {_json_value(r.get(field)) for r in records if r.get(field) not in (None, "")}
        if len(values) > 1:
            add(None, field, "review",
                f"{field.replace('_', ' ').capitalize()} differs between samples; verify that this is intentional.")
    return issues


def summarize_metadata(metadata) -> dict:
    """Serialize every loaded metadata row and validate researcher-facing fields."""
    if metadata is None:
        return {"loaded": False}
    try:
        empty = bool(getattr(metadata, "empty", True))
    except Exception:
        empty = True
    if empty:
        return {"loaded": False}

    columns = _safe(lambda: list(metadata.columns), []) or []
    n_samples = _safe(lambda: int(len(metadata)), 0)

    sample_names = []
    if "sample_name" in columns:
        sample_names = _safe(
            lambda: [str(s) for s in metadata["sample_name"].dropna().unique()], []
        ) or []

    cell_types: dict[str, list] = {}
    try:
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata,
            detect_merged_cell_types_from_metadata,
        )
        cell_types = {
            "immune": list(_safe(lambda: detect_immune_cell_types_from_metadata(metadata), []) or []),
            "organoid": list(_safe(lambda: detect_organoid_types_from_metadata(metadata), []) or []),
            "other": list(_safe(lambda: detect_other_cell_types_from_metadata(metadata), []) or []),
            "merged": list(_safe(lambda: detect_merged_cell_types_from_metadata(metadata), []) or []),
        }
    except Exception:
        cell_types = {}

    records = []
    try:
        for _index, row in metadata.iterrows():
            records.append({str(column): _json_value(row.get(column)) for column in columns})
    except Exception:
        records = []

    return {
        "loaded": True,
        "n_samples": n_samples,
        "sample_names": sample_names,
        "columns": columns,
        "cell_types": cell_types,
        "records": records,
        "validation": validate_metadata_records(records),
    }


def _diff_from_defaults(behav3d_parameters: dict) -> dict:
    """Return {dotted_key: {current, default}} for every parameter whose current
    value differs from the schema default."""
    from behav3d.napari._assistant_schema import flatten_config_to_cards
    from behav3d.napari._assistant_actions import get_by_dotted  # local helper

    diffs: dict[str, dict] = {}
    for card in flatten_config_to_cards():
        key = card["key"]
        if key.startswith("calculated_features."):
            continue
        current = get_by_dotted(behav3d_parameters, key, _MISSING)
        if current is _MISSING:
            continue
        if current != card["default"]:
            diffs[key] = {"current": current, "default": card["default"]}
    return diffs


class _Missing:
    pass


_MISSING = _Missing()


def serialize_queue(queue_panel) -> list[dict]:
    steps = getattr(queue_panel, "_steps", None) or []
    out = []
    for step in steps:
        out.append({
            "type": _safe(lambda s=step: s.step_type.value, "?"),
            "label": _safe(lambda s=step: s.display_label, ""),
            "status": _safe(lambda s=step: s.status.value, "?"),
            "params": _safe(lambda s=step: dict(s.params), {}) or {},
        })
    return out


def _serialize_results(output_dir: str) -> list[dict]:
    if not output_dir:
        return []
    try:
        from behav3d.napari._results_catalog import scan_outputs
        root = Path(output_dir).expanduser()
        return [{
            "id": str(item.path.relative_to(root)),
            "label": item.label,
            "description": item.description,
            "kind": item.kind,
            "category": item.category,
            "subcategory": item.subcategory,
            "cell_type": item.cell_type,
            "viewable": item.is_viewable,
        } for item in scan_outputs(root)]
    except Exception:
        return []


def _active_preview(main_widget) -> dict | None:
    feature_tab = getattr(main_widget, "feature_extraction_tab", None)
    for cell_type, panel in (getattr(feature_tab, "panels", {}) or {}).items():
        if getattr(panel, "_preview_seg_t", None) is not None:
            return {
                "type": "dead_threshold",
                "cell_type": cell_type,
                "sample": _widget_value(getattr(panel, "preview_sample_combo", None)),
                "threshold_percent": _widget_value(getattr(panel, "spin_dead_threshold", None)),
                "overlay": "green cells are alive; red cells are dead",
                "hover": "hover a cell to see its dead-mask percentage",
            }
    return None


def _widget_value(widget):
    """Best-effort JSON-friendly value for a Qt input widget."""
    if widget is None:
        return None
    try:
        from qtpy.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QLineEdit, QSpinBox
        if isinstance(widget, QLineEdit):
            return widget.text()
        if isinstance(widget, QComboBox):
            return widget.currentText()
        if isinstance(widget, QCheckBox):
            return widget.isChecked()
        if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
            return widget.value()
    except Exception:
        pass
    if hasattr(widget, "currentText"):
        return _safe(lambda: widget.currentText(), _MISSING)
    if hasattr(widget, "value"):
        return _safe(lambda: widget.value(), _MISSING)
    if hasattr(widget, "isChecked"):
        return _safe(lambda: widget.isChecked(), _MISSING)
    if hasattr(widget, "text"):
        return _safe(lambda: widget.text(), _MISSING)
    return None


def _clean_widget_values(values: dict) -> dict:
    return {k: v for k, v in values.items() if v is not _MISSING}


def _metadata_builder_state(dp) -> dict:
    """Snapshot the current state of the Metadata Builder form widgets."""
    if dp is None:
        return {"open": False}
    try:
        n_org = len(getattr(dp, "_organoid_name_edits", []))
        n_imm = len(getattr(dp, "_immune_name_edits", []))
        n_oth = len(getattr(dp, "_other_name_edits", []))
        forms = getattr(dp, "_sample_forms", []) or []
        sample_forms = []
        draft_records = []
        for idx, form in enumerate(forms):
            basic = _clean_widget_values({
                k: _widget_value(w) for k, w in (form.get("basic") or {}).items()
            })
            dead = _clean_widget_values({
                k: _widget_value(w) for k, w in (form.get("dead_channel") or {}).items()
            })
            cell_types = {}
            for cell_type, fields in (form.get("cell_types") or {}).items():
                cell_types[cell_type] = _clean_widget_values({
                    k: _widget_value(w) for k, w in fields.items()
                })
            sample_forms.append({
                "index": idx,
                "label": _safe(lambda f=form: f["group"].title(), f"Sample {idx + 1}"),
                "basic": basic,
                "dead_channel": dead,
                "cell_types": cell_types,
            })
            record = dict(basic)
            if "number" in dead:
                record["dead_channel"] = dead["number"]
            if "mask_path" in dead:
                record["dead_mask_path"] = dead["mask_path"]
            if cell_types:
                record["cell_types"] = cell_types
            draft_records.append(record)

        return {
            "open": bool(getattr(dp, "builder_grp", None) and dp.builder_grp.isChecked()),
            "n_samples": _safe(lambda: dp.n_samples_spin.value(), 1),
            "n_organoids": _safe(lambda: dp.n_organoid_spin.value(), 0),
            "n_immune": _safe(lambda: dp.n_immune_spin.value(), 0),
            "n_other": _safe(lambda: dp.n_other_spin.value(), 0),
            "include_dead": _safe(lambda: dp.include_dead_cb.isChecked(), False),
            "cell_types_configured": n_org + n_imm + n_oth > 0,
            "sample_forms_created": bool(forms),
            "sample_form_count": len(forms),
            "organoid_names": [e.text() for e in getattr(dp, "_organoid_name_edits", [])],
            "immune_names": [e.text() for e in getattr(dp, "_immune_name_edits", [])],
            "other_names": [e.text() for e in getattr(dp, "_other_name_edits", [])],
            "sample_forms": sample_forms,
            # These values are the live form draft and can be newer than
            # dp.metadata while an existing CSV is being edited.
            "draft_records": draft_records,
            "draft_validation": validate_metadata_records(draft_records) if draft_records else [],
            "record_source": "metadata_builder_draft" if draft_records else None,
            "save_required": bool(draft_records),
        }
    except Exception:
        return {"open": False}


def _step_readiness(main_widget, ctx: dict) -> dict:
    """Return {step: {ready, blockers}} for every workflow step."""
    md = ctx.get("metadata", {})
    md_loaded = md.get("loaded", False)
    out_set = ctx.get("output_dir_set", False)
    output_dir = ctx.get("output_dir", "")

    def _has_outputs(pattern: str) -> bool:
        """Check if any file matching glob pattern exists in output_dir."""
        import glob
        if not output_dir:
            return False
        return bool(glob.glob(pattern, recursive=True))

    seg_done = _safe(
        lambda: _has_outputs(f"{output_dir}/images/*/*_segments.zarr"), False)
    track_done = _safe(
        lambda: _has_outputs(f"{output_dir}/images/*/*_tracked.zarr"), False)
    feat_done = _safe(
        lambda: _has_outputs(
            f"{output_dir}/analysis/*/track_features/BEHAV3D_*_combined_track_features.csv"),
        False)

    steps = {
        "data_preparation": {"ready": True, "blockers": []},
        "visualization": {
            "ready": md_loaded,
            "blockers": [] if md_loaded else ["metadata not loaded"],
        },
        "segmentation": {
            "ready": md_loaded and out_set,
            "blockers": ([("metadata not loaded" if not md_loaded else None),
                          ("output directory not set" if not out_set else None)])
        },
        "tracking": {
            "ready": md_loaded and out_set and seg_done,
            "blockers": ([("metadata not loaded" if not md_loaded else None),
                          ("output directory not set" if not out_set else None),
                          ("segmentation outputs not found" if not seg_done else None)])
        },
        "feature_extraction": {
            "ready": md_loaded and out_set and track_done,
            "blockers": ([("metadata not loaded" if not md_loaded else None),
                          ("output directory not set" if not out_set else None),
                          ("tracking outputs not found" if not track_done else None)])
        },
        "filtering": {
            "ready": feat_done,
            "blockers": [] if feat_done else ["feature extraction outputs not found"],
        },
        "analysis": {
            "ready": feat_done,
            "blockers": [] if feat_done else ["feature extraction outputs not found"],
        },
    }
    # Clean up None entries in blockers lists
    for v in steps.values():
        v["blockers"] = [b for b in v["blockers"] if b]
    return steps


def _segmentation_state(main_widget) -> dict:
    """Snapshot the visible method and per-cell instance strategies."""
    seg = _safe(lambda: getattr(main_widget, "segmentation_tab", None))
    combo = _safe(lambda: getattr(seg, "method_combo", None))
    if combo is None:
        return {}
    try:
        methods = [combo.itemText(i) for i in range(combo.count())]
    except Exception:
        methods = []
    state = {
        "method": _safe(lambda: combo.currentText(), ""),
        "method_index": _safe(lambda: combo.currentIndex(), 0),
        "available_methods": methods,
    }
    for token, page_attr in (("apoc", "apoc_page"), ("convpaint", "convpaint_page")):
        page = getattr(seg, page_attr, None)
        training = getattr(page, "_training_widget", None) if page is not None else None
        if training is None:
            continue
        global_combo = getattr(training, "strategy_combo", None)
        global_strategy = (
            _safe(lambda c=global_combo: c.currentText(), "") if global_combo is not None else ""
        )
        by_cell_type = {}
        for cell_type, tab in (getattr(training, "tabs", {}) or {}).items():
            per_tab = getattr(tab, "_per_tab_strategy_combo", None)
            by_cell_type[str(cell_type)] = (
                _safe(lambda c=per_tab: c.currentText(), "")
                if per_tab is not None else global_strategy
            )
        state[token] = {
            "global_strategy": global_strategy or None,
            "cell_type_strategies": by_cell_type,
        }
    sam = getattr(seg, "cellpose_sam_page", None)
    if sam is not None:
        all_types = bool(_widget_value(getattr(sam, "check_all_cell_types", None)))
        selected_type = _widget_value(getattr(sam, "cell_type_combo", None))
        state["cellpose_sam"] = {
            "selected_cell_type": selected_type,
            "run_all_cell_types": all_types,
            "process_all_timepoints": _widget_value(
                getattr(sam, "check_process_all", None)
            ),
            "use_3d": _widget_value(getattr(sam, "check_do_3d", None)),
            "force_cpu": _widget_value(getattr(sam, "btn_force_cpu", None)),
            "environment_status": _widget_value(
                getattr(sam, "lbl_env_status", None)
            ),
        }
    return state


def _feature_extraction_state(main_widget) -> dict:
    tab = getattr(main_widget, "feature_extraction_tab", None)
    toggle = getattr(tab, "_ak_toggle_btn", None) if tab is not None else None
    expanded = bool(_widget_value(toggle)) if toggle is not None else False
    active = getattr(tab, "active_killing_panel", None) if tab is not None else None
    if active is None:
        return {"active_killing_open": False}
    targets = []
    target_list = getattr(active, "target_list", None)
    if target_list is not None:
        targets = _safe(
            lambda: [str(item.text()) for item in target_list.selectedItems()],
            [],
        ) or []
    return {
        "active_killing_open": expanded,
        "active_killing": {
            "effector_cell_type": _widget_value(getattr(active, "immune_combo", None)),
            "target_cell_types": targets,
        },
    }


def _analysis_state(main_widget) -> dict:
    analysis = getattr(main_widget, "analysis_tab", None)
    single = getattr(analysis, "single_cell_tab", None) if analysis is not None else None
    if single is None:
        return {}
    outer_tabs = getattr(analysis, "inner_tabs", None)
    outer_index = _safe(outer_tabs.currentIndex, 1) if outer_tabs is not None else 1
    if outer_index != 1:
        return {"view": "death_dynamics"}
    stack = getattr(single, "_stack", None)
    if stack is not None and _safe(stack.currentIndex, 0) == 0:
        view = "single_cell_overview"
    else:
        inner_tabs = getattr(single, "inner_tabs", None)
        inner_index = (
            _safe(inner_tabs.currentIndex, 0)
            if inner_tabs is not None else 0
        )
        view = "behavioral_state" if inner_index == 0 else "state_trajectory"
    return {
        "view": view,
        "selected_cell_type": _widget_value(
            getattr(single, "cell_type_combo", None)
        ),
    }


def _required_params_at_default(
        step: str, step_schema: list, params: dict) -> list:
    """Return up to 6 schema cards for the current step whose value equals the default."""
    from behav3d.napari._assistant_actions import get_by_dotted
    result = []
    for card in (step_schema or []):
        if card.get("type") == "none":
            continue
        key = card.get("key", "")
        default = card.get("default")
        current = get_by_dotted(params, key, None)
        if current is None or current == default:
            result.append({
                "key": key,
                "default": default,
                "description": card.get("description", ""),
            })
        if len(result) >= 6:
            break
    return result


def build_context(main_widget) -> dict:
    """Build the full workflow context dict from a live ``BEHAV3DWidget``."""
    from behav3d.napari._assistant_schema import cards_for_step

    tabs = getattr(main_widget, "tabs", None)
    tab_index = _safe(lambda: tabs.currentIndex(), 0) if tabs is not None else 0
    tab_label = _safe(lambda: tabs.tabText(tab_index), "") if tabs is not None else ""
    step = _TAB_INDEX_TO_STEP.get(tab_index, "general")

    dp = getattr(main_widget, "data_prep_tab", None)
    metadata = getattr(dp, "metadata", None) if dp is not None else None
    output_dir = getattr(dp, "output_dir", "") if dp is not None else ""
    params = getattr(dp, "behav3d_parameters", {}) if dp is not None else {}

    queue_panel = getattr(main_widget, "queue_panel", None)
    queue = serialize_queue(queue_panel) if queue_panel is not None else []

    metadata_summary = summarize_metadata(metadata)
    builder_state = _safe(lambda: _metadata_builder_state(dp), {"open": False}) or {
        "open": False
    }
    draft_records = builder_state.get("draft_records") or []
    if draft_records:
        # Questions about current form entries must use the draft. Otherwise an
        # older loaded DataFrame can make the assistant repeat a completed change.
        metadata_summary["records"] = draft_records
        metadata_summary["validation"] = builder_state.get("draft_validation", [])
        metadata_summary["n_samples"] = len(draft_records)
        metadata_summary["sample_names"] = [
            str(record.get("sample_name")) for record in draft_records
            if record.get("sample_name") not in (None, "")
        ]
        metadata_summary["record_source"] = "metadata_builder_draft"
        metadata_summary["save_required"] = True

    ctx = {
        "current_step": step,
        "current_tab_label": tab_label,
        "tab_index": tab_index,
        "output_dir": output_dir,
        "output_dir_set": bool(output_dir),
        "metadata": metadata_summary,
        "queue": queue,
        "results": _safe(lambda: _serialize_results(output_dir), []) or [],
        "active_preview": _safe(lambda: _active_preview(main_widget), None),
        "parameters": _safe(lambda: _diff_from_defaults(params), {}) or {},
        # The Visualization tab is a viewer with no tunable BEHAV3D parameters that
        # have visible widgets (its schema cards — use_range/start_t/end_t/
        # channel_colors — are phantom controls). Don't feed them to the model, or
        # "Explain this screen" describes parameters the user can't see.
        "step_schema": [] if step == "visualization"
                       else (_safe(lambda: cards_for_step(step), []) or []),
    }
    experiment_reference = _safe(
        lambda: _experiment_reference_context(output_dir, params), None
    )
    if experiment_reference:
        ctx["experiment_reference"] = experiment_reference
    # Keep an open/drafted builder visible even after a tab switch. This avoids
    # regressing to the stale saved DataFrame on the next assistant turn.
    if (step == "data_preparation" or builder_state.get("open")
            or builder_state.get("sample_forms_created")):
        ctx["metadata_builder"] = builder_state
    if step == "segmentation":
        ctx["segmentation"] = _safe(lambda: _segmentation_state(main_widget), {})
    if step == "feature_extraction":
        ctx["feature_extraction"] = _safe(
            lambda: _feature_extraction_state(main_widget), {}
        )
    if step == "analysis":
        ctx["analysis"] = _safe(lambda: _analysis_state(main_widget), {})

    controls = _safe(lambda: control_registry(main_widget), []) or []
    ctx["ui_state"] = {
        "contract_version": CONTROL_CONTRACT_VERSION,
        "controls": controls,
    }

    # Guided-flow enrichments.
    ctx["step_readiness"] = _safe(lambda: _step_readiness(main_widget, ctx), {})
    ctx["active_cell_type"] = _safe(lambda: active_cell_type(main_widget, step), None)
    ctx["required_params_at_default"] = _safe(
        lambda: _required_params_at_default(step, ctx["step_schema"], params), [])
    return ctx


def context_summary_line(ctx: dict) -> str:
    """One-line human summary for the dock's context bar."""
    md = ctx.get("metadata", {})
    n = md.get("n_samples", 0) if md.get("loaded") else 0
    out_ok = "output set" if ctx.get("output_dir_set") else "no output dir"
    nq = len(ctx.get("queue", []))
    step = ctx.get("current_tab_label") or ctx.get("current_step", "")
    details = [str(step), f"{n} sample(s)", out_ok, f"{nq} queued"]
    if ctx.get("experiment_reference"):
        details.append("experiment context")
    active = ctx.get("active_cell_type")
    if active:
        details.insert(1, str(active))
    method = (ctx.get("segmentation") or {}).get("method")
    if method:
        details.insert(1, str(method).split("(")[0].strip())
    return "  |  ".join(details)
