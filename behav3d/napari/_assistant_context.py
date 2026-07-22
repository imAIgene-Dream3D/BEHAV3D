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
      "queue": [{"type": "track", "label": "📍 Batch Tracking", "status": "pending", "params": {...}}],
      "parameters": {                      # only values that differ from defaults
          "tracking.immune.method": {"current": "btrack", "default": "trackpy"},
          ...
      },
      "step_schema": [ <cards for current_step> ],
    }
"""
from __future__ import annotations

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
    return state


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
    # Keep an open/drafted builder visible even after a tab switch. This avoids
    # regressing to the stale saved DataFrame on the next assistant turn.
    if (step == "data_preparation" or builder_state.get("open")
            or builder_state.get("sample_forms_created")):
        ctx["metadata_builder"] = builder_state
    if step == "segmentation":
        ctx["segmentation"] = _safe(lambda: _segmentation_state(main_widget), {})

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
    active = ctx.get("active_cell_type")
    if active:
        details.insert(1, str(active))
    method = (ctx.get("segmentation") or {}).get("method")
    if method:
        details.insert(1, str(method).split("(")[0].strip())
    return "  |  ".join(details)
