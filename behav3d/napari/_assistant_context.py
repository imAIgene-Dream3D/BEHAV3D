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

from typing import Any, Optional

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


def summarize_metadata(metadata) -> dict:
    """Summarise a metadata DataFrame without assuming any column exists."""
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

    return {
        "loaded": True,
        "n_samples": n_samples,
        "sample_names": sample_names,
        "columns": columns,
        "cell_types": cell_types,
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
        for idx, form in enumerate(forms[:5]):
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


def _active_cell_type_tab(main_widget, step: str) -> str | None:
    """Return 'immune', 'organoid', or 'other' for the active cell-type sub-tab."""
    tab_attr = {
        "tracking": "tracking_tab",
        "feature_extraction": "feature_extraction_tab",
        "filtering": "filtering_tab",
    }.get(step)
    if tab_attr is None:
        return None
    tab_widget = _safe(lambda: getattr(main_widget, tab_attr, None))
    if tab_widget is None:
        return None
    panels_tab = _safe(lambda: getattr(tab_widget, "panels_tab", None) or
                       getattr(tab_widget, "tab_widget", None))
    if panels_tab is None:
        return None
    idx = _safe(lambda: panels_tab.currentIndex(), 0)
    label = _safe(lambda: panels_tab.tabText(idx).lower(), "")
    for ct in ("immune", "organoid", "other"):
        if ct in label:
            return ct
    return None


def _segmentation_state(main_widget) -> dict:
    """Snapshot the visible segmentation method selector."""
    seg = _safe(lambda: getattr(main_widget, "segmentation_tab", None))
    combo = _safe(lambda: getattr(seg, "method_combo", None))
    if combo is None:
        return {}
    try:
        methods = [combo.itemText(i) for i in range(combo.count())]
    except Exception:
        methods = []
    return {
        "method": _safe(lambda: combo.currentText(), ""),
        "method_index": _safe(lambda: combo.currentIndex(), 0),
        "available_methods": methods,
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

    ctx = {
        "current_step": step,
        "current_tab_label": tab_label,
        "tab_index": tab_index,
        "output_dir": output_dir,
        "output_dir_set": bool(output_dir),
        "metadata": summarize_metadata(metadata),
        "queue": queue,
        "parameters": _safe(lambda: _diff_from_defaults(params), {}) or {},
        "step_schema": _safe(lambda: cards_for_step(step), []) or [],
    }
    # Include metadata builder state only when on the data_preparation tab.
    if step == "data_preparation":
        ctx["metadata_builder"] = _metadata_builder_state(dp)
    if step == "segmentation":
        ctx["segmentation"] = _safe(lambda: _segmentation_state(main_widget), {})

    # Guided-flow enrichments.
    ctx["step_readiness"] = _safe(lambda: _step_readiness(main_widget, ctx), {})
    ctx["active_cell_type_tab"] = _safe(
        lambda: _active_cell_type_tab(main_widget, step), None)
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
    return f"📍 {step} · {n} sample(s) · {out_ok} · {nq} queued"
