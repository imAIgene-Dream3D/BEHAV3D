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
    return ctx


def context_summary_line(ctx: dict) -> str:
    """One-line human summary for the dock's context bar."""
    md = ctx.get("metadata", {})
    n = md.get("n_samples", 0) if md.get("loaded") else 0
    out_ok = "output set" if ctx.get("output_dir_set") else "no output dir"
    nq = len(ctx.get("queue", []))
    step = ctx.get("current_tab_label") or ctx.get("current_step", "")
    return f"📍 {step} · {n} sample(s) · {out_ok} · {nq} queued"
