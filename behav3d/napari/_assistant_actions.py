"""
BEHAV3D assistant — UI control registry & proposed-action application.

The model never writes to the UI directly. It emits structured **actions**
(``set_parameter``, ``navigate_to_step``, ``add_queue_step``); this module
validates them against the parameter schema, renders a human-readable
*old → new* preview, and — only when the user confirms — applies them.

Apply strategy (robust by design):
  1. Write the value into the live ``behav3d_parameters`` dict (the canonical
     store) via :func:`set_by_dotted`.
  2. Persist to ``behav3d_parameters.yml`` (reusing the tab's own saver when
     available, else a direct yaml dump).
  3. Best-effort push params → widgets via the owning tab's refresh hook and
     briefly pulse the target widget so the user sees where the value landed.

Steps 1–2 are always correct; step 3 is a UX nicety registered per-key where the
widget mapping has been verified. This avoids hard-coding fragile widget
attribute names for all ~195 parameters.
"""
from __future__ import annotations

from typing import Any, Callable, Optional


class _Missing:
    pass


_MISSING = _Missing()


# ---------------------------------------------------------------------------
# Pure dotted-key helpers (no Qt, fully unit-testable)
# ---------------------------------------------------------------------------
def get_by_dotted(d: dict, dotted: str, default: Any = _MISSING) -> Any:
    cur: Any = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def set_by_dotted(d: dict, dotted: str, value: Any) -> None:
    """Set ``d[a][b][c] = value``, creating intermediate dicts as needed."""
    parts = dotted.split(".")
    cur = d
    for part in parts[:-1]:
        nxt = cur.get(part)
        if not isinstance(nxt, dict):
            nxt = {}
            cur[part] = nxt
        cur = nxt
    cur[parts[-1]] = value


# ---------------------------------------------------------------------------
# Validation against schema cards
# ---------------------------------------------------------------------------
def _coerce(card: dict, value: Any) -> Any:
    """Coerce a JSON value (often str/float from the model) to the card's type."""
    typ = card.get("type")
    if typ == "int":
        return int(value)
    if typ == "float":
        return float(value)
    if typ == "bool":
        if isinstance(value, str):
            return value.strip().lower() in ("1", "true", "yes", "on")
        return bool(value)
    if typ == "str":
        return str(value)
    return value


def validate_value(card: dict, value: Any) -> tuple[bool, Any, str]:
    """Return ``(ok, coerced_value, message)``.

    Checks type coercion, enumerated choices, and numeric bounds (min/max are
    optional keys on the card)."""
    try:
        coerced = _coerce(card, value)
    except (ValueError, TypeError):
        return False, value, f"'{value}' is not a valid {card.get('type')} for {card['key']}."

    choices = card.get("choices")
    if choices is not None and coerced not in choices:
        return False, coerced, f"{card['key']} must be one of {choices}; got {coerced!r}."

    lo, hi = card.get("min"), card.get("max")
    if lo is not None and coerced < lo:
        return False, coerced, f"{card['key']} must be ≥ {lo}; got {coerced}."
    if hi is not None and coerced > hi:
        return False, coerced, f"{card['key']} must be ≤ {hi}; got {coerced}."

    return True, coerced, ""


# ---------------------------------------------------------------------------
# Widget-sync registry: dotted key -> applier(main_widget, value) -> widget|None
# The applier performs the actual UI update (so "Fill it in ✨" visibly changes
# the form, not just the in-memory params) and returns the widget to pulse.
# Only verified entries belong here; unmapped keys still update behav3d_parameters
# + the YAML, just without a live widget change.
# ---------------------------------------------------------------------------
WidgetApplier = Callable[[Any, Any], Any]
_WIDGET_APPLIERS: dict[str, WidgetApplier] = {}


def register_applier(dotted: str, applier: WidgetApplier) -> None:
    _WIDGET_APPLIERS[dotted] = applier


def _dp(main_widget):
    return getattr(main_widget, "data_prep_tab", None)


def _apply_metadata_csv(mw, value):
    dp = _dp(mw)
    w = getattr(dp, "csv_path_edit", None) if dp else None
    if w is not None:
        w.setText(str(value))   # populates the Metadata Loader field (user clicks Load)
    return w


def _apply_output_dir(mw, value):
    dp = _dp(mw)
    w = getattr(dp, "output_dir_edit", None) if dp else None
    if w is not None:
        w.setText(str(value))
    if dp is not None:          # programmatic setText doesn't fire editingFinished
        try:
            dp.output_dir = str(value)
        except Exception:
            pass
    return w


def _apply_dim_order(mw, value):
    dp = _dp(mw)
    w = getattr(dp, "dim_apply_all_combo", None) if dp else None
    if w is None:
        return None
    # The combo's options depend on the loaded image's dimensionality, so only
    # set it when the requested order is actually available right now.
    idx = w.findText(str(value))
    if idx < 0:
        return None
    w.setCurrentIndex(idx)
    return w


# Data Preparation tab fields the assistant can fill in visibly.
register_applier("paths.metadata_csv", _apply_metadata_csv)
register_applier("paths.output_dir", _apply_output_dir)
register_applier("dim_order.default_apply_all", _apply_dim_order)


def _push_to_widget(widget, value) -> bool:
    """Best-effort set of a Qt widget's value. Returns True on success."""
    try:
        if hasattr(widget, "setChecked") and isinstance(value, bool):
            widget.setChecked(value)
        elif hasattr(widget, "setValue") and isinstance(value, (int, float)):
            widget.setValue(value)
        elif hasattr(widget, "setCurrentText"):
            widget.setCurrentText(str(value))
        elif hasattr(widget, "setText"):
            widget.setText(str(value))
        else:
            return False
        return True
    except Exception:
        return False


def pulse_widget(widget, msec: int = 1200) -> None:
    """Briefly highlight a widget with a transient stylesheet glow."""
    try:
        from qtpy.QtCore import QTimer
        original = widget.styleSheet()
        widget.setStyleSheet(original + "\nborder: 2px solid #ffd166; border-radius: 4px;")
        QTimer.singleShot(msec, lambda: widget.setStyleSheet(original))
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Action model & application
# ---------------------------------------------------------------------------
class ProposedAction:
    """A single proposed change, ready to preview and (on confirm) apply."""

    def __init__(self, kind: str, **kwargs):
        self.kind = kind          # "set_parameter" | "navigate_to_step" | "add_queue_step"
        self.data = kwargs
        self.ok = True
        self.message = ""
        self.preview = ""


def _card_index(cards: list[dict]) -> dict[str, dict]:
    return {c["key"]: c for c in cards}


def build_actions(raw_tool_calls: list[dict], cards: list[dict], params: dict) -> list[ProposedAction]:
    """Turn raw model tool-calls into validated ProposedAction objects with previews."""
    idx = _card_index(cards)
    actions: list[ProposedAction] = []
    for call in raw_tool_calls or []:
        name = call.get("name")
        args = call.get("arguments", {}) or {}
        if name == "set_parameter":
            key = args.get("key")
            card = idx.get(key)
            act = ProposedAction("set_parameter", key=key, value=args.get("value"))
            if card is None:
                act.ok = False
                act.message = f"Unknown parameter '{key}'."
            else:
                ok, coerced, msg = validate_value(card, args.get("value"))
                act.ok, act.message = ok, msg
                act.data["value"] = coerced
                old = get_by_dotted(params, key, card.get("default"))
                act.preview = f"{key}:  {old!r}  →  {coerced!r}"
                if ok and coerced == old:
                    # No-op: proposed value already equals the current value
                    # (e.g. an empty placeholder '' -> ''). Don't offer to apply.
                    act.ok = False
                    act.message = f"Already set to {coerced!r} — no change needed."
                elif ok and isinstance(coerced, str) and coerced.strip() == "":
                    # Empty string proposal for a value the model doesn't know yet.
                    act.ok = False
                    act.message = "No value provided — set this yourself or tell me the value."
            actions.append(act)
        elif name == "navigate_to_step":
            step = args.get("step")
            act = ProposedAction("navigate_to_step", step=step)
            act.preview = f"Go to the {step} tab"
            actions.append(act)
        elif name == "add_queue_step":
            act = ProposedAction("add_queue_step", step_type=args.get("step_type"),
                                 params=args.get("params", {}))
            act.preview = f"Add '{args.get('step_type')}' to the processing queue"
            actions.append(act)
    return actions


def apply_set_parameter(main_widget, key: str, value: Any):
    """Apply a validated set_parameter. Returns (stored_ok, visible_field_updated)."""
    dp = getattr(main_widget, "data_prep_tab", None)
    if dp is None:
        return False, False
    params = getattr(dp, "behav3d_parameters", None)
    if params is None:
        return False, False

    # 1. canonical store
    set_by_dotted(params, key, value)

    # 2. persist to behav3d_parameters.yml (best-effort)
    _persist_params(dp)

    # 3. live UI update + pulse (for keys with a verified applier)
    widget = None
    applier = _WIDGET_APPLIERS.get(key)
    if applier is not None:
        try:
            widget = applier(main_widget, value)
        except Exception:
            widget = None
        if widget is not None:
            pulse_widget(widget)

    # (stored_ok, visible_field_updated)
    return True, (widget is not None)


def _persist_params(dp) -> None:
    import os
    out_dir = getattr(dp, "output_dir", "") or ""
    params = getattr(dp, "behav3d_parameters", None)
    if not out_dir or params is None:
        return
    try:
        import yaml
        path = os.path.join(out_dir, "behav3d_parameters.yml")
        with open(path, "w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)
    except Exception:
        pass


def apply_navigate(main_widget, step: str) -> bool:
    from behav3d.napari._assistant_context import _TAB_INDEX_TO_STEP
    step_to_index = {v: k for k, v in _TAB_INDEX_TO_STEP.items()}
    idx = step_to_index.get(step)
    tabs = getattr(main_widget, "tabs", None)
    if idx is None or tabs is None:
        return False
    try:
        tabs.setCurrentIndex(idx)
        return True
    except Exception:
        return False


def apply_action(main_widget, action: ProposedAction) -> bool:
    if not action.ok:
        return False
    if action.kind == "set_parameter":
        stored, visible = apply_set_parameter(
            main_widget, action.data["key"], action.data["value"])
        action.data["widget_updated"] = visible   # dock uses this for an honest message
        return stored
    if action.kind == "navigate_to_step":
        return apply_navigate(main_widget, action.data["step"])
    if action.kind == "add_queue_step":
        return _apply_add_queue_step(main_widget, action)
    return False


def _apply_add_queue_step(main_widget, action: ProposedAction) -> bool:
    qp = getattr(main_widget, "queue_panel", None)
    if qp is None:
        return False
    try:
        from behav3d.napari._queue import StepType
        st = StepType(action.data["step_type"])
        qp.add_step(st, params=action.data.get("params", {}))
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Tool schema advertised to the model (kept here so client & server agree)
# ---------------------------------------------------------------------------
TOOL_SCHEMA = [
    {
        "name": "set_parameter",
        "description": "Propose setting a BEHAV3D parameter. Requires user confirmation. "
                       "Use dotted keys from the provided step_schema (e.g. "
                       "'tracking.immune.trackpy.search_range_px').",
        "parameters": {
            "type": "object",
            "properties": {
                "key": {"type": "string"},
                "value": {},
            },
            "required": ["key", "value"],
        },
    },
    {
        "name": "navigate_to_step",
        "description": "Propose switching the active tab to a workflow step.",
        "parameters": {
            "type": "object",
            "properties": {
                "step": {"type": "string",
                         "enum": ["data_preparation", "visualization", "segmentation",
                                  "tracking", "feature_extraction", "filtering", "analysis"]},
            },
            "required": ["step"],
        },
    },
    {
        "name": "add_queue_step",
        "description": "Propose adding a processing step to the queue.",
        "parameters": {
            "type": "object",
            "properties": {
                "step_type": {"type": "string"},
                "params": {"type": "object"},
            },
            "required": ["step_type"],
        },
    },
]
