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


# ---------------------------------------------------------------------------
# Group appliers: handle a whole dotted-key prefix. Per-cell-category params live
# in per-cell-type panels (tab.panels[cell_type]); the config is keyed by category
# (immune/organoid/other), so we fan a category write out to every cell-type panel
# in that category. Binding is one-way at build time, so we set widgets directly.
# ---------------------------------------------------------------------------
_GROUP_APPLIERS: list = []   # (prefix, fn(main_widget, dotted_key, value) -> widget|None)


def register_group_applier(prefix: str, fn: Callable) -> None:
    _GROUP_APPLIERS.append((prefix, fn))


def _category_of(main_widget, cell_type: str) -> Optional[str]:
    dp = _dp(main_widget)
    md = getattr(dp, "metadata", None) if dp else None
    if md is None:
        return None
    try:
        from behav3d.widgets.utils import detect_cell_type_category
        return detect_cell_type_category(cell_type, md)
    except Exception:
        return None


def _set_value(widget, value) -> bool:
    return _push_to_widget(widget, value) if widget is not None else False


def _category_panels(main_widget, tab_attr: str, category: str):
    """Yield (cell_type, panel) for every panel in `category`, force-building the
    tab's per-cell-type panels first if they haven't been created yet (they are
    built lazily on tab-show / metadata_updated)."""
    tab = getattr(main_widget, tab_attr, None)
    if tab is None:
        return []
    panels = getattr(tab, "panels", {}) or {}
    if not panels and hasattr(tab, "_on_metadata_updated"):
        if getattr(_dp(main_widget), "metadata", None) is not None:
            try:
                tab._on_metadata_updated()
            except Exception:
                pass
            panels = getattr(tab, "panels", {}) or {}
    return [(ct, p) for ct, p in panels.items() if _category_of(main_widget, ct) == category]


def _fan_to_panels(main_widget, tab_attr: str, category: str, widget_attr: str, value):
    """Set `panel.<widget_attr> = value` across all panels in `category`."""
    last = None
    for _ct, panel in _category_panels(main_widget, tab_attr, category):
        w = getattr(panel, widget_attr, None)
        if _set_value(w, value):
            last = w
    return last


def _fan_checkset(main_widget, tab_attr, category, dict_attr, selected):
    """For a dict-of-checkboxes (e.g. feature_checks, bt_hyp_checks): check the
    boxes whose key is in `selected`, uncheck the rest."""
    chosen = set(selected) if isinstance(selected, (list, tuple, set)) else set()
    last = None
    for _ct, panel in _category_panels(main_widget, tab_attr, category):
        checks = getattr(panel, dict_attr, None) or {}
        for key, w in checks.items():
            try:
                w.setChecked(key in chosen)
                last = w
            except Exception:
                pass
    return last


# --- Filtering tab (track_filtering.<category>.<field>) --------------------
_FILTER_FIELD = {
    "min_track_length": "spin_min_length",
    "min_track_length_enabled": "en_min_length",
    "max_track_length": "spin_max_length",
    "max_track_length_enabled": "en_max_length",
    "exp_duration": "spin_exp_duration",
    "exp_duration_enabled": "en_exp_duration",
}


def _apply_filtering(main_widget, dotted: str, value):
    parts = dotted.split(".")           # track_filtering.<cat>.<field>
    if len(parts) != 3:
        return None
    _, category, field = parts
    widget_attr = _FILTER_FIELD.get(field)
    if widget_attr is None:
        return None
    return _fan_to_panels(main_widget, "filtering_tab", category, widget_attr, value)


register_group_applier("track_filtering.", _apply_filtering)


# --- Tracking tab (tracking.<category>.{method | lap.* | trackpy.* | btrack.*}) ---
_TRACK_FIELD = {
    "method": "combo_method",
    "lap.track_cost_px": "lap_track_cost",
    "lap.gap_close_cost_px": "lap_gap_cost",
    "lap.gap_close_max_frames": "lap_gap_frames",
    "lap.merging_cost_px": "lap_merge_cost",
    "lap.splitting_cost_px": "lap_split_cost",
    "trackpy.search_range_px": "tp_search_range",
    "trackpy.memory_frames": "tp_memory",
    "trackpy.adaptive_stop": "tp_adaptive_stop",
    "trackpy.adaptive_step": "tp_adaptive_step",
    "btrack.config_preset": "bt_config_preset",
    "btrack.config_path": "bt_config_path",
    "btrack.max_search_radius": "bt_max_search_radius",
    "btrack.update_method": "bt_update_method",
    "btrack.step_size": "bt_step_size",
    "btrack.n_workers": "bt_n_workers",
    "btrack.use_optimize": "bt_use_optimize",
    "btrack.dist_thresh": "bt_dist_thresh",
    "btrack.time_thresh": "bt_time_thresh",
}


def _apply_tracking(main_widget, dotted: str, value):
    parts = dotted.split(".")                 # tracking.<cat>.<rest...>
    if len(parts) < 3:
        return None
    category, rest = parts[1], ".".join(parts[2:])
    if rest == "btrack.hypotheses":
        return _fan_checkset(main_widget, "tracking_tab", category, "bt_hyp_checks", value)
    widget_attr = _TRACK_FIELD.get(rest)
    if widget_attr is None:
        return None
    return _fan_to_panels(main_widget, "tracking_tab", category, widget_attr, value)


register_group_applier("tracking.", _apply_tracking)


# --- Feature Extraction tab (features.<category>.<field>) ------------------
_FEATURE_FIELD = {"contact_threshold": "contact_threshold", "n_workers": "spin_workers"}


def _apply_features(main_widget, dotted: str, value):
    parts = dotted.split(".")                 # features.<cat>.<field>
    if len(parts) != 3:
        return None
    category, field = parts[1], parts[2]
    if field == "features_choice":
        return _fan_checkset(main_widget, "feature_extraction_tab", category, "feature_checks", value)
    widget_attr = _FEATURE_FIELD.get(field)
    if widget_attr is None:
        return None
    return _fan_to_panels(main_widget, "feature_extraction_tab", category, widget_attr, value)


register_group_applier("features.", _apply_features)


# --- Death dynamics (death_dynamics.<category>.dead_perc_threshold) ---------
def _apply_death(main_widget, dotted: str, value):
    parts = dotted.split(".")
    if len(parts) != 3 or parts[2] != "dead_perc_threshold":
        return None
    return _fan_to_panels(main_widget, "feature_extraction_tab", parts[1],
                          "spin_dead_threshold", value)


register_group_applier("death_dynamics.", _apply_death)


def _push_to_widget(widget, value) -> bool:
    """Best-effort set of a Qt widget's value. Returns True on success."""
    try:
        # Combos: items are often labelled (e.g. "btrack (Bayesian)") while the
        # config value is the bare token ("btrack"). Match exact → starts-with →
        # contains so labelled options still resolve; fail if none match.
        from qtpy.QtWidgets import QComboBox
        if isinstance(widget, QComboBox):
            from qtpy.QtCore import Qt
            s = str(value)
            for flag in (Qt.MatchFixedString, Qt.MatchStartsWith, Qt.MatchContains):
                idx = widget.findText(s, flag)
                if idx >= 0:
                    widget.setCurrentIndex(idx)
                    return True
            return False
        if hasattr(widget, "setChecked") and isinstance(value, bool):
            widget.setChecked(value)
        elif hasattr(widget, "setValue") and isinstance(value, (int, float)):
            widget.setValue(value)
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

    # 3. live UI update + pulse — exact-key applier first, then group (prefix) appliers
    widget = None
    applier = _WIDGET_APPLIERS.get(key)
    if applier is not None:
        try:
            widget = applier(main_widget, value)
        except Exception:
            widget = None
    if widget is None:
        for prefix, fn in _GROUP_APPLIERS:
            if key.startswith(prefix):
                try:
                    widget = fn(main_widget, key, value)
                except Exception:
                    widget = None
                break
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
