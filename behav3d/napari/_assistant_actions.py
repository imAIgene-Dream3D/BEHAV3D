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


_CATEGORY_LABELS = {
    "immune": "immune cells",
    "organoid": "organoids",
    "other": "other cells",
    "all_organoids": "all organoids",
}

_LEAF_LABELS = {
    "metadata_csv": "metadata CSV",
    "output_dir": "output directory",
    "default_apply_all": "dimension order",
    "examples_per_sample": "examples per sample",
    "sample_specific_classifier": "sample-specific classifier",
    "workers": "worker count",
    "use_all_timepoints": "use all timepoints",
    "tp_start": "first timepoint",
    "tp_end": "last timepoint",
    "number_of_channels": "number of channels",
    "labels_mode": "channel-label setup",
    "channel_labels": "channel labels",
    "per_sample_channel_labels": "per-sample channel labels",
    "track_organoids_together": "track organoids together",
    "method": "method",
    "overwrite": "overwrite existing outputs",
    "search_range_px": "search range",
    "memory_frames": "memory frames",
    "track_cost_px": "frame-to-frame linking distance",
    "gap_close_cost_px": "gap-closing distance",
    "gap_close_max_frames": "maximum gap length",
    "merging_cost_px": "merging distance",
    "splitting_cost_px": "splitting distance",
    "config_preset": "btrack preset",
    "config_path": "btrack config path",
    "max_search_radius": "maximum search radius",
    "update_method": "btrack update method",
    "step_size": "btrack step size",
    "n_workers": "worker count",
    "use_optimize": "global optimization",
    "use_visual_features": "use visual features",
    "hypotheses": "btrack hypotheses",
    "dist_thresh": "distance threshold",
    "time_thresh": "time threshold",
    "features_choice": "feature groups",
    "contact_threshold": "contact distance",
    "dead_perc_threshold": "dead-pixel threshold",
    "exp_duration": "experiment duration",
    "exp_duration_enabled": "use experiment-duration filter",
    "min_track_length": "minimum track length",
    "min_track_length_enabled": "use minimum-track-length filter",
    "max_track_length": "maximum track length",
    "max_track_length_enabled": "use maximum-track-length filter",
    "umap_min_dist": "UMAP minimum distance",
    "umap_n_neighbors": "UMAP neighbors",
    "nr_of_clusters": "number of clusters",
    "observation_window": "observation window",
    "death_signal_column": "death-signal column",
    "killing_threshold_multiplier": "killing threshold multiplier",
    "min_contact_duration": "minimum contact duration",
    "use_absolute_threshold": "use absolute killing threshold",
    "absolute_killing_threshold": "absolute killing threshold",
    "save_results": "save results",
}

_METHOD_GROUP_LABELS = {
    "lap": "LAP",
    "trackpy": "TrackPy",
    "btrack": "btrack",
}


def humanize_parameter_key(key: str | None) -> str:
    """Return a researcher-facing label for an internal dotted parameter key."""
    if not key:
        return "Parameter"
    parts = str(key).split(".")
    leaf = parts[-1]
    label = _LEAF_LABELS.get(leaf, leaf.replace("_", " "))

    category = next((p for p in parts if p in _CATEGORY_LABELS), None)
    method = next((p for p in parts if p in _METHOD_GROUP_LABELS), None)
    top = parts[0] if parts else ""

    prefixes: list[str] = []
    if category:
        prefixes.append(_CATEGORY_LABELS[category])
    if method and leaf != "method":
        prefixes.append(_METHOD_GROUP_LABELS[method])
    elif top == "pixel_classifier":
        prefixes.append("pixel classifier")
    elif top == "cellpose":
        prefixes.append("Cellpose")
    elif top == "features":
        prefixes.append("feature extraction")
    elif top == "track_filtering":
        prefixes.append("track filtering")
    elif top == "analysis":
        prefixes.append("analysis")
    elif top == "active_killing":
        prefixes.append("active killing")
    elif top == "death_dynamics":
        prefixes.append("death dynamics")

    return ": ".join(prefixes + [label]) if prefixes else label.capitalize()


def _display_value(value: Any) -> str:
    if isinstance(value, str):
        return f'"{value}"'
    return repr(value)


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
        label = humanize_parameter_key(card.get("key"))
        return False, value, f"'{value}' is not a valid {card.get('type')} for {label}."

    choices = card.get("choices")
    if choices is not None and coerced not in choices:
        label = humanize_parameter_key(card.get("key"))
        return False, coerced, f"{label} must be one of {choices}; got {coerced!r}."

    lo, hi = card.get("min"), card.get("max")
    if lo is not None and coerced < lo:
        label = humanize_parameter_key(card.get("key"))
        return False, coerced, f"{label} must be at least {lo}; got {coerced}."
    if hi is not None and coerced > hi:
        label = humanize_parameter_key(card.get("key"))
        return False, coerced, f"{label} must be at most {hi}; got {coerced}."

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


# --- Segmentation tab pages (global, non-cell-type widgets) ----------------
def _seg_page(main_widget, page_attr):
    seg = getattr(main_widget, "segmentation_tab", None)
    return getattr(seg, page_attr, None) if seg is not None else None


_PIXCLF_FIELD = {
    "examples_per_sample": "spin_examples",
    "workers": "spin_workers",
    "use_all_timepoints": "check_process_all",
    "tp_start": "spin_t_start",
    "tp_end": "spin_t_end",
}


def _apply_pixel_classifier(main_widget, dotted: str, value):
    parts = dotted.split(".")                  # pixel_classifier.<field>
    if len(parts) != 2:
        return None
    attr = _PIXCLF_FIELD.get(parts[1])
    page = _seg_page(main_widget, "pixel_classifier_page")
    w = getattr(page, attr, None) if (page and attr) else None
    return w if _set_value(w, value) else None


register_group_applier("pixel_classifier.", _apply_pixel_classifier)


_CELLPOSE_FIELD = {"number_of_channels": "spin_manual_channels"}


def _apply_cellpose(main_widget, dotted: str, value):
    parts = dotted.split(".")                  # cellpose.<field>
    if len(parts) != 2:
        return None
    attr = _CELLPOSE_FIELD.get(parts[1])
    page = _seg_page(main_widget, "cellpose_page")
    w = getattr(page, attr, None) if (page and attr) else None
    return w if _set_value(w, value) else None


register_group_applier("cellpose.", _apply_cellpose)


def _push_to_widget(widget, value) -> bool:
    """Best-effort set of a Qt widget's value. Returns True on success."""
    try:
        # Combos: items are often labelled (e.g. "btrack (Bayesian)") while the
        # config value is the bare token ("btrack"). Match exact → starts-with →
        # contains so labelled options still resolve; fail if none match.
        from qtpy.QtWidgets import QCheckBox, QComboBox, QDoubleSpinBox, QSpinBox
        if isinstance(widget, QComboBox):
            from qtpy.QtCore import Qt
            s = str(value)
            for flag in (Qt.MatchFixedString, Qt.MatchStartsWith, Qt.MatchContains):
                idx = widget.findText(s, flag)
                if idx >= 0:
                    widget.setCurrentIndex(idx)
                    return True
            return False
        if isinstance(widget, QCheckBox):
            widget.setChecked(_coerce_bool(value))
        elif isinstance(widget, QSpinBox):
            widget.setValue(int(value))
        elif isinstance(widget, QDoubleSpinBox):
            widget.setValue(float(value))
        elif hasattr(widget, "setValue") and isinstance(value, (int, float)):
            widget.setValue(value)
        elif hasattr(widget, "setText"):
            widget.setText(str(value))
        else:
            return False
        return True
    except Exception:
        return False


# When True, pulse_widget() is a no-op. Set during batch/bulk applies so a flood
# of fills doesn't schedule hundreds of QTimers + stylesheet thrash (which can
# freeze the GUI). See AssistantDock._on_tool_calls and apply_bulk_fill_metadata.
_SUPPRESS_PULSE = False


def set_pulse_suppressed(flag: bool) -> None:
    global _SUPPRESS_PULSE
    _SUPPRESS_PULSE = bool(flag)


def pulse_widget(widget, msec: int = 1200) -> None:
    """Briefly highlight a widget with a transient stylesheet glow."""
    if _SUPPRESS_PULSE:
        return
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


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def normalize_metadata_line_value(value: Any) -> Any:
    """Use the CSV-safe sentinel for a configured population that is absent."""
    if value is None:
        return "not_added"
    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    if normalized in {
        "none", "null", "n/a", "na", "absent", "not_added", "(not_added)",
    }:
        return "not_added"
    return value


def _coerce_control_value(control: dict, value: Any) -> tuple[bool, Any, str]:
    current = control.get("value")
    try:
        if isinstance(current, bool):
            coerced = _coerce_bool(value)
        elif isinstance(current, int):
            coerced = int(value)
        elif isinstance(current, float):
            coerced = float(value)
        elif isinstance(current, list):
            if not isinstance(value, (list, tuple, set)):
                return False, value, "Choose one or more values from the available options."
            coerced = list(value)
        else:
            coerced = str(value)
    except (TypeError, ValueError):
        return False, value, f"{value!r} is not valid for {control.get('label', 'this control')}."
    choices = control.get("choices")
    if choices:
        requested = coerced if isinstance(coerced, list) else [coerced]
        resolved = []
        for item in requested:
            text = str(item).lower()
            match = next((choice for choice in choices
                          if str(choice).lower() == text
                          or str(choice).lower().startswith(text)), None)
            if match is None:
                return False, coerced, f"Choose a value shown in {control.get('label', 'this control')}."
            resolved.append(match)
        coerced = resolved if isinstance(coerced, list) else resolved[0]
    required_choices = control.get("required_choices") or []
    if isinstance(coerced, list) and not set(required_choices).issubset(set(coerced)):
        required_labels = ", ".join(str(item) for item in required_choices)
        return (
            False,
            coerced,
            f"{control.get('label', 'This selection')} must keep: {required_labels}.",
        )
    minimum, maximum = control.get("minimum"), control.get("maximum")
    if minimum is not None and isinstance(coerced, (int, float)) and coerced < minimum:
        return False, coerced, f"{control.get('label')} must be at least {minimum}."
    if maximum is not None and isinstance(coerced, (int, float)) and coerced > maximum:
        return False, coerced, f"{control.get('label')} must be at most {maximum}."
    return True, coerced, ""


def build_actions(
    raw_tool_calls: list[dict],
    cards: list[dict],
    params: dict,
    controls: list[dict] | None = None,
) -> list[ProposedAction]:
    """Turn raw model tool-calls into validated ProposedAction objects with previews."""
    idx = _card_index(cards)
    actions: list[ProposedAction] = []
    for call in raw_tool_calls or []:
        name = call.get("name")
        args = call.get("arguments", {}) or {}
        if name == "set_ui_value":
            control_id = str(args.get("control_id") or "")
            control = next((c for c in (controls or []) if c.get("id") == control_id), None)
            requested_value = args.get("value")
            if control_id.startswith("metadata.samples.") and control_id.endswith(".line"):
                requested_value = normalize_metadata_line_value(requested_value)
            act = ProposedAction("set_ui_value", control_id=control_id,
                                 value=requested_value)
            if control is None:
                act.ok = False
                act.message = "That field is not available in the current interface."
            elif not control.get("enabled", False):
                act.ok = False
                act.message = f"{control.get('label', 'That field')} is currently disabled."
            elif control.get("method") and not control.get("visible", False):
                act.ok = False
                act.message = "That field belongs to a different method than the one selected."
            else:
                ok, coerced, message = _coerce_control_value(control, requested_value)
                act.ok, act.message = ok, message
                act.data.update(value=coerced, label=control.get("label"),
                                old_value=control.get("value"))
                act.preview = (f"{control.get('label')}: "
                               f"{_display_value(control.get('value'))} -> {_display_value(coerced)}")
                if ok and coerced == control.get("value"):
                    act.ok = False
                    act.data["no_op"] = True
                    act.message = "This value is already set."
            actions.append(act)
        elif name == "set_parameter":
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
                act.data["label"] = humanize_parameter_key(key)
                old = get_by_dotted(params, key, card.get("default"))
                act.preview = (
                    f"{humanize_parameter_key(key)}: "
                    f"{_display_value(old)} -> {_display_value(coerced)}"
                )
                if ok and coerced == old:
                    # No-op: proposed value already equals the current value
                    # (e.g. an empty placeholder '' -> ''). Don't offer to apply.
                    act.ok = False
                    act.data["no_op"] = True
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
        elif name == "fill_metadata_builder":
            field = args.get("field")
            value = args.get("value")
            if field == "cell_line":
                value = normalize_metadata_line_value(value)
            index = _safe_int(args.get("index", 0), 0)
            act = ProposedAction(
                "fill_metadata_builder",
                field=field,
                value=value,
                index=index,
                cell_type=args.get("cell_type"),
            )
            act.preview = _metadata_builder_preview(field, value, index, args.get("cell_type"))
            if field is None:
                act.ok = False
                act.message = "No field specified for fill_metadata_builder."
            elif field in (
                "cell_line", "cell_condition", "cell_segments_image_path",
                "cell_tracks_image_path", "cell_tracks_csv_path",
            ) and not args.get("cell_type"):
                act.ok = False
                act.message = "Choose the visible cell type name for this sample field."
            actions.append(act)
        elif name == "bulk_fill_metadata":
            samples = args.get("samples", []) or []
            for sample in samples:
                if not isinstance(sample, dict):
                    continue
                for values in (sample.get("cell_types") or {}).values():
                    if isinstance(values, dict) and "line" in values:
                        values["line"] = normalize_metadata_line_value(values["line"])
            act = ProposedAction(
                "bulk_fill_metadata",
                n_samples=args.get("n_samples"),
                n_organoids=args.get("n_organoids"),
                n_immune=args.get("n_immune"),
                n_other=args.get("n_other"),
                include_dead=args.get("include_dead"),
                organoid_names=args.get("organoid_names", []) or [],
                immune_names=args.get("immune_names", []) or [],
                other_names=args.get("other_names", []) or [],
                samples=samples,
                sample_count=len(samples),
            )
            act.preview = f"Fill the Metadata Builder ({len(samples)} samples)"
            if not samples:
                act.ok = False
                act.message = "No per-sample data provided for bulk fill."
            actions.append(act)
        elif name == "select_segmentation_method":
            value = args.get("value")
            act = ProposedAction("select_segmentation_method", value=value)
            act.preview = f"Segmentation method → {value}"
            if not value:
                act.ok = False
                act.message = "No segmentation method provided."
            actions.append(act)
        elif name == "show_track_length_distribution":
            cell_type = str(args.get("cell_type") or "")
            act = ProposedAction("show_track_length_distribution", cell_type=cell_type)
            act.preview = f"Show track-length distributions for {cell_type}"
            if not cell_type:
                act.ok = False
                act.message = "Choose a cell type first."
            actions.append(act)
        elif name == "recommend_edt":
            cell_type = str(args.get("cell_type") or "").strip()
            try:
                diameter_um = float(args.get("cell_diameter_um", 10.0))
                cells_across = args.get("organoid_cells_across")
                cells_across = None if cells_across is None else float(cells_across)
            except (TypeError, ValueError):
                diameter_um, cells_across = 10.0, None
            act = ProposedAction(
                "recommend_edt", cell_type=cell_type,
                cell_diameter_um=diameter_um,
                organoid_cells_across=cells_across,
            )
            act.preview = f"Calculate EDT starting values for {cell_type}"
            if not cell_type:
                act.ok = False
                act.message = "Choose a cell type first."
            elif diameter_um <= 0 or (cells_across is not None and cells_across < 1):
                act.ok = False
                act.message = "Object size values must be greater than zero."
            actions.append(act)
        elif name == "summarize_track_counts":
            cell_type = str(args.get("cell_type") or "").strip()
            try:
                position_t = int(args.get("position_t", 0))
                raw_lengths = args.get("minimum_lengths") or [20, 50, 100, 200]
                minimum_lengths = sorted({int(value) for value in raw_lengths})
            except (TypeError, ValueError):
                position_t, minimum_lengths = 0, []
            act = ProposedAction(
                "summarize_track_counts", cell_type=cell_type,
                position_t=position_t, minimum_lengths=minimum_lengths,
            )
            act.preview = (
                f"Count {cell_type} tracks at timepoint {position_t} for minimum "
                f"lengths {minimum_lengths}"
            )
            if not cell_type:
                act.ok = False
                act.message = "Choose a cell type first."
            elif position_t < 0:
                act.ok = False
                act.message = "The timepoint must be zero or greater."
            elif not minimum_lengths or any(value < 1 for value in minimum_lengths):
                act.ok = False
                act.message = "Minimum track lengths must be positive integers."
            elif len(minimum_lengths) > 20:
                act.ok = False
                act.message = "Compare at most 20 minimum track lengths at once."
            actions.append(act)
        elif name == "create_cell_type_group":
            group_name = str(args.get("group_name") or "").strip()
            members = [str(v) for v in (args.get("members") or [])]
            act = ProposedAction("create_cell_type_group", group_name=group_name,
                                 members=members)
            act.preview = f"Create cell-type group '{group_name}' from {', '.join(members)}"
            if not group_name or not members:
                act.ok = False
                act.message = "Provide a group name and at least one cell type."
            elif not all(c.isalnum() or c in "_-" for c in group_name):
                act.ok = False
                act.message = "Group names may only contain letters, numbers, hyphens, and underscores."
            actions.append(act)
        elif name == "create_btrack_config_copy":
            cell_type = str(args.get("cell_type") or "")
            destination = str(args.get("destination") or "").strip()
            act = ProposedAction("create_btrack_config_copy", cell_type=cell_type,
                                 destination=destination)
            act.preview = f"Create a custom btrack configuration for {cell_type} at {destination}"
            if not cell_type or not destination:
                act.ok = False
                act.message = "Choose a cell type and destination file."
            elif not destination.lower().endswith(".json"):
                act.ok = False
                act.message = "The custom btrack configuration must be a JSON file."
            actions.append(act)
        elif name == "open_result":
            result_id = str(args.get("result_id") or "")
            act = ProposedAction("open_result", result_id=result_id)
            act.preview = f"Open result: {result_id}"
            if not result_id:
                act.ok = False
                act.message = "Choose a result first."
            actions.append(act)
        elif name == "save_metadata":
            act = ProposedAction("save_metadata")
            act.preview = "Save and activate the Metadata Builder draft"
            actions.append(act)
        elif name == "load_metadata":
            act = ProposedAction("load_metadata")
            act.preview = "Load the selected metadata CSV"
            actions.append(act)
        elif name == "open_analysis_view":
            view = str(args.get("view") or "")
            act = ProposedAction("open_analysis_view", view=view)
            labels = {
                "death_dynamics": "Death Dynamics",
                "behavioral_state": "Behavioral State",
                "state_trajectory": "State Trajectory",
            }
            act.preview = f"Open {labels.get(view, 'analysis')}"
            if view not in labels:
                act.ok = False
                act.message = "Choose an available analysis view."
            actions.append(act)
    return actions


def _metadata_builder_preview(
    field: Optional[str],
    value: Any,
    index: int,
    cell_type: Optional[str] = None,
) -> str:
    _TRIGGER_LABELS = {
        "open_builder": "Open the Metadata Builder",
        "configure_cell_types": "Configure cell type name fields",
        "create_sample_forms": "Create per-sample input forms",
        "fill_down": "Copy Sample 1 settings to all samples",
    }
    if field in _TRIGGER_LABELS:
        return _TRIGGER_LABELS[field]
    if field == "n_samples":      return f"Number of samples: {value}"
    if field == "n_organoids":    return f"Organoid cell types: {value}"
    if field == "n_immune":       return f"Immune cell types: {value}"
    if field == "n_other":        return f"Other cell types: {value}"
    if field == "include_dead":   return f"Include dead channel: {'yes' if _coerce_bool(value) else 'no'}"
    if field == "immune_multicolor":
        return f"Immune type {index + 1} multicolor: {'yes' if _coerce_bool(value) else 'no'}"
    if field == "immune_multicolor_channels":
        return f"Immune type {index + 1} channels: {value}"
    if field == "organoid_name":  return f"Organoid type {index + 1} name → \"{value}\""
    if field == "immune_name":    return f"Immune type {index + 1} name → \"{value}\""
    if field == "other_name":     return f"Other type {index + 1} name → \"{value}\""
    if field == "sample_name":    return f"Sample {index + 1} name → \"{value}\""
    if field == "raw_image_path": return f"Sample {index + 1} image path → \"{value}\""
    if field == "pixel_distance_xy": return f"Sample {index + 1} pixel XY: {value} μm"
    if field == "pixel_distance_z":  return f"Sample {index + 1} pixel Z: {value} μm"
    if field == "time_interval":  return f"Sample {index + 1} time interval: {value}"
    if field == "time_unit":      return f"Sample {index + 1} time unit: {value}"
    if field == "exp_nr":         return f"Sample {index + 1} experiment #: {value}"
    if field == "well":           return f"Sample {index + 1} well: {value}"
    if field == "dimension_order": return f"Sample {index + 1} dimension order: {value}"
    if field == "dead_channel_number": return f"Sample {index + 1} dead channel #: {value}"
    if field == "dead_mask_path": return f"Sample {index + 1} dead mask path → \"{value}\""
    cell = f" {cell_type}" if cell_type else ""
    if field == "cell_line": return f"Sample {index + 1}{cell} line → \"{value}\""
    if field == "cell_condition": return f"Sample {index + 1}{cell} condition → \"{value}\""
    if field == "cell_segments_image_path":
        return f"Sample {index + 1}{cell} segments path → \"{value}\""
    if field == "cell_tracks_image_path":
        return f"Sample {index + 1}{cell} tracks image path → \"{value}\""
    if field == "cell_tracks_csv_path":
        return f"Sample {index + 1}{cell} tracks CSV path → \"{value}\""
    return f"Metadata builder: {field} = {value!r}"


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


def _values_match(actual: Any, requested: Any) -> bool:
    if isinstance(actual, str) and isinstance(requested, str):
        return actual.strip().lower() == requested.strip().lower() or (
            bool(requested.strip()) and actual.strip().lower().startswith(requested.strip().lower())
        )
    if isinstance(actual, float) or isinstance(requested, float):
        try:
            return abs(float(actual) - float(requested)) < 1e-9
        except (TypeError, ValueError):
            pass
    return actual == requested


def apply_set_ui_value(main_widget, control_id: str, value: Any) -> bool:
    """Set one exact live control and verify the value was actually applied."""
    from behav3d.napari._assistant_controls import find_control

    binding = find_control(main_widget, control_id)
    if binding is None or not binding.get("enabled", False):
        return False
    setter = binding.get("_setter")
    if setter is None or not setter(value):
        return False
    persist = binding.get("_persist")
    if persist is not None:
        try:
            persist()
        except Exception:
            return False
    refreshed = find_control(main_widget, control_id)
    actual = refreshed.get("value") if refreshed is not None else None
    if not _values_match(actual, value):
        return False
    widget = binding.get("_widget")
    if widget is not None:
        pulse_widget(widget)
    dp = getattr(main_widget, "data_prep_tab", None)
    if dp is not None:
        _persist_params(dp)
    return True


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


def _coerce_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "on")
    return bool(value)


def _metadata_cell_widget(form: dict, cell_type: Optional[str], field: str):
    if not cell_type:
        return None
    cell_fields = (form.get("cell_types") or {}).get(str(cell_type))
    if not cell_fields:
        return None
    return cell_fields.get(field)


def apply_fill_metadata_builder(main_widget, action: ProposedAction) -> bool:
    """Apply a fill_metadata_builder action — sets widgets in the Metadata Builder."""
    dp = _dp(main_widget)
    if dp is None:
        return False
    field = action.data.get("field")
    value = action.data.get("value")
    index = int(action.data.get("index", 0))
    cell_type = action.data.get("cell_type")
    try:
        if field == "open_builder":
            dp.builder_grp.setChecked(True)
            # Toggling can open a modal when metadata is already loaded. If the
            # user cancels, the handler restores the unchecked state.
            return bool(dp.builder_grp.isChecked())
        # Silently open the builder if it isn't already (non-destructive).
        if not getattr(dp.builder_grp, "isChecked", lambda: True)():
            dp.builder_grp.setChecked(True)
        if field == "n_samples":
            dp.n_samples_spin.setValue(int(value)); pulse_widget(dp.n_samples_spin); return True
        if field == "n_organoids":
            dp.n_organoid_spin.setValue(int(value)); pulse_widget(dp.n_organoid_spin); return True
        if field == "n_immune":
            dp.n_immune_spin.setValue(int(value)); pulse_widget(dp.n_immune_spin); return True
        if field == "n_other":
            dp.n_other_spin.setValue(int(value)); pulse_widget(dp.n_other_spin); return True
        if field == "include_dead":
            dp.include_dead_cb.setChecked(_coerce_bool(value)); pulse_widget(dp.include_dead_cb); return True
        if field == "configure_cell_types":
            dp._on_configure_cell_types(); return True
        if field == "immune_multicolor":
            flags = getattr(dp, "_immune_multicolor_flags", [])
            if 0 <= index < len(flags):
                flags[index].setChecked(_coerce_bool(value)); pulse_widget(flags[index]); return True
        if field == "immune_multicolor_channels":
            counts = getattr(dp, "_immune_multicolor_counts", [])
            if 0 <= index < len(counts):
                counts[index].setValue(int(value)); pulse_widget(counts[index]); return True
        if field == "organoid_name":
            edits = getattr(dp, "_organoid_name_edits", [])
            if 0 <= index < len(edits):
                edits[index].setText(str(value)); pulse_widget(edits[index]); return True
        if field == "immune_name":
            edits = getattr(dp, "_immune_name_edits", [])
            if 0 <= index < len(edits):
                edits[index].setText(str(value)); pulse_widget(edits[index]); return True
        if field == "other_name":
            edits = getattr(dp, "_other_name_edits", [])
            if 0 <= index < len(edits):
                edits[index].setText(str(value)); pulse_widget(edits[index]); return True
        if field == "create_sample_forms":
            dp._build_sample_forms(); return True
        if field == "fill_down":
            dp._on_fill_down(); return True
        if field in ("sample_name", "exp_nr", "well", "raw_image_path",
                     "dimension_order", "pixel_distance_xy", "pixel_distance_z",
                     "time_interval", "time_unit"):
            forms = getattr(dp, "_sample_forms", [])
            if 0 <= index < len(forms):
                w = forms[index]["basic"].get(field)
                if w is not None and _push_to_widget(w, value):
                    pulse_widget(w); return True
        if field in ("dead_channel_number", "dead_mask_path"):
            forms = getattr(dp, "_sample_forms", [])
            if 0 <= index < len(forms):
                dead_key = "number" if field == "dead_channel_number" else "mask_path"
                w = forms[index].get("dead_channel", {}).get(dead_key)
                if w is not None and _push_to_widget(w, value):
                    pulse_widget(w); return True
        cell_field_map = {
            "cell_line": "line",
            "cell_condition": "condition",
            "cell_segments_image_path": "segments_image_path",
            "cell_tracks_image_path": "tracks_image_path",
            "cell_tracks_csv_path": "tracks_csv_path",
        }
        if field in cell_field_map:
            forms = getattr(dp, "_sample_forms", [])
            if 0 <= index < len(forms):
                w = _metadata_cell_widget(forms[index], cell_type, cell_field_map[field])
                if w is not None and _push_to_widget(w, value):
                    pulse_widget(w); return True
    except Exception:
        pass
    return False


def apply_bulk_fill_metadata(main_widget, action: ProposedAction) -> bool:
    """Fill the ENTIRE Metadata Builder in one guarded pass.

    Sets the population counts and cell-type names, then does exactly ONE
    configure pass and ONE sample-form rebuild, then writes every provided
    per-sample value. Pulses are suppressed so the batch doesn't schedule a
    storm of QTimers (the cause of the GUI freeze on "fill everything")."""
    dp = _dp(main_widget)
    if dp is None:
        return False
    d = action.data
    set_pulse_suppressed(True)
    try:
        if not getattr(dp.builder_grp, "isChecked", lambda: True)():
            dp.builder_grp.setChecked(True)
        # 1. population counts (only when provided)
        for field, attr in (("n_samples", "n_samples_spin"),
                            ("n_organoids", "n_organoid_spin"),
                            ("n_immune", "n_immune_spin"),
                            ("n_other", "n_other_spin")):
            if d.get(field) is not None:
                try:
                    getattr(dp, attr).setValue(int(d[field]))
                except Exception:
                    pass
        if d.get("include_dead") is not None:
            dp.include_dead_cb.setChecked(_coerce_bool(d["include_dead"]))
        # 2. ONE configure pass + cell-type names
        dp._on_configure_cell_types(force=True)
        for attr, key in (("_organoid_name_edits", "organoid_names"),
                          ("_immune_name_edits", "immune_names"),
                          ("_other_name_edits", "other_names")):
            edits = getattr(dp, attr, [])
            for i, name in enumerate(d.get(key, []) or []):
                if i < len(edits):
                    edits[i].setText(str(name))
        # 3. ONE sample-form rebuild (picks up the names just set)
        dp._build_sample_forms(force=True)
        # 4. all per-sample values in a single pass
        forms = getattr(dp, "_sample_forms", [])
        cell_field_map = {
            "line": "line", "cell_line": "line",
            "condition": "condition", "cell_condition": "condition",
            "segments_image_path": "segments_image_path",
            "cell_segments_image_path": "segments_image_path",
            "tracks_image_path": "tracks_image_path",
            "cell_tracks_image_path": "tracks_image_path",
            "tracks_csv_path": "tracks_csv_path",
            "cell_tracks_csv_path": "tracks_csv_path",
        }
        for idx, sample in enumerate(d.get("samples", []) or []):
            if idx >= len(forms) or not isinstance(sample, dict):
                break
            form = forms[idx]
            for k, v in sample.items():
                if k == "cell_types":
                    continue
                if k in ("dead_channel_number", "dead_mask_path"):
                    dk = "number" if k == "dead_channel_number" else "mask_path"
                    w = form.get("dead_channel", {}).get(dk)
                else:
                    w = form.get("basic", {}).get(k)
                if w is not None:
                    _push_to_widget(w, v)
            for ct_name, ct_vals in (sample.get("cell_types") or {}).items():
                ct_fields = (form.get("cell_types") or {}).get(str(ct_name), {})
                if not isinstance(ct_vals, dict):
                    continue
                for fk, fv in ct_vals.items():
                    target = ct_fields.get(cell_field_map.get(fk, fk))
                    if target is not None:
                        _push_to_widget(target, fv)
        return True
    except Exception:
        return False
    finally:
        set_pulse_suppressed(False)


def apply_select_segmentation_method(main_widget, value) -> bool:
    """Select the global Segmentation Method dropdown by its visible label.

    Matches exact → starts-with → contains so a leading label token ('APOC',
    'Pixel Classifier', 'Cellpose', …) resolves to the labelled combo item.
    Starts-with ordering keeps 'Pixel Classifier' from matching
    'ConvPaint (DL pixel classifier)'. Setting the index swaps the param page."""
    seg = getattr(main_widget, "segmentation_tab", None)
    combo = getattr(seg, "method_combo", None) if seg is not None else None
    if combo is None or value is None:
        return False
    try:
        from qtpy.QtCore import Qt
        s = str(value).strip()
        for flag in (Qt.MatchFixedString, Qt.MatchStartsWith, Qt.MatchContains):
            i = combo.findText(s, flag)
            if i >= 0:
                combo.setCurrentIndex(i)
                pulse_widget(combo)
                return True
    except Exception:
        pass
    return False


def apply_action(main_widget, action: ProposedAction) -> bool:
    if not action.ok:
        return False
    if action.kind == "set_parameter":
        stored, visible = apply_set_parameter(
            main_widget, action.data["key"], action.data["value"])
        action.data["widget_updated"] = visible   # dock uses this for an honest message
        return stored
    if action.kind == "set_ui_value":
        ok = apply_set_ui_value(main_widget, action.data["control_id"],
                                action.data["value"])
        action.data["widget_updated"] = ok
        return ok
    if action.kind == "navigate_to_step":
        return apply_navigate(main_widget, action.data["step"])
    if action.kind == "add_queue_step":
        return _apply_add_queue_step(main_widget, action)
    if action.kind == "fill_metadata_builder":
        ok = apply_fill_metadata_builder(main_widget, action)
        action.data["widget_updated"] = ok
        return ok
    if action.kind == "bulk_fill_metadata":
        ok = apply_bulk_fill_metadata(main_widget, action)
        action.data["widget_updated"] = ok
        return ok
    if action.kind == "select_segmentation_method":
        ok = apply_select_segmentation_method(main_widget, action.data.get("value"))
        action.data["widget_updated"] = ok
        return ok
    if action.kind == "show_track_length_distribution":
        return _apply_show_track_length_distribution(main_widget, action)
    if action.kind == "recommend_edt":
        return _apply_recommend_edt(main_widget, action)
    if action.kind == "summarize_track_counts":
        return _apply_summarize_track_counts(main_widget, action)
    if action.kind == "create_cell_type_group":
        return _apply_create_cell_type_group(main_widget, action)
    if action.kind == "create_btrack_config_copy":
        return _apply_create_btrack_config_copy(main_widget, action)
    if action.kind == "open_result":
        return _apply_open_result(main_widget, action)
    if action.kind == "save_metadata":
        dp = _dp(main_widget)
        callback = getattr(dp, "_on_save_metadata", None) if dp is not None else None
        if callback is None or not bool(callback()):
            return False
        action.data["result_markdown"] = (
            "Metadata was saved and activated for the other workflow tabs."
        )
        return True
    if action.kind == "load_metadata":
        dp = _dp(main_widget)
        callback = getattr(dp, "_on_load_metadata", None) if dp is not None else None
        if callback is None or not bool(callback()):
            return False
        action.data["result_markdown"] = (
            "Metadata loading has started. The Metadata Loader status will update "
            "when validation and image inspection finish."
        )
        return True
    if action.kind == "open_analysis_view":
        return _apply_open_analysis_view(main_widget, action.data.get("view"))
    return False


def _apply_open_analysis_view(main_widget, view: str) -> bool:
    if not apply_navigate(main_widget, "analysis"):
        return False
    analysis = getattr(main_widget, "analysis_tab", None)
    tabs = getattr(analysis, "inner_tabs", None) if analysis is not None else None
    if tabs is None:
        return False
    try:
        if view == "death_dynamics":
            tabs.setCurrentIndex(0)
            return True
        single = getattr(analysis, "single_cell_tab", None)
        start = getattr(single, "_on_guided_start", None) if single is not None else None
        if start is None:
            return False
        tabs.setCurrentIndex(1)
        start("state" if view == "behavioral_state" else "track")
        return True
    except Exception:
        return False


def _apply_show_track_length_distribution(main_widget, action: ProposedAction) -> bool:
    tab = getattr(main_widget, "filtering_tab", None)
    panel = (getattr(tab, "panels", {}) or {}).get(action.data.get("cell_type"))
    callback = getattr(panel, "_on_preview_lengths_clicked", None) if panel is not None else None
    if callback is None:
        return False
    callback()
    return True


def _apply_recommend_edt(main_widget, action: ProposedAction) -> bool:
    dp = _dp(main_widget)
    metadata = getattr(dp, "metadata", None) if dp is not None else None
    try:
        from behav3d.napari._assistant_context import _metadata_builder_state
        builder_state = _metadata_builder_state(dp)
    except Exception:
        builder_state = {}
    records = builder_state.get("draft_records") or []
    if not records and metadata is not None and not getattr(metadata, "empty", True):
        records = metadata.to_dict(orient="records")
    if not records:
        action.data["result_markdown"] = (
            "Load or save experiment metadata first so I can read the XY pixel size."
        )
        return True

    cell_type = action.data["cell_type"]
    organoid_names = builder_state.get("organoid_names") or []
    available_types = set(organoid_names)
    available_types.update(builder_state.get("immune_names") or [])
    available_types.update(builder_state.get("other_names") or [])
    is_organoid = cell_type in organoid_names
    try:
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata,
        )
        if metadata is not None and not getattr(metadata, "empty", True):
            metadata_organoids = detect_organoid_types_from_metadata(metadata)
            is_organoid = is_organoid or cell_type in metadata_organoids
            available_types.update(metadata_organoids)
            available_types.update(detect_immune_cell_types_from_metadata(metadata))
            available_types.update(detect_other_cell_types_from_metadata(metadata))
    except Exception:
        pass
    if available_types and cell_type not in available_types:
        action.data["result_markdown"] = (
            f"I could not find **{cell_type}** in the current metadata. Choose one of: "
            + ", ".join(sorted(available_types)) + "."
        )
        return True
    cells_across = action.data.get("organoid_cells_across")
    if is_organoid and cells_across is None:
        action.data["result_markdown"] = (
            f"For **{cell_type}**, approximately how many cell widths span the "
            "organoid diameter? I will use 10 µm per cell unless you provide a "
            "different typical cell diameter."
        )
        return True
    if not is_organoid:
        cells_across = None

    try:
        from behav3d.napari._assistant_recommendations import (
            calculate_edt_recommendations,
            format_edt_recommendations,
        )
        result = calculate_edt_recommendations(
            records,
            cell_diameter_um=action.data.get("cell_diameter_um", 10.0),
            organoid_cells_across=cells_across,
        )
        action.data["recommendation"] = result
        action.data["result_markdown"] = format_edt_recommendations(result, cell_type)
        return True
    except Exception as exc:
        action.data["result_markdown"] = f"Could not calculate EDT values: {exc}"
        return True


def _apply_summarize_track_counts(main_widget, action: ProposedAction) -> bool:
    dp = _dp(main_widget)
    output_dir = getattr(dp, "output_dir", "") if dp is not None else ""
    if not output_dir:
        action.data["result_markdown"] = (
            "Set the output directory first so I can locate the track feature CSV."
        )
        return True
    try:
        from behav3d.analysis.track_counts import (
            format_track_count_summary,
            generate_track_count_summary,
        )
        result = generate_track_count_summary(
            output_dir,
            action.data["cell_type"],
            minimum_lengths=action.data["minimum_lengths"],
            position_t=action.data["position_t"],
        )
        action.data["summary"] = result
        action.data["result_markdown"] = format_track_count_summary(
            result, action.data["cell_type"]
        )
        try:
            from behav3d.napari._results_panel import notify_results_changed
            notify_results_changed(main_widget)
        except Exception:
            pass
        return True
    except Exception as exc:
        action.data["result_markdown"] = f"Could not calculate track counts: {exc}"
        return True


def _apply_create_cell_type_group(main_widget, action: ProposedAction) -> bool:
    """Create the same metadata columns as the Feature Extraction grouping dialog."""
    dp = _dp(main_widget)
    md = getattr(dp, "metadata", None) if dp is not None else None
    if md is None or getattr(md, "empty", True):
        return False
    group_name = action.data["group_name"]
    members = action.data["members"]
    try:
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_merged_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata,
        )
        available = set(detect_organoid_types_from_metadata(md))
        available.update(detect_immune_cell_types_from_metadata(md))
        available.update(detect_other_cell_types_from_metadata(md))
        if not set(members).issubset(available):
            return False
        merged_name = f"{group_name}_merged"
        if group_name in available or merged_name in set(detect_merged_cell_types_from_metadata(md)):
            return False
        organoids = set(detect_organoid_types_from_metadata(md))
        immune = set(detect_immune_cell_types_from_metadata(md))
        prefix = "or_" if any(m in organoids for m in members) else (
            "im_" if any(m in immune for m in members) else "ot_")
        for suffix in ("line_condition", "segments_image_path",
                       "tracks_image_path", "tracks_csv_path"):
            md[f"{prefix}{merged_name}_{suffix}"] = None
        csv_path = (getattr(dp, "metadata_csv_path", None)
                    or getattr(dp, "metadata_csv", None)
                    or (getattr(dp, "behav3d_parameters", {}) or {})
                    .get("paths", {}).get("metadata_csv"))
        if csv_path:
            md.to_csv(csv_path, index=False)
        tab = getattr(main_widget, "feature_extraction_tab", None)
        if tab is not None and hasattr(tab, "_on_metadata_updated"):
            tab._on_metadata_updated()
        return True
    except Exception:
        return False


def _apply_create_btrack_config_copy(main_widget, action: ProposedAction) -> bool:
    from pathlib import Path
    import shutil

    destination = Path(action.data["destination"]).expanduser()
    if destination.exists():
        return False
    source = (Path(__file__).resolve().parents[1] / "preprocessing" / "tracking"
              / "models" / "cell_config.json")
    if not source.exists():
        return False
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    except OSError:
        return False
    control_id = f"tracking.{action.data['cell_type']}.btrack.config_preset"
    # Selecting Custom JSON exposes the path field; both writes are postcondition checked.
    preset_ok = apply_set_ui_value(main_widget, control_id, "Custom JSON")
    path_ok = apply_set_ui_value(
        main_widget,
        f"tracking.{action.data['cell_type']}.btrack.config_path",
        str(destination),
    )
    return preset_ok and path_ok


def _apply_open_result(main_widget, action: ProposedAction) -> bool:
    from pathlib import Path
    from behav3d.napari._results_catalog import scan_outputs

    dp = _dp(main_widget)
    output_dir = getattr(dp, "output_dir", "") if dp is not None else ""
    if not output_dir:
        return False
    results = scan_outputs(Path(output_dir))
    wanted = action.data["result_id"]
    match = next((item for item in results
                  if str(item.path.relative_to(Path(output_dir))) == wanted), None)
    if match is None:
        return False
    panel = getattr(main_widget, "results_panel", None)
    opener = getattr(panel, "_open_in_napari", None) if panel is not None else None
    if opener is None or not match.is_viewable:
        return False
    opener(match)
    return True


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
        "name": "set_ui_value",
        "description": (
            "Set one exact field that exists in ui_state.controls. Use its id exactly. "
            "Never invent an id and never use internal configuration keys. Metadata and "
            "Data Preparation edits always require confirmation, including blank fields. "
            "When the user explicitly asks to fill or change available values, emit the "
            "corresponding calls now instead of only describing the intended changes."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "control_id": {"type": "string"},
                "value": {},
            },
            "required": ["control_id", "value"],
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
                "step_type": {
                    "type": "string",
                    "enum": ["segment", "train", "track",
                             "feature_extract", "filter", "active_killing"],
                },
                "params": {"type": "object"},
            },
            "required": ["step_type"],
        },
    },
    {
        "name": "fill_metadata_builder",
        "description": (
            "Fill a field in the Metadata Builder on the Data Preparation tab. "
            "Use this to guide the user step-by-step through experiment setup. "
            "REQUIRED ORDER: open_builder → n_samples / n_organoids / n_immune / "
            "n_other / include_dead → configure_cell_types → organoid_name / "
            "immune_name / other_name (0-indexed; use immune_multicolor and "
            "immune_multicolor_channels when applicable) → create_sample_forms → "
            "per-sample fields (sample_name, exp_nr, well, raw_image_path, "
            "dimension_order, pixel_distance_xy, pixel_distance_z, time_interval, "
            "time_unit, dead_channel_number, dead_mask_path — sample index 0, 1, …). "
            "For per-sample cell-type rows, use cell_line, cell_condition, "
            "cell_segments_image_path, cell_tracks_image_path, or cell_tracks_csv_path "
            "with the exact visible cell_type name. "
            "After filling Sample 1, call fill_down to copy shared values to all samples."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "field": {
                    "type": "string",
                    "enum": [
                        "open_builder",
                        "n_samples", "n_organoids", "n_immune", "n_other", "include_dead",
                        "configure_cell_types",
                        "organoid_name", "immune_name", "other_name",
                        "immune_multicolor", "immune_multicolor_channels",
                        "create_sample_forms",
                        "sample_name", "exp_nr", "well", "raw_image_path",
                        "dimension_order", "pixel_distance_xy", "pixel_distance_z",
                        "time_interval", "time_unit",
                        "dead_channel_number", "dead_mask_path",
                        "cell_line", "cell_condition", "cell_segments_image_path",
                        "cell_tracks_image_path", "cell_tracks_csv_path",
                        "fill_down",
                    ],
                    "description": "Which field or action to perform.",
                },
                "value": {
                    "description": (
                        "The value to set. Omit for configure_cell_types, "
                        "create_sample_forms, fill_down, open_builder."
                    ),
                },
                "index": {
                    "type": "integer",
                    "description": (
                        "0-based index: cell type (for organoid_name/immune_name/other_name) "
                        "or sample (for per-sample fields)."
                    ),
                    "default": 0,
                },
                "cell_type": {
                    "type": "string",
                    "description": (
                        "Exact visible cell type name for per-sample cell-type fields "
                        "(cell_line, cell_condition, cell_segments_image_path, "
                        "cell_tracks_image_path, cell_tracks_csv_path)."
                    ),
                },
            },
            "required": ["field"],
        },
    },
    {
        "name": "bulk_fill_metadata",
        "description": (
            "Fill the ENTIRE Metadata Builder in ONE call. Use this (instead of "
            "many fill_metadata_builder calls) when the user gives lots of values "
            "at once, describes a new multi-sample experiment, provides a list of image "
            "files, or says 'fill everything'. Unknown fields may be omitted; include one "
            "possibly partial sample object per described movie so the forms are created. "
            "Never guess dimension_order, time_unit, image paths, sample names, or wells. "
            "This action opens the Metadata Builder itself; do not call open_builder first. "
            "The app rebuilds the forms once and applies all values in a single "
            "pass with no per-field confirmation. Provide the population counts, "
            "cell-type name lists, and one dict per sample."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "n_samples": {"type": "integer"},
                "n_organoids": {"type": "integer"},
                "n_immune": {"type": "integer"},
                "n_other": {"type": "integer"},
                "include_dead": {"type": "boolean"},
                "organoid_names": {"type": "array", "items": {"type": "string"}},
                "immune_names": {"type": "array", "items": {"type": "string"}},
                "other_names": {"type": "array", "items": {"type": "string"}},
                "samples": {
                    "type": "array",
                    "description": (
                        "One object per sample, in order. Recognised keys: "
                        "sample_name, exp_nr, well, raw_image_path, dimension_order, "
                        "pixel_distance_xy, pixel_distance_z, time_interval, time_unit, "
                        "dead_channel_number, dead_mask_path, and an optional "
                "'cell_types' object mapping each visible cell-type name to "
                "{line, condition, segments_image_path, tracks_image_path, "
                "tracks_csv_path}. For a configured population confirmed absent "
                "from a sample, use the literal line value 'not_added'."
                    ),
                    "items": {"type": "object"},
                },
            },
            "required": ["samples"],
        },
    },
    {
        "name": "recommend_edt",
        "description": (
            "Calculate EDT/watershed starting values from each sample's metadata XY "
            "resolution. For ordinary cells, use the default 10 um diameter unless "
            "the user provides another value. For organoids, first ask approximately "
            "how many cell widths span the organoid diameter, then pass that number as "
            "organoid_cells_across. This action reports recommendations but does not "
            "change the EDT field."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "cell_type": {"type": "string"},
                "cell_diameter_um": {
                    "type": "number", "default": 10,
                    "description": "Typical diameter of one cell in micrometers.",
                },
                "organoid_cells_across": {
                    "type": "number",
                    "description": (
                        "Approximate number of cell widths across the organoid diameter. "
                        "Required for organoid recommendations."
                    ),
                },
            },
            "required": ["cell_type"],
        },
    },
    {
        "name": "summarize_track_counts",
        "description": (
            "Read the unfiltered combined track-features CSV, count unique tracks "
            "present at one position_t, compare minimum track-length thresholds, and "
            "save the resulting sample-by-threshold table under quality_control. Use "
            "this for both the standard 20/50/100/200 summary and arbitrary interactive "
            "threshold or timepoint questions."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "cell_type": {"type": "string"},
                "position_t": {"type": "integer", "minimum": 0, "default": 0},
                "minimum_lengths": {
                    "type": "array",
                    "items": {"type": "integer", "minimum": 1},
                    "maxItems": 20,
                    "default": [20, 50, 100, 200],
                },
            },
            "required": ["cell_type"],
        },
    },
    {
        "name": "show_track_length_distribution",
        "description": "Open the existing track-length distribution preview for one cell type.",
        "parameters": {
            "type": "object",
            "properties": {"cell_type": {"type": "string"}},
            "required": ["cell_type"],
        },
    },
    {
        "name": "create_cell_type_group",
        "description": (
            "Create a merged cell-type group in metadata. This changes the metadata CSV "
            "and therefore always requires confirmation. It does not require re-tracking."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "group_name": {"type": "string"},
                "members": {"type": "array", "items": {"type": "string"}},
            },
            "required": ["group_name", "members"],
        },
    },
    {
        "name": "create_btrack_config_copy",
        "description": (
            "Copy the bundled cell btrack configuration to a new JSON file and select "
            "it for one cell type. File creation always requires confirmation."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "cell_type": {"type": "string"},
                "destination": {"type": "string"},
            },
            "required": ["cell_type", "destination"],
        },
    },
    {
        "name": "open_result",
        "description": "Open a viewable result listed in context.results using its exact id.",
        "parameters": {
            "type": "object",
            "properties": {"result_id": {"type": "string"}},
            "required": ["result_id"],
        },
    },
    {
        "name": "save_metadata",
        "description": (
            "Save the current valid Metadata Builder draft to metadata.csv and "
            "activate it for all workflow tabs. Use only when this action is "
            "available and the user explicitly asks to save. The client requires "
            "confirmation before writing."
        ),
        "parameters": {"type": "object", "properties": {}},
    },
    {
        "name": "load_metadata",
        "description": (
            "Start loading and validating the metadata CSV already selected in the "
            "Metadata Loader. Use only when this action is available and the user "
            "explicitly asks to load it. Loading completes asynchronously."
        ),
        "parameters": {"type": "object", "properties": {}},
    },
    {
        "name": "open_analysis_view",
        "description": (
            "Open one specific Analysis view. Use this instead of navigating to the "
            "generic Analysis tab when the user names Death Dynamics, Behavioral "
            "State, or State Trajectory."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "view": {
                    "type": "string",
                    "enum": [
                        "death_dynamics",
                        "behavioral_state",
                        "state_trajectory",
                    ],
                },
            },
            "required": ["view"],
        },
    },
]
