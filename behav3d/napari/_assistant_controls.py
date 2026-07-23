"""Live UI control contract used by the BEHAV3D assistant.

The assistant may only edit controls returned by :func:`control_bindings`.
Keeping context serialization and action application on this same registry
prevents stale config keys from being mistaken for fields on the screen.
"""
from __future__ import annotations

from typing import Any, Callable


CONTROL_CONTRACT_VERSION = "2.5"


def _safe(fn: Callable, default=None):
    try:
        return fn()
    except Exception:
        return default


def _widget_value(widget):
    if widget is None:
        return None
    if hasattr(widget, "isChecked"):
        return _safe(widget.isChecked)
    if hasattr(widget, "currentText"):
        return _safe(widget.currentText)
    if hasattr(widget, "value"):
        return _safe(widget.value)
    if hasattr(widget, "text"):
        return _safe(widget.text)
    return None


def _widget_choices(widget) -> list[str] | None:
    if widget is None or not hasattr(widget, "count") or not hasattr(widget, "itemText"):
        return None
    return [_safe(lambda i=i: str(widget.itemText(i)), "") for i in range(widget.count())]


def _coerce_bool(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def set_widget_value(widget, value: Any) -> bool:
    """Set a common Qt input using duck typing so this module stays testable."""
    if widget is None:
        return False
    try:
        if hasattr(widget, "setChecked") and hasattr(widget, "isChecked"):
            widget.setChecked(_coerce_bool(value))
            return True
        if hasattr(widget, "setCurrentIndex") and hasattr(widget, "itemText"):
            wanted = str(value).strip().lower()
            choices = _widget_choices(widget) or []
            matches = [
                i for i, label in enumerate(choices)
                if label.lower() == wanted or label.lower().startswith(wanted)
            ]
            if not matches:
                matches = [i for i, label in enumerate(choices) if wanted in label.lower()]
            if not matches:
                return False
            widget.setCurrentIndex(matches[0])
            return True
        if hasattr(widget, "setValue") and hasattr(widget, "value"):
            current = widget.value()
            widget.setValue(int(value) if isinstance(current, int) and not isinstance(current, bool)
                            else float(value))
            return True
        if hasattr(widget, "setText"):
            widget.setText(str(value))
            return True
    except (TypeError, ValueError):
        return False
    return False


def _is_visible(widget) -> bool:
    if widget is None:
        return False
    return not bool(_safe(widget.isHidden, False)) if hasattr(widget, "isHidden") else True


def _is_enabled(widget) -> bool:
    return bool(_safe(widget.isEnabled, True)) if widget is not None else False


def _binding(
    control_id: str,
    label: str,
    widget,
    *,
    step: str,
    default: Any = None,
    unit: str | None = None,
    method: str | None = None,
    strategy: str | None = None,
    cell_type: str | None = None,
    visible: bool | None = None,
    getter: Callable | None = None,
    setter: Callable | None = None,
    persist: Callable | None = None,
) -> dict:
    get = getter or (lambda w=widget: _widget_value(w))
    set_ = setter or (lambda value, w=widget: set_widget_value(w, value))
    choices = _widget_choices(widget)
    record = {
        "id": control_id,
        "label": label,
        "value": _safe(get),
        "default": default,
        "unit": unit,
        "choices": choices,
        "visible": _is_visible(widget) if visible is None else bool(visible),
        "enabled": _is_enabled(widget),
        "step": step,
        "method": method,
        "strategy": strategy,
        "cell_type": cell_type,
        "_widget": widget,
        "_getter": get,
        "_setter": set_,
        "_persist": persist,
    }
    for name, attr in (("minimum", "minimum"), ("maximum", "maximum")):
        if widget is not None and hasattr(widget, attr):
            record[name] = _safe(getattr(widget, attr))
    return record


def _checkset_binding(control_id, label, checks, **kwargs):
    checks = checks or {}

    def get():
        return [name for name, check in checks.items() if _safe(check.isChecked, False)]

    def set_(value):
        if not isinstance(value, (list, tuple, set)):
            return False
        selected = {str(item) for item in value}
        if not selected.issubset(set(checks)):
            return False
        for name, check in checks.items():
            if _is_enabled(check):
                check.setChecked(name in selected)
        return True

    first = next(iter(checks.values()), None)
    item = _binding(control_id, label, first, getter=get, setter=set_, **kwargs)
    item["choices"] = list(checks)
    item["enabled"] = any(_is_enabled(check) for check in checks.values())
    return item


def _metadata_bindings(main_widget) -> list[dict]:
    dp = getattr(main_widget, "data_prep_tab", None)
    if dp is None:
        return []
    out = []
    def set_output_dir(value):
        widget = getattr(dp, "output_dir_edit", None)
        if not set_widget_value(widget, value):
            return False
        dp.output_dir = str(value)
        params = getattr(dp, "behav3d_parameters", None)
        if isinstance(params, dict):
            params.setdefault("paths", {})["output_dir"] = str(value)
        return True

    globals_ = [
        ("data.output_directory", "Output directory", getattr(dp, "output_dir_edit", None),
         set_output_dir),
        ("data.metadata_csv", "Metadata CSV", getattr(dp, "csv_path_edit", None), None),
        ("data.dimension_order_all", "Dimension order for all samples",
         getattr(dp, "dim_apply_all_combo", None), None),
    ]
    for cid, label, widget, setter in globals_:
        if widget is not None:
            out.append(_binding(cid, label, widget, step="data_preparation", setter=setter))

    builder = [
        ("metadata.number_of_samples", "Number of samples", "n_samples_spin", 1),
        ("metadata.number_of_organoid_types", "Number of organoid types", "n_organoid_spin", 0),
        ("metadata.number_of_immune_types", "Number of immune cell types", "n_immune_spin", 0),
        ("metadata.number_of_other_types", "Number of other cell types", "n_other_spin", 0),
        ("metadata.include_dead_channel", "Include a dead-cell channel", "include_dead_cb", False),
    ]
    for cid, label, attr, default in builder:
        widget = getattr(dp, attr, None)
        if widget is not None:
            out.append(_binding(cid, label, widget, step="data_preparation", default=default))

    basic_labels = {
        "sample_name": "Sample name", "exp_nr": "Experiment number", "well": "Well",
        "raw_image_path": "Raw image path", "dimension_order": "Dimension order",
        "pixel_distance_xy": "Pixel size XY", "pixel_distance_z": "Pixel size Z",
        "time_interval": "Time interval", "time_unit": "Time unit",
    }
    units = {"pixel_distance_xy": "um", "pixel_distance_z": "um"}
    for index, form in enumerate(getattr(dp, "_sample_forms", []) or []):
        for field, widget in (form.get("basic") or {}).items():
            out.append(_binding(
                f"metadata.samples.{index}.{field}",
                f"Sample {index + 1}: {basic_labels.get(field, field.replace('_', ' '))}",
                widget, step="data_preparation", unit=units.get(field),
            ))
        for field, widget in (form.get("dead_channel") or {}).items():
            suffix = "dead_channel" if field == "number" else "dead_mask_path"
            label = "Dead channel number" if field == "number" else "Dead mask path"
            out.append(_binding(f"metadata.samples.{index}.{suffix}",
                                f"Sample {index + 1}: {label}", widget,
                                step="data_preparation"))
        for cell_type, fields in (form.get("cell_types") or {}).items():
            for field, widget in fields.items():
                out.append(_binding(
                    f"metadata.samples.{index}.cell_types.{cell_type}.{field}",
                    f"Sample {index + 1}, {cell_type}: {field.replace('_', ' ')}",
                    widget, step="data_preparation", cell_type=str(cell_type),
                ))
    return out


def _segmentation_bindings(main_widget) -> list[dict]:
    seg = getattr(main_widget, "segmentation_tab", None)
    combo = getattr(seg, "method_combo", None) if seg is not None else None
    if combo is None:
        return []
    current_index = _safe(combo.currentIndex, 0)
    out = [_binding("segmentation.method", "Segmentation method", combo,
                    step="segmentation")]
    pages = [
        (0, "apoc", "APOC", "apoc_page"),
        (1, "convpaint", "ConvPaint", "convpaint_page"),
        (2, "random_forest", "Pixel Classifier (Random Forest)", "pixel_classifier_page"),
    ]
    for index, token, method, attr in pages:
        page = getattr(seg, attr, None)
        widget = getattr(page, "spin_examples", None) if page is not None else None
        if widget is not None:
            persist_name = {
                "apoc": "_save_apoc_params_to_yaml",
                "convpaint": "_save_params_to_yaml",
                "random_forest": "_persist_params",
            }[token]
            out.append(_binding(
                f"segmentation.{token}.examples_per_sample",
                f"{method}: examples per sample", widget,
                step="segmentation", default=3, method=method,
                visible=current_index == index,
                persist=getattr(page, persist_name, None),
            ))

    def instance_strategy(training, tab) -> str | None:
        per_tab = getattr(tab, "_per_tab_strategy_combo", None)
        global_combo = getattr(training, "strategy_combo", None)
        combo = per_tab if per_tab is not None else global_combo
        return _safe(combo.currentText, None) if combo is not None else None

    # Instance-segmentation fields are created dynamically for the selected
    # strategy. Expose only controls that actually exist so advice and edits are
    # grounded in the active method, cell type, and watershed strategy.
    apoc = getattr(seg, "apoc_page", None)
    apoc_training = getattr(apoc, "_training_widget", None) if apoc is not None else None
    for cell_type, tab in (getattr(apoc_training, "tabs", {}) or {}).items():
        def persist_apoc(page=apoc):
            page._save_apoc_params_to_yaml(
                updated_apoc_params=page._collect_apoc_tab_config()
            )

        strategy = instance_strategy(apoc_training, tab)
        strategy_combo = getattr(tab, "_per_tab_strategy_combo", None)
        if strategy_combo is not None:
            out.append(_binding(
                f"segmentation.apoc.{cell_type}.strategy",
                f"{cell_type}: APOC instance segmentation strategy",
                strategy_combo, step="segmentation", method="APOC",
                strategy=strategy, cell_type=str(cell_type),
                visible=current_index == 0, persist=persist_apoc,
            ))
        for suffix, label, attr, unit in (
            ("mask_threshold", "mask threshold", "prob_mask_threshold_spin", None),
            ("seed_threshold", "seed threshold", "prob_seed_threshold_spin", None),
            ("edt_threshold", "EDT threshold", "edt_threshold_spin", "pixels"),
        ):
            widget = getattr(tab, attr, None)
            if widget is None:
                continue
            out.append(_binding(
                f"segmentation.apoc.{cell_type}.{suffix}",
                f"{cell_type}: APOC {label}", widget,
                step="segmentation", unit=unit, method="APOC",
                strategy=strategy, cell_type=str(cell_type),
                visible=current_index == 0, persist=persist_apoc,
            ))

    convpaint = getattr(seg, "convpaint_page", None)
    conv_training = (
        getattr(convpaint, "_training_widget", None) if convpaint is not None else None
    )
    for cell_type, tab in (getattr(conv_training, "tabs", {}) or {}).items():
        strategy = instance_strategy(conv_training, tab)

        def persist_convpaint(page=convpaint, training_tab=tab):
            page._save_params_to_yaml(training_tab.collect_params())

        strategy_combo = getattr(tab, "_per_tab_strategy_combo", None)
        if strategy_combo is not None:
            out.append(_binding(
                f"segmentation.convpaint.{cell_type}.strategy",
                f"{cell_type}: ConvPaint instance segmentation strategy",
                strategy_combo, step="segmentation", method="ConvPaint",
                strategy=strategy, cell_type=str(cell_type),
                visible=current_index == 1, persist=persist_convpaint,
            ))
        for suffix, label, attr, unit in (
            ("mask_threshold", "mask threshold", "prob_mask_threshold_spin", None),
            ("seed_threshold", "seed threshold", "prob_seed_threshold_spin", None),
            ("edt_threshold", "EDT threshold", "edt_threshold_spin", "pixels"),
        ):
            widget = getattr(tab, attr, None)
            if widget is None:
                continue
            out.append(_binding(
                f"segmentation.convpaint.{cell_type}.{suffix}",
                f"{cell_type}: ConvPaint {label}", widget,
                step="segmentation", unit=unit, method="ConvPaint",
                strategy=strategy, cell_type=str(cell_type),
                visible=current_index == 1, persist=persist_convpaint,
            ))

    random_forest = getattr(seg, "pixel_classifier_page", None)
    for cell_type, widgets in (getattr(random_forest, "param_widgets", {}) or {}).items():
        widget = (widgets or {}).get("edt")
        if widget is None:
            continue
        out.append(_binding(
            f"segmentation.random_forest.{cell_type}.edt_threshold",
            f"{cell_type}: Pixel Classifier EDT threshold", widget,
            step="segmentation", unit="pixels",
            method="Pixel Classifier (Random Forest)", cell_type=str(cell_type),
            visible=current_index == 2,
            persist=getattr(random_forest, "_persist_params", None),
        ))
    cellpose = getattr(seg, "cellpose_page", None)
    channels = getattr(cellpose, "spin_manual_channels", None) if cellpose is not None else None
    if channels is not None:
        out.append(_binding("segmentation.cellpose.number_of_channels",
                            "Cellpose: number of channels", channels,
                            step="segmentation", method="Cellpose",
                            visible=current_index == 3,
                            persist=getattr(cellpose, "_persist_channel_config", None)))
    return out


def _tracking_bindings(main_widget) -> list[dict]:
    tab = getattr(main_widget, "tracking_tab", None)
    panels = getattr(tab, "panels", {}) or {} if tab is not None else {}
    out = []
    specs = [
        ("method", "Tracking method", "combo_method", None, None),
        ("lap.frame_link_distance", "LAP frame-to-frame linking distance", "lap_track_cost", "px", "LAP"),
        ("lap.gap_closing_distance", "LAP gap-closing distance", "lap_gap_cost", "px", "LAP"),
        ("lap.maximum_gap", "LAP maximum gap", "lap_gap_frames", "frames", "LAP"),
        ("lap.merging_distance", "LAP merging distance", "lap_merge_cost", "px", "LAP"),
        ("lap.splitting_distance", "LAP splitting distance", "lap_split_cost", "px", "LAP"),
        ("trackpy.search_range", "TrackPy search range", "tp_search_range", "px", "TrackPy"),
        ("trackpy.memory", "TrackPy memory", "tp_memory", "frames", "TrackPy"),
        ("trackpy.adaptive_stop", "TrackPy adaptive stop", "tp_adaptive_stop", "px", "TrackPy"),
        ("trackpy.adaptive_step", "TrackPy adaptive step", "tp_adaptive_step", None, "TrackPy"),
        ("btrack.config_preset", "btrack configuration preset", "bt_config_preset", None, "btrack"),
        ("btrack.config_path", "btrack custom configuration", "bt_config_path", None, "btrack"),
        ("btrack.maximum_search_radius", "btrack maximum search radius", "bt_max_search_radius", "um", "btrack"),
        ("btrack.update_method", "btrack update method", "bt_update_method", None, "btrack"),
        ("btrack.step_size", "btrack step size", "bt_step_size", "frames", "btrack"),
        ("btrack.use_visual_features", "Use visual features", "bt_use_visual_features", None, "btrack"),
        ("btrack.workers", "btrack workers", "bt_n_workers", None, "btrack"),
        ("btrack.use_global_optimization", "Enable global track optimization", "bt_use_optimize", None, "btrack"),
        ("btrack.distance_threshold", "btrack optimizer distance threshold", "bt_dist_thresh", "um", "btrack"),
        ("btrack.time_threshold", "btrack optimizer time threshold", "bt_time_thresh", "frames", "btrack"),
    ]
    for cell_type, panel in panels.items():
        method_index = _safe(panel.combo_method.currentIndex, 0)
        optimizer_enabled = bool(_safe(
            getattr(panel, "bt_use_optimize", None).isChecked, False
        )) if getattr(panel, "bt_use_optimize", None) is not None else False
        for suffix, label, attr, unit, method in specs:
            if suffix in {
                "btrack.distance_threshold", "btrack.time_threshold",
            } and not optimizer_enabled:
                continue
            widget = getattr(panel, attr, None)
            if widget is None:
                continue
            visible = True
            if method == "LAP":
                visible = method_index == 0
            elif method == "TrackPy":
                visible = method_index == 1
            elif method == "btrack":
                visible = method_index == 3
            out.append(_binding(f"tracking.{cell_type}.{suffix}",
                                f"{cell_type}: {label}", widget, step="tracking",
                                unit=unit, method=method, cell_type=str(cell_type),
                                visible=visible, persist=getattr(panel, "_persist", None)))
        checks = getattr(panel, "bt_hyp_checks", {}) or {}
        if checks and optimizer_enabled:
            out.append(_checkset_binding(
                f"tracking.{cell_type}.btrack.hypotheses",
                f"{cell_type}: btrack optimization hypotheses", checks,
                step="tracking", method="btrack", cell_type=str(cell_type),
                visible=method_index == 3,
                persist=getattr(panel, "_persist", None),
            ))
    return out


def _feature_bindings(main_widget) -> list[dict]:
    tab = getattr(main_widget, "feature_extraction_tab", None)
    panels = getattr(tab, "panels", {}) or {} if tab is not None else {}
    out = []
    for cell_type, panel in panels.items():
        checks = getattr(panel, "feature_checks", {}) or {}
        if checks:
            out.append(_checkset_binding(f"features.{cell_type}.feature_groups",
                                         f"{cell_type}: feature groups", checks,
                                         step="feature_extraction", cell_type=str(cell_type),
                                         persist=getattr(panel, "_persist", None)))
        for suffix, label, attr, unit in [
            ("contact_distance", "contact distance", "contact_threshold", "um"),
            ("dead_pixel_threshold", "dead-pixel threshold", "spin_dead_threshold", "%"),
            ("workers", "workers", "spin_workers", None),
        ]:
            widget = getattr(panel, attr, None)
            if widget is not None:
                out.append(_binding(f"features.{cell_type}.{suffix}",
                                    f"{cell_type}: {label}", widget,
                                    step="feature_extraction", unit=unit,
                                    cell_type=str(cell_type), persist=getattr(panel, "_persist", None)))
    return out


def _filter_bindings(main_widget) -> list[dict]:
    tab = getattr(main_widget, "filtering_tab", None)
    panels = getattr(tab, "panels", {}) or {} if tab is not None else {}
    out = []
    specs = [
        ("trim_to_maximum.enabled", "Trim full time series", "en_exp_duration", None),
        ("trim_to_maximum.timepoints", "Maximum timepoints", "spin_exp_duration", "timepoints"),
        ("minimum_length.enabled", "Filter short tracks", "en_min_length", None),
        ("minimum_length.timepoints", "Minimum track length", "spin_min_length", "timepoints"),
        ("maximum_length.enabled", "Trim retained tracks to a common length", "en_max_length", None),
        ("maximum_length.timepoints", "Common output track length", "spin_max_length", "timepoints"),
        ("minimum_initial_size.enabled", "Filter by initial size", "check_filter_min_size", None),
        ("minimum_initial_size.pixels", "Minimum initial size", "spin_min_size_t1", "pixels"),
        ("dead_at_first_timepoint.enabled", "Filter dead cells at first timepoint", "check_filter_dead_t0", None),
        ("time_unit", "Filtering time unit", "combo_time_type", None),
    ]
    for cell_type, panel in panels.items():
        for suffix, label, attr, unit in specs:
            widget = getattr(panel, attr, None)
            if widget is not None:
                out.append(_binding(f"filtering.{cell_type}.{suffix}",
                                    f"{cell_type}: {label}", widget, step="filtering",
                                    unit=unit, cell_type=str(cell_type),
                                    persist=getattr(panel, "_persist", None)))
    return out


def _analysis_bindings(main_widget) -> list[dict]:
    analysis = getattr(main_widget, "analysis_tab", None)
    single = getattr(analysis, "single_cell_tab", None) if analysis is not None else None
    state = getattr(single, "state_tab", None) if single is not None else None
    if state is None:
        return []
    cell_type = _safe(single.cell_type_combo.currentText, "") or None
    out = []
    persist = ((lambda ct=cell_type: state._persist_state_cfg(ct))
               if cell_type else None)
    mode_combo = getattr(state, "combo_hmm_n_states_mode", None)
    sticky_check = getattr(state, "chk_hmm_sticky", None)
    mode = _safe(mode_combo.currentText, "fixed") if mode_combo is not None else "fixed"
    sticky = bool(_safe(sticky_check.isChecked, False)) if sticky_check is not None else False
    specs = [
        ("window_size", "Window size", "spin_window_size", "timepoints", None),
        ("smooth_window", "Feature smoothing window", "spin_hmm_feature_smoothing_window", "timepoints", None),
        ("lower_quantile_cap", "Low percentile cap", "spin_quant_lo", None, None),
        ("upper_quantile_cap", "High percentile cap", "spin_quant_hi", None, None),
        ("state_selection_mode", "State selection mode", "combo_hmm_n_states_mode", None, None),
        ("number_of_states", "Number of states", "spin_hmm_n_states", None, mode == "fixed"),
        ("minimum_states", "Minimum states", "spin_hmm_k_min", None, mode == "auto"),
        ("maximum_states", "Maximum states", "spin_hmm_k_max", None, mode == "auto"),
        ("start_offset", "Start offset", "spin_hmm_start_offset", "timepoints", None),
        ("skipped_frames", "Skipped frames", "combo_hmm_start_offset_fill_mode", None, None),
        ("covariance_type", "Covariance type", "combo_hmm_covariance_type", None, None),
        ("maximum_iterations", "Maximum iterations", "spin_hmm_n_iter", None, None),
        ("convergence_tolerance", "Convergence tolerance", "spin_hmm_tol", None, None),
        ("minimum_covariance", "Minimum covariance", "spin_hmm_min_covar", None, None),
        ("sticky.enabled", "Use sticky HMM", "chk_hmm_sticky", None, None),
        ("sticky.kappa", "Stickiness kappa", "spin_hmm_stickiness_kappa", None, sticky),
        ("sticky.transition_alpha", "Transition matrix alpha", "spin_hmm_transmat_alpha", None, sticky),
        ("random_seed", "Random seed", "spin_seed", None, None),
    ]
    for suffix, label, attr, unit, visible in specs:
        widget = getattr(state, attr, None)
        if widget is not None:
            out.append(_binding(f"analysis.state_classification.{cell_type or 'selected'}.{suffix}",
                                f"{cell_type or 'Selected cell type'}: {label}", widget,
                                step="analysis", method="HMM", cell_type=cell_type,
                                unit=unit, visible=visible, persist=persist))
    window_checks = {
        "net_displacement": getattr(state, "chk_net_disp", None),
        "straightness": getattr(state, "chk_straight", None),
        "mean_square_displacement": getattr(state, "chk_msd", None),
    }
    window_checks = {key: widget for key, widget in window_checks.items()
                     if widget is not None}
    for suffix, label, checks in [
        ("timepoint_features", "Timepoint features", getattr(state, "_timepoint_checkboxes", {})),
        ("log_scaled_features", "Log-scaled features", getattr(state, "_logscale_checkboxes", {})),
        ("binary_feature_groups", "Binary feature groups", getattr(state, "_bingrp_checkboxes", {})),
        ("window_features", "Additional window features", window_checks),
    ]:
        if checks:
            out.append(_checkset_binding(
                f"analysis.state_classification.{cell_type or 'selected'}.{suffix}",
                f"{cell_type or 'Selected cell type'}: {label}", checks,
                step="analysis", method="HMM", cell_type=cell_type,
                persist=persist,
            ))
    return out


def control_bindings(main_widget) -> list[dict]:
    """Return every assistant-editable control currently instantiated in the UI."""
    out = []
    for builder in (_metadata_bindings, _segmentation_bindings, _tracking_bindings,
                    _feature_bindings, _filter_bindings, _analysis_bindings):
        out.extend(_safe(lambda b=builder: b(main_widget), []) or [])
    return out


def control_registry(main_widget) -> list[dict]:
    """Return JSON-safe public records for the model context."""
    return [{k: v for k, v in binding.items() if not k.startswith("_")}
            for binding in control_bindings(main_widget)]


def find_control(main_widget, control_id: str) -> dict | None:
    return next((item for item in control_bindings(main_widget)
                 if item["id"] == control_id), None)


def active_cell_type(main_widget, step: str) -> str | None:
    tab_attr = {
        "tracking": "tracking_tab",
        "feature_extraction": "feature_extraction_tab",
        "filtering": "filtering_tab",
    }.get(step)
    if tab_attr:
        tab = getattr(main_widget, tab_attr, None)
        tabs = getattr(tab, "cell_tabs", None) if tab is not None else None
        panel = _safe(tabs.currentWidget) if tabs is not None else None
        return getattr(panel, "cell_type", None)
    if step == "segmentation":
        segmentation = getattr(main_widget, "segmentation_tab", None)
        method_combo = getattr(segmentation, "method_combo", None)
        method_index = _safe(method_combo.currentIndex, -1) if method_combo is not None else -1
        page_attr = {0: "apoc_page", 1: "convpaint_page"}.get(method_index)
        page = getattr(segmentation, page_attr, None) if page_attr else None
        training = getattr(page, "_training_widget", None) if page is not None else None
        tabs = getattr(training, "tab_widget", None) if training is not None else None
        index = _safe(tabs.currentIndex, -1) if tabs is not None else -1
        cell_types = getattr(training, "_tab_cell_types", []) or []
        if 0 <= index < len(cell_types):
            return str(cell_types[index])
        # The random-forest page has one vertical panel per type rather than an
        # active sub-tab; return the only type when the choice is unambiguous.
        if method_index == 2:
            random_forest = getattr(segmentation, "pixel_classifier_page", None)
            available = list((getattr(random_forest, "param_widgets", {}) or {}).keys())
            if len(available) == 1:
                return str(available[0])
        return None
    if step == "analysis":
        analysis = getattr(main_widget, "analysis_tab", None)
        single = getattr(analysis, "single_cell_tab", None) if analysis is not None else None
        return _safe(single.cell_type_combo.currentText) if single is not None else None
    return None
