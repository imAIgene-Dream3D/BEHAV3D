"""Replay PI feedback scenarios against a deployed BEHAV3D Assistant API."""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path
from typing import Any

import requests

from behav3d.napari._assistant_actions import TOOL_SCHEMA


DATASET_DESCRIPTION = (
    "I have an intravital imaging dataset of a popliteal lymph node and want to set "
    "up a BEHAV3D analysis. I have three movies, 2 immune types and 1 other type, "
    "collagen. CH1 is T cells labelled with CMTMR, CH2 is HIV-GFP-infected T cells, "
    "and CH3 is collagen. Z spacing is 4 um, XY pixel size is 1.15 um, and time-lapses "
    "were acquired every 15, but I do not know whether that was seconds or minutes. "
    "Help me set up the analysis. Where should I start?"
)


def _control(
    control_id: str,
    label: str,
    value: Any,
    *,
    choices: list[str] | None = None,
    method: str | None = None,
    strategy: str | None = None,
    cell_type: str | None = None,
    unit: str | None = None,
    visible: bool = True,
    enabled: bool = True,
    required_choices: list[str] | None = None,
    active: bool | None = None,
) -> dict:
    control = {
        "id": control_id,
        "label": label,
        "value": value,
        "choices": choices,
        "visible": visible,
        "enabled": enabled,
        "method": method,
        "strategy": strategy,
        "cell_type": cell_type,
        "unit": unit,
    }
    if required_choices is not None:
        control["required_choices"] = required_choices
    if active is not None:
        control["active"] = active
    return control


def _metadata_records() -> list[dict]:
    return [
        {
            "sample_name": f"Movie{i}",
            "pixel_distance_xy": 1.15,
            "pixel_distance_z": 4.0,
            "time_interval": 15,
            "time_unit": "s",
        }
        for i in range(1, 4)
    ]


def _context(step: str, controls: list[dict], **extra) -> dict:
    context = {
        "current_step": step,
        "current_tab_label": step.replace("_", " ").title(),
        "assistant_session": {"intent": "free_form"},
        "ui_state": {"controls": controls},
        "metadata": {"loaded": True, "records": _metadata_records(), "validation": []},
    }
    context.update(extra)
    return context


def _metadata_setup_case() -> dict:
    return {
        "name": "metadata_setup",
        "messages": [{"role": "user", "content": DATASET_DESCRIPTION}],
        "context": _context(
            "data_preparation", [],
            metadata={"loaded": False, "records": [], "validation": []},
            metadata_builder={"open": False, "n_samples": 1},
        ),
        "check": _check_metadata_setup,
    }


def _pixel_fill_case() -> dict:
    controls = []
    records = []
    for index in range(3):
        records.append({
            "sample_name": f"Movie{index + 1}",
            "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0,
            "time_interval": 15,
            "time_unit": "s",
        })
        controls.extend([
            _control(
                f"metadata.samples.{index}.pixel_distance_xy",
                f"Sample {index + 1}: XY pixel size", 0.5, unit="um",
            ),
            _control(
                f"metadata.samples.{index}.pixel_distance_z",
                f"Sample {index + 1}: Z pixel size", 2.0, unit="um",
            ),
        ])
    return {
        "name": "fill_pixel_sizes",
        "messages": [
            {"role": "user", "content": DATASET_DESCRIPTION},
            {"role": "assistant", "content": "I can help prepare the metadata."},
            {"role": "user", "content": "Can you fill in the correct pixel size?"},
        ],
        "context": _context(
            "data_preparation", controls,
            metadata={
                "loaded": True,
                "record_source": "metadata_builder_draft",
                "records": records,
                "validation": [],
                "save_required": True,
            },
            metadata_builder={"open": True, "n_samples": 3},
        ),
        "check": _check_pixel_fill,
    }


def _segmentation_case() -> dict:
    strategy = "APOC Probability Map + Watershed"
    controls = [
        _control(
            "segmentation.apoc.Tcell_cmtrm.mask_threshold",
            "Tcell cmtrm: APOC mask threshold", 0.4,
            method="APOC", strategy=strategy, cell_type="Tcell_cmtrm",
        ),
        _control(
            "segmentation.apoc.Tcell_cmtrm.seed_threshold",
            "Tcell cmtrm: APOC seed threshold", 0.7,
            method="APOC", strategy=strategy, cell_type="Tcell_cmtrm",
        ),
    ]
    return {
        "name": "merged_cell_thresholds",
        "messages": [{
            "role": "user",
            "content": "I see a lot of cells that are not split in the T cell CMTMR channel.",
        }],
        "context": _context(
            "segmentation", controls,
            active_cell_type="Tcell_cmtrm",
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {"cell_type_strategies": {"Tcell_cmtrm": strategy}},
            },
        ),
        "check": _check_segmentation,
    }


def _method_requires_signal_context_case() -> dict:
    metadata = {
        "loaded": True,
        "records": [{
            "sample_name": "Movie1",
            "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0,
            "time_interval": 2,
            "time_unit": "min",
        }],
        "image_dimensions": [{
            "sample_name": "Movie1",
            "shape": "(50, 4, 20, 512, 512)",
            "dimension_order": "TCZYX",
            "channel_count": 4,
            "timepoint_count": 50,
        }],
        "validation": [],
    }
    return {
        "name": "segmentation_method_asks_about_signal",
        "messages": [{
            "role": "user",
            "content": "Choose the best segmentation method for this experiment.",
        }],
        "context": _context(
            "segmentation", [], metadata=metadata,
            segmentation={"method": "APOC (GPU)", "method_index": 0},
        ),
        "check": _check_method_requires_signal_context,
    }


def _cellpose_requires_bleed_confirmation_case() -> dict:
    return {
        "name": "cellpose_sam_requires_bleed_confirmation",
        "messages": [{
            "role": "user",
            "content": "Would Cellpose-SAM work for my data?",
        }],
        "context": _context(
            "segmentation", [],
            segmentation={"method": "APOC (GPU)", "method_index": 0},
        ),
        "check": _check_method_requires_signal_context,
    }


def _confirmed_bleed_through_case() -> dict:
    return {
        "name": "confirmed_bleed_through_prefers_apoc",
        "messages": [
            {
                "role": "user",
                "content": "Would Cellpose-SAM work for my data?",
            },
            {
                "role": "assistant",
                "content": (
                    "Is each target isolated in a clean channel, or is signal from "
                    "other cell types visible in the same channel?"
                ),
            },
            {
                "role": "user",
                "content": (
                    "There is signal bleed-through: some channels contain signal "
                    "from more than one cell type."
                ),
            },
        ],
        "context": _context(
            "segmentation", [],
            segmentation={"method": "APOC (GPU)", "method_index": 0},
        ),
        "check": _check_confirmed_bleed_through,
    }


def _apoc_channel_map_case() -> dict:
    choices = ["Channel 0", "Channel 1", "Channel 2", "Channel 3"]
    controls = [
        _control(
            f"segmentation.apoc.{cell_type}.input_channels",
            f"{cell_type}: APOC image channel inputs",
            list(choices),
            choices=choices,
            method="APOC",
            cell_type=cell_type,
        )
        for cell_type in ("27t", "mdo", "tcell", "dead")
    ]
    apoc_cells = {
        cell_type: {
            "available_input_channels": choices,
            "selected_input_channels": choices,
            "channel_controls_ready": True,
            "trained_classifier_found": False,
        }
        for cell_type in ("27t", "mdo", "tcell", "dead")
    }
    return {
        "name": "apoc_applies_channel_map",
        "messages": [{
            "role": "user",
            "content": (
                "Here is the channel map: Channel 0 shows tcell; Channel 1 shows "
                "27t; Channel 2 shows both 27t and mdo; Channel 3 is the dead "
                "signal. For this experiment use Channels 1 and 2 for both 27t "
                "and mdo, Channel 0 for tcell, and Channel 3 for dead. Please set "
                "the APOC image channel inputs."
            ),
        }],
        "context": _context(
            "segmentation", controls, active_cell_type="27t",
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {
                    "training_data_loaded": True,
                    "cell_types": apoc_cells,
                },
            },
        ),
        "check": _check_apoc_channel_map,
    }


def _apoc_channels_wait_for_training_data_case() -> dict:
    return {
        "name": "apoc_channels_wait_for_training_data",
        "messages": [{
            "role": "user",
            "content": "How do I choose the APOC image channel inputs?",
        }],
        "context": _context(
            "segmentation", [], active_cell_type="27t",
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {
                    "training_data_loaded": False,
                    "cell_types": {
                        "27t": {
                            "available_input_channels": [],
                            "selected_input_channels": [],
                            "channel_controls_ready": False,
                            "trained_classifier_found": False,
                        },
                    },
                },
            },
            current_log={
                "source": "segmentation",
                "recent_lines": [
                    "APOC training controls locked. Load training data to enable classifier training.",
                ],
                "errors": [],
                "has_explicit_error": False,
            },
        ),
        "check": _check_apoc_channels_wait_for_training_data,
    }


def _apoc_feature_preset_case() -> dict:
    controls = [_control(
        "segmentation.apoc.tcell.feature_preset",
        "tcell: APOC feature scale preset",
        "Large structures",
        choices=[
            "Small structures", "Medium structures",
            "Large structures", "Custom feature selection",
        ],
        method="APOC",
        cell_type="tcell",
    )]
    return {
        "name": "apoc_tunes_classifier_features",
        "messages": [{
            "role": "user",
            "content": (
                "These are individual T cells. Please tune the APOC features using "
                "the normal preset for small structures. Do not change minimum size "
                "or the segmentation thresholds."
            ),
        }],
        "context": _context(
            "segmentation", controls, active_cell_type="tcell",
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {
                    "training_data_loaded": True,
                    "cell_types": {
                        "tcell": {
                            "feature_preset": "large_preset",
                            "channel_controls_ready": True,
                            "trained_classifier_found": False,
                        },
                    },
                },
            },
        ),
        "check": _check_apoc_feature_preset,
    }


def _swapped_channel_metadata_case() -> dict:
    return {
        "name": "swapped_channels_stay_out_of_metadata",
        "messages": [{
            "role": "user",
            "content": (
                "I have an experiment where two replicates have the immune channels "
                "swapped. Walk me through where I name them and set the channels."
            ),
        }],
        "context": _context(
            "data_preparation", [],
            metadata_builder={"open": True, "sample_forms_created": True},
        ),
        "check": _check_swapped_channel_metadata,
    }


def _apoc_invalid_channel_index_case() -> dict:
    metadata = {
        "loaded": True,
        "records": _metadata_records(),
        "image_dimensions": [{
            "sample_name": "Movie1",
            "shape": "(50, 4, 20, 512, 512)",
            "dimension_order": "TCZYX",
            "channel_count": 4,
            "timepoint_count": 50,
        }],
        "validation": [],
    }
    return {
        "name": "apoc_rejects_invalid_channel_index",
        "messages": [{
            "role": "user",
            "content": (
                "13T is ch1, blue ch0, green ch4. "
                "Which channels do I pick for APOC?"
            ),
        }],
        "context": _context(
            "segmentation", [], metadata=metadata,
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_invalid_channel_index,
    }


def _apoc_dead_channel_case() -> dict:
    return {
        "name": "apoc_does_not_use_dead_as_negative_class",
        "messages": [{
            "role": "user",
            "content": (
                "13T is ch1, blue ch0, green ch3, dead ch2. "
                "Which channels do I pick for APOC?"
            ),
        }],
        "context": _context(
            "segmentation", [],
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_dead_channel,
    }


def _apoc_feature_grid_case() -> dict:
    controls = [
        _control(
            "segmentation.apoc.tcell.feature_scales",
            "tcell: APOC custom feature scales",
            "1, 2, 5, 10",
            method="APOC", cell_type="tcell", unit="pixels",
        ),
        _control(
            "segmentation.apoc.tcell.feature_filters",
            "tcell: APOC custom feature filters",
            ["Gaussian blur at sigma 1 px"],
            choices=[
                "Gaussian blur at sigma 1 px",
                "Difference of Gaussians at sigma 2 px",
                "Laplacian of Gaussian at sigma 2 px",
                "Sobel edge at sigma 2 px",
            ],
            method="APOC", cell_type="tcell",
        ),
    ]
    return {
        "name": "apoc_feature_grid_is_exposed",
        "messages": [{
            "role": "user",
            "content": "Can I tune the actual APOC feature values and filter sigmas?",
        }],
        "context": _context(
            "segmentation", controls, active_cell_type="tcell",
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_feature_grid,
    }


def _apoc_tune_features_explanation_case() -> dict:
    controls = [
        _control(
            "segmentation.apoc.27t.feature_scales",
            "27t: APOC custom feature scales",
            "0.3, 0.5, 1, 2, 3, 4, 5, 10, 15, 25",
            method="APOC", cell_type="27t", unit="pixels",
        ),
        _control(
            "segmentation.apoc.27t.feature_filters",
            "27t: APOC custom feature filters",
            ["Gaussian blur at sigma 1 px"],
            method="APOC", cell_type="27t",
        ),
    ]
    return {
        "name": "apoc_tune_features_stays_in_segmentation",
        "messages": [{
            "role": "user",
            "content": "Explain to me how to tune the features.",
        }],
        "context": _context(
            "segmentation", controls, active_cell_type="27t",
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_tune_features_explanation,
    }


def _apoc_mdo_feature_recommendation_case() -> dict:
    controls = [
        _control(
            "segmentation.apoc.mdo.feature_preset",
            "mdo: APOC feature scale preset",
            "Medium structures",
            method="APOC", cell_type="mdo",
        ),
        _control(
            "segmentation.apoc.mdo.feature_scales",
            "mdo: APOC custom feature scales",
            "0.3, 0.5, 1, 2, 3, 4, 5, 10, 15, 25",
            method="APOC", cell_type="mdo", unit="pixels",
        ),
        _control(
            "segmentation.apoc.mdo.feature_filters",
            "mdo: APOC custom feature filters",
            [],
            method="APOC", cell_type="mdo",
        ),
    ]
    metadata = {
        "loaded": True,
        "records": _metadata_records(),
        "cell_types": {"organoid": ["27t", "mdo"], "immune": ["tcell"]},
        "validation": [],
    }
    return {
        "name": "apoc_recommends_tune_panel_for_mdo",
        "messages": [
            {
                "role": "user",
                "content": "Explain how to tune the APOC features.",
            },
            {
                "role": "assistant",
                "content": "The Tune Features panel controls classifier filters.",
            },
            {
                "role": "user",
                "content": "What do you recommend for the mdo?",
            },
        ],
        "context": _context(
            "segmentation", controls, metadata=metadata, active_cell_type="mdo",
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_mdo_feature_recommendation,
    }


def _apoc_fill_organoid_features_case() -> dict:
    controls = []
    for cell_type in ("27t", "mdo"):
        controls.extend([
            _control(
                f"segmentation.apoc.{cell_type}.feature_preset",
                f"{cell_type}: APOC feature scale preset",
                "Medium structures",
                choices=[
                    "Small structures", "Medium structures",
                    "Large structures", "Custom feature selection",
                ],
                method="APOC", cell_type=cell_type,
            ),
            _control(
                f"segmentation.apoc.{cell_type}.show_feature_tuning",
                f"{cell_type}: show APOC custom feature tuning",
                False,
                method="APOC", cell_type=cell_type,
            ),
        ])
    metadata = {
        "loaded": True,
        "records": _metadata_records(),
        "cell_types": {"organoid": ["27t", "mdo"], "immune": ["tcell"]},
        "validation": [],
    }
    return {
        "name": "apoc_fills_organoid_tune_features",
        "messages": [{
            "role": "user",
            "content": "Fill in the correct features for the MDO and the 27T for me.",
        }],
        "context": _context(
            "segmentation", controls, metadata=metadata, active_cell_type="mdo",
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_apoc_fill_organoid_features,
    }


def _tracking_which_method_case() -> dict:
    return {
        "name": "tracking_which_method_stays_on_tab",
        "messages": [{"role": "user", "content": "Which method?"}],
        "context": _context(
            "tracking", [], active_cell_type="tcell",
        ),
        "check": _check_tracking_which_method,
    }


def _btrack_step2_enable_case() -> dict:
    controls = [
        _control(
            "tracking.tcell.btrack.use_global_optimization",
            "tcell: Enable global track optimization",
            False, method="btrack", cell_type="tcell",
        ),
        _control(
            "tracking.tcell.btrack.distance_threshold",
            "tcell: btrack optimizer distance threshold",
            60, method="btrack", cell_type="tcell", unit="um",
            visible=False, enabled=False,
        ),
        _control(
            "tracking.tcell.btrack.time_threshold",
            "tcell: btrack optimizer time threshold",
            3, method="btrack", cell_type="tcell", unit="frames",
            visible=False, enabled=False,
        ),
    ]
    return {
        "name": "btrack_step2_is_enabled_after_step1",
        "messages": [
            {"role": "user", "content": "Set up tracking for the T cells."},
            {"role": "assistant", "content": "The Step 1 search radius is set."},
            {"role": "user", "content": "I mean the Step 2 of tracking."},
        ],
        "context": _context(
            "tracking", controls, active_cell_type="tcell",
        ),
        "check": _check_btrack_step2_enable,
    }


def _btrack_step2_tune_case() -> dict:
    controls = [
        _control(
            "tracking.tcell.btrack.use_global_optimization",
            "tcell: Enable global track optimization",
            True, method="btrack", cell_type="tcell",
        ),
        _control(
            "tracking.tcell.btrack.distance_threshold",
            "tcell: btrack optimizer distance threshold",
            60, method="btrack", cell_type="tcell", unit="um",
        ),
        _control(
            "tracking.tcell.btrack.time_threshold",
            "tcell: btrack optimizer time threshold",
            3, method="btrack", cell_type="tcell", unit="frames",
        ),
        _control(
            "tracking.tcell.btrack.hypotheses",
            "tcell: btrack optimization hypotheses",
            ["P_FP", "P_init", "P_term", "P_link"],
            choices=[
                "P_FP", "P_init", "P_term", "P_link",
                "P_branch", "P_dead", "P_merge",
            ],
            method="btrack", cell_type="tcell",
        ),
    ]
    return {
        "name": "btrack_step2_uses_measured_gap",
        "messages": [{
            "role": "user",
            "content": "For Step 2, bridge gaps up to 4 frames and 40 um.",
        }],
        "context": _context(
            "tracking", controls, active_cell_type="tcell",
        ),
        "check": _check_btrack_step2_tune,
    }


def _segmentation_minimum_size_case() -> dict:
    controls = [_control(
        "segmentation.apoc.tcell.minimum_size",
        "tcell: APOC minimum object size",
        10,
        method="APOC", cell_type="tcell", unit="voxels",
    )]
    metadata = {
        "loaded": True,
        "records": [{
            "sample_name": "Movie1",
            "pixel_distance_xy": 1.0,
            "pixel_distance_z": 2.0,
            "distance_unit": "um",
            "time_interval": 2,
            "time_unit": "min",
        }],
        "cell_types": {"organoid": [], "immune": ["tcell"]},
        "validation": [],
    }
    return {
        "name": "segmentation_minimum_size_uses_cell_volume",
        "messages": [{
            "role": "user",
            "content": (
                "My T cells are approximately 10 um across. "
                "Set the segmentation Minimum size."
            ),
        }],
        "context": _context(
            "segmentation", controls, metadata=metadata,
            active_cell_type="tcell",
            segmentation={"method": "APOC (GPU)"},
        ),
        "check": _check_segmentation_minimum_size,
    }


def _mask_edt_direction_case() -> dict:
    strategy = "APOC Mask + EDT/Watershed Resegmentation"
    return {
        "name": "mask_edt_direction_and_fallback",
        "messages": [
            {
                "role": "user",
                "content": "I used EDT 50 for organoids before.",
            },
            {
                "role": "assistant",
                "content": "A higher EDT makes splitting harder.",
            },
            {
                "role": "user",
                "content": (
                    "That contradicts the implementation. Explain the EDT direction "
                    "properly."
                ),
            },
        ],
        "context": _context(
            "segmentation", [], active_cell_type="13T",
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {"cell_type_strategies": {"13T": strategy}},
            },
        ),
        "check": _check_mask_edt_direction,
    }


def _contact_and_dead_threshold_case() -> dict:
    controls = [_control(
        "features.tcell.contact_distance",
        "tcell: contact distance",
        1.01,
        cell_type="tcell", unit="um",
    )]
    metadata = {
        "loaded": True,
        "records": [{
            "sample_name": "Movie1",
            "pixel_distance_xy": 1.01,
            "pixel_distance_z": 1.05,
            "time_interval": 2,
            "time_unit": "min",
        }],
        "validation": [],
    }
    return {
        "name": "contact_distance_uses_xy_pixel_scale",
        "messages": [{
            "role": "user",
            "content": (
                "How can I set the contact and dead-mask percentage threshold "
                "correctly? Is 1.01 um strict touching?"
            ),
        }],
        "context": _context(
            "feature_extraction", controls, metadata=metadata,
            active_cell_type="tcell",
        ),
        "check": _check_contact_and_dead_threshold,
    }


def _loaded_metadata_not_unsaved_case() -> dict:
    return {
        "name": "loaded_metadata_is_not_called_unsaved",
        "messages": [{
            "role": "user",
            "content": "The metadata is loaded. Let's do the segmentation.",
        }],
        "context": _context(
            "segmentation", [],
            metadata_builder={
                "open": True,
                "record_source": "loaded_metadata_copy",
                "save_required": False,
            },
            segmentation={"method": "APOC (GPU)", "method_index": 0},
        ),
        "check": _check_loaded_metadata_not_unsaved,
    }


def _external_zarr_reload_case() -> dict:
    return {
        "name": "external_zarr_requires_metadata_reload",
        "messages": [{
            "role": "user",
            "content": (
                "I converted the images to Zarr and updated the metadata file "
                "outside BEHAV3D. Load Dataset still shows nothing. What should I do?"
            ),
        }],
        "context": _context("visualization", []),
        "check": _check_external_zarr_reload,
    }


def _missing_log_error_case() -> dict:
    return {
        "name": "failed_load_requests_log_error",
        "messages": [{
            "role": "user",
            "content": (
                "I clicked Generate Training Data but nothing appeared. What is wrong?"
            ),
        }],
        "context": _context(
            "segmentation", [],
            segmentation={
                "method": "APOC (GPU)",
                "apoc": {"training_data_loaded": False},
            },
            current_log={
                "source": "segmentation",
                "recent_lines": [
                    "APOC training pipeline using GPU device: NVIDIA GPU",
                    "Loading training data...",
                    "Running _load_training_images...",
                ],
                "errors": [],
                "has_explicit_error": False,
            },
        ),
        "check": _check_missing_log_error,
    }


def _tracking_guide_case() -> dict:
    choices = [
        "LAP (laptrack)", "TrackPy", "Propagation",
        "Reporter Propagation", "btrack (Bayesian)",
    ]
    controls = [
        _control(
            f"tracking.{cell_type}.method", f"{cell_type}: Tracking method", method,
            choices=choices, cell_type=cell_type,
        )
        for cell_type, method in (
            ("Tcell_HIV", "btrack (Bayesian)"),
            ("Tcell_cmtrm", "btrack (Bayesian)"),
            ("collagen", "LAP (laptrack)"),
        )
    ]
    return {
        "name": "tracking_asks_about_motion",
        "messages": [{"role": "user", "content": "Guide tracking"}],
        "context": _context("tracking", controls, active_cell_type="Tcell_HIV"),
        "check": _check_tracking_guide,
    }


def _stationary_tracking_case() -> dict:
    controls = [_control(
        "tracking.collagen.method", "collagen: Tracking method", "LAP (laptrack)",
        choices=[
            "LAP (laptrack)", "TrackPy", "Propagation",
            "Reporter Propagation", "btrack (Bayesian)",
        ],
        cell_type="collagen",
    )]
    return {
        "name": "stationary_uses_propagation",
        "messages": [
            {"role": "user", "content": "Guide tracking"},
            {"role": "assistant", "content": "How much does collagen move between frames?"},
            {"role": "user", "content": "Collagen stays stationary and overlaps between frames."},
        ],
        "context": _context("tracking", controls, active_cell_type="collagen"),
        "check": _check_stationary_tracking,
    }


def _tracking_radius_case() -> dict:
    controls = [_control(
        "tracking.Tcell_cmtrm.btrack.maximum_search_radius",
        "Tcell cmtrm: btrack maximum search radius", 100,
        method="btrack", cell_type="Tcell_cmtrm", unit="um",
    )]
    return {
        "name": "tracking_radius_uses_interval",
        "messages": [{
            "role": "user",
            "content": (
                "The fastest Tcell_cmtrm cells move about 60 um per minute. "
                "Please adjust their maximum search radius for our 15 second interval."
            ),
        }],
        "context": _context("tracking", controls, active_cell_type="Tcell_cmtrm"),
        "check": _check_tracking_radius,
    }


def _filtering_case() -> dict:
    controls = [
        _control("filtering.Tcell_HIV.minimum_length.enabled", "Filter short tracks", True),
        _control(
            "filtering.Tcell_HIV.minimum_length.timepoints",
            "Minimum track length", 30, unit="timepoints",
        ),
        _control(
            "filtering.Tcell_HIV.maximum_length.enabled",
            "Trim retained tracks to a common length", True,
        ),
        _control(
            "filtering.Tcell_HIV.maximum_length.timepoints",
            "Common output track length", 30, unit="timepoints",
        ),
    ]
    return {
        "name": "equal_track_lengths_are_valid",
        "messages": [{"role": "user", "content": "Review filters"}],
        "context": _context("filtering", controls, active_cell_type="Tcell_HIV"),
        "check": _check_filtering,
    }


def _reporter_propagation_case() -> dict:
    controls = [_control(
        "tracking.calcium_reporter.method",
        "calcium reporter: Tracking method",
        "btrack (Bayesian)",
        choices=[
            "LAP (laptrack)", "TrackPy", "Propagation",
            "Reporter Propagation", "btrack (Bayesian)",
        ],
        cell_type="calcium_reporter",
    )]
    return {
        "name": "static_reporter_uses_reporter_propagation",
        "messages": [{
            "role": "user",
            "content": (
                "These calcium reporter cells do not move or change shape, but their "
                "fluorescence flickers and some frames cannot be segmented. Please set "
                "the appropriate tracking method."
            ),
        }],
        "context": _context(
            "tracking", controls, active_cell_type="calcium_reporter"
        ),
        "check": _check_reporter_propagation,
    }


def _active_killing_case() -> dict:
    controls = [
        _control(
            "features.active_killing.target_types",
            "Active Killing: Target cell type",
            ["organoid1", "organoid2"],
            choices=["organoid1", "organoid2"],
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.observation_window",
            "Active Killing: Observation window",
            3,
            unit="timepoints",
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.death_signal",
            "Active Killing: Death or reporter signal",
            "Dead-mask percentage",
            choices=[
                "Dead-mask percentage", "Mean dead-dye intensity",
                "Dead-mask pixel count",
            ],
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.use_absolute_threshold",
            "Active Killing: Use an absolute signal-increase threshold",
            False,
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.absolute_threshold",
            "Active Killing: Absolute signal-increase threshold",
            0.0,
            unit="pixels",
            method="Active Killing",
            cell_type="tcell",
            active=False,
        ),
    ]
    metadata = {
        "loaded": True,
        "records": [{
            "sample_name": "Movie1",
            "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0,
            "time_interval": 2,
            "time_unit": "min",
        }],
        "validation": [],
    }
    return {
        "name": "active_killing_uses_cadence",
        "messages": [{
            "role": "user",
            "content": (
                "Configure active killing for tcell against organoid1 only. I expect "
                "killing within 10 minutes and images are every 2 minutes. Use an "
                "absolute threshold of 30 dead pixels."
            ),
        }],
        "context": _context(
            "feature_extraction", controls, metadata=metadata,
            active_cell_type="tcell",
            feature_extraction={"active_killing_open": True},
        ),
        "check": _check_active_killing,
    }


def _feature_group_dead_dye_case() -> dict:
    choices = [
        "movement", "intensity", "morphology", "contact",
        "invasiveness", "death",
    ]
    controls = [_control(
        "features.tcell.feature_groups",
        "tcell: feature groups",
        choices,
        choices=choices,
        cell_type="tcell",
        required_choices=["movement", "intensity", "contact", "death"],
    )]
    return {
        "name": "tcell_features_keep_dead_dye_intensity",
        "messages": [{
            "role": "user",
            "content": "Now adjust T cells. Should I drop intensity?",
        }],
        "context": _context(
            "feature_extraction", controls, active_cell_type="tcell",
        ),
        "check": _check_feature_group_dead_dye,
    }


def _active_killing_complete_acceptance_case() -> dict:
    controls = [
        _control(
            "features.active_killing.target_types",
            "Active Killing: Target cell type",
            ["27t", "mdo"],
            choices=["27t", "mdo"],
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.observation_window",
            "Active Killing: Observation window",
            5,
            unit="timepoints",
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.death_signal",
            "Active Killing: Death or reporter signal",
            "Dead-mask percentage",
            choices=[
                "Dead-mask percentage", "Mean dead-dye intensity",
                "Dead-mask pixel count",
            ],
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.use_absolute_threshold",
            "Active Killing: Use an absolute signal-increase threshold",
            False,
            method="Active Killing",
            cell_type="tcell",
        ),
        _control(
            "features.active_killing.absolute_threshold",
            "Active Killing: Absolute signal-increase threshold",
            0.0,
            unit="pixels",
            method="Active Killing",
            cell_type="tcell",
            active=False,
        ),
        _control(
            "features.active_killing.minimum_contact_duration",
            "Active Killing: Minimum contact duration",
            1,
            unit="timepoints",
            method="Active Killing",
            cell_type="tcell",
        ),
    ]
    return {
        "name": "active_killing_accepts_complete_setup",
        "messages": [
            {
                "role": "assistant",
                "content": (
                    "Active Killing configuration for tcell against 27t and mdo: "
                    "Death signal: Dead-mask pixel count. Absolute threshold: "
                    "30 dead pixels. Observation window: 5 timepoints. "
                    "Minimum contact duration: 1 frame."
                ),
            },
            {"role": "user", "content": "Ok, these settings seem ok"},
        ],
        "context": _context(
            "feature_extraction", controls, active_cell_type="tcell",
            feature_extraction={
                "active_killing_open": True,
                "active_killing": {
                    "setup_ready": True,
                    "setup_issues": [],
                },
            },
        ),
        "check": _check_active_killing_complete_acceptance,
    }


def _hmm_movement_controls() -> list[dict]:
    prefix = "analysis.state_classification.tcell."
    return [
        _control(
            prefix + "timepoint_features",
            "tcell: Timepoint features",
            ["speed"],
            choices=[
                "speed", "displacement", "cumulative_displacement",
                "displacement_from_origin", "directional_persistence",
                "median_turning_angle", "mean_dead_dye",
            ],
            method="HMM",
            cell_type="tcell",
        ),
        _control(
            prefix + "window_features",
            "tcell: Additional window features",
            ["net_displacement"],
            choices=[
                "net_displacement", "straightness",
                "mean_square_displacement",
            ],
            method="HMM",
            cell_type="tcell",
        ),
        _control(
            prefix + "binary_feature_groups",
            "tcell: Binary feature groups",
            ["27t_contact", "mdo_contact"],
            choices=["27t_contact", "mdo_contact"],
            method="HMM",
            cell_type="tcell",
        ),
    ]


def _hmm_movement_options_case() -> dict:
    controls = _hmm_movement_controls()
    return {
        "name": "hmm_lists_all_movement_options",
        "messages": [{
            "role": "user",
            "content": "Only movement features of relevance",
        }],
        "context": _context(
            "analysis", controls, active_cell_type="tcell",
            analysis={"view": "behavioral_state", "selected_cell_type": "tcell"},
        ),
        "check": _check_hmm_movement_options,
    }


def _hmm_apply_all_movement_case() -> dict:
    return {
        "name": "hmm_applies_all_offered_movement_features",
        "messages": [
            {
                "role": "assistant",
                "content": (
                    "Choose the feature names you want, or say use all available "
                    "movement features."
                ),
            },
            {"role": "user", "content": "Use all available movement features"},
        ],
        "context": _context(
            "analysis", _hmm_movement_controls(), active_cell_type="tcell",
            analysis={"view": "behavioral_state", "selected_cell_type": "tcell"},
        ),
        "check": _check_hmm_apply_all_movement,
    }


def _active_killing_zero_threshold_readiness_case() -> dict:
    return {
        "name": "active_killing_zero_threshold_is_not_ready",
        "messages": [
            {"role": "assistant", "content": "Active Killing setup"},
            {"role": "user", "content": "Is it ready?"},
        ],
        "context": _context(
            "feature_extraction", [], active_cell_type="tcell",
            feature_extraction={
                "active_killing_open": True,
                "active_killing": {
                    "setup_ready": False,
                    "setup_issues": [
                        "Absolute signal-increase threshold must be greater than 0."
                    ],
                },
            },
        ),
        "check": _check_active_killing_zero_threshold_readiness,
    }


def _hmm_single_frame_case() -> dict:
    controls = [
        _control(
            "analysis.state_classification.tcell.window_size",
            "tcell: Window size", 5, unit="timepoints", method="HMM",
            cell_type="tcell",
        ),
        _control(
            "analysis.state_classification.tcell.smooth_window",
            "tcell: Feature smoothing window", 5, unit="timepoints",
            method="HMM", cell_type="tcell",
        ),
        _control(
            "analysis.state_classification.tcell.start_offset",
            "tcell: Start offset", 1, unit="timepoints", method="HMM",
            cell_type="tcell",
        ),
    ]
    return {
        "name": "hmm_single_frame_events",
        "messages": [{
            "role": "user",
            "content": (
                "My calcium events are single-timepoint peaks. Please set the HMM "
                "windows appropriately and leave the start offset at its recommendation."
            ),
        }],
        "context": _context(
            "analysis", controls, active_cell_type="tcell",
            analysis={"view": "behavioral_state", "selected_cell_type": "tcell"},
        ),
        "check": _check_hmm_single_frame,
    }


def _trajectory_linkage_case() -> dict:
    controls = [_control(
        "analysis.state_trajectory.tcell.linkage",
        "tcell: Agglomerative linkage",
        "single",
        choices=["average", "complete", "single"],
        method="State Trajectory",
        cell_type="tcell",
    )]
    return {
        "name": "trajectory_uses_average_linkage",
        "messages": [{
            "role": "user",
            "content": "Set the recommended linkage for State Trajectory clustering.",
        }],
        "context": _context(
            "analysis", controls, active_cell_type="tcell",
            analysis={"view": "state_trajectory", "selected_cell_type": "tcell"},
        ),
        "check": _check_trajectory_linkage,
    }


def _functional_experiment_context_case() -> dict:
    reference = {
        "notes": [{
            "source": "README_BEHAV3D_FUNC_MUC1_Exp085.md",
            "text": (
                "Exp085 is an isogenic functional comparison. organoid1 is 27T_V2_KO "
                "(MUC1 knockout) and organoid2 is 27T_V2_OE (MUC1 rescue). Both are "
                "present in each well with the same T-cell product, so the primary "
                "comparison is paired KO versus rescue within each well. NCAM+ has "
                "3 wells and NCAM- has 5, so that axis is exploratory. Invasiveness "
                "was not computed in this experiment."
            ),
        }],
        "saved_configurations": [{
            "source": "behav3d_parameters_FUNC_clean.yml",
            "settings": {
                "features": {
                    "organoid1": {
                        "features_choice": [
                            "intensity", "morphology", "contact", "death",
                        ],
                    },
                    "organoid2": {
                        "features_choice": [
                            "intensity", "morphology", "contact", "death",
                        ],
                    },
                },
            },
        }],
        "configuration_caveat": "Saved settings are not proof that a module ran.",
    }
    return {
        "name": "functional_experiment_uses_paired_design",
        "messages": [{
            "role": "user",
            "content": (
                "What is the cleanest comparison for the MUC1 question, and can I "
                "interpret invasiveness from this experiment?"
            ),
        }],
        "context": _context(
            "analysis", [], experiment_reference=reference, results=[],
        ),
        "check": _check_functional_experiment_context,
    }


def _safety_profiling_context_case() -> dict:
    reference = {
        "notes": [{
            "source": "README_BEHAV3D_SafetyProfiling_Exp010.md",
            "text": (
                "Exp010 is multi-organoid safety profiling. In combined wells, tumor "
                "27T and healthy MDO share the same well with TEG cells. Active killing "
                "is a contact-associated relative rise: percentage_dead_mask reaches "
                "at least 1.5 times baseline within 5 frames (10 minutes), after at "
                "least one frame of contact. There is one well per control and two "
                "combined wells, so comparisons are descriptive/exploratory."
            ),
        }],
        "saved_configurations": [{
            "source": "behav3d_parameters_multiorganoid_clean.yml",
            "settings": {
                "features": {
                    "TEG": {
                        "features_choice": [
                            "movement", "contact", "invasiveness", "death",
                        ],
                    },
                },
            },
        }],
    }
    return {
        "name": "safety_profiling_preserves_operational_definition",
        "messages": [{
            "role": "user",
            "content": (
                "How should I frame the safety comparison, and what exactly does "
                "active killing mean here?"
            ),
        }],
        "context": _context(
            "analysis", [],
            metadata={
                "loaded": True,
                "records": [{
                    "sample_name": "Exp010",
                    "pixel_distance_xy": 0.5,
                    "pixel_distance_z": 2.0,
                    "time_interval": 2,
                    "time_unit": "min",
                }],
                "validation": [],
            },
            experiment_reference=reference, results=[],
        ),
        "check": _check_safety_profiling_context,
    }


def _historical_btrack_examples_case() -> dict:
    controls = [
        _control(
            "tracking.tcell.btrack.maximum_search_radius",
            "tcell: Maximum search radius", 25,
            method="btrack", cell_type="tcell", unit="um",
        ),
        _control(
            "tracking.tcell.btrack.distance_threshold",
            "tcell: Distance threshold", 30,
            method="btrack", cell_type="tcell", unit="um",
        ),
        _control(
            "tracking.tcell.btrack.time_threshold",
            "tcell: Time threshold", 2,
            method="btrack", cell_type="tcell", unit="frames",
        ),
    ]
    return {
        "name": "historical_btrack_values_are_examples_only",
        "messages": [{
            "role": "user",
            "content": (
                "Do you have example values from a previous experiment for btrack "
                "on T cells? Can I copy them into this experiment?"
            ),
        }],
        "context": _context(
            "tracking", controls, active_cell_type="tcell",
            metadata={
                "loaded": True,
                "records": [{
                    "sample_name": "Current01",
                    "pixel_distance_xy": 0.6,
                    "pixel_distance_z": 2.0,
                    "time_interval": 30,
                    "time_unit": "s",
                }],
                "cell_types": {"immune": ["tcell"]},
                "validation": [],
            },
            tracking={"method": "btrack"},
        ),
        "check": _check_historical_btrack_examples,
    }


def _historical_calcium_example_case() -> dict:
    return {
        "name": "historical_calcium_profile_explains_reporter_propagation",
        "messages": [{
            "role": "user",
            "content": (
                "Was there a previous experiment with near-static calcium reporter "
                "cells? What method and example values did it use?"
            ),
        }],
        "context": _context(
            "tracking", [], active_cell_type="islets",
            metadata={
                "loaded": True,
                "records": [{
                    "sample_name": "NewCalcium",
                    "pixel_distance_xy": 0.5,
                    "pixel_distance_z": 2.5,
                    "time_interval": 10,
                    "time_unit": "s",
                }],
                "cell_types": {"other": ["islets"]},
                "validation": [],
            },
        ),
        "check": _check_historical_calcium_example,
    }


def _experiment_reference_metadata_conflict_case() -> dict:
    reference = {
        "notes": [{
            "source": "README_BEHAV3D_IVM_HIV.md",
            "text": "The CS002 frame interval was 10 seconds.",
        }],
        "source_policy": (
            "Use live metadata for acquisition facts. State discrepancies."
        ),
    }
    return {
        "name": "live_metadata_wins_reference_cadence_conflict",
        "messages": [{
            "role": "user",
            "content": (
                "The experiment README says CS002 used 10-second frames. What does "
                "the currently loaded metadata say, and which value should I use?"
            ),
        }],
        "context": _context(
            "data_preparation", [],
            metadata={
                "loaded": True,
                "records": [{
                    "sample_name": "CS002",
                    "pixel_distance_xy": 1.15,
                    "pixel_distance_z": 4.0,
                    "time_interval": 15,
                    "time_unit": "s",
                }],
                "validation": [],
            },
            experiment_reference=reference,
        ),
        "check": _check_experiment_reference_metadata_conflict,
    }


def _organoid_line_grouping_case() -> dict:
    return {
        "name": "organoid_lines_require_processing_choice",
        "messages": [{
            "role": "user",
            "content": (
                "I co-culture T cells and macrophages with three different organoid "
                "lines. Can you help me set up the metadata?"
            ),
        }],
        "context": _context(
            "data_preparation", [],
            metadata={"loaded": False, "records": [], "validation": []},
            metadata_builder={"open": False, "sample_forms_created": False},
        ),
        "check": _check_organoid_line_grouping,
    }


def _metadata_identifier_confirmation_case() -> dict:
    return {
        "name": "metadata_identifiers_require_confirmation",
        "messages": [{
            "role": "user",
            "content": (
                "No well identifier. M21 and M23 are macrophage lines, and "
                "T-cells are GD2_CART."
            ),
        }],
        "context": _context(
            "data_preparation", [],
            metadata={
                "loaded": True,
                "record_source": "metadata_builder_draft",
                "records": [{"sample_name": "Img001_DO7_GD2"}],
                "validation": [],
                "save_required": True,
            },
            metadata_builder={
                "open": True,
                "sample_forms_created": True,
                "actions": {"save_available": False, "load_available": False},
            },
        ),
        "check": _check_metadata_identifier_confirmation,
    }


def _metadata_time_conversion_case() -> dict:
    controls = []
    records = []
    for index in range(8):
        records.append({
            "sample_name": f"Img{index + 1:03d}",
            "time_interval": 2,
            "time_unit": "s",
        })
        controls.extend([
            _control(
                f"metadata.samples.{index}.time_interval",
                f"Sample {index + 1}: Time interval", 2.0,
            ),
            _control(
                f"metadata.samples.{index}.time_unit",
                f"Sample {index + 1}: Time unit", "s",
                choices=["s", "min", "h"],
            ),
        ])
    return {
        "name": "metadata_minutes_are_converted_for_all_samples",
        "messages": [{
            "role": "user",
            "content": "The time unit is s while I set 2 minutes.",
        }],
        "context": _context(
            "data_preparation", controls,
            metadata={
                "loaded": True,
                "record_source": "metadata_builder_draft",
                "records": records,
                "validation": [],
                "save_required": True,
            },
            metadata_builder={
                "open": True,
                "sample_forms_created": True,
                "actions": {"save_available": False, "load_available": False},
            },
        ),
        "check": _check_metadata_time_conversion,
    }


def _metadata_completion_case() -> dict:
    validation = [
        {
            "severity": "error",
            "field": "well",
            "message": "Img001: mandatory well is missing.",
        },
        {
            "severity": "error",
            "field": "T-cells.line",
            "message": "Img001: mandatory T-cells line is missing.",
        },
        {
            "severity": "error",
            "field": "Macrophages.line",
            "message": "Img001: mandatory Macrophages line is missing.",
        },
    ]
    return {
        "name": "metadata_completion_reports_real_mandatory_fields",
        "messages": [{
            "role": "user",
            "content": "Is this all that is needed?",
        }],
        "context": _context(
            "data_preparation", [],
            metadata={
                "loaded": True,
                "record_source": "metadata_builder_draft",
                "records": [{"sample_name": "Img001"}],
                "validation": validation,
                "save_required": True,
            },
            metadata_builder={
                "open": True,
                "sample_forms_created": True,
                "draft_validation": validation,
                "actions": {"save_available": False, "load_available": False},
            },
        ),
        "check": _check_metadata_completion,
    }


def _metadata_save_case() -> dict:
    return {
        "name": "metadata_save_uses_real_action",
        "messages": [{"role": "user", "content": "Yeah save it please"}],
        "context": _context(
            "data_preparation", [],
            metadata={
                "loaded": True,
                "record_source": "metadata_builder_draft",
                "records": [{"sample_name": "Img001"}],
                "validation": [],
                "save_required": True,
            },
            metadata_builder={
                "open": True,
                "sample_forms_created": True,
                "draft_validation": [],
                "actions": {
                    "save_available": True,
                    "save_also_loads": True,
                    "load_available": False,
                },
            },
        ),
        "check": _check_metadata_save,
    }


def _choose_analysis_case() -> dict:
    return {
        "name": "choose_analysis_explains_options",
        "messages": [{"role": "user", "content": "Choose analysis"}],
        "context": _context(
            "analysis", [],
            assistant_session={"intent": "choose_analysis"},
            analysis={"view": "death_dynamics"},
        ),
        "check": _check_choose_analysis,
    }


def _open_death_dynamics_case() -> dict:
    return {
        "name": "death_dynamics_opens_specific_view",
        "messages": [{
            "role": "user",
            "content": "Can you take me to Death Dynamics?",
        }],
        "context": _context(
            "analysis", [],
            analysis={"view": "behavioral_state"},
        ),
        "check": _check_open_death_dynamics,
    }


def _tool_calls(result: dict, name: str) -> list[dict]:
    return [call.get("arguments", {}) for call in result["calls"] if call.get("name") == name]


def _check_organoid_line_grouping(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "separate organoid types" not in text or "one organoid type" not in text:
        errors.append("did not ask how organoid lines should map to processing types")
    if "same movie" not in text or "segmented or tracked separately" not in text:
        errors.append("did not explain when separate organoid types are appropriate")
    if result["calls"]:
        errors.append("attempted metadata edits before the grouping choice was confirmed")
    return errors


def _check_metadata_identifier_confirmation(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in ("well", "mandatory", "condition is optional", "m21/m23", "none"):
        if phrase not in text:
            errors.append(f"missing identifier guidance: {phrase}")
    if result["calls"]:
        errors.append("filled filename-derived values before confirmation")
    return errors


def _check_metadata_time_conversion(result: dict) -> list[str]:
    calls = _tool_calls(result, "set_ui_value")
    changed = {call.get("control_id"): call.get("value") for call in calls}
    errors = []
    if "120 seconds" not in result["text"].lower():
        errors.append("did not explain the 2 minute to 120 second conversion")
    for index in range(8):
        control_id = f"metadata.samples.{index}.time_interval"
        if changed.get(control_id) != 120:
            errors.append(f"missing {control_id} -> 120")
    if any(str(call.get("control_id", "")).startswith("data.output") for call in calls):
        errors.append("emitted an unrelated output-directory action")
    return errors


def _check_metadata_completion(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "mandatory well", "mandatory t-cells line",
        "mandatory macrophages line", "condition", "optional", "none",
    ):
        if phrase not in text:
            errors.append(f"missing completeness detail: {phrase}")
    if _tool_calls(result, "save_metadata"):
        errors.append("offered save despite mandatory validation errors")
    return errors


def _check_metadata_save(result: dict) -> list[str]:
    errors = []
    if _tool_calls(result, "save_metadata") != [{}]:
        errors.append("did not call the real save metadata action")
    if _tool_calls(result, "load_metadata"):
        errors.append("offered a redundant load after save")
    if "confirm" not in result["text"].lower():
        errors.append("did not tell the user that the write requires confirmation")
    return errors


def _check_choose_analysis(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in ("death dynamics", "behavioral state", "state trajectory"):
        if phrase not in text:
            errors.append(f"did not explain {phrase}")
    if result["calls"]:
        errors.append("Choose analysis navigated instead of explaining options")
    return errors


def _check_open_death_dynamics(result: dict) -> list[str]:
    calls = _tool_calls(result, "open_analysis_view")
    if calls != [{"view": "death_dynamics"}]:
        return [f"expected Death Dynamics subview action, got {result['calls']!r}"]
    if _tool_calls(result, "navigate_to_step"):
        return ["used generic Analysis navigation instead of the named subview"]
    return []


def _check_metadata_setup(result: dict) -> list[str]:
    errors = []
    bulk = _tool_calls(result, "bulk_fill_metadata")
    if not bulk:
        return errors + ["did not propose bulk metadata population"]
    args = bulk[0]
    expected = {"n_samples": 3, "n_immune": 2, "n_other": 1}
    for key, value in expected.items():
        if args.get(key) != value:
            errors.append(f"{key} was {args.get(key)!r}, expected {value}")
    samples = args.get("samples") or []
    if len(samples) != 3:
        errors.append(f"created {len(samples)} sample objects instead of 3")
    for index, sample in enumerate(samples):
        for key, value in (("pixel_distance_xy", 1.15), ("pixel_distance_z", 4)):
            if sample.get(key) != value:
                errors.append(f"sample {index + 1} {key} was not populated with {value}")
        if sample.get("time_interval") not in (None, 15):
            errors.append(
                f"sample {index + 1} time_interval was "
                f"{sample.get('time_interval')!r}, expected 15 or blank"
            )
        for unknown in ("dimension_order", "time_unit", "raw_image_path", "sample_name", "well"):
            if sample.get(unknown) not in (None, ""):
                errors.append(f"sample {index + 1} invented {unknown}={sample.get(unknown)!r}")
    immune_names = [str(item).lower() for item in (args.get("immune_names") or [])]
    other_names = [str(item).lower() for item in (args.get("other_names") or [])]
    if len(immune_names) != 2:
        errors.append(f"created {len(immune_names)} immune names instead of 2")
    if not any("hiv" in name and ("t" in name or "cell" in name) for name in immune_names):
        errors.append("immune names did not include an HIV-specific T-cell population")
    if not any(
        ("cmtmr" in name) or ("t" in name and "hiv" not in name)
        for name in immune_names
    ):
        errors.append("immune names did not include the uninfected T-cell population")
    if not any("collagen" in name for name in other_names):
        errors.append("other cell-type names did not include collagen")
    return errors


def _check_pixel_fill(result: dict) -> list[str]:
    calls = _tool_calls(result, "set_ui_value")
    changed = {call.get("control_id"): call.get("value") for call in calls}
    errors = []
    for index in range(3):
        expected = {
            f"metadata.samples.{index}.pixel_distance_xy": 1.15,
            f"metadata.samples.{index}.pixel_distance_z": 4,
        }
        for control_id, value in expected.items():
            if changed.get(control_id) != value:
                errors.append(f"missing {control_id} -> {value}")
    return errors


def _check_segmentation(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    bad_directions = (
        "lower the mask", "lower mask", "decrease the mask", "decrease mask",
        "lower the seed", "lower seed", "decrease the seed", "decrease seed",
    )
    if any(phrase in text for phrase in bad_directions):
        errors.append("response still recommends lowering Mask or Seed threshold")
    if "requiring less confidence" in text:
        errors.append("incorrectly described a higher Seed threshold as requiring less confidence")
    calls = _tool_calls(result, "set_ui_value")
    values = {call.get("control_id"): call.get("value") for call in calls}
    mask = values.get("segmentation.apoc.Tcell_cmtrm.mask_threshold")
    seed = values.get("segmentation.apoc.Tcell_cmtrm.seed_threshold")
    if mask is not None and float(mask) <= 0.4:
        errors.append(f"Mask threshold action did not increase the value: {mask}")
    if seed is not None and float(seed) <= 0.7:
        errors.append(f"Seed threshold action did not increase the value: {seed}")
    if not calls and not any(word in text for word in ("raise", "increase", "higher")):
        errors.append("response does not advise increasing the active thresholds")
    return errors


def _check_method_requires_signal_context(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if not (
        ("bleed" in text or "same channel" in text)
        and ("clean" in text or "isolated" in text)
        and "?" in result["text"]
    ):
        errors.append(
            "did not ask whether the target is isolated/clean or has bleed-through"
        )
    if result["calls"]:
        errors.append("changed a method before the signal context was known")
    if "gfp" in text or "rfp" in text:
        errors.append("asked for irrelevant fluorophore labels instead of visible cell signals")
    if any(phrase in text for phrase in (
        "solid default", "good default", "recommend apoc",
        "best choice is apoc", "apoc is the practical",
    )):
        errors.append("characterized APOC as preferred before signal context was known")
    if any(phrase in text for phrase in (
        "yes, very likely", "best choice is cellpose",
        "recommend cellpose-sam", "my recommendation is cellpose",
    )):
        errors.append("recommended Cellpose-SAM before bleed-through was confirmed")
    return errors


def _check_confirmed_bleed_through(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "apoc" not in text:
        errors.append("did not recommend APOC after bleed-through was confirmed")
    if any(phrase in text for phrase in (
        "recommend cellpose-sam", "cellpose-sam is the best",
        "cellpose-sam is a good fit", "cellpose-sam is the right choice",
    )):
        errors.append("still recommended Cellpose-SAM despite bleed-through")
    return errors


def _check_apoc_channel_map(result: dict) -> list[str]:
    changed = _changed_values(result)
    expected = {
        "segmentation.apoc.27t.input_channels": ["Channel 1", "Channel 2"],
        "segmentation.apoc.mdo.input_channels": ["Channel 1", "Channel 2"],
        "segmentation.apoc.tcell.input_channels": ["Channel 0"],
        "segmentation.apoc.dead.input_channels": ["Channel 3"],
    }
    errors = []
    for control_id, value in expected.items():
        if changed.get(control_id) != value:
            errors.append(
                f"{control_id} was {changed.get(control_id)!r}, expected {value!r}"
            )
    text = result["text"].lower()
    if "apoc does not have" in text or "uses all" in text:
        errors.append("incorrectly denied APOC per-model channel selection")
    return errors


def _check_apoc_channels_wait_for_training_data(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "generate training data" not in text:
        errors.append("did not explain that Generate Training Data unlocks channels")
    if any(phrase in text for phrase in (
        "apoc does not have", "no channel selection", "always uses all",
        "switch to cellpose", "need to switch to cellpose",
    )):
        errors.append("misrepresented unavailable pre-load APOC channel controls")
    if result["calls"]:
        errors.append("attempted to set channel controls before they existed")
    return errors


def _check_apoc_feature_preset(result: dict) -> list[str]:
    changed = _changed_values(result)
    errors = []
    expected_id = "segmentation.apoc.tcell.feature_preset"
    if changed.get(expected_id) != "Small structures":
        errors.append(
            f"feature preset was {changed.get(expected_id)!r}, expected Small structures"
        )
    forbidden = (
        "minimum_size", "minimum size", "mask_threshold", "mask threshold",
        "seed_threshold", "seed threshold",
    )
    for call in _tool_calls(result, "set_ui_value"):
        control_id = str(call.get("control_id") or "").lower()
        if any(token.replace(" ", "_") in control_id for token in forbidden):
            errors.append(f"changed post-processing instead of features: {control_id}")
    return errors


def _check_swapped_channel_metadata(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "metadata builder does not map raw channel indices",
        "line", "processing slots", "segmentation method",
    ):
        if phrase not in text:
            errors.append(f"missing swapped-channel boundary guidance: {phrase}")
    if any(phrase in text for phrase in (
        "each sample form has a channel", "same for all", "per-sample dropdown",
    )):
        errors.append("invented metadata channel-mapping controls")
    if result["calls"]:
        errors.append("attempted edits before the slot workflow was confirmed")
    return errors


def _check_apoc_invalid_channel_index(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "indexed 0-3" not in text or "conflicts with the data" not in text:
        errors.append("did not catch Channel 4 against a four-channel image")
    if any(phrase in text for phrase in (
        "channel 4 for green", "add the dead", "dead channel helps",
    )):
        errors.append("continued with channel advice despite the index conflict")
    if result["calls"]:
        errors.append("attempted channel edits before resolving the conflict")
    return errors


def _check_apoc_dead_channel(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "13t -> channel 1", "blue -> channel 0", "green -> channel 3",
        "do not automatically add", "negative/background class",
        "dead cell is still",
    ):
        if phrase not in text:
            errors.append(f"missing APOC channel guidance: {phrase}")
    if any(phrase in text for phrase in (
        "dead channel helps", "separate class", "exclude dead",
    )):
        errors.append("recommended the dead signal as an APOC negative class")
    if result["calls"]:
        errors.append("attempted channel edits for an advice-only question")
    return errors


def _check_apoc_feature_grid(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "feature scales in pixels", "gaussian blur", "difference of gaussians",
        "laplacian of gaussian", "sobel-of-gaussian",
        "not a structure tensor", "current live apoc controls",
    ):
        if phrase not in text:
            errors.append(f"missing APOC feature-grid detail: {phrase}")
    positive_switch = any(phrase in text for phrase in (
        "need to switch to pixel classifier",
        "should switch to pixel classifier",
        "recommend switching to pixel classifier",
    )) and "no need to switch to pixel classifier" not in text
    if (
        "apoc does not expose" in text
        or "handled internally" in text
        or positive_switch
    ):
        errors.append("denied the APOC feature grid or suggested the wrong method")
    if result["calls"]:
        errors.append("attempted a feature edit without requested values")
    return errors


def _check_apoc_tune_features_explanation(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "segmentation > apoc > tune features",
        "gaussian blur", "difference of gaussians",
        "laplacian of gaussian", "sobel-of-gaussian",
        "small structures", "medium", "large",
        "original intensity", "changes here require retraining",
    ):
        if phrase not in text:
            errors.append(f"missing Tune Features explanation: {phrase}")
    if any(phrase in text for phrase in (
        "feature groups", "contact distance", "dead-pixel threshold",
        "raise the seed threshold",
    )):
        errors.append("drifted from APOC Tune Features into another control group")
    if result["calls"]:
        errors.append("attempted edits for an explanation-only request")
    return errors


def _check_apoc_mdo_feature_recommendation(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "mdo", "organoid", "large structures", "1, 2, 5, 10, and 25 pixels",
        "probability-map preview",
    ):
        if phrase not in text:
            errors.append(f"missing MDO Tune Features recommendation: {phrase}")
    if any(phrase in text for phrase in (
        "raise the seed threshold", "set minimum size", "feature extraction settings",
    )):
        errors.append("answered with post-processing or Feature Extraction")
    if result["calls"]:
        errors.append("attempted an edit for a recommendation-only request")
    return errors


def _check_apoc_fill_organoid_features(result: dict) -> list[str]:
    changed = _changed_values(result)
    expected = {
        "segmentation.apoc.27t.feature_preset": "Large structures",
        "segmentation.apoc.mdo.feature_preset": "Large structures",
        "segmentation.apoc.27t.show_feature_tuning": True,
        "segmentation.apoc.mdo.show_feature_tuning": True,
    }
    errors = []
    for control_id, value in expected.items():
        if changed.get(control_id) != value:
            errors.append(
                f"{control_id} was {changed.get(control_id)!r}, expected {value!r}"
            )
    text = result["text"].lower()
    if "segmentation > apoc > tune features" not in text:
        errors.append("did not identify the requested Segmentation panel")
    if any(str(control_id).startswith("features.") for control_id in changed):
        errors.append("changed Feature Extraction instead of APOC classifier features")
    if any(term in str(control_id) for control_id in changed for term in (
        "minimum_size", "seed_threshold", "mask_threshold",
    )):
        errors.append("changed post-processing instead of the Tune Features panel")
    return errors


def _check_tracking_which_method(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "consecutive frames" not in text or "tcell" not in text:
        errors.append("did not ask for movement for the active tracking cell type")
    if any(term in text for term in (
        "segmentation", "cellpose", "apoc", "pixel classifier",
    )):
        errors.append("answered a Tracking question with segmentation methods")
    if result["calls"]:
        errors.append("changed tracking method before motion was known")
    return errors


def _check_btrack_step2_enable(result: dict) -> list[str]:
    changed = _changed_values(result)
    errors = []
    expected = "tracking.tcell.btrack.use_global_optimization"
    if changed.get(expected) is not True:
        errors.append("did not propose enabling btrack global optimization")
    text = result["text"].lower()
    for phrase in (
        "global hypothesis optimizer", "not the organoid propagation",
        "step 1", "false positive", "initialization", "termination", "linking",
        "largest missing-frame gap",
    ):
        if phrase not in text:
            errors.append(f"missing btrack Step 2 guidance: {phrase}")
    if any(term in text for term in ("segmentation", "cellpose", "apoc")):
        errors.append("left the Tracking module while explaining Step 2")
    return errors


def _check_btrack_step2_tune(result: dict) -> list[str]:
    changed = _changed_values(result)
    expected = {
        "tracking.tcell.btrack.distance_threshold": 40,
        "tracking.tcell.btrack.time_threshold": 4,
    }
    errors = []
    for control_id, value in expected.items():
        if changed.get(control_id) != value:
            errors.append(
                f"{control_id} was {changed.get(control_id)!r}, expected {value!r}"
            )
    text = result["text"].lower()
    if "p_fp, p_init, p_term, p_link" not in text:
        errors.append("did not preserve the normal Step 2 hypotheses")
    if "branching, death, or merging" not in text:
        errors.append("did not scope optional biological hypotheses")
    return errors


def _check_segmentation_minimum_size(result: dict) -> list[str]:
    changed = _changed_values(result)
    errors = []
    if changed.get("segmentation.apoc.tcell.minimum_size") != 131:
        errors.append("did not set the half-volume starting point to 131 voxels")
    text = result["text"].lower()
    for phrase in (
        "10 µm diameter", "524 µm³", "50%", "131 voxels",
        "post-processing exclusion cutoff", "preview",
    ):
        if phrase not in text:
            errors.append(f"missing Minimum size calculation detail: {phrase}")
    if any(phrase in text for phrase in (
        "standard default", "typical minimum", "feature extraction",
    )):
        errors.append("used an ungrounded default or wrong module")
    return errors


def _check_mask_edt_direction(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "higher edt generally means more splitting",
        "raise for more splitting, lower for less",
        "falls back to the original unsplit component",
    ):
        if phrase not in text:
            errors.append(f"missing Mask + EDT direction detail: {phrase}")
    if any(phrase in text for phrase in (
        "higher edt means less splitting", "lower it for more splitting",
        "50 is more conservative", "harder to split",
    )):
        errors.append("repeated the reversed Mask + EDT direction")
    if result["calls"]:
        errors.append("attempted an EDT change during an explanation request")
    return errors


def _check_contact_and_dead_threshold(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "0 µm means strict mask touching", "one xy pixel",
        "not a voxel diagonal", "one-pixel xy gap",
        "green is below", "red is above", "universal numeric range",
    ):
        if phrase not in text:
            errors.append(f"missing feature-threshold guidance: {phrase}")
    if "1.01 µm is strict" in text or "one voxel diagonal" in text:
        errors.append("misstated the 1.01 µm contact scale")
    if result["calls"]:
        errors.append("attempted a threshold edit without a requested value")
    return errors


def _check_loaded_metadata_not_unsaved(result: dict) -> list[str]:
    text = result["text"].lower()
    if any(phrase in text for phrase in (
        "save metadata", "not yet saved", "metadata draft", "save the metadata",
    )):
        return ["incorrectly described loaded metadata as an unsaved draft"]
    return []


def _check_external_zarr_reload(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "load metadata" not in text:
        errors.append("did not tell the user to reload externally updated metadata")
    if "dimension order" in text and "also" not in text:
        errors.append("treated dimension order as the sole cause after external conversion")
    return errors


def _check_missing_log_error(result: dict) -> list[str]:
    text = result["text"].lower()
    asks_for_log = (
        ("copy" in text or "paste" in text)
        and "log" in text
        and "error" in text
    )
    errors = []
    if not asks_for_log:
        errors.append("did not request the exact log error before diagnosing")
    if any(phrase in text for phrase in (
        "the reason is", "this is because", "dimension order is blank",
    )):
        errors.append("claimed an unsupported cause without a log error")
    if any(phrase in text for phrase in (
        "file path", "dimension order", "gpu memory", "insufficient memory",
    )):
        errors.append("speculated about possible causes before reading the error")
    return errors


def _check_tracking_guide(result: dict) -> list[str]:
    text = result["text"].lower()
    asks_motion = any(word in text for word in ("move", "motion", "displacement", "overlap"))
    asks_question = "?" in result["text"] or "how much" in text
    errors = []
    if not (asks_motion and asks_question):
        errors.append("did not ask about frame-to-frame motion before choosing methods")
    if result["calls"]:
        errors.append("changed tracking settings before learning how structures move")
    if any(name in text for name in (
        "btrack", "propagation", "trackpy", "laptrack",
    )):
        errors.append("discussed named tracking methods before learning how the structure moves")
    if any(token in text for token in ("um/min", "µm/min", "micrometres per minute")):
        errors.append("invented a numeric example speed before the user quantified movement")
    return errors


def _check_stationary_tracking(result: dict) -> list[str]:
    text = result["text"].lower()
    method_calls = _tool_calls(result, "set_ui_value")
    proposed = any(
        call.get("control_id") == "tracking.collagen.method"
        and str(call.get("value", "")).lower().startswith("propagation")
        for call in method_calls
    )
    if "propagation" not in text and not proposed:
        return ["did not recommend Propagation for stationary collagen"]
    return []


def _check_tracking_radius(result: dict) -> list[str]:
    calls = _tool_calls(result, "set_ui_value")
    match = next((
        call for call in calls
        if call.get("control_id") == "tracking.Tcell_cmtrm.btrack.maximum_search_radius"
    ), None)
    if match is None:
        return ["did not propose a maximum search-radius action"]
    try:
        value = float(match.get("value"))
    except (TypeError, ValueError):
        return [f"search radius was not numeric: {match.get('value')!r}"]
    if not 16.5 <= value <= 19:
        return [f"search radius was {value:g} um; expected 16.5-19 um from 15 um/frame plus margin"]
    return []


def _check_filtering(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if any(phrase in text for phrase in (
        "are contradictory", "is contradictory", "settings conflict",
        "are conflicting", "is incompatible",
    )):
        errors.append("described equal minimum and maximum track lengths as a conflict")
    if not any(word in text for word in (
        "same length", "equal length", "common length", "comparable", "uniform length",
    )):
        errors.append("did not explain the fixed-length comparison workflow")
    if "reasonable threshold" in text or "reasonable minimum" in text:
        errors.append("endorsed the minimum before reading the track-length distribution")
    if "summarize_track_counts" in text:
        errors.append("exposed the internal track-count tool name")
    return errors


def _changed_values(result: dict) -> dict:
    return {
        call.get("control_id"): call.get("value")
        for call in _tool_calls(result, "set_ui_value")
    }


def _check_reporter_propagation(result: dict) -> list[str]:
    changed = _changed_values(result)
    value = changed.get("tracking.calcium_reporter.method")
    if str(value).lower() != "reporter propagation":
        return [f"tracking method was {value!r}, expected Reporter Propagation"]
    return []


def _check_active_killing(result: dict) -> list[str]:
    changed = _changed_values(result)
    expected = {
        "features.active_killing.target_types": ["organoid1"],
        "features.active_killing.observation_window": 5,
        "features.active_killing.death_signal": "Dead-mask pixel count",
        "features.active_killing.use_absolute_threshold": True,
        "features.active_killing.absolute_threshold": 30,
    }
    errors = []
    for control_id, value in expected.items():
        if changed.get(control_id) != value:
            errors.append(
                f"{control_id} was {changed.get(control_id)!r}, expected {value!r}"
            )
    text = result["text"].lower()
    if any(name in text for name in (
        "nr_dead_mask_pixels", "percentage_dead_mask", "mean_dead_dye",
    )):
        errors.append("exposed an internal Active Killing column name")
    if any(phrase in text for phrase in (
        "i've set", "i have set", "changes are applied", "changes were applied",
    )):
        errors.append("claimed proposed Active Killing changes were already applied")
    return errors


def _check_feature_group_dead_dye(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "required", "intensity", "mean dead-dye intensity",
        "will not suggest removing",
    ):
        if phrase not in text:
            errors.append(f"missing mandatory feature guidance: {phrase}")
    if "drop intensity" in text or "remove intensity" in text:
        errors.append("still suggested removing required T-cell intensity")
    if result["calls"]:
        errors.append("changed feature groups before the optional-group choice")
    return errors


def _check_active_killing_complete_acceptance(result: dict) -> list[str]:
    changed = _changed_values(result)
    expected = {
        "features.active_killing.target_types": ["27t", "mdo"],
        "features.active_killing.observation_window": 5,
        "features.active_killing.death_signal": "Dead-mask pixel count",
        "features.active_killing.use_absolute_threshold": True,
        "features.active_killing.absolute_threshold": 30,
        "features.active_killing.minimum_contact_duration": 1,
    }
    errors = []
    for control_id, value in expected.items():
        if changed.get(control_id) != value:
            errors.append(
                f"{control_id} was {changed.get(control_id)!r}, expected {value!r}"
            )
    text = result["text"].lower()
    for phrase in (
        "complete agreed active killing setup",
        "independently",
        "combined analysis",
        "not ready until every action card",
    ):
        if phrase not in text:
            errors.append(f"missing setup-completeness guidance: {phrase}")
    return errors


def _check_hmm_movement_options(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "speed", "displacement", "cumulative displacement",
        "displacement from origin", "directional persistence",
        "median turning angle", "net displacement", "straightness",
        "mean square displacement", "use all available movement features",
    ):
        if phrase not in text:
            errors.append(f"missing movement option: {phrase}")
    if "mean dead dye" in text:
        errors.append("included a non-movement intensity feature")
    if result["calls"]:
        errors.append("changed HMM inputs before the researcher chose features")
    return errors


def _check_hmm_apply_all_movement(result: dict) -> list[str]:
    changed = _changed_values(result)
    prefix = "analysis.state_classification.tcell."
    errors = []
    expected_timepoint = [
        "speed", "displacement", "cumulative_displacement",
        "displacement_from_origin", "directional_persistence",
        "median_turning_angle",
    ]
    expected_window = [
        "net_displacement", "straightness", "mean_square_displacement",
    ]
    if changed.get(prefix + "timepoint_features") != expected_timepoint:
        errors.append("did not propose all offered timepoint movement features")
    if changed.get(prefix + "window_features") != expected_window:
        errors.append("did not propose all offered window movement features")
    if "mean_dead_dye" in str(changed):
        errors.append("included a non-movement intensity feature")
    if "complete movement-only selection" not in result["text"].lower():
        errors.append("did not describe the two-list selection as complete")
    return errors


def _check_active_killing_zero_threshold_readiness(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "not ready yet" not in text or "greater than 0" not in text:
        errors.append("did not report the zero absolute threshold as incomplete")
    if result["calls"]:
        errors.append("attempted another partial edit instead of reporting readiness")
    return errors


def _check_hmm_single_frame(result: dict) -> list[str]:
    changed = _changed_values(result)
    errors = []
    for suffix in ("window_size", "smooth_window"):
        control_id = f"analysis.state_classification.tcell.{suffix}"
        if changed.get(control_id) != 1:
            errors.append(f"{control_id} was not set to 1")
    if "analysis.state_classification.tcell.start_offset" in changed:
        errors.append("needlessly rewrote the already-correct Start offset")
    return errors


def _check_trajectory_linkage(result: dict) -> list[str]:
    changed = _changed_values(result)
    value = changed.get("analysis.state_trajectory.tcell.linkage")
    if str(value).lower() != "average":
        return [f"trajectory linkage was {value!r}, expected average"]
    return []


def _check_functional_experiment_context(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if not (
        "within" in text and "well" in text
        and ("paired" in text or "same well" in text)
        and "ko" in text and ("rescue" in text or "oe" in text)
    ):
        errors.append("did not prioritize the paired within-well KO/rescue comparison")
    if not any(phrase in text for phrase in (
        "not computed", "not available", "wasn't computed", "was not run",
        "would need", "need to enable",
    )):
        errors.append("did not state that invasiveness is unavailable")
    if result["calls"]:
        errors.append("attempted a UI action for an interpretation-only question")
    if any(name in text for name in (
        "{target}_invasiveness_perc", "is_active_killing",
    )):
        errors.append("exposed an internal analysis column name")
    return errors


def _check_safety_profiling_context(result: dict) -> list[str]:
    text = result["text"].lower()
    plain = text.replace("*", "").replace("`", "")
    errors = []
    if not all(token in text for token in ("27t", "mdo", "1.5")):
        errors.append("lost the tumor/healthy identities or 1.5x definition")
    if not (
        ("5 frame" in text or "five frame" in text)
        and ("10 minute" in text or "ten minute" in text)
    ):
        errors.append("lost the 5-frame / 10-minute observation window")
    if not any(word in text for word in ("exploratory", "descriptive", "small")):
        errors.append("omitted the small-sample interpretation caveat")
    if any(name in text for name in (
        "percentage_dead_mask", "killing_efficiency", "is_active_killing",
    )):
        errors.append("exposed an internal analysis column name")
    if re.search(r"(?:one\s+)?target type\s*(?:is|=|:|\()\s*teg\b", plain):
        errors.append("misidentified the TEG effector as a target type")
    if any(phrase in plain for phrase in (
        "you've configured", "you have configured", "already configured",
        "configured this using", "your setting",
    )):
        errors.append("presented a reference definition as an applied configuration")
    if any(phrase in plain for phrase in (
        "1-3 frame", "1–3 frame", "one to three frame",
    )):
        errors.append("replaced the stated one-frame minimum with a generic range")
    if result["calls"]:
        errors.append("attempted a UI action for an interpretation-only question")
    return errors


def _check_historical_btrack_examples(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in ("ivm", "cd4", "12", "10", "100"):
        if phrase not in text:
            errors.append(f"missing historical btrack context: {phrase}")
    if not any(phrase in text for phrase in (
        "not a default", "not defaults", "do not copy", "shouldn't copy",
        "cannot copy", "historical example",
    )):
        errors.append("presented historical btrack values as reusable defaults")
    if not (
        any(term in text for term in ("per-frame", "per frame", "displacement"))
        and any(term in text for term in ("gap", "preview", "measure"))
    ):
        errors.append("did not explain how to calibrate the current experiment")
    if result["calls"]:
        errors.append("attempted to apply historical btrack values")
    return errors


def _check_historical_calcium_example(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    for phrase in (
        "reporter propagation", "0.33", "2", "5", "32", "100", "10%",
    ):
        if phrase not in text:
            errors.append(f"missing calcium reference detail: {phrase}")
    if not all(term in text for term in ("cellpose", "import")):
        errors.append("omitted the external Cellpose-SAM import workflow")
    if not any(term in text for term in (
        "historical", "example", "not a default", "not universal",
    )):
        errors.append("did not label calcium settings as historical examples")
    if result["calls"]:
        errors.append("attempted to apply the calcium reference settings")
    return errors


def _check_experiment_reference_metadata_conflict(result: dict) -> list[str]:
    text = result["text"].lower()
    errors = []
    if "15" not in text or "10" not in text:
        errors.append("did not report both sides of the cadence discrepancy")
    if not any(term in text for term in (
        "metadata", "currently loaded", "live",
    )):
        errors.append("did not identify live metadata as the current source")
    if not any(term in text for term in (
        "conflict", "disagree", "discrep", "unless", "confirm",
    )):
        errors.append("silently reconciled conflicting reference sources")
    if any(term in text for term in ("time_interval", "time_unit")):
        errors.append("exposed internal metadata field names")
    if result["calls"]:
        errors.append("attempted to edit metadata during a source comparison")
    return errors


SCENARIOS = [
    _organoid_line_grouping_case,
    _metadata_identifier_confirmation_case,
    _metadata_time_conversion_case,
    _metadata_completion_case,
    _metadata_save_case,
    _choose_analysis_case,
    _open_death_dynamics_case,
    _metadata_setup_case,
    _pixel_fill_case,
    _segmentation_case,
    _method_requires_signal_context_case,
    _cellpose_requires_bleed_confirmation_case,
    _confirmed_bleed_through_case,
    _apoc_channel_map_case,
    _apoc_channels_wait_for_training_data_case,
    _apoc_feature_preset_case,
    _swapped_channel_metadata_case,
    _apoc_invalid_channel_index_case,
    _apoc_dead_channel_case,
    _apoc_feature_grid_case,
    _apoc_tune_features_explanation_case,
    _apoc_mdo_feature_recommendation_case,
    _apoc_fill_organoid_features_case,
    _tracking_which_method_case,
    _btrack_step2_enable_case,
    _btrack_step2_tune_case,
    _segmentation_minimum_size_case,
    _mask_edt_direction_case,
    _contact_and_dead_threshold_case,
    _loaded_metadata_not_unsaved_case,
    _external_zarr_reload_case,
    _missing_log_error_case,
    _tracking_guide_case,
    _stationary_tracking_case,
    _tracking_radius_case,
    _filtering_case,
    _reporter_propagation_case,
    _active_killing_case,
    _feature_group_dead_dye_case,
    _active_killing_complete_acceptance_case,
    _hmm_movement_options_case,
    _hmm_apply_all_movement_case,
    _active_killing_zero_threshold_readiness_case,
    _hmm_single_frame_case,
    _trajectory_linkage_case,
    _functional_experiment_context_case,
    _safety_profiling_context_case,
    _historical_btrack_examples_case,
    _historical_calcium_example_case,
    _experiment_reference_metadata_conflict_case,
]


def _default_endpoint() -> str:
    endpoint = os.environ.get("BEHAV3D_ASSISTANT_ENDPOINT", "")
    if endpoint:
        return endpoint.rstrip("/")
    config = Path(__file__).parents[1] / "napari" / "assistant_config.json"
    if config.exists():
        return str(json.loads(config.read_text(encoding="utf-8")).get("endpoint", "")).rstrip("/")
    return ""


def _request(endpoint: str, case: dict, timeout: int) -> dict:
    response = requests.post(
        f"{endpoint}/chat",
        json={"messages": case["messages"], "context": case["context"], "tools": TOOL_SCHEMA},
        headers={"Content-Type": "application/json"},
        stream=True,
        timeout=(10, timeout),
    )
    response.raise_for_status()
    text_parts: list[str] = []
    calls: list[dict] = []
    for raw in response.iter_lines(decode_unicode=True):
        if not raw:
            continue
        line = raw[6:] if raw.startswith("data: ") else raw
        event = json.loads(line)
        if event.get("type") == "token":
            text_parts.append(event.get("text", ""))
        elif event.get("type") == "tool_calls":
            calls.extend(event.get("calls") or [])
        elif event.get("type") == "error":
            raise RuntimeError(event.get("message") or "Unknown assistant API error")
    return {"text": "".join(text_parts).strip(), "calls": calls}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", default=_default_endpoint())
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument(
        "--scenario", action="append",
        choices=[factory()["name"] for factory in SCENARIOS],
        help="Run only this scenario; repeat to select several.",
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    if not args.endpoint:
        parser.error("Provide --endpoint or configure BEHAV3D_ASSISTANT_ENDPOINT")

    selected = [factory() for factory in SCENARIOS]
    if args.scenario:
        selected = [case for case in selected if case["name"] in args.scenario]

    report = []
    for case in selected:
        print(f"\n=== {case['name']} ===", flush=True)
        try:
            result = _request(args.endpoint, case, args.timeout)
            errors = case["check"](result)
        except Exception as exc:
            result, errors = {"text": "", "calls": []}, [str(exc)]
        report.append({"scenario": case["name"], **result, "errors": errors})
        print(result["text"] or "(no visible text)")
        print("Tool calls:", json.dumps(result["calls"], indent=2, ensure_ascii=False))
        print("PASS" if not errors else "FAIL: " + "; ".join(errors))

    if args.json_out:
        args.json_out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    failures = sum(bool(item["errors"]) for item in report)
    print(f"\n{len(report) - failures}/{len(report)} scenarios passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
