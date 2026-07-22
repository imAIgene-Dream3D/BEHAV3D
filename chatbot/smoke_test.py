"""Replay PI feedback scenarios against a deployed BEHAV3D Assistant API."""
from __future__ import annotations

import argparse
import json
import os
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
) -> dict:
    return {
        "id": control_id,
        "label": label,
        "value": value,
        "choices": choices,
        "visible": True,
        "enabled": True,
        "method": method,
        "strategy": strategy,
        "cell_type": cell_type,
        "unit": unit,
    }


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


def _tracking_guide_case() -> dict:
    choices = ["LAP (laptrack)", "TrackPy", "Propagation", "btrack (Bayesian)"]
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
        choices=["LAP (laptrack)", "TrackPy", "Propagation", "btrack (Bayesian)"],
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


def _tool_calls(result: dict, name: str) -> list[dict]:
    return [call.get("arguments", {}) for call in result["calls"] if call.get("name") == name]


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


def _check_tracking_guide(result: dict) -> list[str]:
    text = result["text"].lower()
    asks_motion = any(word in text for word in ("move", "motion", "displacement", "overlap"))
    asks_question = "?" in result["text"] or "how much" in text
    errors = []
    if not (asks_motion and asks_question):
        errors.append("did not ask about frame-to-frame motion before choosing methods")
    if result["calls"]:
        errors.append("changed tracking settings before learning how structures move")
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
    if any(phrase in text for phrase in (
        "are contradictory", "is contradictory", "settings conflict",
        "are conflicting", "is incompatible",
    )):
        return ["described equal minimum and maximum track lengths as a conflict"]
    if not any(word in text for word in ("same length", "equal length", "common length", "comparable")):
        return ["did not explain the fixed-length comparison workflow"]
    if "reasonable threshold" in text or "reasonable minimum" in text:
        return ["endorsed the minimum before reading the track-length distribution"]
    return []


SCENARIOS = [
    _metadata_setup_case,
    _pixel_fill_case,
    _segmentation_case,
    _tracking_guide_case,
    _stationary_tracking_case,
    _tracking_radius_case,
    _filtering_case,
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
