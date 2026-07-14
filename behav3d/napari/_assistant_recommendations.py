"""Deterministic calculations used by the BEHAV3D assistant."""
from __future__ import annotations

import math
from statistics import median
from typing import Any


EDT_DIAMETER_FRACTIONS = (0.20, 0.25, 0.30)


def _distance_um(value: Any, unit: Any) -> float:
    from behav3d.core.utils import convert_distance

    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError("Pixel size must be greater than zero.")
    normalized_unit = str(unit or "um").strip() or "um"
    if normalized_unit in {"um", "μm", "µm"}:
        normalized_unit = "um"
    return float(convert_distance(number, normalized_unit))


def _half_pixel(value: float) -> float:
    return max(0.5, round(float(value) * 2.0) / 2.0)


def calculate_edt_recommendations(
    records: list[dict],
    *,
    cell_diameter_um: float = 10.0,
    organoid_cells_across: float | None = None,
) -> dict:
    """Convert metadata resolution into transparent EDT starting values.

    EDT in the current segmentation implementation is measured in voxel/pixel
    units. Candidate thresholds are 20%, 25%, and 30% of the expected XY object
    diameter. These are starting points for previewing, not fitted parameters.

    ``organoid_cells_across`` means the approximate number of 10-um cell widths
    spanning the organoid diameter. This direct geometric input avoids pretending
    that a total cell count uniquely determines organoid diameter.
    """
    diameter_um = float(cell_diameter_um)
    if not math.isfinite(diameter_um) or diameter_um <= 0:
        raise ValueError("Cell diameter must be greater than zero.")
    if organoid_cells_across is not None:
        cells_across = float(organoid_cells_across)
        if not math.isfinite(cells_across) or cells_across < 1:
            raise ValueError("Organoid size must be at least one cell width.")
        object_diameter_um = diameter_um * cells_across
    else:
        cells_across = None
        object_diameter_um = diameter_um

    rows = []
    skipped = []
    for index, record in enumerate(records or []):
        sample = str(record.get("sample_name") or f"Sample {index + 1}")
        try:
            xy_um = _distance_um(
                record.get("pixel_distance_xy"), record.get("distance_unit", "um")
            )
        except (AssertionError, TypeError, ValueError) as exc:
            skipped.append({"sample_name": sample, "reason": str(exc)})
            continue
        diameter_px = object_diameter_um / xy_um
        candidates = [_half_pixel(diameter_px * fraction)
                      for fraction in EDT_DIAMETER_FRACTIONS]
        candidates = list(dict.fromkeys(candidates))
        rows.append({
            "sample_name": sample,
            "pixel_size_xy_um": xy_um,
            "object_diameter_um": object_diameter_um,
            "object_diameter_px": diameter_px,
            "object_radius_px": diameter_px / 2.0,
            "edt_candidates_px": candidates,
            "edt_start_px": _half_pixel(diameter_px * 0.25),
        })

    if not rows:
        raise ValueError("No sample has a valid XY pixel size in the metadata.")

    return {
        "cell_diameter_um": diameter_um,
        "organoid_cells_across": cells_across,
        "object_diameter_um": object_diameter_um,
        "rows": rows,
        "skipped": skipped,
        "global_start_px": _half_pixel(median(row["edt_start_px"] for row in rows)),
        "fractions": list(EDT_DIAMETER_FRACTIONS),
    }


def format_edt_recommendations(result: dict, cell_type: str) -> str:
    rows = result["rows"]
    object_um = result["object_diameter_um"]
    cells_across = result.get("organoid_cells_across")
    assumption = f"a {object_um:g} µm object diameter"
    if cells_across is not None:
        assumption += (
            f" ({cells_across:g} cell widths across, using "
            f"{result['cell_diameter_um']:g} µm per cell)"
        )

    lines = [
        f"**EDT starting points for {cell_type}**",
        "",
        f"Assumption: {assumption}.",
        "",
        "| Sample | XY resolution (µm/pixel) | Diameter (pixels) | EDT values to try |",
        "|---|---:|---:|---:|",
    ]
    for row in rows:
        trials = ", ".join(f"{value:g}" for value in row["edt_candidates_px"])
        lines.append(
            f"| {row['sample_name']} | {row['pixel_size_xy_um']:.4g} | "
            f"{row['object_diameter_px']:.2f} | {trials} |"
        )
    lines.extend([
        "",
        f"A practical global starting value is **{result['global_start_px']:g} pixels** "
        "(the median middle candidate across samples). Lower EDT values split touching "
        "objects more aggressively; higher values keep more objects merged.",
        "",
        "These candidates are 20%, 25%, and 30% of the expected XY diameter. "
        "They are XY-based starting points, so preview segmentation around them before "
        "running the batch.",
    ])
    if any(max(row["edt_candidates_px"]) > 50 for row in rows):
        lines.append(
            "Some candidates exceed the current EDT control maximum of 50 pixels; "
            "review the object-size assumption before applying them."
        )
    if result.get("skipped"):
        samples = ", ".join(item["sample_name"] for item in result["skipped"])
        lines.append(f"Skipped samples with invalid resolution metadata: {samples}.")
    return "\n".join(lines)
