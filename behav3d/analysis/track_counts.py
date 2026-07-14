"""Track-count summaries for evaluating minimum-length filters."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_MINIMUM_LENGTHS = (20, 50, 100, 200)
MAXIMUM_LENGTH_QUERIES = 20


def _track_id_column(frame: pd.DataFrame) -> str:
    for name in ("TrackID", "track_id"):
        if name in frame.columns:
            return name
    raise ValueError("Track feature CSV must contain TrackID or track_id.")


def _normalize_thresholds(values: Iterable[int] | None) -> list[int]:
    source = DEFAULT_MINIMUM_LENGTHS if values is None else values
    thresholds = []
    for value in source:
        number = int(value)
        if number < 1:
            raise ValueError("Minimum track lengths must be positive integers.")
        if number not in thresholds:
            thresholds.append(number)
    if not thresholds:
        raise ValueError("Provide at least one minimum track length.")
    if len(thresholds) > MAXIMUM_LENGTH_QUERIES:
        raise ValueError(
            f"Compare at most {MAXIMUM_LENGTH_QUERIES} minimum track lengths at once."
        )
    return sorted(thresholds)


def calculate_track_count_table(
    frame: pd.DataFrame,
    *,
    minimum_lengths: Iterable[int] | None = None,
    position_t: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Count tracks present at one timepoint before/after length filters.

    Track length is the number of unique ``position_t`` values for each
    ``(sample_name, TrackID)`` pair. A filtered count includes only tracks whose
    full length meets the threshold and that have a row at the requested timepoint.
    """
    thresholds = _normalize_thresholds(minimum_lengths)
    query_t = int(position_t)
    if query_t < 0:
        raise ValueError("Timepoint must be zero or greater.")
    required = {"sample_name", "position_t"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Track feature CSV is missing: {', '.join(sorted(missing))}.")
    track_col = _track_id_column(frame)

    clean = frame[["sample_name", track_col, "position_t"]].copy()
    clean["position_t"] = pd.to_numeric(clean["position_t"], errors="coerce")
    clean = clean.dropna(subset=["sample_name", track_col, "position_t"])
    clean = clean.drop_duplicates(["sample_name", track_col, "position_t"])
    if clean.empty:
        raise ValueError("Track feature CSV has no valid sample, track, and timepoint rows.")

    samples = [str(value) for value in clean["sample_name"].drop_duplicates()]
    clean["sample_name"] = clean["sample_name"].astype(str)
    lengths = (
        clean.groupby(["sample_name", track_col], sort=False)["position_t"]
        .nunique()
        .rename("track_length")
        .reset_index()
    )
    present = clean.loc[clean["position_t"] == query_t, ["sample_name", track_col]]
    present = present.drop_duplicates()

    base_label = f"Cell count at timepoint {query_t} (not filtered)"
    rows = []
    for sample in samples:
        sample_present = present.loc[present["sample_name"] == sample, [track_col]]
        row = {"sample_name": sample, base_label: int(sample_present[track_col].nunique())}
        sample_lengths = lengths.loc[lengths["sample_name"] == sample]
        for threshold in thresholds:
            eligible = sample_lengths.loc[
                sample_lengths["track_length"] >= threshold, track_col
            ]
            row[f">={threshold} timepoints"] = int(
                sample_present.loc[sample_present[track_col].isin(eligible), track_col].nunique()
            )
        rows.append(row)

    table = pd.DataFrame(rows)
    numeric_columns = [column for column in table.columns if column != "sample_name"]
    mean_row = {"sample_name": "Mean"}
    mean_row.update({column: float(table[column].mean()) for column in numeric_columns})
    table = pd.concat([table, pd.DataFrame([mean_row])], ignore_index=True)

    observed = sorted(clean["position_t"].unique().tolist())
    details = {
        "position_t": query_t,
        "minimum_lengths": thresholds,
        "observed_min_t": observed[0],
        "observed_max_t": observed[-1],
        "position_present": query_t in set(observed),
        "track_id_column": track_col,
    }
    return table, details


def track_features_path(output_dir: str | Path, cell_type: str) -> Path:
    return (
        Path(output_dir) / "analysis" / str(cell_type) / "track_features"
        / f"BEHAV3D_{cell_type}_combined_track_features.csv"
    )


def generate_track_count_summary(
    output_dir: str | Path,
    cell_type: str,
    *,
    minimum_lengths: Iterable[int] | None = None,
    position_t: int = 0,
) -> dict:
    source = track_features_path(output_dir, cell_type)
    if not source.exists():
        raise FileNotFoundError(f"Unfiltered track features not found: {source}")
    frame = pd.read_csv(source, low_memory=False)
    table, details = calculate_track_count_table(
        frame, minimum_lengths=minimum_lengths, position_t=position_t
    )
    qc_dir = Path(output_dir) / "analysis" / str(cell_type) / "quality_control"
    if details["position_t"] == 0 and details["minimum_lengths"] == list(
        DEFAULT_MINIMUM_LENGTHS
    ):
        filename = f"BEHAV3D_{cell_type}_track_count_summary.csv"
    else:
        threshold_slug = "-".join(str(value) for value in details["minimum_lengths"])
        filename = (
            f"BEHAV3D_{cell_type}_track_count_query_t{details['position_t']}_"
            f"min_{threshold_slug}.csv"
        )
    destination = qc_dir / filename
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, index=False)
    return {
        "table": table,
        "source_path": source,
        "output_path": destination,
        **details,
    }


def _format_number(value) -> str:
    number = float(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:.2f}".rstrip("0").rstrip(".")


def format_track_count_summary(result: dict, cell_type: str) -> str:
    table = result["table"]
    columns = list(table.columns)
    headers = ["Sample", "Unfiltered"]
    headers.extend(
        column.replace(" timepoints", "").replace(">=", "≥")
        for column in columns[2:]
    )
    lines = [
        f"**Track counts for {cell_type} at timepoint {result['position_t']}**",
        "",
        "| " + " | ".join(headers) + " |",
        "|" + "|".join("---" if i == 0 else "---:" for i in range(len(headers))) + "|",
    ]
    for _, row in table.iterrows():
        values = [str(row[columns[0]])]
        values.extend(_format_number(row[column]) for column in columns[1:])
        lines.append("| " + " | ".join(values) + " |")
    lines.extend([
        "",
        f"Saved under **Filtering results** as `{result['output_path'].name}`.",
    ])
    if not result["position_present"]:
        lines.append(
            f"No rows were recorded at timepoint {result['position_t']}; the observed "
            f"range is {result['observed_min_t']:g} to {result['observed_max_t']:g}."
        )
    return "\n".join(lines)
