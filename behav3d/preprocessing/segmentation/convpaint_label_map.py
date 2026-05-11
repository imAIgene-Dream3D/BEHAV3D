"""
Label-map helpers for the unified multi-class ConvPaint classifier.

A unified ConvPaint model classifies every voxel into one of N+1 classes:

* class index 1 = background
* class indices 2..N+1 = one per cell type

The mapping between class indices and cell-type names is the single source of
truth shared between the training widget (``convpaint_train.py``) and the
batch inference queue (``convpaint_segment.py``). It is persisted as a JSON
sidecar next to the trained ``ConvPaintModel_Unified.pkl`` file so that
inference can decode predictions without re-deriving the order from the
current metadata (which may have changed since training).

Schema (``schema_version: 1``)::

    {
        "schema_version": 1,
        "background_label": 1,
        "cell_type_order": ["organoid_a", "tcells", "tcells_1_multicolor"],
        "label_to_celltype": {"2": "organoid_a", "3": "tcells", ...},
        "celltype_to_label": {"organoid_a": 2, "tcells": 3, ...}
    }

JSON object keys must be strings, so ``label_to_celltype`` keys are stringified
integers. Use :func:`celltype_for_index` / :func:`index_for_celltype` to read
the map without worrying about the key type.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Optional


SCHEMA_VERSION = 1
BACKGROUND_LABEL = 1
FIRST_CELL_TYPE_LABEL = BACKGROUND_LABEL + 1

UNIFIED_MODEL_FILENAME = "ConvPaintModel_Unified.pkl"
UNIFIED_LABEL_MAP_FILENAME = "ConvPaintModel_Unified_LabelMap.json"
UNIFIED_INPUT_CHANNELS_FILENAME = "ConvPaintModel_Unified_InputChannels.json"
DEATH_INPUT_CHANNELS_FILENAME = "ConvPaintModel_Dead_InputChannels.json"
UNIFIED_USER_LABELS_FILENAME = "PixelClassifier_UserUnifiedLabels.zarr"
UNIFIED_PREDICTED_LABELS_FILENAME = "ConvPaintModel_Unified_PredictedLabels.zarr"


# ---------------------------------------------------------------------------
# Canonical paths
# ---------------------------------------------------------------------------

def unified_model_path(pixel_class_outdir) -> Path:
    """Path to the unified ConvPaint model file."""
    return Path(pixel_class_outdir) / UNIFIED_MODEL_FILENAME


def unified_label_map_path(pixel_class_outdir) -> Path:
    """Path to the JSON label-map sidecar that accompanies the unified model."""
    return Path(pixel_class_outdir) / UNIFIED_LABEL_MAP_FILENAME


def unified_input_channels_path(pixel_class_outdir) -> Path:
    """Path to the unified model input-channel sidecar."""
    return Path(pixel_class_outdir) / UNIFIED_INPUT_CHANNELS_FILENAME


def death_input_channels_path(pixel_class_outdir) -> Path:
    """Path to the death model input-channel sidecar."""
    return Path(pixel_class_outdir) / DEATH_INPUT_CHANNELS_FILENAME


def unified_user_labels_path(pixel_class_outdir) -> Path:
    """Path to the unified ``User Provided Labels`` zarr cache."""
    return Path(pixel_class_outdir) / UNIFIED_USER_LABELS_FILENAME


def unified_predicted_labels_path(pixel_class_outdir) -> Path:
    """Path to the cached multi-class predicted-labels volume (training stack)."""
    return Path(pixel_class_outdir) / UNIFIED_PREDICTED_LABELS_FILENAME


# ---------------------------------------------------------------------------
# Build / save / load
# ---------------------------------------------------------------------------

def build_label_map(cell_types: Iterable[str]) -> dict:
    """Build a deterministic label map for the given cell-type order.

    The order of ``cell_types`` is preserved (e.g. organoid → immune → other,
    matching the panel's ``all_cell_types`` ordering). Duplicates are dropped
    while keeping the first occurrence. Empty strings and ``None`` are ignored.

    Parameters
    ----------
    cell_types : iterable of str
        Cell-type names to assign indices to. Must NOT include ``"background"``
        (background is always class 1 and is implicit).

    Returns
    -------
    dict
        Label-map dictionary matching the documented schema.

    Raises
    ------
    ValueError
        If ``cell_types`` is empty after filtering, or contains ``"background"``.
    """
    seen: set[str] = set()
    ordered: list[str] = []
    for ct in cell_types:
        if not ct:
            continue
        ct_str = str(ct).strip()
        if not ct_str:
            continue
        if ct_str.lower() == "background":
            raise ValueError(
                "'background' is reserved as class 1 and cannot appear in "
                "the unified cell-type list."
            )
        if ct_str in seen:
            continue
        seen.add(ct_str)
        ordered.append(ct_str)

    if not ordered:
        raise ValueError(
            "build_label_map() requires at least one cell type; got empty list."
        )

    label_to_celltype: dict[str, str] = {}
    celltype_to_label: dict[str, int] = {}
    for offset, ct in enumerate(ordered):
        idx = FIRST_CELL_TYPE_LABEL + offset
        label_to_celltype[str(idx)] = ct
        celltype_to_label[ct] = idx

    return {
        "schema_version": SCHEMA_VERSION,
        "background_label": BACKGROUND_LABEL,
        "cell_type_order": ordered,
        "label_to_celltype": label_to_celltype,
        "celltype_to_label": celltype_to_label,
    }


def save_label_map(path, label_map: dict) -> None:
    """Persist ``label_map`` as JSON, overwriting any existing file."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as fh:
        json.dump(label_map, fh, indent=2, sort_keys=False)


def load_label_map(path) -> dict:
    """Load a label map from JSON and validate the schema version.

    Raises
    ------
    FileNotFoundError
        If the sidecar does not exist.
    ValueError
        If the schema version is unknown or the file is malformed.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Label-map sidecar not found: {p}")
    with p.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError(f"Label-map at {p} is not a JSON object.")

    schema = int(data.get("schema_version", 0))
    if schema != SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported label-map schema_version={schema} (expected {SCHEMA_VERSION}). "
            f"File: {p}"
        )

    required = {"background_label", "cell_type_order",
                "label_to_celltype", "celltype_to_label"}
    missing = required - set(data.keys())
    if missing:
        raise ValueError(
            f"Label-map at {p} is missing required keys: {sorted(missing)}"
        )

    # Normalize: ensure label_to_celltype keys are strings (json already enforces this)
    # and celltype_to_label values are ints (json may store as int already, but be safe).
    data["celltype_to_label"] = {
        str(k): int(v) for k, v in data["celltype_to_label"].items()
    }
    data["label_to_celltype"] = {
        str(k): str(v) for k, v in data["label_to_celltype"].items()
    }
    return data


def save_input_channels(path, channel_indices, *, model_name=None) -> None:
    """Persist the raw channel indices used to train a ConvPaint model.

    This sidecar lets batch inference feed the same channel subset to a newly
    trained model while old models without a sidecar keep using all channels.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    channels = []
    for idx in channel_indices or []:
        try:
            idx = int(idx)
        except (TypeError, ValueError):
            continue
        if idx >= 0 and idx not in channels:
            channels.append(idx)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "model_name": str(model_name or ""),
        "input_channel_indices": channels,
    }
    with p.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)


def load_input_channels(path) -> Optional[list[int]]:
    """Load a ConvPaint input-channel sidecar.

    Returns ``None`` when the sidecar is missing or malformed so callers can
    safely fall back to the legacy all-channel behavior.
    """
    p = Path(path)
    if not p.exists():
        return None
    try:
        with p.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception:
        return None
    if not isinstance(data, dict):
        return None
    try:
        schema = int(data.get("schema_version", 0))
    except (TypeError, ValueError):
        return None
    if schema != SCHEMA_VERSION:
        return None
    raw_channels = data.get("input_channel_indices")
    if not isinstance(raw_channels, list):
        return None
    channels = []
    for idx in raw_channels:
        try:
            idx = int(idx)
        except (TypeError, ValueError):
            continue
        if idx >= 0 and idx not in channels:
            channels.append(idx)
    return channels or None


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def cell_type_order(label_map: dict) -> list[str]:
    """Return the ordered list of cell-type names (excluding background)."""
    return list(label_map.get("cell_type_order", []))


def n_classes(label_map: dict) -> int:
    """Total number of classes including background."""
    return 1 + len(cell_type_order(label_map))


def index_for_celltype(label_map: dict, cell_type: str) -> Optional[int]:
    """Return the class index for ``cell_type`` (or ``None`` if unknown)."""
    if not cell_type:
        return None
    raw = label_map.get("celltype_to_label", {}).get(str(cell_type))
    if raw is None:
        return None
    return int(raw)


def celltype_for_index(label_map: dict, idx: int) -> Optional[str]:
    """Return the cell-type name for class index ``idx`` (or ``None``).

    Returns ``None`` for the background label so callers always receive a real
    cell-type name (or nothing).
    """
    if idx is None:
        return None
    if int(idx) == int(label_map.get("background_label", BACKGROUND_LABEL)):
        return None
    return label_map.get("label_to_celltype", {}).get(str(int(idx)))


def label_map_matches(label_map: dict, cell_types: Iterable[str]) -> bool:
    """Return ``True`` iff ``cell_types`` (in order) matches the saved map.

    Useful at training time to decide whether the cached unified labels can be
    reused or must be rebuilt because the cell-type set has changed.
    """
    expected = list(cell_type_order(label_map))
    given = [str(ct).strip() for ct in cell_types if ct]
    return expected == given
