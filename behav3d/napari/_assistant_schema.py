"""
BEHAV3D assistant — parameter schema cards.

Flattens the canonical ``_DEFAULT_CONFIG`` (and ``behav3d_calculated_features``)
from :mod:`behav3d.widgets.utils` into a flat list of per-parameter "cards".

A card is a small, self-describing record for one tunable parameter:

    {
        "key": "tracking.immune.trackpy.search_range_px",  # dotted path
        "title": "search_range_px",
        "step": "tracking",          # workflow step it belongs to
        "category": "immune",        # cell-type category, or None
        "type": "int",               # int|float|bool|str|list|none
        "default": 31,
        "choices": None,             # list of allowed values, or None
        "description": "...",        # human-readable guidance for the bot
    }

These cards are the single source of truth shared by:
  * the context serializer (``_assistant_context``) — to describe the current
    step's parameters to the model, and
  * the Modal ingestion job — emitted to ``schema_cards.json`` and indexed for RAG.

The module imports ``_DEFAULT_CONFIG`` *lazily* so it stays cheap to import and
unit-testable (callers may inject a synthetic config instead).
"""
from __future__ import annotations

from typing import Any, Optional

# ---------------------------------------------------------------------------
# Top-level config key  ->  workflow step label
# ---------------------------------------------------------------------------
STEP_MAP = {
    "seed": "general",
    "paths": "data_preparation",
    "dim_order": "data_preparation",
    "signal_unmixing": "preprocessing",
    "pixel_classifier": "segmentation",
    "cellpose": "segmentation",
    "tracking": "tracking",
    "tracking_visualization": "visualization",
    "features": "feature_extraction",
    "track_filtering": "filtering",
    "analysis": "analysis",
    "backprojection": "analysis",
    "active_killing": "feature_extraction",
    "death_dynamics": "feature_extraction",
}

CELL_CATEGORIES = ("immune", "organoid", "other")

# ---------------------------------------------------------------------------
# Known enumerated choices, keyed by the *leaf* parameter name.
# ---------------------------------------------------------------------------
CHOICES: dict[str, list] = {
    "method": ["trackpy", "lap", "btrack", "propagation", "propagation_all_organoids"],
    "update_method": ["EXACT", "APPROXIMATE"],
    "config_preset": ["cell", "particle"],
    "mode": ["mean", "max", "median", "sum"],
    # Must match behav3d.widgets.utils.DIM_ORDER_OPTIONS exactly (these are the
    # actual items in the "Set all" dim-order combo; setCurrentText no-ops otherwise).
    "default_apply_all": [
        "TCZYX", "TZCYX", "CTZYX", "CZTYX", "ZCTYX", "ZTCYX",
        "CZYX", "ZCYX", "TZYX", "ZTYX", "TCYX", "CTYX",
        "ZYX", "CYX", "TYX", "YX",
    ],
    "labels_mode": ["same_for_all", "per_sample"],
    "death_signal_column": ["percentage_dead_mask", "mean_dead_dye"],
}

# ---------------------------------------------------------------------------
# Human-readable descriptions, keyed by leaf parameter name. These give the RAG
# index meaningful per-parameter text. Generic fallback is synthesised for keys
# not listed here.
# ---------------------------------------------------------------------------
_DESCRIPTIONS: dict[str, str] = {
    "seed": "Random seed for reproducible clustering / UMAP / sampling.",
    "metadata_csv": "Path to the BEHAV3D metadata.csv describing every sample.",
    "output_dir": "Directory where all BEHAV3D outputs (zarr, csv, results) are written.",
    "default_apply_all": "Dimension order of the raw images (e.g. TCZYX = time, channel, z, y, x). Used to interpret every input image unless overridden per sample.",
    # pixel classifier
    "examples_per_sample": "How many annotated example timepoints to draw per sample when training the pixel classifier. More examples = better but slower training.",
    "sample_specific_classifier": "Train a separate classifier per sample instead of one shared classifier. Use when samples differ a lot in appearance.",
    "workers": "Number of parallel CPU worker processes. Higher = faster but more RAM/CPU.",
    "use_all_timepoints": "Process every timepoint. Disable to restrict to a [tp_start, tp_end] range.",
    # cellpose
    "number_of_channels": "Number of image channels passed to Cellpose.",
    "labels_mode": "'same_for_all' applies one channel→cell-type label map to every sample; 'per_sample' lets each sample differ.",
    "channel_labels": "Mapping of channel index → cell-type name (e.g. {0: 'organoid', 1: 'tcell'}).",
    # tracking — general
    "track_organoids_together": "Track all organoids as a single merged object rather than individually.",
    "method": "Tracking algorithm: 'trackpy' (linking by proximity), 'lap' (linear assignment with gap-closing), 'btrack' (Bayesian, best for crowded/dividing cells), or 'propagation' (label propagation for organoids).",
    "overwrite": "Re-run and overwrite existing outputs instead of skipping completed samples.",
    # tracking — trackpy
    "search_range_px": "Max distance (pixels) a cell may move between consecutive frames. Set just above the fastest expected displacement; too small drops links, too large makes spurious ones.",
    "memory_frames": "How many frames a cell may disappear (occlusion/missed detection) and still be re-linked to the same track.",
    "adaptive_stop": "trackpy adaptive-search lower bound; linking subdivides the search range down to this value when ambiguous.",
    "adaptive_step": "Factor by which trackpy shrinks the search range on each adaptive subdivision (0–1).",
    # tracking — lap
    "track_cost_px": "LAP frame-to-frame linking cost cutoff in pixels (akin to a max move distance).",
    "gap_close_cost_px": "LAP cost cutoff for bridging gaps where a cell was briefly undetected.",
    "gap_close_max_frames": "Maximum gap length (frames) LAP will bridge when gap-closing.",
    "merging_cost_px": "LAP cost for allowing two tracks to merge (0 disables merging).",
    "splitting_cost_px": "LAP cost for allowing a track to split (0 disables splitting).",
    # tracking — btrack
    "config_preset": "btrack motion-model preset: 'cell' (typical migrating cells) or 'particle' (fast, near-random motion).",
    "config_path": "Optional path to a custom btrack config JSON; overrides the preset.",
    "max_search_radius": "btrack maximum linking radius in pixels.",
    "update_method": "btrack hypothesis update: 'EXACT' (accurate, slower) or 'APPROXIMATE' (faster for large datasets).",
    "step_size": "btrack optimiser batch/step size; larger trades memory for speed.",
    "n_workers": "Parallel workers for this step.",
    "use_optimize": "Run btrack's global optimisation pass (better tracks, slower).",
    "hypotheses": "btrack hypotheses to evaluate (P_FP false-positive, P_init initialisation, P_term termination, P_link linking, P_branch division).",
    "dist_thresh": "btrack distance threshold (pixels) for candidate links.",
    "time_thresh": "btrack time threshold (frames) for candidate links.",
    # features
    "features_choice": "Which feature groups to compute for this cell category: movement, intensity, contact, death, morphology.",
    "contact_threshold": "Distance (pixels) within which two cells count as in contact. 0 = touching only.",
    # filtering
    "exp_duration": "Total experiment duration (hours); used to convert frame counts to time.",
    "exp_duration_enabled": "Whether to apply the experiment-duration based filter.",
    "min_track_length": "Discard tracks shorter than this many frames.",
    "min_track_length_enabled": "Whether the minimum-track-length filter is active.",
    "max_track_length": "Discard tracks longer than this many frames.",
    "max_track_length_enabled": "Whether the maximum-track-length filter is active.",
    # analysis
    "umap_min_dist": "UMAP min_dist: lower packs clusters tighter, higher spreads points out (0–1).",
    "umap_n_neighbors": "UMAP neighbourhood size: low = local structure, high = global structure.",
    "nr_of_clusters": "Number of behavioural clusters to extract.",
    # active killing
    "observation_window": "Frames over which a death-signal increase is assessed for active killing.",
    "death_signal_column": "Feature column used as the death/viability signal.",
    "killing_threshold_multiplier": "Multiplier over baseline death signal that flags an active-killing event.",
    "min_contact_duration": "Minimum frames of immune–target contact required to call active killing.",
    # death dynamics
    "dead_perc_threshold": "Fraction of dead-mask pixels above which a cell is considered dead.",
}


def _infer_type(value: Any) -> str:
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    if isinstance(value, (list, tuple)):
        return "list"
    if value is None:
        return "none"
    return type(value).__name__


def _describe(leaf: str, dotted: str, default: Any, typ: str) -> str:
    if leaf in _DESCRIPTIONS:
        return _DESCRIPTIONS[leaf]
    return f"BEHAV3D parameter '{dotted}' (type {typ}, default {default!r})."


def _step_and_category(parts: list[str]) -> tuple[str, Optional[str]]:
    top = parts[0]
    step = STEP_MAP.get(top, top)
    category = parts[1] if len(parts) > 1 and parts[1] in CELL_CATEGORIES else None
    return step, category


def flatten_config_to_cards(
    config: Optional[dict] = None,
    calculated_features: Optional[dict] = None,
) -> list[dict]:
    """Flatten a nested BEHAV3D config dict into a list of parameter cards.

    If *config* is None, ``_DEFAULT_CONFIG`` is imported lazily from
    :mod:`behav3d.widgets.utils`. Pass a synthetic dict in tests to avoid the
    heavy import.
    """
    if config is None:
        from behav3d.widgets.utils import _DEFAULT_CONFIG  # lazy, heavy import
        config = _DEFAULT_CONFIG
    if calculated_features is None:
        try:
            from behav3d.widgets.utils import behav3d_calculated_features
            calculated_features = behav3d_calculated_features
        except Exception:
            calculated_features = {}

    cards: list[dict] = []

    def is_param_group(node: Any) -> bool:
        # Recurse into a dict only if it groups named sub-parameters (non-empty,
        # string keys). Empty dicts and int-keyed value maps (e.g. channel_labels
        # = {0: 'organoid'}) are treated as leaves, not recursed into.
        return (
            isinstance(node, dict)
            and len(node) > 0
            and all(isinstance(k, str) for k in node)
        )

    def walk(node: Any, parts: list[str]) -> None:
        if is_param_group(node):
            for k, v in node.items():
                walk(v, parts + [str(k)])
            return
        # leaf (scalar, list, empty dict, or an int-keyed value mapping)
        leaf = parts[-1]
        dotted = ".".join(parts)
        typ = _infer_type(node)
        step, category = _step_and_category(parts)
        cards.append({
            "key": dotted,
            "title": leaf,
            "step": step,
            "category": category,
            "type": typ,
            "default": node,
            "choices": CHOICES.get(leaf),
            "description": _describe(leaf, dotted, node, typ),
        })

    for top_key, sub in config.items():
        walk(sub, [top_key])

    # Append the computed-feature catalogue as reference cards (not tunable).
    for group, feats in (calculated_features or {}).items():
        cards.append({
            "key": f"calculated_features.{group}",
            "title": group,
            "step": "feature_extraction",
            "category": None,
            "type": "list",
            "default": list(feats),
            "choices": None,
            "description": (
                f"Computed '{group}' feature columns produced by BEHAV3D "
                f"(reference, not user-tunable): {', '.join(feats)}."
            ),
        })

    return cards


def cards_for_step(step: str, config: Optional[dict] = None) -> list[dict]:
    """Return only the cards belonging to *step* (e.g. 'tracking')."""
    return [c for c in flatten_config_to_cards(config) if c["step"] == step]


def dump_cards_json(path: str, config: Optional[dict] = None) -> int:
    """Write all cards to *path* as JSON. Returns the number of cards.

    Used by the Modal ingestion job (run once in the behav3d env) so the cloud
    side never needs to import the heavy ``behav3d.widgets.utils`` module.
    """
    import json
    cards = flatten_config_to_cards(config)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cards, f, indent=2, default=str)
    return len(cards)
