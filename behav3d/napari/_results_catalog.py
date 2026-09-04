"""
BEHAV3D napari plugin – Results catalog.

Walks ``output_dir/analysis/`` once and returns a flat list of
:class:`ResultFile` entries grouped by *category* (Filtering, Feature
Extraction, Analysis) and *subcategory* (Death Dynamics, Single Cell).

Categorisation is rule-based:

1. ``quality_control/`` in path parts -> ``filtering``.
2. ``track_features/`` in path parts -> ``feature_extraction``.
3. ``behavioral_states`` / ``track_classification`` /
   ``state_classification`` / ``behavioral_clustering`` in path parts,
   OR filename matches single-cell clustering patterns
   (``BEHAV3D_*_UMAP_*``, ``state_*``, ``behavioral_clustering_*``,
   ``track_class_*``, ``clustering_diagnostics_*``,
   ``example_tracks_overview*``, ``*_confusion_matrices*``,
   ``transition_matrix_heatmap``, ``sankey_all_pairs``) ->
   ``analysis / single_cell``.
4. Everything else under ``analysis/`` -> ``analysis / death_dynamics``.

The "cell type" bucket is simply the first path component under
``analysis/`` (e.g. ``tcell``, ``organoid1``, ``morpho-dead``,
``multi_organoid_comparison``).

Each :class:`ResultFile` carries a human-readable ``description`` looked
up in :data:`FILE_CATALOG`. The catalog is intentionally a list of
``(glob_pattern, description)`` tuples so order matters: the first
match wins and we list specific patterns before generic fallbacks.
"""
from __future__ import annotations

import fnmatch
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


# ---------------------------------------------------------------------------
# Catalog
# ---------------------------------------------------------------------------
# Filename-pattern -> human description. Matched in order; first match wins.
FILE_CATALOG: list[tuple[str, str]] = [
    # ── Filtering / QC ────────────────────────────────────────────────────
    (
        "BEHAV3D_filter_counts.pdf",
        "Bar plots of track counts after each filtering step (per sample).",
    ),
    (
        "BEHAV3D_*_track_count_summary.csv",
        "Standard per-sample cell-count comparison at timepoint 0 for minimum "
        "track lengths of 20, 50, 100, and 200.",
    ),
    (
        "BEHAV3D_*_track_count_query_t*_min_*.csv",
        "Per-sample cell counts at one timepoint before and after minimum "
        "track-length thresholds.",
    ),
    (
        "BEHAV3D_*_touching_distribution.pdf",
        "Distribution of contacts between tracks and the given target cell"
        " type, used to tune contact thresholds.",
    ),
    (
        "BEHAV3D_dead_dye_distribution.pdf",
        "Distribution of dead-dye intensity at all timepoints and at t=1,"
        " used to tune the dead-dye threshold.",
    ),

    # ── Feature Extraction ────────────────────────────────────────────────
    (
        "BEHAV3D_*_combined_track_features.csv",
        "Per-timepoint features for every track, combined across samples"
        " (raw, pre-filter).",
    ),
    (
        "BEHAV3D_*_combined_track_features_filtered.csv",
        "Per-timepoint features after filtering.",
    ),
    (
        "BEHAV3D_*_combined_track_features_summarized.csv",
        "Track-level summary statistics computed from timepoint features.",
    ),
    (
        "BEHAV3D_*_combined_track_features_clustered.csv",
        "Track-level features with UMAP cluster assignment.",
    ),
    (
        "BEHAV3D_*_advanced_track_features.csv",
        "Advanced (active-killing-aware) per-track features.",
    ),
    (
        "*_intensity.csv",
        "Per-timepoint intensity features for a single sample (intermediate).",
    ),
    (
        "*_contact.csv",
        "Per-timepoint contact features for a single sample (intermediate).",
    ),
    (
        "*_dead_mask.csv",
        "Per-timepoint dead-mask features for a single sample (intermediate).",
    ),
    (
        "*_morphology.csv",
        "Per-timepoint morphology features for a single sample (intermediate).",
    ),
    (
        "*_track_features.csv",
        "Per-track features for a single sample (intermediate).",
    ),

    # ── Analysis: Invasiveness ────────────────────────────────────────────
    (
        "invasiveness_analysis_*.pdf",
        "Invasiveness analysis for an immune cell type vs one or more"
        " targets: fraction / percentage over time + per-movie summary.",
    ),
    (
        "invasiveness_fraction_over_time_*.csv",
        "Per-(sample, timepoint) fraction of invasive immune cells per target.",
    ),
    (
        "invasiveness_perc_over_time_*.csv",
        "Per-(sample, timepoint) mean/median invasiveness percentage per target.",
    ),
    (
        "invasiveness_per_movie_summary_*.csv",
        "Per-movie invasiveness summary (one row per sample x target).",
    ),

    # ── Analysis: Death Dynamics ──────────────────────────────────────────
    (
        "interaction_analysis_*_vs_*.pdf",
        "Interaction analysis between a target cell type and an interaction"
        " cell type: contact / fate / cumulative-to-death plots.",
    ),
    (
        "multi_organoid_interaction_comparison*.pdf",
        "Interaction Overview: violin, before-death curves, and active-killing "
        "dashboard (one or more organoid types).",
    ),
    (
        "*_organoid_analysis.pdf",
        "Per-sample organoid death-dynamics analysis.",
    ),
    (
        "combined_general_*_dynamics_analysis.pdf",
        "Combined general organoid dynamics analysis for one organoid type.",
    ),
    (
        "multi_organoid_death_dynamics_comparison.pdf",
        "Cross-organoid comparison of death dynamics.",
    ),
    (
        "morpho_dead_umap.pdf",
        "UMAP of organoid morpho-dead features colored by cluster +"
        " per-feature panels.",
    ),
    (
        "morpho_dead_heatmap.pdf",
        "Cluster x feature heatmap + violin plots for the morpho-dead"
        " clustering.",
    ),
    (
        "morpho_dead_dynamics.pdf",
        "Death dynamics per organoid colored by morpho-dead cluster.",
    ),
    (
        "morpho_dead_snapshots.pdf",
        "Image crops per cluster at window start/end for the morpho-dead"
        " clustering.",
    ),
    (
        "morpho_dead_composition.pdf",
        "Stacked composition bar plots of morpho-dead clusters per condition.",
    ),
    (
        "morpho_dead_*.csv",
        "Per-organoid morpho-dead results and cluster summaries.",
    ),

    # ── Analysis: Single Cell (clustering / states / tracks) ─────────────
    (
        "BEHAV3D_*_UMAP_clusters.pdf",
        "Clustered UMAP with backprojected per-track features.",
    ),
    (
        "BEHAV3D_*_UMAP_cluster_feature_heatmap.pdf",
        "Heatmap of mean feature values per UMAP cluster.",
    ),
    (
        "BEHAV3D_*_UMAP_cluster_percentages*.pdf",
        "Barplots of cluster percentages per sample / per condition.",
    ),
    (
        "BEHAV3D_*_UMAP_clusters.csv",
        "UMAP coordinates and cluster assignment for each track.",
    ),
    (
        "behavioral_clustering_diagnostics.pdf",
        "Diagnostic plots for the behavioral-state clustering (UMAP,"
        " silhouette, etc.).",
    ),
    (
        "behavioral_clustering_*_proportions.pdf",
        "Cluster vs. binary-group proportions for the behavioral-state"
        " clustering.",
    ),
    (
        "state_classification_diagnostics.pdf",
        "Diagnostics of the state classification model.",
    ),
    (
        "state_composition_report_*.pdf",
        "State composition report for a given grouping.",
    ),
    (
        "transition_matrix_heatmap.pdf",
        "Heatmap of state-to-state transition probabilities.",
    ),
    (
        "sankey_all_pairs.pdf",
        "Sankey diagrams of state transitions for all pairs.",
    ),
    (
        "*_confusion_matrices.pdf",
        "Confusion matrices for a classifier.",
    ),
    (
        "track_class_proportions_by_sample_*.pdf",
        "Per-sample proportions of track classes.",
    ),
    (
        "clustering_diagnostics_*.pdf",
        "Clustering diagnostics for a given cluster key.",
    ),
    (
        "example_tracks_overview*.pdf",
        "Exemplar track previews per cluster (renamed/classified variants).",
    ),

    # ── Generic fallbacks (always last) ───────────────────────────────────
    ("*.pdf", "PDF output (no specific description registered)."),
    ("*.csv", "Tabular data (CSV)."),
    ("*.tsv", "Tabular data (TSV)."),
    ("*.png", "Static image output."),
    ("*.jpg", "Static image output."),
    ("*.jpeg", "Static image output."),
    ("*.tif", "Static image output (TIFF)."),
    ("*.tiff", "Static image output (TIFF)."),
    ("*.html", "Interactive HTML output (opens externally)."),
    ("*.json", "JSON data."),
    ("*.yaml", "YAML configuration."),
    ("*.yml", "YAML configuration."),
    ("*.zarr", "Zarr image store."),
]


def describe(path: Path) -> str:
    """Return the catalog description for ``path``'s filename."""
    name = path.name
    for pattern, desc in FILE_CATALOG:
        if fnmatch.fnmatch(name, pattern):
            return desc
    return ""


# ---------------------------------------------------------------------------
# Result file model
# ---------------------------------------------------------------------------
_VIEWABLE_KINDS = {"pdf", "image"}

_IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}


def _kind_for(path: Path) -> str:
    ext = path.suffix.lower()
    if ext == ".pdf":
        return "pdf"
    if ext in _IMAGE_EXTS:
        return "image"
    if ext in {".csv", ".tsv"}:
        return "csv"
    if ext in {".html", ".htm"}:
        return "html"
    if ext == ".zarr" or (path.is_dir() and ext == ".zarr"):
        return "zarr"
    return "other"


@dataclass(frozen=True)
class ResultFile:
    path: Path
    label: str          # short, e.g. "interaction_analysis_tcell_vs_organoid1.pdf"
    description: str    # tooltip, from FILE_CATALOG (or "")
    kind: str           # "pdf" | "image" | "csv" | "html" | "zarr" | "other"
    category: str       # "filtering" | "feature_extraction" | "analysis"
    subcategory: Optional[str]  # "death_dynamics" | "single_cell" | None
    cell_type: Optional[str]    # first folder under analysis/, e.g. "tcell"

    @property
    def is_viewable(self) -> bool:
        return self.kind in _VIEWABLE_KINDS


# ---------------------------------------------------------------------------
# Categorisation rules
# ---------------------------------------------------------------------------
_SINGLE_CELL_DIRS = {
    "behavioral_states",
    "behavioral_clustering",
    "track_classification",
    "state_classification",
}
_SINGLE_CELL_NAME_PATTERNS = (
    "BEHAV3D_*_UMAP_*",
    "state_*",
    "behavioral_clustering_*",
    "track_class_*",
    "clustering_diagnostics_*",
    "example_tracks_overview*",
    "*_confusion_matrices*",
    "transition_matrix_heatmap.*",
    "sankey_all_pairs.*",
)


def _classify(path: Path, analysis_root: Path) -> tuple[str, Optional[str], Optional[str]]:
    """Return ``(category, subcategory, cell_type)`` for ``path``."""
    try:
        rel = path.relative_to(analysis_root)
    except ValueError:
        return ("analysis", "death_dynamics", None)

    parts = rel.parts
    cell_type = parts[0] if parts else None

    in_sc_dir = any(p in _SINGLE_CELL_DIRS for p in parts)
    name = path.name
    matches_sc_name = any(
        fnmatch.fnmatch(name, pat) for pat in _SINGLE_CELL_NAME_PATTERNS
    )
    if in_sc_dir or matches_sc_name:
        return ("analysis", "single_cell", cell_type)

    if "quality_control" in parts:
        return ("filtering", None, cell_type)
    if "track_features" in parts:
        return ("feature_extraction", None, cell_type)

    return ("analysis", "death_dynamics", cell_type)


# ---------------------------------------------------------------------------
# Walking
# ---------------------------------------------------------------------------
_SKIP_DIR_NAMES = {
    "__pycache__",
    ".cache",
    ".ipynb_checkpoints",
}

# Extensions we care about. Other files (logs, temp, etc.) are skipped.
_ALLOWED_EXTS = {
    ".pdf",
    ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp",
    ".csv", ".tsv",
    ".html", ".htm",
    ".json", ".yaml", ".yml",
    ".txt", ".md",
}


def scan_outputs(output_dir: Path) -> list[ResultFile]:
    """Walk ``output_dir/analysis/`` and return discovered files.

    Returns an empty list when ``output_dir`` does not contain an
    ``analysis`` subfolder yet.
    """
    output_dir = Path(output_dir).expanduser()
    analysis_root = output_dir / "analysis"
    if not analysis_root.is_dir():
        return []

    results: list[ResultFile] = []
    zarr_dirs: list[Path] = []
    for root, dirs, files in os.walk(analysis_root):
        # Prune skip dirs and don't descend into .zarr stores (huge, opaque),
        # collecting the stores as we go. The walk already visits every .zarr
        # by name here, so a second ``rglob("*.zarr")`` pass would re-traverse
        # the whole tree -- and, unlike this loop, would descend *into* every
        # store and stat each of its chunk files.
        keep = []
        for d in dirs:
            if d in _SKIP_DIR_NAMES:
                continue
            if d.lower().endswith(".zarr"):
                zarr_dirs.append(Path(root, d))
            else:
                keep.append(d)
        dirs[:] = keep
        for fname in files:
            ext = os.path.splitext(fname)[1].lower()
            if ext not in _ALLOWED_EXTS:
                continue
            p = Path(root, fname)
            category, subcategory, cell_type = _classify(p, analysis_root)
            results.append(
                ResultFile(
                    path=p,
                    label=fname,
                    description=describe(p),
                    kind=_kind_for(p),
                    category=category,
                    subcategory=subcategory,
                    cell_type=cell_type,
                )
            )

    # Also surface .zarr stores at directory level (for "Reveal in folder"),
    # from the directories the walk above already pruned.
    for child in zarr_dirs:
        category, subcategory, cell_type = _classify(child, analysis_root)
        results.append(
            ResultFile(
                path=child,
                label=child.name,
                description=describe(child),
                kind="zarr",
                category=category,
                subcategory=subcategory,
                cell_type=cell_type,
            )
        )

    results.sort(key=lambda r: (str(r.path).lower(),))
    return results


# ---------------------------------------------------------------------------
# Convenience: friendly subcategory / category labels
# ---------------------------------------------------------------------------
CATEGORY_LABELS = {
    "filtering": "Filtering",
    "feature_extraction": "Feature Extraction",
    "analysis": "Analysis",
}

SUBCATEGORY_LABELS = {
    "death_dynamics": "Death Dynamics",
    "single_cell": "Single Cell",
}


def group_by_tree(
    files: Iterable[ResultFile],
) -> dict[str, dict[Optional[str], dict[Optional[str], list[ResultFile]]]]:
    """Group files as ``{category: {subcategory: {cell_type: [files]}}}``.

    ``None`` keys are used where the dimension does not apply
    (e.g. ``subcategory`` is ``None`` for filtering / feature
    extraction).
    """
    tree: dict[str, dict[Optional[str], dict[Optional[str], list[ResultFile]]]] = {}
    for f in files:
        cat = tree.setdefault(f.category, {})
        sub = cat.setdefault(f.subcategory, {})
        bucket = sub.setdefault(f.cell_type, [])
        bucket.append(f)
    return tree


__all__ = [
    "FILE_CATALOG",
    "CATEGORY_LABELS",
    "SUBCATEGORY_LABELS",
    "ResultFile",
    "describe",
    "scan_outputs",
    "group_by_tree",
]
