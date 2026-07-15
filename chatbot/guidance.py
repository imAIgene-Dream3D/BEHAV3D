"""Versioned, feedback-grounded guidance for the BEHAV3D assistant."""
from __future__ import annotations


KNOWLEDGE_VERSION = "2026.07.14.2"


GUIDANCE_CARDS = {
    "metadata": (
        "Loaded metadata is authoritative. Read every record and validation note before asking "
        "for a value. Never ask the user to repeat a value that is present. A folder of TIFFs is "
        "not one supported movie input: each sample needs one multidimensional image file or "
        "hyperstack. Valid 5D dimension orders are TCZYX, TZCYX, ZCTYX, ZTCYX, CZTYX, and CTZYX. "
        "When asked what analysis to run, first ask the research question; do not infer it from "
        "metadata alone."
    ),
    "segmentation": (
        "Method ladder: start with APOC for a fast GPU-assisted interactive baseline. Use ConvPaint "
        "for more complex appearance when extra training and compute are acceptable. Use Pixel "
        "Classifier (Random Forest) when no suitable GPU is available. Use a custom-trained Cellpose "
        "model when repeatable accuracy across experiments justifies model training. Import accepts "
        "existing TIFF or zarr segmentations. After a method is selected, discuss only controls whose "
        "method matches that selection. APOC is not the Random Forest pixel classifier. For EDT advice, "
        "calculate the expected XY diameter from each sample's pixel_distance_xy. A 10 um cell is the "
        "default reference. Try 20%, 25%, and 30% of the diameter in pixels as transparent preview "
        "starting points. For organoids, first ask how many cell widths span the diameter."
    ),
    "tracking": (
        "Read the current values for the exact cell type before recommending changes. btrack maximum "
        "search radius and optimizer distance threshold are physical distances after metadata scaling, "
        "normally micrometres, not pixels. Use visual features adds image-derived measurements to "
        "linking; global optimization is a separate second stage. Discuss P_branch only when global "
        "optimization is already enabled and divisions matter. For customization, copy the bundled "
        "cell_config.json and edit the copy. Describe tracking failures as hypotheses to check, not "
        "certain diagnoses."
    ),
    "feature_extraction": (
        "Cell-type grouping is available in Feature Extraction and creates a merged metadata type; "
        "existing tracks do not need to be recomputed. In the dead-threshold preview, green means alive "
        "and red means dead; hovering shows the dead-mask percentage. If the preview is ambiguous, ask "
        "what looks wrong before recommending a threshold."
    ),
    "filtering": (
        "Use the track-length distribution preview before choosing a minimum length. A higher minimum "
        "removes short unreliable tracks but can exclude real short-lived behavior. Visual accuracy and "
        "the downstream analysis matter more than forcing equal track lengths. Trim to equal lengths only "
        "when the selected downstream workflow requires comparable windows. For count questions, use "
        "the unfiltered combined track-features CSV. Track length is the number of unique position_t "
        "values per sample and TrackID. Count only qualifying tracks that are present at the requested "
        "position_t; never estimate the result."
    ),
    "analysis": (
        "State Classification uses an HMM on per-timepoint and rolling-window features. Track "
        "Classification contains whole-trajectory clustering, including legacy DTW; do not conflate them. "
        "Begin with a small set of biologically interpretable features. Treat 3D morphology cautiously "
        "when Z sampling is coarse. Window size controls rolling feature calculation; Smooth window "
        "controls feature smoothing. All continuous features are standardized; log1p is optional for "
        "strongly right-skewed non-negative features. HMM diagnostics and the state-feature heatmap are in "
        "analysis/<cell_type>/behavioral_states/processing/hmm_behavioral_classification/quality_control/"
        "raw/behavioral_clustering_diagnostics.pdf."
    ),
}


def select_guidance_cards(context: dict, user_message: str = "", intent: str | None = None) -> list[dict]:
    """Select deterministic cards before generic vector retrieval."""
    step = str(context.get("current_step") or "")
    query = f"{user_message} {intent or ''}".lower()
    selected: list[str] = []
    if step in GUIDANCE_CARDS:
        selected.append(step)
    keyword_map = {
        "metadata": ("metadata", "sample", "dimension", "image path", "tiff", "voxel"),
        "segmentation": ("segment", "apoc", "convpaint", "cellpose", "random forest", "edt"),
        "tracking": ("track", "btrack", "branch", "visual feature", "search radius"),
        "feature_extraction": ("feature", "group", "death", "dead threshold", "preview"),
        "filtering": ("filter", "track length", "short track", "distribution", "cell count", "timepoint"),
        "analysis": ("analysis", "hmm", "dtw", "state", "cluster", "log", "heatmap"),
    }
    for key, words in keyword_map.items():
        if key not in selected and any(word in query for word in words):
            selected.append(key)
    return [{"id": key, "text": GUIDANCE_CARDS[key]} for key in selected]
