"""Versioned, feedback-grounded guidance for the BEHAV3D assistant."""
from __future__ import annotations


KNOWLEDGE_VERSION = "2026.07.27.27"


GUIDANCE_CARDS = {
    "metadata": (
        "Loaded metadata is authoritative. Read every record and validation note before asking "
        "for a value. Never ask the user to repeat a value that is present. A folder of TIFFs is "
        "not one supported movie input: each sample needs one multidimensional image file or "
        "hyperstack. Valid 5D dimension orders are TCZYX, TZCYX, ZCTYX, ZTCYX, CZTYX, and CTZYX. "
        "When a researcher describes several movies, cell types, and acquisition settings before "
        "metadata exists, open the Metadata Builder and populate all known values immediately; omit "
        "unknown fields and ask only for the one ambiguity that matters next. If asked to correct a "
        "shared acquisition value, propose that correction for every sample. When asked what "
        "analysis to run, first ask the research question; do not infer it from metadata alone. "
        "An open builder containing an unchanged copy of loaded metadata does not need to be saved "
        "again. After an external Zarr conversion or external metadata-path edit, use Load Metadata "
        "to propagate the new paths. The in-app Zarr converter reloads all tabs automatically. "
        "The Metadata Builder does not map raw channel indices to cell types. Channel-input or "
        "channel-label controls belong to the selected segmentation method. For swapped-channel "
        "replicates, generic processing population slots with the true identity recorded in each "
        "sample's Line field are a valid metadata structure; confirm that this slot-based workflow "
        "matches the intended downstream processing before editing it. Organoid lines are not "
        "necessarily separate processing types: ask whether they coexist in "
        "one movie and need separate segmentation/tracking, or should share one organoid type with "
        "line identity recorded per sample. Well and each configured population's line are mandatory; "
        "condition is optional. When no physical wells exist, propose one shared identifier such as "
        "'1'. When a configured population is confirmed absent, describe it as not added and use "
        "the literal CSV-safe line value 'not_added'; never use None. Metadata clarification "
        "prompts must stay experiment-neutral: do not introduce example strain, line, or population "
        "names that were not read from the live context, and do not interrupt an informational "
        "analysis question with a metadata-building clarification."
    ),
    "experiment_design": (
        "Use the loaded experiment reference only for the current dataset. Preserve explicit "
        "population identities, assay purpose, operational definitions, scope exclusions, and "
        "caveats. Apply this source hierarchy: live metadata is authoritative for acquisition facts "
        "and configured populations; the experiment README describes study intent and operational "
        "definitions; YAML describes saved configuration; discovered output files prove execution. "
        "When sources disagree, name the discrepancy and use live metadata unless the researcher "
        "confirms a correction. When two populations share a well and the reference identifies a paired design, "
        "prioritize the within-well comparison because dose, timing, and imaging conditions are "
        "shared. Distinguish an isogenic functional comparison from tumor-versus-healthy safety "
        "profiling; they answer different questions. Call small, unequal, or unreplicated condition "
        "comparisons descriptive or exploratory when the reference says so. Matched segmentation, "
        "feature, and filtering settings strengthen a symmetric comparison but do not by themselves "
        "prove a biological effect. Keep per-object death thresholds, population death-dynamics "
        "thresholds, and contact-associated Active Killing definitions separate. A saved parameter "
        "block shows configuration, not execution. Before saying invasiveness, Active Killing, HMM, "
        "trajectory clustering, or another readout is available, require a corresponding discovered "
        "result; otherwise say it was configured or described but no result was found."
    ),
    "reference_examples": (
        "Historical experiment profiles are provenance-labeled examples, never defaults. Mention "
        "the source profile, compare its image spacing, cadence, object scale, motion, signal quality, "
        "and method with the current experiment, and never edit a live control from an example alone. "
        "Older YAMLs may use legacy labels or units. IVM HIV example: CSV 1.15 um XY, 4 um Z, 15 s; "
        "ConvPaint Probability + Watershed; historical btrack maximum radii 12/10, optimizer distance "
        "26, time 6/4 frames; fixed 15-frame tracks and a three-state motility HMM. Its README says "
        "one sample was 10 s, conflicting with the CSV, so use live metadata. CD4/CD8-13T example: "
        "1.01 um XY, 1.05 um Z, 2 min; APOC Mask + EDT; organoid propagation; matched T-cell btrack "
        "100 maximum radius, 60 optimizer distance, 3 frames; Active Killing 7 frames, relative 2x "
        "rise, minimum 3 contact frames. Its YAML contains multiple segmentation snapshots, so those "
        "post-processing values are not a final preset. Multi-organoid safety example: 1.77 um "
        "isotropic, 2 min; paired 27T-versus-MDO targets with TEG effectors; organoid propagation; "
        "README Active Killing definition 5 frames, relative 1.5x rise, one contact frame; configured "
        "does not mean executed. Calcium reporter example: 0.33 um XY, 2 um Z, 5 s, 32 frames; external "
        "Cellpose-SAM import; Reporter Propagation for near-static intermittent detections; historical "
        "100-voxel noise and 10% overlap values; a five-state HMM on one reporter-intensity fold-change "
        "feature with smoothing 1 and a full 32-frame trajectory. Use it to explain intensity-driven "
        "behavior, not as a motility template. For any spatial btrack value, calculate or measure "
        "per-frame displacement; for optimizer values, ask for the largest spatial and missing-frame "
        "gap; for segmentation thresholds, use the current resolution and preview."
    ),
    "segmentation": (
        "Choose the segmentation method from the image signal, resolution, and available compute; "
        "the method currently shown in the dropdown may only be the UI default. Signal bleed-through "
        "is not metadata and only the researcher can confirm it, so ask whether each target is isolated "
        "in a clean channel before comparing or recommending any method. When that information is "
        "missing, ask the channel-cleanliness question and stop rather than calling the selected method "
        "a good default. Cellpose-SAM is the accuracy-first "
        "zero-shot choice only for confirmed clean, high-resolution, low-bleed-through inputs when "
        "a strong GPU or HPC is available and the number of movies/timepoints is manageable. APOC "
        "is the normal-workstation default for lower-resolution live imaging, bleed-through, or many "
        "similar experiments; its training can be reused. Try ConvPaint when APOC misses complex "
        "structures. Retrained classic Cellpose is the longest setup and is reserved for complex "
        "heterogeneous data where accuracy justifies creating ground-truth masks. The CPU Random "
        "Forest pixel classifier remains a fallback when a suitable GPU is unavailable. Import accepts "
        "existing TIFF or zarr segmentations. Image shape can provide the channel count but cannot say "
        "what biological signal is visible in each channel; ask for a channel-by-channel map and never "
        "infer signal identities or absence of bleed-through from filenames. Fluorophore names such as "
        "GFP/RFP are not needed, so do not ask for them or use them in examples. Discuss only controls for the "
        "selected method and exact cell type. For EDT advice, calculate the expected XY diameter from each "
        "sample's XY pixel size, using a 10 um cell by default and 20%, 25%, and 30% of its pixel "
        "diameter as transparent preview starting points. For organoids, first ask how many cell "
        "widths span the diameter. For a single-cell Minimum size starting point, ask for the "
        "estimated cell diameter, approximate a full spherical volume, and start at 50% of that "
        "volume so incomplete or dim segmentations are not immediately discarded. Convert to voxels "
        "using XY pixel size squared times Z pixel size when the control is in voxels, and label the "
        "result as a preview value rather than ground truth."
    ),
    "apoc": (
        "APOC has per-cell-type Image Channel Inputs. After Generate Training Data finishes, select "
        "only the channels the researcher identifies as informative for that model, including shared "
        "channels when they contain useful target signal or context. Do not add a dead channel merely "
        "as a negative class, to exclude dead cells, or to identify background: a dead target is still "
        "the same target population. Include it only when the researcher confirms that its signal is "
        "genuinely informative for that target model. If channel checkboxes are not ready, that means "
        "training data is not "
        "loaded; it does not mean APOC always uses every channel. Configure channels and classifier "
        "features before training, even when saved annotations already exist. Small structures is the "
        "normal feature preset for single cells and Large structures for organoid-scale objects. "
        "Tune Features refers to custom feature scales in pixels and a filter grid containing Gaussian "
        "blur, Difference of Gaussians, Laplacian of Gaussian, and Sobel-of-Gaussian (SoG). Small "
        "structures selects all four filters at sigma 1, 2, and 5 pixels plus original intensity; "
        "Medium structures uses 1, 2, 5, and 15; Large structures uses 1, 2, 5, 10, and 25. The "
        "original-intensity checkbox adds raw pixel intensity as a feature, so never claim APOC "
        "cannot learn from raw voxels. Tune Features never means Minimum size, EDT, Mask threshold, "
        "Seed threshold, or Feature Extraction. A classifier is trained only when a classifier "
        "file is actually found. APOC "
        "direct instance segmentation is only for sparse, non-touching objects where every "
        "connected region should remain one instance. Mask + EDT is preferred for similarly sized "
        "touching objects. Probability Map + Watershed is the default and handles heterogeneous "
        "sizes. In probability mode, Mask threshold defines the foreground contour; Seed threshold "
        "is the main splitting lever. Raise Seed threshold in small increments to split merged "
        "objects, but watch for cells disappearing if it becomes too strict. Raising Mask threshold "
        "can tighten borders and can add splitting when combined with Seed threshold. If tuning is "
        "not enough, annotate more background at touching-cell boundaries and retrain. In Mask + "
        "EDT, higher EDT generally splits more; in Peak EDT, lower EDT splits more because it is a "
        "peak-height filter. Use small opening values only and preview them; Fill holes is useful "
        "for hollow or cytoplasmic objects. Minimum size is post-processing: it excludes objects "
        "below the threshold and never merges a small object into a larger one. For organoids, prefer the initial-size "
        "filter in Filtering so segmentation and propagation retain small fragments."
    ),
    "convpaint": (
        "ConvPaint uses the same Mask + EDT and Probability Map + Watershed instance strategies as "
        "APOC. VGG16 is the sensible feature-model starting point for single cells; heavier DINO "
        "models cost more compute. Choose only channels relevant to the target. Normalization mode "
        "1 leaves intensities unchanged, mode 2 normalizes the whole stack consistently and preserves "
        "relative temporal changes, and mode 3 normalizes each image independently to correct "
        "non-biological frame-to-frame drift. Keep downsampling at 1 for single-cell precision; use "
        "it mainly for high-resolution objects roughly 100-500 pixels across. Use smoothing around "
        "1-3 only for salt-and-pepper labels. Keep iterations and tree depth at defaults unless "
        "convergence requires tuning; extra depth can overfit small annotation sets. Prediction "
        "tiling is for memory-limited large images. Dask only parallelizes tiled prediction and is "
        "still a beta path; annotation tiling is a separate training optimization. Disabled training "
        "options before training data exists are expected."
    ),
    "cellpose_sam": (
        "Cellpose-SAM requires its one-time sidecar environment, then should normally use the default "
        "GPU. CPU-only is much slower. It is zero-shot: cpsam and cpsam_v2 need no annotations, so "
        "they can be compared directly. Select the target cell type and one to three input channels; "
        "one type is segmented per run unless Run all cell types in one batch is enabled. Process all "
        "timepoints normally; a partial inclusive range preserves existing results outside the range. "
        "Leave Diameter at 0 for automatic sizing unless objects fall outside roughly 7.5-120 pixels. "
        "Flow threshold defaults to 0.4 and is ignored in 3D; raise it to keep more imperfect masks "
        "or lower it for cleaner shapes. Lower Cell probability threshold to detect more/larger "
        "objects and raise it to reject dim detections. 3D is most accurate but around 10-15 times "
        "slower; use 2D plus Stitch threshold for speed or strongly anisotropic Z data. Reduce Batch "
        "size first for GPU out-of-memory errors. Keep Remove flat segments on unless true objects "
        "occupy one slice. Size filters exclude objects and do not merge them."
    ),
    "tracking": (
        "Read current values for the exact cell type and establish observed movement between "
        "consecutive frames. For a generic tracking guide when movement has not been supplied, ask "
        "only how far the active structure moves or whether it remains overlapping, then stop without "
        "listing or evaluating any named method. Once movement is known, use Propagation for slow, "
        "non-dividing, non-touching structures that "
        "remain spatially overlapping. Use Reporter Propagation for genuinely static objects whose "
        "reporter signal flickers or disappears, such as calcium traces; it reuses one reference "
        "shape and is inappropriate for real movement or shape change. btrack is the routine default "
        "for motile cells and the only moving-object tracker in regular use. LAP and TrackPy have no "
        "identified routine advantage and should not be suggested without a specific reason. Set "
        "btrack maximum search radius from the fastest plausible one-frame displacement plus only a "
        "10-25% margin, using metadata time interval and units. Never invent a typical speed. Global "
        "optimization is Step 2, a separate refinement stage for reconnecting fragments and resolving "
        "track hypotheses. After Step 1 has been previewed and its search radius looks correct, offer "
        "to enable Step 2 as part of a complete btrack setup. Its normal starting hypotheses are false "
        "positive, initialization, termination, and linking; enable branching, death, or merging only "
        "when those events should be modeled. Distance threshold limits spatial gaps and Time threshold "
        "limits missing frames, so ask for the largest gap to bridge rather than copying defaults. "
        "Change Step size only after an "
        "out-of-memory error, lowering it to reduce RAM. When more than one organoid type exists, "
        "track all organoid types together with Propagation while preserving their origin types."
    ),
    "feature_extraction": (
        "Select feature families for the biological question. Morphology is the most expensive and "
        "should be computed only when shape matters. Movement is fast and is the default family for "
        "motility, but is usually unnecessary for organoids. The live UI makes Intensity and Contact "
        "mandatory for every population, Movement mandatory for immune/other populations, and Death "
        "mandatory when a dead channel is configured. Intensity includes mean dead-dye and channel "
        "intensities, so never suggest dropping it from T cells with a configured dead channel. Only "
        "optional groups such as Morphology and Invasiveness should be presented as removable. "
        "Contact distance 0 means strict mask "
        "touching; increase it only when proximity should count as contact. A positive distance equal "
        "to the XY pixel size permits a one-XY-pixel gap; it is not strict touching and is not a voxel "
        "diagonal. Any contact-distance "
        "change requires feature extraction to be run again, so plan it as a batch decision rather "
        "than an interactive sweep. In the dead-threshold preview, green means alive and red means "
        "dead; hovering shows the dead-mask percentage. If ambiguous, ask what looks wrong before "
        "recommending a threshold. Cell-type grouping creates a merged metadata type without "
        "requiring tracks to be recomputed."
    ),
    "active_killing": (
        "Active Killing accepts multiple selected targets in one setup. The implementation runs each "
        "target independently and also creates a combined analysis when more than one target is selected. "
        "Observation window counts forward from contact and must "
        "be calibrated to expected killing delay and the metadata time interval, not copied as a "
        "universal number. Dead-mask pixel count with an absolute increase threshold is the general "
        "default. Percentage is reasonable only when target sizes are comparable; mean dye intensity "
        "is useful for diffuse reporters or when no dead-mask segmentation exists. A multiplier is "
        "reserved for a single target line or heterogeneous baselines within one well and can be "
        "biased across target lines with different baselines. Calibrate an absolute dead-mask voxel "
        "threshold from cell size and XY and Z sampling, then inspect it visually; do not reuse 20-30 "
        "pixels blindly. Absolute-threshold mode is incomplete while its value is 0. When the user "
        "accepts a proposed setup, include targets, signal, threshold mode and value, observation "
        "window, and minimum contact duration in one proposal. Minimum contact duration must reflect "
        "imaging cadence and plausible biology."
    ),
    "filtering": (
        "Filtering must be run even when every filter is disabled: it creates the downstream CSV and "
        "interpolates missing timepoints. Minimum track length is optional because state analysis "
        "supports unequal lengths, but removing short tracks can reduce noise and computation. Read "
        "the track-length distribution before endorsing a number. Trim to a common maximum only for "
        "analyses requiring uniform windows. Equal minimum and maximum lengths are valid: the minimum "
        "discards shorter tracks and the maximum trims retained tracks. Splitting divides a long "
        "track into full-size chunks and discards the remainder. For organoids, filtering by initial "
        "size here is preferred to segmentation-time size removal. For count questions, use the "
        "unfiltered combined track-features CSV: track length is unique timepoints per sample and "
        "track ID, and only qualifying tracks present at the requested timepoint are counted. Never "
        "estimate the result."
    ),
    "analysis": (
        "First identify the research question and the active analysis view. A general analysis overview "
        "must include Death Dynamics, Interaction Analysis, Invasiveness Analysis, Active Killing, "
        "Behavioral State, State Trajectory, Contact-Based Grouping, Contact State-Shift Analysis, "
        "and Backprojection, then use the live metadata to suggest "
        "dataset-specific questions and a sensible sequence. It must not repeatedly navigate to the "
        "Analysis tab. Open a named "
        "view directly when the researcher asks. Behavioral State assigns "
        "an HMM state per timepoint. State Trajectory clusters whole trajectories using dynamic time "
        "warping. They answer different questions and their controls must not be mixed. Treat 3D "
        "morphology cautiously when Z sampling is coarse. Group cells together when populations from "
        "different channels should be analyzed jointly while retaining split results. Death Dynamics "
        "uses an organoid type plus a condition/line grouping and currently reports the mixed "
        "percentage. Interaction analysis separates live and dead targets; its overview window looks "
        "back before death, with 60 minutes as an example rather than a universal default, and an "
        "optional time range can restrict the experiment. For invasiveness, keep both the per-timepoint "
        "trace and a per-movie summary: Mean averages timepoints, while AUC integrates the trace and "
        "normalizes by elapsed time."
    ),
    "hmm": (
        "For Behavioral State, always read and state the live selected cell type before giving setup or "
        "binary-group advice. Interpret a contact group from the selected cell's perspective: for example, "
        "Macrophages contact while T-cells are selected means a T-cell touching a macrophage. Begin with "
        "a small biologically interpretable feature set. Window size "
        "5 is the default for rolling features; use 1 for genuinely single-timepoint events such as "
        "calcium peaks. Smooth window usually matches Window size and suppresses one-frame "
        "segmentation errors. Continuous features are standardized. Apply log scaling selectively "
        "only after inspecting a heavily right-skewed non-negative distribution; speed is a common "
        "zero-inflated example. Do not use percentile clipping routinely because it has degraded HMM "
        "results. Binary groups are applied post hoc to split HMM motion states; they are not HMM "
        "inputs. For movement-only analysis, enumerate all movement columns available in the live "
        "Timepoint features control and all rolling choices in Additional window features. These can "
        "include speed, displacement, cumulative displacement, displacement from origin, directional "
        "persistence, turning/reversal measures, net displacement, straightness, and mean square "
        "displacement, depending on the loaded CSV. Never imply that the two currently selected "
        "features are the complete list; offer the available choices before changing them. Use fixed "
        "state selection because automatic selection has not performed well. Keep "
        "Start offset at 1 so the undefined first-frame speed does not make every track begin static. "
        "Name, order, and color states from the feature heatmap and per-state distributions in "
        "behavioral_clustering_diagnostics.pdf. In Step 2, assigning the same biological name to "
        "multiple primary clusters merges them. Then review and rename the full behavioral clusters "
        "formed with binary groups, which can collapse unwanted state-plus-binary combinations. "
        "Those labels propagate to later reports."
    ),
    "trajectory": (
        "State Trajectory clusters whole tracks. Trajectory size cannot exceed the trim length already "
        "set in Filtering. Average linkage is the default; Complete is a reasonable comparison, while "
        "Single linkage performs poorly. If Behavioral State should remain unfiltered, trim or divide "
        "tracks here instead. Original BEHAV3D feature-based mode is deprecated and requires equal "
        "track lengths. The UMAP may look poor even when agglomerative clusters are sensible. "
        "Categorical DTW supports Contact-Based Grouping, which compares trajectory clusters for "
        "tracks with and without a sufficiently long contact bout. It also supports Contact State-Shift "
        "Analysis, which compares behavioral-state composition before and after contact against "
        "timing-matched no-contact tracks and therefore additionally requires Behavioral State results. "
        "Do not repeat the obsolete claim that these current contact analyses inherently produce empty "
        "output."
    ),
}


def select_guidance_cards(
    context: dict,
    user_message: str = "",
    intent: str | None = None,
) -> list[dict]:
    """Select deterministic cards before generic vector retrieval."""
    step = str(context.get("current_step") or "")
    query = f"{user_message} {intent or ''}".lower()
    selected: list[str] = []
    if step in GUIDANCE_CARDS:
        selected.append(step)
    if context.get("experiment_reference"):
        selected.append("experiment_design")
    example_request = any(phrase in query for phrase in (
        "example value", "example setting", "example configuration",
        "previous experiment", "past experiment", "historical value",
        "historical setting", "reference profile", "reference configuration",
        "similar experiment", "similar dataset", "values used before",
        "what did you use", "used previously", "prior experiment",
    ))
    if example_request:
        selected.append("reference_examples")

    segmentation_method = str(
        (context.get("segmentation") or {}).get("method") or ""
    ).lower()
    for token, card_id in (
        ("apoc", "apoc"),
        ("convpaint", "convpaint"),
        ("cellpose-sam", "cellpose_sam"),
    ):
        if token in segmentation_method and card_id not in selected:
            selected.append(card_id)

    special_keywords = {
        "apoc": (
            "apoc", "probability map", "mask threshold", "seed threshold",
            "opening", "fill holes",
        ),
        "convpaint": ("convpaint", "vgg16", "dinov2", "dask"),
        "cellpose_sam": ("cellpose-sam", "cellpose sam", "cpsam", "zero-shot"),
        "active_killing": (
            "active killing", "observation window", "top killers",
            "killing threshold", "effector",
        ),
        "hmm": (
            "hmm", "behavioral state", "behavioural state", "state selection",
            "smooth window", "log scaling", "percentile clipping",
        ),
        "trajectory": (
            "state trajectory", "trajectory clustering", "linkage",
            "divide long tracks", "original behav3d", "exemplar tracks",
        ),
    }
    for key, words in special_keywords.items():
        if key not in selected and any(word in query for word in words):
            selected.append(key)

    view = str((context.get("analysis") or {}).get("view") or "")
    if view == "behavioral_state" and "hmm" not in selected:
        selected.append("hmm")
    elif view == "state_trajectory" and "trajectory" not in selected:
        selected.append("trajectory")
    if (
        (context.get("feature_extraction") or {}).get("active_killing_open")
        and "active_killing" not in selected
    ):
        selected.append("active_killing")

    keyword_map = {
        "metadata": ("metadata", "sample", "dimension", "image path", "tiff", "voxel"),
        "segmentation": ("segment", "cellpose", "random forest"),
        "tracking": (
            "tracking", "btrack", "branch", "visual feature", "search radius",
            "propagation", "movement", "motion",
        ),
        "feature_extraction": (
            "feature extraction", "feature group", "contact threshold",
            "death preview", "dead threshold",
        ),
        "filtering": (
            "filter", "track length", "short track", "distribution",
            "cell count", "timepoint",
        ),
        "analysis": ("analysis", "cluster", "heatmap"),
    }
    for key, words in keyword_map.items():
        if key not in selected and any(word in query for word in words):
            selected.append(key)
    return [{"id": key, "text": GUIDANCE_CARDS[key]} for key in selected]
