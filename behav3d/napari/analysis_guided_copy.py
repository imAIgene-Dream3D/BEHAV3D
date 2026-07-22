"""Plain-language, biologist-facing copy for the Analysis tab's Guided mode.

Guided mode is an *explanation-first on-ramp*, not an auto-solver: it teaches
what each analysis does (especially the concept of clustering) and is honest
about the data-dependent choices the user must make, then hands them into the
real settings form.

Keeping the copy here (separate from layout code in ``_guided.py`` and the two
sub-tab modules) means the wording can be edited without touching any widget
code. The vocabulary is deliberately **cell-type-neutral** — "target" and
"interacting" cells are roles, mapping to whatever cell types the user defined;
they are not hardcoded to immune-vs-organoid.

Each analysis spec is a dict with:

- ``id``        — stable key (not shown; used for Start wiring)
- ``title``     — the real analysis name, unchanged (shown as the heading)
- ``subtitle``  — the plain-language question (shown under the title)
- ``color``     — accent hex for the leading marker
- ``what_does`` — one short paragraph: what the analysis does
- ``concept``   — optional {"term", "text"} callout explaining a jargon term
- ``what_get``  — one line: what the user receives
- ``decide``    — list of {"label", "tag", "kind"} the user will choose
- ``seed``      — the question sent to the assistant by "Ask the assistant"
- ``start_label`` — optional override for the Start button text (defaults to
  "Start — open the settings  ▸" in AnalysisExplainer)
- ``has_params``  — optional; only meaningful for the nested pipeline cards
  consumed by StateClassificationSubTab / TrackClassificationSubTab. False
  means the pipeline has no configurable settings, so its Start button
  triggers the existing run button directly instead of opening an isolated
  settings view.

``kind`` is one of ``required`` / ``suggested`` / ``estimated`` / ``default``
and controls how honestly each choice is labelled and coloured.
"""
from __future__ import annotations


# Decision-tier tags reused across analyses so the honesty labels stay
# consistent everywhere.
_ONLY_YOU = {"tag": "only you know this", "kind": "required"}
_SUGGESTED = {"tag": "suggested set, adjustable", "kind": "suggested"}
_ESTIMATED = {"tag": "can be estimated", "kind": "estimated"}
_DEFAULTS = {"tag": "tested defaults, optional", "kind": "default"}


def _decide(label: str, base: dict) -> dict:
    return {"label": label, **base}


DEATH_DYNAMICS = {
    "id": "death_dynamics",
    "title": "Death Dynamics",
    "subtitle": "How and when do target cells die?",
    "color": "#d85a30",
    "what_does": (
        "Tracks the fraction of your target population that dies over the "
        "course of each movie, so you can compare killing between conditions."
    ),
    "concept": None,
    "what_get": "A death-curve PDF per target (and a combined view across targets).",
    "decide": [
        _decide("Which target cell type(s)", _ONLY_YOU),
        _decide("Time range to analyse", _ESTIMATED),
    ],
    "seed": (
        "Explain BEHAV3D Death Dynamics analysis for a biologist: what it "
        "measures, what the death curve means, and how the dead-cell threshold "
        "is decided."
    ),
}

INTERACTION = {
    "id": "interaction",
    "title": "Interaction Analysis",
    "subtitle": "Do interacting cells contact targets before death?",
    "color": "#534ab7",
    "what_does": (
        "Measures contact between an interacting cell type and target cells in "
        "the window before each target death — i.e. whether killing is contact "
        "dependent."
    ),
    "concept": None,
    "what_get": "An interaction PDF per target, plus a cross-condition overview.",
    "decide": [
        _decide("Which target and interacting cell types", _ONLY_YOU),
        _decide("Time window before death", _ESTIMATED),
        _decide("How to group conditions", _SUGGESTED),
    ],
    "seed": (
        "Explain BEHAV3D Interaction Analysis for a biologist: what contact/"
        "interaction it measures, and how to choose the time window before "
        "death."
    ),
}

INVASIVENESS = {
    "id": "invasiveness",
    "title": "Invasiveness Analysis",
    "subtitle": "How deep do interacting cells penetrate?",
    "color": "#0f6e56",
    "what_does": (
        "Scores how far an interacting cell type penetrates into the target "
        "structure over the movie, summarised per movie."
    ),
    "concept": None,
    "what_get": "An invasiveness PDF summarising penetration per condition.",
    "decide": [
        _decide("Which interacting and target cell types", _ONLY_YOU),
        _decide("Per-movie summary statistic", _SUGGESTED),
        _decide("Time range to analyse", _ESTIMATED),
    ],
    "seed": (
        "Explain BEHAV3D Invasiveness Analysis for a biologist: what "
        "penetration it measures and which summary statistic to pick."
    ),
}

BEHAVIORAL_STATE = {
    "id": "state",
    "title": "Behavioral State",
    "subtitle": "What behaviour states do cells pass through?",
    "color": "#993556",
    "what_does": (
        "Looks at how each cell behaves moment-to-moment and finds a few "
        "recurring \"states\" it switches between over time — for example a "
        "resting, a moving, and an interacting state. You don't define the "
        "states by hand; the analysis discovers them from the data."
    ),
    "concept": {
        "term": "\"Clustering\" means",
        "text": (
            "the computer sorts many cells into a handful of look-alike groups, "
            "so instead of eyeballing thousands of tracks you read a few "
            "well-defined behaviour types."
        ),
    },
    "what_get": (
        "State maps overlaid on your cells, and a breakdown of how much time "
        "each condition spends in each state."
    ),
    "decide": [
        _decide("Which cell type", _ONLY_YOU),
        _decide("Which features describe behaviour", _SUGGESTED),
        _decide("How many states", _ESTIMATED),
        _decide("Expert statistics (covariance, iterations, seed…)", _DEFAULTS),
    ],
    "seed": (
        "Explain BEHAV3D behaviour-state analysis for a biologist: what it "
        "does, which features it uses, and how to choose the number of states."
    ),
}

STATE_TRAJECTORY = {
    "id": "track",
    "title": "State Trajectory",
    "subtitle": "What movement patterns exist across whole tracks?",
    "color": "#0f6e56",
    "what_does": (
        "Groups whole cell tracks by the shape of their behaviour over time — "
        "instead of moment-to-moment states, it sorts each cell's entire path "
        "into a few \"movement types\"."
    ),
    "concept": {
        "term": "Same idea, whole track",
        "text": (
            "it clusters each cell's complete trajectory, so look-alike paths "
            "end up in the same group and you read a few movement types instead "
            "of every individual track."
        ),
    },
    "what_get": "Cluster and exemplar-track plots, and per-condition composition.",
    "decide": [
        _decide("Which cell type", _ONLY_YOU),
        _decide("How many clusters", _ESTIMATED),
        _decide("Expert statistics (UMAP, linkage, classifier…)", _DEFAULTS),
    ],
    "seed": (
        "Explain BEHAV3D track/state-trajectory analysis for a biologist: what "
        "it groups, what UMAP and clustering do here, and how to choose the "
        "number of clusters."
    ),
}


# ── Single-cell "extendable plotting" sub-sections ──────────────────────────
# These are nested one level deeper than the top-level cards above: an entry
# card (Level 0) leads to a list of per-pipeline cards (Level 1); pipeline
# cards with configurable parameters lead to that pipeline's isolated
# settings (Level 2), while parameter-free pipelines run immediately from
# Level 1. See `has_params` below — StateClassificationSubTab and
# TrackClassificationSubTab read it to decide whether Start opens settings
# or triggers the existing run button directly.

STATE_REPORTS_ENTRY = {
    "id": "state_reports",
    "title": "Reports & Plots",
    "subtitle": "Turn your classified behavioral states into shareable figures.",
    "color": "#4a90d9",
    "show_explainer": False,
    "what_does": (
        "State composition, state transition, and condition comparison "
        "reports are all built from the behavioral states you classified in "
        "Step 2 — pick which one(s) you need."
    ),
    "concept": None,
    "what_get": (
        "A choice of report PDFs: state composition, state transitions, "
        "and/or a two-condition statistical comparison."
    ),
    "decide": [],
    "start_label": "Generate analysis and plots  ▸",
    "seed": (
        "Explain what reports are available after BEHAV3D behavioral state "
        "classification and when to use each one."
    ),
}

STATE_COMPOSITION_REPORT = {
    "id": "state_composition",
    "title": "State Composition Report",
    "subtitle": "How much time does each condition spend in each state?",
    "color": "#4a90d9",
    "what_does": (
        "Plots the proportion of each behavioral state per sample, "
        "optionally grouped by one or more metadata columns."
    ),
    "concept": None,
    "what_get": "A composition PDF showing the state breakdown per sample/group.",
    "decide": [
        _decide("Which column(s) to group composition plots by", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain the BEHAV3D State Composition Report: what it plots and "
        "how to choose grouping columns."
    ),
}

STATE_TRANSITION_REPORT = {
    "id": "state_transition",
    "title": "State Transition Report",
    "subtitle": "Which states do cells switch between, and how often?",
    "color": "#4a90d9",
    "what_does": (
        "Plots state-to-state transition probabilities, showing which "
        "behavioral states cells move into and out of."
    ),
    "concept": None,
    "what_get": "A transition-probability plot — no further configuration needed.",
    "decide": [],
    "has_params": False,
    "start_label": "Generate State Transition plots  ▸",
    "seed": (
        "Explain the BEHAV3D State Transition Report: what a state-to-state "
        "transition probability plot shows."
    ),
}

STATE_COMPARISON_REPORT = {
    "id": "state_comparison",
    "title": "Condition Comparison Report",
    "subtitle": "Does state proportion differ significantly between two conditions?",
    "color": "#4a90d9",
    "what_does": (
        "Compares overall state proportions between two levels of a "
        "condition column using Welch's t-test, optionally faceted by "
        "another column."
    ),
    "concept": {
        "term": "Welch's t-test",
        "text": (
            "compares two group means without assuming equal variances — "
            "safer than a standard t-test when group sizes or spreads differ."
        ),
    },
    "what_get": "A condition-comparison PDF with significance annotations.",
    "decide": [
        _decide("Which condition column and levels to compare", _ONLY_YOU),
        _decide("Column to facet by", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain the BEHAV3D state Condition Comparison Report: what "
        "Welch's t-test compares here and how to choose condition/levels/facet."
    ),
}

TRACK_PLOTS_ENTRY = {
    "id": "track_plots",
    "title": "Create Plots",
    "subtitle": "Turn your track clustering into diagnostics and figures.",
    "color": "#c98a2c",
    "show_explainer": False,
    "what_does": (
        "Diagnostics, track proportions, condition comparison, "
        "contact-based grouping, and exemplar tracks are all built from the "
        "track clustering above — pick which one(s) you need."
    ),
    "concept": None,
    "what_get": (
        "A choice of diagnostic and summary PDFs, plus optional "
        "exemplar-track PDFs/MP4s, for your track clustering."
    ),
    "decide": [],
    "start_label": "Generate analysis and plots  ▸",
    "seed": (
        "Explain what plots and reports are available after BEHAV3D "
        "track/state-trajectory clustering and when to use each one."
    ),
}

TRACK_DIAGNOSTICS = {
    "id": "track_diagnostics",
    "title": "Diagnostic",
    "subtitle": "How well-separated are your trajectory clusters?",
    "color": "#c98a2c",
    "what_does": (
        "Runs quality-control diagnostics on the trajectory clustering "
        "(silhouette/distance plots), so you can judge whether the clusters "
        "are well separated."
    ),
    "concept": None,
    "what_get": "Diagnostic PDF(s) with silhouette and distance plots.",
    "decide": [],
    "has_params": False,
    "start_label": "Generate Diagnostics plots  ▸",
    "seed": (
        "Explain the BEHAV3D track-clustering diagnostics: what a "
        "silhouette plot tells you about cluster quality."
    ),
}

TRACK_PROPORTIONS = {
    "id": "track_proportions",
    "title": "Track Proportions",
    "subtitle": "How do movement-type proportions vary across samples?",
    "color": "#c98a2c",
    "what_does": (
        "Plots how track cluster proportions vary across samples, "
        "optionally grouped by one or more metadata columns."
    ),
    "concept": None,
    "what_get": "A track-proportion PDF per grouping.",
    "decide": [
        _decide("Which column(s) to group proportion plots by", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain the BEHAV3D Track Proportions plot: what it shows and how "
        "to choose grouping columns."
    ),
}

TRACK_COMPARISON_REPORT = {
    "id": "track_comparison",
    "title": "Condition Comparison Report",
    "subtitle": "Does track-cluster proportion differ significantly between two conditions?",
    "color": "#c98a2c",
    "what_does": (
        "Compares overall track cluster proportions between two levels of "
        "a condition column using Welch's t-test, optionally faceted by "
        "another column."
    ),
    "concept": {
        "term": "Welch's t-test",
        "text": (
            "compares two group means without assuming equal variances — "
            "safer than a standard t-test when group sizes or spreads differ."
        ),
    },
    "what_get": "A condition-comparison PDF with significance annotations.",
    "decide": [
        _decide("Which condition column and levels to compare", _ONLY_YOU),
        _decide("Column to facet by", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain the BEHAV3D track Condition Comparison Report: what "
        "Welch's t-test compares here and how to choose condition/levels/facet."
    ),
}

TRACK_CONTACT_GROUPING = {
    "id": "track_contact",
    "title": "Contact-based Grouping",
    "subtitle": "Do tracks with a sustained contact behave differently?",
    "color": "#c98a2c",
    "what_does": (
        "Splits tracks into 'contact' vs 'no_contact' groups based on "
        "whether each track has at least one sufficiently long unbroken run "
        "of contact timepoints, then compares them (optionally split "
        "further by other columns)."
    ),
    "concept": {
        "term": "Contiguous bout",
        "text": (
            "a track counts as 'contact' only if it has an unbroken run of "
            "consecutive contact timepoints at least as long as the minimum "
            "you set — brief, isolated contacts don't count."
        ),
    },
    "what_get": "A contact-vs-no-contact comparison PDF.",
    "decide": [
        _decide("Which column marks contact", _ONLY_YOU),
        _decide("Minimum contiguous contact bout length", _ESTIMATED),
        _decide("Additional column(s) to split by", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain the BEHAV3D Contact-based Grouping analysis: how the "
        "contact/no-contact split works and how to choose the minimum bout "
        "length."
    ),
}

TRACK_EXEMPLAR_TRACKS = {
    "id": "track_exemplars",
    "title": "Exemplar Tracks",
    "subtitle": "What does a typical track from each cluster actually look like?",
    "color": "#c98a2c",
    "what_does": (
        "Selects and renders representative example tracks per cluster, so "
        "you can sanity-check what each movement-type cluster looks like in "
        "practice."
    ),
    "concept": None,
    "what_get": (
        "A statebar overview PDF and/or backprojection PDFs/MP4s of "
        "exemplar tracks, depending on which outputs you enable."
    ),
    "decide": [
        _decide("Number of exemplars per cluster", _DEFAULTS),
        _decide("Which output(s) to generate (statebars / PDFs / MP4)", _SUGGESTED),
    ],
    "has_params": True,
    "seed": (
        "Explain BEHAV3D Exemplar Tracks: how exemplar tracks are chosen "
        "per cluster and what the statebar/backprojection outputs show."
    ),
}

# Level-1 pipeline lists, in display order (matches existing widget order).
STATE_REPORT_PIPELINES = [
    STATE_COMPOSITION_REPORT, STATE_TRANSITION_REPORT, STATE_COMPARISON_REPORT,
]
TRACK_PLOT_PIPELINES = [
    TRACK_DIAGNOSTICS, TRACK_PROPORTIONS, TRACK_COMPARISON_REPORT,
    TRACK_CONTACT_GROUPING, TRACK_EXEMPLAR_TRACKS,
]


# Grouped per sub-tab, in display order.
DEATH_DYNAMICS_ANALYSES = [DEATH_DYNAMICS, INTERACTION, INVASIVENESS]
SINGLE_CELL_ANALYSES = [BEHAVIORAL_STATE, STATE_TRAJECTORY]


# Short intro shown at the top of each Guided page.
GUIDED_INTRO = (
    "New to the analysis? Open \"See what this does\" under any analysis for a "
    "plain-language description. Guided mode explains and suggests — it never "
    "runs anything on its own; you confirm every choice in the settings."
)
