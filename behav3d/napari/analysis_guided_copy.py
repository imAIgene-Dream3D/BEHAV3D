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


# Grouped per sub-tab, in display order.
DEATH_DYNAMICS_ANALYSES = [DEATH_DYNAMICS, INTERACTION, INVASIVENESS]
SINGLE_CELL_ANALYSES = [BEHAVIORAL_STATE, STATE_TRAJECTORY]


# Short intro shown at the top of each Guided page.
GUIDED_INTRO = (
    "New to the analysis? Open \"See what this does\" under any analysis for a "
    "plain-language description. Guided mode explains and suggests — it never "
    "runs anything on its own; you confirm every choice in the settings."
)
