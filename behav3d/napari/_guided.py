"""Reusable widgets for the Analysis tab's Guided on-ramp.

This module is deliberately free of any analysis logic — it only renders the
explanation-first "Guided" overview that each sub-tab shows by default.

Pieces
------
- :class:`AnalysisExplainer` — one analysis block: real name + plain-language
  subtitle, a collapsed explainer behind a tinted "See what this does" cue, an
  honest "what you'll decide" list, an "Ask the assistant" button, and an
  always-visible "Start" button.
- :class:`GuidedPanel`     — the scrollable list of explainers for one sub-tab.
- :func:`make_back_header` — the "‹ Back to overview" + title bar shown above
  an isolated analysis's settings once Start has been pressed.

Copy lives in :mod:`behav3d.napari.analysis_guided_copy`; the two sub-tab
modules own the Start wiring and the QStackedWidget that holds the Guided
overview (page 0) and, once Start is pressed, that one analysis's isolated
settings (page 1+), with a Back button to return to page 0.
"""
from __future__ import annotations

from typing import Callable

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QToolButton,
    QFrame, QScrollArea, QSizePolicy,
)
from qtpy.QtCore import Qt


# ── Assistant hand-off ──────────────────────────────────────────────────────
def find_main_widget(w):
    """Walk the Qt parent chain to the BEHAV3DWidget (holds ``assistant``)."""
    node = w
    while node is not None:
        if hasattr(node, "assistant") and hasattr(node, "_assistant_dock"):
            return node
        node = node.parent() if hasattr(node, "parent") else None
    return None


def ask_assistant(origin: QWidget, prompt: str) -> None:
    """Reveal the assistant dock and send ``prompt`` on the user's behalf."""
    main = find_main_widget(origin)
    if main is None or getattr(main, "assistant", None) is None:
        return
    try:
        main.ask_assistant(prompt)
    except Exception:
        pass


# ── Tag colouring for the "what you'll decide" list ─────────────────────────
_TAG_COLORS = {
    "required": "#ff8a80",    # only you know this
    "suggested": "#90caf9",   # suggested set, adjustable
    "estimated": "#90caf9",   # can be estimated
    "default": "#888888",     # tested defaults, optional
}


# ── One analysis explainer ──────────────────────────────────────────────────
class AnalysisExplainer(QWidget):
    """A single analysis block for the Guided page.

    ``spec`` is one of the dicts in :mod:`analysis_guided_copy`. ``on_start`` is
    called (no args) when the user presses Start.
    """

    def __init__(self, spec: dict, on_start: Callable[[], None], parent=None):
        super().__init__(parent)
        self._spec = spec
        self._on_start = on_start
        self._cue = None
        self._body = None
        show_explainer = spec.get("show_explainer", True)

        card = QFrame(self)
        card.setObjectName("guidedCard")
        card.setStyleSheet(
            "#guidedCard { background:#2a2a2e; border:1px solid #3a3a40; "
            "border-radius:10px; }"
        )
        wrap = QVBoxLayout(self)
        wrap.setContentsMargins(0, 0, 0, 6)
        wrap.addWidget(card)

        lay = QVBoxLayout(card)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(8)

        # Header: colored marker + real name, plain-language subtitle underneath.
        title = QLabel(
            f"<span style='color:{spec['color']}; font-size:15px;'>●</span>"
            f"&nbsp;<span style='font-weight:bold; font-size:14px; color:#eaeaea;'>"
            f"{spec['title']}</span>"
        )
        title.setTextFormat(Qt.RichText)
        lay.addWidget(title)
        subtitle = QLabel(spec["subtitle"])
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color:#9aa0a6; font-size:12px; margin-left:16px;")
        lay.addWidget(subtitle)

        if show_explainer:
            # Tinted "See what this does" cue that toggles the explainer body.
            self._cue = QToolButton()
            self._cue.setCheckable(True)
            self._cue.setChecked(False)
            self._cue.setCursor(Qt.PointingHandCursor)
            self._cue.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
            self._cue.setArrowType(Qt.RightArrow)
            self._cue.setText("ⓘ  New to this? See what this analysis does")
            self._cue.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            self._cue.setStyleSheet(
                "QToolButton { background:#1f3a5f; color:#bbdefb; border:1px solid "
                "#3f6ea5; border-radius:6px; padding:6px 10px; text-align:left; "
                "font-size:12px; font-weight:bold; }"
                "QToolButton:hover { background:#274a78; }"
            )
            self._cue.clicked.connect(self._toggle_body)
            lay.addWidget(self._cue)

            # Explainer body (hidden by default).
            self._body = QFrame()
            self._body.setVisible(False)
            body = QVBoxLayout(self._body)
            body.setContentsMargins(4, 4, 4, 2)
            body.setSpacing(8)

            what = QLabel(spec["what_does"])
            what.setWordWrap(True)
            what.setStyleSheet("color:#dcdcdc; font-size:12px;")
            body.addWidget(what)

            if spec.get("concept"):
                c = spec["concept"]
                callout = QLabel(
                    f"<span style='color:#ffd54f;'>💡</span> "
                    f"<b style='color:#eaeaea;'>{c['term']}</b> "
                    f"<span style='color:#c5cae9;'>{c['text']}</span>"
                )
                callout.setTextFormat(Qt.RichText)
                callout.setWordWrap(True)
                callout.setStyleSheet(
                    "background:#20304a; border-radius:6px; padding:8px 10px; "
                    "font-size:12px;"
                )
                body.addWidget(callout)

            get = QLabel(f"<b>What you'll get:</b>&nbsp;{spec['what_get']}")
            get.setTextFormat(Qt.RichText)
            get.setWordWrap(True)
            get.setStyleSheet("color:#c8c8c8; font-size:12px;")
            body.addWidget(get)

            decide_hdr = QLabel("WHAT YOU'LL DECIDE")
            decide_hdr.setStyleSheet(
                "color:#8a8f98; font-size:10px; font-weight:bold; letter-spacing:1px;"
            )
            body.addWidget(decide_hdr)
            hint = QLabel(
                "These depend on your biology — no setting can pick them perfectly. "
                "The assistant can go deeper on any of them."
            )
            hint.setWordWrap(True)
            hint.setStyleSheet("color:#7f858d; font-size:11px;")
            body.addWidget(hint)

            for item in spec.get("decide", []):
                color = _TAG_COLORS.get(item.get("kind", "default"), "#888888")
                row = QLabel(
                    f"<span style='color:#cfcfcf; font-size:12px;'>• {item['label']}</span> "
                    f"<span style='color:{color}; font-size:11px;'>{item['tag']}</span>"
                )
                row.setTextFormat(Qt.RichText)
                row.setWordWrap(True)
                body.addWidget(row)

            ask = QPushButton("💬  Ask the assistant")
            ask.setCursor(Qt.PointingHandCursor)
            ask.setStyleSheet(
                "QPushButton { background:transparent; color:#bbdefb; border:1px "
                "solid #3f6ea5; border-radius:6px; padding:5px 12px; font-size:12px; }"
                "QPushButton:hover { background:#20304a; }"
            )
            ask.clicked.connect(self._on_ask)
            ask_row = QHBoxLayout()
            ask_row.addWidget(ask)
            ask_row.addStretch()
            body.addLayout(ask_row)

            lay.addWidget(self._body)

        # Always-visible primary action.
        start = QPushButton(spec.get("start_label", "Start — open the settings  ▸"))
        start.setCursor(Qt.PointingHandCursor)
        start.setStyleSheet(
            "QPushButton { background:#1e88e5; color:white; font-weight:bold; "
            "border:none; border-radius:6px; padding:7px 16px; font-size:13px; }"
            "QPushButton:hover { background:#1976d2; }"
        )
        start.clicked.connect(lambda: self._on_start())
        start_row = QHBoxLayout()
        start_row.addWidget(start)
        start_row.addStretch()
        lay.addLayout(start_row)

    def _toggle_body(self):
        show = self._cue.isChecked()
        self._body.setVisible(show)
        self._cue.setArrowType(Qt.DownArrow if show else Qt.RightArrow)

    def _on_ask(self):
        ask_assistant(self, self._spec.get("seed", ""))


# ── Guided page for one sub-tab ─────────────────────────────────────────────
class GuidedPanel(QScrollArea):
    """Scrollable list of :class:`AnalysisExplainer` for one sub-tab.

    ``analyses`` is a list of copy specs; ``start_cb`` is called with the
    analysis ``id`` when its Start button is pressed.
    """

    def __init__(self, analyses: list[dict], start_cb: Callable[[str], None],
                 intro: str = "", parent=None):
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setFrameShape(QFrame.NoFrame)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        if intro:
            intro_lbl = QLabel(intro)
            intro_lbl.setWordWrap(True)
            intro_lbl.setStyleSheet(
                "color:#9aa0a6; font-size:11px; padding:2px 4px 6px 4px;"
            )
            lay.addWidget(intro_lbl)

        for spec in analyses:
            aid = spec["id"]
            lay.addWidget(
                AnalysisExplainer(spec, on_start=lambda a=aid: start_cb(a))
            )

        lay.addStretch()
        self.setWidget(content)


# ── Back header for an isolated analysis view ───────────────────────────────
def make_nav_button(label: str, on_click: Callable[[], None]) -> QPushButton:
    """A blue nav button ("‹ Back", "‹ Back to overview") shared by every
    isolated-view header."""
    btn = QPushButton(label)
    btn.setCursor(Qt.PointingHandCursor)
    btn.setStyleSheet(
        "QPushButton { background:#1e88e5; color:white; border:1px solid #1565c0; "
        "border-radius:6px; padding:5px 12px; font-size:12px; font-weight:bold; }"
        "QPushButton:hover { background:#1976d2; border-color:#1565c0; color:white; }"
        "QPushButton:pressed { background:#1565c0; }"
    )
    btn.clicked.connect(on_click)
    return btn


def make_back_header(on_back: Callable[[], None], label: str = "‹ Back to overview") -> "tuple[QWidget, QLabel]":
    """Return a (bar, title_label) pair: a Back button plus a title label the
    caller updates per-analysis (e.g. via ``setText``)."""
    bar = QWidget()
    lay = QHBoxLayout(bar)
    lay.setContentsMargins(4, 4, 4, 4)
    lay.setSpacing(10)

    lay.addWidget(make_nav_button(label, on_back))

    title = QLabel("")
    title.setStyleSheet("color:#eaeaea; font-size:13px; font-weight:bold;")
    lay.addWidget(title)
    lay.addStretch()

    return bar, title
