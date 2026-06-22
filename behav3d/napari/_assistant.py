"""
BEHAV3D assistant — the right-side co-pilot dock widget.

``AssistantDock`` is added as a *separate* napari dock (right area) and holds a
reference to the main ``BEHAV3DWidget`` so it can read live workflow context
(:func:`_assistant_context.build_context`) and apply proposed parameter changes
(:mod:`_assistant_actions`) after the user confirms.

Layout (top → bottom):
  • context bar          — 📍 step · N samples · output dir · M queued (live)
  • chat transcript      — QTextBrowser, markdown, streamed assistant tokens
  • proposed-action tray — confirm cards: "old -> new  [Fill it in] [Dismiss]"
  • quick actions        — Explain this screen · Tell me about my data
  • input row            — text box + Send

All network I/O happens on a QThread worker (``ChatWorker``); the GUI never blocks.
"""
from __future__ import annotations

from typing import Optional

from qtpy.QtCore import Qt, QThread, QTimer, Signal
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit,
    QTextBrowser, QFrame, QSizePolicy,
)

from behav3d.napari._assistant_context import build_context, context_summary_line
from behav3d.napari._assistant_actions import (
    build_actions, apply_action, TOOL_SCHEMA, ProposedAction, humanize_parameter_key,
)
from behav3d.napari._assistant_client import ChatWorker, load_client_config


_UNSET = object()  # sentinel for "parameter not yet in config"
_METADATA_FIELD_DEFAULTS = {
    "n_samples": 1,
    "n_organoids": 0,
    "n_immune": 0,
    "n_other": 0,
    "include_dead": False,
    "exp_nr": 1,
    "pixel_distance_xy": 0.5,
    "pixel_distance_z": 2.0,
    "time_interval": 1.0,
    "time_unit": "s",
}


def _coerce_bool(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "on")
    return bool(value)


def _metadata_field_has_value(dp, field: str, index: int, cell_type: str | None = None) -> bool:
    """Return True if the metadata builder field already has a non-default value set."""
    if dp is None:
        return False
    try:
        from qtpy.QtWidgets import QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox
        if field == "n_samples":
            return dp.n_samples_spin.value() != _METADATA_FIELD_DEFAULTS["n_samples"]
        if field == "n_organoids":
            return dp.n_organoid_spin.value() != _METADATA_FIELD_DEFAULTS["n_organoids"]
        if field == "n_immune":
            return dp.n_immune_spin.value() != _METADATA_FIELD_DEFAULTS["n_immune"]
        if field == "n_other":
            return dp.n_other_spin.value() != _METADATA_FIELD_DEFAULTS["n_other"]
        if field == "include_dead":
            return dp.include_dead_cb.isChecked() != _METADATA_FIELD_DEFAULTS["include_dead"]
        if field in ("immune_multicolor", "immune_multicolor_channels"):
            flags = getattr(dp, "_immune_multicolor_flags", [])
            counts = getattr(dp, "_immune_multicolor_counts", [])
            if field == "immune_multicolor":
                return index < len(flags) and flags[index].isChecked()
            return index < len(counts) and counts[index].value() != 2
        for attr, cat in (("_organoid_name_edits", "organoid_name"),
                          ("_immune_name_edits",   "immune_name"),
                          ("_other_name_edits",    "other_name")):
            if field == cat:
                edits = getattr(dp, attr, [])
                return index < len(edits) and bool(edits[index].text().strip())
        # Per-sample fields
        forms = getattr(dp, "_sample_forms", [])
        if index >= len(forms):
            return False
        if field in ("dead_channel_number", "dead_mask_path"):
            dead_key = "number" if field == "dead_channel_number" else "mask_path"
            w = forms[index].get("dead_channel", {}).get(dead_key)
        elif field in (
            "cell_line", "cell_condition", "cell_segments_image_path",
            "cell_tracks_image_path", "cell_tracks_csv_path",
        ):
            field_map = {
                "cell_line": "line",
                "cell_condition": "condition",
                "cell_segments_image_path": "segments_image_path",
                "cell_tracks_image_path": "tracks_image_path",
                "cell_tracks_csv_path": "tracks_csv_path",
            }
            cell_fields = forms[index].get("cell_types", {}).get(str(cell_type), {})
            w = cell_fields.get(field_map[field])
        else:
            w = forms[index]["basic"].get(field)
        if w is None:
            return False
        if isinstance(w, QLineEdit):
            return bool(w.text().strip())
        if isinstance(w, (QSpinBox, QDoubleSpinBox)):
            val = w.value()
            if field in _METADATA_FIELD_DEFAULTS:
                return val != _METADATA_FIELD_DEFAULTS[field]
            return val != 0 and val != w.minimum()
        if isinstance(w, QComboBox):
            default = _METADATA_FIELD_DEFAULTS.get(field)
            if default is not None:
                return w.currentText() != str(default)
            return w.currentIndex() > 0
    except Exception:
        pass
    return False


def _fill_value_matches_current(
    dp,
    field: str,
    index: int,
    proposed,
    cell_type: str | None = None,
) -> bool:
    """Return True if the proposed fill value is the same as what's already in the widget.
    Used to auto-apply no-op repeated calls without showing a confirm card."""
    if dp is None:
        return False
    try:
        from qtpy.QtWidgets import QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox
        if field == "n_samples":
            return dp.n_samples_spin.value() == int(proposed)
        if field == "n_organoids":
            return dp.n_organoid_spin.value() == int(proposed)
        if field == "n_immune":
            return dp.n_immune_spin.value() == int(proposed)
        if field == "n_other":
            return dp.n_other_spin.value() == int(proposed)
        if field == "include_dead":
            return dp.include_dead_cb.isChecked() == _coerce_bool(proposed)
        if field in ("immune_multicolor", "immune_multicolor_channels"):
            flags = getattr(dp, "_immune_multicolor_flags", [])
            counts = getattr(dp, "_immune_multicolor_counts", [])
            if field == "immune_multicolor":
                return index < len(flags) and flags[index].isChecked() == _coerce_bool(proposed)
            return index < len(counts) and counts[index].value() == int(proposed)
        for attr, cat in (("_organoid_name_edits", "organoid_name"),
                          ("_immune_name_edits",   "immune_name"),
                          ("_other_name_edits",    "other_name")):
            if field == cat:
                edits = getattr(dp, attr, [])
                return (index < len(edits) and
                        edits[index].text().strip() == str(proposed).strip())
        forms = getattr(dp, "_sample_forms", [])
        if index >= len(forms):
            return False
        if field in ("dead_channel_number", "dead_mask_path"):
            dead_key = "number" if field == "dead_channel_number" else "mask_path"
            w = forms[index].get("dead_channel", {}).get(dead_key)
        elif field in (
            "cell_line", "cell_condition", "cell_segments_image_path",
            "cell_tracks_image_path", "cell_tracks_csv_path",
        ):
            field_map = {
                "cell_line": "line",
                "cell_condition": "condition",
                "cell_segments_image_path": "segments_image_path",
                "cell_tracks_image_path": "tracks_image_path",
                "cell_tracks_csv_path": "tracks_csv_path",
            }
            cell_fields = forms[index].get("cell_types", {}).get(str(cell_type), {})
            w = cell_fields.get(field_map[field])
        else:
            w = forms[index]["basic"].get(field)
        if w is None:
            return False
        if isinstance(w, QLineEdit):
            return w.text().strip() == str(proposed).strip()
        if isinstance(w, (QSpinBox, QDoubleSpinBox)):
            return w.value() == float(proposed)
        if isinstance(w, QComboBox):
            return w.currentText() == str(proposed)
    except Exception:
        pass
    return False


_SYSTEM_PRIMER = (
    "You are the BEHAV3D co-pilot, embedded in the napari plugin for analysing "
    "cell behaviour in 3D fluorescent imaging. Help the user choose methods and "
    "parameter values for their data. Ground answers in the provided step_schema "
    "and retrieved docs. When you recommend a concrete value, also emit a "
    "set_parameter tool call (the user confirms before it applies).\n\n"
    "FORMATTING — keep responses easy to read in a narrow side panel:\n"
    "- Keep it short: lead with a one-sentence answer, then detail only if useful.\n"
    "- Use short paragraphs (1–2 sentences) separated by a blank line.\n"
    "- Use bullet points or a numbered list for any set of items or steps.\n"
    "- Use researcher-facing labels, not internal variable names.\n"
    "- Avoid long dense paragraphs and walls of text."
)


class _ActionCard(QFrame):
    """A confirm card for a single proposed action."""

    confirmed = Signal(object)   # ProposedAction
    dismissed = Signal(object)

    def __init__(self, action: ProposedAction, parent=None):
        super().__init__(parent)
        self.action = action
        self.setStyleSheet(
            "QFrame { background:#20303a; border:1px solid #3a5a6a; border-radius:4px; }"
        )
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(4)

        title = QLabel(action.preview or action.kind)
        title.setWordWrap(True)
        title.setStyleSheet("color:#e0f0ff; font-size:11px;")
        lay.addWidget(title)

        if not action.ok:
            warn = QLabel(action.message)
            warn.setWordWrap(True)
            warn.setStyleSheet("color:#ffb3b3; font-size:10px;")
            lay.addWidget(warn)

        row = QHBoxLayout()
        row.addStretch(1)
        if action.ok:
            btn_apply = QPushButton("Fill it in")
            btn_apply.setStyleSheet("font-size:10px;")
            btn_apply.clicked.connect(lambda: self.confirmed.emit(self.action))
            row.addWidget(btn_apply)
        btn_dismiss = QPushButton("Dismiss")
        btn_dismiss.setStyleSheet("font-size:10px;")
        btn_dismiss.clicked.connect(lambda: self.dismissed.emit(self.action))
        row.addWidget(btn_dismiss)
        lay.addLayout(row)


class AssistantDock(QWidget):
    """Right-side co-pilot chat dock for BEHAV3D."""

    def __init__(self, main_widget, parent=None):
        super().__init__(parent)
        self.main_widget = main_widget
        self.setMinimumWidth(300)
        self._history: list[dict] = []        # [{role, content}] sent to the model
        self._threads: list[tuple] = []       # keep (thread, worker) refs alive
        self._md_log: list[str] = []          # finalised transcript blocks (markdown)
        self._streaming_text = None           # in-progress assistant text, or None
        self._greeted_tabs: set = set()       # tabs greeted this session (no repeat proactive)
        self._guided_flow_active: bool = False  # True while a step-by-step guide is running
        self._pending_auto_continue: list | None = None  # deferred after _on_finished
        self._auto_continue_turns: int = 0
        self._last_auto_continue_signature: tuple | None = None
        # Progress-based auto-continue guard: keep going as long as each turn applies a
        # NEW distinct value; only stop after consecutive stalled turns (no new progress).
        self._auto_continue_seen_signatures: set[tuple] = set()
        self._auto_continue_stall_count: int = 0
        # The metadata-builder field applied on the current turn — drives the
        # deterministic next-question logic in _auto_continue.
        self._last_applied_md_field: str | None = None
        # Deterministic metadata wizard state (no LLM): _md_current is the question
        # awaiting an answer; _md_queue holds the remaining questions in this phase.
        self._md_current: dict | None = None
        self._md_queue: list[dict] = []
        self._md_phase: str | None = None
        self._confirmed_parameter_keys: set[str] = set()
        # Coalesce streamed-token re-renders. Re-rendering the whole transcript on
        # every token saturates the GUI thread on long answers (the napari window
        # freezes until the stream ends). This single-shot timer batches tokens into
        # at most one render per interval; the expensive full-document style walk is
        # applied only on the final render of the turn.
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(90)
        self._render_timer.timeout.connect(self._render_streaming)

        # Cohesive dark-theme styling for the whole dock (matches napari surfaces).
        self.setStyleSheet("""
            AssistantDock { background-color: #2a2d30; }
            QTextBrowser {
                background-color: #232629;
                color: #e6e6e6;
                border: 1px solid #3a3f44;
                border-radius: 8px;
                padding: 12px;
                font-family: -apple-system, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
                font-size: 13px;
                selection-background-color: #3d6b8a;
            }
            QLineEdit {
                background-color: #2f3338; color: #e6e6e6;
                border: 1px solid #3a3f44; border-radius: 8px;
                padding: 7px 10px; font-size: 13px;
            }
            QLineEdit:focus { border: 1px solid #5285a6; }
            QPushButton {
                background-color: #353a3f; color: #e6e6e6;
                border: 1px solid #454b50; border-radius: 8px;
                padding: 7px 10px; font-size: 12px;
            }
            QPushButton:hover { background-color: #3f464c; }
            QPushButton:pressed { background-color: #2c3035; }
        """)

        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # --- context bar --------------------------------------------------
        self.context_bar = QLabel("BEHAV3D Assistant")
        self.context_bar.setStyleSheet(
            "background:#2f3338; color:#9fc6e0; padding:7px 10px;"
            "border:1px solid #3a3f44; border-radius:8px; font-size:11px; font-weight:600;"
        )
        self.context_bar.setWordWrap(True)
        root.addWidget(self.context_bar)

        # --- transcript ---------------------------------------------------
        self.transcript = QTextBrowser()
        self.transcript.setOpenExternalLinks(True)
        self.transcript.setFrameShape(QFrame.NoFrame)
        self.transcript.document().setDocumentMargin(2)
        root.addWidget(self.transcript, stretch=1)

        # --- proposed-action tray ----------------------------------------
        self.action_tray = QWidget()
        self.action_tray_layout = QVBoxLayout(self.action_tray)
        self.action_tray_layout.setContentsMargins(0, 0, 0, 0)
        self.action_tray_layout.setSpacing(4)
        root.addWidget(self.action_tray)

        # --- quick actions (rebuilt on tab switch) ------------------------
        self.quick_widget = QWidget()
        self._quick_layout = QVBoxLayout(self.quick_widget)
        self._quick_layout.setContentsMargins(0, 0, 0, 0)
        self._quick_layout.setSpacing(4)
        root.addWidget(self.quick_widget)

        # --- input row ----------------------------------------------------
        input_row = QHBoxLayout()
        self.input = QLineEdit()
        self.input.setPlaceholderText("Ask about parameters, methods, your data…")
        self.input.returnPressed.connect(self._on_send)
        input_row.addWidget(self.input, stretch=1)
        self.btn_send = QPushButton("Send")
        self.btn_send.clicked.connect(self._on_send)
        input_row.addWidget(self.btn_send)
        root.addLayout(input_row)

        self.refresh_context_bar()
        self._greet()
        try:
            initial_idx = self.main_widget.tabs.currentIndex()
            self.main_widget.tabs.currentChanged.connect(self._on_tab_changed)
        except Exception:
            initial_idx = 0
        self._greeted_tabs.add(initial_idx)
        self._set_quick_buttons(initial_idx)

    # ------------------------------------------------------------------
    # Context bar
    # ------------------------------------------------------------------
    def refresh_context_bar(self):
        try:
            ctx = build_context(self.main_widget)
            self.context_bar.setText(context_summary_line(ctx))
        except Exception:
            self.context_bar.setText("BEHAV3D Assistant")

    def _greet(self):
        cfg = load_client_config()
        mode = "" if cfg.get("endpoint") else "  *(offline mode — no server configured)*"
        self._append_md(
            f"**BEHAV3D Assistant**{mode}\n\nAsk me which method or parameter "
            "values suit your data. I can also fill the forms in for you "
            "(you confirm first). Try *“Explain this screen”* to start."
        )

    # ------------------------------------------------------------------
    # Transcript helpers
    # ------------------------------------------------------------------
    def _render(self, style: bool = True):
        """Re-render the whole transcript as markdown.

        QTextCursor.insertMarkdown does not exist in PyQt5, so we keep a log of
        markdown blocks and render the document with QTextEdit.setMarkdown (which
        does render bold/italic/lists/code), then pin the scrollbar to the bottom.

        ``style=False`` skips the full-document style walk — used for the throttled
        renders during token streaming so a long answer doesn't freeze the GUI.
        """
        # Any direct/full render supersedes a pending throttled one.
        self._render_timer.stop()
        blocks = list(self._md_log)
        if self._streaming_text is not None:
            blocks.append(f"**BEHAV3D Assistant**\n\n{self._streaming_text or '...'}")
        try:
            # A thin rule between turns gives clear visual separation.
            self.transcript.setMarkdown("\n\n---\n\n".join(blocks))
            if style:
                self._style_blocks()
        except Exception:
            # Fallback for any Qt build lacking setMarkdown.
            self.transcript.setPlainText("\n\n".join(blocks))
        sb = self.transcript.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _render_streaming(self):
        """Throttled render during streaming — skips the expensive style walk,
        which is applied once when the turn finishes (_on_finished → _render)."""
        self._render(style=False)

    def _schedule_render(self):
        """Request a render soon, coalescing bursts of streamed tokens into at
        most one render per timer interval."""
        if not self._render_timer.isActive():
            self._render_timer.start()

    def _style_blocks(self):
        """setMarkdown renders paragraphs/lists very tightly. Walk the document
        and add line spacing + paragraph spacing so responses are readable."""
        try:
            from qtpy.QtGui import QTextCursor, QTextBlockFormat
        except Exception:
            return
        doc = self.transcript.document()
        block = doc.begin()
        while block.isValid():
            cur = QTextCursor(block)
            bf = cur.blockFormat()
            bf.setLineHeight(140, QTextBlockFormat.ProportionalHeight)  # 1.4× line spacing
            bf.setTopMargin(2.0)
            bf.setBottomMargin(9.0)                                     # gap between paragraphs
            cur.setBlockFormat(bf)
            block = block.next()

    def _append_md(self, markdown: str):
        self._md_log.append(markdown)
        self._render()

    def _append_user(self, text: str):
        self._append_md(f"**You**\n\n{text}")

    # ------------------------------------------------------------------
    # Sending
    # ------------------------------------------------------------------
    def _on_send(self):
        text = self.input.text().strip()
        if not text:
            return
        self.input.clear()
        # Deterministic metadata wizard owns the input while it's running (the
        # questions and answers never touch the LLM, so it cannot stall or loop).
        if self._md_current is not None:
            self._md_handle_answer(text)
            return
        self._guided_flow_active = False   # free-form input exits any guided flow
        self._send_message(text)

    def _send_message(
        self,
        text: str,
        *,
        show_as_user: bool = True,
        display_text: str | None = None,
        reset_auto_loop: bool = True,
    ):
        if reset_auto_loop:
            self._auto_continue_turns = 0
            self._last_auto_continue_signature = None
            self._auto_continue_seen_signatures = set()
            self._auto_continue_stall_count = 0
            self._last_applied_md_field = None
            # Any LLM-bound send (free-form or another guide button) cancels the
            # deterministic metadata wizard.
            self._md_current = None
            self._md_queue = []
            self._md_phase = None
        if show_as_user:
            self._append_user(display_text or text)
        self._history.append({"role": "user", "content": text})
        self._set_busy(True)
        self._streaming_text = ""        # opens a live "Assistant:" block
        self._render()

        ctx = {}
        try:
            ctx = build_context(self.main_widget)
            ctx["assistant_session"] = {
                "confirmed_parameter_keys": sorted(self._confirmed_parameter_keys),
            }
        except Exception:
            ctx = {}

        messages = [{"role": "system", "content": _SYSTEM_PRIMER}] + self._history
        worker = ChatWorker(messages=messages, context=ctx, tools=TOOL_SCHEMA)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.token.connect(self._on_token)
        worker.tool_calls.connect(self._on_tool_calls)
        worker.degraded.connect(self._on_degraded)
        worker.error.connect(self._on_error)
        worker.finished.connect(thread.quit)
        worker.finished.connect(self._on_finished)
        thread.finished.connect(lambda: self._cleanup_thread(thread, worker))
        self._threads.append((thread, worker))
        thread.start()

    def _set_busy(self, busy: bool):
        self.btn_send.setEnabled(not busy)
        # Leave the input box editable while a turn is in flight so the user can
        # type their next answer without waiting for the stream to finish.
        for i in range(self._quick_layout.count()):
            row = self._quick_layout.itemAt(i)
            sub = row.layout() if row else None
            if sub:
                for j in range(sub.count()):
                    w = sub.itemAt(j).widget() if sub.itemAt(j) else None
                    if w:
                        w.setEnabled(not busy)

    # ------------------------------------------------------------------
    # Worker signal handlers
    # ------------------------------------------------------------------
    def _on_token(self, chunk: str):
        if self._streaming_text is None:
            self._streaming_text = ""
        self._streaming_text += chunk
        self._schedule_render()

    def _on_degraded(self, full_text: str):
        self._streaming_text = full_text
        self._render()

    def _on_error(self, message: str):
        # close any open streaming block, then show the error
        self._finalize_streaming()
        self._append_md(message)
        # Defensive: ensure the dock never stays stuck "thinking" if signal
        # ordering ever changes (worker.finished still clears busy too).
        self._set_busy(False)

    def _on_tool_calls(self, calls: list):
        try:
            ctx = build_context(self.main_widget)
            cards = ctx.get("step_schema", [])
            # include all cards so cross-step suggestions still validate
            from behav3d.napari._assistant_schema import flatten_config_to_cards
            cards = flatten_config_to_cards()
            params = getattr(self.main_widget.data_prep_tab, "behav3d_parameters", {})
            actions = build_actions(calls, cards, params)
        except Exception:
            actions = []
        auto_previews = []
        noop_previews = []
        had_auto = False
        # Suppress per-widget pulse glows when applying a batch — a flood of fills
        # would otherwise schedule hundreds of QTimers and thrash stylesheets,
        # which can freeze the GUI. One context-bar refresh happens at the end.
        from behav3d.napari._assistant_actions import set_pulse_suppressed
        batch = len(actions) > 1
        if batch:
            set_pulse_suppressed(True)
        try:
            for act in actions:
                if self._is_noop_action(act):
                    self._record_noop_action(act)
                    noop_previews.append(self._noop_preview(act))
                    continue
                if self._should_auto_apply(act):
                    if act.ok:
                        # Classify BEFORE applying (widget still holds pre-apply state).
                        _silent_auto = self._is_silent_auto_apply(act)
                        try:
                            ok = apply_action(self.main_widget, act)
                        except Exception:
                            ok = False
                        if ok:
                            had_auto = True
                            if act.kind == "bulk_fill_metadata":
                                n = (act.data.get("sample_count")
                                     or len(act.data.get("samples", []) or []))
                                self._append_md(
                                    f"Filled the Metadata Builder for {n} "
                                    f"sample{'s' if n != 1 else ''}. Review the values "
                                    "and tell me anything you'd like to change."
                                )
                                # Bulk fill completes the form in one pass — do NOT add
                                # to auto_previews, so it doesn't trigger an auto-continue.
                            # Only count as a meaningful change if it wasn't a structural
                            # trigger or a same-value no-op re-call from the model.
                            elif not _silent_auto:
                                auto_previews.append(act.preview)
                                # Remember the metadata field just filled so the
                                # deterministic counts driver can ask the next one.
                                if act.kind == "fill_metadata_builder":
                                    self._last_applied_md_field = act.data.get("field")
                            if not batch:
                                self.refresh_context_bar()
                else:
                    self._add_action_card(act)
        finally:
            if batch:
                set_pulse_suppressed(False)
                self.refresh_context_bar()

        if auto_previews:
            self._append_md("Filled in: " + " · ".join(auto_previews))
        # A no-op means the proposed value was already present. Treat that as
        # confirmation from the user's previous answer and move forward silently
        # instead of rendering a scary "Already set" card.
        if (
            noop_previews and not auto_previews and self.action_tray_layout.count() == 0
            and not self._text_asks_question(self._streaming_text)
        ):
            self._pending_auto_continue = noop_previews
        # Defer auto-continue until _on_finished so ChatWorker1's _finalize_streaming()
        # runs first — otherwise it wipes _streaming_text that ChatWorker2 just started.
        # Only fire for value-setting actions (non-empty previews). Structural-only actions
        # (open_builder, configure_cell_types, create_sample_forms, fill_down) must not
        # trigger a new LLM turn — the bot should ask the next question in its current
        # streaming response, otherwise the loop never terminates.
        if (
            had_auto and auto_previews and self.action_tray_layout.count() == 0
            and not self._text_asks_question(self._streaming_text)
        ):
            self._pending_auto_continue = auto_previews

    @staticmethod
    def _text_asks_question(text: str | None) -> bool:
        """True if the streamed assistant text already poses a follow-up question.

        The model usually streams a one-line acknowledgement ("Got it — 22 samples.")
        alongside its tool call; that must NOT suppress auto-continue, otherwise the
        guided flow stalls. Only an actual question (a '?' at/near the end) means the
        bot has already prompted the user, so we should not auto-continue."""
        t = (text or "").strip()
        if not t:
            return False
        tail = t.rstrip(" *_`)\n\t")
        if tail.endswith("?"):
            return True
        # A '?' within the last stretch of text also counts (e.g. trailing emoji/note).
        return "?" in t[-120:]

    def _is_noop_action(self, action: ProposedAction) -> bool:
        return bool(action.data.get("no_op"))

    def _record_noop_action(self, action: ProposedAction) -> None:
        if action.kind == "set_parameter":
            key = action.data.get("key")
            if key:
                self._confirmed_parameter_keys.add(str(key))

    def _noop_preview(self, action: ProposedAction) -> str:
        if action.kind == "set_parameter":
            label = action.data.get("label") or humanize_parameter_key(action.data.get("key"))
            return f"{label} already set to {action.data.get('value')!r}"
        return action.preview or "Value already set"

    def _should_auto_apply(self, action: ProposedAction) -> bool:
        """Auto-apply by default; only show a confirm card if a value is already set."""
        if not action.ok:
            return False

        if action.kind == "navigate_to_step":
            return True  # navigation never overwrites anything

        if action.kind == "add_queue_step":
            return False  # adding to the queue is a deliberate action

        if action.kind == "bulk_fill_metadata":
            return True  # deterministic bulk fill applies in one guarded pass

        if action.kind == "select_segmentation_method":
            # Apply only when it changes the current selection.
            seg = getattr(self.main_widget, "segmentation_tab", None)
            combo = getattr(seg, "method_combo", None) if seg is not None else None
            if combo is None:
                return False
            cur = combo.currentText()
            val = str(action.data.get("value", "")).strip()
            return not (cur == val or cur.startswith(val) or (val and val in cur))

        if action.kind == "fill_metadata_builder":
            field = action.data.get("field")
            # Structural triggers never overwrite anything.
            if field in ("open_builder", "configure_cell_types",
                         "create_sample_forms", "fill_down"):
                return True
            dp = getattr(self.main_widget, "data_prep_tab", None)
            idx = int(action.data.get("index", 0))
            cell_type = action.data.get("cell_type")
            if not _metadata_field_has_value(dp, field, idx, cell_type):
                return True   # blank field — auto-apply
            # Field already has a value; auto-apply only if the proposed value is the
            # same (bot redundantly repeating itself). Different value → confirm card.
            return _fill_value_matches_current(
                dp, field, idx, action.data.get("value"), cell_type
            )

        if action.kind == "set_parameter":
            try:
                from behav3d.napari._assistant_actions import get_by_dotted
                from behav3d.napari._assistant_schema import flatten_config_to_cards
                key = action.data.get("key")
                params = getattr(
                    self.main_widget.data_prep_tab, "behav3d_parameters", {})
                current = get_by_dotted(params, key, _UNSET)
                if current is _UNSET:
                    return True  # never been set
                for card in flatten_config_to_cards():
                    if card["key"] == key:
                        return current == card["default"]
            except Exception:
                pass
            return True

        return True  # unknown action kind: auto-apply

    def _is_silent_auto_apply(self, action: ProposedAction) -> bool:
        """True when the action auto-applies but should NOT count as a new user answer.

        Covers two cases:
        - Structural triggers (open_builder, configure_cell_types, etc.) — no user data.
        - Same-value no-ops — model re-confirming a value that is already set.
          Detected BEFORE applying, while the widget still holds its current value.
        """
        if action.kind != "fill_metadata_builder":
            return False
        field = action.data.get("field")
        if field in ("open_builder", "configure_cell_types",
                     "create_sample_forms", "fill_down"):
            return True
        # Same-value no-op: proposed equals current, checked pre-apply.
        # This covers both "field at a meaningful value" and "field still at default"
        # (e.g. n_samples=1 → 1==1) — both are no-ops that should not trigger auto-continue.
        dp = getattr(self.main_widget, "data_prep_tab", None)
        idx = int(action.data.get("index", 0))
        return _fill_value_matches_current(
            dp, field, idx, action.data.get("value"), action.data.get("cell_type")
        )

    def _finalize_streaming(self):
        """Move the in-progress assistant text into the finalised transcript log."""
        if self._streaming_text:
            self._md_log.append(f"**BEHAV3D Assistant**\n\n{self._streaming_text}")
            self._history.append({"role": "assistant", "content": self._streaming_text})
        self._streaming_text = None

    def _on_finished(self):
        had_text = bool(self._streaming_text)
        self._finalize_streaming()
        self._render()
        self._set_busy(False)
        if had_text:
            self._append_log()
        # Fire any deferred auto-continue AFTER this turn is fully finalised.
        pending = self._pending_auto_continue
        if pending is not None:
            self._pending_auto_continue = None
            self._auto_continue(pending)

    def _cleanup_thread(self, thread, worker):
        self._threads = [(t, w) for (t, w) in self._threads if t is not thread]
        worker.deleteLater()
        thread.deleteLater()

    # ------------------------------------------------------------------
    # Proposed-action tray
    # ------------------------------------------------------------------
    def _add_action_card(self, action: ProposedAction):
        card = _ActionCard(action)
        card.confirmed.connect(self._apply_action)
        card.dismissed.connect(lambda a, c=card: self._remove_card(c))
        self.action_tray_layout.addWidget(card)

    def _remove_card(self, card: _ActionCard):
        card.setParent(None)
        card.deleteLater()

    def _apply_action(self, action: ProposedAction):
        ok = False
        try:
            ok = apply_action(self.main_widget, action)
        except Exception as e:
            self._append_md(f"\nCould not apply: {e}")
        if ok:
            if action.kind == "set_parameter" and not action.data.get("widget_updated"):
                # stored in config, but no matching field exists on the current screen
                label = action.data.get("label") or humanize_parameter_key(action.data.get("key"))
                self._append_md(
                    f"Saved to config: **{label}** = "
                    f"`{action.data.get('value')!r}`\n\n*This step has no form field "
                    "here yet — it'll be used when you reach that part of the pipeline.*"
                )
            else:
                self._append_md(f"Filled in: {action.preview}")
            self.refresh_context_bar()
        else:
            self._append_md(f"Could not apply: {action.preview}")
        # remove all cards for this action
        for i in reversed(range(self.action_tray_layout.count())):
            w = self.action_tray_layout.itemAt(i).widget()
            if isinstance(w, _ActionCard) and w.action is action:
                self._remove_card(w)
        if ok and action.kind == "fill_metadata_builder":
            self._last_applied_md_field = action.data.get("field")
        # Auto-continue: once the tray is empty, prompt the bot for the next step.
        if ok and self.action_tray_layout.count() == 0:
            self._auto_continue([action.preview] if action.preview else [])

    def _auto_continue(self, applied: list[str] | None = None):
        """Silently prompt the bot to continue, naming what was just applied.

        Continues as long as each turn makes NEW progress (a distinct applied value);
        only stops after two consecutive stalled turns (a repeated or empty signature),
        so legitimate multi-field guided flows run to completion instead of cutting off
        after a fixed number of turns."""
        signature = tuple(applied or [])
        made_progress = bool(signature) and signature not in self._auto_continue_seen_signatures
        if made_progress:
            self._auto_continue_seen_signatures.add(signature)
            self._auto_continue_stall_count = 0
        else:
            self._auto_continue_stall_count += 1
        # Stop only when we've stalled (no NEW distinct value) twice in a row.
        if self._auto_continue_stall_count >= 2:
            self._append_md(
                "I paused because I wasn't making new progress. Tell me the next value "
                "you want to set, or use one of the buttons below."
            )
            return
        # Hard safety ceiling against a runaway loop in pathological cases.
        if self._auto_continue_turns >= 40:
            return
        # --- Deterministic metadata guidance -------------------------------
        # The model reliably FILLS an answer and ASKS a given question, but is
        # unreliable at DECIDING the next field on its own (it loops on the field
        # it just set). So during the metadata guided flow we compute the next
        # question and tell it to ask exactly that. Verified against the live model:
        # explicit next-question = 6/6 pass; generic "ask the next field" = 1/6.
        md_field = self._last_applied_md_field
        self._last_applied_md_field = None
        if self._guided_flow_active and md_field in self._MD_COUNT_ORDER:
            i = self._MD_COUNT_ORDER.index(md_field)
            nxt = self._MD_COUNT_ORDER[i + 1] if i + 1 < len(self._MD_COUNT_ORDER) else None
            if nxt is None:
                # All four counts captured → build the cell-type fields and the
                # per-sample forms (structural + idempotent), then hand off.
                self._md_apply_structural("configure_cell_types")
                self._md_apply_structural("create_sample_forms")
                self.refresh_context_bar()
                self._append_md(
                    "All set — the cell-type fields and per-sample forms are ready below. "
                    "For each sample, fill in the **image path** and acquisition settings; "
                    "use **Fill All from Sample 1** to copy shared values (pixel size, time "
                    "interval), rename the cell types if you like, then **Save Metadata CSV** "
                    "and click **Load Metadata**. Ask me about any field you're unsure of."
                )
                self._guided_flow_active = False
                return
            self._auto_continue_turns += 1
            self._last_auto_continue_signature = signature
            label = "; ".join(applied) if applied else "That"
            self._send_message(
                f"{label} is set — do NOT call any tool now. Ask me exactly this next "
                f'question: "{self._MD_COUNT_Q[nxt]}"',
                show_as_user=False, reset_auto_loop=False,
            )
            return

        # --- Generic fallback (set_parameter flows, etc.) ------------------
        self._auto_continue_turns += 1
        self._last_auto_continue_signature = signature
        if applied:
            msg = (
                f"Applied: {'; '.join(applied)}. Those values are now set — do NOT set "
                "them again. Ask me the next field that is still missing in the builder; "
                "only call a tool after I answer."
            )
        else:
            msg = ("Ask me the next field that is still missing; only call a tool after "
                   "I answer.")
        self._send_message(msg, show_as_user=False, reset_auto_loop=False)

    # Canonical order + questions for the deterministic metadata counts phase.
    _MD_COUNT_ORDER = ["n_samples", "n_organoids", "n_immune", "n_other"]
    _MD_COUNT_Q = {
        "n_organoids": "How many organoid types do you have?",
        "n_immune": "How many immune cell types do you have?",
        "n_other": "How many other (non-organoid, non-immune) cell types do you have?",
    }

    def _md_apply_structural(self, field: str):
        """Apply a structural metadata-builder step (configure_cell_types /
        create_sample_forms) directly, without an LLM round-trip."""
        try:
            act = ProposedAction("fill_metadata_builder", field=field, value=None, index=0)
            apply_action(self.main_widget, act)
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Quick actions
    # ------------------------------------------------------------------
    def _explain_screen(self):
        ctx = {}
        try:
            ctx = build_context(self.main_widget)
        except Exception:
            pass
        step = ctx.get("current_tab_label") or ctx.get("current_step", "this step")
        if ctx.get("current_step") == "visualization":
            # The Visualization tab has no tunable parameters — describe only the
            # controls that are actually on screen, not config-only fields.
            self._send_message(
                "Explain the Visualization tab. It has no tunable parameters — only "
                "these on-screen controls: the Dataset section (Sample selector, "
                "'Clear existing layers before loading', and 'Load Dataset into Napari'), "
                "the Visibility Toggles (Raw, Segments, Tracked Segments, Tracks), and "
                "the Manual Edition section (pick tracked segments and 'Edit tracked "
                "segments'). Describe what each of these controls does — do not mention "
                "any other parameters.",
                display_text="Explain this screen",
            )
            return
        self._send_message(
            f"Explain the {step} tab: what does it do, and what do its key "
            "parameters mean for someone configuring it for the first time?",
            display_text="Explain this screen",
        )

    def _start_interview(self):
        self._send_message(
            "I'd like help configuring BEHAV3D for my data. Ask me a short series "
            "of questions about my experiment (imaging modality, cell types, frame "
            "interval, pixel size, expected cell diameter, motility) and then "
            "recommend sensible default parameters for the current step.",
            display_text="Tell me about my data",
        )

    def _start_metadata_guide(self):
        self._append_user("Build metadata")
        self._begin_metadata_wizard()

    # ------------------------------------------------------------------
    # Deterministic metadata wizard (no LLM — cannot stall or loop)
    # ------------------------------------------------------------------
    _MD_COUNT_STEPS = [
        {"kind": "count", "field": "n_samples",
         "question": "How many samples (fields of view) do you have?"},
        {"kind": "count", "field": "n_organoids",
         "question": "How many organoid types do you have?"},
        {"kind": "count", "field": "n_immune",
         "question": "How many immune cell types do you have?"},
        {"kind": "count", "field": "n_other",
         "question": "How many other (non-organoid, non-immune) cell types do you have?"},
    ]

    def _begin_metadata_wizard(self):
        """Drive the structural metadata setup (sample/cell-type counts and names)
        entirely on the frontend. Questions are posted directly and the typed answer
        is applied immediately — the LLM is never in the loop, so it can't stall,
        loop, or skip a field."""
        self._guided_flow_active = True
        self._md_apply_structural("open_builder")
        self._append_md("Opened the Metadata Builder")
        self.refresh_context_bar()
        self._md_phase = "counts"
        self._md_queue = [dict(s) for s in self._MD_COUNT_STEPS]
        self._md_current = None
        self._md_ask_next()

    def _md_ask_next(self):
        """Post the next question in the current phase, or advance phases."""
        if self._md_queue:
            self._md_current = self._md_queue.pop(0)
            self._append_md(f"**BEHAV3D Assistant**\n\n{self._md_current['question']}")
            return
        if self._md_phase == "counts":
            # Counts done → build the cell-type name fields, then ask each name.
            self._md_apply_structural("configure_cell_types")
            self.refresh_context_bar()
            self._md_phase = "names"
            self._md_queue = self._md_build_name_steps()
            self._md_ask_next()
            return
        if self._md_phase == "names":
            # Names done → build the per-sample forms and hand off.
            self._md_apply_structural("create_sample_forms")
            self.refresh_context_bar()
            self._md_finish()
            return
        self._md_finish()

    def _md_build_name_steps(self) -> list[dict]:
        dp = getattr(self.main_widget, "data_prep_tab", None)
        steps: list[dict] = []
        if dp is None:
            return steps
        spec = [
            ("n_organoid_spin", "organoid_name", "organoid", 'e.g. "tumor"'),
            ("n_immune_spin", "immune_name", "immune cell", 'e.g. "CD8"'),
            ("n_other_spin", "other_name", "other cell", ""),
        ]
        for spin_attr, field, label, hint in spec:
            spin = getattr(dp, spin_attr, None)
            n = spin.value() if spin is not None else 0
            for i in range(n):
                eg = f" ({hint})" if hint else ""
                steps.append({"kind": "name", "field": field, "index": i,
                              "question": f"What is the name of {label} type {i + 1}?{eg}"})
        return steps

    def _md_handle_answer(self, text: str):
        """Apply the typed answer to the current wizard question and advance."""
        item = self._md_current
        self._append_user(text)
        if item["kind"] == "count":
            try:
                value = int("".join(ch for ch in text if (ch.isdigit() or ch == "-")))
            except ValueError:
                # Not a number — treat as a free-form question: exit the wizard and
                # hand the message to the assistant instead of forcing a number.
                self._md_current = None
                self._md_queue = []
                self._md_phase = None
                self._send_message(text, show_as_user=False)
                return
        else:  # name
            value = text.strip()
            if not value:
                self._append_md(f"**BEHAV3D Assistant**\n\n{item['question']}")
                return  # re-ask; keep _md_current
        self._md_current = None
        self._md_apply_value(item["field"], value, item.get("index", 0))
        self.refresh_context_bar()
        self._md_ask_next()

    def _md_apply_value(self, field: str, value, index: int = 0):
        from behav3d.napari._assistant_actions import _metadata_builder_preview
        act = ProposedAction("fill_metadata_builder", field=field, value=value, index=index)
        try:
            ok = apply_action(self.main_widget, act)
        except Exception:
            ok = False
        if ok:
            try:
                preview = _metadata_builder_preview(field, value, index)
            except Exception:
                preview = f"{field} = {value}"
            self._append_md(f"Filled in: {preview}")

    def _md_finish(self):
        self._md_current = None
        self._md_queue = []
        self._md_phase = None
        self._guided_flow_active = False
        self._append_md(
            "All set — the cell-type fields and per-sample forms are ready below. "
            "For each sample, fill in the **image path** and acquisition settings; use "
            "**Fill All from Sample 1** to copy shared values (pixel size, time interval), "
            "then **Save Metadata CSV** and click **Load Metadata**. Ask me about any "
            "field you're unsure of."
        )

    # ------------------------------------------------------------------
    # Tab-change awareness and dynamic quick buttons
    # ------------------------------------------------------------------
    def _on_tab_changed(self, index: int):
        """Called whenever the user switches tabs in the main widget."""
        self.refresh_context_bar()
        self._set_quick_buttons(index)
        if index not in self._greeted_tabs and not self._guided_flow_active:
            self._greeted_tabs.add(index)
            self._dispatch_proactive(index)

    def _dispatch_proactive(self, index: int):
        dispatch = {
            0: self._proactive_data_prep,
            2: self._proactive_segmentation,
            3: self._proactive_tracking,
            4: self._proactive_feature_extraction,
        }
        fn = dispatch.get(index)
        if fn:
            fn()

    def _set_quick_buttons(self, tab_index: int):
        """Clear and repopulate quick-action buttons for the given tab."""
        lay = self._quick_layout
        while lay.count():
            item = lay.takeAt(0)
            sub = item.layout()
            if sub:
                while sub.count():
                    si = sub.takeAt(0)
                    w = si.widget()
                    if w:
                        w.setParent(None)
                        w.deleteLater()

        def _btn(label, handler, tip=""):
            b = QPushButton(label)
            b.setStyleSheet("font-size:10px;")
            if tip:
                b.setToolTip(tip)
            b.clicked.connect(handler)
            return b

        def _row(*btns):
            r = QHBoxLayout()
            for b in btns:
                r.addWidget(b)
            lay.addLayout(r)

        if tab_index == 0:   # data_preparation
            _row(_btn("Walk through setup", self._start_setup_guide),
                 _btn("Check what's missing", self._check_data_prereqs))
            _row(_btn("Build metadata", self._start_metadata_guide))
        elif tab_index == 2:  # segmentation
            _row(_btn("Guide segmentation", self._start_segmentation_guide),
                 _btn("Choose a method", self._explain_seg_methods))
        elif tab_index == 3:  # tracking
            _row(_btn("Guide tracking", self._start_tracking_guide),
                 _btn("Which method?", self._explain_tracking_methods))
        elif tab_index == 4:  # feature_extraction
            _row(_btn("Guide setup", self._start_feature_guide),
                 _btn("Check prerequisites", self._check_feature_prereqs))
        else:
            _row(_btn("Explain this screen", self._explain_screen),
                 _btn("Tell me about my data", self._start_interview))

    def _proactive_message(self, prompt: str):
        """Trigger a proactive LLM check; displays as the bot's initiative, not a user message."""
        self._append_md("*Checking your setup...*")
        self._send_message(prompt, show_as_user=False)

    # ------------------------------------------------------------------
    # Proactive per-tab openers (fire once per tab per session)
    # ------------------------------------------------------------------
    def _proactive_data_prep(self):
        self._proactive_message(
            "I am on the Data Preparation tab. "
            "Check the context: is output_dir_set? Is metadata loaded (metadata.loaded)? "
            "If output_dir is not set, mention that first and offer to help. "
            "If output_dir is set but metadata is not loaded, offer to walk through the Metadata Builder. "
            "If both are ready, confirm briefly and suggest moving to the next step. "
            "Keep it to 2–3 sentences maximum."
        )

    def _proactive_segmentation(self):
        self._proactive_message(
            "I just switched to the Segmentation tab. "
            "Check context: is metadata loaded and output dir set? "
            "If not, state which prereq is missing and offer to call navigate_to_step to Data Preparation. "
            "If prerequisites are met, mention the current visible Method dropdown choice and offer "
            "a guided walkthrough. Treat APOC, ConvPaint, Pixel Classifier, Cellpose, and Import "
            "segmentation as distinct methods. "
            "Keep it to 2–3 sentences."
        )

    def _proactive_tracking(self):
        self._proactive_message(
            "I just switched to the Tracking tab. "
            "Check context: is metadata loaded? What cell types are listed in metadata.cell_types? "
            "If metadata is not loaded, redirect with navigate_to_step to Data Preparation. "
            "If loaded, state which cell types you see (immune, organoid, other) and ask if "
            "I'd like a guided walkthrough to configure tracking methods and parameters. "
            "Keep it to 2–3 sentences."
        )

    def _proactive_feature_extraction(self):
        self._proactive_message(
            "I just switched to the Feature Extraction tab. "
            "Check context: is metadata loaded? What cell types are present? "
            "Briefly say what this step does and mention whether key prereqs are satisfied. "
            "Offer to guide me through setup. Keep it to 2–3 sentences."
        )

    # ------------------------------------------------------------------
    # Guided-flow handlers for per-tab quick buttons
    # ------------------------------------------------------------------
    def _start_setup_guide(self):
        self._guided_flow_active = True
        ctx = {}
        try:
            ctx = build_context(self.main_widget)
        except Exception:
            pass
        if ctx.get("output_dir_set") and not ctx.get("metadata", {}).get("loaded"):
            # Output dir is already set → run the deterministic metadata wizard
            # (same robust path as the "Build metadata" button).
            self._append_user("Walk through setup")
            self._begin_metadata_wizard()
        else:
            self._send_message(
                "I need help setting up Data Preparation from scratch. "
                "Check context — is output_dir_set? Is metadata loaded? "
                "If output dir is not set, ask me for the path. "
                "Otherwise open the Metadata Builder and ask me 'How many samples?' — "
                "one question at a time.",
                display_text="Walk through setup",
            )

    def _check_data_prereqs(self):
        self._send_message(
            "Look at the current context (output_dir_set, metadata.loaded, metadata_builder state, "
            "metadata.n_samples, metadata.cell_types) and list exactly what is missing or "
            "incomplete on the Data Preparation tab. If everything is complete, say so and suggest "
            "the next workflow step.",
            display_text="Check what's missing",
        )

    def _start_segmentation_guide(self):
        self._guided_flow_active = True
        self._send_message(
            "I'd like a guided walkthrough of Segmentation setup. "
            "First check that metadata is loaded and output dir is set — if not, redirect me. "
            "Then ask which visible segmentation Method fits my data. The available choices are "
            "APOC (GPU), ConvPaint (DL pixel classifier), Pixel Classifier (Random Forest), "
            "Cellpose (Deep Learning), and Import segmentation. APOC is its own method, not "
            "Pixel Classifier. One question at a time — just ask for method first.",
            display_text="Guide segmentation",
        )

    def _explain_seg_methods(self):
        self._send_message(
            "Explain the visible segmentation methods in BEHAV3D: APOC (GPU), ConvPaint, "
            "Pixel Classifier (Random Forest), Cellpose, and Import segmentation. "
            "Recommend one for each of my cell types when that makes sense. "
            "Read my cell types from the context and be concrete in your recommendation.",
            display_text="Choose a segmentation method",
        )

    def _start_tracking_guide(self):
        self._guided_flow_active = True
        self._send_message(
            "I'd like a guided walkthrough of Tracking setup. "
            "Read my cell types from context and start by listing which types you see. "
            "For each type ask which tracking method to use — state the recommended default "
            "(btrack for immune, propagation for organoid, lap for other) and ask for confirmation. "
            "Then ask the key parameter for that method. One question at a time.",
            display_text="Guide tracking",
        )

    def _explain_tracking_methods(self):
        self._send_message(
            "Looking at my cell types in context, which tracking method do you recommend for each type "
            "and why? Compare btrack, LAP, propagation, and trackpy. Give a concrete recommendation.",
            display_text="Which tracking method?",
        )

    def _start_feature_guide(self):
        self._guided_flow_active = True
        self._send_message(
            "I'd like a guided walkthrough of Feature Extraction setup. "
            "Check context: are metadata loaded and output dir set? What cell types are present? "
            "If prereqs are missing, redirect me. Otherwise ask one question at a time, "
            "starting with: which features do I want to extract for each cell type?",
            display_text="Guide feature extraction",
        )

    def _check_feature_prereqs(self):
        self._send_message(
            "Check whether the prerequisites for Feature Extraction are satisfied: "
            "metadata loaded, output dir set, segmentation masks present, tracking outputs present. "
            "Tell me specifically what is ready and what is missing. "
            "If anything is missing, offer to navigate me to the right tab to fix it.",
            display_text="Check prerequisites",
        )

    # ------------------------------------------------------------------
    # Session replay log
    # ------------------------------------------------------------------
    def _append_log(self):
        import os
        dp = getattr(self.main_widget, "data_prep_tab", None)
        out_dir = getattr(dp, "output_dir", "") if dp else ""
        if not out_dir or not self._history:
            return
        try:
            path = os.path.join(out_dir, "assistant_log.md")
            last = self._history[-1]
            with open(path, "a", encoding="utf-8") as f:
                f.write(f"\n**{last['role']}:** {last['content']}\n")
        except Exception:
            pass
