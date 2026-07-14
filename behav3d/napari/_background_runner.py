"""Shared background-execution helpers for the BEHAV3D napari plugin.

Mirrors the pattern established by
:mod:`behav3d.napari._segment_editor` (see ``_EditWorker`` /
``_run_operation_async``) and makes it reusable from the Segmentation,
Tracking, Feature Extraction, and Filtering tabs.

Three concerns are packaged here:

1. :class:`_AsyncWorker` – a generic ``QObject`` that executes an arbitrary
   callable in a ``QThread`` and exposes ``done`` / ``failed`` / ``finished``
   / ``progress`` Qt signals.  The callable is invoked with an optional
   ``progress_cb`` keyword which, when called by the backend, drives both
   the in-tab progress widget and the napari activity-dock entry.

2. :class:`BackgroundOperation` – per-owner-widget orchestrator.  Wires
   the worker up to a :class:`ProgressBarRow`, manages the napari
   activity-dock entry, disables/restores a set of buttons while running,
   and routes ``done`` / ``failed`` back to the Qt thread.

3. :class:`ProgressBarRow` – the tiny ``QProgressBar`` + status ``QLabel``
   row each tab embeds once and reuses across operations.  Has explicit
   ``start`` / ``set_busy`` / ``update`` / ``finish`` helpers so calling
   code stays terse.

A :class:`ThreadSafeLogger` helper is included so backend functions that
accept a ``log_callback=`` can keep writing to the GUI ``QTextEdit`` from
the worker thread without touching Qt widgets directly.

This module is imported only by the Qt-based napari plugin tabs.  It is
intentionally **not** referenced by ``behav3d/widgets/*`` (the notebook
ipywidgets layer) so notebook behaviour is unaffected.
"""
from __future__ import annotations

import time
import traceback
from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence

from qtpy.QtCore import QObject, QThread, QTimer, Qt, Signal
from qtpy.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QWidget,
)

try:
    from napari.utils import progress as _NapariProgress
except Exception:  # pragma: no cover - napari < 0.5 or import error fallback
    _NapariProgress = None  # type: ignore[assignment]


__all__ = [
    "ProgressBarRow",
    "BackgroundOperation",
    "ThreadSafeLogger",
    "make_activity_progress",
    "update_activity_progress",
    "close_activity_progress",
    "show_napari_activity_panel",
    "hide_napari_activity_panel",
    "fire_extra_callback",
]


def fire_extra_callback(extra_callbacks, key: str, *args) -> None:
    """Best-effort fire of ``extra_callbacks[key](*args)``.

    The processing queue plumbs ``extra_callbacks={"on_done": cb,
    "on_failed": cb, "progress": cb}`` through every Run-batch method.
    Each method invokes this helper at the same point as its existing
    on_done / on_failed wrappers.  Default ``None`` keeps GUI / notebook
    behaviour untouched; bad callbacks never propagate.
    """
    if not extra_callbacks:
        return
    cb = extra_callbacks.get(key)
    if cb is None:
        return
    try:
        cb(*args)
    except Exception:
        traceback.print_exc()


# ---------------------------------------------------------------------------
# Tab-side progress widget
# ---------------------------------------------------------------------------
class ProgressBarRow(QWidget):
    """Tiny progress bar + status label embedded in a tab.

    The widget is hidden by default; call :meth:`start` to make it visible
    in busy/indeterminate mode, :meth:`update` to switch to determinate
    mode, and :meth:`finish` to hide it again.
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        self.lbl = QLabel("")
        self.lbl.setStyleSheet("color:#555; font-size:11px; min-width:140px;")
        self.bar = QProgressBar()
        self.bar.setRange(0, 0)  # indeterminate by default
        self.bar.setMaximumHeight(14)
        self.bar.setTextVisible(True)
        lay.addWidget(self.lbl)
        lay.addWidget(self.bar, stretch=1)
        self.setVisible(False)
        # Wall-clock when the current run started (set in start/set_busy).
        # Used to derive an ETA shown on the label during update().
        self._start_time: Optional[float] = None

    # -- public API ----------------------------------------------------
    def start(self, desc: str = "Running…") -> None:
        """Show in busy/indeterminate mode."""
        self._start_time = time.time()
        self.lbl.setText(desc)
        self.bar.setRange(0, 0)
        self.bar.setValue(0)
        self.setVisible(True)

    def set_busy(self, desc: Optional[str] = None) -> None:
        """Switch back to indeterminate (busy) mode."""
        if self._start_time is None:
            self._start_time = time.time()
        if desc is not None:
            self.lbl.setText(desc)
        self.bar.setRange(0, 0)
        self.bar.setValue(0)
        self.setVisible(True)

    def update(self, current: int, total: int, label: str = "") -> None:
        """Switch to determinate mode and update the bar.

        If ``total <= 0`` the bar stays indeterminate (only the label is
        updated).  The label always gets a ``" — NN% — ETA H:MM:SS"`` suffix
        when both ``current > 0`` and ``total > 0`` so the user sees live
        percentage and estimated time even when the backend only emits a
        bare description.
        """
        if total > 0:
            if self.bar.maximum() != total:
                self.bar.setRange(0, total)
            self.bar.setValue(max(0, min(current, total)))
        else:
            self.bar.setRange(0, 0)

        base_label = label if label else self.lbl.text().split(" — ")[0]
        suffix = self._format_progress_suffix(current, total)
        full_label = f"{base_label}{suffix}" if suffix else base_label
        if full_label:
            self.lbl.setText(full_label)
        self.setVisible(True)

    def finish(self) -> None:
        """Hide the row."""
        self.setVisible(False)
        self.lbl.setText("")
        self.bar.setRange(0, 0)
        self.bar.setValue(0)
        self._start_time = None

    # -- helpers -------------------------------------------------------
    def _format_progress_suffix(self, current: int, total: int) -> str:
        """Return ``" — NN% — ETA H:MM:SS"`` (or just ``" — NN%"``)."""
        if total <= 0 or current < 0:
            return ""
        pct = int(round(100.0 * current / total))
        pct = max(0, min(pct, 100))
        if self._start_time is None or current <= 0:
            return f" — {pct}%"
        elapsed = max(0.0, time.time() - self._start_time)
        if elapsed <= 0:
            return f" — {pct}%"
        rate = current / elapsed
        if rate <= 0:
            return f" — {pct}%"
        remaining_s = max(0.0, (total - current) / rate)
        return f" — {pct}% — ETA {self._format_eta(remaining_s)}"

    @staticmethod
    def _format_eta(seconds: float) -> str:
        """Format seconds as ``H:MM:SS`` (or ``M:SS`` when under one hour)."""
        s = int(round(seconds))
        h, rem = divmod(s, 3600)
        m, sec = divmod(rem, 60)
        if h > 0:
            return f"{h}:{m:02d}:{sec:02d}"
        return f"{m}:{sec:02d}"


# ---------------------------------------------------------------------------
# Thread-safe log adapter
# ---------------------------------------------------------------------------
class ThreadSafeLogger(QObject):
    """Bridge a worker-thread ``log_callback`` to a Qt-thread sink.

    Backend functions that accept ``log_callback=`` may emit log lines from
    the worker thread; :class:`~qtpy.QtWidgets.QTextEdit` is not
    thread-safe so we must marshal those calls back to the Qt thread via a
    signal/slot.

    The instance itself is callable, mirroring the original
    ``log_callback(msg)`` interface so it is a drop-in replacement:

    .. code-block:: python

        logger = ThreadSafeLogger(self.log)
        run_pixel_classifier_segmentation(..., log_callback=logger)
    """

    log = Signal(str)

    def __init__(self, sink: Callable[[str], None]) -> None:
        super().__init__()
        self._sink = sink
        # Queue-connect so calls cross-thread; if the sink is on the Qt
        # thread the slot fires there.
        self.log.connect(self._dispatch, Qt.QueuedConnection)

    def _dispatch(self, msg: str) -> None:
        try:
            self._sink(msg)
        except Exception:
            # Never let a broken log handler take down a worker thread.
            traceback.print_exc()

    def __call__(self, msg: str) -> None:
        try:
            self.log.emit(str(msg))
        except Exception:
            # Fallback: write directly if the signal failed for some reason.
            try:
                self._sink(str(msg))
            except Exception:
                traceback.print_exc()


# ---------------------------------------------------------------------------
# napari activity-dock helpers (lifted from _segment_editor.py)
# ---------------------------------------------------------------------------
def make_activity_progress(viewer, desc: str = "Running…"):
    """Create an indeterminate napari activity-dock progress bar.

    ``total=0`` initially renders as a spinning busy indicator until the
    real total is set on the first progress tick via
    :func:`update_activity_progress`.

    Returns ``None`` when napari's ``progress`` helper is unavailable, so
    callers can no-op gracefully.
    """
    if _NapariProgress is None:
        return None
    try:
        pbr = _NapariProgress(total=0, desc=desc)
        if viewer is not None:
            show_napari_activity_panel(viewer)
        return pbr
    except Exception:
        return None


def show_napari_activity_panel(viewer) -> None:
    """Best-effort attempt to make the napari activity panel visible.

    The internal attribute name for the activity dock varies across napari
    versions; we try each in turn and fall back to a generic dock-search.
    """
    try:
        win = viewer.window._qt_window
    except Exception:
        return
    for attr in (
        "_activity_dialog",
        "_qt_activity",
        "_activity_dock",
    ):
        obj = getattr(win, attr, None)
        if obj is not None and hasattr(obj, "show"):
            try:
                obj.show()
                obj.raise_()
            except Exception:
                pass
            return
    try:
        from qtpy.QtWidgets import QDialog, QDockWidget
        for child in win.findChildren((QDockWidget, QDialog)):
            if "activity" in (child.objectName() or "").lower():
                child.show()
                child.raise_()
                return
    except Exception:
        pass


def hide_napari_activity_panel(viewer) -> None:
    """Best-effort attempt to hide/collapse the napari activity panel."""
    try:
        win = viewer.window._qt_window
    except Exception:
        return
    for attr in (
        "_activity_dialog",
        "_qt_activity",
        "_activity_dock",
    ):
        obj = getattr(win, attr, None)
        if obj is not None and hasattr(obj, "hide"):
            try:
                obj.hide()
            except Exception:
                pass
            return
    try:
        from qtpy.QtWidgets import QDialog, QDockWidget
        for child in win.findChildren((QDockWidget, QDialog)):
            if "activity" in (child.objectName() or "").lower():
                child.hide()
                return
    except Exception:
        pass


def update_activity_progress(pbr, current: int, total: int) -> None:
    """Update an existing activity-dock progress entry.

    First call re-initialises the bar with the correct ``total`` via
    ``reset()`` so the Qt widget's maximum is set properly and the
    percentage / ETA is meaningful.  Subsequent calls use ``update(delta)``
    (the canonical tqdm API) which internally refreshes the Qt widget.
    """
    if pbr is None:
        return
    try:
        if total > 0 and pbr.total != total:
            try:
                pbr.reset(total=total)
            except Exception:
                pbr.total = total
                pbr.n = 0
        delta = current - pbr.n
        if delta > 0:
            pbr.update(delta)
        elif delta < 0:
            pbr.n = current
            pbr.refresh()
    except Exception:
        pass


def close_activity_progress(pbr, viewer=None, *, delay_ms: int = 1000) -> None:
    """Fill the bar to 100 %, wait ``delay_ms`` ms, then close the entry."""
    if pbr is None:
        return
    try:
        if pbr.total and pbr.total > 0 and pbr.n < pbr.total:
            pbr.update(pbr.total - pbr.n)
    except Exception:
        pass

    def _do_close():
        try:
            pbr.close()
        except Exception:
            pass
        if viewer is not None:
            hide_napari_activity_panel(viewer)

    QTimer.singleShot(max(0, int(delay_ms)), _do_close)


# ---------------------------------------------------------------------------
# Generic worker
# ---------------------------------------------------------------------------
class _AsyncWorker(QObject):
    """Run ``fn(*args, **kwargs)`` in a background ``QThread``.

    Emits, in order:
        - ``progress(current, total, label)`` zero or more times,
        - either ``done(result)`` or ``failed(error_msg)`` exactly once,
        - ``finished()`` exactly once (no arguments) so the owning
          ``QThread`` can be told to quit via a clean signal connection.

    ``progress_cb`` is injected into ``kwargs`` so callers can either:
        - accept it and emit it themselves, or
        - have the backend ignore it (it defaults to ``None`` everywhere
          else, so the call is safe even if the backend doesn't know
          about it — see :class:`BackgroundOperation` for the inject-only
          path that omits the kwarg).
    """

    done = Signal(object)
    failed = Signal(str)
    finished = Signal()
    progress = Signal(int, int, str)

    def __init__(self, fn: Callable[..., Any], args: tuple, kwargs: dict,
                 inject_progress: bool = True) -> None:
        super().__init__()
        self._fn = fn
        self._args = tuple(args or ())
        self._kwargs = dict(kwargs or {})
        if inject_progress:
            # Bind ``progress_cb`` to the signal's emit method so the
            # backend can call it freely from the worker thread.  Qt's
            # auto-connection to the Qt thread is handled at the
            # ``progress.connect(...)`` site in :class:`BackgroundOperation`.
            self._kwargs.setdefault("progress_cb", self._emit_progress)

    def _emit_progress(self, current: int, total: int, label: str = "") -> None:
        try:
            self.progress.emit(int(current), int(total), str(label))
        except Exception:
            # Never let a progress emit kill the worker.
            pass

    def run(self) -> None:
        try:
            result = self._fn(*self._args, **self._kwargs)
            self.done.emit(result)
        except Exception as exc:
            traceback.print_exc()
            self.failed.emit(str(exc))
        finally:
            self.finished.emit()


# ---------------------------------------------------------------------------
# Background-operation orchestrator
# ---------------------------------------------------------------------------
@dataclass
class _RunState:
    thread: QThread
    worker: _AsyncWorker
    on_done: Optional[Callable[[Any], None]]
    on_failed: Optional[Callable[[str], None]]
    activity_progress: Any
    buttons: List[QWidget]
    progress_row: Optional[ProgressBarRow]
    finish_delay_ms: int
    viewer: Any = None


class BackgroundOperation(QObject):
    """Run a backend callable in the background with progress feedback.

    One instance is intended to be owned by a single widget (tab or sub-
    panel) and re-used across runs.  Only one operation may be active at a
    time on a given instance; calling :meth:`run` again while a previous
    operation is in flight raises ``RuntimeError`` so callers get a clear
    signal rather than silently dropping the request.

    Typical usage from a Qt button slot::

        self._bg.run(
            fn=self._run_tracking_for,
            args=(self.cell_type,),
            kwargs={"overwrite": True},
            desc=f"Tracking {self.cell_type}…",
            progress_row=self.tab.progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=self._on_tracking_done,
            on_failed=self._on_tracking_failed,
        )

    ``fn`` is invoked with ``progress_cb`` automatically injected (unless
    ``inject_progress=False``).  ``progress_cb(current, total, label)``
    drives both the in-tab progress row and the napari activity-dock
    entry.

    External observers (e.g. the processing queue) can subscribe to the
    public :attr:`progress` Qt signal to receive the same
    ``(current, total, label)`` tuples without owning the underlying
    worker.
    """

    progress = Signal(int, int, str)

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._state: Optional[_RunState] = None

    # ------------------------------------------------------------------
    def is_running(self) -> bool:
        st = self._state
        return st is not None and st.thread is not None

    # ------------------------------------------------------------------
    def run(
        self,
        fn: Callable[..., Any],
        args: Sequence[Any] = (),
        kwargs: Optional[Dict[str, Any]] = None,
        *,
        desc: str = "Running…",
        progress_row: Optional[ProgressBarRow] = None,
        buttons: Optional[Iterable[QWidget]] = None,
        viewer=None,
        on_done: Optional[Callable[[Any], None]] = None,
        on_failed: Optional[Callable[[str], None]] = None,
        inject_progress: bool = True,
        indeterminate: bool = False,
        finish_delay_ms: int = 600,
    ) -> None:
        """Kick off a background run.

        Parameters
        ----------
        fn:
            Callable executed on the worker thread.  It receives
            ``*args, **kwargs``; when ``inject_progress=True`` (default) a
            ``progress_cb`` kwarg is added so the backend can report
            progress.
        desc:
            Short label shown in the in-tab progress row and the napari
            activity-dock entry.
        progress_row:
            Optional :class:`ProgressBarRow` to drive.  Tabs typically
            instantiate one of these once and pass it here.
        buttons:
            Iterable of ``QWidget`` to disable while the operation runs.
            They are re-enabled on completion / failure (including via
            crash).
        viewer:
            ``napari.Viewer`` used to surface / hide the activity dock.
        on_done / on_failed:
            Callbacks invoked on the Qt thread when the operation
            finishes.  ``on_done`` receives the value returned by ``fn``.
        inject_progress:
            When ``False``, ``progress_cb`` is not added to ``kwargs``
            (use this for backends that take other progress keywords or
            none at all).
        indeterminate:
            When ``True``, force the in-tab progress row into busy mode
            even when the worker emits progress (useful for backends that
            don't expose a per-sample loop).
        finish_delay_ms:
            How long the activity-dock bar lingers at 100 % before
            disappearing.
        """
        if self.is_running():
            print("BackgroundOperation already running, ignoring request.", flush=True)
            return

        kwargs = dict(kwargs or {})
        btn_list = [b for b in (buttons or []) if b is not None]

        # Disable buttons up-front on the Qt thread.
        for b in btn_list:
            try:
                b.setEnabled(False)
            except Exception:
                pass

        # Tab-side progress row.
        if progress_row is not None:
            if indeterminate:
                progress_row.set_busy(desc)
            else:
                progress_row.start(desc)

        # napari activity-dock entry.
        activity_progress = make_activity_progress(viewer, desc=desc)

        # Worker + thread.
        worker = _AsyncWorker(fn, tuple(args), kwargs,
                              inject_progress=inject_progress)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(thread.quit)

        state = _RunState(
            thread=thread,
            worker=worker,
            on_done=on_done,
            on_failed=on_failed,
            activity_progress=activity_progress,
            buttons=btn_list,
            progress_row=progress_row,
            finish_delay_ms=finish_delay_ms,
            viewer=viewer,
        )
        self._state = state

        # Progress routing — both the tab row and the activity dock.
        if not indeterminate:
            worker.progress.connect(self._on_progress, Qt.QueuedConnection)
        worker.done.connect(self._on_done, Qt.QueuedConnection)
        worker.failed.connect(self._on_failed, Qt.QueuedConnection)
        thread.finished.connect(self._on_thread_finished, Qt.QueuedConnection)

        thread.start()

    # ------------------------------------------------------------------
    def _on_progress(self, current: int, total: int, label: str) -> None:
        st = self._state
        if st is None:
            return
        # Tab row
        if st.progress_row is not None:
            try:
                st.progress_row.update(int(current), int(total), str(label))
            except Exception:
                pass
        # Activity dock
        try:
            update_activity_progress(st.activity_progress, int(current), int(total))
        except Exception:
            pass
        # Re-emit on the public signal so external observers (e.g. the
        # processing queue) receive the same tuple.
        try:
            self.progress.emit(int(current), int(total), str(label))
        except Exception:
            pass

    # ------------------------------------------------------------------
    def _on_done(self, result: Any) -> None:
        st = self._state
        if st is None:
            return
        cb = st.on_done
        # Close UI feedback before calling user callback so the callback's
        # own dialogs / dom updates don't overlap with the activity bar
        # being torn down.
        self._cleanup_ui()
        if cb is not None:
            try:
                cb(result)
            except Exception:
                traceback.print_exc()

    # ------------------------------------------------------------------
    def _on_failed(self, error: str) -> None:
        st = self._state
        if st is None:
            return
        cb = st.on_failed
        self._cleanup_ui()
        if cb is not None:
            try:
                cb(error)
            except Exception:
                traceback.print_exc()

    # ------------------------------------------------------------------
    def _cleanup_ui(self) -> None:
        st = self._state
        if st is None:
            return
        # Re-enable buttons.
        for b in st.buttons:
            try:
                b.setEnabled(True)
            except Exception:
                pass
        # Finish the in-tab progress row.
        if st.progress_row is not None:
            try:
                st.progress_row.finish()
            except Exception:
                pass
        # Close the activity-dock entry AND hide the dock itself so the
        # panel collapses after the run finishes.
        try:
            close_activity_progress(
                st.activity_progress,
                viewer=st.viewer,
                delay_ms=st.finish_delay_ms,
            )
        except Exception:
            pass
        st.activity_progress = None

    # ------------------------------------------------------------------
    def _on_thread_finished(self) -> None:
        """Called once the worker's QThread has fully stopped.

        Safe place to drop the Python references; clearing earlier (inside
        ``_on_done``/``_on_failed``) would release the last reference to
        ``QThread`` while the OS thread is still alive — a recipe for
        ``QThread: Destroyed while thread is still running`` and the
        resulting hard crash.
        """
        self._state = None

    # ------------------------------------------------------------------
    def wait(self, timeout_ms: int = -1) -> bool:
        """Block the current thread until the operation completes.

        Returns ``True`` if the worker thread terminated within the
        timeout, ``False`` otherwise.  Useful only in tests / shutdown
        paths — interactive code should rely on the ``on_done`` callback.
        """
        st = self._state
        if st is None or st.thread is None:
            return True
        try:
            return bool(st.thread.wait(timeout_ms if timeout_ms >= 0 else 0x7FFFFFFF))
        except Exception:
            return False
