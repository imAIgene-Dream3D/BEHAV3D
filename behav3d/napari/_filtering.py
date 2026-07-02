"""
BEHAV3D napari plugin – Filtering Tab.

Provides per-cell-type sub-tabs with track-length filters, dead-at-t0 filter,
time-unit toggle, and batch-run capability including summarization.
Mirrors the architecture of _tracking.py / _feature_extraction.py.

Run buttons execute the backend in a background ``QThread`` (via
:class:`behav3d.napari._background_runner.BackgroundOperation`) so napari
stays responsive.  ``filter_tracks`` / ``summarize_track_features`` have
no per-sample loop accessible, so the per-cell-type Run uses an
indeterminate (busy) progress bar.  The batch Run shows determinate
progress at cell-type granularity.  Queue compatibility is preserved by
defaulting ``block=True``.
"""
import sys
import traceback
from pathlib import Path

import yaml
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QTabWidget, QTextEdit, QCheckBox,
    QSpinBox, QGroupBox, QComboBox, QMessageBox, QScrollArea,
    QSplitter,
)
from qtpy.QtCore import Qt, Signal

from behav3d.core.qt_help import make_help_row
from behav3d.napari._results_panel import (
    ResultsPanel,
    notify_results_changed,
)
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    fire_extra_callback,
)


# ═══════════════════════════════════════════════════════════════════════════
# Per-cell-type panel
# ═══════════════════════════════════════════════════════════════════════════
class CellTypeFilterPanel(QWidget):
    """Filtering controls for one cell type."""

    def __init__(
        self,
        cell_type: str,
        category: str,
        metadata_loader,
        all_cell_types: list,
        category_types: list,
        log_callback=None,
        parent=None,
        tab_progress_row=None,
        viewer=None,
    ):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types
        self.log = log_callback or (lambda m: None)
        self.viewer = viewer

        # Background-execution infrastructure.
        self.tab_progress_row = tab_progress_row
        self._bg = BackgroundOperation(self)

        # Detect dead channel
        md = metadata_loader.metadata
        self._has_dead = False
        if md is not None:
            self._has_dead = "dead_channel" in md.columns and md["dead_channel"].notna().any()

        # Read saved config
        params = self.metadata_loader.behav3d_parameters
        cfg = params.get("track_filtering", {}).get(self.cell_type, {}) or {}

        self._init_ui(cfg)

    # ---- UI ---------------------------------------------------------------
    def _init_ui(self, cfg):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # ── Preview ────────────────────────────────────────────────────────
        self.btn_preview_lengths = QPushButton("📊 See Track Length Distributions")
        self.btn_preview_lengths.setStyleSheet(
            "QPushButton { background: #17a2b8; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px; font-size: 12px; } "
            "QPushButton:hover { background: #138496; }"
        )
        self.btn_preview_lengths.clicked.connect(self._on_preview_lengths_clicked)
        layout.addWidget(self.btn_preview_lengths)

        # ── Experiment duration ────────────────────────────────────────
        self.en_exp_duration = QCheckBox("Trim full time series to max timepoints")
        self.en_exp_duration.setChecked(bool(cfg.get("exp_duration_enabled", True)))
        layout.addWidget(self.en_exp_duration)

        self.spin_exp_duration = QSpinBox()
        self.spin_exp_duration.setRange(1, 99999)
        self.spin_exp_duration.setValue(int(cfg.get("exp_duration", 350)))
        self.spin_exp_duration.setMaximumWidth(100)
        dur_form = QFormLayout()
        dur_form.setContentsMargins(20, 0, 0, 0)
        dur_form.addRow("Max timepoints:", make_help_row(
            self.spin_exp_duration,
            "Max Experiment Duration",
            "Cut off all data after this many timepoints/hours.\n\n"
            "Timepoints beyond this are removed."
        ))
        self.dur_widget = QWidget()
        self.dur_widget.setLayout(dur_form)
        layout.addWidget(self.dur_widget)

        def _toggle_dur(state):
            self.dur_widget.setVisible(self.en_exp_duration.isChecked())
        self.en_exp_duration.stateChanged.connect(_toggle_dur)
        _toggle_dur(None)

        # ── Min track length ──────────────────────────────────────────
        self.en_min_length = QCheckBox("Filter tracks shorter than minimal length")
        self.en_min_length.setChecked(bool(cfg.get("min_length_enabled", True)))
        layout.addWidget(self.en_min_length)

        self.spin_min_length = QSpinBox()
        self.spin_min_length.setRange(1, 99999)
        self.spin_min_length.setValue(int(cfg.get("min_track_length", 30)))
        self.spin_min_length.setMaximumWidth(100)
        min_form = QFormLayout()
        min_form.setContentsMargins(20, 0, 0, 0)
        min_form.addRow("Min length:", make_help_row(
            self.spin_min_length,
            "Minimum Track Length",
            "Tracks with fewer timepoints than this are removed.\n\n"
            "Helps eliminate short, unreliable tracks."
        ))
        self.min_widget = QWidget()
        self.min_widget.setLayout(min_form)
        layout.addWidget(self.min_widget)

        def _toggle_min(state):
            self.min_widget.setVisible(self.en_min_length.isChecked())
        self.en_min_length.stateChanged.connect(_toggle_min)
        _toggle_min(None)

        # ── Max track length ──────────────────────────────────────────
        self.en_max_length = QCheckBox("Trim tracks to maximum length")
        self.en_max_length.setChecked(bool(cfg.get("max_length_enabled", True)))
        layout.addWidget(self.en_max_length)

        self.spin_max_length = QSpinBox()
        self.spin_max_length.setRange(1, 99999)
        self.spin_max_length.setValue(int(cfg.get("max_track_length", 30)))
        self.spin_max_length.setMaximumWidth(100)
        max_form = QFormLayout()
        max_form.setContentsMargins(20, 0, 0, 0)
        max_form.addRow("Max length:", make_help_row(
            self.spin_max_length,
            "Maximum Track Length",
            "Tracks longer than this are trimmed to this length.\n\n"
            "Useful for DTW which works best with equal-length tracks."
        ))
        self.max_widget = QWidget()
        self.max_widget.setLayout(max_form)
        layout.addWidget(self.max_widget)

        def _toggle_max(state):
            self.max_widget.setVisible(self.en_max_length.isChecked())
        self.en_max_length.stateChanged.connect(_toggle_max)
        _toggle_max(None)

        # ── Min size at first timepoint filter ───────────────────────────
        self.check_filter_min_size = QCheckBox("Filter by minimal size at first timepoint")
        self.check_filter_min_size.setChecked(bool(cfg.get("filter_min_size_t1", False)))
        layout.addWidget(self.check_filter_min_size)

        self.spin_min_size_t1 = QSpinBox()
        self.spin_min_size_t1.setRange(1, 999999)
        self.spin_min_size_t1.setValue(int(cfg.get("min_size_t1", 1000)))
        self.spin_min_size_t1.setMaximumWidth(100)
        size_form = QFormLayout()
        size_form.setContentsMargins(20, 0, 0, 0)
        size_form.addRow("Min size (px):", self.spin_min_size_t1)
        self.size_widget = QWidget()
        self.size_widget.setLayout(size_form)
        layout.addWidget(self.size_widget)

        def _toggle_size(state):
            self.size_widget.setVisible(self.check_filter_min_size.isChecked())
        self.check_filter_min_size.stateChanged.connect(_toggle_size)
        _toggle_size(None)

        # ── Filter dead at first timepoint ────────────────────────────────
        self.check_filter_dead_t0 = QCheckBox("Filter dead cells at first timepoint")
        self.check_filter_dead_t0.setChecked(bool(cfg.get("filter_t0_dead", False)))
        self.check_filter_dead_t0.setToolTip(
            "When enabled, tracks where the cell is already marked as 'dead'\n"
            "at the first timepoint (relative_time == 1) are removed.\n\n"
            "Requires the 'dead' column to have been created during\n"
            "Feature Extraction with a dead threshold."
        )
        self.check_filter_dead_t0.setVisible(self._has_dead)
        layout.addWidget(self.check_filter_dead_t0)



        # ── Time unit ─────────────────────────────────────────────────
        unit_label = QLabel("Unit for time-based filters:")
        unit_label.setStyleSheet("font-weight: bold; margin-top: 6px;")
        layout.addWidget(unit_label)

        self.combo_time_type = QComboBox()
        self.combo_time_type.addItems(["frames", "hours"])
        self.combo_time_type.setCurrentText(cfg.get("time_type", "frames"))
        self.combo_time_type.setMaximumWidth(120)
        layout.addWidget(self.combo_time_type)

        # ── Apply settings buttons ────────────────────────────────────
        cat_label = self.category.capitalize() + "s" if self.category != "other" else "Other types"
        sync_row = QHBoxLayout()
        btn_apply_cat = QPushButton(f"Apply to all {cat_label}")
        btn_apply_cat.clicked.connect(lambda: self._apply_to_others(category_only=True))
        btn_apply_all = QPushButton("Apply to all")
        btn_apply_all.clicked.connect(lambda: self._apply_to_others(category_only=False))
        sync_row.addWidget(btn_apply_cat)
        sync_row.addWidget(btn_apply_all)
        layout.addLayout(sync_row)

        # ── Run button ────────────────────────────────────────────────
        self.btn_run = QPushButton(f"Filter {self.cell_type.capitalize()} Tracks & Summarize")
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        # ``clicked`` emits a ``bool`` (checked state) that would override the
        # ``interactive=True`` default of ``_on_run_clicked``; wrap in a lambda
        # so the overwrite prompt is not silently skipped on button presses.
        self.btn_run.clicked.connect(lambda: self._on_run_clicked(interactive=True))
        layout.addWidget(self.btn_run)

        layout.addStretch()

    # ---- Helpers ----------------------------------------------------------
    def _on_preview_lengths_clicked(self):
        if self.viewer is not None and self.viewer.layers:
            reply = QMessageBox.question(
                self, "Warning",
                "This will remove current layers in the viewer. Continue?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if reply == QMessageBox.No:
                return

        ct = self.cell_type
        out = Path(self.metadata_loader.output_dir) if self.metadata_loader.output_dir else None
        if not out:
            self.log(f"⚠️ Metadata output dir not set.")
            return
            
        csv_path = out / "analysis" / ct / "track_features" / f"BEHAV3D_{ct}_combined_track_features.csv"
        if not csv_path.exists():
            self.log(f"⚠️ Unfiltered track features not found:\n  {csv_path}")
            return
            
        try:
            import pandas as pd
            import matplotlib.pyplot as plt
            import tempfile
            from qtpy.QtGui import QDesktopServices
            from qtpy.QtCore import QUrl
            import os
            
            df = pd.read_csv(csv_path)
            if "sample_name" not in df.columns or "TrackID" not in df.columns:
                self.log(f"⚠️ CSV missing required columns (sample_name, TrackID).")
                return
                
            samples = df["sample_name"].unique()
            n_samples = len(samples)
            if n_samples == 0:
                self.log(f"⚠️ No samples found in CSV.")
                return
                
            cols = min(3, n_samples)
            rows = (n_samples + cols - 1) // cols
            fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), squeeze=False)
            axes = axes.flatten()
            
            for i, sample in enumerate(samples):
                ax = axes[i]
                sample_df = df[df["sample_name"] == sample]
                track_lengths = sample_df.groupby("TrackID").size()
                
                ax.hist(track_lengths, bins=50, color='skyblue', edgecolor='black')
                ax.set_title(sample)
                ax.set_xlabel("Track Length (timepoints)")
                ax.set_ylabel("Frequency")
                
            for j in range(i + 1, len(axes)):
                fig.delaxes(axes[j])
                
            plt.tight_layout()
            
            fd, temp_path = tempfile.mkstemp(suffix=".png", prefix=f"track_lengths_{ct}_")
            os.close(fd)
            plt.savefig(temp_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            QDesktopServices.openUrl(QUrl.fromLocalFile(temp_path))
            self.log(f"📊 Track length distribution generated.")
            
            # Actually clear layers if the user consented to this preview.
            if self.viewer is not None:
                self.viewer.layers.clear()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.log(f"❌ Error generating plot: {e}")

    def _collect_params(self) -> dict:
        d = {
            "exp_duration_enabled": self.en_exp_duration.isChecked(),
            "exp_duration": int(self.spin_exp_duration.value()),
            "min_length_enabled": self.en_min_length.isChecked(),
            "min_track_length": int(self.spin_min_length.value()),
            "max_length_enabled": self.en_max_length.isChecked(),
            "max_track_length": int(self.spin_max_length.value()),
            "filter_min_size_t1": self.check_filter_min_size.isChecked(),
            "min_size_t1": int(self.spin_min_size_t1.value()),
            "filter_t0_dead": self.check_filter_dead_t0.isChecked(),
            "time_type": self.combo_time_type.currentText(),
        }
        return d

    def _apply_to_others(self, category_only=False):
        parent_tab = self.parent()
        while parent_tab and not hasattr(parent_tab, "panels"):
            parent_tab = parent_tab.parent()
        if not parent_tab:
            return

        settings = self._collect_params()
        targets = self.category_types if category_only else self.all_cell_types
        count = 0
        for ct in targets:
            if ct == self.cell_type:
                continue
            if ct in parent_tab.panels:
                p = parent_tab.panels[ct]
                p.en_exp_duration.setChecked(settings["exp_duration_enabled"])
                p.spin_exp_duration.setValue(settings["exp_duration"])
                p.en_min_length.setChecked(settings["min_length_enabled"])
                p.spin_min_length.setValue(settings["min_track_length"])
                p.en_max_length.setChecked(settings["max_length_enabled"])
                p.spin_max_length.setValue(settings["max_track_length"])
                p.check_filter_min_size.setChecked(settings["filter_min_size_t1"])
                p.spin_min_size_t1.setValue(settings["min_size_t1"])
                p.check_filter_dead_t0.setChecked(settings.get("filter_t0_dead", False))
                p.combo_time_type.setCurrentText(settings["time_type"])
                count += 1
        scope = "category" if category_only else "all"
        self.log(f"Applied filter settings to {count} other cell types ({scope}).")

    def _persist(self):
        params = self.metadata_loader.behav3d_parameters
        filtering = params.setdefault("track_filtering", {})
        filtering[self.cell_type] = self._collect_params()

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    def collect_runtime_params(self) -> dict:
        """Snapshot Qt widget values for thread-safe handoff to a worker."""
        return {
            "exp_duration_enabled": bool(self.en_exp_duration.isChecked()),
            "exp_duration": int(self.spin_exp_duration.value()),
            "min_length_enabled": bool(self.en_min_length.isChecked()),
            "min_track_length": int(self.spin_min_length.value()),
            "max_length_enabled": bool(self.en_max_length.isChecked()),
            "max_track_length": int(self.spin_max_length.value()),
            "filter_min_size": bool(self.check_filter_min_size.isChecked()),
            "min_size_t1": int(self.spin_min_size_t1.value()),
            "filter_t0_dead": bool(self.check_filter_dead_t0.isChecked()),
            "time_type": self.combo_time_type.currentText(),
        }

    def _run_filtering_for(self, cell_type: str, overwrite: bool = False,
                           params: dict = None):
        """Run track filtering + summarization for one cell type.

        Pass ``params`` when calling from a background worker (must be
        a snapshot from :meth:`collect_runtime_params`).  When ``None``,
        widget values are read directly — only safe on the Qt thread.
        """
        from behav3d.analysis.filtering import filter_tracks
        from behav3d.analysis import summarize_track_features

        if params is None:
            params = self.collect_runtime_params()

        print(f"\n{'='*50}", file=sys.stderr)
        print(f"  Filtering: {cell_type}", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        # Check for advanced features CSV
        active_killing_dir = Path(out_dir) / "analysis" / cell_type / "active_killing"
        adv_path = active_killing_dir / f"BEHAV3D_{cell_type}_advanced_track_features.csv"
        df_input_path = str(adv_path) if adv_path.exists() else None

        filter_kwargs = {
            "metadata": self.metadata_loader.metadata,
            "output_dir": out_dir,
            "cell_type": cell_type,
            "exp_duration": (int(params["exp_duration"]) if params["exp_duration_enabled"] else None),
            "min_track_length": (int(params["min_track_length"]) if params["min_length_enabled"] else None),
            "max_track_length": (int(params["max_track_length"]) if params["max_length_enabled"] else None),
            "df_input_path": df_input_path,
            "time_type": params["time_type"],
            "plot_results": True,
            "filter_t0_dead": bool(params["filter_t0_dead"]),
        }

        if params["filter_min_size"]:
            filter_kwargs["min_size"] = int(params["min_size_t1"])

        filter_tracks(**filter_kwargs)

        print(f"  Summarizing {cell_type} tracks…", file=sys.stderr)
        summarize_track_features(output_dir=out_dir, cell_type=cell_type)
        print(f"✅ {cell_type} filtering & summarization finished.", file=sys.stderr)

    def _check_existing_filtering(self, cell_types: list) -> list:
        warnings = []
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in cell_types:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            filtered = feat_dir / f"BEHAV3D_{ct}_filtered_track_features.csv"
            if filtered.exists():
                warnings.append(f"{ct} filtered data ({filtered.name})")
        return warnings

    def _on_run_clicked(self, interactive=True):
        """Run filtering for this cell type in the background.

        ``filter_tracks`` has no per-sample loop accessible from the
        top-level call, so the in-tab progress row stays in indeterminate
        (busy) mode for the duration of the run.
        """
        if self._bg.is_running():
            self.log("⚠️ A filtering run is already in progress for this panel.")
            return

        self._persist()
        self.log(f"Running filtering for: {self.cell_type}")

        overwrite = False
        existing = self._check_existing_filtering([self.cell_type])
        if existing:
            if interactive:
                from behav3d.napari._overwrite_prompt import prompt_overwrite_single
                choice = prompt_overwrite_single(
                    self,
                    "Overwrite Existing Filtered Data?",
                    existing,
                )
                if choice != "overwrite":
                    self.log(f"Filtering for {self.cell_type} cancelled.")
                    return
                overwrite = True
            else:
                overwrite = True

        params_snapshot = self.collect_runtime_params()
        cell_type = self.cell_type

        def _do_filtering(progress_cb=None):
            # filter_tracks has no per-sample loop, so we don't drive
            # progress_cb here — the in-tab bar remains indeterminate.
            return self._run_filtering_for(
                cell_type, overwrite=overwrite, params=params_snapshot,
            )

        def _on_done(_result):
            self.log(f"✅ {cell_type} filtering finished.")
            notify_results_changed(self)

        def _on_failed(err: str):
            self.log(f"Error during filtering: {err}")
            notify_results_changed(self)

        self._bg.run(
            fn=_do_filtering,
            desc=f"Filtering {cell_type}\u2026",
            progress_row=self.tab_progress_row,
            buttons=[self.btn_run],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
            inject_progress=False,
            indeterminate=True,
        )


# ═══════════════════════════════════════════════════════════════════════════
# FilteringTab — main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class FilteringTab(QWidget):
    # Worker threads call ``self._log(...)`` during batch filtering; QTextEdit
    # is not thread-safe and direct ``append`` calls from a non-GUI thread
    # raise the "Cannot queue arguments of type 'QTextCursor'" warning.  We
    # route every log message through this queue-connected signal so the
    # actual ``QTextEdit.append`` always runs on the GUI thread.
    _log_message = Signal(str)

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeFilterPanel] = {}

        # Background-execution infrastructure for batch runs.
        self._bg = BackgroundOperation(self)

        self._init_ui()

        # Queue-connect so cross-thread emits get marshalled to the GUI thread.
        self._log_message.connect(self._append_log_line, Qt.QueuedConnection)

        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Vertical splitter: tab content on top, shared Results panel
        # underneath (visible from this tab and re-scanning the same
        # output directory used by Analysis / Feature Extraction).
        self.splitter = QSplitter(Qt.Vertical, self)
        self.splitter.setChildrenCollapsible(True)
        outer.addWidget(self.splitter)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.splitter.addWidget(scroll)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        scroll.setWidget(content)

        self.results_panel = ResultsPanel(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            parent=self,
        )
        self.splitter.addWidget(self.results_panel)
        self.splitter.setStretchFactor(0, 3)
        self.splitter.setStretchFactor(1, 2)
        self.splitter.setCollapsible(0, False)
        self.splitter.setCollapsible(1, True)
        self.splitter.setSizes([600, 400])

        self.cell_tabs = QTabWidget()
        self.cell_tabs.setTabPosition(QTabWidget.West)
        layout.addWidget(self.cell_tabs)

        self.btn_run_batch = QPushButton("Run Batch Filtering (All Cell Types)")
        self.btn_run_batch.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_batch.clicked.connect(self._on_run_batch_clicked)

        self.btn_queue_filter = QPushButton("+🛒")
        self.btn_queue_filter.setFixedSize(36, 32)
        self.btn_queue_filter.setToolTip("Add Filtering to Processing Queue")
        self.btn_queue_filter.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_row = QHBoxLayout()
        batch_row.setSpacing(4)
        batch_row.addWidget(self.btn_run_batch, stretch=1)
        batch_row.addWidget(self.btn_queue_filter)
        layout.addLayout(batch_row)

        # Hidden until metadata is loaded
        self.btn_run_batch.setVisible(False)
        self.btn_queue_filter.setVisible(False)

        # Shared progress row — fed by every Run button on this tab.
        self.progress_row = ProgressBarRow()
        layout.addWidget(self.progress_row)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(120)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        self._placeholder = QLabel("Load metadata in the Data Preparation tab to see filtering options.")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.cell_tabs.addTab(self._placeholder, "—")

    def _log(self, msg):
        """Thread-safe log entrypoint.

        Emits via ``_log_message`` (queue-connected) so calls from worker
        threads — e.g. ``_do_batch`` running on a ``QThread`` — never touch
        the ``QTextEdit`` directly.  Calls from the GUI thread still work
        because Qt processes queued events between worker steps.
        """
        try:
            self._log_message.emit(str(msg))
        except Exception:
            # As a last resort, fall back to a direct write; this only
            # matters on the GUI thread and avoids losing log lines if the
            # signal infrastructure is somehow torn down.
            try:
                self._append_log_line(str(msg))
            except Exception:
                pass

    def _append_log_line(self, msg: str):
        """Slot — always invoked on the GUI thread."""
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing filtering tabs…")
        self._rebuild_tabs()

    def request_tab_exit(self) -> bool:
        """Block tab switching while a background filtering run is in flight."""
        from qtpy.QtWidgets import QMessageBox

        if self._bg.is_running():
            QMessageBox.information(
                self,
                "Operation in progress",
                "A filtering run is still in progress. Please wait for "
                "it to finish before switching tabs.",
            )
            return False
        for panel in self.panels.values():
            bg = getattr(panel, "_bg", None)
            if bg is not None and bg.is_running():
                QMessageBox.information(
                    self,
                    "Operation in progress",
                    "A filtering run is still in progress. Please wait "
                    "for it to finish before switching tabs.",
                )
                return False
        return True

    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            filter_multicolor_inputs,
        )
        md = self.metadata_loader.metadata
        if md is None:
            return [], [], []
        return (
            filter_multicolor_inputs(detect_organoid_types_from_metadata(md)),
            filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)),
            filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)),
        )

    def _rebuild_tabs(self):
        self.cell_tabs.clear()
        self.panels.clear()

        org, imm, oth = self._detect_cell_types()
        all_types = org + imm + oth

        if not all_types:
            self.cell_tabs.addTab(self._placeholder, "—")
            self.btn_run_batch.setVisible(False)
            self.btn_queue_filter.setVisible(False)
            return

        self.btn_run_batch.setVisible(True)
        self.btn_queue_filter.setVisible(True)

        color_map = {"organoid": "🟣", "immune": "🔵", "other": "🟡"}

        for ct in org:
            self._add_panel(ct, "organoid", all_types, org, color_map)
        for ct in imm:
            self._add_panel(ct, "immune", all_types, imm, color_map)
        for ct in oth:
            self._add_panel(ct, "other", all_types, oth, color_map)

    def _add_panel(self, ct, category, all_types, cat_types, color_map):
        panel = CellTypeFilterPanel(
            cell_type=ct, category=category,
            metadata_loader=self.metadata_loader,
            all_cell_types=all_types, category_types=cat_types,
            log_callback=self._log,
            tab_progress_row=self.progress_row,
            viewer=self.viewer,
        )
        self.panels[ct] = panel
        icon = color_map.get(category, "")
        self.cell_tabs.addTab(panel, f"{icon} {ct}")

    def _on_run_batch_clicked(self):
        """User-triggered batch run — asynchronous."""
        self.run_batch_filtering(interactive=True, block=False)

    def run_batch_filtering(self, interactive=True, skip_existing=False, block=True,
                            extra_callbacks=None):
        """Run filtering for all cell types sequentially.

        ``block=True`` (default, queue path) runs synchronously on the
        caller's thread.  ``block=False`` (GUI button) runs in a background
        worker with a determinate progress bar at cell-type granularity.

        ``extra_callbacks`` is the queue's chaining hook
        (``{"on_done": cb, "on_failed": cb}``).
        """
        if not self.panels:
            self._log("No cell type panels available.")
            fire_extra_callback(extra_callbacks, "on_failed", "no cell type panels")
            return

        if not block and self._bg.is_running():
            self._log("⚠️ A batch filtering run is already in progress.")
            fire_extra_callback(extra_callbacks, "on_failed", "already running")
            return

        total = len(self.panels)
        self._log(f"Starting batch filtering for {total} cell types…")

        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Running batch filtering for {total} cell types", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        # Overwrite check
        all_cts = list(self.panels.keys())
        existing = []
        existing_cts = set()
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in all_cts:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            filtered = feat_dir / f"BEHAV3D_{ct}_filtered_track_features.csv"
            if filtered.exists():
                existing.append(f"{ct} filtered data ({filtered.name})")
                existing_cts.add(ct)

        skip_existing_flag = skip_existing
        overwrite = not skip_existing
        if existing:
            if interactive:
                from behav3d.napari._overwrite_prompt import prompt_overwrite_batch
                choice = prompt_overwrite_batch(
                    self,
                    "Overwrite Existing Filtered Data?",
                    existing,
                )
                if choice == "cancel":
                    self._log("Batch filtering cancelled.")
                    fire_extra_callback(extra_callbacks, "on_failed", "cancelled")
                    return
                skip_existing_flag = (choice == "skip")
                overwrite = not skip_existing_flag
            else:
                overwrite = True

        # Persist + snapshot every panel's Qt widget state on the Qt
        # thread — the worker must not read widgets.
        panel_params: dict[str, dict] = {}
        for ct, panel in self.panels.items():
            panel._persist()
            panel_params[ct] = panel.collect_runtime_params()

        def _do_batch(progress_cb=None):
            for i, (ct, panel) in enumerate(self.panels.items()):
                if progress_cb is not None:
                    try:
                        progress_cb(i, total, f"Filtering: {ct}")
                    except Exception:
                        pass
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- [{i + 1}/{total}] Skipping {ct} (existing data) ---")
                    continue
                print(f"\n▶ [{i + 1}/{total}] Filtering: {ct}…", file=sys.stderr)
                self._log(f"--- [{i + 1}/{total}] Filtering: {ct} ---")
                panel._run_filtering_for(
                    ct, overwrite=overwrite,
                    params=panel_params[ct],
                )
                self._log(f"Done: {ct}")
            if progress_cb is not None:
                try:
                    progress_cb(total, total, "done")
                except Exception:
                    pass

            self._log("✅ Batch filtering finished.")
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"  ✅ Batch filtering complete", file=sys.stderr)
            print(f"{'='*60}\n", file=sys.stderr)
            return None

        if block:
            try:
                result = _do_batch(progress_cb=None)
                fire_extra_callback(extra_callbacks, "on_done", result)
            except Exception as e:
                traceback.print_exc()
                self._log(f"❌ Batch filtering error: {e}")
                fire_extra_callback(extra_callbacks, "on_failed", str(e))
            finally:
                try:
                    self.results_panel.refresh()
                except Exception:
                    traceback.print_exc()
            return

        def _on_done(result):
            try:
                self.results_panel.refresh()
            except Exception:
                traceback.print_exc()
            fire_extra_callback(extra_callbacks, "on_done", result)

        def _on_failed(err: str):
            self._log(f"❌ Batch filtering error: {err}")
            try:
                self.results_panel.refresh()
            except Exception:
                traceback.print_exc()
            fire_extra_callback(extra_callbacks, "on_failed", err)

        self._bg.run(
            fn=_do_batch,
            desc=f"Batch filtering ({total} cell types)\u2026",
            progress_row=self.progress_row,
            buttons=[self.btn_run_batch],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
        )
