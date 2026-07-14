"""Global workers controller for the BEHAV3D napari plugin.

A single :class:`GlobalWorkersController` instance is created once and shared
across all pipeline tabs.  Every ``spin_workers`` spinbox in the UI is
*linked* to it: changing any one spinbox updates all the others and persists
the value to the top-level ``n_workers`` key in ``behav3d_parameters.yml``.

Usage
-----
Create the controller in the main plugin widget::

    from behav3d.napari._global_workers import GlobalWorkersController
    self.workers_ctrl = GlobalWorkersController(metadata_loader=self.data_prep_tab)

Link each tab's spinbox during tab construction::

    self.workers_ctrl.link(self.segmentation_tab.spin_workers)
    self.workers_ctrl.link(self.feature_extraction_tab.spin_workers)
    self.workers_ctrl.link(self.edition_tab.spin_workers)
    # …etc.

The controller automatically clamps each linked spinbox to its own range, so
a spinbox that only goes up to 8 will be set to ``min(n_workers, 8)`` when
the global value changes.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

from qtpy.QtWidgets import QSpinBox

if TYPE_CHECKING:
    pass


class GlobalWorkersController:
    """Observable controller that keeps all ``spin_workers`` spinboxes in sync.

    Parameters
    ----------
    metadata_loader:
        The :class:`~behav3d.napari._data_preparation.DataPreparationTab`
        instance.  Used to access ``output_dir`` and
        ``behav3d_parameters`` (the shared in-memory YAML dict).
    default_workers:
        Initial value used when no saved value is found.
    """

    _KEY = "n_workers"

    def __init__(self, metadata_loader, default_workers: int | None = None) -> None:
        self._loader = metadata_loader
        if default_workers is None:
            import multiprocessing as _mp
            n_cores = max(1, _mp.cpu_count())
            default_workers = max(1, n_cores // 2)
        self._value: int = default_workers
        self._spinboxes: list[QSpinBox] = []
        self._updating: bool = False  # re-entrancy guard

        # Try to read saved value from YAML right away.
        self._reload_from_params()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    @property
    def value(self) -> int:
        return self._value

    def link(self, spinbox: QSpinBox) -> None:
        """Register *spinbox* and keep it in sync with the global value.

        The spinbox's current value is set to ``clamp(global_value, min, max)``.
        When the spinbox changes, the global value and all other spinboxes are
        updated; the new value is also saved to YAML.
        """
        if spinbox in self._spinboxes:
            return
        self._spinboxes.append(spinbox)
        self._apply_to_spinbox(spinbox, self._value)
        spinbox.valueChanged.connect(lambda v, s=spinbox: self._on_spinbox_changed(v, s))

    def set_value(self, value: int, save: bool = True) -> None:
        """Set the global worker count programmatically and propagate it."""
        self._value = max(1, int(value))
        self._propagate(source=None)
        if save:
            self._save_to_params()

    def reload(self) -> None:
        """Re-read ``n_workers`` from the shared parameter dict and propagate."""
        self._reload_from_params()
        self._propagate(source=None)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _on_spinbox_changed(self, value: int, source: QSpinBox) -> None:
        if self._updating:
            return
        self._value = max(1, int(value))
        self._propagate(source=source)
        self._save_to_params()

    def _propagate(self, source: QSpinBox | None) -> None:
        self._updating = True
        try:
            for sb in self._spinboxes:
                if sb is source:
                    continue  # already has the right value
                self._apply_to_spinbox(sb, self._value)
        finally:
            self._updating = False

    @staticmethod
    def _apply_to_spinbox(sb: QSpinBox, value: int) -> None:
        clamped = max(sb.minimum(), min(sb.maximum(), value))
        if sb.value() != clamped:
            sb.setValue(clamped)

    def _save_to_params(self) -> None:
        """Persist *n_workers* to the top-level key of ``behav3d_parameters``."""
        params = getattr(self._loader, "behav3d_parameters", None)
        if params is None:
            return
        params[self._KEY] = self._value
        out_dir = getattr(self._loader, "output_dir", "")
        if not out_dir:
            return
        params_path = Path(out_dir) / "behav3d_parameters.yml"
        try:
            import yaml
            with open(params_path, "w") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        except Exception:
            pass  # best-effort

    def _reload_from_params(self) -> None:
        """Read saved *n_workers* from the shared parameter dict.

        First looks for the top-level ``n_workers`` key (written by this
        controller).  If absent, falls back to the first non-``None`` value
        found in the older per-tab keys so that existing parameter files from
        before the global workers feature was introduced still work correctly.
        """
        params = getattr(self._loader, "behav3d_parameters", None)
        if not params:
            return
        saved = params.get(self._KEY)
        if saved is None:
            # Fallback: check legacy per-tab locations in priority order.
            for _key_path in (
                ("pixel_classifier", "workers"),
                ("feature_extraction", "n_workers"),
                ("apoc", "workers"),
            ):
                section = params
                for k in _key_path:
                    if not isinstance(section, dict):
                        section = None
                        break
                    section = section.get(k)
                if section is not None:
                    saved = section
                    break
        if saved is not None:
            try:
                self._value = max(1, int(saved))
            except (TypeError, ValueError):
                pass
