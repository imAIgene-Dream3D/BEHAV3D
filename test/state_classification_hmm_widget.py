"""Compatibility wrapper for the promoted production HMM state widget."""

from behav3d.widgets.state_classification import (
    StateClassificationHMMDeploymentPanel,
    StateClassificationHMMPanel,
    StateClassificationPanel,
    _winfo,
)

__all__ = [
    "StateClassificationHMMDeploymentPanel",
    "StateClassificationHMMPanel",
    "StateClassificationPanel",
    "_winfo",
]
