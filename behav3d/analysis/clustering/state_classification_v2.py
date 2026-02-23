"""Compatibility wrappers for legacy v2 API.

Use `behav3d.analysis.clustering.state_classification` instead.
"""

from behav3d.analysis.clustering.state_classification import (
    run_state_classification,
    apply_classifier,
    load_state_classifier_artifact,
    build_state_preprocessing_artifact,
    plot_binary_group_behavioral_cluster_grid,
    plot_behavioral_cluster_binary_group_grid,
)


def run_state_classification_v2(*args, **kwargs):
    return run_state_classification(*args, **kwargs)


def load_state_v2_classifier_artifact(*args, **kwargs):
    return load_state_classifier_artifact(*args, **kwargs)


def build_state_v2_preprocessing_artifact(*args, **kwargs):
    return build_state_preprocessing_artifact(*args, **kwargs)
