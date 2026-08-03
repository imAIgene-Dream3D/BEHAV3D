from pathlib import Path

import pandas as pd

from behav3d.analysis.organoid_analysis import resolve_organoid_seg_path


def _sample_metadata(**cols):
    return pd.DataFrame([cols])


def test_resolves_from_metadata_when_available(tmp_path):
    seg_path = "/data/exp84/images/S1/S1_organoid_tracked.zarr"
    sample_metadata = _sample_metadata(
        sample_name="S1", or_organoid_tracks_image_path=seg_path
    )
    img_outdir = tmp_path / "combined_analysis" / "images" / "S1"

    resolved = resolve_organoid_seg_path(sample_metadata, "organoid", img_outdir, "S1")

    assert resolved == Path(seg_path)


def test_resolves_immune_cell_type_via_im_prefix(tmp_path):
    seg_path = "/data/exp84/images/S1/S1_macrophage_tracked.zarr"
    sample_metadata = _sample_metadata(
        sample_name="S1", im_macrophage_tracks_image_path=seg_path
    )
    img_outdir = tmp_path / "combined_analysis" / "images" / "S1"

    resolved = resolve_organoid_seg_path(sample_metadata, "macrophage", img_outdir, "S1")

    assert resolved == Path(seg_path)


def test_falls_back_to_constructed_path_when_metadata_empty(tmp_path):
    img_outdir = tmp_path / "images" / "S1"

    resolved = resolve_organoid_seg_path(pd.DataFrame(), "organoid", img_outdir, "S1")

    assert resolved == img_outdir / "S1_organoid_tracked.zarr"


def test_falls_back_to_constructed_path_when_column_missing(tmp_path):
    sample_metadata = _sample_metadata(sample_name="S1", raw_image_path="/data/S1.czi")
    img_outdir = tmp_path / "images" / "S1"

    resolved = resolve_organoid_seg_path(sample_metadata, "organoid", img_outdir, "S1")

    assert resolved == img_outdir / "S1_organoid_tracked.zarr"


def test_falls_back_to_constructed_path_when_metadata_value_blank(tmp_path):
    sample_metadata = _sample_metadata(sample_name="S1", or_organoid_tracks_image_path="")
    img_outdir = tmp_path / "images" / "S1"

    resolved = resolve_organoid_seg_path(sample_metadata, "organoid", img_outdir, "S1")

    assert resolved == img_outdir / "S1_organoid_tracked.zarr"
