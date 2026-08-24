"""Tests for the shared segmentation progress journal.

These cover the bookkeeping only - no GPU, no classifiers, no images. The rules
being pinned down here are the ones that decide whether a run reuses frames on disk
or recomputes them, so getting them wrong silently corrupts outputs.
"""
from pathlib import Path
import shutil
import uuid

import numpy as np
import pytest
import zarr

from behav3d.preprocessing.segmentation.segment_journal import (
    JOURNAL_VERSION,
    check_fingerprint,
    describe_state,
    file_fingerprint,
    journal_path,
    journal_state,
    mark_done,
    new_journal,
    open_output_zarr,
    params_fingerprint,
    plan_output,
    preflight_conflicts,
    raise_for_conflicts,
    read_journal,
    should_recreate,
    write_journal,
)


@pytest.fixture
def case_dir():
    root = Path(__file__).resolve().parent / ".tmp_segment_journal"
    root.mkdir(exist_ok=True)
    path = root / uuid.uuid4().hex
    path.mkdir()
    yield path
    shutil.rmtree(path, ignore_errors=True)


def _make_zarr(path, shape=(5, 2, 4, 4), dtype="uint16"):
    arr = zarr.open(str(path), mode="w", shape=shape, dtype=dtype,
                    chunks=(1,) + tuple(shape[1:]))
    arr[:] = 1
    return arr


# ── round-trip and corruption ────────────────────────────────────────────────

def test_journal_round_trip(case_dir):
    p = case_dir / "a.progress.json"
    j = new_journal("fp1", (10, 2, 3, 3), done=[0, 2, 1])
    write_journal(p, j)
    loaded = read_journal(p)
    assert loaded["fingerprint"] == "fp1"
    assert loaded["done"] == [0, 1, 2]
    assert loaded["shape"] == [10, 2, 3, 3]


def test_truncated_journal_reads_as_none(case_dir):
    """A half-written journal is the expected casualty of a power cut."""
    p = case_dir / "a.progress.json"
    write_journal(p, new_journal("fp1", (10, 2, 3, 3), done=[0, 1]))
    text = p.read_text(encoding="utf-8")
    p.write_text(text[: len(text) // 2], encoding="utf-8")
    assert read_journal(p) is None


def test_version_mismatch_reads_as_none(case_dir):
    p = case_dir / "a.progress.json"
    j = new_journal("fp1", (10, 2, 3, 3))
    j["version"] = JOURNAL_VERSION + 1
    write_journal(p, j)
    assert read_journal(p) is None


def test_missing_journal_reads_as_none(case_dir):
    assert read_journal(case_dir / "nope.progress.json") is None


def test_mark_done_persists_every_frame(case_dir):
    p = case_dir / "a.progress.json"
    j = new_journal("fp1", (3, 1, 2, 2))
    for t in range(3):
        mark_done(j, p, t)
        assert read_journal(p)["done"] == list(range(t + 1))


# ── fingerprints ─────────────────────────────────────────────────────────────

def test_fingerprint_is_key_order_insensitive():
    a = params_fingerprint({"edt_threshold": 1.0, "size": 10})
    b = params_fingerprint({"size": 10, "edt_threshold": 1.0})
    assert a == b


def test_fingerprint_changes_with_a_segmentation_param():
    a = params_fingerprint({"edt_threshold": 1.0, "size": 10})
    b = params_fingerprint({"edt_threshold": 1.8, "size": 10})
    assert a != b


def test_fingerprint_ignores_run_bookkeeping():
    """n_workers and friends never reach the payload, so they cannot invalidate."""
    payload = {"edt_threshold": 1.0}
    assert params_fingerprint(payload) == params_fingerprint(dict(payload))


def test_file_fingerprint_tracks_content(case_dir):
    f = case_dir / "clf.cl"
    f.write_text("alpha", encoding="utf-8")
    first = file_fingerprint(f)
    f.write_text("beta", encoding="utf-8")
    assert file_fingerprint(f) != first
    assert file_fingerprint(case_dir / "absent.cl") is None


# ── check_fingerprint ────────────────────────────────────────────────────────

def test_check_fingerprint_passes_on_match():
    check_fingerprint(new_journal("fp1", (3, 1, 1, 1)), "fp1", "x")


def test_check_fingerprint_raises_on_mismatch():
    with pytest.raises(RuntimeError, match="settings have changed"):
        check_fingerprint(new_journal("fp1", (3, 1, 1, 1)), "fp2", "tcell segments for S1")


def test_check_fingerprint_ignores_missing_journal():
    """Legacy arrays predate journalling and must never raise."""
    check_fingerprint(None, "fp1", "x")


# ── journal_state ────────────────────────────────────────────────────────────

def test_state_missing(case_dir):
    assert journal_state(case_dir / "absent.zarr")[0] == "missing"


def test_state_unknown_without_journal(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z)
    assert journal_state(z) == ("unknown", 0, 0)


def test_state_partial_then_complete(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    j = new_journal("fp1", (4, 1, 2, 2), done=[0, 1])
    write_journal(journal_path(z), j)
    assert journal_state(z) == ("partial", 2, 4)

    mark_done(j, journal_path(z), 2)
    mark_done(j, journal_path(z), 3)
    assert journal_state(z) == ("complete", 4, 4)


def test_describe_state_wording():
    assert describe_state("partial", 47, 120) == "INCOMPLETE (47/120)"
    assert describe_state("complete", 120, 120) == "complete (120/120)"
    assert describe_state("unknown", 0, 0) == "unknown (no journal)"


# ── should_recreate ──────────────────────────────────────────────────────────

def test_recreate_when_absent(case_dir):
    assert should_recreate(case_dir / "absent.zarr", (4, 1, 2, 2), "uint16",
                           overwrite_existing=False, timepoint_range=None)


def test_recreate_on_full_range_overwrite(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    assert should_recreate(z, (4, 1, 2, 2), "uint16",
                           overwrite_existing=True, timepoint_range=None)


def test_never_recreate_on_sub_range_overwrite(case_dir):
    """Wiping the array for a sub-range run would zero every frame outside it."""
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    assert not should_recreate(z, (4, 1, 2, 2), "uint16",
                               overwrite_existing=True, timepoint_range=(0, 1))


def test_recreate_on_shape_mismatch(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    assert should_recreate(z, (9, 1, 2, 2), "uint16",
                           overwrite_existing=False, timepoint_range=None)


# ── plan_output ──────────────────────────────────────────────────────────────

def _plan(z, shape, fp, **kw):
    kw.setdefault("overwrite_existing", False)
    kw.setdefault("timepoint_range", None)
    kw.setdefault("requested_timepoints", list(range(shape[0])))
    return plan_output(z, shape, "uint16", fp, "tcell segments for S1", **kw)


def test_plan_legacy_array_counts_as_complete(case_dir):
    """No journal + full range: skip, exactly as before journalling existed."""
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    plan = _plan(z, (4, 1, 2, 2), "fp1")
    assert plan.legacy and plan.complete and not plan.recreate


def test_plan_resumes_a_partial_journal(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp1", (4, 1, 2, 2), done=[0, 1]))
    plan = _plan(z, (4, 1, 2, 2), "fp1")
    assert plan.done == {0, 1}
    assert not plan.complete and not plan.recreate


def test_plan_complete_journal_skips(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp1", (4, 1, 2, 2), done=[0, 1, 2, 3]))
    assert _plan(z, (4, 1, 2, 2), "fp1").complete


def test_plan_raises_on_changed_settings(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp_old", (4, 1, 2, 2), done=[0, 1]))
    with pytest.raises(RuntimeError, match="settings have changed"):
        _plan(z, (4, 1, 2, 2), "fp_new")


def test_plan_overwrite_ignores_the_journal_entirely(case_dir):
    """Overwrite means redo; a stale fingerprint is not an error there."""
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp_old", (4, 1, 2, 2), done=[0, 1, 2, 3]))
    plan = _plan(z, (4, 1, 2, 2), "fp_new", overwrite_existing=True)
    assert plan.done == set() and not plan.complete and plan.recreate


def test_plan_sub_range_overwrite_keeps_the_array(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    plan = _plan(z, (4, 1, 2, 2), "fp1", overwrite_existing=True,
                 timepoint_range=(0, 1), requested_timepoints=[0, 1])
    assert not plan.recreate


# ── the defect this replaces ─────────────────────────────────────────────────

def test_sub_range_overwrite_preserves_out_of_range_frames(case_dir):
    """The old marker scheme wiped the whole array and kept claiming it was done."""
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))  # every frame filled with 1s

    recreate = should_recreate(z, (4, 1, 2, 2), "uint16",
                               overwrite_existing=True, timepoint_range=(0, 1))
    arr = open_output_zarr(z, "uint16", (4, 1, 2, 2), recreate)
    arr[0] = 7
    arr[1] = 7

    assert np.all(np.asarray(arr[2]) == 1), "frames outside the range were destroyed"
    assert np.all(np.asarray(arr[3]) == 1)


def test_recreating_an_array_discards_its_journal(case_dir):
    """A journal that outlives its data is the bug this module exists to end."""
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp1", (4, 1, 2, 2), done=[0, 1, 2, 3]))
    open_output_zarr(z, "uint16", (4, 1, 2, 2), recreate=True)
    assert not journal_path(z).exists()


# ── preflight ────────────────────────────────────────────────────────────────

def test_preflight_collects_every_conflict(case_dir):
    entries = []
    for name in ("a", "b"):
        z = case_dir / f"{name}.zarr"
        _make_zarr(z, shape=(4, 1, 2, 2))
        write_journal(journal_path(z), new_journal("fp_old", (4, 1, 2, 2), done=[0]))
        entries.append((z, (4, 1, 2, 2), "uint16", "fp_new", f"{name} segments for S1"))

    conflicts = preflight_conflicts(entries, overwrite_existing=False, timepoint_range=None)
    assert len(conflicts) == 2

    with pytest.raises(RuntimeError) as exc:
        raise_for_conflicts(conflicts, "APOC")
    msg = str(exc.value)
    assert "a segments for S1" in msg and "b segments for S1" in msg
    assert "Nothing has been modified." in msg


def test_preflight_is_silent_under_overwrite(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    write_journal(journal_path(z), new_journal("fp_old", (4, 1, 2, 2), done=[0]))
    entries = [(z, (4, 1, 2, 2), "uint16", "fp_new", "a segments for S1")]
    assert preflight_conflicts(entries, overwrite_existing=True, timepoint_range=None) == []


def test_preflight_ignores_legacy_arrays(case_dir):
    z = case_dir / "a.zarr"
    _make_zarr(z, shape=(4, 1, 2, 2))
    entries = [(z, (4, 1, 2, 2), "uint16", "fp_new", "a segments for S1")]
    assert preflight_conflicts(entries, overwrite_existing=False, timepoint_range=None) == []


def test_raise_for_conflicts_is_a_noop_when_clean():
    raise_for_conflicts([], "APOC")
