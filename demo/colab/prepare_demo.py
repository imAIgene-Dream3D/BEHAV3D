"""Repoint a BEHAV3D demo bundle at the machine it was just downloaded onto.

BEHAV3D stores **absolute** paths: every ``*_path`` column of ``metadata.csv``
(``raw_image_path``, ``im_<type>_segments_image_path``, ``*_tracks_csv_path``,
...) and the ``paths.metadata_csv`` / ``paths.output_dir`` keys of
``<output_dir>/behav3d_parameters.yml``. Those paths point at whatever machine
produced the bundle, so they must be rewritten before the GUI can open it.

The script finds the longest common directory of all recorded paths (that is
the bundle root on the producing machine) and swaps it for the local root. A
hand-authored bundle may instead use the literal token ``{DEMO_ROOT}``, which
is substituted directly.

Usage::

    python prepare_demo.py --root /content/behav3d_demo
    python prepare_demo.py --root /content/behav3d_demo --dry-run
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path, PurePosixPath

import pandas as pd
import yaml

TOKEN = "{DEMO_ROOT}"
PATH_SUFFIX = "_path"


def _norm(value: str) -> str:
    """Windows-produced bundles carry backslashes; POSIX-ify them."""
    return str(value).replace("\\", "/")


def _is_pathlike(value) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False
    text = _norm(value)
    return "/" in text or text.startswith(TOKEN)


def _path_columns(frame: pd.DataFrame):
    return [c for c in frame.columns if c.endswith(PATH_SUFFIX)]


def _collect(frame: pd.DataFrame):
    values = []
    for col in _path_columns(frame):
        values += [_norm(v) for v in frame[col].dropna().tolist() if _is_pathlike(v)]
    return values


def _detect_old_root(values) -> str | None:
    """Longest common parent directory of all recorded paths."""
    absolute = [v for v in values if not v.startswith(TOKEN)]
    if not absolute:
        return None
    if len(absolute) == 1:
        return str(PurePosixPath(absolute[0]).parent)
    try:
        return os.path.commonpath(absolute).replace("\\", "/")
    except ValueError:      # mixed drive letters / mixed absolute-relative
        return None


def _rewrite(value, old_root: str | None, new_root: str):
    if not _is_pathlike(value):
        return value, False
    text = _norm(value)
    if TOKEN in text:
        return text.replace(TOKEN, new_root).replace("//", "/"), True
    if old_root and text.startswith(old_root):
        return new_root + text[len(old_root):], True
    return text, False


def rewrite_metadata(csv_path: Path, new_root: str, dry_run=False):
    frame = pd.read_csv(csv_path, low_memory=False)
    old_root = _detect_old_root(_collect(frame))
    if old_root == new_root:
        print(f"  metadata.csv already points at {new_root} - unchanged")
        return 0, old_root

    changed = 0
    for col in _path_columns(frame):
        new_values = []
        for value in frame[col]:
            rewritten, did = _rewrite(value, old_root, new_root)
            changed += bool(did)
            new_values.append(rewritten)
        frame[col] = new_values

    if not dry_run and changed:
        backup = csv_path.with_suffix(".original.csv")
        if not backup.exists():
            csv_path.replace(backup)
        frame.to_csv(csv_path, index=False)
    return changed, old_root


def rewrite_parameters(yml_path: Path, new_root: str, csv_path: Path, dry_run=False):
    """Point ``behav3d_parameters.yml`` at the local metadata CSV and output dir."""
    config = yaml.safe_load(yml_path.read_text()) or {} if yml_path.exists() else {}
    config.setdefault("paths", {})
    config["paths"]["metadata_csv"] = _norm(csv_path)
    config["paths"]["output_dir"] = _norm(yml_path.parent)
    if not dry_run:
        yml_path.parent.mkdir(parents=True, exist_ok=True)
        yml_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return config["paths"]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True,
                        help="local folder holding raw/, metadata.csv and output/")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would change without writing")
    args = parser.parse_args(argv)

    root = Path(args.root).resolve()
    new_root = _norm(root)
    csv_path = root / "metadata.csv"
    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found - is this a BEHAV3D demo bundle?", file=sys.stderr)
        return 1

    print(f"Preparing demo bundle at {new_root}")
    changed, old_root = rewrite_metadata(csv_path, new_root, args.dry_run)
    if old_root and old_root != new_root:
        print(f"  rewrote {changed} path(s): {old_root} -> {new_root}")
    elif not old_root:
        print("  no absolute paths detected in metadata.csv (nothing to rewrite)")

    output_dir = root / "output"
    paths = rewrite_parameters(output_dir / "behav3d_parameters.yml",
                               new_root, csv_path, args.dry_run)
    print(f"  behav3d_parameters.yml -> metadata_csv={paths['metadata_csv']}")

    missing = _check_referenced_files(csv_path)
    if missing:
        print(f"  WARNING: {len(missing)} referenced file(s) are missing, e.g. {missing[0]}")
    else:
        print("  all referenced files exist")

    print("\nOpen these in the BEHAV3D Explorer 'Data Preparation' tab:")
    print(f"  Metadata CSV : {csv_path}")
    print(f"  Output folder: {output_dir}")
    return 0


def _check_referenced_files(csv_path: Path, limit=5):
    frame = pd.read_csv(csv_path, low_memory=False)
    missing = []
    for value in _collect(frame):
        if not Path(value).exists():
            missing.append(value)
            if len(missing) >= limit:
                break
    return missing


if __name__ == "__main__":
    raise SystemExit(main())
