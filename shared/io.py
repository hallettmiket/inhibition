"""
Purpose: Reading and writing approach frames, with the contract enforced at the door.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: candidate frames from any approach
Output: validated parquet/CSV under the governed append-only tree, plus manifests

WHY WRITES GO THROUGH HERE. Two rules are easy to violate by accident and
expensive to discover late:

1. **Append-only means append-only.** Outputs under
   `/data/lab_vm/append_only/inhibition/` are never overwritten. A re-run writes
   the next integer version; the previous one is retired through
   `data/ready_to_delete.md`, not deleted in place. That is what makes an old
   result still checkable against the file it actually used.

2. **A frame that fails the schema must not reach disk.** Validating on read is
   too late — by then the GUI is already pooling something malformed, and the
   approach that produced it may have finished hours ago.

Every write also emits a manifest, so "what produced this file" is answerable
without asking whoever ran it.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from . import schema
from .manifest import Manifest

log = logging.getLogger(__name__)

APPEND_ONLY_ROOT = Path("/data/lab_vm/append_only/inhibition")
_VERSIONED = re.compile(r"^(?P<stem>.+)_(?P<version>\d+)$")


class IOError_(RuntimeError):
    """A read or write violated the storage rules."""


def next_version(directory: Path, stem: str, suffix: str) -> Path:
    """Return the next unused integer-versioned path, e.g. ``D4_2.parquet``.

    Integer versioning per the lab data rule: the largest integer is newest, and
    old versions are retired via ready_to_delete.md rather than overwritten.
    """
    directory.mkdir(parents=True, exist_ok=True)
    existing = []
    for p in directory.glob(f"{stem}_*{suffix}"):
        m = _VERSIONED.match(p.stem)
        if m and m.group("stem") == stem:
            try:
                existing.append(int(m.group("version")))
            except ValueError:
                continue
    return directory / f"{stem}_{max(existing, default=0) + 1}{suffix}"


def approach_dir(approach: str, experiment: str) -> Path:
    """Governed output directory for one approach's experiment."""
    return APPEND_ONLY_ROOT / experiment


def write_full_frame(df: pd.DataFrame, *, approach: str, experiment: str,
                     stage: str, params: dict | None = None,
                     inputs: dict[str, Path] | None = None) -> Path:
    """Validate and write an approach's FULL candidate frame as parquet.

    Raises
    ------
    SchemaError
        If the frame violates the D^i contract — before anything is written.
    """
    report = schema.validate_full(df, approach)
    for w in report.warnings:
        log.warning("[%s] %s", approach, w)
    report.raise_if_bad(f"{approach} full frame")

    d = approach_dir(approach, experiment)
    out = next_version(d, f"D{approach[-1]}", ".parquet")
    df.to_parquet(out, index=False)

    mf = Manifest(stage=stage, approach=approach, params=params or {})
    for name, p in (inputs or {}).items():
        mf.add_input(name, p)
    mf.add_output("full_frame", out)
    mf.note(f"{len(df)} candidates; "
            f"{int(df['rejected_at'].notna().sum())} stamped rejected and RETAINED")
    mf.write(d, filename=f"{out.stem}_manifest.json")
    log.info("[%s] wrote full frame: %s (%d rows)", approach, out, len(df))
    return out


def write_top10(df: pd.DataFrame, *, approach: str, experiment: str,
                stage: str = "top10", params: dict | None = None) -> Path:
    """Validate and write the GUI hand-off frame as CSV."""
    report = schema.validate_top10(df, approach)
    for w in report.warnings:
        log.warning("[%s] %s", approach, w)
    report.raise_if_bad(f"{approach} top-10")

    d = approach_dir(approach, experiment)
    out = next_version(d, f"D{approach[-1]}_top10", ".csv")
    df.to_csv(out, index=False)

    (Manifest(stage=stage, approach=approach, params=params or {})
     .add_output("top10", out)
     .note(f"{len(df)} candidates ranked on "
           f"{df['rank_metric_name'].iloc[0] if len(df) else 'n/a'} "
           f"({df['rank_metric_direction'].iloc[0] if len(df) else 'n/a'})")
     .write(d, filename=f"{out.stem}_manifest.json"))
    log.info("[%s] wrote top-10: %s", approach, out)
    return out


def latest(directory: Path, stem: str, suffix: str) -> Path | None:
    """Highest-numbered version of a versioned artifact, or None."""
    best, best_n = None, -1
    for p in Path(directory).glob(f"{stem}_*{suffix}"):
        m = _VERSIONED.match(p.stem)
        if m and m.group("stem") == stem:
            try:
                n = int(m.group("version"))
            except ValueError:
                continue
            if n > best_n:
                best, best_n = p, n
    return best


def read_top10(approach: str, experiment: str) -> pd.DataFrame:
    """Read an approach's latest top-10 hand-off, validating it on the way in.

    The GUI uses this. Validating on read as well as on write is deliberate: a
    file can be edited by hand between the two, and a malformed hand-off should
    fail at the point of use rather than render as a plausible-looking table.
    """
    d = approach_dir(approach, experiment)
    p = latest(d, f"D{approach[-1]}_top10", ".csv")
    if p is None:
        raise IOError_(f"no top-10 hand-off for {approach} under {d}")
    df = pd.read_csv(p)
    schema.validate_top10(df, approach).raise_if_bad(f"{p.name}")
    return df


def pool_top10(experiment: str, approaches: tuple[str, ...] = ("t1", "t2", "t3", "t4")
               ) -> pd.DataFrame:
    """Concatenate the available top-10 frames for the integration phase.

    Concatenation only — NO numeric join, and no cross-approach normalisation of
    rank metrics. The approaches rank on different quantities in different units
    and directions; the GUI presents them side by side and the human adjudicates
    (Rev 3 section 7). A missing approach is reported and skipped rather than
    blocking the pool.
    """
    frames = []
    for a in approaches:
        try:
            frames.append(read_top10(a, experiment))
        except IOError_ as exc:
            log.warning("pool: %s", exc)
    if not frames:
        raise IOError_(f"no approach hand-offs found for experiment {experiment!r}")
    pooled = pd.concat(frames, ignore_index=True)
    log.info("pooled %d candidates from %d approach(es): %s",
             len(pooled), len(frames), sorted(pooled["approach"].unique()))
    return pooled
