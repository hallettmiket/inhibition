"""
Purpose: Load and validate decision records; the feed for the GUI's Decisions pane.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: decisions/*.md — YAML frontmatter plus prose
Output: validated records, queryable by approach/status; an integrity check

Files are the source of truth. This module is how the GUI (and anyone else)
reads them without re-implementing the format, and how CI catches a malformed
or dangling record before it reaches a reader.

The validation matters more than it looks. A decision log people stop trusting
is worse than no decision log, and the fastest way to lose trust is a record
that says 'supersedes D0003' when D0003 says nothing about being superseded.
"""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DIR = _REPO_ROOT / "decisions"

VALID_STATUS = {"proposed", "accepted", "superseded", "rejected"}
VALID_APPROACH = {"shared", "t1", "t2", "t3", "t4", "integration"}
VALID_ORIGIN = {"spec", "adversary", "implementation", "user"}
REQUIRED_FIELDS = {"id", "title", "date", "status", "approach"}

_FRONTMATTER = re.compile(r"\A---\n(.*?)\n---\n(.*)\Z", re.DOTALL)


class DecisionError(RuntimeError):
    """A decision record is malformed, duplicated, or dangling."""


@dataclass
class Decision:
    """One decision record."""

    id: str
    title: str
    date: str
    status: str
    approach: str
    path: Path
    body: str = ""
    decided_by: str = ""
    origin: str = "implementation"
    supersedes: list[str] = field(default_factory=list)
    superseded_by: str | None = None
    affects: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    runbook: str | None = None

    def section(self, name: str) -> str:
        """Extract one '## Name' section from the prose body, or ''."""
        m = re.search(rf"^## {re.escape(name)}\n(.*?)(?=^## |\Z)",
                      self.body, re.DOTALL | re.MULTILINE)
        return m.group(1).strip() if m else ""

    def to_dict(self) -> dict:
        """Flat dict for the GUI — frontmatter plus the three prose sections."""
        return {
            "id": self.id, "title": self.title, "date": self.date,
            "status": self.status, "approach": self.approach,
            "decided_by": self.decided_by, "origin": self.origin,
            "supersedes": self.supersedes, "superseded_by": self.superseded_by,
            "affects": self.affects, "evidence": self.evidence,
            "runbook": self.runbook, "path": str(self.path),
            "context": self.section("Context"),
            "decision": self.section("Decision"),
            "consequences": self.section("Consequences"),
        }


def _parse(path: Path) -> Decision:
    text = path.read_text(encoding="utf-8")
    m = _FRONTMATTER.match(text)
    if not m:
        raise DecisionError(f"{path.name}: no YAML frontmatter block")
    meta = yaml.safe_load(m.group(1)) or {}
    missing = REQUIRED_FIELDS - set(meta)
    if missing:
        raise DecisionError(f"{path.name}: missing field(s) {sorted(missing)}")
    if meta["status"] not in VALID_STATUS:
        raise DecisionError(f"{path.name}: status {meta['status']!r} not in {sorted(VALID_STATUS)}")
    if meta["approach"] not in VALID_APPROACH:
        raise DecisionError(f"{path.name}: approach {meta['approach']!r} not in {sorted(VALID_APPROACH)}")
    if meta.get("origin", "implementation") not in VALID_ORIGIN:
        raise DecisionError(f"{path.name}: origin {meta['origin']!r} not in {sorted(VALID_ORIGIN)}")
    return Decision(
        id=str(meta["id"]), title=str(meta["title"]), date=str(meta["date"]),
        status=meta["status"], approach=meta["approach"], path=path,
        body=m.group(2), decided_by=str(meta.get("decided_by", "")),
        origin=meta.get("origin", "implementation"),
        supersedes=list(meta.get("supersedes") or []),
        superseded_by=meta.get("superseded_by"),
        affects=list(meta.get("affects") or []),
        evidence=[str(e) for e in (meta.get("evidence") or [])],
        runbook=meta.get("runbook"),
    )


def load(directory: Path | str | None = None) -> list[Decision]:
    """Load every decision record, sorted by id.

    Raises
    ------
    DecisionError
        On a malformed record, a duplicate id, a dangling ``supersedes``
        reference, or an asymmetric supersede link.
    """
    d = Path(directory) if directory else DEFAULT_DIR
    if not d.is_dir():
        raise DecisionError(f"decisions directory not found: {d}")

    records = [_parse(p) for p in sorted(d.glob("D*.md"))]

    seen: dict[str, Path] = {}
    for r in records:
        if r.id in seen:
            raise DecisionError(f"duplicate id {r.id}: {seen[r.id].name} and {r.path.name}")
        seen[r.id] = r.path

    by_id = {r.id: r for r in records}
    for r in records:
        for sid in r.supersedes:
            target = by_id.get(sid)
            if target is None:
                raise DecisionError(f"{r.id} supersedes unknown id {sid}")
            # Asymmetry is the failure that quietly rots a decision log: the new
            # record claims to replace the old one, the old one still reads as
            # current, and whoever finds the old one first believes it.
            if target.status != "superseded" or target.superseded_by != r.id:
                raise DecisionError(
                    f"{r.id} supersedes {sid}, but {sid} is status={target.status!r} "
                    f"superseded_by={target.superseded_by!r}; set them to "
                    f"'superseded' and {r.id!r}."
                )
    return records


def by_approach(approach: str, directory: Path | str | None = None) -> list[Decision]:
    """Records governing one approach, plus the shared ones it inherits."""
    return [r for r in load(directory) if r.approach in (approach, "shared")]


def affecting(path_fragment: str, directory: Path | str | None = None) -> list[Decision]:
    """Records whose ``affects`` mentions a path — 'why is this file like this?'"""
    return [r for r in load(directory)
            if any(path_fragment in a for a in r.affects)]


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate and query decision records.")
    ap.add_argument("action", choices=["check", "list"])
    ap.add_argument("--approach", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    records = by_approach(args.approach) if args.approach else load()
    if args.action == "check":
        active = [r for r in records if r.status == "accepted"]
        print(f"OK — {len(records)} record(s), {len(active)} accepted, no duplicate "
              "ids, no dangling or asymmetric supersede links")
        no_ev = [r.id for r in records if not r.evidence]
        if no_ev:
            print(f"note: {len(no_ev)} record(s) carry no evidence: {', '.join(no_ev)}")
        return
    for r in records:
        mark = {"accepted": " ", "superseded": "~", "rejected": "x", "proposed": "?"}[r.status]
        print(f"[{mark}] {r.id}  {r.approach:11s} {r.title}")
        if r.evidence:
            print(f"       evidence: {r.evidence[0][:88]}")


if __name__ == "__main__":
    main()
