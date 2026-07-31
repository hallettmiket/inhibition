"""
Purpose: Every decision record must parse, so the GUI cannot be broken by one.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: decisions/*.md
Output: pass/fail

WHY THIS EXISTS. Two records were written with `origin: audit`, which is not in
the vocabulary. `shared.decisions.load()` correctly refused the whole set --
and nothing ran it until a user clicked the Decisions panel, which then threw a
DecisionError instead of rendering. The validator was right; what was missing
was anything exercising it before a person did.

Writing a record is the most common way to touch this repo without running any
code, so it is the change most likely to reach a user unexercised.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from shared import decisions as dec

DECISIONS = Path(__file__).resolve().parent.parent / "decisions"
FILES = sorted(DECISIONS.glob("D*.md"))


def test_there_are_records_to_check():
    assert FILES, "no decision records found — the glob or the path is wrong"


def test_the_whole_set_loads():
    """The GUI calls this; if it raises, every decisions panel is dead."""
    records = dec.load()
    assert len(records) == len(FILES)


@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name.split("-")[0])
def test_each_record_parses_individually(path):
    """Named per record, so a failure says WHICH one rather than 'the set'."""
    try:
        dec._parse(path)
    except dec.DecisionError as exc:
        pytest.fail(f"{path.name}: {exc}")


def test_every_origin_is_in_the_vocabulary():
    bad = []
    for p in FILES:
        for line in p.read_text().splitlines()[:30]:
            if line.startswith("origin:"):
                val = line.split(":", 1)[1].split("#")[0].strip()
                if val and val not in dec.VALID_ORIGIN:
                    bad.append((p.name, val))
                break
    assert not bad, (
        f"origin values outside {sorted(dec.VALID_ORIGIN)}: {bad}. Either use "
        "an existing value or widen the vocabulary deliberately — do not "
        "invent one per record.")


def test_every_status_is_in_the_vocabulary():
    valid = getattr(dec, "VALID_STATUS", None)
    if valid is None:
        pytest.skip("no status vocabulary declared")
    bad = []
    for p in FILES:
        for line in p.read_text().splitlines()[:30]:
            if line.startswith("status:"):
                val = line.split(":", 1)[1].split("#")[0].strip()
                if val and val not in valid:
                    bad.append((p.name, val))
                break
    assert not bad, f"status values outside {sorted(valid)}: {bad}"
