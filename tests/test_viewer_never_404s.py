#!/usr/bin/env python3
"""
Purpose: no rail row can point the viewer at a page that does not exist.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17

@tt8804 sent a screenshot of the MD results page with the server's 404 rendered
inside the viewer pane: "Error code: 404 / File not found".

TWO FAULTS, ONE SYMPTOM. `sweep_combine` and `mdprio_combine` only ever LINK
`sweep_pages/<ident>.html`; the page is written by `sweep_report.py --ident`,
which nothing in the refresh loop ran -- so a mode swept after the last manual
build had no page (39 of 129 were missing). And `show()` fell back to
`<ident>.html` when a row carried no `data-src`, which for a queued molecule is
nothing at all. A row that cannot be viewed is the same failure as a page that
does not render, and neither builder checked its own links.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                # noqa: E402

COMBINED = rp.reports_dir() / "combined.html"
pytestmark = pytest.mark.skipif(
    not COMBINED.is_file(), reason="no combined.html built for this topic yet")


def _rows() -> list[tuple[str, str]]:
    """(ident, viewer target) for every rail row, mirroring show()'s rule:
    the row's data-src if it has one, else <ident>.html."""
    h = COMBINED.read_text()
    out = []
    for m in re.finditer(r"id='b_([^']+)'", h):
        ident = m.group(1)
        seg = h[max(0, m.start() - 600):m.start()]
        ds = re.findall(r"data-src='([^']*)'", seg)
        out.append((ident, ds[-1] if ds else f"{ident}.html"))
    return out


def test_there_are_rows_at_all():
    assert _rows(), "combined.html has no rail rows — the parse below is vacuous"


def test_every_row_points_at_a_file_that_exists():
    missing = [(i, t) for i, t in _rows()
               if not (rp.reports_dir() / t).is_file()]
    assert not missing, "rows whose viewer target is absent (these render a 404 "
    f"inside the pane): {missing}"


def test_a_row_may_only_fall_through_if_it_has_its_own_report():
    """Falling through to `<ident>.html` is correct for a finished 100 ns run --
    that page exists. It is the QUEUED rows that must carry a data-src, because
    for them the bare ident resolves to nothing. That asymmetry is the bug."""
    offenders = [i for i, t in _rows()
                 if t == f"{i}.html" and not (rp.reports_dir() / t).is_file()]
    assert not offenders, (
        f"rows with no data-src and no report page of their own: {offenders}")


def test_every_swept_mode_has_a_sweep_page():
    """sweep_report is what writes these, and it must actually be run."""
    import csv
    import glob
    import os
    import time
    pages = rp.reports_dir() / "sweep_pages"
    if not pages.is_dir():
        pytest.skip("no sweep pages for this topic yet")
    # A mode that finished in the last quarter hour has not had a refresh cycle
    # yet, and failing on that would make this test track the clock rather than
    # the defect. The defect was 39 modes, some of them days old.
    GRACE_S = 900
    now = time.time()
    missing = set()
    for f in glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")):
        if now - os.path.getmtime(f) < GRACE_S:
            continue
        with open(f) as fh:
            for r in csv.DictReader(fh):
                i = r.get("ident")
                if r.get("status") == "ok" and i and not (pages / f"{i}.html").is_file():
                    missing.add(i)
    assert not missing, (
        f"{len(missing)} swept mode(s) have no sweep_pages/<ident>.html; "
        f"run scripts/sweep_report.py --ident for each: {sorted(missing)[:5]}")
