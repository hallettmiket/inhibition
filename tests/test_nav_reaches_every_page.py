#!/usr/bin/env python3
"""
Purpose: every page the GUI builds is reachable from every other page.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17

@tt8804: "lets see the how it works slides. not working with the link".

The deck was built, served with the right bytes, and rendered correctly. It had
no link outside the MD results toolbar, so from Home, Ranking or Sweep there was
no way in -- and a page nobody can navigate to is indistinguishable from one
that is broken. Nothing tested reachability, so nothing caught it.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gui_shell as gs                # noqa: E402


def test_the_explainer_is_linked_from_the_nav():
    assert 'href="pipeline.html"' in gs.nav("index.html")


def test_the_explainer_is_reachable_from_every_step():
    """Not just from MD results, which is where the only link used to live."""
    for href, _, _ in gs.STEPS:
        assert 'href="pipeline.html"' in gs.nav(href), f"unreachable from {href}"


def test_the_explainer_is_not_a_funnel_step():
    """The stepper draws a chevron between consecutive items, so an explainer
    inside STEPS reads as the pipeline's last stage."""
    assert "pipeline.html" not in [s[0] for s in gs.STEPS]
    assert gs.ASIDE[0] == "pipeline.html"


def test_the_aside_takes_no_chevron():
    assert "#steps a.aside::before{display:none}" in gs.CSS


def test_the_explainer_highlights_when_it_is_the_current_page():
    h = gs.nav("pipeline.html")
    assert re.search(r'href="pipeline\.html" class="aside on"', h)


def test_every_nav_target_is_a_page_some_builder_writes():
    """A nav entry pointing at a filename nothing produces is a dead link the
    moment the topic changes."""
    built = set()
    for p in (REPO / "scripts").glob("*.py"):
        src = p.read_text()
        built |= set(re.findall(r'["\'](\w+\.html)["\']', src))
    for href, label, _ in list(gs.STEPS) + [gs.ASIDE]:
        assert href in built, f"nav links {href} ({label}) but no builder writes it"
