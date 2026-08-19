#!/usr/bin/env python3
"""
Purpose: the three-tier production verdict has one home, and controls are
         run-scoped like every other topic.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17

@tt8804: "so there should be held in pocket, held but not optimal and below
max .35 is optimal", and "get rid of the yellow controls, they dont belong in
this version".
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import residence_tier as rt          # noqa: E402
from shared import run_paths as rp               # noqa: E402
from shared import target_config as tc           # noqa: E402


def test_four_tiers_and_nothing_else():
    assert [k for k, _, _ in rt.TIERS] == ["optimal", "held", "unstable", "left"]


def test_optimal_bar_is_the_production_one_not_the_sweep_bar():
    """The production bar is its own config key.

    It was `md_survivor_rmsd_nm()` -- the 8 ns triage bar -- and sharing one
    number across both timescales made the tier unreachable by construction: a
    ligand explores more in 100 ns than in 8, and 13 finished runs bottomed out
    at 0.410 nm, so 'optimal' stayed empty while nine runs never left the
    pocket. @tt8804 set production to 0.45; the sweep bar stays 0.35.
    """
    assert rt.optimal_nm() == float(tc.md_production_optimal_rmsd_nm())
    assert rt.optimal_nm() != float(tc.md_survivor_rmsd_nm())


def test_the_two_bars_have_not_been_re_merged():
    assert float(tc.md_survivor_rmsd_nm()) == 0.35
    assert float(tc.md_production_optimal_rmsd_nm()) == 0.45


@pytest.mark.parametrize("mx,dissoc,res,want", [
    (0.10, False, 1.00, "optimal"),
    (0.449, False, 1.00, "optimal"),
    (0.45, False, 1.00, "held"),      # the bar is strict: 0.45 is NOT optimal
    (0.466, False, 1.00, "held"),
    (1.140, False, 1.00, "held"),     # brushed past 1.0 nm but came back
    # Came back and stayed, so not dissociated -- but 23% of the run was spent
    # outside the pocket. It is not the same object as a run that never moved.
    (5.765, False, 0.769, "unstable"),
    (0.40, False, 0.80, "unstable"),  # tight max cannot rescue poor residence
    (2.140, True, 0.558, "left"),
    (0.20, True, 1.00, "left"),       # dissociation outranks everything
])
def test_tiers(mx, dissoc, res, want):
    assert rt.tier(mx, dissoc, res) == want


def test_unmeasured_raises_rather_than_passing():
    """An unmeasured value entering a table as a passing one is this project's
    recurring defect -- 'value taken by label not identity'."""
    with pytest.raises(ValueError):
        rt.tier(None, False, 1.0)
    with pytest.raises(ValueError):
        rt.tier(0.2, None, 1.0)
    # An unknown residence must not quietly pass as a clean hold.
    with pytest.raises(ValueError):
        rt.tier(0.2, False, None)


def test_the_two_green_tiers_are_told_apart_by_their_words():
    """optimal and held are both green on purpose -- a tighter reading of the
    same good outcome -- so the WORD has to carry the distinction."""
    greens = [k for k, _, c in rt.TIERS if c == "good"]
    assert greens == ["optimal", "held"]
    assert len({rt.label(k) for k in greens}) == 2


def test_residence_floor_is_configured_not_hardcoded():
    assert rt.residence_floor() == float(tc.md_held_residence_floor())


def test_every_tier_carries_a_label_not_just_a_colour():
    """report_theme's own note: these tables are read by people with red-green
    deficiency, so colour can only ever be the second signal."""
    for key, label, colour in rt.TIERS:
        assert label and label != colour
        assert colour in ("good", "warn", "bad")
        assert label in rt.badge(key)


def test_controls_topic_is_run_scoped():
    assert rp.controls_topic("nac_v9") == "crystal_controls_nac_v9"
    assert rp.controls_topic() .endswith(rp.topic())


def test_no_script_reads_an_unscoped_controls_directory():
    """The leak this closes: mdprio_combine read a flat `crystal_controls`
    topic, so nac_v5 displayed eight controls it never produced -- including
    xtal_6VAJ, built against the receptor 3IKD replaced."""
    bad = []
    for p in sorted((REPO / "scripts").glob("*.py")) + sorted((REPO / "shared").glob("*.py")):
        if p.name == "run_paths.py":
            continue
        src = p.read_text()
        for m in re.finditer(r'Topic\(\s*["\'][^"\']+["\']\s*,\s*["\']([^"\']+)["\']', src):
            topic = m.group(1)
            if "control" in topic and not topic.endswith("_"):
                bad.append(f"{p.name}: Topic(..., {topic!r}) is not run-scoped")
    assert not bad, "\n".join(bad)
