"""
Purpose: Tests for T_4 step 8 — within-class ranking and the quota shortlist.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: synthetic candidate frames
Output: pytest pass/fail

The central claim of step 8 is that ranking within warhead class produces a
different — and more useful — shortlist than a global sort. `test_global_sort_
would_starve_a_class` states that claim as an executable assertion rather than a
comment, using a frame built so that one class's raw affinities dominate.

The other tests pin the properties that are easy to invert by accident: which
direction is "better", that undocked rows stay null rather than sorting to the
bottom, and that flags survive into the shortlist reason.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

_spec = importlib.util.spec_from_file_location(
    "t4_rank", REPO / "approaches" / "t4_combinatorial" / "04_rank_within_class.py")
rank_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rank_mod)


def _frame() -> pd.DataFrame:
    """Three classes with deliberately unequal sizes and affinity ranges.

    Class C has the best raw affinities but only 4 successful docks; class B is
    small and mediocre; class A is large and middling. A global sort therefore
    excludes B entirely, which is the outcome the design exists to prevent.
    """
    rows = []
    for i in range(30):
        rows.append(dict(warhead_class="A", affinity_kcal=-9.0 + i * 0.1, HAC=30,
                         reactivity_flag="IN_WINDOW", candidate_id=f"a{i}"))
    for i in range(5):
        rows.append(dict(warhead_class="B", affinity_kcal=-7.0 + i * 0.1, HAC=20,
                         reactivity_flag="IN_WINDOW", candidate_id=f"b{i}"))
    for i in range(4):
        rows.append(dict(warhead_class="C", affinity_kcal=-11.0 + i * 0.1, HAC=44,
                         reactivity_flag="OUTSIDE_WINDOW", candidate_id=f"c{i}"))
    for i in range(3):          # enumerated, never docked
        rows.append(dict(warhead_class="A", affinity_kcal=None, HAC=30,
                         reactivity_flag="IN_WINDOW", candidate_id=f"x{i}"))
    return pd.DataFrame(rows)


@pytest.fixture
def ranked():
    return rank_mod.rank_within_class(_frame(), min_docked=20)


def test_lower_affinity_ranks_first(ranked):
    """affinity_kcal is lower-is-better (D0015); rank 1 must be the minimum."""
    a = ranked[ranked.warhead_class == "A"].dropna(subset=["class_rank"])
    assert a.loc[a.class_rank == 1, "affinity_kcal"].iloc[0] == a.affinity_kcal.min()


def test_percentile_spans_the_class(ranked):
    assert ranked.loc[ranked.candidate_id == "a0", "class_percentile"].iloc[0] == 100.0
    assert ranked.loc[ranked.candidate_id == "a29", "class_percentile"].iloc[0] == 0.0


def test_undocked_rows_keep_null_rank(ranked):
    """"Did not dock" and "docked badly" are different facts."""
    x = ranked[ranked.candidate_id.str.startswith("x")]
    assert x.class_rank.isna().all()
    assert x.class_percentile.isna().all()


def test_small_classes_are_flagged_not_selective(ranked):
    assert bool(ranked.loc[ranked.candidate_id == "a0", "rank_is_selective"].iloc[0])
    assert not bool(ranked.loc[ranked.candidate_id == "b0", "rank_is_selective"].iloc[0])


def test_ligand_efficiency_is_affinity_per_heavy_atom(ranked):
    r = ranked.loc[ranked.candidate_id == "a0"].iloc[0]
    assert abs(r.ligand_efficiency - (-r.affinity_kcal / r.HAC)) < 1e-6


def test_every_class_gets_its_quota(ranked):
    s = rank_mod.build_shortlist(ranked, quota=3, include_flagged=True)
    short = s[s.shortlist]
    assert len(short) == 9
    assert dict(short.warhead_class.value_counts()) == {"A": 3, "B": 3, "C": 3}


def test_global_sort_would_starve_a_class(ranked):
    """The design claim, as a test rather than a comment."""
    df = _frame().dropna(subset=["affinity_kcal"])
    global_top = set(df.nsmallest(9, "affinity_kcal").warhead_class)
    s = rank_mod.build_shortlist(ranked, quota=3, include_flagged=True)
    within_top = set(s[s.shortlist].warhead_class)
    assert "B" not in global_top, "fixture no longer demonstrates the failure mode"
    assert within_top == {"A", "B", "C"}


def test_flags_travel_into_the_shortlist_reason(ranked):
    s = rank_mod.build_shortlist(ranked, quota=3, include_flagged=True)
    short = s[s.shortlist]
    assert "OUTSIDE_WINDOW" in str(short[short.warhead_class == "C"].shortlist_reason.iloc[0])
    assert "not selective" in str(short[short.warhead_class == "B"].shortlist_reason.iloc[0])
    assert s.loc[~s.shortlist, "shortlist_reason"].isna().all()


def test_flagged_classes_can_be_excluded(ranked):
    s = rank_mod.build_shortlist(ranked, quota=3, include_flagged=False)
    assert set(s[s.shortlist].warhead_class) == {"A", "B"}


def test_quota_larger_than_class_takes_the_whole_class(ranked):
    s = rank_mod.build_shortlist(ranked, quota=50, include_flagged=True)
    short = s[s.shortlist]
    assert (short.warhead_class == "C").sum() == 4      # class C only has 4 docked
    assert short.affinity_kcal.notna().all()            # never an undocked row
