"""
Purpose: the GUI must show the SYNTHESIZABLE shortlist by default, not the raw one.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: integration.app.data
Output: pass/fail

WHY THIS EXISTS. `scripts/reshortlist_synthesizable.py` built `shortlist_synth`
-- a different SET of candidates, with rule failures dropped and the next-best
passers promoted into their slots. The GUI kept reading the raw `shortlist`
column and merely sorted failures to the bottom of it. Both lists existed, both
looked plausible, and the panel reported "N candidates now occupy the top slots
that failures held" while in fact promoting nobody: the replacements were never
in the frame it was reading.

That is this project's recurring defect in its purest form -- a value selected
by the name that was already there rather than by the one that answers the
question, failing silently because both are populated. The user caught it by
noticing failures still on screen.
"""

from __future__ import annotations

import pandas as pd
import pytest

from integration.app import data as D


def _frame(with_synth: bool = True) -> pd.DataFrame:
    """Five candidates; the raw top-3 contains one rule failure."""
    df = pd.DataFrame({
        "candidate_id": ["a", "b", "c", "d", "e"],
        "canonical_smiles": ["C"] * 5,
        "approach": ["t1"] * 5,
        "rejected_at": [pd.NA] * 5,
        "rank": [1, 2, 3, 4, 5],
        "shortlist": [True, True, True, False, False],
        "vina_affinity": [-9.0, -8.5, -8.0, -7.5, -7.0],
    })
    if with_synth:
        # 'c' failed a rule; 'd' was promoted in its place.
        df["shortlist_synth"] = [True, True, False, True, False]
        df["rank_synth"] = [1, 2, pd.NA, 3, pd.NA]
    return df


def test_default_prefers_the_synthesizable_column(monkeypatch):
    monkeypatch.setattr(D, "load_frame", lambda a: (_frame(), "D1_test.parquet"))
    s = D.shortlist("t1")
    ids = set(s["candidate_id"])
    assert ids == {"a", "b", "d"}, (
        f"default shortlist returned {sorted(ids)}; it must read "
        f"{D.SYNTH_COLUMN}, dropping the rule failure 'c' and promoting 'd'")
    assert (s["shortlist_column"] == D.SYNTH_COLUMN).all()


def test_the_raw_list_is_still_reachable(monkeypatch):
    monkeypatch.setattr(D, "load_frame", lambda a: (_frame(), "D1_test.parquet"))
    s = D.shortlist("t1", prefer_synth=False)
    assert set(s["candidate_id"]) == {"a", "b", "c"}
    assert (s["shortlist_column"] == D.BASE_COLUMN).all()


def test_falls_back_when_the_rebuild_is_absent(monkeypatch):
    """A frame predating the rebuild must render, not vanish."""
    monkeypatch.setattr(D, "load_frame",
                        lambda a: (_frame(with_synth=False), "D1_old.parquet"))
    s = D.shortlist("t1")
    assert set(s["candidate_id"]) == {"a", "b", "c"}
    assert (s["shortlist_column"] == D.BASE_COLUMN).all(), (
        "with no shortlist_synth the caller must be told it got the raw list, "
        "so the panel cannot claim a filter was applied")


def test_delta_reports_what_changed(monkeypatch):
    monkeypatch.setattr(D, "load_frame", lambda a: (_frame(), "D1_test.parquet"))
    d = D.shortlist_delta("t1")
    assert d == {"available": True, "promoted": 1, "dropped": 1, "kept": 2}


def test_delta_is_unavailable_without_the_rebuild(monkeypatch):
    monkeypatch.setattr(D, "load_frame",
                        lambda a: (_frame(with_synth=False), "D1_old.parquet"))
    assert D.shortlist_delta("t1")["available"] is False


def test_display_rank_uses_the_synth_ordering(monkeypatch):
    monkeypatch.setattr(D, "load_frame", lambda a: (_frame(), "D1_test.parquet"))
    s = D.shortlist("t1").sort_values("display_rank")
    assert list(s["candidate_id"]) == ["a", "b", "d"]
    # The raw rank survives as provenance, with a gap where 'c' was removed.
    assert list(s["rank"]) == [1, 2, 4]


@pytest.mark.parametrize("approach", sorted(D.APPROACHES))
def test_live_frames_do_not_regress(approach):
    """On the real data, the default list must contain no rule failure."""
    from shared import synthesizability as syn
    s = D.shortlist(approach)
    if s.empty or "shortlist_column" not in s.columns:
        pytest.skip(f"{approach}: no shortlist on disk")
    if s["shortlist_column"].iloc[0] != D.SYNTH_COLUMN:
        pytest.skip(f"{approach}: frame predates the rebuild")
    bad = [r.candidate_id for r in s.itertuples()
           if syn.violations(str(r.canonical_smiles))]
    assert not bad, (
        f"{approach}: the synthesizable shortlist still contains rule "
        f"failures {bad} — rerun scripts/reshortlist_synthesizable.py")
