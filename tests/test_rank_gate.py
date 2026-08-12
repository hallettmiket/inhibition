"""
Purpose: prove the rank gate counts poses, not a fraction of the cloud.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: none — synthetic mode tables with a known answer
Output: pass/fail

WHY THIS EXISTS (#65). The gate deciding which modes may hold a rank was
`consensus >= 0.05`. `consensus` is exactly `mode_size / n_poses`, so on the
3.0.0 run — where every cloud holds 500 poses — that rule was precisely
`n_poses_mode >= 25`. Nothing measured 25. It is what 0.05 of 500 comes to.

The defect is not the number, it is the DIVISION. A size floor written as a
fraction of the workload moves when the workload moves: double
`docking.n_runs` and the same 0.05 becomes "at least 50 poses", roughly halving
what is ranked, with nothing in any output announcing that the gate changed.
That is `how_this_project_breaks.md` disguise #3 — a constant sized against a
workload that has since grown — and the test below is written so it fails if
the gate ever goes back to depending on cloud size.

The checks:

  1. the gate admits on ABSOLUTE population, and is unchanged by cloud size
  2. the superseded fraction gate IS affected by cloud size — asserted so the
     difference between the two is demonstrated rather than claimed
  3. the two config keys that must agree do agree (`ranking.mode_gate.min_poses`
     and `sweep_rule.min_mode_poses`), by assertion rather than by aliasing
  4. the gate REFUSES a frame that cannot carry it instead of falling back
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import target_config as tc                # noqa: E402
import rank_v2 as rv                                  # noqa: E402


def _modes(sizes, n_poses):
    """One warhead class, one mode per entry in `sizes`, all equally scored."""
    return pd.DataFrame([
        {"ident": f"x_m{i}", "warhead_class": "chloroacetamide",
         "n_poses_mode": s, "n_poses": n_poses, "consensus": s / n_poses,
         "score": 10.0 - i}
        for i, s in enumerate(sizes)])


def test_gate_counts_poses_and_ignores_cloud_size():
    """A 12-pose mode is admitted whether the cloud held 500 poses or 5,000."""
    m = tc.rank_min_mode_poses()
    sizes = [m - 1, m, m + 50]
    for cloud in (500, 5000):
        out = rv.filter_and_rank(_modes(sizes, cloud), "score", "consensus",
                                 quota=1.0, floor=0.05, gate="mode_poses")
        out = out.sort_values("n_poses_mode")
        assert list(out.passes) == [False, True, True], (
            f"the gate changed with a cloud of {cloud} poses; it must not")
        assert set(out.rank_gate) == {"mode_poses"}
        assert set(out.rank_gate_min) == {m}


def test_the_superseded_fraction_gate_moves_with_the_cloud():
    """The defect itself, demonstrated — so the fix is measured against it.

    The SAME three modes, the SAME threshold, a different cloud size, and a
    different answer. Nothing in the output of either run says the rule moved.
    """
    sizes = [30, 60, 120]
    small = rv.filter_and_rank(_modes(sizes, 500), "score", "consensus",
                               quota=1.0, floor=0.05, gate="consensus_fraction")
    big = rv.filter_and_rank(_modes(sizes, 5000), "score", "consensus",
                             quota=1.0, floor=0.05, gate="consensus_fraction")
    assert list(small.sort_values("n_poses_mode").passes) == [True, True, True]
    assert list(big.sort_values("n_poses_mode").passes) == [False, False, False]


def test_ranking_and_sweep_floors_agree():
    """They are the same estimability question and must not drift apart silently.

    Kept as two config keys ON PURPOSE, so one could be changed without the
    other -- but a change to either should be a decision, and this fails until
    it is made in both places or this assertion is deliberately relaxed.
    """
    assert tc.rank_min_mode_poses() == tc.sweep_min_mode_poses()
    assert tc.rank_min_mode_poses() == int(tc.get("splitting.stage2.min_mode_size"))


def test_gate_refuses_a_frame_it_cannot_evaluate():
    """No silent fallback to the fraction on a frame with no per-mode population.

    A 2.1.0 aggregate has no `n_poses_mode`. Falling back would rank it by a
    different rule under the same filename, which is D0080's defect exactly.
    """
    df = _modes([12, 40], 500).drop(columns=["n_poses_mode"])
    with pytest.raises(SystemExit, match="consensus_fraction"):
        rv.filter_and_rank(df, "score", "consensus", quota=1.0, floor=0.05,
                           gate="mode_poses")


def test_unknown_gate_is_refused():
    with pytest.raises(ValueError, match="unknown rank gate"):
        rv.filter_and_rank(_modes([40], 500), "score", "consensus",
                           quota=1.0, floor=0.05, gate="whatever")
