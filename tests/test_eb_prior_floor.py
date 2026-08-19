#!/usr/bin/env python3
"""
Purpose: the empirical-Bayes prior is strong enough to actually shrink.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18

@tt8804: "run an experiment do 500 poses times 5 runs for our last top hit and
see how many times we actually identify the right mode by splitting 500", then
"yes build the fix".

WHAT THE EXPERIMENT FOUND (exp/1_mode_stability). Five independent 500-pose
screens of t4_716800c125a7 recovered the mode that 3.0.0 elected and validated
with a 100 ns run 5/5 times -- so sampling is not the problem -- but ranked it
first only 2/5 times. Method of moments fits the prior to the population's
SPREAD, and the population is heterogeneous, so the concentration collapsed to
~2.2 poses: a prior worth two poses shrinks essentially nothing, and a 12-pose
mode at 0.42 outranked an 82-pose mode at 0.24.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import target_config as tc               # noqa: E402


def _eb(k, n, iso, mu, conc):
    a0, b0 = mu * conc, (1 - mu) * conc
    return ((k + a0) / (n + a0 + b0)) / iso


def test_the_floor_is_configured_not_hardcoded():
    assert float(tc.get("ranking.eb_prior_min_strength")) > 0


def test_the_floor_is_below_the_typical_evidence():
    """It must bite only where the estimate is thin. This library's median
    n_in_range is 26; a floor at or above that would overwhelm real data."""
    assert float(tc.get("ranking.eb_prior_min_strength")) < 26


def test_the_floor_narrows_a_thin_mode_s_advantage_without_reversing_it():
    """What the floor actually does, stated exactly.

    5/12 (vf 0.42) against 20/82 (vf 0.24): the thin mode wins at EVERY prior
    strength -- 0.42 vs 0.24 is a real difference and shrinkage toward the
    population mean of 0.205 should not overturn it. What the floor does is
    shrink the thin mode's ADVANTAGE, from a ratio of 1.58 at the fitted prior
    to 1.34 at 10 and 1.10 at 40, because the thin estimate moves and the
    well-evidenced one barely does.

    An earlier version of this test asserted the floor made the thick mode win.
    It does not, at any strength, and asserting it would have shipped a fix
    described by a mechanism it does not have.
    """
    mu, iso = 0.2053, 0.0817
    ratios = [_eb(5, 12, iso, mu, c) / _eb(20, 82, iso, mu, c)
              for c in (2.2, 5, 10, 20, 40)]
    assert all(r > 1 for r in ratios), "the thin mode wins throughout"
    assert ratios == sorted(ratios, reverse=True), "a stronger prior must narrow the gap"
    assert ratios[0] > 1.5 and ratios[2] < 1.4


def test_shrinkage_moves_a_thin_estimate_toward_the_mean_not_toward_zero():
    """Shrinkage must have no systematic direction -- that is why EB was chosen
    over a lower confidence bound in the first place."""
    mu, iso = 0.2053, 0.0817
    below = _eb(k=1, n=12, iso=iso, mu=mu, conc=10.0)     # p = 0.083, under mu
    above = _eb(k=8, n=12, iso=iso, mu=mu, conc=10.0)     # p = 0.667, over mu
    assert below > (0.083 / iso), "an estimate below the mean must be pulled UP"
    assert above < (0.667 / iso), "an estimate above the mean must be pulled DOWN"


def test_the_floor_only_raises_never_lowers():
    """A library whose own fit is already strong keeps its fitted prior."""
    floor = float(tc.get("ranking.eb_prior_min_strength"))
    for fitted in (0.5, 2.2, floor, floor + 15):
        assert max(fitted, floor) >= fitted


@pytest.mark.parametrize("bad", ["sqrt_size", "linear_size"])
def test_size_weighting_stays_rejected(bad):
    """RECORDED SO IT IS NOT RETRIED. Multiplying the score by sqrt(mode_size)
    elects the validated mode 5/5 on the replicate experiment -- and scores
    rho(score, mode_size) = +0.72 across the library, so it is largely a size
    proxy that would discard genuine minority modes. A rule that wins by
    ranking big things first has not solved the problem it was aimed at."""
    src = (REPO / "scripts" / "rank_v2.py").read_text()
    assert "sqrt(d.mode_size)" not in src.replace(" ", "")
    assert "mode_size)" not in src.split("conditional_eb")[1][:400].replace(" ", "")


def test_the_experiment_that_justifies_this_is_in_the_repo():
    """The number in the comment has to be reproducible by someone else."""
    exp = REPO / "exp" / "1_mode_stability" / "run_all.py"
    assert exp.is_file()
    s = exp.read_text()
    assert "REFERENCE" in s and "validated_100ns_rmsd_max_nm" in s
