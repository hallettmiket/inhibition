"""
Purpose: the fast ROC-AUC must be the SAME statistic as the definition it replaced.
Author: @tt8804 (with Claude Code)
Date: 2026-09-09

WHY THIS EXISTS. On 2026-09-09 `roc_auc` was changed from the pairwise
definition -- O(actives x decoys) -- to a sort-and-binary-search form, because
the gate's population grew from ~1,700 property-matched decoys to 361,354
ASSAYED inactives and the pairwise form stopped finishing (24.6 billion
comparisons across a 2,000-replicate bootstrap).

A speedup that quietly changes the number is worse than a slow gate, because
the verdict would still be a plausible float in the right range. So the old
implementation is kept as `_roc_auc_pairwise` and this asserts the new one
agrees with it EXACTLY, including the tie handling, which is the part most
likely to differ: the definition weights a tie 0.5, and a searchsorted form
gets that right only if it takes the gap between the "left" and "right"
insertion points.

This is the habit `how_this_project_breaks.md` recommends and scores as the
route that caught #33 and #34: run the new thing on something whose answer is
already written down. Here the written-down answer is the previous code.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import enrichment_gate as eg      # noqa: E402


def _population(rng, n_pos, n_neg, *, decimals=None):
    """Scores for a labelled population; `decimals` forces ties by rounding."""
    scores = rng.normal(size=n_pos + n_neg)
    if decimals is not None:
        scores = np.round(scores, decimals)
    labels = [1] * n_pos + [0] * n_neg
    return list(scores), labels


@pytest.mark.parametrize("higher_is_better", [True, False])
@pytest.mark.parametrize("n_pos,n_neg", [(3, 7), (10, 200), (34, 500)])
def test_fast_auc_equals_the_pairwise_definition(n_pos, n_neg, higher_is_better):
    rng = np.random.default_rng(0)
    scores, labels = _population(rng, n_pos, n_neg)
    fast = eg.roc_auc(scores, labels, higher_is_better=higher_is_better)
    slow = eg._roc_auc_pairwise(scores, labels, higher_is_better=higher_is_better)
    assert fast == pytest.approx(slow, abs=1e-12)


@pytest.mark.parametrize("higher_is_better", [True, False])
def test_ties_are_weighted_a_half_exactly_as_before(higher_is_better):
    """THE PART MOST LIKELY TO DIVERGE.

    Rounding to one decimal over a few hundred draws guarantees many exact
    ties between actives and decoys, which is also the realistic case: a
    Tanimoto over a small model set takes few distinct values.
    """
    rng = np.random.default_rng(7)
    scores, labels = _population(rng, 34, 800, decimals=1)
    assert len(set(scores)) < len(scores), "fixture failed to produce ties"
    fast = eg.roc_auc(scores, labels, higher_is_better=higher_is_better)
    slow = eg._roc_auc_pairwise(scores, labels, higher_is_better=higher_is_better)
    assert fast == pytest.approx(slow, abs=1e-12)


def test_all_tied_scores_give_exactly_one_half():
    """A scorer that separates nothing must read 0.5, not 0.0 or 1.0."""
    scores = [1.0] * 50
    labels = [1] * 10 + [0] * 40
    assert eg.roc_auc(scores, labels, higher_is_better=True) == 0.5
    assert eg.roc_auc(scores, labels, higher_is_better=False) == 0.5


def test_perfect_and_inverted_separation():
    scores = [1.0, 2.0, 3.0] + [-1.0, -2.0, -3.0]
    labels = [1, 1, 1, 0, 0, 0]
    assert eg.roc_auc(scores, labels, higher_is_better=True) == 1.0
    assert eg.roc_auc(scores, labels, higher_is_better=False) == 0.0


def test_direction_flip_is_the_complement():
    rng = np.random.default_rng(3)
    scores, labels = _population(rng, 20, 300)
    up = eg.roc_auc(scores, labels, higher_is_better=True)
    down = eg.roc_auc(scores, labels, higher_is_better=False)
    assert up + down == pytest.approx(1.0, abs=1e-12)


def test_auc_still_refuses_a_population_with_no_decoys():
    with pytest.raises(eg.EnrichmentGateError):
        eg.roc_auc([1.0, 2.0], [1, 1], higher_is_better=True)
    with pytest.raises(eg.EnrichmentGateError):
        eg.roc_auc([1.0, 2.0], [0, 0], higher_is_better=True)


def _bootstrap_ci_old(scores, labels, *, higher_is_better, n_boot, seed=42):
    """The pre-2026-09-09 bootstrap, verbatim, for the agreement check."""
    import random
    rng = random.Random(seed)
    ai = [i for i, l in enumerate(labels) if l == 1]
    di = [i for i, l in enumerate(labels) if l == 0]
    vals = []
    for _ in range(n_boot):
        sa = [rng.choice(ai) for _ in ai]
        sd = [rng.choice(di) for _ in di]
        vals.append(eg._roc_auc_pairwise(
            [scores[i] for i in sa + sd],
            [1] * len(sa) + [0] * len(sd), higher_is_better=higher_is_better))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1])


def test_bootstrap_agrees_with_the_old_sampler_within_monte_carlo_error():
    """The ESTIMATOR is preserved; the RNG stream deliberately is not.

    Vectorising the resampling changed which indices are drawn, so the
    endpoints cannot match to the digit and asserting that they do would be
    asserting the wrong thing. What must hold is that both samplers estimate
    the same interval. Checked on a population small enough for the old
    implementation to finish.
    """
    rng = np.random.default_rng(11)
    scores, labels = _population(rng, 15, 400)
    new = eg.bootstrap_ci(scores, labels, higher_is_better=True, n_boot=400)
    old = _bootstrap_ci_old(scores, labels, higher_is_better=True, n_boot=400)
    assert new[0] == pytest.approx(old[0], abs=0.05)
    assert new[1] == pytest.approx(old[1], abs=0.05)
    assert new[0] < new[1]


def test_bootstrap_brackets_the_point_estimate():
    rng = np.random.default_rng(5)
    scores, labels = _population(rng, 30, 900)
    point = eg.roc_auc(scores, labels, higher_is_better=True)
    lo, hi = eg.bootstrap_ci(scores, labels, higher_is_better=True, n_boot=500)
    assert lo <= point <= hi


def test_fast_path_is_actually_fast_at_the_measured_gate_size():
    """The regression this whole change exists to prevent.

    34 actives against 361,354 decoys is the real gate population. The pairwise
    form needs 12.3M comparisons for ONE call; this asserts the shipped path
    does a single call in well under a second, so a 2,000-replicate bootstrap
    is minutes rather than days.
    """
    import time
    rng = np.random.default_rng(1)
    scores, labels = _population(rng, 34, 361_354)
    t = time.perf_counter()
    auc = eg.roc_auc(scores, labels, higher_is_better=True)
    dt = time.perf_counter() - t
    assert 0.0 <= auc <= 1.0
    assert dt < 1.0, f"one AUC at gate scale took {dt:.2f}s; the bootstrap "\
                     f"calls it 2,000 times"
