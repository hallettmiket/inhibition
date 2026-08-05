"""
Purpose: the degree-2 sampler returns exactly `target` items, uniformly, with no population estimate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: the reservoir draw used by scripts/sample_t2_degree2.py
Output: pass/fail

WHY THIS EXISTS. The sampler previously kept each child with
`p = target / (n_parents * MEAN_CHILDREN_PER_PARENT)` -- a LINEAR extrapolation
from a 30-parent probe, applied to a deduplicated union that grows SUBLINEARLY
because parents' children overlap. Measured on the ATRA run: estimated
7,800,890 against a realised 4,063,427, so p was half what it should have been
and the frame kept 15,653 against a target of 30,000.

The sample was still unbiased -- every molecule had the same p -- so nothing
about the CHEMISTRY was wrong. Only the SIZE was, and nothing flagged it: the
docstring claimed the realised size "varies slightly around target", and a
48% shortfall reads as slightly if nobody checks. Worse, the constant was
measured on ATRA and would have been applied unchanged to four seeds whose
degree-1 neighbourhoods differ by up to 9x.

Reservoir sampling needs no estimate at all, which deletes the class rather
than retuning the constant. These tests pin the two properties that matter:
the size is EXACT, and the draw is UNIFORM over everything seen.
"""

from __future__ import annotations

import random
from collections import Counter


def reservoir(stream, target: int, rng: random.Random) -> list:
    """Algorithm R — the draw as implemented in the sampler's main loop."""
    kept: list = []
    n_eligible = 0
    for item in stream:
        n_eligible += 1
        if len(kept) < target:
            kept.append(item)
        else:
            j = rng.randrange(n_eligible)
            if j < target:
                kept[j] = item
    return kept


def test_the_sample_is_exactly_the_target_size():
    """The property the Bernoulli draw could not give.

    ATRA asked for 30,000 and got 15,653 because the population was guessed.
    """
    rng = random.Random(20260804)
    # Every population here EXCEEDS the target; the smaller-than-target case
    # is a different contract and has its own test below.
    for pop in (5_001, 10_000, 250_000):
        assert len(reservoir(range(pop), 5_000, rng)) == 5_000


def test_a_population_smaller_than_the_target_returns_all_of_it():
    """The honest answer, not a scaled-down one."""
    rng = random.Random(20260804)
    got = reservoir(range(120), 5_000, rng)
    assert len(got) == 120
    assert sorted(got) == list(range(120))


def test_the_draw_is_uniform_over_the_whole_stream():
    """Every eligible item equally likely — including the ones seen LAST.

    The failure this catches is a draw that fills up early and then rejects
    everything after, which would sample the first parents' children and none
    of the rest. That is exactly the `frontier_cap` truncation-in-parent-order
    the sampler exists to avoid, and it would look like a valid sample.
    """
    pop, target, trials = 200, 20, 4000
    counts = Counter()
    rng = random.Random(20260804)
    for _ in range(trials):
        counts.update(reservoir(range(pop), target, rng))

    assert len(counts) == pop, "some items were never sampled at all"
    expected = trials * target / pop
    worst = max(abs(c - expected) / expected for c in counts.values())
    assert worst < 0.25, (
        f"selection frequency deviates by {worst:.0%} from uniform; the draw "
        "favours part of the stream")


def test_the_tail_of_the_stream_is_not_starved():
    """Sharper version of the above, aimed at the specific failure mode."""
    pop, target, trials = 200, 20, 4000
    rng = random.Random(20260804)
    first_decile = last_decile = 0
    for _ in range(trials):
        got = reservoir(range(pop), target, rng)
        first_decile += sum(1 for x in got if x < pop // 10)
        last_decile += sum(1 for x in got if x >= pop - pop // 10)
    ratio = last_decile / first_decile
    assert 0.8 < ratio < 1.25, (
        f"last decile sampled {ratio:.2f}x the first — the draw is "
        "position-biased, which is the defect a uniform sample exists to avoid")
