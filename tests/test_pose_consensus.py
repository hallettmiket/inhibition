"""
Purpose: pose consensus reads the reactive region, states its N, and cannot call a bimodal molecule confident.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: shared/pose_consensus.py
Output: pass/fail

EVERY FIXTURE HERE IS HAND-BUILT AND ITS ANSWER IS KNOWN BEFORE THE MODULE RUNS.
None of them is docked output. A docked pose set has no independent "correct"
consensus — whatever the module returns IS the answer — so a test written
against one would assert that the code does what the code does. Testing D0055's
generated numbers taught the same lesson from the other side.

The four cases the design has to survive, in the order they matter:

* **identical poses** — the ceiling. If perfect agreement does not read 1.0,
  nothing below it means anything.
* **maximally scattered poses** — the floor, and the case that must be
  distinguishable from "could not be measured".
* **two tight clusters** — THE interesting one. Each cluster is internally
  perfect, so any statistic computed WITHIN a mode (largest-cluster fraction,
  the medoid's own neighbourhood, the mean of a cluster) reports a molecule
  that has confidently done two incompatible things as confident. The test
  asserts both that the fixture really is bimodal and that the headline number
  falls below half, so it documents the difference rather than merely
  forbidding the wrong answer.
* **too few poses** — must RAISE. D0068's whole shape is a number that means
  something other than what it appears to; a molecule that produced two poses
  scoring 0.0 would be indistinguishable from one whose poses genuinely
  disagree, and the first says fix the docking while the second says drop the
  molecule.

Plus the two parameter guards the module exists to enforce: the top-N is stated
in every result, and results measured at different N refuse to be compared.

MUTATION-TESTED 2026-08-06. `DEFAULT_TOLERANCE_A` 1.0 -> 6.5 fails only
`test_the_default_tolerance_is_the_measured_placement_radius`; every functional
test passes its tolerance explicitly so that none of them silently depends on
the default. `MIN_POSES_FOR_CONSENSUS` 3 -> 2 fails only
`test_too_few_poses_raise_rather_than_scoring_zero`;
`test_a_top_n_below_the_minimum_is_refused` and
`test_zero_and_unmeasurable_are_different_answers` derive their arguments FROM
the constant and so stay true under the mutation, which is the point — a test
that fails for every mutation localises nothing.

    Worth recording from mutation B: with the floor at 2 the code does not just
    return a debatable number, it computes `np.mean` of an empty jackknife
    slice and reports `nan` with a RuntimeWarning. The floor is load-bearing,
    not decorative.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_consensus as pc          # noqa: E402

TOL = 1.0          # passed explicitly everywhere, never inherited from the module
ATOMS = (7, 8, 9)  # a plausible reactive_atom_smarts match: C, LG, neighbour


def _warhead(offset: float | tuple[float, float, float]) -> np.ndarray:
    """Three reactive atoms in a FIXED internal geometry, translated to `offset`.

    Internal geometry is constant across every fixture, so every RMSD in this
    file is a statement about WHERE the reactive region sits in the site and
    never about the ligand's own conformation. That is the D0062 distinction
    the module is built on.
    """
    base = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.4, 0.0]])
    shift = np.array([offset, 0.0, 0.0]) if np.isscalar(offset) else np.asarray(offset)
    return base + np.asarray(shift, dtype=float)


def _pose(energy: float, offset, atom_ids: tuple[int, ...] = ATOMS) -> pc.ReactivePose:
    return pc.ReactivePose(energy=energy, reactive_xyz=_warhead(offset),
                           atom_ids=atom_ids)


def _tight(n: int, centre: float, energy0: float = -9.0) -> list[pc.ReactivePose]:
    """`n` poses within 0.25 A of `centre` — one mode, comfortably inside TOL."""
    return [_pose(energy0 + 0.1 * j, centre + 0.05 * j) for j in range(n)]


# --- the ceiling ----------------------------------------------------------

def test_identical_poses_are_perfect_consensus():
    """If this is not 1.0, no number below it is interpretable."""
    poses = [_pose(-9.0 - 0.1 * j, 0.0) for j in range(5)]
    r = pc.consensus(poses, top_n=5, tolerance_a=TOL)

    assert r.agreement == pytest.approx(1.0)
    assert r.agreement_jackknife == (pytest.approx(1.0), pytest.approx(1.0))
    assert r.median_rmsd == pytest.approx(0.0)
    assert r.max_rmsd == pytest.approx(0.0)
    assert r.n_modes == 1
    assert r.dominant_mode_fraction == pytest.approx(1.0)
    assert not r.is_multimodal


# --- the floor ------------------------------------------------------------

def test_maximally_scattered_poses_have_no_consensus():
    """Every pair further apart than tolerance: agreement 0.0, every pose its own mode."""
    poses = [_pose(-9.0 + j, 8.0 * j) for j in range(5)]
    r = pc.consensus(poses, top_n=5, tolerance_a=TOL)

    assert r.agreement == pytest.approx(0.0)
    assert r.agreement_jackknife == (pytest.approx(0.0), pytest.approx(0.0))
    assert r.n_modes == 5
    assert r.dominant_mode_fraction == pytest.approx(1 / 5)
    assert r.median_rmsd > TOL


def test_zero_and_unmeasurable_are_different_answers():
    """The floor is reachable, so it must not also mean 'not measured'.

    Three scattered poses score exactly 0.0. Two poses RAISE. A caller reading
    a column of consensus values can therefore tell 'this molecule's poses
    disagree' from 'this molecule was never measured' — which is the whole
    reason the raise exists.
    """
    scattered = [_pose(-9.0 + j, 8.0 * j) for j in range(3)]
    assert pc.consensus(scattered, top_n=3, tolerance_a=TOL).agreement == 0.0

    # Derived from the constant rather than written as a literal 2, so this
    # test asserts the DISTINCTION and leaves the floor's VALUE to the one test
    # that names it.
    too_few = scattered[:pc.MIN_POSES_FOR_CONSENSUS - 1]
    with pytest.raises(pc.ConsensusError):
        pc.consensus(too_few, top_n=3, tolerance_a=TOL)


# --- the interesting one --------------------------------------------------

def test_two_tight_clusters_do_not_read_as_high_consensus():
    """A molecule that has confidently done two incompatible things.

    Each mode is internally perfect, so every within-mode statistic says
    'confident'. The headline counts pose PAIRS, and in a balanced bimodal set
    the cross-mode pairs are the majority (25 of 45 here), so it lands below
    half. The largest-cluster fraction is asserted alongside at 0.5 to record
    exactly what the rejected design would have reported.
    """
    poses = _tight(5, 0.0) + _tight(5, 6.0, energy0=-8.5)
    r = pc.consensus(poses, top_n=10, tolerance_a=TOL)

    # The fixture has to actually be bimodal or this test proves nothing.
    assert r.n_modes == 2, "fixture is not bimodal; the test would be vacuous"
    assert r.dominant_mode_fraction == pytest.approx(0.5)

    within = 2 * (5 * 4 // 2)
    assert r.agreement == pytest.approx(within / 45)
    assert r.agreement < 0.5, (
        "a balanced bimodal molecule must not read as agreeing")
    assert r.is_multimodal

    # The distribution says the same thing independently: the median pair is a
    # CROSS-mode pair, so the central value is far outside tolerance even
    # though most individual poses have a near-perfect neighbour.
    assert r.median_rmsd > 5.0
    assert r.iqr_rmsd[0] < TOL < r.iqr_rmsd[1], (
        "the IQR should straddle the tolerance for a bimodal set")


def test_the_bimodal_representative_is_a_real_pose_and_not_the_gap():
    """`pose_vector.representative` returns an index; the mean would be fiction.

    The mean reactive-region position of a bimodal set sits in the trough
    between the modes, where no pose ever was. Both are asserted so the test
    documents the difference rather than just forbidding one.
    """
    poses = _tight(5, 0.0) + _tight(5, 6.0, energy0=-8.5)
    r = pc.consensus(poses, top_n=10, tolerance_a=TOL)

    rep = poses[r.representative_index].reactive_xyz
    assert any(np.allclose(rep, p.reactive_xyz) for p in poses)

    mean_xyz = np.mean([p.reactive_xyz for p in poses], axis=0)
    assert not any(np.allclose(mean_xyz, p.reactive_xyz, atol=TOL)
                   for p in poses), (
        "the mean coincided with a real pose, so this fixture no longer "
        "demonstrates the hazard it exists to demonstrate")


# --- could-not-be-measured ------------------------------------------------

def test_too_few_poses_raise_rather_than_scoring_zero():
    """Two poses give one pairwise distance: a central value with no spread.

    Names `MIN_POSES_FOR_CONSENSUS` so that changing it fails here and only
    here.
    """
    assert pc.MIN_POSES_FOR_CONSENSUS == 3, (
        "two poses yield one pairwise distance, so there is no IQR and no "
        "jackknife; the floor exists to stop a bare central value being "
        "reported as if it had a spread")

    with pytest.raises(pc.ConsensusError, match="NOT been measured"):
        pc.consensus([_pose(-9.0, 0.0), _pose(-8.0, 0.0)], top_n=5,
                     tolerance_a=TOL)


def test_a_top_n_below_the_minimum_is_refused():
    """Derived FROM the constant, so it stays true under a mutation of it."""
    poses = _tight(5, 0.0)
    with pytest.raises(pc.ConsensusError, match="no spread"):
        pc.consensus(poses, top_n=pc.MIN_POSES_FOR_CONSENSUS - 1,
                     tolerance_a=TOL)


def test_an_empty_pose_list_is_not_a_zero():
    with pytest.raises(pc.ConsensusError):
        pc.consensus([], top_n=5, tolerance_a=TOL)


# --- N is part of the metric ----------------------------------------------

def test_top_n_has_no_default_and_must_be_named():
    """D0068: a number whose defining parameter is implicit is the failure mode."""
    with pytest.raises(TypeError):
        pc.consensus(_tight(5, 0.0))                     # type: ignore[call-arg]
    with pytest.raises(TypeError):
        pc.consensus(_tight(5, 0.0), 3)                  # type: ignore[misc]


def test_the_result_states_the_n_it_was_measured_at():
    poses = _tight(8, 0.0)
    r = pc.consensus(poses, top_n=5, tolerance_a=TOL)
    assert r.top_n_requested == 5
    assert r.n_poses == 5
    assert "N=5" in str(r)


def test_a_short_pose_set_is_measured_at_what_exists_and_says_so(caplog):
    poses = _tight(4, 0.0)
    with caplog.at_level("WARNING"):
        r = pc.consensus(poses, top_n=10, tolerance_a=TOL)
    assert (r.top_n_requested, r.n_poses) == (10, 4)
    assert any("not comparable" in rec.message for rec in caplog.records)


def test_results_at_different_n_refuse_to_be_compared():
    """Both values are populated, plausible and in [0, 1]; only the guard notices."""
    a = pc.consensus(_tight(10, 0.0), top_n=10, tolerance_a=TOL)
    b = pc.consensus(_tight(10, 0.0), top_n=4, tolerance_a=TOL)
    with pytest.raises(pc.ConsensusError, match="different N"):
        pc.require_same_n([a, b])

    c = pc.consensus(_tight(10, 0.0), top_n=10, tolerance_a=2.0)
    with pytest.raises(pc.ConsensusError, match="different tolerances"):
        pc.require_same_n([a, c])

    assert pc.require_same_n([a, a]) == (10, TOL)


# --- the tolerance --------------------------------------------------------

def test_the_default_tolerance_is_the_measured_placement_radius():
    """1.0 A is D0062's reactive-region radius, reused rather than re-chosen.

    The only test that depends on the default's VALUE. Everything else states
    its tolerance, so a change to this constant localises here.
    """
    assert pc.DEFAULT_TOLERANCE_A == 1.0


def test_the_default_is_what_is_used_when_the_tolerance_is_omitted():
    poses = _tight(5, 0.0) + _tight(5, 6.0, energy0=-8.5)
    assert (pc.consensus(poses, top_n=10)
            == pc.consensus(poses, top_n=10,
                            tolerance_a=pc.DEFAULT_TOLERANCE_A))


def test_a_zero_tolerance_is_refused():
    """Every pose becomes its own mode by construction, which measures nothing."""
    with pytest.raises(pc.ConsensusError, match="must be positive"):
        pc.consensus(_tight(5, 0.0), top_n=5, tolerance_a=0.0)


# --- what is compared -----------------------------------------------------

def test_translated_poses_do_not_agree_even_with_identical_internal_geometry():
    """No superposition — the question is where the warhead sits IN THE SITE.

    These two poses have exactly the same internal geometry, so a superposing
    RMSD would call them identical. They are 5 A apart in the pocket.
    """
    a, b = _pose(-9.0, 0.0), _pose(-9.0, 5.0)
    assert pc.reactive_rmsd(a, b) == pytest.approx(5.0)

    r = pc.consensus([a, b, _pose(-9.0, 10.0)], top_n=3, tolerance_a=TOL)
    assert r.agreement == pytest.approx(0.0)


def test_poses_describing_different_atoms_refuse_to_compare():
    """A symmetric warhead matched at two centres must not be averaged over."""
    poses = [_pose(-9.0, 0.0), _pose(-8.9, 0.0),
             _pose(-8.8, 0.0, atom_ids=(11, 12, 13))]
    with pytest.raises(pc.ConsensusError, match="not the same atoms"):
        pc.consensus(poses, top_n=3, tolerance_a=TOL)


def test_a_pose_with_no_reactive_atoms_is_an_error():
    with pytest.raises(pc.ConsensusError, match="no reactive atoms"):
        pc.ReactivePose(-9.0, np.empty((0, 3)), ())


def test_a_pose_with_no_energy_cannot_enter_a_top_n():
    with pytest.raises(pc.ConsensusError, match="cannot be placed"):
        pc.ReactivePose(float("nan"), _warhead(0.0), ATOMS)


# --- energy selects, and does not weight ----------------------------------

def test_energy_chooses_which_poses_are_compared():
    """The top-3 are one tight mode; widening to 6 pulls in the scattered tail."""
    poses = _tight(3, 0.0) + [_pose(-5.0 + j, 8.0 * (j + 1)) for j in range(3)]
    assert pc.consensus(poses, top_n=3, tolerance_a=TOL).agreement == pytest.approx(1.0)
    assert pc.consensus(poses, top_n=6, tolerance_a=TOL).agreement < 0.5


def test_energy_values_do_not_enter_the_number():
    """Only the ORDER matters; the score contributes no weight.

    Consensus is meant to be information the docking score does not already
    carry (D0041 and four other levels of theory have failed on that score), so
    rescaling every energy must leave the number untouched — and reversing the
    order must too, when the whole set is selected either way.
    """
    coords = [0.0, 0.05, 0.1, 6.0, 6.05, 6.1]
    base = [_pose(-9.0 + j, c) for j, c in enumerate(coords)]
    rescaled = [_pose(10.0 * (-9.0 + j), c) for j, c in enumerate(coords)]
    reversed_order = [_pose(-(-9.0 + j), c) for j, c in enumerate(coords)]

    r = pc.consensus(base, top_n=6, tolerance_a=TOL)
    assert pc.consensus(rescaled, top_n=6, tolerance_a=TOL) == r
    assert pc.consensus(reversed_order, top_n=6, tolerance_a=TOL).agreement == r.agreement


def test_energy_ties_are_broken_deterministically():
    """A re-run must not silently swap a tied pair and move the answer.

    The tie is placed ON the top-N boundary so it actually decides the answer:
    poses 2 and 3 both score -8.0, only one fits in the top 3, and they sit in
    different modes. Earlier position wins, so the tight pose is selected and
    agreement is 1.0. If the rule were unstable this would flip to 1/3.
    """
    poses = [_pose(-9.0, 0.0), _pose(-8.5, 0.05),
             _pose(-8.0, 0.1), _pose(-8.0, 20.0), _pose(-7.0, 0.15)]
    first = pc.consensus(poses, top_n=3, tolerance_a=TOL)

    assert first == pc.consensus(poses, top_n=3, tolerance_a=TOL)
    assert first.agreement == pytest.approx(1.0)
    assert first.n_modes == 1


# --- never a bare central value -------------------------------------------

def test_no_central_value_is_reported_without_its_spread():
    """Three tight poses and one outlier: the headline rests on which pose left.

    Dropping the outlier gives 1.0; dropping a tight pose gives 1/3. A bare
    0.5 would hide that entirely, which is D0065's four-seed lesson (AUC 0.872,
    0.881, 0.852, 0.822 — quote the spread) arriving in a different metric.
    """
    poses = _tight(3, 0.0) + [_pose(-7.0, 20.0)]
    r = pc.consensus(poses, top_n=4, tolerance_a=TOL)

    assert r.agreement == pytest.approx(3 / 6)
    lo, hi = r.agreement_jackknife
    assert lo == pytest.approx(1 / 3)
    assert hi == pytest.approx(1.0)
    assert lo < r.agreement < hi
    # And the RMSD side says the same: the quartiles sit on opposite sides of
    # the tolerance, so no single central value describes this set either.
    assert r.iqr_rmsd[0] < TOL < r.iqr_rmsd[1]
    assert r.iqr_rmsd[0] < r.median_rmsd < r.iqr_rmsd[1]


def test_the_rmsd_distribution_is_reported_so_the_tolerance_can_be_redrawn():
    """Raw distances survive alongside the verdict, as `NACResult`'s do."""
    poses = _tight(3, 0.0) + [_pose(-7.0, 3.0)]
    r = pc.consensus(poses, top_n=4, tolerance_a=TOL)
    assert r.max_rmsd == pytest.approx(3.0)
    assert r.iqr_rmsd[0] <= r.median_rmsd <= r.iqr_rmsd[1]
    # Re-drawing the tolerance at 3.5 A would merge everything, and the reported
    # maximum is enough to know that without re-docking.
    assert pc.consensus(poses, top_n=4, tolerance_a=3.5).agreement == pytest.approx(1.0)


# --- the representative is always real ------------------------------------

def test_the_representative_indexes_the_callers_own_pose_list():
    """Not the top-N slice: an index that addresses the wrong list is D0047's shape."""
    poses = [_pose(-5.0, 20.0)] + _tight(4, 0.0)      # worst pose FIRST
    r = pc.consensus(poses, top_n=4, tolerance_a=TOL)
    assert 1 <= r.representative_index <= 4, (
        "the representative must be one of the top-4, addressed in the "
        "caller's indexing")
    assert np.allclose(poses[r.representative_index].reactive_xyz,
                       _warhead(0.0), atol=0.3)
