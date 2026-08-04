"""
Purpose: the pose vector is comparable across molecule sizes, and no representative is ever synthesised.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: shared/pose_vector.py
Output: pass/fail

WHY THIS EXISTS. Issue #14 proposes ranking on a "weighted average pose". The
averaging is the one part that cannot be built as described: if one pose puts a
warhead in the Cys113 sub-pocket and another has it flipped, the average puts it
between them -- somewhere neither pose ever was, quite possibly inside the
protein -- and its fit score is a perfectly plausible number computed from a
structure that does not exist. That is the failure shape the whole of
`docs/how_this_project_breaks.md` is about.

So the module returns an INDEX into real poses, and these tests pin that:

* `test_representative_is_always_one_of_the_inputs` — the returned pose is one
  that a program actually generated, on random inputs, not a blend.
* `test_a_flipped_pose_pair_never_averages_into_the_gap` — the concrete failure
  case, built deliberately bimodal. The medoid sits on a real mode; the mean
  sits in the trough where nothing is. Both are asserted, so the test documents
  the difference rather than just forbidding one of them.
* `test_weights_bias_the_choice_and_never_blend` — "weighted average pose" done
  honestly: confidence changes WHICH real pose is picked, never the geometry.

Plus the property the whole idea depends on: a 6-atom fragment and a 40-atom
derivative must produce vectors of the same shape, or nothing downstream can
compare them.
"""

from __future__ import annotations

import numpy as np
import pytest

from shared import pose_vector as pv


def _residues(n: int = 6, spacing: float = 5.0) -> dict[int, np.ndarray]:
    """n single-atom residues on a line, `spacing` apart."""
    return {100 + i: np.array([[i * spacing, 0.0, 0.0]]) for i in range(n)}


def _at(x: float, n_atoms: int = 1) -> np.ndarray:
    """A ligand of `n_atoms` heavy atoms sitting at x."""
    return np.array([[x, 0.0, 0.0]] * n_atoms, dtype=float)


# --- the vector -----------------------------------------------------------

def test_vector_length_is_independent_of_molecule_size():
    """The property everything downstream rests on.

    T_1 spans fragment and lead space (12 to 45 heavy atoms); if the vector
    grew with the molecule, no two arms could be compared at all.
    """
    res = _residues()
    small = pv.contact_vector(_at(0.0, n_atoms=6), res)
    large = pv.contact_vector(_at(0.0, n_atoms=40), res)
    assert small.values.shape == large.values.shape == (6,)
    assert small.resi == large.resi


def test_distances_are_clipped_at_the_ceiling():
    """Whether a residue is 14 A or 22 A away is not information.

    Unclipped, far residues dominate a Euclidean comparison purely by spread.
    """
    res = _residues(n=3, spacing=50.0)
    v = pv.contact_vector(_at(0.0), res)
    assert v.values[0] == pytest.approx(0.0)
    assert v.values[1] == pv.CONTACT_CEILING_A
    assert v.values[2] == pv.CONTACT_CEILING_A


def test_a_residue_with_no_atoms_is_recorded_not_dropped():
    """Dropping it would shorten one pose's vector and leave it looking valid."""
    res = _residues(n=3)
    res[101] = np.empty((0, 3))
    v = pv.contact_vector(_at(0.0), res)
    assert len(v.values) == 3
    assert v.values[1] == pv.CONTACT_CEILING_A


def test_vectors_on_different_bases_refuse_to_compare():
    a = pv.contact_vector(_at(0.0), _residues(n=4))
    b = pv.contact_vector(_at(0.0), _residues(n=5))
    with pytest.raises(ValueError, match="different residue bases"):
        pv.require_same_basis([a, b])


def test_an_empty_pose_is_an_error_not_a_zero_vector():
    with pytest.raises(ValueError, match="no heavy atoms"):
        pv.contact_vector(np.empty((0, 3)), _residues())


# --- the representative: the part #14 asked for ---------------------------

def test_representative_is_always_one_of_the_inputs():
    rng = np.random.default_rng(20260804)
    res = _residues(n=5)
    for _ in range(25):
        vs = [pv.contact_vector(_at(float(x)), res)
              for x in rng.uniform(0, 20, size=7)]
        i = pv.representative(vs)
        assert 0 <= i < len(vs)
        assert isinstance(i, int)


def test_a_flipped_pose_pair_never_averages_into_the_gap():
    """The concrete failure #14's 'weighted average pose' would produce.

    Two binding modes at opposite ends of the pocket. The medoid is one of the
    real poses. The MEAN sits between them, touching nothing either mode
    touches -- a plausible vector describing a structure that never existed.
    """
    res = _residues(n=5, spacing=5.0)
    left = [pv.contact_vector(_at(0.0 + j * 0.1), res) for j in range(3)]
    right = [pv.contact_vector(_at(20.0 + j * 0.1), res) for j in range(3)]
    poses = left + right

    assert pv.is_multimodal(poses, threshold=2.0), (
        "the fixture is not actually bimodal; this test would prove nothing")

    rep = pv.representative(poses)
    assert any(np.allclose(poses[rep].values, p.values) for p in poses)

    mean = pv.mean_vector(poses)
    assert not any(np.allclose(mean, p.values, atol=0.5) for p in poses), (
        "the mean coincided with a real pose, so this fixture no longer "
        "demonstrates the hazard it exists to demonstrate")


def test_mean_vector_is_not_a_pose_type():
    """It must not be passable anywhere a pose is expected."""
    res = _residues()
    vs = [pv.contact_vector(_at(float(x)), res) for x in (0.0, 1.0, 2.0)]
    assert not isinstance(pv.mean_vector(vs), pv.PoseVector)


def test_weights_bias_the_choice_and_never_blend():
    """'Weighted average pose', done honestly: weight WHICH, not the geometry."""
    res = _residues(n=5, spacing=5.0)
    poses = [pv.contact_vector(_at(x), res) for x in (0.0, 10.0, 20.0)]

    assert pv.representative(poses) == 1, "unweighted medoid should be the middle"

    heavy_on_first = np.array([50.0, 1.0, 1.0])
    chosen = pv.representative(poses, weights=heavy_on_first)
    assert chosen == 0
    # Still a real input pose, not a blend of pose 0 and the medoid.
    assert np.allclose(poses[chosen].values, poses[0].values)


def test_zero_weights_everywhere_is_refused():
    res = _residues()
    poses = [pv.contact_vector(_at(x), res) for x in (0.0, 5.0)]
    with pytest.raises(ValueError, match="all weights are zero"):
        pv.representative(poses, weights=np.zeros(2))


# --- clustering + fit score -----------------------------------------------

def test_clustering_separates_two_modes_and_merges_one():
    res = _residues(n=5, spacing=5.0)
    two = [pv.contact_vector(_at(x), res) for x in (0.0, 0.2, 20.0, 20.2)]
    assert len(set(pv.cluster(two, threshold=2.0))) == 2
    one = [pv.contact_vector(_at(x), res) for x in (0.0, 0.2, 0.4)]
    assert len(set(pv.cluster(one, threshold=2.0))) == 1


def test_fit_score_is_lower_for_a_pose_matching_the_reference():
    """LOWER IS BETTER, and it must actually discriminate."""
    res = _residues(n=5, spacing=5.0)
    reference = [pv.contact_vector(_at(0.0 + j * 0.1), res) for j in range(5)]
    med, mad = pv.reference_profile(reference, threshold=2.0)

    like_reference = pv.contact_vector(_at(0.05), res)
    unlike = pv.contact_vector(_at(20.0), res)
    assert pv.fit_score(like_reference, med, mad) < pv.fit_score(unlike, med, mad)


def test_reference_profile_warns_when_the_references_are_multimodal(caplog):
    """A single median over two modes describes neither."""
    res = _residues(n=5, spacing=5.0)
    mixed = ([pv.contact_vector(_at(0.0), res)] * 3
             + [pv.contact_vector(_at(20.0), res)] * 3)
    with caplog.at_level("WARNING"):
        pv.reference_profile(mixed, threshold=2.0)
    assert any("MULTIMODAL" in r.message for r in caplog.records)
