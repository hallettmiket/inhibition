"""`pose_contacts.split_poses` -- the production entry point, and its guarantees.

The two rules this replaces both failed the same way: they produced groups whose
members did not resemble each other, and nothing in the code checked. So the
tests here are about the PROPERTIES a caller relies on, not about reproducing
numbers:

  * every pose gets a label and NO pose is dropped -- there is no noise class,
    which is how `<topic>_allposes` came to hold 79% of its cloud (D0093);
  * within-group distance is bounded by the tolerance STRUCTURALLY;
  * a looser tolerance can only merge, never split (monotonicity) -- a rule that
    fails this has no length scale, which is D0090's finding about HDBSCAN;
  * the landmark set is resolved by glob and RAISES on a missing residue rather
    than quietly measuring in fewer dimensions.
"""

from __future__ import annotations

import numpy as np
import pytest

from shared import pose_contacts as pc


# --------------------------------------------------------------------------- #
#  landmark resolution
# --------------------------------------------------------------------------- #
def test_landmarks_resolve_by_glob_and_exclude_water():
    names = pc.landmark_residues(15)
    assert len(names) == 15
    assert len(set(names)) == 15, "a landmark is listed twice"
    assert not any(n.endswith(":HOH") for n in names), (
        "water is a landmark modelled inconsistently between structures; it is "
        "excluded in the reference file, not by the caller")
    assert all(n.count(":") == 2 for n in names)


def test_landmarks_refuse_more_than_exist():
    with pytest.raises(ValueError):
        pc.landmark_residues(10_000)


def test_missing_landmark_raises_rather_than_shrinking_the_metric(tmp_path):
    """A dropped landmark changes the metric's dimension without changing its
    name, and every tolerance measured under the old dimension still looks fine."""
    pdb = tmp_path / "tiny.pdb"
    pdb.write_text(
        "ATOM      1  CA  SER A 114      10.000  10.000  10.000  1.00  0.00           C\n")
    with pytest.raises(ValueError, match="no heavy atoms"):
        pc.receptor_landmark_coords(["A:114:SER", "A:68:ARG"], receptor_pdb=pdb)


# --------------------------------------------------------------------------- #
#  the guarantees, on synthetic clouds
# --------------------------------------------------------------------------- #
def _cloud(n=80, atoms=12, res=6, seed=0, spread=0.5):
    rng = np.random.default_rng(seed)
    centres = rng.uniform(2.0, 8.0, (4, atoms, res))
    T = np.concatenate([c + rng.normal(0, spread, (n // 4, atoms, res))
                        for c in centres]).astype(np.float32)
    return np.clip(T, 0.5, pc.CAP_A), pc.atom_weights(rng.uniform(0.3, 1.6, atoms))


@pytest.mark.parametrize("tol", [0.3, 0.8, 1.5, 3.0])
def test_every_pose_is_labelled_and_nothing_is_dropped(tol):
    T, w = _cloud()
    lab = pc.group(pc.pose_distances(T, w), tol)
    assert len(lab) == len(T)
    assert lab.min() >= 0, "there must be no noise label; a lone pose is a group of one"
    assert set(lab.tolist()) == set(range(lab.max() + 1)), "labels must be dense"
    assert np.bincount(lab).sum() == len(T)


@pytest.mark.parametrize("tol", [0.3, 0.8, 1.5, 3.0])
def test_within_group_distance_is_bounded_by_the_tolerance(tol):
    T, w = _cloud(seed=3)
    D = pc.pose_distances(T, w)
    assert pc.within_group_max(D, pc.group(D, tol)) <= tol + 1e-9


def test_looser_tolerance_only_merges():
    """Monotonicity. A rule without a length scale fails this (D0090)."""
    T, w = _cloud(seed=5)
    D = pc.pose_distances(T, w)
    counts = [pc.group(D, t).max() + 1 for t in (0.3, 0.6, 1.0, 1.6, 2.4, 3.2)]
    assert counts == sorted(counts, reverse=True), (
        f"group count must fall monotonically as the cut loosens, got {counts}")


def test_two_separated_populations_are_not_merged():
    """The bound is worth nothing if the rule cannot separate the obvious case."""
    rng = np.random.default_rng(11)
    a = np.full((30, 10, 5), 3.0) + rng.normal(0, 0.05, (30, 10, 5))
    b = np.full((30, 10, 5), 8.0) + rng.normal(0, 0.05, (30, 10, 5))
    T = np.concatenate([a, b]).astype(np.float32)
    lab = pc.group(pc.pose_distances(T, np.ones(10)), 1.0)
    assert lab[:30].std() == 0 and lab[30:].std() == 0
    assert lab[0] != lab[30], "two populations 5 A apart were put in one group"


def test_medoids_are_real_poses():
    T, w = _cloud(seed=7)
    D = pc.pose_distances(T, w)
    lab = pc.group(D, 0.9)
    med = pc.medoids(D, lab)
    assert set(med) == set(range(lab.max() + 1))
    for k, i in med.items():
        assert lab[i] == k, "a medoid must belong to the group it represents"
        assert 0 <= i < len(T)


# --------------------------------------------------------------------------- #
#  the tolerance, and the weakness it carries
# --------------------------------------------------------------------------- #
def test_tolerance_is_the_calibrated_median():
    rmsf = np.array([1.0, 2.0, 3.0, 4.0])
    assert pc.tolerance_for(rmsf) == pytest.approx(2.5 / pc.RMSF_CALIBRATION)


def test_descriptor_tolerance_refuses_unfitted_coefficients():
    """D0094: the coefficients are deliberately not pinned in source."""
    with pytest.raises(ValueError, match="fitted coefficients"):
        pc.tolerance_from_descriptors(None)


def test_atom_weights_are_clipped_and_mean_one():
    w = pc.atom_weights(np.array([0.001, 1.0, 1.0, 1.0, 50.0]))
    assert w.max() <= pc.MAX_WEIGHT_RATIO + 1e-9
    assert w.min() >= 1.0 / pc.MAX_WEIGHT_RATIO - 1e-9
