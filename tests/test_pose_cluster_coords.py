"""HDBSCAN in 3N coordinate space: the properties `cluster_coords` must hold.

The claim in the docstring -- that clustering the raw coordinates gives the SAME
partition as clustering the precomputed RMSD matrix -- is an argument about a
constant scale factor. Arguments about scale factors are exactly the kind of
thing this project has been wrong about before, so it is asserted on real
clouds rather than reasoned about in a comment.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_cluster as pclust                   # noqa: E402

pytest.importorskip("sklearn", reason="HDBSCAN lives in scikit-learn")


def _cloud(groups, n_atoms=12, spread=0.2, seed=0):
    """Poses in `len(groups)` separated clusters, `groups[i]` poses each."""
    rng = np.random.default_rng(seed)
    out = []
    for g, n in enumerate(groups):
        centre = rng.normal(0, 8.0, (n_atoms, 3)) + g * 25.0
        for _ in range(n):
            out.append(centre + rng.normal(0, spread, (n_atoms, 3)))
    return np.array(out)


def _same_partition(a: np.ndarray, b: np.ndarray) -> bool:
    """Do two labellings group the poses identically? Ids need not match."""
    if len(a) != len(b):
        return False
    if ((a < 0) != (b < 0)).any():
        return False
    pairs_a = {(i, j) for i in range(len(a)) for j in range(i + 1, len(a))
               if a[i] >= 0 and a[i] == a[j]}
    pairs_b = {(i, j) for i in range(len(b)) for j in range(i + 1, len(b))
               if b[i] >= 0 and b[i] == b[j]}
    return pairs_a == pairs_b


@pytest.mark.parametrize("seed", range(5))
def test_coords_and_rmsd_matrix_agree(seed):
    """The scale factor is `sqrt(n_atoms)`, and HDBSCAN is blind to it."""
    coords = _cloud([12, 9, 6, 4], seed=seed)
    assert _same_partition(pclust.cluster_coords(coords),
                           pclust.cluster(coords))


def test_finds_the_planted_groups():
    coords = _cloud([20, 15, 8])
    lab = pclust.cluster_coords(coords)
    sizes = sorted((int((lab == c).sum()) for c in set(lab) - {-1}),
                   reverse=True)
    assert sizes == [20, 15, 8]
    assert (lab < 0).sum() == 0


def test_a_lone_pose_is_an_orphan_not_a_mode():
    """A mode is a REPEATED arrangement. One draw with a label is not one."""
    coords = np.concatenate([_cloud([12, 9]), _cloud([1], seed=99) + 200.0])
    lab = pclust.cluster_coords(coords)
    assert lab[-1] == -1
    assert (lab >= 0).sum() == 21


def test_one_arrangement_only_comes_back_all_orphan():
    """KNOWN BEHAVIOUR, asserted so it is not rediscovered as a bug.

    sklearn's `allow_single_cluster` defaults to False, so a cloud holding ONE
    arrangement has no split to make and every pose is labelled noise. Harmless
    on a 400-pose docking cloud, which yields 35-81 modes; it would bite a
    caller that fed in a handful of near-identical poses and read the empty
    result as "no modes found" rather than "one mode, unsplittable".
    """
    assert (pclust.cluster_coords(_cloud([10])) == -1).all()


def test_modes_are_renumbered_largest_first():
    """Ids must not follow scan order, or two runs disagree on 'mode 0'."""
    coords = np.concatenate([_cloud([5]), _cloud([14], seed=7) + 70.0])
    lab = pclust.cluster_coords(coords)
    assert int((lab == 0).sum()) == 14
    assert int((lab == 1).sum()) == 5


def test_too_few_poses_is_all_orphan_not_a_crash():
    coords = _cloud([2])
    assert (pclust.cluster_coords(coords, min_cluster_size=3) == -1).all()
    assert len(pclust.cluster_coords(np.empty((0, 5, 3)))) == 0


def test_min_cluster_size_is_a_count_not_a_distance():
    """Scaling every coordinate must not change the partition."""
    coords = _cloud([10, 7])
    assert _same_partition(pclust.cluster_coords(coords),
                           pclust.cluster_coords(coords * 3.0))
