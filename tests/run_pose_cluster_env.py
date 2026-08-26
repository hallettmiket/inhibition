"""
Purpose: run the coordinate-space HDBSCAN checks in the env that actually has scikit-learn.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-19
Input: shared/pose_cluster.py
Output: pass/fail per property, non-zero exit on any failure

WHY THIS EXISTS. `tests/test_pose_cluster_coords.py` asserts the property the
new splitter's docstring claims -- that clustering raw 3N coordinates gives the
same partition as clustering the precomputed RMSD matrix -- and it needs
scikit-learn. sklearn lives in `dwi_admet`; `dwi_admet` has no pytest; the
suite's own interpreter (`dwi_cheminf`) has pytest and no sklearn. So under the
suite it `importorskip`s and reports as a skip, which reads as "covered".

Same shape and same reasoning as `run_app_renders_gui_env.py` (#45): a skipped
test that reads as covered is worse than no test. A PLAIN SCRIPT, not a pytest
file, so it runs under an interpreter with no pytest in it -- which is why the
assertions are mirrored here rather than imported. When either env gains the
other's dependency this can be deleted in favour of the pytest file.

    /data/lab_vm/envs/dwi_admet/bin/python tests/run_pose_cluster_env.py
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_cluster as pclust                    # noqa: E402


def cloud(groups, n_atoms=12, spread=0.2, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for g, n in enumerate(groups):
        centre = rng.normal(0, 8.0, (n_atoms, 3)) + g * 25.0
        for _ in range(n):
            out.append(centre + rng.normal(0, spread, (n_atoms, 3)))
    return np.array(out)


def same_partition(a, b) -> bool:
    """Do two labellings group the poses identically? Ids need not match."""
    if len(a) != len(b) or ((a < 0) != (b < 0)).any():
        return False
    def pairs(l):
        return {(i, j) for i in range(len(l)) for j in range(i + 1, len(l))
                if l[i] >= 0 and l[i] == l[j]}
    return pairs(a) == pairs(b)


def coords_and_rmsd_matrix_agree(seed):
    c = cloud([12, 9, 6, 4], seed=seed)
    assert same_partition(pclust.cluster_coords(c), pclust.cluster(c))


def finds_the_planted_groups():
    lab = pclust.cluster_coords(cloud([20, 15, 8]))
    assert sorted((int((lab == c).sum()) for c in set(lab) - {-1}),
                  reverse=True) == [20, 15, 8]
    assert (lab < 0).sum() == 0


def a_lone_pose_is_an_orphan_not_a_mode():
    lab = pclust.cluster_coords(
        np.concatenate([cloud([12, 9]), cloud([1], seed=99) + 200.0]))
    assert lab[-1] == -1 and (lab >= 0).sum() == 21


def one_arrangement_only_comes_back_all_orphan():
    """KNOWN BEHAVIOUR, asserted so it is not rediscovered as a bug.

    sklearn's `allow_single_cluster` defaults to False, so a cloud holding ONE
    arrangement has no split to make and every pose is labelled noise. Harmless
    on a 400-pose docking cloud, which yields 35-81 modes; it would bite a
    caller that fed in a handful of near-identical poses and read the empty
    result as "no modes found" rather than "one mode, unsplittable".
    """
    lab = pclust.cluster_coords(cloud([10]))
    assert (lab == -1).all()


def modes_are_renumbered_largest_first():
    lab = pclust.cluster_coords(
        np.concatenate([cloud([5]), cloud([14], seed=7) + 70.0]))
    assert int((lab == 0).sum()) == 14 and int((lab == 1).sum()) == 5


def too_few_poses_is_all_orphan_not_a_crash():
    assert (pclust.cluster_coords(cloud([2]), min_cluster_size=3) == -1).all()
    assert len(pclust.cluster_coords(np.empty((0, 5, 3)))) == 0


def min_cluster_size_is_a_count_not_a_distance():
    c = cloud([10, 7])
    assert same_partition(pclust.cluster_coords(c),
                          pclust.cluster_coords(c * 3.0))


def main() -> int:
    cases = [(f"coords_and_rmsd_matrix_agree[{s}]",
              (lambda s=s: coords_and_rmsd_matrix_agree(s))) for s in range(5)]
    cases += [(f.__name__, f) for f in (
        finds_the_planted_groups, a_lone_pose_is_an_orphan_not_a_mode,
        one_arrangement_only_comes_back_all_orphan,
        modes_are_renumbered_largest_first,
        too_few_poses_is_all_orphan_not_a_crash,
        min_cluster_size_is_a_count_not_a_distance)]
    rc = 0
    for name, fn in cases:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:                                    # noqa: BLE001
            rc = 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print("\n  all passed" if rc == 0 else "\n  FAILURES", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
