"""The second-stage split (#61), and the properties it must not break."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_subsplit as psub                   # noqa: E402


def _cloud(groups, n_atoms=8, spread=0.15, seed=0):
    """Poses in `len(groups)` well-separated clusters, `groups[i]` poses each."""
    rng = np.random.default_rng(seed)
    out = []
    for g, n in enumerate(groups):
        centre = np.zeros((n_atoms, 3)) + g * 10.0
        for _ in range(n):
            out.append(centre + rng.normal(0, spread, (n_atoms, 3)))
    return np.array(out)


def test_a_mode_holding_two_separated_clouds_is_split():
    coords = _cloud([30, 30])
    labels = np.zeros(60, dtype=int)                    # one mode, two clouds
    new, info = psub.subdivide(labels, coords, max_sub=5)
    assert len(set(new)) == 2, "two well-separated clusters must become two modes"
    assert info["subdivided"] == 1


def test_a_tight_mode_is_left_alone():
    coords = _cloud([40], spread=0.05)
    labels = np.zeros(40, dtype=int)
    new, _ = psub.subdivide(labels, coords, max_sub=5)
    assert len(set(new)) <= 5
    assert (new >= 0).all()


def test_no_pose_is_ever_dropped():
    """Splitting may only redistribute. Losing a pose here loses it everywhere:
    the exported representatives are all anything downstream sees."""
    coords = _cloud([25, 25, 25])
    labels = np.zeros(75, dtype=int)
    new, _ = psub.subdivide(labels, coords, max_sub=5)
    assert (new >= 0).sum() == 75


def test_unassigned_poses_stay_unassigned():
    coords = _cloud([20, 20])
    labels = np.array([0] * 20 + [-1] * 20)
    new, _ = psub.subdivide(labels, coords, max_sub=5)
    assert (new[20:] == -1).all(), "-1 means the clustering rejected it; that survives"


def test_small_modes_are_not_subdivided():
    """Sub-clusters of two or three poses are noise with a row and a sweep."""
    coords = _cloud([4, 4])
    labels = np.zeros(8, dtype=int)
    new, info = psub.subdivide(labels, coords, max_sub=5, min_size=12)
    assert len(set(new)) == 1
    assert info["subdivided"] == 0


def test_max_sub_of_one_is_a_no_op():
    coords = _cloud([30, 30])
    labels = np.zeros(60, dtype=int)
    new, info = psub.subdivide(labels, coords, max_sub=1)
    assert len(set(new)) == 1 and info["subdivided"] == 0


def test_labels_are_renumbered_contiguously_from_zero():
    """Downstream treats these as ordinary modes; a gap in the numbering would
    make `_m<k>` idents non-contiguous and the rank tables harder to join."""
    coords = _cloud([30, 30, 30])
    labels = np.array([0] * 30 + [1] * 30 + [2] * 30)
    new, _ = psub.subdivide(labels, coords, max_sub=5)
    seen = sorted(set(int(x) for x in new) - {-1})
    assert seen == list(range(len(seen)))
