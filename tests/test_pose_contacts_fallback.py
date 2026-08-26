"""The numpy fallbacks in `pose_contacts` must agree with scipy, exactly.

The GUI environment has no scipy, so `pose_distances` and `group` each carry a
numpy path. Two implementations of one metric is the near-miss-neighbour shape
this project keeps failing on -- both are populated, both are plausible, and a
divergence would show up as a slightly different grouping that nothing flags. So
they are compared here rather than trusted.

PARTITIONS, NOT LABEL INTEGERS. scipy numbers clusters by its own dendrogram
order; agreement means "the same poses are grouped together", so the comparison
is on the induced partition. Comparing integers would fail on a correct
implementation and pass on nothing useful.

THE COMPARISONS ARE CHECKED TO BE ABLE TO FAIL. `test_partition_check_can_fail`
feeds the same comparison a SINGLE-linkage clustering, which is the specific
wrong answer complete linkage protects against, and asserts it is rejected.
"""

from __future__ import annotations

import numpy as np
import pytest

from shared import pose_contacts as pc

scipy = pytest.importorskip("scipy", reason="the reference side of the comparison")


def _partition(labels) -> set:
    return {frozenset(np.flatnonzero(np.asarray(labels) == k).tolist())
            for k in set(np.asarray(labels).tolist())}


def _cloud(n=90, atoms=14, res=6, seed=0):
    rng = np.random.default_rng(seed)
    # three loose blobs, so the cut has something to separate
    centres = rng.uniform(2.0, 9.0, (3, atoms, res))
    T = np.concatenate([c + rng.normal(0, 0.55, (n // 3, atoms, res))
                        for c in centres]).astype(np.float32)
    return np.clip(T, 0.5, pc.CAP_A), pc.atom_weights(rng.uniform(0.3, 1.8, atoms))


def test_euclidean_fallback_matches_scipy():
    from scipy.spatial.distance import pdist, squareform
    rng = np.random.default_rng(3)
    U = rng.normal(0, 3.0, (120, 40))
    assert np.abs(_ref := squareform(pdist(U)) - pc._euclidean_numpy(U)).max() < 1e-9


def test_euclidean_fallback_is_a_metric_shape():
    U = np.random.default_rng(4).normal(0, 2.0, (40, 12))
    D = pc._euclidean_numpy(U)
    assert np.allclose(np.diag(D), 0.0)
    assert np.allclose(D, D.T)
    assert (D >= 0).all()


@pytest.mark.parametrize("tol", [0.4, 0.9, 1.6, 3.0])
def test_complete_linkage_fallback_matches_scipy(tol):
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    T, w = _cloud()
    D = pc.pose_distances(T, w)
    ref = fcluster(linkage(squareform(D, checks=False), method="complete"),
                   t=tol, criterion="distance") - 1
    assert _partition(pc._complete_linkage_numpy(D, tol)) == _partition(ref)


@pytest.mark.parametrize("tol", [0.4, 0.9, 1.6, 3.0])
def test_fallback_honours_the_guarantee(tol):
    """The whole reason complete linkage was chosen; it must hold on this path too."""
    T, w = _cloud(seed=5)
    D = pc.pose_distances(T, w)
    lab = pc._complete_linkage_numpy(D, tol)
    assert pc.within_group_max(D, lab) <= tol + 1e-9


def _chain(n=40, atoms=14, res=6, step=0.2):
    """Poses evenly spaced along one direction -- the case that separates the two.

    EVERY coordinate is shifted by `i * step`, which makes the contact distance
    exactly `|i - j| * step`: the atom weights sum to 1 and the residue average
    divides by the same n_r the sum introduces, so the metric's normalisation
    cancels and the fixture's geometry is exact rather than approximate.

    Single linkage then walks the chain into one cluster at any cut above `step`,
    while complete linkage cuts it into segments of diameter <= tol. The chain
    must SPAN more than tol or both are right and the fixture proves nothing --
    the first two versions of this failed exactly there, once with well-separated
    blobs and once with a chain shorter than the cut.
    """
    base = np.random.default_rng(2).uniform(0.6, 1.0, (atoms, res))
    T = np.stack([base + i * step for i in range(n)]).astype(np.float32)
    assert T.max() < pc.CAP_A, "the cap would flatten the chain"
    return T, np.ones(atoms)


def test_partition_check_can_fail():
    """The comparison must reject the wrong clustering, or it proves nothing.

    Single linkage is the specific failure complete linkage exists to prevent --
    it chains, which is D0088's defect. If the partition comparison passed it,
    every test above would be vacuous.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    T, w = _chain()
    D = pc.pose_distances(T, w)
    tol = 1.2
    single = fcluster(linkage(squareform(D, checks=False), method="single"),
                      t=tol, criterion="distance") - 1
    ours = pc._complete_linkage_numpy(D, tol)
    assert _partition(single) != _partition(ours), (
        "single and complete linkage agree on this cloud, so it cannot "
        "discriminate -- change the fixture, not this assertion")
    assert pc.within_group_max(D, single) > tol, (
        "single linkage did not violate the diameter bound here")


def test_group_uses_the_fallback_when_scipy_is_absent(monkeypatch):
    """Simulate the GUI environment and check the same answer comes back."""
    import builtins
    T, w = _cloud(seed=11)
    D = pc.pose_distances(T, w)
    want = pc.group(D, 1.0)
    real = builtins.__import__

    def no_scipy(name, *a, **k):
        if name.startswith("scipy"):
            raise ImportError("scipy is not available in dwi_gui")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_scipy)
    with pytest.raises(ImportError):
        __import__("scipy.cluster.hierarchy")
    assert _partition(pc.group(D, 1.0)) == _partition(want)
