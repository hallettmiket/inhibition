"""Grouping poses by what they touch, at the molecule's own scale.

WHY NOT CARTESIAN RMSD, AND WHY NOT A SPATIAL PARTITION. Both were tried and
recorded. RMSD collapses a pose into one scalar, so it cannot say WHAT moved --
a rotated tail and a bodily shift score alike. A partition of 3D space cannot
bound a distance in configuration space at all: two poses may share a centroid
and be 180 degree flips, and measured, a 3 A partition put 1,758 of 6,000 poses
in one cell spanning 9.13 A (D0091).

WHAT A POSE IS HERE. A table of distances, ligand atom by key residue, capped.
That is a description of where the molecule sits RELATIVE TO THE POCKET rather
than in the box's coordinates -- so it is invariant to any rigid motion of the
complex, and orientation is carried implicitly: flip the molecule and a
different atom is the one near a given residue.

THE ATOM WEIGHTS ARE THE MOLECULE'S OWN; THE TOLERANCE, AS SHIPPED, IS NOT.
The distinction is D0094 and it was got wrong here first. Each atom's
contribution is weighted by the inverse of its predicted fluctuation, and that is
validated: `exp/15` ranks per-atom flexibility WITHIN a molecule at rho = 0.657
over 147 modes, 100% positive. A pose differing only in a wagging tail is the
same pose, because the tail wags anyway; the same deviation in a rigid core is a
different pose.

The TOLERANCE asks a different question -- how one molecule's scale compares to
another's -- and the ensemble does not answer it. ACROSS molecules the same
prediction correlates with measurement at rho = +0.112, CI [-0.06, +0.27],
crossing zero; the prediction varies at CV 0.15 where the truth varies at 0.45.
Dividing it by 2.21 does not beat writing one number down for every molecule
(Wilcoxon p = 0.515). Rotatable-bond count does, at rho = +0.523 -- see
`tolerance_from_descriptors`, which is measured but NOT adopted, because adopting
it re-groups every cloud and invalidates the measurements made under this one.

COMPLETE LINKAGE, BECAUSE IT IS THE ONLY ONE THAT GUARANTEES THE THING WE NEED.
Group distance is defined by the FARTHEST pair, so cutting the tree at `tol`
bounds every within-group distance by `tol` structurally, not on average and not
usually. DBSCAN's eps is a neighbour radius and chains (D0088); HDBSCAN has no
length scale at all and its cluster count grows linearly with sampling (D0090).

Cys113 IS NOT PRIVILEGED. It is one landmark among the rest, or absent. The
warhead's geometry against it is how poses are RANKED, one stage later; using it
to form the groups would grade them on the axis that built them (D0088).
"""

from __future__ import annotations

import numpy as np

#: Beyond this, "far" is "far": an uncapped tail lets one remote residue dominate.
CAP_A = 10.0
#: Predicted RMSF runs 2.21x the value MD measures (exp/15, n=147; 95% CI
#: [1.95, 2.51] over a cluster bootstrap by ident, and per-molecule ratios span
#: 0.90-6.80 with only 35% inside +-25% of it -- D0094). Applied so
#: the tolerance is expressed on the scale the trajectories actually show.
RMSF_CALIBRATION = 2.21
#: No atom weight may exceed this multiple of the median, or one near-rigid atom
#: dominates the metric and the others stop mattering.
MAX_WEIGHT_RATIO = 4.0


def contact_tensor(xyz: np.ndarray, residues: list, cap: float = CAP_A) -> np.ndarray:
    """(poses, atoms, residues): capped min distance from each atom to each residue.

    `residues` is a list of (n_res_atoms, 3) coordinate arrays.
    """
    n_p, n_a = xyz.shape[0], xyz.shape[1]
    out = np.empty((n_p, n_a, len(residues)), dtype=np.float32)
    flat = xyz.reshape(-1, 3)
    for j, r in enumerate(residues):
        d = np.sqrt(((flat[:, None, :] - r[None, :, :]) ** 2).sum(-1)).min(1)
        out[:, :, j] = np.minimum(d, cap).reshape(n_p, n_a)
    return out


def atom_weights(rmsf: np.ndarray, calibration: float = RMSF_CALIBRATION,
                 max_ratio: float = MAX_WEIGHT_RATIO) -> np.ndarray:
    """1/fluctuation per atom, normalised to mean 1 and clipped.

    Clipped because the weight is a reciprocal: an atom the ensemble happens to
    place almost identically in every conformer would otherwise carry unbounded
    weight and the metric would be about that one atom.
    """
    f = np.asarray(rmsf, dtype=float) / calibration
    f = np.maximum(f, 1e-3)
    w = 1.0 / f
    w = w / w.mean()
    return np.clip(w, 1.0 / max_ratio, max_ratio)


def pose_distances(T: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(poses, poses) weighted RMS contact deviation, in angstroms.

    Weighted over ATOMS and averaged over residues, so the number stays on the
    same scale as a distance however many residues are carried -- adding a
    landmark must not silently change what `tol` means.
    """
    n_p, n_a, n_r = T.shape
    # THE METRIC IS EUCLIDEAN IN A RESCALED SPACE, so it is computed as one.
    #   d(i,j)^2 = mean_r sum_a w_a (T[i,a,r] - T[j,a,r])^2
    #            = || U[i] - U[j] ||^2 / n_r      with U = sqrt(w_a) * T
    # Written as an explicit double loop this was O(n^2) python; as a pdist it
    # is the same number to floating point and reaches 6,000 poses.
    W = np.sqrt(w / w.sum())[None, :, None]
    U = (T * W).reshape(n_p, -1).astype(np.float64)
    try:
        from scipy.spatial.distance import pdist, squareform
        D = squareform(pdist(U))
    except ImportError:
        D = _euclidean_numpy(U)
    return D / np.sqrt(n_r)


def _euclidean_numpy(U: np.ndarray) -> np.ndarray:
    """Pairwise Euclidean distance without scipy, agreeing with pdist to 1e-9.

    THE GUI ENVIRONMENT HAS NO SCIPY. Rather than stand up a fourth environment
    to run a viewer, both scipy-dependent steps carry a numpy path, and
    `tests/test_pose_contacts_fallback.py` asserts the two agree -- a fallback
    nobody compares is a second implementation of the metric, which is exactly
    the near-miss-neighbour shape this project keeps failing on.

    The Gram identity loses precision when the subtraction cancels, so the
    diagonal is zeroed and negatives are clipped before the square root.
    """
    g = U @ U.T
    sq = np.diag(g)
    d2 = sq[:, None] + sq[None, :] - 2.0 * g
    np.maximum(d2, 0.0, out=d2)
    np.fill_diagonal(d2, 0.0)
    D = np.sqrt(d2)
    return (D + D.T) / 2.0


def group(D: np.ndarray, tol: float) -> np.ndarray:
    """Complete-linkage labels at a `tol` cut. Every within-group pair is <= tol."""
    n = len(D)
    if n < 2:
        return np.zeros(n, dtype=int)
    try:
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform
        Z = linkage(squareform(D, checks=False), method="complete")
        return fcluster(Z, t=tol, criterion="distance").astype(int) - 1
    except ImportError:
        return _complete_linkage_numpy(D, tol)


def _complete_linkage_numpy(D: np.ndarray, tol: float) -> np.ndarray:
    """Complete linkage cut at `tol`, without scipy.

    Lance-Williams for complete linkage is `d(A u B, C) = max(d(A,C), d(B,C))`,
    so merging is a row-wise maximum and no dendrogram has to be stored. Merging
    stops when the closest pair exceeds `tol`, which IS the cut -- the guarantee
    that every within-group pair is <= tol falls out of the update rule rather
    than being checked afterwards.

    Ties are broken by the lowest index pair, which scipy also does, so the two
    agree on label content; the LABEL NUMBERS need not match and the test
    compares partitions rather than integers.
    """
    n = len(D)
    W = np.array(D, dtype=float, copy=True)
    np.fill_diagonal(W, np.inf)
    alive = np.ones(n, dtype=bool)
    member = [[i] for i in range(n)]
    n_alive = n
    while n_alive > 1:
        # ARGMIN OVER THE WHOLE MATRIX, not over a gathered sub-block. Dead rows
        # and columns are already +inf, so the answer is the same, and gathering
        # `W[np.ix_(idx, idx)]` copied the full matrix on every one of ~n
        # iterations -- the copy, not the search, was the cost.
        flat = int(np.argmin(W))
        a, b = divmod(flat, n)
        if W[a, b] > tol:
            break
        a, b = (a, b) if a < b else (b, a)
        np.maximum(W[a], W[b], out=W[a])
        W[:, a] = W[a]
        W[a, a] = np.inf
        alive[b] = False
        W[b, :] = np.inf
        W[:, b] = np.inf
        member[a].extend(member[b])
        member[b] = []
        n_alive -= 1
    lab = np.empty(n, dtype=int)
    for k, a in enumerate(np.flatnonzero(alive)):
        lab[member[a]] = k
    return lab


def medoids(D: np.ndarray, labels: np.ndarray) -> dict:
    """{label: pose index} -- the most central member, never a synthetic mean."""
    out = {}
    for k in sorted(set(labels.tolist())):
        idx = np.flatnonzero(labels == k)
        out[int(k)] = int(idx[np.argmin(D[np.ix_(idx, idx)].sum(axis=1))])
    return out


def within_group_max(D: np.ndarray, labels: np.ndarray) -> float:
    """The largest within-group distance. MUST be <= tol; the guarantee is the point."""
    worst = 0.0
    for k in set(labels.tolist()):
        idx = np.flatnonzero(labels == k)
        if len(idx) > 1:
            worst = max(worst, float(D[np.ix_(idx, idx)].max()))
    return worst


#: Out-of-sample coefficients for log(RMSF) ~ a + b*rotatable_bonds, from
#: 20x5-fold grouped CV over 119 molecules (exp/18). PROPOSED, NOT ADOPTED.
DESCRIPTOR_TOLERANCE = {"intercept": None, "rotb": None}


def tolerance_from_descriptors(mol, ensemble_rmsf=None, coeffs=None) -> float:
    """A tolerance from rotatable-bond count, refitted rather than hardcoded.

    MEASURED AND NOT ADOPTED. `exp/18/tolerance_model.py` scores it out of sample
    at 26.2% median relative error against the shipped 32.2% and a flat
    constant's 33.1%, and 24.9% when the ensemble is carried alongside -- which is
    the measurement's OWN reproducibility (same molecule, different trajectory,
    CV 0.24), so it is at the ceiling this data can support. Adopting it re-groups
    every cloud, so it stays proposed until the re-screen (#79) that D0092's
    conclusions would otherwise have to be re-derived against.

    THE COEFFICIENTS ARE NOT BAKED IN. `DESCRIPTOR_TOLERANCE` is deliberately
    None: a fitted constant pinned in source is exactly the stale-pin family this
    project has hit five times (catalogue disguise #3), and this one would go
    stale the moment the training set grows. The caller passes coefficients from
    the fit artefact, or this raises.
    """
    c = coeffs or DESCRIPTOR_TOLERANCE
    if c.get("intercept") is None or c.get("rotb") is None:
        raise ValueError(
            "tolerance_from_descriptors needs fitted coefficients; run "
            "exp/18_rmsf_calibration/tolerance_model.py and pass its fit. They "
            "are not pinned here on purpose (D0094).")
    from rdkit.Chem import rdMolDescriptors
    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    v = float(np.exp(c["intercept"] + c["rotb"] * n_rot))
    if c.get("ens") is not None:
        if ensemble_rmsf is None:
            raise ValueError("these coefficients use the ensemble term; pass it")
        v = float(np.exp(c["intercept"] + c["rotb"] * n_rot
                         + c["ens"] * float(np.median(ensemble_rmsf))))
    return v
