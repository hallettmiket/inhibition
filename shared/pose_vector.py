"""
Purpose: turn a pose into a fixed-length pocket-contact vector; cluster poses; pick a REAL representative.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: a pose's heavy-atom coordinates + the receptor's pocket residues
Output: an (n_residues,) contact vector, cluster labels, and a medoid INDEX

Issue #14. The ranking question is "does this molecule sit in the pocket the way
things that bind sit in it", and that is a question about GEOMETRY, not about
the docking score -- which matters here because D0041 and D0046 have measured
that the score does not discriminate on this target.

## The vector

For each residue lining the site, the **minimum distance from any ligand heavy
atom to any heavy atom of that residue**. One number per residue, so the vector
is a fixed length regardless of how big the molecule is -- a 12-heavy-atom T_1
fragment and a 45-atom liu derivative produce vectors of the same shape and are
directly comparable. That is the property the whole idea needs, and it is why
this is a per-RESIDUE reduction rather than a per-ATOM one.

**It is agnostic by construction.** Nobody declares that the warhead belongs
near Cys113 or the R-group in the proline pocket. The vector records where the
molecule actually is; what counts as a good arrangement is learned separately,
from poses that were experimentally observed (see `reference_profile`). T_3's
prior knowledge falls out of that automatically if it is right, and is
contradicted by the data if it is not.

**It is symmetry-invariant for free.** Two poses related by flipping a
symmetric ring have identical contacts and land on the same vector. `redock_04`
needed `symmetric_rmsd` and a graph match to handle that explicitly; here it
costs nothing, because a contact profile cannot see which of two equivalent
atoms is which. Where RMSD asks "did the atoms land in the same places", this
asks "did the molecule touch the same things" -- which is nearer the question a
chemist reading a campaign is actually asking.

## Why there is no `mean_pose` in this module

Averaging conformations is the one operation that must not exist here. If one
pose places a warhead in the Cys113 sub-pocket and another has it flipped, the
average places it BETWEEN them -- somewhere neither pose ever visited, quite
possibly inside the protein. The resulting fit score is a perfectly plausible
number computed from a structure that does not exist, which is the exact failure
shape catalogued in `docs/how_this_project_breaks.md`.

So `representative()` returns an **INDEX into the poses you passed in**. It is
always a real, physically realisable pose that a program actually generated.
The medoid (minimum summed distance to the rest of its cluster) is the closest
honest analogue of "the average pose", and `mean_vector()` exists only for
summarising a set of ALREADY-CLUSTERED vectors -- it is documented as not being
a pose and cannot be handed back as one.

## What this module does NOT claim

A contact vector says where a molecule sits. It says nothing about whether it
STAYS there -- that is what BPMD is for (#14) -- and nothing about affinity.
Per #13, no score derived from this may label anything downstream until it has
been scored on D0046's 80 redock cases against crystal ground truth. Top-1 is
22.5% and best-of-9 is 55%; those are the numbers to beat and the ceiling.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

log = logging.getLogger(__name__)

# Distance beyond which a residue is "not touched". Contacts saturate: whether
# a residue is 14 A or 22 A away is not information about the binding mode, and
# leaving it unclipped lets far-away residues dominate a Euclidean comparison
# purely through their spread. 8.0 A matches the shell `pocket_residues.json`
# was measured at, so the clip and the residue list agree about what "lining
# the site" means.
CONTACT_CEILING_A = 8.0


@dataclass(frozen=True)
class PoseVector:
    """A pose's contact profile, with the residue order it was built against.

    `resi` travels WITH the vector rather than being assumed. Two vectors built
    against different residue orders are not comparable, and a bare ndarray
    gives a caller no way to notice -- it would subtract them happily and
    return a plausible number. `require_same_basis` makes that an error.
    """

    resi: tuple[int, ...]
    values: np.ndarray

    def __post_init__(self) -> None:
        if len(self.resi) != len(self.values):
            raise ValueError(
                f"{len(self.resi)} residues but {len(self.values)} values")


def require_same_basis(vectors: list[PoseVector]) -> tuple[int, ...]:
    """Every vector must be built against the same residues, in the same order."""
    if not vectors:
        raise ValueError("no vectors")
    basis = vectors[0].resi
    for v in vectors[1:]:
        if v.resi != basis:
            raise ValueError(
                "pose vectors were built against different residue bases and "
                "cannot be compared; rebuild them against one receptor")
    return basis


def contact_vector(ligand_xyz: np.ndarray,
                   residue_xyz: dict[int, np.ndarray],
                   resi: tuple[int, ...] | None = None,
                   ceiling: float = CONTACT_CEILING_A) -> PoseVector:
    """Minimum ligand-heavy-atom distance to each pocket residue, clipped.

    `ligand_xyz` is (n_atoms, 3) heavy atoms only. `residue_xyz` maps residue
    number -> (n_atoms, 3) for that residue's heavy atoms.

    A residue with no coordinates is recorded at the ceiling rather than
    dropped: dropping it would silently shorten the vector for one pose and
    leave it comparable-looking against full-length ones.
    """
    if ligand_xyz.ndim != 2 or ligand_xyz.shape[1] != 3:
        raise ValueError(f"ligand_xyz must be (n, 3), got {ligand_xyz.shape}")
    if ligand_xyz.shape[0] == 0:
        raise ValueError("pose has no heavy atoms")

    order = tuple(resi) if resi is not None else tuple(sorted(residue_xyz))
    out = np.full(len(order), ceiling, dtype=float)
    for i, r in enumerate(order):
        rx = residue_xyz.get(r)
        if rx is None or len(rx) == 0:
            continue
        d = np.linalg.norm(ligand_xyz[:, None, :] - rx[None, :, :], axis=2)
        out[i] = min(float(d.min()), ceiling)
    return PoseVector(order, out)


def distance_matrix(vectors: list[PoseVector]) -> np.ndarray:
    """Pairwise Euclidean distance between contact profiles."""
    require_same_basis(vectors)
    m = np.vstack([v.values for v in vectors])
    return np.linalg.norm(m[:, None, :] - m[None, :, :], axis=2)


def cluster(vectors: list[PoseVector], threshold: float) -> list[int]:
    """Single-linkage clustering of poses by contact profile.

    Returns a label per pose. Single-linkage on purpose: two poses that differ
    by a small rotation form a chain of near-neighbours and belong together,
    and complete-linkage would split them on the widest pair.

    THIS IS NOT A CONFIDENCE MEASURE. Cluster population over Vina's nine modes
    counts the output diversification filter, not agreement -- Vina spreads its
    reported modes with a minimum-RMSD floor, and our own median
    nearest-neighbour RMSD of 1.176 A sits on it. Population here is only
    meaningful across INDEPENDENT runs (different seeds) or when the poses came
    from somewhere other than one Vina invocation.
    """
    n = len(vectors)
    if n == 0:
        return []
    d = distance_matrix(vectors)
    label = list(range(n))
    changed = True
    while changed:                       # tiny n (9 modes); clarity over speed
        changed = False
        for i in range(n):
            for j in range(i + 1, n):
                if d[i, j] <= threshold and label[i] != label[j]:
                    lo, hi = sorted((label[i], label[j]))
                    label = [lo if x == hi else x for x in label]
                    changed = True
    remap = {v: k for k, v in enumerate(sorted(set(label)))}
    return [remap[x] for x in label]


def representative(vectors: list[PoseVector],
                   members: list[int] | None = None,
                   weights: np.ndarray | None = None) -> int:
    """The medoid's INDEX — a real pose, never a synthesised one.

    Returns an index into `vectors`, so the caller gets back a conformation
    some program actually generated. This is the whole reason the function
    returns an int rather than a vector: an index cannot be mistaken for, or
    quietly turned into, an averaged structure.

    `weights` (e.g. BPMD confidence per pose) biases WHICH REAL POSE is chosen
    -- a high-confidence pose is more likely to be the representative -- and
    never blends coordinates. That is the honest reading of "a weighted average
    pose": weight the choice, not the geometry.
    """
    idx = list(range(len(vectors))) if members is None else list(members)
    if not idx:
        raise ValueError("no members to choose a representative from")
    if len(idx) == 1:
        return idx[0]
    d = distance_matrix([vectors[i] for i in idx])
    cost = d.sum(axis=1)
    if weights is not None:
        w = np.asarray(weights, dtype=float)[idx]
        if np.any(w < 0):
            raise ValueError("weights must be non-negative")
        if not np.any(w > 0):
            raise ValueError("all weights are zero; no pose can be preferred")
        # Divide, so HIGHER confidence lowers the cost. Guarded against zero.
        cost = cost / np.clip(w, 1e-9, None)
    return idx[int(np.argmin(cost))]


def mean_vector(vectors: list[PoseVector]) -> np.ndarray:
    """Elementwise mean of contact profiles. **NOT A POSE.**

    Deliberately returns a bare ndarray and not a `PoseVector`, so it cannot be
    passed anywhere a pose is expected. Legitimate for summarising a set of
    vectors that are already known to be one cluster -- e.g. building a
    reference profile from many crystal structures. NOT legitimate as a stand-in
    for a conformation: the mean of a bimodal set lies in the trough, where no
    pose ever was. Use `representative()` for that.
    """
    require_same_basis(vectors)
    return np.vstack([v.values for v in vectors]).mean(axis=0)


def is_multimodal(vectors: list[PoseVector], threshold: float) -> bool:
    """Does this set contain more than one binding mode at `threshold`?

    The guard to call BEFORE summarising a set with `mean_vector`. If it is
    True, the mean describes nothing real.
    """
    return len(set(cluster(vectors, threshold))) > 1


def reference_profile(vectors: list[PoseVector],
                      threshold: float) -> tuple[np.ndarray, np.ndarray]:
    """(median, MAD) contact profile over EXPERIMENTALLY OBSERVED poses.

    This is what makes the fit score agnostic: the expected arrangement is
    learned from crystal structures rather than declared by hand. Median and
    MAD rather than mean and SD because the input is a small heterogeneous set
    of real ligands and one unusual binder should not move the target.

    Warns rather than raising on a multimodal input -- Pin1's site genuinely
    binds more than one way, and a caller may legitimately want a profile per
    mode. It must be a decision, not an accident.

    LEAKAGE IS THE CALLER'S PROBLEM AND IT IS A REAL ONE: fitting this on the
    same actives the enrichment gate scores means training on the test set,
    which is exactly why Boltz-2 was ruled out (state_of_the_project section 4).
    Hold the gate's actives out.
    """
    require_same_basis(vectors)
    m = np.vstack([v.values for v in vectors])
    med = np.median(m, axis=0)
    mad = np.median(np.abs(m - med), axis=0)
    if is_multimodal(vectors, threshold):
        log.warning(
            "the reference poses are MULTIMODAL at threshold %.2f — a single "
            "median profile describes none of the modes well. Consider one "
            "profile per mode.", threshold)
    return med, mad


def fit_score(vector: PoseVector, median: np.ndarray, mad: np.ndarray,
              floor: float = 0.5) -> float:
    """How far this pose's contacts sit from the observed-binder profile.

    LOWER IS BETTER — a robust z-like distance in units of MAD, so a residue
    whose contact distance varies a lot among real binders contributes less
    than one that is always the same. `floor` keeps a residue with near-zero
    MAD from dominating through division by a tiny number.

    Registered in `rank_shortlist.LOWER_IS_BETTER` before anything sorts on it;
    that registry is what caught the cnn_affinity/affinity_kcal mix-up.
    """
    if len(vector.values) != len(median) or len(median) != len(mad):
        raise ValueError("vector, median and mad have different lengths")
    z = (vector.values - median) / np.clip(mad, floor, None)
    return float(np.sqrt(np.mean(z ** 2)))
