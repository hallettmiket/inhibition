"""Second-stage splitting: subdivide a binding mode on whole-molecule geometry.

WHY (#61). The first stage clusters on the reactive atom's position and the
direction the warhead faces. That is what makes a mode mechanistically meaningful
and it is deliberately blind to the rest of the molecule -- so two poses that
place the warhead identically and hang the scaffold differently are ONE mode, and
one representative has to stand for both. Sulfopin's 456 poses formed a single
mode containing a pose 1.43 A from its crystal structure and the pose the
pipeline kept at 5.06 A.

HOW FINELY, MEASURED RATHER THAN CHOSEN. On the 82-case pose-recovery benchmark
(real Pin1 crystal ligands docked into our 3IKD, crystal pose known for each),
recovery within 2 A against the number of representatives carried:

    k = 1 (before)   22.0%      k = 4   39.0%
    k = 2            29.3%      k = 5   40.2%
    k = 3            36.6%      k = 9   41.5%   <- the ceiling, all poses

k=1 -> k=4 gains 14 cases and loses 0 (McNemar p = 1.2e-4). k=5 -> k=9 gains 1
and loses 0 (p = 1.00). Five representatives capture 97% of the ceiling, and
splitting cannot lose a case by construction: more representatives only add
coverage. Clustering also beats taking the top-k by energy at every k, by up to
7.3 points -- WHICH poses are kept matters, not just how many.

The 456-pose Sulfopin cloud agrees independently: k = 4 is the first cut whose
medoid lands inside 2 A (1.52 A), and k = 8..40 does no better.

WHAT THIS CHANGES DOWNSTREAM, AND IT IS NOT COSMETIC. Subdividing a mode changes
`consensus`, which is mode_size / n_poses, and `conditional_eb` is computed from
it. A library screened with sub-splitting is NOT score-comparable with one
screened without it. That is why this is a parameter with an explicit default
rather than a silent change, and why the sub-split count is recorded on every
row.
"""

from __future__ import annotations

import numpy as np

#: Representatives per mode. 5 sits at the measured plateau: 97% of the ceiling
#: on the benchmark, and indistinguishable from carrying every pose (p = 1.00).
DEFAULT_MAX_SUB = 5

#: A mode smaller than this is not subdivided. Sub-clusters of one or two poses
#: are not binding modes, they are noise given its own row, and each one costs a
#: sweep.
MIN_MODE_SIZE = 12

#: Cut diameter, Angstrom. A single representative cannot stand for poses further
#: apart than the accuracy being claimed, and 2 A is the bar this field uses to
#: call a docked pose correct -- the same bar the recovery numbers above are
#: measured at. Tying the cut to it means the resolution of the split and the
#: resolution of the claim are one number, not two.
CUT_A = 2.0


def _pairwise_rmsd(coords: np.ndarray) -> np.ndarray:
    """In-place heavy-atom RMSD between every pair of poses.

    NO SUPERPOSITION. Every pose is already in the receptor's frame, and fitting
    them onto each other would ask whether they are the same SHAPE rather than
    whether they are in the same PLACE -- and place is the whole question.
    """
    n = len(coords)
    d = np.zeros((n, n))
    for i in range(n):
        d[i] = np.sqrt(((coords - coords[i]) ** 2).sum(axis=2).mean(axis=1))
    np.fill_diagonal(d, 0.0)
    return d


def subdivide(labels: np.ndarray, coords: np.ndarray,
              max_sub: int = DEFAULT_MAX_SUB,
              min_size: int = MIN_MODE_SIZE) -> tuple[np.ndarray, dict]:
    """Split each mode into at most `max_sub` sub-modes on whole-molecule RMSD.

    `labels` is the first-stage mode assignment (-1 = unassigned, preserved).
    `coords` is (n_poses, n_heavy_atoms, 3) in the receptor frame.

    Returns (new_labels, info). New labels are renumbered contiguously from 0, so
    everything downstream keeps treating them as ordinary modes -- the mode
    abstraction gets finer, nothing else changes shape. Unassigned poses stay -1.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    labels = np.asarray(labels)
    out = np.full(len(labels), -1, dtype=int)
    info: dict = {"modes_in": 0, "modes_out": 0, "subdivided": 0}
    nxt = 0
    for m in sorted(set(int(x) for x in labels) - {-1}):
        idx = np.flatnonzero(labels == m)
        info["modes_in"] += 1
        if len(idx) < max(min_size, 2 * 2) or max_sub <= 1:
            out[idx] = nxt; nxt += 1; info["modes_out"] += 1
            continue
        d = _pairwise_rmsd(coords[idx])
        z = linkage(squareform(d, checks=False), method="average")
        # CUT BY DIAMETER, THEN CAP. `maxclust` always returns k clusters, so a
        # mode that is genuinely one tight cluster would be cut into `max_sub`
        # pieces and cost `max_sub` sweeps for nothing. A distance cut asks the
        # question that matters instead -- is anything in here further apart than
        # the accuracy we claim -- and leaves a tight mode alone. The cap then
        # bounds the cost when a mode really is that heterogeneous.
        sub = fcluster(z, CUT_A, criterion="distance")
        if len(np.unique(sub)) > max_sub:
            sub = fcluster(z, max_sub, criterion="maxclust")
        got = 0
        for s in np.unique(sub):
            members = idx[sub == s]
            # A sub-cluster of one pose is not a mode. Fold singletons back into
            # the nearest surviving sub-cluster rather than giving noise a row
            # and a sweep.
            if len(members) < 2 and len(np.unique(sub)) > 1:
                continue
            out[members] = nxt; nxt += 1; got += 1
        stray = idx[out[idx] < 0]
        if len(stray):
            # Whatever was folded out goes to the largest sub-cluster of its own
            # mode, chosen by size, so no pose is silently dropped.
            sizes = [(int((out[idx] == v).sum()), v) for v in set(out[idx]) if v >= 0]
            out[stray] = max(sizes)[1] if sizes else nxt
            if not sizes:
                nxt += 1; got += 1
        info["modes_out"] += got
        if got > 1:
            info["subdivided"] += 1
    info["max_sub"] = max_sub
    info["min_size"] = min_size
    return out, info
