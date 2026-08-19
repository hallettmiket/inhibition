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

#: Poses a SUB-cluster needs to earn its own row (@tt8804: "why would we even
#: have a mode cap. we should just establish a min consensus score for a mode
#: (maybe as low as 3 poses)").
#:
#: The right control on a split is how much evidence a sub-mode carries, not how
#: many sub-modes came out. `max_sub` bounds the count, so once a parent holds
#: more distinct poses than the cap the surplus is merged back into whichever
#: cluster is nearest or largest -- widening exactly the modes the cut had just
#: separated. A minimum size bounds the same cost without that side effect:
#: clusters below it are folded, everything above it keeps its own row however
#: many there are.
#:
#: Set `max_sub=None` to let the diameter cut govern with this as the only
#: guard. The default stays 5 until the pairing is measured -- see
#: exp/2_mode_homogeneity.
MIN_SUB_SIZE = 3

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
              max_sub: int | None = DEFAULT_MAX_SUB,
              min_size: int = MIN_MODE_SIZE,
              min_sub_size: int = MIN_SUB_SIZE,
              cut_a: float = CUT_A) -> tuple[np.ndarray, dict]:
    """Split each mode into at most `max_sub` sub-modes on whole-molecule RMSD.

    `labels` is the first-stage mode assignment (-1 = unassigned, preserved).
    `coords` is (n_poses, n_heavy_atoms, 3) in the receptor frame.

    Returns (new_labels, info). New labels are renumbered contiguously from 0, so
    everything downstream keeps treating them as ordinary modes -- the mode
    abstraction gets finer, nothing else changes shape. Unassigned poses stay -1.

    `info["parent"]` maps each new mode to (first_stage_mode, sub_index), and
    `info["label"]` to its display name: `1a`, `1b`. WITHOUT THAT MAPPING THE
    PROVENANCE IS GONE -- a molecule showing m0..m4 would be indistinguishable
    from one with five genuine first-stage modes, when it may be one mode split
    five ways, and those are different claims about the pose cloud. The ident
    stays numeric (`_m3`) so `mode_key` and every join keep working; the letter
    is a display label carried alongside.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    labels = np.asarray(labels)
    out = np.full(len(labels), -1, dtype=int)
    info: dict = {"modes_in": 0, "modes_out": 0, "subdivided": 0,
                  "parent": {}, "label": {}}
    nxt = 0

    def _claim(parent: int, sub_i: int) -> int:
        nonlocal nxt
        new = nxt; nxt += 1
        info["parent"][new] = (int(parent), int(sub_i))
        info["label"][new] = f"{parent}{chr(ord('a') + sub_i)}" if sub_i >= 0 \
            else str(parent)
        return new
    for m in sorted(set(int(x) for x in labels) - {-1}):
        idx = np.flatnonzero(labels == m)
        info["modes_in"] += 1
        if len(idx) < max(min_size, 2 * 2) or (max_sub is not None and max_sub <= 1):
            out[idx] = _claim(m, -1); info["modes_out"] += 1
            continue
        d = _pairwise_rmsd(coords[idx])
        z = linkage(squareform(d, checks=False), method="average")
        # CUT BY DIAMETER, THEN CAP. `maxclust` always returns k clusters, so a
        # mode that is genuinely one tight cluster would be cut into `max_sub`
        # pieces and cost `max_sub` sweeps for nothing. A distance cut asks the
        # question that matters instead -- is anything in here further apart than
        # the accuracy we claim -- and leaves a tight mode alone. The cap then
        # bounds the cost when a mode really is that heterogeneous.
        sub = fcluster(z, cut_a, criterion="distance")
        # THE CAP IS OPTIONAL. With max_sub set this re-cuts at k=max_sub, which
        # discards the diameter answer and merges poses the cut separated; that
        # is the behaviour `min_sub_size` exists to replace. With max_sub=None
        # the cut stands and size is the only guard.
        if max_sub is not None and len(np.unique(sub)) > max_sub:
            sub = fcluster(z, max_sub, criterion="maxclust")
        got = 0
        keep_n = {s: int((sub == s).sum()) for s in np.unique(sub)}
        floor = max(2, min_sub_size) if max_sub is None else 2
        for s in np.unique(sub):
            members = idx[sub == s]
            # A sub-cluster below the evidence floor is not a mode. Fold it back
            # rather than giving noise its own row and its own sweep.
            if keep_n[s] < floor and len(np.unique(sub)) > 1:
                continue
            out[members] = _claim(m, got); got += 1
        stray = idx[out[idx] < 0]
        if len(stray):
            # Whatever was folded out goes to the largest sub-cluster of its own
            # mode, chosen by size, so no pose is silently dropped.
            sizes = [(int((out[idx] == v).sum()), v) for v in set(out[idx]) if v >= 0]
            out[stray] = max(sizes)[1] if sizes else _claim(m, got)
            if not sizes:
                got += 1
        info["modes_out"] += got
        if got > 1:
            info["subdivided"] += 1
        elif got == 1:
            only = int(out[idx][0])
            info["parent"][only] = (int(m), -1)
            info["label"][only] = str(m)
    info["max_sub"] = max_sub
    info["min_sub_size"] = min_sub_size
    info["cut_a"] = cut_a
    info["min_size"] = min_size
    return out, info
