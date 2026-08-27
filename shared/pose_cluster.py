"""Modes by POSE SIMILARITY alone, ranked afterwards by attack geometry.

@tt8804: *"we first cluster into groups to make modes and then we rank the
modes"*, and *"there is only stage 2 -- what you are describing as stage 1
happens later"*, and *"use HDBSCAN"*.

WHAT THIS REPLACES, AND WHY THE OLD ORDER WAS CIRCULAR. `pose_modes.split`
clusters on the reactive atom's POSITION and the direction its warhead faces,
and `pose_subsplit.subdivide` then subdivides those on whole-molecule RMSD. The
comment on the first step says it clusters "never on the NAC geometry itself,
which is the score" -- but the score's distance term IS the distance from that
atom to Cys113's SG, which is fixed within a run, and its angle term IS that
direction. So the pipeline clustered on the score and then scored the clusters.

That circularity is what made a mode's `viable_fraction` behave like a mixing
ratio rather than a property: groups were formed along the very axis they were
later graded on, so a group's grade described where the cut happened to fall.

Here there is ONE clustering step and it knows nothing about the anchor: poses
are grouped by how similar the whole molecule's placement is, full stop. Attack
geometry is used later, to rank the groups (`mode_angle_score`), never to form
them.

WHY HDBSCAN AND NOT DBSCAN. DBSCAN takes a single radius `eps` and links any two
poses inside it, so a chain of poses each within `eps` of the next becomes one
group however wide it grows -- measured on this library, a group formed at
eps 3.0 spans 4.22 A. It also needs that radius chosen in advance, and the right
value moves with how densely the cloud was sampled: at a fixed `eps` the shipped
splitter collapsed a 2,000-pose cloud into a SINGLE group of 1,173 poses.

HDBSCAN asks for a minimum group SIZE instead of a radius. It builds the
hierarchy over all densities at once and keeps the groups that persist across
it, so a tight group inside a diffuse halo survives as its own group rather than
being absorbed, and nothing has to be retuned when sampling depth changes. The
halo is labelled noise (-1), which is the honest answer for a pose that belongs
to no repeated arrangement.

STATUS: SUPERSEDED. NO PRODUCTION CALLER -- only exp/4,5,6,7,8,9,10 and their
tests, which are the record of why it was not adopted: it discards 29% of a cloud
as noise, lost the MD-validated pose in 3 of 30 replicates, kept 1 of 3 modes
across an independent draw (#78), and its cluster count grows linearly with
sampling because HDBSCAN has no length scale (D0090). Replaced by
`shared/pose_contacts`. Listed in data/ready_to_delete.md; delete when those
experiments are archived. Do not wire this into anything new.
"""

from __future__ import annotations

import numpy as np

#: Poses a group needs to exist. @tt8804: "maybe as low as 3 poses" -- a mode is
#: a REPEATED arrangement, and three independent docking runs landing in the same
#: place is the smallest claim worth making. Below this a "mode" is one lucky
#: draw with a label.
MIN_CLUSTER_SIZE = 3

#: How conservative the density estimate is. Left at `min_cluster_size` (the
#: HDBSCAN default behaviour) so there is one knob rather than two, and the
#: number that decides what counts as a mode is the same number a reader sees.
MIN_SAMPLES = None

#: Which groups to keep from the hierarchy. MEASURED on a real 393-pose cloud
#: (t4_716800c125a7), median / widest pair inside a mode and the largest mode:
#:
#:   eom,  min_samples 3    49 modes  1.56 / 5.43 A  largest 20
#:   eom,  min_samples 10   12 modes  3.68 / 4.72 A  largest 64
#:   leaf, min_samples 3    55 modes  1.54 / 4.05 A  largest 14   <- this
#:
#: `leaf` takes the finest groups in the hierarchy rather than the most
#: persistent, which is what "smaller more homogeneous clusters of very similar
#: poses" asks for (@tt8804). `eom` at a high min_samples reproduces the old
#: failure -- a 64-pose group nearly 4 A across.
SELECTION = "leaf"


def rmsd_matrix(coords: np.ndarray) -> np.ndarray:
    """Pairwise heavy-atom RMSD, (n_poses, n_atoms, 3) -> (n_poses, n_poses).

    NO SUPERPOSITION, deliberately: every pose is already in the receptor's
    frame, and fitting them onto each other would ask whether they are the same
    SHAPE rather than whether they sit in the same PLACE. Place is the question.
    """
    n = len(coords)
    d = np.zeros((n, n), dtype=float)
    for i in range(n):
        d[i] = np.sqrt(((coords - coords[i]) ** 2).sum(axis=2).mean(axis=1))
    np.fill_diagonal(d, 0.0)
    return (d + d.T) / 2.0            # exact symmetry; HDBSCAN requires it


def cluster(coords: np.ndarray,
            min_cluster_size: int = MIN_CLUSTER_SIZE,
            min_samples: int | None = MIN_SAMPLES,
            selection: str = SELECTION) -> np.ndarray:
    """Mode label per pose; -1 is noise. Groups are renumbered by size.

    `selection="eom"` (excess of mass) keeps the most persistent groups in the
    hierarchy; `"leaf"` keeps the finest ones instead, which is the setting to
    reach for if modes come out too coarse. Neither takes a radius.

    Renumbered largest-first so mode 0 means the same thing across runs -- ids
    that follow scan order make two runs of one molecule disagree about which
    group is "mode 0" for no reason anyone can see.
    """
    from sklearn.cluster import HDBSCAN

    n = len(coords)
    if n == 0:
        return np.empty(0, dtype=int)
    if n < max(2, min_cluster_size):
        return np.full(n, -1, dtype=int)

    d = rmsd_matrix(coords)
    # copy=True: HDBSCAN mutates a precomputed matrix in place otherwise, and
    # callers reuse it (the widest-pair diagnostics below, for one).
    lab = HDBSCAN(min_cluster_size=int(min_cluster_size),
                  min_samples=(int(min_samples) if min_samples else None),
                  metric="precomputed",
                  cluster_selection_method=selection,
                  copy=True).fit_predict(d)

    order = [c for c, _ in sorted(
        ((c, int((lab == c).sum())) for c in set(lab) - {-1}),
        key=lambda kv: -kv[1])]
    remap = {c: i for i, c in enumerate(order)}
    return np.array([remap.get(int(x), -1) for x in lab], dtype=int)


def mode_angle_score(labels: np.ndarray, angle_deg: np.ndarray,
                     in_range: np.ndarray, angle_max_deg: float,
                     consensus_w: float = 0.0) -> dict[int, float]:
    """Rank a mode by ATTACK ANGLE. Higher is better; 1.0 is dead-on.

    @tt8804: *"we just rank by attack angle ... and we heavily decrease how much
    we rank by consensus."*

    The angle is a property of the pose, so it neither rewards a mode for being
    large nor punishes it for being small. That matters more than it sounds:
    mode count grows with how deeply the cloud was sampled, so `consensus`
    (mode_size / n_poses) is partly a statement about `docking.n_runs`. Ranking
    on it makes the shortlist depend on the sampling budget.

    Only poses at a reactable DISTANCE contribute -- the approach angle of a pose
    8 A away is not an attack angle. A mode with none of those is unscorable and
    is absent from the result rather than scored 0, because "never got close" and
    "got close at a bad angle" are different facts.

    `consensus_w` keeps population available as a dial: 0 ignores it entirely.
    """
    labels = np.asarray(labels)
    ang = np.asarray(angle_deg, dtype=float)
    inr = np.asarray(in_range, dtype=bool)
    if not (len(labels) == len(ang) == len(inr)):
        raise ValueError("labels, angle_deg and in_range must describe the same "
                         f"poses; got {len(labels)}, {len(ang)}, {len(inr)}")
    n_tot = max(len(labels), 1)
    out: dict[int, float] = {}
    for c in sorted({int(x) for x in labels if x >= 0}):
        m = labels == c
        sel = m & inr
        if not sel.any():
            continue
        med = float(np.median(ang[sel]))
        quality = max(0.0, (angle_max_deg - med) / angle_max_deg)
        if consensus_w:
            quality *= (float(m.sum()) / n_tot) ** consensus_w
        out[c] = quality
    return out


def _renumber(lab: np.ndarray) -> np.ndarray:
    """Largest group becomes 0, next 1, ...; -1 stays noise.

    Ids that follow scan order make two runs of one molecule disagree about
    which group is "mode 0" for no reason a reader can see.
    """
    order = [c for c, _ in sorted(
        ((c, int((lab == c).sum())) for c in set(int(x) for x in lab) - {-1}),
        key=lambda kv: -kv[1])]
    remap = {c: i for i, c in enumerate(order)}
    return np.array([remap.get(int(x), -1) for x in lab], dtype=int)


def cluster_coords(coords: np.ndarray,
                   min_cluster_size: int = MIN_CLUSTER_SIZE,
                   min_samples: int | None = MIN_SAMPLES,
                   selection: str = SELECTION) -> np.ndarray:
    """HDBSCAN on the raw coordinates: (n_poses, n_atoms, 3) -> label per pose.

    @tt8804: *"use HDBSCAN on only the molecules dimensions in 3d space
    (3 x atoms dimensions) so that we generate clusters that are essentially the
    same poses being recreated."*

    Each pose is ONE POINT in 3N-dimensional space -- x1,y1,z1,x2,y2,z2,... for
    every heavy atom, in the receptor's frame, with no superposition and nothing
    about the anchor. Two poses are near each other exactly when every atom of
    one sits near the same atom of the other, which is what "the same pose being
    recreated" means. HDBSCAN then keeps the groups that persist across
    densities and labels the rest noise, so a pose that belongs to no repeated
    arrangement is an ORPHAN rather than being forced into the nearest mode.

    SAME PARTITION AS `cluster()`, AND NOT BY COINCIDENCE. Euclidean distance in
    3N space is `sqrt(n_atoms)` times the un-superposed RMSD `cluster()` feeds in
    precomputed -- one constant factor for a given molecule. HDBSCAN's hierarchy,
    its mutual-reachability graph and its excess-of-mass / leaf selection are all
    invariant to a global rescale of the distances (nothing here sets
    `cluster_selection_epsilon`, which is the one parameter that would not be),
    so the labels are identical. `tests/test_pose_cluster_coords.py` asserts it
    on random clouds rather than leaving it as an argument.

    THE DIFFERENCE IS COST, NOT ANSWER. `cluster()` materialises an
    (n_poses x n_poses) matrix built in Python over every atom; this hands
    sklearn the (n_poses x 3*n_atoms) array and lets it index. On a 500-pose,
    50-heavy-atom cloud that is ~1,500 floats per pose against a 250,000-cell
    matrix, and it is what makes a 561-molecule library tractable in one pass.

    UNITS. Distances here are Angstrom-of-whole-molecule-displacement, NOT RMSD:
    a group whose Euclidean diameter is `d` has RMSD diameter `d / sqrt(n_atoms)`.
    `min_cluster_size` is a COUNT, so it is unaffected -- but anything that reads
    a distance out of this space has to divide, and `exp/7_coord_modes` reports
    widths in RMSD for exactly that reason.
    """
    from sklearn.cluster import HDBSCAN

    n = len(coords)
    if n == 0:
        return np.empty(0, dtype=int)
    if n < max(2, min_cluster_size):
        return np.full(n, -1, dtype=int)
    x = np.asarray(coords, dtype=float).reshape(n, -1)
    lab = HDBSCAN(min_cluster_size=int(min_cluster_size),
                  min_samples=(int(min_samples) if min_samples else None),
                  metric="euclidean",
                  # Explicit because sklearn 1.10 flips the default and warns
                  # until it does. `x` is this function's own reshape, so a copy
                  # costs nothing and no caller's array can be mutated.
                  copy=True,
                  cluster_selection_method=selection).fit_predict(x)
    return _renumber(lab)
