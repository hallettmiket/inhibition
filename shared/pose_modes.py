"""
Purpose: split a molecule's docked poses into binding modes, each of which becomes a candidate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: a meeko-rebuilt multi-conformer mol + its reactive SMARTS match
Output: a mode label per pose, plus a geometric identity per mode

The 2.2.0 core feature (`docs/build_plan_2.2.0.md` §1.2). @tt8804: *"within the
200 to split/group poses by consensus and we treat all those split poses as
separate candidates and we filter by our current methods."*

WHY THIS EXISTS. #30 measured that the crystallographically correct pose is
present in a 200-run docking 93.3% of the time and reaches the scoring window
33.3% of the time, because the window is the 20 lowest-energy poses and energy
places the correct pose at a rank indistinguishable from uniform (KS p = 0.666).
The information is generated and then discarded. Splitting replaces "keep the 20
best-scoring poses" with "keep a representative of every distinct thing the
molecule does".

WHAT IS CLUSTERED, AND WHAT IS DELIBERATELY NOT.

  * Clustered on the REACTIVE ATOM'S POSITION and the direction its warhead
    faces. D0062: whole-molecule RMSD is the wrong endpoint for a covalent
    question -- two poses that place the warhead identically and differ in a
    distal ring are one mode for our purposes, and RMSD would call them two.
  * NOT clustered on docking energy, at all. #23/#30. Using energy to define or
    order modes would re-import the exact defect this version exists to remove.
    Energy is reported per mode; it never shapes one.
  * NOT clustered on distance-to-SG or the NAC angle. Those are the SCORE. A
    mode defined partly by its own score cannot then be scored honestly -- it
    would be guaranteed to look internally consistent.

EPS IS A RESOLUTION, AND IT CARRIES ITS JUSTIFICATION (D0068). Every number in
this project has to say where it came from. `eps` does not set the number of
modes -- that is measured -- but it does set what counts as "the same place",
so it is calibrated against crystal ground truth and reported with its
stability, never chosen for appearance.

STATUS: LIVE, and the incumbent -- `config/target.yaml: splitting.method` is
`warhead_dbscan`, which is this module plus `pose_subsplit`. It carries a known
defect (D0088): it clusters on the reactive atom's POSITION and the direction its
warhead faces, then each group is graded by how often it reaches attack geometry
-- which is position and direction. It forms groups along the axis it grades them
on. The median mode spans 3.51 A and 42% hold two populations under one label.
`shared/pose_contacts.split_poses` is the built replacement; the switch waits on
the re-screen (#79). Five modules in this repo group poses -- see
docs/pose_frameworks.md before adding a sixth.
"""

from __future__ import annotations

import numpy as np

#: Angstroms of positional difference treated as equivalent to one radian of
#: orientational difference. Two poses whose warheads sit 1 A apart but point
#: 30 degrees differently are then about equally distinct either way. Without a
#: conversion the position term dominates and orientation stops mattering --
#: which would merge a productive approach with an inverted one sitting in the
#: same spot, and those are chemically opposite.
ANGSTROM_PER_RADIAN = 2.0

#: A mode holding fewer than this fraction of a molecule's poses is LABELLED
#: noise, not deleted. Same rule the consensus floor already follows: a rare
#: geometry might be a real minor mode, and silently dropping it would make the
#: mode count depend on how hard we happened to sample.
#:
#: CALIBRATED, NOT CHOSEN (D0068). Swept against 15 crystal complexes docked
#: TWICE at 500 runs each, scoring reproducibility across the two independent
#: dockings:
#:
#:   min_pop  modes/mol  count reproduces  mode-0 reproduces  crystal in mode 0
#:      0.02        2.9               47%                80%                93%
#:      0.05        2.1               73%                87%                87%
#:      0.10        1.4               87%                73%                87%
#:      0.20        1.1               93%                60%                60%
#:
#: 0.05 is the knee. Below it the mode COUNT is barely reproducible (47%), which
#: is the documented failure condition -- "a mode is then an artefact of the
#: clustering, not a property of the ligand". Above it the clustering collapses
#: toward one mode per molecule and stops being pose splitting at all; by 0.20
#: it is degenerate and accuracy collapses with it. The one point of crystal
#: accuracy given up between 0.02 and 0.05 buys 26 points of count
#: reproducibility, and a mode that does not reproduce is not worth scoring.
MIN_POPULATION_FRAC = 0.05

#: Clustering resolution, in the same units as `distances` (Angstrom-equivalents).
#:
#: CALIBRATED, NOT CHOSEN (D0068). Swept against crystal ground truth on the same
#: 15 complexes at 500 runs:
#:
#:   eps  modes/mol  crystal in a NAMED mode  crystal in the TOP mode
#:   1.0        5.5                     60%                      33%
#:   2.0        3.9                    100%                      87%
#:   3.0        2.9                    100%                      93%
#:   4.0        1.5                    100%                     100%
#:
#: 3.0 is adopted. The apparent perfection at 4.0-5.0 is DEGENERATE: at 1.1-1.5
#: modes per molecule there is essentially one cluster, so "the crystal pose is
#: in the top mode" is a tautology rather than a result. 3.0 keeps genuine
#: multi-modality (2.9 modes/molecule) while placing the crystal pose in a named
#: mode every time.
DEFAULT_EPS = 3.0


def features(mol, smarts_match: tuple[int, ...]) -> np.ndarray:
    """(n_poses, 6): reactive-atom xyz, then the unit vector its warhead faces.

    The direction runs from the reactive atom toward the centroid of the rest of
    the SMARTS match -- the leaving group for SN2, the sp2 system for a Michael
    acceptor. It is mechanism-agnostic on purpose: it asks "which way is the
    warhead pointing", which is the thing that distinguishes a productive
    approach from an inverted one, without needing a per-mechanism special case
    that could silently disagree with `nac_criterion`.
    """
    rx = int(smarts_match[0])
    rest = [int(i) for i in smarts_match[1:]]
    out = np.empty((mol.GetNumConformers(), 6))
    for cid in range(mol.GetNumConformers()):
        pos = mol.GetConformer(cid).GetPositions()
        p = pos[rx]
        v = (pos[rest].mean(axis=0) - p) if rest else np.array([0.0, 0.0, 1.0])
        n = np.linalg.norm(v)
        out[cid, :3] = p
        out[cid, 3:] = v / n if n > 1e-9 else np.array([0.0, 0.0, 1.0])
    return out


def distances(feat: np.ndarray) -> np.ndarray:
    """Pairwise pose distance: positional separation + scaled angular difference."""
    pos, dirs = feat[:, :3], feat[:, 3:]
    dpos = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
    cos = np.clip(dirs @ dirs.T, -1.0, 1.0)
    dang = np.arccos(cos)
    return dpos + ANGSTROM_PER_RADIAN * dang


def _dbscan(dist: np.ndarray, eps: float, min_samples: int) -> np.ndarray:
    """DBSCAN on a precomputed distance matrix. -1 is noise.

    Written out rather than imported from scikit-learn ON PURPOSE. This module is
    called by `nac_screen_v2`, which runs in the `dwi_reactive` environment --
    and sklearn is not installed there. A shared module on the screen's hot path
    must not depend on a package the screen's own interpreter lacks; that is an
    ImportError three hours into a library re-dock. The algorithm is thirty
    lines and has no tuning of its own beyond the two arguments already being
    justified above.
    """
    n = len(dist)
    neigh = [np.flatnonzero(dist[i] <= eps) for i in range(n)]
    core = np.array([len(neigh[i]) >= min_samples for i in range(n)])
    labels = np.full(n, -1, dtype=int)
    cid = 0
    for i in range(n):
        if labels[i] != -1 or not core[i]:
            continue
        labels[i] = cid
        queue = list(neigh[i])
        while queue:
            j = queue.pop()
            if labels[j] == -1:
                labels[j] = cid          # reachable: joins, core or not
                if core[j]:
                    queue.extend(neigh[j])
        cid += 1
    return labels


def _complete_linkage(dist: np.ndarray, diameter: float,
                      min_samples: int) -> np.ndarray:
    """Cluster so NO two members are further apart than `diameter`. -1 is noise.

    The alternative to DBSCAN, and the difference is the whole point (D0086).
    DBSCAN bounds the LINK: A-B-C-D each within `eps` of the next is one cluster
    however wide the chain grows. Complete linkage bounds the DIAMETER: the
    merge criterion is the FURTHEST pair, so a cluster is only formed if every
    member is within the tolerance of every other. Same number, enforced.

    Clusters below `min_samples` are labelled noise, matching DBSCAN's treatment
    so the two are swappable without changing what "unassigned" means.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform
    n = len(dist)
    if n < 2:
        return np.zeros(n, dtype=int) if n else np.empty(0, dtype=int)
    z = linkage(squareform(dist, checks=False), method="complete")
    lab = fcluster(z, diameter, criterion="distance")
    out = np.full(n, -1, dtype=int)
    for c in np.unique(lab):
        m = lab == c
        if int(m.sum()) >= min_samples:
            out[m] = int(c)
    return out


def split(feat: np.ndarray, eps: float = DEFAULT_EPS,
          min_population_frac: float = MIN_POPULATION_FRAC,
          method: str = "dbscan") -> np.ndarray:
    """Mode label per pose; -1 means noise. Count is measured, not requested.

    Density-based, so the number of modes falls out of the pose cloud rather
    than being asked for. `min_samples` scales with how many poses the molecule
    has, so the same rule applies whether it was docked 200 times or 2,000 -- a
    fixed count would make a molecule look more multi-modal simply for having
    been sampled harder.

    `method="complete"` swaps the density rule for complete linkage at the same
    `eps`, read as a DIAMETER rather than a neighbour radius -- see D0086 and
    `_complete_linkage`. The default is unchanged until the comparison in
    exp/3_linkage is decided.
    """
    n = len(feat)
    if n == 0:
        return np.empty(0, dtype=int)
    min_samples = max(3, int(round(min_population_frac * n)))
    if method == "complete":
        lab = _complete_linkage(distances(feat), eps, min_samples)
    elif method == "dbscan":
        lab = _dbscan(distances(feat), eps, min_samples)
    else:
        raise ValueError(f"unknown split method {method!r}; "
                         "expected 'dbscan' or 'complete'")
    # Relabel so 0 is the most populated mode, 1 the next, and so on. Without
    # this, mode ids are an artefact of scan order and two runs of the same
    # molecule would disagree about which mode is "mode 0".
    order = [c for c, _ in sorted(
        ((c, int((lab == c).sum())) for c in set(lab) - {-1}),
        key=lambda kv: -kv[1])]
    remap = {c: i for i, c in enumerate(order)}
    return np.array([remap.get(l, -1) for l in lab])


def identity(feat: np.ndarray, labels: np.ndarray, mode: int) -> dict:
    """A mode's geometric fingerprint — what makes it recognisable across runs.

    Named by WHERE the warhead sits and WHICH WAY it points, never by rank or
    by index. A mode named by its rank is a mode that a re-dock silently
    redefines, and every stability claim about modes then means nothing.
    """
    m = feat[labels == mode]
    if not len(m):
        return {}
    c = m[:, :3].mean(axis=0)
    d = m[:, 3:].mean(axis=0)
    nd = np.linalg.norm(d)
    d = d / nd if nd > 1e-9 else np.array([0.0, 0.0, 1.0])
    # NOT "n_poses". The screen merges this dict into a row that already has
    # `n_poses` meaning "poses this molecule produced"; the same key here means
    # "poses in this mode", and update() silently overwrote the total with the
    # mode count. One key, two meanings, no error -- the aggregate row reported
    # 468 of 468 for a molecule with 500 poses.
    return {"mode_size": int(len(m)),
            "centroid_x": float(c[0]), "centroid_y": float(c[1]),
            "centroid_z": float(c[2]),
            "dir_x": float(d[0]), "dir_y": float(d[1]), "dir_z": float(d[2]),
            "spread_a": float(np.linalg.norm(m[:, :3] - c, axis=1).mean()),
            # How tightly the mode agrees on orientation: the resultant length
            # of the mean of unit vectors, 1.0 = perfectly aligned, 0 = uniformly
            # scattered. NOT divided by n -- `nd` is already the norm of a MEAN,
            # so dividing again drove every mode to ~0.00 and made the field
            # look measured while carrying nothing.
            "dir_coherence": float(nd)}


def match_modes(id_a: dict, id_b: dict, tol_a: float = 2.0,
                tol_deg: float = 40.0) -> bool:
    """Are two modes, from independent dockings, the same mode?

    The basis of the stability test that `docs/build_plan_2.2.0.md` §1.2 makes an
    acceptance criterion: if mode membership is not reproducible across re-docks,
    mode-level scores inherit that noise and a "mode" is an artefact of the
    clustering rather than a property of the ligand.
    """
    if not id_a or not id_b:
        return False
    ca = np.array([id_a["centroid_x"], id_a["centroid_y"], id_a["centroid_z"]])
    cb = np.array([id_b["centroid_x"], id_b["centroid_y"], id_b["centroid_z"]])
    da = np.array([id_a["dir_x"], id_a["dir_y"], id_a["dir_z"]])
    db = np.array([id_b["dir_x"], id_b["dir_y"], id_b["dir_z"]])
    ang = np.degrees(np.arccos(np.clip(float(da @ db), -1.0, 1.0)))
    return bool(np.linalg.norm(ca - cb) <= tol_a and ang <= tol_deg)
