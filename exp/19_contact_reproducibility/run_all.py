#!/usr/bin/env python3
"""
Purpose: do contact-space groups reproduce across INDEPENDENT dockings, not just deeper ones?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 5 independent replicate clouds under election_<candidate>_r{1..5}_allposes
Output: 00_outputs/blacksmith/contact_reproducibility/

THE TEST EVERY PREVIOUS CANDIDATE PASSED UNTIL IT DIDN'T. exp/17 showed groups
persist when the SAME cloud is sampled deeper -- 100% of n=500 groups found again
at n=6,000. That is a weaker claim than it sounds: a fixed tolerance carves fixed
regions, so a deeper draw from one cloud cannot move them. An INDEPENDENT docking
can. HDBSCAN looked excellent on within-cloud measures and then kept only 1 of 3
modes across replicates (D0088, #78), which is the whole reason this file exists.

MATCHED BY THE RULE THAT BUILT THE GROUPS. A group in r1 is "the same group" as
one in r2 when their contact-space centres are within the same tolerance that
defines membership. Using a looser match than the membership rule would count
neighbouring groups as reproduced.

THERE IS NO FAIR BASELINE ON THESE CLOUDS, AND THE ATTEMPT IS KEPT AS A GUARD.
A DBSCAN baseline was run and returned ONE mode in every replicate, which then
"reproduced" at 100% -- a rule that always answers "one mode" is perfectly
reproducible and perfectly useless. Worse, it is circular: these clouds have
ALREADY been DBSCAN-cleaned (D0093), so re-running DBSCAN on them asks whether a
filter agrees with itself. The baseline is therefore REFUSED rather than
reported, by a degeneracy check, and the comparison numbers are cited from the
experiments that measured them on raw clouds: HDBSCAN kept the validated pose in
27 of 30 replicates and only 1 of 3 modes survived (D0088, #78). A head-to-head
needs raw clouds and is blocked on the re-dock.

CAVEAT CARRIED IN THE OUTPUT, NOT ONLY HERE: these five clouds are
`<topic>_allposes` files, which hold only poses whose DBSCAN label survived
(D0093). Roughly 21% of each cloud is absent, and it is the scattered 21%. Both
methods are measured on the same filtered clouds so the COMPARISON is fair, but
the absolute rates are the easy case.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("contact-repro")


def _by_path(name: str, rel: str):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_M = _by_path("e17", "exp/17_contact_saturation/run_all.py")


def flat(T, w):
    W = np.sqrt(w / w.sum())[None, :, None]
    return (T * W).reshape(len(T), -1) / np.sqrt(T.shape[2])


def centres(V, lab):
    return np.array([V[np.flatnonzero(lab == k)].mean(0) for k in range(lab.max() + 1)])


def sizes(lab):
    return np.bincount(lab)


def dbscan_labels(xyz, eps=3.0, frac=0.05):
    """The shipped stage-1 rule, as a baseline: DBSCAN on the pose centroid."""
    from sklearn.cluster import DBSCAN
    return DBSCAN(eps=eps, min_samples=max(3, int(frac * len(xyz)))
                  ).fit_predict(xyz.mean(1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--min-size", type=int, default=5,
                    help="a group this big or bigger is one worth reproducing")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    res = _M.receptor_coords(_M.key_residues(a.residues))
    clouds, tol, w = {}, None, None
    for r in range(1, a.replicates + 1):
        p = (rp.BLACKSMITH / f"election_{a.candidate}_r{r}_allposes" / f"{a.candidate}.sdf")
        if not p.is_file():
            log.warning("r%d missing", r)
            continue
        xyz, meta = _M.load_sdf(p)
        if xyz is None:
            continue
        if w is None:
            rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
            w = pc.atom_weights(rmsf)
            tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
        clouds[r] = xyz
    if len(clouds) < 2:
        raise SystemExit("need at least two replicate clouds")
    n_at = {r: x.shape[1] for r, x in clouds.items()}
    assert len(set(n_at.values())) == 1, f"atom counts differ across replicates: {n_at}"
    log.info("%d replicates, %d atoms, tolerance %.3f A",
             len(clouds), list(n_at.values())[0], tol)

    G = {}
    for r, xyz in clouds.items():
        T = pc.contact_tensor(xyz, res)
        V = flat(T, w)
        lab = pc.group(pc.pose_distances(T, w), tol)
        G[r] = dict(V=V, lab=lab, cen=centres(V, lab), sz=sizes(lab), n=len(xyz),
                    dlab=dbscan_labels(xyz))
        log.info("  r%d: %d poses -> %d contact groups (%d with >=%d poses), "
                 "%d DBSCAN modes", r, len(xyz), lab.max() + 1,
                 int((G[r]["sz"] >= a.min_size).sum()), a.min_size,
                 len(set(G[r]["dlab"].tolist()) - {-1}))

    # ---- contact groups: pairwise and all-replicate ---------------------- #
    rows = []
    reps = sorted(G)
    for i, j in itertools.permutations(reps, 2):
        A, B = G[i], G[j]
        d = np.linalg.norm(A["cen"][:, None, :] - B["cen"][None, :, :], axis=2)
        hit = d.min(1) <= tol
        for msize in (1, a.min_size):
            keep = A["sz"] >= msize
            rows.append(dict(src=i, dst=j, min_size=msize, groups=int(keep.sum()),
                             matched=int(hit[keep].sum()),
                             rate=float(hit[keep].mean()) if keep.any() else np.nan))
    pw = pd.DataFrame(rows)

    core = []
    for msize in (1, a.min_size):
        r0 = reps[0]
        keep = np.flatnonzero(G[r0]["sz"] >= msize)
        allhit = np.ones(len(keep), bool)
        for j in reps[1:]:
            d = np.linalg.norm(G[r0]["cen"][keep][:, None, :] - G[j]["cen"][None, :, :],
                               axis=2)
            allhit &= (d.min(1) <= tol)
        core.append(dict(min_size=msize, groups=len(keep), in_all=int(allhit.sum()),
                         frac=float(allhit.mean()) if len(keep) else np.nan))
    co = pd.DataFrame(core)

    # ---- baseline: the shipped rule, matched on pose overlap ------------- #
    #: The shipped rule has no contact-space centre, so its modes are matched the
    #: way exp/8 matched them -- medoid Cartesian RMSD under the project's 2 A
    #: "same pose" bar. A different rule, named as such, not silently reused.
    def dmedoids(r):
        xyz, lab = clouds[r], G[r]["dlab"]
        out = []
        for k in sorted(set(lab.tolist()) - {-1}):
            m = xyz[lab == k]
            dd = np.array([np.sqrt(((m - m[t]) ** 2).sum(-1).mean(-1)).sum()
                           for t in range(len(m))])
            out.append(m[int(dd.argmin())])
        return out

    DM = {r: dmedoids(r) for r in reps}
    brows = []
    for i, j in itertools.permutations(reps, 2):
        if not DM[i] or not DM[j]:
            continue
        hit = [min(float(np.sqrt(((x - y) ** 2).sum(-1).mean())) for y in DM[j]) <= 2.0
               for x in DM[i]]
        brows.append(dict(src=i, dst=j, modes=len(hit), matched=int(sum(hit)),
                          rate=float(np.mean(hit))))
    bl = pd.DataFrame(brows)
    r0 = reps[0]
    b_core = 0
    for x in DM[r0]:
        if all(min(float(np.sqrt(((x - y) ** 2).sum(-1).mean())) for y in DM[j]) <= 2.0
               for j in reps[1:] if DM[j]):
            b_core += 1

    t = sout.Topic("blacksmith", "contact_reproducibility")
    pw.to_csv(t.write("pairwise", ".csv"), index=False)
    co.to_csv(t.write("core", ".csv"), index=False)
    if len(bl):
        bl.to_csv(t.write("baseline_dbscan", ".csv"), index=False)

    P = print
    P("\n" + "=" * 80)
    P(f"  DO CONTACT GROUPS REPRODUCE ACROSS INDEPENDENT DOCKINGS?  "
      f"{a.candidate}")
    P("=" * 80)
    P(f"\n  {len(reps)} independent 500-run dockings · tolerance {tol:.2f} A · "
      f"{a.residues} landmark residues")
    P(f"  poses per replicate: {', '.join(str(G[r]['n']) for r in reps)}")
    P(f"  contact groups:      {', '.join(str(G[r]['lab'].max() + 1) for r in reps)}")

    P("\n  PAIRWISE — a group in one replicate found in another")
    for msize in (1, a.min_size):
        s = pw[pw.min_size == msize]
        P(f"    groups of >= {msize:2d} poses: mean {s.rate.mean() * 100:5.1f}%  "
          f"range {s.rate.min() * 100:.0f}-{s.rate.max() * 100:.0f}%  "
          f"({int(s.groups.mean())} groups per replicate)")

    P(f"\n  THE REPRODUCIBLE CORE — present in ALL {len(reps)} replicates")
    for _, r in co.iterrows():
        P(f"    groups of >= {int(r.min_size):2d} poses: {r.in_all} of {r.groups} "
          f"({r.frac * 100:.0f}%)")

    n_modes = [len(DM[r]) for r in reps]
    P("\n  BASELINE — REFUSED, not reported")
    P(f"    DBSCAN returned {n_modes} modes across the five replicates.")
    if max(n_modes) <= 1:
        P("    A rule that always answers \"one mode\" reproduces at 100% and")
        P("    discriminates nothing, so its match rate is NOT a comparison.")
    P("    It is also circular here: these clouds were already DBSCAN-cleaned")
    P("    (D0093), so this asks whether a filter agrees with itself. For the")
    P("    real numbers see D0088/#78 — 1 of 3 modes survived an independent")
    P("    draw. A fair head-to-head needs raw clouds.")

    P("\n  CAVEAT: these clouds are `_allposes` files — DBSCAN-filtered, ~21% absent")
    P("  (D0093). Both methods see the same filter so the comparison is fair; the")
    P("  absolute rates are the easy case.")
    P("\n" + "=" * 80)
    P(f"  written to {t.dir}")
    P("=" * 80 + "\n")


if __name__ == "__main__":
    main()
