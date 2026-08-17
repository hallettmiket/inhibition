#!/usr/bin/env python3
"""
Purpose: write each mode's representative pose from a run's persisted cloud.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: --topic <nac_v4> [--cloud-dir <dir>] [--jobs N]
Output: <topic>_poses/<ident>.sdf, one conformer per mode

WHY THIS EXISTS (2026-08-12). `nac_screen_v2` derived its table topic from
`--topic` but hardcoded its two pose directories. The 3.0.0 run therefore wrote
tables to `nac_v4` and aimed its poses at `nac_v3_poses`, where every molecule
already had a file from Aug 07; the `if not sdf.exists()` append-only guard --
behaving correctly -- skipped all of them. The screen exited 0 having written no
representatives. The cause is fixed at source (`topic_paths`, tests in
tests/test_topic_paths.py); this recovers the run that already happened.

RECOVERY IS POSSIBLE BECAUSE #44 HELD. `persist_all_poses` wrote the full cloud
for 5,420 molecules, median 431 poses, each stamped with the `mode` the screen
assigned it -- post-sub-split, so `1a`/`1b` are already separate labels. Those
stamps agree with the tables for 5,336 of 5,336 molecules checked and disagree
for none. The representative is a pure function of the cloud plus the per-pose
geometry, both of which are on disk, so NO RE-DOCKING IS NEEDED.

THE SELECTION RULE IS IMPORTED, NOT REIMPLEMENTED. `representative_indices` in
`nac_screen_v2` is the one definition of "the typical pose of the well-anchored
quartile"; a copy here would be a second definition free to drift, which is the
same class of defect that caused the problem being repaired.

NOTHING IS RE-MEASURED. The first attempt recomputed each pose's approach
geometry from the cloud coordinates, which needs the position of Cys113's SG --
and that is NOT recoverable. The sidechain is flexible during docking, so
`sg_position` read the sulfur from the DLG's first docked model, and the DLG dies
with its work directory. Reading SG from the rigid receptor instead measures the
approach to where the sulfur STARTED, which is a different number.

So the write permutation is undone instead (`cloud_order`). The screen wrote the
cloud in `argsort(labels)` order and recorded every pose's geometry against its
`pose_idx`; reconstructing that permutation pairs each cloud pose with the row
the screen wrote for it. Every per-pose number used here is therefore the
screen's own, and the reconstruction is PROVEN per pose, not assumed: the mode
the recovered order implies must equal the `mode` each pose actually carries.

THE REACTIVE ATOM IS PROVEN TOO. Clustering features are built on the
reactive-atom index, and where a molecule has several SMARTS matches the screen
picked the docked one from the DLG. With the DLG gone, each candidate is solved
for the sulfur position its recorded distances imply (`trilaterate`): the true
atom fits one point, a wrong centre fits none. A molecule that cannot be pinned
to exactly one match is refused and counted, never written from a guess.
"""

from __future__ import annotations

import argparse
import csv
import glob
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import target_config as tc        # noqa: E402

log = logging.getLogger("rebuild-reps")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
#: Worst allowed residual, Angstrom, when a candidate reactive atom's recorded
#: distances are solved for a single sulfur position. The docked atom fits far
#: inside this; a wrong reactive centre has no consistent solution at all, so the
#: test does not depend on the exact value.
SG_FIT_TOL = 0.05


def per_pose_table(topic: str) -> dict:
    """{ident: {pose_idx: (mode, distance, angle)}} -- the screen's own geometry.

    `pose_idx` IS the conformer index of the mol the screen measured, so this
    table addresses poses the same way the cloud does once the write permutation
    is undone (see `cloud_order`). Nothing here is recomputed.

    THE JOIN KEY IS NOT `energy_rank`. `write_sdf` stamps that property with the
    WRITE SEQUENCE -- 1..n in mode order -- so every cloud file carries 1,2,3,...
    and the property identifies nothing despite its name. Joining on it silently
    pairs each pose with an unrelated row.
    """
    out: dict = defaultdict(dict)
    files = sorted(glob.glob(str(B / topic / "poses_s*_*.csv")))
    if not files:
        raise SystemExit(f"no poses_s*.csv in {B / topic}")
    for f in files:
        with open(f) as fh:
            for r in csv.DictReader(fh):
                try:
                    a = r["angle"]
                    out[r["ident"]][int(float(r["pose_idx"]))] = (
                        int(float(r["mode"])), float(r["distance"]),
                        float(a) if a not in ("", "None", "nan") else np.nan)
                except (ValueError, KeyError):
                    continue
    log.info("per-pose geometry for %d molecules from %d files", len(out), len(files))
    return out


def cloud_order(labels: np.ndarray, mode_ids: list) -> np.ndarray:
    """The pose indices the cloud holds, in the order it holds them.

    Reproduces `write_sdf`'s argument exactly -- the screen wrote

        order = [i for i in np.argsort(labels, kind="stable") if labels[i] in mode_ids]

    so cloud conformer j is original pose `order[j]`. Recovering this is what
    makes the repair exact rather than approximate: with it, every per-pose
    number comes from the screen's own measurement instead of a re-measurement
    that would need the flexible Cys113 SG, and THAT is unrecoverable -- the
    sidechain moves during docking, `sg_position` read it from the DLG, and the
    DLG is deleted with its work directory.
    """
    return np.array([i for i in np.argsort(labels, kind="stable")
                     if labels[i] in mode_ids], dtype=int)


def trilaterate(pts: np.ndarray, d: np.ndarray):
    """(the point at distance d[i] from pts[i], worst residual).

    Used only to decide WHICH reactive atom was docked when a molecule has
    several SMARTS matches. The screen resolved that from the DLG; with the DLG
    gone, the test is that an atom's positions and the recorded distances agree
    on a single sulfur location. The true reactive atom fits to ~0; a wrong one
    is angstroms out and has no consistent solution.
    """
    p0, d0 = pts[0], d[0]
    A = 2.0 * (p0[None, :] - pts[1:])
    b = ((p0 ** 2).sum() - d0 ** 2) - ((pts[1:] ** 2).sum(axis=1) - d[1:] ** 2)
    s, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = float(np.nanmax(np.abs(np.linalg.norm(pts - s, axis=1) - d)))
    return s, resid


def table_modes(topic: str) -> dict:
    """{ident: [mode, ...]} exactly as the run's aggregate rows describe them.

    The aggregate is the authority on which modes exist, because it is what the
    ranking reads. A cloud that offers a different set is a cloud from a
    different run, and writing a representative from it would recreate the very
    table/pose mismatch this script exists to repair.
    """
    out: dict = defaultdict(list)
    for f in sorted(glob.glob(str(B / topic / "agg_s*_*.csv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh):
                if r.get("status") != "ok":
                    continue
                try:
                    out[r["parent_ident"]].append(int(float(r["mode"])))
                except (ValueError, KeyError):
                    continue
    return {k: sorted(v) for k, v in out.items()}


def load_cloud(path: Path):
    """(mol with every cloud pose as a conformer, [energy_rank], [mode]).

    SANITIZED ONCE, ON THE CONTAINER, AND THAT IS NOT OPTIONAL. Poses are read
    with `sanitize=False` -- docked geometry regularly trips valence perception
    on the way in -- but the molecule is then sanitized before any substructure
    search. Without it, aromaticity is never perceived, and an aromatic reactive
    SMARTS matches NOTHING: `[c]([Cl])[n]` found 0 matches on every
    snar_chloroazine cloud and 2 on the same molecule once sanitized. That cost
    165 molecules to a "reactive SMARTS does not match" refusal that was a
    property of how the file was read, not of the molecule.

    Sanitizing does not reorder atoms, so the SMARTS indices still address the
    atoms the conformers hold -- which is the invariant the whole rebuild rests
    on.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    supp = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    base, eranks, modes = None, [], []
    for m in supp:
        if m is None:
            continue
        if not (m.HasProp("energy_rank") and m.HasProp("mode")):
            continue
        eranks.append(int(m.GetProp("energy_rank")))
        modes.append(int(m.GetProp("mode")))
        if base is None:
            base = Chem.Mol(m)
            base.RemoveAllConformers()
        conf = Chem.Conformer(m.GetConformer())
        base.AddConformer(conf, assignId=True)
    if base is None:
        raise ValueError("no usable poses")
    try:
        Chem.SanitizeMol(base)
    except Exception as e:                                 # noqa: BLE001
        raise ValueError(f"cloud will not sanitize: {str(e)[:40]}") from e
    return base, np.array(eranks), np.array(modes)


def resolve_match(mol, smarts: str, dist: np.ndarray):
    """The SMARTS match that was actually docked, proven against the record.

    Where a molecule has one match this confirms it; where it has several -- two
    reactive centres, only one of which was docked -- the screen resolved it from
    the DLG, which no longer exists. So each candidate is TESTED by
    `trilaterate`: does this atom's trajectory across the cloud, together with
    the distances the screen recorded, place the sulfur somewhere consistent?

    A molecule that cannot be pinned to exactly one match is refused and counted,
    never written from the first match -- taking the first is precisely what
    silently picks the wrong carbon.
    """
    from rdkit import Chem
    patt = Chem.MolFromSmarts(smarts)
    if patt is None:
        raise ValueError(f"unparseable reactive SMARTS {smarts!r}")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        raise ValueError("reactive SMARTS does not match the cloud mol")
    pos = np.array([mol.GetConformer(c).GetPositions()
                    for c in range(mol.GetNumConformers())])
    ok = np.isfinite(dist)
    if ok.sum() < 6:
        raise ValueError("too few measured poses to pin the reactive atom")
    good = []
    for mt in matches:
        _, resid = trilaterate(pos[ok, int(mt[0]), :], dist[ok])
        if resid <= SG_FIT_TOL:
            good.append((resid, mt))
    if not good:
        raise ValueError(f"no match fits the recorded distances "
                         f"({len(matches)} candidates)")
    if len(good) == 1:
        return good[0][1]
    # ORDER MATTERS: the screen takes `hits[0]` in GetSubstructMatches order, so
    # this must too. Sorting by fit residual first would pick a different member
    # of an equally-valid set and put the representative in a different place.
    good = [(r, mt) for mt in matches for r, m2 in good if m2 == mt]

    # SEVERAL MATCHES SHARING ONE REACTIVE ATOM IS NOT AMBIGUITY -- the same rule
    # `nac_screen.rebuild_and_match` applies, and it has to be the same rule or
    # this rebuild disagrees with the run it is repairing. A chloroazine's ipso
    # carbon sits between two ring nitrogens, so `[c]([Cl])[n]` matches twice
    # with the same attacked atom and a different nitrogen; both triples are
    # coplanar with the ring, so they define the same plane and the same
    # criterion. Refusing them cost 170 molecules here and both SNAr positives on
    # an earlier run.
    #
    # Verified, not assumed. If two matches at one reactive atom ever define
    # genuinely different planes, the criterion means two different things and
    # picking either would be arbitrary -- so that still refuses.
    hits = [mt for _r, mt in good if int(mt[0]) == int(good[0][1][0])]
    if len(hits) != len(good):
        raise ValueError(f"{len(good)} matches fit at different atoms; ambiguous")
    pos0 = pos[0]
    normals = []
    for h in hits:
        if len(h) < 3:
            continue
        a, b = pos0[h[1]] - pos0[h[0]], pos0[h[2]] - pos0[h[0]]
        n = np.cross(a, b)
        if np.linalg.norm(n) > 1e-6:
            normals.append(n / np.linalg.norm(n))
    for other in normals[1:]:
        if abs(float(normals[0] @ other)) < 0.98:              # ~11 degrees
            raise ValueError("matches at one reactive atom define "
                             "non-parallel planes")
    return hits[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=tc.topic())
    ap.add_argument("--cloud-dir", default=None,
                    help="where the persisted clouds are; defaults to "
                         "<topic>_allposes, but the 3.0.0 run misfiled its "
                         "clouds into nac_v3_allposes and this is how they are "
                         "read from where they actually landed")
    ap.add_argument("--min-mtime", default=None,
                    help="ISO timestamp; refuse clouds written before it, so a "
                         "file left by an earlier run cannot be mistaken for "
                         "this run's output")
    ap.add_argument("--out-dir", default=None,
                    help="write here instead of <topic>_poses; used to verify "
                         "the rebuild against a screen that wrote its own")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import datetime as _dt

    import nac_rank as nr
    import nac_screen_v2 as nv

    _, pose_dir, default_cloud = nv.topic_paths(args.topic)
    cloud_dir = Path(args.cloud_dir) if args.cloud_dir else default_cloud
    if args.out_dir:
        pose_dir = Path(args.out_dir)
    cutoff = (_dt.datetime.fromisoformat(args.min_mtime).timestamp()
              if args.min_mtime else None)
    log.info("clouds  %s", cloud_dir)
    log.info("writing %s", pose_dir)

    per_pose = per_pose_table(args.topic)
    tbl_modes = table_modes(args.topic)
    cands = {c.ident: c for c in nr.load_candidates()}

    n_ok = n_written = 0
    skip: dict = defaultdict(int)
    idents = sorted(per_pose)
    if args.limit:
        idents = idents[:args.limit]

    for i, ident in enumerate(idents, 1):
        if i % 250 == 0:
            log.info("  %d/%d  written=%d  skipped=%d",
                     i, len(idents), n_written, sum(skip.values()))
        # CHEAPEST REFUSAL FIRST. An already-written molecule is settled, and
        # deciding that costs one stat() -- so it happens before the ~430-pose
        # cloud is parsed, not after. Ordering it the other way made a resumed
        # pass re-read the entire library to write nothing.
        if not args.dry_run and (pose_dir / f"{ident}.sdf").exists():
            skip["already written"] += 1
            continue
        cand = cands.get(ident)
        if cand is None:
            skip["no candidate record"] += 1
            continue
        f = cloud_dir / f"{ident}.sdf"
        if not f.is_file():
            skip["no persisted cloud"] += 1
            continue
        # A CLOUD FROM AN EARLIER RUN IS NOT THIS RUN'S CLOUD. This is the exact
        # confusion being repaired, so it is refused rather than tolerated.
        if cutoff is not None and f.stat().st_mtime < cutoff:
            skip["cloud pre-dates this run"] += 1
            continue
        try:
            mol, eranks, modes = load_cloud(f)
        except Exception as e:                             # noqa: BLE001
            skip[f"unreadable cloud: {e}"] += 1
            continue

        from shared import nac_criterion as nac
        from shared import pose_modes as pmod

        g = per_pose[ident]
        n_pose = max(g) + 1
        labels = np.full(n_pose, -1, dtype=int)
        for pi, (md, _d, _a) in g.items():
            labels[pi] = md
        mode_ids = sorted(set(int(x) for x in labels) - {-1})
        if mode_ids != sorted(tbl_modes.get(ident, [])):
            skip["per-pose modes disagree with the aggregate"] += 1
            continue

        # UNDO THE WRITE PERMUTATION, THEN PROVE IT. If the reconstructed order
        # is right, the mode it implies for every cloud pose equals the `mode`
        # that pose actually carries. This is checked per pose, not per set: an
        # order that is wrong only for a few poses would still produce the right
        # set of modes and put those representatives in the wrong cluster.
        order = cloud_order(labels, mode_ids)
        if len(order) != len(modes) or not np.array_equal(labels[order], modes):
            skip["cloud order could not be reconstructed"] += 1
            continue

        dist = np.array([g[int(i)][1] for i in order], dtype=float)
        ang = np.array([g[int(i)][2] for i in order], dtype=float)
        try:
            match = resolve_match(mol, cand.reactive_smarts, dist)
        except Exception as e:                             # noqa: BLE001
            skip[f"reactive atom: {str(e)[:44]}"] += 1
            continue

        # The screen's own numbers, not a re-measurement -- so `anchor` here is
        # the `anchor` that chose the representative there.
        anchor = np.array([nac.anchor_quality(d, a, cand.mechanism)
                           for d, a in zip(dist, ang)])
        dmat = pmod.distances(pmod.features(mol, match))
        reps = nv.representative_indices(modes, anchor, dmat, mode_ids)

        n_ok += 1
        if args.dry_run:
            continue
        pose_dir.mkdir(parents=True, exist_ok=True)
        dest = pose_dir / f"{ident}.sdf"
        if dest.exists():                       # append_only: never overwrite
            skip["already written"] += 1
            continue
        nv.write_sdf(mol, reps, dest, modes=mode_ids)
        n_written += 1

    print(f"\n  {n_ok} molecules resolved, {n_written} written -> {pose_dir}")
    if skip:
        print("  skipped:")
        for k, v in sorted(skip.items(), key=lambda kv: -kv[1]):
            print(f"    {v:6d}  {k}")


if __name__ == "__main__":
    main()
