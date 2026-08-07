"""
Purpose: automatic selection — walk the ranked list in order, re-measure the pose being elevated, and emit an elevation queue.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: rank_v2/rank_v2_<TIER>_<N>.csv + the persisted top-20 poses
Output: 00_outputs/blacksmith/elevation_queue/queue_<N>.csv

@tt8804, #22: "The default mode is automatic that just selects molecules for
elevation in order and double checking that the interaction distances and angles
are appropriate."

THE RE-CHECK MEASURES THE POSE, NOT THE SCORE. The ranking score is an aggregate
over a window of poses; what gets elevated is ONE specific pose. Those are
different objects and can disagree — a molecule can rank well on a window while
the single pose handed to MD is not in attack geometry. Re-reading the score
would launch a 4 GPU-hour trajectory on an assumption. So the pose is re-opened,
its reactive atom re-matched by SMARTS, and its distance and angle re-measured
against Cys113's SG in the receptor the MD will actually use.

A molecule whose elevated pose fails the geometry check is NOT silently dropped:
it is queued with `geometry_ok = False` and a reason, because "the ranking liked
it but its best pose is not reaction-competent" is a finding about the ranking,
not a filtering step to hide.

WHY THE QUEUE IS PER CLASS. The ranking is per warhead class (D0073), so taking
"the top N overall" would re-import the cross-class rigidity bias the per-class
quota exists to remove. Selection takes the top `--per-class` from each class.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import nac_criterion as nac           # noqa: E402
from shared import outputs as sout                # noqa: E402
from shared import receptors as R                 # noqa: E402

log = logging.getLogger("select")
OUT = sout.Topic("blacksmith", "elevation_queue")
DATA = Path("/data/lab_vm/append_only/inhibition")
RANK = DATA / "00_outputs/blacksmith/rank_v2"
POSES = DATA / "00_outputs/blacksmith/nac_v2_poses"


def sg_from_receptor() -> np.ndarray:
    """Cys113 SG, identified by its RECORDED COORDINATE, not by a residue number.

    The residue number is not stable across the pipeline: the prepared receptor
    keeps Pin1's own numbering (residues 51-163, so Cys113 is 113), while tleap
    renumbers to 1-115 for MD and Cys113 becomes 63. Hard-coding either one works
    on one file and silently finds nothing -- or worse, finds Cys57 -- on the
    other.

    `prepare_3ikd_receptor` recorded the SG coordinate at preparation time. That
    is the identity. This finds the CYS SG nearest it and refuses if the match is
    not essentially exact, so a receptor swap cannot pass unnoticed.
    """
    import json
    rec = R.resolve_3ikd_ian(noligand=True)
    meta = json.loads(sout.latest_path("blacksmith", "receptor_3ikd",
                                       "prep_3IKD", ".json").read_text())
    want = meta["verified"]["cys113_sg"]
    ref = np.array([want["x"], want["y"], want["z"]])

    best, best_d, best_resi = None, np.inf, None
    for line in rec.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != "CYS" or line[12:16].strip() != "SG":
            continue
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        d = float(np.linalg.norm(xyz - ref))
        if d < best_d:
            best, best_d, best_resi = xyz, d, line[22:26].strip()
    if best is None or best_d > 0.05:
        raise SystemExit(
            f"Cys113 SG not found at its recorded position in {rec.name}. "
            f"Nearest CYS SG is residue {best_resi} at {best_d:.3f} A from "
            f"{ref}. Either the receptor changed or the record is stale; "
            f"refusing to elevate against an unverified sulfur.")
    log.info("Cys113 SG = residue %s at %s (%.3f A from the recorded value)",
             best_resi, np.round(best, 3), best_d)
    return best


def recheck(ident: str, smarts: str, mechanism: str, sg: np.ndarray,
            policy: str = "best_viable") -> dict:
    """Re-measure the poses and choose WHICH ONE to elevate.

    The persisted SDF is ordered by docking energy, and the lowest-energy pose is
    frequently not the reaction-competent one -- measured here, two of the top
    three T_3 molecules had a rank-1 pose outside attack geometry, one of them
    8.0 A from the sulfur. That is the same thing the redock benchmark measured
    directly (AutoDock's own ordering puts the right pose first 18.3% of the
    time) and what D0068 found (more search finds lower-energy, less
    reaction-competent poses).

    So `best_viable` walks the poses in energy order and elevates the FIRST that
    passes the geometry check, recording which rank that was. Starting a 4
    GPU-hour trajectory from a pose 8 A off the sulfur would measure nothing
    about the reaction.

    `lowest_energy` is kept for comparison, not because it is defensible here.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = {"geometry_ok": False, "reason": "", "pose_rank": None,
           "n_poses_scanned": 0, "n_poses_viable": 0}
    sdf = POSES / f"{ident}.sdf"
    if not sdf.is_file():
        out["reason"] = "no persisted pose"
        return out
    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
    if not mols:
        out["reason"] = "no readable poses"
        return out
    patt = Chem.MolFromSmarts(smarts) if smarts else None
    if patt is None:
        out["reason"] = "no reactive SMARTS for this class"
        return out

    measured = []
    for k, mol in enumerate(mols, 1):
        match = mol.GetSubstructMatch(patt)
        if not match:
            continue
        conf = mol.GetConformer()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in match])
        try:
            # keywords, not positions: the signature is
            # measure(mechanism, coords, sg), and a silent transposition here
            # would measure geometry against the wrong atom set.
            measured.append((k, nac.measure(mechanism=mechanism, coords=coords, sg=sg)))
        except Exception:                          # noqa: BLE001
            continue
    if not measured:
        out["reason"] = "reactive SMARTS matched no persisted pose"
        return out

    out["n_poses_scanned"] = len(measured)
    out["n_poses_viable"] = sum(1 for _, r in measured if r.viable)

    if policy == "best_viable":
        pick = next(((k, r) for k, r in measured if r.viable), None)
        if pick is None:
            k, r = measured[0]
            out.update(_pose_fields(k, r, sdf))
            out["reason"] = (f"no viable pose among {len(measured)} persisted; "
                             f"best-energy pose is d={r.distance:.2f} A, "
                             f"{r.angle_kind}={r.angle:.1f} deg")
            return out
        k, r = pick
    else:
        k, r = measured[0]

    out.update(_pose_fields(k, r, sdf))
    if not r.viable:
        out["reason"] = (f"pose not in attack geometry (d={r.distance:.2f} A, "
                         f"{r.angle_kind}={r.angle:.1f} deg)")
    return out


def _pose_fields(rank: int, r, sdf: Path) -> dict:
    return {"pose_rank": rank, "pose_distance_A": r.distance,
            "pose_angle_deg": r.angle, "pose_angle_kind": r.angle_kind,
            "pose_approach_deg": r.approach_angle,
            "geometry_ok": bool(r.viable), "sdf": str(sdf)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--tier", default="both", choices=("T3", "T4", "both"),
                    help="'both' writes ONE queue. Writing a queue per tier and "
                         "letting downstream take the newest is how an empty T3 "
                         "queue stranded 17 good T4 molecules overnight.")
    ap.add_argument("--per-class", type=int, default=2)
    ap.add_argument("--score", default="weighted_score",
                    help="which ranking to select from, BY NAME. Taking the "
                         "newest file instead is how the overnight run selected "
                         "off enrichment_joint by accident.")
    ap.add_argument("--classes", nargs="*", default=None,
                    help="restrict to these warhead classes")
    ap.add_argument("--pose-policy", default="best_viable",
                    choices=("best_viable", "lowest_energy"),
                    help="which persisted pose gets elevated")
    ap.add_argument("--require-geometry", action="store_true",
                    help="queue only poses that pass the re-check")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    tiers = ["T3", "T4"] if args.tier == "both" else [args.tier]
    frames = []
    for t in tiers:
        fs = sorted(glob.glob(str(RANK / f"rank_v2_{t}_{args.score}_*.csv")),
                    key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            log.warning("no %s ranking for %s under %s", args.score, t, RANK)
            continue
        d = pd.read_csv(fs[-1])
        d["tier"] = t
        frames.append(d)
        log.info("ranking %s (%d rows)", Path(fs[-1]).name, len(d))
    if not frames:
        raise SystemExit(f"no ranking for {tiers} under {RANK}")
    df = pd.concat(frames, ignore_index=True)

    surv = df[df.passes & df.class_rank.notna()].copy()
    if args.classes:
        surv = surv[surv.warhead_class.isin(args.classes)]

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    smarts = dict(zip(wh.class_id, wh.reactive_atom_smarts))
    mech = dict(zip(wh.class_id, wh.mechanism))
    sg = sg_from_receptor()
    log.info("Cys113 SG at %s (from the receptor the MD will use)", np.round(sg, 3))

    rows = []
    for (tier, cls), g in surv.groupby(["tier", "warhead_class"]):
        for r in g.nsmallest(args.per_class, "class_rank").itertuples():
            rec = {"ident": r.ident, "warhead_class": cls,
                   "mechanism": mech.get(cls), "class_rank": int(r.class_rank),
                   "tier": args.tier,
                   "enrichment_conditional": getattr(r, "enrichment_conditional", np.nan),
                   "enrichment_joint": getattr(r, "enrichment_joint", np.nan),
                   "consensus_gnina": getattr(r, "consensus_gnina", np.nan),
                   "consensus_autodock": getattr(r, "consensus_autodock", np.nan),
                   "QED": getattr(r, "QED", np.nan),
                   "smiles": getattr(r, "smiles", None),
                   "selected_by": "automatic", "selected_reason":
                   f"rank {int(r.class_rank)} of {len(g)} in {cls}"}
            rec.update(recheck(r.ident, smarts.get(cls, ""), mech.get(cls, ""),
                               sg, args.pose_policy))
            rows.append(rec)

    q = pd.DataFrame(rows).sort_values(["warhead_class", "class_rank"])
    if args.require_geometry:
        dropped = (~q.geometry_ok).sum()
        q = q[q.geometry_ok]
        log.info("--require-geometry dropped %d", dropped)

    dest = OUT.write("queue", ".csv")
    q.to_csv(dest, index=False)

    print(f"\n=== elevation queue: {len(q)} molecules, top {args.per_class} per class ===\n")
    print(f"  {'ident':<20}{'class':<22}{'rk':>3}{'cond':>7}{'cons':>6}"
          f"{'d(A)':>7}{'angle':>7}{'pose':>5}  geom")
    for r in q.itertuples():
        d = getattr(r, "pose_distance_A", np.nan)
        a = getattr(r, "pose_angle_deg", np.nan)
        print(f"  {r.ident:<20}{r.warhead_class:<22}{r.class_rank:>3}"
              f"{r.enrichment_conditional:>7.2f}"
              f"{r.consensus_gnina if r.consensus_gnina == r.consensus_gnina else float('nan'):>6.2f}"
              f"{d:>7.2f}{a:>7.1f}{r.pose_rank if r.pose_rank==r.pose_rank else -1:>5}"
              f"  {'OK' if r.geometry_ok else 'FAIL'}")
    bad = q[~q.geometry_ok]
    if len(bad):
        print(f"\n  {len(bad)} queued with a FAILING geometry re-check — kept, not "
              f"hidden.\n  The ranking liked them; their elevated pose is not "
              f"reaction-competent.")
        for r in bad.itertuples():
            print(f"    {r.ident}: {r.reason}")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
