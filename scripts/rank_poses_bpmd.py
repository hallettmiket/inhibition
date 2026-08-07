"""
Purpose: rank the POSES within one molecule by BPMD confidence, so 100 ns MD starts from the pose that earned it.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/elevation_queue/queue_<N>.csv + the persisted v2 poses
Output: 00_outputs/blacksmith/pose_rank_bpmd/pose_rank_<N>.csv

@tt8804: "selecting a molecule doesnt mean we just send its lowest energy pose,
we need to rerun the poses through bpmd first to get a confidence score per pose
and from there we run MD on the most confident best angled pose. We need to
sperate that we are ranking molecules using pose consensus filtering and
enrichment/anchor scoring and then within each molecule during elevation we rank
the poses."

TWO DIFFERENT RANKING PROBLEMS, AND THEY WERE BEING CONFLATED.

  BETWEEN molecules (stage 2)  consensus filter + the weighted anchoring score.
                               Answers: which molecules are worth the compute.
  WITHIN a molecule (stage 4)  THIS. Answers: of this molecule's poses, which
                               one is real enough to spend 4 GPU-hours on.

The first cut of selection picked the elevated pose on geometry alone -- the
first pose in energy order that cleared the near-attack window. That is a
statement about whether a pose COULD react, not about whether it is physically
there. Docking energy is a poor witness (AutoDock's own ordering puts the
right pose first 18.3% of the time on this receptor), and geometry is a
threshold, not evidence of stability. BPMD is the per-pose measurement that
arbitrates, and it is cheap next to what it protects: ~1 GPU-hour to choose the
pose, against 4 to run the wrong one.

THE READOUT IS OCCUPANCY, NOT ESCAPE COST. The completed elevation run measured
`bias_at_exit` separating nothing (all p >= 0.08) while tracking occupancy at
rho = 0.974 -- the escape-cost term is nearly inert at 3 ns. `frac_in_window` is
what discriminated crystallographic positives from candidates (p = 0.005-0.021),
so it is what ranks poses here.

"MOST CONFIDENT, BEST ANGLED" IS TWO CRITERIA AND THEY ARE KEPT SEPARATE. A pose
that holds the window under bias but sits at a poor approach angle is not the
same as one that holds it in near-ideal geometry. Both are reported, and the
winner is chosen on occupancy among poses that clear the geometry criterion --
so stability decides among reaction-competent poses rather than overriding them.
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
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac           # noqa: E402
from shared import outputs as sout                # noqa: E402
import nac_screen as ns                           # noqa: E402
import bpmd_run as br                             # noqa: E402

log = logging.getLogger("pose-rank")
OUT = sout.Topic("blacksmith", "pose_rank_bpmd")
DATA = Path("/data/lab_vm/append_only/inhibition")
QUEUE = DATA / "00_outputs/blacksmith/elevation_queue"
POSES = DATA / "00_outputs/blacksmith/nac_v2_poses"


def viable_ranks(ident: str, smarts: str, mechanism: str, sg: np.ndarray,
                 max_poses: int) -> list[dict]:
    """The pose ranks worth spending BPMD on: reaction-competent ones, best first."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sdf = POSES / f"{ident}.sdf"
    if not sdf.is_file():
        return []
    patt = Chem.MolFromSmarts(smarts) if smarts else None
    if patt is None:
        return []
    out = []
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False):
        if m is None or not m.HasProp("pose_rank"):
            continue
        match = m.GetSubstructMatch(patt)
        if not match:
            continue
        conf = m.GetConformer()
        coords = np.array([list(conf.GetAtomPosition(i)) for i in match])
        try:
            r = nac.measure(mechanism=mechanism, coords=coords, sg=sg)
        except Exception:                          # noqa: BLE001
            continue
        if r.viable:
            out.append({"pose_rank": int(m.GetProp("pose_rank")),
                        "distance_A": r.distance, "angle_deg": r.angle,
                        "angle_kind": r.angle_kind,
                        "approach_deg": r.approach_angle})
    return out[:max_poses]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--max-poses", type=int, default=3,
                    help="candidate poses to test per molecule")
    ap.add_argument("--replicates", type=int, default=2)
    ap.add_argument("--production-ps", type=float, default=3000.0)
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--molecules", nargs="*", default=None)
    ap.add_argument("--queue", default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    q = Path(args.queue) if args.queue else None
    if q is None:
        fs = sorted(glob.glob(str(QUEUE / "queue_*.csv")),
                    key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        if not fs:
            raise SystemExit(f"no queue under {QUEUE}")
        q = Path(fs[-1])
    qdf = pd.read_csv(q)
    if args.molecules:
        qdf = qdf[qdf.ident.isin(args.molecules)]
    qdf = qdf[qdf.geometry_ok] if "geometry_ok" in qdf.columns else qdf
    log.info("queue %s: %d molecules to rank poses for", q.name, len(qdf))

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    smarts = dict(zip(wh.class_id, wh.reactive_atom_smarts))
    mech = dict(zip(wh.class_id, wh.mechanism))
    sg = br.sg_in_pose_frame(br.CYS_PDB, br.CYX_INDEX) if hasattr(br, "CYS_PDB") \
        else None
    if sg is None:
        import select_elevate as se
        sg = se.sg_from_receptor()

    rows = []
    for r in qdf.itertuples():
        cls = r.warhead_class
        cands = viable_ranks(r.ident, smarts.get(cls, ""), mech.get(cls, ""),
                             sg, args.max_poses)
        if not cands:
            log.warning("%s: no reaction-competent pose to test", r.ident)
            continue
        log.info("%s: testing %d poses (ranks %s)", r.ident, len(cands),
                 [c["pose_rank"] for c in cands])
        if args.dry_run:
            for c in cands:
                rows.append({"ident": r.ident, "warhead_class": cls, **c,
                             "status": "dry-run"})
            continue

        cand = ns.Candidate(ident=r.ident, smiles=r.smiles, warhead_class=cls,
                            mechanism=mech.get(cls, ""),
                            reactive_smarts=smarts.get(cls, ""), label="candidate")
        for c in cands:
            got = []
            try:
                got = br.run_pose(
                    cand, replicates=args.replicates,
                    production_ps=args.production_ps, gpu=args.gpu,
                    threads=args.threads, nrun=200, dock_gpu=str(args.gpu),
                    allow_redock=False, on_row=lambda row: None,
                    reuse_equilibration=True,
                    pose_rank=c["pose_rank"], sdf=POSES / f"{r.ident}.sdf")
            except Exception as exc:               # noqa: BLE001
                log.warning("%s pose %d: %s", r.ident, c["pose_rank"],
                            str(exc)[:140])
            ok = [g for g in got if g.get("status") == "ok"]
            rec = {"ident": r.ident, "warhead_class": cls, **c,
                   "n_replicates": len(ok),
                   "status": "ok" if ok else "failed"}
            if ok:
                rec["frac_in_window"] = float(np.mean([g["frac_in_window"] for g in ok]))
                rec["frac_in_window_sd"] = float(np.std([g["frac_in_window"] for g in ok]))
                rec["max_cv_nm"] = float(np.mean([g["max_cv_nm"] for g in ok]))
            rows.append(rec)
            log.info("  pose %d: %s", c["pose_rank"],
                     f"occupancy {rec.get('frac_in_window', float('nan')):.3f}"
                     if ok else "failed")

    if not rows:
        raise SystemExit("nothing measured")
    df = pd.DataFrame(rows)
    dest = OUT.write("pose_rank", ".csv")

    # the winner per molecule: highest occupancy among reaction-competent poses
    if "frac_in_window" in df.columns:
        df["is_winner"] = False
        for ident, g in df[df.status == "ok"].groupby("ident"):
            df.loc[g.frac_in_window.idxmax(), "is_winner"] = True
    df.to_csv(dest, index=False)

    print(f"\n=== pose ranking within each molecule ===\n")
    print(f"  {'ident':<22}{'pose':>5}{'d(A)':>7}{'angle':>7}"
          f"{'occupancy':>11}{'reps':>6}  winner")
    for r in df.itertuples():
        occ = getattr(r, "frac_in_window", float("nan"))
        print(f"  {r.ident:<22}{r.pose_rank:>5}{r.distance_A:>7.2f}"
              f"{r.angle_deg:>7.1f}{occ:>11.3f}{getattr(r, 'n_replicates', 0):>6}"
              f"  {'<--' if getattr(r, 'is_winner', False) else ''}")
    print(f"\n  -> {dest}")
    print("\n  The winner is the most STABLE pose among the reaction-competent "
          "ones —\n  stability decides among viable poses rather than overriding "
          "the geometry.\n  Occupancy is the readout because bias-at-exit "
          "separated nothing at 3 ns.")


if __name__ == "__main__":
    main()
