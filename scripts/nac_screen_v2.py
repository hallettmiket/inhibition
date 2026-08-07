"""
Purpose: re-run the geometric screen PERSISTING per-pose geometry and gnina scores, so the score can be repaired without docking again.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: all warhead-bearing candidates (T_3 + T_4), the chemist's 3IKD
Output: 00_outputs/blacksmith/nac_v2/{poses,agg}_s<shard>_<N>.csv + top-20 poses on disk

WHY THIS RE-RUN EXISTS. The 2.0.0 screen computed per-pose distances and angles,
used them to make one summary number per molecule, and threw the poses away. 79
molecules out of 5,765 have poses on disk. So two known repairs to the ranking
score are both blocked on the same missing data, and neither can be tested
without docking again:

  1. THE SCORE IS PARTLY MEASURING OUR OWN THUMB ON THE SCALE.
     `viable_fraction` is the JOINT rate P(distance in window AND angle in
     window), while `isotropic_null` is purely orientational, P(angle). The
     numerator therefore carries a distance hit rate -- and distance is exactly
     what the reactive potential is designed to bias (D0064: the reactive
     potential is a SAMPLER, not a criterion). Part of the score is how well the
     sampler worked on a molecule rather than the molecule's own ability to
     orient. The conditional form, P(angle | distance in window), cancels it:
     the bias applies equally to the numerator and to the conditioning set.

  2. WE ARE SCORING THE WRONG TEN POSES. Measured on this receptor tonight over
     82 Pin1 crystal ligands, AutoDock's own ordering puts a sub-2 A pose first
     18.3% of the time and gnina's CNN gets 26.8%, with sampling held fixed. Both
     consensus and the geometric score are computed on AutoDock's top ten, so a
     better ordering changes which ten they see.

WHAT IS PERSISTED, chosen so no third re-run is needed:

  poses_*.csv   one row per pose for the TOP 20 BY ENERGY: energy, distance,
                angle, approach, viable, plus gnina Affinity/CNNscore/
                CNNaffinity. Twenty is past every window worth testing (D0068
                argues for top-N; the redock benchmark found the right pose
                always inside the top 10 when present).
  agg_*.csv     counts over ALL poses -- n_poses, n_in_range, n_viable,
                n_viable_given_in_range -- so the conditional score can also be
                computed over the whole population, not only the retained window.
  poses/        the top-20 poses themselves, as SDF. The thing 2.0.0 discarded.

POSE ORDER IS THE JOIN KEY between geometry and gnina scores, so it is asserted
rather than assumed at every hop, exactly as `rescore_benchmark` does: counts
must match, and SDF coordinates are checked against the conformer they came from.
This project's defect catalogue is values taken by position instead of identity.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac           # noqa: E402
from shared import outputs as sout                # noqa: E402
from shared import covalent_protocol as cp        # noqa: E402
from shared import receptors as R                 # noqa: E402
import nac_screen as ns                           # noqa: E402
import nac_rank as nr                             # noqa: E402

log = logging.getLogger("nac-v2")
OUT = sout.Topic("blacksmith", "nac_v2")
POSE_DIR = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/nac_v2_poses")
KEEP_TOP = 20
_COORD_TOL = 0.05


def write_sdf(mol, order: list[int], dest: Path) -> int:
    """Write the top-`order` conformers to one SDF, stamping `pose_rank`.

    The rank is written as a PROPERTY, not left to file position, because
    `bpmd_run.read_pose` selects a pose by its own `pose_rank` and refuses to
    take one by index -- a pose identified by where it sits in a file is a pose
    that a re-sort silently redefines.
    """
    from rdkit import Chem
    w = Chem.SDWriter(str(dest))
    n = 0
    for rank, i in enumerate(order, 1):
        mol.SetProp("pose_rank", str(rank))
        mol.SetProp("energy_rank", str(rank))
        w.write(mol, confId=i)
        n += 1
    w.close()
    return n


def gnina_scores(receptor: Path, sdf: Path, gpu: str) -> pd.DataFrame:
    env = cp.gnina_env()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    r = subprocess.run(
        [str(cp.GNINA_BIN), "--receptor", str(receptor), "--ligand", str(sdf),
         "--score_only", "--cnn_scoring", "rescore", "--seed", "42"],
        capture_output=True, text=True, env=env, timeout=1800)
    if r.returncode != 0:
        raise RuntimeError(f"gnina rc={r.returncode}: {r.stderr[-200:]}")
    rows, cur = [], {}
    for line in r.stdout.splitlines():
        m = re.match(r"^(Affinity|CNNscore|CNNaffinity|CNNvariance):\s+(-?[\d.]+)", line)
        if not m:
            continue
        if m.group(1) == "Affinity" and cur:
            rows.append(cur)
            cur = {}
        cur[m.group(1)] = float(m.group(2))
    if cur:
        rows.append(cur)
    return pd.DataFrame(rows)


def one(cand, rec_dir: Path, plain_rec: Path, nrun: int, gpu: str,
        do_gnina: bool) -> tuple[pd.DataFrame, dict]:
    """Dock one candidate; return (per-pose rows, aggregate row)."""
    from rdkit import Chem
    work = Path(tempfile.mkdtemp(prefix="nacv2_"))
    agg = {"ident": cand.ident, "warhead_class": cand.warhead_class,
           "mechanism": cand.mechanism, "nrun": nrun}
    try:
        best = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res = ns.measure_dlg(dlg, cand)
            en = ns.pose_energies(dlg)
            if len(en) != len(res):
                raise ValueError(f"{len(en)} energies vs {len(res)} poses")
            frac = nac.viable_fraction(res)
            if best is None or frac > best[0]:
                best = (frac, res, en, dlg)
        if best is None:
            raise ValueError("no usable poses")
        frac, res, en, dlg = best

        in_rng = [nac.NAC_DIST_MIN <= r.distance <= nac.NAC_DIST_MAX for r in res]
        agg.update({
            "n_poses": len(res),
            "n_in_range": int(sum(in_rng)),
            "n_viable": int(sum(r.viable for r in res)),
            "n_viable_given_in_range": int(sum(r.viable and i for r, i in zip(res, in_rng))),
            "viable_fraction": frac,
            "enrichment": frac / nac.isotropic_null(cand.mechanism),
            "isotropic_null": nac.isotropic_null(cand.mechanism),
            "status": "ok",
        })

        order = sorted((i for i, e in enumerate(en) if not np.isnan(e)),
                       key=lambda i: en[i])[:KEEP_TOP]
        rows = pd.DataFrame([{
            "ident": cand.ident, "warhead_class": cand.warhead_class,
            "mechanism": cand.mechanism, "energy_rank": k + 1, "pose_idx": i,
            "energy": en[i], "distance": res[i].distance, "angle": res[i].angle,
            "approach_angle": res[i].approach_angle,   # None for SN2 by design
            "angle_kind": res[i].angle_kind,
            "viable": bool(res[i].viable),
            "in_range": bool(in_rng[i]),
        } for k, i in enumerate(order)])

        # persist the poses themselves -- the thing 2.0.0 threw away
        from meeko import PDBQTMolecule, RDKitMolCreate
        pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
        mol = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None][0]
        POSE_DIR.mkdir(parents=True, exist_ok=True)
        sdf = POSE_DIR / f"{cand.ident}.sdf"
        if not sdf.exists():                       # append_only: never overwrite
            write_sdf(mol, order, sdf)

        if do_gnina:
            tmp = work / "top.sdf"
            write_sdf(mol, order, tmp)
            g = gnina_scores(plain_rec, tmp, gpu)
            if len(g) != len(rows):
                raise ValueError(f"gnina returned {len(g)} scores for {len(rows)} poses")
            # verify the conformers gnina saw are the ones we measured
            supp = Chem.SDMolSupplier(str(tmp), removeHs=False)
            for k, (m, i) in enumerate(zip(supp, order)):
                if m is None:
                    raise ValueError(f"pose {k} unreadable in the SDF handed to gnina")
                a = mol.GetConformer(i).GetPositions()
                b = m.GetConformer().GetPositions()
                if a.shape != b.shape or float(np.abs(a - b).max()) > _COORD_TOL:
                    raise ValueError(f"pose {k}: SDF is not the conformer measured")
            for c in ("Affinity", "CNNscore", "CNNaffinity", "CNNvariance"):
                rows[c] = g[c].values if c in g else np.nan
        return rows, agg
    except Exception as exc:                       # noqa: BLE001
        agg["status"] = f"failed: {str(exc)[:160]}"
        return pd.DataFrame(), agg
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-gnina", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [s{args.shard}] %(message)s")
    os.nice(19)

    R.resolve_3ikd_ian()
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    plain_rec = sout.latest_path("blacksmith", "receptor_3ikd", "3IKD_prepared", ".pdbqt")
    log.info("reactive receptor ready; plain %s (3IKD_ian verified)", Path(plain_rec).name)

    cands = nr.load_candidates()
    cands = [c for i, c in enumerate(cands) if i % args.n_shards == args.shard]
    if args.limit:
        cands = cands[:args.limit]

    # Resume. A shard that dies at 90% must not restart from zero overnight.
    # Two rules, both learned the hard way:
    #   - only `ok` rows count as done. nac_rank.refine() counts `failed:` rows
    #     as done and therefore never retries a transient failure.
    #   - the run count is part of the identity. Skipping on ident alone would
    #     seat a 200-run measurement inside a re-run at a different effort,
    #     which is the same defect the BPMD resume had (catalogue #22).
    done = set()
    for f in sorted(glob.glob(str(OUT.dir / f"agg_s{args.shard}_*.csv"))):
        try:
            d = pd.read_csv(f)
        except Exception:                          # noqa: BLE001
            continue
        if {"ident", "status", "nrun"} <= set(d.columns):
            done |= set(d[(d.status == "ok") & (d.nrun == args.nrun)].ident)
    if done:
        before = len(cands)
        cands = [c for c in cands if c.ident not in done]
        log.info("resume: %d already complete at nrun=%d, %d remain of %d",
                 len(done), args.nrun, len(cands), before)
    log.info("%d candidates on this shard", len(cands))

    pose_buf, agg_buf, done = [], [], 0
    for i, c in enumerate(cands, 1):
        rows, agg = one(c, rec_dir, Path(plain_rec), args.nrun, args.gpu,
                        not args.no_gnina)
        agg_buf.append(agg)
        if not rows.empty:
            pose_buf.append(rows)
        done += 1
        if agg["status"] != "ok":
            log.warning("%s: %s", c.ident, agg["status"])
        if done % args.chunk == 0 or i == len(cands):
            if pose_buf:
                pd.concat(pose_buf, ignore_index=True).to_csv(
                    OUT.write(f"poses_s{args.shard}", ".csv"), index=False)
            pd.DataFrame(agg_buf).to_csv(
                OUT.write(f"agg_s{args.shard}", ".csv"), index=False)
            ok = sum(a["status"] == "ok" for a in agg_buf)
            log.info("[%d/%d] flushed; %d ok", i, len(cands), ok)
            pose_buf, agg_buf = [], []

    log.info("shard %d complete", args.shard)


if __name__ == "__main__":
    main()
