"""
Purpose: save the ACTUAL poses the near-attack score was computed from, so the GUI can show them.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: the shortlist_nac candidates from the latest T_3/T_4 frames
Output: append_only/.../nac_poses/<candidate_id>.sdf + nac_pose_path in new frames

WHY THIS EXISTS, AND WHY THE POSES ALREADY ON DISK WOULD BE THE WRONG ANSWER.

`scripts/nac_rank.py` docked every candidate into the **reactive 3IKD** receptor
and threw the poses away with the temporary directory. The frames still carry
`pose_path` from the earlier production docking — a different receptor and a
different protocol — and `pose3d.find_pose` resolves those happily.

So a pose viewer wired naively into the near-attack panel would show a
confident, correctly-rendered pose that **is not the pose the number describes**.
That is this project's signature defect exactly: a value taken by the nearest
plausible source rather than by identity, and `find_pose`'s own docstring already
records a version of it ("overlap is exactly zero across 4080 files").

This re-docks the shortlist under the SAME protocol that produced the score and
writes the pose beside it, so the two cannot drift apart.

WHICH POSE IS SAVED. The best NAC-viable pose — the lowest-energy pose that
CLEARS the geometric criterion — not the lowest-energy pose overall. The panel is
showing "can this molecule present its warhead", and the pose that answers that
question is the one that does so, not the one that binds best while pointing the
wrong way (D0062: 52.7% of poses get the reactive region right while failing a
whole-molecule RMSD test).

When no pose clears the criterion the molecule is written with
`nac_pose_viable = False` and its best pose anyway, rather than being skipped.
"Nothing to show" and "nothing worked" are different facts, and a gap in the
viewer would read as the former.
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import nac_criterion as nac          # noqa: E402
from shared import pose_consensus as pc          # noqa: E402
import nac_screen as ns                          # noqa: E402
import nac_rank as nr                            # noqa: E402

log = logging.getLogger("nac-poses")

DATA = Path("/data/lab_vm/append_only/inhibition")
POSES = DATA / "00_outputs" / "blacksmith" / "nac_poses"
FRAMES = {"T_3": ("03_t3_reinvent", "D3"), "T_4": ("04_t4_combinatorial", "D4")}


def latest_frame(subdir: str, stem: str) -> Path:
    fs = list((DATA / subdir).glob(f"{stem}_*.parquet"))
    return max(fs, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1)))


def next_version(subdir: str, stem: str) -> int:
    fs = list((DATA / subdir).glob(f"{stem}_*.parquet"))
    return max(int(re.search(r"_(\d+)\.parquet$", p.name).group(1)) for p in fs) + 1


def export_one(cand: ns.Candidate, rec_dir: Path, nrun: int, gpu: str,
               top_n: int = 10) -> dict:
    """Re-dock under the scoring protocol and write the TOP-N poses.

    ALL TOP-N, NOT JUST THE BEST. Two reasons, and the first is the panel's
    whole point: consensus is a statement about whether the top poses AGREE
    (`shared/pose_consensus.py`), and a single pose makes that unreadable — a
    molecule whose top ten land in one place and one whose top ten scatter across
    the pocket look identical when you only draw one of them.

    The second is that the viewer already supports every mode and selects one by
    default; showing a single pose was throwing that away.

    ORDERED WITH THE BEST NAC-VIABLE POSE FIRST, so the mode shown by default is
    the one the enrichment is about, not merely the lowest-energy one. The rest
    follow by energy, each carrying whether it cleared the criterion.
    """
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    work = Path(tempfile.mkdtemp(prefix="nacpose_"))
    try:
        best = None
        for j, lig in enumerate(ns.prepare_ligand(cand, work / "lig.pdbqt")):
            dlg = ns.dock(lig, rec_dir, work / f"c{j}", nrun, gpu)
            res = ns.measure_dlg(dlg, cand)
            energies = ns.pose_energies(dlg)
            if len(energies) != len(res):
                raise ValueError("energy/geometry length mismatch")
            pm = PDBQTMolecule.from_file(str(dlg), is_dlg=True, skip_typing=True)
            mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
            if not mols:
                continue
            mol = mols[0]
            # Prefer the lowest-energy pose that CLEARS the criterion.
            viable = [(e, i) for i, (r, e) in enumerate(zip(res, energies))
                      if r.viable and not np.isnan(e)]
            match = mol.GetSubstructMatches(Chem.MolFromSmarts(cand.reactive_smarts))
            if not match:
                continue
            rx_idx = list(match[0])
            if viable:
                e, i = min(viable)
                cand_best = (True, e, i, mol, res[i], res, energies)
            else:
                fin = [(e, i) for i, e in enumerate(energies) if not np.isnan(e)]
                if not fin:
                    continue
                e, i = min(fin)
                cand_best = (False, e, i, mol, res[i], res, energies)
            if best is None or (cand_best[0], -cand_best[1]) > (best[0], -best[1]):
                best = cand_best

        if best is None:
            raise ValueError("no usable pose")
        viable, energy, idx, mol, geom, res, energies = best

        # Rank: best NAC-viable first, then the rest by energy. `order` indexes
        # into the docked conformers.
        finite = [i for i, e in enumerate(energies) if not np.isnan(e)]
        viable_i = sorted((i for i in finite if res[i].viable),
                          key=lambda i: energies[i])
        other_i = sorted((i for i in finite if not res[i].viable),
                         key=lambda i: energies[i])
        order = (viable_i + other_i)[:top_n]

        POSES.mkdir(parents=True, exist_ok=True)
        out = POSES / f"{cand.ident}.sdf"
        w = Chem.SDWriter(str(out))
        for rank, i in enumerate(order, 1):
            one = Chem.Mol(mol)
            one.RemoveAllConformers()
            one.AddConformer(mol.GetConformer(i), assignId=True)
            one.SetProp("_Name", f"{cand.ident}_pose{rank}")
            one.SetProp("pose_rank", str(rank))
            one.SetProp("nac_viable", str(res[i].viable))
            one.SetProp("nac_energy_kcal", f"{energies[i]:.3f}")
            one.SetProp("nac_distance_A", f"{res[i].distance:.3f}")
            one.SetProp("nac_angle_deg", f"{res[i].angle:.2f}")
            one.SetProp("nac_angle_kind", res[i].angle_kind)
            one.SetProp("receptor", "3IKD reactive (D0059/D0064)")
            w.write(one)
        w.close()

        # Consensus over the SAME poses, so what the viewer draws and what the
        # number describes cannot drift apart.
        cons = None
        try:
            rp = [pc.ReactivePose(energy=float(energies[i]),
                                  reactive_xyz=np.asarray(
                                      mol.GetConformer(i).GetPositions()[list(rx_idx)]),
                                  atom_ids=tuple(rx_idx))
                  for i in finite]
            if len(rp) >= pc.MIN_POSES_FOR_CONSENSUS:
                cons = pc.consensus(rp, top_n=min(top_n, len(rp)))
        except Exception as exc:                        # noqa: BLE001
            log.debug("%s: consensus unavailable: %s", cand.ident, exc)

        return {"ident": cand.ident, "nac_pose_path": str(out),
                "nac_pose_viable": viable, "nac_pose_energy": energy,
                "nac_pose_distance": geom.distance, "nac_pose_angle": geom.angle,
                "nac_n_poses_shown": len(order),
                "nac_n_viable_shown": sum(1 for i in order if res[i].viable),
                "nac_consensus": cons.agreement if cons else np.nan,
                "nac_consensus_n": cons.n_poses if cons else np.nan,
                "nac_consensus_median_rmsd": cons.median_rmsd if cons else np.nan,
                "nac_consensus_modes": cons.n_modes if cons else np.nan,
                "status": "ok"}
    except Exception as exc:                            # noqa: BLE001
        return {"ident": cand.ident, "nac_pose_path": None,
                "nac_pose_viable": None, "status": f"failed: {str(exc)[:120]}"}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    by_id = {c.ident: c for c in nr.load_candidates()}
    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)

    for approach, (subdir, stem) in FRAMES.items():
        src = latest_frame(subdir, stem)
        df = pd.read_parquet(src)
        if "shortlist_nac" not in df.columns:
            log.warning("%s: no shortlist_nac in %s", approach, src.name)
            continue
        short = df[df.shortlist_nac].nlargest(args.top, "nac_enrichment")
        log.info("%s: exporting %d poses", approach, len(short))

        rows = []
        for i, cid in enumerate(short.candidate_id.astype(str), 1):
            c = by_id.get(cid)
            if c is None:
                rows.append({"ident": cid, "nac_pose_path": None,
                             "status": "not in candidate set"})
                continue
            r = export_one(c, rec, args.nrun, args.gpu)
            rows.append(r)
            log.info("  [%d/%d] %s %s", i, len(short), cid,
                     "viable" if r.get("nac_pose_viable") else r["status"])

        poses = pd.DataFrame(rows)
        merged = df.merge(poses.drop(columns=["status"]), how="left",
                          left_on=df.candidate_id.astype(str), right_on="ident")
        merged = merged.drop(columns=[c for c in ("key_0", "ident") if c in merged])
        dest = DATA / subdir / f"{stem}_{next_version(subdir, stem)}.parquet"
        merged.to_parquet(dest, index=False)
        ok = int(poses.nac_pose_path.notna().sum())
        log.info("  %d/%d poses written -> %s", ok, len(poses), dest.name)


if __name__ == "__main__":
    main()
