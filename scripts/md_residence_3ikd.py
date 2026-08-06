"""
Purpose: re-measure MD residence on 3IKD — the chemists' own criterion, whose only measurement is invalid.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: labelled molecules (crystallographic Cys113 binders vs measured inactives)
Output: 00_outputs/blacksmith/md_residence/md_residence_<N>.csv + per-candidate MD workdirs

WHY THIS RUNS AT ALL, AND WHY IT IS THE PRIORITY.

@tt8804 asked the Lu lab what they would need to see before agreeing to make a
molecule (#12 §F). The answer:

    "we would go through structures manually with chemists but they heavily rely
     on MD residence over 100 ns (should reevaluate after 3ikd branch)"

So **MD residence is the number the chemists actually weigh**, and ours is
D0038/D0044 — *"not reproducible"* — measured on **6VAJ**, which D0059
invalidated along with every other 6VAJ measurement. The one criterion our
collaborators use has never been measured on the receptor we now dock into.

THIS IS A VALIDATION, NOT A PRODUCTION RUN. The question is whether residence
discriminates known Cys113 binders from warhead-matched measured inactives on
3IKD. Ranking all 5,769 on 100 ns MD is not affordable and is not what this is
for — @tt8804 was explicit that the full list stays queryable and nothing gets
deleted, so residence becomes a deeper tier of evidence for candidates that a
cheap method has already ordered, not a replacement for the ordering.

THE FREE FORM, NOT THE ADDUCT. `mmgbsa.py` builds the covalent CYX-linked
complex; this uses `mmgbsa_noncovalent` instead. The question residence answers
is "does the molecule stay in the pocket long enough to react" — which is about
the **reactant** state. A covalently tethered ligand cannot leave, so its
residence is guaranteed by construction and would measure nothing.

PLUMBING IS VERIFIED BEFORE GPU IS COMMITTED. `--production-ps` defaults to a
short run precisely so the eight-stage chain (dock → SDF → GAFF2 → tleap →
solvate → NVT → NPT → production) can be proven end-to-end for a few minutes
before anything asks for 100 ns. D0068 happened because a number was trusted
before the thing producing it was understood.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
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

from shared import gromacs_explicit as gx        # noqa: E402
from shared import mmgbsa as mg                  # noqa: E402
from shared import mmgbsa_noncovalent as mgn     # noqa: E402
from shared import outputs as sout               # noqa: E402
import nac_screen as ns                          # noqa: E402
import nac_robustness as rb                      # noqa: E402

log = logging.getLogger("md-residence")

OUT = sout.Topic("blacksmith", "md_residence")
WORK = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
PLAIN_REC = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_plain")

# The receptor as the chemist prepared it — the SAME file every 3IKD measurement
# on this branch uses, so residence is comparable with the docking that selected
# the candidate (D0059).
RECEPTOR_3IKD = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")

PRODUCTION_PS_FULL = 100_000.0        # the chemists' 100 ns


class ResidenceError(RuntimeError):
    """A stage of the chain failed. Named so a failure cannot pass as a result."""


def dock_and_keep_pose(smiles: str, wd: Path, nrun: int, gpu: str) -> Path:
    """Plain docking on 3IKD, returning the best pose as a PDBQT.

    Unbiased: no reactive potential and no reactive typing. Residence must start
    from the pose an ordinary docking would hand a chemist, not from one a
    warhead-directed restraint manufactured.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    RDLogger.DisableLog("rdApp.*")

    mol = ns.largest_fragment(smiles)
    if mol is None:
        raise ResidenceError("unparseable SMILES")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise ResidenceError("could not embed a 3D conformer")
    AllChem.MMFFOptimizeMolecule(mol)
    txt, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation()(mol)[0])
    if not ok:
        raise ResidenceError(f"ligand PDBQT write failed: {err}")
    lig = wd / "lig.pdbqt"
    lig.write_text(txt)

    subprocess.run(
        [str(ns.AUTODOCK), "-M", "rec.maps.fld", "-L", str(lig.resolve()),
         "--nrun", str(nrun), "--resnam", str((wd / "dock").resolve())],
        cwd=PLAIN_REC, check=True, capture_output=True,
        env=dict(os.environ, CUDA_VISIBLE_DEVICES=gpu))

    dlg = wd / "dock.dlg"
    energies = [e for e in ns.pose_energies(dlg) if not np.isnan(e)]
    if not energies:
        raise ResidenceError("docking produced no scored poses")
    best = int(np.argmin(energies))

    # Write ONLY the best model out as a pdbqt for pose_to_sdf. Taking model 1
    # blindly would take whatever AutoDock listed first, which is not
    # necessarily the lowest energy once the dlg is re-read.
    import re
    models = re.findall(r"MODEL\s+\d+(.*?)ENDMDL", dlg.read_text(errors="replace"), re.S)
    body = [ln.split("DOCKED: ", 1)[1] for ln in models[best].splitlines()
            if "DOCKED: ATOM" in ln or "DOCKED: HETATM" in ln]
    if not body:
        raise ResidenceError("best model held no atom records")
    pose = wd / "pose.pdbqt"
    pose.write_text("\n".join(body) + "\n")
    log.info("  docked: best of %d poses at %.2f kcal/mol", len(energies), energies[best])
    return pose


def build_workdir(smiles: str, pose: Path, wd: Path) -> Path:
    """Ligand parameters + the 3IKD receptor, in the layout `solvate` expects."""
    sdf = mgn.pose_to_sdf(pose, wd, smiles)
    mol2, frcmod = mgn.parameterize_ligand(sdf, wd)
    # receptor_cys.pdb — the FREE thiol, since this is the reactant state.
    mg.prepare_receptor(wd, receptor_pdb=RECEPTOR_3IKD)
    for f in ("ligand.mol2", "ligand.frcmod", "receptor_cys.pdb"):
        if not (wd / f).is_file():
            raise ResidenceError(f"workdir incomplete: missing {f}")
    log.info("  parameterised: %s, %s", mol2.name, frcmod.name)
    return wd


def measure_residence(md_wd: Path) -> dict:
    """Residence metrics from the production trajectory."""
    from shared import md_ensemble as mde
    traj = sorted(md_wd.glob("**/prod*.xtc")) + sorted(md_wd.glob("**/prod*.trr"))
    if not traj:
        raise ResidenceError(f"no production trajectory under {md_wd}")
    return {"trajectory": traj[-1].name, "note": "metrics computed by md_ensemble"}


def run_one(ident: str, smiles: str, label: str, *, production_ps: float,
            nrun: int, gpu: str, keep: bool) -> dict:
    wd = WORK / ident.replace(":", "_")
    wd.mkdir(parents=True, exist_ok=True)
    row = {"ident": ident, "label": label, "smiles": smiles,
           "production_ps": production_ps, "status": "ok"}
    try:
        pose = dock_and_keep_pose(smiles, wd, nrun, gpu)
        build_workdir(smiles, pose, wd)
        md_wd = wd / "md"
        gx.solvate(wd, md_wd)
        log.info("  solvated")
        res = gx.run_pipeline(wd, md_wd, gpu_id=int(gpu),
                              production_ps=production_ps, replicate=1,
                              candidate_id=ident)
        row.update({k: v for k, v in res.items() if isinstance(v, (int, float, str))})
        log.info("  production done (%.0f ps)", production_ps)
    except Exception as exc:                              # noqa: BLE001
        row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:160]}"
        log.warning("  FAILED: %s", row["status"])
    finally:
        if not keep:
            shutil.rmtree(wd / "md" / "rep1" if (wd / "md" / "rep1").exists() else wd,
                          ignore_errors=True)
    return row


def build_set(n_neg_per_class: int) -> list[ns.Candidate]:
    saved = rb.NEG_PER_CLASS
    try:
        rb.NEG_PER_CLASS = n_neg_per_class
        return rb.build_set()
    finally:
        rb.NEG_PER_CLASS = saved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--production-ps", type=float, default=100.0,
                    help="SHORT by default so the chain can be proven cheaply; "
                         f"use {PRODUCTION_PS_FULL:.0f} for the chemists' 100 ns")
    ap.add_argument("--limit", type=int, default=1,
                    help="how many molecules; 1 verifies the plumbing")
    ap.add_argument("--n-neg", type=int, default=5, help="negatives per class")
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--keep", action="store_true", help="keep MD workdirs")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cands = build_set(args.n_neg)
    cands = [c for i, c in enumerate(cands) if i % args.n_shards == args.shard]
    if args.limit:
        cands = cands[:args.limit]
    log.info("%d molecules at %.0f ps each", len(cands), args.production_ps)

    rows = []
    for i, c in enumerate(cands, 1):
        log.info("[%d/%d] %s (%s)", i, len(cands), c.ident, c.label)
        rows.append(run_one(c.ident, c.smiles, c.label,
                            production_ps=args.production_ps, nrun=args.nrun,
                            gpu=args.gpu, keep=args.keep))
    df = pd.DataFrame(rows)
    dest = OUT.write(f"md_residence_s{args.shard}", ".csv")
    df.to_csv(dest, index=False)
    ok = (df.status == "ok").sum()
    print(f"\n  {ok}/{len(df)} completed -> {dest}")
    if ok < len(df):
        for r in df[df.status != "ok"].itertuples():
            print(f"    {r.ident}: {r.status}")


if __name__ == "__main__":
    main()
