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

ONE NAMED CANDIDATE, FROM A POSE THAT ALREADY EXISTS (`--candidate`/`--pose`).
The validation set above is the reason this script was written, but the chemists
also ask for a deep workup on a SINGLE molecule, and re-docking it would answer
a different question than the one the ranking asked. `--pose` therefore starts
the chain from a saved pose instead of `dock_and_keep_pose`.

    THE TWO ENTRY POINTS DO NOT SAMPLE THE SAME THING, and the output records
    which one ran (`pose_source`). `dock_and_keep_pose` is deliberately
    unbiased — plain AutoDock-GPU, no reactive potential. The poses under
    `export_nac_poses.poses_dir()` come from the REACTIVE receptor, where a
    warhead-directed potential steered the search (D0063/D0064). Residence
    measured from a reactive pose is conditioned on the warhead already being
    presented; residence measured from a plain pose is not. Neither is wrong,
    they are different questions, and reporting one as the other would be the
    project's signature defect in yet another dress.
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


def saved_pose_to_sdf(pose_sdf: Path, wd: Path, smiles: str,
                      pose_rank: int = 1) -> Path:
    """A pose from an existing multi-model SDF, in `pose_to_sdf`'s output form.

    Deliberately mirrors `mmgbsa_noncovalent.pose_to_sdf` step for step —
    hydrogens stripped by hand, bond orders imposed from the candidate's own
    SMILES, hydrogens re-added at that protonation state — so the two entry
    points hand antechamber the same KIND of file and a difference downstream
    cannot be an artefact of how the pose was read.

    `pose_rank` is matched against the SDF's own `pose_rank` property rather
    than taken by position. The export writes them in rank order today, which
    is exactly the assumption that makes taking element [0] look correct
    forever and fail silently the day the order changes.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mols = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m]
    if not mols:
        raise ResidenceError(f"no readable pose in {pose_sdf}")

    ranked = [m for m in mols
              if m.HasProp("pose_rank") and int(m.GetProp("pose_rank")) == pose_rank]
    if not ranked:
        raise ResidenceError(
            f"{pose_sdf.name} carries no pose with pose_rank={pose_rank} "
            f"(has {[m.GetProp('pose_rank') for m in mols if m.HasProp('pose_rank')]})")
    pose = ranked[0]

    # STRIP HYDROGENS, THEN GIVE THE HEAVY ATOMS THEIR IMPLICIT VALENCE BACK.
    #
    # `mmgbsa_noncovalent.pose_to_sdf` removes hydrogens the same way, but its
    # input is an UNSANITIZED PDB where the heavy atoms never had `noImplicit`
    # set. An SDF read through SDMolSupplier is sanitized with explicit
    # hydrogens, so every heavy atom carries `noImplicit=True`. Deleting the
    # hydrogens then leaves carbons that RDKit believes are complete at valence
    # 2, `Chem.AddHs` finds nothing to add, and the file that reaches
    # antechamber is heavy-atoms-only.
    #
    # THAT IS NOT A CRASH-ON-CONTACT BUG. RDKit sanitizes it, writes it, and
    # reports the correct SMILES for it, because a heavy-atom skeleton with
    # consistent bond orders is a valid molecule — just not this one. It got as
    # far as antechamber's acdoctor ("Weird atomic valence (2) for atom C1")
    # only because acdoctor checks valences that RDKit had been told not to.
    #
    # COUNT THE HYDROGENS BEFORE DELETING THEM, AND RESTATE THE COUNT.
    #
    # The first version of this fix cleared `noImplicit` and zeroed
    # `numExplicitHs`, letting RDKit re-derive the hydrogen count from valence.
    # That works for carbon and it is wrong for an aromatic nitrogen, where the
    # hydrogen is not implied by the valence — it is the whole difference
    # between a pyrrole-type N and a pyridine-type one. Imidazole has one of
    # each; drop both hydrogens and the ring can no longer be kekulized, so
    # sanitisation throws and the molecule never reaches MD at all.
    #
    # Measured over all 5,790 persisted poses: re-deriving reproduces the pose's
    # own formula for 94.1% of them and LOSES 339 (5.9%), including the
    # reference BJP-06-005-3. Restating the observed count reproduces 100.0% and
    # breaks none of the ones that already worked.
    n_hydrogens = {a.GetIdx(): sum(1 for nb in a.GetNeighbors()
                                   if nb.GetAtomicNum() == 1)
                   for a in pose.GetAtoms() if a.GetAtomicNum() > 1}
    rw = Chem.RWMol(pose)
    for idx in sorted([a.GetIdx() for a in pose.GetAtoms()
                       if a.GetAtomicNum() == 1], reverse=True):
        rw.RemoveAtom(idx)
    for a in rw.GetAtoms():
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(n_hydrogens[a.GetIdx()])
    Chem.SanitizeMol(rw)
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ResidenceError(f"unparseable candidate SMILES {smiles!r}")
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(template, rw.GetMol())
    except Exception as exc:                              # noqa: BLE001
        raise ResidenceError(
            f"saved pose does not match the frame's SMILES ({exc})") from exc

    fixed = Chem.AddHs(fixed, addCoords=True)

    # The molecule written out must be the molecule the frame names, hydrogens
    # included. Comparing FORMULAE rather than heavy-atom counts is the point:
    # the bug above passed a heavy-atom check with room to spare and was only
    # visible in the hydrogens.
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    want = CalcMolFormula(Chem.AddHs(template))
    got = CalcMolFormula(fixed)
    if want != got:
        raise ResidenceError(
            f"pose rebuilt as {got} but the frame's SMILES is {want}; the "
            "structure handed to antechamber is not this candidate")

    out = wd / "ligand_pose.sdf"
    w = Chem.SDWriter(str(out))
    w.write(fixed)
    w.close()
    log.info("  pose: rank %d of %d from %s", pose_rank, len(mols), pose_sdf.name)
    return out


def build_workdir(smiles: str, pose: Path, wd: Path, *,
                  net_charge: int = 0, pose_rank: int = 1) -> Path:
    """Ligand parameters + the 3IKD receptor, in the layout `solvate` expects."""
    if pose.suffix.lower() == ".sdf":
        sdf = saved_pose_to_sdf(pose, wd, smiles, pose_rank)
    else:
        sdf = mgn.pose_to_sdf(pose, wd, smiles)

    # THE CHARGE ANTECHAMBER IS TOLD MUST BE THE CHARGE THE STRUCTURE HAS.
    # `-nc` sets the total for AM1-BCC; hand it +1 for a structure whose formal
    # charges sum to 0 and sqm solves an open-shell radical cation, then
    # distributes an entire spurious electron across the molecule — populated,
    # plausible, wrong. The frame's `charge_ph74` describes the species at pH
    # 7.4; the docked pose is whatever the docking built. When they disagree
    # that is a real finding about the pipeline, not something to paper over,
    # so it is raised rather than silently coerced.
    from rdkit import Chem as _Chem
    built = [m for m in _Chem.SDMolSupplier(str(sdf), removeHs=False) if m][0]
    structural_q = _Chem.GetFormalCharge(built)
    if structural_q != net_charge:
        raise ResidenceError(
            f"net_charge={net_charge:+d} was requested but the pose's formal "
            f"charges sum to {structural_q:+d}. Protonate the pose to match, or "
            f"pass --net-charge {structural_q:+d} and record that the simulated "
            "species is not the pH-7.4 species.")
    mol2, frcmod = mgn.parameterize_ligand(sdf, wd, net_charge=net_charge)
    # receptor_cys.pdb — the FREE thiol, since this is the reactant state.
    mg.prepare_receptor(wd, receptor_pdb=RECEPTOR_3IKD)
    for f in ("ligand.mol2", "ligand.frcmod", "receptor_cys.pdb"):
        if not (wd / f).is_file():
            raise ResidenceError(f"workdir incomplete: missing {f}")
    log.info("  parameterised: %s, %s", mol2.name, frcmod.name)
    return wd


def measure_residence(rep_wd: Path) -> dict:
    """Residence metrics from the production trajectory.

    THIS USED TO RETURN A PLACEHOLDER. The previous body located the trajectory
    and returned `{"trajectory": ..., "note": "metrics computed by md_ensemble"}`
    — a dict that reads like a result, contains no measurement, and was never
    called by `run_one` anyway. So the "verified end to end" chain produced a
    completed 100 ns trajectory and NO residence number, and the CSV looked
    populated because the pipeline metadata (atoms, waters, ns/day) filled it.
    A missing metric that renders as a full row is exactly the disguise this
    project keeps re-finding.

    `md_ensemble.residence_metrics` cannot be used here and that is not a
    detail: it takes an (n_frames, n_atoms, 3) array from an OpenMM IMPLICIT
    trajectory. This tier is explicit-solvent GROMACS, whose XTC has no such
    reader in this environment (see gromacs_analysis's own docstring). The two
    modules measure ligand RMSD the same way — nanometres of ligand
    displacement after superposing on the protein — so `explicit_ligand_rmsd_*`
    IS comparable with `ligand_rmsd_*`. The CONTACT COUNTS ARE NOT: `gmx
    mindist -on` and the implicit tier's heavy-atom pair count are different
    definitions and differ several-fold on the same complex. The keys stay
    distinctly named so the two cannot be silently pooled.
    """
    from shared import gromacs_analysis as ga
    traj = rep_wd / "prod.xtc"
    if not traj.is_file() or traj.stat().st_size == 0:
        raise ResidenceError(f"no production trajectory at {traj}")
    out = ga.analyse(rep_wd)
    out.update(residence_uncertainty(rep_wd))
    return out


def residence_uncertainty(rep_wd: Path) -> dict:
    """A correlation-corrected error bar for a SINGLE trajectory.

    `gromacs_analysis.analyse` returns means with no uncertainty, and a mean
    over 10,000 correlated frames quoted bare invites the reader to treat it as
    if it had 10,000 independent samples. Consecutive MD frames are nothing of
    the kind.

    `md_ensemble.statistical_inefficiency` (g = 1 + 2*sum of the
    autocorrelation) is reused rather than reimplemented, so this number and
    the implicit-solvent tier's are corrected the same way. The SEM is widened
    by sqrt(g), and `n_independent` is reported so a reader can see how few
    samples 100 ns actually contains.

    THIS IS A WITHIN-RUN ERROR BAR AND NOT A REPLICATE SPREAD. It describes the
    precision of the mean of ONE trajectory. It says nothing about whether a
    second trajectory from different velocities would land anywhere near it —
    which, on this project's own evidence, is the question that matters
    (D0044: explicit-solvent residence "is not reproducible either"). Quoting
    it as though it bounded run-to-run variation would understate the real
    uncertainty by whatever that variation happens to be.
    """
    from shared import md_ensemble as mde
    rmsd = ga_read(rep_wd / "rmsd.xvg")
    cont = ga_read(rep_wd / "numcont.xvg")
    out = {}
    for name, arr in (("explicit_ligand_rmsd_nm", rmsd), ("gmx_contacts", cont)):
        y = arr[:, 1]
        g = mde.statistical_inefficiency(y)
        n_ind = max(1.0, y.size / g)
        out[f"{name}_sd"] = round(float(y.std(ddof=1)), 4)
        out[f"{name}_sem"] = round(float(y.std(ddof=1) / np.sqrt(n_ind)), 4)
        out[f"{name}_stat_inefficiency"] = round(float(g), 2)
        out[f"{name}_n_independent"] = round(float(n_ind), 1)
    out["uncertainty_scope"] = ("within one trajectory, autocorrelation-"
                                "corrected; NOT a replicate spread (D0044)")
    return out


def ga_read(path: Path) -> np.ndarray:
    """The XVG `gromacs_analysis` already wrote, re-read rather than recomputed."""
    from shared import gromacs_analysis as ga
    return ga._read_xvg(path)


def run_one(ident: str, smiles: str, label: str, *, production_ps: float,
            nrun: int, gpu: str, keep: bool, pose: Path | None = None,
            pose_rank: int = 1, net_charge: int = 0,
            work_root: Path | None = None) -> dict:
    wd = (work_root or WORK) / ident.replace(":", "_")
    wd.mkdir(parents=True, exist_ok=True)
    row = {"ident": ident, "label": label, "smiles": smiles,
           "production_ps": production_ps, "net_charge": net_charge,
           "pose_source": "saved" if pose else "plain redock",
           "pose_path": str(pose) if pose else "",
           "status": "ok"}
    try:
        p = pose if pose else dock_and_keep_pose(smiles, wd, nrun, gpu)
        build_workdir(smiles, p, wd, net_charge=net_charge, pose_rank=pose_rank)
        md_wd = wd / "md"
        gx.solvate(wd, md_wd)
        log.info("  solvated")
        res = gx.run_pipeline(wd, md_wd, gpu_id=int(gpu),
                              production_ps=production_ps, replicate=1,
                              candidate_id=ident)
        row.update({k: v for k, v in res.items() if isinstance(v, (int, float, str))})
        log.info("  production done (%.0f ps)", production_ps)
        # Measure, and let a measurement failure FAIL THE ROW rather than
        # leaving a complete-looking row with no residence in it.
        row.update(measure_residence(md_wd / "rep1"))
        log.info("  residence: ligand RMSD %.3f nm mean, engaged in %.0f%% of frames",
                 row["explicit_ligand_rmsd_nm_mean"],
                 100 * row["explicit_frac_frames_engaged"])
    except Exception as exc:                              # noqa: BLE001
        row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:160]}"
        log.warning("  FAILED: %s", row["status"])
    finally:
        if not keep:
            shutil.rmtree(wd / "md" / "rep1" if (wd / "md" / "rep1").exists() else wd,
                          ignore_errors=True)
    return row


def candidate_from_frames(cid: str) -> tuple[str, int]:
    """(canonical SMILES, formal charge at pH 7.4) for a named candidate.

    Read from the newest T_3/T_4 frame rather than passed on the command line,
    so the molecule simulated is the molecule the ranking ranked. `charge_ph74`
    is returned alongside because antechamber's `-nc` must match the species,
    and the non-covalent path's default of 0 is a guess (D0058).
    """
    import pandas as _pd
    import export_nac_poses as enp
    for subdir, stem in enp.FRAMES.values():
        f = enp.latest_frame(subdir, stem)
        df = _pd.read_parquet(f)
        hit = df[df["candidate_id"].astype(str) == cid]
        if len(hit):
            r = hit.iloc[0]
            q = r.get("charge_ph74", 0)
            log.info("  %s found in %s: charge_ph74=%s", cid, f.name, q)
            return str(r["canonical_smiles"]), int(0 if _pd.isna(q) else q)
    raise ResidenceError(f"{cid} is in no current frame")


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
    ap.add_argument("--candidate", help="run ONE named candidate instead of the "
                                        "validation set")
    ap.add_argument("--pose", help="start from this saved pose (SDF or PDBQT) "
                                   "instead of re-docking; implies --candidate")
    ap.add_argument("--pose-rank", type=int, default=1,
                    help="which pose_rank to take from a multi-model SDF")
    ap.add_argument("--net-charge", type=int, default=None,
                    help="ligand formal charge for antechamber; default is the "
                         "frame's charge_ph74 in --candidate mode, else 0")
    ap.add_argument("--tag", default=None, help="output stem suffix")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.candidate:
        smiles, q74 = candidate_from_frames(args.candidate)
        pose = Path(args.pose) if args.pose else None
        if pose and not pose.is_file():
            raise SystemExit(f"no pose at {pose}")
        nc = args.net_charge if args.net_charge is not None else q74
        if nc != q74:
            log.warning("net charge %+d OVERRIDES the frame's charge_ph74 %+d", nc, q74)
        log.info("single candidate %s at %.0f ps, net charge %+d, pose %s",
                 args.candidate, args.production_ps, nc, pose or "(re-dock)")
        row = run_one(args.candidate, smiles, "candidate",
                      production_ps=args.production_ps, nrun=args.nrun,
                      gpu=args.gpu, keep=args.keep, pose=pose,
                      pose_rank=args.pose_rank, net_charge=nc)
        df = pd.DataFrame([row])
        dest = OUT.write(f"md_residence_{args.tag or args.candidate}", ".csv")
        df.to_csv(dest, index=False)
        print(f"\n  {row['status']} -> {dest}")
        for k, v in row.items():
            print(f"    {k} = {v}")
        return

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
