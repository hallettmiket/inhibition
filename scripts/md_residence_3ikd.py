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
import json
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
from shared import mode_key                      # noqa: E402
from shared import target_config as tc           # noqa: E402
from shared import mmgbsa as mg                  # noqa: E402
from shared import mmgbsa_noncovalent as mgn     # noqa: E402
from shared import outputs as sout               # noqa: E402
from shared import run_paths as rp        # noqa: E402
import nac_screen as ns                          # noqa: E402
import nac_robustness as rb                      # noqa: E402

log = logging.getLogger("md-residence")

OUT = sout.Topic("blacksmith", rp.residence_topic())
WORK = rp.residence_work()
PLAIN_REC = rp.receptor_plain()

# The receptor as the chemist prepared it — the SAME file every 3IKD measurement
# on this branch uses, so residence is comparable with the docking that selected
# the candidate (D0059).
RECEPTOR_3IKD = rp.receptor_prep()
#: SMILES + pH 7.4 charge for poses that are not library candidates (crystal
#: controls, references). One JSON per ident, written beside the pose.
POSE_SIDECARS = rp.sidecars()

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
        # A CRYSTAL-DERIVED POSE HAS NO HYDROGENS AND ONLY CONNECTIVITY, so a
        # sanitising read assigns implicit hydrogens from single-bond valence --
        # a sulfone oxygen becomes a hydroxyl, and imposing the template's double
        # bond then exceeds its valence. Re-read unsanitised so bond orders come
        # from the template rather than being argued with.
        mols = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False,
                                              sanitize=False) if m]
    if not mols:
        raise ResidenceError(f"no readable pose in {pose_sdf}")
    if not any(a.GetAtomicNum() == 1 for a in mols[0].GetAtoms()):
        mols = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False,
                                              sanitize=False) if m]

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
    # A CRYSTAL POSE CARRIES NO HYDROGENS AT ALL, and restating "zero observed"
    # for every heavy atom is then a statement about the deposition, not about
    # the molecule -- the rebuild comes back as C11NO3S against the template's
    # C11H21NO3S and the formula guard correctly refuses it. Where the pose has
    # no hydrogens anywhere, let the template supply them instead: the count is
    # then derived from the SMILES we are matching against, which is the only
    # source of truth available. Poses that DO carry hydrogens keep the restated
    # counts, because that is what protects aromatic nitrogens (above).
    crystal_like = not any(a.GetAtomicNum() == 1 for a in pose.GetAtoms())
    if crystal_like:
        log.info("  pose carries no hydrogens (crystal-derived); "
                 "taking hydrogen counts from the template")
    rw = Chem.RWMol(pose)
    for idx in sorted([a.GetIdx() for a in pose.GetAtoms()
                       if a.GetAtomicNum() == 1], reverse=True):
        rw.RemoveAtom(idx)
    if crystal_like:
        # The SDF write stamps `noImplicit`, so AddHs later adds nothing and the
        # rebuild comes back as C11NO3S. Clearing it lets valence supply the
        # hydrogens, which for a hydrogen-free deposition is the only honest
        # source -- and the formula guard below still has to agree with the
        # template before anything reaches antechamber.
        for a in rw.GetAtoms():
            a.SetNoImplicit(False)
            a.SetNumExplicitHs(0)
        rw.UpdatePropertyCache(strict=False)
    else:
        for a in rw.GetAtoms():
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(n_hydrogens[a.GetIdx()])
    Chem.SanitizeMol(rw)
    # THE TEMPLATE MUST BE THE SPECIES THAT WAS DOCKED, NOT THE NEUTRAL ONE.
    #
    # Since D0074 the screen protonates at pH 7.4 before docking, so the saved
    # pose is the ionised form. The frame's `canonical_smiles` is still the
    # neutral molecule, and comparing the two made the formula guard below fire
    # on every ionisable candidate: "pose rebuilt as C18H22BrN4O3S+ but the
    # frame's SMILES is C18H21BrN4O3S". That is a real mismatch and the guard was
    # right to refuse -- but the fix is to compare against the right species, not
    # to weaken the check.
    #
    # It cost 69 of 131 sweeps in the 2026-08-08 run, and not at random: 85% of
    # bdhi_c5 and 68% of bdhi_c4 against 15% of acrylamide, because the
    # bromo-dihydroisoxazoles carry an ionisable centre. Two of the three
    # priority warhead classes were nearly erased by a units mismatch.
    #
    # Protonated with the SAME function the screen used, so the two cannot drift.
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ResidenceError(f"unparseable candidate SMILES {smiles!r}")
    pose_charge = Chem.GetFormalCharge(rw.GetMol())
    if pose_charge != Chem.GetFormalCharge(template):
        try:
            from shared import ionisation as ion
            prot = ion.protonate({"x": smiles}).get("x")
        except Exception as exc:                          # noqa: BLE001
            raise ResidenceError(
                f"pose carries charge {pose_charge:+d} but the frame's SMILES is "
                f"neutral, and it could not be protonated to compare ({exc})"
            ) from exc
        cand = prot and Chem.MolFromSmiles(prot)
        if cand is None or Chem.GetFormalCharge(cand) != pose_charge:
            raise ResidenceError(
                f"pose carries charge {pose_charge:+d}; protonating the frame's "
                f"SMILES at pH 7.4 gave {prot!r}, which does not match. The pose "
                "and the frame disagree about which species this is.")
        template = cand
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(template, rw.GetMol())
    except Exception as exc:                              # noqa: BLE001
        raise ResidenceError(
            f"saved pose does not match the frame's SMILES ({exc})") from exc

    fixed.UpdatePropertyCache(strict=False)
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
    # WHICH MODE WAS SIMULATED, written beside the trajectory it produced.
    #
    # The result row records `pose_path` but not which pose inside that file was
    # taken, and a v3 pose SDF holds one pose PER MODE. So a finished 100 ns run
    # could not be joined to the 10 ns sweep, which is scored per mode -- the
    # bare ident in the sweep table is mode 0, and at least one molecule reached
    # 100 ns on a minority mode (#46). Matching on the bare ident silently
    # compares one mode's trajectory against a different mode's sweep.
    #
    # This blocked #35's re-derivation and #36's false-negative rate, both for
    # the same reason: the run did not say which pose it ran. It says so now.
    (wd / "pose_provenance.json").write_text(json.dumps({
        "pose_sdf": str(pose_sdf),
        "pose_rank": int(pose_rank),
        "mode": (int(pose.GetProp("mode")) if pose.HasProp("mode") else None),
        "n_poses_in_file": len(mols),
    }, indent=2), encoding="utf-8")
    log.info("  pose: rank %d of %d from %s (mode %s)", pose_rank, len(mols),
             pose_sdf.name,
             pose.GetProp("mode") if pose.HasProp("mode") else "unrecorded")
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

    # THE MOL2 MUST HOLD THE POSE WE ASKED FOR.
    #
    # `tleap` takes the ligand's starting COORDINATES from `ligand.mol2`
    # (`LIG = loadmol2`), and antechamber writes that mol2 from whatever SDF it
    # was handed. So anything that supplies a mol2 from elsewhere -- a cache
    # keyed on the molecule rather than the pose, a stale workdir, a copied
    # directory -- starts the simulation from a DIFFERENT pose of the right
    # molecule and reports a perfectly plausible stability for it.
    #
    # Nothing checked this. It is checked now, cheaply, against the SDF this
    # build was given, because the failure is silent by construction: both
    # poses are real, both are of this molecule, and every downstream artefact
    # would look normal.
    _assert_mol2_matches_pose(wd / "ligand.mol2", sdf)
    log.info("  parameterised: %s, %s", mol2.name, frcmod.name)
    return wd


def _assert_mol2_matches_pose(mol2: Path, sdf: Path, tol_a: float = 0.05) -> None:
    """Heavy-atom coordinates in the mol2 must be the pose's, within `tol_a`.

    Compared as SORTED coordinate rows rather than by atom index: antechamber
    may reorder atoms, and an index-wise comparison would fail for a correct
    build. What must hold is that the SET of heavy-atom positions is the pose's.

    Raises rather than warning. A wrong pose here is not a degraded result, it
    is a different experiment wearing this one's name.
    """
    from rdkit import Chem as _Chem
    import numpy as _np

    mols = [m for m in _Chem.SDMolSupplier(str(sdf), removeHs=False) if m]
    if not mols:
        raise ResidenceError(f"cannot verify pose: no molecule in {sdf.name}")
    conf = mols[0].GetConformer()
    want = _np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                       conf.GetAtomPosition(i).z]
                      for i, a in enumerate(mols[0].GetAtoms())
                      if a.GetAtomicNum() > 1])

    got, in_atoms = [], False
    for ln in mol2.read_text(errors="replace").splitlines():
        if ln.startswith("@<TRIPOS>ATOM"):
            in_atoms = True; continue
        if ln.startswith("@<TRIPOS>"):
            in_atoms = False; continue
        if in_atoms and len(ln.split()) >= 6:
            f = ln.split()
            if not f[5].upper().startswith("H"):          # gaff2 type, not element
                got.append([float(f[2]), float(f[3]), float(f[4])])
    got = _np.array(got)

    if len(got) != len(want):
        raise ResidenceError(
            f"{mol2.name} has {len(got)} heavy atoms, the pose {sdf.name} has "
            f"{len(want)} — the mol2 was not built from this pose")
    a = want[_np.lexsort(want.T)]
    b = got[_np.lexsort(got.T)]
    worst = float(_np.abs(a - b).max())
    if worst > tol_a:
        raise ResidenceError(
            f"{mol2.name} coordinates differ from {sdf.name} by up to "
            f"{worst:.2f} A — the simulation would start from a DIFFERENT pose "
            f"of this molecule than the one requested. Refusing to run.")


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
            work_root: Path | None = None, replicate: int = 1,
            stop_after: str | None = None,
            reuse_equilibration: bool = False) -> dict:
    wd = (work_root or WORK) / ident.replace(":", "_")
    wd.mkdir(parents=True, exist_ok=True)
    # THE ROW CARRIES THE PAIR, NOT JUST THE LABEL. `ident` may now be
    # `<parent>_m<mode>`, and every table that joins these rows to the sweep or
    # the ranking joins on (parent_ident, mode) -- see shared/mode_key, which
    # exists because a merge on `ident` silently dropped every mode-0 row. Being
    # written here means a consumer never has to re-parse the label to find out
    # what the run was of.
    _parent, _mode = mode_key.split_ident(ident)
    row = {"ident": ident, "parent_ident": _parent, "mode": _mode,
           "label": label, "smiles": smiles,
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
                              production_ps=production_ps, replicate=replicate,
                              candidate_id=ident,
                              stop_after=stop_after,
                              reuse_equilibration=reuse_equilibration)
        if stop_after:
            # NO ROW. A row means a measured result, and there is nothing
            # measured here -- writing one would be the `stage0 only` defect
            # again, where a placeholder in the results table was counted as a
            # finished sweep.
            log.info("  stopped after %s; no row written (caller decides)",
                     stop_after)
            return {"status": f"stopped after {stop_after}",
                    "equilibration_dir": res.get("equilibration_dir")}
        row.update({k: v for k, v in res.items() if isinstance(v, (int, float, str))})
        log.info("  production done (%.0f ps)", production_ps)
        # Measure, and let a measurement failure FAIL THE ROW rather than
        # leaving a complete-looking row with no residence in it.
        #
        # THE RUN SAYS WHERE IT WROTE. run_pipeline puts replicate N in rep<N>,
        # but this read "rep1" regardless: replicate 2 ran its full 100 ns and
        # then failed measurement against a directory it had never written, so a
        # finished trajectory produced a row saying "failed". Same defect class
        # as the rest of this project -- a path taken by DEFAULT rather than by
        # identity -- and it arrived with the --replicate flag that fixed the
        # seed. Fixing one half of "replicate" and not the other is how.
        rep_dir = Path(res.get("equilibration_dir") or (md_wd / f"rep{replicate}"))
        # Carry the pose provenance onto the ROW, not only into the workdir --
        # the row is what survives into 00_outputs and what any later join reads.
        prov = wd / "pose_provenance.json"
        if prov.is_file():
            row.update({f"pose_{k}" if k != "pose_sdf" else k: v
                        for k, v in json.loads(prov.read_text()).items()})
        row.update(measure_residence(rep_dir))
        log.info("  residence: ligand RMSD %.3f nm mean, engaged in %.0f%% of frames",
                 row["explicit_ligand_rmsd_nm_mean"],
                 100 * row["explicit_frac_frames_engaged"])
    except Exception as exc:                              # noqa: BLE001
        row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:160]}"
        log.warning("  FAILED: %s", row["status"])
    finally:
        if not keep:
            r = wd / "md" / f"rep{replicate}"
            shutil.rmtree(r if r.exists() else wd, ignore_errors=True)
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
    # REFERENCE MOLECULES ARE NOT CANDIDATES. Crystal controls (#48) and any
    # other externally-sourced pose live in no T_3/T_4 frame by construction, so
    # the lookup that keeps a CANDIDATE honest -- simulate what the ranking
    # ranked -- has nothing to find. A sidecar written beside the pose carries
    # the SMILES and pH 7.4 charge for exactly these, and is consulted only after
    # the frames have been tried, so it can never shadow a real candidate.
    side = _sidecar(cid)
    if side is not None:
        log.info("  %s from pose sidecar: charge_ph74=%s", cid, side[1])
        return side
    raise ResidenceError(f"{cid} is in no current frame and has no pose sidecar")


def _sidecar(cid: str):
    """(smiles, charge) for a non-candidate pose, or None."""
    import json as _json
    for d in (POSE_SIDECARS,):
        f = d / f"{cid}.json"
        if f.is_file():
            try:
                j = _json.loads(f.read_text())
                return str(j["canonical_smiles"]), int(j.get("charge_ph74", 0))
            except Exception as exc:                       # noqa: BLE001
                raise ResidenceError(f"{cid}: unreadable sidecar {f} ({exc})")
    return None


def build_set(n_neg_per_class: int) -> list[ns.Candidate]:
    saved = rb.NEG_PER_CLASS
    try:
        rb.NEG_PER_CLASS = n_neg_per_class
        return rb.build_set()
    finally:
        rb.NEG_PER_CLASS = saved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    # THE SPEC LIVES IN config/target.yaml, NOT IN THIS DEFAULT. It read 100.0
    # ps -- a thousandth of `md.production_ps` -- while the only reader of that
    # key was `pipeline_schematic`, the DIAGRAM. So the GUI stated 100 ns, this
    # stated 0.1 ns, and the number that actually ran came from a flag typed
    # into a scratch shell script. The flag still wins when passed; what changes
    # is that omitting it now means the spec rather than a placeholder.
    ap.add_argument("--production-ps", type=float, default=tc.md_production_ps(),
                    help="production length in ps; defaults to md.production_ps "
                         "in config/target.yaml. Pass a small value explicitly "
                         "to prove the chain cheaply -- the point is that a "
                         "short run is now something you ASK for rather than "
                         "something you get by forgetting a flag.")
    ap.add_argument("--limit", type=int, default=1,
                    help="how many molecules; 1 verifies the plumbing")
    ap.add_argument("--n-neg", type=int, default=5, help="negatives per class")
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--keep", action="store_true", help="keep MD workdirs")
    # THE SWEEP AND THE FULL RUN MUST NOT SHARE A DIRECTORY.
    #
    # The workdir is <root>/<candidate>/md/rep1 regardless of --tag, so a 10 ns
    # sweep and a 100 ns run of the SAME molecule would land on top of each
    # other -- and the second would find a prod.xtc already there. `work_root`
    # existed in run_one and was never exposed, so there was no way to say it.
    ap.add_argument("--work-root", default=None,
                    help="override the MD working root (keeps a 10 ns sweep and "
                         "a 100 ns run of the same molecule apart)")
    ap.add_argument("--candidate", help="run ONE named candidate instead of the "
                                        "validation set")
    ap.add_argument("--pose", help="start from this saved pose (SDF or PDBQT) "
                                   "instead of re-docking; implies --candidate")
    ap.add_argument("--pose-rank", type=int, default=1,
                    help="which pose_rank to take from a multi-model SDF")
    # POSE RANK AND MODE ARE NOT THE SAME NUMBER and neither implies the other:
    # a mode is a cluster of poses, so mode 3 of one molecule may be pose_rank 4
    # and of another pose_rank 9. Deriving one from the other is a guess; the
    # caller knows both, so it states both.
    ap.add_argument("--mode", type=int, default=None,
                    help="binding mode id; makes the run ident <candidate>_m<N> "
                         "so two modes of one molecule get separate workdirs, "
                         "rows and reports. Omit for molecule-level behaviour.")
    ap.add_argument("--net-charge", type=int, default=None,
                    help="ligand formal charge for antechamber; default is the "
                         "frame's charge_ph74 in --candidate mode, else 0")
    # REPLICATE IS A REAL PARAMETER, NOT A DIRECTORY. gromacs_explicit derives a
    # distinct, reproducible velocity seed per (candidate, replicate); the runner
    # hardcoded replicate=1, so re-running into a different --work-root would have
    # reproduced the SAME trajectory bit for bit and called it a replicate. A
    # sibling that cannot differ is not a replicate, it is a copy.
    ap.add_argument("--replicate", type=int, default=tc.md_replicates(),
                    help="replicate index; selects the velocity seed (D0038)")
    ap.add_argument("--tag", default=None, help="output stem suffix")
    # TWO-PHASE MD, so a caller can look at the pose after equilibration and
    # decide whether production is worth running (2026-09-02). Both default off,
    # so every existing caller behaves exactly as before.
    #
    # WHY THIS IS WORTH A FLAG. Measured over the first nac_v8 sweeps, the
    # DOCKED geometry does not predict the outcome -- five modes all at
    # 2.78-2.97 A gave 0.000 to 0.926 attack-ready -- while the distance AFTER
    # the 300 ps unrestrained equilibration does (3.58 A -> 0.926, 6.46 A ->
    # 0.000). Equilibration is 300 ps of the 1,500 ps a sweep runs, so a pose
    # that has already left can be dropped for a fifth of the cost.
    ap.add_argument("--stop-after", default=None, choices=("npt",),
                    help="end after equilibration and write NO row; the caller "
                         "inspects the frame and decides about production")
    ap.add_argument("--reuse-equilibration", action="store_true",
                    help="continue from an existing npt frame instead of "
                         "re-equilibrating (pair with an earlier --stop-after)")
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
        # THE RUN IS IDENTIFIED BY THE MODE IT RAN, NOT BY THE MOLECULE.
        #
        # `ident` reaches three places that all assumed one run per molecule: the
        # workdir (<work_root>/<ident>), the row written to md_residence, and the
        # report filename. So two modes of one molecule shared a directory --
        # and build_workdir rebuilds in place, so the second overwrote the
        # first's finished trajectory while its row survived. t4_c8c3aec07421
        # (_m1 and _m5) was queued to do exactly that.
        #
        # `--mode` makes the mode part of the identity instead of a column
        # beside it, which is what `shared/mode_key` says an ident is for: the
        # display label of a (parent, mode) pair. The workdir separates itself,
        # the rows key per mode, and the results rail can be per mode rather
        # than one row per molecule showing its best. @tt8804: "it doesnt show
        # sub modes only molecules".
        #
        # Omitting it keeps the old molecule-level behaviour, so every existing
        # invocation and all 63 legacy rows still mean what they meant.
        ident = args.candidate
        if args.mode is not None:
            ident = f"{args.candidate}_m{int(args.mode)}"
        log.info("single candidate %s at %.0f ps, net charge %+d, pose %s",
                 ident, args.production_ps, nc, pose or "(re-dock)")
        row = run_one(ident, smiles, "candidate",
                      production_ps=args.production_ps, nrun=args.nrun,
                      gpu=args.gpu, keep=args.keep, pose=pose,
                      pose_rank=args.pose_rank, net_charge=nc,
                      work_root=Path(args.work_root) if args.work_root else None,
                      replicate=args.replicate,
                      stop_after=args.stop_after,
                      reuse_equilibration=args.reuse_equilibration)
        if args.stop_after:
            # NO ROW FOR A HALF-RUN. `run_one` already declined to measure; if
            # this wrote the row anyway the residence table would carry a
            # `stopped after npt` entry with no residence in it, and anything
            # resuming on ident would count it as a finished 100 ns run. That is
            # the `stage0 only` defect exactly, one table over.
            print(f"\n  {row['status']}; no row written")
            print(f"  equilibration at {row.get('equilibration_dir')}")
            return
        df = pd.DataFrame([row])
        dest = OUT.write(f"md_residence_{args.tag or ident}", ".csv")
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
