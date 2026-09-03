"""
Purpose: Non-covalent MM-GBSA rescoring — the T5 tier for T_1 and T_2.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: a docked Vina-GPU pose (PDBQT) + the prepared receptor
Output: dG_bind and its per-term breakdown, per candidate

THE ADAPTER T_2's PLAN NAMED. T_2's documented pipeline lists "T5 · Physics
rescoring on survivors" as step 8, with the MM-GBSA adapter "specified with its
interface fixed but not yet implemented". This is that adapter. It is
deliberately the SAME three-leg protocol, the same implicit-solvent model and
the same minimiser as the covalent version, so a T_2 number and a T_4 number
differ only where the chemistry differs.

SIMPLER THAN THE COVALENT CASE, NOT HARDER. Everything that makes
`shared.mmgbsa` intricate is covalent-specific: the CYX residue, the link-atom
cut, the hand-derived GAFF2 junction parameters, the capping hydrogen and the
charge correction. None of it applies here. What is reused verbatim is the part
that must not diverge — receptor preparation (including the histidine states
`reduce` assigned), the minimisation input, the energy parsing and dG.

ONE REAL ADVANTAGE OVER T_4. The covalent dG carries a constant term for the
bond formed and the bonds broken, which cancels only within a warhead class
(D0020, D0023). There is no such term here, so a non-covalent dG IS comparable
across the whole approach. T_2 may rank on it globally; T_4 may not.

WHAT THIS IS NOT. Single-structure MM-GBSA on a minimised geometry — no
ensemble, no explicit solvent, no entropy. The plan's "ΔG_bind ± uncertainty"
needs the short explicit-solvent MD listed beside it in the same step, which
does not exist in this repo. A single minimisation has no error bar, so this
module reports no uncertainty rather than inventing one.

READ IT WITH D0031. The candidates reaching this tier were selected by docking,
and on class-matched decoys docking sits at chance on this receptor. So this is
not confirmation of a docking ranking; it is an INDEPENDENT estimate that may
disagree with it, and the disagreement is the informative part.
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from . import mmgbsa as mg

log = logging.getLogger(__name__)

OBABEL = "/data/lab_vm/envs/dwi_cheminf/bin/obabel"


def pose_to_sdf(pose_pdbqt: Path, workdir: Path, smiles: str) -> Path:
    """First (best) pose from a Vina-GPU PDBQT, as a chemically correct molecule.

    COORDINATES FROM THE POSE, BOND ORDERS FROM THE SMILES. PDBQT does not
    record bond orders or formal charges — it carries atom types for scoring, not
    chemistry. Letting the converter re-perceive them produced aromatic carbons
    with three connections and no hydrogen, and antechamber refused the file:
    "Weird atomic valence (3) for atom C1". Any parameterisation built from a
    PDBQT alone is guessing at the molecule it is parameterising.

    So the docked heavy-atom geometry is kept and the connectivity is imposed
    from the candidate's own canonical SMILES, which the frame already carries.
    Hydrogens are then added at that geometry rather than by re-embedding, which
    would discard the pose this stage exists to score.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # VIA PDB, NOT SDF. obabel's SDF writer commits to bond orders it inferred
    # from the PDBQT and gets them wrong — it produced a three-valent oxygen,
    # which RDKit rejects before the template can correct anything. A PDB carries
    # connectivity without bond orders, which is exactly the input
    # AssignBondOrdersFromTemplate is designed for. Same pose either way; the
    # difference is whether the converter is allowed to guess.
    raw = workdir / "pose_raw.pdb"
    cmd = [OBABEL, "-ipdbqt", str(pose_pdbqt), "-opdb", "-O", str(raw),
           "-f", "1", "-l", "1"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if not raw.is_file() or raw.stat().st_size == 0:
        raise mg.MMGBSAError(
            f"could not convert {pose_pdbqt.name} to PDB: "
            f"{(proc.stderr or proc.stdout)[-300:]}")

    pose = Chem.MolFromPDBFile(str(raw), sanitize=False, removeHs=False)
    if pose is None:
        raise mg.MMGBSAError(f"RDKit could not read the converted pose {raw.name}")

    # STRIP THE POSE'S HYDROGENS BEFORE MATCHING, and do it explicitly.
    #
    # Vina keeps polar hydrogens, so a carboxylic acid pose carries its -OH
    # hydrogen. The template is the species at pH 7.4 — a carboxylATE, whose
    # oxygen has no hydrogen and a -1 charge — and imposing that on an oxygen
    # that still holds an H is a valence contradiction: "Explicit valence for
    # atom # 23 O, 2, is greater than permitted". Every T_1 and T_2 candidate
    # failed this way the moment protonation was switched on.
    #
    # `removeHs=True` on the reader does NOT help: it is a no-op under
    # sanitize=False, which is why the atom count was unchanged. The hydrogens
    # are therefore deleted by hand, the heavy-atom skeleton is matched to the
    # template, and hydrogens are put back afterwards at the template's
    # protonation state.
    rw = Chem.RWMol(pose)
    for idx in sorted([a.GetIdx() for a in pose.GetAtoms()
                       if a.GetAtomicNum() == 1], reverse=True):
        rw.RemoveAtom(idx)
    pose = rw.GetMol()
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise mg.MMGBSAError(f"unparseable candidate SMILES {smiles!r}")
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(template, pose)
    except Exception as exc:  # noqa: BLE001
        raise mg.MMGBSAError(
            f"could not map the docked pose onto its own SMILES ({exc}); the "
            "pose and the frame disagree about the molecule") from exc

    fixed = Chem.AddHs(fixed, addCoords=True)
    out = workdir / "ligand_pose.sdf"
    w = Chem.SDWriter(str(out))
    w.write(fixed)
    w.close()
    return out


def parameterize_ligand(pose_sdf: Path, workdir: Path,
                        net_charge: int = 0) -> tuple[Path, Path]:
    """AM1-BCC/GAFF2 parameters for the docked pose. No cap, no correction.

    `net_charge` MUST be the ligand's formal charge. antechamber assigns AM1-BCC
    charges under the total it is given, so a carboxylate passed as neutral
    silently redistributes a whole electron across the molecule.
    """
    mol2 = workdir / "ligand.mol2"

    # ---- AM1-BCC CHARGES ARE CACHED PER MOLECULE; COORDINATES NEVER ARE ----
    #
    # MEASURED 2026-09-02: this call is 98-134 s of a 256-286 s sweep, ~40%.
    # It runs once per (molecule, MODE) and the campaign has ~9 modes per
    # molecule, so the AM1 semi-empirical step is being redone ~9x for an
    # answer that does not depend on which mode is being simulated.
    #
    # WHAT IS SAFE TO REUSE AND WHAT IS NOT. `tleap` takes the ligand's
    # starting COORDINATES from this mol2 (`LIG = loadmol2`), so caching the
    # mol2 itself would start every mode of a molecule from whichever mode
    # populated the cache -- a different pose of the right molecule, silently,
    # with a plausible stability reported for it. Only the CHARGES are reused:
    # `-c rc` reads them from a file while the geometry still comes from THIS
    # pose's SDF.
    #
    # THE KEY IS THE MOLECULE'S IDENTITY AND ITS CHARGE STATE, and what it
    # omits is stated because this repo has been bitten five times by caches
    # keyed on less than their inputs (catalogue 8, 9, 18, 22, and the
    # `lru_cache` on nothing at all). It omits the CONFORMER, deliberately:
    # AM1-BCC charges are conformation-dependent in principle, and reusing one
    # conformer's charges across a molecule's poses is the standard practice
    # this trades a second-order effect for a 40% saving. It does NOT omit the
    # net charge, which changes the charges outright.
    #
    # And the element sequence is verified on every hit, because a charge file
    # is a bare list in atom order: if the order differed, every atom would get
    # someone else's charge and nothing would raise.
    key = _charge_cache_key(pose_sdf, net_charge)
    cached = _charge_cache_dir() / key
    charges = cached / "charges.dat"
    elements = cached / "elements.txt"
    want_elems = _element_sequence(pose_sdf)

    hit = (charges.is_file() and elements.is_file()
           and elements.read_text().split() == want_elems)
    if hit:
        log.info("charges from cache %s (skipping AM1-BCC)", key[:12])
        mg._run([mg._amber("antechamber"),
                 "-i", str(pose_sdf), "-fi", "sdf",
                 "-o", str(mol2), "-fo", "mol2",
                 "-c", "rc", "-cf", str(charges),
                 "-nc", str(int(net_charge)), "-s", "2",
                 "-at", "gaff2", "-pf", "y"],
                workdir, "antechamber.log", timeout=3600)
    else:
        if elements.is_file():
            log.warning("charge cache %s has a different atom order; "
                        "recomputing rather than mis-assigning", key[:12])
        mg._run([mg._amber("antechamber"),
                 "-i", str(pose_sdf), "-fi", "sdf",
                 "-o", str(mol2), "-fo", "mol2",
                 "-c", "bcc", "-nc", str(int(net_charge)), "-s", "2",
                 "-at", "gaff2", "-pf", "y"],
                workdir, "antechamber.log", timeout=3600)
    if not mol2.is_file() or mol2.stat().st_size == 0:
        raise mg.MMGBSAError(
            f"antechamber produced no mol2; see {workdir/'antechamber.log'}")

    if not hit:
        # Populate the cache only from a real AM1-BCC run, never from a `-c rc`
        # one -- that would re-store the charges it just read, and a corrupted
        # entry would then look self-consistent forever.
        _store_charges(mol2, pose_sdf, net_charge)

    frcmod = workdir / "ligand.frcmod"
    mg._run([mg._amber("parmchk2"), "-i", str(mol2), "-f", "mol2",
             "-o", str(frcmod), "-s", "gaff2"],
            workdir, "parmchk2.log")
    log.info("ligand parameterised (non-covalent), net charge %+d", net_charge)
    return mol2, frcmod


def _charge_cache_dir():
    """Scratch, never the append-only root -- this is a derived convenience."""
    from pathlib import Path as _P
    d = _P("/data/lab_vm/modifiable/inhibition/ligand_param_cache")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _element_sequence(sdf) -> list[str]:
    """Element symbols in the SDF's atom order -- what a charge file indexes."""
    from rdkit import Chem as _C
    mols = [m for m in _C.SDMolSupplier(str(sdf), removeHs=False) if m]
    if not mols:
        raise mg.MMGBSAError(f"no readable molecule in {sdf}")
    return [a.GetSymbol() for a in mols[0].GetAtoms()]


def _charge_cache_key(sdf, net_charge: int) -> str:
    """Identity of the MOLECULE plus its charge state -- never the pose.

    From the canonical SMILES, so two conformers of one molecule key the same
    and two different molecules cannot collide however similar their poses.
    """
    import hashlib
    from rdkit import Chem as _C
    mols = [m for m in _C.SDMolSupplier(str(sdf), removeHs=False) if m]
    if not mols:
        raise mg.MMGBSAError(f"no readable molecule in {sdf}")
    smi = _C.MolToSmiles(_C.RemoveHs(mols[0]))
    return hashlib.sha256(f"{smi}|{int(net_charge)}".encode()).hexdigest()[:24]


def _store_charges(mol2, sdf, net_charge: int) -> None:
    """Write this molecule's AM1-BCC charges to the cache, in atom order."""
    key = _charge_cache_key(sdf, net_charge)
    d = _charge_cache_dir() / key
    d.mkdir(parents=True, exist_ok=True)
    q, in_atoms = [], False
    for ln in mol2.read_text(errors="replace").splitlines():
        if ln.startswith("@<TRIPOS>ATOM"):
            in_atoms = True; continue
        if ln.startswith("@<TRIPOS>"):
            in_atoms = False; continue
        f = ln.split()
        if in_atoms and len(f) >= 9:
            q.append(float(f[8]))
    if not q:
        return
    # antechamber's -cf format: whitespace-separated, 8 per line.
    (d / "charges.dat").write_text(
        "\n".join(" ".join(f"{x:10.6f}" for x in q[i:i + 8])
                   for i in range(0, len(q), 8)) + "\n")
    (d / "elements.txt").write_text(" ".join(_element_sequence(sdf)))


def build_topologies(workdir: Path, mol2: Path, frcmod: Path, receptor_pdb: Path
                     ) -> dict[str, tuple[Path, Path]]:
    """Complex, receptor and ligand topologies in one tleap pass.

    The complex is a plain `combine` — no `bond`, because nothing is bonded.
    That single missing line is the whole difference from the covalent build,
    and with it go the junction parameters and the capping logic.
    """
    script = f"""source {mg.PROTEIN_FF}
source {mg.LIGAND_FF}
loadamberparams {frcmod.name}
set default PBRadii {mg.PB_RADII}

LIG = loadmol2 {mol2.name}
saveamberparm LIG ligand.prmtop ligand.inpcrd

rec = loadpdb {receptor_pdb.name}
saveamberparm rec receptor.prmtop receptor.inpcrd

cpx = combine {{rec LIG}}
saveamberparm cpx complex.prmtop complex.inpcrd
quit
"""
    (workdir / "build.leap").write_text(script, encoding="utf-8")
    mg._run([mg._amber("tleap"), "-f", "build.leap"], workdir, "tleap.log")

    legs = {}
    for leg in ("complex", "receptor", "ligand"):
        top, crd = workdir / f"{leg}.prmtop", workdir / f"{leg}.inpcrd"
        if not top.is_file() or top.stat().st_size == 0:
            raise mg.MMGBSAError(
                f"tleap produced no usable {leg} topology; "
                f"see {workdir/'tleap.log'}")
        legs[leg] = (top, crd)

    missing = [l for l in (workdir / "tleap.log").read_text().splitlines()
               if "Could not find" in l]
    if missing:
        raise mg.MMGBSAError(
            f"{len(missing)} missing force-field parameter(s), first: "
            f"{missing[0].strip()}")
    return legs


def verify_complex(legs: dict[str, tuple[Path, Path]]) -> dict:
    """Atom-count bookkeeping: the complex must be receptor plus ligand exactly.

    The covalent build asserts a specific SG-C bond exists. There is no bond to
    check here, so what is checked instead is that nothing was lost or
    duplicated by `combine` — a silently truncated complex would still minimise
    and would still return a plausible dG.
    """
    import parmed

    n = {leg: len(parmed.load_file(str(top)).atoms)
         for leg, (top, _) in legs.items()}
    if n["complex"] != n["receptor"] + n["ligand"]:
        raise mg.MMGBSAError(
            f"complex has {n['complex']} atoms but receptor + ligand is "
            f"{n['receptor']} + {n['ligand']} = {n['receptor'] + n['ligand']}; "
            "combine lost or duplicated atoms")
    info = {"n_complex": n["complex"], "n_receptor": n["receptor"],
            "n_ligand": n["ligand"]}
    log.info("complex verified (non-covalent): %s", info)
    return info


def delta_g(legs: dict[str, mg.LegEnergies]) -> dict:
    """dG = G_complex - G_receptor - G_ligand, with the per-term breakdown."""
    out = mg.delta_g(legs)
    # The covalent scheme's caveats do not apply, and leaving them attached
    # would understate what this number supports.
    out["scheme"] = "non-covalent 3-leg, no link atom"
    out["comparable"] = ("across the whole approach — there is no per-class "
                         "constant bond term to cancel (contrast D0020)")
    # No junction exists here, so the junction parameters are not an input to
    # this number and must not sit in its cache fingerprint — otherwise every
    # non-covalent candidate rescores whenever the covalent frcmod changes. Set
    # to None rather than deleted, so the key is still PRESENT: in
    # mg.cached_result_is_current a MISSING marker means "written before
    # fingerprinting existed, treat as stale", and an absent key here would
    # make every non-covalent result permanently un-cacheable.
    out["junction_frcmod"] = None
    return out
