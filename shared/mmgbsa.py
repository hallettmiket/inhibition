"""
Purpose: Three-leg MM-GBSA on the TRUE covalent adduct (Cys113-S bonded to the ligand).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: a docked adduct pose (SDF) + the prepared receptor
Output: complex/receptor/ligand topologies, minimised structures, and dG

WHAT dG MEANS HERE. A covalent complex has no ligand to remove, so the usual
G_complex - G_receptor - G_ligand is undefined until you say where to cut. This
module cuts at the Cys113 SG-C bond and caps both sides with hydrogen — the
standard link-atom scheme:

    G_complex   the adduct, SG bonded to the ligand, Cys113 built as CYX
    G_receptor  apo protein, Cys113 a normal CYS with its HG back
    G_ligand    the ligand alone, attachment atom capped by H
                (this is EXACTLY the adduct form that was docked)

    dG = G_complex - G_receptor - G_ligand

dG therefore carries a constant chemical term — one S-C bond formed, one S-H and
one C-H broken. That constant is identical for every candidate in a warhead
class and cancels when they are compared, which is the only comparison D0020
licenses anyway. Across classes it does NOT cancel, and dG must not be used to
rank one warhead chemistry against another.

G_receptor is the same apo protein for every candidate, so it is computed once
and cached. It cancels entirely in any ranking; it is kept in dG only so the
number has its conventional meaning.

WHAT THIS IS NOT. Single-structure MM-GBSA on minimised geometries, not
ensemble-averaged over an MD trajectory. That is the standard cheap protocol for
ranking and it is what the plan budgeted. It gives no entropy, and its absolute
values are not free energies of binding.

FOUR THINGS THE RECEPTOR PREP HAD TO FIX, each found by tleap refusing to build:

1. **Histidine protonation was being thrown away.** `reduce` had assigned HIS27
   as doubly-protonated (HIP, +1), HIS59 as HID and HIS64/157 as HIE. tleap
   defaults every HIS to HIE, silently discarding a formal charge and two
   tautomer assignments. Protonation states are read back from which hydrogens
   `reduce` actually placed.
2. **Cys113 must be CYX, not CYS.** CYX is the residue that exists precisely for
   an SG bonded to something else; as CYS, tleap adds an HG and the sulfur ends
   up with an extra bond.
3. **tleap renumbers residues.** The 6VAJ construct starts at 6 and is missing
   40-47, and tleap closes that gap — so Cys113 becomes residue 105 and a bond
   command written against 113 silently targets the wrong residue. Residues are
   renumbered contiguously here and the bond target is VERIFIED by name before
   the topology is saved.
4. **igb=8 requires mbondi3 radii**, and the same radii must be set on all three
   legs. Mixing radii between legs makes the difference meaningless while every
   individual calculation still looks fine.

THE JUNCTION HAS NO PARAMETERS. ff19SB's `S` meets GAFF2's `c3` at the covalent
bond and no force field covers that pair. `data/params/cys_gaff2_junction_2.frcmod`
supplies the missing bond, angle and dihedral terms. Every one is the GAFF2
parameter for the same geometry with GAFF2's thioether sulfur `ss` in place of
the protein's `S` — an analogue rather than a convenience, since a Cys SG that
has bonded a ligand IS a thioether sulfur.

v1 of that file hand-picked analogues from parm19 and covered only sp3 carbon.
Six of the nine warhead classes attach through sp2 or aromatic carbon, so 18 of
27 MM-GBSA builds failed on missing `S-c2`, `S-cc` and `S-ca` terms. Deriving
every term from one source closed the gap and removed the inconsistency of
mixing parameter sets between classes.

This remains an approximation at the one bond the whole calculation is about,
and it is the largest modelling assumption in this module.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent

AMBER_ENV = Path("/data/lab_vm/envs/dwi_amber_md")
RECEPTOR_PDB = Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb")
JUNCTION_FRCMOD = _REPO_ROOT / "data" / "params" / "cys_gaff2_junction_2.frcmod"

COVALENT_RESNUM = "113"          # Cys113, in the ORIGINAL PDB numbering
COVALENT_RESNAME = "CYS"
LIGAND_RESNAME = "LIG"

PROTEIN_FF = "leaprc.protein.ff19SB"
LIGAND_FF = "leaprc.gaff2"
PB_RADII = "mbondi3"             # required by igb=8
IGB = 8

ENERGY_TERMS = ("BOND", "ANGLE", "DIHED", "VDWAALS", "EEL", "EGB", "ESURF",
                "1-4VDW", "1-4EEL")


class MMGBSAError(RuntimeError):
    """A leg could not be built, minimised, or scored."""


@dataclass
class LegEnergies:
    """Parsed sander output for one leg."""

    terms: dict[str, float] = field(default_factory=dict)

    @property
    def total(self) -> float:
        return sum(self.terms.get(k, 0.0) for k in ENERGY_TERMS)


def _amber(tool: str) -> str:
    p = AMBER_ENV / "bin" / tool
    if not p.is_file():
        raise MMGBSAError(f"{tool} not found at {p}; is the dwi_amber_md env built?")
    return str(p)


def _run(cmd: list[str], cwd: Path, log_name: str, timeout: int = 3600) -> None:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    (cwd / log_name).write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise MMGBSAError(
            f"{cmd[0].split('/')[-1]} failed ({proc.returncode}); see {cwd/log_name}\n"
            + proc.stderr[-600:])


# --------------------------------------------------------------------------
# receptor preparation
# --------------------------------------------------------------------------

def _histidine_state(atom_names: set[str]) -> str:
    """HID / HIE / HIP from which hydrogens `reduce` actually placed."""
    hs = tuple(sorted(n for n in atom_names if n in ("HD1", "HE2")))
    return {("HD1",): "HID", ("HE2",): "HIE", ("HD1", "HE2"): "HIP"}.get(hs, "HIE")


def prepare_receptor(workdir: Path, receptor_pdb: Path | None = None
                     ) -> tuple[Path, Path, int, int]:
    """Write the CYX (covalent) and CYS (apo) receptor PDBs, renumbered.

    Returns
    -------
    (cyx_pdb, cys_pdb, cyx_index, n_residues)
        ``cyx_index`` is the 1-based residue index tleap will use, which is only
        trustworthy because the residues were renumbered contiguously first.
    """
    from collections import OrderedDict, defaultdict

    src = Path(receptor_pdb or RECEPTOR_PDB)
    if not src.is_file():
        raise MMGBSAError(f"prepared receptor not found: {src}")
    lines = [l for l in src.read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]

    order: OrderedDict = OrderedDict()
    atoms: dict = defaultdict(set)
    for l in lines:
        key = (l[21], l[22:27])
        order.setdefault(key, l[17:20])
        atoms[key].add(l[12:16].strip())

    renumber = {k: i + 1 for i, k in enumerate(order)}
    his = {k: _histidine_state(atoms[k]) for k in order}

    targets = [k for k in order
               if k[1].strip() == COVALENT_RESNUM and order[k] == COVALENT_RESNAME]
    if len(targets) != 1:
        raise MMGBSAError(
            f"expected exactly one {COVALENT_RESNAME}{COVALENT_RESNUM} in {src.name}, "
            f"found {len(targets)}")
    target = targets[0]

    def emit(covalent: bool) -> str:
        out = []
        for l in lines:
            key = (l[21], l[22:27])
            resname, atom = l[17:20], l[12:16].strip()
            if resname in ("HIS", "HIE", "HID", "HIP"):
                l = l[:17] + his[key] + l[20:]
            elif covalent and key == target:
                if atom == "HG":
                    continue          # CYX has no HG; SG bonds the ligand
                l = l[:17] + "CYX" + l[20:]
            out.append(f"{l[:22]}{renumber[key]:4d} {l[27:]}")
        return "\n".join(out) + "\nTER\nEND\n"

    cyx, cys = workdir / "receptor_cyx.pdb", workdir / "receptor_cys.pdb"
    cyx.write_text(emit(True), encoding="utf-8")
    cys.write_text(emit(False), encoding="utf-8")

    states = {v for k, v in his.items() if order[k] in ("HIS", "HIE", "HID", "HIP")}
    log.info("receptor prepared: %d residues, Cys%s -> index %d, histidine states %s",
             len(order), COVALENT_RESNUM, renumber[target], sorted(states) or "none")
    return cyx, cys, renumber[target], len(order)


# --------------------------------------------------------------------------
# ligand parameterisation
# --------------------------------------------------------------------------

def parameterize_ligand(pose_sdf: Path, workdir: Path, attachment_smarts: str,
                        net_charge: int = 0) -> tuple[Path, Path, str, str, float]:
    """AM1-BCC/GAFF2 parameters for the docked adduct pose.

    Returns
    -------
    (mol2, frcmod, attachment_atom_name, cap_hydrogen_name, corrected_charge)
        The cap hydrogen is removed in tleap and its charge folded onto the
        attachment atom, so the ligand stays neutral. Without that the complex
        picks up a spurious fractional charge that GB electrostatics will
        happily and silently use.
    """
    mols = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m]
    if not mols:
        raise MMGBSAError(f"no readable pose in {pose_sdf}")
    mol = Chem.AddHs(mols[0], addCoords=True)

    patt = Chem.MolFromSmarts(attachment_smarts)
    matches = mol.GetSubstructMatches(patt) if patt else []
    if not matches:
        raise MMGBSAError(
            f"attachment SMARTS {attachment_smarts!r} does not match the pose in "
            f"{pose_sdf.name}; the pose and the warhead class disagree")
    attach_idx = min(m[0] for m in matches)

    lig_sdf = workdir / "ligand_h.sdf"
    w = Chem.SDWriter(str(lig_sdf))
    w.write(mol)
    w.close()

    mol2 = workdir / "ligand.mol2"
    _run([_amber("antechamber"), "-i", str(lig_sdf), "-fi", "sdf",
          "-o", str(mol2), "-fo", "mol2", "-c", "bcc", "-nc", str(net_charge),
          "-at", "gaff2", "-rn", LIGAND_RESNAME, "-s", "0", "-pf", "y"],
         workdir, "antechamber.log", timeout=1800)
    if not mol2.is_file():
        raise MMGBSAError("antechamber produced no mol2 (AM1-BCC likely failed)")

    frcmod = workdir / "ligand.frcmod"
    _run([_amber("parmchk2"), "-i", str(mol2), "-f", "mol2",
          "-o", str(frcmod), "-s", "gaff2"], workdir, "parmchk2.log")

    name, cap, corrected = _cap_hydrogen(mol2, attach_idx)
    log.info("ligand parameterised: attachment %s, cap H %s, corrected charge %+.4f",
             name, cap, corrected)
    return mol2, frcmod, name, cap, corrected


def _cap_hydrogen(mol2: Path, attach_idx: int) -> tuple[str, str, float]:
    """Name the attachment atom and one of its hydrogens; fold that H's charge in.

    antechamber preserves input atom order, so ``attach_idx`` (0-based, from the
    SDF) is mol2 atom ``attach_idx + 1``. That is asserted rather than assumed:
    the attachment atom must not be a hydrogen.
    """
    lines = mol2.read_text().splitlines()

    def section(tag: str) -> list[str]:
        i = lines.index(tag)
        out = []
        for l in lines[i + 1:]:
            if l.startswith("@<TRIPOS>"):
                break
            if l.strip():
                out.append(l)
        return out

    atoms = []
    for l in section("@<TRIPOS>ATOM"):
        p = l.split()
        atoms.append({"id": int(p[0]), "name": p[1], "type": p[5], "charge": float(p[8])})
    bonds = [(int(p[1]), int(p[2]))
             for p in (l.split() for l in section("@<TRIPOS>BOND"))]

    if attach_idx >= len(atoms):
        raise MMGBSAError(f"attachment index {attach_idx} outside the mol2")
    att = atoms[attach_idx]
    if att["type"].lower().startswith("h"):
        raise MMGBSAError(
            f"attachment atom {att['name']} is a hydrogen — atom order between "
            "the SDF and the mol2 does not correspond")

    hydrogens = [a for a in atoms if a["type"].lower().startswith("h")
                 and any({att["id"], a["id"]} == {x, y} for x, y in bonds)]
    if not hydrogens:
        raise MMGBSAError(
            f"attachment atom {att['name']} has no hydrogen to replace; gnina "
            "formed the bond by displacing one, so the pose and the topology "
            "disagree about what is bonded there")
    cap = hydrogens[0]
    return att["name"], cap["name"], att["charge"] + cap["charge"]


# --------------------------------------------------------------------------
# topology construction
# --------------------------------------------------------------------------

def build_topologies(workdir: Path, mol2: Path, frcmod: Path, cyx_pdb: Path,
                     cys_pdb: Path, cyx_index: int, ligand_index: int,
                     attach_name: str, cap_name: str, corrected_charge: float
                     ) -> dict[str, tuple[Path, Path]]:
    """Build complex, receptor and ligand topologies in one tleap pass."""
    if not JUNCTION_FRCMOD.is_file():
        raise MMGBSAError(f"junction parameters missing: {JUNCTION_FRCMOD}")
    shutil.copy(JUNCTION_FRCMOD, workdir / "junction.frcmod")

    script = f"""source {PROTEIN_FF}
source {LIGAND_FF}
loadamberparams {frcmod.name}
loadamberparams junction.frcmod
set default PBRadii {PB_RADII}

LIGF = loadmol2 {mol2.name}
saveamberparm LIGF ligand.prmtop ligand.inpcrd

rec = loadpdb {cys_pdb.name}
saveamberparm rec receptor.prmtop receptor.inpcrd

LIGX = loadmol2 {mol2.name}
remove LIGX LIGX.1.{cap_name}
set LIGX.1.{attach_name} charge {corrected_charge:.6f}
prot = loadpdb {cyx_pdb.name}
cpx = combine {{prot LIGX}}
bond cpx.{cyx_index}.SG cpx.{ligand_index}.{attach_name}
saveamberparm cpx complex.prmtop complex.inpcrd
quit
"""
    (workdir / "build.leap").write_text(script, encoding="utf-8")
    _run([_amber("tleap"), "-f", "build.leap"], workdir, "tleap.log")

    legs = {}
    for leg in ("complex", "receptor", "ligand"):
        top, crd = workdir / f"{leg}.prmtop", workdir / f"{leg}.inpcrd"
        if not top.is_file() or top.stat().st_size == 0:
            raise MMGBSAError(
                f"tleap produced no usable {leg} topology; see {workdir/'tleap.log'}")
        legs[leg] = (top, crd)

    missing = [l for l in (workdir / "tleap.log").read_text().splitlines()
               if "Could not find" in l]
    if missing:
        raise MMGBSAError(
            f"{len(missing)} missing force-field parameter(s), first: "
            f"{missing[0].strip()}. Add it to {JUNCTION_FRCMOD.name} with its "
            "analogue cited, rather than letting tleap substitute silently.")
    verify_complex(legs["complex"][0], cyx_index, attach_name)
    return legs


def verify_complex(prmtop: Path, cyx_index: int, attach_name: str) -> dict:
    """Assert the covalent bond exists and the ligand kept its charge.

    tleap reports Errors = 0 for a complex where the `bond` command targeted the
    wrong residue, so the topology is checked rather than trusted.
    """
    try:
        import parmed as pmd
    except ImportError as exc:
        raise MMGBSAError("parmed is needed to verify the topology") from exc

    p = pmd.load_file(str(prmtop))
    cyx = [r for r in p.residues if r.name == "CYX"]
    if not cyx:
        raise MMGBSAError("no CYX residue in the complex — Cys113 was not converted")
    sgs = [a for r in cyx for a in r.atoms if a.name == "SG"]
    bonded_to_ligand = []
    for sg in sgs:
        for b in sg.bonds:
            other = b.atom1 if b.atom2 is sg else b.atom2
            if other.residue.name == LIGAND_RESNAME:
                bonded_to_ligand.append((sg.residue.idx + 1, other.name))
    if not bonded_to_ligand:
        raise MMGBSAError(
            "no bond between any CYX SG and the ligand: tleap's `bond` command "
            "did not do what it was asked, and it does not report that as an error")
    if len(bonded_to_ligand) > 1:
        raise MMGBSAError(f"more than one SG-ligand bond: {bonded_to_ligand}")

    res_idx, atom_name = bonded_to_ligand[0]
    if atom_name != attach_name:
        raise MMGBSAError(
            f"SG is bonded to ligand atom {atom_name}, not the intended "
            f"{attach_name}")

    lig = [r for r in p.residues if r.name == LIGAND_RESNAME][0]
    att = [a for a in lig.atoms if a.name == attach_name][0]
    if len(att.bonds) != 4:
        log.warning("attachment atom %s has %d bonds, expected 4",
                    attach_name, len(att.bonds))
    info = {"cyx_residue_index": res_idx, "attachment_atom": atom_name,
            "ligand_atoms": len(lig.atoms),
            "ligand_charge": round(sum(a.charge for a in lig.atoms), 4),
            "complex_charge": round(sum(a.charge for a in p.atoms), 4),
            "attachment_bonds": len(att.bonds)}
    log.info("complex verified: %s", info)
    return info


# --------------------------------------------------------------------------
# minimisation + energies
# --------------------------------------------------------------------------

MIN_INPUT = f"""minimise, implicit solvent (GBn2)
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  ntb=0, igb={IGB}, gbsa=1, cut=999.0,
  ntpr=200, ntr=0,
 /
"""


def minimize_and_score(workdir: Path, leg: str) -> LegEnergies:
    """Minimise one leg and parse its final GB energy terms."""
    (workdir / "min.in").write_text(MIN_INPUT, encoding="utf-8")
    _run([_amber("sander"), "-O", "-i", "min.in",
          "-p", f"{leg}.prmtop", "-c", f"{leg}.inpcrd",
          "-o", f"{leg}.min.out", "-r", f"{leg}.min.rst"],
         workdir, f"sander_{leg}.log", timeout=7200)

    out = (workdir / f"{leg}.min.out").read_text()
    if "FINAL RESULTS" not in out:
        raise MMGBSAError(f"{leg}: minimisation did not converge to FINAL RESULTS")
    block = out.split("FINAL RESULTS")[-1]
    terms: dict[str, float] = {}
    for key, val in re.findall(r"([A-Z0-9\-]+)\s*=\s*(-?\d+\.\d+)", block):
        terms.setdefault(key, float(val))
    e = LegEnergies(terms=terms)
    log.info("%s minimised: G = %.2f kcal/mol", leg, e.total)
    return e


def delta_g(legs: dict[str, LegEnergies]) -> dict:
    """dG = G_complex - G_receptor - G_ligand, with the per-term breakdown."""
    for name in ("complex", "receptor", "ligand"):
        if name not in legs:
            raise MMGBSAError(f"missing leg {name!r}")
    dg = legs["complex"].total - legs["receptor"].total - legs["ligand"].total
    per_term = {
        t: round(legs["complex"].terms.get(t, 0.0)
                 - legs["receptor"].terms.get(t, 0.0)
                 - legs["ligand"].terms.get(t, 0.0), 4)
        for t in ENERGY_TERMS
    }
    return {
        "dG_kcal": round(dg, 3),
        "G_complex": round(legs["complex"].total, 3),
        "G_receptor": round(legs["receptor"].total, 3),
        "G_ligand": round(legs["ligand"].total, 3),
        "per_term": per_term,
        "scheme": "link-atom 3-leg, cut at Cys113 SG-C, both sides H-capped",
        "comparable": "within warhead class only (D0020); the constant bond "
                      "term does not cancel across classes",
        "igb": IGB,
        "pb_radii": PB_RADII,
        "ensemble_averaged": False,
    }
