"""
Purpose: the full covalent workup for ONE named candidate — adduct, covalent docking and MM-GBSA — on each receptor separately.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: a candidate_id present in a current T_3/T_4 frame
Output: 00_outputs/blacksmith/covalent_workup/covalent_workup_<tag>_<N>.json

WHY THIS EXISTS. Every covalent path in this repo is hard-wired to **6VAJ**:
`covalent_protocol.dock` defaults to `6VAJ_prepared.pdbqt`, `mmgbsa.RECEPTOR_PDB`
is `6VAJ_prepared.pdb`, and `config/receptor.yaml` still pins `pdb_id: 6VAJ`.
D0059 replaced the receptor with **3IKD** and invalidated every 6VAJ
measurement, but only the NON-covalent side was re-pointed. So the covalent
numbers this project can produce today describe a receptor it no longer uses,
and there is no CLI that docks a single named molecule at all.

Both entry points take `receptor_pdbqt=` / `receptor_pdb=` overrides that
nothing in the repo passes. This script passes them, and runs each receptor as
its own labelled leg.

THE TWO LEGS ARE NOT COMPARABLE AND ARE NOT COMBINED. 6VAJ and 3IKD place the
pocket 48.6 A apart (D0059); their boxes differ in centre and in size (20 A vs
26 A); and gnina's affinity is not receptor-transferable. Reporting a mean of
the two, or quietly preferring whichever is better, would manufacture a number
neither receptor produced. Each leg carries its own `receptor` field and they
are written side by side so the reader does the comparing.

WHAT THE SCORES ARE WORTH. `cnn_affinity`/`cnn_score` are advisory: gnina emits
"CNN scoring not yet calibrated for covalent docking" on every covalent run and
D0011 demoted them accordingly — the flag is carried through as
`cnn_uncalibrated_for_covalent` rather than left in a log. `affinity_kcal` is
the affinity-best mode over all modes, not row 0 (D0047). None of these
validate a ranking: the gate verdict on this stratum is UNDERPOWERED.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import pandas as pd                                  # noqa: E402

from shared import covalent_adduct as cad            # noqa: E402
from shared import covalent_protocol as cp           # noqa: E402
from shared import mmgbsa as mg                      # noqa: E402
from shared import outputs as sout                   # noqa: E402
from shared import warhead_library as wl             # noqa: E402

log = logging.getLogger("covalent-workup")

OUT = sout.Topic("blacksmith", "covalent_workup")
WORK = Path("/data/lab_vm/modifiable/inhibition/covalent_workup")

# Each receptor is a (pdbqt for gnina, box, pdb for tleap) triple. Keeping them
# in one table is the point: a leg cannot be assembled from one receptor's
# docking grid and another's topology, which is exactly how a 48.6 A error
# would enter without raising anything.
RECEPTORS = {
    "6VAJ": {
        "pdbqt": Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdbqt"),
        "box": Path("/data/lab_vm/immutable/inhibition/receptor/box.json"),
        "pdb": Path("/data/lab_vm/immutable/inhibition/receptor/6VAJ_prepared.pdb"),
        "status": "SUPERSEDED by D0059; kept only for continuity with the "
                  "affinity_kcal/cnn_* already on the frame",
        "sg_xyz": (-12.530, -35.870, 8.190),
    },
    "3IKD": {
        "pdbqt": Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                      "receptor_3ikd/3IKD_prepared_1.pdbqt"),
        "box": Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                    "receptor_3ikd/box_3IKD_1.json"),
        "pdb": Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/"
                    "3IKD_noligand.pdb"),
        "status": "current receptor (D0059), chemist-prepared",
        # Cys113 SG, from this receptor's own prepared coordinates. The two
        # receptors' pockets are 48.6 A apart, so this single triple
        # distinguishes them beyond any doubt.
        "sg_xyz": (13.385, 3.989, -2.040),
    },
}

# How far a built receptor's SG may sit from the expected triple before the leg
# is refused. Not zero: prepare_receptor renumbers, re-protonates and may add a
# missing heavy atom, so bitwise equality is the wrong test. 0.5 A is far
# tighter than the 48.6 A that separates the two receptors and far looser than
# any preparation step moves a sulfur.
SG_TOLERANCE_A = 0.5


class WorkupError(RuntimeError):
    """A leg failed. Named so a failure cannot be read as a result."""


def assert_receptor_identity(cyx_pdb: Path, rec_name: str, cyx_index: int) -> dict:
    """The built receptor IS the named one — checked, not inherited.

    `mmgbsa.prepare_receptor(workdir, receptor_pdb=None)` defaults to
    **6VAJ**. Every covalent path in this repo takes that default. A leg that
    forgets the override builds a topology against a pocket 48.6 A away from
    the one the non-covalent legs used, minimises happily, and returns a dG
    with nothing wrong on its face.

    Passing the override is not the same as knowing it took effect, so the
    receptor's own Cys113 SG is read back out of the structure tleap will be
    handed and compared with the coordinate that identifies it. This turns the
    receptor from an inherited default into a checked fact recorded in the
    output.

    THE RESIDUE IS FOUND BY `cyx_index`, NOT BY THE NUMBER 113.
    `prepare_receptor` renumbers residues contiguously so that the index it
    returns is the one tleap will use, which means Cys113 is residue 63 in the
    prepared 3IKD and 100 in the prepared 6VAJ. Searching for residue "113"
    finds nothing in either, and an identity check that cannot locate its own
    atom is worse than none — it fails on correct input and would have to be
    disabled to make anything run.
    """
    want = RECEPTORS[rec_name]["sg_xyz"]
    got = None
    for line in cyx_pdb.read_text().splitlines():
        if line[:6] in ("ATOM  ", "HETATM") and line[12:16].strip() == "SG" \
                and line[17:20].strip() in ("CYS", "CYX") \
                and line[22:26].strip() == str(cyx_index):
            got = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            break
    if got is None:
        raise WorkupError(
            f"no CYS/CYX SG at renumbered residue {cyx_index} in {cyx_pdb}")

    d = sum((a - b) ** 2 for a, b in zip(got, want)) ** 0.5
    other = {n: sum((a - b) ** 2 for a, b in
                    zip(got, RECEPTORS[n]["sg_xyz"])) ** 0.5
             for n in RECEPTORS if n != rec_name}
    if d > SG_TOLERANCE_A:
        nearest = min(other, key=other.get) if other else "?"
        raise WorkupError(
            f"leg claims receptor {rec_name} but its Cys113 SG is at "
            f"{got}, {d:.2f} A from {rec_name}'s {want}"
            + (f" and {other[nearest]:.2f} A from {nearest}'s — this is "
               f"{nearest}, not {rec_name}" if other else ""))
    log.info("  receptor identity: Cys113 SG %.2f A from %s (%s)",
             d, rec_name, ", ".join(f"{n} {v:.1f} A" for n, v in other.items()))
    return {"receptor_verified": rec_name, "sg_xyz_built": list(got),
            "sg_xyz_expected": list(want), "sg_offset_a": round(d, 3),
            "sg_offset_to_other_receptors_a": {k: round(v, 1)
                                               for k, v in other.items()}}


def candidate_row(cid: str) -> pd.Series:
    """The candidate's row from the newest frame that holds it."""
    import export_nac_poses as enp
    for subdir, stem in enp.FRAMES.values():
        f = enp.latest_frame(subdir, stem)
        df = pd.read_parquet(f)
        hit = df[df["candidate_id"].astype(str) == cid]
        if len(hit):
            r = hit.iloc[0].copy()
            r["_frame"] = f.name
            return r
    raise WorkupError(f"{cid} is in no current frame")


def stereo_centres(mol) -> list[tuple[int, str]]:
    """(atom index, CIP label) for every stereocentre, '?' when unassigned."""
    from rdkit import Chem
    m = Chem.Mol(mol)
    if m.GetNumConformers():
        Chem.AssignStereochemistryFrom3D(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                     useLegacyImplementation=False)


def adduct_from_pose(pose_sdf: Path, adduct_smiles: str, leaving_group: str,
                     dest: Path, pose_rank: int = 1) -> tuple[Path, dict]:
    """The adduct built FROM the ranked pose, not re-embedded from its SMILES.

    THE CANDIDATE'S WARHEAD STEREOCENTRE IS UNDEFINED, AND THAT MAKES
    RE-EMBEDDING NON-DETERMINISTIC ACROSS LEGS.
    `FindMolChiralCenters(..., includeUnassigned=True)` on the frame's own
    SMILES returns `[(5, 'R'), (13, '?')]`: the sulfolane carbon is specified
    and the 4,5-dihydroisoxazole ring carbon is not. Whatever runs, runs on ONE
    arbitrary configuration, chosen by whichever embedding happened first.

    They did not choose the same one. Measured on this candidate:

        NAC pose rank 1 (what MD residence and BPMD used)      (5R, 13R)
        adduct embedded from SMILES with randomSeed=42          (5R, 13S)

    So the covalent leg was scoring the OPPOSITE diastereomer from the
    non-covalent legs, and every number would have completed, agreed with
    itself, and been quietly incomparable. Nothing raises: both are valid
    molecules and the SMILES does not claim otherwise.

    Deleting the leaving group from the POSE and reading stereochemistry back
    off the 3D coordinates keeps the configuration the ranking actually saw.
    Verified for this warhead: the CIP label at the ring carbon is unchanged by
    the loss of bromide (R before, R after), so the transfer is a fact about
    this structure and not an assumption about CIP priorities.

    THIS FIXES CONSISTENCY, NOT THE UNDERLYING AMBIGUITY. The other enantiomer
    remains unsimulated, and BDHI enantiomers are reported to differ ~50-fold
    in potency on TG2. The configuration is returned so the report can name it
    rather than imply the molecule has only one.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mols = [m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m]
    ranked = [m for m in mols if m.HasProp("pose_rank")
              and int(m.GetProp("pose_rank")) == pose_rank] or mols[:1]
    pose = Chem.RemoveHs(ranked[0])
    Chem.AssignStereochemistryFrom3D(pose)
    free_centres = stereo_centres(pose)

    lg = Chem.MolFromSmiles(leaving_group) if leaving_group else None
    if lg is None or lg.GetNumAtoms() != 1:
        raise WorkupError(
            f"leaving group {leaving_group!r} is not a single atom; building the "
            "adduct from the pose is only defined for a single-atom leaving "
            "group. Fall back to embedding and record the configuration.")
    sym = lg.GetAtomWithIdx(0).GetSymbol()

    rw = Chem.RWMol(pose)
    victims = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == sym]
    if len(victims) != 1:
        raise WorkupError(f"{len(victims)} {sym} atoms in the pose; the leaving "
                          "group is not identified")
    rw.RemoveAtom(victims[0])
    built = rw.GetMol()
    Chem.SanitizeMol(built)
    Chem.AssignStereochemistryFrom3D(built)

    # The molecule built from the pose must BE the adduct the warhead library
    # defines. Constitution is compared with stereochemistry stripped, because
    # that is exactly the part the two disagree about and the part being
    # transferred.
    want = Chem.MolFromSmiles(adduct_smiles)
    flat = lambda m: Chem.MolToSmiles(Chem.MolFromSmiles(  # noqa: E731
        Chem.MolToSmiles(m, isomericSmiles=False)))
    if flat(built) != flat(want):
        raise WorkupError(
            f"the pose minus {sym} is {flat(built)} but the warhead library's "
            f"adduct is {flat(want)}; they are not the same molecule")

    mol = Chem.AddHs(built, addCoords=True)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise WorkupError("could not embed the adduct")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    Chem.AssignStereochemistryFrom3D(mol)

    final = stereo_centres(Chem.RemoveHs(mol))
    if any(lab == "?" for _, lab in final):
        raise WorkupError(f"adduct still has an unassigned centre: {final}")

    w = Chem.SDWriter(str(dest))
    w.write(mol)
    w.close()
    return dest, {
        "built_from": f"{pose_sdf.name} pose_rank {pose_rank}, minus {sym}",
        "pose_free_centres": [[i, l] for i, l in free_centres],
        "adduct_centres": [[i, l] for i, l in final],
        "unassigned_in_smiles": [i for i, l in
                                 stereo_centres(Chem.MolFromSmiles(adduct_smiles))
                                 if l == "?"],
        "caveat": "the opposite configuration at the warhead ring carbon is "
                  "UNSIMULATED; BDHI enantiomers differ ~50-fold on TG2",
    }


def junction_coverage(mol2: Path, attach_name: str) -> dict:
    """Do the junction parameters cover THIS warhead's attachment carbon?

    The junction frcmod is organised by GAFF2 atom type, and its own header
    says it covers "sp3/sp2/aromatic carbon". That is a claim about types, not
    about this molecule: the type antechamber actually assigned is the only
    thing that decides whether a term exists. Reading the assigned type out of
    the mol2 and listing the terms the bond will need turns "it should be
    covered" into something a reader can check — and, when it is not covered,
    names the missing term instead of leaving tleap to substitute.

    Returns the assigned type, the neighbour types, and the required terms with
    a present/absent verdict on each.
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

    atoms = {}
    for l in section("@<TRIPOS>ATOM"):
        p = l.split()
        atoms[int(p[0])] = {"name": p[1], "type": p[5]}
    by_name = {a["name"]: i for i, a in atoms.items()}
    if attach_name not in by_name:
        raise WorkupError(f"attachment atom {attach_name} not in {mol2.name}")
    aid = by_name[attach_name]
    at = atoms[aid]["type"]

    nbrs = []
    for l in section("@<TRIPOS>BOND"):
        p = l.split()
        x, y = int(p[1]), int(p[2])
        if x == aid:
            nbrs.append(atoms[y])
        elif y == aid:
            nbrs.append(atoms[x])

    # The cap hydrogen is removed and replaced by SG, so its own angle term is
    # not needed; every other neighbour's is.
    heavy = [n for n in nbrs if not n["type"].startswith("h")]
    required = {
        "BOND": [("S", at)],
        "ANGLE": [("2C", "S", at)] + [(n["type"], at, "S") for n in heavy],
        "DIHE": [("X", at, "S", "X")],
    }

    have = _frcmod_terms(mg.JUNCTION_FRCMOD)
    checks = {}
    for section, terms in required.items():
        for t in terms:
            # AMBER terms are direction-agnostic: `S -c2` and `c2-S ` are the
            # same bond. Both orientations are tried rather than assuming the
            # file happens to be written the way this function spells it.
            ok = (t in have[section]) or (tuple(reversed(t)) in have[section])
            checks["%s %s" % (section, "-".join(t))] = ok

    return {
        "attachment_atom": attach_name,
        "attachment_gaff2_type": at,
        "neighbour_types": sorted(n["type"] for n in nbrs),
        "junction_frcmod": mg.JUNCTION_FRCMOD.name,
        "required_terms": checks,
        "all_present": all(checks.values()),
        "missing": [t for t, ok in checks.items() if not ok],
    }


# frcmod atom-type fields are FIXED WIDTH, two characters each, separated by
# '-': `S -c2`, `2C-S -c2`, `X -c2-S -X `. Splitting on whitespace would merge
# `S -c2` into one token and splitting on '-' alone would not strip the padding,
# so both are done and each field is stripped. A looser parser here would report
# a term as present because some other line happened to contain its letters,
# which is the failure mode this whole check exists to prevent.
_WIDTHS = {"BOND": 5, "ANGLE": 8, "DIHE": 11, "IMPROPER": 11}


def _frcmod_terms(path: Path) -> dict[str, set[tuple[str, ...]]]:
    """Every term key defined in an frcmod, by section."""
    out: dict[str, set[tuple[str, ...]]] = {k: set() for k in _WIDTHS}
    section = None
    for line in path.read_text().splitlines():
        head = line.split()[0] if line.split() else ""
        if head in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
            section = head
            continue
        if section not in _WIDTHS or not line.strip():
            continue
        key = line[:_WIDTHS[section]]
        if "-" not in key:
            continue
        out[section].add(tuple(p.strip() for p in key.split("-")))
    return out


def dock_leg(row: pd.Series, adduct, rec_name: str, wd: Path, gpu: int,
             pose_sdf: Path) -> dict:
    """Covalent docking on ONE named receptor."""
    r = RECEPTORS[rec_name]
    for k in ("pdbqt", "box", "pdb"):
        if not r[k].is_file():
            raise WorkupError(f"{rec_name}: missing {k} at {r[k]}")
    wd.mkdir(parents=True, exist_ok=True)
    lig, stereo = adduct_from_pose(pose_sdf, adduct.adduct_smiles,
                                   adduct.leaving_group_smiles, wd / "adduct.sdf")
    log.info("  adduct stereocentres %s (from the ranked pose, not re-embedded)",
             stereo["adduct_centres"])
    out = wd / "adduct_docked.sdf"
    import os
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
    t0 = time.time()
    res = cp.dock(lig, out, str(row["warhead_class"]),
                  receptor_pdbqt=r["pdbqt"], box=r["box"], timeout_s=1800)
    res.update({
        "receptor": rec_name,
        "receptor_status": r["status"],
        "receptor_pdbqt": str(r["pdbqt"]),
        "box": str(r["box"]),
        "wall_seconds": round(time.time() - t0, 1),
        "pose_path": str(out),
        "stereochemistry": stereo,
    })
    return res


def mmgbsa_leg(row: pd.Series, rec_name: str, pose: Path, wd: Path) -> dict:
    """Covalent MM-GBSA on ONE named receptor, from that receptor's own pose."""
    r = RECEPTORS[rec_name]
    wd.mkdir(parents=True, exist_ok=True)
    smarts = cad.adduct_attachment_smarts(str(row["warhead_class"]), library=wl.load())
    cyx, cys, cyx_idx, n_res = mg.prepare_receptor(wd, receptor_pdb=r["pdb"])
    identity = assert_receptor_identity(cyx, rec_name, cyx_idx)
    mol2, frcmod, att, cap, q = mg.parameterize_ligand(pose, wd, smarts, net_charge=0)

    # Checked BEFORE tleap runs, so a missing junction term is named rather
    # than discovered as "tleap produced no usable complex topology".
    cover = junction_coverage(mol2, att)
    log.info("  junction: attachment %s is GAFF2 %r; all terms present = %s",
             att, cover["attachment_gaff2_type"], cover["all_present"])
    if not cover["all_present"]:
        raise WorkupError(
            f"{mg.JUNCTION_FRCMOD.name} does not cover attachment type "
            f"{cover['attachment_gaff2_type']!r}: missing {cover['missing']}. "
            "Add the terms with their gaff2 analogue cited rather than letting "
            "tleap substitute silently.")

    legs = mg.build_topologies(wd, mol2, frcmod, cyx, cys, cyx_idx, n_res + 1,
                               att, cap, q)
    verified = mg.verify_complex(legs["complex"][0], cyx_idx, att)
    energies = {leg: mg.minimize_and_score(wd, leg)
                for leg in ("complex", "receptor", "ligand")}
    # Crystallographic waters ride along in the complex and receptor legs alike
    # (3IKD keeps 2; 6VAJ has none). Counted rather than assumed to cancel.
    n_wat = sum(1 for line in cyx.read_text().splitlines()
                if line[:6] in ("ATOM  ", "HETATM")
                and line[17:20].strip() in mg.SOLVENT_RESNAMES
                and line[12:16].strip() == "O")
    return {"receptor": rec_name, "receptor_status": r["status"],
            "receptor_pdb": str(r["pdb"]), "cyx_index": cyx_idx,
            "n_residues": n_res, "crystallographic_waters": n_wat,
            "attachment_atom": att, "cap_hydrogen": cap,
            "corrected_charge": q, "junction_coverage": cover,
            "receptor_identity": identity,
            "verified": verified, **mg.delta_g(energies)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--receptors", nargs="+", default=["3IKD", "6VAJ"],
                    choices=sorted(RECEPTORS))
    ap.add_argument("--gpu", type=int, default=7)
    ap.add_argument("--pose", default=None,
                    help="ranked pose SDF; default is the newest export")
    ap.add_argument("--skip-mmgbsa", action="store_true")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import export_nac_poses as enp
    pose_sdf = Path(args.pose) if args.pose else \
        enp.poses_dir() / f"{args.candidate}.sdf"
    if not pose_sdf.is_file():
        raise SystemExit(
            f"no ranked pose at {pose_sdf}. The adduct's warhead stereocentre "
            "is undefined, so without the pose there is nothing to inherit the "
            "configuration from and re-embedding would pick one arbitrarily.")

    row = candidate_row(args.candidate)
    lib = wl.load()
    adduct = cad.to_adduct_form(str(row["canonical_smiles"]),
                               str(row["warhead_class"]), library=lib)
    log.info("%s (%s) from %s", args.candidate, row["warhead_class"], row["_frame"])
    log.info("  free    %s", row["canonical_smiles"])
    log.info("  adduct  %s", adduct.adduct_smiles)

    report = {
        "candidate_id": args.candidate,
        "frame": row["_frame"],
        "warhead_class": str(row["warhead_class"]),
        "mechanism": str(row.get("warhead_mechanism", "")),
        "canonical_smiles": str(row["canonical_smiles"]),
        "adduct": adduct.as_dict(),
        "legs": {},
        "not_a_ranking": "covalent scores for one molecule; the gate verdict on "
                         "this stratum is UNDERPOWERED and rank_validated is False",
    }

    for rec in args.receptors:
        wd = WORK / args.candidate / rec
        leg: dict = {}
        log.info("[%s] covalent docking", rec)
        try:
            leg["docking"] = dock_leg(row, adduct, rec, wd, args.gpu, pose_sdf)
            log.info("  affinity %.3f kcal/mol (mode %s of %s), CNNaffinity %.3f",
                     leg["docking"]["affinity_kcal"],
                     leg["docking"]["selected_mode"], leg["docking"]["n_modes"],
                     leg["docking"]["cnn_affinity"])
        except Exception as exc:                          # noqa: BLE001
            leg["docking"] = {"receptor": rec,
                              "failed": f"{type(exc).__name__}: {str(exc)[:400]}"}
            log.warning("  DOCKING FAILED: %s", leg["docking"]["failed"])

        if not args.skip_mmgbsa and "failed" not in leg["docking"]:
            log.info("[%s] MM-GBSA", rec)
            try:
                leg["mmgbsa"] = mmgbsa_leg(row, rec, Path(leg["docking"]["pose_path"]),
                                           wd / "mmgbsa")
                log.info("  dG %.2f kcal/mol", leg["mmgbsa"]["dG_kcal"])
            except Exception as exc:                      # noqa: BLE001
                leg["mmgbsa"] = {"receptor": rec,
                                 "failed": f"{type(exc).__name__}: {str(exc)[:400]}"}
                log.warning("  MM-GBSA FAILED: %s", leg["mmgbsa"]["failed"])
        report["legs"][rec] = leg

    dest = OUT.write(f"covalent_workup_{args.tag or args.candidate}", ".json")
    dest.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(f"\n  -> {dest}")
    print(json.dumps(report["legs"], indent=2, default=str)[:4000])


if __name__ == "__main__":
    main()
