#!/usr/bin/env python3
"""
Purpose: for a candidate that has a crystal control, decide whether the screen
         FAILED TO GENERATE the crystallographic pose or generated it and then
         FAILED TO SELECT it.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-12
Input: the run's persisted pose cloud (`<run.topic>_allposes/<ident>.sdf`) and
       the crystal control in the production receptor's frame
       (`crystal_control_poses/<xtal>.sdf`, written by `crystal_controls.py`)
Output: 00_outputs/blacksmith/crystal_pose_audit/<ident>_poses_<N>.csv (per pose)
        + <ident>_summary_<N>.json

THE QUESTION, AND WHY ONE MEASUREMENT ANSWERS IT (#64). Sulfopin -- the known
covalent Pin1 inhibitor, carried through the screen as an ordinary candidate --
ranks in the bottom 7% of its own warhead class. Two readings fit that equally
well and they need OPPOSITE fixes:

  (a) GENERATION. Docking never produces the crystallographic pose, so no
      selection rule could have picked it and a better scorer changes nothing.
  (b) SELECTION. The pose is in the cloud and the pipeline exports something
      else, in which case replacing the pose generator is effort spent on the
      wrong stage.

Scoring the WHOLE cloud against the crystal separates them: if the best pose in
the cloud is within 2 A, the pose was generated and (a) is refuted.

TWO ENDPOINTS, BECAUSE ONE OF THEM IS THE WRONG QUESTION (D0062). Whole-molecule
RMSD is the field's convention and #64 asks for it, but D0062 measured that it
correlates only rho = +0.433 with error in placing the REACTIVE region, and that
52.7% of poses place the reactive region correctly while failing the 2 A RMSD
test. Our criterion asks whether the warhead can reach Cys113, not whether the
tail is where the crystallographer found it. Both are reported, per pose, and
neither is called "the" recovery number.

THE CRYSTAL LIGAND IS THE ADDUCT AND THE DOCKED LIGAND IS THE FREE FORM.
6VAJ deposits QT7 with its chlorine ALREADY DISPLACED -- 16 modelled atoms, C10
bonded to Cys113 SG at 1.78 A (the LINK record; D0001) -- while `FORMUL` still
lists CL because that field comes from the chemical component definition, which
describes the free ligand. Taking the formula rather than the modelled atoms is
catalogue #21 exactly. So the comparison is over the 16-atom COMMON SUBSTRUCTURE
and the docked chlorine has no counterpart to be scored against. That is stated
in the output rather than left for a reader to infer from an atom count.

THE FRAME IS CHECKED BEFORE ANYTHING IS MEASURED. #64 was blocked because the
only sulfopin SDF the author found was the ORIGIN-FRAMED docking input --
centroid 13.5 A from Cys113 SG -- which gives a best RMSD of 12.07 A and 0 of
455 poses within 2.5 A. That reads exactly like a catastrophic docking failure
and is a coordinate-frame mismatch. `_check_frame` re-derives the distance from
the crystal's reactive carbon to the receptor's own SG and REFUSES to continue
if it is not a bond length, so this cannot be reported as a docking result.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from shared import covalent_protocol as cp    # noqa: E402
from shared import nac_criterion as nac       # noqa: E402
from shared import outputs as sout            # noqa: E402
from shared import receptors as srec          # noqa: E402
from shared import target_config as tc        # noqa: E402

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("crystal-pose-audit")

B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
XTAL_DIR = B / "crystal_control_poses"

#: The redocking success criterion (Astex/PoseBusters), and the stricter line
#: the pose-prediction literature also reports. Same constants as
#: `redock_04_rmsd.py`; they are the convention, not a choice made here.
SUCCESS_A = 2.0
TIGHT_A = 1.0

#: The crystal's reactive carbon is COVALENTLY BONDED to SG, so its distance to
#: SG is a bond length. Anything beyond this is a different coordinate frame,
#: not a worse pose. Generous by design -- the failure it must catch is 13.5 A,
#: and the superposition itself already costs ~0.2 A (1.78 A in 6VAJ's own
#: frame becomes 1.98 A after the fit onto 3IKD).
FRAME_MAX_A = 3.0

#: candidate ident -> the crystal control that IS that molecule. Never inferred
#: from a name match: `t4_sulfopin` and `xtal_6VAJ` share no substring, and a
#: pairing that guessed would be free to guess wrong.
PAIRS = {"t4_sulfopin": "xtal_6VAJ"}


class FrameError(RuntimeError):
    """The crystal and the poses are not in the same coordinate frame."""


# --------------------------------------------------------------------------
# graph matching
#
# Lifted deliberately from `redock_04_rmsd.symmetric_rmsd` rather than imported:
# that module resolves its own paths and topic at import time. The docstrings
# there explain why the graph is flattened; the behaviour must not diverge, and
# `tests/test_crystal_pose_audit.py` asserts the two agree on a shared case.
# --------------------------------------------------------------------------
def _heavy_graph(m: Chem.Mol) -> Chem.Mol:
    """A bond-order- and charge-agnostic heavy-atom copy, for matching only."""
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(0)
    out = rw.GetMol()
    Chem.FastFindRings(out)
    return out


def _coords(m: Chem.Mol) -> np.ndarray:
    c = m.GetConformer()
    return np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])


def _automorphisms(graph: Chem.Mol) -> list[tuple[int, ...]]:
    """The template's symmetry group -- the swaps X-ray density cannot resolve.

    For QT7 there are 12: the three neopentyl methyls (3! = 6) times the two
    sulfone oxygens. A symmetry-naive RMSD would score a correct pose with the
    methyls permuted as an error of a bond length.
    """
    autos = graph.GetSubstructMatches(graph, uniquify=False, maxMatches=50000)
    return list(autos) if autos else [tuple(range(graph.GetNumAtoms()))]


def reactive_atoms(xtal: Chem.Mol) -> tuple[int, list[int]]:
    """(reactive carbon, the acetamide region) as indices into the crystal.

    DERIVED FROM CONNECTIVITY, NOT FROM ATOM NAMES OR POSITIONS. `C10` is the
    name the depositor happened to give the carbon that bonds to SG; another
    entry calls the equivalent atom `C19`, `C14`, `C24`, `C12`, `C3`
    (`crystal_controls_2.csv`). Selecting it by name would work on 6VAJ and
    silently pick the wrong atom on the next structure -- catalogue disguise #1.

    The acetamide adduct is the project's own `adduct_attachment_smarts` for
    chloroacetamide, `[CH3][CX3](=O)[NX3]`: reactive C, carbonyl C, carbonyl O,
    amide N. It is matched here by connectivity for the same reason -- a
    hydrogen-count SMARTS cannot match the docked FREE ligand, where that carbon
    carries a chlorine and two hydrogens rather than three.
    """
    ns = [a.GetIdx() for a in xtal.GetAtoms() if a.GetSymbol() == "N"]
    if len(ns) != 1:
        raise ValueError(f"expected exactly one nitrogen, found {len(ns)}")
    n = ns[0]
    carbonyl = [
        a.GetIdx() for a in xtal.GetAtomWithIdx(n).GetNeighbors()
        if a.GetSymbol() == "C"
        and any(x.GetSymbol() == "O" and x.GetDegree() == 1
                for x in a.GetNeighbors())
    ]
    if len(carbonyl) != 1:
        raise ValueError(f"expected one amide carbonyl, found {len(carbonyl)}")
    c = carbonyl[0]
    oxy = [a.GetIdx() for a in xtal.GetAtomWithIdx(c).GetNeighbors()
           if a.GetSymbol() == "O"][0]
    reactive = [
        a.GetIdx() for a in xtal.GetAtomWithIdx(c).GetNeighbors()
        if a.GetSymbol() == "C" and a.GetDegree() == 1
    ]
    if len(reactive) != 1:
        raise ValueError(
            f"expected one terminal carbon on the carbonyl, found "
            f"{len(reactive)} -- this ligand is not a simple acetamide adduct "
            f"and the reactive atom must be established for it explicitly"
        )
    return reactive[0], [reactive[0], c, oxy, n]


def anchor_sg() -> np.ndarray:
    """The anchor atom's coordinate in the PRODUCTION receptor.

    Residue and atom come from `config/target.yaml`, not from a literal: the
    tool is target-agnostic and `Cys113`/`SG` is a fact about Pin1. The residue
    number is the CRYSTAL numbering -- the MD system renumbers from 1 and calls
    this residue 63 (`PIN1_OFFSET = 50`), which is the confusion that once drew
    a glutamate labelled as the target cysteine.
    """
    rec = srec.resolve_3ikd_ian(noligand=True)
    anchor = str(tc.get("target.anchor"))            # "Cys113"
    atom = str(tc.get("target.anchor_atom"))         # "SG"
    resnum = "".join(ch for ch in anchor if ch.isdigit())
    resname = "".join(ch for ch in anchor if ch.isalpha()).upper()[:3]
    for line in rec.read_text().splitlines():
        if (line.startswith("ATOM") and line[12:16].strip() == atom
                and line[17:20].strip().upper() == resname
                and line[22:26].strip() == resnum):
            return np.array([float(line[30:38]), float(line[38:46]),
                             float(line[46:54])])
    raise FrameError(f"{anchor} {atom} not found in {rec}")


def _warhead_pattern(warhead_class: str) -> tuple[str, str]:
    """(`reactive_atom_smarts`, mechanism) for a class, from the library.

    `reactive_atom_smarts`, NOT `adduct_attachment_smarts` (D0022). These poses
    are the UNREACTED ligand -- non-covalent docking of the free form -- so the
    pattern that names the leaving group is the one that matches, and the
    leaving group is precisely what the SN2 angle is measured against. The
    adduct pattern would match too, on a different set of atoms, and give a
    plausible angle about the wrong geometry: catalogue #21.
    """
    from shared import warhead_library as wl

    df = wl.load()
    row = df[df["class_id"].astype(str) == warhead_class]
    if len(row) != 1:
        raise ValueError(
            f"{len(row)} rows for warhead class {warhead_class!r} in the "
            f"library; expected exactly one")
    r = row.iloc[0]
    return str(r["reactive_atom_smarts"]), str(r["mechanism"])


def _rx(mol: Chem.Mol, patt: Chem.Mol) -> tuple[int, ...]:
    """The reactive-atom match on one persisted pose.

    Sulfopin has one CH2-Cl and this is unambiguous. A molecule with several
    reactive centres is REFUSED rather than resolved by taking the first match:
    the screen breaks that tie by asking the docking output which atom it
    retyped (`nac_screen.rebuild_and_match`), and that information is not in the
    persisted SDF. Guessing here would measure the approach to a different
    carbon than the screen did, and the angle would look entirely reasonable.
    """
    ms = mol.GetSubstructMatches(patt)
    if len(ms) != 1:
        raise ValueError(
            f"{len(ms)} reactive-atom matches on this pose; this audit can "
            f"only recompute the criterion where the centre is unambiguous")
    return ms[0]


def _check_frame(xtal: Chem.Mol, reactive_idx: int, sg: np.ndarray) -> float:
    """Distance from the crystal's reactive carbon to the receptor's Cys113 SG.

    THE GUARD THAT #64 NEEDED. It can fail: point it at the origin-framed
    `m2_covalent_smoke/sulfopin.sdf` and it raises at 13.5 A instead of
    reporting 0 of 455 poses recovered.
    """
    anchor = str(tc.get("target.anchor"))
    atom = str(tc.get("target.anchor_atom"))
    d = float(np.linalg.norm(_coords(xtal)[reactive_idx] - sg))
    if d > FRAME_MAX_A:
        raise FrameError(
            f"the crystal ligand's reactive carbon is {d:.2f} A from "
            f"{anchor} {atom}, and it is COVALENTLY BONDED to it -- so this "
            f"structure is not in the receptor's frame. Measuring poses "
            f"against it would report a frame mismatch as a docking failure. "
            f"Use the crystal control from crystal_controls.py, which "
            f"superposes onto the production receptor and records the fit."
        )
    return d


def _join_screen_rows(poses: pd.DataFrame, allposes: pd.DataFrame,
                      ident: str) -> pd.DataFrame:
    """Attach the screen's own `distance`/`angle`/`viable` to each cloud pose.

    THE JOIN IS BY ORDER, WHICH IS THE DISGUISE THIS PROJECT BREAKS ON, so what
    guarantees the order is written down here and checked in `_check_join`.

    `nac_screen_v2` writes the cloud as
    `argsort(labels, kind="stable")` filtered to the real modes -- i.e. grouped
    by mode, and within a mode in conformer order, which is `pose_idx` order.
    That is the only correspondence available: the SDF carries `pose_rank` and
    `energy_rank`, and BOTH are the write-order counter 1..N stamped by
    `write_sdf`, not the pose's rank by energy. `energy_rank` in a persisted SDF
    is therefore a name that describes a quantity the value is not.

    So: the screen's rows for this molecule, noise mode dropped, sorted by
    (mode, pose_idx). Anything else -- including sorting by the SDF's own
    `energy_rank` -- gives a scrambled pairing that still produces a full table
    of plausible numbers.
    """
    if allposes.empty:
        raise ValueError("rank_v2 returned no per-pose rows")
    t = allposes[allposes["ident"].astype(str) == ident]
    t = t[t["mode"] != -1].sort_values(["mode", "pose_idx"]).reset_index(drop=True)
    if len(t) != len(poses):
        raise ValueError(
            f"{len(poses)} poses in the cloud but {len(t)} scored rows for "
            f"{ident}: the SDF and the screen's table describe different runs")
    out = poses.sort_values(["mode", "pose_rank"]).reset_index(drop=True)
    if not (out["mode"].astype(int).values == t["mode"].astype(int).values).all():
        raise ValueError("mode sequence differs between the cloud and the table")
    for c in ("pose_idx", "energy", "energy_rank", "distance", "angle",
              "viable", "in_range"):
        out[c] = t[c].values
    return out


def _check_join(poses: pd.DataFrame, agg: pd.DataFrame) -> dict:
    """Two independent tests that the order-based join paired the right rows.

    1. The screen's per-mode counts must fall out of the joined table. This
       catches a pairing that crossed mode boundaries.
    2. The geometry must agree pose by pose. Every pose is measured here against
       the RIGID receptor's SG and by the screen against the FLEXIBLE one, so
       the two distances differ by however far that sidechain moved -- a small
       offset, not an arbitrary one. Under a scrambled pairing the correlation
       collapses, because it would be comparing one pose's distance with
       another's. This is the test that a count-based check cannot do.
    """
    got = poses.groupby("mode").agg(n_poses_mode=("mode", "size"),
                                    n_in_range=("in_range", "sum"),
                                    n_viable=("viable", "sum"))
    want = agg.set_index("mode")[["n_poses_mode", "n_in_range", "n_viable"]]
    j = want.join(got, how="outer", lsuffix="_screen", rsuffix="_joined")
    counts_ok = all(
        j[f"{c}_screen"].fillna(-1).astype(int).equals(
            j[f"{c}_joined"].fillna(-1).astype(int))
        for c in ("n_poses_mode", "n_in_range", "n_viable"))
    r = float(np.corrcoef(poses["distance"], poses["distance_rigid_sg"])[0, 1])
    dev = float(np.abs(poses["distance"] - poses["distance_rigid_sg"]).median())
    if not counts_ok:
        raise ValueError(
            f"the join does not reproduce the screen's per-mode counts:\n"
            f"{j.to_string()}")
    if r < 0.9:
        raise ValueError(
            f"joined poses correlate at r = {r:.3f} between the flexible-SG "
            f"distance the screen measured and the rigid-SG distance measured "
            f"here. Below 0.9 the rows are not paired correctly.")
    return {"join_counts_reproduced": counts_ok,
            "join_distance_r": round(r, 4),
            "median_flex_vs_rigid_sg_shift_a": round(dev, 3)}


# --------------------------------------------------------------------------
def audit(ident: str, xtal_ident: str) -> tuple[pd.DataFrame, dict]:
    """Score every persisted pose of `ident` against its crystal control."""
    topic = str(tc.get("run.topic"))
    cloud = B / f"{topic}_allposes" / f"{ident}.sdf"
    reps = B / f"{topic}_poses" / f"{ident}.sdf"
    xpath = XTAL_DIR / f"{xtal_ident}.sdf"
    for p in (cloud, xpath):
        if not p.is_file():
            raise FileNotFoundError(
                f"{p} does not exist. The pose cloud must come from the SAME "
                f"run that produced the scores (#44) -- re-docking gives a "
                f"different cloud, because docking is stochastic with no seed."
            )

    # The run's OWN record of what class it screened this molecule as, through
    # rank_v2's loader so the molecule-level supersession rule is the same one
    # the ranking used. This matters here: two runs of `t4_sulfopin` sit in the
    # aggregates, with `mode` 0-4 in both, and the older one is not the cloud
    # being audited.
    import rank_v2 as rv                                          # noqa: PLC0415

    # `rank_v2.SRC` is a module global that its own `main()` rebinds from
    # `--topic`; imported rather than run, it still points at the nac_v2
    # default. Bound here from `run.topic` for the same reason D0080 exists --
    # the default is a plausible directory that would load a different run's
    # aggregates and report them under this one's name.
    rv.SRC = B / topic
    if not rv.SRC.is_dir():
        raise FileNotFoundError(f"no such topic directory: {rv.SRC}")
    agg, allposes = rv.load_v2()
    mine = agg[agg["parent_ident"].astype(str) == ident]
    if mine.empty:
        raise ValueError(f"{ident} has no aggregate row under topic {topic}")
    classes = sorted(set(mine["warhead_class"].astype(str)))
    if len(classes) != 1:
        raise ValueError(f"{ident} is recorded under {classes}; expected one")
    warhead_class = classes[0]

    xtal = next(iter(Chem.SDMolSupplier(str(xpath), removeHs=False,
                                        sanitize=False)))
    if xtal is None:
        raise ValueError(f"{xpath} unreadable")
    r_idx, region = reactive_atoms(xtal)
    sg = anchor_sg()
    sg_dist = _check_frame(xtal, r_idx, sg)

    # THE CRITERION IS TAKEN FROM THE SCREEN, NOT RECOMPUTED -- and the first
    # version of this script got that wrong in a way worth recording.
    #
    # `nac_criterion.measure` needs the position of Cys113 SG, and the screen
    # takes it from the FLEXIBLE sidechain in the docked model being measured:
    # the sulfur moves during docking, so a distance to the rigid receptor's SG
    # measures the approach to where the sulfur STARTED
    # (`nac_screen.sg_position`). The persisted cloud holds the ligand only, so
    # that per-pose sulfur position is NOT in the archive and the criterion
    # cannot be reproduced from it. Recomputing against the rigid SG looked
    # right -- mode sizes matched exactly, 2/202/41/210 -- and moved 6 poses of
    # mode 3 across the viability bar (10 -> 16). Small, plausible, and wrong.
    #
    # So the screen's own per-pose rows are joined in, and `nac.measure` is
    # still run as an INDEPENDENT CHECK on that join rather than as the answer.
    lib_smarts, mechanism = _warhead_pattern(warhead_class)
    patt = Chem.MolFromSmarts(lib_smarts)

    qt = _heavy_graph(xtal)
    autos = _automorphisms(qt)
    xyz_x = _coords(xtal)

    def score(mol: Chem.Mol) -> tuple[float, float, float] | None:
        heavy = Chem.RemoveHs(mol, sanitize=False)
        qd = _heavy_graph(heavy)
        m = qd.GetSubstructMatch(qt)
        if not m:
            return None
        xd = _coords(heavy)[list(m)]
        # The minimum is over the TEMPLATE's automorphisms, applied to the
        # docked coordinates -- the crystal stays put. Both molecules are
        # already in the receptor's frame, so nothing is superposed: aligning
        # them would erase the displacement being measured (redock_04_rmsd).
        best, best_a = None, None
        for a in autos:
            v = float(np.sqrt(((xd[list(a)] - xyz_x) ** 2).sum(1).mean()))
            if best is None or v < best:
                best, best_a = v, a
        xd_b = xd[list(best_a)]
        reg = float(np.sqrt(((xd_b[region] - xyz_x[region]) ** 2).sum(1).mean()))
        rea = float(np.linalg.norm(xd_b[r_idx] - xyz_x[r_idx]))
        return best, reg, rea

    rows, unmatched = [], 0
    for mol in Chem.SDMolSupplier(str(cloud), removeHs=False, sanitize=False):
        if mol is None:
            unmatched += 1
            continue
        s = score(mol)
        if s is None:
            unmatched += 1
            continue
        props = mol.GetPropsAsDict()
        m = nac.measure(mechanism,
                        mol.GetConformer().GetPositions()[list(_rx(mol, patt))],
                        sg)
        rows.append({
            "ident": ident,
            "mode": props.get("mode"),
            "pose_rank": props.get("pose_rank"),
            "energy_rank": props.get("energy_rank"),
            "rmsd_a": s[0],
            "region_rmsd_a": s[1],
            "reactive_atom_err_a": s[2],
            "distance_rigid_sg": m.distance,
            "angle_rigid_sg": m.angle,
        })
    poses = pd.DataFrame(rows)
    if poses.empty:
        raise ValueError(f"no pose in {cloud} matched the crystal graph")

    poses = _join_screen_rows(poses, allposes, ident)
    check = _check_join(poses, mine)

    # The exported representatives -- what the pipeline actually carries
    # forward. Scored through the SAME function, so a difference between them
    # and the cloud is a difference in the poses and not in the metric.
    rep_rows = []
    if reps.is_file():
        for mol in Chem.SDMolSupplier(str(reps), removeHs=False,
                                      sanitize=False):
            if mol is None:
                continue
            s = score(mol)
            if s is None:
                continue
            props = mol.GetPropsAsDict()
            rep_rows.append({
                "mode": props.get("mode"),
                "mode_label": props.get("mode_label"),
                "rmsd_a": s[0],
                "region_rmsd_a": s[1],
                "reactive_atom_err_a": s[2],
            })
    reps_df = pd.DataFrame(rep_rows)

    near = poses[poses["rmsd_a"] <= SUCCESS_A]
    null = nac.isotropic_null(mechanism)
    best_i = int(poses["rmsd_a"].idxmin())
    summary = {
        "ident": ident,
        "crystal": xtal_ident,
        "run_topic": topic,
        "warhead_class": warhead_class,
        "mechanism": mechanism,
        "reactive_atom_smarts": lib_smarts,
        **check,
        "crystal_atoms": int(xtal.GetNumAtoms()),
        "crystal_has_halogen": any(a.GetSymbol() in ("Cl", "Br", "F", "I")
                                   for a in xtal.GetAtoms()),
        "compared_atoms": int(qt.GetNumAtoms()),
        "automorphisms": len(autos),
        "crystal_reactive_to_sg_a": round(sg_dist, 3),
        "n_poses": int(len(poses)),
        "n_unmatched": int(unmatched),
        "best_rmsd_a": round(float(poses["rmsd_a"].min()), 3),
        "best_rmsd_mode": poses.loc[best_i, "mode"],
        "median_rmsd_a": round(float(poses["rmsd_a"].median()), 3),
        f"n_within_{SUCCESS_A}a": int((poses["rmsd_a"] <= SUCCESS_A).sum()),
        f"frac_within_{SUCCESS_A}a": round(
            float((poses["rmsd_a"] <= SUCCESS_A).mean()), 4),
        f"n_within_{TIGHT_A}a": int((poses["rmsd_a"] <= TIGHT_A).sum()),
        "best_region_rmsd_a": round(float(poses["region_rmsd_a"].min()), 3),
        "n_region_within_1a": int((poses["region_rmsd_a"] <= 1.0).sum()),
        "best_reactive_err_a": round(
            float(poses["reactive_atom_err_a"].min()), 3),
        "n_representatives": int(len(reps_df)),
        "best_representative_rmsd_a": (
            round(float(reps_df["rmsd_a"].min()), 3) if len(reps_df) else None),

        # ---- what the CRITERION does to the poses that ARE the crystal pose --
        # The crystal ligand is the ADDUCT: its reactive carbon is bonded to SG,
        # so it sits BELOW the 2.8 A floor of the near-attack window and cannot
        # be attack-ready in its own frame (D0077). A pose that reproduces it
        # closely inherits that. These three numbers say whether that is a
        # theoretical worry or what is actually happening.
        "n_near_crystal": int(len(near)),
        "n_near_crystal_in_range": int(near["in_range"].sum()),
        "n_near_crystal_viable": int(near["viable"].sum()),
        "median_distance_near_crystal_a": (
            round(float(near["distance"].median()), 3) if len(near) else None),
        "median_distance_all_a": round(float(poses["distance"].median()), 3),
        "median_rmsd_of_viable_a": (
            round(float(poses.loc[poses["viable"], "rmsd_a"].median()), 3)
            if poses["viable"].any() else None),
        "median_rmsd_of_nonviable_a": (
            round(float(poses.loc[~poses["viable"], "rmsd_a"].median()), 3)
            if (~poses["viable"]).any() else None),
        "n_viable": int(poses["viable"].sum()),

        # ---- the criterion's score AS A FUNCTION OF BEING RIGHT --------------
        # `enrichment` = viable_fraction / isotropic_null is what ranks a mode
        # and what the sweep's floor is set on. Recomputing it over successively
        # tighter shells around the crystal pose asks the one question the class
        # rank cannot: what does this criterion award a pose that IS the answer?
        # If enrichment stays below the floor even there, no improvement to
        # docking or to pose selection can put this molecule through the gate,
        # and that is a fact about the criterion rather than about the molecule.
        "isotropic_null": round(null, 5),
        "enrichment_by_shell": {
            k: {"n": int(len(v)), "n_viable": int(v["viable"].sum()),
                "viable_fraction": round(float(v["viable"].mean()), 4),
                "enrichment": round(float(v["viable"].mean() / null), 3)}
            for k, v in (("all_poses", poses),
                         ("within_2.0a_of_crystal", poses[poses.rmsd_a <= 2.0]),
                         ("within_1.5a_of_crystal", poses[poses.rmsd_a <= 1.5]))
            if len(v)
        },
        # Of the poses that reproduce the crystal, what does the criterion
        # reject them ON? Distance and angle are separate bars and they have
        # different remedies.
        "near_crystal_fail_distance": int((~near["in_range"]).sum()),
        "near_crystal_fail_angle": int(
            (near["in_range"] & (near["angle"] < nac.SN2_ANGLE_MIN)).sum())
        if mechanism == "sn2_displacement" else None,
        "near_crystal_median_angle": round(float(near["angle"].median()), 1),
        "all_median_angle": round(float(poses["angle"].median()), 1),
    }
    if len(reps_df) and summary[f"n_within_{SUCCESS_A}a"]:
        # The number the generation-vs-selection question turns on: how much of
        # the cloud is closer to the crystal than the best thing the pipeline
        # chose to carry.
        summary["cloud_better_than_best_rep_frac"] = round(float(
            (poses["rmsd_a"] < reps_df["rmsd_a"].min()).mean()), 4)
    return poses, reps_df, summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate", default="t4_sulfopin")
    ap.add_argument("--crystal", default=None,
                    help="crystal control ident; defaults to the PAIRS entry")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    xt = a.crystal or PAIRS.get(a.candidate)
    if xt is None:
        raise SystemExit(
            f"no crystal control is recorded for {a.candidate}. Add it to "
            f"PAIRS only if that molecule IS the deposited ligand -- a "
            f"chemically similar one is not a recovery reference."
        )

    poses, reps, summary = audit(a.candidate, xt)
    t = sout.Topic("blacksmith", "crystal_pose_audit")
    poses.to_csv(t.write(f"{a.candidate}_poses", ".csv"), index=False)
    t.write(f"{a.candidate}_summary", ".json").write_text(
        json.dumps(summary, indent=2, default=str))

    print("\n" + "=" * 76)
    print(f"  CRYSTAL POSE AUDIT — {a.candidate} vs {xt} ({summary['run_topic']})")
    print("=" * 76)
    print(f"\n  crystal ligand: {summary['crystal_atoms']} atoms, "
          f"halogen present: {summary['crystal_has_halogen']} — compared over "
          f"{summary['compared_atoms']} common heavy atoms, "
          f"{summary['automorphisms']} automorphisms")
    print(f"  frame check: reactive C to Cys113 SG = "
          f"{summary['crystal_reactive_to_sg_a']} A (a bond length)\n")
    print(f"  join to the screen's rows: counts reproduced="
          f"{summary['join_counts_reproduced']}, "
          f"r={summary['join_distance_r']}, "
          f"flexible-vs-rigid SG shift="
          f"{summary['median_flex_vs_rigid_sg_shift_a']} A\n")
    for k in ("n_poses", "best_rmsd_a", "median_rmsd_a",
              f"n_within_{SUCCESS_A}a", f"frac_within_{SUCCESS_A}a",
              f"n_within_{TIGHT_A}a", "best_region_rmsd_a",
              "n_region_within_1a", "best_reactive_err_a",
              "best_representative_rmsd_a", "cloud_better_than_best_rep_frac",
              "n_viable", "n_near_crystal", "n_near_crystal_in_range",
              "n_near_crystal_viable", "median_distance_near_crystal_a",
              "median_distance_all_a", "median_rmsd_of_viable_a",
              "median_rmsd_of_nonviable_a"):
        if k in summary:
            print(f"  {k:34s} {summary[k]}")
    floor = tc.sweep_budget_floor()
    print("\n  enrichment = viable_fraction / isotropic_null "
          f"({summary['isotropic_null']}), by distance from the crystal pose:")
    for k, v in summary["enrichment_by_shell"].items():
        mark = "" if v["enrichment"] >= floor else "   < budget_floor"
        print(f"    {k:26s} n={v['n']:4d}  viable={v['n_viable']:3d}  "
              f"enrichment={v['enrichment']}{mark}")
    print(f"    the sweep's budget_floor is {floor}")
    if len(reps):
        print("\n  exported representatives (what the pipeline carries):")
        print(reps.to_string(index=False))
    print()


if __name__ == "__main__":
    main()
