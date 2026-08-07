"""
Purpose: run binding-pose metadynamics — the execution driver `shared/bpmd.py` was written for.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: an exported near-attack pose (SDF) + the 3IKD receptor the pose was docked into
Output: 00_outputs/blacksmith/bpmd/bpmd_s<shard>_<N>.csv (per replicate) and,
        for --convergence, 00_outputs/blacksmith/bpmd_convergence/*.csv

WHY THIS SCRIPT DID NOT EXIST AND HAD TO. `shared/bpmd.py` is the scoring layer
(15 tests, all on synthetic trajectories with arithmetic answers) and
`shared/gromacs_explicit.py` carries the `plumed=` argument that applies a bias
to production. Neither of them runs anything. BPMD has therefore never been
measured, which is why D0068's correction -- "docking generates and filters,
binding-pose metadynamics ranks" -- is still an intention rather than a result.

THE FIRST DELIVERABLE IS A CONVERGENCE TEST, NOT A RANKING.

`--convergence` runs ONE known-active pose at the full protocol (10 replicas x
10 ns) and then reports what the stability score would have been at 2, 4, 6, 8
and 10 replicas and at 1, 2, 5 and 10 ns. The point is to find the cheapest
protocol that still gives the same answer, by MEASURING it.

That ordering is not fussiness. D0068 is exactly what happens when a number is
believed before the thing producing it is understood: the enrichment metric was
run across 5,769 candidates and only later found to depend on how hard the
search looked, invalidating the shortlist it produced. BPMD's equivalent knob is
replica count and trajectory length. If the score at 4 replicas x 2 ns is within
the replica-to-replica spread of the score at 10 x 10, then the cheap protocol
IS the protocol and the shortlist is affordable. If it is not, we have to know
that before spending a fortnight of GPU on the wrong one.

TWO THINGS THIS DRIVER REFUSES TO GUESS.

1. **The atom indices.** PLUMED counts from 1, Python from 0, and an off-by-one
   here biases a NEIGHBOURING atom -- the run completes, the COLVAR is
   well-formed, and the reported stability is plausible and wrong. Conversion
   goes through `bpmd.bpmd_atom_indices` and nothing else, and the warhead atom
   is located in the solvated system by INTERNAL GEOMETRY (see
   `_locate_ligand_atom`), never by assuming antechamber and tleap preserved the
   SDF's atom order. The identity assumption is checked rather than trusted, and
   it is reported per pose which way it resolved.

2. **Whether the pose is still a near-attack conformation in the MD system.**
   The docking that produced these poses ran with Cys113 as a FLEXIBLE residue,
   so the recorded `nac_distance_A` is the distance to a rotated sidechain. The
   MD system rebuilds the receptor exactly as the chemist prepared it (D0059),
   in its crystallographic rotamer, so the starting CV is a DIFFERENT number.
   Both are recorded side by side, and a pose that starts past the unbound
   threshold is refused rather than scored -- a stability for a pose that was
   never in the window is a number about nothing.

A FAILED REPLICATE IS A ROW, NOT A GAP. Every replicate writes a record with its
status. A silently missing replicate would shrink `n_replicas` and narrow the
reported spread, i.e. it would make the answer look MORE certain for having gone
wrong -- the precise inversion this project keeps catching.

WHICH ENV, BECAUSE ONLY ONE OF THEM WORKS. This driver spans two stacks that
live in different environments: meeko + AutoDock (to build a pose that has not
been exported) and parmed + the Amber/GROMACS binaries (to build and run the
system). `dwi_cheminf`'s meeko is broken (`import meeko` -> no module named
`gemmi`) and `dwi_amber_md` has no meeko at all. **`~/.micromamba/envs/dwi_reactive`
is the one env that carries both**, and it is already `nac_screen.RX_ENV`. The
GROMACS and Amber executables are invoked by absolute path, so the running
interpreter does not need them.

FAIR USE. One GPU, and not 0 or 4 (other people's work). Launch long runs under
`nice -n 19` in tmux:

    tmux new -s bpmd -d '~/.micromamba/envs/dwi_reactive/bin/python \\
        scripts/bpmd_run.py --convergence --gpu 2 2>&1 | tee bpmd.log'
"""

from __future__ import annotations

import argparse
import glob
import itertools
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

from shared import bpmd                            # noqa: E402
from shared import gromacs_explicit as gx          # noqa: E402
from shared import mmgbsa as mg                    # noqa: E402
from shared import mmgbsa_noncovalent as mgn       # noqa: E402
from shared import outputs as sout                 # noqa: E402
import nac_screen as ns                            # noqa: E402
import nac_rank as nr                              # noqa: E402
import export_nac_poses as enp                     # noqa: E402

log = logging.getLogger("bpmd-run")

OUT = sout.Topic("blacksmith", "bpmd")
CONV = sout.Topic("blacksmith", "bpmd_convergence")

WORK = Path("/data/lab_vm/modifiable/inhibition/bpmd_3ikd")
# The CURRENT pose set. `export_nac_poses` versions these (nac_poses/vN) because
# the directory lives under append_only and a re-run must not overwrite, so this
# resolves the newest rather than naming a fixed path.
POSES = enp.poses_dir()

# The receptor as the chemist prepared it — the same file the docking that
# produced these poses used, so the pose and the protein are in one frame (D0059).
RECEPTOR_3IKD = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_prep/3IKD_noligand.pdb")

# The full protocol the convergence test measures DOWN from. Not asserted to be
# necessary — that is the question.
FULL_REPLICAS = 10
FULL_PRODUCTION_PS = 10_000.0

# What --convergence sweeps. Truncations are in ps and must be <= the production
# length; subsample sizes must be <= the replica count.
TRUNCATIONS_PS = (1_000.0, 2_000.0, 5_000.0, 10_000.0)
SUBSAMPLES = (2, 4, 6, 8, 10)
SUBSET_DRAWS = 25                  # per subsample size, when C(k, n) exceeds it
SUBSET_SEED = 0xB1AD5              # fixed; a verdict must not move on re-analysis

# Coordinates survive antechamber and tleap unchanged, so an atom's distances to
# its own molecule's other atoms are identical in the SDF and in the solvated
# system to file precision (gro is 0.001 nm = 0.01 A). These tolerances are about
# rounding, not about chemistry.
MAP_RMS_TOL_A = 0.02               # a genuine match is this good
# A SECOND atom this good is not a near miss, it is a real internal symmetry --
# no rigid transform can tell the two apart. Anything worse than MAP_RMS_TOL_A is
# distinguishable geometry and must NOT be swept into a tie: a threshold loose
# enough to catch "similar" atoms silently hands the tie-break rule two atoms
# that are simply different, and it picks one.
MAP_TIE_TOL_A = MAP_RMS_TOL_A

# What the ligand is NOT. Used to find it by exclusion, which is the only handle
# that survives tleap moving residues around (see `resolve_system`).
SOLVENT_NAMES = {"WAT", "HOH", "Na+", "Cl-", "NA", "CL", "K+", "IP", "IM"}
PROTEIN_NAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "CYM", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH", "ACE", "NME", "NHE",
}


class BPMDRunError(RuntimeError):
    """A stage of the BPMD chain failed. Named so a failure cannot pass as a result."""


# --------------------------------------------------------------------------
# candidates and poses
# --------------------------------------------------------------------------

def candidate_index(with_positives: bool = True) -> dict[str, ns.Candidate]:
    """Every ident BPMD can be asked for: generated candidates and crystal positives.

    The crystallographic positives are included because they are the only
    molecules on this project whose reactivity with Cys113 is a structural fact
    rather than a score, which is what a convergence test needs underneath it.
    """
    by_id = {c.ident: c for c in nr.load_candidates()}
    if with_positives:
        wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
        meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts)
                for r in wh.itertuples()}
        for c in ns.crystal_positives(meta, None):
            by_id.setdefault(c.ident, c)
    return by_id


def resolve_pose(cand: ns.Candidate, *, nrun: int, gpu: str,
                 allow_redock: bool, sdf: Path | None = None
                 ) -> tuple[Path, str]:
    """The SDF BPMD starts from — preferring the pose the score was computed from.

    `export_nac_poses` already wrote the shortlists' poses under the SAME
    protocol that scored them, and starting from those keeps the stability and
    the near-attack number describing one pose rather than two. When the ident
    has no exported pose — every crystallographic positive, since the export ran
    over T_3/T_4 shortlists only — it is re-docked through `export_nac_poses`'s
    own function rather than a second implementation of it.
    """
    if sdf is not None:
        if not Path(sdf).is_file():
            raise BPMDRunError(f"{cand.ident}: pose file {sdf} does not exist")
        return Path(sdf), "explicit"
    p = POSES / f"{cand.ident}.sdf"
    if p.is_file():
        return p, "exported"
    if not allow_redock:
        raise BPMDRunError(
            f"{cand.ident}: no exported pose at {p} and --no-redock was given")
    log.info("%s: no exported pose; re-docking under the scoring protocol",
             cand.ident)
    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    r = enp.export_one(cand, rec, nrun, gpu)
    if r.get("status") != "ok" or not r.get("nac_pose_path"):
        raise BPMDRunError(f"{cand.ident}: pose export failed: {r.get('status')}")
    return Path(r["nac_pose_path"]), "redocked"


def read_pose(sdf: Path, pose_rank: int = 1):
    """The pose as a complete molecule, with its recorded near-attack geometry.

    ONE POSE PER FILE USED TO BE THE EXPORT'S CONTRACT AND IS NOT ANY MORE.
    `export_nac_poses` now writes the TOP-N poses per molecule so that pose
    consensus is visible, and this reader refused every file it produced:

        t4_72f5671e89cb.sdf holds 10 molecules; one pose per file is the
        export's contract

    which is BPMD declining to run on the current pose set while reporting it
    as a malformed input. The contract moved; this reader did not.

    The pose is selected by its own `pose_rank` property, never by position.
    Taking `mols[0]` would work today — the export happens to write them in
    rank order — and would begin biasing a different pose the day that
    changes. This module is more exposed to that than most: BPMD on the wrong
    pose completes normally and reports a perfectly plausible stability.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]
    if not mols:
        raise BPMDRunError(f"no readable molecule in {sdf}")
    if len(mols) > 1:
        ranked = [m for m in mols if m.HasProp("pose_rank")
                  and int(m.GetProp("pose_rank")) == pose_rank]
        if len(ranked) != 1:
            raise BPMDRunError(
                f"{sdf.name} holds {len(mols)} poses and {len(ranked)} carry "
                f"pose_rank={pose_rank}; the rank does not identify one pose. "
                f"Ranks present: "
                f"{sorted(int(m.GetProp('pose_rank')) for m in mols if m.HasProp('pose_rank'))}")
        log.info("  pose_rank %d of %d poses in %s", pose_rank, len(mols), sdf.name)
        mols = ranked
    props = mols[0].GetPropsAsDict()
    # Hydrogens at the pose geometry, not by re-embedding: re-embedding would
    # discard the pose this whole run exists to test.
    mol = Chem.AddHs(mols[0], addCoords=True)
    return mol, props


def sg_in_pose_frame(cys_pdb: Path, cyx_index: int) -> np.ndarray:
    """Cys113's SG from the PREPARED RECEPTOR PDB — the pose's own coordinate frame.

    THIS IS NOT THE SG USED TO BUILD THE CV, and the difference is the point.
    `solvatebox` places the solute in a box, so the solvated system's coordinates
    are the pose's under a rigid motion. Comparing a pose coordinate against a
    solvated-system coordinate therefore measures the translation between two
    frames and calls it a distance.

    The prepared receptor PDB is written from 3IKD_noligand.pdb with coordinates
    untouched, so it is in the same frame as the docked pose, and it is what any
    pose-space geometry must be measured against.
    """
    for l in cys_pdb.read_text().splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[12:16].strip() == "SG" and l[17:20].strip() in ("CYS", "CYX") \
                and l[22:27].strip() == str(cyx_index):
            return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    raise BPMDRunError(f"no CYS{cyx_index} SG in {cys_pdb.name}")


def warhead_atom_in_pose(mol, cand: ns.Candidate, sg_xyz: np.ndarray) -> dict:
    """Which atom of the pose is the electrophile the sulfur attacks.

    `sg_xyz` must be in the POSE's frame (see `sg_in_pose_frame`).

    Index 0 of the reactive SMARTS match is the attacked atom for every mechanism
    in the library (see `nac_criterion.measure`). Several matches can share one
    attacked atom — a chloroazine's ipso carbon sits between two ring nitrogens —
    so matches are reduced to their distinct centres first.

    WHEN GENUINELY DISTINCT CENTRES REMAIN (a fumarate has two), the one nearest
    the sulfur is taken, and the margin to the runner-up is recorded. That is a
    chemical rule and not a default: the centre being presented to Cys113 is the
    one BPMD is being asked about, and for a symmetric acceptor the two are
    equivalent anyway. The margin travels with the result so a near-tie is
    visible rather than absorbed.
    """
    from rdkit import Chem

    patt = Chem.MolFromSmarts(cand.reactive_smarts)
    if patt is None:
        raise BPMDRunError(f"{cand.ident}: unparseable reactive SMARTS "
                           f"{cand.reactive_smarts!r}")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        raise BPMDRunError(f"{cand.ident}: reactive SMARTS does not match its own "
                           "exported pose — the pose and the frame disagree "
                           "about the molecule")
    centres = sorted({m[0] for m in matches})
    pos = mol.GetConformer(0).GetPositions()
    d = [float(np.linalg.norm(pos[i] - sg_xyz)) for i in centres]
    order = np.argsort(d)
    chosen = centres[int(order[0])]
    margin = float(d[int(order[1])] - d[int(order[0])]) if len(centres) > 1 else float("nan")
    return {"warhead_pose_idx": chosen, "n_reactive_centres": len(centres),
            "centre_margin_A": margin}


# --------------------------------------------------------------------------
# building the system
# --------------------------------------------------------------------------

def build_workdir(cand: ns.Candidate, mol, wd: Path) -> tuple[int, int]:
    """GAFF2 ligand parameters + the 3IKD receptor, in the layout `solvate` expects.

    Returns (cys113_residue_index_1based, receptor_residue_count), both over the
    PREPARED PDB's residues. The first is what `resolve_system` converts into a
    position in the built system; the second is reported only. It is explicitly
    NOT where the ligand lands — tleap moves waters during solvation, so the
    ligand's index is not `n_residues` and must be found by identity instead.
    """
    from rdkit import Chem

    wd.mkdir(parents=True, exist_ok=True)
    sdf = wd / "ligand_pose.sdf"
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()

    # THE FORMAL CHARGE OF THE MOLECULE IN THE FILE, not an assumed zero.
    # antechamber distributes AM1-BCC charges under the total it is handed, so a
    # carboxylate declared neutral silently smears a whole electron over the
    # ligand (D0058 is the same mistake from the descriptor side).
    charge = Chem.GetFormalCharge(mol)
    mgn.parameterize_ligand(sdf, wd, net_charge=charge)

    # receptor_cys.pdb — the FREE thiol. BPMD measures the REACTANT state: a
    # covalently tethered ligand cannot leave, so its stability is guaranteed by
    # construction and would measure nothing.
    _, _, cyx_index, n_res = mg.prepare_receptor(wd, receptor_pdb=RECEPTOR_3IKD)
    for f in ("ligand.mol2", "ligand.frcmod", "receptor_cys.pdb"):
        if not (wd / f).is_file():
            raise BPMDRunError(f"workdir incomplete: missing {f}")
    log.info("  parameterised (net charge %+d); Cys113 is residue %d of %d",
             charge, cyx_index, n_res)
    return cyx_index, n_res


def _heavy(struct_atoms) -> tuple[list[int], list[int]]:
    """(indices, atomic numbers) of the heavy atoms, in file order."""
    idx = [i for i, a in enumerate(struct_atoms) if a[1] > 1]
    return idx, [a[1] for a in struct_atoms if a[1] > 1]


def _locate_ligand_atom(ref_xyz: np.ndarray, ref_z: list[int], target: int,
                        sys_xyz: np.ndarray, sys_z: list[int],
                        sg_xyz: np.ndarray) -> dict:
    """Find `target` (an index into the pose) among the solvated system's ligand atoms.

    ORDER IS NOT ASSUMED. antechamber and tleap are free to renumber, and if they
    do, taking the pose's index straight across biases whatever atom happens to
    sit there — a mistake that produces a complete, plausible trajectory.

    What IS invariant is the molecule's internal geometry: neither program moves
    an atom, and `solvatebox` only places the solute in a box. So each atom is
    fingerprinted by the sorted vector of its distances to the ligand's other
    heavy atoms, which is unchanged by any rigid motion. A match must be good to
    within file precision (MAP_RMS_TOL_A) — nothing looser, because a threshold
    that admits merely SIMILAR atoms hands the tie-break below two atoms that are
    just different, and it then picks one.

    A true internal symmetry — the four corners of a ring, the two ends of a
    fumarate — leaves a genuine tie that no rigid transform can break. Those atoms
    are chemically equivalent, so the tie is resolved the way
    `warhead_atom_in_pose` resolves its own: nearest the sulfur. That rule is
    invariant under the frame change, so the two cannot disagree. Reported, never
    silent.
    """
    if len(ref_xyz) != len(sys_xyz):
        raise BPMDRunError(
            f"ligand heavy-atom count changed: pose has {len(ref_xyz)}, the "
            f"solvated system has {len(sys_xyz)}")

    ref_d = np.linalg.norm(ref_xyz[:, None, :] - ref_xyz[None, :, :], axis=-1)
    sys_d = np.linalg.norm(sys_xyz[:, None, :] - sys_xyz[None, :, :], axis=-1)

    # Fast path with a proof attached: if the two distance matrices agree
    # elementwise under the identity mapping and the elements line up, the order
    # was preserved and there is nothing to search.
    order_preserved = (ref_z == sys_z
                       and float(np.abs(ref_d - sys_d).max()) < MAP_RMS_TOL_A)
    if order_preserved:
        return {"sys_heavy_idx": target, "atom_map": "order preserved (verified)",
                "map_rms_A": float(np.abs(ref_d - sys_d).max()),
                "map_runner_up_A": float("nan"), "map_ties": 1}

    prof_ref = np.sort(ref_d[target])
    prof_sys = np.sort(sys_d, axis=1)
    rms = np.sqrt(((prof_sys - prof_ref) ** 2).mean(axis=1))
    rms = np.where(np.asarray(sys_z) == ref_z[target], rms, np.inf)
    ranked = np.argsort(rms)
    best = int(ranked[0])
    if not np.isfinite(rms[best]) or rms[best] > MAP_RMS_TOL_A:
        raise BPMDRunError(
            f"no atom of the solvated ligand matches the pose's warhead "
            f"(best RMS {rms[best]:.3f} A); the built system is not this pose")

    tied = [int(i) for i in ranked if rms[i] <= MAP_TIE_TOL_A]
    if len(tied) > 1:
        # Internal symmetry. The atoms are equivalent, so take the one facing the
        # sulfur — the SAME rule `warhead_atom_in_pose` used to choose among
        # reactive centres. Nearest-to-SG is invariant under the rigid motion
        # that separates the two frames, so the two choices cannot disagree.
        tied.sort(key=lambda i: float(np.linalg.norm(sys_xyz[i] - sg_xyz)))
    runner = float(rms[ranked[1]]) if len(ranked) > 1 else float("nan")
    return {"sys_heavy_idx": tied[0], "atom_map": "matched by internal geometry",
            "map_rms_A": float(rms[best]), "map_runner_up_A": runner,
            "map_ties": len(tied)}


def _pdb_residue_names(pdb: Path) -> list[str]:
    """Residue names of a prepared receptor PDB, in file order, one per residue."""
    from collections import OrderedDict
    o: OrderedDict = OrderedDict()
    for l in pdb.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            o.setdefault((l[21], l[22:27]), l[17:20].strip())
    return list(o.values())


def resolve_system(md_wd: Path, cys_pdb: Path, cyx_index: int) -> dict:
    """Cys113's SG and the ligand residue in the solvated system, located BY IDENTITY.

    TLEAP DOES NOT PRESERVE THE LAYOUT ITS OWN SCRIPT IMPLIES. `combine {rec LIG}`
    does put the ligand after the receptor, but `solvatebox` then moves every
    water — including the receptor's own crystallographic waters — to the end, so
    the ligand does NOT land at `n_residues`. Indexing there hit a WAT, which is
    the lucky version of this bug: 3IKD's receptor_cys.pdb carries 2
    crystallographic waters, and had they been ordered differently the same
    arithmetic would have landed on a real residue and produced a CV between two
    arbitrary atoms, in a run that would have completed and reported a number.

    So neither atom is taken positionally:

      * the LIGAND is the one residue that is neither protein nor solvent;
      * the SULFUR is Cys113's, found by its position in the PROTEIN-ONLY
        sequence, and that whole sequence must match the PDB tleap was handed
        before the index is used. Pin1 has two cysteines, so "it is a CYS" is not
        on its own evidence that it is the right one.
    """
    import parmed

    st = parmed.load_file(str(md_wd / "solv.prmtop"),
                          xyz=str(md_wd / "solv.inpcrd"))

    ligs = [r for r in st.residues
            if r.name not in PROTEIN_NAMES and r.name not in SOLVENT_NAMES]
    if len(ligs) != 1:
        raise BPMDRunError(
            f"expected exactly one non-protein, non-solvent residue (the "
            f"ligand); found {len(ligs)}: {sorted({r.name for r in ligs})}")
    lig = ligs[0]

    sys_prot = [r for r in st.residues if r.name in PROTEIN_NAMES]
    pdb_names = _pdb_residue_names(cys_pdb)
    pdb_prot = [n for n in pdb_names if n in PROTEIN_NAMES]
    if [r.name for r in sys_prot] != pdb_prot:
        raise BPMDRunError(
            f"the built system's protein sequence ({len(sys_prot)} residues) is "
            f"not the one in {cys_pdb.name} ({len(pdb_prot)}); tleap did not "
            "build the receptor that was prepared")

    # `cyx_index` is 1-based over ALL of the PDB's residues; convert it to a
    # position in the protein-only sequence, which is the ordering that survived.
    target_pos = sum(1 for n in pdb_names[:cyx_index - 1] if n in PROTEIN_NAMES)
    if pdb_names[cyx_index - 1] not in ("CYS", "CYX"):
        raise BPMDRunError(f"{cys_pdb.name} residue {cyx_index} is "
                           f"{pdb_names[cyx_index - 1]}, not a cysteine")
    cys = sys_prot[target_pos]
    if cys.name not in ("CYS", "CYX"):
        raise BPMDRunError(f"protein residue {target_pos + 1} is {cys.name}, "
                           "not the catalytic cysteine")
    sgs = [a for a in cys.atoms if a.name == "SG"]
    if len(sgs) != 1:
        raise BPMDRunError(f"{len(sgs)} atoms named SG in {cys.name}")

    sg = sgs[0]
    return {"struct": st, "sg": sg, "lig": lig,
            "sg_xyz": np.array([sg.xx, sg.xy, sg.xz]),
            "sg_label": f"{cys.name}{cyx_index}:SG", "lig_resname": lig.name}


def cv_atoms(md_wd: Path, sysres: dict, mol, warhead_pose_idx: int) -> dict:
    """The CV's two atoms, as PLUMED 1-based indices, plus the starting distance."""
    sg, lig = sysres["sg"], sysres["lig"]
    sg_xyz = sysres["sg_xyz"]

    sys_atoms = [(a.idx, a.atomic_number, (a.xx, a.xy, a.xz), a.name)
                 for a in lig.atoms]
    sys_heavy = [t for t in sys_atoms if t[1] > 1]
    sys_xyz = np.array([t[2] for t in sys_heavy], dtype=float)
    sys_z = [t[1] for t in sys_heavy]

    conf = mol.GetConformer(0).GetPositions()
    ref_heavy = [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() > 1]
    ref_xyz = conf[ref_heavy]
    ref_z = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in ref_heavy]
    if warhead_pose_idx not in ref_heavy:
        raise BPMDRunError("the warhead atom is a hydrogen — the reactive SMARTS "
                           "matched something that cannot be the electrophile")
    target = ref_heavy.index(warhead_pose_idx)

    found = _locate_ligand_atom(ref_xyz, ref_z, target, sys_xyz, sys_z, sg_xyz)
    warhead = sys_heavy[found["sys_heavy_idx"]]
    warhead_serial, sg_serial = warhead[0], sg.idx        # 0-based, parmed order

    # THE ONLY SUPPORTED CONVERSION. Bounds-checked against the system it will
    # actually be run on, and 1-based on the way out.
    gro = md_wd / "sys.gro"
    w1, s1 = bpmd.bpmd_atom_indices(gro, warhead_serial, sg_serial)

    d0_a = float(np.linalg.norm(np.array(warhead[2]) - sg_xyz))
    return {**found,
            "warhead_serial0": warhead_serial, "sg_serial0": sg_serial,
            # Carried so a LATER reader of a coordinate file can check that the
            # serial still names the atom it was derived for, instead of
            # trusting an integer across a file format boundary.
            "warhead_atom_name": warhead[3], "sg_atom_name": sg.name,
            "warhead_plumed": w1, "sg_plumed": s1,
            "sg_atom": sysres["sg_label"], "lig_resname": sysres["lig_resname"],
            "start_distance_A": round(d0_a, 3),
            "start_distance_nm": round(d0_a / 10.0, 4)}


# --------------------------------------------------------------------------
# running
# --------------------------------------------------------------------------

def run_replicate(cand_id: str, src: Path, md_wd: Path, plumed_txt: str,
                  replicate: int, *, gpu: int, threads: int,
                  production_ps: float,
                  reuse_equilibration: bool = False) -> tuple[Path, dict]:
    """One replicate, or the COLVAR a previous run already produced.

    `prod.gro` is written only when mdrun finishes (including a PLUMED COMMITTOR
    stop), so it — not the COLVAR, which exists from the first deposition — is
    what says a replicate completed. Reusing a half-written COLVAR would report a
    pose as stable because the job died before it could escape.
    """
    rep = md_wd / f"rep{replicate}"
    colvar, done = rep / "COLVAR", rep / "prod.gro"
    if colvar.is_file() and done.is_file():
        ran = _mdp_production_ps(rep / "prod.mdp")
        if ran is not None and abs(ran - production_ps) > 1e-6:
            # A 100 ps verification replicate and a 10 ns production replicate
            # both leave a COLVAR and a prod.gro. Pooling them would put a
            # replica that was never given time to escape into the same mean as
            # ones that were, and the convergence test would then be measuring
            # its own bookkeeping.
            log.warning("  rep%d: on disk at %.0f ps but %.0f ps was asked for; "
                        "re-running rather than pooling two protocols",
                        replicate, ran, production_ps)
        else:
            log.info("  rep%d: reusing the completed trajectory already on disk",
                     replicate)
            return colvar, {"reused": True, "velocity_seed":
                            gx.replicate_seed(cand_id, replicate)}
    res = gx.run_pipeline(src, md_wd, gpu_id=gpu, threads=threads,
                          production_ps=production_ps, replicate=replicate,
                          candidate_id=cand_id, plumed=plumed_txt,
                          reuse_equilibration=reuse_equilibration)
    if not res.get("plumed"):
        raise BPMDRunError("run_pipeline reports no PLUMED on production — the "
                           "trajectory is unbiased and is not a BPMD replicate")
    if not colvar.is_file():
        raise BPMDRunError(f"no COLVAR at {colvar}: mdrun accepted -plumed but "
                           "PLUMED wrote nothing")
    return colvar, {"reused": False,
                    "velocity_seed": res.get("velocity_seed"),
                    "ns_per_day": (res.get("stages", {}).get("prod", {})
                                   .get("ns_per_day"))}


def prepare_pose(cand: ns.Candidate, *, nrun: int, gpu: str, allow_redock: bool,
                 pose_rank: int = 1, sdf: Path | None = None) -> dict:
    """Everything up to the first replicate: pose, parameters, solvation, CV atoms.

    Solvation happens ONCE per pose and is shared by every replicate, so the
    replicates differ only in their velocity seed — which is what makes them
    replicates of one system rather than k loosely related simulations.
    """
    sdf, source = resolve_pose(cand, nrun=nrun, gpu=gpu,
                               allow_redock=allow_redock, sdf=sdf)
    mol, props = read_pose(sdf, pose_rank)

    # Each pose rank gets its OWN workdir. Sharing one would let pose 3's
    # solvated system be reused for pose 1 by the resume logic, and the run
    # would look complete while measuring the wrong pose.
    stem = cand.ident.replace(":", "_")
    wd = WORK / (stem if pose_rank == 1 else f"{stem}__p{pose_rank}")
    md_wd = wd / "md"
    cyx_index, _ = build_workdir(cand, mol, wd)
    if not (md_wd / "sys.gro").is_file():
        info = gx.solvate(wd, md_wd)
        log.info("  solvated: %d atoms, %d waters, %d ions, net charge %+.3f",
                 info["atoms"], info["waters"], info["ions"], info["net_charge"])

    # Which reactive centre is being presented is decided in the POSE's frame,
    # against the receptor PDB the pose was docked into. Using the solvated
    # system's sulfur here would compare coordinates across a rigid motion.
    cys_pdb = wd / "receptor_cys.pdb"
    sysres = resolve_system(md_wd, cys_pdb, cyx_index)
    centre = warhead_atom_in_pose(mol, cand, sg_in_pose_frame(cys_pdb, cyx_index))
    cv = cv_atoms(md_wd, sysres, mol, centre["warhead_pose_idx"])

    # THE POSE MUST STILL BE A NEAR-ATTACK CONFORMATION IN THE SYSTEM BEING RUN.
    # It is scored against a different sulfur position than the docking used, and
    # a pose that already sits past the unbound threshold has nothing for
    # metadynamics to push out of.
    if cv["start_distance_nm"] >= bpmd.UNBOUND_NM:
        raise BPMDRunError(
            f"the warhead starts {cv['start_distance_nm']:.2f} nm from SG, at or "
            f"past the unbound threshold {bpmd.UNBOUND_NM} nm — this pose is not "
            "in the pocket once the Cys113 sidechain is back in its "
            "crystallographic rotamer")

    docked = props.get("nac_distance_A")
    log.info("  CV: atoms %d,%d (PLUMED 1-based); start %.2f A vs %s A docked "
             "against the flexible sidechain; %s",
             cv["warhead_plumed"], cv["sg_plumed"], cv["start_distance_A"],
             docked, cv["atom_map"])

    plumed_txt = bpmd.plumed_input(cv["warhead_plumed"], cv["sg_plumed"])
    return {"cand": cand, "sdf": sdf, "pose_source": source, "mol": mol,
            "wd": wd, "md_wd": md_wd, "plumed": plumed_txt,
            "docked_distance_A": docked,
            "nac_viable": props.get("nac_viable"), **centre, **cv}


def run_pose(cand: ns.Candidate, *, replicates: int, production_ps: float,
             gpu: int, threads: int, nrun: int, dock_gpu: str,
             allow_redock: bool, on_row,
             reuse_equilibration: bool = False,
             pose_rank: int = 1, sdf: Path | None = None) -> list[dict]:
    """`replicates` replicates of one pose. Every replicate produces a row."""
    base = {"ident": cand.ident, "warhead_class": cand.warhead_class,
            "mechanism": cand.mechanism, "approach": cand.label,
            "production_ps": production_ps, "pose_rank": pose_rank}
    try:
        prep = prepare_pose(cand, nrun=nrun, gpu=dock_gpu,
                            allow_redock=allow_redock,
                            pose_rank=pose_rank, sdf=sdf)
    except Exception as exc:                              # noqa: BLE001
        # A pose that could not be SET UP is recorded once, for every replicate
        # it owed, so the shortfall is visible in the same table as the results.
        row = {**base, "replicate": 0, "status":
               f"setup failed: {type(exc).__name__}: {str(exc)[:200]}"}
        log.warning("  SETUP FAILED: %s", row["status"])
        on_row(row)
        return [row]

    meta = {k: v for k, v in prep.items()
            if isinstance(v, (int, float, str, bool, type(None)))}
    rows = []
    for k in range(1, replicates + 1):
        row = {**base, **meta, "replicate": k, "status": "ok"}
        try:
            colvar, info = run_replicate(cand.ident, prep["wd"], prep["md_wd"],
                                         prep["plumed"], k, gpu=gpu,
                                         threads=threads,
                                         production_ps=production_ps,
                                         reuse_equilibration=reuse_equilibration)
            r = bpmd.analyse_replica(colvar, k)
            t = _colvar_times(colvar)
            row.update({"colvar": str(colvar), "covered_ps": float(t[-1]),
                        "max_cv_nm": r.max_cv_nm,
                        "frac_in_window": r.frac_in_window,
                        "bias_at_exit_kj": r.bias_at_exit_kj,
                        "escaped": r.escaped, **info})
            log.info("  rep%d: frac_in_window %.3f, max %.2f nm, %s (%.0f ps)",
                     k, r.frac_in_window, r.max_cv_nm,
                     "ESCAPED" if r.escaped else "held", t[-1])
        except Exception as exc:                          # noqa: BLE001
            row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:200]}"
            log.warning("  rep%d FAILED: %s", k, row["status"])
        rows.append(row)
        on_row(row)
    return rows


def _mdp_production_ps(mdp: Path) -> float | None:
    """How long a replicate on disk was actually run for, from its own .mdp.

    `nsteps` is what mdrun was given, so it says what was asked for even when a
    PLUMED COMMITTOR stopped the run early — which is the case the COLVAR's last
    timestamp cannot distinguish from a job that died.
    """
    if not mdp.is_file():
        return None
    for line in mdp.read_text().splitlines():
        k, _, v = line.partition("=")
        if k.strip() == "nsteps":
            try:
                return int(v.strip()) / 500.0        # dt = 0.002 ps
            except ValueError:
                return None
    return None


def _colvar_times(path: Path) -> np.ndarray:
    ts = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        try:
            ts.append(float(line.split()[0]))
        except (IndexError, ValueError):
            continue
    if not ts:
        raise BPMDRunError(f"{path} carries no readable time column")
    return np.array(ts)


# --------------------------------------------------------------------------
# convergence
# --------------------------------------------------------------------------

def truncate_colvar(src: Path, t_ps: float, dest: Path) -> float:
    """Write the frames at time <= `t_ps`; return the last time actually present.

    Truncation is done on the FILE and handed back to `bpmd.analyse_replica`, so
    the shortened protocol is scored by exactly the tested code that scores the
    full one. Re-deriving the reduction here would let the two drift apart.
    """
    keep, last = [], 0.0
    for line in src.read_text().splitlines():
        if line.startswith("#"):
            keep.append(line)
            continue
        if not line.strip():
            continue
        try:
            t = float(line.split()[0])
        except (IndexError, ValueError):
            continue
        if t <= t_ps:
            keep.append(line)
            last = t
    if last == 0.0:
        raise BPMDRunError(f"{src} has no frames at or before {t_ps} ps")
    dest.write_text("\n".join(keep) + "\n")
    return last


def draw_subsets(k: int, n: int, draws: int, rng: np.random.Generator
                 ) -> list[tuple[int, ...]]:
    """Up to `draws` distinct n-subsets of k replicas — all of them when there are few.

    Enumerated exhaustively where possible so the reported spread is the true
    spread over every way the smaller protocol could have come out, rather than a
    sample of it. The seed is fixed: a convergence verdict must not change
    because the analysis was re-run.
    """
    all_subsets = list(itertools.combinations(range(k), n))
    if len(all_subsets) <= draws:
        return all_subsets
    idx = rng.choice(len(all_subsets), size=draws, replace=False)
    return [all_subsets[int(i)] for i in sorted(idx)]


def convergence_report(rows: list[dict], truncations, subsamples) -> pd.DataFrame:
    """How the stability score moves as replicas are dropped and length is cut."""
    ok = [r for r in rows if r["status"] == "ok" and r.get("colvar")]
    if len(ok) < 2:
        raise BPMDRunError(f"only {len(ok)} usable replicas — nothing to converge")
    colvars = [Path(r["colvar"]) for r in ok]
    k = len(colvars)
    rng = np.random.default_rng(SUBSET_SEED)

    out = []
    scratch = Path(tempfile.mkdtemp(prefix="bpmd_trunc_"))
    try:
        for t in truncations:
            reps, covered = [], []
            for i, cv in enumerate(colvars):
                dest = scratch / f"COLVAR_{i}_{int(t)}"
                covered.append(truncate_colvar(cv, t, dest))
                reps.append(bpmd.analyse_replica(dest, i))
            full = bpmd.combine(reps)
            for n in [s for s in subsamples if s <= k] + ([k] if k not in subsamples else []):
                scores = []
                for sub in draw_subsets(k, n, SUBSET_DRAWS, rng):
                    s = bpmd.combine([reps[i] for i in sub])
                    scores.append(s.score)
                scores = np.array(scores)
                out.append({
                    "truncation_ps": t,
                    "median_covered_ps": float(np.median(covered)),
                    "n_replicas": n,
                    "n_subsets": len(scores),
                    "score_mean": float(scores.mean()),
                    "score_sd": float(scores.std(ddof=1)) if len(scores) > 1 else 0.0,
                    "score_min": float(scores.min()),
                    "score_max": float(scores.max()),
                    "full_score": full.score,
                    "full_mean_frac_in_window": full.mean_frac_in_window,
                    "full_spread_frac": full.spread_frac,
                    "full_median_bias_kj": full.median_bias_to_escape_kj,
                    "n_escaped": full.n_escaped,
                })
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return pd.DataFrame(out)


def print_convergence(df: pd.DataFrame, ident: str) -> None:
    ref = df[(df.truncation_ps == df.truncation_ps.max())
             & (df.n_replicas == df.n_replicas.max())]
    full = float(ref.full_score.iloc[0]) if len(ref) else float("nan")
    print(f"\n=== BPMD convergence: {ident} ===")
    print(f"  reference protocol: {int(df.n_replicas.max())} replicas x "
          f"{df.truncation_ps.max()/1000:.0f} ns  ->  score {full:.4f}")
    print(f"  (mean fraction in the near-attack window "
          f"{ref.full_mean_frac_in_window.iloc[0]:.3f} +/- "
          f"{ref.full_spread_frac.iloc[0]:.3f} across replicas; "
          f"{int(ref.n_escaped.iloc[0])} of {int(df.n_replicas.max())} escaped)\n")
    print(f"  score at each (replicas x length), mean over subsets "
          f"[min, max]; deviation from the reference in the last column")
    print(f"  {'ns':>5}  {'k':>3}  {'mean':>8}  {'sd':>7}  "
          f"{'min':>8}  {'max':>8}  {'vs ref':>8}")
    for t in sorted(df.truncation_ps.unique()):
        for _, r in df[df.truncation_ps == t].sort_values("n_replicas").iterrows():
            dev = (r.score_mean - full) / full if full else float("nan")
            print(f"  {t/1000:>5.0f}  {int(r.n_replicas):>3}  {r.score_mean:>8.4f}  "
                  f"{r.score_sd:>7.4f}  {r.score_min:>8.4f}  {r.score_max:>8.4f}  "
                  f"{dev:>+7.1%}")
    print("\n  READ THIS AS: the cheapest (k, length) whose spread still contains "
          "the\n  reference score is the cheapest DEFENSIBLE protocol. A cell "
          "whose min-max\n  range excludes the reference is a protocol that "
          "would have given a\n  different answer, however plausible it looks on "
          "its own.")


# --------------------------------------------------------------------------
# resumability
# --------------------------------------------------------------------------

def _chunk_files(d: Path, pattern: str) -> list[str]:
    """Chunk files in VERSION order — `_10.csv` sorts before `_2.csv` lexically."""
    def key(f: str) -> tuple[int, int]:
        m = re.search(r"_s(\d+)_(\d+)\.csv$", f)
        return (int(m.group(1)), int(m.group(2))) if m else (0, 0)
    return sorted(glob.glob(str(d / pattern)), key=key)


def already_done() -> set[tuple[str, int]]:
    """(ident, replicate) pairs that have already landed OK in a chunk file."""
    done = set()
    for f in _chunk_files(OUT.dir, "bpmd_s*.csv"):
        try:
            df = pd.read_csv(f)
        except Exception as exc:                          # noqa: BLE001
            log.warning("unreadable chunk %s: %s", Path(f).name, exc)
            continue
        if not {"ident", "replicate", "status"} <= set(df.columns):
            continue
        ok = df[df.status == "ok"]
        done.update(zip(ok.ident.astype(str), ok.replicate.astype(int)))
    return done


def load_results() -> pd.DataFrame:
    fs = _chunk_files(OUT.dir, "bpmd_s*.csv")
    if not fs:
        return pd.DataFrame()
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return df.drop_duplicates(["ident", "replicate"], keep="last")


def report() -> None:
    df = load_results()
    if df.empty:
        raise SystemExit("nothing written yet")
    ok = df[df.status == "ok"]
    print(f"\n=== BPMD: {len(df)} replicate records, {len(ok)} ok, "
          f"{len(df) - len(ok)} failed, over {df.ident.nunique()} poses ===")
    if len(ok) != len(df):
        print("\n  failures (recorded, never dropped):")
        print(df[df.status != "ok"].status.str.slice(0, 70)
              .value_counts().head(10).to_string())
    if ok.empty:
        return
    print(f"\n  {'pose':<28}{'k':>3}  {'score':>8}  {'frac':>6}  {'spread':>7}  "
          f"{'esc':>4}")
    stats = []
    for ident, g in ok.groupby("ident"):
        res = [bpmd.ReplicaResult(int(r.replicate), float(r.max_cv_nm),
                                  float(r.frac_in_window),
                                  float(r.bias_at_exit_kj), bool(r.escaped))
               for r in g.itertuples()]
        s = bpmd.combine(res)
        stats.append((s.score, ident, s))
    for score, ident, s in sorted(stats, reverse=True):
        print(f"  {ident[:27]:<28}{s.n_replicas:>3}  {score:>8.4f}  "
              f"{s.mean_frac_in_window:>6.3f}  {s.spread_frac:>7.3f}  "
              f"{s.n_escaped:>4}")
    print("\n  The spread IS the uncertainty on the score. Do not read an "
          "ordering\n  finer than it, and do not quote a score without the "
          "replica count and\n  trajectory length that produced it (D0068).")


# --------------------------------------------------------------------------

def pose_idents_on_disk() -> list[str]:
    """Every exported pose, in a deterministic order so shards are stable."""
    return sorted(p.stem for p in POSES.glob("*.sdf"))


def choose_convergence_pose(by_id: dict[str, ns.Candidate]) -> ns.Candidate:
    """A crystallographic positive — a molecule KNOWN to react with Cys113.

    The convergence test asks how much protocol is needed to reproduce an answer,
    so the pose underneath it has to be one where "stable" is the expected answer
    for an independent reason. `xtal:` idents are ligands with an OBSERVED
    covalent bond to Cys113 at 1.6-1.9 A, which is the only such reason this
    project has.

    NOTE, because it cost time to find: the 50 poses under `nac_poses/` are ALL
    T_3/T_4 generated candidates. `export_nac_poses` walks the shortlists, and no
    crystallographic positive is on a shortlist, so none was ever exported. The
    positive is therefore re-docked here under the same protocol rather than
    quietly substituting a generated molecule for a known active.

    RESTRICTED TO CHLOROACETAMIDE, and for a measured reason. It is the one
    warhead class whose near-attack criterion still discriminates once the
    docking search has converged: AUC 0.756 (p = 0.016) at 2,000 runs, against
    Michael's 0.750 (p = 0.065, not significant) and SNAr's 0.575 (n.s.) — D0068,
    on the same 75 molecules. It also survived ten disjoint draws of negatives at
    0.908 (`nac_robustness`). Resting a convergence test on a class whose
    geometry we cannot reliably score would confound two questions at once.

    A NOTE ON THE IDENT, because it reads wrong at a glance. The first
    chloroacetamide positive by ident is `xtal:6VAJ:QT7`, and D0059 invalidated
    **6VAJ as a receptor**. This uses QT7 as a MOLECULE — a ligand observed
    covalently bonded to Cys113 — and docks it into 3IKD like everything else on
    this branch. Nothing about the 6VAJ structure enters the calculation.

    Deterministic by ident so a re-run picks the same molecule. If its pose fails
    to set up (the driver refuses a pose that is not in the pocket once Cys113 is
    back in its crystallographic rotamer), that is recorded as a failure and the
    next one is chosen by hand with `--pose`, rather than the code quietly
    walking down the list until something works.
    """
    validated = {"sn2_displacement"}
    xtal = sorted((c for i, c in by_id.items()
                   if i.startswith("xtal:") and c.mechanism in validated),
                  key=lambda c: c.ident)
    if not xtal:
        raise BPMDRunError(
            "no crystallographic positive in a validated mechanism class; pass "
            "--pose explicitly and say why in the run's notes")
    return xtal[0]


class _ChunkWriter:
    """Buffered, versioned chunk files — never appended to, never overwritten.

    A replicate costs ~20 minutes of GPU, so a run that dies must not cost the
    replicates that already finished. Each flush allocates the next version under
    the append-only root, which is the same contract `nac_rank` runs under.
    """

    def __init__(self, topic: sout.Topic, stem: str, size: int) -> None:
        self.topic, self.stem, self.size, self.buf = topic, stem, max(1, size), []

    def add(self, row: dict) -> None:
        self.buf.append(row)
        if len(self.buf) >= self.size:
            self.flush()

    def flush(self) -> None:
        if not self.buf:
            return
        dest = self.topic.write(self.stem, ".csv")
        pd.DataFrame(self.buf).to_csv(dest, index=False)
        ok = sum(r.get("status") == "ok" for r in self.buf)
        log.info("wrote %d records (%d ok) -> %s", len(self.buf), ok, dest.name)
        self.buf = []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--convergence", action="store_true",
                    help="measure how the score depends on replica count and "
                         "length, on ONE known-active pose, before any "
                         "production run commits to a protocol")
    ap.add_argument("--pose", default=None, metavar="IDENT",
                    help="a single pose to run (default for --convergence: a "
                         "crystallographic positive)")
    ap.add_argument("--replicates", type=int, default=None,
                    help=f"replicates per pose (default {FULL_REPLICAS} under "
                         "--convergence, else 3)")
    ap.add_argument("--production-ps", type=float, default=None,
                    help="SHORT by default (100 ps) so the chain can be proven "
                         f"cheaply; --convergence defaults to "
                         f"{FULL_PRODUCTION_PS:.0f}")
    ap.add_argument("--gpu", type=int, default=2,
                    help="ONE GPU, and not 0 or 4 — other people's jobs")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dock-gpu", default=None,
                    help="GPU for the re-dock, if a pose has to be built "
                         "(defaults to --gpu)")
    ap.add_argument("--nrun", type=int, default=200,
                    help="docking runs, if a pose has to be re-docked")
    ap.add_argument("--no-redock", action="store_true",
                    help="refuse to build a missing pose rather than docking it")
    ap.add_argument("--limit", type=int, default=None,
                    help="how many poses; 1 verifies the plumbing")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=5,
                    help="replicate records per chunk file")
    ap.add_argument("--report", action="store_true",
                    help="summarise what has actually landed, and exit")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [s{args.shard}] %(message)s")

    if args.report:
        report()
        return

    replicates = args.replicates or (FULL_REPLICAS if args.convergence else 3)
    production_ps = args.production_ps or (FULL_PRODUCTION_PS if args.convergence
                                           else 100.0)
    dock_gpu = args.dock_gpu or str(args.gpu)

    # Resolved before any GPU time: "PLUMED is not installed" must not surface
    # three hours into a run that was never actually biased.
    log.info("PLUMED kernel: %s", gx.plumed_kernel())

    by_id = candidate_index()

    if args.convergence:
        cand = by_id[args.pose] if args.pose else choose_convergence_pose(by_id)
        if args.pose and args.pose not in by_id:
            raise SystemExit(f"unknown ident {args.pose}")
        log.info("convergence on %s (%s, %s): %d replicas x %.0f ps",
                 cand.ident, cand.label, cand.mechanism, replicates,
                 production_ps)
        chunk = _ChunkWriter(OUT, f"bpmd_s{args.shard}", args.chunk)
        rows = run_pose(cand, replicates=replicates, production_ps=production_ps,
                        gpu=args.gpu, threads=args.threads, nrun=args.nrun,
                        dock_gpu=dock_gpu, allow_redock=not args.no_redock,
                        on_row=chunk.add)
        chunk.flush()
        pd.DataFrame(rows).to_csv(
            CONV.write(f"bpmd_conv_replicas_{cand.ident.replace(':', '_')}", ".csv"),
            index=False)
        trunc = tuple(t for t in TRUNCATIONS_PS if t <= production_ps) or (production_ps,)
        df = convergence_report(rows, trunc, SUBSAMPLES)
        dest = CONV.write(f"bpmd_conv_{cand.ident.replace(':', '_')}", ".csv")
        df.to_csv(dest, index=False)
        print_convergence(df, cand.ident)
        print(f"\n  -> {dest}")
        return

    idents = [args.pose] if args.pose else pose_idents_on_disk()
    idents = [i for k, i in enumerate(idents) if k % args.n_shards == args.shard]
    if args.limit:
        idents = idents[:args.limit]
    done = already_done()
    log.info("%d poses assigned to shard %d, %d replicates each at %.0f ps",
             len(idents), args.shard, replicates, production_ps)

    chunk = _ChunkWriter(OUT, f"bpmd_s{args.shard}", args.chunk)
    for i, ident in enumerate(idents, 1):
        cand = by_id.get(ident)
        if cand is None:
            chunk.add({"ident": ident, "replicate": 0,
                       "status": "failed: not in the candidate set"})
            continue
        todo = [k for k in range(1, replicates + 1) if (ident, k) not in done]
        if not todo:
            log.info("[%d/%d] %s: all %d replicates already done", i, len(idents),
                     ident, replicates)
            continue
        log.info("[%d/%d] %s (%s): %d replicates to do", i, len(idents), ident,
                 cand.mechanism, len(todo))
        run_pose(cand, replicates=replicates, production_ps=production_ps,
                 gpu=args.gpu, threads=args.threads, nrun=args.nrun,
                 dock_gpu=dock_gpu, allow_redock=not args.no_redock,
                 on_row=chunk.add)
    chunk.flush()
    report()


if __name__ == "__main__":
    main()
