#!/usr/bin/env python3
"""
Purpose: typed non-covalent interactions per frame, per residue, over an MD run.
Author: Timothy Wu (with Claude Code)
Date: 2026-09-09
Input: a PBC-corrected multi-model PDB (protein + MOL) and the ligand's SDF
Output: a tidy frame x residue x interaction-type table

@twu383, 2026-09-09: *"I have many preferential factors that I can determine
visually but need to find ways to quantitatively describe them for automated
screening. for example we need to profile non-cov interactions of poses at
time x for md runs"*.

WHAT ALREADY EXISTED AND WHY IT IS NOT ENOUGH. `shortlist_report.contacts()`
gives per-residue contact FREQUENCY across frames, and its docstring makes the
right argument -- one frame is close to arbitrary for a molecule that moves.
But it is distance-only: a residue is "in contact", and polar contacts are split
off by element. That cannot express "the aryl stacks on His59" or "it bridges
two of the basic triad", which are the things a chemist reads off the screen.
This types the interaction instead of just measuring proximity.

NO HYDROGENS. `gmx trjconv` writes heavy atoms only, so every criterion here is
a heavy-atom one: an H-bond is a donor-acceptor distance plus an angle taken
through the donor's own root atom, which is the standard treatment for
structures without explicit H. Getting this wrong in the other direction --
assuming hydrogens and silently finding none -- would report zero H-bonds for
every pose, which reads as a chemistry result rather than a parsing fault
(catalogue #31).

TYPING IS BY IDENTITY, NEVER BY POSITION OR ELEMENT ALONE. Protein donors,
acceptors, charges and rings come from a per-residue table keyed on
(resname, atom name); an unknown residue RAISES rather than being silently
treated as apolar. The ligand's are derived from its own SDF through RDKit, so
they follow the molecule rather than a guess from the PDB record.

The thresholds are conventional and are stated as module constants so a caller
can see and change them, not buried in the comparisons.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

# --- criteria, all heavy-atom -------------------------------------------------
#: Apolar C/S within this distance counts as a hydrophobic contact.
HYDROPHOBIC_A = 4.0
#: Donor heavy atom to acceptor heavy atom.
HBOND_A = 3.5
#: Angle at the donor: root-donor...acceptor. Loose because there is no H to
#: place; 90 deg only excludes an acceptor sitting behind the donor's root.
HBOND_MIN_ANGLE = 90.0
#: Charged-group centroid separation for a salt bridge.
SALT_BRIDGE_A = 4.0
#: Aromatic ring centroid separation.
PI_STACK_A = 5.5
#: Interplanar angle below this is face-to-face; above (90 - this) is T-shaped.
PI_PARALLEL_MAX = 30.0
#: Cation to aromatic ring centroid.
PI_CATION_A = 6.0
#: Halogen ... acceptor, with a near-linear C-X...A.
HALOGEN_A = 3.5
HALOGEN_MIN_ANGLE = 140.0

TYPES = ("hydrophobic", "hbond_protein_donor", "hbond_protein_acceptor",
         "salt_bridge", "pi_stack", "pi_cation", "halogen")

#: Protein H-bond donors and acceptors by (resname, atom). Backbone N donates
#: and backbone O accepts for every residue, so they are added separately.
#:
#: HID / HIE ARE DISTINCT AND THAT IS THE POINT. Amber writes the protonation
#: state into the residue name: HID has H on ND1 (ND1 donates, NE2 accepts) and
#: HIE has it on NE2 (the reverse). Treating both as "HIS" would invert the
#: donor/acceptor assignment on half the histidines in the structure -- and
#: His59 and His157 are the catalytic pair here, so that is not a rounding
#: error. `md_movie`/`elevation_report` already carry HIS_FORMS for the same
#: reason.
SIDECHAIN_DONORS = {
    "ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "TRP": ("NE1",),
    "ASN": ("ND2",), "GLN": ("NE2",), "HID": ("ND1",), "HIE": ("NE2",),
    "HIP": ("ND1", "NE2"), "SER": ("OG",), "THR": ("OG1",), "TYR": ("OH",),
    "CYS": ("SG",),
}
SIDECHAIN_ACCEPTORS = {
    "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "ASN": ("OD1",),
    "GLN": ("OE1",), "HID": ("NE2",), "HIE": ("ND1",),
    "SER": ("OG",), "THR": ("OG1",), "TYR": ("OH",), "MET": ("SD",),
    "CYS": ("SG",),
}
#: Formally charged sidechain groups, as the atoms whose centroid is used.
CATIONIC = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "HIP": ("ND1", "NE2")}
ANIONIC = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
#: Aromatic sidechain rings, as the atoms defining the ring plane.
AROMATIC_RINGS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HID": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIE": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "HIP": ("CG", "ND1", "CD2", "CE1", "NE2"),
}
#: Apolar sidechain carbons/sulfur, for hydrophobic contact.
APOLAR_ELEMENTS = ("C", "S")

#: Every residue name this module is prepared to type. An unknown one RAISES:
#: treating it as apolar would silently drop every polar interaction it makes,
#: and a missing interaction is invisible in a way a wrong one is not.
KNOWN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY",
    "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}


class ProfileError(RuntimeError):
    """Raised when a frame cannot be typed, rather than typed wrongly."""


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Angle a-b-c in degrees."""
    v1, v2 = a - b, c - b
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 0.0
    return float(np.degrees(np.arccos(np.clip(v1.dot(v2) / (n1 * n2), -1, 1))))


def _plane_normal(pts: np.ndarray) -> np.ndarray:
    """Best-fit plane normal for a ring, by SVD on the centred points."""
    c = pts - pts.mean(axis=0)
    return np.linalg.svd(c)[2][2]


def ligand_chemistry(sdf: Path) -> dict:
    """Donor / acceptor / apolar / aromatic / halogen / charged atom indices.

    FROM THE MOLECULE, NOT FROM THE PDB RECORD. The movie's MOL residue carries
    element letters and nothing else -- no bond orders, no aromaticity, no formal
    charges -- so typing the ligand from it would miss every aromatic ring and
    every charge. The SDF is the same file the MD was parameterised from, and its
    heavy-atom ORDER is the order the MOL residue appears in (the invariant
    `mdprio_report.nac_series` also relies on), which is what lets these indices
    address frame coordinates.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Lipinski
    RDLogger.DisableLog("rdApp.*")

    mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    if mol is None:
        raise ProfileError(f"RDKit could not read {sdf}")
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {old: new for new, old in enumerate(heavy)}      # sdf idx -> frame idx

    don = {pos[i] for (i,) in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[$([N;!H0]),$([O;H1]),$([S;H1])]")) if i in pos}
    acc = {pos[i] for (i,) in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[$([O]),$([N;!$(N-[!#6;!#1]);!$([N+])]),$([S;X2])]"))
        if i in pos}
    hal = {pos[i] for (i,) in mol.GetSubstructMatches(
        Chem.MolFromSmarts("[F,Cl,Br,I]")) if i in pos}
    apolar = {pos[a.GetIdx()] for a in mol.GetAtoms()
              if a.GetIdx() in pos and a.GetSymbol() in APOLAR_ELEMENTS
              and all(n.GetSymbol() in ("C", "S", "H") for n in a.GetNeighbors())}
    cat = {pos[a.GetIdx()] for a in mol.GetAtoms()
           if a.GetIdx() in pos and a.GetFormalCharge() > 0}
    ani = {pos[a.GetIdx()] for a in mol.GetAtoms()
           if a.GetIdx() in pos and a.GetFormalCharge() < 0}
    rings = [[pos[i] for i in r] for r in mol.GetRingInfo().AtomRings()
             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
             and all(i in pos for i in r)]
    # the halogen's own root carbon, for the C-X...A angle
    hal_root = {}
    for i in hal:
        old = heavy[i]
        nb = [n.GetIdx() for n in mol.GetAtomWithIdx(old).GetNeighbors()
              if n.GetIdx() in pos]
        if nb:
            hal_root[i] = pos[nb[0]]
    return {"n_heavy": len(heavy), "donors": don, "acceptors": acc,
            "apolar": apolar, "halogens": hal, "halogen_roots": hal_root,
            "cations": cat, "anions": ani, "aromatic_rings": rings}


def _read_frames(pdb: Path, stride: int = 1):
    """Yield (frame_index, ligand_xyz, protein_atoms) per MODEL.

    `protein_atoms` is a list of (resname, resid, atom_name, element, xyz).
    """
    lig, prot, idx = [], [], 0
    for ln in Path(pdb).read_text(errors="replace").splitlines():
        if ln.startswith("MODEL"):
            lig, prot = [], []
        elif ln.startswith("ENDMDL"):
            if lig and prot:
                if idx % stride == 0:
                    yield idx, np.array(lig), prot
                idx += 1
        elif ln.startswith(("ATOM", "HETATM")):
            el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if el == "H":
                continue
            try:
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                ri = int(ln[22:26])
            except ValueError:
                continue
            rn, nm = ln[17:20].strip(), ln[12:16].strip()
            if rn == "MOL":
                lig.append(xyz)
            else:
                prot.append((rn, ri, nm, el, xyz))


def profile(movie_pdb: Path, ligand_sdf: Path, *, stride: int = 1,
            total_ps: float | None = None, offset: int = 50):
    """Typed interactions for every frame. Returns a tidy DataFrame.

    Columns: frame, time_ns (when `total_ps` is given), resname, resid_pdb,
    interaction, distance_a, detail.

    `resid_pdb` is the PIN1 number (MD resid + `offset`), so a row can be read
    against the literature without the reader having to know the renumbering --
    the same +50 every other module applies through `md_movie.PIN1_OFFSET`.
    """
    import pandas as pd

    chem = ligand_chemistry(Path(ligand_sdf))
    rows = []
    n_frames = 0
    for fi, L, prot in _read_frames(Path(movie_pdb), stride):
        n_frames += 1
        if L.shape[0] != chem["n_heavy"]:
            raise ProfileError(
                f"{Path(movie_pdb).name} frame {fi} has {L.shape[0]} MOL atoms "
                f"but {Path(ligand_sdf).name} has {chem['n_heavy']} heavy atoms; "
                f"these are not the same molecule and the indices would address "
                f"the wrong atoms")
        unknown = {r[0] for r in prot} - KNOWN_RESIDUES
        if unknown:
            raise ProfileError(
                f"untypable residue(s) {sorted(unknown)}: treating them as "
                f"apolar would silently drop every polar interaction they make")

        P = np.array([r[4] for r in prot])
        d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2)
        by_res: dict[tuple, list[int]] = {}
        for i, (rn, ri, nm, el, _) in enumerate(prot):
            by_res.setdefault((rn, ri), []).append(i)

        for (rn, ri), idxs in by_res.items():
            if d[idxs].min() > PI_CATION_A:          # nothing can reach
                continue
            name = {i: prot[i][2] for i in idxs}
            elem = {i: prot[i][3] for i in idxs}
            hit = []

            # --- hydrophobic ------------------------------------------------
            ap = [i for i in idxs if elem[i] in APOLAR_ELEMENTS]
            if ap and chem["apolar"]:
                sub = d[np.ix_(ap, sorted(chem["apolar"]))]
                if sub.size and sub.min() <= HYDROPHOBIC_A:
                    hit.append(("hydrophobic", float(sub.min()), ""))

            # --- H-bonds, both directions -----------------------------------
            don_at = [i for i in idxs
                      if name[i] == "N" or name[i] in SIDECHAIN_DONORS.get(rn, ())]
            acc_at = [i for i in idxs
                      if name[i] == "O" or name[i] in SIDECHAIN_ACCEPTORS.get(rn, ())]
            for atoms, partner, label in (
                    (don_at, chem["acceptors"], "hbond_protein_donor"),
                    (acc_at, chem["donors"], "hbond_protein_acceptor")):
                if not atoms or not partner:
                    continue
                sub = d[np.ix_(atoms, sorted(partner))]
                if sub.size and sub.min() <= HBOND_A:
                    a, b = np.unravel_index(sub.argmin(), sub.shape)
                    hit.append((label, float(sub.min()),
                                f"{name[atoms[a]]}...lig{sorted(partner)[b]}"))

            # --- salt bridges -----------------------------------------------
            for group, lig_group, lab in ((CATIONIC.get(rn), chem["anions"], "+/-"),
                                          (ANIONIC.get(rn), chem["cations"], "-/+")):
                if not group or not lig_group:
                    continue
                g = [i for i in idxs if name[i] in group]
                if not g:
                    continue
                c1 = P[g].mean(axis=0)
                c2 = L[sorted(lig_group)].mean(axis=0)
                sep = float(np.linalg.norm(c1 - c2))
                if sep <= SALT_BRIDGE_A:
                    hit.append(("salt_bridge", sep, lab))

            # --- pi stacking and pi-cation ----------------------------------
            ring_names = AROMATIC_RINGS.get(rn)
            if ring_names:
                g = [i for i in idxs if name[i] in ring_names]
                if len(g) >= 5:
                    rc, rn_norm = P[g].mean(axis=0), _plane_normal(P[g])
                    for lr in chem["aromatic_rings"]:
                        lc = L[lr].mean(axis=0)
                        sep = float(np.linalg.norm(rc - lc))
                        if sep <= PI_STACK_A:
                            ang = _angle(rc + rn_norm, rc, lc + _plane_normal(L[lr]))
                            ang = min(ang, 180 - ang)
                            kind = ("face-to-face" if ang <= PI_PARALLEL_MAX
                                    else "T-shaped" if ang >= 90 - PI_PARALLEL_MAX
                                    else "offset")
                            hit.append(("pi_stack", sep, f"{ang:.0f}deg {kind}"))
                    for ci in chem["cations"]:
                        sep = float(np.linalg.norm(rc - L[ci]))
                        if sep <= PI_CATION_A:
                            hit.append(("pi_cation", sep, "lig cation"))
            # protein cation onto a ligand ring
            cat_names = CATIONIC.get(rn)
            if cat_names and chem["aromatic_rings"]:
                g = [i for i in idxs if name[i] in cat_names]
                if g:
                    cc = P[g].mean(axis=0)
                    for lr in chem["aromatic_rings"]:
                        sep = float(np.linalg.norm(cc - L[lr].mean(axis=0)))
                        if sep <= PI_CATION_A:
                            hit.append(("pi_cation", sep, "protein cation"))

            # --- halogen bonds ----------------------------------------------
            if chem["halogens"] and acc_at:
                for xi in chem["halogens"]:
                    root = chem["halogen_roots"].get(xi)
                    if root is None:
                        continue
                    for ai in acc_at:
                        sep = float(d[ai, xi])
                        if sep <= HALOGEN_A:
                            ang = _angle(L[root], L[xi], P[ai])
                            if ang >= HALOGEN_MIN_ANGLE:
                                hit.append(("halogen", sep, f"{ang:.0f}deg"))

            for kind, dist, detail in hit:
                rows.append(dict(frame=fi, resname=rn, resid_md=ri,
                                 resid_pdb=ri + offset, interaction=kind,
                                 distance_a=round(dist, 2), detail=detail))

    # A RESULT WITH NO ROWS STILL HAS A SCHEMA. Returning a bare
    # `DataFrame([])` gives a frame with no COLUMNS, so `df.interaction` raises
    # AttributeError on the perfectly ordinary case of a pose that makes no
    # typed interaction -- and a caller then has to guard emptiness before it
    # can ask any question. "None found" and "cannot be asked" are different
    # answers and this returns the first.
    COLUMNS = ["frame", "resname", "resid_md", "resid_pdb",
               "interaction", "distance_a", "detail"]
    df = pd.DataFrame(rows, columns=COLUMNS) if rows else pd.DataFrame(
        {c: pd.Series(dtype="object") for c in COLUMNS})
    if total_ps and n_frames > 1 and len(df):
        df["time_ns"] = df.frame * (total_ps / 1000.0) / (n_frames - 1)
    if not len(df):
        df["time_ns"] = pd.Series(dtype=float)
    df.attrs["n_frames"] = n_frames
    return df


def occupancy(df, n_frames: int | None = None):
    """Fraction of frames each (residue, interaction) is present.

    THE SCREENABLE FORM. A single frame says which residues happened to be near
    the ligand at one instant; the fraction over a run is what a rule can be
    written against -- the same argument `shortlist_report.contacts` makes for
    plain proximity, applied to typed interactions.
    """
    import pandas as pd
    n = n_frames or df.attrs.get("n_frames")
    if not n:
        raise ProfileError("occupancy needs the frame count; pass n_frames")
    if not len(df):
        return pd.DataFrame(columns=["resname", "resid_pdb", "interaction",
                                     "occupancy", "median_distance_a"])
    g = (df.drop_duplicates(["frame", "resid_pdb", "interaction"])
           .groupby(["resname", "resid_pdb", "interaction"])
           .agg(n_frames_present=("frame", "size"),
                median_distance_a=("distance_a", "median"))
           .reset_index())
    g["occupancy"] = g.n_frames_present / n
    return g.sort_values("occupancy", ascending=False).reset_index(drop=True)
