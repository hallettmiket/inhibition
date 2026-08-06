"""
Purpose: the internal strain a docked pose imposes on the ligand, in kcal/mol.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: a pose (coordinates) + the ligand's SMILES
Output: strain energy, and a pose-selection rule built on it

WHY STRAIN AND NOT MORE GEOMETRY. Every pose-selection rule tested in D0061 was,
in one form or another, a re-reading of the same ensemble Vina produced: contact
profiles, cluster membership, centroid position. None beat random, and the
rationale doc's stated reading is that poses inside one ensemble are not
separable by cheap geometry.

Strain is different in kind. It is a property of the **ligand alone** -- the
energetic cost of holding it in the docked conformation rather than a relaxed
one -- so it is orthogonal to the intermolecular score by construction. Docking
optimises the protein-ligand interaction and pays for internal distortion only
weakly; a pose that fits beautifully by burying surface while twisting a bond
into an eclipsed conformation is exactly the artefact this catches, and exactly
the artefact a contact profile cannot see.

PRECEDENT. Strain filtering is established rather than novel: POSIT terminates a
pose when strain exceeds ~10 kcal/mol, and the strain-energy literature it cites
uses thresholds in that region. That threshold is a REFERENCE POINT, not a
transplant -- ours is measured on our own benchmark before use, because a cutoff
imported unmeasured is the shape this project keeps writing (catalogue #19).

HOW IT IS COMPUTED, AND THE HONEST LIMITS.

    strain = E(MMFF, pose coordinates, bonds relaxed) - E(MMFF, fully relaxed)

The pose's heavy-atom positions are held while hydrogens and bond lengths relax,
because docked poses carry no hydrogens and their bond lengths come from a
force field that is not MMFF -- comparing raw MMFF energies of a docked pose
against a minimised one would measure the difference between two force fields,
not the strain. Only the TORSIONAL and steric penalty of the docked arrangement
survives that, which is the part that means anything.

It is a **local** relaxation, not a global conformational search, so this is an
UPPER BOUND on the true relaxed energy and therefore a LOWER BOUND on strain.
Reported as such; a molecule with a genuinely better conformer elsewhere in
torsion space will look less strained than it is.
"""

from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)

# POSIT's reference point, recorded so the number's origin is traceable. NOT
# adopted as our threshold -- see the module docstring.
POSIT_REFERENCE_KCAL = 10.0


class StrainError(RuntimeError):
    """The pose could not be typed well enough to carry a force field."""


def _mmff(mol):
    from rdkit.Chem import AllChem
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise StrainError("MMFF cannot type this molecule")
    return props


def strain_kcal(mol, max_iters: int = 500) -> float:
    """Strain of `mol`'s current conformer, in kcal/mol. LOWER is more relaxed.

    `mol` must carry explicit hydrogens and one conformer holding the pose.

    Returns `E(pose, hydrogens relaxed) - E(locally minimised)`. Both terms come
    from the SAME force field on the SAME molecule, so the force-field offset
    cancels and what remains is the strain the pose imposes.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    work = Chem.Mol(mol)
    props = _mmff(work)

    # 1. Relax hydrogens only, with heavy atoms pinned. Docked poses have no
    #    hydrogens; obabel adds them in arbitrary positions, and their clashes
    #    would otherwise be counted as ligand strain.
    ff = AllChem.MMFFGetMoleculeForceField(work, props)
    if ff is None:
        raise StrainError("could not build an MMFF force field")
    for a in work.GetAtoms():
        if a.GetAtomicNum() > 1:
            ff.AddFixedPoint(a.GetIdx())
    ff.Minimize(maxIts=max_iters)
    e_pose = ff.CalcEnergy()

    # 2. Relax everything from the pose. Local, deliberately: a global search
    #    would find a conformer the molecule may never adopt in the pocket, and
    #    would make strain a property of the search rather than the pose.
    ff2 = AllChem.MMFFGetMoleculeForceField(Chem.Mol(work), props)
    if ff2 is None:
        raise StrainError("could not build an MMFF force field for relaxation")
    ff2.Minimize(maxIts=max_iters * 4)
    e_relaxed = ff2.CalcEnergy()

    return float(e_pose - e_relaxed)


def pose_strains(mol_with_confs) -> list[float]:
    """Strain for every conformer on a molecule, in conformer order."""
    from rdkit import Chem

    out = []
    for cid in range(mol_with_confs.GetNumConformers()):
        one = Chem.Mol(mol_with_confs)
        one.RemoveAllConformers()
        one.AddConformer(mol_with_confs.GetConformer(cid), assignId=True)
        try:
            out.append(strain_kcal(one))
        except Exception as exc:            # noqa: BLE001
            log.debug("strain failed for conformer %d: %s", cid, exc)
            out.append(float("nan"))
    return out


def least_strained(strains: list[float]) -> int:
    """Index of the least-strained pose.

    NaNs are poses the force field could not type. They are ranked LAST rather
    than treated as zero-strain -- "could not be measured" and "is not strained"
    are different facts, and conflating them would silently promote exactly the
    molecules the force field found unreasonable.
    """
    arr = np.asarray(strains, dtype=float)
    if np.all(np.isnan(arr)):
        raise StrainError("no conformer could be typed")
    return int(np.nanargmin(arr))
