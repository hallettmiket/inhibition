"""
Purpose: Turn a pre-reaction candidate into the adduct-form ligand gnina must dock.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: candidate SMILES + warhead class
Output: adduct-form SMILES and the SMARTS locating its attachment atom

WHY (D0022). gnina's `--covalent_lig_atom_pattern` bonds a ligand atom to the
receptor atom. It does not perform reaction chemistry and it does not remove a
leaving group — it replaces an **implicit hydrogen** on the matched atom. Docking
the intact pre-reaction molecule therefore produced 1,683 complexes in which the
reactive carbon still carried the leaving group it is supposed to have lost, and
for leaving groups fixed on a rigid sp2 or aromatic carbon it drove the halogen
into the sulfur: 68% of `bdhi_c5` poses clashed, with
contacts as short as 0.89 Å.

So the ligand handed to gnina must be the **adduct form**: the post-reaction
molecule with the leaving group gone and an implicit hydrogen at the attachment
atom for gnina to replace.

HOW THE LEAVING GROUP IS FOUND. Not by a new hand-written pattern. Every
displacement class already carries a `reactive_atom_smarts` whose first two atoms
are, by construction, the attachment atom and the leaving group's anchor —
`[CH2][Cl]`, `[c]([Cl])[n]`, `[CH2][OX2][SX4](=O)(=O)`. Breaking the bond
between match[0] and match[1] and keeping the fragment that holds match[0]
removes the leaving group whatever its size, from a single chlorine to the
twelve-heavy-atom sulfamate. Reusing the already-validated pattern means the two
cannot drift apart.

MICHAEL ACCEPTORS ARE NOT TRANSFORMED. Nothing leaves in a Michael addition, and
the alkene carbon already carries the implicit hydrogen gnina needs — verified in
the existing poses, where all three Michael classes attached at a carbon that
ended up with four bonds. The true adduct saturates the C=C and adds a hydrogen to the
alpha carbon; the docking input keeps the alkene. That approximation is recorded
in `adduct_approximation` rather than hidden, and it is one hydrogen on a
flexible position, against the twelve atoms it replaces elsewhere.

UNIQUENESS IS A GATE, NOT AN ASSUMPTION. The adduct attachment SMARTS must match
exactly one atom in the product. An R-group that happens to contain its own
acetamide would give a second match, and gnina would bond whichever it found
first — silently docking through the wrong atom. Candidates whose attachment is
ambiguous are rejected here rather than docked.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger

from . import smiles as smi
from . import warhead_library as wl

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)


class AdductError(RuntimeError):
    """A candidate could not be converted to a well-defined adduct form."""


@dataclass(frozen=True)
class AdductForm:
    """The docking-ready form of one candidate."""

    candidate_smiles: str          # what was enumerated (pre-reaction)
    adduct_smiles: str             # what gnina should dock
    attachment_smarts: str         # locates the atom gnina bonds to Cys SG
    attachment_idx: int            # that atom's index in adduct_smiles
    leaving_group_smiles: str | None
    n_atoms_removed: int
    approximation: str | None      # stated where the model is not the true adduct
    degenerate_attachment: str | None = None   # >1 equivalent attack position

    def as_dict(self) -> dict:
        return {
            "adduct_smiles": self.adduct_smiles,
            "adduct_attachment_smarts": self.attachment_smarts,
            "adduct_attachment_idx": self.attachment_idx,
            "leaving_group_smiles": self.leaving_group_smiles,
            "adduct_atoms_removed": self.n_atoms_removed,
            "adduct_approximation": self.approximation,
            "adduct_degenerate_attachment": self.degenerate_attachment,
        }


def _strip_leaving_group(mol: Chem.Mol, reactive_smarts: str) -> tuple[Chem.Mol, str, int]:
    """Break attachment–leaving bond, keep the attachment fragment.

    Returns (adduct_mol, leaving_group_smiles, n_atoms_removed).
    """
    patt = Chem.MolFromSmarts(reactive_smarts)
    if patt is None:
        raise AdductError(f"unparseable reactive_atom_smarts {reactive_smarts!r}")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        raise AdductError(
            f"reactive pattern {reactive_smarts!r} does not match the candidate; "
            "the warhead is absent or was destroyed upstream")
    if len({m[0] for m in matches}) > 1:
        raise AdductError(
            f"reactive pattern {reactive_smarts!r} matches {len(matches)} distinct "
            "attachment atoms; which one reacts is ambiguous")

    attach, leaving = matches[0][0], matches[0][1]
    bond = mol.GetBondBetweenAtoms(attach, leaving)
    if bond is None:
        raise AdductError(
            f"no bond between match atoms {attach} and {leaving} — the SMARTS "
            "does not have the attachment atom and the leaving-group anchor "
            "adjacent, which this transform requires")

    frags = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
    pieces = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=False)
    keep, gone = None, None
    for piece, idxs in zip(pieces, Chem.GetMolFrags(frags)):
        (keep, gone) = (piece, gone) if attach in idxs else (keep, piece)
    if keep is None:
        raise AdductError("fragmentation lost the attachment atom")

    try:
        Chem.SanitizeMol(keep)
    except Exception as exc:  # noqa: BLE001
        raise AdductError(f"adduct fragment does not sanitize: {exc}") from exc

    lg = None
    if gone is not None:
        try:
            Chem.SanitizeMol(gone)
            lg = Chem.MolToSmiles(gone)
        except Exception:  # noqa: BLE001 - the leaving group is reported, not used
            lg = None
    removed = mol.GetNumAtoms() - keep.GetNumAtoms()
    return keep, lg, removed


def to_adduct_form(candidate_smiles: str, warhead_class: str,
                   library=None) -> AdductForm:
    """Build the docking-ready adduct form of one candidate.

    Raises
    ------
    AdductError
        If the warhead is absent, the leaving group cannot be identified, the
        product does not sanitize, or the attachment atom is ambiguous.
    """
    df = library if library is not None else wl.load()
    rows = df[df["class_id"] == warhead_class]
    if rows.empty:
        raise AdductError(f"unknown warhead class {warhead_class!r}")
    row = rows.iloc[0]

    mol = smi.to_mol(candidate_smiles)
    if mol is None:
        raise AdductError(f"unparseable candidate SMILES {candidate_smiles!r}")

    approximation = None
    if bool(row.get("has_leaving_group", True)):
        adduct, lg, removed = _strip_leaving_group(
            mol, str(row["reactive_atom_smarts"]))
    else:
        # Michael addition: nothing leaves and the alkene carbon already carries
        # the implicit hydrogen gnina replaces.
        adduct, lg, removed = mol, None, 0
        approximation = (
            "Michael acceptor docked as the alkene. The true adduct saturates "
            "the C=C and gains one hydrogen on the alpha carbon; that hydrogen "
            "is not modelled in the docked pose.")

    adduct_smiles = smi.canonical(Chem.MolToSmiles(adduct))
    if adduct_smiles is None:
        raise AdductError("adduct does not canonicalize")

    att_smarts = str(row["adduct_attachment_smarts"])
    patt = Chem.MolFromSmarts(att_smarts)
    if patt is None:
        raise AdductError(f"unparseable adduct_attachment_smarts {att_smarts!r}")
    final = smi.to_mol(adduct_smiles)
    matches = final.GetSubstructMatches(patt)
    attach_atoms = {m[0] for m in matches}
    if not attach_atoms:
        raise AdductError(
            f"adduct attachment pattern {att_smarts!r} does not match the "
            f"product {adduct_smiles!r}; gnina would find no reactive atom and "
            "the dock would silently produce nothing")
    degenerate = None
    if len(attach_atoms) > 1:
        # Two different reasons a pattern can match twice, and they need
        # opposite treatments.
        #
        # REAL DEGENERACY: an unsubstituted quinone ring offers C2 and C3, both
        # genuine Michael acceptors related by the ring's near-symmetry. Picking
        # one is chemically legitimate; refusing to pick would discard a whole
        # class over a distinction with no consequence.
        #
        # SPURIOUS AMBIGUITY: an R-group that happens to contain its own
        # matching motif — a pyrimidine on the R-group matched
        # `snar_chloroazine`'s attachment pattern in 30 of 198 candidates. Here
        # the two atoms are in completely different parts of the molecule and
        # gnina would bond whichever it found first, docking through an atom
        # nobody chose.
        #
        # The discriminator is whether the candidates share a ring: warhead
        # attack positions do, an R-group coincidence does not.
        ring_info = final.GetRingInfo()
        shared_ring = any(attach_atoms <= set(ring) for ring in ring_info.AtomRings())
        if not shared_ring:
            raise AdductError(
                f"adduct attachment pattern {att_smarts!r} matches "
                f"{len(attach_atoms)} atoms in different parts of "
                f"{adduct_smiles!r}; gnina bonds the first it finds, so this "
                "candidate would dock through an undetermined atom")
        degenerate = (
            f"{len(attach_atoms)} equivalent attack positions in one ring "
            f"(atoms {sorted(attach_atoms)}); the lowest-indexed was used")
        log.debug("degenerate attachment in %s: %s", adduct_smiles, degenerate)

    idx = min(attach_atoms)
    atom = final.GetAtomWithIdx(idx)
    if atom.GetTotalNumHs() < 1:
        raise AdductError(
            f"attachment atom {idx} ({atom.GetSymbol()}) carries no implicit "
            "hydrogen; gnina replaces a hydrogen to form the bond, so bonding "
            "here would again produce an over-valent atom")

    return AdductForm(
        candidate_smiles=candidate_smiles,
        adduct_smiles=adduct_smiles,
        attachment_smarts=att_smarts,
        attachment_idx=idx,
        leaving_group_smiles=lg,
        n_atoms_removed=removed,
        approximation=approximation,
        degenerate_attachment=degenerate,
    )


def adduct_attachment_smarts(warhead_class: str, library=None) -> str:
    """The SMARTS gnina should be given for this class's adduct form."""
    df = library if library is not None else wl.load()
    rows = df[df["class_id"] == warhead_class]
    if rows.empty:
        raise AdductError(f"unknown warhead class {warhead_class!r}")
    return str(rows.iloc[0]["adduct_attachment_smarts"])
