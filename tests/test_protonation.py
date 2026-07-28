"""
Purpose: Tests for physiological-pH protonation of ligands before scoring.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: representative candidate SMILES
Output: pytest pass/fail

Generators emit neutral SMILES. At pH 7.4 a carboxylic acid is a carboxylate,
and 14 of T_2's 25 shortlisted candidates carry one against a cationic pocket.
Scoring the neutral form is scoring a species that is not there — the same error
as D0022 and D0030, in a third place.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import protonation as prot          # noqa: E402


def test_carboxylic_acid_becomes_a_carboxylate():
    """The defect that motivated this module, stated directly."""
    d = prot.dominant_state("CC(=O)O")
    assert d["protonated_charge"] == -1
    assert d["charge_changed"]
    assert "[O-]" in d["protonated_smiles"]


def test_aliphatic_amine_becomes_an_ammonium():
    d = prot.dominant_state("CCN")
    assert d["protonated_charge"] == +1
    assert d["charge_changed"]


def test_a_neutral_molecule_is_left_alone():
    d = prot.dominant_state("c1ccccc1")
    assert d["protonated_charge"] == 0
    assert not d["charge_changed"]


def test_atra_the_t2_seed_is_an_anion_at_physiological_ph():
    """T_2's whole neighbourhood descends from a carboxylic acid."""
    atra = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)/C)/C"
    d = prot.dominant_state(atra)
    assert d["protonated_charge"] == -1, (
        "ATRA is a retinoic ACID; if it scores as neutral, so does most of T_2")


def test_borderline_groups_are_flagged_not_guessed():
    """An imidazole near pH 7.4 is a judgement call, and must say so."""
    ok, groups = prot.is_confident("c1cnc[nH]1")
    assert not ok
    assert "imidazole" in groups

    ok, groups = prot.is_confident("CC(=O)O")
    assert ok and not groups


def test_heavy_atom_skeleton_is_unchanged():
    """Protonation moves hydrogens only, so a docked pose stays usable."""
    from rdkit import Chem

    smi = "CC(C=Cc1ccc2cc(C)ccc2n1)=CC=CC(C)=CC(=O)O"
    before = Chem.MolFromSmiles(smi).GetNumAtoms()
    after = Chem.MolFromSmiles(prot.dominant_state(smi)["protonated_smiles"]).GetNumAtoms()
    assert before == after


def test_unparseable_smiles_raises_rather_than_returning_a_guess():
    with pytest.raises(ValueError):
        prot.dominant_state("not a molecule")


def test_a_pose_docked_neutral_still_maps_onto_its_anion():
    """Regression: protonation broke every T_1 and T_2 candidate at once.

    Poses were docked in the neutral form, so a carboxylic acid pose carries its
    -OH hydrogen. Matching it against the pH-7.4 template — a carboxylate whose
    oxygen has no hydrogen and a -1 charge — is a valence contradiction unless
    the pose's hydrogens are stripped first. `removeHs=True` on the RDKit PDB
    reader is a no-op under `sanitize=False`, which is what hid this.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    acid = "CC(=O)O"
    anion = prot.dominant_state(acid)["protonated_smiles"]

    # A pose-like molecule: correct heavy atoms, explicit acidic H, no bond orders.
    posed = Chem.AddHs(Chem.MolFromSmiles(acid))
    AllChem.EmbedMolecule(posed, randomSeed=42)

    rw = Chem.RWMol(posed)
    for idx in sorted([a.GetIdx() for a in posed.GetAtoms()
                       if a.GetAtomicNum() == 1], reverse=True):
        rw.RemoveAtom(idx)
    heavy = rw.GetMol()

    fixed = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(anion), heavy)
    assert Chem.GetFormalCharge(fixed) == -1
