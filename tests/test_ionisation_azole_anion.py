"""obabel writes azole anions RDKit cannot read, and `protonate` never checked.

THE DEFECT. At pH 7.4 a tetrazole (pKa ~4.9) is deprotonated, and obabel gets
that right. It then serialises the anion it just made in a Kekule form that puts
three bonds on the anionic nitrogen -- `[N-]1=C(NN=N1)R` -- which RDKit rejects
with "Explicit valence for atom # 0 N, 3, is greater than permitted".

It does this only to an anion it CREATES. Handed the aromatic anion
`c1nnn[n-]1` it round-trips it perfectly. So the bug is invisible to every
molecule that is neutral at the azole, which is every molecule screened to date
-- and it surfaced the first time an ionisable heteroaryl went through, as a
failure two stages downstream in `nac_screen.prepare_ligand` that read like a
property of the molecule.

WHAT `protonate` GUARANTEED, AND WHAT IT DID NOT. It guaranteed IDENTITY: the
right string matched back to the right id, with the recursive-split machinery
above it built precisely so a short return could not slide results onto the
wrong candidates. It said nothing about whether the string was a MOLECULE. That
gap is this project's shape #4 -- a guard scoped to one property while the
neighbouring one goes unchecked.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import ionisation as ion            # noqa: E402

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem                          # noqa: E402

#: Exactly what obabel emitted for issue #81's molecule.
OBABEL_BAD = "O=C(C=Cc1ccccc1CCC1=[N-]N=NN1)n1c(c2ccccc2)cc2ccc(Br)cc12"
#: The species it should have written: the aromatic tetrazolate.
EXPECTED = "O=C(C=Cc1ccccc1CCc1nnn[n-]1)n1c(-c2ccccc2)cc2ccc(Br)cc21"


def test_the_input_really_is_unparseable():
    """The premise. If RDKit ever accepts this, the repair is dead code."""
    assert Chem.MolFromSmiles(OBABEL_BAD) is None


def test_repair_produces_the_aromatic_tetrazolate():
    got = ion._repair_azole_anion(OBABEL_BAD)
    assert got is not None
    assert Chem.CanonSmiles(got) == Chem.CanonSmiles(EXPECTED)


def test_repair_preserves_obabels_charge():
    """The repair fixes VALENCE. Changing the protonation state would be D0074."""
    m = Chem.MolFromSmiles(OBABEL_BAD, sanitize=False)
    m.UpdatePropertyCache(strict=False)
    before = sum(a.GetFormalCharge() for a in m.GetAtoms())
    after = Chem.GetFormalCharge(Chem.MolFromSmiles(ion._repair_azole_anion(OBABEL_BAD)))
    assert before == after == -1


def test_repair_preserves_the_heavy_atom_skeleton():
    """A repair that edited the molecule would be worse than the failure."""
    bad = Chem.MolFromSmiles(OBABEL_BAD, sanitize=False)
    bad.UpdatePropertyCache(strict=False)
    fixed = Chem.MolFromSmiles(ion._repair_azole_anion(OBABEL_BAD))
    assert bad.GetNumAtoms() == fixed.GetNumAtoms()
    assert (sorted(a.GetSymbol() for a in bad.GetAtoms())
            == sorted(a.GetSymbol() for a in fixed.GetAtoms()))


@pytest.mark.parametrize("smi", [
    "CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl",     # Sulfopin -- neutral, fine
    "c1nnn[n-]1",                              # already correct
    "CC(=O)[O-]",                              # an anion that is not this shape
])
def test_repair_declines_anything_that_is_not_this_defect(smi):
    """It must not be a general-purpose molecule rewriter."""
    assert ion._repair_azole_anion(smi) is None


def test_validated_passes_good_smiles_through_untouched():
    good = "CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl"
    assert ion._validated("x", good) == good


def test_validated_drops_what_it_cannot_repair(caplog):
    """ABSENT, never substituted. The caller's contract is that a missing id
    means no species could be built -- docking the neutral form in its place is
    exactly the substitution D0074 exists to prevent."""
    assert ion._validated("x", "not a smiles at all [[[") is None


def test_protonate_returns_a_parseable_species_for_the_issue81_control():
    """End to end, through the same call `nac_screen.prepare_ligand` makes."""
    neutral = "O=C(C=Cc1ccccc1CCc1nnn[nH]1)n1c(-c2ccccc2)cc2ccc(Br)cc21"
    got = ion.protonate({"ian_ctrl_issue81": neutral})
    assert "ian_ctrl_issue81" in got, "the control was dropped again"
    m = Chem.MolFromSmiles(got["ian_ctrl_issue81"])
    assert m is not None
    assert Chem.GetFormalCharge(m) == -1, "the tetrazolate is the pH 7.4 species"


def test_the_michael_acceptor_survives_the_repair():
    """The warhead is what the screen measures; the repair must not disturb it.

    `prepare_ligand` raises if the reactive SMARTS stops matching, so a repair
    that perturbed the enone would surface as 'this molecule has no reactive
    centre' rather than as a repair bug.
    """
    m = Chem.MolFromSmiles(ion._repair_azole_anion(OBABEL_BAD))
    patt = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=O")
    assert len(m.GetSubstructMatches(patt)) == 1
