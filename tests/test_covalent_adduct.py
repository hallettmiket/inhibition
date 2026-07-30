"""
Purpose: Tests for the pre-reaction -> adduct transform that feeds covalent docking.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the shipped warhead library + representative candidate SMILES
Output: pytest pass/fail

D0022 cost a full 1,683-ligand docking run because nothing checked what was
actually handed to gnina. These tests pin the properties that make the adduct
form correct, and the two that would have caught the original defect outright:
the leaving group must be gone, and the attachment atom must have a hydrogen for
gnina to replace.

`test_sn2_acetamides_converge` is the load-bearing one. Chloroacetamide,
sulfamate and sulfonate are SN2 at the same CH2 and differ only in what leaves,
so their adducts are one molecule. If that ever stops being true, either the
transform or the library has drifted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from rdkit import Chem, RDLogger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_adduct as ca          # noqa: E402
from shared import warhead_library as wl          # noqa: E402

RDLogger.DisableLog("rdApp.*")

# One representative candidate per class, all on the same phenyl R-group so the
# comparisons below are like-for-like.
CANDIDATES = {
    "chloroacetamide":
        "O=C(CCl)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1",
    "sulfamate_acetamide":
        "O=C(COS(=O)(=O)Nc1ccccc1)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1",
    "sulfonate_acetamide":
        "CS(=O)(=O)OCC(=O)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1",
    "snar_chloroazine":
        "O=[N+]([O-])c1cnc(Cl)nc1N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1",
    "acrylamide":
        "C=CC(=O)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1",
    "bdhi_c5":
        "O=S1(=O)CC[C@@H](N(c2ccccc2)C2CC(Br)=NO2)C1",
    "naphthoquinone_c2":
        "O=C1c2ccccc2C(=O)C(c2ccccc2)=C1",
    "naphthoquinone_benzo":
        "O=C1C=CC(=O)c2cc(-c3ccccc3)ccc21",
}

HALOGENS = {"F", "Cl", "Br", "I"}


@pytest.fixture(scope="module")
def library():
    return wl.load()


@pytest.mark.parametrize("cls", sorted(CANDIDATES))
def test_adduct_is_a_valid_molecule(library, cls):
    a = ca.to_adduct_form(CANDIDATES[cls], cls, library=library)
    m = Chem.MolFromSmiles(a.adduct_smiles)
    assert m is not None, f"{cls}: adduct {a.adduct_smiles!r} does not sanitize"
    assert not any(x.GetAtomicNum() == 0 for x in m.GetAtoms())


@pytest.mark.parametrize("cls", sorted(CANDIDATES))
def test_attachment_atom_has_a_hydrogen_to_give_up(library, cls):
    """gnina forms the bond by replacing an implicit H — there must be one."""
    a = ca.to_adduct_form(CANDIDATES[cls], cls, library=library)
    m = Chem.MolFromSmiles(a.adduct_smiles)
    atom = m.GetAtomWithIdx(a.attachment_idx)
    assert atom.GetTotalNumHs() >= 1, (
        f"{cls}: attachment atom {atom.GetSymbol()} has no hydrogen, so bonding "
        "Cys SG here would over-fill its valence")


@pytest.mark.parametrize("cls", ["chloroacetamide", "snar_chloroazine", "bdhi_c5"])
def test_halogen_leaving_group_is_removed(library, cls):
    """The defect that started D0022: docking a ligand that kept its halogen."""
    pre = Chem.MolFromSmiles(CANDIDATES[cls])
    assert any(x.GetSymbol() in HALOGENS for x in pre.GetAtoms()), "fixture has no halogen"
    a = ca.to_adduct_form(CANDIDATES[cls], cls, library=library)
    post = Chem.MolFromSmiles(a.adduct_smiles)
    assert not any(x.GetSymbol() in HALOGENS for x in post.GetAtoms()), (
        f"{cls}: adduct {a.adduct_smiles!r} still carries its leaving group")
    assert a.n_atoms_removed == 1


def test_sulfamate_sheds_its_whole_leaving_group(library):
    """Eleven heavy atoms, not one — the transform must not stop at the first bond."""
    a = ca.to_adduct_form(CANDIDATES["sulfamate_acetamide"],
                          "sulfamate_acetamide", library=library)
    assert a.n_atoms_removed == 11, (
        f"expected the whole OS(=O)(=O)NPh to leave, removed {a.n_atoms_removed}")
    assert "S(=O)(=O)N" not in a.adduct_smiles


def test_sn2_acetamides_converge(library):
    """The three SN2 acetamides differ only in what leaves, so share one adduct."""
    got = {cls: ca.to_adduct_form(CANDIDATES[cls], cls, library=library).adduct_smiles
           for cls in ("chloroacetamide", "sulfamate_acetamide", "sulfonate_acetamide")}
    assert len(set(got.values())) == 1, (
        "the three SN2 acetamides must give one identical adduct; got "
        + "; ".join(f"{k}={v}" for k, v in got.items()))


def test_acrylamide_adduct_is_saturated(library):
    """D0030. Acrylamide's adduct is the SATURATED thioether, not the alkene.

    Leaving the C=C in place hands gnina an sp2 carbon and produces a vinyl
    thioether — planar and rigid where the real linker rotates freely. T_3's
    only warhead is acrylamide, so this was about to bias its entire search
    against decorations that need to bend.
    """
    a = ca.to_adduct_form(CANDIDATES["acrylamide"], "acrylamide", library=library)
    assert a.n_atoms_removed == 0, "nothing leaves in a Michael addition"
    assert a.leaving_group_smiles is None

    m = Chem.MolFromSmiles(a.adduct_smiles)
    assert m is not None
    assert not m.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[CX3][CX3]=O")), (
        "the acceptor alkene survived the transform; the docked ligand would be "
        "a vinyl thioether")
    # The product is a propanamide whose terminal carbon carries the hydrogen
    # gnina replaces — exactly parallel to the acetamide classes' CH3.
    assert m.HasSubstructMatch(Chem.MolFromSmarts("[CH3][CH2][CX3](=O)[NX3]"))
    att = m.GetAtomWithIdx(a.attachment_idx)
    assert att.GetTotalNumHs() == 3, "attachment carbon should be a CH3"
    assert str(att.GetHybridization()) == "SP3"


def test_quinone_adduct_keeps_its_alkene(library):
    """D0030, the other half: quinones re-aromatize, so sp2 is CORRECT.

    Thiol addition to a 1,4-naphthoquinone gives a hydroquinone that re-oxidizes
    to the 2-thio-quinone. Saturating it here would model a species that does
    not persist — the opposite error to the acrylamide one.
    """
    for cls in ("naphthoquinone_c2", "naphthoquinone_benzo"):
        a = ca.to_adduct_form(CANDIDATES[cls], cls, library=library)
        m = Chem.MolFromSmiles(a.adduct_smiles)
        assert m.GetAtomWithIdx(a.attachment_idx).GetTotalNumHs() == 1
        assert str(m.GetAtomWithIdx(a.attachment_idx).GetHybridization()) == "SP2"
        assert a.approximation and "re-aromatized" in a.approximation


def test_saturation_is_driven_by_the_library_not_the_mechanism(library):
    """`michael_addition` alone must never decide whether to saturate.

    Both acrylamide and the quinones are `michael_addition`; they need opposite
    treatments. The discriminator is re-aromatizability, which is chemistry a
    SMARTS cannot infer, so it lives in the library as a declared column.
    """
    michael = library[library.mechanism == "michael_addition"]
    # NOT a count. This asserted len(michael) == 3 and failed the moment a
    # legitimate fourth Michael acceptor (cinnamamide) was added -- a passing
    # test became a blocker for correct work. What the test is FOR is that
    # `mechanism` does not determine saturation, so assert that directly: the
    # same mechanism must carry both answers.
    assert len(michael) >= 3, "need several michael_addition classes to test this"
    assert set(michael["adduct_saturates_alkene"]) == {True, False}, (
        "michael_addition classes must disagree on saturation, or the column is "
        "redundant with mechanism and this guard proves nothing")
    saturating = set(michael.loc[michael.adduct_saturates_alkene, "class_id"])
    assert saturating == {"acrylamide", "cinnamamide"}, (
        "acrylamide and cinnamamide saturate — no re-aromatization is available; "
        "the quinones do not, because the hydroquinone re-oxidizes")


def test_ambiguous_attachment_from_the_rgroup_is_rejected(library):
    """An R-group carrying its own pyrimidine must not dock through a coin flip.

    Thirty of 198 `snar_chloroazine` candidates hit this with the first,
    looser attachment SMARTS. The tightened pattern resolves it; if it ever
    loosens again, this fails rather than silently docking the wrong atom.
    """
    tricky = "O=[N+]([O-])c1cnc(Cl)nc1N(c1ncncc1)[C@@H]1CCS(=O)(=O)C1"
    a = ca.to_adduct_form(tricky, "snar_chloroazine", library=library)
    m = Chem.MolFromSmiles(a.adduct_smiles)
    patt = Chem.MolFromSmarts(a.attachment_smarts)
    assert len({x[0] for x in m.GetSubstructMatches(patt)}) == 1


def test_degenerate_ring_attachment_is_allowed_but_recorded(library):
    """Quinone C2/C3 are both genuine acceptors — pick one, say so."""
    benzo = "O=C1C=CC(=O)c2cc(N(c3ccccc3)[C@@H]3CCS(=O)(=O)C3)ccc21"
    a = ca.to_adduct_form(benzo, "naphthoquinone_benzo", library=library)
    assert a.degenerate_attachment is not None
    assert "one ring" in a.degenerate_attachment


def test_missing_warhead_is_an_error(library):
    """A molecule with no warhead must not silently produce an adduct."""
    with pytest.raises(ca.AdductError):
        ca.to_adduct_form("c1ccccc1C(=O)NC", "chloroacetamide", library=library)


def test_unknown_class_is_an_error(library):
    with pytest.raises(ca.AdductError):
        ca.to_adduct_form(CANDIDATES["chloroacetamide"], "not_a_class", library=library)


def test_protocol_uses_the_adduct_smarts_not_the_reactive_one(library):
    """Regression for the root cause: the protocol must not hand gnina a pattern
    that names a leaving group the docked ligand no longer has."""
    from shared import covalent_protocol as cp
    smarts = cp.load_warhead_smarts()
    for cls, s in smarts.items():
        row = library[library.class_id == cls].iloc[0]
        assert s == row["adduct_attachment_smarts"], f"{cls} using the wrong SMARTS"
        if row["has_leaving_group"]:
            assert s != row["reactive_atom_smarts"]
