"""
Purpose: Validation suite for the hand-authored warhead-class library.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: data/reference/warhead_classes_3.csv, config/seeds.yaml
Output: pytest pass/fail

WHY THIS FILE EXISTS. The warhead library is hand-written chemistry data, and it
has produced three separate defects, each caught downstream after wasting
compute rather than at the source:

1. `bdhi_*` was written WITHOUT its bromine — a ring-opening warhead with no
   leaving group.
2. `snar_chloroazine` was written WITHOUT its chlorine. The 5b validity gate
   caught it only after 198 products had been enumerated (0/198 intact).
3. `sulfamate_acetamide` carried TWO attachment points, so all 198 of its
   products kept a dangling dummy atom. It passed core verification, the alert
   gate AND the warhead gate, and reached covalent docking as a non-molecule.

Every defect was a property of one CSV row and testable in milliseconds. These
tests assert those properties directly, and each historical defect gets an
explicit regression test that feeds the broken value back in and demands the
guard fire.

The tests run against the SHIPPED library, not a fixture — a fixture would pass
happily while the real data rotted.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml
from rdkit import Chem, RDLogger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import smiles as smi                    # noqa: E402
from shared import warhead_library as wl            # noqa: E402

RDLogger.DisableLog("rdApp.*")

# Substructure that MUST be present for each declared leaving group. Michael
# acceptors add rather than displace, so they legitimately have none.
LEAVING_GROUP_SMARTS = {
    "chloride": "[Cl]",
    "bromide": "[Br]",
    "iodide": "[I]",
    "fluoride": "[F]",
    "sulfamate": "[OX2][SX4](=O)(=O)[NX3]",
    "sulfonate": "[OX2][SX4](=O)(=O)[CX4]",
    "none_addition": None,
}

CORE_SMILES = yaml.safe_load(
    (REPO / "config" / "seeds.yaml").read_text(encoding="utf-8")
)["seeds"]["sulfopin"]["protected_core"]["smiles"]

# A neutral R-group stand-in, so coupling tests exercise the warhead only.
DUMMY_RGROUP = "[*]c1ccccc1"


@pytest.fixture(scope="module")
def library():
    return wl.load()


@pytest.fixture(scope="module")
def rows(library):
    return [r for _, r in library.iterrows()
            if str(r["warhead_fragment_smiles"]) != wl.UNRESOLVED]


def _couple(warhead_frag: str, rgroup_frag: str = DUMMY_RGROUP) -> str | None:
    """Mirror of 01_enumerate.couple, imported by path (the module is numbered)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "t4_enumerate", REPO / "approaches" / "t4_combinatorial" / "01_enumerate.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.couple(CORE_SMILES, warhead_frag, rgroup_frag)


def test_library_loads(library):
    """The shipped library passes its own loader."""
    assert len(library) > 0
    assert library["class_id"].is_unique


def test_every_fragment_has_exactly_one_attachment_point(rows):
    """Regression: sulfamate_acetamide shipped with two, producing 198 non-molecules."""
    for r in rows:
        n = str(r["warhead_fragment_smiles"]).count("*")
        assert n == 1, (
            f"{r['class_id']}: {n} attachment points. The core coupling fills "
            "exactly one; any other survives as a dangling dummy atom.")


def test_every_fragment_sanitizes(rows):
    """A fragment that only parses unsanitized will fail later, in docking."""
    for r in rows:
        frag = str(r["warhead_fragment_smiles"])
        m = Chem.MolFromSmiles(frag)
        assert m is not None, f"{r['class_id']}: {frag!r} does not sanitize"


def test_reactive_atom_smarts_matches_its_own_fragment(rows):
    """The reactive atom the docking protocol targets must exist in the fragment.

    This is the check that would have caught both missing-halogen defects at
    source: each SMARTS names the halogen, so a fragment without it cannot
    match. gnina was instead asked to constrain an atom pattern that was not
    there, and 168 of 192 docks failed with no useful message.
    """
    for r in rows:
        frag = str(r["warhead_fragment_smiles"])
        patt = Chem.MolFromSmarts(str(r["reactive_atom_smarts"]))
        assert patt is not None, f"{r['class_id']}: unparseable reactive_atom_smarts"
        m = Chem.MolFromSmiles(frag)
        assert m is not None and m.HasSubstructMatch(patt), (
            f"{r['class_id']}: reactive_atom_smarts "
            f"{r['reactive_atom_smarts']!r} does not match its own fragment "
            f"{frag!r} — the warhead is mis-specified, not merely unusual.")


def test_declared_leaving_group_is_present(rows):
    """Regression: bdhi_* shipped without Br, snar_chloroazine without Cl."""
    for r in rows:
        lg = str(r.get("leaving_group", "")).strip()
        assert lg in LEAVING_GROUP_SMARTS, (
            f"{r['class_id']}: unknown leaving_group {lg!r}; add it to "
            "LEAVING_GROUP_SMARTS with the substructure it implies.")
        smarts = LEAVING_GROUP_SMARTS[lg]
        if smarts is None:
            continue
        m = Chem.MolFromSmiles(str(r["warhead_fragment_smiles"]))
        patt = Chem.MolFromSmarts(smarts)
        assert m.HasSubstructMatch(patt), (
            f"{r['class_id']}: declares leaving_group={lg!r} but the fragment "
            f"{r['warhead_fragment_smiles']!r} contains no such group. A "
            "displacement warhead with nothing to displace is not a warhead.")


def test_coupling_to_core_yields_a_complete_molecule(rows):
    """Every class must couple to the sulfopin core and leave no dummy atoms.

    Core verification alone is not enough: the malformed sulfamate contained the
    core perfectly and was still not a molecule.
    """
    for r in rows:
        product = _couple(str(r["warhead_fragment_smiles"]))
        assert product is not None, f"{r['class_id']}: coupling to the core failed"
        m = Chem.MolFromSmiles(product)
        assert m is not None, f"{r['class_id']}: product {product!r} does not sanitize"
        dummies = [a for a in m.GetAtoms() if a.GetAtomicNum() == 0]
        assert not dummies, (
            f"{r['class_id']}: product carries {len(dummies)} dangling dummy "
            f"atom(s): {product}")
        assert smi.candidate_id(product, prefix="t4") is not None, (
            f"{r['class_id']}: product has no InChIKey, so its candidate_id "
            "would collide with every other keyless product.")


def test_warhead_survives_coupling(rows):
    """Gate 5b as a unit test: is it still an electrophile once attached?"""
    for r in rows:
        product = _couple(str(r["warhead_fragment_smiles"]))
        patt = Chem.MolFromSmarts(str(r["reactive_atom_smarts"]))
        m = Chem.MolFromSmiles(product)
        assert m.HasSubstructMatch(patt), (
            f"{r['class_id']}: the reactive atom is destroyed by coupling to the "
            "core — enumerating this class would produce 198 dead products.")


def test_statuses_are_known(library):
    assert set(library["structure_status"]) <= set(wl.STATUS_TIERS)


# --- Regression tests: feed the historical defects back in ------------------

def test_two_attachment_points_is_rejected(tmp_path, library):
    """The sulfamate defect, verbatim."""
    broken = library.copy()
    broken.loc[broken["class_id"] == "sulfamate_acetamide",
               "warhead_fragment_smiles"] = "[*]C(=O)COS(=O)(=O)N[*]"
    p = tmp_path / "broken_two_attachments.csv"
    broken.to_csv(p, index=False)
    with pytest.raises(wl.WarheadLibraryError, match="attachment points"):
        wl.load(p)


def test_zero_attachment_points_is_rejected(tmp_path, library):
    broken = library.copy()
    broken.loc[broken["class_id"] == "chloroacetamide",
               "warhead_fragment_smiles"] = "CC(=O)CCl"
    p = tmp_path / "broken_no_attachment.csv"
    broken.to_csv(p, index=False)
    with pytest.raises(wl.WarheadLibraryError, match="no attachment point"):
        wl.load(p)


@pytest.mark.parametrize("class_id,broken_smiles,defect", [
    ("snar_chloroazine", "[*]c1ncc([N+](=O)[O-])cn1", "chlorine dropped"),
    ("bdhi_c5", "[*]C1=NOCC1", "bromine dropped"),
])
def test_missing_leaving_group_is_caught(library, class_id, broken_smiles, defect):
    """Both halogen defects, as they were actually written.

    These reach `load()` cleanly — they parse and have one attachment point — so
    only the leaving-group and reactive-atom checks catch them. That is why
    those checks exist.
    """
    row = library[library["class_id"] == class_id].iloc[0]
    m = Chem.MolFromSmiles(broken_smiles)
    patt = Chem.MolFromSmarts(str(row["reactive_atom_smarts"]))
    lg = LEAVING_GROUP_SMARTS[str(row["leaving_group"])]
    assert not m.HasSubstructMatch(patt), f"{class_id}: {defect} should not match"
    assert not m.HasSubstructMatch(Chem.MolFromSmarts(lg)), (
        f"{class_id}: {defect} should leave no leaving group")


def test_unknown_status_is_rejected(tmp_path, library):
    broken = library.copy()
    broken.loc[0, "structure_status"] = "PROBABLY_FINE"
    p = tmp_path / "broken_status.csv"
    broken.to_csv(p, index=False)
    with pytest.raises(wl.WarheadLibraryError, match="unknown structure_status"):
        wl.load(p)
