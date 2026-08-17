"""Building the pH 7.4 species: the +-1 ceiling that deleted every dication.

`stamp_identity.protonate()` applies the charge the frame already recorded in
`charge_ph74`. It used to protonate exactly ONE site -- `hits[0][0]` -- and then
test `GetFormalCharge(out) == want`; worse, the working molecule was rebuilt from
the input INSIDE the pattern loop, so nothing accumulated and each pattern
produced a fresh +-1.

A molecule needing +2 therefore made +1, failed the equality, tried the next
pattern, made +1 again, and returned None -- however many basic sites it had. It
was never a chemistry failure and never a hard case: it was arithmetic.

WHAT IT COST. Every dication in the library was stamped `docked_species_ok =
False` and dropped from screening: all 60 of D4's failures and 7 of D3's. The 60
are ALL BDHI (30 bdhi_c4, 30 bdhi_c5, zero acrylamide), so both BDHI arms entered
the screen 15% short against a full acrylamide arm -- and the only symptom would
have been BDHI appearing to underperform. @tt8804: "we have gone over this
goddamn protonation issue so many times".
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))


def _si():
    spec = importlib.util.spec_from_file_location(
        "stamp_identity", REPO / "scripts" / "stamp_identity.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _charge(smi):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smi)
    return None if m is None else Chem.GetFormalCharge(m)


#: A real D4 bdhi_c4 row that the ceiling deleted: a tertiary amine plus an
#: N-aryl piperazine, recorded charge +2.
DICATION = "O=S1(=O)CC[C@@H](N(Cc2ccc(N3CCNCC3)cc2)C2CC(Br)=NO2)C1"


def test_a_dication_is_built_not_refused():
    """THE REGRESSION, on one of the 60."""
    out = _si().protonate(DICATION, 2)
    assert out is not None, "returned None for a molecule with two basic sites"
    assert _charge(out) == 2


def test_the_protons_land_on_the_basic_nitrogens_not_the_anilinic_one():
    """A piperazine's two N are not equivalent: the aryl-attached N is anilinic
    (pKa ~4-5, neutral at 7.4) and only the distal amine is basic. Protonating
    the anilinic N would invent a species that does not exist at this pH -- the
    `!$(N-a)` exclusion in CATION_SITES is what prevents it, and it must survive
    the change to accumulate across sites."""
    out = _si().protonate(DICATION, 2)
    # the aryl-attached ring N stays neutral: an aromatic carbon bonded to a
    # neutral N3, not to [N+]
    assert "cc2)C2CC" in out or "c2ccc(" in out
    assert out.count("[NH2+]") + out.count("[NH+]") == 2


@pytest.mark.parametrize("smi,want", [
    ("CCN", 1),                       # single aliphatic amine
    ("NCCN", 2),                       # diamine -- needs accumulation
    ("NCCNCCN", 3),                    # triamine
    ("CC(=O)O", -1),                   # carboxylic acid
    ("OC(=O)CCC(=O)O", -2),            # diacid -- needs accumulation
])
def test_charges_beyond_one_are_reachable_in_both_directions(smi, want):
    out = _si().protonate(smi, want)
    assert out is not None, f"{smi}: could not build {want:+d}"
    assert _charge(out) == want


def test_a_neutral_molecule_is_returned_unchanged():
    out = _si().protonate("c1ccccc1", 0)
    assert out is not None and _charge(out) == 0


def test_a_charge_the_structure_cannot_carry_is_still_refused():
    """STILL AN EQUALITY, NOT A BEST EFFORT. Benzene has no basic site, so a
    recorded +1 is a real disagreement between the ionisation model and the
    structure. Docking it at the wrong charge is the silent substitution this
    script exists to prevent, so it stays stamped rather than being approximated."""
    assert _si().protonate("c1ccccc1", 1) is None


def test_a_missing_charge_is_not_read_as_zero():
    """pandas NA is not falsy. A missing charge means nobody computed one, which
    is a different claim from neutral."""
    import numpy as np
    assert _si().protonate("CCN", float("nan")) is None
    assert _si().protonate("CCN", None) is None


def test_a_molecule_already_carrying_the_charge_is_not_protonated_again():
    """135 of an earlier pass's 234 failures were molecules whose SMILES already
    carried the charge. Protonating one of those would produce a +2."""
    out = _si().protonate("CC[NH3+]", 1)
    assert out is not None and _charge(out) == 1
