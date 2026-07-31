"""
Purpose: The pocket size ceiling must never reject a known Pin1 binder.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: data/reference/pin1_reference_binders_3.csv
Output: pass/fail

A generation filter that removes a molecule known to bind the target is not a
filter, it is a bug that looks like tidiness. This is the test that sets the
ceiling: raise it if a future reference binder exceeds it, never lower it to
make a generated set look cleaner.
"""

from __future__ import annotations

import pandas as pd
import pytest
from rdkit import Chem

from shared import pocket_size as ps
from shared import reference_set as rs

# Peptidic macrocycles are a different modality; none of T_1..T_4 generates
# them, and they are excluded from the gate for want of property-matched
# decoys. They are exempt from the ceiling by design, named individually so
# the exemption cannot silently widen.
PEPTIDIC = {
    "Wildemann-macrocyclic-peptide",
    "Liu-Pei-cyclic-peptide",
    "Jiang-Pei-bicyclic-CPP-peptide",
}


@pytest.fixture(scope="module")
def binders() -> pd.DataFrame:
    return rs.load().master


def test_no_known_nonpeptidic_binder_is_rejected(binders):
    rejected = []
    for _, r in binders.iterrows():
        if r["name"] in PEPTIDIC:
            continue
        m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
        if m is None:
            continue
        if not ps.fits_pocket(str(r["canonical_smiles"])):
            rejected.append((r["name"], m.GetNumHeavyAtoms()))
    assert not rejected, (
        f"the ceiling of {ps.MAX_HEAVY_ATOMS} rejects known binder(s): "
        f"{rejected}. Raise MAX_HEAVY_ATOMS — never keep a filter that removes "
        "a molecule the target is known to bind.")


def test_the_peptides_are_above_the_ceiling(binders):
    """The exemption must be doing real work, not covering for a loose limit."""
    for name in PEPTIDIC:
        row = binders[binders["name"] == name]
        if row.empty:
            continue
        m = Chem.MolFromSmiles(str(row.iloc[0]["canonical_smiles"]))
        assert m.GetNumHeavyAtoms() > ps.MAX_HEAVY_ATOMS, (
            f"{name} is below the ceiling, so exempting it is meaningless")


def test_ceiling_is_generous_relative_to_the_structural_estimate():
    """The brief was a loose limit, and 'loose' should be checkable."""
    assert ps.MAX_HEAVY_ATOMS > ps.HEAVY_ATOMS_TIGHT_PACKING
    assert ps.MAX_HEAVY_ATOMS >= 1.5 * ps.HEAVY_ATOMS_TIGHT_PACKING
    assert ps.MAX_HEAVY_ATOMS > ps.LARGEST_KNOWN_NONPEPTIDIC_HEAVY_ATOMS


def test_obviously_oversized_molecules_are_pruned():
    huge = "C" * 200
    assert not ps.fits_pocket(huge)
    assert "heavy atoms" in (ps.why(huge) or "")


def test_a_real_binder_passes_and_reports_no_reason():
    sulfopin = "O=C(CCl)N(Cc1ccccc1)C1CCS(=O)(=O)C1"
    assert ps.fits_pocket(sulfopin)
    assert ps.why(sulfopin) is None


def test_unparseable_is_pruned_not_crashed():
    assert not ps.fits_pocket("not a molecule")
    assert ps.why("not a molecule") == "unparseable"
