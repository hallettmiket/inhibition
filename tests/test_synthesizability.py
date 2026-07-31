"""
Purpose: A synthesizability rule that rejects a real molecule is wrong.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: the reference binders + hand-written positive and negative cases
Output: pass/fail

THIS TEST FOUND TWO BROKEN RULES BEFORE THEY WERE USED. The first version of
`adjacent_quaternary_ring_carbons` counted a ring carbon's OWN ring bonds as
substituents, so every ordinary fused bicyclic matched -- 47% of T_1, plus
EGCG, a molecule found in tea. A `three_fused_rings_sharing_one_atom` rule
built on `[R3]` rejected PiB, a clinical PET tracer. Both were removed or
corrected.

The cost of a false positive here is invisible and permanent: a good candidate
is dropped and nobody ever learns it existed. So the bar is that no molecule
anyone has actually made may be rejected.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from shared import reference_set as rs
from shared import synthesizability as syn

# Its SMILES is the literal string UNVERIFIED in the reference set, so it is
# unparseable by construction. Named explicitly rather than filtered by a
# predicate, so that if it is ever resolved this test starts checking it.
UNRESOLVED_BINDERS = {"Byun-BDHI-fragment"}


def test_no_known_binder_is_called_unsynthesizable():
    """Every one of these compounds exists. Rejecting one is a rule bug."""
    rejected = []
    for _, r in rs.load().master.iterrows():
        if r["name"] in UNRESOLVED_BINDERS:
            continue
        v = syn.violations(str(r["canonical_smiles"]))
        if v:
            rejected.append((r["name"], [x.name for x in v]))
    assert not rejected, (
        f"rules reject molecules that have been made: {rejected}. The compound "
        "exists, so the rule is wrong — fix or remove the rule, never the "
        "expectation.")


@pytest.mark.parametrize("name,smiles", [
    # Real drugs and natural products spanning the structural motifs the rules
    # brush against: fused polycyclics, quaternary centres, dense heteroatoms.
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("morphine-like fused system", "CN1CCC23c4c5ccc(O)c4OC2C(O)C=CC3C1C5"),
    ("cholesterol core", "CC(C)CCCC(C)C1CCC2C1(C)CCC1C2CC=C2CC(O)CCC12C"),
    ("penicillin core", "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"),
    ("EGCG", "O[C@@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1"),
])
def test_real_molecules_pass(name, smiles):
    assert Chem.MolFromSmiles(smiles) is not None, f"{name}: bad test SMILES"
    assert syn.is_plausible(smiles), f"{name} rejected: {syn.explain(smiles)}"


@pytest.mark.parametrize("smiles,expect", [
    ("CC1(C)C(C)(C)CCCC1", "adjacent_quaternary_ring_carbons"),
    ("CCOOCC", "peroxide_or_higher"),
    ("OC(O)CCC", "geminal_diol_or_hemiketal"),
    ("CNNNC", "nitrogen_nitrogen_nitrogen_chain"),
])
def test_genuinely_bad_motifs_are_caught(smiles, expect):
    names = [r.name for r in syn.violations(smiles)]
    assert expect in names, f"{smiles} should trip {expect}, got {names}"


def test_unparseable_is_reported_not_crashed():
    v = syn.violations("this is not a molecule")
    assert [r.name for r in v] == ["unparseable"]


def test_every_rule_has_a_rationale_and_valid_smarts():
    """A rule a chemist cannot argue with is a rule they cannot correct."""
    for rule in syn.RULES:
        assert rule.why.strip(), f"{rule.name} has no rationale"
        assert Chem.MolFromSmarts(rule.smarts) is not None, \
            f"{rule.name} has unparseable SMARTS {rule.smarts!r}"
