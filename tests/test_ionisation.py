"""
Purpose: charge at pH 7.4 must match what was docked, and must survive obabel's failure modes.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: shared/ionisation.py
Output: pass/fail

WHY THIS EXISTS. #6 item 7 asks for charge-stratified ranking. The obvious
column to stratify on already exists -- `descriptors.formal_charge` -- and it
is **0 for essentially every molecule in the project**, because it is the
charge of the NEUTRAL canonical SMILES while docking protonated for pH 7.4.
Stratifying on it would produce one stratum and look like it worked.

The five T_2 seeds are the sharpest demonstration: all five are
`formal_charge = 0`, and at pH 7.4 they span four distinct charge states.

Two failure modes of the obabel call are pinned here because both are silent:

* **One bad molecule halts the whole stream.** Measured on T_1: obabel
  reported "143 molecules converted" from a 4,803-line file, wrote nothing to
  stderr and exited 0. The 144th is `OC[P@TB14](O)(O)(O)CO` -- pentavalent
  phosphorus with trigonal-bipyramidal stereochemistry, a DiffSBDD artefact.
  97% of the arm vanished with no error.
* **Results must be matched by identity, not position.** A short return with
  positional matching shifts every subsequent charge onto the wrong candidate,
  and nothing about the output looks wrong.
"""

from __future__ import annotations

import pytest

from shared import ionisation as ion

# All five are formal_charge = 0 on the neutral SMILES.
SEEDS = {
    "atra": ("CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C(=O)O)/C)/C", -1),
    "du_xu": ("O=C(N[C@H](C/C=C/c1cccc(F)c1)C(=O)O)c1ccc2ccccc2c1", -1),
    "guo_pfizer": ("O=C(N[C@@H](COP(=O)(O)O)Cc1cccc(F)c1)c1cc2ccccc2s1", -2),
    "potter_astex": ("COC(=O)[C@@H](Cc1nc2ccccc2[nH]1)NC(=O)c1cc(-c2ccc(CN)cc2)oc1", 1),
    "sulfopin": ("CC(C)(C)CN(C1CCS(=O)(=O)C1)C(=O)CCl", 0),
}

#: The molecule that silently halted a 4,803-molecule obabel stream.
POISON = "OC[P@TB14](O)(O)(O)CO"


def test_the_seeds_span_four_charge_states_where_formal_charge_reports_one():
    """The whole reason this module exists."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    got = ion.charge_at_ph({k: v[0] for k, v in SEEDS.items()})
    for name, (smiles, expected) in SEEDS.items():
        assert got[name] == expected, f"{name}: {got[name]} != {expected}"
        assert Chem.GetFormalCharge(Chem.MolFromSmiles(smiles)) == 0, (
            f"{name} is no longer neutral as drawn; this fixture's premise is "
            "gone and the test no longer demonstrates anything")
    assert len(set(got.values())) == 4, "the seeds should span four states"


def test_one_poison_molecule_does_not_take_the_batch_with_it():
    """The measured T_1 failure, as a fixture.

    Without the split-retry, obabel converts everything up to the poison and
    silently stops -- so the molecules AFTER it come back missing.
    """
    good = {f"good_{i}": s for i, (s, _) in enumerate(SEEDS.values())}
    # The poison goes FIRST, so a stream that halts on it loses every good
    # molecule after it — which is exactly the measured T_1 failure.
    ordered = {"poison": POISON, **good}

    got = ion.charge_at_ph(ordered)
    recovered = [k for k in good if got.get(k) is not None]
    assert len(recovered) == len(good), (
        f"only {len(recovered)}/{len(good)} good molecules survived a poisoned "
        "batch; the split-retry is not isolating the failure")
    assert got["poison"] is None, "the poison itself should be unconvertible"


def test_an_unconvertible_molecule_is_None_and_not_zero():
    """None and 0 are different facts.

    Folding failures into `neutral` would put them in the largest stratum.
    """
    got = ion.charge_at_ph({"poison": POISON})
    assert got["poison"] is None
    assert ion.charge_class(None) == "unknown"
    assert ion.charge_class(0) == "neutral"


def test_phosphate_is_detected_on_the_one_seed_that_has_it():
    for name, (smiles, _) in SEEDS.items():
        assert ion.has_phosphate(smiles) is (name == "guo_pfizer"), name


def test_charge_class_is_coarse_on_purpose():
    """Vina has no electrostatic term, so -1 vs -2 is not a resolution we have."""
    assert ion.charge_class(-1) == ion.charge_class(-2) == "anion"
    assert ion.charge_class(1) == ion.charge_class(2) == "cation"


def test_an_id_that_would_break_the_title_field_is_refused():
    with pytest.raises(ValueError, match="tab or newline"):
        ion.protonate({"bad\tid": "CCO"})
