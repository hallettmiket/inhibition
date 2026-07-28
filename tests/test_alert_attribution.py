"""
Purpose: Tests for attributing structural alerts without cutting the molecule.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: representative covalent candidates
Output: pytest pass/fail

D0025. The old two-tier gate excised the R-group and capped the severed bond,
which invented functional groups: an amide cut from its nitrogen and H-capped is
a formamide (BRENK: aldehyde), a thioether is a thiol. Those two accounted for
3,014 of T_3's 5,270 rejections.

The replacement screens the intact molecule and attributes each alert by the
atoms it matched. These tests pin the properties that make that correct, and the
two regressions found while building it — both of which reintroduced the exact
false positive the two-tier design exists to prevent, through new doors.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from rdkit import RDLogger

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import alerts                          # noqa: E402

RDLogger.DisableLog("rdApp.*")

T4_SCOPE = "N[CH]1CCS(=O)(=O)C1"                   # sulfolane + N (warhead is outside)
T3_SCOPE = "N(C(=O)C=C)C1CCS(=O)(=O)C1"            # includes the acrylamide

CHLOROACETAMIDE = "O=C(CCl)N(c1ccccc1)[C@@H]1CCS(=O)(=O)C1"
UREA_DECORATION = "C=CC(=O)N(C(=O)NCc1ccccc1OC)C1CCS(=O)(=O)C1"


def test_no_aldehyde_is_invented_from_an_amide():
    """The defect that started D0025, stated directly."""
    att = alerts.attribute_alerts(UREA_DECORATION, T3_SCOPE)
    assert "aldehyde" not in att.rgroup + att.boundary + att.core, (
        "an aldehyde was attributed to a molecule that contains none; the "
        "screen is cutting the molecule again")


def test_warhead_alerts_are_excused_as_core():
    """A chloroacetamide's own reactivity is the mechanism, not a liability."""
    att = alerts.attribute_alerts(CHLOROACETAMIDE, T4_SCOPE)
    assert att.core, "the warhead's alerts should be attributed to the core"
    assert "alkyl_halide" in att.core


def test_alpha_halo_carbonyl_does_not_straddle_the_boundary():
    """Regression: excusing only the reactive ATOMS is not enough.

    chloroacetamide's reactive-atom SMARTS is `[CH2][Cl]` — two atoms. Excusing
    just those leaves the carbonyl exposed, so BRENK's `alpha_halo_carbonyl`
    spans core and decoration and is charged to the decoration. That rejected
    all 1,683 T_4 survivors. The whole warhead group must be excused.
    """
    att = alerts.attribute_alerts(CHLOROACETAMIDE, T4_SCOPE)
    assert "alpha_halo_carbonyl" not in att.boundary, (
        "the alpha-halo-carbonyl IS the chloroacetamide warhead; it must not be "
        "charged to the decoration")
    assert att.attributable == 0, f"unexpected decoration alerts: {att.rgroup + att.boundary}"


def test_whole_t4_library_still_passes_broadly():
    """Regression: attributing by core alone rejected 1683/1683.

    A sanity bound rather than an exact count — the point is that a plain
    sulfopin+aryl candidate is not a liability.
    """
    for r in ("c1ccccc1", "Cc1ccccc1", "CCc1ccccc1"):
        smiles = f"O=C(CCl)N({r})[C@@H]1CCS(=O)(=O)C1"
        att = alerts.attribute_alerts(smiles, T4_SCOPE)
        assert att.attributable == 0, f"{smiles} -> {att.rgroup + att.boundary}"


def test_a_real_decoration_alert_is_still_caught():
    """The gate must not have been softened into uselessness."""
    nitro = "O=C(CCl)N(c1ccc([N+](=O)[O-])cc1O)[C@@H]1CCS(=O)(=O)C1"
    att = alerts.attribute_alerts(nitro, T4_SCOPE)
    assert att.attributable > 0, "a nitro/phenol decoration should still flag"


def test_boundary_alerts_are_reported_separately():
    """A spanning match is charged to the decoration but stays visible.

    The T_3 scaffold nitrogen already carries the acrylamide carbonyl; an acyl
    decoration makes it a genuine acyclic imide. Real chemistry, and a judgement
    call about attribution, so it must not be silently folded in with clean hits.
    """
    acyl = "C=CC(=O)N(C(=O)c1ccccc1)C1CCS(=O)(=O)C1"
    att = alerts.attribute_alerts(acyl, T3_SCOPE)
    assert "acyclic_imide" in att.boundary
    assert att.attributable >= 1
    cols = alerts.two_tier(acyl, T3_SCOPE).to_columns()
    assert cols["boundary_alert_total"] >= 1
    assert "acyclic_imide" in cols["boundary_alert_names"]


def test_missing_core_is_not_a_pass():
    att = alerts.attribute_alerts("c1ccccc1C(=O)NC", T4_SCOPE)
    assert not att.core_found
    assert not alerts.two_tier("c1ccccc1C(=O)NC", T4_SCOPE).passes_gate


def test_gate_and_stamp_agree():
    """A frame must not disagree with itself.

    Re-running annotation over its own output inherited the previous run's
    rejections, leaving 1,114 T_3 rows stamped `alerts` while carrying
    `alert_gate_pass = True`.
    """
    import pandas as pd

    from shared import annotate as ann

    df = pd.DataFrame({"canonical_smiles": [CHLOROACETAMIDE, UREA_DECORATION],
                       "rejected_at": ["alerts", "alerts"]})   # stale stamps
    out = ann.annotate(df, approach="t4", core_smarts=T4_SCOPE)
    disagree = out["rejected_at"].notna() & out["alert_gate_pass"].eq(True)
    assert not disagree.any(), "stamped rejected while the gate says pass"


def test_excused_alert_passes_but_is_carried(): 
    """D0026: a named excusal is a decision, and it must stay visible."""
    acyl = "C=CC(=O)N(C(=O)c1ccccc1)C1CCS(=O)(=O)C1"
    strict = alerts.two_tier(acyl, T3_SCOPE)
    assert not strict.passes_gate

    lenient = alerts.two_tier(acyl, T3_SCOPE, excused_alerts=("acyclic_imide",))
    assert lenient.passes_gate
    assert "acyclic_imide" in lenient.excused, "the excusal must be recorded"
    cols = lenient.to_columns()
    assert cols["excused_alert_names"] == "acyclic_imide"
    assert cols["excused_alert_total"] == 1


def test_excusing_one_alert_does_not_admit_others():
    """The point of naming: excusing imides must not let a thioester through."""
    both = "C=CC(=O)N(C(=O)CSC(=O)C)C1CCS(=O)(=O)C1"
    r = alerts.two_tier(both, T3_SCOPE, excused_alerts=("acyclic_imide",))
    other = [n for n in r.attributed.rgroup + r.attributed.boundary
             if n != "acyclic_imide"]
    if other:
        assert not r.passes_gate, f"excusal leaked past the named alert: {other}"
