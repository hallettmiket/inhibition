"""
Purpose: covalent chemotypes are counted by warhead class, non-covalent by structure (D0045).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: shared.enrichment_gate, shared.warhead_library
Output: pass/fail

WHY THIS EXISTS. The gate matched decoys by warhead class (D0031) and then
counted chemotypes by whole-molecule similarity — one definition to build the
comparison, another to size it. On the real actives the two are nearly
orthogonal and give 4 vs 6 against a floor of 6, so the choice decides whether
a verdict is possible at all.

The counting rule is therefore load-bearing, and it is the kind of rule that
would be easy to "simplify" back into a single clustering call by someone who
did not know it was decided. Pinned here with the reason.
"""

from __future__ import annotations

import pandas as pd
import pytest

from shared import enrichment_gate as eg
from shared import reference_set as rs
from shared import warhead_library as wl


def _actives(warheads, smiles=None):
    n = len(warheads)
    return pd.DataFrame({
        "canonical_smiles": smiles or [f"C{'C' * i}O" for i in range(n)],
        "warhead_class": warheads,
        "label": [1] * n,
    })


# --------------------------------------------------------------------------
# canonical_class: the prose -> vocabulary map
# --------------------------------------------------------------------------

def test_the_same_warhead_written_two_ways_is_one_class():
    """The overcount that would have cleared the floor by accident."""
    a = wl.canonical_class("chloroacetamide")
    b = wl.canonical_class("chloroacetamide (N-methyl peptidomimetic)")
    assert a == b == "chloroacetamide"


def test_absent_or_unestablished_warheads_map_to_none():
    for v in (None, "", "  ", "UNVERIFIED", "unverified", "nan"):
        assert wl.canonical_class(v) is None


def test_an_unknown_warhead_raises_rather_than_becoming_its_own_class():
    """A new active must not silently inflate the denominator.

    This is the failure mode the whole decision is about: a value that becomes
    its own category by default, in the direction that certifies a verdict.
    """
    with pytest.raises(wl.WarheadClassError, match="canonical class_id"):
        wl.canonical_class("some brand new electrophile nobody mapped")


def test_every_warhead_in_the_reference_set_maps_or_is_none():
    """The live data must not trip the raise above."""
    m = pd.read_csv(rs.DEFAULT_MASTER)
    cov = m[m["mechanism"] == "covalent_cys113"]
    for _, r in cov.iterrows():
        wl.canonical_class(r["warhead_class"])  # must not raise


# --------------------------------------------------------------------------
# chemotype_ids: stratum-dependent counting
# --------------------------------------------------------------------------

def test_covalent_counts_by_warhead_not_by_structure():
    # Three molecules, near-identical structurally, three DIFFERENT warheads.
    a = _actives(
        ["chloroacetamide", "sulfamate acetamide",
         "cinnamamide (aryl Michael acceptor; acrylamide-class)"],
        smiles=["CCO", "CCO", "CCO"])
    ids, method = eg.chemotype_ids(a, "covalent")
    assert method == "warhead_class"
    assert eg.n_chemotypes(ids) == 3, (
        "identical structures with distinct warheads must count as three "
        "covalent chemotypes")


def test_covalent_collapses_one_warhead_written_two_ways():
    a = _actives(["chloroacetamide", "chloroacetamide (N-methyl peptidomimetic)"],
                 smiles=["CCO", "c1ccccc1CCN"])
    ids, _ = eg.chemotype_ids(a, "covalent")
    assert eg.n_chemotypes(ids) == 1, (
        "structurally unrelated molecules sharing a warhead are ONE covalent "
        "chemotype")


def test_an_unestablished_warhead_adds_no_chemotype():
    a = _actives(["chloroacetamide", "UNVERIFIED"])
    ids, _ = eg.chemotype_ids(a, "covalent")
    assert eg.n_chemotypes(ids) == 1
    assert -1 in ids, "the unestablished active must be marked, not dropped"


def test_non_covalent_still_uses_structure():
    a = _actives([None, None, None],
                 smiles=["CCCCCCCCO", "c1ccccc1C(=O)NCCN", "OC(=O)Cc1ccncc1"])
    ids, method = eg.chemotype_ids(a, "non_covalent")
    assert method.startswith("ecfp4"), (
        "a reversible binder has no warhead; structure is all there is")
    assert eg.n_chemotypes(ids) >= 1


def test_covalent_without_warhead_data_refuses_rather_than_falling_back():
    """A silent fallback would reintroduce exactly the defect D0045 removed."""
    a = pd.DataFrame({"canonical_smiles": ["CCO"], "label": [1]})
    with pytest.raises(eg.EnrichmentGateError, match="warhead_class"):
        eg.chemotype_ids(a, "covalent")


# --------------------------------------------------------------------------
# The live numbers this decision turns on
# --------------------------------------------------------------------------

def test_the_covalent_stratum_is_underpowered_under_the_chosen_definition():
    """Records the unfavourable outcome so it cannot drift unnoticed.

    D0045 was decided BEFORE these counts were taken, and the honest result is
    that warhead-class counting leaves the covalent gate below its floor of 6.
    If this ever passes 6, that should be because chemotypes were ADDED — not
    because the definition loosened.
    """
    m = pd.read_csv(rs.DEFAULT_MASTER)
    lead = m[(m.mechanism == "covalent_cys113")
             & (m.tier == "lead")
             & (m.canonical_smiles != rs.UNVERIFIED)]
    ids, method = eg.chemotype_ids(lead, "covalent")
    assert method == "warhead_class"
    n = eg.n_chemotypes(ids)
    assert n == 4, (
        f"covalent lead-tier chemotypes changed from 4 to {n}. If chemotypes "
        "were added, update this number and say which. If the DEFINITION "
        "changed, that needs a decision record superseding D0045.")
