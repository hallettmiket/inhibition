"""
Purpose: Pin the D0033 defect -- a leg total that silently dropped three of
         sander's energy terms and still looked like a plausible dG.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: a synthetic sander FINAL RESULTS block
Output: pass/fail

The bug was not that a number was wrong by a little. It was that
`sum(terms.get(k, 0.0) for k in ENERGY_TERMS)` treats a key that does not
exist as a term worth zero, so a parse failure and a genuinely-zero term are
indistinguishable. These tests hold both ends: the parser must read the
multi-word labels, and the total must refuse to sum a partial energy.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import mmgbsa  # noqa: E402

# A real block, trimmed. Note "1-4 VDW" and "1-4 EEL" carry a space, which is
# exactly what the original token regex could not survive.
BLOCK = """
                    FINAL RESULTS

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
   1000      -7.2863E+03     2.4170E-01     5.4982E+00     CE1       371

 BOND    =      100.4672  ANGLE   =      344.1302  DIHED      =      633.9605
 VDWAALS =    -1304.1752  EEL     =    -9419.9331  EGB        =    -2651.2718
 1-4 VDW =      448.2115  1-4 EEL =     4403.4840  RESTRAINT  =        0.0000
 ESURF   =       39.5085
 CMAP    =      119.3561
"""


def test_multiword_terms_are_parsed_under_their_own_names():
    terms = mmgbsa.parse_energy_block(BLOCK)
    assert terms["1-4 VDW"] == pytest.approx(448.2115)
    assert terms["1-4 EEL"] == pytest.approx(4403.4840)
    assert terms["CMAP"] == pytest.approx(119.3561)


def test_one_four_eel_does_not_collide_with_plain_eel():
    """The original parser stored 1-4 EEL's value under 'EEL' or dropped it."""
    terms = mmgbsa.parse_energy_block(BLOCK)
    assert terms["EEL"] == pytest.approx(-9419.9331)
    assert terms["1-4 EEL"] != terms["EEL"]


def test_total_matches_sander_own_total():
    """The summed terms must equal the ENERGY sander prints on the NSTEP line.

    This is the check that would have caught the bug on day one: sander states
    its own total, and the code that re-derives it can be compared against it.
    RESTRAINT is 0.0 here so its exclusion does not affect the comparison.
    """
    terms = mmgbsa.parse_energy_block(BLOCK)
    total = mmgbsa.LegEnergies(terms=terms).total
    assert total == pytest.approx(-7286.3, abs=0.1)


def test_partial_parse_raises_instead_of_summing_to_something_plausible():
    """A missing term must fail loudly, not quietly shrink the total."""
    terms = mmgbsa.parse_energy_block(BLOCK)
    del terms["1-4 EEL"]
    with pytest.raises(mmgbsa.MMGBSAError, match="missing"):
        _ = mmgbsa.LegEnergies(terms=terms).total


def test_ligand_leg_may_omit_cmap():
    """A ligand-only leg has no protein backbone, so CMAP is legitimately absent."""
    terms = mmgbsa.parse_energy_block(BLOCK)
    del terms["CMAP"]
    total = mmgbsa.LegEnergies(terms=terms).total
    assert total == pytest.approx(-7286.3 - 119.3561, abs=0.1)
