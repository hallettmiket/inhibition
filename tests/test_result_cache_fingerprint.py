"""
Purpose: Pin the cache fingerprint that stops superseded results being reused.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: synthetic result dicts
Output: pass/fail

This project has now been bitten four times by a cache keyed on nothing but the
result file's existence:

  1. the MD cache returned a 40 ps smoke-test trajectory for a 2 ns request;
  2-4. three separate copies of the MM-GBSA result cache returned pre-D0033
     energies, leaving 11 of 27 T_4 rows wrong by up to 28 kcal/mol and
     inverting the chloroacetamide ordering.

Every one was silent: the run reported success and the number looked plausible.
The rule these tests hold is that a cache entry must carry the parameters that
produced it, and that a MISSING marker means stale -- anything written before
fingerprinting existed cannot be assumed current.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import mmgbsa as mg  # noqa: E402


def _current(**over) -> dict:
    d = {"dG_kcal": -15.0,
         "energy_terms": list(mg.ENERGY_TERMS),
         "igb": mg.IGB,
         "pb_radii": mg.PB_RADII,
         "junction_frcmod": mg.JUNCTION_FRCMOD.name}
    d.update(over)
    return d


def test_current_covalent_result_is_reused():
    assert mg.cached_result_is_current(_current()) is True


def test_non_covalent_result_records_none_and_is_current():
    """No junction exists, so None is correct -- and must not force a rescore."""
    assert mg.cached_result_is_current(_current(junction_frcmod=None)) is True


def test_missing_marker_means_stale():
    """The D0033 case: results written before fingerprinting must not be trusted."""
    assert mg.cached_result_is_current({"dG_kcal": -15.0}) is False


def test_absent_junction_key_is_stale_but_explicit_none_is_not():
    """The distinction the whole fingerprint turns on."""
    d = _current()
    del d["junction_frcmod"]
    assert mg.cached_result_is_current(d) is False
    assert mg.cached_result_is_current(_current(junction_frcmod=None)) is True


def test_superseded_energy_terms_are_stale():
    assert mg.cached_result_is_current(
        _current(energy_terms=["BOND", "ANGLE"])) is False


def test_dropping_one_term_is_detected():
    """D0033 dropped exactly three terms and the total still looked plausible."""
    partial = [t for t in mg.ENERGY_TERMS if t not in ("1-4 EEL",)]
    assert mg.cached_result_is_current(_current(energy_terms=partial)) is False


def test_superseded_junction_parameters_are_stale():
    assert mg.cached_result_is_current(
        _current(junction_frcmod="cys_gaff2_junction_2.frcmod")) is False


def test_changed_solvent_model_is_stale():
    assert mg.cached_result_is_current(_current(igb=5)) is False
    assert mg.cached_result_is_current(_current(pb_radii="mbondi2")) is False


def test_recorded_failure_is_a_valid_answer():
    """A candidate that genuinely cannot be built should not be retried forever."""
    assert mg.cached_result_is_current({"mmgbsa_error": "tleap failed"}) is True


def test_non_dict_is_never_current():
    for bad in (None, [], "result", 3):
        assert mg.cached_result_is_current(bad) is False


def test_delta_g_stamps_what_the_check_reads():
    """The producer and the validator must agree, or nothing is ever cacheable."""
    legs = {name: mg.LegEnergies(terms={t: 1.0 for t in mg.ENERGY_TERMS})
            for name in ("complex", "receptor", "ligand")}
    out = mg.delta_g(legs)
    assert mg.cached_result_is_current(out) is True
