"""
Purpose: Pin the ensemble module's uncertainty maths and the hybridisation-aware
         valence check at the covalent junction.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: synthetic series and GAFF type strings
Output: pass/fail

The valence check previously hardcoded "expected 4", which fired on every
naphthoquinone -- output that is correct by design (D0030 keeps the quinone
adduct sp2). A check that cries wolf on correct output is worse than no check,
so the expectation is now derived from the atom type and pinned here.

The statistical inefficiency is what turns a standard error computed from 90
correlated frames into an honest one. If it silently returned 1.0 the SEM would
be too narrow by sqrt(g) -- up to 5x on the trajectories measured here -- and
nothing downstream would look wrong.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import mmgbsa as mg  # noqa: E402
from shared import mmgbsa_ensemble as me  # noqa: E402


# --------------------------------------------------------------------------
# hybridisation-aware valence


@pytest.mark.parametrize("gaff_type,expected", [
    ("c3", 4), ("cx", 4), ("cy", 4),          # sp3 -> saturating adducts
    ("c2", 3), ("ca", 3), ("cc", 3), ("cd", 3),  # sp2 -> quinone stays sp2
    ("c", 3),                                  # carbonyl carbon
    ("c1", 2),                                 # sp
])
def test_expected_valence_by_gaff_type(gaff_type, expected):
    assert mg._expected_valence(gaff_type) == expected


def test_unknown_type_returns_none_rather_than_guessing():
    """An unknown type must not silently default to 4 and reject good output."""
    assert mg._expected_valence("zz") is None
    assert mg._expected_valence("") is None
    assert mg._expected_valence(None) is None


def test_sp2_attachment_is_not_treated_as_an_error():
    """The D0030 case: a re-aromatized quinone adduct carries THREE bonds."""
    assert mg._expected_valence("cc") == 3
    assert mg._expected_valence("cc") != mg._expected_valence("c3")


# --------------------------------------------------------------------------
# statistical inefficiency


def test_independent_series_has_inefficiency_near_one():
    rng = np.random.default_rng(0)
    g = me.statistical_inefficiency(rng.normal(size=4000))
    assert 0.8 < g < 1.6


def test_correlated_series_has_inefficiency_above_one():
    """A random walk is strongly autocorrelated and must be penalised."""
    rng = np.random.default_rng(1)
    walk = np.cumsum(rng.normal(size=2000))
    assert me.statistical_inefficiency(walk) > 5.0


def test_inefficiency_never_below_one():
    """g < 1 would NARROW the error bar, which it must never do."""
    rng = np.random.default_rng(2)
    for n in (4, 10, 100):
        assert me.statistical_inefficiency(rng.normal(size=n)) >= 1.0


def test_constant_series_does_not_divide_by_zero():
    """A fully dissociated ligand gives a near-constant dG (measured: SD 0.04)."""
    assert me.statistical_inefficiency(np.full(50, -7.5)) == 1.0


def test_widening_uses_sqrt_g_not_g():
    """Pin the relationship the reported SEM depends on."""
    rng = np.random.default_rng(3)
    x = np.cumsum(rng.normal(size=500))
    g = me.statistical_inefficiency(x)
    n_eff = x.size / g
    naive = x.std(ddof=1) / np.sqrt(x.size)
    corrected = x.std(ddof=1) / np.sqrt(n_eff)
    assert corrected == pytest.approx(naive * np.sqrt(g), rel=1e-9)
    assert corrected > naive


# --------------------------------------------------------------------------
# cap placement geometry


def test_cap_is_placed_along_the_severed_bond_at_the_stated_length():
    """The link-atom construction, checked on coordinates rather than trusted."""
    frame = np.array([[0.0, 0.0, 0.0],     # 0: anchor (ligand C)
                      [3.0, 0.0, 0.0],     # 1: severed partner (protein S)
                      [0.0, 5.0, 0.0]],    # 2: unrelated
                     dtype=float)
    m = me.LegMap(leg="ligand", complex_index=np.array([0, 2, -1]),
                  cap_index=2, cap_anchor_complex=0,
                  cap_direction_complex=1, cap_bond_length_a=1.09)
    out = me.slice_frame(frame, m)
    assert out[2] == pytest.approx([1.09, 0.0, 0.0])
    assert np.linalg.norm(out[2] - frame[0]) == pytest.approx(1.09)


def test_coincident_atoms_raise_rather_than_produce_a_nan_cap():
    frame = np.zeros((2, 3), dtype=float)
    m = me.LegMap(leg="ligand", complex_index=np.array([0, -1]),
                  cap_index=1, cap_anchor_complex=0,
                  cap_direction_complex=0, cap_bond_length_a=1.09)
    with pytest.raises(me.EnsembleError, match="coincident"):
        me.slice_frame(frame, m)


def test_non_covalent_leg_places_no_cap():
    """cap_index == -1 means every atom comes from the complex, untouched."""
    frame = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=float)
    m = me.LegMap(leg="ligand", complex_index=np.array([1, 0]), cap_index=-1,
                  cap_anchor_complex=-1, cap_direction_complex=-1,
                  cap_bond_length_a=0.0)
    out = me.slice_frame(frame, m)
    assert out[0] == pytest.approx([4.0, 5.0, 6.0])
    assert out[1] == pytest.approx([1.0, 2.0, 3.0])
