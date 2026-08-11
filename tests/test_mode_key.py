"""The mode key, and the collision it exists to stop (#53)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import mode_key as mk                          # noqa: E402


def test_a_bare_ident_reports_no_mode_rather_than_mode_zero():
    """None, not 0. Reading a bare ident as mode 0 is the assumption that let
    the sweep/rank collision go unnoticed."""
    assert mk.split_ident("t4_x") == ("t4_x", None)
    assert mk.split_ident("t4_x_m0") == ("t4_x", 0)
    assert mk.split_ident("t4_x_m12") == ("t4_x", 12)


def test_a_molecule_name_containing_m_digits_is_not_truncated():
    assert mk.split_ident("t4_m1abc_m2") == ("t4_m1abc", 2)


def test_the_naive_ident_join_drops_every_mode_zero_row():
    """The failure this module exists to prevent, asserted rather than described.

    rank_v2 writes `_m0`; attack_sweep wrote the bare ident for mode 0. A merge
    on `ident` therefore loses exactly the rows that were simulated, silently.
    """
    rank = pd.DataFrame({"ident": ["t4_a_m0", "t4_a_m1"],
                         "parent_ident": ["t4_a", "t4_a"],
                         "mode": [0, 1], "class_rank": [56, 58]})
    sweep = pd.DataFrame({"ident": ["t4_a"], "parent_ident": ["t4_a"],
                          "mode": [0], "frac_attack_ready": [0.78]})

    naive = rank.merge(sweep[["ident", "frac_attack_ready"]], on="ident", how="left")
    assert naive.frac_attack_ready.notna().sum() == 0, (
        "the naive join is supposed to fail; if it now works the fixture is wrong")

    good = mk.join(rank, sweep[["ident", "frac_attack_ready"]],
                   right_bare_is_mode_zero=True)
    got = good[good.mode_key == "t4_a|0"].frac_attack_ready
    assert len(got) == 1 and got.iloc[0] == pytest.approx(0.78)
    assert good[good.mode_key == "t4_a|1"].frac_attack_ready.isna().all()


def test_bare_is_mode_zero_must_be_asked_for():
    """Off by default: in rank_v2 a bare ident is a molecule-level row, so
    silently calling it mode 0 would invent a mode that was never scored."""
    sweep = pd.DataFrame({"ident": ["t4_a"]})
    assert mk.add_key(sweep)["mode"].isna().all()
    assert mk.add_key(sweep, bare_is_mode_zero=True)["mode"].tolist() == [0]


def test_existing_columns_win_over_parsing_the_label():
    """`ident` is a display label. If the frame states the pair, believe it."""
    d = pd.DataFrame({"ident": ["t4_a"], "parent_ident": ["t4_a"], "mode": [3]})
    assert mk.add_key(d)["mode_key"].tolist() == ["t4_a|3"]
