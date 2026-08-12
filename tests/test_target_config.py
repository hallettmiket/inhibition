"""The target config, and the one thing it must refuse to do."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import target_config as tc                     # noqa: E402


def test_the_config_exists_and_names_the_target():
    c = tc.load()
    assert tc.get("target.name", c) == "Pin1"
    assert tc.get("target.pdb", c) == "3IKD"
    assert tc.get("target.anchor", c) == "Cys113"


def test_the_sweep_floor_refuses_rather_than_guessing():
    """It is measured per target, never inherited. A plausible default would be
    indistinguishable from a measured one in every artefact downstream."""
    with pytest.raises(tc.ConfigError) as e:
        tc.sweep_floor()
    assert "measured, not" in str(e.value)
    assert "pilot" in str(e.value)


def test_a_missing_key_raises_instead_of_returning_none():
    with pytest.raises(tc.ConfigError):
        tc.get("md.no_such_setting")
    assert tc.get("md.no_such_setting", default=7) == 7


def test_the_decisions_this_release_rests_on_are_recorded():
    c = tc.load()
    assert tc.get("splitting.stage2.enabled", c) is True          # #61
    assert tc.get("splitting.stage2.cut_diameter_a", c) == 2.0
    assert tc.get("md.salt_molar", c) == 0.15                     # #57
    assert tc.get("sweep_rule.select_by", c) == "per_mode"        # #53
    assert tc.get("chemistry.docked_species", c) == "ph_7.4"      # #58
    assert tc.get("docking.persist_all_poses", c) is True         # #44


def test_summary_states_that_the_floor_is_unset():
    assert "UNSET" in tc.summary()
