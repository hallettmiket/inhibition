"""
Purpose: The curator filter must not silently mis-parse a chemist's constraint.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: small in-memory frames
Output: pass/fail

The failure this guards against is specific and was real: lowercasing the whole
constraint line turned `no [Cl]` into `[cl]`, which is valid SMARTS matching an
aromatic chlorine and therefore matches nothing. The filter reported "removed
0", which is exactly what a working filter reports when no molecule carries the
group. A chemist would have concluded the shortlist was chlorine-free.
"""

from __future__ import annotations

import pandas as pd
import pytest

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "integration" / "app"))
import curate  # noqa: E402


@pytest.fixture
def frame() -> pd.DataFrame:
    return pd.DataFrame({
        "candidate_id": ["a", "b", "c"],
        "canonical_smiles": [
            "Clc1ccccc1",                    # chlorobenzene
            "c1ccccc1",                      # benzene
            "O=S(=O)(N)c1ccccc1",            # benzenesulfonamide
        ],
        "rank": [1, 2, 3],
    })


def test_explicit_smarts_case_is_preserved(frame):
    """`[Cl]` must not be lowercased into `[cl]`, which matches nothing."""
    kept, rules = curate.apply(frame, "no [Cl]")
    assert rules[0].removed == 1, "chlorobenzene should have been removed"
    assert set(kept["candidate_id"]) == {"b", "c"}


def test_smarts_and_the_named_group_agree(frame):
    by_smarts = curate.apply(frame, "no [Cl]")[0]
    by_name = curate.apply(frame, "no chlorine")[0]
    assert list(by_smarts["candidate_id"]) == list(by_name["candidate_id"])


@pytest.mark.parametrize("spec", ["no chlorine", "NO CHLORINE", "No Chlorine"])
def test_group_names_are_case_insensitive(frame, spec):
    assert curate.apply(frame, spec)[0]["candidate_id"].tolist() == ["b", "c"]


def test_require_keeps_only_matches(frame):
    kept, _ = curate.apply(frame, "require sulfonamide")
    assert kept["candidate_id"].tolist() == ["c"]


def test_numeric_bound(frame):
    kept, rules = curate.apply(frame, "heavy_atoms <= 6")
    assert kept["candidate_id"].tolist() == ["b"]
    assert rules[0].removed == 2


def test_order_is_preserved_not_rescored(frame):
    """Filtering removes rows; it must never reorder them (D0041, D0043)."""
    kept, _ = curate.apply(frame, "heavy_atoms <= 12")
    assert kept["rank"].tolist() == sorted(kept["rank"].tolist())


def test_each_rule_reports_its_own_count(frame):
    _, rules = curate.apply(frame, "no chlorine\nheavy_atoms <= 6")
    assert [r.removed for r in rules] == [1, 1]


@pytest.mark.parametrize("bad", ["no chlorien", "mw ~ 450", "no [Cl", "wibble < 3"])
def test_unparseable_constraints_are_refused_not_ignored(frame, bad):
    """A typo must raise, never quietly filter nothing."""
    with pytest.raises(curate.ConstraintError):
        curate.apply(frame, bad)


def test_empty_spec_is_a_no_op(frame):
    kept, rules = curate.apply(frame, "")
    assert len(kept) == len(frame) and rules == []
