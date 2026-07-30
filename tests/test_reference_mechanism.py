"""
Purpose: The reference set must reject a mechanism value no stratum selects on.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: a temporary reference CSV
Output: pass/fail

WHY THIS TEST EXISTS. Liu-2024-C3 was added with mechanism
"non_covalent_active_site", which is descriptive and correct English and
matches nothing. The gate selects strata with `master.mechanism ==
"non_covalent"`, so the active was dropped with no error, no warning, and a log
line reading "5 actives" that looked exactly like the previous run. It cost the
sixth independent chemotype -- the gate's own verdict floor -- and was found
only because the count was one lower than expected.
"""

from __future__ import annotations

import pandas as pd
import pytest

from shared import reference_set as rs


def _write(tmp_path, mechanism: str):
    master = tmp_path / "master.csv"
    pd.DataFrame([{
        "name": "Test-Binder",
        "canonical_smiles": "c1ccccc1",
        "mechanism": mechanism,
    }]).to_csv(master, index=False)
    return master


def test_unrecognised_mechanism_raises(tmp_path):
    master = _write(tmp_path, "non_covalent_active_site")
    with pytest.raises(rs.ReferenceSetError) as exc:
        rs.load(master_path=master)
    msg = str(exc.value)
    assert "non_covalent_active_site" in msg
    assert "Test-Binder" in msg, "the offending row must be named, not just the value"


@pytest.mark.parametrize("mechanism", sorted(rs.VALID_MECHANISMS))
def test_every_valid_mechanism_is_accepted(tmp_path, mechanism):
    master = _write(tmp_path, mechanism)
    rs.load(master_path=master)


def test_the_shipped_reference_set_validates():
    """The file the pipeline actually defaults to must pass its own check."""
    ref = rs.load()
    assert set(ref.master["mechanism"].dropna()) <= rs.VALID_MECHANISMS


def test_liu_2024_c3_is_selectable_as_non_covalent():
    """The specific active whose mechanism typo cost the sixth chemotype."""
    master = rs.load().master
    row = master[master["name"] == "Liu-2024-C3"]
    assert len(row) == 1, "Liu-2024-C3 missing from the reference set"
    assert row.iloc[0]["mechanism"] == "non_covalent"
