"""
Purpose: Tests for parsing ALL of Vina's reported modes, not just the first.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-01
Input: synthetic Vina-GPU output PDBQTs
Output: pytest pass/fail

Issue #10. `collect_scores` regex-searched for the FIRST `REMARK VINA RESULT`
line and returned one float per ligand, so 8 of every 9 poses we computed were
written to disk and never read. These tests pin the full parse and, more
importantly, pin the two things that are easy to get wrong about it:

* mode 1's RMSDs are 0.0 BY CONSTRUCTION (Vina reports RMSD relative to the
  best mode), so a nearest-neighbour spread that includes mode 1 is always 0
  and would read as perfect convergence;
* a single-mode ligand has no neighbour at all, and must report NaN rather
  than 0.0 for the same reason.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import noncovalent_dock_run as ncd          # noqa: E402


def _pdbqt(modes: list[tuple[float, float, float]]) -> str:
    out = []
    for i, (aff, lb, ub) in enumerate(modes, 1):
        out += [f"MODEL {i}",
                f"REMARK VINA RESULT:    {aff:.1f}      {lb:.3f}      {ub:.3f}",
                "ATOM      1  C   LIG A   1       0.000   0.000   0.000",
                "ENDMDL"]
    return "\n".join(out) + "\n"


def _write(tmp_path: Path, name: str, modes) -> Path:
    p = tmp_path / f"{name}_out.pdbqt"
    p.write_text(_pdbqt(modes))
    return p


def test_every_mode_is_parsed_not_just_the_first(tmp_path):
    _write(tmp_path, "lig1", [(-8.3, 0.0, 0.0), (-8.2, 1.1, 2.2),
                              (-7.9, 3.4, 5.6)])
    df = ncd.collect_modes(tmp_path)
    assert len(df) == 1
    r = df.iloc[0]
    assert r["candidate_id"] == "lig1"
    assert r["vina_n_modes"] == 3
    assert r["vina_affinity"] == pytest.approx(-8.3)


def test_the_best_affinity_is_unchanged_from_the_old_single_value_parse(tmp_path):
    """Whatever else changes, the ranked column must not move."""
    _write(tmp_path, "lig1", [(-8.3, 0.0, 0.0), (-8.2, 1.1, 2.2)])
    _write(tmp_path, "lig2", [(-6.0, 0.0, 0.0)])
    assert ncd.collect_scores(tmp_path) == {"lig1": -8.3, "lig2": -6.0}


def test_mode2_gap_is_the_margin_the_top_pose_was_chosen_by(tmp_path):
    """#10's point, as an assertion: the margin is small, and now visible."""
    _write(tmp_path, "lig1", [(-8.30, 0.0, 0.0), (-8.20, 1.1, 2.2)])
    r = ncd.collect_modes(tmp_path).iloc[0]
    assert r["vina_mode2_gap"] == pytest.approx(0.10)
    assert r["vina_affinity_spread"] == pytest.approx(0.10)


def test_nearest_neighbour_rmsd_excludes_mode_1(tmp_path):
    """Mode 1's RMSD to itself is 0.0 by construction and means nothing.

    Including it would report 0.0 for every ligand ever docked, which reads as
    perfect pose convergence -- a populated, plausible, entirely wrong number.
    """
    _write(tmp_path, "lig1", [(-8.3, 0.0, 0.0), (-8.2, 1.5, 2.0),
                              (-7.9, 1.1, 1.8)])
    r = ncd.collect_modes(tmp_path).iloc[0]
    assert r["vina_mode_rmsd_nn"] == pytest.approx(1.1)


def test_a_single_mode_ligand_reports_no_neighbour_rather_than_zero(tmp_path):
    _write(tmp_path, "lonely", [(-7.0, 0.0, 0.0)])
    r = ncd.collect_modes(tmp_path).iloc[0]
    assert r["vina_n_modes"] == 1
    assert pd.isna(r["vina_mode2_gap"])
    assert pd.isna(r["vina_mode_rmsd_nn"])


def test_an_empty_or_unparseable_file_is_skipped_not_zero_filled(tmp_path):
    (tmp_path / "broken_out.pdbqt").write_text("MODEL 1\nENDMDL\n")
    _write(tmp_path, "good", [(-8.0, 0.0, 0.0), (-7.5, 2.0, 3.0)])
    df = ncd.collect_modes(tmp_path)
    assert set(df["candidate_id"]) == {"good"}


def test_the_declared_columns_are_the_columns_produced(tmp_path):
    """MODE_COLS drives the merge's drop list; drift makes it silently stale."""
    _write(tmp_path, "lig1", [(-8.0, 0.0, 0.0), (-7.5, 2.0, 3.0)])
    df = ncd.collect_modes(tmp_path)
    assert list(df.columns) == ["candidate_id", *ncd.MODE_COLS]


def test_no_modes_at_all_yields_the_declared_empty_shape(tmp_path):
    df = ncd.collect_modes(tmp_path)
    assert df.empty
    assert list(df.columns) == ["candidate_id", *ncd.MODE_COLS]
