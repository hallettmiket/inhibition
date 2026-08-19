#!/usr/bin/env python3
"""
Purpose: the election benchmark is anchored on a validated POSE, never a mixture's mean.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18

WHY THIS TEST EXISTS. `exp/4_election` first matched a replicate's modes against
`dir_x/y/z` from the rank table -- the MEAN warhead direction of the validated
mode. That mode is a mixture (D0086) whose poses span about 83 degrees, so its
mean direction is not a stable quantity, and matching a fresh mode's mean against
a mixture's mean rewards whichever rule produces mixtures. It reported DBSCAN
9/10 and complete linkage 4/10; re-anchored on the pose that was actually
elevated to 100 ns, the same experiment reports 6/20 and 11/20. The measurement
had been confirming its own assumption.

A benchmark that can be wrong in that direction is worse than no benchmark, so
the anchor is asserted here rather than left to a comment.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

EXP = REPO / "exp" / "4_election" / "run_all.py"
pytestmark = pytest.mark.skipif(not EXP.is_file(), reason="election experiment absent")


def _mod():
    spec = importlib.util.spec_from_file_location("election_exp", EXP)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_every_reference_names_the_pose_it_was_validated_on():
    """A reference is a mode a run elected AND a trajectory confirmed. Without
    the pose_rank there is no single geometry to anchor on."""
    m = _mod()
    assert m.REFERENCES, "at least one validated molecule must be registered"
    for cand, meta in m.REFERENCES.items():
        assert "pose_rank" in meta, f"{cand}: no pose_rank — the anchor is ambiguous"
        assert "validated" in meta and meta["validated"], \
            f"{cand}: must record what the 100 ns run measured"
        assert "topic" in meta and "mode" in meta


def test_the_reference_is_not_read_from_the_rank_table_means():
    """`dir_x`/`centroid_x` are per-MODE averages. Anchoring on them is the
    retracted measurement."""
    src = EXP.read_text()
    resolver = src.split("def reference_for(")[1].split("\ndef ")[0]
    for col in ("dir_x", "dir_y", "dir_z", "centroid_x", "centroid_y", "centroid_z"):
        assert f"r.{col}" not in resolver, \
            f"reference_for reads {col} — that is a mixture's mean, not a pose"
    assert "pose_rank" in resolver and "features" in resolver


def test_the_resolved_reference_is_a_unit_direction_and_a_point():
    m = _mod()
    for cand in m.REFERENCES:
        ref = m.reference_for(cand)
        assert ref["centroid"].shape == (3,)
        assert ref["direction"].shape == (3,)
        assert np.isclose(np.linalg.norm(ref["direction"]), 1.0, atol=1e-6), \
            "the warhead direction must be a unit vector"
        assert cand in ref["note"] or ref["candidate"] == cand


def test_an_unvalidated_molecule_is_refused_not_guessed():
    """Scoring against a mode nothing confirmed would manufacture a benchmark."""
    m = _mod()
    with pytest.raises(SystemExit):
        m.reference_for("t4_not_a_real_molecule")


def test_both_validated_molecules_resolve_to_different_geometries():
    """A resolver that silently returned the same anchor for everything would
    pass every other test here."""
    m = _mod()
    if len(m.REFERENCES) < 2:
        pytest.skip("only one validated molecule registered")
    refs = [m.reference_for(c) for c in m.REFERENCES]
    assert np.linalg.norm(refs[0]["centroid"] - refs[1]["centroid"]) > 0.5
