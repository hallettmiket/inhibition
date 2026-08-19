#!/usr/bin/env python3
"""
Purpose: a persisted pose can be traced to its own measurements, and a mode's width is bounded.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18

D0086 steps 1 and 2, plus the seed (#77). Three defects that all made the screen
unauditable rather than merely imprecise:

  #76  the all-poses SDF numbered poses by POSITION and was never rewritten, so
       the stored cloud could not be joined to the per-pose table beside it --
       and after a re-screen it described a different run entirely
  #77  AutoDock-GPU was invoked without --seed, so no run could be reproduced
  D0086 stage 1 clustered with DBSCAN, which bounds the LINK, so a mode's
       diameter was unbounded and a "mode" could hold two populations
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_modes as pmod                # noqa: E402
from shared import target_config as tc               # noqa: E402


# --------------------------------------------------------------- provenance --
def test_every_persisted_pose_carries_its_conformer_id():
    """`pose_idx` is the key the per-pose table uses; without it on the SDF the
    cloud and its measurements are two lists at the same offset."""
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    w = src.split("def write_sdf")[1].split("\ndef ")[0]
    assert 'SetProp("pose_idx"' in w
    assert re.search(r'SetProp\("pose_idx",\s*str\(int\(i\)\)\)', w), \
        "pose_idx must be the conformer id i, not the enumeration position"


def test_the_cloud_is_rewritten_with_its_run():
    """`if not adest.exists()` kept the PREVIOUS run's cloud next to the current
    run's table -- config's `persist_all_poses` calls that out as a rule (#44)."""
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    blk = src.split("ALL_POSE_DIR.mkdir")[1][:1400]
    # CODE ONLY. The comment explaining the fix necessarily quotes the guard it
    # removed, and a naive substring check matches its own rationale.
    code = "\n".join(ln for ln in blk.splitlines()
                     if not ln.lstrip().startswith("#"))
    assert "if not adest.exists()" not in code


# --------------------------------------------------------------------- seed --
def test_docking_can_be_seeded_and_is_by_default():
    src = (REPO / "scripts" / "nac_screen.py").read_text()
    d = src.split("def dock(")[1].split("\ndef ")[0]
    assert '"--seed"' in d, "the autodock argv must be able to carry a seed"
    assert "if seed is not None" in d, "seed=None must keep clock behaviour"
    assert tc.get("docking.seed", default=None) is not None


def test_the_screen_passes_the_configured_seed():
    """Threaded as a PARAMETER, not read from `args` inside the worker.

    `one()` has no `args` in scope, and reading it there made every molecule fail
    with `name 'args' is not defined` while the screen reported "0 modes from 1
    molecules" and exited 0 -- a broken run that looked like an empty one.
    """
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    assert "seed=_seed_for(args)" in src, "main() must resolve the seed"
    worker = src.split("def one(")[1].split("\ndef ")[0]
    assert "seed=seed" in worker, "one() must pass its seed parameter through"
    assert "_seed_for(args)" not in worker, "one() has no args in scope"


def test_a_negative_seed_means_clock_behaviour():
    """A replicate experiment needs INDEPENDENT draws; with `docking.seed` pinned
    it would otherwise screen five identical clouds and answer nothing."""
    src = (REPO / "scripts" / "nac_screen_v2.py").read_text()
    f = src.split("def _seed_for(")[1].split("\ndef ")[0]
    assert "int(v) < 0" in f and "return None" in f


# ------------------------------------------------------------------ linkage --
def _chain(step: float, n: int) -> np.ndarray:
    """n poses in a line, each `step` from the next: the shape DBSCAN merges."""
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            d[i, j] = abs(i - j) * step
    return d


def test_dbscan_merges_a_chain_wider_than_its_own_tolerance():
    """The defect, stated as a test: eps bounds the NEAREST neighbour, so the
    group's width grows with the chain. This is why a stage-1 parent spans
    4.22 A under a 3 A rule."""
    lab = pmod._dbscan(_chain(1.0, 10), eps=1.5, min_samples=3)
    assert len(set(lab[lab >= 0])) == 1, "DBSCAN chains the whole line into one"


def test_complete_linkage_bounds_the_diameter():
    """The fix: no two members further apart than the tolerance."""
    d = _chain(1.0, 10)
    lab = pmod._complete_linkage(d, diameter=1.5, min_samples=3)
    for c in set(lab[lab >= 0]):
        m = lab == c
        assert d[np.ix_(m, m)].max() <= 1.5 + 1e-9


def test_both_rules_call_a_tight_cloud_one_mode():
    """The fix must not fragment a genuinely tight cluster."""
    d = _chain(0.05, 12)
    for lab in (pmod._dbscan(d, 1.5, 3), pmod._complete_linkage(d, 1.5, 3)):
        assert len(set(lab[lab >= 0])) == 1


def test_small_clusters_are_noise_under_both_rules():
    """`-1` must mean the same thing either way, or the two are not swappable."""
    d = np.full((4, 4), 9.0)
    np.fill_diagonal(d, 0.0)
    for lab in (pmod._dbscan(d, 1.0, 3), pmod._complete_linkage(d, 1.0, 3)):
        assert (lab == -1).all()


def test_the_default_is_unchanged_until_the_comparison_is_decided():
    """D0086 step 3: measure before switching. exp/3_linkage is that measurement."""
    src = (REPO / "shared" / "pose_modes.py").read_text()
    assert 'method: str = "dbscan"' in src
    assert (REPO / "exp" / "3_linkage" / "run_all.py").is_file()


def test_an_unknown_method_raises_rather_than_silently_picking_one():
    with pytest.raises(ValueError):
        pmod.split(np.zeros((5, 6)), method="kmeans")
