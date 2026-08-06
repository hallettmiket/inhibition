"""
Purpose: the BPMD scoring layer, tested on synthetic trajectories with known answers.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06

Trajectories are hand-built so the right answer is arithmetic, not whatever the
module says. A stable pose, a pose that escapes, and a pose that was never in the
window are each constructed explicitly.

The index test carries the most weight. PLUMED atom indices are 1-based; passing
a Python 0-based index biases a NEIGHBOURING atom, and the simulation runs to
completion and reports a perfectly plausible stability. That is this project's
signature failure shape, and it is worth a guard that cannot be bypassed.
"""

from __future__ import annotations

import numpy as np
import pytest

from shared import bpmd


def write_colvar(tmp_path, distances, biases, name="COLVAR"):
    p = tmp_path / name
    lines = ["#! FIELDS time d metad.bias"]
    lines += [f"{i*0.5:.1f} {d:.4f} {b:.4f}" for i, (d, b) in enumerate(zip(distances, biases))]
    p.write_text("\n".join(lines) + "\n")
    return p


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def test_read_colvar_skips_comments_and_blanks(tmp_path):
    p = tmp_path / "COLVAR"
    p.write_text("#! FIELDS time d metad.bias\n\n0.0 0.35 0.0\n# mid-file comment\n0.5 0.36 1.2\n")
    d, b = bpmd.read_colvar(p)
    assert list(d) == pytest.approx([0.35, 0.36])
    assert list(b) == pytest.approx([0.0, 1.2])


def test_missing_colvar_raises(tmp_path):
    with pytest.raises(bpmd.BPMDError, match="no COLVAR"):
        bpmd.read_colvar(tmp_path / "nope")


def test_empty_colvar_raises(tmp_path):
    p = tmp_path / "COLVAR"; p.write_text("#! FIELDS time d metad.bias\n")
    with pytest.raises(bpmd.BPMDError, match="no readable frames"):
        bpmd.read_colvar(p)


# --------------------------------------------------------------------------
# per-replica reduction
# --------------------------------------------------------------------------

def test_stable_pose_never_escapes(tmp_path):
    """Sits at 0.35 nm, inside the window, for the whole run."""
    n = 200
    c = write_colvar(tmp_path, [0.35] * n, np.linspace(0, 50, n))
    r = bpmd.analyse_replica(c, 0)
    assert r.frac_in_window == pytest.approx(1.0)
    assert not r.escaped
    assert np.isnan(r.bias_at_exit_kj), "a pose that never left has no escape cost"


def test_escaping_pose_records_the_bias_it_took(tmp_path):
    """In the window for half the run, then pushed out past 1.0 nm."""
    d = [0.35] * 100 + list(np.linspace(0.45, 1.4, 100))
    b = list(np.linspace(0, 40, 200))
    r = bpmd.analyse_replica(write_colvar(tmp_path, d, b), 1)
    assert r.escaped
    assert r.frac_in_window == pytest.approx(0.5, abs=0.02)
    # the bias standing at the first frame past 1.0 nm
    first = next(i for i, x in enumerate(d) if x >= bpmd.UNBOUND_NM)
    # 1e-3, because the fixture writes COLVAR at 4 decimal places — tightening
    # this further would test the fixture's formatting, not the module.
    assert r.bias_at_exit_kj == pytest.approx(b[first], abs=1e-3)


def test_pose_outside_the_window_scores_zero_occupancy(tmp_path):
    """Never near-attack at all — 0.6 nm throughout, well past the 0.42 nm edge."""
    r = bpmd.analyse_replica(write_colvar(tmp_path, [0.6] * 50, [0.0] * 50), 2)
    assert r.frac_in_window == 0.0
    assert not r.escaped


def test_window_edges_are_inclusive(tmp_path):
    d = [bpmd.NAC_MIN_NM, bpmd.NAC_MAX_NM, bpmd.NAC_MAX_NM + 0.01]
    r = bpmd.analyse_replica(write_colvar(tmp_path, d, [0.0] * 3), 3)
    assert r.frac_in_window == pytest.approx(2 / 3)


# --------------------------------------------------------------------------
# pooling
# --------------------------------------------------------------------------

def test_non_escapers_are_not_dropped_from_the_median(tmp_path):
    """Dropping them would bias the escape cost toward the WEAK poses.

    Only poses that left have a finite cost, so a median over escapers alone
    describes the ones that failed — exactly backwards for a stability ranking.
    """
    res = [bpmd.ReplicaResult(0, 0.4, 1.0, float("nan"), False),
           bpmd.ReplicaResult(1, 1.2, 0.5, 10.0, True),
           bpmd.ReplicaResult(2, 1.3, 0.5, 30.0, True)]
    s = bpmd.combine(res)
    assert s.n_escaped == 2
    assert s.median_bias_to_escape_kj == pytest.approx(20.0)
    assert s.mean_frac_in_window == pytest.approx(2 / 3, abs=1e-6)


def test_spread_is_reported(tmp_path):
    res = [bpmd.ReplicaResult(i, 0.4, f, float("nan"), False)
           for i, f in enumerate([0.2, 0.8])]
    s = bpmd.combine(res)
    assert s.spread_frac > 0, "replica variation IS the uncertainty and must surface"


def test_score_requires_both_occupancy_and_cost():
    """Multiplicative on purpose — one without the other is not stability."""
    stable = bpmd.PoseStability(5, 0.9, 0.05, 30.0, 5)
    resists_but_absent = bpmd.PoseStability(5, 0.05, 0.01, 30.0, 5)
    present_but_free = bpmd.PoseStability(5, 0.9, 0.05, 0.5, 5)
    assert stable.score > resists_but_absent.score
    assert stable.score > present_but_free.score


def test_score_is_finite_when_nothing_escaped():
    s = bpmd.PoseStability(3, 0.8, 0.1, float("nan"), 0)
    assert np.isfinite(s.score) and s.score > 0


def test_combine_refuses_an_empty_set():
    with pytest.raises(bpmd.BPMDError):
        bpmd.combine([])


# --------------------------------------------------------------------------
# the index guard
# --------------------------------------------------------------------------

def test_plumed_input_rejects_zero_based_indices():
    """A 0-based index biases the neighbouring atom and the run still succeeds."""
    with pytest.raises(bpmd.BPMDError, match="1-based"):
        bpmd.plumed_input(0, 42)
    with pytest.raises(bpmd.BPMDError, match="1-based"):
        bpmd.plumed_input(42, 0)


def test_atom_indices_convert_and_bounds_check(tmp_path):
    gro = tmp_path / "sys.gro"
    gro.write_text("title\n  100\n" + "\n" * 100 + "1.0 1.0 1.0\n")
    assert bpmd.bpmd_atom_indices(gro, 0, 99) == (1, 100)
    with pytest.raises(bpmd.BPMDError, match="outside"):
        bpmd.bpmd_atom_indices(gro, 100, 5)
    with pytest.raises(bpmd.BPMDError, match="outside"):
        bpmd.bpmd_atom_indices(gro, 5, -1)


def test_plumed_input_biases_the_distance_not_rmsd():
    """The CV choice is the point — D0062 measured RMSD as the wrong endpoint."""
    txt = bpmd.plumed_input(10, 20)
    assert "DISTANCE ATOMS=10,20" in txt
    # Directives only — the header COMMENT mentions RMSD precisely to explain why
    # it is not the CV, and testing the prose would forbid documenting the choice.
    directives = [l for l in txt.splitlines() if l.strip() and not l.lstrip().startswith("#")]
    assert not any("RMSD" in l for l in directives), \
        f"no PLUMED action may use RMSD: {[l for l in directives if 'RMSD' in l]}"
    assert "BIASFACTOR" in txt, "must be well-tempered so the estimate converges"
