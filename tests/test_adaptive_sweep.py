"""Variable-length sweeps must stay comparable with fixed-length ones.

@twu383, 2026-09-03: *"we can just make the runs longer until the mol leaves or
reaches 10 ns max"*.

THE HAZARD IS NOT THE MD, IT IS THE DENOMINATOR. `frac_attack_ready` is a
fraction OF THE RUN. Once runs have different lengths it silently stops being
one quantity: a mode that left at 2 ns and one capped at 10 ns do not share a
window, and the 676 rows already collected are all 1.2 ns. That is this
project's signature defect with a new cause -- two populations under one column
name, both populated and plausible.

`frac_attack_ready_common` is the answer: always the first COMMON_WINDOW_PS, on
every row, whatever the run length. For a fixed 1.2 ns sweep it is identical to
`frac_attack_ready`, so nothing already written changes meaning.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))


def _asw():
    try:
        import attack_sweep as m
    except Exception as exc:                              # noqa: BLE001
        pytest.skip(f"attack_sweep not importable here: {exc}")
    return m


def test_a_fixed_length_sweep_is_unchanged_by_the_new_column():
    """The 676 rows already collected must keep meaning what they meant."""
    a = _asw()
    n = 100
    frame_ps = a.COMMON_WINDOW_PS / n              # exactly one common window
    dist = np.linspace(2.9, 3.4, n)
    ang = np.full(n, 5.0)
    st = a.geometry_stats(dist, ang, "off-normal", frame_ps=frame_ps)
    assert st["frac_attack_ready_common"] == pytest.approx(
        st["frac_attack_ready"]), (
        "on a run exactly one common-window long the two figures must be equal; "
        "if they are not, every existing row changed meaning")


def test_the_common_window_ignores_everything_past_it():
    """A long tail must not move the comparable figure."""
    a = _asw()
    n_win = 100
    frame_ps = a.COMMON_WINDOW_PS / n_win
    good = np.full(n_win, 3.0)                     # engaged through the window
    tail = np.full(400, 9.0)                       # long, and gone
    ang = np.full(n_win + 400, 5.0)
    st = a.geometry_stats(np.concatenate([good, tail]), ang, "off-normal",
                          frame_ps=frame_ps)
    assert st["frac_attack_ready_common"] == pytest.approx(1.0), (
        "the tail leaked into the common window")
    assert st["frac_attack_ready"] < 0.3, (
        "the full-run figure should be dragged down by the tail")
    # and the row must say which window each was taken over
    assert st["common_window_ps"] == a.COMMON_WINDOW_PS
    assert st["frame_ps"] == pytest.approx(frame_ps)


def test_the_two_figures_are_reported_side_by_side():
    """Either alone is a trap: one is comparable, the other is what happened."""
    a = _asw()
    st = a.geometry_stats(np.full(50, 3.0), np.full(50, 5.0), "off-normal",
                          frame_ps=24.0)
    for c in ("frac_attack_ready", "frac_attack_ready_common",
              "common_window_ps", "frame_ps"):
        assert c in st, f"{c} missing — the window a fraction was taken over "\
                        f"has to travel with it"


def test_extension_refuses_without_a_checkpoint(tmp_path):
    """A restart from prod.gro would reset the thermostat mid-trajectory.

    That produces a settling transient at every join, which reads as the
    molecule moving. Refusing is the only safe answer.
    """
    from shared import gromacs_explicit as ge
    (tmp_path / "prod.tpr").write_bytes(b"not really a tpr")
    with pytest.raises(ge.GromacsError, match="prod.cpt"):
        ge.extend_production(tmp_path, 1000.0)


def test_an_unmeasurable_frame_extends_rather_than_stopping(monkeypatch, tmp_path):
    """Fail-safe direction: stopping early throws the molecule away.

    `_equil_distance` returns None when it cannot be certain (missing frame,
    triclinic box, ambiguous sulfur). Treating None as "left" would silently
    truncate runs for an infrastructure reason.
    """
    a = _asw()
    monkeypatch.setattr(a, "_equil_distance", lambda *args, **kw: None)
    calls = []

    class _FakeGE:
        @staticmethod
        def extend_production(rep, add, gpu_id=None, threads=8):
            calls.append(add)
            return 1200.0 + sum(calls)

    monkeypatch.setattr(a, "_ge", lambda: _FakeGE)
    out = a.adaptive_extend("c", tmp_path, tmp_path / "p.sdf", 1, 0,
                            start_ps=1200.0, max_ps=5200.0, chunk_ps=2000.0,
                            leave_a=6.0)
    assert out["left"] is False
    assert out["total_ps"] == pytest.approx(5200.0)
    assert calls, "an unmeasurable frame stopped the run instead of extending"


def test_a_departed_molecule_stops_immediately(monkeypatch, tmp_path):
    a = _asw()
    monkeypatch.setattr(a, "_equil_distance", lambda *args, **kw: 9.0)

    class _FakeGE:
        @staticmethod
        def extend_production(rep, add, gpu_id=None, threads=8):
            raise AssertionError("extended a molecule that had already left")

    monkeypatch.setattr(a, "_ge", lambda: _FakeGE)
    out = a.adaptive_extend("c", tmp_path, tmp_path / "p.sdf", 1, 0,
                            start_ps=1200.0, max_ps=10000.0, chunk_ps=2000.0,
                            leave_a=6.0)
    assert out["left"] is True
    assert out["left_at_ps"] == pytest.approx(1200.0)
    assert out["extensions"] == 0


def test_a_failed_extension_is_not_reported_as_a_departure(monkeypatch, tmp_path):
    """Inventing a result is worse than a short run."""
    a = _asw()
    monkeypatch.setattr(a, "_equil_distance", lambda *args, **kw: 3.0)

    class _FakeGE:
        @staticmethod
        def extend_production(rep, add, gpu_id=None, threads=8):
            raise RuntimeError("gmx fell over")

    monkeypatch.setattr(a, "_ge", lambda: _FakeGE)
    out = a.adaptive_extend("c", tmp_path, tmp_path / "p.sdf", 1, 0,
                            start_ps=1200.0, max_ps=10000.0, chunk_ps=2000.0,
                            leave_a=6.0)
    assert out["left"] is False, "a GROMACS failure was recorded as the "\
                                 "molecule leaving"
    assert out["total_ps"] == pytest.approx(1200.0)


def test_adaptive_is_off_by_default():
    """Turning it on changes what every fraction is a fraction of."""
    import subprocess
    h = subprocess.run([sys.executable, str(REPO / "scripts/attack_sweep.py"),
                        "--help"], capture_output=True, text=True)
    assert "--adaptive-max-ps" in h.stdout
    assert "0 = off" in h.stdout, "the default must be stated as off"


# --------------------- the analysis files must track the trajectory ---------
def test_the_rmsd_window_is_read_back_not_assumed(tmp_path):
    """Half the gate comes from `rmsd.xvg`; the row must say what it covers.

    THE DEFECT (found 2026-09-03, by @twu383 asking whether the viewer showed
    how long a run held). `rmsd.xvg` is written when `md_residence` returns,
    which is BEFORE any adaptive extension. A 5.2 ns run therefore kept a 1.2 ns
    RMSD trace, so `frac_attack_ready` covered the full run while `rmsd_max_a`
    beside it covered the first quarter -- two windows inside one row, with
    nothing saying so. It also made every adaptive plot stop at 1.2 ns.

    `t4_b00da4134a24_m50` read rmsd_max 0.224 nm on the stale trace and
    0.810 nm over its real 5.2 ns, which flips the pose-held half of the gate.
    """
    a = _asw()
    (tmp_path / "rmsd.xvg").write_text(
        "# c\n@ t\n0.0 0.10\n2.5 0.30\n5.2 0.81\n")     # gmx writes ns
    assert a._rmsd_window_ps(tmp_path) == pytest.approx(5200.0), (
        "the window is not read back from the trace itself")
    assert a._rmsd_window_ps(tmp_path / "nope") is None


def test_the_protein_rmsd_cache_is_keyed_on_currency(tmp_path):
    """A cache keyed on EXISTENCE cannot notice a longer trajectory.

    `protein_rmsd` returned its output whenever the file merely existed, so
    after an extension rewrote `whole.xtc` the protein trace still covered the
    first 1.2 ns while the ligand trace covered the whole run -- one line
    stopping a quarter of the way across the plot, which reads as a deliberate
    choice rather than a stale file.
    """
    import inspect
    from shared import gromacs_analysis as ga
    src = inspect.getsource(ga.protein_rmsd)
    assert "st_mtime" in src, (
        "protein_rmsd does not compare its output against the trajectory; a "
        "stale trace will be served after every extension")
    assert "whole" in src


def test_a_stale_analysis_is_refreshed_after_an_extension():
    """The extension path must redo the analysis, not just the movie."""
    src = (REPO / "scripts" / "attack_sweep.py").read_text()
    i = src.find("adaptive_extend(cand, rep")
    assert i > 0, "the adaptive call site moved"
    seg = src[i:i + 1500]
    assert "gromacs_analysis" in seg and "analyse" in seg, (
        "attack_sweep does not re-analyse after extending; rmsd.xvg and "
        "mindist.xvg would keep covering only the first chunk")
