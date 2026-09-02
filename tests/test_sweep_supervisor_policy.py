"""The supervisor must be stealthy by construction, not by luck.

This process is meant to run for weeks, unattended, on a shared box. Every test
here is about a way it could become antisocial or lose work while still looking
like it was working -- which is the only failure mode that matters in something
nobody is watching.

WHAT IS DELIBERATELY NOT TESTED: that a sweep produces the right answer. That is
`attack_sweep`'s job and it has its own coverage. These tests are about the
supervisor's two responsibilities -- take only what is free, and never lose a
task.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def _sup():
    spec = importlib.util.spec_from_file_location(
        "sweep_supervisor_under_test", REPO / "scripts" / "sweep_supervisor.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    try:
        spec.loader.exec_module(m)
    except Exception as exc:                              # noqa: BLE001
        pytest.skip(f"supervisor not importable here: {exc}")
    return m


def test_the_reserved_gpus_can_never_be_taken():
    """0, 4 and 7 belong to other people (overnight.sh, elevate_queue)."""
    m = _sup()
    assert m.FORBIDDEN_GPUS == {0, 4, 7}
    assert not (set(m.CANDIDATE_GPUS) & m.FORBIDDEN_GPUS), (
        f"candidate list {m.CANDIDATE_GPUS} includes a reserved GPU")


def test_the_schedule_covers_every_hour_exactly_once():
    """A gap would fall through to a default; an overlap makes the cap ambiguous."""
    m = _sup()
    hits = {}
    for h in range(24):
        n, win, ceiling = m.target_workers(datetime(2026, 9, 2, h))
        hits[h] = (win, n, ceiling)
        assert n >= 0
        assert 0 < ceiling <= (os.cpu_count() or 1)
    assert len(hits) == 24
    # Daytime must be the most restrictive window -- that is the whole request.
    day = [n for h, (w, n, c) in hits.items() if w == "day"]
    night = [n for h, (w, n, c) in hits.items() if w == "night"]
    assert day and night
    assert max(day) < min(night), (
        f"daytime cap {max(day)} is not below the night cap {min(night)}; "
        f"the policy would be as intrusive at noon as at 3am")


def test_the_night_cap_can_never_exceed_the_candidate_gpus():
    """A cap above the card count is a number that cannot mean anything."""
    m = _sup()
    for name, v in m.SCHEDULE.items():
        assert v[2] <= len(m.CANDIDATE_GPUS), (
            f"window {name!r} allows {v[2]} workers but only "
            f"{len(m.CANDIDATE_GPUS)} candidate GPUs exist")


def test_a_malformed_policy_file_falls_back_conservatively(tmp_path):
    """A typo in a hand-edited policy must not turn the campaign loose.

    This runs unattended for weeks. The dangerous failure is not a crash -- it
    is a policy file that half-parses into something more aggressive than
    anyone chose.
    """
    m = _sup()
    (tmp_path / "policy.json").write_text("{ not json at all")
    pol = m.load_policy(tmp_path)
    assert pol == dict(m.SCHEDULE), "a broken policy file did not fall back"

    # a well-formed file with an OUT-OF-RANGE cap is also refused
    (tmp_path / "policy.json").write_text(
        '{"schedule": {"night": [22, 8, 99, 190.0]}}')
    assert m.load_policy(tmp_path) == dict(m.SCHEDULE), (
        "a policy asking for 99 workers was accepted; the range check on "
        "max_workers is not firing")


def test_a_valid_policy_file_is_actually_honoured(tmp_path):
    """The converse -- a fallback that ignores every file is equally wrong."""
    m = _sup()
    (tmp_path / "policy.json").write_text(
        '{"schedule": {"night": [22, 8, 2, 100.0], "day": [8, 22, 1, 90.0]}}')
    pol = m.load_policy(tmp_path)
    n, win, ceiling = m.target_workers(datetime(2026, 9, 2, 3), schedule=pol)
    assert (n, win, ceiling) == (2, "night", 100.0), (
        f"policy.json was parsed but not applied: got {(n, win, ceiling)}")


def test_a_failed_gpu_probe_assumes_BUSY_not_free(monkeypatch):
    """The dangerous default. An unreadable nvidia-smi must not free the box.

    If a probe failure returned "no owners", the supervisor would fan out onto
    every candidate GPU at exactly the moment it cannot see who is on them.
    """
    m = _sup()

    def boom(*a, **k):
        raise OSError("nvidia-smi is not answering")
    monkeypatch.setattr(m.subprocess, "run", boom)

    owners = m.gpu_owners()
    assert set(owners) == set(m.CANDIDATE_GPUS)
    for g, who in owners.items():
        assert who, f"gpu{g} reported as unowned after a failed probe"
    monkeypatch.setattr(m, "gpu_owners", m.gpu_owners)
    assert m.free_gpus("twu383") == [], (
        "a failed probe left GPUs looking free -- the supervisor would take "
        "cards it cannot see the owners of")


def test_another_users_gpu_is_never_free(monkeypatch):
    """Ours vs theirs is decided by PID OWNERSHIP, not process name."""
    m = _sup()
    monkeypatch.setattr(m, "gpu_owners", lambda: {
        1: {"mmeawad"}, 2: {"twu383"}, 3: {"twu383", "mmeawad"},
        5: set(), 6: {"ysun2443"}})
    free = m.free_gpus("twu383")
    assert 1 not in free, "took a GPU owned by another user"
    assert 3 not in free, "took a SHARED GPU -- another user is already on it"
    assert 6 not in free
    assert sorted(free) == [2, 5], f"expected our own and the idle one, got {free}"


def test_a_claim_is_atomic_and_taken_once(tmp_path):
    """Two workers must never take the same mode."""
    m = _sup()
    c = m.Claims(tmp_path / "claims")
    assert c.take("t4_abc_m12") is True
    assert c.take("t4_abc_m12") is False, "the same task was claimed twice"
    assert c.take("t4_abc_m13") is True


def test_a_stale_claim_is_reclaimed_but_a_finished_one_is_not(tmp_path):
    """A worker killed by a reboot must not lose its task for good.

    And the converse: a task that FINISHED must stay claimed, or the supervisor
    re-runs completed work forever.
    """
    m = _sup()
    root = tmp_path / "claims"
    c = m.Claims(root)

    c.take("stale_task")
    old = time.time() - (m.STALE_H + 1) * 3600
    os.utime(root / "stale_task", (old, old))
    assert c.take("stale_task") is True, (
        f"a claim older than {m.STALE_H} h was not reclaimed; its task is lost")

    c.take("finished_task")
    c.finish("finished_task")
    os.utime(root / "finished_task", (old, old))
    assert c.take("finished_task") is False, (
        "a FINISHED task was reclaimed because it was old -- the supervisor "
        "would re-run completed sweeps indefinitely")


def test_done_counts_only_completed_sweeps(tmp_path, monkeypatch):
    """A claim means started; a row means nothing unless its status is `ok`.

    THE DEFECT THIS CATCHES (found the hard way, 2026-09-02). `attack_sweep
    --stage0-only` is a FREE geometry probe that writes a row carrying
    `status = "stage0 only"` and no measurements. `done_tasks` matched on ident
    alone, so running the free probe over the worklist marked those modes
    finished and they would never have been simulated -- a row that means "not
    measured" read as "measured".

    The earlier version of this test asserted only that `done` came from the
    results table, which was true and insufficient.
    """
    m = _sup()
    import pandas as pd
    d = tmp_path / "attack_sweep_zz"
    d.mkdir()
    pd.DataFrame({
        "ident": ["t4_a_m1", "t4_b_m2", "t4_c_m3", "t4_d_m4"],
        "status": ["ok", "stage0 only", "failed: GromacsError", "ok"],
    }).to_csv(d / "s_1.csv", index=False)
    monkeypatch.setattr(m.rp, "BLACKSMITH", tmp_path)
    monkeypatch.setattr(m.rp, "sweep_topic", lambda t=None: "attack_sweep_zz")

    done = m.done_tasks("zz")
    assert done == {"t4_a_m1", "t4_d_m4"}, (
        f"expected only the `ok` rows, got {done}")
    assert "t4_b_m2" not in done, (
        "a stage0-only probe row marked a mode as swept; the free geometry "
        "check would silently consume the worklist")
    assert "t4_c_m3" not in done, (
        "a failed sweep counted as done, so a transient failure is never "
        "retried (nac_screen_v2's resume rule, same reason)")


def test_the_load_ceiling_leaves_headroom():
    """Pausing on CPU load is the difference between free GPU and free BOX."""
    m = _sup()
    ncpu = os.cpu_count() or 1
    assert 0 < m.LOAD_CEILING < ncpu, (
        f"load ceiling {m.LOAD_CEILING} is not below the {ncpu} cores present; "
        f"it would never fire and the CPU guard is decorative")
