"""BPMD must run the pose it was asked for, and must not think it already has.

THE DEFECT (D0105). `read_pose(sdf, pose_rank)` has selected by the `pose_rank`
PROPERTY since it was written, `run_pose` accepts the argument, `prepare_pose`
honours it, and the workdir is named `<stem>__p<rank>` for ranks other than 1.
The whole mechanism existed. Both call sites in `main()` simply never passed it,
and there was no CLI flag to supply one -- so every BPMD replicate this project
has run was on pose_rank 1, whatever pose the ranking had actually chosen.

Nothing could have caught it by looking at output. pose_rank 1 is a real pose of
the right molecule; it parameterises, biases and reports an ordinary stability
score. The mode representatives are multi-pose (165 poses in one file for
`t4_80fbed3bdf1e`), so the wrong one is always available and always plausible.

THE SECOND HALF is the resume key. `already_done()` returned `(ident,
replicate)`, so a finished run of pose_rank 1 marked pose_rank 11 as done -- the
molecule matches, the replicate matches -- and the second pose would never be
simulated while the table said it had been. Adding the flag WITHOUT fixing the
key installs a worse bug than it removes. That is `how_this_project_breaks` #22
(this same function, keyed on trajectory length) meeting #23 (an asset matched
on the molecule where the runner writes one per (molecule, pose_rank)).

WHAT THESE TESTS ARE CAREFUL NOT TO DO. Asserting that `--pose-rank` parses
would pass while the value went nowhere -- the vacuous-guard trap this repo has
already been caught by twice (`how_this_project_breaks`, disguise #4). So the
first test drives `main()` with the runner replaced and asserts the value
ARRIVES; the second asserts the resume key separates two poses of one molecule.
Both fail if the plumbing is removed.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]


def _load():
    """Import `scripts/bpmd_run.py` by path -- it is a script, not a package."""
    spec = importlib.util.spec_from_file_location(
        "bpmd_run_under_test", REPO / "scripts" / "bpmd_run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:                                  # noqa: BLE001
        pytest.skip(f"bpmd_run is not importable in this env: {exc}")
    return mod


def test_pose_rank_reaches_the_runner(monkeypatch, tmp_path):
    """`--pose-rank 11` must arrive at `run_pose`, not stop at argparse.

    THE ASSERTION IS ON THE VALUE THE RUNNER RECEIVED. A test that only checked
    `args.pose_rank == 11` would have passed throughout the entire period the
    defect existed, because argparse was never the broken part.
    """
    m = _load()
    seen = {}

    def fake_run_pose(cand, **kw):
        seen.update(kw)
        seen["ident"] = cand.ident
        return []

    class Cand:
        ident = "t4_80fbed3bdf1e"
        warhead_class = "bdhi_c4"
        mechanism = "sn2_ring_opening"
        label = "T_4"

    monkeypatch.setattr(m, "run_pose", fake_run_pose)
    monkeypatch.setattr(m, "candidate_index", lambda: {"t4_80fbed3bdf1e": Cand()})
    monkeypatch.setattr(m, "already_done", lambda: set())
    monkeypatch.setattr(m, "report", lambda: None)
    monkeypatch.setattr(m.gx, "plumed_kernel", lambda: "stub")
    monkeypatch.setattr(m, "set_poses_dir", lambda d: None, raising=False)

    class _Chunk:
        def __init__(self, *a, **k): pass
        def add(self, row): pass
        def flush(self): pass
    monkeypatch.setattr(m, "_ChunkWriter", _Chunk)

    monkeypatch.setattr(sys, "argv", [
        "bpmd_run.py", "--pose", "t4_80fbed3bdf1e", "--pose-rank", "11",
        "--replicates", "3", "--production-ps", "10000", "--gpu", "3",
        "--no-redock"])
    m.main()

    assert seen, "run_pose was never called — the test drove nothing"
    assert seen.get("pose_rank") == 11, (
        "run_pose received pose_rank=%r, not 11. The CLI value is not reaching "
        "the runner, so BPMD is simulating a pose nobody asked for (D0105)."
        % seen.get("pose_rank"))


def test_pose_rank_is_part_of_the_resume_key(monkeypatch, tmp_path):
    """A finished pose_rank 1 must NOT mark pose_rank 11 as already done."""
    m = _load()
    chunk = tmp_path / "bpmd_s0_1.csv"
    pd.DataFrame([
        {"ident": "t4_80fbed3bdf1e", "pose_rank": 1, "replicate": 1, "status": "ok"},
        {"ident": "t4_80fbed3bdf1e", "pose_rank": 1, "replicate": 2, "status": "ok"},
    ]).to_csv(chunk, index=False)

    monkeypatch.setattr(m, "_chunk_files", lambda d, pat: [chunk])
    done = m.already_done()

    assert ("t4_80fbed3bdf1e", 1, 1) in done, "the finished pose is not recorded"
    assert ("t4_80fbed3bdf1e", 11, 1) not in done, (
        "pose_rank 11 is reported as done on the strength of pose_rank 1's run. "
        "The resume key does not distinguish two poses of one molecule, so the "
        "pose the ranking actually chose would be silently skipped (D0105).")


def test_rows_without_pose_rank_read_as_rank_1(monkeypatch, tmp_path):
    """Pre-flag rows have no `pose_rank` column; they were all rank 1.

    Not a compatibility shim for its own sake: reading them as anything else --
    or dropping them -- would re-run six replicates that are already on disk.
    """
    m = _load()
    chunk = tmp_path / "bpmd_s0_1.csv"
    pd.DataFrame([
        {"ident": "t3_00055649a545", "replicate": 1, "status": "ok"},
    ]).to_csv(chunk, index=False)

    monkeypatch.setattr(m, "_chunk_files", lambda d, pat: [chunk])
    assert ("t3_00055649a545", 1, 1) in m.already_done()
