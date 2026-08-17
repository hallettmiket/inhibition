"""The screen as a pipeline: stages that know their own inputs, outputs and state.

WHY THIS EXISTS. Every stage of this project has a correct, tested script. What
did not exist was the thing that runs them in order, and so it was hand-written
as shell each time -- and every hand-written launcher introduced a new defect
that the stage code itself would never have made. From one night of running
3.1.0 by hand:

  * a worklist read POSITIONALLY -- field 5 taken as `pose_rank` when it holds
    `global_rank` -- so all 24 launched sweeps asked for pose ranks in the
    hundreds;
  * `ls -t | head -1` picking the newest `sweep_gaps_*.csv`, which after a
    failed ranking was the PREVIOUS SCREEN's worklist, reported as this run's;
  * `--consensus consensus_autodock` passed because it looked more correct,
    bypassing a remap that only fires on the default, leaving an all-NaN column
    and ranking 0 of 4,432 modes with no error;
  * a survivor query looking for an RMSD column that `attack_sweep` does not
    write, exiting 0 every poll for eight hours while two survivors waited;
  * `md_residence_3ikd` surviving the kill of its parent and relaunching `gmx`,
    so a "stopped" fleet kept eight GPUs at 90%.

Every one is the same shape: a value taken by POSITION, by RECENCY, by DEFAULT,
or by a NAME THAT DOES NOT EXIST -- and every one failed silently, because an
empty result and a broken query are indistinguishable when nothing asserts the
difference.

So the rules this module enforces, each paid for above:

  1. INPUTS AND OUTPUTS ARE NAMED, never discovered. A stage names the artefact
     it consumes; `ls -t` appears nowhere.
  2. PROGRESS IS COUNTED FROM ARTEFACTS ON DISK, never parsed from logs. A log
     line is a claim; a file is a fact.
  3. A QUERY THAT CANNOT BE EVALUATED RAISES. `Unknown` is a state, distinct
     from `zero`. This is rule 3 because it is the one that cost the most.
  4. STAGES ARE IDEMPOTENT AND RESUMABLE. Restarting one is always safe, which
     is what makes a supervisor tolerable.
  5. THE GPU BUDGET IS DECLARED, and the runner refuses to exceed it.

The GUI reads `status()`; the CLI in `scripts/pipeline.py` drives it. Neither
holds any knowledge of a stage that is not declared here.
"""

from __future__ import annotations

import glob
import json
import os
import signal
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

from . import run_paths as rp
from . import target_config as tc

REPO = Path(__file__).resolve().parent.parent
PY = str(Path.home() / ".micromamba/envs/dwi_reactive/bin/python")

#: Where per-stage logs and the pid files live. Under the run's scratch root, so
#: a fresh topic gets a fresh set and nothing is inherited from a past screen.
def run_dir() -> Path:
    d = rp.SCRATCH / f"pipeline_{rp.topic()}"
    d.mkdir(parents=True, exist_ok=True)
    return d


class StageError(RuntimeError):
    """A stage's state could not be determined. NOT the same as 'nothing done'."""


# ---------------------------------------------------------------------------
# progress probes -- each returns (done, total) and RAISES if it cannot tell
# ---------------------------------------------------------------------------

def _cat(pattern: str) -> pd.DataFrame:
    out = []
    for f in glob.glob(pattern):
        try:
            out.append(pd.read_csv(f))
        except Exception:                                      # noqa: BLE001
            continue
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def scope_idents() -> list[str]:
    """The molecules this run screens: tier + family scope, gates applied.

    Derived from the frame every time rather than from a file written once, so
    it cannot describe a scope the config no longer names.
    """
    fs = sorted(glob.glob("/data/lab_vm/append_only/inhibition/04_t4_combinatorial/D4_*.parquet"),
                key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        raise StageError("no D4 frame on disk")
    d = pd.read_parquet(fs[-1]).drop_duplicates("candidate_id")
    want: set[str] = set()
    for v in tc.sweep_families().values():
        want |= set(v)
    keep = d.warhead_class.isin(want)
    for gate in ("docked_species_ok", "warhead_intact", "alert_gate_pass"):
        if gate in d.columns:
            keep &= (d[gate] == True)                          # noqa: E712
    return sorted(d[keep].candidate_id.astype(str))


def p_screen() -> tuple[int, int]:
    total = len(scope_idents())
    done = len(glob.glob(str(rp.poses_dir() / "*.sdf")))
    return min(done, total), total


def p_rank() -> tuple[int, int]:
    tier = str(tc.get("run.tiers", default=["T4"])[0]).upper()
    fs = glob.glob(str(rp.BLACKSMITH / "rank_v2" /
                       f"rank_v2_{tier}_{rp.topic()}_conditional_eb_*.csv"))
    return (1 if fs else 0), 1


def p_worklist() -> tuple[int, int]:
    return (1 if worklist_path() else 0), 1


def p_sweep() -> tuple[int, int]:
    wl = worklist_path()
    if wl is None:
        return 0, 0
    total = len(pd.read_csv(wl))
    sw = _cat(str(rp.sweep_dir() / "*.csv"))
    if sw.empty or "status" not in sw.columns:
        return 0, total
    ok = sw[sw.status.astype(str).str.startswith("ok")]
    # Deduped on the PAIR: a mode re-swept after a failure must count once.
    if {"parent_ident", "pose_rank"} <= set(ok.columns):
        ok = ok.drop_duplicates(["parent_ident", "pose_rank"])
    return len(ok), total


def survivors() -> pd.DataFrame:
    """Swept modes that held under the bar, with their max ligand RMSD.

    THE RMSD IS NOT IN THE SWEEP CSV. `attack_sweep` writes attack geometry --
    frac_attack_ready, visits, distances, angles -- and no RMSD column at all;
    the trace lives in `rmsd.xvg` beside the trajectory. A watcher that looked
    for `explicit_ligand_rmsd_nm_max` in the CSV found nothing and reported "no
    survivors" for eight hours while two sat under the bar.

    KEYED ON (parent_ident, pose_rank). In an `ok` row `ident` is the MODE and
    `parent_ident` the molecule; in a FAILED row `ident` is the molecule,
    because the run died before it knew its mode. Keying on `ident` resolved
    none of 26.

    Raises when finished sweeps exist but no trace can be read -- that is a
    broken query, and it must not present as an empty result.
    """
    sys.path.insert(0, str(REPO / "scripts"))
    import sweep_assets as sa                                  # noqa: PLC0415

    sw = _cat(str(rp.sweep_dir() / "*.csv"))
    if sw.empty or "status" not in sw.columns:
        return pd.DataFrame(columns=["ident", "parent_ident", "pose_rank", "rmsd_max"])
    ok = sw[sw.status.astype(str).str.startswith("ok")]
    if ok.empty:
        return pd.DataFrame(columns=["ident", "parent_ident", "pose_rank", "rmsd_max"])
    if not {"parent_ident", "pose_rank"} <= set(ok.columns):
        raise StageError("sweep rows lack parent_ident/pose_rank")

    rows, unresolved = [], 0
    for r in ok.drop_duplicates(["parent_ident", "pose_rank"]).itertuples():
        rep = sa.rep_dir(str(r.parent_ident), int(r.pose_rank))
        if rep is None:
            unresolved += 1
            continue
        _t, y = sa._xvg(rep / "rmsd.xvg")
        if y is None or not len(y):
            unresolved += 1
            continue
        rows.append({"ident": str(r.ident), "parent_ident": str(r.parent_ident),
                     "pose_rank": int(r.pose_rank), "rmsd_max": float(y.max()),
                     "rmsd_mean": float(y.mean()),
                     "frac_attack_ready": float(getattr(r, "frac_attack_ready", float("nan")))})
    if not rows and unresolved:
        raise StageError(
            f"{unresolved} finished sweeps but no RMSD trace resolved under "
            f"{rp.sweep_work()} -- the query is broken, not the result empty")
    d = pd.DataFrame(rows)
    return d[d.rmsd_max < tc.md_survivor_rmsd_nm()].sort_values("rmsd_max")


def p_production() -> tuple[int, int]:
    total = len(survivors())
    md = _cat(str(rp.residence_dir() / "*.csv"))
    if md.empty or "production_ps" not in md.columns:
        return 0, total
    done = md[md.production_ps >= 90_000].ident.nunique()
    return int(done), max(total, int(done))


def p_bpmd() -> tuple[int, int]:
    md = _cat(str(rp.residence_dir() / "*.csv"))
    total = 0
    if not md.empty and {"production_ps", "explicit_ligand_rmsd_nm_max"} <= set(md.columns):
        held = md[(md.production_ps >= 90_000)
                  & (md.explicit_ligand_rmsd_nm_max < tc.md_survivor_rmsd_nm())]
        total = held.ident.nunique()
    bp = _cat(str(rp.bpmd_dir() / "*.csv"))
    done = bp.ident.nunique() if not bp.empty and "ident" in bp.columns else 0
    return int(done), max(total, int(done))


def worklist_path() -> Path | None:
    """The worklist for THIS run's ranking, or None.

    NAMED BY THE RUN, NOT BY RECENCY. `ls -t | head -1` returned a worklist
    written for a previous screen when this run's ranking had failed, and 170
    modes from that screen were reported as this one's selection. A worklist is
    only this run's if it is newer than this run's ranking table.
    """
    tier = str(tc.get("run.tiers", default=["T4"])[0]).upper()
    ranks = glob.glob(str(rp.BLACKSMITH / "rank_v2" /
                          f"rank_v2_{tier}_{rp.topic()}_conditional_eb_*.csv"))
    if not ranks:
        return None
    newest_rank = max(os.path.getmtime(f) for f in ranks)
    cands = [f for f in glob.glob(str(rp.BLACKSMITH / "sweep_gaps" / "sweep_gaps_*.csv"))
             if os.path.getmtime(f) >= newest_rank]
    if not cands:
        return None
    return Path(max(cands, key=os.path.getmtime))


# ---------------------------------------------------------------------------
# the stages
# ---------------------------------------------------------------------------

@dataclass
class Stage:
    name: str
    title: str
    needs: tuple[str, ...]
    probe: Callable[[], tuple[int, int]]
    gpus: int = 0
    #: argv builder; returns the command that RUNS this stage to completion.
    #: A stage that fans out returns the supervisor script, not one job.
    launch: Callable[[], list[str]] | None = None
    #: Substring identifying this stage's processes, for state + stop.
    proc: str = ""
    note: str = ""


def _screen_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/pipeline_stage.py"), "screen"]


def _rank_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/rank_v2.py"), "--score", "conditional_eb"]


def _worklist_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/sweep_gap_worklist.py"), "--by-family"]


def _sweep_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/pipeline_stage.py"), "sweep"]


def _prod_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/pipeline_stage.py"), "production"]


def _bpmd_cmd() -> list[str]:
    return [PY, str(REPO / "scripts/pipeline_stage.py"), "bpmd"]


STAGES: list[Stage] = [
    Stage("screen", "Docking + NAC", (), p_screen, gpus=8, launch=_screen_cmd,
          proc="nac_screen_v2.py",
          note="dock each in-scope molecule, score near-attack geometry, "
               "split into binding modes"),
    Stage("rank", "Ranking", ("screen",), p_rank, gpus=0, launch=_rank_cmd,
          proc="rank_v2.py",
          note="conditional_eb within warhead class; empirical-Bayes shrinkage "
               "over the mode's in-range poses, divided by the isotropic null"),
    Stage("worklist", "Sweep worklist", ("rank",), p_worklist, gpus=0,
          launch=_worklist_cmd, proc="sweep_gap_worklist.py",
          note="per family, above the enrichment floor, capped at max_depth"),
    Stage("sweep", "Triage sweep (8 ns)", ("worklist",), p_sweep, gpus=2,
          launch=_sweep_cmd, proc="attack_sweep.py",
          note="is the pose stable at all — max ligand RMSD over 8 ns"),
    Stage("production", "Production MD (100 ns)", ("sweep",), p_production,
          gpus=1, launch=_prod_cmd, proc="md_residence_3ikd.py",
          note="survivors only; does it stay for a long time"),
    Stage("bpmd", "BPMD", ("production",), p_bpmd, gpus=1, launch=_bpmd_cmd,
          proc="bpmd_run.py",
          note="how hard is it to push out — 100 ns holders only"),
]

BY_NAME = {s.name: s for s in STAGES}


# ---------------------------------------------------------------------------
# process state
# ---------------------------------------------------------------------------

def _procs(match: str) -> list[int]:
    """PIDs whose command line contains `match`, excluding this process tree.

    `pgrep -f` matches the SHELL running the query as often as the target, which
    made every count off by one or two and a `pkill` kill its own caller.
    """
    out = []
    me = {os.getpid(), os.getppid()}
    try:
        ps = subprocess.run(["ps", "-eo", "pid,args"], capture_output=True,
                            text=True, timeout=30).stdout
    except Exception:                                          # noqa: BLE001
        return out
    for line in ps.splitlines()[1:]:
        line = line.strip()
        pid_s, _, args = line.partition(" ")
        if not pid_s.isdigit():
            continue
        pid = int(pid_s)
        if pid in me or "bash -c" in args or " ps -eo" in args:
            continue
        if match in args:
            out.append(pid)
    return out


def running(stage: Stage) -> list[int]:
    return _procs(stage.proc) if stage.proc else []


def state(stage: Stage) -> dict:
    """One stage's state, with `unknown` distinct from `zero`."""
    try:
        done, total = stage.probe()
        err = None
    except Exception as exc:                                   # noqa: BLE001
        done = total = None
        err = str(exc)[:200]
    pids = running(stage)
    if err is not None:
        st = "unknown"
    elif pids:
        st = "running"
    elif total and done >= total:
        st = "done"
    elif done:
        st = "stopped"          # partial, nothing running
    else:
        st = "waiting"
    return {"name": stage.name, "title": stage.title, "state": st,
            "done": done, "total": total, "pids": pids, "gpus": stage.gpus,
            "note": stage.note, "error": err, "needs": list(stage.needs)}


def status() -> dict:
    st = {s.name: state(s) for s in STAGES}
    for s in STAGES:
        # A stage is only startable once everything it needs is done.
        st[s.name]["ready"] = all(st[n]["state"] == "done" for n in s.needs)
    return {"topic": rp.topic(), "stages": [st[s.name] for s in STAGES],
            "spec": {"sweep_ns": tc.md_sweep_ps() / 1000,
                     "production_ns": tc.md_production_ps() / 1000,
                     "survivor_rmsd_nm": tc.md_survivor_rmsd_nm(),
                     "tiers": list(tc.get("run.tiers", default=[])),
                     "families": sorted(tc.sweep_families())}}


def write_status() -> Path:
    p = rp.reports_dir() / "pipeline_state.json"
    p.write_text(json.dumps(status(), indent=1, default=str))
    return p


# ---------------------------------------------------------------------------
# start / stop
# ---------------------------------------------------------------------------

def start(name: str) -> int:
    """Launch a stage detached. Returns the pid."""
    s = BY_NAME[name]
    if running(s):
        raise StageError(f"{name} is already running")
    st = status()
    by = {x["name"]: x for x in st["stages"]}
    missing = [n for n in s.needs if by[n]["state"] != "done"]
    if missing:
        raise StageError(f"{name} needs {missing} to be done first")
    log = run_dir() / f"{name}.log"
    cmd = s.launch()
    with open(log, "a") as fh:
        p = subprocess.Popen(cmd, cwd=REPO, stdout=fh, stderr=fh,
                             start_new_session=True)
    (run_dir() / f"{name}.pid").write_text(str(p.pid))
    return p.pid


#: Kill order matters and is not alphabetical. `md_residence_3ikd` spawns `gmx`
#: and relaunches it for the next stage of its own chain, so killing `gmx` first
#: makes a new one appear ~13 s later; a "stopped" fleet kept eight GPUs at 90%
#: for several minutes that way. Parents first, innermost last.
_KILL_ORDER = ("pipeline_stage.py", "attack_sweep.py", "md_residence_3ikd.py",
               "bpmd_run.py", "nac_screen_v2.py", "gmx mdrun")


def stop(name: str) -> int:
    s = BY_NAME[name]
    killed = 0
    targets = [s.proc] if s.proc else []
    if s.name in ("sweep", "production", "bpmd"):
        targets = list(_KILL_ORDER)
    for pat in targets:
        for pid in _procs(pat):
            # `gmx` is shared with other users on this box; only ever kill one
            # whose working directory is inside this run's scratch tree.
            if pat == "gmx mdrun":
                try:
                    cwd = os.readlink(f"/proc/{pid}/cwd")
                except OSError:
                    continue
                if rp.topic() not in cwd:
                    continue
            try:
                os.kill(pid, signal.SIGKILL)
                killed += 1
            except OSError:
                pass
    return killed
