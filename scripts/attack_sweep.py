"""
Purpose: the cheap attack-geometry gate — free pose check, then a 10 ns sweep.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: candidates + their elevated poses
Output: 00_outputs/blacksmith/attack_sweep/attack_sweep_<N>.csv, ranked

Implements #32 and `docs/prereg_attack_sweep.md`. @tt8804: *"100% residence can
still result in a warhead that does not stay near the attack distance or angle...
we should do 10 ns sweeps to first filter for that ability to start at cys113 and
stay near before comitting to 100 ns MD."*

WHY THE FUNNEL NEEDED THIS. Residence and attack geometry are nearly independent
across the six molecules measured so far -- rho = +0.319, p = 0.538. The two with
PERFECT residence are attack-ready 0.4% and 1.7% of frames; the one that
dissociated at 81 ns is attack-ready 55.2%. Ranking on residence selects against
the thing we want.

TWO STAGES, AND THE FIRST IS FREE.

  STAGE 0  the elevated pose's OWN geometry, before any simulation. Only 1 of 6
           poses elevated in the bornite cohort started inside the attack window
           at all, and it is the only one that got anywhere. If that separation
           holds it is a filter that costs nothing, and no sweep should be run
           on a pose that starts 9 A away pointing the wrong direction.

  STAGE 1  10 ns of unbiased MD, then attack geometry per frame. ~0.4 GPU-h
           against ~4 for a full run, so it pays for itself on one rejection.

BOTH OBSERVABLES ARE REPORTED, BECAUSE THEY MAY DIVERGE.
  frac_attack_ready  how LONG the molecule is competent
  n_visits           how MANY independent excursions it makes into attack
                     geometry -- and a covalent reaction needs ONE good approach,
                     not sustained occupancy, so this is the more mechanistically
                     honest quantity. On the six existing trajectories the two
                     rank identically (rho = +1.000); that may not survive a
                     cohort with more range, which is why both are kept.

THIS RANKS, IT DOES NOT THRESHOLD. There is no evidence for a cut-off value --
one molecule above 5% and five below is not a threshold, it is one molecule. The
output is an ordering; how deep to elevate is a separate decision with its own
cost argument.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import md_movie as mov                  # noqa: E402
from shared import nac_criterion as nac             # noqa: E402
from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("attack-sweep")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
#: The sweep writes to its OWN root. The workdir is <root>/<candidate>/md/rep1
#: regardless of tag, so a 10 ns sweep and a later 100 ns run of the same
#: molecule would otherwise collide -- and the 100 ns run would find a finished
#: 10 ns prod.xtc sitting there and skip itself.
SWEEP_ROOT = Path("/data/lab_vm/modifiable/inhibition/attack_sweep_10ns")
POSES = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
OUT = sout.Topic("blacksmith", "attack_sweep")
PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"

SWEEP_PS = 10_000.0        # 10 ns
FRAMES = 500               # 20 ps resolution over the sweep

#: A visit has to LAST to count as a visit (#34, adversary audit).
#:
#: Measured on the six 100 ns trajectories, the median attack-ready episode is
#: one or two frames for five of the six -- 100% of two molecules' episodes are
#: <=200 ps. Those are boundary recrossings, not resolved approaches. With a mean
#: episode of ~1 frame the visit count is algebraically almost
#: `frac_attack_ready x n_frames`, which is why the two ranked identically
#: (rho = +1.000) and why calling them separate observables was not defensible.
#:
#: IN PICOSECONDS, NOT FRAMES, because a frame rule makes the count a function of
#: the save interval: one molecule gave 57/26/14/7/3 visits at
#: 100 ps/200 ps/500 ps/1 ns/2 ns. The sweep saves every 20 ps and the validation
#: ran at 100 ps, so a frame-based count would not be comparable between them.
MIN_DWELL_PS = 100.0


def _mp():
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = m
    spec.loader.exec_module(m)
    return m


def competent(angle: np.ndarray, kind: str) -> np.ndarray:
    """Angular competence by the mechanism's own criterion, never one constant."""
    if "anti" in (kind or ""):
        return angle >= nac.SN2_ANGLE_MIN
    return angle <= nac.PERPENDICULAR_MAX_OFF_NORMAL


def _episodes(ready: np.ndarray) -> list[tuple[int, int]]:
    """(start, length) of every contiguous run of True."""
    out, i, n = [], 0, len(ready)
    while i < n:
        if ready[i]:
            j = i
            while j < n and ready[j]:
                j += 1
            out.append((i, j - i))
            i = j
        else:
            i += 1
    return out


def geometry_stats(dist: np.ndarray, angle: np.ndarray, kind: str,
                   frame_ps: float) -> dict:
    """Geometry readings for one trajectory.

    `frame_ps` is REQUIRED rather than defaulted: `n_visits` is meaningless
    without it (see MIN_DWELL_PS), and a default would silently produce a number
    that is not comparable to the one it is being validated against.
    """
    if not len(dist):
        # nac_criterion raises rather than returning a false verdict when it
        # cannot measure; this did the opposite and raised IndexError on
        # `ready[0]` from deep inside a worker, where it reads as a crash rather
        # than as "there is no trajectory here".
        raise ValueError("empty trajectory: no frames to score")
    inw = (dist >= nac.NAC_DIST_MIN) & (dist <= nac.NAC_DIST_MAX)
    ready = inw & competent(angle, kind)
    # An excursion is a rising edge: not-ready -> ready. The first frame counts
    # as a visit if it is already ready, otherwise a pose that starts competent
    # and never leaves would be recorded as zero visits.
    raw_visits = int(np.sum(np.diff(ready.astype(int)) == 1) + (1 if ready[0] else 0))
    min_frames = max(1, int(round(MIN_DWELL_PS / frame_ps)))
    eps = _episodes(ready)
    visits = sum(1 for _, ln in eps if ln >= min_frames)
    return {
        "frac_in_window": float(inw.mean()),
        "frac_attack_ready": float(ready.mean()),
        "n_visits": visits,
        # Both are kept so the debounce is inspectable rather than assumed. A
        # large gap between them means the molecule is skimming the boundary,
        # not approaching.
        "n_visits_raw": raw_visits,
        "min_dwell_ps": MIN_DWELL_PS,
        "frame_ps": float(frame_ps),
        "median_episode_ps": float(np.median([ln for _, ln in eps]) * frame_ps)
                             if eps else 0.0,
        "start_dist_a": float(dist[0]),
        "start_angle_deg": float(angle[0]),
        "start_attack_ready": bool(ready[0]),
        "median_dist_a": float(np.median(dist)),
        "median_angle_deg": float(np.median(angle)),
        "min_dist_a": float(dist.min()),
    }


def pose_mode(pose_sdf: Path, pose_rank: int) -> int | None:
    """The `mode` property of the pose with this `pose_rank`, or None.

    Read by identity: the pose whose own `pose_rank` matches, and then that
    pose's own `mode`. Nothing here counts positions in the file.
    """
    if not pose_sdf.is_file():
        return None
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=False):
            if m is None or not m.HasProp("pose_rank"):
                continue
            if int(m.GetProp("pose_rank")) == pose_rank:
                return int(m.GetProp("mode")) if m.HasProp("mode") else None
    except Exception as exc:                                # noqa: BLE001
        log.warning("%s: could not read mode (%s)", pose_sdf.name, exc)
    return None


class SweepError(RuntimeError):
    """A sweep could not be run or reused. Named so it cannot pass as a result."""


def _finished(rep: Path) -> bool:
    """True only if mdrun reported it finished this trajectory.

    Existence of `prod.xtc` is not completion -- see run_sweep. GROMACS writes
    "Finished mdrun" to the log on a clean exit, and that is the only marker
    here that distinguishes a finished run from one in progress.
    """
    log_f, xtc = rep / "prod.log", rep / "prod.xtc"
    if not (log_f.is_file() and xtc.is_file()):
        return False
    try:
        return "Finished mdrun" in log_f.read_text(errors="replace")
    except OSError:
        return False


def run_sweep(cand: str, pose: Path, pose_rank: int, gpu: int,
              ps: float, net_charge: int | None) -> Path | None:
    """10 ns of MD through the production script, so the physics is identical."""
    # ONE WORKDIR PER (MOLECULE, MODE). md_residence names the directory after
    # the candidate, so sweeping two modes of one molecule would put them in the
    # same place and the second would find the first's finished trajectory and
    # skip itself -- reporting mode 0's result as mode 4's.
    # THE LENGTH IS PART OF THE PATH, not just the tag.
    #
    # The resume guard checked only that `prod.xtc` existed, so a finished 200 ps
    # trajectory silently satisfied a 10 ns request -- caught tonight when a
    # smoke test's 200 ps run was reused as if it were the real sweep. That is
    # the same defect class the per-mode workdir fixed for pose_rank, with the
    # length dimension left open. A 10 ns answer read off 200 ps of dynamics
    # would be indistinguishable from a real one in every artefact downstream.
    #
    # AND THE GUARD MUST CHECK COMPLETION, NOT EXISTENCE. `prod.xtc` exists the
    # moment mdrun starts writing it, so an IN-PROGRESS run satisfied this test
    # and a second process read the partial trajectory as a finished one. Caught
    # 2026-08-11 the first time two workers were run: `t4_e0b03662d460_m1` was
    # 3.4 ns into its 10 ns when a second invocation analysed the same directory
    # and wrote `status: ok, frac_attack_ready 0.0` in four seconds. A partial
    # answer stamped ok is indistinguishable from a real one downstream -- the
    # same defect this docstring already records for the 200 ps case, with the
    # completeness dimension left open instead of the length one.
    root = SWEEP_ROOT / f"rank{pose_rank}_{int(ps)}ps"
    rep = root / cand / "md" / "rep1"
    if _finished(rep):
        log.info("%s: %d ps trajectory already complete, not re-running",
                 cand, int(ps))
        return rep
    if (rep / "prod.xtc").is_file():
        # Present but unfinished: either another process is running it right now,
        # or one died partway. Neither is a result. Refuse rather than analyse it.
        raise SweepError(
            f"{cand}: an unfinished {int(ps)} ps trajectory is already in "
            f"{rep} — another worker may be running it. Not analysing a partial "
            f"run as a complete one.")
    cmd = [str(PY), str(REPO / "scripts/md_residence_3ikd.py"),
           "--candidate", cand, "--pose", str(pose),
           "--pose-rank", str(pose_rank), "--production-ps", str(int(ps)),
           "--gpu", str(gpu), "--keep", "--tag", f"sweep_r{pose_rank}",
           "--work-root", str(root)]
    if net_charge is not None:
        cmd += ["--net-charge", str(net_charge)]
    log.info("%s: %.0f ps sweep on GPU %d", cand, ps, gpu)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
    if r.returncode != 0:
        log.warning("%s: sweep failed rc=%d %s", cand, r.returncode, r.stderr[-300:])
        return None
    return rep if (rep / "prod.xtc").is_file() else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--pose-dir", default=str(POSES / "nac_v3_poses"))
    ap.add_argument("--pose-rank", type=int, default=1)
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--sweep-ps", type=float, default=SWEEP_PS)
    ap.add_argument("--stage0-only", action="store_true",
                    help="report starting geometry only — costs no GPU at all")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    mp = _mp()

    rows = []
    for cand in args.candidates:
        pose = Path(args.pose_dir) / f"{cand}.sdf"
        # THE MODE IS READ FROM THE POSE, NOT DERIVED FROM ITS RANK.
        #
        # This was `mode = pose_rank - 1`, which assumes the exported poses are
        # in mode order. `export_nac_poses` sorts NAC-viable first by energy and
        # then the rest by energy, so that ordering is not guaranteed by
        # construction -- it happens to hold for all 1,751 poses on disk today,
        # which is exactly the kind of invariant that holds until it does not.
        # The SDF carries an explicit `mode` property; read it (#53).
        mode = pose_mode(pose, args.pose_rank)
        # `ident` carries the MODE, so two modes of one molecule are two rows
        # rather than one overwriting the other downstream.
        #
        # ALWAYS `_m<mode>`, INCLUDING MODE 0. It used to write the bare ident
        # for pose_rank 1, so mode 0 was `t4_x` in this table and `t4_x_m0` in
        # the rank table. Every obvious join on `ident` therefore dropped exactly
        # the modes that were simulated, which is how #53 stayed invisible and
        # how a wrong correlation reached #36. Join on (parent_ident, mode).
        rec = {"ident": f"{cand}_m{mode}" if mode is not None else cand,
               "parent_ident": cand, "pose_rank": args.pose_rank,
               "mode": mode, "sweep_ps": args.sweep_ps}
        if not pose.is_file():
            rec["status"] = f"no pose at {pose}"
            rows.append(rec); log.warning("%s: %s", cand, rec["status"]); continue

        if args.stage0_only:
            rec["status"] = "stage0 only"
            rows.append(rec); continue

        try:
            rep = run_sweep(cand, pose, args.pose_rank, args.gpu, args.sweep_ps, None)
        except SweepError as exc:
            # Recorded, never dropped, and never as a number: the row says why
            # there is no reading rather than carrying a value from a partial run.
            rec["status"] = f"skipped: {exc}"
            rows.append(rec); log.warning("%s: %s", cand, exc); continue
        if rep is None:
            rec["status"] = "sweep failed"
            rows.append(rec); continue

        dense = rep / "sweep_dense.pdb"
        if not dense.is_file():
            mov.build_movie_pdb(rep, dense, n_frames=FRAMES)
        s = mp.nac_series(cand, rep, dense)
        if s is None:
            rec["status"] = "no attack-geometry series"
            rows.append(rec); continue

        rec["status"] = "ok"
        rec["mechanism"] = s["mechanism"]
        # The frame spacing is DERIVED from the run, not assumed: the movie is
        # built with n_frames=FRAMES over `ps` picoseconds, and `n_visits` is
        # meaningless without it (see MIN_DWELL_PS). Taking the actual series
        # length rather than FRAMES, because build_movie_pdb can return fewer
        # frames than asked for and a stale divisor would silently rescale every
        # dwell time.
        n_fr = max(1, len(s["dist"]))
        rec.update(geometry_stats(s["dist"], s["angle"], s["kind"],
                                  frame_ps=float(args.sweep_ps) / n_fr))
        rows.append(rec)
        log.info("%s: attack-ready %.3f over %d visits (start ready=%s)",
                 cand, rec["frac_attack_ready"], rec["n_visits"],
                 rec["start_attack_ready"])

    t = pd.DataFrame(rows)
    ok = t[t.status == "ok"] if "status" in t.columns else t
    if not ok.empty:
        t = t.sort_values("frac_attack_ready", ascending=False, na_position="last")
    dest = OUT.write("attack_sweep", ".csv")
    t.to_csv(dest, index=False)

    print(f"\nAttack-geometry sweep — {args.sweep_ps/1000:.0f} ns, ranked\n")
    if ok.empty:
        print("  no molecules completed"); print(f"\n  -> {dest}"); return
    cols = ["ident", "start_attack_ready", "frac_in_window",
            "frac_attack_ready", "n_visits", "min_dist_a"]
    print(t[[c for c in cols if c in t.columns]].round(4).to_string(index=False))
    # REJECT-ONLY, NOT A RANKING (#34).
    #
    # Agreement between the 10 ns sweep and the 100 ns run is rho = +0.60 on
    # `frac_attack_ready` -- the quantity printed here -- against +0.83 for mere
    # proximity, which is the number older docs quote. The pre-registered reading
    # table licenses +0.60 for discarding the bottom and explicitly NOT for
    # ordering the middle, and this printout used to announce the opposite.
    #
    # The order below is still by attack-readiness, because something has to be
    # printed in some order and the elevation queue has to take a best-first
    # view when capacity binds. What has changed is the claim attached to it.
    print("\n  ORDERED, NOT RANKED. rho(10 ns, 100 ns) = +0.60 on this reading —")
    print("  enough to reject the bottom, NOT enough to order the middle. Treat")
    print("  the ordering as a queue, not as a statement that #2 beats #5.")
    print("  A survivor is a molecule with a SUSTAINED episode (n_visits > 0);")
    print("  a single 20 ps touch is not evidence of anything.")
    print(f"\n  readings fixed in advance: docs/prereg_attack_sweep.md")
    print(f"  -> {dest}")


if __name__ == "__main__":
    main()
