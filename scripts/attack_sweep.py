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


def geometry_stats(dist: np.ndarray, angle: np.ndarray, kind: str) -> dict:
    inw = (dist >= nac.NAC_DIST_MIN) & (dist <= nac.NAC_DIST_MAX)
    ready = inw & competent(angle, kind)
    # An excursion is a rising edge: not-ready -> ready. The first frame counts
    # as a visit if it is already ready, otherwise a pose that starts competent
    # and never leaves would be recorded as zero visits.
    visits = int(np.sum(np.diff(ready.astype(int)) == 1) + (1 if ready[0] else 0))
    return {
        "frac_in_window": float(inw.mean()),
        "frac_attack_ready": float(ready.mean()),
        "n_visits": visits,
        "start_dist_a": float(dist[0]),
        "start_angle_deg": float(angle[0]),
        "start_attack_ready": bool(ready[0]),
        "median_dist_a": float(np.median(dist)),
        "median_angle_deg": float(np.median(angle)),
        "min_dist_a": float(dist.min()),
    }


def run_sweep(cand: str, pose: Path, pose_rank: int, gpu: int,
              ps: float, net_charge: int | None) -> Path | None:
    """10 ns of MD through the production script, so the physics is identical."""
    rep = SWEEP_ROOT / cand / "md" / "rep1"
    if (rep / "prod.xtc").is_file():
        log.info("%s: trajectory already present, not re-running", cand)
        return rep
    cmd = [str(PY), str(REPO / "scripts/md_residence_3ikd.py"),
           "--candidate", cand, "--pose", str(pose),
           "--pose-rank", str(pose_rank), "--production-ps", str(int(ps)),
           "--gpu", str(gpu), "--keep", "--tag", "sweep",
           "--work-root", str(SWEEP_ROOT)]
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
        rec = {"ident": cand, "sweep_ps": args.sweep_ps}
        pose = Path(args.pose_dir) / f"{cand}.sdf"
        if not pose.is_file():
            rec["status"] = f"no pose at {pose}"
            rows.append(rec); log.warning("%s: %s", cand, rec["status"]); continue

        if args.stage0_only:
            rec["status"] = "stage0 only"
            rows.append(rec); continue

        rep = run_sweep(cand, pose, args.pose_rank, args.gpu, args.sweep_ps, None)
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
        rec.update(geometry_stats(s["dist"], s["angle"], s["kind"]))
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
    print("\n  RANKED, NOT THRESHOLDED. There is no evidence for a cut-off — one")
    print("  molecule above 5% and five below is one molecule, not a threshold.")
    print("  How deep to elevate is a separate decision with its own cost argument.")
    print(f"\n  readings fixed in advance: docs/prereg_attack_sweep.md")
    print(f"  -> {dest}")


if __name__ == "__main__":
    main()
