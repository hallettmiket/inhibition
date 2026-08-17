"""
Purpose: #32 — would a 10 ns sweep have predicted the 100 ns attack geometry?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: the completed 100 ns MD-priority trajectories
Output: 00_outputs/blacksmith/mdprio_reports/attack_sweep_<N>.csv

@tt8804 in #32: *"100% residence can still result in a warhead that does not stay
near the attack distance or angle... we should do 10 ns sweeps to first filter for
that ability to start at cys113 and stay near before comitting to 100 ns MD."*

The observation is already confirmed: both molecules that held for the full 100 ns
are inside the attack window AND angle-competent in under 1% of frames. Residence
and reaction competence are different readings.

THE PROPOSAL IS A CHEAP PRE-FILTER, AND ITS ONLY QUESTION IS WHETHER 10 ns SEES
WHAT 100 ns SEES. That is answerable without running anything: five 100 ns
trajectories exist, and their first 10 ns is a 10 ns sweep. This truncates each
one and asks whether the early window ranks the molecules the way the full run
does.

WHAT WOULD KILL THE IDEA. If a molecule's attack geometry in 0-10 ns is
uncorrelated with its 0-100 ns geometry, the sweep filters on noise and would
discard good molecules 10x more cheaply than before. If it tracks, the filter is
free: 10 ns is ~0.4 GPU-h against ~4 GPU-h, so ten molecules screen for the price
of one elevated.

THE EARLY WINDOW IS NOT NEUTRAL, and that is the interesting risk rather than a
technicality. Every run STARTS from the docked pose, which was chosen for its
anchoring geometry -- so the first nanoseconds are biased toward looking
competent, and a sweep could pass everything. The equilibration-only fraction is
reported separately so that bias is visible rather than assumed away.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402
from shared import md_movie as mov                  # noqa: E402
from shared import nac_criterion as nac             # noqa: E402
from shared import run_paths as rp            # noqa: E402

log = logging.getLogger("attack-sweep")
MD = Path("/data/lab_vm/modifiable/inhibition/md_residence_3ikd")
OUT = sout.Topic("blacksmith", rp.reports_topic())   # scoped to this run (#74)

#: 100 ps resolution over 100 ns. Enough that a 10 ns window still holds 100
#: points, which is what makes the truncation comparison meaningful at all.
FRAMES = 1000
SWEEP_NS = 10.0


def _mp():
    spec = importlib.util.spec_from_file_location(
        "mdprio_report", REPO / "scripts" / "mdprio_report.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["mdprio_report"] = m
    spec.loader.exec_module(m)
    return m


def competent(angle: np.ndarray, kind: str) -> np.ndarray:
    """Per-frame angular competence, by the mechanism's own criterion."""
    if "anti" in (kind or ""):
        return angle >= nac.SN2_ANGLE_MIN
    return angle <= nac.PERPENDICULAR_MAX_OFF_NORMAL


def stats(dist: np.ndarray, angle: np.ndarray, kind: str) -> dict:
    inw = (dist >= nac.NAC_DIST_MIN) & (dist <= nac.NAC_DIST_MAX)
    ok = competent(angle, kind)
    return {"frac_in_window": float(inw.mean()),
            "frac_attack_ready": float((inw & ok).mean()),
            "median_dist_a": float(np.median(dist)),
            "median_angle_deg": float(np.median(angle))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", nargs="+", required=True)
    ap.add_argument("--sweep-ns", type=float, default=SWEEP_NS)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    mp = _mp()

    rows = []
    for cand in args.candidates:
        rep = MD / cand / "md" / "rep1"
        if not (rep / "prod.xtc").is_file():
            log.warning("%s: no trajectory", cand)
            continue
        dense = rep / "sweep_dense.pdb"
        if not dense.is_file():
            mov.build_movie_pdb(rep, dense, n_frames=FRAMES)
        s = mp.nac_series(cand, rep, dense)
        if s is None:
            log.warning("%s: no attack-geometry series", cand)
            continue
        t, d, a, kind = s["t"], s["dist"], s["angle"], s["kind"]
        early = t <= args.sweep_ns
        rec = {"ident": cand, "n_frames": len(t), "mechanism": s["mechanism"]}
        for tag, sel in (("full", np.ones_like(early)), ("sweep", early)):
            for k, v in stats(d[sel], a[sel], kind).items():
                rec[f"{tag}_{k}"] = v
        rows.append(rec)
        log.info("%s: sweep ready %.3f -> full %.3f", cand,
                 rec["sweep_frac_attack_ready"], rec["full_frac_attack_ready"])

    t = pd.DataFrame(rows)
    if t.empty:
        print("no trajectories analysed")
        return
    dest = OUT.write("attack_sweep", ".csv")
    t.to_csv(dest, index=False)

    print(f"\n#32 — would a {args.sweep_ns:.0f} ns sweep have seen what 100 ns saw?"
          f"   (n={len(t)})\n")
    cols = ["ident", "sweep_frac_in_window", "full_frac_in_window",
            "sweep_frac_attack_ready", "full_frac_attack_ready"]
    print(t[cols].round(4).to_string(index=False))

    if len(t) >= 3:
        from scipy.stats import spearmanr
        for m in ("frac_in_window", "frac_attack_ready"):
            r, p = spearmanr(t[f"sweep_{m}"], t[f"full_{m}"])
            print(f"\n  rho(sweep, full) for {m}: {r:+.3f}  (p={p:.3f}, n={len(t)})")
    print(f"\n  A {args.sweep_ns:.0f} ns sweep costs ~0.4 GPU-h against ~4 for a "
          f"full run, so it pays for itself if it ORDERS molecules correctly —")
    print("  it does not need to reproduce the values.")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
