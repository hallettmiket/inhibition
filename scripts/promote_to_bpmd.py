#!/usr/bin/env python3
"""
Purpose: promote a 100 ns holder to BPMD, automatically and idempotently.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-16
Input: --gpu N [--watch] [--bar 0.35] [--replicates 3] [--production-ps 10000]
Output: launches scripts/bpmd_run.py per qualifying molecule; prints a manifest

@tt8804: "now that 100 ns MD is much tighter and we expect even less candidates,
lets automatically run BPMD on any mols that keep below 0.35 nm rmsd max for the
whole 100 ns".

THE FOURTH STAGE OF THE CASCADE, and the first one that is not a stability test.
Docking ranks the best case; the 8 ns sweep asks whether the pose is stable at
all; 100 ns asks whether it is stable for a long time. BPMD asks something
different -- how much bias does it take to push the ligand OUT -- so it separates
poses that merely survived an unperturbed trajectory from poses that are actually
held. A molecule reaching this stage has already earned ~4.5 GPU-hours, so the
population is small by construction and the per-molecule cost is affordable.

THE BAR IS THE SAME NUMBER AS THE SWEEP'S, DELIBERATELY: max ligand RMSD < 0.35
nm, now over the full 100 ns rather than 8. Reading it from
`md.sweep_survivor_rmsd_nm` keeps one definition of "did not move"; a second
constant here would drift the moment either is retuned.

IDEMPOTENT, BECAUSE IT IS MEANT TO BE LEFT RUNNING. Every launch checks the bpmd
topic for an existing row on that ident. `--watch` polls, so a 100 ns run that
lands at 3am is promoted without anybody waiting for it.

MAX RMSD IS A NOISY STATISTIC AND THAT IS RECORDED HERE. It is an extreme value
over 10,001 frames, so it is set by the single largest excursion. Three
replicates of t4_716800c125a7 gave 0.275, 0.317 and 0.625 nm while their MEANS
were 0.148 and 0.237 -- the max varies 2.3x where the mean varies 1.6x. The
promotion therefore also records mean RMSD and engaged fraction, so a later
decision can be re-derived on a steadier statistic without re-running anything.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("promote-bpmd")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
PY = Path.home() / ".micromamba/envs/dwi_reactive/bin/python"

RMSD_MAX = "explicit_ligand_rmsd_nm_max"
RMSD_MEAN = "explicit_ligand_rmsd_nm_mean"
ENGAGED = "explicit_frac_frames_engaged"


def holders(bar: float, min_ps: float = 90_000.0) -> pd.DataFrame:
    """Molecules whose 100 ns run stayed under `bar`, newest reading per ident.

    `min_ps` guards against promoting on a partial run: a 50 ns trajectory has
    had half the chance to wander, so its max RMSD is not comparable with a
    100 ns one and would clear the bar too easily.
    """
    rows = []
    for f in glob.glob(str(rp.residence_dir() / "*.csv")):
        try:
            d = pd.read_csv(f)
        except Exception:                                  # noqa: BLE001
            continue
        if "production_ps" not in d.columns or RMSD_MAX not in d.columns:
            continue
        d = d[(d.production_ps >= min_ps) & d[RMSD_MAX].notna()]
        if len(d):
            d = d.copy(); d["_t"] = os.path.getmtime(f); rows.append(d)
    if not rows:
        return pd.DataFrame()
    d = pd.concat(rows, ignore_index=True).sort_values("_t")
    d = d.drop_duplicates("ident", keep="last")
    return d[d[RMSD_MAX] < bar]


def already_run() -> set:
    out = set()
    for f in glob.glob(str(rp.bpmd_dir() / "*.csv")):
        try:
            d = pd.read_csv(f, usecols=["ident"])
        except Exception:                                  # noqa: BLE001
            continue
        out |= set(d.ident.astype(str))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--gpu", type=int, required=True)
    ap.add_argument("--bar", type=float, default=None,
                    help="max ligand RMSD over the 100 ns; defaults to "
                         "md.sweep_survivor_rmsd_nm so there is one definition")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--production-ps", type=float, default=10_000)
    ap.add_argument("--watch", action="store_true",
                    help="keep polling; promote each 100 ns run as it lands")
    ap.add_argument("--poll-s", type=int, default=900)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.bar is None:
        from shared import target_config as tc
        args.bar = float(tc.get("md.sweep_survivor_rmsd_nm", default=0.35))
    log.info("promoting 100 ns runs with max ligand RMSD < %.2f nm to BPMD "
             "(%d replicates x %.0f ps) on GPU %d",
             args.bar, args.replicates, args.production_ps, args.gpu)

    seen: set = set()
    while True:
        done = already_run()
        h = holders(args.bar)
        todo = [r for r in h.to_dict("records")
                if str(r["ident"]) not in done and str(r["ident"]) not in seen]
        if h.empty:
            log.info("no 100 ns run has cleared %.2f nm yet", args.bar)
        for r in todo:
            ident = str(r["ident"])
            log.info("PROMOTE %s  max %.3f nm  mean %s  engaged %s",
                     ident, r[RMSD_MAX],
                     f"{r.get(RMSD_MEAN):.3f}" if pd.notna(r.get(RMSD_MEAN)) else "—",
                     f"{r.get(ENGAGED):.3f}" if pd.notna(r.get(ENGAGED)) else "—")
            if args.dry_run:
                seen.add(ident); continue
            cmd = [str(PY), str(REPO / "scripts/bpmd_run.py"),
                   "--pose", ident, "--gpu", str(args.gpu),
                   "--replicates", str(args.replicates),
                   "--production-ps", str(args.production_ps)]
            rc = subprocess.run(cmd, capture_output=True, text=True)
            if rc.returncode != 0:
                # Recorded, not silently retried forever: a molecule that cannot
                # be set up will fail identically on every poll.
                log.warning("%s: bpmd failed rc=%d %s", ident, rc.returncode,
                            (rc.stderr or rc.stdout)[-200:])
            else:
                log.info("%s: bpmd done", ident)
            seen.add(ident)
        if not args.watch:
            break
        time.sleep(args.poll_s)


if __name__ == "__main__":
    main()
