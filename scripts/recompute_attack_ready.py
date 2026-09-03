#!/usr/bin/env python3
"""
Purpose: re-score finished sweeps under the current attack-ready definition, without re-running MD.
Author: @twu383 (with Claude Code)
Date: 2026-09-02
Input: the run's sweep rows + each mode's persisted `sweep_dense.pdb`
Output: a NEW versioned attack_sweep_<N>.csv carrying the recomputed readings

@twu383, 2026-09-02: *"we need to clearly establish that attack ready means
within angle and under 3 A at any time"*.

`frac_attack_ready` used `nac_criterion.NAC_DIST_MAX` (4.2 A), the near-attack
WINDOW, so a mode whose warhead sat at a trajectory median of 3.6 A scored 93%
attack-ready. The worklist selects modes at < 3.0 A, so the selection and the
readout disagreed about what "close" meant, and the readout was the looser of
the two.

WHY THIS DOES NOT NEED THE GPU. The per-frame geometry comes from
`sweep_dense.pdb`, which `attack_sweep` already writes and keeps for every
completed sweep. Re-scoring is arithmetic over frames that are already on disk:
94 trajectories re-read in seconds against ~6 GPU-hours to re-simulate them.

THE OLD ROWS ARE NOT EDITED. The outputs root is append-only, and a rewritten
row would make two different definitions indistinguishable in the same file.
This writes a new versioned table; every row in it carries
`attack_ready_max_a`, so which definition produced a number is a property of
the row and not of when you happened to read it.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                        # noqa: E402
from shared import run_paths as rp                        # noqa: E402

log = logging.getLogger("recompute-ar")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    import attack_sweep as asw
    # THE SAME LOADER `attack_sweep` USES, so the geometry comes from one
    # implementation. A second copy of "the attack angle" is how a plot and a
    # ranking come to disagree while both look right (mdprio_report.nac_series).
    mp = asw._mp()

    topic = args.topic or rp.topic()
    hi = asw.attack_ready_max_a()
    log.info("attack ready := %.1f <= d < %.1f A AND angle within the "
             "mechanism's own criterion", asw.nac.NAC_DIST_MIN, hi)

    fs = sorted(glob.glob(str(rp.BLACKSMITH / rp.sweep_topic(topic) / "*.csv")))
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d.drop_duplicates("ident", keep="last")
    ok = d[d.status.astype(str) == "ok"].copy()
    log.info("%d completed sweeps to re-score", len(ok))

    import sweep_assets as sa
    wl = sorted((rp.BLACKSMITH / f"sweep_worklist_{topic}").glob("worklist_*.csv"),
                key=lambda f: int(f.stem.split("_")[-1]))
    prank = {}
    if wl:
        w = pd.read_csv(wl[-1])
        key = "task_id" if "task_id" in w.columns else "ident"
        prank = dict(zip(w[key].astype(str), w.pose_rank.astype(int)))

    rows, missed = [], 0
    for r in ok.itertuples():
        ident = str(r.ident)
        parent = ident.rsplit("_m", 1)[0]
        rep = sa.rep_dir(parent, prank.get(ident))
        dense = (rep / "sweep_dense.pdb") if rep else None
        if dense is None or not dense.is_file():
            missed += 1
            continue
        try:
            s = mp.nac_series(parent, rep, dense)
            if s is None:
                missed += 1
                continue
            dist = np.asarray(s["dist"], dtype=float)
            angle = np.asarray(s["angle"], dtype=float)
            kind = str(getattr(r, "angle_kind", "off-normal"))
            frame_ps = float(getattr(r, "frame_ps", 0) or 0) or (
                float(getattr(r, "sweep_ps", 1200.0)) / max(1, len(dist)))
            st = asw.geometry_stats(dist, angle, kind, frame_ps)
            # The RMSD half of the 100 ns gate, from the same `rmsd.xvg` the
            # run already wrote. Re-scored here too, so an old row and a new one
            # carry the verdict under the SAME definition rather than one having
            # it and the other not.
            st.update(asw.rmsd_stats(rep))
        except Exception as exc:                          # noqa: BLE001
            log.warning("%s: %s", ident, exc)
            missed += 1
            continue
        base = {c: getattr(r, c, None) for c in ok.columns if c != "Index"}
        base.update(st)
        base.update(asw.elevation_verdict(base))
        base["rescored_from"] = "sweep_dense.pdb"
        rows.append(base)

    if missed:
        log.warning("%d sweep(s) had no readable trajectory and keep their "
                    "original reading", missed)
    if not rows:
        raise SystemExit("nothing re-scored")

    new = pd.DataFrame(rows)
    old_ar = ok.set_index("ident").frac_attack_ready
    new["frac_attack_ready_old"] = new.ident.map(old_ar)
    moved = (new.frac_attack_ready - new.frac_attack_ready_old).abs() > 1e-9
    log.info("re-scored %d; %d changed", len(new), int(moved.sum()))
    ang = " + angle" if asw.attack_ready_use_angle() else " (distance only)"
    th = asw.elevation_thresholds()
    print(f"\n  attack ready := {asw.nac.NAC_DIST_MIN:.1f} <= d < {hi:.1f} A{ang}")
    print(f"  elevate      := engaged >= {th['occupancy_min']*100:.0f}%  AND  "
          f"(rmsd max < {th['rmsd_max_a']:.1f} A or mean < {th['rmsd_mean_a']:.1f} A)")
    print(f"  {'':22s}{'was':>8}{'now':>8}")
    for _, x in new.nlargest(8, "frac_attack_ready_old").iterrows():
        print(f"  {x.ident:<22}{x.frac_attack_ready_old*100:7.1f}%"
              f"{x.frac_attack_ready*100:7.1f}%")
    print(f"\n  modes >1% engaged: "
          f"{int((new.frac_attack_ready > 0.01).sum())} of {len(new)} "
          f"(was {int((new.frac_attack_ready_old > 0.01).sum())})")
    if "elevate" in new.columns:
        held = int(new.pose_held.fillna(False).astype(bool).sum())
        eng = int(new.warhead_engaged.fillna(False).astype(bool).sum())
        print(f"  pose held (rmsd):   {held} of {len(new)}")
        print(f"  warhead engaged:    {eng} of {len(new)}")
        print(f"  ELEVATE to 100 ns:  {int(new.elevate.fillna(False).sum())} "
              f"of {len(new)}")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return
    OUT = sout.Topic("blacksmith", rp.sweep_topic(topic))
    dest = OUT.write("attack_sweep", ".csv")
    new.drop(columns=["frac_attack_ready_old"]).to_csv(dest, index=False)
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
