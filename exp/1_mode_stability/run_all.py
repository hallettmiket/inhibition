#!/usr/bin/env python3
"""
Purpose: how often does a fresh 500-pose docking cloud elect the SAME binding mode?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18
Input: --candidate (default t4_716800c125a7), --replicates 5, --nrun 500, --gpu
Output: append_only/00_outputs/blacksmith/mode_stability_<candidate>/ + a summary table

@tt8804: "mol generation is 500 poses split into modes. run an experiment do 500
poses times 5 runs for our last top hit and see how many times we actually
identify the right mode by splitting 500".

WHY. AutoDock-GPU's Lamarckian GA is seeded from the clock -- `nac_screen.dock`
passes `--nrun` and no `--seed` -- so each screen draws a different cloud. The
mode split and the ranking are deterministic GIVEN the poses, so all downstream
instability enters here. Across 3.0.0 and 3.1.0 the same 504 molecules ranked at
rho = +0.43 and only 22.6% kept the same winning sub-mode; t4_716800c125a7 went
from class_rank 26 (a mode that then held 100 ns at 0.317 nm) to 89 (a mode that
sweeps at 0.854 nm). This measures that instability directly, on one molecule,
without the confound of two different code versions.

WHAT "THE RIGHT MODE" MEANS HERE. Not an oracle -- there is no ground truth for
which mode is correct. It is the mode 3.0.0 elected and then VALIDATED with a
100 ns run that held: centroid and warhead direction recorded in REFERENCE below.
A replicate "recovers" it if some mode of that replicate sits within MATCH_A of
that centroid and within MATCH_DEG of that direction; it "elects" it if that mode
is also the replicate's top-ranked one. Recovery and election are reported
separately, because a screen that finds the mode but ranks it 4th is a different
failure from one that never samples it.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                       # noqa: E402
from shared import target_config as tc                   # noqa: E402

log = logging.getLogger("mode-stability")

#: The mode 3.0.0 elected AND validated at 100 ns (max ligand RMSD 0.317 nm).
#: Taken from rank_v2_T4_nac_v4_enrichment_conditional_4.csv, mode 0.
REFERENCE = {
    "candidate": "t4_716800c125a7",
    "centroid": np.array([14.194164419827397, 3.632780812374533, 1.0417123323638144]),
    "direction": np.array([0.7756940765993148, -0.3236334203074604, -0.5418118758285304]),
    "mode_size": 73, "viable_fraction": 0.3835616438356164,
    "enrichment": 4.698542083349057, "class_rank": 26.0,
    "validated_100ns_rmsd_max_nm": 0.317,
}

#: A mode matches the reference if its warhead cluster sits this close in space
#: and points this nearly the same way. 2.0 A is the stage-2 cut diameter -- two
#: clusters further apart than that are, by the screen's own definition,
#: different modes. 45 deg is half the angular width the stage-1 direction
#: clustering tolerates.
MATCH_A = 2.0
MATCH_DEG = 45.0


def screen_once(cand: str, topic: str, nrun: int, gpu: str, sub_split: int) -> Path:
    """One independent screen of one molecule into its own topic.

    Its own topic per replicate, because rank_v2 concatenates every agg file in a
    topic: two replicates sharing one would be read as one molecule with twice
    the modes.
    """
    only = Path(f"/tmp/mode_stability_{topic}.txt")
    only.write_text(cand + "\n")
    cmd = [sys.executable, str(REPO / "scripts/nac_screen_v2.py"),
           "--only", str(only), "--topic", topic, "--nrun", str(nrun),
           "--gpu", gpu, "--sub-split", str(sub_split), "--all-poses"]
    log.info("replicate -> topic %s (nrun=%d, gpu=%s)", topic, nrun, gpu)
    r = subprocess.run(["nice", "-n", "19"] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        log.error("screen failed for %s:\n%s", topic, r.stderr[-2000:])
        raise RuntimeError(f"screen failed for {topic}")
    return rp.BLACKSMITH / topic


def modes_of(topic_dir: Path) -> pd.DataFrame:
    """Every mode this replicate produced, with its geometry and its score."""
    fs = sorted(topic_dir.glob("agg_s*_*.csv"))
    if not fs:
        raise FileNotFoundError(f"no agg_s*.csv under {topic_dir}")
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d


def match(d: pd.DataFrame) -> pd.DataFrame:
    """Distance and angle from each mode to the reference."""
    c = d[["centroid_x", "centroid_y", "centroid_z"]].to_numpy(float)
    v = d[["dir_x", "dir_y", "dir_z"]].to_numpy(float)
    d = d.copy()
    d["dist_a"] = np.linalg.norm(c - REFERENCE["centroid"], axis=1)
    n = np.linalg.norm(v, axis=1)
    n[n == 0] = np.nan
    cos = np.clip((v @ REFERENCE["direction"]) / n, -1.0, 1.0)
    d["angle_deg"] = np.degrees(np.arccos(cos))
    d["matches_ref"] = (d.dist_a <= MATCH_A) & (d.angle_deg <= MATCH_DEG)
    return d


def rank_within(d: pd.DataFrame) -> pd.DataFrame:
    """Order this replicate's modes the way the screen would.

    `viable_fraction / isotropic_null` is `enrichment`, which is what the sweep
    rule ranks on. Modes with nothing in range are unrankable and are left NaN
    rather than given a 0 -- a mode that never reached the window is not a mode
    that reached it badly.
    """
    d = d.copy()
    key = "enrichment" if "enrichment" in d.columns else "viable_fraction"
    d.loc[d.get("n_in_range", 1) == 0, key] = np.nan
    d["rank_in_replicate"] = d[key].rank(ascending=False, method="min")
    return d.sort_values(key, ascending=False, na_position="last")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default=REFERENCE["candidate"])
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--nrun", type=int, default=int(tc.get("docking.n_runs", default=500)))
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--sub-split", type=int, default=5)
    ap.add_argument("--skip-screen", action="store_true",
                    help="reuse replicate topics already on disk")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    rows, per_replicate = [], []
    for i in range(1, args.replicates + 1):
        topic = f"mode_stability_{args.candidate}_r{i}"
        tdir = rp.BLACKSMITH / topic
        if not args.skip_screen or not tdir.is_dir():
            screen_once(args.candidate, topic, args.nrun, args.gpu, args.sub_split)
        d = rank_within(match(modes_of(tdir)))
        d["replicate"] = i
        per_replicate.append(d)
        hit = d[d.matches_ref]
        top = d.iloc[0] if len(d) else None
        rows.append({
            "replicate": i,
            "n_modes": len(d),
            "cloud": int(d.n_poses.iloc[0]) if "n_poses" in d.columns and len(d) else np.nan,
            "assigned": int(d.mode_size.sum()) if "mode_size" in d.columns else np.nan,
            "recovered": bool(len(hit)),
            "ref_rank": float(hit.rank_in_replicate.min()) if len(hit) else np.nan,
            "ref_viable_fraction": float(hit.viable_fraction.max()) if len(hit) else np.nan,
            "elected": bool(len(hit) and hit.rank_in_replicate.min() == 1),
            "top_mode_dist_a": float(top.dist_a) if top is not None else np.nan,
            "top_mode_vf": float(top.viable_fraction) if top is not None else np.nan,
        })

    out = rp.BLACKSMITH / f"mode_stability_{args.candidate}"
    out.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "summary_1.csv", index=False)
    pd.concat(per_replicate, ignore_index=True).to_csv(out / "all_modes_1.csv", index=False)

    n = len(summary)
    rec, ele = int(summary.recovered.sum()), int(summary.elected.sum())
    print(f"\n  {args.candidate}: {n} independent {args.nrun}-pose screens")
    print(f"  reference = the mode 3.0.0 elected and validated at 100 ns "
          f"(max RMSD {REFERENCE['validated_100ns_rmsd_max_nm']} nm)\n")
    print(summary.to_string(index=False))
    print(f"\n  RECOVERED (mode sampled at all) : {rec}/{n}")
    print(f"  ELECTED   (and ranked first)    : {ele}/{n}")
    if rec:
        print(f"  when recovered, its rank among that replicate's modes: "
              f"{sorted(summary.ref_rank.dropna().astype(int).tolist())}")
    (out / "summary_1.json").write_text(json.dumps(
        {"candidate": args.candidate, "replicates": n, "nrun": args.nrun,
         "recovered": rec, "elected": ele,
         "match_a": MATCH_A, "match_deg": MATCH_DEG}, indent=2))
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
