"""
Purpose: tier-1 warhead drift for an ARBITRARY set of molecules, with the pose named rather than defaulted.
Author: @twu383 (with Claude Code)
Date: 2026-08-31
Input: a run topic with a ranked screen + a list of candidate idents
Output: 00_outputs/blacksmith/panel_tier1_<topic>/panel_t1_<N>.csv

WHY NOT `elevation_run.py --tier 1`. That runs `cohort()` -- the pre-registered
stratified cohort from `elevation_cohort/` -- and `--only` filters WITHIN it. A
molecule that is not in that CSV cannot be measured by it at all, which is the
right behaviour for a pre-registered experiment and useless for asking "does a
VALIDATED warhead on this exact scaffold behave differently from an unvalidated
one". This runs the identical measurement on a named set.

WHAT IS MEASURED, and it is the only readout on this project that has passed its
own validation: the warhead-to-SG distance BEFORE and AFTER 300 ps of
unrestrained equilibration (`gromacs_explicit.NVT_PS` 100 + `NPT_PS` 200, no
position restraints), per replicate. D0071 measured that it separates
crystallographic Cys113 binders from generated candidates at p = 0.007, Cliff
delta -0.781; D0072 rested a NO GO on it.

THE POSE IS NAMED, NOT DEFAULTED. `elevation_run.tier1_pose` calls
`prepare_pose` without `pose_rank`, so it takes rank 1 -- the same defect D0105
records in `bpmd_run`. Here the rank is chosen per molecule as the TOP mode by
`engagement` (`ranking.score_by_tier` for T_4, D0098) and written into every
row, so the comparison is "each molecule's best mode" rather than "whatever
sorted first".

COMPARABILITY IS THE WHOLE POINT, so the things that would break it are held
fixed and recorded: one screen (one topic, one seed, one receptor, one splitter)
produced every pose, and every molecule shares an R-group by construction. What
differs between rows is the warhead, and the run records `docked_distance_nm` so
a reader can check that the starting distances really were comparable rather
than assuming it.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import bpmd_run as br                                  # noqa: E402
from shared import bpmd                                # noqa: E402
from shared import gromacs_explicit as gx              # noqa: E402
from shared import outputs as sout                     # noqa: E402
from shared import run_paths as rp                     # noqa: E402


def _elevation_run():
    """`elevation_run.distance_nm`, imported by path.

    Its module name starts with a digit-free word but it is a SCRIPT, and the
    measurement must be the same function tier 1 used -- a second
    implementation of "distance across a periodic box" is how two tables that
    look comparable stop being comparable. It handles the PBC wrap that a naive
    norm silently gets wrong when the ligand crosses a box face.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "elevation_run_for_panel", REPO / "scripts" / "elevation_run.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


log = logging.getLogger("panel-t1")


def top_mode_by_engagement(topic: str) -> dict[str, int]:
    """{ident: pose_rank} for each molecule's highest-engagement mode.

    `pose_rank` is `mode + 1` in the representative SDF, but that is an
    observation about how the export happens to be ordered, not a contract --
    so the rank is READ from the pose file's own `pose_rank` property, keyed on
    the mode. Deriving it arithmetically is how the wrong pose gets simulated
    the day the export changes (D0105).
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_engagement_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no engagement ranking for topic {topic!r}")
    d = pd.read_csv(fs[-1])
    log.info("ranking %s: %d modes over %d molecules",
             Path(fs[-1]).name, len(d), d.parent_ident.nunique())

    out = {}
    for ident, g in d.groupby("parent_ident"):
        best = g.loc[g.engagement.idxmax()]
        mode = int(best["mode"])
        sdf = rp.BLACKSMITH / f"{topic}_poses" / f"{ident}.sdf"
        rank = None
        for m in Chem.SDMolSupplier(str(sdf), removeHs=False):
            if m is None or not m.HasProp("mode"):
                continue
            if int(m.GetProp("mode")) == mode:
                rank = int(m.GetProp("pose_rank"))
                break
        if rank is None:
            raise SystemExit(f"{ident}: mode {mode} has no representative in "
                             f"{sdf.name}; refusing to guess a pose_rank")
        out[str(ident)] = rank
        log.info("  %-18s mode %-4d -> pose_rank %-4d engagement %.4f "
                 "(mode_size %d)", ident, mode, rank, best.engagement,
                 int(best.mode_size))
    return out


def measure(cand, pose_rank: int, sdf: Path, *, replicates: int, gpu: int,
            threads: int, on_row) -> None:
    """Equilibrate `replicates` replicas and record the warhead's drift."""
    base = {"ident": cand.ident, "warhead_class": cand.warhead_class,
            "mechanism": cand.mechanism, "pose_rank": pose_rank}
    try:
        prep = br.prepare_pose(cand, nrun=0, gpu=str(gpu), allow_redock=False,
                               pose_rank=pose_rank, sdf=sdf)
    except Exception as exc:                                  # noqa: BLE001
        row = {**base, "replicate": 0,
               "status": f"setup failed: {type(exc).__name__}: {str(exc)[:200]}"}
        log.warning("  SETUP FAILED: %s", row["status"])
        on_row(row)
        return

    meta = {k: v for k, v in prep.items()
            if isinstance(v, (int, float, str, bool, type(None)))}
    # The docked distance IN THE MD FRAME -- not the docking-time value, which
    # was measured against a flexible Cys113 sidechain and is a distance to a
    # different sulfur position. `elevation_run.tier1_pose` makes the same
    # distinction and for the same reason.
    d_dock = float(prep["start_distance_nm"])
    want = (prep["warhead_atom_name"], prep["sg_atom_name"])

    for k in range(1, replicates + 1):
        row = {**base, **meta, "replicate": k, "status": "ok",
               "docked_distance_nm": round(d_dock, 4)}
        try:
            t0 = time.time()
            res = gx.run_pipeline(prep["wd"], prep["md_wd"], gpu_id=gpu,
                                  threads=threads, replicate=k,
                                  candidate_id=cand.ident, stop_after="npt",
                                  reuse_equilibration=True)
            eq = Path(res["equilibration_dir"])
            d = ER.distance_nm(eq / "npt.gro", prep["warhead_serial0"],
                               prep["sg_serial0"], want)
            delta = d["distance_nm"] - d_dock
            row.update({
                "npt_distance_nm": round(d["distance_nm"], 4),
                "delta_nm": round(delta, 4),
                "abs_delta_nm": round(abs(delta), 4),
                "pbc_wrapped": d["pbc_wrapped"],
                "start_in_window": bool(bpmd.NAC_MIN_NM <= d_dock <= bpmd.NAC_MAX_NM),
                "npt_in_window": bool(bpmd.NAC_MIN_NM <= d["distance_nm"]
                                      <= bpmd.NAC_MAX_NM),
                "equilibration_ps": gx.NVT_PS + gx.NPT_PS,
                "velocity_seed": res["velocity_seed"],
                "seconds": round(time.time() - t0, 1)})
            log.info("    rep%d: %.3f -> %.3f nm (delta %+.3f), %s",
                     k, d_dock, d["distance_nm"], delta,
                     "still in the NAC window" if row["npt_in_window"] else "outside")
        except Exception as exc:                              # noqa: BLE001
            row["status"] = f"failed: {type(exc).__name__}: {str(exc)[:200]}"
            log.warning("    rep%d FAILED: %s", k, row["status"])
        on_row(row)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", required=True)
    ap.add_argument("--idents", required=True, help="file of candidate idents")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--gpu", type=int, required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [s{args.shard}] %(message)s")
    if args.gpu in (0, 4, 7):
        raise SystemExit(f"GPU {args.gpu} is reserved for other users' jobs")
    os.nice(19)

    want = [l.strip() for l in Path(args.idents).read_text().splitlines() if l.strip()]
    ranks = top_mode_by_engagement(args.topic)
    by_id = br.candidate_index()

    todo = [i for i in want if i in ranks]
    missing = [i for i in want if i not in ranks]
    if missing:
        log.warning("%d ident(s) absent from the ranking and SKIPPED: %s",
                    len(missing), missing)
    todo = [t for k, t in enumerate(todo) if k % args.n_shards == args.shard]
    log.info("%d molecules on this shard, %d replicates each",
             len(todo), args.replicates)

    global ER
    ER = _elevation_run()
    OUT = sout.Topic("blacksmith", f"panel_tier1_{args.topic}")
    chunk = br._ChunkWriter(OUT, f"panel_t1_s{args.shard}", 50)
    pose_dir = rp.BLACKSMITH / f"{args.topic}_poses"
    for i, ident in enumerate(todo, 1):
        cand = by_id.get(ident)
        if cand is None:
            chunk.add({"ident": ident, "replicate": 0,
                       "status": "failed: not in the candidate set"})
            continue
        log.info("[%d/%d] %s (%s) pose_rank %d",
                 i, len(todo), ident, cand.warhead_class, ranks[ident])
        measure(cand, ranks[ident], pose_dir / f"{ident}.sdf",
                replicates=args.replicates, gpu=args.gpu,
                threads=args.threads, on_row=chunk.add)
    chunk.flush()
    print(f"\n  -> {OUT.dir}")


if __name__ == "__main__":
    main()
