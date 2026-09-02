#!/usr/bin/env python3
"""
Purpose: dock a reference compound through the screen's own path, keeping every pose.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-11
Input: --ident (a pose_sidecars entry) or --smiles, plus its warhead class
Output: 00_outputs/blacksmith/<topic>/{agg,poses}_*.csv + nac_v3_poses-style SDF

DOES THE SCREEN EVEN GENERATE THE CRYSTAL POSE? (@tt8804, #56/#39). Sulfopin is
the one molecule whose true binding position is known, and it is **not in the
library** -- 0 rows in the ranking -- so nothing in this project has ever docked
it through the same path a candidate takes. `nac_poses/xtal:6VAJ:QT7.sdf` is the
crystal pose put through the exporter, not a docking of it.

This runs one named compound through `nac_screen_v2.one`, the identical function
every candidate goes through: same receptor, same reactive potential, same
`--nrun`, same mode clustering. Anything else would answer a different question
than "would our screen have found it".

WRITTEN WITH --all-poses ALWAYS. The whole point is to compare the pose CLOUD
against a known answer, and this repo has already lost one cloud to a temporary
directory (#44, now a rule in CLAUDE.md).

ITS OWN TOPIC, NEVER `nac_v3`. `rank_v2` concatenates every `agg_s*_*.csv` in a
topic, so writing a reference compound into the production topic would put it in
the ranking as if it were a library candidate.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import nac_screen as ns                                    # noqa: E402
import nac_screen_v2 as nsv                                # noqa: E402
from shared import outputs as sout                         # noqa: E402
from shared import reference_set as rs                     # noqa: E402
from shared import pose_subsplit as psub                   # noqa: E402

log = logging.getLogger("dock-ref")
B = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
SIDECARS = B / "pose_sidecars"


def warhead_row(class_id: str) -> pd.Series:
    """The named class, from the HIGHEST-NUMBERED warhead library.

    RESOLVED NUMERICALLY, via the project's own resolver. This was
    `sorted(glob(...))[-1]`, which sorts LEXICALLY -- so once the library
    reached `_10`, `"..._9.csv"` sorted last and the newest file became
    invisible. `_10` adds exactly one class, `cinnamamide`, so every caller
    here silently saw a 10-class library instead of an 11-class one.

    This is the stale-pin defect (catalogue #6, #7) wearing a glob's clothes:
    the code looks like it resolves dynamically, and it does -- to the wrong
    file, for as long as the version count has two digits. `latest_reference`
    compares the integer, so it cannot do this.
    """
    return rs.load_warhead_row(class_id)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--ident", required=True,
                    help="name for the run; if a pose_sidecars/<ident>.json "
                         "exists its canonical_smiles is used")
    ap.add_argument("--smiles", default=None, help="overrides the sidecar")
    ap.add_argument("--warhead-class", required=True)
    ap.add_argument("--nrun", type=int, default=500,
                    help="MUST match the library's, or the pose cloud is not "
                         "comparable to a candidate's")
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--topic", default="ref_modes")
    # THE SEED IS THE LIBRARY'S UNLESS SAID OTHERWISE. #77 fixed the screen's
    # seed so a run is reproducible; this path passed nothing, so `one()` fell
    # to `seed=None` and drew from the clock -- a control that cannot be re-run
    # to the same cloud, which is the one property a control needs most.
    ap.add_argument("--seed", type=int, default=nsv._cfg("docking.seed", None),
                    help="docking seed; defaults to docking.seed from the config")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    smi = args.smiles
    if smi is None:
        f = SIDECARS / f"{args.ident}.json"
        if not f.is_file():
            raise SystemExit(f"no --smiles and no sidecar at {f}")
        smi = json.loads(f.read_text())["canonical_smiles"]
    w = warhead_row(args.warhead_class)
    log.info("%s  %s", args.ident, smi)
    log.info("  class %s / %s, reactive %s",
             w.class_id, w.mechanism, w.reactive_atom_smarts)

    cand = ns.Candidate(ident=args.ident, smiles=smi,
                        warhead_class=str(w.class_id),
                        mechanism=str(w.mechanism),
                        reactive_smarts=str(w.reactive_atom_smarts),
                        label="positive")

    # ALL FOUR OUTPUT PATHS MOVE TOGETHER, via the screen's own function.
    # This used to be `nsv.OUT = sout.Topic(...)` and nothing else, so the tables
    # landed in `--topic` while `one()` -- which reads POSE_DIR and ALL_POSE_DIR
    # off the module at write time -- put this compound's representative AND its
    # whole 500-pose cloud into the PRODUCTION run's directories. A reference
    # compound is not in the library, so anything globbing `<run.topic>_allposes`
    # would have picked it up as a candidate. Same defect as 2.2.0's, one caller
    # out from where the test for it looks.
    OUT, _pose_dir, _all_pose_dir = nsv.use_topic(args.topic)
    log.info("topic %s -> poses %s, clouds %s",
             args.topic, _pose_dir.name, _all_pose_dir.name)
    # The SAME receptors the production screen builds, resolved the same way --
    # a reference docked into a differently-prepared receptor answers a different
    # question than "would our screen have found it" (D0059).
    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    plain_rec = sout.latest_path("blacksmith", "receptor_3ikd",
                                 "3IKD_prepared", ".pdbqt")
    # SPLITTING IS PASSED EXPLICITLY, NOT LEFT TO THE SIGNATURE DEFAULT. The
    # point of this script is a number comparable to the library's, and
    # `consensus` = mode_size/n_poses -- which `conditional_eb` is computed from
    # -- changes with the sub-split. `one()`'s default is psub.DEFAULT_MAX_SUB;
    # the screen's default is the config. Those agree today only because
    # `contact_linkage` forces the sub-split to 1 regardless. Reading the config
    # here says what was run and keeps agreeing if the method changes.
    split_method = nsv._cfg("splitting.method", "warhead_dbscan")
    sub_split = (nsv._cfg("splitting.stage2.max_sub", psub.DEFAULT_MAX_SUB)
                 if nsv._cfg("splitting.stage2.enabled", True) else 1)
    log.info("  splitting %s, sub_split %s, seed %s, nrun %d",
             split_method, sub_split, args.seed, args.nrun)
    poses, aggs = nsv.one(cand, rec_dir, plain_rec, args.nrun,
                          args.gpu, do_gnina=False, all_poses=True,
                          sub_split=sub_split, seed=args.seed,
                          split_method=split_method)

    pd.DataFrame(poses).to_csv(OUT.write(f"poses_{args.ident}", ".csv"), index=False)
    a = pd.DataFrame(aggs)
    dest = OUT.write(f"agg_{args.ident}", ".csv")
    a.to_csv(dest, index=False)
    print(f"\n  {len(a)} modes over {args.nrun} poses -> {dest}")
    cols = [c for c in ("ident", "mode", "n_poses_mode", "consensus",
                        "viable_fraction", "enrichment", "mean_energy")
            if c in a.columns]
    print(a[cols].to_string(index=False))


if __name__ == "__main__":
    main()
