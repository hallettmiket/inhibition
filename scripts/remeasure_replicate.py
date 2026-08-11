#!/usr/bin/env python3
"""Re-measure a replicate whose trajectory finished but whose row said "failed".

Why this exists. `run_one` measured residence at `md/rep1` regardless of which
replicate had been asked for, while `gromacs_explicit.run_pipeline` writes
replicate N to `md/rep<N>`. Replicate 2 of t4_716800c125a7 therefore ran its full
100 ns, wrote 962 MB of trajectory, and produced a row reading

    failed: ResidenceError: no production trajectory at .../md/rep1/prod.xtc

The run was fine. The path was wrong. That is fixed in `md_residence_3ikd.py`
(the directory is now taken from what run_pipeline reports it wrote), but a run
already in flight has the old code loaded, and the trajectories are on disk
because `--keep` was passed. This measures them where they actually are, so a
completed 100 ns is not thrown away over a filename.

Usage:

    python3 scripts/remeasure_replicate.py \\
        --dir /data/lab_vm/modifiable/inhibition/md_residence_3ikd_rep2/\\
t4_716800c125a7/md/rep2 \\
        --candidate t4_716800c125a7 --replicate 2

Writes an integer-versioned CSV alongside every other md_residence row. It does
NOT overwrite the failed row: that row is the honest record that the run was
mismeasured, and the versioning convention is to add, not to rewrite.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import pandas as pd                                        # noqa: E402

import md_residence_3ikd as mr                             # noqa: E402

log = logging.getLogger("remeasure")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", required=True,
                    help="the rep<N> directory the run actually wrote")
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--replicate", type=int, required=True)
    ap.add_argument("--production-ps", type=float, default=100000.0)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    d = Path(args.dir)
    traj = d / "prod.xtc"
    if not traj.is_file():
        raise SystemExit(f"no trajectory at {traj}")

    row = {"ident": args.candidate, "label": "candidate",
           "replicate": args.replicate, "production_ps": args.production_ps,
           "status": "ok",
           # Say plainly where the numbers came from. A row that was rescued
           # from a mismeasured run should not be indistinguishable from one
           # that was measured correctly the first time.
           "remeasured_from": str(d),
           "remeasure_reason": "run_one measured md/rep1 regardless of "
                               "--replicate; trajectory was intact"}
    row.update(mr.measure_residence(d))

    dest = mr.OUT.write(
        f"md_residence_rep{args.replicate}_{args.candidate}_remeasured", ".csv")
    pd.DataFrame([row]).to_csv(dest, index=False)
    print(f"\n  ok -> {dest}")
    for k in ("ns_analysed", "explicit_frac_frames_engaged",
              "explicit_ligand_rmsd_nm_mean", "explicit_ligand_rmsd_nm_max",
              "explicit_ligand_rmsd_nm_final"):
        if k in row:
            print(f"    {k} = {row[k]}")


if __name__ == "__main__":
    main()
