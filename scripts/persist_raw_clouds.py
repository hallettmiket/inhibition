#!/usr/bin/env python3
"""
Purpose: dock molecules and persist the RAW pose cloud, unfiltered by any clustering
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: --molecules (default: exactly the 20 exp/16 sampled), the reactive 3IKD receptor
Output: 00_outputs/blacksmith/raw_cloud_<ident>/cloud_*.sdf

D0093's FIRST FIX. `<topic>_allposes` holds only poses whose DBSCAN label is in
`mode_ids` (nac_screen_v2.py:501), so ~21% of every production cloud is absent --
and it is the scattered 21%, the poses hardest to group. exp/14, exp/15 and
exp/16 all read that path, which means every candidate replacement for DBSCAN so
far was evaluated on clouds DBSCAN had already cleaned.

This writes what the name `allposes` promises: every pose the docking returned,
with nothing removed. Same receptor, same box, same run count as production -- the
ONLY difference is that no clustering runs before the write, so the comparison
against the filtered clouds isolates the filter.

THE SAME MOLECULES exp/16 SAMPLED, by reproducing its selection rather than
drawing fresh ones. A re-measurement on a different sample would confound the
filter with the draw, and the whole point is to isolate one of them.

DOCKING IS STOCHASTIC AND SEEDED. `docking.seed` is set (#77), so a re-dock of the
same molecule is reproducible; but a cloud persisted here is NOT the cloud the
existing scores were computed from, and must not be shown beside them (#44). It
exists to re-run the CLUSTERING experiments, which need only the geometry.
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import nac_screen as ns                             # noqa: E402
import nac_rank as nr                               # noqa: E402
from shared import outputs as sout                  # noqa: E402
from shared import run_paths as rp                  # noqa: E402
from shared import target_config as tc              # noqa: E402

log = logging.getLogger("raw-clouds")


def exp16_selection(n: int, seed: int) -> list[str]:
    """The molecules exp/16 used, by re-running its draw -- not a fresh sample."""
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{rp.topic()}_conditional_eb_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit("no rank_v2 table; cannot reproduce exp/16's selection")
    rk = pd.read_csv(fs[-1])
    mols = [m for m in rk.parent_ident.dropna().unique()
            if (rp.allposes_dir() / f"{m}.sdf").is_file()]
    rng = np.random.default_rng(seed)
    return [str(x) for x in rng.choice(mols, size=min(n, len(mols)), replace=False)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--molecules", default="", help="comma-separated idents")
    ap.add_argument("--n-molecules", type=int, default=20)
    ap.add_argument("--n-runs", type=int, default=0, help="default: docking.n_runs")
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--skip-existing", action="store_true", default=True)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    n_runs = a.n_runs or int(tc.get("docking.n_runs"))
    want = ([x.strip() for x in a.molecules.split(",") if x.strip()]
            if a.molecules else exp16_selection(a.n_molecules, a.seed))
    log.info("%d molecules, %d runs each, GPU %s", len(want), n_runs, a.gpu)

    cands = {c.ident: c for c in nr.load_candidates()}
    missing = [m for m in want if m not in cands]
    if missing:
        log.warning("%d not in the candidate table, skipped: %s",
                    len(missing), ", ".join(missing[:4]))
    want = [m for m in want if m in cands]

    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    from rdkit import Chem

    done: list = []
    failed: list = []
    for i, ident in enumerate(want, 1):
        outdir = rp.BLACKSMITH / f"raw_cloud_{ident}"
        if a.skip_existing and list(outdir.glob("cloud_*.sdf")):
            log.info("[%2d/%d] %s already persisted", i, len(want), ident)
            done.append(ident)
            continue
        work = Path(tempfile.mkdtemp(prefix=f"raw_{ident[:12]}_"))
        try:
            ligs = list(ns.prepare_ligand(cands[ident], work / "lig.pdbqt"))
            if not ligs:
                raise RuntimeError("ligand preparation produced nothing")
            dlg = ns.dock(ligs[0], rec_dir, work / "d", n_runs, a.gpu, seed=a.seed)
            mol, match = ns.rebuild_and_match(dlg, cands[ident])
            n_conf = mol.GetNumConformers()
            if n_conf < n_runs * 0.5:
                raise RuntimeError(f"only {n_conf} of {n_runs} poses came back")
            tp = sout.Topic("blacksmith", f"raw_cloud_{ident}")
            f = tp.write("cloud", ".sdf")
            w = Chem.SDWriter(str(f))
            for cid in range(n_conf):
                # NO FILTER OF ANY KIND between the dock and this write.
                mol.SetProp("pose_index", str(cid))
                w.write(mol, confId=cid)
            w.close()
            log.info("[%2d/%d] %s: %d poses -> %s", i, len(want), ident, n_conf, f.name)
            done.append(ident)
        except Exception as exc:                                # noqa: BLE001
            log.error("[%2d/%d] %s FAILED: %s", i, len(want), ident, str(exc)[:140])
            failed.append((ident, str(exc)[:140]))
        finally:
            import shutil
            shutil.rmtree(work, ignore_errors=True)

    print("\n" + "=" * 74)
    print(f"  RAW CLOUDS PERSISTED — {len(done)} ok, {len(failed)} failed")
    print("=" * 74)
    for ident, why in failed:
        print(f"    {ident}: {why}")
    print(f"\n  {rp.BLACKSMITH}/raw_cloud_<ident>/cloud_*.sdf")
    print("  These are UNFILTERED. Do not show them beside scores computed from")
    print("  the production clouds (#44) -- they are a different docking run.\n")
    # A WHOLESALE FAILURE MUST NOT EXIT 0. The first run of this hit
    # `No module named 'gemmi'` on all 20 molecules -- the wrong environment --
    # printed "0 ok, 20 failed", and RETURNED SUCCESS, so the chain step after it
    # ran against nothing and also reported success. A downstream stage cannot
    # tell an empty result from a clean one unless this says so.
    if not done:
        raise SystemExit(
            f"no cloud was persisted ({len(failed)} failures). This needs the "
            "docking environment: ~/.micromamba/envs/dwi_reactive/bin/python")


if __name__ == "__main__":
    main()
