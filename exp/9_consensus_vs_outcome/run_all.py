#!/usr/bin/env python3
"""
Purpose: does a pose that shows up MANY TIMES actually do better in dynamics?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-17
Input: the nac_v5 sweep + 100 ns results, and each swept molecule's pose cloud
Output: 00_outputs/blacksmith/consensus_vs_outcome_nac_v5/

@tt8804: "I feel like a pose that shows up many times is more likely to be the
real pose but I guess not ... maybe rerun a test with our tighter pose
splitting."

THE INTUITION IS THE STANDARD ONE and it is what `consensus` encodes: a docking
that keeps returning to the same place is telling you something. The project has
measured against it twice already -- D0071 (neither ranking metric predicts pose
stability) and D0073 (consensus DEPLETES validated mechanisms) -- but both were
measured under the SHIPPED clustering, where a "mode" routinely holds 65-86% of
the cloud. A mode that large has a consensus near 1 whatever the molecule does,
so the test may have been measuring the clustering rather than the intuition.

THIS RE-ASKS IT UNDER TIGHT CLUSTERING. For every swept mode, the pose that was
actually simulated is located in an HDBSCAN re-clustering of its own cloud, and
that group's SIZE is correlated against what the trajectory then did. Loose and
tight are computed on the same molecules and the same outcomes, so the
comparison isolates the clustering.

NOISE IS A SINGLETON, NOT A DISCARD (@tt8804: "treat the noise as singleton
poses and reframe clustering as a cost saving mechanism to collapse the number
of poses"). A pose HDBSCAN calls noise is a group of one -- size 1 -- not a pose
that failed to exist. Under that reading clustering stops being an ontology of
binding modes and becomes de-duplication for a simulation budget, and "how many
poses agree with this one" stays defined for every pose in the cloud.

OUTCOMES, both pre-existing:
  frac_attack_ready   fraction of the 5 ns triage in attack geometry
  explicit_ligand_rmsd_nm_max   how far the ligand wandered over 100 ns
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import pose_cluster as pc               # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("consensus-vs-outcome")


def sweep_table() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.sweep_dir() / "attack_sweep_*.csv")),
                key=os.path.getmtime)
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d = d.drop_duplicates("ident", keep="last")
    return d[d.status.astype(str).str.startswith("ok")]


def md_table() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.residence_dir() / "*.csv")), key=os.path.getmtime)
    if not fs:
        return pd.DataFrame()
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    keep = [c for c in ("ident", "parent_ident", "production_ps",
                        "explicit_ligand_rmsd_nm_max",
                        "explicit_ligand_rmsd_nm_mean") if c in d.columns]
    d = d[keep]
    d = d[d.production_ps >= 50_000] if "production_ps" in d else d
    return d.drop_duplicates("ident", keep="last")


def cloud_and_rep(ident: str, pose_rank: int):
    """(cloud coords, the simulated pose's coords) — both heavy atoms only.

    The simulated pose is taken from the REPRESENTATIVES file by `pose_rank`,
    which is how the sweep asked for it, then located inside the full cloud by
    coordinates. Located by geometry rather than by index because the two files
    are written in different orders and an index would silently pair the wrong
    pose (`build_plan_next.md` §1.6a).
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    cloud_f = rp.allposes_dir() / f"{ident}.sdf"
    rep_f = rp.poses_dir() / f"{ident}.sdf"
    if not cloud_f.is_file() or not rep_f.is_file():
        return None, None
    ms = [m for m in Chem.SDMolSupplier(str(cloud_f), removeHs=False,
                                        sanitize=False) if m is not None]
    if not ms:
        return None, None
    heavy = [i for i, a in enumerate(ms[0].GetAtoms()) if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[heavy] for m in ms])
    want = None
    for m in Chem.SDMolSupplier(str(rep_f), removeHs=False, sanitize=False):
        if m is None or not m.HasProp("pose_rank"):
            continue
        if int(m.GetProp("pose_rank")) == int(pose_rank):
            want = m.GetConformer().GetPositions()[heavy]
            break
    return xyz, want


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--limit", type=int, default=0, help="0 = all swept modes")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    sw = sweep_table()
    md = md_table()
    log.info("swept modes with an ok result: %d;  100 ns rows: %d", len(sw), len(md))

    rows = []
    items = list(sw.itertuples())
    if a.limit:
        items = items[:a.limit]
    for n, r in enumerate(items, 1):
        parent = str(r.parent_ident)
        xyz, want = cloud_and_rep(parent, r.pose_rank)
        if xyz is None or want is None:
            continue
        # which cloud pose IS the simulated one
        d0 = np.sqrt(((xyz - want) ** 2).sum(axis=2).mean(axis=1))
        j = int(np.argmin(d0))
        if d0[j] > 0.05:                       # not the same pose; refuse to guess
            continue
        lab = pc.cluster(xyz)
        k = int(lab[j])
        # NOISE IS A SINGLETON. size 1, not "absent".
        tight = 1 if k == -1 else int((lab == k).sum())
        rows.append(dict(
            ident=str(r.ident), parent_ident=parent, mode=r.mode,
            n_poses=len(xyz),
            tight_size=tight, tight_consensus=tight / len(xyz),
            was_noise=(k == -1),
            frac_attack_ready=float(r.frac_attack_ready),
            n_visits=float(r.n_visits)))
        if n % 25 == 0:
            log.info("  %d/%d", n, len(items))

    d = pd.DataFrame(rows)
    if d.empty:
        raise SystemExit("nothing joined")

    # the SHIPPED consensus, from the run's own aggregates
    agg = pd.concat([pd.read_csv(f) for f in
                     sorted(glob.glob(str(rp.BLACKSMITH / rp.topic() / "agg_s*.csv")))],
                    ignore_index=True)
    agg = agg.drop_duplicates("ident", keep="last")[["ident", "consensus",
                                                     "n_poses_mode"]]
    d = d.merge(agg.rename(columns={"consensus": "loose_consensus",
                                    "n_poses_mode": "loose_size"}),
                on="ident", how="left")
    if not md.empty:
        d = d.merge(md, on="ident", how="left")

    t = sout.Topic("blacksmith", "consensus_vs_outcome_nac_v5")
    d.to_csv(t.write("joined", ".csv"), index=False)

    from scipy.stats import spearmanr
    print("\n" + "=" * 74)
    print("  DOES A POSE THAT SHOWS UP MANY TIMES DO BETTER IN DYNAMICS?")
    print("=" * 74)
    print(f"\n  swept modes joined: {len(d)}   "
          f"of which HDBSCAN called the simulated pose noise: "
          f"{int(d.was_noise.sum())} ({d.was_noise.mean()*100:.0f}%)")
    print(f"  group size holding the simulated pose — "
          f"loose median {d.loose_size.median():.0f}, tight median {d.tight_size.median():.0f}")

    for outcome, label in (("frac_attack_ready", "5 ns attack-ready fraction"),
                           ("explicit_ligand_rmsd_nm_max", "100 ns max ligand RMSD")):
        if outcome not in d.columns:
            continue
        sub = d[d[outcome].notna()]
        if len(sub) < 8:
            print(f"\n  {label}: only {len(sub)} rows, skipped")
            continue
        print(f"\n  vs {label}  (n = {len(sub)})")
        for col, name in (("loose_consensus", "consensus, SHIPPED clustering"),
                          ("tight_consensus", "consensus, HDBSCAN (tight)  "),
                          ("tight_size", "group size, HDBSCAN (tight) ")):
            if col not in sub or sub[col].notna().sum() < 8:
                continue
            rho, p = spearmanr(sub[col], sub[outcome], nan_policy="omit")
            flag = "" if p >= 0.05 else "   *"
            print(f"    {name}  rho = {rho:+.3f}   p = {p:.3f}{flag}")
    print()


if __name__ == "__main__":
    main()
