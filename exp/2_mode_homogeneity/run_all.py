#!/usr/bin/env python3
"""
Purpose: does letting the 2 A cut govern actually make a mode one pose?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18
Input: --candidate, --topic (reads the screen's own pose cloud; no GPU, no docking)
Output: append_only/00_outputs/blacksmith/mode_homogeneity_<candidate>/

@tt8804: "why would only 20 out of 82 poses be good?? that suggests that they
arent the same poses. I dont think you are splitting correctly, modes should be
essentially the same pose within a few A."

THE CLAIM UNDER TEST. If a mode is one pose, its poses either reach attack
geometry or they do not -- viable fraction sits near 0 or near 1. Measured on
nac_v5 it sits between 0.1 and 0.9 for 42% of modes, and the median mode spans
3.51 A in warhead-anchor distance. So "fraction viable" has been reporting the
mixing ratio of two populations rather than a property of a pose.

WHAT CHANGED. `pose_subsplit.subdivide` cut by diameter and then, whenever that
asked for more than `max_sub` clusters, threw the answer away and re-cut with
`maxclust` at max_sub -- merging poses the cut had just separated. It also folded
strays into the LARGEST sub-cluster rather than the nearest. Both are fixed; this
measures whether the fix does what it claims, and picks `max_sub` on evidence
rather than on the recovery benchmark, which asked a different question.

This re-splits poses that are already on disk. It does not re-dock, so the only
thing varying is the splitting rule.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_subsplit as psub                # noqa: E402
from shared import run_paths as rp                      # noqa: E402
from shared import target_config as tc                  # noqa: E402

log = logging.getLogger("mode-homogeneity")


def heavy_coords(sdf: Path) -> dict[int, np.ndarray]:
    """Heavy-atom coordinates per pose, keyed by 0-based pose index."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = {}
    for i, m in enumerate(Chem.SDMolSupplier(str(sdf), removeHs=True, sanitize=False)):
        if m is None:
            continue
        c = m.GetConformer()
        out[i] = np.array([list(c.GetAtomPosition(a.GetIdx()))
                           for a in m.GetAtoms() if a.GetAtomicNum() > 1], float)
    return out


def poses_table(topic: str, cand: str) -> pd.DataFrame:
    fs = glob.glob(str(rp.BLACKSMITH / topic / "poses_s*_*.csv"))
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d[d.ident == cand].sort_values("pose_idx").reset_index(drop=True)


def parent_of(topic: str, cand: str) -> dict[int, int]:
    """final mode -> stage-1 parent mode, from the rank table."""
    fs = sorted(glob.glob(str(rp.BLACKSMITH / f"rank_v2/rank_v2_*_{topic}_*.csv")))
    if not fs:
        return {}
    d = pd.read_csv(fs[-1])
    d = d[d.parent_ident == cand]
    if "parent_mode" not in d.columns:
        return {}
    return {int(r["mode"]): int(r.parent_mode) for _, r in d.iterrows()
            if pd.notna(r.get("parent_mode"))}


def homogeneity(groups: dict[int, np.ndarray], p: pd.DataFrame) -> dict:
    """Span, and how bimodal viability is, over modes with >= MIN poses."""
    spans, vfs, sizes = [], [], []
    for _, rows in groups.items():
        if len(rows) < 12:
            continue
        g = p.iloc[rows]
        spans.append(float(g.distance.max() - g.distance.min()))
        vfs.append(float(g.viable.mean()))
        sizes.append(len(g))
    if not spans:
        return {"modes": 0}
    spans, vfs = np.array(spans), np.array(vfs)
    return {"modes": len(spans),
            "median_span_a": float(np.median(spans)),
            "frac_span_over_2a": float((spans > 2).mean()),
            "frac_mixed_vf": float(((vfs > 0.1) & (vfs < 0.9)).mean()),
            "median_size": float(np.median(sizes))}


def split_with(parents: dict[int, list[int]], coords: dict[int, np.ndarray],
               max_sub: int) -> dict[int, np.ndarray]:
    """Re-split every stage-1 parent at a given cap; return final groups."""
    out: dict[int, np.ndarray] = {}
    nxt = 0
    for pm, rows in parents.items():
        rows = np.array(sorted(rows))
        c = np.array([coords[i] for i in rows])
        lab = np.zeros(len(rows), int)
        sub, _ = psub.subdivide(lab, c, max_sub=max_sub,
                               min_size=psub.MIN_MODE_SIZE)
        for s in np.unique(sub):
            if s < 0:
                continue
            out[nxt] = rows[sub == s]
            nxt += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--topic", default=tc.topic())
    ap.add_argument("--caps", default="5,8,12,20")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    p = poses_table(args.topic, args.candidate)
    sdf = rp.BLACKSMITH / f"{args.topic}_allposes" / f"{args.candidate}.sdf"
    coords = heavy_coords(sdf)
    n = min(len(p), len(coords))
    p = p.iloc[:n]
    log.info("%s: %d poses with coordinates", args.candidate, n)

    # AS SHIPPED: whatever the screen actually produced, read back off disk.
    shipped = {int(m): np.flatnonzero((p["mode"] == m).to_numpy())
               for m in sorted(p["mode"].unique()) if m >= 0}
    rows = [dict(rule="as shipped (nac_v5)", **homogeneity(shipped, p))]

    # Stage-1 parents, reconstructed so the re-split starts where the screen did.
    pm = parent_of(args.topic, args.candidate)
    parents: dict[int, list[int]] = {}
    for m, idxs in shipped.items():
        parents.setdefault(pm.get(m, m), []).extend(idxs.tolist())
    log.info("%d stage-1 parent mode(s), sizes %s", len(parents),
             sorted(len(v) for v in parents.values()))

    for cap in [int(x) for x in args.caps.split(",")]:
        g = split_with(parents, coords, cap)
        rows.append(dict(rule=f"cut governs, cap {cap}", **homogeneity(g, p)))

    out = rp.BLACKSMITH / f"mode_homogeneity_{args.candidate}"
    out.mkdir(parents=True, exist_ok=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / "homogeneity_1.csv", index=False)
    print(f"\n  {args.candidate}, {n} poses, cut = {psub.CUT_A} A\n")
    print(t.to_string(index=False))
    print("\n  frac_mixed_vf is the headline: the share of modes whose viable")
    print("  fraction is neither near 0 nor near 1 -- i.e. still a mixture.")
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
