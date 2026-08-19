#!/usr/bin/env python3
"""
Purpose: does bounding the mode's DIAMETER instead of its LINK give modes that behave as one pose?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-18
Input: --candidates (screened fresh into their own topic), --eps
Output: append_only/00_outputs/blacksmith/linkage_compare/

D0086 step 2. Stage 1 is DBSCAN at eps 3.0, which bounds the distance to the
NEAREST neighbour: A-B-C-D each within 3 A of the next is one mode however wide
the chain grows, and measured on nac_v5 a stage-1 parent spans 4.22 A (median)
and 83 deg. Complete linkage at the same 3 A bounds the FURTHEST pair instead,
so no two poses in a mode exceed the tolerance already chosen.

THE TEST. A mode that is one pose has poses that either reach attack geometry or
do not, so its viable fraction sits near 0 or near 1. `pure` is the share of
modes where it does. Reported beside the cost -- modes per molecule, and what
share of the cloud stays assigned -- because a rule that purifies by discarding
the cloud has not helped (a 0.1 A RMSD cut scores 100% pure on 4% of the poses).

WHY IT RE-SCREENS. The stored clouds cannot be joined to their own per-pose
tables (#76): the SDF numbered poses by position and was not rewritten on a
re-run. Both are fixed, so a fresh screen produces a cloud that carries
`pose_idx` and matches the table beside it. ~30 s per molecule on one GPU.
"""

from __future__ import annotations

import argparse
import glob
import logging
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_modes as pmod                  # noqa: E402
from shared import run_paths as rp                     # noqa: E402

log = logging.getLogger("linkage")


def screen(cands: list[str], topic: str, gpu: str, nrun: int) -> None:
    only = Path(f"/tmp/linkage_{topic}.txt")
    only.write_text("\n".join(cands) + "\n")
    cmd = [sys.executable, str(REPO / "scripts/nac_screen_v2.py"),
           "--only", str(only), "--topic", topic, "--nrun", str(nrun),
           "--gpu", gpu, "--all-poses"]
    r = subprocess.run(["nice", "-n", "19"] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])


def cloud(topic: str, cand: str):
    """(feature matrix, per-pose table) aligned on pose_idx.

    Alignment is the whole reason this experiment can run: the SDF now carries
    `pose_idx`, so a pose in the cloud and its row in the table are the same
    object rather than two things at the same offset.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    fs = glob.glob(str(rp.BLACKSMITH / topic / "poses_s*_*.csv"))
    p = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    p = p[p.ident == cand].set_index("pose_idx")
    sdf = rp.BLACKSMITH / f"{topic}_allposes" / f"{cand}.sdf"
    idx, pos = [], []
    for m in Chem.SDMolSupplier(str(sdf), removeHs=True, sanitize=False):
        if m is None or not m.HasProp("pose_idx"):
            continue
        i = int(m.GetProp("pose_idx"))
        if i not in p.index:
            continue
        idx.append(i)
        pos.append(np.array(m.GetConformer().GetPositions()))
    if not idx:
        raise RuntimeError(f"{cand}: no pose in the cloud carries a usable pose_idx")
    return np.array(idx), np.array(pos), p.loc[idx]


def score(labels: np.ndarray, tab: pd.DataFrame) -> dict | None:
    keep = [c for c in np.unique(labels) if c >= 0]
    groups = [np.flatnonzero(labels == c) for c in keep]
    groups = [g for g in groups if len(g) >= 3]
    if not groups:
        return None
    d = tab.distance.to_numpy()
    v = tab.viable.to_numpy().astype(float)
    span = np.array([d[g].max() - d[g].min() for g in groups])
    vf = np.array([v[g].mean() for g in groups])
    return {"modes": len(groups),
            "assigned": sum(len(g) for g in groups) / len(labels),
            "median_span_a": float(np.median(span)),
            "worst_span_a": float(span.max()),
            "pure": float(((vf < 0.1) | (vf > 0.9)).mean()),
            "largest": int(max(len(g) for g in groups))}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidates", default=None,
                    help="comma-separated; default: a spread of nac_v5 molecules")
    ap.add_argument("--topic", default="linkage_compare")
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--nrun", type=int, default=500)
    ap.add_argument("--eps", type=float, default=pmod.DEFAULT_EPS)
    ap.add_argument("--skip-screen", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    if args.candidates:
        cands = args.candidates.split(",")
    else:
        fs = glob.glob(str(rp.BLACKSMITH / "nac_v5" / "poses_s*_*.csv"))
        allc = sorted(pd.concat([pd.read_csv(f, usecols=["ident"]) for f in fs],
                                ignore_index=True).ident.unique())
        cands = ["t4_716800c125a7"] + [c for c in allc[::37] if c != "t4_716800c125a7"][:11]

    if not args.skip_screen:
        log.info("screening %d molecule(s) into %s", len(cands), args.topic)
        screen(cands, args.topic, args.gpu, args.nrun)

    rows = []
    for c in cands:
        try:
            _, pos, tab = cloud(args.topic, c)
        except Exception as exc:                          # noqa: BLE001
            log.warning("%s: %s", c, exc)
            continue
        # The stage-1 feature needs the reactive-atom match, which the screen
        # already resolved; re-deriving it here would be a second answer to the
        # same question. Cluster the whole-pose coordinates' warhead proxy
        # instead: the per-pose distance and angle ARE the criterion, so they
        # cannot be used; what is available and independent is the pose geometry
        # itself, which is what `distances` consumes via features. So compare on
        # the ONE thing that differs between the two rules -- the linkage -- by
        # feeding both the identical distance matrix.
        n = len(pos)
        flat = pos.reshape(n, -1)
        dist = np.sqrt(((flat[:, None, :] - flat[None, :, :]) ** 2)
                       .reshape(n, n, -1).mean(axis=2))
        ms = max(3, int(round(pmod.MIN_POPULATION_FRAC * n)))
        for name, lab in (("dbscan (link)", pmod._dbscan(dist, args.eps, ms)),
                          ("complete (diameter)",
                           pmod._complete_linkage(dist, args.eps, ms))):
            s = score(lab, tab)
            if s:
                rows.append(dict(candidate=c, rule=name, **s))

    if not rows:
        log.error("nothing scored")
        return 1
    d = pd.DataFrame(rows)
    out = rp.BLACKSMITH / "linkage_compare"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "per_molecule_1.csv", index=False)
    agg = d.groupby("rule").agg(molecules=("candidate", "nunique"),
                                modes=("modes", "mean"),
                                assigned=("assigned", "mean"),
                                median_span=("median_span_a", "mean"),
                                worst_span=("worst_span_a", "mean"),
                                largest=("largest", "mean"),
                                pure=("pure", "mean"))
    print(f"\n  eps = {args.eps} (a neighbour radius for dbscan, a diameter for complete)\n")
    print(agg.to_string())
    print("\n  pure = share of modes whose viable fraction is <0.1 or >0.9")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
