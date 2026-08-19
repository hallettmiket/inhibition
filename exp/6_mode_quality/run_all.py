#!/usr/bin/env python3
"""
Purpose: are the modes RIGHT? Width, purity, and whether the validated pose sits in a clean one.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-19
Input: --molecules, --sets (reuses the seeded replicate screens already on disk)
Output: append_only/00_outputs/blacksmith/mode_quality/

@tt8804: "can you stop talking about the ranking and focus on getting the right
modes".

NO SCORING HERE, DELIBERATELY. Every earlier comparison mixed clustering and
ranking together and could not say which one was responsible for what. This asks
only whether a clustering rule produces MODES -- groups of poses that are the
same pose -- and it does so with three measurements that do not involve a score:

  WIDTH     the widest heavy-atom RMSD between any two poses in a mode. A mode
            is supposed to be "essentially the same pose within a few A"
            (@tt8804), so this is the definition made numeric.

  PURITY    a pose either reaches attack geometry or does not. If a mode is one
            pose, its poses agree -- the fraction that are viable sits near 0 or
            near 1. A mode at 0.4 is two populations wearing one label, and its
            viable fraction is a mixing ratio rather than a property.
            Viability is used as a LABEL to test homogeneity against, never as
            an input to the clustering.

  HOME      the mode that contains the pose a 100 ns trajectory validated: how
            wide it is, how big, and whether it is pure. A rule can look good on
            averages and still put the one pose we know about in a bag.

The methods differ ONLY in how poses are grouped. Nothing here ranks anything.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import pose_cluster as pclust                # noqa: E402
from shared import pose_modes as pmod                    # noqa: E402
from shared import pose_subsplit as psub                 # noqa: E402
from shared import run_paths as rp                       # noqa: E402

log = logging.getLogger("mode-quality")

sys.path.insert(0, str(REPO / "exp" / "4_election"))


def _election_mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "el4", REPO / "exp" / "4_election" / "run_all.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def group(feat, heavy, method: str) -> np.ndarray:
    """Labels per pose. The ONLY thing that differs between the arms."""
    if method == "hdbscan":
        return pclust.cluster(heavy)
    lab = pmod.split(feat, method=("complete" if method == "complete" else "dbscan"))
    if (lab >= 0).sum() == 0:
        return lab
    kw = dict(max_sub=None, min_sub_size=3, cut_a=1.0) if method == "fine" else {}
    sub, _ = psub.subdivide(lab, heavy, **kw)
    return sub


def quality(lab, heavy, tab, ref_idx: int | None) -> dict:
    gs = sorted({int(x) for x in lab if x >= 0})
    if not gs:
        return {"modes": 0}
    via = tab.viable.to_numpy().astype(bool)
    widths, purity, sizes = [], [], []
    for c in gs:
        m = lab == c
        widths.append(float(pclust.rmsd_matrix(heavy[m]).max()))
        f = float(via[m].mean())
        purity.append(f < 0.1 or f > 0.9)
        sizes.append(int(m.sum()))
    out = {"modes": len(gs), "noise": float((lab == -1).mean()),
           "width_med": float(np.median(widths)),
           "width_p90": float(np.quantile(widths, 0.9)),
           "width_max": float(max(widths)),
           "pure_frac": float(np.mean(purity)),
           "largest": int(max(sizes))}
    if ref_idx is not None and lab[ref_idx] >= 0:
        m = lab == lab[ref_idx]
        f = float(via[m].mean())
        out.update({"home_n": int(m.sum()),
                    "home_width": float(pclust.rmsd_matrix(heavy[m]).max()),
                    "home_viable_frac": f,
                    "home_pure": bool(f < 0.1 or f > 0.9)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--molecules",
                    default="t4_716800c125a7,t4_80fbed3bdf1e,t4_b306425b6a73")
    ap.add_argument("--sets", default="1000,2000|3000,4000|5000,6000",
                    help="seed0 groups, one per molecule, '|'-separated")
    ap.add_argument("--methods", default="dbscan,complete,fine,hdbscan")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")
    el = _election_mod()

    rows = []
    for mol, sets in zip(args.molecules.split(","), args.sets.split("|")):
        ref = el.reference_for(mol)
        for s0 in sets.split(","):
            for r in range(1, 6):
                topic = f"election_{mol}_s{s0}_r{r}"
                try:
                    smarts = el.smarts_for(mol, topic)
                    _, feat, heavy, tab = el.cloud(topic, mol, smarts)
                except Exception as exc:                     # noqa: BLE001
                    log.warning("%s: %s", topic, exc)
                    continue
                # WHICH POSE IS THE VALIDATED ONE, in this replicate: the pose
                # closest to the geometry the trajectory confirmed. Matched on
                # the warhead's position and direction, which is what the
                # reference records -- and used ONLY to ask which mode it landed
                # in, never to form modes.
                dpos = np.linalg.norm(feat[:, :3] - ref["centroid"], axis=1)
                cosa = np.clip(feat[:, 3:] @ ref["direction"], -1, 1)
                score = dpos + 2.0 * np.arccos(cosa)
                ref_idx = int(np.argmin(score)) if score.min() < 2.5 else None
                for meth in args.methods.split(","):
                    q = quality(group(feat, heavy, meth), heavy, tab, ref_idx)
                    rows.append(dict(molecule=mol, set=s0, rep=r, method=meth,
                                     ref_found=ref_idx is not None, **q))
        log.info("%s done", mol)

    d = pd.DataFrame(rows)
    out = rp.BLACKSMITH / "mode_quality"
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "quality_1.csv", index=False)

    print(f"\n  {len(d)} (replicate x method) rows over "
          f"{d.molecule.nunique()} molecules\n")
    print(f"  {'method':<10}{'modes':>7}{'noise':>7}{'width med':>11}{'width p90':>11}"
          f"{'widest':>8}{'pure':>7}{'largest':>9}")
    for m, g in d.groupby("method"):
        print(f"  {m:<10}{g.modes.mean():>7.1f}{100*g.noise.mean():>6.0f}%"
              f"{g.width_med.mean():>11.2f}{g.width_p90.mean():>11.2f}"
              f"{g.width_max.mean():>8.2f}{100*g.pure_frac.mean():>6.0f}%"
              f"{g.largest.mean():>9.0f}")
    print(f"\n  THE MODE HOLDING THE VALIDATED POSE")
    print(f"  {'method':<10}{'found':>8}{'size':>7}{'width':>8}{'viable frac':>13}{'pure':>7}")
    for m, g in d.groupby("method"):
        h = g.dropna(subset=["home_n"])
        if h.empty:
            print(f"  {m:<10}{'0':>8}"); continue
        print(f"  {m:<10}{len(h):>4}/{len(g):<3}{h.home_n.mean():>7.0f}"
              f"{h.home_width.mean():>8.2f}{h.home_viable_frac.mean():>13.2f}"
              f"{100*h.home_pure.mean():>6.0f}%")
    (out / "quality_1.json").write_text(json.dumps(
        {"rows": len(d), "molecules": sorted(d.molecule.unique())}, indent=2))
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
