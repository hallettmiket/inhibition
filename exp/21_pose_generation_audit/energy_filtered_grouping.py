#!/usr/bin/env python3
"""
Purpose: do the grouping and saturation conclusions survive an energy filter?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the RAW clouds re-persisted WITH per-pose energies
Output: 00_outputs/blacksmith/pose_generation_audit/

EVERY CLUSTERING RESULT SO FAR WAS MEASURED ON UNFILTERED CLOUDS. exp/16, 17, 19
and 20 read pose sets that carried no energies at all, so the best-scoring pose
and the 500th counted equally. exp/21 shows the tail is real docking output that
the SCORE ranks correctly -- poses with >30% of atoms uncontacted sit at the 88th
energy percentile. The question this asks is whether including them is what
produced the two headline findings:

  * "no molecule has a consensus pose" (exp/20: top group holds 1.4-7.2%)
  * "the group count never saturates" (exp/17: b = +0.69)

Both would read very differently if they are properties of the tail rather than
of the cloud.

THE SIZE-MATCHED CONTROL IS THE WHOLE EXPERIMENT. Keeping the best 25% leaves
125 poses instead of 500, and a smaller sample concentrates its top group for
purely arithmetic reasons -- so "the top-1 share rose" proves nothing on its own.
The first version of this script printed exactly that under a heading claiming
concentration was "not just a smaller n", which it had not measured. Every
filtered number is therefore compared against a RANDOM subset of identical size
from the same cloud.

THE FILTER IS A RANK, NOT A THRESHOLD. An absolute kcal/mol cut would keep 90% of
one molecule and 5% of another, since the energy scale shifts with size -- so
"the best N%" is used, and the count of poses kept is reported beside every
number so a change in concentration cannot be confused with a change in n.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("energy-filter")
KEEP = [1.00, 0.50, 0.25, 0.10]


def _by_path(name, rel):
    import importlib.util
    sp = importlib.util.spec_from_file_location(name, REPO / rel)
    m = importlib.util.module_from_spec(sp)
    sys.modules[name] = m
    sp.loader.exec_module(m)
    return m


_M = _by_path("e17", "exp/17_contact_saturation/run_all.py")


def load_with_energy(f: Path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ms = [m for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=True)
          if m is not None]
    if len(ms) < 30 or not ms[0].HasProp("free_energy_kcal"):
        return None, None, None
    hv = [a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[hv] for m in ms])
    en = np.array([float(m.GetProp("free_energy_kcal")) for m in ms])
    return xyz, en, (ms[0], hv)


def stats(D, tol, xyz):
    lab = pc.group(D, tol)
    sz = np.bincount(lab).astype(float)
    p = sz / sz.sum()
    worst = 0.0
    for k in range(len(sz)):
        idx = np.flatnonzero(lab == k)
        if len(idx) < 2:
            continue
        s = xyz[idx] if len(idx) <= 60 else xyz[idx[:60]]
        dd = np.array([np.sqrt(((s - s[i]) ** 2).sum(-1).mean(-1)) for i in range(len(s))])
        worst = max(worst, float(dd[np.triu_indices(len(s), 1)].max()))
    return dict(groups=len(sz), top1=float(p.max() * 100),
                inv_simpson=float(1 / np.sum(p ** 2)),
                singles=float((sz == 1).mean() * 100), rmsd_max=worst)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--control-draws", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    res = _M.receptor_coords(_M.key_residues(a.residues))
    rows, ladder = [], []
    for n, d in enumerate(sorted(rp.BLACKSMITH.glob("raw_cloud_*")), 1):
        ident = d.name[len("raw_cloud_"):]
        fs = sorted(d.glob("cloud_*.sdf"), key=os.path.getmtime)
        xyz, en, meta = load_with_energy(fs[-1])
        if xyz is None:
            continue
        rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
        w = pc.atom_weights(rmsf)
        tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
        order = np.argsort(en)
        for kf in KEEP:
            k = max(20, int(round(len(xyz) * kf)))
            sel = np.sort(order[:k])
            sub = xyz[sel]
            D = pc.pose_distances(pc.contact_tensor(sub, res), w)
            rows.append(dict(ident=ident, keep=kf, poses=len(sub), tol=tol,
                             e_cut=float(en[order[k - 1]]),
                             **stats(D, tol, sub)))
        # the saturation ladder, best-25% vs all, on the same molecule
        if len(xyz) >= 1500:
            for kf in (1.00, 0.25):
                k = int(round(len(xyz) * kf))
                pool = np.sort(order[:k])
                for m_ in (100, 200, 350, 500, 750, 1000):
                    if m_ > len(pool):
                        continue
                    idx = np.sort(np.random.default_rng(a.seed).choice(pool, m_, False))
                    s = xyz[idx]
                    D = pc.pose_distances(pc.contact_tensor(s, res), w)
                    ladder.append(dict(ident=ident, keep=kf, poses=m_,
                                       groups=int(pc.group(D, tol).max() + 1)))
        log.info("  %2d %s (%d poses)", n, ident, len(xyz))

    # ---- the size-matched control ---- #
    ctl = []
    rngc = np.random.default_rng(a.seed + 11)
    for n, dd in enumerate(sorted(rp.BLACKSMITH.glob("raw_cloud_*")), 1):
        ident = dd.name[len("raw_cloud_"):]
        fs = sorted(dd.glob("cloud_*.sdf"), key=os.path.getmtime)
        xyz, en, meta = load_with_energy(fs[-1])
        if xyz is None or len(xyz) < 400:
            continue
        rmsf = _M.predict_rmsf(meta[0], meta[1], 50, a.seed)
        w = pc.atom_weights(rmsf)
        tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
        for kf in [k for k in KEEP if k < 1.0]:
            k = max(20, int(round(len(xyz) * kf)))
            for r in range(a.control_draws):
                sel = np.sort(rngc.choice(len(xyz), k, replace=False))
                sub = xyz[sel]
                D = pc.pose_distances(pc.contact_tensor(sub, res), w)
                ctl.append(dict(ident=ident, keep=kf, draw=r, **stats(D, tol, sub)))
    C = pd.DataFrame(ctl)

    d = pd.DataFrame(rows)
    t = sout.Topic("blacksmith", "pose_generation_audit")
    if len(C):
        C.to_csv(t.write("energy_filtered_control", ".csv"), index=False)
    d.to_csv(t.write("energy_filtered", ".csv"), index=False)
    L = pd.DataFrame(ladder)
    if len(L):
        L.to_csv(t.write("energy_filtered_ladder", ".csv"), index=False)

    P = print
    P("\n" + "=" * 80)
    P("  DOES AN ENERGY FILTER CHANGE THE GROUPING CONCLUSIONS?")
    P("=" * 80)
    P(f"\n  {d.ident.nunique()} molecules with per-pose energies\n")
    P(f"    {'keep':>6}{'poses':>8}{'groups':>8}{'top-1':>9}{'eff # poses':>13}"
      f"{'singletons':>12}{'worst RMSD':>12}")
    for kf in KEEP:
        s = d[d.keep == kf]
        P(f"    {kf * 100:5.0f}%{s.poses.median():8.0f}{s.groups.median():8.0f}"
          f"{s.top1.median():8.1f}%{s.inv_simpson.median():13.0f}"
          f"{s.singles.median():11.0f}%{s.rmsd_max.max():11.2f}A")
    P("\n  'eff # poses' is the inverse Simpson index — the effective number of")
    P("  distinct poses, dominated by the big groups rather than the singleton tail.")

    P("\n  THE SIZE-MATCHED CONTROL — best N% against a RANDOM N% of the same")
    P("  cloud. Without this, a rising top-1 share is just a smaller sample.\n")
    if len(C):
        from scipy.stats import wilcoxon
        P(f"    {'keep':>6}{'top-1 best':>12}{'top-1 random':>14}{'ratio':>8}"
          f"{'eff# best':>11}{'eff# random':>13}")
        for kf in [k for k in KEEP if k < 1.0]:
            b = d[d.keep == kf]
            r = C[C.keep == kf]
            P(f"    {kf * 100:5.0f}%{b.top1.median():11.1f}%{r.top1.median():13.1f}%"
              f"{b.top1.median() / r.top1.median():8.2f}x"
              f"{b.inv_simpson.median():11.0f}{r.inv_simpson.median():13.0f}")
        for kf in [k for k in KEEP if k < 1.0]:
            bb = d[d.keep == kf].set_index("ident").top1
            rr = C[C.keep == kf].groupby("ident").top1.median()
            j = bb.index.intersection(rr.index)
            try:
                pv = wilcoxon(bb[j], rr[j]).pvalue
            except ValueError:
                pv = float("nan")
            P(f"    keep {kf * 100:3.0f}%: energy-selected is more concentrated in "
              f"{int((bb[j] > rr[j]).sum())} of {len(j)} molecules  (p = {pv:.2g})")

    if len(L):
        P("\n  DOES THE COUNT STILL CLIMB WITHIN THE GOOD-ENERGY POSES?\n")
        for ident, g in L.groupby("ident"):
            for kf in sorted(g.keep.unique()):
                s = g[g.keep == kf].groupby("poses").groups_ if False else \
                    g[g.keep == kf].groupby("poses")["groups"].mean()
                b, _ = _M.fit_power(s.index.values, s.values)
                P(f"    {ident}  keep {kf * 100:3.0f}%:  b = {b:+.3f}   "
                  f"({int(s.iloc[0])} groups at n={s.index[0]} -> "
                  f"{int(s.iloc[-1])} at n={s.index[-1]})")

    P("\n" + "=" * 80)
    P(f"  written to {t.dir}")
    P("=" * 80 + "\n")


if __name__ == "__main__":
    main()
