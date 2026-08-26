#!/usr/bin/env python3
"""
Purpose: does the number of contact-space groups taper off as the pose cloud deepens?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 6,000-pose deep cloud, the 5 independent replicate clouds, production clouds
Output: 00_outputs/blacksmith/contact_saturation/

exp/16 SHOWED THE GROUPS ARE TIGHT AND THERE ARE ~174 OF THEM PER MOLECULE. This
asks the only question that decides whether 174 is acceptable: is it a number, or
is it a point on a line that keeps climbing? If groups keep accumulating with
docking depth the count is an artefact of how long we docked and no shortlist
built on it means anything.

WHAT IS GUARANTEED AND THEREFORE NOT A FINDING. Complete linkage at a fixed
absolute tolerance cannot produce more groups than the covering number of the
occupied region at that scale, so the count is bounded a priori. D0091 is the
record of reporting exactly that kind of apparatus-entailed bound as a result,
and this docstring exists so it is not done twice. The empirical content here is
the RATE -- where on the curve 500 poses sits, and what the plateau is -- never
the fact that a plateau exists.

THE TOLERANCE IS HELD FIXED ALONG EACH LADDER. It is the molecule's own predicted
RMSF (exp/15) and does not depend on n, so a rung's count changes for one reason
only. A rule whose length scale moved with the sample -- the 5%-of-sample density
threshold behind D0088 -- is why mode counts never saturated before.

SUBSAMPLING IS CHECKED, NOT ASSUMED. AutoDock's LGA runs are independent, so an
n-subset of a 6,000-run cloud should be distributed like an n-run docking. Should
be. Panel B measures it: five INDEPENDENT dockings against subsamples of the deep
cloud at the same n. If independent dockings find more groups, every rung below
6,000 understates and the taper is an artefact of the pool being finite.
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
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("contact-saturation")


def _by_path(name: str, path: Path):
    """Import a module by FILE, never by name.

    exp/15 and exp/16 both expose `run_all`, so `import run_all` returns whichever
    sys.path entry won -- the plainest possible instance of selection by name over
    identity, in a repo whose whole failure catalogue is that shape. Each helper is
    loaded from its own file and bound to its own alias.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_RMSF = _by_path("exp15_rmsf", REPO / "exp" / "15_rmsf_predictor" / "run_all.py")
_CONTACT = _by_path("exp16_contact", REPO / "exp" / "16_contact_clustering" / "run_all.py")
predict_rmsf = _RMSF.predict_rmsf
key_residues = _CONTACT.key_residues
receptor_coords = _CONTACT.receptor_coords

LADDER = [100, 200, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000]


def draws_for(n: int) -> int:
    """Repeats at each rung. More where the scatter is largest, not uniformly."""
    return 5 if n <= 1000 else 3 if n <= 3000 else 2 if n < 6000 else 1


# --------------------------------------------------------------------------- #
#  fits
# --------------------------------------------------------------------------- #
def fit_power(n, y) -> tuple:
    """(b, rss) for y = a*n^b. b=1 is linear growth, b=0 is a flat plateau."""
    n, y = np.asarray(n, float), np.asarray(y, float)
    ok = (n > 0) & (y > 0)
    if ok.sum() < 3:
        return float("nan"), float("nan")
    A = np.vstack([np.ones(ok.sum()), np.log(n[ok])]).T
    coef, *_ = np.linalg.lstsq(A, np.log(y[ok]), rcond=None)
    pred = np.exp(A @ coef)
    return float(coef[1]), float(((y[ok] - pred) ** 2).sum())


def fit_saturating(n, y) -> tuple:
    """(G_inf, K, n95, rss) for the species-accumulation shape y = G*n/(n+K).

    Chosen over an exponential because it is the standard rarefaction curve and
    because 1/y = (K/G)(1/n) + 1/G is linear, so the fit has no starting guess to
    get wrong. n95 -- the depth reaching 95% of the plateau -- is 19K, and it is
    the number the campaign actually needs.
    """
    n, y = np.asarray(n, float), np.asarray(y, float)
    ok = (n > 0) & (y > 0)
    if ok.sum() < 3:
        return (float("nan"),) * 4
    A = np.vstack([np.ones(ok.sum()), 1.0 / n[ok]]).T
    coef, *_ = np.linalg.lstsq(A, 1.0 / y[ok], rcond=None)
    if coef[0] <= 0:                       # no finite plateau implied by the data
        return float("inf"), float("nan"), float("inf"), float("nan")
    G = 1.0 / coef[0]
    K = coef[1] * G
    pred = G * n[ok] / (n[ok] + K)
    return G, K, 19.0 * K, float(((y[ok] - pred) ** 2).sum())


# --------------------------------------------------------------------------- #
#  one clustering
# --------------------------------------------------------------------------- #
def cluster(xyz, res, w, tol) -> dict:
    """Group one pose set at a FIXED tolerance and describe the result."""
    T = pc.contact_tensor(xyz, res)
    D = pc.pose_distances(T, w)
    lab = pc.group(D, tol)
    worst = pc.within_group_max(D, lab)
    assert worst <= tol + 1e-6, f"linkage broke its guarantee: {worst} > {tol}"
    sizes = np.bincount(lab)
    return dict(groups=len(sizes), largest=int(sizes.max()),
                singletons=int((sizes == 1).sum()),
                non_singletons=int((sizes > 1).sum()),
                worst_contact_a=worst)


def cart_max(xyz, lab_sizes_idx) -> float:
    """Worst unaligned heavy-atom RMSD inside any group -- does quality hold up?"""
    out = 0.0
    for idx in lab_sizes_idx:
        if len(idx) < 2:
            continue
        s = xyz[idx] if len(idx) <= 60 else xyz[idx[:60]]
        d = np.array([np.sqrt(((s - s[i]) ** 2).sum(-1).mean(-1)) for i in range(len(s))])
        out = max(out, float(d[np.triu_indices(len(s), 1)].max()))
    return out


def cluster_full(xyz, res, w, tol) -> dict:
    """`cluster`, plus the Cartesian-width check, on the same labelling."""
    T = pc.contact_tensor(xyz, res)
    D = pc.pose_distances(T, w)
    lab = pc.group(D, tol)
    worst = pc.within_group_max(D, lab)
    assert worst <= tol + 1e-6, f"linkage broke its guarantee: {worst} > {tol}"
    sizes = np.bincount(lab)
    idxs = [np.flatnonzero(lab == k) for k in range(len(sizes))]
    return dict(groups=len(sizes), largest=int(sizes.max()),
                singletons=int((sizes == 1).sum()),
                non_singletons=int((sizes > 1).sum()),
                worst_contact_a=worst, cart_max=cart_max(xyz, idxs))


# --------------------------------------------------------------------------- #
def load_sdf(path: Path):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    ms = [m for m in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
          if m is not None]
    if len(ms) < 30:
        return None, None
    hv = [a.GetIdx() for a in ms[0].GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.array([m.GetConformer().GetPositions()[hv] for m in ms])
    return xyz, (ms[0], hv)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--candidate", default="t4_716800c125a7")
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--conformers", type=int, default=50)
    ap.add_argument("--shallow-molecules", type=int, default=12)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    names = key_residues(a.residues)
    res = receptor_coords(names)
    log.info("landmarks (%d): %s", len(names), ", ".join(names))

    rng = np.random.default_rng(a.seed)
    t = sout.Topic("blacksmith", "contact_saturation")

    # ---------------- Panel A: the deep ladder ---------------- #
    deep = rp.BLACKSMITH / f"deep_cloud_{a.candidate}" / "cloud_1.sdf"
    xyz, meta = load_sdf(deep)
    if xyz is None:
        raise SystemExit(f"no deep cloud at {deep}")
    tmpl, hv = meta
    rmsf = predict_rmsf(tmpl, hv, a.conformers, a.seed)
    w = pc.atom_weights(rmsf)
    tol = float(np.median(rmsf) / pc.RMSF_CALIBRATION)
    log.info("deep cloud: %d poses x %d heavy atoms, tolerance %.3f A (FIXED)",
             len(xyz), xyz.shape[1], tol)

    rows = []
    for k in [x for x in LADDER if x <= len(xyz)]:
        for d in range(draws_for(k)):
            idx = (np.arange(len(xyz)) if k == len(xyz)
                   else rng.choice(len(xyz), size=k, replace=False))
            r = cluster_full(xyz[idx], res, w, tol)
            r.update(poses=k, draw=d, source="subsample")
            rows.append(r)
        m = np.mean([r["groups"] for r in rows if r["poses"] == k])
        log.info("  n=%5d -> %6.1f groups (mean of %d)", k, m, draws_for(k))
    A = pd.DataFrame(rows)
    A.to_csv(t.write("deep_ladder", ".csv"), index=False)

    agg = A.groupby("poses")["groups"].mean()
    b, rss_p = fit_power(agg.index.values, agg.values)
    G, K, n95, rss_s = fit_saturating(agg.index.values, agg.values)

    # ---------------- Panel B: is subsampling honest? ---------------- #
    brows = []
    for r in range(1, 6):
        p = rp.BLACKSMITH / f"election_{a.candidate}_r{r}_allposes" / f"{a.candidate}.sdf"
        if not p.is_file():
            continue
        rx, _ = load_sdf(p)
        if rx is None or rx.shape[1] != xyz.shape[1]:
            log.warning("  r%d: skipped (atom count %s)", r,
                        None if rx is None else rx.shape[1])
            continue
        c = cluster(rx, res, w, tol)
        brows.append(dict(replicate=r, poses=len(rx), groups=c["groups"],
                          source="independent"))
        # the matched subsample of the deep cloud, same n
        for d in range(3):
            sub = rng.choice(len(xyz), size=len(rx), replace=False)
            c2 = cluster(xyz[sub], res, w, tol)
            brows.append(dict(replicate=r, poses=len(rx), groups=c2["groups"],
                              source="subsample"))
    B = pd.DataFrame(brows)
    if len(B):
        B.to_csv(t.write("independence_check", ".csv"), index=False)

    # ---------------- Panel C: is the shape molecule-specific? ---------------- #
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{rp.topic()}_conditional_eb_*.csv")),
                key=os.path.getmtime)
    mols = [m for m in pd.read_csv(fs[-1]).parent_ident.dropna().unique()
            if (rp.allposes_dir() / f"{m}.sdf").is_file()]
    mols = list(np.random.default_rng(a.seed).choice(
        mols, size=min(a.shallow_molecules, len(mols)), replace=False))
    crows = []
    for i, ident in enumerate(mols, 1):
        cx, cm = load_sdf(rp.allposes_dir() / f"{ident}.sdf")
        if cx is None:
            continue
        try:
            rm = predict_rmsf(cm[0], cm[1], a.conformers, a.seed)
        except Exception as exc:                            # noqa: BLE001
            log.warning("  %s: no RMSF (%s)", ident, str(exc)[:50])
            continue
        ww, tt = pc.atom_weights(rm), float(np.median(rm) / pc.RMSF_CALIBRATION)
        pts = []
        for k in [x for x in LADDER if x <= len(cx)] + [len(cx)]:
            if k in pts:
                continue
            pts.append(k)
            for d in range(2 if k < len(cx) else 1):
                sub = (np.arange(len(cx)) if k == len(cx)
                       else rng.choice(len(cx), size=k, replace=False))
                c = cluster(cx[sub], res, ww, tt)
                crows.append(dict(ident=ident, poses=k, draw=d,
                                  groups=c["groups"], tol_a=tt))
        log.info("  %2d/%d %s: ladder to %d poses", i, len(mols), ident, len(cx))
    C = pd.DataFrame(crows)
    if len(C):
        C.to_csv(t.write("shallow_ladders", ".csv"), index=False)

    # ---------------- report ---------------- #
    print("\n" + "=" * 78)
    print("  DOES THE GROUP COUNT TAPER?  contact-space grouping, fixed tolerance")
    print("=" * 78)
    print(f"\n  A. DEEP LADDER — {a.candidate}, tolerance {tol:.2f} A held fixed\n")
    print(f"    {'poses':>7} {'groups':>9} {'per 1k new':>11} {'largest':>8} "
          f"{'singl.':>7} {'worst RMSD':>11}")
    prev = None
    for k, g in agg.items():
        sub = A[A.poses == k]
        rate = "" if prev is None else f"{(g - prev[1]) / (k - prev[0]) * 1000:9.1f}"
        print(f"    {k:7,} {g:9.1f} {rate:>11} {sub.largest.mean():8.1f} "
              f"{sub.singletons.mean() / g * 100:6.0f}% {sub.cart_max.max():10.2f} A")
        prev = (k, g)
    print(f"\n    power-law exponent b (groups ~ n^b): {b:+.3f}   "
          f"(1.0 = linear, 0.0 = flat)")
    if np.isfinite(G):
        print(f"    saturating fit: plateau {G:,.0f} groups, "
              f"95% of it at n = {n95:,.0f} poses")
        print(f"    fit residual  saturating {rss_s:,.0f}  vs  power law {rss_p:,.0f}"
              f"   -> {'saturating' if rss_s < rss_p else 'power law'} fits better")
        print(f"    at 500 poses we are at {agg.get(500, float('nan')) / G * 100:.0f}%"
              f" of the plateau; at 6,000, {agg.iloc[-1] / G * 100:.0f}%")
    else:
        print("    saturating fit: NO finite plateau implied — the count is still "
              "climbing linearly")

    if len(B):
        print(f"\n  B. IS SUBSAMPLING HONEST?  independent dockings vs subsamples "
              f"at the same n\n")
        for n_, sub in B.groupby("poses"):
            ind = sub[sub.source == "independent"].groups
            ss = sub[sub.source == "subsample"].groups
            print(f"    n={n_:5,}  independent {ind.mean():6.1f} "
                  f"(n={len(ind)})   subsample {ss.mean():6.1f} (n={len(ss)})   "
                  f"ratio {ind.mean() / ss.mean():.3f}")
        ind = B[B.source == "independent"].groups.mean()
        ss = B[B.source == "subsample"].groups.mean()
        r = ind / ss
        # TWO-SIDED, because the first version was not. It read
        # `"UNDERSTATES" if r > 1.05 else "fair stand-in"`, so a ratio of 0.461 --
        # a factor of 2.2 in the other direction -- printed "fair stand-in". A
        # guard that can only fail one way cannot fail the way this one did.
        verdict = ("Subsampling UNDERSTATES; the ladder is biased toward taper."
                   if r > 1.05 else
                   "Subsampling OVERSTATES; the two clouds are not the same "
                   "population (see D0093 -- `allposes` is DBSCAN-cleaned)."
                   if r < 0.95 else
                   "Subsampling is a fair stand-in for an independent docking.")
        print(f"\n    overall independent/subsample = {r:.3f}. {verdict}")

    if len(C):
        print(f"\n  C. IS THE SHAPE MOLECULE-SPECIFIC?  {C.ident.nunique()} molecules, "
              f"ladders to their own cloud depth\n")
        bs = []
        for ident, sub in C.groupby("ident"):
            m = sub.groupby("poses")["groups"].mean()
            bb, _ = fit_power(m.index.values, m.values)
            bs.append(bb)
        bs = np.array([x for x in bs if np.isfinite(x)])
        print(f"    exponent b per molecule: median {np.median(bs):+.3f}   "
              f"range {bs.min():+.3f} to {bs.max():+.3f}")
        print(f"    molecules with b > 0.9 (essentially linear): "
              f"{(bs > 0.9).sum()} of {len(bs)}")

    print("\n" + "=" * 78)
    print(f"  written to {t.dir}")
    print("=" * 78 + "\n")


if __name__ == "__main__":
    main()
