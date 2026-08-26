#!/usr/bin/env python3
"""
Purpose: does any molecule's pose cloud actually converge, at any scale?
Author: Timothy Wu (with Claude Code)
Date: 2026-08-26
Input: the 20 RAW clouds from scripts/persist_raw_clouds.py
Output: 00_outputs/blacksmith/cloud_concentration/

@tt8804: "are there any molecules that have truley consensus poses? or very
little? whats the spread on pose number". At the shipped 0.6 A tolerance the
answer is no for all 20 -- the largest group holds 1.4-7.2% of its cloud. But
that is one tolerance, and "is there a consensus pose" is a question about the
molecule, not about a cut. So it is asked across scales.

CONCENTRATION, NOT COUNT. A group count answers "how finely did we cut", which is
a property of the tolerance (D0092). What a consensus pose would look like is a
cloud whose MASS is concentrated: one group holding a large share. So the
statistics here are share-based -- top-1 share, and the inverse Simpson index,
which is the EFFECTIVE number of distinct poses and is dominated by the big
groups rather than by the singleton tail. A cloud of 500 poses with one true pose
would score near 1; a uniformly scattered cloud scores near the group count.

THE NULL IS COMPUTED, NOT ASSUMED. A cloud of poses scattered at random still
produces a largest group, and at a tight tolerance that group is small for
reasons that have nothing to do with binding. So each molecule is compared
against ITS OWN shuffle: the same poses with contact vectors resampled per
coordinate, which destroys pose identity while preserving each coordinate's
marginal distribution. Without that, "the top group holds 3%" is unreadable.
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
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                  # noqa: E402
from shared import pose_contacts as pc              # noqa: E402
from shared import run_paths as rp                  # noqa: E402

log = logging.getLogger("cloud-concentration")
TOLS = [0.6, 1.0, 1.5, 2.0, 2.5, 3.0]


def _by_path(name: str, rel: str):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_M = _by_path("e17", "exp/17_contact_saturation/run_all.py")


def concentration(lab: np.ndarray) -> dict:
    sz = np.bincount(lab).astype(float)
    p = sz / sz.sum()
    return dict(groups=len(sz), largest=int(sz.max()),
                top1=float(p.max() * 100),
                top5=float(np.sort(p)[::-1][:5].sum() * 100),
                inv_simpson=float(1.0 / np.sum(p ** 2)),
                singles=float((sz == 1).mean() * 100))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--residues", type=int, default=15)
    ap.add_argument("--null-draws", type=int, default=3)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr
    from rdkit.Chem import rdMolDescriptors

    res = _M.receptor_coords(_M.key_residues(a.residues))
    dirs = sorted(rp.BLACKSMITH.glob("raw_cloud_*"))
    if not dirs:
        raise SystemExit("no raw clouds; run scripts/persist_raw_clouds.py first")

    rng = np.random.default_rng(a.seed)
    rows, nulls = [], []
    for n, d in enumerate(dirs, 1):
        ident = d.name[len("raw_cloud_"):]
        fs = sorted(d.glob("cloud_*.sdf"), key=os.path.getmtime)
        xyz, meta = _M.load_sdf(fs[-1])
        if xyz is None:
            continue
        tmpl, hv = meta
        rmsf = _M.predict_rmsf(tmpl, hv, 50, a.seed)
        w = pc.atom_weights(rmsf)
        T = pc.contact_tensor(xyz, res)
        D = pc.pose_distances(T, w)
        rotb = rdMolDescriptors.CalcNumRotatableBonds(tmpl)
        for t_ in TOLS:
            r = concentration(pc.group(D, t_))
            rows.append(dict(ident=ident, tol=t_, poses=len(xyz),
                             rotb=rotb, heavy=xyz.shape[1], **r))
        # the null: same marginals, pose identity destroyed
        for k in range(a.null_draws):
            Ts = T.copy()
            for c in range(Ts.shape[1]):
                for j in range(Ts.shape[2]):
                    Ts[:, c, j] = Ts[rng.permutation(len(Ts)), c, j]
            Dn = pc.pose_distances(Ts, w)
            for t_ in TOLS:
                nulls.append(dict(ident=ident, tol=t_, draw=k,
                                  **concentration(pc.group(Dn, t_))))
        log.info("  %2d/%d %s (%d rot. bonds)", n, len(dirs), ident, rotb)

    d = pd.DataFrame(rows)
    nl = pd.DataFrame(nulls)
    t = sout.Topic("blacksmith", "cloud_concentration")
    d.to_csv(t.write("per_molecule_per_tolerance", ".csv"), index=False)
    nl.to_csv(t.write("shuffled_null", ".csv"), index=False)

    P = print
    P("\n" + "=" * 82)
    P("  DOES ANY MOLECULE'S POSE CLOUD CONVERGE?   "
      f"{d.ident.nunique()} molecules, 500 poses each")
    P("=" * 82)
    P("\n  SHARE OF THE CLOUD IN THE SINGLE LARGEST GROUP\n")
    P(f"    {'tolerance':>10} {'median':>8} {'best mol':>10} {'worst':>8}"
      f"   {'shuffled null':>14}   {'effective # poses':>18}")
    for t_ in TOLS:
        s = d[d.tol == t_]
        ns = nl[nl.tol == t_]
        P(f"    {t_:9.1f}A {s.top1.median():7.1f}% {s.top1.max():9.1f}% "
          f"{s.top1.min():7.1f}%   {ns.top1.median():13.1f}%   "
          f"{s.inv_simpson.median():17.0f}")
    P("\n  'effective # poses' is the inverse Simpson index: 1 would mean one pose,")
    P("  500 would mean every pose distinct. It counts the big groups, not the tail.")

    P("\n  THE MOST AND LEAST CONVERGED MOLECULES (at 2.0 A)\n")
    s = d[d.tol == 2.0].sort_values("top1", ascending=False)
    P(f"    {'molecule':<20}{'top-1':>8}{'top-5':>8}{'eff #':>8}{'rot.bonds':>11}")
    for _, r in pd.concat([s.head(3), s.tail(3)]).iterrows():
        P(f"    {r.ident:<20}{r.top1:7.1f}%{r.top5:7.1f}%{r.inv_simpson:8.0f}"
          f"{int(r.rotb):11d}")

    P("\n  DOES FLEXIBILITY EXPLAIN THE SPREAD?\n")
    for t_ in (1.0, 2.0, 3.0):
        s = d[d.tol == t_]
        rho_t = spearmanr(s.rotb, s.top1)[0]
        rho_e = spearmanr(s.rotb, s.inv_simpson)[0]
        P(f"    tol {t_:.1f} A:  rho(rotatable bonds, top-1 share) = {rho_t:+.3f}   "
          f"rho(rotatable bonds, effective # poses) = {rho_e:+.3f}")

    P("\n  AGAINST THE NULL — is the cloud more concentrated than chance?\n")
    for t_ in TOLS:
        s = d[d.tol == t_].set_index("ident").top1
        ns = nl[nl.tol == t_].groupby("ident").top1.median()
        j = s.index.intersection(ns.index)
        # THE NULL DIES AT LOOSE TOLERANCES AND MUST SAY SO. Shuffling each
        # contact coordinate independently pulls every pose toward the marginal
        # mean, so above ~2 A the shuffled cloud collapses into ONE group and
        # scores 100%. The resulting "real / shuffled = 0.4x" is not evidence
        # that real poses are less concentrated than chance -- it is the null
        # being degenerate. Reporting the ratio there would be a comparison
        # against an apparatus artefact, which is D0091's mistake.
        if float(ns[j].median()) > 50.0:
            P(f"    tol {t_:.1f} A: NULL DEGENERATE — the shuffled cloud collapses "
              f"to one group ({ns[j].median():.0f}% top-1); no comparison possible")
            continue
        ratio = (s[j] / ns[j])
        P(f"    tol {t_:.1f} A: real / shuffled top-1 share = "
          f"{ratio.median():.2f}x   ({int((ratio > 1.5).sum())} of {len(j)} "
          f"molecules above 1.5x)")

    P("\n" + "=" * 82)
    P(f"  written to {t.dir}")
    P("=" * 82 + "\n")


if __name__ == "__main__":
    main()
