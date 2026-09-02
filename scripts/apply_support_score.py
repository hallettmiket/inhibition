#!/usr/bin/env python3
"""
Purpose: compute the per-mode rank score (engagement x support) for a whole run.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-29
Input: the run's mode representatives + rank_v2 engagement table
Output: 00_outputs/blacksmith/rank_score_<topic>/rank_score_<N>.csv

@tt8804: "upgrade the ranking to consider more residues than cys113 but in a
conservative manner ... only adding basic derivative/supporting rules that would
be essential for cys113 engagement".

WHAT IS SCORED, AND WHY SO LITTLE. Six residues sit inside the 4.2 A near-attack
window on the prepared 3IKD. Only SER114 and SER115 earn credit here -- see
`engagement_rank.SUPPORT_ATOMS` for why His59, the leucines and the Arg loop are
all excluded. Each of those would express a preference about what KIND of
molecule wins, which is the influence being ruled out.

THE SCORE CANNOT PROMOTE A POSE THAT CANNOT REACT. `rank_score = anchor_quality
* support`, support in [1.0, 1.15]. A pose with no attack geometry scores zero
however well supported, and an unsupported pose keeps its geometry score exactly.

MEASURED ON THE REPRESENTATIVE, which is the pose that would actually be
simulated -- not on the mode's mean. D0098: the representative's geometry
predicts the MD outcome at rho = +0.652 where every mode-level aggregate managed
~0.11, because a shipped mode spans 0.776 of the anchor scale.
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

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import engagement_rank as er           # noqa: E402
from shared import outputs as sout                 # noqa: E402
from shared import run_paths as rp                 # noqa: E402

log = logging.getLogger("support-score")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default=None)
    ap.add_argument("--max-support", type=float, default=er.MAX_SUPPORT)
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    topic = a.topic or rp.topic()
    fs = sorted(glob.glob(str(rp.BLACKSMITH /
                              f"rank_v2/rank_v2_T4_{topic}_engagement_*.csv")),
                key=os.path.getmtime)
    if not fs:
        raise SystemExit(f"no engagement ranking for {topic!r}")
    rk = pd.read_csv(fs[-1])
    sup = er.receptor_support_atoms()
    log.info("%s: %d modes; support atoms %s", topic, len(rk),
             [f"{rn}{ri}:{an}" for _, ri, rn, an in er.SUPPORT_ATOMS])

    rows, missing = [], 0
    poses = rp.BLACKSMITH / f"{topic}_poses"
    for n, f in enumerate(sorted(poses.glob("*.sdf")), 1):
        ident = f.stem
        for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=False):
            if m is None or not m.HasProp("mode"):
                missing += 1
                continue
            hv = [i for i, at in enumerate(m.GetAtoms()) if at.GetAtomicNum() > 1]
            P = m.GetConformer().GetPositions()[hv]
            els = [m.GetAtomWithIdx(i).GetSymbol() for i in hv]
            rows.append(dict(parent_ident=ident, mode=int(m.GetProp("mode")),
                             support=er.support_factor(P, els, sup, a.max_support)))
        if n % 400 == 0:
            log.info("  %d/%d molecules", n, len(list(poses.glob('*.sdf'))))
    s = pd.DataFrame(rows)
    log.info("support computed for %d modes (%d unreadable)", len(s), missing)

    d = rk.merge(s, on=["parent_ident", "mode"], how="left")
    n_nosup = int(d.support.isna().sum())
    # A MODE WITH NO REPRESENTATIVE GETS support = 1.0, NOT a dropped row: it
    # keeps its geometry score untouched, which is the conservative direction.
    d["support"] = d.support.fillna(1.0)
    d["rank_score"] = [er.rank_score(x, y) for x, y in zip(d.engagement, d.support)]

    t = sout.Topic("blacksmith", f"rank_score_{topic}")
    out = t.write("rank_score", ".csv")
    d.to_csv(out, index=False)

    sup_any = (d.support > 1.0)
    P = print
    P("\n" + "=" * 78)
    P("  PER-MODE RANK SCORE = engagement x support")
    P("=" * 78)
    P(f"\n  {len(d):,} modes · {d.parent_ident.nunique():,} molecules"
      f" · {n_nosup:,} without a representative (support = 1.0)\n")
    P(f"  supported by at least one serine : {int(sup_any.sum()):,} "
      f"({sup_any.mean()*100:.1f}%)")
    P(f"  supported by both                : {int((d.support >= 1 + a.max_support - 1e-9).sum()):,}")
    P(f"  support factor range             : {d.support.min():.3f} - {d.support.max():.3f}")
    P(f"\n  engagement  median {d.engagement.median():.4f}  max {d.engagement.max():.4f}")
    P(f"  rank_score  median {d.rank_score.median():.4f}  max {d.rank_score.max():.4f}")
    from scipy.stats import spearmanr
    ok = d.engagement.notna() & d.rank_score.notna()
    P(f"\n  rho(engagement, rank_score) = {spearmanr(d.engagement[ok], d.rank_score[ok])[0]:+.4f}")
    moved = (d.sort_values('rank_score', ascending=False).head(450).parent_ident.nunique())
    base = (d.sort_values('engagement', ascending=False).head(450).parent_ident.nunique())
    top_e = set(d.sort_values("engagement", ascending=False).head(450).index)
    top_r = set(d.sort_values("rank_score", ascending=False).head(450).index)
    P(f"  top-450 overlap with the geometry-only ordering: {len(top_e & top_r)}/450")
    P(f"  molecules represented there: {base} -> {moved}")
    P("\n" + "=" * 78)
    P(f"  written to {out}")
    P("=" * 78 + "\n")


if __name__ == "__main__":
    main()
