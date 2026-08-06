"""
Purpose: pick the stratified cohort that tests WHICH ranking metric selects for physical stability.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: the consensus screen (all 5,769) + the candidate frames
Output: 00_outputs/blacksmith/elevation_cohort/elevation_cohort_<N>.csv

THE QUESTION, in @tt8804's words: take a few from each group — high-enrichment /
high-consensus, low-enrichment / high-consensus, high-enrichment / low-consensus
— and see which poses survive. That directly tests which metric is selecting for
something physical.

THREE CONFOUNDS FORCED THE DESIGN AWAY FROM A CLEAN 2x2, and each is measured
rather than assumed:

1. **Consensus is partly rigidity.** Spearman -0.259 against rotatable bonds
   (p = 4e-89). Rigid molecules have fewer ways to sit AND are trivially more
   stable under dynamics, so an unmatched comparison would find "high consensus
   is stable" for a reason that has nothing to do with the pocket. Groups are
   therefore **matched on rotatable-bond count**, not merely balanced in size.

2. **Enrichment and warhead class are confounded.** The high-enrichment cells are
   BDHI-dominated (10 molecules, 8 BDHI) while the low-enrichment cells are
   acrylamide-dominated (340, of which 240 acrylamide). "Does enrichment predict
   stability" and "does BDHI predict stability" cannot be separated across the
   whole set — only WITHIN BDHI. So the core comparison is BDHI-only, and that
   is a limit on what the result can claim, stated here rather than discovered
   later.

3. **High enrichment is nearly a SUBSET of high consensus.** The
   high-enrichment/low-consensus cell holds 5 molecules. Enrichment is not an
   independent axis; it is a stricter filter sitting inside consensus. So the
   informative contrast is A vs B — does the extra strictness buy anything — and
   not a four-cell factorial that the data cannot populate.

THE ANCHOR IS WHAT MAKES THIS INTERPRETABLE. Crystallographic positives are the
only molecules we KNOW react with Cys113. Without them, "group A is more stable
than group D" is a statement about two arbitrary groups. With them, the question
becomes "which group behaves like a molecule that actually reacts" — which is
the question worth answering.

TWO-TIER READOUT, and tier 1 is free:

    tier 1  warhead displacement across 300 ps of UNRESTRAINED equilibration.
            `gromacs_explicit` applies no position restraints during NVT/NPT, so
            every run already measures whether the docked pose survives plain
            dynamics before any bias is applied. That is a physical stability
            signal we were discarding.
    tier 2  BPMD: how hard the warhead is to push out of the near-attack window.

Reporting tier 1 also keeps tier 2 honest — BPMD starts from the POST-
equilibration pose, so a molecule that drifted 3 A during equilibration is not
having its docked pose tested at all.
"""

from __future__ import annotations

import argparse
import glob
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout               # noqa: E402

log = logging.getLogger("elevation-cohort")
OUT = sout.Topic("blacksmith", "elevation_cohort")
DATA = Path("/data/lab_vm/append_only/inhibition")

# Thresholds, fixed before selection. HI_CONS is the level the >5.7 band showed
# (median 0.933); LO_CONS is comfortably below the pool median (~0.29).
HI_CONS, LO_CONS = 0.90, 0.50
HI_ENR, LO_ENR = 5.0, 3.0
PER_GROUP = 8


def load() -> pd.DataFrame:
    fs = glob.glob(str(OUT.dir.parent / "nac_consensus" / "*.csv"))
    c = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    c = c.drop_duplicates("ident")
    c = c[(c.status == "ok") & c.consensus.notna()].copy()

    smi, extra = {}, {}
    for sub, stem in (("03_t3_reinvent", "D3"), ("04_t4_combinatorial", "D4")):
        g = list((DATA / sub).glob(f"{stem}_*.parquet"))
        d = pd.read_parquet(max(g, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1))))
        smi.update(dict(zip(d.candidate_id.astype(str), d.canonical_smiles)))
        for col in ("QED", "nac_pose_path"):
            if col in d.columns:
                extra.setdefault(col, {}).update(dict(zip(d.candidate_id.astype(str), d[col])))
    c["smiles"] = c.ident.map(smi)
    for col, m in extra.items():
        c[col] = c.ident.map(m)
    c = c[c.smiles.notna()].copy()

    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as rdd
    RDLogger.DisableLog("rdApp.*")
    def rotb(s):
        m = Chem.MolFromSmiles(s)
        return np.nan if m is None else rdd.CalcNumRotatableBonds(m)
    c["rotb"] = c.smiles.map(rotb)
    return c.dropna(subset=["rotb"])


def match_rotb(pool: pd.DataFrame, target: pd.Series, n: int,
               rng: np.random.Generator) -> pd.DataFrame:
    """Sample `n` from `pool` whose rotatable-bond distribution matches `target`.

    Nearest-neighbour on rotb, sampled without replacement. Matching rather than
    randomising because the confound is known and directional: consensus already
    correlates with rigidity at rho = -0.259, so a random draw would reproduce it.
    """
    avail, picked = pool.copy(), []
    for t in sorted(target.values)[:n]:
        if avail.empty:
            break
        d = (avail.rotb - t).abs()
        best = d[d == d.min()].index
        i = rng.choice(best)
        picked.append(i)
        avail = avail.drop(i)
    return pool.loc[picked]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--per-group", type=int, default=PER_GROUP)
    ap.add_argument("--seed", type=int, default=0xE1E7A7)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    rng = np.random.default_rng(args.seed)

    c = load()
    log.info("%d molecules with consensus + properties", len(c))

    bdhi = c[c.warhead_class.str.startswith("bdhi")]
    A = bdhi[(bdhi.enrichment >= HI_ENR) & (bdhi.consensus >= HI_CONS)]
    Bp = bdhi[(bdhi.enrichment < LO_ENR) & (bdhi.consensus >= HI_CONS)]
    Dp = bdhi[(bdhi.enrichment < LO_ENR) & (bdhi.consensus < LO_CONS)]

    n = min(args.per_group, len(A))
    A = A.nlargest(n, "consensus")
    B = match_rotb(Bp, A.rotb, n, rng)
    D = match_rotb(Dp, A.rotb, n, rng)

    # The validated warhead class, high consensus. Cannot be rotb-matched to the
    # BDHI groups and is not meant to be — it answers a different question:
    # does consensus-selection inside the ONE class with AUC 0.908 find molecules
    # that behave like known binders?
    V = c[(c.warhead_class == "chloroacetamide") & (c.consensus >= HI_CONS)]
    V = V.nlargest(min(args.per_group, len(V)), "consensus")

    groups = {"A_hiEnr_hiCons_bdhi": A, "B_loEnr_hiCons_bdhi": B,
              "D_loEnr_loCons_bdhi": D, "V_hiCons_chloroacetamide": V}

    rows = []
    print(f"\n=== elevation cohort ===")
    print(f"  {'group':<28}{'n':>3}{'enrich':>9}{'consens':>9}{'rotb':>7}{'QED':>7}")
    for g, df in groups.items():
        if df.empty:
            log.warning("%s is empty", g); continue
        print(f"  {g:<28}{len(df):>3}{df.enrichment.median():>9.2f}"
              f"{df.consensus.median():>9.3f}{df.rotb.median():>7.1f}"
              f"{df.QED.median() if 'QED' in df else float('nan'):>7.3f}")
        for r in df.itertuples():
            rows.append({"group": g, "ident": r.ident, "smiles": r.smiles,
                         "warhead_class": r.warhead_class,
                         "enrichment": r.enrichment, "consensus": r.consensus,
                         "consensus_modes": r.consensus_modes, "rotb": r.rotb,
                         "QED": getattr(r, "QED", np.nan),
                         "nac_pose_path": getattr(r, "nac_pose_path", None)})

    out = pd.DataFrame(rows)
    dest = OUT.write("elevation_cohort", ".csv")
    out.to_csv(dest, index=False)
    print(f"\n  {len(out)} molecules -> {dest}")
    print(f"  rotb matched across the three BDHI groups: "
          f"{[float(groups[g].rotb.median()) for g in list(groups)[:3] if not groups[g].empty]}")
    print("\n  Crystallographic positives are added by the RUNNER as the anchor —")
    print("  they are the only molecules known to react, and without them a")
    print("  between-group difference has nothing to be measured against.")


if __name__ == "__main__":
    main()
