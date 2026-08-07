"""
Purpose: the 2.1.0 ranking — consensus filtered WITHIN warhead class, then ranked on anchoring geometry.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: the consensus screen (all 5,765 measured) + candidate frames for properties
Output: 00_outputs/blacksmith/rank_2_1/rank_2_1_<TIER>_<N>.csv, one ranked list per warhead class

@tt8804's framework (#18): consensus as a FILTER, then rank the surviving poses
on the anchoring score, then selection, then the elevation suite.

WHAT CHANGED FROM 2.0.0, AND WHY

1. THE FILTER IS PER WARHEAD CLASS, NOT LIBRARY-WIDE. D0073 measured a single
   0.90 bar promoting rigid chemistry and demoting flexible chemistry: pass rates
   ran BDHI 16.6%, SNAr 13.9%, Michael 6.3%, chloroacetamide 2.9%, while pass rate
   is monotone in rotatable-bond count (19.5% at 0-2 down to 2.9% at 7+). Pose
   agreement is partly a statement about how many ways a molecule CAN sit, so
   comparing a fused ring system's consensus to an acyclic chain's compares their
   ring counts. Taking a fixed FRACTION within each class removes the
   cross-class part of that, because every molecule then competes only against
   molecules with the same warhead.

   It does NOT remove the within-class part. Rigidity still varies inside a class
   and still helps. That is measured and reported per class rather than claimed
   to be solved.

2. NO VALIDATION GATE. Earlier drafts restricted the shortlist to warhead classes
   with crystallographic positives. @tt8804 ruled that out: 15 depositions across
   three classes is not enough to decide which chemistry to pursue, and the
   chloroacetamide series those positives mostly come from is not of interest.
   Classes are therefore ranked on their own terms, and no class is promoted or
   demoted for having, or lacking, an external validation set.

   (The crystallographic anchor keeps its OTHER job: calibrating the elevation
   stability assay, where 8 molecules make one comparison readable. That is not
   the same as using them to choose chemistry.)

3. NO SINGLE CROSS-CLASS TOP-N. The output is one ranked list per class. A merged
   "top 50 overall" would re-import exactly the bias point 1 removes, and D0073's
   phrasing stands: a top-N-overall list built on consensus is a rigidity ranking
   wearing a geometry label.

WHAT THIS RANKING IS NOT YET

The rank-within-filter component is `enrichment`, and it carries two open
problems that are stated here rather than discovered later:

  - D0071 showed enrichment does not predict pose stability on a pre-registered
    cohort. It has not been shown to predict anything physical.
  - `viable_fraction` is the JOINT rate P(distance AND angle) while
    `isotropic_null` is purely orientational P(angle), so the score conflates a
    distance hit rate that the reactive potential deliberately biases (D0064)
    with the orientation term it was meant to measure. The conditional repair is
    designed in `docs/ranking_2.1.0_design.md` §3.1 and is NOT applied here,
    because it needs per-pose distances and angles that were never persisted.

So this is the framework's SHAPE with 2.0.0's score still inside it. The shape is
the part that is defensible today; the score is the part under repair.
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

from shared import outputs as sout               # noqa: E402

log = logging.getLogger("rank-2-1")
OUT = sout.Topic("blacksmith", "rank_2_1")
DATA = Path("/data/lab_vm/append_only/inhibition")
CONS = DATA / "00_outputs/blacksmith/nac_consensus"

#: Fraction of each warhead class that survives the consensus filter. A quota
#: rather than an absolute bar, because an absolute bar is what let rigidity
#: decide which chemistry survived (D0073). Provisional: @tt8804 asked for tests
#: to set this quantitatively, and those have not been run.
DEFAULT_QUOTA = 0.20

#: Below this, "agreement among the top poses" is not a meaningful statement --
#: the poses disagree. A class whose quota reaches down here is telling you the
#: class has no well-determined poses, not that its 20th percentile is good.
CONSENSUS_FLOOR = 0.50


def load() -> pd.DataFrame:
    fs = sorted(glob.glob(str(CONS / "*.csv")))
    if not fs:
        raise SystemExit(f"no consensus shards under {CONS}")
    c = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates("ident")
    c = c[(c.status == "ok") & c.consensus.notna()].copy()

    smi, extra = {}, {}
    for sub, stem in (("03_t3_reinvent", "D3"), ("04_t4_combinatorial", "D4")):
        g = list((DATA / sub).glob(f"{stem}_*.parquet"))
        if not g:
            continue
        d = pd.read_parquet(max(g, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p.name).group(1))))
        smi.update(dict(zip(d.candidate_id.astype(str), d.canonical_smiles)))
        for col in ("QED", "nac_pose_path"):
            if col in d.columns:
                extra.setdefault(col, {}).update(dict(zip(d.candidate_id.astype(str), d[col])))
    c["smiles"] = c.ident.map(smi)
    for col, m in extra.items():
        c[col] = c.ident.map(m)

    c["tier"] = c.ident.str.split("_").str[0].str.upper()
    return c


def add_rotb(df: pd.DataFrame) -> pd.DataFrame:
    """Rotatable-bond count, kept as a REPORTED covariate rather than a correction.

    D0073 measured consensus correlating with rigidity at rho = -0.259. Silently
    residualising consensus on rotb would hide how much of each class's ranking
    is flexibility; showing the correlation per class lets a reader see it.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors as rdd
    RDLogger.DisableLog("rdApp.*")

    def rotb(s):
        if not isinstance(s, str):
            return np.nan
        m = Chem.MolFromSmiles(s)
        return np.nan if m is None else rdd.CalcNumRotatableBonds(m)

    df["rotb"] = df.smiles.map(rotb)
    return df


def rank_within_class(df: pd.DataFrame, quota: float, floor: float) -> pd.DataFrame:
    """Filter on consensus within each warhead class, then rank on enrichment.

    Returns every molecule with its filter verdict and within-class rank, not
    just the survivors -- @tt8804: "we dont want to delete any candidates".
    """
    out = []
    for cls, g in df.groupby("warhead_class"):
        g = g.copy()
        k = max(1, int(round(len(g) * quota)))
        cut = g.consensus.nlargest(k).min()
        g["class_n"] = len(g)
        g["consensus_cut"] = cut
        g["passes_consensus"] = (g.consensus >= cut) & (g.consensus >= floor)
        g["consensus_pct"] = g.consensus.rank(pct=True) * 100
        # rank is assigned ONLY among survivors; non-survivors get NaN rather
        # than a rank that would invite reading them as merely lower-placed
        s = g[g.passes_consensus].copy()
        s["class_rank"] = s.enrichment.rank(ascending=False, method="min")
        g = g.merge(s[["ident", "class_rank"]], on="ident", how="left")
        out.append(g)
    return pd.concat(out, ignore_index=True)


def class_report(df: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import spearmanr
    rows = []
    for cls, g in df.groupby("warhead_class"):
        s = g[g.passes_consensus]
        m = g.dropna(subset=["rotb", "consensus"])
        rho = spearmanr(m.rotb, m.consensus).statistic if len(m) > 10 else np.nan
        rows.append({
            "warhead_class": cls, "n": len(g), "survivors": len(s),
            "consensus_cut": g.consensus_cut.iloc[0],
            "median_consensus_kept": s.consensus.median() if len(s) else np.nan,
            "median_enrichment_kept": s.enrichment.median() if len(s) else np.nan,
            "median_rotb_kept": s.rotb.median() if len(s) else np.nan,
            "median_rotb_all": g.rotb.median(),
            "rho_rotb_consensus": rho,
        })
    return pd.DataFrame(rows).sort_values("survivors", ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--tier", choices=("T3", "T4", "both"), default="both")
    ap.add_argument("--quota", type=float, default=DEFAULT_QUOTA)
    ap.add_argument("--floor", type=float, default=CONSENSUS_FLOOR)
    ap.add_argument("--top", type=int, default=10, help="rows to print per class")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    c = add_rotb(load())
    log.info("%d molecules with consensus; T3 %d, T4 %d",
             len(c), (c.tier == "T3").sum(), (c.tier == "T4").sum())

    tiers = ["T3", "T4"] if args.tier == "both" else [args.tier]
    for tier in tiers:
        g = c[c.tier == tier].copy()
        if g.empty:
            log.warning("%s empty", tier)
            continue
        ranked = rank_within_class(g, args.quota, args.floor)
        dest = OUT.write(f"rank_2_1_{tier}", ".csv")
        ranked.sort_values(["warhead_class", "class_rank"]).to_csv(dest, index=False)

        rep = class_report(ranked)
        print(f"\n{'='*78}\n{tier}  —  {len(ranked)} molecules, quota {args.quota:.0%} "
              f"within class, consensus floor {args.floor}\n{'='*78}")
        print(f"  {'class':<24}{'n':>6}{'kept':>6}{'cut':>7}{'cons':>7}"
              f"{'enrich':>8}{'rotb':>7}{'rho(rotb,cons)':>16}")
        for r in rep.itertuples():
            print(f"  {r.warhead_class:<24}{r.n:>6}{r.survivors:>6}{r.consensus_cut:>7.2f}"
                  f"{r.median_consensus_kept:>7.2f}{r.median_enrichment_kept:>8.2f}"
                  f"{r.median_rotb_kept:>7.1f}{r.rho_rotb_consensus:>16.3f}")

        print(f"\n  top {args.top} per class, ranked on enrichment among survivors:")
        for cls, g2 in ranked[ranked.passes_consensus].groupby("warhead_class"):
            t = g2.nsmallest(args.top, "class_rank")
            print(f"\n   --- {cls} ---")
            print(f"   {'rank':>5} {'ident':<20}{'cons':>6}{'enrich':>8}{'QED':>7}{'rotb':>6}")
            for r in t.itertuples():
                q = getattr(r, "QED", np.nan)
                print(f"   {int(r.class_rank):>5} {r.ident:<20}{r.consensus:>6.2f}"
                      f"{r.enrichment:>8.2f}{q if q == q else float('nan'):>7.3f}{r.rotb:>6.0f}")
        print(f"\n  -> {dest}")

    print("\n  NOTE: one ranked list PER CLASS by design. A merged top-N would "
          "re-import\n  the cross-class rigidity bias the per-class quota exists "
          "to remove (D0073).")
    print("  The FILTER shape is defensible today; the RANKING score inside it is "
          "still\n  enrichment, which D0071 showed does not predict stability. "
          "See docs/ranking_2.1.0_design.md §3.1.")


if __name__ == "__main__":
    main()
