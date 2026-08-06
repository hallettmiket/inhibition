"""
Purpose: does the criterion measure Cys113 recognition, or just "can a warhead point at a sulfur"?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: the 15 crystallographic positives + matched negatives, scored at a DECOY cysteine
Output: 00_outputs/blacksmith/nac_decoy/nac_decoy_s<shard>_<N>.csv + a report

THE CONTROL THAT SEPARATES TWO EXPLANATIONS OF THE SAME RESULT.

Crystallographic Cys113 binders enrich 2.39x over chance at Cys113 while
warhead-matched measured inactives sit at 0.82-1.14x (D0065). Two readings:

  (a) the criterion measures RECOGNITION -- these molecules are shaped to bind
      this pocket with the warhead presented to this cysteine; or
  (b) the criterion measures GENERIC warhead accessibility -- they simply have
      exposed, unhindered warheads that could point at any sulfur anywhere.

Only (a) supports using this to rank candidates for Pin1. Nothing measured so far
distinguishes them, because everything has been scored at one site.

THE DECOY. Pin1's other cysteine, **Cys57**, 12.4 A from Cys113 in the same
chain. Same protein, same file, same preparation, same reactive parameterisation,
same box size -- only the site differs (box centre 12.97/14.65/0.40 against
14.04/7.19/-2.11). Every alternative explanation involving preparation,
protonation or parameters is held fixed by construction.

    PRE-REGISTERED READING. If the positives enrich at Cys113 and NOT at Cys57,
    reading (a) holds and the signal is site-specific. If they enrich equally at
    both, reading (b) holds, the ranking is measuring warhead exposure rather
    than Pin1 recognition, and D0065's interpretation must be withdrawn even
    though its AUCs stand.

WHAT THIS IS NOT. Not a claim that Cys57 is undruggable or that nothing binds
there. The comparison is strictly within-molecule: the same molecules, the same
criterion, two sites.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout              # noqa: E402
import nac_screen as ns                         # noqa: E402
import nac_rank as nr                           # noqa: E402
import nac_robustness as rb                     # noqa: E402

log = logging.getLogger("nac-decoy")

OUT = sout.Topic("blacksmith", "nac_decoy")
DECOY = Path("/data/lab_vm/modifiable/inhibition/receptor_3ikd_decoy_cys57")
N_NEG_PER_CLASS = 40          # enough to place the positives, not a validation set


def build_set() -> list[ns.Candidate]:
    """The crystallographic positives plus a modest matched negative sample.

    The negatives are here to give the positives a reference distribution AT THE
    DECOY SITE, not to re-validate anything -- the question is whether the
    positives' advantage survives the move, so they are the ones that matter.
    """
    saved = rb.NEG_PER_CLASS
    try:
        rb.NEG_PER_CLASS = N_NEG_PER_CLASS
        return rb.build_set()
    finally:
        rb.NEG_PER_CLASS = saved


def run_shard(shard: int, n_shards: int, nrun: int, gpu: str, chunk: int) -> None:
    if not (DECOY / "rec.reactive_config").is_file():
        raise SystemExit(f"decoy receptor not built at {DECOY}")
    mine = [c for i, c in enumerate(build_set()) if i % n_shards == shard]
    done = set()
    for f in glob.glob(str(OUT.dir / "nac_decoy_s*.csv")):
        try:
            done.update(pd.read_csv(f).ident.astype(str))
        except Exception:                              # noqa: BLE001
            pass
    todo = [c for c in mine if c.ident not in done]
    log.info("shard %d/%d: %d assigned, %d to do", shard, n_shards, len(mine), len(todo))

    buf = []
    for k, c in enumerate(todo, 1):
        # DECOY receptor dir — the only thing that differs from the main run.
        row = nr.score_one(c, DECOY, nrun, gpu)
        row["label"] = c.label
        buf.append(row)
        if len(buf) >= chunk or k == len(todo):
            dest = OUT.write(f"nac_decoy_s{shard}", ".csv")
            pd.DataFrame(buf).to_csv(dest, index=False)
            log.info("shard %d: %d/%d -> %s", shard, k, len(todo), dest.name)
            buf = []


def report() -> None:
    """Compare each molecule's enrichment at the decoy site with its value at Cys113."""
    dfs = []
    for f in glob.glob(str(OUT.dir / "nac_decoy_s*.csv")):
        dfs.append(pd.read_csv(f))
    if not dfs:
        raise SystemExit("no decoy output yet")
    decoy = pd.concat(dfs, ignore_index=True).drop_duplicates("ident")
    decoy = decoy[decoy.status == "ok"]

    main = []
    for f in glob.glob(str(sout.Topic("blacksmith", "nac_robust").dir / "nac_robust_s*.csv")):
        main.append(pd.read_csv(f))
    if not main:
        raise SystemExit("no Cys113 comparator (run nac_robustness first)")
    cys113 = pd.concat(main, ignore_index=True).drop_duplicates("ident")
    cys113 = cys113[cys113.status == "ok"]

    m = decoy.merge(cys113[["ident", "enrichment"]], on="ident",
                    suffixes=("_decoy", "_cys113"))
    print(f"\n=== decoy-site control: {len(m)} molecules scored at BOTH sites ===")
    if m.empty:
        return
    for lab, g in m.groupby("label"):
        print(f"\n  {lab}  (n={len(g)})")
        print(f"    Cys113 median enrichment  {g.enrichment_cys113.median():>5.2f}x")
        print(f"    Cys57  median enrichment  {g.enrichment_decoy.median():>5.2f}x")
        print(f"    per-molecule ratio Cys113/Cys57  median "
              f"{(g.enrichment_cys113 / g.enrichment_decoy.replace(0, np.nan)).median():>5.2f}")

    pos, neg = m[m.label == "positive"], m[m.label == "negative"]
    if len(pos) >= 2 and len(neg) >= 5:
        from scipy.stats import mannwhitneyu, wilcoxon
        for site, col in (("Cys113", "enrichment_cys113"), ("Cys57 (decoy)", "enrichment_decoy")):
            u, p = mannwhitneyu(pos[col], neg[col], alternative="greater")
            print(f"\n  positives vs negatives at {site:<14} "
                  f"AUC {u/(len(pos)*len(neg)):.3f}  p={p:.4f}")
        try:
            _, pw = wilcoxon(pos.enrichment_cys113, pos.enrichment_decoy,
                             alternative="greater")
            print(f"\n  positives enrich MORE at Cys113 than at Cys57 "
                  f"(paired Wilcoxon): p={pw:.4f}")
        except ValueError as exc:
            print(f"\n  paired test unavailable: {exc}")
    print("\n  READING, fixed in advance: site-specific signal means positives separate")
    print("  at Cys113 and NOT at Cys57. Equal separation at both means the criterion")
    print("  measures warhead exposure, not Pin1 recognition — and D0065's")
    print("  interpretation would have to be withdrawn even though its AUCs stand.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--chunk", type=int, default=40)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [d{args.shard}] %(message)s")
    if args.report:
        report()
        return
    run_shard(args.shard, args.n_shards, args.nrun, args.gpu, args.chunk)


if __name__ == "__main__":
    main()
