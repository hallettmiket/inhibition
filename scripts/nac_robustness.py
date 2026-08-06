"""
Purpose: is the validation AUC real, or an artefact of which 30 negatives were drawn?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: the 15 crystallographic positives + a large warhead-matched negative pool
Output: 00_outputs/blacksmith/nac_robust/nac_robust_s<shard>_<N>.csv + a report

THE CHECK THE HEADLINE RESULT NEEDS BEFORE ANYONE BUILDS ON IT.

`nac_screen` reported chloroacetamide AUC 0.822 (p=0.0020) against **30**
warhead-matched measured inactives, drawn once. That is a single sample from a
pool of 642, and a single sample can be lucky. The number has already been shown
to move: four independent runs gave 0.872, 0.881, 0.852, 0.822 -- but those
varied the DOCKING seed while holding the negatives fixed, so they measure only
one of the two sources of variance.

This measures the other, and does it by MEASUREMENT rather than by resampling
theory:

  1. Score a large negative pool -- 300 per class, ten times the original draw,
     shuffled with a DIFFERENT seed so it is not a superset of a set already
     known to work.
  2. Split it into TEN DISJOINT subsets of 30 and compute the AUC of each.
     Disjoint, not bootstrapped: ten independent answers to the actual question
     asked, "what would we have concluded had we drawn a different 30".
  3. Report the full-pool AUC with a bootstrap CI as well, since 300 negatives
     is the better point estimate regardless.

WHAT WOULD FALSIFY THE RESULT. If the ten disjoint AUCs straddle 0.5, the
headline was a draw artefact and the claim must be withdrawn. If they cluster
above 0.7, the choice of negatives was not what produced it.

WHY DISJOINT SUBSETS AND NOT ONLY A BOOTSTRAP. A bootstrap resamples the pool it
is given and therefore inherits whatever is peculiar about that pool. Ten
disjoint draws of genuinely different compounds is the stronger evidence, and it
is the design that answers the question in the form it was asked.
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

from shared import nac_criterion as nac         # noqa: E402
from shared import outputs as sout              # noqa: E402
import nac_screen as ns                         # noqa: E402
import nac_rank as nr                           # noqa: E402

log = logging.getLogger("nac-robust")

OUT = sout.Topic("blacksmith", "nac_robust")
NEG_PER_CLASS = 300
SUBSET_SIZE = 30                 # the size the headline result used
# Deliberately NOT nac_screen's 0xC0FFEE: reusing that seed would make the pool
# a superset of a draw already known to give a positive answer.
POOL_SEED = 0xBEEF11


def build_set() -> list[ns.Candidate]:
    """The 15 crystallographic positives plus a large, independently-drawn negative pool."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")
    meta = {r.class_id: (r.mechanism, r.reactive_atom_smarts) for r in wh.itertuples()}
    pos = ns.crystal_positives(meta, None)
    want = {c.warhead_class for c in pos}
    log.info("positives: %d across %s", len(pos), sorted(want))

    inact = pd.read_csv("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith/"
                        "measured_inactives/aid504891_inactives_1.csv")
    pool = inact.dropna(subset=["canonical_smiles"]).sample(frac=1.0,
                                                            random_state=POOL_SEED)
    neg = []
    for cls in sorted(want):
        mech, smarts = meta[cls]
        patt = Chem.MolFromSmarts(smarts)
        taken = 0
        for r in pool.itertuples():
            if taken >= NEG_PER_CLASS:
                break
            m = ns.largest_fragment(r.canonical_smiles)
            if m is None or not m.HasSubstructMatch(patt):
                continue
            neg.append(ns.Candidate(f"neg:{cls}:{int(r.PUBCHEM_CID)}",
                                    r.canonical_smiles, cls, mech, smarts, "negative"))
            taken += 1
        log.info("  %s: %d negatives", cls, taken)

    allc = pos + neg
    allc.sort(key=lambda c: c.ident)
    return allc


def already_done() -> set[str]:
    ids = set()
    for f in glob.glob(str(OUT.dir / "nac_robust_s*.csv")):
        try:
            ids.update(pd.read_csv(f).ident.astype(str))
        except Exception as exc:                       # noqa: BLE001
            log.warning("unreadable chunk %s: %s", Path(f).name, exc)
    return ids


def run_shard(shard: int, n_shards: int, nrun: int, gpu: str, chunk: int) -> None:
    rec = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    mine = [c for i, c in enumerate(build_set()) if i % n_shards == shard]
    done = already_done()
    todo = [c for c in mine if c.ident not in done]
    log.info("shard %d/%d: %d assigned, %d to do", shard, n_shards, len(mine), len(todo))

    buf = []
    for k, c in enumerate(todo, 1):
        row = nr.score_one(c, rec, nrun, gpu)
        row["label"] = c.label
        buf.append(row)
        if len(buf) >= chunk or k == len(todo):
            dest = OUT.write(f"nac_robust_s{shard}", ".csv")
            pd.DataFrame(buf).to_csv(dest, index=False)
            log.info("shard %d: %d/%d, wrote %d -> %s", shard, k, len(todo),
                     len(buf), dest.name)
            buf = []


# --------------------------------------------------------------------------

def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC, ties counted as half — the standard convention."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def report() -> None:
    fs = sorted(glob.glob(str(OUT.dir / "nac_robust_s*.csv")))
    if not fs:
        raise SystemExit("nothing written yet")
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df = df.drop_duplicates("ident", keep="first")
    ok = df[df.status == "ok"].copy()
    print(f"\n=== robustness: {len(df)} scored, {len(ok)} ok, "
          f"{len(df)-len(ok)} failed ===")
    if ok.empty:
        return

    rng = np.random.default_rng(0xD15EA5E)
    for cls, g in ok.groupby("warhead_class"):
        p = g[g.label == "positive"].enrichment.values
        n = g[g.label == "negative"].enrichment.values
        if len(p) < 2 or len(n) < SUBSET_SIZE:
            print(f"\n  {cls}: {len(p)} pos / {len(n)} neg — too few to split")
            continue
        full = auc(p, n)

        # TEN DISJOINT SUBSETS: the direct answer to "what if we had drawn a
        # different 30". Not resampled — genuinely different compounds each time.
        idx = rng.permutation(len(n))
        k = min(10, len(n) // SUBSET_SIZE)
        subs = [auc(p, n[idx[i*SUBSET_SIZE:(i+1)*SUBSET_SIZE]]) for i in range(k)]

        boot = [auc(p[rng.integers(0, len(p), len(p))],
                    n[rng.integers(0, len(n), len(n))]) for _ in range(2000)]
        lo, hi = np.percentile(boot, [2.5, 97.5])

        print(f"\n  {cls}  ({len(p)} positives, {len(n)} negatives)")
        print(f"    full-pool AUC          {full:.3f}   bootstrap 95% CI [{lo:.3f}, {hi:.3f}]")
        print(f"    {k} DISJOINT subsets of {SUBSET_SIZE}:")
        print(f"      {'  '.join(f'{s:.3f}' for s in subs)}")
        print(f"      min {min(subs):.3f}   median {np.median(subs):.3f}   max {max(subs):.3f}")
        verdict = ("ROBUST — every disjoint draw agrees" if min(subs) > 0.6 else
                   "FRAGILE — some draws straddle chance" if min(subs) < 0.5 else
                   "MIXED — directionally consistent, spread is wide")
        print(f"      {verdict}")

    p = ok[ok.label == "positive"].enrichment.values
    n = ok[ok.label == "negative"].enrichment.values
    print(f"\n  POOLED on enrichment: AUC {auc(p, n):.3f}  "
          f"({len(p)} positives vs {len(n)} negatives)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--gpu", default="1")
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format=f"%(levelname)s [r{args.shard}] %(message)s")
    if args.report:
        report()
        return
    run_shard(args.shard, args.n_shards, args.nrun, args.gpu, args.chunk)


if __name__ == "__main__":
    main()
