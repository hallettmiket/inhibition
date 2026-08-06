"""
Purpose: does the energy of a near-attack pose add anything over how often one is reached?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: 00_outputs/blacksmith/nac_robust/*.csv (labelled, carries per-pose energies)
Output: a report — no new artefact, this decides whether stage 4 exists at all

PRE-REGISTERED. Written and committed BEFORE the robustness run produced the data
it reads, so the candidate rules below were fixed without knowing which wins.
That is the whole point: with five rules and one labelled set, picking the rule
after seeing the answer would manufacture a result, and D0045 exists because that
has happened on this project before.

THE QUESTION. `enrichment` measures how OFTEN a molecule reaches a
mechanism-appropriate near-attack conformation. It says nothing about whether it
gets there in a good pose or a strained one. `docs/ranking_rationale.md` stage 4
proposed ranking the survivors by the stability of the near-attack geometry;
per-pose binding energies are the cheap version of that question.

WHAT WOULD MAKE STAGE 4 REAL. A rule combining geometry with energy must beat
enrichment alone, on the same labelled set, by more than the noise. If none does,
stage 4 as specified does not exist and enrichment is the whole ranking -- which
is a perfectly good outcome and must be reported as such rather than quietly
dropped.

THE ENERGY IS NOT TRUSTED TO RANK MOLECULES. Five measurements say the docking
score carries no signal on this target (D0041, D0046, D0036, D0038/D0044, D0061),
and rule E1 below is included precisely to re-confirm that on this set rather
than assume it. The rules that matter (C1-C3) ask it only to compare poses OF ONE
MOLECULE that have already cleared the geometric gate -- a much weaker question
than "which molecule binds better".
"""

from __future__ import annotations

import glob
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout              # noqa: E402

log = logging.getLogger("nac-stage4")
ROBUST = sout.Topic("blacksmith", "nac_robust")

# --------------------------------------------------------------------------
# THE PRE-REGISTERED RULES. Fixed before the data existed. Higher = better in
# every case, so all are scored the same way and none needs a sign flip decided
# after the fact.
# --------------------------------------------------------------------------

RULES = {
    # The incumbent. Everything else has to beat this.
    "G1_enrichment":
        lambda d: d.enrichment,

    # Energy alone, as a control. Expected to fail -- five prior measurements say
    # the score carries no signal here. If it SUCCEEDS, something is wrong with
    # the earlier work or with this set, and that matters more than stage 4.
    "E1_best_viable_dg_alone":
        lambda d: -d.best_viable_dg,

    # C1: gate on geometry, then rank by the best near-attack pose's energy.
    # The literal reading of ranking_rationale stage 4 -- "binary gate, then
    # continuous rank".
    "C1_gate_then_energy":
        lambda d: np.where(d.enrichment >= 1.0, -d.best_viable_dg, -np.inf),

    # C2: reward reaching a NAC often AND well. Multiplicative, so a molecule
    # must do both; a good energy cannot rescue a molecule that never gets there.
    "C2_enrichment_x_energy":
        lambda d: d.enrichment * (-d.best_viable_dg),

    # C3: does the molecule PAY to satisfy the geometry? The gap between its best
    # pose overall and its best pose that clears the gate. Small gap = the
    # near-attack geometry is nearly free; large gap = it must contort to react.
    # This is the only rule here that is not a re-reading of the score's level,
    # and the one most likely to carry information enrichment does not.
    "C3_nac_energy_penalty":
        lambda d: -(d.best_viable_dg - d.best_dg),
}


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC, ties as half. NaNs rank last, never dropped."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    worst = np.nanmin([np.nanmin(pos), np.nanmin(neg)]) - 1e9
    pos = np.where(np.isnan(pos), worst, pos)
    neg = np.where(np.isnan(neg), worst, neg)
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def bootstrap_delta(p1, n1, p0, n0, reps: int = 4000) -> tuple[float, float, float]:
    """CI on the AUC DIFFERENCE between a rule and the incumbent.

    Resampled on the SAME molecules for both rules, so the paired structure is
    kept. Comparing two independent CIs instead would overstate the uncertainty
    of the difference and hide a real improvement -- or invent one.
    """
    rng = np.random.default_rng(0x5747E4)
    d = []
    for _ in range(reps):
        i = rng.integers(0, len(p1), len(p1))
        j = rng.integers(0, len(n1), len(n1))
        d.append(auc(p1[i], n1[j]) - auc(p0[i], n0[j]))
    d = np.array(d)
    lo, hi = np.percentile(d, [2.5, 97.5])
    return float(np.mean(d)), float(lo), float(hi)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fs = sorted(glob.glob(str(ROBUST.dir / "nac_robust_s*.csv")))
    if not fs:
        raise SystemExit("no robustness output yet — stage 4 needs its labelled set")
    df = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    df = df.drop_duplicates("ident", keep="first")
    ok = df[df.status == "ok"].copy()
    if "best_viable_dg" not in ok.columns:
        raise SystemExit("this run predates per-pose energy capture")

    print(f"\n=== stage 4: does pose energy add to pose frequency? ===")
    print(f"  {len(ok)} molecules "
          f"({(ok.label == 'positive').sum()} positive, "
          f"{(ok.label == 'negative').sum()} negative)")
    print(f"  molecules with NO viable pose (best_viable_dg is NaN): "
          f"{int(ok.best_viable_dg.isna().sum())} — ranked last, never dropped")

    # Per class, because enrichment is only comparable within a mechanism and
    # the energy scale differs between chemotypes too.
    for cls, g in ok.groupby("warhead_class"):
        p, n = g[g.label == "positive"], g[g.label == "negative"]
        if len(p) < 2 or len(n) < 5:
            print(f"\n  {cls}: {len(p)} pos / {len(n)} neg — too few, skipped")
            continue
        print(f"\n  {cls}  ({len(p)} positives, {len(n)} negatives)")
        base_p, base_n = RULES["G1_enrichment"](p).values, RULES["G1_enrichment"](n).values
        base = auc(base_p, base_n)
        rows = []
        for name, fn in RULES.items():
            vp, vn = np.asarray(fn(p), float), np.asarray(fn(n), float)
            a = auc(vp, vn)
            if name == "G1_enrichment":
                rows.append((name, a, 0.0, 0.0, 0.0))
                continue
            m, lo, hi = bootstrap_delta(vp, vn, base_p, base_n)
            rows.append((name, a, m, lo, hi))
        for name, a, m, lo, hi in rows:
            if name == "G1_enrichment":
                print(f"    {name:<26} AUC {a:.3f}   (incumbent)")
            else:
                verdict = ("BEATS enrichment" if lo > 0 else
                           "worse than enrichment" if hi < 0 else
                           "indistinguishable")
                print(f"    {name:<26} AUC {a:.3f}   delta {m:+.3f} "
                      f"[{lo:+.3f}, {hi:+.3f}]   {verdict}")

    print("\n  READING, fixed in advance: a rule earns stage 4 only if its delta CI")
    print("  excludes zero. If none does, enrichment IS the ranking and stage 4 as")
    print("  specified does not exist — which is a result, not a gap.")


if __name__ == "__main__":
    main()
