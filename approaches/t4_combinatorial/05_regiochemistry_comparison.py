"""
Purpose: T_4 — decide the untested BDHI and naphthoquinone attachment regiochemistries.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the latest D4 frame (post within-class ranking)
Output: a paired comparison per regiochemistry pair, written beside the frame

THE QUESTION THIS ANSWERS. Two warhead chemotypes reached enumeration with a
VERIFIED class but an UNTESTED attachment point:

  BDHI            — C3 bears the Br and is where Cys attacks, so the core must
                    attach at C4 or C5. Nobody knew which.
  1,4-naphthoquinone — attach at C2 (adjacent to the Michael acceptor) or on the
                    benzo ring (further away, electronically milder).

Rather than choose by intuition, `t4_combinatorial.yaml` enumerated all four as
separate classes and declared `discriminated_by: [warhead_validity_gate_5b,
covalent_docking_geometry, lumo_window]`. Both survived 5b and the LUMO window
did not separate them, so the discriminating evidence is the docking geometry —
which is what this script reads out.

WHY PAIRED. Each regiochemistry was enumerated against the SAME 187 R-groups, so
the two samples are matched pair-for-pair on the only other varying factor. A
Wilcoxon signed-rank test on the within-pair differences is therefore the right
test and is far more powerful than comparing two independent distributions. It
also makes the comparison robust to the R-group library's own composition: an
R-group that docks well will do so in both arms and cancels.

WHY THE MEDIAN, NOT THE BEST. Best-in-class is one lucky R-group and is a poor
estimate of whether a geometry works. The question here is whether Cys113 can
reach the reactive atom with the core in the way — a property of the attachment
geometry that should show across the whole R-group series, not in its tail.

WHAT THIS IS NOT. A docking-score difference is not a free-energy difference,
and a consistent difference could in principle be a consistent docking artifact.
This decides which regiochemistry to CARRY FORWARD, on the best evidence
available at this stage. It does not assert a measured potency difference.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                      # noqa: E402
from shared.manifest import Manifest              # noqa: E402

log = logging.getLogger("t4-regio")

EXPERIMENT = "04_t4_combinatorial"
METRIC = "affinity_kcal"          # lower is better (D0015)

# Competing attachment points for one chemotype. Order is arbitrary; the test is
# two-sided and the winner is read off the sign of the median difference.
REGIO_PAIRS = [
    ("bdhi", "bdhi_c4", "bdhi_c5"),
    ("naphthoquinone", "naphthoquinone_c2", "naphthoquinone_benzo"),
]

# Above this no-pose fraction in either arm, most of the affinity column is not
# a binding energy, and a rank test on it compares noise.
HEAVY_CENSORING = 0.5

# A Vina-style score at or above zero means the search found no favourable pose
# at all. Those are not "weak binders" on a continuum — they are failures to
# place the ligand, and counting them is more informative than averaging them.
NO_POSE_THRESHOLD = 0.0


def compare_pair(df: pd.DataFrame, label: str, a: str, b: str) -> dict | None:
    """Paired comparison of two regiochemistries across the shared R-group series."""
    da = df[df["warhead_class"] == a].set_index("rgroup_id")
    db = df[df["warhead_class"] == b].set_index("rgroup_id")
    if da.empty or db.empty:
        log.warning("%s: one arm is absent (%s=%d, %s=%d) — skipping",
                    label, a, len(da), b, len(db))
        return None

    shared = da.index.intersection(db.index)
    pair = pd.DataFrame({a: da.loc[shared, METRIC], b: db.loc[shared, METRIC]}).dropna()
    n = len(pair)
    if n < 3:
        log.warning("%s: only %d paired observations — not comparing", label, n)
        return None

    diff = pair[a] - pair[b]          # negative => a is better (lower is better)
    try:
        stat, p = stats.wilcoxon(pair[a], pair[b])
    except ValueError as exc:         # all differences zero
        log.warning("%s: Wilcoxon undefined (%s)", label, exc)
        stat, p = float("nan"), float("nan")

    # Matched-pairs rank-biserial correlation, from the SIGNED RANK SUMS —
    # r = (W+ - W-) / (W+ + W-). Computing it from win counts instead would
    # throw away the magnitudes, which is the whole reason to use a signed-rank
    # test: here one arm wins only slightly more often but wins by a wide
    # margin when it does, and a count-based effect size reports that as noise.
    nz = diff[diff != 0]
    if len(nz):
        ranks = stats.rankdata(nz.abs())
        w_pos = float(ranks[(nz > 0).to_numpy()].sum())   # b better
        w_neg = float(ranks[(nz < 0).to_numpy()].sum())   # a better
        effect = (w_pos - w_neg) / (w_pos + w_neg)
    else:
        effect = float("nan")

    n_a_better = int((diff < 0).sum())
    n_b_better = int((diff > 0).sum())
    n_tied = int((diff == 0).sum())

    # PRIMARY ENDPOINT when poses fail often. A score at or above zero is not a
    # weak binding energy, it is "no favourable pose found", so averaging it
    # with real energies is meaningless. Whether a geometry can be docked AT ALL
    # is a paired binary outcome, and McNemar is the paired test for that.
    fail_a = (pair[a] >= NO_POSE_THRESHOLD).to_numpy()
    fail_b = (pair[b] >= NO_POSE_THRESHOLD).to_numpy()
    b_only = int((fail_a & ~fail_b).sum())    # a fails, b succeeds
    a_only = int((~fail_a & fail_b).sum())    # b fails, a succeeds
    if b_only + a_only > 0:
        p_mcnemar = float(stats.binomtest(b_only, b_only + a_only, 0.5).pvalue)
    else:
        p_mcnemar = float("nan")

    winner = a if diff.median() < 0 else b
    loser = b if winner == a else a

    res = {
        "chemotype": label,
        "arm_a": a, "arm_b": b,
        "n_paired": n,
        "median_a": round(float(pair[a].median()), 3),
        "median_b": round(float(pair[b].median()), 3),
        "best_a": round(float(pair[a].min()), 3),
        "best_b": round(float(pair[b].min()), 3),
        "median_paired_difference": round(float(diff.median()), 3),
        "wilcoxon_statistic": None if np.isnan(stat) else float(stat),
        "p_value": None if np.isnan(p) else float(p),
        "rank_biserial_effect": round(float(effect), 3),
        "n_pairs_a_better": n_a_better,
        "n_pairs_b_better": n_b_better,
        "n_pairs_tied": n_tied,
        "no_pose_fraction_a": round(float(fail_a.mean()), 3),
        "no_pose_fraction_b": round(float(fail_b.mean()), 3),
        "n_discordant_b_poses_a_does_not": b_only,
        "n_discordant_a_poses_b_does_not": a_only,
        "p_value_mcnemar_pose_success": None if np.isnan(p_mcnemar) else p_mcnemar,
        "heavy_censoring": bool(max(fail_a.mean(), fail_b.mean()) > HEAVY_CENSORING),
        "primary_endpoint": ("pose_success"
                             if max(fail_a.mean(), fail_b.mean()) > HEAVY_CENSORING
                             else "affinity"),
        "winner": winner,
        "loser": loser,
    }
    return res


def verdict(r: dict) -> str:
    """Graded verdict, in the vocabulary the gates already use.

    Which endpoint is primary depends on how much of the affinity column is
    real. This rule is stated in terms of the DATA STRUCTURE, not the outcome:

    - **Light censoring** — both endpoints must agree for STRONG: the winner
      docks more often (McNemar on paired pose success) *and* scores better
      when both dock (signed-rank on affinity). Agreement between a binary and
      a continuous read of the same geometry is harder to produce by artifact
      than either alone.
    - **Heavy censoring** (>50% no-pose in either arm) — pose success is
      primary and the affinity test is reported but not graded on. You cannot
      compare the magnitudes of numbers that are not binding energies.

    DISCLOSURE: the censoring branch was added after seeing that
    `naphthoquinone_c2` fails to pose 97% of the time — i.e. after seeing the
    data, which is exactly the circumstance in which a threshold change can
    become a way of choosing the answer. It is recorded here because it does
    not change any conclusion: both endpoints favour the same arm in both
    pairs, so only the confidence LABEL moves, never the winner. The rule was
    fixed before the thresholds were set, and the alternative — grading the
    naphthoquinone comparison on a rank test over 97%-censored data — would be
    wrong on its own terms whichever arm it happened to favour.
    """
    p = r["p_value"]
    p_pose = r["p_value_mcnemar_pose_success"]
    eff = abs(r["rank_biserial_effect"])
    gap = abs(r["median_paired_difference"])
    if p is None:
        return "UNDECIDED"

    heavy = max(r["no_pose_fraction_a"], r["no_pose_fraction_b"]) > HEAVY_CENSORING
    pose_favours_winner = (
        (r["no_pose_fraction_a"] > r["no_pose_fraction_b"]) == (r["winner"] == r["arm_b"])
    )
    pose_agrees = p_pose is not None and p_pose < 0.05 and pose_favours_winner

    if heavy:
        if p_pose is None or not pose_favours_winner:
            return "UNDECIDED"
        disc_win = r["n_discordant_b_poses_a_does_not"]
        disc_lose = r["n_discordant_a_poses_b_does_not"]
        if r["winner"] == r["arm_a"]:
            disc_win, disc_lose = disc_lose, disc_win
        if p_pose < 0.001 and disc_win >= 5 * max(disc_lose, 1):
            return "STRONG"
        if p_pose < 0.05:
            return "WEAK"
        return "UNDERPOWERED"

    if p < 0.001 and eff >= 0.5 and gap >= 1.0 and pose_agrees:
        return "STRONG"
    if p < 0.001 and eff >= 0.3 and gap >= 1.0:
        return "WEAK"
    if p < 0.05 and eff >= 0.3:
        return "WEAK"
    if p < 0.05:
        return "UNDERPOWERED"
    return "NO_DIFFERENCE"


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Decide the untested BDHI / naphthoquinone regiochemistries.")
    ap.add_argument("--out-name", default="regiochemistry_comparison",
                    help="stem for the JSON written beside the frame")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    df, frame_path = dio.latest_frame(EXPERIMENT, "t4")
    log.info("loaded %s (%d rows)", frame_path.name, len(df))
    if METRIC not in df.columns:
        raise SystemExit(f"frame has no {METRIC!r} — run 03_covalent_dock.py first")

    results = [r for r in
               (compare_pair(df, label, a, b) for label, a, b in REGIO_PAIRS)
               if r is not None]
    if not results:
        raise SystemExit("no regiochemistry pair could be compared")

    for r in results:
        r["verdict"] = verdict(r)

    out_dir = dio.approach_dir("t4", EXPERIMENT)
    out = dio.next_version(out_dir, args.out_name, ".json")
    out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    (Manifest(stage="t4_regiochemistry_comparison", approach="t4",
              params={"test": "Wilcoxon signed-rank, paired on rgroup_id",
                      "metric": f"{METRIC} (lower better, D0015)",
                      "effect_size": "matched-pairs rank-biserial correlation",
                      "no_pose_threshold_kcal": NO_POSE_THRESHOLD})
     .add_input("d4_frame", frame_path)
     .add_output("comparison", out)
     .note("; ".join(f"{r['chemotype']}: {r['winner']} ({r['verdict']})"
                     for r in results))
     .write(out_dir, filename=f"{out.stem}_manifest.json"))

    print(f"\nT_4 regiochemistry comparison -> {out}\n")
    for r in results:
        print(f"  {r['chemotype'].upper()}  ({r['n_paired']} paired R-groups)")
        print(f"    {r['arm_a']:22s} median {r['median_a']:6.2f}  best {r['best_a']:6.2f}  "
              f"no pose {r['no_pose_fraction_a']:.0%}")
        print(f"    {r['arm_b']:22s} median {r['median_b']:6.2f}  best {r['best_b']:6.2f}  "
              f"no pose {r['no_pose_fraction_b']:.0%}")
        p = r["p_value"]
        print(f"    paired difference     median {r['median_paired_difference']:+.2f} kcal/mol, "
              f"p = {p:.3g}" if p is not None else "    p undefined")
        print(f"    effect (rank-biserial) {r['rank_biserial_effect']:+.3f}   "
              f"{r['arm_b']} better in {r['n_pairs_b_better']}/{r['n_paired']} pairs")
        pp = r["p_value_mcnemar_pose_success"]
        print(f"    pose success (McNemar) {r['n_discordant_b_poses_a_does_not']} pairs "
              f"where only {r['arm_b']} poses vs "
              f"{r['n_discordant_a_poses_b_does_not']} where only {r['arm_a']} does"
              + (f", p = {pp:.3g}" if pp is not None else ""))
        print(f"    -> {r['verdict']}: carry {r['winner']}, drop {r['loser']}\n")

    print("  A docking-score difference is not a measured potency difference. This")
    print("  chooses which attachment to carry forward, on the evidence available.")


if __name__ == "__main__":
    main()
