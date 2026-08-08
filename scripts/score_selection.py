"""
Purpose: choose the 2.2.0 ranking score against criteria fixed before the numbers.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: nac_v3 per-mode rows, the screened references, and (optionally) a
       convergence replicate at a different sampling depth
Output: 00_outputs/blacksmith/score_selection/score_selection_<N>.csv

Scores `docs/prereg_score_selection.md`. Six candidate scores against three
tests, with the thresholds written down first.

WHY A SCRIPT RATHER THAN A PARAGRAPH. Comparing six scores on one validation set
and reporting the best is how a score gets chosen by noise -- and this project
has already shipped a ranking that put the crystallographically-confirmed parent
compound dead last of 5,765. The reading table is transcribed here as a literal
and the verdict is LOOKED UP, so the answer cannot be adjusted after seeing it.

T1 CONVERGENCE IS DISQUALIFYING AND IT COMES FIRST, because a score that does
not reproduce cannot be validated by anything else. A good AUC on an
irreproducible score is a good measurement of that day's noise. D0068 has asked
for this since the beginning and nothing has yet passed it.

T2 USES A VALIDATION SET THAT HAS NEVER BEEN USED. The 22 reference molecules
carry potency annotations -- nanomolar covalent leads on one side, historical
promiscuous compounds (KPT-6566, Juglone, ATRA) on the other. Only 10 carry a
warhead this criterion can score, and the promiscuous three are all
naphthoquinone-type, so the comparison is confounded with warhead class. Both the
pooled AUC and a within-class version are reported; neither is treated as more
than the weak instrument it is at n = 10.
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

from shared import outputs as sout                  # noqa: E402

log = logging.getLogger("score-sel")
DATA = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith")
OUT = sout.Topic("blacksmith", "score_selection")

SCORES = ["viable_fraction", "enrichment_joint", "enrichment_conditional",
          "conditional_x_consensus", "conditional_lcb", "anchor_quality_max"]

#: Reference potency classes, from reference_assignment_1.csv's own annotations.
PROMISCUOUS = {"KPT-6566", "Juglone", "ATRA", "EGCG", "PiB"}

#: The reading table, transcribed from the pre-registration. The verdict is
#: looked up in this, never argued.
READINGS = {
    "t1_sampling":  [(0.0, 0.5, "DISQUALIFIED"), (0.5, 0.7, "acceptable"), (0.7, 1.01, "good")],
    "t1_rerun":     [(0.0, 0.6, "DISQUALIFIED"), (0.6, 0.8, "acceptable"), (0.8, 1.01, "good")],
    "t2_sulfopin":  [(0.0, 25.0, "DISQUALIFIED"), (25.0, 60.0, "acceptable"), (60.0, 101.0, "good")],
    "t2_auc":       [(0.0, 0.5, "DISQUALIFIED"), (0.5, 0.7, "acceptable"), (0.7, 1.01, "good")],
    "t3_energy":    [(0.0, 0.2, "good"), (0.2, 0.4, "acceptable"), (0.4, 1.01, "DISQUALIFIED")],
}


def verdict(metric: str, value: float) -> str:
    if value != value:
        return "not measured"
    for lo, hi, txt in READINGS[metric]:
        if lo <= value < hi:
            return txt
    return "?"


def load(topic: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(DATA / topic / "agg_s*_*.csv")))
    fs += sorted(glob.glob(str(DATA / topic / "agg_ref_*.csv")))
    fs += sorted(glob.glob(str(DATA / topic / "refs_*.csv")))
    if not fs:
        raise SystemExit(f"no aggregates under {DATA / topic}")
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return d[d.status == "ok"].copy()


def add_scores(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["enrichment_joint"] = d.get("enrichment", np.nan)
    n = d.n_in_range.astype(float).replace(0, np.nan)
    phat = d.n_viable_given_in_range / n
    d["enrichment_conditional"] = phat / d.isotropic_null
    if "consensus" in d.columns:
        d["conditional_x_consensus"] = d.enrichment_conditional * d.consensus
    z = 1.96
    den = 1.0 + z * z / n
    centre = (phat + z * z / (2 * n)) / den
    half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / den
    d["conditional_lcb"] = (centre - half).clip(lower=0.0) / d.isotropic_null
    return d


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC. Ties count a half, which is what makes it an AUC."""
    pos = pos[~np.isnan(pos)]; neg = neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    gt = sum((p > neg).sum() + 0.5 * (p == neg).sum() for p in pos)
    return float(gt / (len(pos) * len(neg)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--topic", default="nac_v3")
    ap.add_argument("--converge-dir", default=None,
                    help="a replicate at a different sampling depth, for T1")
    ap.add_argument("--rerun-dir", default=None,
                    help="an independent re-run at the SAME depth, for T1")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from scipy.stats import spearmanr

    d = add_scores(load(args.topic))
    # one row per molecule, on its dominant mode, so a molecule counts once
    key = "parent_ident" if "parent_ident" in d.columns else "ident"
    d[key] = d[key].fillna(d["ident"])
    dom = d.sort_values("consensus", ascending=False).drop_duplicates(key, keep="first")
    refs = dom[dom.get("is_reference", False) == True]           # noqa: E712
    lib = dom[dom.get("is_reference", False) != True]            # noqa: E712
    log.info("%d library molecules, %d reference rows", len(lib), len(refs))

    def refname(i):
        m = re.match(r"ref_(.+?)__", str(i))
        return m.group(1).replace("-", " ") if m else str(i)
    if len(refs):
        refs = refs.copy()
        refs["name"] = refs.ident.map(refname)
        refs["promiscuous"] = refs.name.str.replace(" ", "-").isin(PROMISCUOUS)

    rows = []
    for sc in SCORES:
        if sc not in dom.columns:
            log.warning("%s: absent, skipping", sc)
            continue
        rec = {"score": sc}
        # T2 — Sulfopin percentile and lead-vs-promiscuous AUC
        if len(refs):
            sp = refs[refs.name.str.contains("Sulfopin", case=False, na=False)]
            if len(sp) and lib[sc].notna().any():
                rec["t2_sulfopin"] = float((lib[sc] < sp[sc].max()).mean() * 100)
            leads = refs[~refs.promiscuous][sc].values
            prom = refs[refs.promiscuous][sc].values
            rec["t2_auc"] = auc(leads, prom)
            rec["n_leads"], rec["n_promiscuous"] = len(leads), len(prom)
        # T3 — independence from docking energy
        if "mean_energy" in lib.columns:
            r, _ = spearmanr(lib[sc], lib.mean_energy, nan_policy="omit")
            rec["t3_energy"] = abs(float(r))
        # T1 — convergence, if replicates were supplied
        for tag, path in (("t1_sampling", args.converge_dir), ("t1_rerun", args.rerun_dir)):
            if not path:
                continue
            try:
                o = add_scores(load(Path(path).name))
            except SystemExit:
                continue
            ok2 = o.sort_values("consensus", ascending=False).drop_duplicates(
                "parent_ident" if "parent_ident" in o.columns else "ident", keep="first")
            j = dom.merge(ok2, on=key, suffixes=("_a", "_b"))
            if len(j) > 10 and f"{sc}_a" in j and f"{sc}_b" in j:
                r, _ = spearmanr(j[f"{sc}_a"], j[f"{sc}_b"], nan_policy="omit")
                rec[tag] = float(r)
        rows.append(rec)

    t = pd.DataFrame(rows)
    for m in READINGS:
        if m in t.columns:
            t[f"{m}_verdict"] = t[m].map(lambda v: verdict(m, v))
    dest = OUT.write("score_selection", ".csv")
    t.to_csv(dest, index=False)

    print("\n" + "=" * 78)
    print("  SCORE SELECTION — docs/prereg_score_selection.md, as written")
    print("=" * 78 + "\n")
    show = ["score"] + [c for c in ("t1_sampling", "t1_rerun", "t2_sulfopin",
                                    "t2_auc", "t3_energy") if c in t.columns]
    print(t[show].round(3).to_string(index=False))
    dq = [c for c in t.columns if c.endswith("_verdict")]
    if dq:
        print("\n  verdicts, looked up in the pre-registered table:")
        for r in t.itertuples():
            bad = [c.replace("_verdict", "") for c in dq
                   if getattr(r, c, "") == "DISQUALIFIED"]
            state = f"DISQUALIFIED on {', '.join(bad)}" if bad else "eligible"
            print(f"    {r.score:<26} {state}")
    if "t1_sampling" not in t.columns:
        print("\n  T1 NOT MEASURED — no convergence replicate supplied, and T1 is")
        print("  DISQUALIFYING and ranked first. No score can be chosen until it runs.")
    print(f"\n  -> {dest}")


if __name__ == "__main__":
    main()
