"""
Purpose: score the pre-registered elevation experiment against the readings fixed before it ran.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: 00_outputs/blacksmith/elevation_tier1/*.csv and .../bpmd/bpmd_s*.csv
Output: 00_outputs/blacksmith/elevation_analysis/elevation_<tier|readings>_<N>.csv

THE READINGS WERE FIXED IN ADVANCE AND THIS SCRIPT EVALUATES THEM, IT DOES NOT
CHOOSE THEM. `docs/elevation_prereg.md` names four mutually-exclusive patterns
over the BDHI groups (B ~ A > D; A > B > D; A ~ B ~ D; D >= A,B) and one
descriptive statement about the validated class (V ~ REF). Whichever pattern the
data shows, its conclusion is the conclusion -- including "report as a failure,
do not reinterpret".

ONE THING THE PREREG DID NOT FIX, AND IT IS NAMED RATHER THAN SMUGGLED. It fixed
the readings but not the numeric line between "~" and ">". A rule is therefore
chosen HERE, at analysis time, and it is declared post-hoc wherever it is
reported:

    ">"  means the Mann-Whitney test is significant after Holm correction across
         the three pre-registered contrasts AND |Cliff's delta| >= 0.5
    "~"  means anything else

Both halves are required because either alone is misreadable at n = 8: a p-value
because n = 8 makes significance itself a coin-flip on one or two molecules, and
an effect size because a large delta on eight points is routine under the null.

NO SIGNIFICANCE CLAIM IS COMPUTED FOR GROUP V. n = 5, the prereg forbids it, and
the way that rule survives contact with a tempting result is by the code refusing
to produce the number at all.

DIRECTION, BECAUSE THE TWO TIERS POINT OPPOSITE WAYS. Tier 1 is a displacement:
SMALLER is more stable. Tier 2 is a stability score: LARGER is more stable. Every
comparison is expressed as "more stable", and the sign flip is applied once, in
`STABILITY_SIGN`, rather than in each caller.
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

from shared import bpmd                            # noqa: E402
from shared import outputs as sout                 # noqa: E402
import elevation_run as er                         # noqa: E402

log = logging.getLogger("elevation-analysis")
OUT = sout.Topic("blacksmith", "elevation_analysis")

GROUPS = ["A_hiEnr_hiCons_bdhi", "B_loEnr_hiCons_bdhi", "D_loEnr_loCons_bdhi",
          "V_hiCons_chloroacetamide", er.REF_GROUP]
SHORT = {g: g[0] for g in GROUPS}

# The three contrasts the prereg's readings are written over. REF comparisons are
# reported too -- the anchor exists to be compared against -- but they are marked
# as anchor comparisons, not as one of the three the correction is applied over.
PREREG_CONTRASTS = [("A_hiEnr_hiCons_bdhi", "B_loEnr_hiCons_bdhi"),
                    ("B_loEnr_hiCons_bdhi", "D_loEnr_loCons_bdhi"),
                    ("A_hiEnr_hiCons_bdhi", "D_loEnr_loCons_bdhi")]
ANCHOR_CONTRASTS = [(g, er.REF_GROUP) for g in GROUPS[:4]]

# n = 5. The prereg forbids a significance claim from it, so no p-value is
# computed for any contrast involving it.
NO_INFERENCE = {"V_hiCons_chloroacetamide"}

# Tier 1 is a displacement (smaller = more stable); tier 2 is a stability score
# (larger = more stable). Applied once so no caller has to remember.
STABILITY_SIGN = {"tier1": -1.0, "tier2": +1.0}

DELTA_LARGE = 0.5
ALPHA = 0.05


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def _read(d: Path, pattern: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(d / pattern)),
                key=lambda f: [int(x) for x in re.findall(r"_s?(\d+)", Path(f).stem)])
    if not fs:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def tier1_per_molecule() -> pd.DataFrame:
    """One row per molecule: the warhead's drift across unrestrained equilibration.

    `abs_delta_nm` is the pre-registered readout -- the prereg says "the
    warhead's displacement from its docked position", and a displacement is a
    magnitude. The signed mean is carried beside it because a molecule that moved
    0.2 nm CLOSER is not the same event as one that moved 0.2 nm away, and only
    the sign distinguishes them.
    """
    df = _read(er.T1.dir, "elevation_t1_s*.csv")
    if df.empty:
        raise SystemExit("no tier-1 records")
    df = df.drop_duplicates(["ident", "replicate"], keep="last")
    ok = df[df.status == "ok"].copy()
    rows = []
    for (g, ident), d in ok.groupby(["group", "ident"]):
        rows.append({
            "group": g, "ident": ident,
            "warhead_class": d.warhead_class.iloc[0],
            "n_replicas": len(d),
            "docked_nm": float(d.docked_distance_nm.iloc[0]),
            "tier1": float(d.abs_delta_nm.mean()),
            "tier1_sd": float(d.abs_delta_nm.std(ddof=1)) if len(d) > 1 else 0.0,
            "signed_delta_nm": float(d.delta_nm.mean()),
            "npt_nm": float(d.npt_distance_nm.mean()),
            "npt_nm_sd": float(d.npt_distance_nm.std(ddof=1)) if len(d) > 1 else 0.0,
            "frac_replicas_in_window": float(d.npt_in_window.mean()),
            "any_pbc_wrapped": bool(d.pbc_wrapped.any()),
        })
    return pd.DataFrame(rows), df


def minimisation_split(raw: pd.DataFrame) -> pd.DataFrame:
    """POST-HOC, and labelled as such wherever it is reported.

    The pre-registered tier-1 window runs from the docked pose to the start of
    production, which contains an energy minimisation as well as the 300 ps of
    NVT/NPT. Splitting the two was not registered and does not replace the
    registered readout -- but `min.gro` is already on disk, minimisation is
    shared by every replica of a molecule, and the split answers a question the
    combined number cannot: whether a pose was already relaxing out of the
    pocket before any thermal motion was applied.
    """
    rows = []
    for ident, d in raw[raw.status == "ok"].groupby("ident"):
        r = d.iloc[0]
        gro = er.br.WORK / str(ident).replace(":", "_") / "md" / "min.gro"
        if not gro.is_file():
            continue
        try:
            got = er.distance_nm(gro, int(r.warhead_serial0), int(r.sg_serial0),
                                 (r.warhead_atom_name, r.sg_atom_name))
        except er.ElevationError as exc:
            log.warning("%s: %s", ident, exc)
            continue
        rows.append({"group": r.group, "ident": ident,
                     "docked_nm": float(r.docked_distance_nm),
                     "after_min_nm": round(got["distance_nm"], 4),
                     "min_delta_nm": round(got["distance_nm"] - float(r.docked_distance_nm), 4),
                     "npt_delta_nm": round(float(d.npt_distance_nm.mean())
                                           - got["distance_nm"], 4)})
    return pd.DataFrame(rows)


def tier2_per_molecule() -> pd.DataFrame:
    """One row per molecule: the BPMD stability score and the spread behind it.

    Read from the elevation topic and NOT from `bpmd/`, which still holds ok
    replicates for two of these molecules at 300 ps and 10,000 ps from earlier
    protocol work. Those are not this experiment's protocol, and a between-group
    comparison that mixed them would be measuring its own bookkeeping.
    """
    df = _read(er.T2.dir, "elevation_t2_s*.csv")
    if df.empty:
        return pd.DataFrame(), df
    df = df.drop_duplicates(["ident", "replicate"], keep="last")
    if "group" not in df.columns:
        return pd.DataFrame(), df
    ok = df[(df.status == "ok") & df.group.notna()].copy()
    lens = sorted(ok.production_ps.unique())
    if len(lens) > 1:
        raise SystemExit(
            f"tier-2 records span more than one trajectory length {lens}; a "
            "between-group comparison needs protocol consistency above all else")
    rows = []
    for (g, ident), d in ok.groupby(["group", "ident"]):
        res = [bpmd.ReplicaResult(int(r.replicate), float(r.max_cv_nm),
                                  float(r.frac_in_window),
                                  float(r.bias_at_exit_kj), bool(r.escaped))
               for r in d.itertuples()]
        s = bpmd.combine(res)
        rows.append({
            "group": g, "ident": ident,
            "warhead_class": d.warhead_class.iloc[0],
            "n_replicas": s.n_replicas,
            "tier2": s.score,
            "frac_in_window": s.mean_frac_in_window,
            "tier2_spread": s.spread_frac,
            "median_bias_kj": s.median_bias_to_escape_kj,
            "n_escaped": s.n_escaped,
            "production_ps": float(d.production_ps.iloc[0]),
            "max_cv_nm": float(d.max_cv_nm.max()),
        })
    return pd.DataFrame(rows), df


# --------------------------------------------------------------------------
# statistics
# --------------------------------------------------------------------------

def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """P(x > y) - P(x < y). Identical to the rank-biserial correlation from U.

    Ties are counted as neither, which is what makes it readable on a small
    sample: it is a statement about pairs, and n = 8 vs 8 is 64 pairs.
    """
    gt = sum((a > b) for a in x for b in y)
    lt = sum((a < b) for a in x for b in y)
    return (gt - lt) / (len(x) * len(y))


def compare(df: pd.DataFrame, col: str, g1: str, g2: str, tier: str) -> dict:
    """One group contrast, expressed as "more stable" regardless of tier direction."""
    from scipy import stats

    x = df.loc[df.group == g1, col].to_numpy(float)
    y = df.loc[df.group == g2, col].to_numpy(float)
    row = {"tier": tier, "metric": col, "group_1": SHORT[g1], "group_2": SHORT[g2],
           "n_1": len(x), "n_2": len(y),
           "median_1": float(np.median(x)) if len(x) else np.nan,
           "median_2": float(np.median(y)) if len(y) else np.nan}
    # The effect-size column is present on EVERY row, populated or not. A column
    # that appears only when the comparison succeeded makes every downstream
    # reader conditional on the data, and `readings` then fails on an empty
    # group instead of reporting one.
    row["cliffs_delta_more_stable"] = np.nan
    if len(x) < 2 or len(y) < 2:
        return {**row, "p": np.nan, "note": "too few molecules to compare"}

    s = STABILITY_SIGN[tier]
    row["cliffs_delta_more_stable"] = round(cliffs_delta(s * x, s * y), 3)

    # THE PREREG FORBIDS A SIGNIFICANCE CLAIM FROM n = 5, so the p-value is not
    # computed rather than computed and then asked to be ignored.
    if g1 in NO_INFERENCE or g2 in NO_INFERENCE:
        return {**row, "p": np.nan,
                "note": "descriptive only — the prereg forbids inference at n = 5"}
    u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    return {**row, "U": float(u), "p": float(p), "note": ""}


def holm(ps: list[float]) -> list[float]:
    """Holm-Bonferroni, applied only across the three pre-registered contrasts."""
    idx = [i for i, p in enumerate(ps) if np.isfinite(p)]
    out = [np.nan] * len(ps)
    m, run = len(idx), 0.0
    for rank, i in enumerate(sorted(idx, key=lambda i: ps[i])):
        run = max(run, (m - rank) * ps[i])
        out[i] = min(1.0, run)
    return out


def contrast_table(df: pd.DataFrame, col: str, tier: str) -> pd.DataFrame:
    rows = [{**compare(df, col, a, b, tier), "family": "pre-registered"}
            for a, b in PREREG_CONTRASTS]
    ps = [r.get("p", np.nan) for r in rows]
    for r, p in zip(rows, holm(ps)):
        r["p_holm"] = p
        r["verdict"] = (">" if (np.isfinite(p) and p < ALPHA
                                and abs(r.get("cliffs_delta_more_stable", 0)) >= DELTA_LARGE)
                        else "~")
    for a, b in ANCHOR_CONTRASTS:
        r = compare(df, col, a, b, tier)
        r.update({"family": "anchor", "p_holm": np.nan,
                  "verdict": ("" if not np.isfinite(r.get("p", np.nan))
                              else (">" if (r["p"] < ALPHA
                                            and abs(r.get("cliffs_delta_more_stable", 0)) >= DELTA_LARGE)
                                    else "~"))})
        rows.append(r)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# the pre-registered readings
# --------------------------------------------------------------------------

def readings(con: pd.DataFrame, per: pd.DataFrame, col: str, tier: str) -> dict:
    """Which of the prereg's four mutually-exclusive patterns the data shows."""
    def v(a, b):
        r = con[(con.group_1 == a) & (con.group_2 == b) & (con.family == "pre-registered")]
        if r.empty:
            return "~", 0.0
        d = float(r.cliffs_delta_more_stable.iloc[0])
        # A contrast that could not be computed is "~" with no effect, not a
        # NaN propagating into a direction test that would then read as False
        # for both directions and silently land on "A ~ B ~ D".
        return r.verdict.iloc[0], (d if np.isfinite(d) else 0.0)

    # AN INCOMPLETE RUN MUST NOT PRODUCE A READING. With a group still empty
    # every contrast is "~", which lands on "A ~ B ~ D" -- "neither metric
    # predicts stability" -- and that is a substantive conclusion the prereg
    # attaches a real consequence to. Reporting it because the jobs had not
    # finished would be the most expensive failure available here, so it is
    # refused by name rather than caveated in prose.
    thin = {SHORT[g]: int((per.group == g).sum()) for g in GROUPS[:3]
            if (per.group == g).sum() < 2}
    if thin:
        return {"tier": tier, "metric": col,
                "reading": "NOT READABLE — the run is incomplete",
                "conclusion": f"groups {thin} have too few molecules with "
                              "results; no pre-registered reading is emitted",
                "A_vs_B": "", "B_vs_D": "", "A_vs_D": "",
                "delta_A_B": np.nan, "delta_B_D": np.nan, "delta_A_D": np.nan,
                "V_median": np.nan, "REF_median": np.nan,
                "V_vs_REF_delta": np.nan, "V_vs_REF": ""}

    ab, d_ab = v("A", "B")
    bd, d_bd = v("B", "D")
    ad, d_ad = v("A", "D")

    # "A > D" / "B > D" require the direction as well as the magnitude: a
    # significant contrast in the WRONG direction is the prereg's fourth reading,
    # not its first.
    a_beats_d = ad == ">" and d_ad > 0
    b_beats_d = bd == ">" and d_bd > 0
    d_beats = (ad == ">" and d_ad < 0) or (bd == ">" and d_bd < 0)
    a_beats_b = ab == ">" and d_ab > 0

    if d_beats:
        reading, conclusion = (
            "D >= A, B",
            "Something is wrong with the design or the metric direction. The "
            "prereg says: report as a failure, do not reinterpret.")
    elif a_beats_b and (a_beats_d or b_beats_d):
        reading, conclusion = (
            "A > B, both > D",
            "Enrichment adds something real beyond consensus. Both belong in "
            "the ranking.")
    elif (not a_beats_b) and (a_beats_d and b_beats_d):
        reading, conclusion = (
            "B ~ A, both > D",
            "Consensus is the filter; enrichment adds nothing. The shortlist "
            "should be drawn on consensus, and the 397 are the real candidate "
            "pool.")
    else:
        reading, conclusion = (
            "A ~ B ~ D",
            "Neither metric predicts stability. This tier is measuring "
            "something orthogonal, and the ranking has no physical support "
            "from this experiment.")

    # V ~ REF is descriptive by construction: n = 5.
    def med(g):
        s = per.loc[per.group == g, col]
        return float(s.median()) if len(s) else np.nan
    mv, mr = med("V_hiCons_chloroacetamide"), med(er.REF_GROUP)
    dv = con[(con.group_1 == "V") & (con.group_2 == "R")]
    v_ref = (float(dv.cliffs_delta_more_stable.iloc[0])
             if len(dv) and np.isfinite(dv.cliffs_delta_more_stable.iloc[0]) else np.nan)

    return {"tier": tier, "metric": col, "reading": reading,
            "conclusion": conclusion,
            "A_vs_B": ab, "B_vs_D": bd, "A_vs_D": ad,
            "delta_A_B": d_ab, "delta_B_D": d_bd, "delta_A_D": d_ad,
            "V_median": mv, "REF_median": mr, "V_vs_REF_delta": v_ref,
            "V_vs_REF": "descriptive only (n = 5; no significance claim)"}


# --------------------------------------------------------------------------

def group_summary(per: pd.DataFrame, col: str, sd_col: str | None) -> pd.DataFrame:
    rows = []
    for g in GROUPS:
        d = per[per.group == g]
        if d.empty:
            continue
        s = d[col]
        r = {"group": g, "n_molecules": len(d),
             "median": float(s.median()),
             "q1": float(s.quantile(0.25)), "q3": float(s.quantile(0.75)),
             "min": float(s.min()), "max": float(s.max()),
             "mean_replicas": float(d.n_replicas.mean())}
        if sd_col and sd_col in d:
            r["median_replica_spread"] = float(d[sd_col].median())
        rows.append(r)
    return pd.DataFrame(rows)


def show(title: str, df: pd.DataFrame) -> None:
    print(f"\n=== {title} ===")
    print(df.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--write", action="store_true", help="persist the tables")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    t1, t1_raw = tier1_per_molecule()
    show("tier 1 — |delta d| across 300 ps of unrestrained equilibration (nm; "
         "SMALLER is more stable)", group_summary(t1, "tier1", "tier1_sd"))
    c1 = contrast_table(t1, "tier1", "tier1")
    show("tier 1 contrasts", c1)
    r1 = readings(c1, t1, "tier1", "tier1")
    print(f"\n  TIER-1 READING: {r1['reading']}\n  {r1['conclusion']}")

    frames = {"tier1_per_molecule": t1, "tier1_contrasts": c1,
              "tier1_replicates": t1_raw}
    reads = [r1]

    ms = minimisation_split(t1_raw)
    if not ms.empty:
        show("POST-HOC (not pre-registered) — where the tier-1 drift happened: "
             "energy minimisation vs the 300 ps of dynamics (nm)",
             ms.groupby("group")[["min_delta_nm", "npt_delta_nm"]]
             .median().round(3).reset_index())
        frames["tier1_minimisation_split"] = ms

    t2, t2_raw = tier2_per_molecule()
    if t2.empty:
        print("\n  tier 2: nothing recorded yet")
    else:
        show("tier 2 — BPMD stability score (LARGER is more stable)",
             group_summary(t2, "tier2", "tier2_spread"))
        show("tier 2 — mean fraction of time in the near-attack window",
             group_summary(t2, "frac_in_window", "tier2_spread"))
        c2 = contrast_table(t2, "tier2", "tier2")
        show("tier 2 contrasts", c2)
        r2 = readings(c2, t2, "tier2", "tier2")
        print(f"\n  TIER-2 READING: {r2['reading']}\n  {r2['conclusion']}")
        frames.update({"tier2_per_molecule": t2, "tier2_contrasts": c2,
                       "tier2_replicates": t2_raw})
        reads.append(r2)

    show("pre-registered readings", pd.DataFrame(reads)[
        ["tier", "reading", "A_vs_B", "B_vs_D", "A_vs_D", "V_median",
         "REF_median", "V_vs_REF_delta"]])

    if args.write:
        for name, df in frames.items():
            dest = OUT.write(f"elevation_{name}", ".csv")
            df.to_csv(dest, index=False)
            log.info("%-24s %3d rows -> %s", name, len(df), dest.name)
        dest = OUT.write("elevation_readings", ".csv")
        pd.DataFrame(reads).to_csv(dest, index=False)
        log.info("readings -> %s", dest.name)


if __name__ == "__main__":
    main()
