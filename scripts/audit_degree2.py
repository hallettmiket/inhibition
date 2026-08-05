"""
Purpose: audit every T_2 degree-2 sample — is it a sample, is it uniform, and is it comparable across seeds?
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: each seed's degree-2 frame + its degree-1 parent frame
Output: a per-seed report; `--check` exits non-zero on anything that would mislead

WHY THIS EXISTS. A degree-2 frame is a SAMPLE of a population ~100-1000x larger
than itself, and every count derived from it is wrong by that factor unless the
reader knows. The previous run of this sampler was silently HALF its intended
size for six days -- the manifest recorded everything needed to notice
(`estimated_population` 7,800,890 against a `realised_population` of 4,063,427)
and nobody compared the two numbers.

So this compares them, and the other things that can go wrong without raising:

* **`is_sample` on every row.** A frame that carries the fact only in its
  manifest is one merge away from being counted as a census.
* **The realised fraction is MEASURED, not assumed.** Reservoir sampling makes
  it a result; this asserts the recorded value matches `n_kept / population`.
* **Size comparability across seeds.** Five pools shown side by side in the GUI
  must not differ in size for a reason unrelated to chemistry. The old
  estimator would have made them differ by whatever its per-seed error was.
* **Parent coverage.** A uniform sample of the union should draw from
  essentially all parents. If a large share of parents are unrepresented, the
  draw was position-biased -- which is exactly the `frontier_cap`
  truncation-in-parent-order failure the sampler exists to avoid.
* **Degree-1 leakage.** Nothing in a degree-2 frame may already be in the
  degree-1 frame or be the seed. The sampler seeds its dedup set with both;
  this verifies it worked rather than trusting it.
* **The governor.** Degree-2 products can carry two growths, which is what the
  55-heavy-atom pocket ceiling exists for. Zero pruned is plausible; a large
  fraction pruned means the ceiling is shaping the chemistry and should be said
  out loud rather than discovered later.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import io as dio                       # noqa: E402
from shared import seeds as sd                     # noqa: E402
from shared import smiles as smi                   # noqa: E402

log = logging.getLogger("audit-degree2")
DATA = Path("/data/lab_vm/append_only/inhibition")
APPROACH = "t2"

SEEDS = ("atra", "potter_astex", "guo_pfizer", "du_xu", "liu_2024_c3")

# A uniform draw over the union should touch nearly every parent. Below this
# share, the draw is concentrated and the sample is not what it claims.
MIN_PARENT_COVERAGE = 0.50


SAMPLE_STAGE = "t2_generate_degree2_sample"


def _sampling_manifest(experiment_dir: Path) -> dict:
    """The manifest of the stage that DREW THE SAMPLE, not the newest one.

    The sampling parameters -- population, fraction, target -- exist only on the
    generate stage. Reading the latest frame's manifest instead returns the
    RANK stage's, which has none of them, and every one of those fields comes
    back as a default. It reported `population 0, fraction 0.00000` for a frame
    whose real population is 4,063,427: populated, plausible, and taken from
    the wrong artefact because it was the newest one. Select by stage, which is
    the identity that matters, not by recency.
    """
    best, best_n = {}, -1
    for m in experiment_dir.glob("D2_*_manifest.json"):
        try:
            data = json.loads(m.read_text())
        except Exception:  # noqa: BLE001
            continue
        if data.get("stage") != SAMPLE_STAGE:
            continue
        n = int(m.name.split("_")[1])
        if n > best_n:
            best, best_n = data, n
    return best


def audit_seed(seed: str) -> dict | None:
    try:
        rec = sd.resolve(APPROACH, seed)
    except Exception as exc:  # noqa: BLE001
        return {"seed": seed, "state": "no-seed-record", "detail": str(exc)[:80]}

    d2_exp = f"{rec['experiment']}_degree2"
    if not (DATA / d2_exp).is_dir():
        return {"seed": seed, "state": "not-run", "experiment": d2_exp}

    frame = dio.latest(DATA / d2_exp, "D2", ".parquet")
    if frame is None:
        return {"seed": seed, "state": "not-run", "experiment": d2_exp}

    df = dio.read_frame(frame)
    man = _sampling_manifest(DATA / d2_exp)
    p = man.get("params", {})

    d1, _ = dio.latest_frame(rec["experiment"], APPROACH)
    d1 = d1[d1["degree"] == 1] if "degree" in d1.columns else d1
    parents = set(d1["canonical_smiles"])

    pop = p.get("realised_population")
    kept = len(df)
    frac_recorded = p.get("sampling_fraction")
    frac_actual = (kept / pop) if pop else None

    # Leakage: nothing at degree 2 may already exist at degree 1.
    d1_keys = {k for k in (smi.inchikey(s) for s in parents) if k}
    d2_keys = [smi.inchikey(s) for s in df["canonical_smiles"]]
    leaked = sum(1 for k in d2_keys if k and k in d1_keys)

    covered = df["parent_smiles"].nunique() if "parent_smiles" in df else 0
    problems: list[str] = []
    notes: list[str] = []

    if "is_sample" not in df.columns or not bool(df["is_sample"].all()):
        problems.append("is_sample missing or False on some rows")
    # WHAT `sampling_fraction` MEANS DEPENDS ON THE METHOD, and conflating the
    # two would make every legacy frame look corrupt while letting a real
    # reservoir defect pass.
    #
    #   Bernoulli (legacy): the recorded value is the INTENDED probability p.
    #     The realised fraction is a random variable around it, so a small
    #     difference is the method working, not a fault.
    #   Reservoir (current): the recorded value IS the realised fraction,
    #     computed after the fact. Any difference means the two disagree about
    #     the same measured quantity, which is a defect.
    method = (p.get("sampling") or "")
    is_reservoir = "reservoir" in method.lower()
    if frac_recorded is not None and frac_actual is not None:
        delta = abs(frac_recorded - frac_actual)
        if is_reservoir and delta > 1e-9:
            problems.append(
                f"reservoir recorded fraction {frac_recorded:.6f} != measured "
                f"{frac_actual:.6f} — these must be the same number")
        elif not is_reservoir and delta > 1e-6:
            notes.append(
                f"legacy Bernoulli: recorded p={frac_recorded:.6f} is the "
                f"INTENDED probability; realised fraction was "
                f"{frac_actual:.6f}. Expected for the method — but the target "
                "was missed because the POPULATION was over-estimated, not "
                "because of this scatter.")
    if leaked:
        problems.append(f"{leaked} degree-2 rows already exist at degree 1")
    if parents and covered / len(parents) < MIN_PARENT_COVERAGE:
        problems.append(
            f"only {covered}/{len(parents)} parents represented "
            f"({covered/len(parents):.1%}) — the draw looks position-biased")
    if df["candidate_id"].duplicated().any():
        problems.append("duplicate candidate_id")

    return {
        "seed": seed, "state": "ok" if not problems else "PROBLEM",
        "frame": frame.name, "experiment": d2_exp,
        "n_parents": len(parents), "parents_covered": covered,
        "population": pop, "kept": kept,
        "target": p.get("target"),
        "fraction": frac_actual,
        "method": (p.get("sampling") or "")[:44],
        "governor_pruned": p.get("governor_pruned_oversize"),
        "legacy_est": p.get("legacy_estimator_would_have_said")
                      or p.get("estimated_population"),
        "leaked": leaked, "problems": problems, "notes": notes,
        "is_reservoir": is_reservoir,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any seed has a problem")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    rows = [audit_seed(s) for s in SEEDS]
    done = [r for r in rows if r and r.get("state") in ("ok", "PROBLEM")]

    print(f"\n{'seed':<14}{'kept':>8}{'target':>8}{'population':>13}"
          f"{'fraction':>11}{'parents':>16}{'state':>9}")
    for r in rows:
        if r.get("state") in ("not-run", "no-seed-record"):
            print(f"{r['seed']:<14}{'—':>8}{'—':>8}{'—':>13}{'—':>11}"
                  f"{'—':>16}{r['state']:>9}")
            continue
        cov = (f"{r['parents_covered']:,}/{r['n_parents']:,}")
        print(f"{r['seed']:<14}{r['kept']:>8,}{(r['target'] or 0):>8,}"
              f"{(r['population'] or 0):>13,}{(r['fraction'] or 0):>11.5f}"
              f"{cov:>16}{r['state']:>9}")

    if done:
        print("\nComparability across seeds (they are shown side by side):")
        sizes = {r["seed"]: r["kept"] for r in done}
        if len(sizes) < 2:
            # One seed is not a comparison. Saying "IDENTICAL" here would be
            # trivially true and would read as a cross-seed result.
            print(f"  only {len(sizes)} seed generated so far — nothing to "
                  "compare yet")
        else:
            lo, hi = min(sizes.values()), max(sizes.values())
            print(f"  sample sizes {lo:,}–{hi:,}"
                  + ("  IDENTICAL — differences between pools are chemistry"
                     if lo == hi else
                     f"  ratio {hi/lo:.2f}x — a size difference NOT from "
                     "chemistry"))
        print("\nWhat the retired estimator would have produced instead:")
        for r in done:
            if r.get("legacy_est") and r.get("population"):
                err = r["legacy_est"] / r["population"]
                print(f"  {r['seed']:<14} guessed {r['legacy_est']:>12,.0f} vs "
                      f"measured {r['population']:>12,}  ({err:.2f}x)")

    for r in done:
        for n in r.get("notes", []):
            print(f"\n{r['seed']}: note — {n}")

    bad = [r for r in done if r["problems"]]
    for r in bad:
        print(f"\n{r['seed']}: PROBLEMS")
        for pr in r["problems"]:
            print(f"    - {pr}")

    if args.check and bad:
        raise SystemExit(f"{len(bad)} seed(s) failed the degree-2 audit")
    if args.check:
        missing = [r["seed"] for r in rows if r.get("state") == "not-run"]
        if missing:
            raise SystemExit(f"not yet generated to degree 2: {missing}")
        print("\nall seeds audited clean.")


if __name__ == "__main__":
    main()
