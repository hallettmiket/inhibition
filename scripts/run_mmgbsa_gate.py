"""
Purpose: Put MM-GBSA through the SAME enrichment gate that docking failed.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the class-matched gate poses under 00_shared_substrate/m3_enrichment/covalent
Output: an `mmgbsa_dG` verdict merged into the enrichment-gate token

Run:  python scripts/run_mmgbsa_gate.py [--workers 6] [--limit N]

WHY THIS IS THE RIGHT NEXT EXPERIMENT. D0031 measured docking against
class-matched decoys and found it indistinguishable from chance: covalent
ROC-AUC 0.537, non-covalent 0.535, EF1% 0.0 for both. The build therefore has no
working discriminator, and its four shortlists are orderings nothing validates.

The instinct is to treat MM-GBSA as a refinement of the docking ranking — a
higher rung on a fidelity ladder whose lower rung selects what climbs it. But if
the lower rung is at chance, refinement is the wrong frame entirely. MM-GBSA is
an INDEPENDENT estimator, and the question worth asking is not "does it polish
the shortlist" but **"can it do what docking could not"**.

That question is testable with what already exists: the same actives, the same
class-matched decoys (D0031), the same graded gate. Nothing new has to be
measured or assumed.

EITHER ANSWER IS WORTH HAVING. If MM-GBSA separates actives from same-chemotype
decoys, the build acquires the validated ranking it currently lacks. If it does
not, that is a stronger negative result than docking's alone — two independent
methods failing on the same shallow, solvent-exposed pocket says something about
the pocket rather than about one scoring function.

SCOPE, STATED UP FRONT. The covalent junction parameters cover an sp3 attachment
carbon only; an sp2 attachment fails with a missing `2C - S - cc` angle. So the
naphthoquinone chemotype cannot be scored until those parameters are extended,
and this gate covers the chemotypes that build. Which ones those are is reported
rather than quietly dropped, because a gate silently narrowed to one chemotype
and a gate designed for one are very different claims.

DIRECTION. dG is kcal/mol and LOWER is better, like every other rank metric in
this project. Recorded explicitly so it cannot be inverted by assumption.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import covalent_adduct as cad          # noqa: E402
from shared import enrichment_gate as eg           # noqa: E402
from shared import mmgbsa as mg                    # noqa: E402
from shared import warhead_library as wl           # noqa: E402

log = logging.getLogger("mmgbsa-gate")

GATE_DIR = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"
                "/m3_enrichment/covalent")
RESULTS = GATE_DIR / "results_adduct_classmatched.jsonl"
WORK_ROOT = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"
                 "/mmgbsa_gate")
NICE = 19

# One thread per worker. antechamber's AM1-BCC step shells out to sqm, which is
# OpenMP-threaded and will otherwise take every core on the machine.
_SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                  "OPENBLAS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}


def score_one(job: dict) -> dict:
    """MM-GBSA on one gate ligand, from the pose the gate already docked."""
    os.nice(NICE)
    os.environ.update(_SINGLE_THREAD)
    wd = WORK_ROOT / job["id"]
    wd.mkdir(parents=True, exist_ok=True)

    done = wd / "result.json"
    if done.is_file():
        return json.loads(done.read_text())

    try:
        pose = GATE_DIR / f"{job['id']}_docked.sdf"
        if not pose.is_file():
            raise mg.MMGBSAError(f"no docked pose at {pose}")

        lib = wl.load()
        smarts = cad.adduct_attachment_smarts(job["warhead_class"], library=lib)
        cyx, cys, cyx_idx, n_res = mg.prepare_receptor(wd)
        mol2, frcmod, att, cap, q = mg.parameterize_ligand(
            pose, wd, smarts, net_charge=0)
        legs = mg.build_topologies(wd, mol2, frcmod, cyx, cys, cyx_idx,
                                   n_res + 1, att, cap, q)
        mg.verify_complex(legs["complex"][0], cyx_idx, att)
        energies = {leg: mg.minimize_and_score(wd, leg)
                    for leg in ("complex", "receptor", "ligand")}
        out = {**job, **mg.delta_g(energies)}
        out.pop("per_term", None)
        done.write_text(json.dumps(out, indent=2), encoding="utf-8")
        return out
    except Exception as exc:  # noqa: BLE001 - one failure must not end the run
        (wd / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return {**job, "mmgbsa_error": str(exc)[:300]}


def load_jobs() -> list[dict]:
    """The class-matched gate ligands, from the run that produced D0031."""
    if not RESULTS.is_file():
        raise SystemExit(f"no class-matched gate results at {RESULTS}; run "
                         "`scripts/run_enrichment_gate.py covalent` first")
    jobs = []
    for line in RESULTS.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        if r.get("error") or r.get("affinity_kcal") is None:
            continue
        jobs.append({"id": r["id"], "name": r.get("name"),
                     "label": int(r["label"]),
                     "warhead_class": r["warhead_class"],
                     "canonical_smiles": r.get("adduct_smiles") or r.get("smiles"),
                     "affinity_kcal": r.get("affinity_kcal")})
    return jobs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    WORK_ROOT.mkdir(parents=True, exist_ok=True)

    jobs = load_jobs()
    if args.limit:
        jobs = jobs[:args.limit]
    n_act = sum(1 for j in jobs if j["label"] == 1)
    log.info("MM-GBSA gate on %d ligands (%d actives, %d decoys), %d workers",
             len(jobs), n_act, len(jobs) - n_act, args.workers)

    scored, failed = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(score_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            (failed if "mmgbsa_error" in r else scored).append(r)
            if i % 10 == 0 or "mmgbsa_error" in r:
                log.info("[%d/%d] %s %s", i, len(jobs), r["id"],
                         f"dG {r['dG_kcal']:.2f}" if "dG_kcal" in r
                         else f"FAILED {r['mmgbsa_error'][:70]}")

    if not scored:
        raise SystemExit("no ligand produced a dG; nothing to grade")

    df = pd.DataFrame(scored)
    # WHICH CHEMOTYPES SURVIVED IS PART OF THE RESULT, not a footnote.
    by_class = df.groupby("warhead_class")["label"].agg(["size", "sum"])
    fail_by_class = (pd.DataFrame(failed).groupby("warhead_class").size()
                     if failed else pd.Series(dtype=int))
    print("\n=== coverage ===")
    print(f"{'chemotype':24s} {'scored':>7s} {'actives':>8s} {'failed':>7s}")
    for cls in sorted(set(by_class.index) | set(fail_by_class.index)):
        print(f"{cls:24s} {int(by_class['size'].get(cls, 0)):7d} "
              f"{int(by_class['sum'].get(cls, 0)):8d} "
              f"{int(fail_by_class.get(cls, 0)):7d}")

    usable = df[df["label"].isin([0, 1])].copy()
    act = usable[usable.label == 1]
    if act.empty:
        raise SystemExit(
            "every ACTIVE failed to parameterise, so there is nothing to "
            "enrich for. This is a coverage failure, not a negative result.")

    res = eg.evaluate(usable, metric="dG_kcal", stratum="covalent",
                      higher_is_better=False)
    print(f"\n[covalent/mmgbsa_dG] {res.verdict}")
    print(f"  ROC-AUC {res.roc_auc:.3f}  CI[{res.roc_auc_ci[0]:.3f},"
          f"{res.roc_auc_ci[1]:.3f}]  EF1% {res.ef_1pct:.1f}  "
          f"BEDROC {res.bedroc:.3f}")
    print(f"  {res.n_actives} actives / {res.n_decoys} decoys / "
          f"{res.n_chemotypes} chemotypes")
    for r in res.reasons:
        print(f"  - {r}")

    # The comparison that motivated the run: same ligands, same gate, two methods.
    dock = eg.evaluate(usable, metric="affinity_kcal", stratum="covalent",
                       higher_is_better=False)
    print(f"\n  docking on the SAME {len(usable)} ligands: "
          f"ROC-AUC {dock.roc_auc:.3f}, EF1% {dock.ef_1pct:.1f}")
    delta = res.roc_auc - dock.roc_auc
    print(f"  MM-GBSA - docking = {delta:+.3f} ROC-AUC")

    res.metric = "mmgbsa_dG"
    eg.write_token([res])
    log.info("merged mmgbsa_dG into the enrichment-gate token")

    if failed:
        print(f"\n  {len(failed)} ligand(s) could not be parameterised. The "
              "junction covers an sp3 attachment carbon only; an sp2 attachment "
              "needs the frcmod extended (see D0027's method).")


if __name__ == "__main__":
    main()
