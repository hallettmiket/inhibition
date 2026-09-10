"""
Purpose: how much did each CReM edit move a T_2 degree-2 molecule -- structurally, and in its score?
Author: @tt8804 (with Claude Code)
Date: 2026-09-10
Input: each degree-2 pool's ranked frame (parent_smiles + ph2d_max_similarity
       already on it, from D0117) + its distinct parents
Output: a summary JSON under append_only, printed to stdout

WHY "VS PARENT" IS A DIFFERENT QUESTION FROM D0117. D0117 asked whether a
molecule resembles a KNOWN Pin1 binder. This asks whether one more CReM edit
moved the molecule TOWARD that resemblance or away from it -- the "does a
second edit help" question `integration/app/data.py` already names but never
answers. It needs the PARENT's own score under the identical held-out model
set, which nothing before this computed -- only the degree-2 CHILDREN were
ever scored.

NOT THE SAME COMPUTATION AS "NOVELTY", AND NOT A SUBSTITUTE FOR IT.
`docs/shared/gates.md` control B4 forbids computing NOVELTY against the seed
-- seed-relative novelty mechanically rewards T_1 (no seed) and penalises T_2
for doing what it was asked. This is not that: the seed stays excluded from
the MODEL SET exactly as D0117 left it (`hold_out=seed_smiles`, unchanged),
and what is compared here is the immediate one-edit-away PARENT, not the seed.
`novelty_external` (vs the frozen reference set) is untouched by this script
and remains the only axis anything is ranked on.

TWO FINGERPRINTS, EACH REUSED FROM WHERE THIS PROJECT ALREADY USES IT, NOT A
THIRD CHOICE INVENTED HERE:

  * STRUCTURAL similarity to parent -- ECFP4 (Morgan r=2), the SAME
    definition `t2_generate_degree2_sample`'s own manifest already names for
    novelty ("1 - max Tanimoto (ECFP4) vs the external set") and
    `enrichment_gate` uses for chemotype clustering. Answers "how big was
    the edit".
  * The SCORE delta -- the Gobbi 2D pharmacophore Tanimoto D0117 measured
    and gated, applied to the parent under the identical held-out model set
    the child was scored against. Answers "did the edit move the molecule
    toward known actives".

SCOPE, STATED RATHER THAN SILENTLY WIDENED. The two degree-2 pools ranked in
D0117 (60,000 molecules), each versus its own immediate parent (already on
the frame as `parent_smiles`). NOT extended to the other four T_2 seeds'
degree-1 pools, which have never been pharmacophore-scored at all -- that is
a bigger, separate decision (scoring five more pools against their own
held-out model sets) and is not silently taken here.

THIS IS DIAGNOSTIC, NOT A RANK METRIC. Neither column is registered in
`rank_shortlist.RANK_DIRECTION` and neither frame is re-versioned: ranking
candidates by closeness to their own parent would be backwards (CReM's whole
point is producing something different), and this is written to
`00_outputs/`, not merged onto the candidate frame, so it cannot be mistaken
for one.
"""

from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator as fpg

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
RDLogger.DisableLog("rdApp.*")

from shared import compute as cmp                    # noqa: E402
from shared import io as dio                          # noqa: E402
from shared import outputs as out                     # noqa: E402
from shared import pharmacophore as ph                # noqa: E402
from shared.manifest import Manifest                  # noqa: E402
from scripts.pharmacophore_gate import score_parallel  # noqa: E402

log = logging.getLogger("pharmacophore-parent-drift")

DATA = Path("/data/lab_vm/append_only/inhibition")
AGENT = "blacksmith"
TOPIC = "pharmacophore_gate"

# The exact ECFP4 definition this project already uses for structural
# similarity -- not a new choice made here.
_ECFP4 = fpg.GetMorganGenerator(radius=2, fpSize=2048)

POOLS = {
    "atra": "02_t2_atra_crem_degree2",
    "guo_pfizer": "02e_t2_guo_crem_degree2",
}


def _ecfp4(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    return _ECFP4.GetFingerprint(mol) if mol else None


def run_pool(seed_name: str, experiment: str, workers: int) -> tuple[dict, pd.DataFrame]:
    t0 = time.time()
    frame_path = dio.latest(DATA / experiment, "D2", ".parquet")
    df = dio.read_frame(frame_path)
    for col in ("ph2d_max_similarity", "parent_smiles"):
        if col not in df.columns:
            raise SystemExit(f"{frame_path} has no {col!r} -- run "
                             "scripts/rank_pharmacophore.py first")
    log.info("[%s] %s: %d rows", seed_name, frame_path.name, len(df))

    # --- structural similarity to the immediate parent (ECFP4) -----------
    uniq = pd.Index(df.canonical_smiles).append(pd.Index(df.parent_smiles)).unique()
    fps = {s: _ecfp4(s) for s in uniq}
    n_bad = sum(v is None for v in fps.values())
    if n_bad:
        log.warning("[%s] %d of %d unique SMILES did not parse for ECFP4",
                    seed_name, n_bad, len(uniq))
    sim_to_parent = np.array([
        DataStructs.TanimotoSimilarity(fps[c], fps[p])
        if fps.get(c) is not None and fps.get(p) is not None else np.nan
        for c, p in zip(df.canonical_smiles, df.parent_smiles)])
    df = df.assign(ecfp4_similarity_to_parent=sim_to_parent)

    # --- the parent's OWN pharmacophore score, same held-out model set ---
    seed_smiles = yaml.safe_load(
        (REPO / "config" / "seeds.yaml").read_text())["seeds"][seed_name][
        "canonical_smiles"]
    model = ph.build_model_set(hold_out=seed_smiles)
    distinct_parents = df.parent_smiles.dropna().unique().tolist()
    parent_scored = score_parallel(distinct_parents, model, workers)
    parent_lookup = pd.Series(parent_scored[ph.METRIC].values,
                              index=distinct_parents)
    df["parent_ph2d_max_similarity"] = df.parent_smiles.map(parent_lookup)
    df["ph2d_delta_vs_parent"] = (
        df[ph.METRIC].astype(float) - df["parent_ph2d_max_similarity"].astype(float))

    elapsed = time.time() - t0
    delta = df.ph2d_delta_vs_parent
    n_scored = int(delta.notna().sum())
    n_improved = int((delta > 0).sum())
    n_worsened = int((delta < 0).sum())
    summary = {
        "seed_name": seed_name, "experiment": experiment,
        "frame": frame_path.name, "n_rows": len(df),
        "n_distinct_parents": len(distinct_parents),
        "n_unique_smiles_fingerprinted": int(len(uniq)),
        "elapsed_s": round(elapsed, 1),
        "ecfp4_similarity_to_parent": {
            "median": float(np.nanmedian(sim_to_parent)),
            "mean": float(np.nanmean(sim_to_parent)),
            "p10": float(np.nanpercentile(sim_to_parent, 10)),
            "p90": float(np.nanpercentile(sim_to_parent, 90)),
        },
        "ph2d_delta_vs_parent": {
            "median": float(delta.median()), "mean": float(delta.mean()),
            "n_scored": n_scored, "n_improved": n_improved,
            "n_worsened": n_worsened,
            "frac_improved": round(n_improved / n_scored, 4) if n_scored else None,
        },
    }
    log.info("[%s] done in %.1fs — median ECFP4-to-parent %.3f, "
             "median score delta %+.4f, %d/%d improved (%d worsened)",
             seed_name, elapsed, summary["ecfp4_similarity_to_parent"]["median"],
             summary["ph2d_delta_vs_parent"]["median"], n_improved, n_scored,
             n_worsened)
    return summary, df


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    workers = cmp.available_workers(None)
    log.info("using %d CPU workers", workers)
    t0 = time.time()
    summaries = {}
    for seed_name, experiment in POOLS.items():
        summaries[seed_name], _ = run_pool(seed_name, experiment, workers)
    total = time.time() - t0

    print("\n" + "=" * 78)
    print("EACH MOLECULE vs ITS PARENT -- structural drift + score delta")
    print("=" * 78)
    for name, s in summaries.items():
        print(f"\n[{name}]  {s['n_rows']:,} molecules, "
              f"{s['n_distinct_parents']:,} distinct parents, {s['elapsed_s']}s")
        e = s["ecfp4_similarity_to_parent"]
        print(f"  ECFP4 similarity to parent : median {e['median']:.3f} "
              f"(p10 {e['p10']:.3f}, p90 {e['p90']:.3f})")
        d = s["ph2d_delta_vs_parent"]
        fi = 100 * d["frac_improved"] if d["frac_improved"] is not None else float("nan")
        print(f"  ph2d score delta vs parent : median {d['median']:+.4f}  "
              f"({d['n_improved']}/{d['n_scored']} improved = {fi:.1f}%, "
              f"{d['n_worsened']} worsened)")
    print(f"\ntotal wall-clock: {total:.1f}s at {workers} workers "
          f"(no GPU used)")

    p = out.write_path(AGENT, TOPIC, "pharmacophore_parent_drift", ".json")
    p.write_text(json.dumps(summaries, indent=2), encoding="utf-8")
    man = Manifest(
        stage="pharmacophore_parent_drift", approach="t2",
        params={"pools": list(POOLS), "total_elapsed_s": round(total, 1),
                "workers": workers,
                "note": "diagnostic only -- not a rank metric, not merged "
                        "onto the candidate frame"})
    man.add_output("summary", p)
    man.write(p.parent, filename=f"{p.stem}_manifest.json")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
