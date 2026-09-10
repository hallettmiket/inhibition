"""
Purpose: run the pharmacophore scorer through the enrichment gate on MEASURED inactives.
Author: @tt8804 (with Claude Code)
Date: 2026-09-09
Input: AID 504891 actives + inactives; the frozen reference binder set
Output: a GateResult per model set, written to append_only + logged

WHY THIS RUNS BEFORE ANY RANKING. #4 Phase 2.2 asks for an orthogonal scorer
and says to MEASURE IT FIRST. Ranking 60,000 molecules on an unmeasured scorer
produces another ordering stamped `rank_validated = False`, which this project
already has ~53,500 of. The gate costs CPU minutes on data already on disk.

WHY THIS IS NOT A FORK OF `run_enrichment_gate.py`. That script DOCKS actives
and decoys through Vina and gnina on up to six GPUs -- it is the gate for a
docking score. The reusable part is `enrichment_gate.evaluate()`, which is
scorer-agnostic: it takes any metric column with its direction and does the
bootstrap CI, the leave-one-chemotype-out and the graded verdict. Calling it
directly is the whole integration.

THE TRAIN/TEST SPLIT IS REAL, AND IT WAS NOT ARRANGED. The model set is the
published Pin1 binders; the gate actives are AID 504891's qHTS hits.
`ingest_measured_inactives` recorded "reference binders also called ACTIVE
here: 0", so the two sets are disjoint as a matter of record rather than by
construction here. That is what makes this a recovery test rather than a
restatement: nothing in the model set can be found by looking for itself.

THE NEGATIVES WERE ASSAYED, NOT ASSUMED. 361,354 measured inactives, against
the property-matched decoys D0041's docking gate had to build. A scorer that
enriches against assayed negatives has cleared a bar docking never did here.

TWO MODEL SETS ARE SCORED, AND THE SECOND IS A SENSITIVITY CHECK, NOT A
RETRY. The declared rule (see `shared/pharmacophore`) excludes a molecule with
4+ amide-SMARTS matches as peptide-like. PiB is a bis-IMIDE: each imide
nitrogen matches twice, so it is excluded for a reason that names the wrong
chemistry. The rule was written down before it was applied and is NOT being
edited to taste (D0045's discipline); instead both model sets are scored and
both verdicts reported. If the verdict moves, that is a finding about the
rule's sensitivity and belongs in the record. If it does not, the point is
closed honestly and cheaply.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import compute as cmp                    # noqa: E402
from shared import enrichment_gate as eg             # noqa: E402
from shared import outputs as out                    # noqa: E402
from shared import pharmacophore as ph               # noqa: E402
from shared.manifest import Manifest                 # noqa: E402

log = logging.getLogger("pharmacophore-gate")

MEASURED = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith"
                "/measured_inactives")
STRATUM = "pharmacophore_aid504891"
AGENT = "blacksmith"
TOPIC = "pharmacophore_gate"

_MODEL_FPS: list = []
_MODEL_NAMES: list = []


def _init_worker(smiles: list[str], names: list[str]) -> None:
    """Build the model fingerprints once per worker, not once per candidate."""
    cmp.pin_to_one_thread()
    global _MODEL_FPS, _MODEL_NAMES
    _MODEL_FPS = [ph.fingerprint(s) for s in smiles]
    _MODEL_NAMES = names


def _score_chunk(chunk: list[str]) -> list[tuple]:
    from rdkit import DataStructs
    rows = []
    for smi in chunk:
        fp = ph.fingerprint(smi) if isinstance(smi, str) else None
        if fp is None:
            rows.append((np.nan, np.nan, None))
            continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, _MODEL_FPS))
        j = int(sims.argmax())
        rows.append((float(sims[j]), float(sims.mean()), _MODEL_NAMES[j]))
    return rows


def score_parallel(smiles: list[str], model: pd.DataFrame,
                   workers: int) -> pd.DataFrame:
    """`pharmacophore.score`, chunked across processes.

    The serial path is 4.87 ms/molecule -- 29 minutes for this gate's 361,388
    rows -- so this is a convenience, not a necessity. It calls the SAME
    fingerprint and the same Tanimoto as `pharmacophore.score`; a second
    scoring implementation free to drift from the first is how this project
    acquired two definitions of a merge.
    """
    chunks = [smiles[i:i + 2000] for i in range(0, len(smiles), 2000)]
    rows: list[tuple] = []
    with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker,
            initargs=(list(model.canonical_smiles), list(model.name))) as ex:
        for i, res in enumerate(ex.map(_score_chunk, chunks), start=1):
            rows.extend(res)
            if i % 25 == 0:
                log.info("  scored %d/%d", len(rows), len(smiles))
    return pd.DataFrame(rows, columns=[ph.METRIC, "ph2d_mean_similarity",
                                       "ph2d_nearest_active"])


def load_gate_population() -> pd.DataFrame:
    """The 34 assayed actives and 361,354 assayed inactives, labelled."""
    act = pd.read_csv(MEASURED / "aid504891_actives_1.csv")
    ina = pd.read_csv(MEASURED / "aid504891_inactives_1.csv")
    act = act[act.canonical_smiles.notna()].assign(label=1)
    ina = ina[ina.canonical_smiles.notna()].assign(label=0)
    df = pd.concat([act, ina], ignore_index=True)[["canonical_smiles", "label"]]
    log.info("gate population: %d actives, %d measured inactives",
             int(df.label.sum()), int((df.label == 0).sum()))
    return df


def run_one(pop: pd.DataFrame, model: pd.DataFrame, label: str,
            workers: int) -> eg.GateResult:
    log.info("[%s] model set: %d molecules — %s", label, len(model),
             ", ".join(model.name))
    scored = pop.copy()
    scored = pd.concat(
        [scored.reset_index(drop=True),
         score_parallel(list(pop.canonical_smiles), model, workers)], axis=1)
    n_unscored = int(scored[ph.METRIC].isna().sum())
    if n_unscored:
        log.warning("[%s] %d of %d rows could not be fingerprinted and are "
                    "dropped by evaluate()", label, n_unscored, len(scored))
    res = eg.evaluate(scored, metric=ph.METRIC, stratum=STRATUM,
                      higher_is_better=not ph.LOWER_IS_BETTER)
    # Provenance of the score itself, not just of the verdict: which model
    # molecule each active was nearest to. A gate that passes because every
    # active is nearest to ONE model molecule is a different claim from one
    # that passes on the set, and only this makes the difference visible.
    near = scored[scored.label == 1]["ph2d_nearest_active"].value_counts()
    log.info("[%s] actives' nearest model molecule: %s", label, near.to_dict())
    return res, scored


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limit-inactives", type=int, default=None,
                    help="subsample the measured inactives (smoke tests only; "
                         "a subsampled gate is NOT the gate)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    workers = cmp.available_workers(args.workers)
    log.info("using %d CPU workers", workers)

    pop = load_gate_population()
    if args.limit_inactives:
        act = pop[pop.label == 1]
        ina = pop[pop.label == 0].sample(args.limit_inactives, random_state=0)
        pop = pd.concat([act, ina], ignore_index=True)
        log.warning("SUBSAMPLED to %d inactives — this is a smoke test, not "
                    "the gate", args.limit_inactives)

    # The gate actives are AID hits and the model set is the literature
    # binders, so there is no seed to hold out here. Passed explicitly because
    # `build_model_set` refuses to default it.
    declared = ph.build_model_set(hold_out=None)

    results = {}
    res, scored = run_one(pop, declared, "declared rule", workers)
    results["declared"] = res

    # Sensitivity: PiB restored. See the module docstring.
    import pandas as _pd
    from shared import reference_set as rs
    ref = _pd.read_csv(rs.latest_reference("pin1_reference_binders"))
    pib = ref[ref.name == "PiB"]
    if len(pib) and "PiB" not in set(declared.name):
        extra = declared.copy()
        row = {"name": "PiB", "canonical_smiles": pib.iloc[0].canonical_smiles,
               "tier": pib.iloc[0].get("tier"), "n_amide": 4, "n_heavy": 32}
        row["fp"] = ph.fingerprint(row["canonical_smiles"])
        extra = _pd.concat([extra, _pd.DataFrame([row])], ignore_index=True)
        results["with_pib"], _ = run_one(pop, extra, "sensitivity: +PiB",
                                         workers)

    # ---- report -------------------------------------------------------
    print("\n" + "=" * 78)
    print(f"PHARMACOPHORE ENRICHMENT GATE — stratum {STRATUM}")
    print("=" * 78)
    for k, r in results.items():
        print(f"\n[{k}]  VERDICT: {r.verdict}")
        print(f"  ROC-AUC   {r.roc_auc:.4f}  CI[{r.roc_auc_ci[0]:.3f}, "
              f"{r.roc_auc_ci[1]:.3f}]")
        print(f"  EF 1%     {r.ef_1pct:.2f}")
        print(f"  BEDROC    {r.bedroc:.4f}")
        print(f"  actives {r.n_actives}  decoys {r.n_decoys}  "
              f"chemotypes {r.n_chemotypes} ({r.chemotype_method})")
        if r.per_chemotype_auc:
            v = list(r.per_chemotype_auc.values())
            print(f"  leave-one-chemotype-out AUC: min {min(v):.3f} "
                  f"max {max(v):.3f}")
        for reason in r.reasons:
            print(f"    - {reason}")

    # THE TOKEN, SO THE RANKING CAN READ THIS VERDICT. `rank_shortlist.
    # attach_gate` reads the shared token and stamps UNGATED when a stratum is
    # missing -- which reads as "no gate was ever run" rather than "the gate
    # ran and its result went somewhere else".
    #
    # ONLY THE DECLARED RULE IS WRITTEN. The sensitivity arm is a check on the
    # model-set rule, not a second verdict about the scorer; writing both would
    # let a reader pick the friendlier one. `write_token` merges per stratum
    # AND per metric, so this adds a new stratum without touching the
    # non_covalent / covalent verdicts every other approach ranks on.
    eg.write_token([results["declared"]])

    payload = {k: r.to_dict() for k, r in results.items()}
    p = out.write_path(AGENT, TOPIC, "pharmacophore_gate", ".json")
    p.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    man = Manifest(
        stage="pharmacophore_gate", approach="t2",
        params={"stratum": STRATUM, "metric": ph.METRIC,
                "higher_is_better": True,
                "model_rule": {"max_amide_bonds": ph.MAX_AMIDE_BONDS,
                               "max_heavy_atoms": ph.MAX_HEAVY_ATOMS},
                "model_set_declared": list(declared.name),
                "n_actives": int(pop.label.sum()),
                "n_measured_inactives": int((pop.label == 0).sum()),
                "subsampled": bool(args.limit_inactives),
                "verdicts": {k: r.verdict for k, r in results.items()}})
    man.add_input("aid504891_actives", MEASURED / "aid504891_actives_1.csv")
    man.add_input("aid504891_inactives", MEASURED / "aid504891_inactives_1.csv")
    man.add_output("gate", p)
    man.write(p.parent, filename=f"{p.stem}_manifest.json")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
