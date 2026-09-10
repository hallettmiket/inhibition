"""
Purpose: can the scorer recover a KNOWN Pin1 binder it was not shown, out of 361,354 measured inactives?
Author: @tt8804 (with Claude Code)
Date: 2026-09-09
Input: the frozen reference binder set + AID 504891's assayed inactives
Output: each held-out binder's percentile among the inactives, and the median

WHY THIS IS THE STRONGER TEST. D0117 graded the scorer on recovering AID
504891's qHTS hits, which are hits of unverified mechanism from a different
population. This asks a strictly easier and more direct question:

    hide ONE published Pin1 binder from the model set, score it against
    361,354 molecules that were ASSAYED and found inactive, and see where
    it lands.

If a ligand-based scorer cannot rank a real Pin1 binder above measured
inactives when FOURTEEN of that binder's peers are in the model set, it is not
going to rank a novel CReM product. And unlike the AID gate, a failure here
cannot be explained away by the actives being a different population: these are
the reference set's own molecules.

This is the habit `how_this_project_breaks.md` scores as the best-performing
route of all -- run the pipeline on something whose answer is already written
down. The answer written down here is "these 15 molecules bind Pin1".

WHY IT IS CHEAP DESPITE 15 FOLDS. The inactives' fingerprints do not depend on
the model set, so each worker computes a molecule's fingerprint ONCE and
evaluates all 15 folds against it. The cost is one pass over the inactives, not
fifteen.

READ THE PERCENTILE, NOT THE VERDICT. There is no gate token written here and
no STRONG/WEAK grade: this is a recovery measurement, not a second enrichment
gate, and dressing it as one would put a second verdict on the same scorer for
a reader to choose between.
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
from shared import outputs as out                    # noqa: E402
from shared import pharmacophore as ph               # noqa: E402
from shared.manifest import Manifest                 # noqa: E402

log = logging.getLogger("pharmacophore-loo")

MEASURED = Path("/data/lab_vm/append_only/inhibition/00_outputs/blacksmith"
                "/measured_inactives")
AGENT = "blacksmith"
TOPIC = "pharmacophore_gate"

_FOLDS: list[list] = []          # per fold: the 14 model fingerprints


def _init(fold_smiles: list[list[str]]) -> None:
    cmp.pin_to_one_thread()
    global _FOLDS
    _FOLDS = [[ph.fingerprint(s) for s in fold] for fold in fold_smiles]


def _score_chunk(chunk: list[str]) -> np.ndarray:
    """max-Tanimoto for every fold, for every molecule in the chunk."""
    from rdkit import DataStructs
    rows = np.full((len(chunk), len(_FOLDS)), np.nan)
    for i, smi in enumerate(chunk):
        fp = ph.fingerprint(smi) if isinstance(smi, str) else None
        if fp is None:
            continue
        for j, model in enumerate(_FOLDS):
            rows[i, j] = max(DataStructs.BulkTanimotoSimilarity(fp, model))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limit-inactives", type=int, default=None,
                    help="smoke tests only; a subsampled background is not the test")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")
    workers = cmp.available_workers(args.workers)

    model = ph.build_model_set(hold_out=None)
    names = list(model.name)
    smiles = list(model.canonical_smiles)
    n = len(names)
    log.info("leave-one-out over %d model molecules, %d workers", n, workers)

    # Fold j holds out molecule j. The held-out molecule is scored against the
    # other n-1 -- never against itself, which would return 1.0 by definition.
    folds = [[s for k, s in enumerate(smiles) if k != j] for j in range(n)]

    ina = pd.read_csv(MEASURED / "aid504891_inactives_1.csv")
    ina = ina[ina.canonical_smiles.notna()]
    if args.limit_inactives:
        ina = ina.sample(args.limit_inactives, random_state=0)
        log.warning("SUBSAMPLED to %d inactives — smoke test, not the measurement",
                    len(ina))
    background = list(ina.canonical_smiles)
    log.info("background: %d measured inactives", len(background))

    # The held-out binders go through the SAME code path as the background, in
    # the same call, so a difference between them cannot come from a difference
    # in how they were scored.
    everything = smiles + background
    chunks = [everything[i:i + 2000] for i in range(0, len(everything), 2000)]
    parts = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init,
                             initargs=(folds,)) as ex:
        for i, res in enumerate(ex.map(_score_chunk, chunks), start=1):
            parts.append(res)
            if i % 25 == 0:
                log.info("  scored %d/%d", sum(len(p) for p in parts),
                         len(everything))
    mat = np.vstack(parts)                       # (n + n_background, n_folds)
    held = mat[:n, :]
    bg = mat[n:, :]

    rows = []
    for j, name in enumerate(names):
        score = held[j, j]                       # molecule j under fold j
        col = bg[:, j]
        col = col[~np.isnan(col)]
        better = int((col >= score).sum())
        pct = 100.0 * better / len(col)          # smaller = ranked higher
        rows.append({"held_out": name, "score": float(score),
                     "n_inactives_scoring_at_least_as_high": better,
                     "percentile_from_top": round(pct, 4),
                     "in_top_1pct": bool(pct <= 1.0)})
    res = pd.DataFrame(rows).sort_values("percentile_from_top")

    print("\n" + "=" * 84)
    print("LEAVE-ONE-OUT RECOVERY — a known Pin1 binder vs 361,354 MEASURED inactives")
    print("=" * 84)
    print(res.to_string(index=False))
    med = res.percentile_from_top.median()
    n_top1 = int(res.in_top_1pct.sum())
    print(f"\nmedian percentile from the top : {med:.3f}%")
    print(f"held-out binders in the top 1% : {n_top1}/{n}")
    print(f"expected in top 1% by chance   : {0.01*n:.2f}")
    print("\nA scorer that ranks these molecules is one whose held-out binders sit")
    print("near the top. Chance puts the median at 50%.")

    p = out.write_path(AGENT, TOPIC, "pharmacophore_loo_recovery", ".json")
    p.write_text(json.dumps({"median_percentile": med, "n_top_1pct": n_top1,
                             "n_folds": n, "n_background": len(bg),
                             "per_fold": rows}, indent=2), encoding="utf-8")
    man = Manifest(
        stage="pharmacophore_loo_recovery", approach="t2",
        params={"metric": ph.METRIC, "n_folds": n,
                "n_background_inactives": len(bg),
                "subsampled": bool(args.limit_inactives),
                "median_percentile_from_top": med, "n_top_1pct": n_top1,
                "model_rule": {"max_amide_bonds": ph.MAX_AMIDE_BONDS,
                               "max_heavy_atoms": ph.MAX_HEAVY_ATOMS},
                "note": "a recovery measurement, NOT a second enrichment gate; "
                        "no verdict is written and no gate token is touched"})
    man.add_input("inactives", MEASURED / "aid504891_inactives_1.csv")
    man.add_output("recovery", p)
    man.write(p.parent, filename=f"{p.stem}_manifest.json")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
