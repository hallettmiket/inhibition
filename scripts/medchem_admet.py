"""
Purpose: run ADMET-AI over a small set of SMILES and emit one tidy CSV per input.
Author: @tt8804 (with Claude Code)
Date: 2026-08-06
Input: a .smi file (one SMILES per line) and an output directory
Output: <outdir>/admet_ai_<N>.csv -- one row per molecule, one column per endpoint,
        plus the DrugBank-approved percentile for every endpoint

WHY A DRIVER AND NOT `admet_predict`. The env's console script dies on import:
`chemprop 2.3.0` does a hard `import cuik_molmaker` at
`chemprop/featurizers/molgraph/molecule.py:4`, and the installed
`cuik_molmaker_pin` wheel is linked against the auditwheel-mangled RDKit
libraries of the *pip* rdkit wheel (`libRDKitAbbreviations-a55b7b38.so.1`),
which this env does not have -- its RDKit is the conda build. The import
therefore raises `ImportError` before any model is loaded.

THE STUB IS SAFE, AND THAT IS A CLAIM WITH A REASON. `cuik_molmaker` is
referenced in exactly one class, `CuikmolmakerMolGraphFeaturizer` (lines
198-225 of that file), which is an *alternative* batched featurizer. ADMET-AI
loads its checkpoints through `chemprop.models.load_model`, and those
checkpoints carry `SimpleMoleculeMolGraphFeaturizer`. So the module object is
imported and never called. `--assert-featurizer` (default on) checks that at
run time rather than trusting this paragraph: it walks every loaded model and
fails if any featurizer is the cuik one. A stub that silently changed the
featurisation would change every number in this file, so the check is a
required-fail probe and not a nicety.

NOT PATCHING THE ENV. Per D0010, pip installs stay behind the target env's own
bin, and this env belongs to another user. The stub lives in this process only.
"""

from __future__ import annotations

import argparse
import logging
import sys
import types
from pathlib import Path

log = logging.getLogger("medchem-admet")


def _stub_cuik_molmaker() -> None:
    """Pre-empt `import cuik_molmaker` with an empty module. See the header."""
    if "cuik_molmaker" not in sys.modules:
        sys.modules["cuik_molmaker"] = types.ModuleType("cuik_molmaker")


def _assert_no_cuik_featurizer(model) -> None:
    """Fail loudly if any loaded model actually wants the stubbed featurizer.

    The stub is only sound because nothing calls it. This is the guard that can
    fail, rather than a comment asserting it cannot.
    """
    bad = []
    for task_models in model.model_lists:
        for m in task_models:
            f = getattr(getattr(m, "message_passing", None), "featurizer", None)
            f = f or getattr(m, "featurizer", None)
            if f is not None and "Cuik" in type(f).__name__:
                bad.append(type(f).__name__)
    if bad:
        raise SystemExit(
            "ABORT: a loaded chemprop model uses "
            f"{sorted(set(bad))}, which this process stubbed out. "
            "Every prediction would be featurised by an empty module."
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smiles", required=True, type=Path,
                    help="file with one SMILES per line")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--version", type=int, default=1,
                    help="integer suffix; never overwrite an existing version")
    ap.add_argument("--no-assert-featurizer", dest="assert_featurizer",
                    action="store_false")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    _stub_cuik_molmaker()

    import pandas as pd
    from admet_ai import ADMETModel

    smis = [s.strip() for s in args.smiles.read_text().splitlines() if s.strip()]
    log.info("%d SMILES from %s", len(smis), args.smiles)

    model = ADMETModel()
    if args.assert_featurizer:
        _assert_no_cuik_featurizer(model)
        log.info("featurizer check passed: no model uses the stubbed module")

    preds: pd.DataFrame = model.predict(smiles=smis)
    preds.index.name = "smiles"

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"admet_ai_{args.version}.csv"
    if out.exists():  # append-only: version, do not overwrite
        raise SystemExit(f"{out} exists -- bump --version rather than overwriting")
    preds.to_csv(out)
    log.info("wrote %s  (%d rows x %d cols)", out, len(preds), preds.shape[1])


if __name__ == "__main__":
    main()
