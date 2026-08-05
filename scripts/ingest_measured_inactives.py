"""
Purpose: turn PubChem AID 504891 into an ASSAYED active/inactive set for the enrichment gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: the AID 504891 datatable download (387,572 rows, SMILES included)
Output: 00_outputs/blacksmith/measured_inactives/{actives,inactives}_<N>.csv

#4 Phase 1. The gate's weakest component is its NEGATIVES: they are
property-matched decoys ASSUMED to be inactive, and D0041's null could in
principle be an artefact of how they were built. This dataset replaces the
assumption with a measurement.

AID 504891 is *qHTS Assay to Find Inhibitors of Pin1* (NCGC, PI Kun Ping Lu) --
the same Pin1, and the same lab #12's questions are addressed to.

## Inconclusive is neither, and is dropped

    Active         34
    Inactive  364,905
    Inconclusive  22,628

The 22,628 inconclusive compounds are NOT inactives. Folding them in would add
6% of the pool as negatives that the assay explicitly declined to call, and
every one that is really an active would push enrichment DOWN -- making the
gate look worse for a reason that is a data-handling choice rather than a
property of docking. They are written to their own file so the decision is
visible and reversible, never merged.

## Leakage is checked, not assumed

Our reference actives came from ChEMBL2288. If any of them also appear in this
assay's INACTIVE list, we would be scoring a molecule as a decoy that we
elsewhere call a known binder. The overlap is computed on InChIKey and
reported; overlapping compounds are removed from the inactive pool, because a
compound cannot be both.

That check is the reason this is a script and not a one-line download: the same
reasoning that ruled out Boltz-2 (`state_of_the_project` section 4) applies
here in reverse.

## What this does NOT do

It does not run the gate. Building the pool and measuring enrichment against it
are separate steps on purpose -- the pool should be inspectable before anything
is concluded from it.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                 # noqa: E402
from shared import reference_set as rs             # noqa: E402
from shared import smiles as smi                   # noqa: E402

log = logging.getLogger("ingest-measured")

OUT = sout.Topic("blacksmith", "measured_inactives")
SOURCE = Path("/data/lab_vm/modifiable/inhibition/measured_inactives"
              "/aid504891_full.csv")
AID = 504891

SMILES_COL = "PUBCHEM_EXT_DATASOURCE_SMILES"
OUTCOME_COL = "PUBCHEM_ACTIVITY_OUTCOME"
KEEP = ["PUBCHEM_CID", "PUBCHEM_SID", SMILES_COL, OUTCOME_COL,
        "PUBCHEM_ACTIVITY_SCORE", "Potency", "Efficacy", "Phenotype"]


def load() -> pd.DataFrame:
    if not SOURCE.is_file():
        raise SystemExit(
            f"{SOURCE} not found. Download it first:\n"
            "  curl -L -o aid504891_full.csv 'https://pubchem.ncbi.nlm.nih.gov"
            "/assay/pcget.cgi?query=download&record_type=datatable&actvty=all"
            f"&response_type=save&aid={AID}'")
    df = pd.read_csv(SOURCE, low_memory=False,
                     usecols=lambda c: c in KEEP)
    log.info("loaded %d rows from %s", len(df), SOURCE.name)
    return df


def canonicalise(df: pd.DataFrame) -> pd.DataFrame:
    out = df[df[SMILES_COL].notna()].copy()
    log.info("%d rows carry a SMILES", len(out))
    out["canonical_smiles"] = [smi.canonical(s) for s in out[SMILES_COL]]
    bad = out["canonical_smiles"].isna().sum()
    if bad:
        log.warning("%d SMILES did not canonicalise and are dropped", bad)
    out = out[out["canonical_smiles"].notna()].copy()
    out["inchikey"] = [smi.inchikey(s) for s in out["canonical_smiles"]]
    out = out[out["inchikey"].notna()]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    raw = load()
    counts = raw[OUTCOME_COL].value_counts().to_dict()
    log.info("outcomes: %s", counts)

    df = canonicalise(raw)
    actives = df[df[OUTCOME_COL] == "Active"].drop_duplicates("inchikey")
    inactives = df[df[OUTCOME_COL] == "Inactive"].drop_duplicates("inchikey")
    incon = df[df[OUTCOME_COL] == "Inconclusive"].drop_duplicates("inchikey")

    # LEAKAGE: a molecule cannot be both a known binder and a decoy.
    ref = pd.read_csv(rs.DEFAULT_MASTER)
    ref_keys = {k for k in (smi.inchikey(s) for s in ref["canonical_smiles"])
                if k}
    overlap = inactives[inactives["inchikey"].isin(ref_keys)]
    if len(overlap):
        log.warning("%d reference binder(s) appear in this assay's INACTIVE "
                    "list and are removed from the decoy pool: %s",
                    len(overlap), list(overlap["PUBCHEM_CID"])[:10])
    inactives = inactives[~inactives["inchikey"].isin(ref_keys)]

    active_overlap = actives[actives["inchikey"].isin(ref_keys)]

    print(f"\nAID {AID} — qHTS Assay to Find Inhibitors of Pin1 (NCGC, Lu)")
    print(f"  raw rows                {len(raw):,}")
    print(f"  assayed ACTIVE          {len(actives):,}")
    print(f"  assayed INACTIVE        {len(inactives):,}   "
          "<- measured negatives, the point of this")
    print(f"  inconclusive (dropped)  {len(incon):,}")
    print(f"  reference binders also called inactive here   {len(overlap)}")
    print(f"  reference binders also called ACTIVE here     {len(active_overlap)}"
          "   <- independent corroboration")

    if args.dry_run:
        print("\n  (dry run — nothing written)")
        return

    cols = ["PUBCHEM_CID", "PUBCHEM_SID", "canonical_smiles", "inchikey",
            OUTCOME_COL, "PUBCHEM_ACTIVITY_SCORE", "Potency", "Efficacy"]
    cols = [c for c in cols if c in df.columns]
    for name, frame in (("actives", actives), ("inactives", inactives),
                        ("inconclusive", incon)):
        dest = OUT.write(f"aid{AID}_{name}", ".csv")
        frame[cols].to_csv(dest, index=False)
        print(f"  wrote {len(frame):>7,} -> {dest.name}")

    print("\n  These negatives were ASSAYED, not assumed. Phase 2.1 measures "
          "enrichment\n  against them; on D0046's evidence it should still "
          "fail, and that is the test.")


if __name__ == "__main__":
    main()
