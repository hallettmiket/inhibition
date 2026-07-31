"""
Purpose: Rebuild the non-covalent decoy set so every non-covalent active has
         enough property-matched decoys to enter the enrichment gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: data/reference/pin1_reference_binders_<latest>.csv (resolved by glob),
       the cached ChEMBL pool
Output: append_only/.../decoys/decoys_non_covalent_2.csv + a per-active report

Run:  python scripts/build_noncovalent_decoys_2.py [--per-active 60] [--dry-run]

WHY THIS EXISTS. Liu-2024-C3 is the sixth independent chemotype on the
non-covalent side, and six is the gate's floor -- with it, T_1 and T_2 can
return a verdict other than UNDERPOWERED for the first time. But it has
MW 547 and the existing decoy set (decoys_non_covalent_1.csv) tops out at
MW 478, so it has ZERO property-matched decoys and the gate's
min_decoys_per_active = 10 rule would exclude it, exactly as it excluded EGCG
and the peptidic macrocycles. The chemotype would be present in the reference
file and absent from the gate.

The cached ChEMBL pool already holds 887 molecules within MW +/-50 and
logP +/-1.5 of Liu-2024-C3. Nothing needs fetching; the previous build simply
never reached that far up the mass range.

WHY append_only AND NOT immutable. immutable/ is read-only by project rule.
The earlier pools were written there before that was enforced; this one is not.
The loader is pointed here explicitly rather than by search order, so which file
a run consumed is recorded rather than inferred.

WHAT A DECOY IS HERE. A ChEMBL molecule matched to one active on molecular
weight and logP, and NOT similar to any known active (Tanimoto < 0.4 on Morgan2
against every entry in the reference set). Property-matched decoys are presumed
inactive, not known inactive -- standard for enrichment benchmarking, and the
similarity cut is what stops a genuine binder being labelled a decoy.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors
from rdkit.Chem import rdFingerprintGenerator as fpg

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("build-noncov-decoys")

# Resolved by glob, never pinned by hand — a literal version here goes stale
# silently, because the old file still exists and still parses. See
# shared.reference_set.latest_reference for the three times that has happened.
from shared.reference_set import latest_reference  # noqa: E402

ACTIVES = latest_reference("pin1_reference_binders")
POOL = Path("/data/lab_vm/immutable/inhibition/decoys/chembl_pool.csv")
OUT_DIR = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/decoys")

# The non-covalent actives the gate can actually use. The peptidic macrocycles
# (Wildemann, Liu-Pei, Jiang-Pei) are deliberately absent: they sit outside the
# ChEMBL pool's chemical space entirely, and matching them on MW/logP alone
# would pair a bicyclic peptide with small-molecule decoys, which is not a
# property match in any sense that makes the comparison mean anything.
GATE_ACTIVES = [
    "ATRA",
    "PiB",
    "Du-Xu-naphthalenecarboxamide",
    "Guo-Pfizer-benzothiophene-phosphonate",
    "Potter-Astex-indole-furancarboxamide",
    "Liu-2024-C3",
]

MW_WINDOW = 50.0
LOGP_WINDOW = 1.5
MAX_SIMILARITY_TO_ANY_ACTIVE = 0.4

_gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)


def _props(smiles: str):
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    return m, Descriptors.MolWt(m), Crippen.MolLogP(m)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-active", type=int, default=60,
                    help="decoys to draw per active (gate minimum is 10)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    ref = pd.read_csv(ACTIVES)
    pool = pd.read_csv(POOL)

    # Every active in the reference set, not just the gate ones: a decoy must be
    # dissimilar to ALL known binders, including the covalent ones and the
    # macrocycles, or we risk labelling a real binder as a negative.
    active_fps = []
    for s in ref["canonical_smiles"].dropna():
        m = Chem.MolFromSmiles(str(s))
        if m is not None:
            active_fps.append(_gen.GetFingerprint(m))
    log.info("similarity screen against %d known binders", len(active_fps))

    cand = []
    for cid, smi in zip(pool["chembl_id"], pool["smiles"]):
        p = _props(smi)
        if p is None:
            continue
        m, mw, logp = p
        cand.append((cid, smi, mw, logp, _gen.GetFingerprint(m)))
    log.info("candidate pool: %d parseable molecules", len(cand))

    rows, report = [], []
    used: set[str] = set()
    for name in GATE_ACTIVES:
        hit = ref[ref["name"] == name]
        if hit.empty:
            log.error("%s not in %s; skipping", name, ACTIVES.name)
            continue
        p = _props(hit.iloc[0]["canonical_smiles"])
        if p is None:
            log.error("%s: unparseable SMILES; skipping", name)
            continue
        _, a_mw, a_logp = p

        matched = []
        for cid, smi, mw, logp, fp in cand:
            if cid in used:
                continue
            if abs(mw - a_mw) > MW_WINDOW or abs(logp - a_logp) > LOGP_WINDOW:
                continue
            sim = max(DataStructs.BulkTanimotoSimilarity(fp, active_fps))
            if sim >= MAX_SIMILARITY_TO_ANY_ACTIVE:
                continue
            matched.append((cid, smi))
            if len(matched) >= args.per_active:
                break

        for cid, smi in matched:
            used.add(cid)
            rows.append({"chembl_id": cid, "smiles": smi,
                         "matched_active": name, "stratum": "non_covalent",
                         "label": 0})
        status = "OK" if len(matched) >= 10 else "BELOW GATE MINIMUM"
        report.append((name, round(a_mw), round(a_logp, 2), len(matched), status))
        log.info("%-40s MW %4.0f logP %5.2f -> %3d decoys  %s",
                 name, a_mw, a_logp, len(matched), status)

    out = pd.DataFrame(rows)
    print()
    print(f"{'active':<40}{'MW':>6}{'logP':>7}{'decoys':>8}  status")
    print("-" * 78)
    for r in report:
        print(f"{r[0]:<40}{r[1]:>6}{r[2]:>7}{r[3]:>8}  {r[4]}")
    print(f"\ntotal decoys: {len(out)}   distinct: {out.chembl_id.nunique()}")

    short = [r[0] for r in report if r[3] < 10]
    if short:
        print(f"\nSTILL BELOW THE GATE MINIMUM: {short}")
        print("These actives will be excluded by the gate exactly as before.")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUT_DIR / "decoys_non_covalent_2.csv"
    out.to_csv(dest, index=False)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
