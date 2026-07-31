"""
Purpose: Deepen a thin decoy chemotype from PubChem when ChEMBL cannot supply it.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-30
Input: a chemotype from decoy_chemotypes_4.csv
Output: append_only/.../decoys/class_pools_4/<class>.csv, merged with ChEMBL

Run:  python scripts/fetch_pubchem_class_pool.py --chemotype sulfamate_acetamide

WHY A SECOND SOURCE. `sulfamate_acetamide` holds SIX molecules in ChEMBL, all
six adduct-valid, and property matching to Reddi-4d and 4g yields 0 and 1
against a gate minimum of 10. The chemotype is not rare in chemistry; it is
rare in a BIOACTIVITY database, because ChEMBL only contains compounds somebody
assayed. PubChem's substructure search returns 310 for the same pattern.

WHY NOT ZINC, WHICH IS WHAT WAS ASKED FOR. ZINC20's documented REST endpoints
redirect to HTML rather than returning data, and ZINC22/CartBlanche is built for
lookup by ZINC ID or exact SMILES, not substructure search over the catalogue.
PubChem is the source that actually supports the query we need. Recorded here so
the next person does not repeat the attempt.

THE BIAS THIS PARTLY FIXES. Every decoy drawn from ChEMBL is a molecule with
measured bioactivity against something, so our "presumed inactives" come from
the most biologically active chemical space there is — enriched for privileged
scaffolds and promiscuous binders. PubChem is broader and includes vendor
catalogue compounds that were never assayed, which is a better prior for
"probably not a binder". It is not a complete fix and is not claimed as one.

MEMBERSHIP IS STILL DECIDED LOCALLY. PubChem's fastsubstructure is the
RETRIEVAL step and is deliberately loose. Every molecule is then put through the
same `verify_class` the ChEMBL path uses — whole-group SMARTS plus a successful
adduct transform — so a change of source cannot change what counts as class
membership.
"""

from __future__ import annotations

import argparse
import io
import logging
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("pubchem-pool")

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BATCH = 200


def _get(url: str, retries: int = 4, timeout: int = 300) -> bytes:
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"PubChem request failed after {retries}: {last}")


def cids_for(query_smiles: str, max_records: int) -> list[int]:
    q = urllib.parse.quote(query_smiles)
    url = f"{PUG}/compound/fastsubstructure/smiles/{q}/cids/TXT?MaxRecords={max_records}"
    txt = _get(url).decode(errors="ignore")
    return [int(x) for x in txt.split() if x.strip().isdigit()]


def smiles_for(cids: list[int]) -> pd.DataFrame:
    rows = []
    for i in range(0, len(cids), BATCH):
        chunk = ",".join(str(c) for c in cids[i:i + BATCH])
        url = f"{PUG}/compound/cid/{chunk}/property/CanonicalSMILES/CSV"
        try:
            df = pd.read_csv(io.BytesIO(_get(url)))
            rows.append(df)
        except Exception as exc:  # noqa: BLE001
            log.warning("batch %d failed: %s", i, str(exc)[:120])
        time.sleep(0.25)          # PubChem asks for <= 5 requests/second
        log.info("fetched %d/%d SMILES", min(i + BATCH, len(cids)), len(cids))
    if not rows:
        return pd.DataFrame(columns=["CID", "CanonicalSMILES"])
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chemotype", required=True)
    ap.add_argument("--max-records", type=int, default=5000)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from shared import decoys_classmatched as dcm
    from shared import warhead_library as wl

    lib = wl.load()
    chemo = dcm.load_chemotypes()
    row = chemo[chemo["chemotype"] == args.chemotype]
    if row.empty:
        raise SystemExit(f"no chemotype {args.chemotype!r} in {dcm.CHEMOTYPES.name}")
    row = row.iloc[0]
    rep = str(row["representative_class"])
    query = str(row["chembl_query"])
    log.info("chemotype %s: query %r -> class test %s", args.chemotype, query, rep)

    cids = cids_for(query, args.max_records)
    log.info("PubChem returned %d CIDs", len(cids))
    if not cids:
        raise SystemExit("PubChem returned nothing for this query")

    df = smiles_for(cids)
    # PubChem renames this column between API revisions -- it returned
    # `CanonicalSMILES` and now returns `ConnectivitySMILES`. Hard-coding the
    # name produced a KeyError; picking whichever SMILES-like column came back
    # keeps working, and raising when there is none beats guessing.
    smi_cols = [c for c in df.columns if "SMILES" in c.upper()]
    if not smi_cols:
        raise SystemExit(f"no SMILES column in PubChem response; got {list(df.columns)}")
    if len(smi_cols) > 1:
        log.warning("several SMILES columns %s; using %s", smi_cols, smi_cols[0])
    df = df.rename(columns={"CID": "cid", smi_cols[0]: "smiles"}).dropna()
    log.info("resolved %d SMILES from column %r", len(df), smi_cols[0])

    patt = dcm.warhead_group_pattern(row)
    keep = []
    for cid, s in zip(df["cid"], df["smiles"]):
        if dcm.verify_class(str(s), rep, patt, lib):
            keep.append({"chembl_id": f"PUBCHEM{int(cid)}", "smiles": str(s)})
    log.info("verify_class passes: %d of %d", len(keep), len(df))

    cache = dcm.CLASS_POOL_DIR / f"{args.chemotype}.csv"
    existing = pd.DataFrame(columns=["chembl_id", "smiles"])
    if cache.is_file():
        try:
            existing = pd.read_csv(cache)
        except pd.errors.EmptyDataError:
            pass
    merged = (pd.concat([existing, pd.DataFrame(keep)], ignore_index=True)
              .drop_duplicates("smiles"))
    print(f"\n{args.chemotype}: ChEMBL {len(existing)} + PubChem {len(keep)} "
          f"-> {len(merged)} unique (gate minimum 10 AFTER property matching)")

    if args.dry_run:
        print("--dry-run: nothing written")
        return
    dcm.CLASS_POOL_DIR.mkdir(parents=True, exist_ok=True)
    merged.to_csv(cache, index=False)
    print(f"wrote {cache}")


if __name__ == "__main__":
    main()
