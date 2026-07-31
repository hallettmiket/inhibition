"""
Purpose: Build the class-matched covalent decoy set (D0031).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: data/reference/pin1_reference_binders_2.csv (frozen actives)
Output: append_only/.../decoys/decoys_covalent_10.csv + provenance

Run:  python scripts/build_covalent_decoys.py [--force-fetch] [--n-per-active 50]

Every covalent active is matched to decoys drawn from ITS OWN warhead class.
Where a class has no usable decoys the shortfall is REPORTED and the active is
marked untestable, never topped up from another class — that top-up is what
D0028 found the previous set doing, and it is what made the gate partly a
chemotype comparison instead of a binding one.
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

from shared import decoys_classmatched as dcm      # noqa: E402
from shared import smiles as smi                   # noqa: E402
from shared import warhead_library as wl           # noqa: E402
from shared.manifest import Manifest               # noqa: E402

log = logging.getLogger("build-decoys")

ACTIVES = REPO / "data" / "reference" / "pin1_reference_binders_3.csv"
# append_only, not immutable: immutable/ is read-only by project rule, and the
# earlier versions were written there before that was enforced.
OUT = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate/decoys")
OUT_NAME = "decoys_covalent_10.csv"

# The covalent actives, mapped to the library class each one actually belongs
# to. Kept here rather than parsed from the free-text `warhead_class` column,
# which holds prose like "1;4-naphthoquinone (Michael acceptor)".
# Mapped to CHEMOTYPE (data/reference/decoy_chemotypes_1.csv), not to a warhead
# library class. The library's classes are T_4's enumeration attachment points:
# `naphthoquinone_c2` and `_benzo` are two positions on one chemistry, and
# Juglone (5-hydroxy) sits at neither, so classifying it by an enumeration unit
# rejected it from its own chemotype.
ACTIVE_CHEMOTYPE = {
    "Sulfopin": "chloroacetamide",
    "BJP-06-005-3": "chloroacetamide",
    "Reddi-2023-4d": "sulfamate_acetamide",
    "Reddi-2023-4g": "sulfamate_acetamide",
    "Juglone": "naphthoquinone",
    "Tian-chloropyrimidine-covalent-6a": "snar_chloroazine",
    "Ieda-2019-(S)-2": "cinnamamide",
}

# KPT-6566 IS DELIBERATELY ABSENT. It is a self-immolative aryl-sulfonyl-acetate:
# the species that alkylates Cys113 is a naphthoquinone RELEASED from it, not the
# deposited molecule. Treating the parent as a naphthoquinone would dock a
# structure that never reaches the cysteine — the same "which species is being
# scored" error as D0022 and D0030. Scoring it needs the released fragment as its
# own entry, which the actives table does not currently carry.
EXCLUDED_ACTIVES = {
    "KPT-6566": "self-immolative prodrug; the covalent species is a released "
                "naphthoquinone, not this molecule",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-per-active", type=int, default=50)
    ap.add_argument("--force-fetch", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    lib = wl.load()
    chemo = dcm.load_chemotypes()
    actives = pd.read_csv(ACTIVES)
    for name, why in EXCLUDED_ACTIVES.items():
        log.warning("excluding active %s: %s", name, why)
    actives = actives[actives["name"].isin(ACTIVE_CHEMOTYPE)].copy()
    actives["library_class"] = actives["name"].map(ACTIVE_CHEMOTYPE)
    actives["canonical_smiles"] = actives["canonical_smiles"].map(smi.canonical)
    bad = actives[actives["canonical_smiles"].isna()]
    if len(bad):
        log.warning("dropping %d active(s) with unusable SMILES: %s",
                    len(bad), list(bad["name"]))
        actives = actives.dropna(subset=["canonical_smiles"])
    log.info("%d covalent actives across %d chemotypes", len(actives),
             actives["library_class"].nunique())

    # SELF-CHECK: every active must satisfy its OWN class test.
    #
    # The class test decides which decoys count as same-class, so if an active
    # fails it the test is wrong, not the active. Built without this, a SMILES/
    # SMARTS mix-up made the whole-group pattern match nothing at all and every
    # class was reported "chemically unavailable" while the pools held thousands
    # of molecules. A control whose own positives fail it is not a control.
    failures = []
    for _, a in actives.iterrows():
        cls = a["library_class"]
        row = chemo[chemo["chemotype"] == cls].iloc[0]
        patt = dcm.warhead_group_pattern(row)
        if not dcm.verify_class(a["canonical_smiles"],
                                str(row["representative_class"]), patt, lib):
            failures.append(f"{a['name']} ({cls})")
    if failures:
        raise SystemExit(
            "the class test rejects its own actives, so it cannot decide class "
            "membership for decoys either: " + "; ".join(failures))
    log.info("self-check: all %d actives satisfy their own class test", len(actives))

    active_fps = [dcm._fp(s) for s in actives["canonical_smiles"]]
    active_fps = [f for f in active_fps if f is not None]

    pools: dict[str, dcm.ClassPool] = {}
    for cls in sorted(actives["library_class"].unique()):
        log.info("--- building chemotype pool: %s ---", cls)
        pools[cls] = dcm.build_class_pool(cls, lib, chemo, force=args.force_fetch)
        p = pools[cls]
        log.info("[%s] retrieved %d -> whole-group %d -> adduct-valid %d",
                 cls, p.n_retrieved, p.n_group_match, p.n_adduct_ok)
        if p.unavailable_reason:
            log.warning("[%s] UNUSABLE: %s", cls, p.unavailable_reason)

    frames, report, used = [], [], set()
    for _, a in actives.iterrows():
        cls = a["library_class"]
        got = dcm.match_within_class(a, pools[cls], active_fps,
                                     n_per_active=args.n_per_active, used=used)
        if len(got):
            used.update(got["chembl_id"])
            got = got.assign(active=a["name"], warhead_class=cls,
                             class_matched=True)
            frames.append(got)
        report.append({
            "active": a["name"], "warhead_class": cls,
            "n_decoys": int(len(got)),
            "target": args.n_per_active,
            "shortfall": int(args.n_per_active - len(got)),
            "testable": bool(len(got) >= 10),
            "class_pool_adduct_valid": pools[cls].n_adduct_ok,
            "unavailable_reason": pools[cls].unavailable_reason,
        })
        log.info("[%s] %s -> %d decoys", cls, a["name"], len(got))

    rep = pd.DataFrame(report)
    print("\n=== per-active decoy yield ===")
    print(rep.to_string(index=False))

    if not frames:
        raise SystemExit("no class-matched decoys could be built at all")

    ds = pd.concat(frames, ignore_index=True).drop_duplicates("chembl_id")
    OUT.mkdir(parents=True, exist_ok=True)
    out_csv = OUT / OUT_NAME
    ds.to_csv(out_csv, index=False)

    man = Manifest(stage="build_covalent_decoys",
                   params={"n_per_active": args.n_per_active,
                           "match_tolerance": dcm.MATCH_TOLERANCE,
                           "max_similarity_to_active": dcm.MAX_SIMILARITY_TO_ACTIVE,
                           "class_matched": True,
                           "cross_class_topup": False},
                   inputs={"actives": ACTIVES})
    (OUT / OUT_NAME.replace(".csv", "_report.json")).write_text(
        json.dumps({"per_active": report,
                    "class_pools": {c: {"query": p.query,
                                        "retrieved": p.n_retrieved,
                                        "whole_group": p.n_group_match,
                                        "adduct_valid": p.n_adduct_ok,
                                        "unavailable_reason": p.unavailable_reason}
                                    for c, p in pools.items()}},
                   indent=2), encoding="utf-8")
    try:
        man.write(out_csv)
    except Exception as exc:  # noqa: BLE001 - the CSV is the deliverable
        log.warning("manifest not written: %s", exc)

    print(f"\nwrote {len(ds)} class-matched decoys -> {out_csv}")
    print(f"  classes covered: {sorted(ds['warhead_class'].unique())}")
    untestable = rep[~rep["testable"]]
    if len(untestable):
        print("\n  UNTESTABLE ACTIVES (fewer than 10 same-class decoys):")
        for _, r in untestable.iterrows():
            print(f"    {r['active']} ({r['warhead_class']}): {r['n_decoys']} "
                  f"— {r['unavailable_reason'] or 'property matching found too few'}")
        print("\n  These are NOT topped up from another class. Doing so is what")
        print("  made the previous gate partly a chemotype comparison (D0028).")


if __name__ == "__main__":
    main()
