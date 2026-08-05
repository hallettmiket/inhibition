"""
Purpose: add pH-7.4 charge and a phosphate label to every arm's frame.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: each approach's latest frame
Output: the frame with `charge_ph74`, `charge_class`, `has_phosphate` merged in

#6 items 5 and 7, both decided and neither implemented. They land together
because "label the phosphate rather than protecting it" is only honest once
phosphate-free molecules are actually evaluated rather than being filtered out
by charge stratification.

WHAT STRATIFYING ON CHARGE IS AND IS NOT. Vina carries no electrostatic term --
verified in `noncovalent_dock_run`, which also notes obabel writes all-zero
partial charges either way. So a charge stratum is NOT a claim that the score
models electrostatics. It is the narrower statement that comparing a dianion's
score with a cation's compares two different physical situations, and that a
ranking which mixes them is partly ordering by charge state. Ranking WITHIN a
stratum removes that; it does not make the score valid.

Once these columns exist, stratified ranking needs no new machinery:
`rank_shortlist.rank(..., group_col="charge_class")` already ranks within
groups -- the same mechanism T_4 uses for warhead class.

THE COLUMNS ARE MERGED ON candidate_id, NOT ASSIGNED BY POSITION, and stale
copies are dropped before the merge so a re-run cannot produce `_x`/`_y`. That
is catalogue entry #5's fix applied here rather than rediscovered.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "integration" / "app"))

from shared import io as dio                       # noqa: E402
from shared import ionisation as ion               # noqa: E402

log = logging.getLogger("annotate-ionisation")

COLS = ("charge_ph74", "charge_class", "has_phosphate")

EXPERIMENTS = {
    "t1": ("01_t1_de_novo", "D1"),
    "t2": ("02_t2_atra_crem", "D2"),
    "t3": ("03_t3_reinvent", "D3"),
    "t4": ("04_t4_combinatorial", "D4"),
}


def annotate(approach: str, experiment: str, prefix: str,
             dry_run: bool = False) -> dict:
    frame_path = dio.latest(dio.DATA_ROOT / experiment
                            if hasattr(dio, "DATA_ROOT")
                            else Path("/data/lab_vm/append_only/inhibition")
                            / experiment, prefix, ".parquet")
    if frame_path is None:
        return {"approach": approach, "state": "no-frame"}
    df = dio.read_frame(frame_path)

    smiles_by_id = dict(zip(df["candidate_id"], df["canonical_smiles"]))
    charges = ion.charge_at_ph(smiles_by_id)

    add = pd.DataFrame({
        "candidate_id": list(smiles_by_id),
        "charge_ph74": [charges[c] for c in smiles_by_id],
        "charge_class": [ion.charge_class(charges[c]) for c in smiles_by_id],
        "has_phosphate": [ion.has_phosphate(smiles_by_id[c])
                          for c in smiles_by_id],
    })
    # Nullable Int64: a molecule obabel could not convert has NO charge, and
    # float64 would render it as `-1.0` in the GUI and any CSV export.
    add["charge_ph74"] = add["charge_ph74"].astype("Int64")

    stale = [c for c in (*COLS, *(f"{c}{s}" for c in COLS for s in ("_x", "_y")))
             if c in df.columns]
    if stale:
        log.info("%s: dropping %d stale column(s) before merge", approach,
                 len(stale))
        df = df.drop(columns=stale)
    merged = df.merge(add, on="candidate_id", how="left")
    if len(merged) != len(df):
        raise RuntimeError(f"{approach}: merge changed row count "
                           f"{len(df)} -> {len(merged)}")
    suffixed = [c for c in merged.columns if c.endswith(("_x", "_y"))]
    if suffixed:
        raise RuntimeError(f"{approach}: merge produced {suffixed}")

    counts = merged["charge_class"].value_counts().to_dict()
    n_phos = int(merged["has_phosphate"].fillna(False).sum())
    out = None
    if not dry_run:
        out = dio.write_full_frame(
            merged, approach=approach, experiment=experiment,
            stage=f"{approach}_ionisation",
            params={"ph": ion.DEFAULT_PH,
                    "tool": "obabel -p (the same call docking used)",
                    "phosphate_smarts": ion.PHOSPHATE_SMARTS,
                    "charge_classes": counts,
                    "n_phosphate": n_phos,
                    "n_unknown_charge": int(counts.get("unknown", 0)),
                    "why": "#6 items 5+7 — phosphate is LABELLED not filtered; "
                           "charge_class exists so rank(group_col=...) can "
                           "stratify. Vina has no electrostatic term, so a "
                           "stratum is not a claim about modelled energy."},
            inputs={"frame": frame_path})
    return {"approach": approach, "state": "ok", "frame": frame_path.name,
            "out": out.name if out else "(dry run)", "rows": len(merged),
            "classes": counts, "n_phosphate": n_phos}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--approach", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    todo = ([args.approach] if args.approach else list(EXPERIMENTS))
    print(f"\n{'arm':<5}{'rows':>8}{'anion':>8}{'neutral':>9}{'cation':>8}"
          f"{'unknown':>9}{'phosphate':>11}  frame")
    for a in todo:
        experiment, prefix = EXPERIMENTS[a]
        r = annotate(a, experiment, prefix, dry_run=args.dry_run)
        if r["state"] != "ok":
            print(f"{a:<5}  {r['state']}")
            continue
        c = r["classes"]
        print(f"{a:<5}{r['rows']:>8,}{c.get('anion',0):>8,}"
              f"{c.get('neutral',0):>9,}{c.get('cation',0):>8,}"
              f"{c.get('unknown',0):>9,}{r['n_phosphate']:>11,}  {r['out']}")

    print("\n  charge_class is a LABEL and a stratification key. It is not a "
          "filter,\n  and it is not a claim that the score models "
          "electrostatics — Vina has no\n  electrostatic term. Stratified "
          "ranking: rank(group_col='charge_class').")


if __name__ == "__main__":
    main()
