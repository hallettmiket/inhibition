"""
Purpose: dock the non-covalent shortlists into all four ensemble receptors and combine on the MEDIAN.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: T_1 + every T_2 seed variant's shortlist; the four prepared receptors
Output: 00_outputs/blacksmith/ensemble_shortlist/ensemble_shortlist_<N>.csv

#6 item 6 / D0052, at the scope decided on 2026-08-04: **shortlists first**.
Mike left the choice open -- full pools (~265 GPU-hours) or shortlists only
(minutes) -- and observed that the shortlist run answers the question actually
worth buying: *does the ensemble change who is on top?* The full re-dock is
only worth it if the answer is yes.

## The combination rule was pre-registered, and is not chosen here

`vina_affinity_ensemble_median` is the rank metric. That was fixed in D0052
BEFORE any ensemble result existed, which is the D0045 discipline: choosing the
rule after seeing which one improves the ranking is choosing the answer.

Best-across is computed and carried, never sorted on. It is a maximum over four
correlated draws, so its upward bias grows with the width of a ligand's score
distribution -- and that width scales with conformational flexibility.
liu_2024_c3 averages 10.65 rotatable bonds against du_xu's 4.81, so ranking on
best-across would hand the flexible pool an advantage unrelated to binding,
reintroducing the artefact class D0049 removed and doing it invisibly, because
"we docked into an ensemble" reads as a refinement.

## Every receptor is docked FRESH, including 6VAJ

The frames already hold a 6VAJ `vina_affinity`, and reusing it would save a few
minutes. It is not done: those scores came from a different invocation with a
different random seed, and a median over three fresh draws plus one old one is
not a median over four comparable numbers. The whole point of an ensemble is
that the receptor is the only thing that varies.

## What this does NOT do

It does not rank anything. A median over four receptors is a NEW metric and the
enrichment gate has never been run on it; D0016/D0041's verdict is a 6VAJ
number measured with `box_expanded.json`. Since D0051 an unknown metric returns
UNGATED and fails closed, so this is safe rather than silent -- but the gate
must run on the ensemble metric before it orders anything, and that is a
separate step.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "integration" / "app"))

from shared import noncovalent_dock_run as ncd     # noqa: E402
from shared import outputs as sout                 # noqa: E402

log = logging.getLogger("ensemble-shortlist")

OUT = sout.Topic("blacksmith", "ensemble_shortlist")
WORK = Path("/data/lab_vm/append_only/inhibition/06_ensemble_shortlist/docking")


def collect_shortlist() -> pd.DataFrame:
    """Every non-covalent shortlisted candidate, across T_1 and all T_2 seeds.

    T_3/T_4 are excluded: they are covalent and dock through gnina with a
    bonded restraint, so a Vina score on them is not the same quantity.
    """
    import data as appdata                          # noqa: PLC0415
    rows = []
    df1, _ = appdata.load_frame("t1")
    col = appdata.shortlist_column(df1)
    s = df1[df1[col] == True]                       # noqa: E712
    rows.append(s.assign(pool="t1"))
    for key, v in appdata.variants("t2").items():
        d, _ = appdata.load_variant_frame("t2", key)
        if d is None:
            continue
        c = appdata.shortlist_column(d)
        if c not in d.columns:
            continue
        rows.append(d[d[c] == True].assign(pool=f"t2:{key}"))  # noqa: E712
    out = pd.concat(rows, ignore_index=True)
    # One molecule can be shortlisted in more than one T_2 seed pool. Dock it
    # ONCE -- the pose is a property of the molecule and the receptor, not of
    # which pool surfaced it -- and carry the pools it came from.
    pools = (out.groupby("candidate_id")["pool"]
             .apply(lambda s: ",".join(sorted(set(s)))).rename("pools"))
    out = out.drop_duplicates("candidate_id").merge(pools, on="candidate_id")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--gpus", default="1,2,3,5",
                    help="one GPU per receptor, in ENSEMBLE order")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    receptors = list(ncd.ENSEMBLE)
    missing = [r.tag for r in receptors
               if not (r.pdbqt.is_file() and r.box.is_file())]
    if missing:
        raise SystemExit(
            f"receptors not prepared: {missing}. Run "
            "scripts/prepare_ensemble_receptors.py first.")
    log.info("ensemble: %s", ", ".join(r.tag for r in receptors))

    shortlist = collect_shortlist()
    if args.limit:
        shortlist = shortlist.head(args.limit)
    log.info("%d distinct shortlisted molecules across T_1 + every T_2 seed",
             len(shortlist))

    WORK.mkdir(parents=True, exist_ok=True)
    ligand_dir = WORK / f"ligands_{ncd.LIGAND_PREP_TAG}"
    prep = ncd.prepare_ligands(shortlist, ligand_dir)
    n_ready = sum(1 for r in prep if r["ok"])
    log.info("%d/%d ligands prepared", n_ready, len(prep))
    if not n_ready:
        raise SystemExit("no ligand survived preparation")

    gpus = [int(g) for g in args.gpus.split(",")]
    per_receptor: dict[str, pd.DataFrame] = {}
    for i, rec in enumerate(receptors):
        gpu = gpus[i % len(gpus)]
        out_dir = ncd.pose_dir(WORK, rec)
        log.info("docking %d ligands into %s on GPU %d -> %s",
                 n_ready, rec.tag, gpu, out_dir.name)
        t0 = time.time()
        ncd.run_vina_gpu(ligand_dir, out_dir, gpu, rec)
        modes = ncd.collect_modes(out_dir, rec)
        log.info("  %s: %d ligands scored in %.1f min",
                 rec.tag, len(modes), (time.time() - t0) / 60)
        per_receptor[rec.tag] = modes

    combined = ncd.combine_ensemble(per_receptor)
    combined = combined.merge(
        shortlist[["candidate_id", "pools", "canonical_smiles"]],
        on="candidate_id", how="left")
    # The existing 6VAJ score from the production frames, carried for
    # comparison ONLY -- it is a different invocation with a different seed and
    # is not part of the median.
    combined = combined.merge(
        shortlist[["candidate_id", "vina_affinity"]]
        .rename(columns={"vina_affinity": "vina_affinity_production_6VAJ"}),
        on="candidate_id", how="left")

    dest = OUT.write("ensemble_shortlist", ".csv")
    combined.to_csv(dest, index=False)

    print(f"\nensemble shortlist -> {dest}")
    print(f"  molecules            {len(combined)}")
    print(f"  receptors            {', '.join(per_receptor)}")
    full = int((combined[ncd.ENSEMBLE_N] == len(per_receptor)).sum())
    print(f"  scored on all {len(per_receptor)}      {full}")
    print(f"\n  {ncd.ENSEMBLE_MEDIAN} is the rank metric (D0052, "
          "pre-registered).")
    print("  It is UNGATED: the enrichment gate has never been run on a "
          "median-over-four\n  metric, and D0016/D0041's verdict is a 6VAJ "
          "number. Nothing may rank on it yet.")


if __name__ == "__main__":
    main()
