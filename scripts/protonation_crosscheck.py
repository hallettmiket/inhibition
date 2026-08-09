"""
Purpose: check the pH 7.4 species we dock against a literature protonation model.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-08
Input: candidate SMILES from the current T_3/T_4 frames
Output: 00_outputs/blacksmith/protonation/crosscheck_<N>.csv

@tt8804 chose the cross-check over replacing the model, and that is the right way
round. `obabel -p 7.4` stays operative because it is what actually protonated the
structure that was DOCKED -- swap it and every pose on disk describes a different
species from the one being ranked. This does not change the pipeline's chemistry;
it says where that chemistry is likely wrong.

WHY DIMORPHITE CANNOT SIMPLY REPLACE IT. Dimorphite-DL (Ropp, Kaminsky,
Yablonski & Durrant, *J Cheminform* 2019;11:14) is the standard tool, but it is a
state ENUMERATOR: at pH 7.4 it returns every microstate whose predicted pKa
interval spans the pH, not one answer. Measured here, `t4_8b12474d07bd` comes back
as four states spanning charge -1 to +1, and Sulfopin comes back with its AMIDE
nitrogen protonated -- a species with a conjugate-acid pKa near -0.5 that does not
exist at physiological pH. A pipeline needs exactly one species per molecule; an
enumerator hands it a set, and picking from that set by rule is how a protonation
model quietly becomes a tuning knob.

WHAT THIS REPORTS

  agrees        obabel's species is among dimorphite's states, and its charge
                matches dimorphite's most common charge. Nothing to do.
  charge_differs obabel's charge is not the one dimorphite favours. The molecule
                was docked, scored and simulated as a species the literature
                model considers less likely. These are the rows worth a chemist's
                eye -- especially where |charge| >= 1, since charge drives both
                the docking electrostatics and the ion count in the MD box.
  absent        obabel's species is not among dimorphite's states at all.

NO VERDICT IS COMPUTED. Disagreement is not proof either tool is right; both are
predictors. The output is a worklist for a human, not a gate.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import ionisation as ion                  # noqa: E402
from shared import outputs as sout                    # noqa: E402

log = logging.getLogger("protonation-xcheck")
FRAMES = {
    "T3": "/data/lab_vm/append_only/inhibition/03_t3_reinvent",
    "T4": "/data/lab_vm/append_only/inhibition/04_t4_combinatorial",
}


def _latest(d: str, stem: str) -> Path | None:
    import glob
    import re
    fs = glob.glob(f"{d}/{stem}_*.parquet")
    if not fs:
        return None
    return Path(max(fs, key=lambda p: int(re.search(r"_(\d+)\.parquet$", p).group(1))))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default=None,
                    help="file of idents, one per line — e.g. the re-sweep list")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from dimorphite_dl import protonate_smiles

    frames = []
    for tier, d in FRAMES.items():
        f = _latest(d, "D3" if tier == "T3" else "D4")
        if f is None:
            continue
        df = pd.read_parquet(f)[["candidate_id", "canonical_smiles"]].copy()
        df["tier"] = tier
        frames.append(df)
        log.info("%s: %s (%d rows)", tier, f.name, len(df))
    d = pd.concat(frames, ignore_index=True).dropna(subset=["canonical_smiles"])
    d = d.drop_duplicates("candidate_id")

    if args.only:
        want = {ln.strip().split()[0] for ln in Path(args.only).read_text().splitlines()
                if ln.strip()}
        d = d[d.candidate_id.isin(want)]
        log.info("--only: %d of %d requested idents present", len(d), len(want))
    if args.limit:
        d = d.head(args.limit)
    log.info("cross-checking %d molecules", len(d))

    # obabel in ONE batch -- it is a subprocess per call otherwise, and 5,700
    # process spawns dominate the runtime of a job that is otherwise arithmetic.
    obab = ion.protonate(dict(zip(d.candidate_id, d.canonical_smiles)))

    rows = []
    for i, r in enumerate(d.itertuples(), 1):
        if i % 500 == 0:
            log.info("  %d/%d", i, len(d))
        o = obab.get(r.candidate_id)
        om = Chem.MolFromSmiles(o) if isinstance(o, str) else None
        if om is None:
            rows.append({"candidate_id": r.candidate_id, "tier": r.tier,
                         "verdict": "obabel_failed"})
            continue
        oq = Chem.GetFormalCharge(om)
        ocan = Chem.MolToSmiles(om)
        try:
            states = protonate_smiles(r.canonical_smiles, ph_min=7.4, ph_max=7.4)
        except Exception:                              # noqa: BLE001
            rows.append({"candidate_id": r.candidate_id, "tier": r.tier,
                         "obabel_charge": oq, "verdict": "dimorphite_failed"})
            continue
        qs, cans = [], set()
        for s in states:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            qs.append(Chem.GetFormalCharge(m))
            cans.add(Chem.MolToSmiles(m))
        if not qs:
            rows.append({"candidate_id": r.candidate_id, "tier": r.tier,
                         "obabel_charge": oq, "verdict": "dimorphite_failed"})
            continue
        fav = max(set(qs), key=qs.count)               # dimorphite's modal charge
        verdict = ("agrees" if ocan in cans and oq == fav
                   else "charge_differs" if ocan in cans or oq != fav
                   else "absent")
        rows.append({"candidate_id": r.candidate_id, "tier": r.tier,
                     "obabel_charge": oq, "dimorphite_modal_charge": fav,
                     "dimorphite_n_states": len(qs),
                     "dimorphite_charges": ",".join(map(str, sorted(set(qs)))),
                     "species_in_dimorphite_set": ocan in cans,
                     "verdict": verdict})

    t = pd.DataFrame(rows)
    dest = sout.Topic("blacksmith", "protonation").write("crosscheck", ".csv")
    t.to_csv(dest, index=False)

    print("\n" + "=" * 70)
    print("  pH 7.4 species: obabel (operative) vs Dimorphite-DL (reference)")
    print("=" * 70)
    print(t.verdict.value_counts().to_string())
    if "obabel_charge" in t.columns:
        bad = t[t.verdict == "charge_differs"]
        print(f"\n  charge disagreements: {len(bad)}")
        if len(bad):
            print(bad.groupby(["obabel_charge", "dimorphite_modal_charge"])
                    .size().rename("n").to_string())
    print(f"\n  -> {dest}")
    print("  No verdict is applied. Both are predictors; this is a worklist.\n")


if __name__ == "__main__":
    main()
