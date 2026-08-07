"""
Purpose: put the known Pin1 binders through the SAME criterion as the candidates, so the GUI can compare against them.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: data/reference/pin1_reference_binders_4.csv + the v2 screen machinery
Output: 00_outputs/blacksmith/nac_v2/{agg,poses}_ref_<N>.csv + poses, marked is_reference

@tt8804, #19: "parent molecules need to be carried through all stages
automatically and highlighted on the gui for easy comparison."

WARHEAD CLASS IS ASSIGNED BY SMARTS, NOT BY THE NAME IN THE FILE. The reference
table's `warhead_class` column is free prose written for a human -- "1;4-
naphthoquinone (Michael acceptor)", "2-chloro-5-nitropyrimidine (SNAr
electrophile...)" -- and only 2 of 22 rows happen to match a canonical class_id
as a string. Hand-mapping that prose onto class IDs would be a value taken by a
label, which is the defect class this project keeps paying for (D0067 mapped a
mechanism NAME onto the wrong geometry and read 374 candidates as dead).

So each reference molecule is tested against every class's
`reactive_atom_smarts` and assigned the classes that actually MATCH its
structure. A molecule matching several is reported as such rather than silently
resolved, and a molecule matching none is reported as unscoreable rather than
forced into a class.

THESE ARE NOT A VALIDATION SET AND ARE NOT USED AS ONE. @tt8804 has ruled that
the known Pin1 binders are too few and too poor to decide which chemistry to
pursue, and no gate here is calibrated against them. They are on the screen for
one reason: so a chemist looking at a candidate can see what the incumbent
compound scores on the identical measurement. Sulfopin is the one that matters --
it is the parent.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from shared import outputs as sout                # noqa: E402
from shared import receptors as R                 # noqa: E402
import nac_screen as ns                           # noqa: E402
import nac_screen_v2 as v2                        # noqa: E402

log = logging.getLogger("screen-refs")
OUT = sout.Topic("blacksmith", "nac_v2")
REF = REPO / "data/reference/pin1_reference_binders_4.csv"


def assign_classes(smiles: str, wh: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Every warhead class whose reactive SMARTS matches this molecule."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = ns.largest_fragment(smiles)
    if mol is None:
        return []
    hits = []
    for r in wh.itertuples():
        patt = Chem.MolFromSmarts(r.reactive_atom_smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            hits.append((r.class_id, r.mechanism, r.reactive_atom_smarts))
    return hits


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--gpu", default="6")
    ap.add_argument("--nrun", type=int, default=200)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    R.resolve_3ikd_ian()
    refs = pd.read_csv(REF)
    wh = pd.read_csv(REPO / "data/reference/warhead_classes_10.csv")

    rows, cands = [], []
    for r in refs.itertuples():
        smi = getattr(r, "canonical_smiles", None)
        rec = {"name": r.name, "tier": r.tier, "potency": r.potency,
               "mechanism_declared": r.mechanism, "pdb": getattr(r, "pdb", None),
               "warhead_text": r.warhead_class}
        if not isinstance(smi, str) or smi == "UNVERIFIED":
            rec["assigned"] = ""
            rec["note"] = "no verified SMILES"
            rows.append(rec)
            continue
        hits = assign_classes(smi, wh)
        rec["assigned"] = "|".join(h[0] for h in hits)
        rec["n_classes_matched"] = len(hits)
        if not hits:
            rec["note"] = "no warhead SMARTS matches — not scoreable by this criterion"
        elif len(hits) > 1:
            rec["note"] = f"matches {len(hits)} classes; scored under each"
        else:
            rec["note"] = ""
        rows.append(rec)
        for cls, mech, smarts in hits:
            ident = f"ref_{r.name.replace(' ', '-')}__{cls}"
            cands.append(ns.Candidate(ident=ident, smiles=smi, warhead_class=cls,
                                      mechanism=mech, reactive_smarts=smarts,
                                      label="positive"))

    tbl = pd.DataFrame(rows)
    print("\n=== warhead classes assigned BY SMARTS, not by the file's prose ===\n")
    print(f"  {'name':<38}{'declared (prose)':<30}{'assigned (SMARTS)':<26}")
    for r in tbl.itertuples():
        txt = str(r.warhead_text)[:28]
        print(f"  {r.name[:37]:<38}{txt:<30}{r.assigned or '—':<26}")
    n_ok = sum(1 for r in tbl.itertuples() if r.assigned)
    print(f"\n  {n_ok} of {len(tbl)} reference molecules carry a warhead this "
          f"criterion can score\n  ({len(cands)} molecule×class jobs)")
    for r in tbl.itertuples():
        if getattr(r, "note", ""):
            print(f"    {r.name}: {r.note}")

    if args.dry_run or not cands:
        return

    rec_dir = ns.build_reactive_receptor(ns.RX_RECEPTOR)
    plain = sout.latest_path("blacksmith", "receptor_3ikd", "3IKD_prepared", ".pdbqt")
    log.info("screening %d reference jobs on GPU %s", len(cands), args.gpu)

    pose_buf, agg_buf = [], []
    for i, c in enumerate(cands, 1):
        prow, agg = v2.one(c, rec_dir, Path(plain), args.nrun, args.gpu, True)
        agg["is_reference"] = True
        agg_buf.append(agg)
        if not prow.empty:
            pose_buf.append(prow)
        log.info("[%d/%d] %s: %s", i, len(cands), c.ident, agg["status"])

    if pose_buf:
        pd.concat(pose_buf, ignore_index=True).to_csv(
            OUT.write("poses_ref", ".csv"), index=False)
    a = pd.DataFrame(agg_buf)
    a.to_csv(OUT.write("agg_ref", ".csv"), index=False)

    meta = OUT.write("reference_assignment", ".csv")
    tbl.to_csv(meta, index=False)
    ok = a[a.status == "ok"]
    print(f"\n=== {len(ok)} of {len(a)} scored ===")
    if len(ok):
        show = ok[["ident", "warhead_class", "n_in_range", "n_viable", "enrichment"]]
        print(show.to_string(index=False))
    print(f"\n  -> {meta}")
    print("\n  NOT a validation set. On the screen so a chemist can see what the "
          "incumbent\n  scores on the identical measurement — Sulfopin above all, "
          "it is the parent.")


if __name__ == "__main__":
    main()
