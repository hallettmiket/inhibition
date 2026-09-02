"""
Purpose: add a HAND-DESIGNED analogue to a candidate library, annotated by the same code the enumeration uses.
Author: @twu383 (with Claude Code)
Date: 2026-08-31
Input: a parent candidate_id already in the frame + the analogue's SMILES
Output: the next D<N> frame, parent frame + one row, with its manifest

WHY THIS EXISTS. T_4 is a combinatorial enumeration: every candidate in it is a
warhead crossed with an R-group from `rgroup_library`. A medicinal chemist
proposing "same molecule, sulfone swapped for gem-difluoro" is proposing a
molecule the enumeration cannot produce, and there was no way to get it into the
screen except by hand-editing a parquet -- which is how a row with a plausible
but uncomputed number gets onto disk.

THE ROW STARTS EMPTY, AND THAT IS THE WHOLE DESIGN. The obvious implementation
is `row = parent_row.copy()` then overwrite what changed. That is exactly the
defect this repo catalogues: `affinity_kcal`, `engagement`, `class_rank`,
`nac3_*` and forty other columns would silently carry the PARENT's docking
result on the ANALOGUE's row, every one of them populated and plausible, and
nothing downstream would ever disagree. So the row is built as all-NA over the
parent frame's schema and each column is filled only by code that actually
computed it for THIS molecule. A column this script forgets is NA, which reads
as "not measured" -- the failure mode that announces itself.

WHAT IS INHERITED, EXPLICITLY. Only facts about the WARHEAD, which is unchanged
by construction and is verified unchanged before anything is written:
`warhead_class`, `warhead_mechanism`, `warhead_status`, `adduct_class`,
`adduct_attachment_smarts`, `leaving_group_smiles`. If the analogue's warhead
does not match the parent's, this refuses to run rather than annotating a
molecule as something it is not.

WHAT IS RECOMPUTED: candidate_id (from the InChIKey, never assigned), every
descriptor, novelty, both alert tiers, the pH 7.4 docked species, the covalent
adduct, the synthesizability rules, and the xTB LUMO -- because LUMO is a
whole-molecule property and swapping a sulfone for CF2 moves it, so inheriting
the parent's `reactivity_in_window` would be a reactivity claim nobody made.

WHAT IS LEFT NA: everything produced by docking, ranking, MMGBSA, NAC screening
or MD. Those are properties of a RUN, and this molecule has not been in one.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "approaches" / "t4_combinatorial"))

from shared import alerts as al                      # noqa: E402
from shared import covalent_adduct as ca             # noqa: E402
from shared import descriptors as desc               # noqa: E402
from shared import io as dio                         # noqa: E402
from shared import ionisation as ion                 # noqa: E402
from shared import novelty as nov                    # noqa: E402
from shared import smiles as smi                     # noqa: E402
from shared import synthesizability as syn           # noqa: E402
from shared import warhead_library as wl             # noqa: E402

log = logging.getLogger("add-analogue")

#: Facts about the warhead, which this script requires to be UNCHANGED. Anything
#: not listed here is either recomputed or left NA.
WARHEAD_INHERITED = (
    "warhead_class", "warhead_mechanism", "warhead_status",
    "adduct_class", "adduct_attachment_smarts", "leaving_group_smiles",
)

#: Where the frame for an approach lives, and the stem its versions carry.
FRAME_DIR = {"t4": ("04_t4_combinatorial", "D4"),
             "t3": ("03_t3_reinvent", "D3")}


def _warhead_match_count(smiles: str, reactive_smarts: str) -> int:
    """How many times the warhead's reactive group appears.

    A count rather than a boolean: an analogue that accidentally introduces a
    SECOND copy of the reactive group is a different chemical proposition from
    the one the chemist drew, and `to_adduct_form` would pick one of them by
    position.
    """
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    q = Chem.MolFromSmarts(reactive_smarts)
    if m is None or q is None:
        raise SystemExit(f"unparseable molecule or warhead SMARTS: {smiles!r}")
    return len(m.GetSubstructMatches(q))


def build_row(parent: pd.Series, analogue_smiles: str, *, note: str,
              core_smarts: str) -> dict:
    """One fully annotated candidate row, built from nothing but the SMILES."""
    from rdkit import Chem

    m = Chem.MolFromSmiles(analogue_smiles)
    if m is None:
        raise SystemExit(f"RDKit cannot parse the analogue: {analogue_smiles!r}")
    canon = Chem.MolToSmiles(m)

    cls = str(parent["warhead_class"])
    lib = wl.load()
    meta = lib[lib.class_id == cls]
    if meta.empty:
        raise SystemExit(f"warhead class {cls!r} is not in {wl.DEFAULT_LIBRARY}")
    reactive = str(meta.iloc[0].reactive_atom_smarts)

    # THE WARHEAD MUST SURVIVE THE EDIT. This is the one thing the caller is
    # asserting by naming a parent at all, so it is checked rather than trusted.
    n_par = _warhead_match_count(str(parent["canonical_smiles"]), reactive)
    n_ana = _warhead_match_count(canon, reactive)
    if n_ana != n_par:
        raise SystemExit(
            f"warhead count changed: parent has {n_par} match(es) of "
            f"{reactive!r}, analogue has {n_ana}. This is not the same warhead, "
            f"so it cannot inherit {cls!r}. Refusing.")

    cid = smi.candidate_id(canon, prefix="t4")
    if cid is None:
        raise SystemExit("analogue has no InChIKey and cannot be keyed")

    row: dict = {
        "canonical_smiles": canon,
        "candidate_id": cid,
        "approach": str(parent["approach"]),
        "inchikey": smi.inchikey(canon),
        "inchikey_ok": True,
        "core_verified": bool(smi.has_substructure(canon, core_smarts)),
        "molecule_complete": True,
        "warhead_intact": n_ana == n_par,
        # PROVENANCE, in the columns a reader of this frame already looks at.
        # `rgroup_id` is where a T_4 row records where its R-group came from,
        # and "this one came from a person" belongs in the same place rather
        # than in a column nothing reads.
        "rgroup_id": f"designed:{parent['candidate_id']}",
        "linker_id": str(parent["linker_id"]),
        "designed_note": note,
        "designed_parent": str(parent["candidate_id"]),
    }
    for c in WARHEAD_INHERITED:
        if c in parent.index:
            row[c] = parent[c]

    log.info("descriptors")
    row.update(desc.compute(canon))

    log.info("novelty against the frozen external set")
    row["novelty_external"] = nov.novelty(canon)

    log.info("two-tier structural alerts")
    # `to_columns()` rather than a hand-written mapping: the frame's alert
    # columns are DEFINED by this method, so a hand-rolled dict would be a
    # second, silently-drifting definition of the same schema.
    tt = al.two_tier(canon, core_smarts)
    row.update(tt.to_columns())

    # THE GATE CANNOT BE SCOPED TO A MOLECULE THAT NO LONGER HAS THE CORE, AND
    # `False` IS THE WRONG WAY TO SAY SO.
    #
    # `two_tier` returns `passes_gate=False, reason="core not found"` when the
    # core is absent, and for its own caller that is right -- an ENUMERATION
    # product that lost the core is a coupling bug. A hand-designed analogue
    # that modifies the core is a different thing entirely, and the two arrive
    # at the same `False`. Left alone, `screen_frame`'s rule
    # (`alert_gate_pass == False` -> `rejected_at = "alerts"`) would stamp this
    # molecule as carrying a dirty decoration, and a reader comparing it with
    # its parent would read a chemistry difference where there is none: the
    # WHOLE-MOLECULE alerts are computed on both and are directly comparable.
    #
    # So: NA, not False -- the convention this repo already adopted for T_1's
    # advisory alerts (`shared/alerts.screen_frame`, and state_of_the_project
    # section 6). `False` means tested and failed; nothing here was tested.
    if not row["core_verified"]:
        log.warning("core absent: this analogue MODIFIES the core, so the "
                    "two-tier alert gate is not applicable -- recording NA, "
                    "not False. Whole-molecule alerts are still computed.")
        row["alert_gate_pass"] = pd.NA
        row["alert_gate_applied"] = False
        row["alert_gate_reason"] = (
            "core modified by design; the two-tier gate cannot be scoped. "
            "Whole-molecule alerts are reported and ARE comparable.")
        # Attribution columns are NA, not 0. `core_alert_total = 0` would read
        # as "no alerts on the core"; nothing was attributed to anything.
        for c in ("core_alert_total", "core_alert_names", "rgroup_alert_total",
                  "rgroup_alert_names", "boundary_alert_total",
                  "boundary_alert_names", "rgroup_smiles"):
            row[c] = pd.NA

    log.info("pH 7.4 docked species")
    prot = ion.protonate({cid: canon})
    dsmi = prot.get(cid)
    row["docked_smiles"] = dsmi
    row["docked_species_ok"] = dsmi is not None
    row["charge_ph74"] = ion.charge_from_smiles(dsmi) if dsmi else None
    row["charge_class"] = ion.charge_class(row["charge_ph74"])
    row["has_phosphate"] = ion.has_phosphate(canon)

    log.info("covalent adduct")
    try:
        add = ca.to_adduct_form(canon, cls)
        row.update(add.as_dict())
    except ca.AdductError as e:                        # noqa: BLE001
        # Loud and recorded, never silently skipped: an analogue whose adduct
        # cannot be built is one the covalent docking cannot score.
        log.error("adduct could not be built: %s", e)
        row["adduct_smiles"] = None

    log.info("synthesizability rules")
    viol = syn.violations(canon)
    row["synth_fail"] = bool(viol)
    row["synth_violations"] = "|".join(v.name for v in viol) if viol else None

    log.info("xTB LUMO (whole-molecule -- NOT inherited from the parent)")
    row.update(_reactivity(canon, int(row["charge_ph74"] or 0)))

    row["rejected_at"] = None
    return row


def _reactivity(canon: str, charge: int) -> dict:
    """LUMO and the window verdict, computed for THIS molecule."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "t4_reactivity", REPO / "approaches/t4_combinatorial/02_reactivity_triage.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    r = mod.lumo_ev(canon, charge=charge)
    if r is None:
        log.error("xTB failed -- reactivity left NA rather than assumed in-window")
        return {}
    win = mod.anchor_window()
    inside = win["window_lo"] <= r["lumo_ev"] <= win["window_hi"]
    log.info("LUMO %+.3f eV; window [%.3f, %.3f] -> %s",
             r["lumo_ev"], win["window_lo"], win["window_hi"],
             "in window" if inside else "OUTSIDE")
    return {"lumo_ev": r["lumo_ev"], "homo_lumo_gap_ev": r["gap_ev"],
            "reactivity_in_window": inside,
            "reactivity_flag": "in_window" if inside else "outside_window"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--approach", default="t4", choices=sorted(FRAME_DIR))
    ap.add_argument("--parent", required=True,
                    help="candidate_id whose warhead the analogue keeps")
    ap.add_argument("--smiles", required=True, help="the analogue")
    ap.add_argument("--note", required=True,
                    help="why this molecule exists, in one line")
    ap.add_argument("--dry-run", action="store_true",
                    help="annotate and print; write nothing")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    subdir, stem = FRAME_DIR[args.approach]
    df, src = dio.latest_frame(subdir, args.approach)
    log.info("parent frame %s (%d rows, %d columns)", src.name, len(df), df.shape[1])

    par = df[df.candidate_id == args.parent]
    if par.empty:
        raise SystemExit(f"parent {args.parent!r} is not in {src.name}")
    parent = par.iloc[0]

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "t4_enum", REPO / "approaches/t4_combinatorial/01_enumerate.py")
    enum = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(enum)
    _, core_smarts = enum.load_core()

    row = build_row(parent, args.smiles, note=args.note, core_smarts=core_smarts)

    if row["candidate_id"] in set(df.candidate_id):
        raise SystemExit(f"{row['candidate_id']} is already in {src.name} -- "
                         "this molecule is not new")

    # THE ALL-NA BASE. Every column the parent frame has, none of its values.
    blank = {c: pd.NA for c in df.columns}
    blank.update({k: v for k, v in row.items()})
    new = pd.DataFrame([blank])
    for c in df.columns:                       # keep the parent frame's dtypes
        if c in new and df[c].dtype != object:
            new[c] = new[c].astype(df[c].dtype, errors="ignore")

    filled = [c for c in df.columns if c in row and row[c] is not None]
    left_na = [c for c in df.columns if c not in row]
    print(f"\n  analogue  {row['candidate_id']}")
    print(f"  smiles    {row['canonical_smiles']}")
    print(f"  parent    {args.parent}  ({parent['canonical_smiles']})")
    print(f"  warhead   {row.get('warhead_class')} (verified unchanged)")
    print(f"  filled    {len(filled)} columns computed for this molecule")
    print(f"  left NA   {len(left_na)} columns -- every docking / ranking / MD "
          f"quantity, because it has not been in a run")
    for c in ("MW", "HAC", "cLogP", "TPSA", "QED", "SAscore", "lumo_ev",
              "reactivity_in_window", "alert_gate_pass", "synth_fail",
              "charge_ph74", "docked_species_ok"):
        if c in row:
            print(f"    {c:24s} {row[c]}")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return

    out = pd.concat([df, new], ignore_index=True)
    dest = dio.write_full_frame(
        out, approach=args.approach, experiment=subdir, stage="designed_analogue",
        params={"parent": args.parent, "analogue": row["canonical_smiles"],
                "candidate_id": row["candidate_id"], "note": args.note,
                "columns_filled": len(filled), "columns_left_na": len(left_na)},
        inputs={"parent_frame": src, "warhead_library": wl.DEFAULT_LIBRARY})
    print(f"\n  wrote {dest} ({len(out)} rows, was {len(df)})")


if __name__ == "__main__":
    main()
