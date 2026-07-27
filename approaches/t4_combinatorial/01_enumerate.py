"""
Purpose: T_4 stages 1-5b — enumerate warhead x R-group, verify, label, gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: seeds.yaml core, warhead_classes_3.csv, rgroup_library_1.csv
Output: D4.parquet (full frame, rejects RETAINED) + manifest

Implements Rev 3 section 6 steps 1-5b:

  1  define the core and its two labelled attachment points
  2  enumerate combinatorially ON THE MOLECULAR GRAPH
  3  verify the intact core in EVERY product
  4  RDKit descriptors
  5  two-tier structural alerts (whole = advisory, R-group = gating)
  5b warhead-electrophilicity validity gate

GRAPH, NEVER STRING. Coupling by splicing SMILES text silently produces
malformed molecules whenever a fragment contains a ring, because SMILES
ring-closure digits are positional rather than compositional. This was a real
bug in the prior project run, not a hypothetical. Everything here goes through
the RDKit reaction engine.

VERIFY EVERY ROW, NOT A SAMPLE. A reaction definition correct in general can
still misfire on an unusual substrate, and a product that lost the core is not a
member of this library at all.

STAMP, DO NOT DELETE. Failing a gate sets `rejected_at` and skips the expensive
downstream tiers. The row stays. Only step 3 (core verification) removes
anything, because a product without the core is not a candidate.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import alerts as al                      # noqa: E402
from shared import descriptors as desc               # noqa: E402
from shared import io as dio                         # noqa: E402
from shared import novelty as nov                    # noqa: E402
from shared import smiles as smi                     # noqa: E402
from shared import warhead_library as wl             # noqa: E402

log = logging.getLogger("t4-enumerate")

EXPERIMENT = "04_t4_combinatorial"
RGROUPS = REPO / "data" / "reference" / "rgroup_library_1.csv"

# The coupling reaction from the spec: join an attachment point on one fragment
# to an attachment point on the other, dropping both dummies.
COUPLE = rdChemReactions.ReactionFromSmarts("[*:1]-[#0].[#0]-[*:2]>>[*:1]-[*:2]")


def load_core() -> tuple[str, str]:
    """Return (core SMILES with labelled dummies, alert-scoping SMARTS)."""
    cfg = yaml.safe_load((REPO / "config" / "seeds.yaml").read_text(encoding="utf-8"))
    pc = cfg["seeds"]["sulfopin"]["protected_core"]
    return pc["smiles"], pc["alert_scoping_smarts"]


def couple(core_smiles: str, warhead_frag: str, rgroup_frag: str) -> str | None:
    """Attach warhead at [1*] and R-group at [2*], on the graph.

    Returns
    -------
    str or None
        Canonical SMILES of the product, or None if either coupling failed.
    """
    core = Chem.MolFromSmiles(core_smiles)
    wh = Chem.MolFromSmiles(warhead_frag)
    rg = Chem.MolFromSmiles(rgroup_frag)
    if core is None or wh is None or rg is None:
        return None

    # Attach at the isotope-labelled dummies in order: 1 = warhead, 2 = R-group.
    for isotope, frag in ((1, wh), (2, rg)):
        target = None
        for atom in core.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetIsotope() == isotope:
                target = atom.GetIdx()
                break
        if target is None:
            return None
        # Zero the isotope so the generic reaction template matches this dummy
        # and not the other one still waiting its turn.
        rw = Chem.RWMol(core)
        rw.GetAtomWithIdx(target).SetIsotope(0)
        try:
            products = COUPLE.RunReactants((rw.GetMol(), frag))
        except Exception:  # noqa: BLE001 - a failed coupling is a None, not a crash
            return None
        if not products:
            return None
        core = products[0][0]
        try:
            Chem.SanitizeMol(core)
        except Exception:  # noqa: BLE001
            return None
    return smi.canonical(Chem.MolToSmiles(core))


def enumerate_library(allow_statuses: tuple[str, ...]) -> pd.DataFrame:
    """Full cross-product of enumerable warheads x R-groups."""
    core_smiles, _ = load_core()
    warheads = wl.enumerable(allow_statuses=allow_statuses)
    rgroups = pd.read_csv(RGROUPS)
    log.info("enumerating %d warhead classes x %d R-groups = %d products",
             len(warheads), len(rgroups), len(warheads) * len(rgroups))

    rows, failed = [], 0
    for _, w in warheads.iterrows():
        for _, r in rgroups.iterrows():
            product = couple(core_smiles, str(w["warhead_fragment_smiles"]),
                             str(r["fragment_smiles"]))
            if product is None:
                failed += 1
                continue
            rows.append({
                "canonical_smiles": product,
                "candidate_id": smi.candidate_id(product, prefix="t4"),
                "approach": "t4",
                "warhead_class": w["class_id"],
                "warhead_status": w["structure_status"],
                "warhead_mechanism": w["mechanism"],
                "rgroup_id": r["rgroup_id"],
                "linker_id": r["linker_id"],
                "aryl_smiles": r["aryl_smiles"],
            })
    if failed:
        log.warning("%d coupling(s) failed and were not enumerated", failed)
    df = pd.DataFrame(rows)
    before = len(df)
    df = df.drop_duplicates("canonical_smiles").reset_index(drop=True)
    if before != len(df):
        log.info("dropped %d duplicate product(s) — different (warhead, R-group) "
                 "pairs can give the same molecule", before - len(df))
    return df


def verify_core(df: pd.DataFrame, core_smarts: str) -> pd.DataFrame:
    """Step 3: confirm the intact core in EVERY product; drop any that lost it.

    This is the one step that removes rows. A product missing the core is not a
    member of this library, so retaining it would corrupt the frame rather than
    preserve information.
    """
    keep = [smi.has_substructure(s, core_smarts) for s in df["canonical_smiles"]]
    lost = len(df) - sum(keep)
    if lost:
        log.error("%d product(s) LOST THE CORE and were removed — check the "
                  "coupling reaction, not the data", lost)
    else:
        log.info("core verified intact in all %d products", len(df))
    out = df[pd.Series(keep, index=df.index)].reset_index(drop=True)
    out["core_verified"] = True
    return out


def warhead_validity_gate(df: pd.DataFrame) -> pd.DataFrame:
    """Step 5b: is the ATTACHED warhead still a reactive electrophile?

    The failure this blocks is invisible: in the prior run 6 of 16 warhead
    classes collapsed to an inert formamide or sulfonamide once coupled to the
    core, and docking happily returned plausible scores for molecules that
    could not react at all.
    """
    lib = wl.load()
    patterns = {str(r["class_id"]): Chem.MolFromSmarts(str(r["reactive_atom_smarts"]))
                for _, r in lib.iterrows()}
    ok = []
    for _, row in df.iterrows():
        patt = patterns.get(row["warhead_class"])
        m = smi.to_mol(row["canonical_smiles"])
        ok.append(bool(patt is not None and m is not None and m.HasSubstructMatch(patt)))
    out = df.copy()
    out["warhead_intact"] = ok

    dead = out.loc[~out["warhead_intact"], "warhead_class"].value_counts().to_dict()
    if dead:
        log.warning("warhead-validity gate: %d candidate(s) lost their electrophile "
                    "on coupling, by class: %s", int((~out["warhead_intact"]).sum()), dead)
    per_class = out.groupby("warhead_class")["warhead_intact"].mean()
    fully_dead = [c for c, frac in per_class.items() if frac == 0.0]
    if fully_dead:
        log.error("warhead class(es) DEAD in every product: %s. That is a "
                  "chemistry problem with the attachment, not a data problem.",
                  fully_dead)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-statuses", default=None,
                    help="comma-separated; defaults to the t4 config")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    cfg = yaml.safe_load(
        (REPO / "config" / "approaches" / "t4_combinatorial.yaml").read_text(encoding="utf-8"))
    allow = tuple(args.allow_statuses.split(",") if args.allow_statuses
                  else cfg["warheads"]["allow_statuses"])
    _, alert_core = load_core()

    df = enumerate_library(allow)
    df = verify_core(df, alert_core)

    log.info("computing descriptors")
    df = desc.compute_frame(df)
    log.info("computing novelty against the frozen external set")
    df = nov.novelty_frame(df)

    log.info("two-tier structural alerts (whole advisory, R-group gating)")
    df = al.screen_frame(df, core_smarts=alert_core)

    df = warhead_validity_gate(df)
    # 5b stamps only rows not already stamped by the alert gate — a candidate
    # keeps the reason it was FIRST set aside.
    stamp = (~df["warhead_intact"]) & df["rejected_at"].isna()
    df.loc[stamp, "rejected_at"] = "warhead_validity"

    survivors = df["rejected_at"].isna().sum()
    log.info("enumerated %d candidates; %d survive to covalent docking",
             len(df), int(survivors))

    out = dio.write_full_frame(
        df, approach="t4", experiment=EXPERIMENT, stage="t4_enumerate",
        params={"allow_statuses": list(allow), "library_size": len(df),
                "n_warhead_classes": int(df["warhead_class"].nunique()),
                "n_rgroups": int(df["rgroup_id"].nunique())},
        inputs={"warhead_library": REPO / "data/reference/warhead_classes_3.csv",
                "rgroup_library": RGROUPS})

    print(f"\nT_4 enumeration -> {out}")
    print(f"  library size        {len(df)}")
    print(f"  warhead classes     {df['warhead_class'].nunique()}")
    print(f"  R-groups            {df['rgroup_id'].nunique()}")
    print(f"  core verified       {int(df['core_verified'].sum())}/{len(df)}")
    print(f"  warhead intact      {int(df['warhead_intact'].sum())}/{len(df)}")
    print(f"  survive to docking  {int(survivors)}")
    print("\n  stamped rejections (retained in the frame):")
    for reason, n in df["rejected_at"].value_counts().items():
        print(f"    {reason:22s} {n}")
    print("\n  per warhead class:")
    for cls, g in df.groupby("warhead_class"):
        print(f"    {cls:22s} n={len(g):4d}  intact={int(g['warhead_intact'].sum()):4d}"
              f"  surviving={int(g['rejected_at'].isna().sum()):4d}")


if __name__ == "__main__":
    main()
