"""
Purpose: verify which Pin1 PDB entries are covalent, TO WHAT, and with what warhead.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-01
Input: pdb_entry_survey_<latest>.csv + pdb_ligand_survey_<latest>.csv (blacksmith)
Output: 00_outputs/blacksmith/pdb_covalent/covalent_links_<N>.csv + a printed report

#4 Phase 0.3a. The plan's wording is the point: **RCSB says "covalent" but not
to WHAT.** Pin1 has a second ligandable cysteine (Cys57,
10.1016/j.chembiol.2023.11.015), so an entry annotated covalent may not be a
Cys113 adduct at all, and counting it as one would inflate the chemotype
denominator with molecules that do not attack the site T_3 and T_4 target.

This resolves that from `_struct_conn` in each entry's mmCIF header, which
records every covalent linkage as an explicit atom-to-atom pair.

WHY THE ADDUCT PATTERN, NOT THE REACTIVE ONE. A deposited covalent ligand is
the PRODUCT: the warhead has already reacted and its leaving group is gone. A
chloroacetamide adduct has no chlorine left, so matching
`reactive_atom_smarts` against a deposited SMILES would fail on exactly the
molecules we are trying to classify. `warhead_classes` carries
`adduct_attachment_smarts` for this, and it is what D0022 ("dock the adduct,
not the pre-reaction ligand") and D0030 (the adduct saturates for acrylamide
but not for quinone) already established as the right frame.

WHAT THIS DOES NOT DO. It does not invent a class for a ligand that matches
nothing. `canonical_class()` is a PROSE canonicaliser that raises on anything
unmapped, deliberately, and there is no structural classifier for chemistries
the library has never seen — which is precisely the four the PDB survey
suggested (aryl aldehyde, maleate ester, SuFEx, Mannich). Unmatched adducts are
reported as UNCLASSIFIED with their SMILES so a chemist can name them. Guessing
here would manufacture the chemotype count this exercise exists to measure
honestly, which is the D0045 failure in a new costume.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import outputs as sout                  # noqa: E402
from shared import warhead_library as wl            # noqa: E402

log = logging.getLogger("curate-covalent")

SURVEY = sout.Topic("blacksmith", "redock_pin1")
OUT = sout.Topic("blacksmith", "pdb_covalent")

HEADER_URL = "https://files.rcsb.org/header/{pdb}.cif"
CATALYTIC_CYS = 113
SECOND_CYS = 57
MAX_WORKERS = 8          # polite against RCSB; the header files are small

# Peptide capping and blocking groups. These ARE covalently bonded to the
# chain and do appear as `covale` links, but they are chemistry the
# crystallographer added to terminate a peptide, not a candidate warhead.
# Counting ACE as a covalent chemotype would inflate the denominator D0045
# exists to keep honest.
CAPPING_COMPS = {"ACE", "NH2", "NME", "FOR", "TFA"}


def fetch_header(pdb: str, retries: int = 3) -> str | None:
    """The entry's mmCIF HEADER, which carries `_struct_conn` without coordinates."""
    req = urllib.request.Request(HEADER_URL.format(pdb=pdb.lower()),
                                 headers={"User-Agent": "inhibition/curate"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as fh:
                return fh.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:          # entry has no header; not an error
                return None
            log.debug("%s attempt %d: HTTP %s", pdb, attempt + 1, exc.code)
        except Exception as exc:  # noqa: BLE001 - transient network
            log.debug("%s attempt %d: %s", pdb, attempt + 1, exc)
        time.sleep(1.5 * (attempt + 1))
    log.warning("%s: header unavailable after %d attempts", pdb, retries)
    return None


def parse_struct_conn(cif: str) -> list[dict]:
    """Every `covale` linkage in the header, as partner pairs.

    TWO SERIALISATIONS, BOTH REQUIRED. mmCIF writes a category with several
    records as a `loop_`, and a category with EXACTLY ONE record as plain
    `_tag value` pairs. Handling only the loop form silently returns zero
    links for every single-adduct entry — which is most of them, including
    6VAJ, this project's own receptor. The first version of this function did
    exactly that and reported 0 covalent entries across the survey, which looks
    like a finding rather than a parse failure.

    Tags are read and rows mapped onto them by NAME. Reading by position would
    break the moment RCSB reorders a column — catalogue entry #2 in
    `docs/how_this_project_breaks.md`.
    """
    lines = cif.splitlines()
    out: list[dict] = []

    # --- loop_ form: several records ---------------------------------------
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        j, tags = i + 1, []
        while j < len(lines) and lines[j].strip().startswith("_struct_conn."):
            tags.append(lines[j].strip().split(".", 1)[1])
            j += 1
        if not tags:
            i += 1
            continue
        while j < len(lines):
            row = lines[j].strip()
            if not row or row.startswith(("#", "loop_", "_")):
                break
            vals = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", row)
            if len(vals) == len(tags):
                out.append(dict(zip(tags, (v.strip("'\"") for v in vals))))
            j += 1
        i = j

    # --- key-value form: exactly one record --------------------------------
    if not out:
        rec: dict[str, str] = {}
        for line in lines:
            s = line.strip()
            if not s.startswith("_struct_conn."):
                continue
            parts = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", s)
            if len(parts) >= 2:
                rec[parts[0].split(".", 1)[1]] = parts[1].strip("'\"")
        if rec:
            out.append(rec)

    return [r for r in out if r.get("conn_type_id", "").lower() == "covale"]


def ligand_cys_links(rows: list[dict], pdb: str) -> list[dict]:
    """Covalent links joining a non-polymer component to a cysteine."""
    links = []
    for r in rows:
        for a, b in (("ptnr1", "ptnr2"), ("ptnr2", "ptnr1")):
            comp_a = r.get(f"{a}_label_comp_id", "")
            comp_b = r.get(f"{b}_label_comp_id", "")
            if comp_b != "CYS" or comp_a in ("CYS", "", "?"):
                continue
            if comp_a in CAPPING_COMPS:      # crystallographer's cap, not a warhead
                continue
            seq = r.get(f"{b}_auth_seq_id") or r.get(f"{b}_label_seq_id") or "?"
            links.append({
                "pdb_id": pdb.upper(),
                "comp_id": comp_a,
                "ligand_atom": r.get(f"{a}_label_atom_id", "?"),
                "cys_seq_id": seq,
                "cys_atom": r.get(f"{b}_label_atom_id", "?"),
                "dist_a": r.get("pdbx_dist_value", ""),
            })
            break
    return links


def classify_adduct(smiles: str, lib: pd.DataFrame) -> tuple[str | None, str]:
    """Match a deposited ligand against the library's warhead patterns.

    THE PDB COMPONENT SMILES IS THE FREE LIGAND, NOT THE ADDUCT. This function
    first matched `adduct_attachment_smarts`, on the reasoning that a
    crystallised covalent ligand is the reacted product. That is true of the
    COORDINATES and false of the SMILES: the chemical component dictionary
    defines each component standalone, and the bond to Cys is recorded only in
    `_struct_conn`. So QT7 — sulfopin, this project's own reference ligand —
    is deposited as `...C(=O)CCl` with its chlorine intact, and came back
    UNCLASSIFIED against the adduct pattern that expects the chlorine gone.

    That is worth stating plainly because it is the failure this repo
    catalogues: a populated, plausible result (11 "novel" chemistries!) that
    was computed from the wrong column, and which looked like a discovery.

    So `reactive_atom_smarts` is matched FIRST — the warhead as drawn — and the
    adduct pattern only as a fallback, for any source that does supply reacted
    structures. Both are reported so the basis of a call is never implicit.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    m = Chem.MolFromSmiles(str(smiles)) if smiles and str(smiles) != "nan" else None
    if m is None:
        return None, "unparseable or absent SMILES"

    for col, note in (("reactive_atom_smarts", "warhead-as-drawn"),
                      ("adduct_attachment_smarts", "adduct pattern")):
        hits = []
        for _, row in lib.iterrows():
            pat = row.get(col)
            if not pat or str(pat) in ("nan", wl.UNRESOLVED):
                continue
            q = Chem.MolFromSmarts(str(pat))
            if q is not None and m.HasSubstructMatch(q):
                hits.append(str(row["class_id"]))
        if hits:
            # Several classes share a pattern BY DESIGN (D0029: three warhead
            # classes give one identical adduct), so >1 hit is not ambiguity.
            return hits[0], f"{note}: {'; '.join(hits)}" if len(hits) > 1 else note
    return None, "no library warhead or adduct pattern matches"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None,
                    help="only the first N entries (smoke testing)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s %(name)s: %(message)s")

    entries = pd.read_csv(SURVEY.latest("pdb_entry_survey", ".csv"))
    ligands = pd.read_csv(SURVEY.latest("pdb_ligand_survey", ".csv"))
    pdbs = list(entries["pdb_id"].astype(str).unique())
    if args.limit:
        pdbs = pdbs[:args.limit]
    log.info("resolving _struct_conn for %d Pin1 entries", len(pdbs))

    links: list[dict] = []
    missing = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for pdb, cif in zip(pdbs, ex.map(fetch_header, pdbs)):
            if cif is None:
                missing.append(pdb)
                continue
            links.extend(ligand_cys_links(parse_struct_conn(cif), pdb))
    log.info("%d ligand-cysteine covalent link(s) across %d entries; "
             "%d header(s) unavailable",
             len(links), len({l["pdb_id"] for l in links}), len(missing))

    df = pd.DataFrame(links)
    if df.empty:
        raise SystemExit("no ligand-cysteine covalent links found")

    smi = (ligands.dropna(subset=["smiles"])
           .drop_duplicates("comp_id").set_index("comp_id")["smiles"])
    df["smiles"] = df["comp_id"].map(smi)
    df["cys_seq_id"] = pd.to_numeric(df["cys_seq_id"], errors="coerce").astype("Int64")
    df["site"] = df["cys_seq_id"].map(
        lambda s: "Cys113 (catalytic)" if s == CATALYTIC_CYS
        else (f"Cys{SECOND_CYS} (second site)" if s == SECOND_CYS
              else (f"Cys{s}" if pd.notna(s) else "unknown")))

    lib = wl.load()
    cls = df["smiles"].map(lambda s: classify_adduct(s, lib))
    df["warhead_class"] = [c for c, _ in cls]
    df["classify_note"] = [n for _, n in cls]

    dest = OUT.write("covalent_links", ".csv")
    df.sort_values(["site", "pdb_id"]).to_csv(dest, index=False)

    # ---- the report -------------------------------------------------------
    print(f"\nCOVALENT LINKAGES IN THE Pin1 PDB  ({len(pdbs)} entries surveyed)")
    print(f"  headers unavailable: {len(missing)}"
          + (f" — {', '.join(missing[:8])}" if missing else ""))
    print(f"\n  ligand-cysteine covalent links: {len(df)} "
          f"across {df['pdb_id'].nunique()} entries, "
          f"{df['comp_id'].nunique()} distinct ligands\n")
    print("  BY SITE — this is what 'covalent' alone does not tell you:")
    for site, g in df.groupby("site"):
        print(f"    {site:24s} {len(g):3d} link(s), "
              f"{g['comp_id'].nunique():3d} distinct ligand(s), "
              f"{g['pdb_id'].nunique():3d} entries")

    cat = df[df["cys_seq_id"] == CATALYTIC_CYS]
    print(f"\n  CHEMOTYPES AT Cys113, counted by WARHEAD CLASS (D0045):")
    known = cat[cat["warhead_class"].notna()]
    n_classes = known["warhead_class"].nunique()
    for c, g in known.groupby("warhead_class"):
        print(f"    {c:26s} {g['comp_id'].nunique():3d} distinct ligand(s)")
    unk = cat[cat["warhead_class"].isna()]
    print(f"\n    mapped to a library class : {n_classes} class(es), "
          f"{known['comp_id'].nunique()} ligand(s)")
    print(f"    UNCLASSIFIED              : {unk['comp_id'].nunique()} "
          "distinct ligand(s) — the library has no adduct pattern for these")
    print(f"\n  VERDICT: {n_classes} chemotype(s) against D0045's floor of 6 — "
          f"{'CLEARS' if n_classes >= 6 else 'UNDERPOWERED'}")
    if len(unk):
        print("\n  Unclassified adducts, for a chemist to name (comp_id: SMILES):")
        for cid, g in unk.groupby("comp_id"):
            print(f"    {cid:6s} {str(g['smiles'].iloc[0])[:88]}")
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
