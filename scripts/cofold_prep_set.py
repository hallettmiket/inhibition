"""
Purpose: build the Boltz-2 benchmark set from deposited Pin1 complexes.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-07
Input: PDB IDs (fetched as mmCIF from RCSB) + a held-out/in-training split
Output: a CSV for cofold_bench.py + one truth PDB per entry

Turns the era split in `docs/prereg_cofolding.md` into files. Three things here
are decisions, not plumbing:

WHICH LIGAND IS THE INHIBITOR. Half these entries carry citrate, PEG or a second
copy of a fragment alongside the compound of interest. Picking "the HETATM group"
would benchmark Boltz on citrate. The inhibitor is chosen as the ligand COVALENTLY
BONDED to a cysteine where struct_conn records one -- these are covalent binders
and that bond is the ground truth -- and otherwise the largest non-additive by
heavy-atom count. The rule used is recorded per row so a wrong pick is visible
rather than silent.

THE SMILES IS THE DEPOSITED CHEMISTRY, POST-REACTION. A covalent ligand is
deposited as the adduct, still carrying the bond to Cys113 and usually missing its
leaving group. Boltz predicts NON-covalent complexes, so it is given the ligand as
deposited and asked to place it; the comparison is of position, not of chemistry.
This is the hard limit the prereg states -- co-folding cannot speak to reaction.

THE ERA SPLIT IS DATA, NOT ASSUMPTION. `held_out` comes from the deposition date
actually parsed out of the file, not from a hand-written list, and the date is
carried into the CSV so any row can be re-judged when Boltz-2's real cutoff is
confirmed.
"""

from __future__ import annotations

import argparse
import logging
import sys
import urllib.request
from pathlib import Path

import gemmi
import pandas as pd


def _is_aa(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def _is_water(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_water())

log = logging.getLogger("prep")

#: Crystallisation additives and cryoprotectants -- never the compound of interest.
ADDITIVES = {
    "HOH", "CIT", "1PG", "PEG", "PG4", "EDO", "GOL", "SO4", "PO4", "ACT",
    "DMS", "MPD", "TRS", "IMD", "CL", "NA", "MG", "ZN", "CA", "K", "IOD",
    "FLC", "MES", "EPE", "NO3", "ACY", "FMT", "P6G", "2PE", "12P", "PGE",
}

#: The 15 Pin1 complexes, split by deposition era (docs/prereg_cofolding.md).
HELD_OUT = ["9INN", "9INO", "9INP", "9INQ", "9JF6", "9JFH",
            "9V6G", "9V6I", "9V6P", "9V6W"]
IN_TRAINING = ["7EFJ", "7EFX", "7EKV", "7F0M", "6VAJ"]


def fetch(pdb_id: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{pdb_id}.cif"
    if not p.is_file() or p.stat().st_size < 1000:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        log.info("fetching %s", pdb_id)
        urllib.request.urlretrieve(url, p)
    return p


def deposition_date(doc) -> str:
    """The date the entry entered the PDB — the only basis for the era split.

    Read from the file rather than assumed from the ID, so a row can be
    re-judged against Boltz-2's real cutoff instead of against a guess.
    """
    v = doc.sole_block().find_value(
        "_pdbx_database_status.recvd_initial_deposition_date")
    return (v or "").strip('"\'')


def construct_sequence(doc) -> str:
    """The polymer sequence the entry actually crystallised, one-letter, gapless.

    Boltz is given THIS rather than one canonical Pin1 sequence for all fifteen,
    because the set is not one construct: nine of the ten held-out entries are
    the isolated PPIase domain (113 residues) while four of the five in-training
    entries are full-length with the WW domain. Handing every entry the PPIase
    domain would hand the contamination control a different, harder problem than
    the primary test -- and T2 only means something if it is measured the same
    way as T1.

    Taken from `_entity_poly`, the deposited construct, not from the modelled CA
    trace, which is missing every disordered residue.
    """
    blk = doc.sole_block()
    tab = blk.find("_entity_poly.", ["type", "pdbx_seq_one_letter_code_can"])
    for row in tab:
        if "polypeptide" in row.str(0):
            return "".join(row.str(1).split())
    return ""


def covalent_partners(doc) -> set[str]:
    """comp_ids that struct_conn records as covalently bonded to a residue."""
    out = set()
    blk = doc.sole_block()
    tab = blk.find("_struct_conn.", ["conn_type_id", "ptnr1_label_comp_id",
                                     "ptnr2_label_comp_id"])
    for row in tab:
        if row[0].strip('"\'') != "covale":
            continue
        for c in (row[1], row[2]):
            c = c.strip('"\'')
            if c not in ADDITIVES and len(c) > 2:
                out.add(c)
    return out


def smiles_for(comp_id: str, cache: Path) -> str | None:
    """Isomeric SMILES from the PDB chemical-component dictionary.

    Entry files carry only atom records for their ligands, not chemistry, so the
    SMILES has to come from the component definition. Preference order is
    OpenEye SMILES_CANONICAL (isomeric, and what the RCSB front page shows),
    then any SMILES — never a guess reconstructed from coordinates.
    """
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{comp_id}.cif"
    if not p.is_file() or p.stat().st_size < 200:
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/ligands/download/{comp_id}.cif", p)
        except Exception as e:                       # noqa: BLE001
            log.warning("%s: component fetch failed (%s)", comp_id, e)
            return None
    blk = gemmi.cif.read(str(p)).sole_block()
    tab = blk.find("_pdbx_chem_comp_descriptor.",
                   ["comp_id", "type", "program", "descriptor"])
    best = None
    for row in tab:
        typ = row.str(1).upper()
        d = row.str(3)
        if typ == "SMILES_CANONICAL" and "OpenEye" in row.str(2):
            return d
        if typ.startswith("SMILES"):
            best = best or d
    return best


def pick_ligand(st, doc) -> tuple[str | None, str]:
    """(comp_id, why) — covalent partner first, else largest non-additive."""
    counts: dict[str, int] = {}
    for ch in st[0]:
        for res in ch:
            if res.name in ADDITIVES or _is_aa(res.name) or _is_water(res.name):
                continue
            counts[res.name] = max(counts.get(res.name, 0),
                                   sum(1 for a in res if a.element.name != "H"))
    if not counts:
        return None, "no non-additive ligand"
    cov = covalent_partners(doc) & set(counts)
    if len(cov) == 1:
        c = cov.pop()
        return c, f"covalent to protein ({counts[c]} heavy atoms)"
    pool = cov or set(counts)
    c = max(pool, key=lambda k: counts[k])
    return c, ("largest of %d covalent partners" % len(cov) if cov
               else f"largest non-additive ({counts[c]} heavy atoms), no covale record")


def write_truth(st, comp_id: str, out: Path) -> tuple[int, int]:
    """ONE protein chain + the single ligand copy bound to it. → (n_copies_seen, n_res)

    Boltz is asked for a monomer, so the truth must be a monomer too. Most of
    these entries have two protein chains in the asymmetric unit and therefore
    two copies of the ligand; superposing a monomer prediction onto a dimer and
    comparing to whichever copy parsed first would produce an RMSD that depends
    on file order. The first protein chain is kept, and with it the ligand copy
    physically closest to it — the one actually in ITS pocket.
    """
    st = st.clone()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    st.setup_entities()

    prot, ligs = None, []
    for ch in st[0]:
        aas = [r for r in ch if _is_aa(r.name)]
        if aas and prot is None:
            prot = (ch.name, aas)
        ligs += [r for r in ch if r.name == comp_id]
    if prot is None or not ligs:
        return len(ligs), 0

    ref = [a.pos for r in prot[1] for a in r]

    def near(res) -> float:
        return min(a.pos.dist(p) for a in res for p in ref)

    best = min(ligs, key=near)
    keep = gemmi.Structure()
    keep.spacegroup_hm, keep.cell = st.spacegroup_hm, st.cell
    m = gemmi.Model("1")
    pc = gemmi.Chain(prot[0])
    for r in prot[1]:
        pc.add_residue(r)
    lc = gemmi.Chain("B" if prot[0] != "B" else "L")
    lc.add_residue(best)
    m.add_chain(pc)
    m.add_chain(lc)
    keep.add_model(m)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep.write_pdb(str(out))
    return len(ligs), len(prot[1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cache", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    outdir = Path(args.outdir)
    cache = Path(args.cache) if args.cache else outdir / "cif"
    rows = []
    for pid in HELD_OUT + IN_TRAINING:
        cif = fetch(pid, cache)
        doc = gemmi.cif.read(str(cif))
        st = gemmi.read_structure(str(cif))
        st.setup_entities()
        comp, why = pick_ligand(st, doc)
        rec = {"ident": pid, "held_out": pid in HELD_OUT,
               "deposited": deposition_date(doc),
               "comp_id": comp, "ligand_rule": why,
               "sequence": construct_sequence(doc)}
        if comp is None:
            rec["status"] = "no ligand"
            rows.append(rec)
            continue
        smi = smiles_for(comp, cache / "components")
        truth = outdir / "truth" / f"{pid}.pdb"
        rec["n_copies"], rec["n_res"] = write_truth(st, comp, truth)
        rec["smiles"] = smi
        rec["truth"] = str(truth)
        rec["status"] = "ok" if smi else "no SMILES in chem_comp"
        rows.append(rec)
        log.info("%s  %-7s %s copies, %s res  %s", pid, comp, rec["n_copies"], rec["n_res"], why)

    d = pd.DataFrame(rows)
    d["seq_len"] = d.sequence.fillna("").str.len()
    p = outdir / "cofold_set.csv"
    d.to_csv(p, index=False)
    print(d[["ident", "held_out", "deposited", "comp_id", "seq_len", "n_res",
             "status"]].to_string(index=False))
    print(f"\n{d.sequence.nunique()} distinct constructs → that many MSAs, "
          "not fifteen")
    print(f"\n{p}")
    bad = d[d.status != "ok"]
    if len(bad):
        print(f"\n{len(bad)} entries unusable — they are dropped, not guessed at:")
        print(bad[["ident", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
