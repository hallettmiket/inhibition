#!/usr/bin/env python3
"""
Purpose: stamp InChIKey and the docked species onto a candidate frame.
Author: Timothy Wu (with Claude Code)
Date: 2026-08-11
Input: --frame <D3|D4> (defaults to the newest of each)
Output: a new integer-versioned parquet with identity columns added

WHY (#58). Two identity defects, one cause: the pipeline carries SMILES and
nothing else, and a SMILES is not an identifier.

INCHIKEY. **1,782 of 1,782** T4 candidates carry a stereocentre and the frame had
no InChIKey column. The 21b memo named this as the single most likely way to be
wrong without anyone noticing -- "a SMILES that loses the wedge is
indistinguishable from the wrong enantiomer" -- and it has already bitten: the
first covalent workup refused because the adduct's warhead stereocentre was
undefined. The InChIKey is the only identifier here that survives a lost wedge.

DOCKED SPECIES. @tt8804, 2026-08-11: dock the **pH 7.4 form**, and correct the
warhead library to match. The covalent workup was refusing on t4_716800c125a7
because the pose carried the protonated amine (D4: `charge_ph74 = +1`) while the
library's adduct template carried the neutral form -- two stages holding
different species of one molecule. 726 of 1,782 T4 candidates are cations, so
this decides what 41% of the library docks AS.

`docked_smiles` is now written explicitly rather than inferred at dock time, so
the species that was docked is a recorded fact rather than a reconstruction.
"""

from __future__ import annotations

import argparse
import glob
import logging
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

log = logging.getLogger("stamp-identity")
DATA = Path("/data/lab_vm/append_only/inhibition")
FRAMES = {"D3": "03_t3_reinvent", "D4": "04_t4_combinatorial"}


def latest(stem: str) -> Path:
    fs = sorted(glob.glob(str(DATA / FRAMES[stem] / f"{stem}_*.parquet")),
                key=lambda q: int(q.rsplit("_", 1)[1].split(".")[0]))
    if not fs:
        raise SystemExit(f"no {stem} frame")
    return Path(fs[-1])


def protonate(smiles: str, target_charge: float):
    """The species at pH 7.4, as SMILES, or None if it cannot be built.

    Deliberately narrow: it neutralises nothing and invents nothing. It only
    applies the charge the frame ALREADY recorded in `charge_ph74`, which was
    computed upstream, so this cannot disagree with the column it is keyed on.
    A molecule whose recorded charge is 0 is returned unchanged.
    """
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    # pandas NA is not falsy -- `if target_charge` raises on it. A missing charge
    # is NOT zero: it means nobody computed one, and docking such a molecule as
    # neutral would be the same silent substitution this whole script exists to
    # stop. It returns None and gets stamped.
    import pandas as _pd
    if target_charge is None or _pd.isna(target_charge):
        return None
    want = int(target_charge)
    if want == 0:
        return Chem.MolToSmiles(m)
    # ALREADY THE RIGHT SPECIES. 135 of the 234 first-pass failures were
    # molecules whose SMILES already carried the charge -- the frame and the
    # structure agreed and the function refused anyway. Protonating one of these
    # would have produced a +2.
    if Chem.GetFormalCharge(m) == want:
        return Chem.MolToSmiles(m)

    # In basicity order, most basic first, so the site that is actually
    # protonated at 7.4 is the one that gets the proton. Amides and anilines stay
    # excluded: they are not basic here, and protonating one invents a species
    # that does not exist.
    CATION_SITES = (
        "[NX3;$(N=C(N)N);!+]",                       # guanidine
        "[NX2;$(N=C[NX3]);!+]",                      # amidine, the imine N
        "[NX3;H0,H1,H2;!$(N[C,S]=[O,S,N]);!$(N-a);!$(N=*);!+]",   # aliphatic amine
        "[nX2;!+]",                                  # pyridine-type aromatic N
    )
    ANION_SITES = (
        "[CX3](=O)[OX2H1]",                          # carboxylic acid
        "[SX4](=O)(=O)[OX2H1]",                      # sulfonic acid
        "[PX4](=O)([OX2H1])",                        # phosphate/phosphonate
        "[nX3H1]1nnnc1", "[nX3H1]1ncnn1",            # tetrazole
    )
    for sma in (CATION_SITES if want > 0 else ANION_SITES):
        patt = Chem.MolFromSmarts(sma)
        if patt is None:
            continue
        hits = m.GetSubstructMatches(patt)
        if not hits:
            continue
        em = Chem.RWMol(m)
        if want > 0:
            a = em.GetAtomWithIdx(hits[0][0])
            a.SetFormalCharge(1)
            a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
            a.SetNoImplicit(True)
        else:
            # the acidic O/N is the last atom of every ANION pattern above
            a = em.GetAtomWithIdx(hits[0][-1])
            a.SetFormalCharge(-1)
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(True)
        try:
            out = em.GetMol(); Chem.SanitizeMol(out)
        except Exception:                                  # noqa: BLE001
            continue
        if Chem.GetFormalCharge(out) == want:
            return Chem.MolToSmiles(out)
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--frame", choices=("D3", "D4"), action="append")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    for stem in (args.frame or ["D3", "D4"]):
        f = latest(stem)
        d = pd.read_parquet(f)
        log.info("%s: %d rows", f.name, len(d))

        keys, mismatch = [], 0
        for smi in d.canonical_smiles.astype(str):
            m = Chem.MolFromSmiles(smi)
            keys.append(Chem.MolToInchiKey(m) if m is not None else None)
        d["inchikey"] = keys
        d["inchikey_ok"] = d.inchikey.notna()

        # THE SPECIES THAT WILL BE DOCKED, WRITTEN DOWN. @tt8804 chose the pH 7.4
        # form; the frame already carries the charge, so this only realises it.
        if "charge_ph74" in d.columns:
            sp, failed = [], 0
            for smi, q in zip(d.canonical_smiles.astype(str), d.charge_ph74):
                out = protonate(smi, q)
                if out is None:
                    failed += 1
                sp.append(out)
            d["docked_smiles"] = sp
            d["docked_charge"] = d.charge_ph74
            # A molecule whose pH-7.4 form could not be built is STAMPED, not
            # silently docked as the neutral: that is exactly the substitution
            # #58 is about.
            d["docked_species_ok"] = d.docked_smiles.notna()
            log.info("  pH 7.4 species: %d built, %d could not be "
                     "(stamped docked_species_ok = False)", len(d) - failed, failed)
            # A build failure is NOT a difference. Counting NaN as "differs"
            # reported 726 changed when 492 had, which reads as full coverage.
            built = d[d.docked_smiles.notna()]
            mismatch = int((built.docked_smiles != built.canonical_smiles).sum())
            log.info("  species differs from the neutral SMILES: %d of %d built",
                     mismatch, len(built))

        dupes = d[d.inchikey.notna()].inchikey.duplicated().sum()
        log.info("  InChIKeys: %d written, %d duplicate",
                 int(d.inchikey.notna().sum()), int(dupes))

        if args.dry_run:
            continue
        n = int(f.stem.rsplit("_", 1)[1]) + 1
        dest = f.parent / f"{stem}_{n}.parquet"
        d.to_parquet(dest, index=False)
        print(f"  -> {dest}")


if __name__ == "__main__":
    main()
