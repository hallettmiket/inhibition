"""
Purpose: Build T_4's R-group library from frequency in a pinned ChEMBL pool.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the cached ChEMBL pool (immutable/inhibition/decoys/chembl_pool.csv)
Output: data/reference/rgroup_library_1.csv + a provenance sidecar

WHY DERIVE RATHER THAN HAND-PICK. D0004 rejected inheriting the prior run's
444-member library because it predated the reference set and so could not have
been grounded in it. Replacing it with substituents I choose by intuition would
be no better grounded — just differently arbitrary, and much harder for a
reviewer to argue with.

So the aryl/heteroaryl set is taken by FREQUENCY from the ChEMBL pool already
pinned in sources.lock.json. "These are the substituents medicinal chemists
actually use" is a claim anyone can check by re-running this against the same
hash-pinned file.

The LINKER set is small and explicit, because it is genuinely a design choice
rather than a frequency question: it is what connects the sulfolane nitrogen to
the aryl group, and the verified anchors are informative about it (Sulfopin uses
neopentyl, Reddi 4g cyclohexylmethyl — both CH2-linked).

CARTESIAN SIZE. |warheads| x |linkers| x |aryls|. With 8 enumerable warhead
classes (D0013) this grows fast, so the aryl count is capped and the resulting
library size is recorded, not assumed — gates.yaml deliberately holds
library_size: null until enumeration pins it.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared.manifest import Manifest          # noqa: E402

log = logging.getLogger("t4-rgroups")

POOL = Path("/data/lab_vm/immutable/inhibition/decoys/chembl_pool.csv")
OUT_CSV = REPO / "data" / "reference" / "rgroup_library_1.csv"
OUT_PROV = REPO / "data" / "reference" / "rgroup_library_1_provenance.md"

# Linkers: what joins the core nitrogen to the aryl group. Explicit and small —
# this is a design decision, and the anchors are informative (Sulfopin's
# neopentyl and Reddi 4g's cyclohexylmethyl are both CH2-linked).
LINKERS: dict[str, str] = {
    "direct": "",              # aryl bonded straight to N
    "1C": "C",                 # -CH2-  (as in Sulfopin, Reddi 4g)
    "2C": "CC",                # -CH2CH2-
    "3C": "CCC",               # -CH2CH2CH2-
    "alpha_Me": "C(C)",        # -CH(CH3)-  branch at the alpha position
    "gem_diMe": "C(C)(C)",     # -C(CH3)2-  quaternary, conformationally rigid
}

# An aryl fragment is only useful if it is small enough to be decoration rather
# than a second scaffold.
MAX_ARYL_HEAVY_ATOMS = 12
MIN_ARYL_HEAVY_ATOMS = 5


def extract_ring_systems(smiles_list: list[str], *, max_n: int) -> list[tuple[str, int]]:
    """Most frequent small aromatic ring systems in the pool, by Murcko scaffold.

    Returns
    -------
    list of (smiles, count)
        Ordered by descending frequency.
    """
    counts: Counter[str] = Counter()
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        try:
            scaf = MurckoScaffold.GetScaffoldForMol(m)
        except Exception:  # noqa: BLE001 - odd valences are skipped, not fatal
            continue
        if scaf is None or scaf.GetNumAtoms() == 0:
            continue
        # Break the scaffold into its individual ring systems so fused pairs and
        # linked biaryls both contribute their parts.
        for frag in Chem.GetMolFrags(scaf, asMols=True, sanitizeFrags=False):
            try:
                Chem.SanitizeMol(frag)
            except Exception:  # noqa: BLE001
                continue
            n = frag.GetNumHeavyAtoms()
            if not (MIN_ARYL_HEAVY_ATOMS <= n <= MAX_ARYL_HEAVY_ATOMS):
                continue
            if not any(a.GetIsAromatic() for a in frag.GetAtoms()):
                continue
            counts[Chem.MolToSmiles(frag)] += 1
    return counts.most_common(max_n)


def build(max_aryls: int = 40) -> pd.DataFrame:
    """Assemble the R-group library as linker x aryl."""
    if not POOL.is_file():
        raise SystemExit(f"ChEMBL pool not found: {POOL} — run shared.sources first")
    pool = pd.read_csv(POOL)
    log.info("scanning %d pooled molecules for frequent ring systems", len(pool))
    rings = extract_ring_systems(pool["smiles"].tolist(), max_n=max_aryls)
    log.info("kept %d aryl/heteroaryl ring systems", len(rings))

    rows = []
    for aryl_smiles, freq in rings:
        m = Chem.MolFromSmiles(aryl_smiles)
        if m is None:
            continue
        aryl_id = f"ar{len(rows) // max(1, len(LINKERS)):03d}"
        for lk_id, lk in LINKERS.items():
            # The fragment as attached to the core nitrogen: [*] marks the bond
            # to N, then the linker, then the aryl.
            frag = f"[*]{lk}{aryl_smiles}" if lk else f"[*]{aryl_smiles}"
            if Chem.MolFromSmiles(frag.replace("[*]", "C")) is None:
                continue
            rows.append({
                "rgroup_id": f"{lk_id}_{Chem.MolToSmiles(m)[:20]}",
                "linker_id": lk_id,
                "linker_smiles": lk,
                "aryl_smiles": aryl_smiles,
                "fragment_smiles": frag,
                "aryl_heavy_atoms": m.GetNumHeavyAtoms(),
                "aryl_mw": round(Descriptors.MolWt(m), 2),
                "chembl_pool_frequency": freq,
            })
    df = pd.DataFrame(rows).drop_duplicates("fragment_smiles").reset_index(drop=True)
    df["rgroup_id"] = [f"rg{i:04d}" for i in range(len(df))]
    return df


def write_provenance(df: pd.DataFrame, max_aryls: int) -> None:
    n_aryl = df["aryl_smiles"].nunique()
    OUT_PROV.write_text(f"""# R-group library — provenance

`rgroup_library_1.csv`, built by
`approaches/t4_combinatorial/build_rgroup_library.py`.

## How it was derived

**Aryl / heteroaryl groups: by frequency, not by taste.** The {n_aryl} ring
systems are the most common small aromatic Murcko ring systems in the ChEMBL
pool pinned in `config/sources.lock.json`
(`{POOL.name}`, {len(pd.read_csv(POOL))} molecules). Ring systems were kept when
they had {MIN_ARYL_HEAVY_ATOMS}-{MAX_ARYL_HEAVY_ATOMS} heavy atoms and at least
one aromatic atom — small enough to be decoration rather than a second scaffold.

This is deliberate. D0004 rejected inheriting the prior run's 444-member library
because it predated the reference set. Replacing it with substituents chosen by
intuition would be no better grounded, only differently arbitrary. "These are
the substituents medicinal chemists actually use" is a claim a reviewer can
check by re-running this against the same hash-pinned file.

**Linkers: explicit, and a design choice.** {len(LINKERS)} linkers connect the
sulfolane nitrogen to the aryl group. This is not a frequency question, and the
verified anchors are informative: Sulfopin uses neopentyl and Reddi 4g uses
cyclohexylmethyl — both CH2-linked, which is why `1C` is included and why direct
attachment and longer chains bracket it.

| id | SMILES | rationale |
|---|---|---|
| `direct` | (none) | aryl bonded straight to N |
| `1C` | `C` | the anchors' own linker |
| `2C` | `CC` | one atom longer |
| `3C` | `CCC` | reach into a deeper sub-pocket |
| `alpha_Me` | `C(C)` | branch at the alpha carbon |
| `gem_diMe` | `C(C)(C)` | quaternary, conformationally restricted |

## Size

{len(df)} R-groups = {n_aryl} aryls x {len(LINKERS)} linkers (minus any that
failed to construct). The enumerated library is this multiplied by the enumerable
warhead classes, so `gates.yaml` keeps `library_size: null` until the
enumeration stage pins the real number.

## Known limitations

- **Frequency is not suitability.** A substituent common in ChEMBL is common
  across all targets; nothing here is Pin1-specific. The pocket is shallow and
  solvent-exposed, and a chemist may well want groups this method will not
  surface.
- **Murcko scaffolds discard substituents on the ring**, so `phenyl` stands in
  for every substituted phenyl. That keeps the set small and orthogonal to the
  linker axis, but it means fluorinated and methylated variants are absent.
  Adding them is a CSV edit.
- The pool was filtered to MW 150-700 when fetched, so ring systems appearing
  only in larger molecules are underrepresented.
""", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-aryls", type=int, default=40)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    df = build(max_aryls=args.max_aryls)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    write_provenance(df, args.max_aryls)

    (Manifest(stage="build_rgroup_library", approach="t4",
              params={"max_aryls": args.max_aryls, "linkers": list(LINKERS),
                      "aryl_heavy_atom_range": [MIN_ARYL_HEAVY_ATOMS,
                                                MAX_ARYL_HEAVY_ATOMS]})
     .add_input("chembl_pool", POOL)
     .add_output("rgroup_library", OUT_CSV)
     .note(f"{len(df)} R-groups from {df['aryl_smiles'].nunique()} aryls x "
           f"{len(LINKERS)} linkers, derived by frequency not by hand")
     .write(Path("/data/lab_vm/append_only/inhibition/04_t4_combinatorial"),
            filename="rgroup_library_manifest.json"))

    print(f"{len(df)} R-groups -> {OUT_CSV}")
    print(f"  {df['aryl_smiles'].nunique()} distinct aryls x {len(LINKERS)} linkers")
    print("\nmost frequent aryls:")
    top = df.drop_duplicates("aryl_smiles").nlargest(8, "chembl_pool_frequency")
    for _, r in top.iterrows():
        print(f"  {r.aryl_smiles:28s} n={r.chembl_pool_frequency:5d}  "
              f"heavy={r.aryl_heavy_atoms}")


if __name__ == "__main__":
    main()
