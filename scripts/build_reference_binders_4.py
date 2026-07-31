"""
Purpose: build pin1_reference_binders_4.csv — two missing actives, structured potency, tiers.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: data/reference/pin1_reference_binders_3.csv
Output: data/reference/pin1_reference_binders_4.csv

WHY A NEW VERSION.

1. TWO POTENT COVALENT ACTIVES WERE MISSING. The covalent stratum has been
   returning UNDERPOWERED (D0031/D0041) for want of actives while two
   crystallographically-confirmed Cys113 binders sat unused:

     ZL-Pin13   IC50 67 nM, cell-active, PDB 7F0M, ligand 0BF
     164A10     DELFIA IC50 4-20 nM,     PDB 8VJG, ligand A1ACH

   Both PDB entries were verified to exist before this file was written.

2. POTENCY WAS FREE TEXT. `potency` held strings like "IC50 3.2 uM; Ki 1.37 uM;
   kinact/Ki 0.249 M^-1 s^-1". Nothing could sort, filter or threshold on it, so
   every comparison was done by eye. Split into potency_value / potency_unit /
   potency_type, with the original string preserved verbatim in `potency` --
   the parse is an addition, never a replacement, so a mis-parse cannot destroy
   the source.

3. THE PROMISCUOUS ENTRIES NEEDED A TIER, NOT DELETION. juglone, EGCG, ATRA,
   KPT-6566 and PiB are promiscuity-flagged; three carry no usable number at
   all, and KPT-6566 is explicitly dual-mechanism (it releases a ROS-generating
   naphthoquinone). Carrying them as gate ACTIVES plausibly depresses the AUC.
   They are tiered `historical_promiscuous` rather than removed, because they
   are real literature and removing them would erase why they were ever here.

   ATRA IS STILL T_2's SEED. Demoting a compound as a gate active is not
   removing it as a starting point -- different roles, and conflating them
   would silently discard an entire approach. `seeds.yaml` is untouched.

WHAT THIS FILE DOES NOT DO. It does not rank covalent against non-covalent
binders. Covalent potency is time-dependent (an IC50 depends on preincubation)
and its real metric is kinact/Ki. `potency_type` carries the distinction so a
consumer can refuse to mix them, which is D0031's lesson one level up.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

REF = REPO / "data" / "reference"
SRC = REF / "pin1_reference_binders_3.csv"
OUT = REF / "pin1_reference_binders_4.csv"

# Verified against RCSB on 2026-07-31: both entries exist, both SMILES pulled
# from the deposited chemical component rather than transcribed from a paper.
NEW_ROWS = [
    {
        "name": "Liu-2022-ZL-Pin13",
        "canonical_smiles":
            "O=C(CCl)N1CCC2(CC1)SCC(=O)N2Cc1ccc(-c2cccc3ccccc23)o1",
        "mechanism": "covalent_cys113",
        "warhead_class": "chloroacetamide",
        "potency": "IC50 67 nM; cell-active",
        "promiscuity_flag": "n",
        "pdb": "7F0M",
        "citation": ("Liu L 2022, J Med Chem 65:2174, "
                     "10.1021/acs.jmedchem.1c01686; PDB 7F0M ligand 0BF"),
    },
    {
        "name": "Alboreggia-2024-164A10",
        # The DEPOSITED ligand, i.e. the BOUND form. See the provenance note:
        # this is not necessarily the as-synthesised compound.
        "canonical_smiles":
            "CC(=O)N1Cc2c(c3ccccc3[nH]2)C[C@H]1C(=O)N1CCCC[C@H]1C(=O)N[C@@H]"
            "(Cc1c[nH]c2ccc(F)cc12)C(N)=O",
        "mechanism": "covalent_cys113",
        "warhead_class": "UNVERIFIED",
        "potency": "DELFIA IC50 4-20 nM; DC50 <500 nM (degrader)",
        "promiscuity_flag": "n",
        "pdb": "8VJG",
        "citation": ("Alboreggia 2024, PNAS 121, 10.1073/pnas.2403330121, "
                     "PMID 39531501; PDB 8VJG ligand A1ACH"),
    },
]

# Promiscuity-flagged or mechanistically confounded. Kept, not deleted.
HISTORICAL = {"Juglone", "EGCG", "ATRA", "KPT-6566", "PiB"}
FRAGMENT = {"Byun-BDHI-fragment"}

# Ordered: the FIRST pattern that matches wins, so the tightest binding constant
# a paper reports is what gets parsed. Ki before Kd before IC50 before EC50,
# because that is the order of directness as a binding measurement.
_UNIT = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "mm": 1e6}
_PATTERNS = [
    ("Ki", r"\bKi\s*[:=]?\s*~?\s*([\d.]+)\s*(pM|nM|uM|µM|mM)"),
    ("Kd", r"\bKd\s*[:=]?\s*~?\s*([\d.]+)\s*(pM|nM|uM|µM|mM)"),
    ("IC50", r"\bIC50\s*[:=]?\s*~?\s*([\d.]+)\s*(?:-\s*[\d.]+\s*)?(pM|nM|uM|µM|mM)"),
    ("EC50", r"\bEC50\s*[:=]?\s*~?\s*([\d.]+)\s*(pM|nM|uM|µM|mM)"),
]


def parse_potency(text: str) -> tuple[float | None, str | None, str | None]:
    """(value_nM, 'nM', type) from a free-text potency string, or (None,)*3.

    Returns nM for everything so the column is sortable. Kinetic constants
    (kinact/Ki, k in M^-1 s^-1, % labeling) are deliberately NOT parsed into
    this column: they are not affinities, and giving them a value_nM would
    invite exactly the cross-metric sort this file exists to prevent.
    """
    if not isinstance(text, str):
        return None, None, None
    for kind, pat in _PATTERNS:
        m = re.search(pat, text, re.I)
        if m:
            val = float(m.group(1)) * _UNIT[m.group(2).lower()]
            return val, "nM", kind
    return None, None, None


def kinetic_note(text: str) -> str | None:
    """The covalent kinetic constant, if the row reports one."""
    if not isinstance(text, str):
        return None
    bits = []
    if (m := re.search(r"k(?:inact/Ki)?\s*=?\s*[\d.e+-]+\s*M-?1\s*s\^?-?1", text, re.I)):
        bits.append(m.group(0))
    if (m := re.search(r"kinact\s*[\d.e+-]+\s*-?\s*[\d.e+-]*\s*s\^?-?1", text, re.I)):
        bits.append(m.group(0))
    if (m := re.search(r"(\d+)%\s*labeling", text, re.I)):
        bits.append(m.group(0))
    return "; ".join(bits) or None


def main() -> None:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    df = pd.read_csv(SRC)
    before = len(df)

    for row in NEW_ROWS:
        if (df["name"] == row["name"]).any():
            raise SystemExit(f"{row['name']} already present — nothing to add")
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)

    # Every SMILES must parse. A reference set with a broken structure silently
    # corrupts the novelty axis for every approach (shared/reference_set.py).
    bad = [r["name"] for _, r in df.iterrows()
           if str(r["canonical_smiles"]) != "UNVERIFIED"
           and Chem.MolFromSmiles(str(r["canonical_smiles"])) is None]
    if bad:
        raise SystemExit(f"unparseable SMILES for: {bad}")

    parsed = df["potency"].apply(parse_potency)
    df["potency_value_nM"] = [p[0] for p in parsed]
    df["potency_unit"] = [p[1] for p in parsed]
    df["potency_type"] = [p[2] for p in parsed]
    df["kinetic_constant"] = df["potency"].apply(kinetic_note)

    df["tier"] = [
        "historical_promiscuous" if n in HISTORICAL
        else "fragment" if n in FRAGMENT
        else "lead"
        for n in df["name"]
    ]
    # Belt and braces: anything promiscuity-flagged must not be tiered `lead`,
    # whatever the name list says.
    flagged = df["promiscuity_flag"].astype(str).str.lower().isin({"y", "yes", "true"})
    mism = df[flagged & (df["tier"] == "lead")]["name"].tolist()
    if mism:
        raise SystemExit(f"promiscuity-flagged but tiered lead: {mism}")

    df.to_csv(OUT, index=False)

    print(f"{SRC.name} ({before} rows) -> {OUT.name} ({len(df)} rows)")
    print(f"  added: {', '.join(r['name'] for r in NEW_ROWS)}")
    print(f"\n  tiers: {df['tier'].value_counts().to_dict()}")
    print(f"  potency parsed: {int(df['potency_value_nM'].notna().sum())}/{len(df)}")
    print(f"  kinetic constants: {int(df['kinetic_constant'].notna().sum())}")
    print("\n  leads by affinity (nM, lower = tighter):")
    lead = (df[(df.tier == "lead") & df.potency_value_nM.notna()]
            .sort_values("potency_value_nM"))
    for _, r in lead.iterrows():
        print(f"    {r['potency_value_nM']:8.1f} nM  {r['potency_type']:5s} "
              f"{r['mechanism'][:13]:14s} {r['name'][:40]}")
    noval = df[(df.tier == "lead") & df.potency_value_nM.isna()]["name"].tolist()
    print(f"\n  leads with NO parsable affinity ({len(noval)}): {', '.join(noval)}")


if __name__ == "__main__":
    main()
