"""
Purpose: Survey every Pin1 (UniProt Q13526) PDB entry and build the redocking
         benchmark's candidate table -- entry metadata + non-polymer ligands,
         with crystallisation additives filtered out.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-31
Input: RCSB search API (exact_match on the UniProt accession) + RCSB GraphQL
Output: outputs/blacksmith/redock_pin1/pdb_ligand_survey_1.csv
        outputs/blacksmith/redock_pin1/pdb_entry_survey_1.csv

WHY THE FILTER IS BY NAME *AND* MASS, NOT FREQUENCY. The naive "most common
HET code" list for this target is dominated by PEG fragments, glycerol and
sulfate -- cryoprotectants and buffer, present in most crystals and binding
nothing. Ranking ligands by frequency therefore returns the additives, not the
chemistry. Membership here requires passing a mass window AND not matching an
additive name/code pattern, and every rejection is written to the CSV with its
reason so the drop list is auditable rather than implied.
"""

from __future__ import annotations

import json
import logging
import re
import time
import urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

log = logging.getLogger(__name__)

OUT_DIR = Path(__file__).resolve().parent.parent / "outputs" / "blacksmith" / "redock_pin1"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
UNIPROT = "Q13526"

# Mass window. The lower bound removes ions, single glycols and DMSO; the upper
# bound is generous (peptidomimetics run large).
#
# 150 Da IS NOT THE FLOOR, BECAUSE THE FRAGMENTS ARE REAL. A 150 Da cut removes
# 22 genuine Xiao-series fragments (1H-indol-4-ylmethanol, quinolin-5-amine,
# ...) alongside the buffer. Those are deposited ligands with real density and
# dropping them silently would be exactly the filtered-subset dishonesty this
# benchmark is meant to avoid. They are KEPT and tiered instead: `drug_like`
# (>= 150 Da) carries the headline number, `fragment` (< 150 Da) is reported
# separately, because recovering a 7-heavy-atom fragment within 2 A is a much
# weaker claim than recovering a 30-atom inhibitor and pooling them would
# flatter the result.
MIN_MW = 100.0
FRAGMENT_MW = 150.0
MAX_MW = 900.0

# Explicit additive codes. This is receptor_prep.STRIPPABLE_HET plus the wider
# set of PEG/cryo/buffer codes that appear across 190 entries rather than in
# 6VAJ alone.
ADDITIVE_CODES = {
    "HOH", "WAT", "DOD",
    "GOL", "EDO", "MPD", "PGO", "PDO", "BU3", "MRD", "TRE", "SUC", "GLC",
    "PEG", "PG4", "PGE", "1PE", "P6G", "2PE", "XPE", "7PE", "12P", "15P",
    "PE4", "PE8", "P33", "P4G", "PG5", "PG6", "M2M", "TOE", "SPD", "SPM",
    "SO4", "PO4", "NO3", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "CD", "NI",
    "CU", "FE", "IOD", "BR", "F", "ACY", "CIT", "FLC", "TLA", "MLA", "MES",
    "EPE", "TRS", "IMD", "BTB", "CAC", "HEP", "MOP", "BIS", "NHE",
    "DMS", "ACT", "FMT", "OXL", "GLY", "BME", "DTT", "DTU", "MRC",
    "IPA", "MOH", "EOH", "ACN", "ACE", "NH4", "AZI", "SCN", "PER",
    "LDA", "C8E", "BOG", "LMT", "DDQ", "OCT", "HEX", "MYR", "PLM", "OLA",
    "UNX", "UNL", "UNK",
    # Longer PEGs that clear the 100 Da floor on mass alone. PE3 (634 Da,
    # "tridecaoxahentetracontane-1,41-diol") passed both the mass window and
    # the name regex on the first pass and was about to be redocked as a
    # ligand -- the exact failure mode the brief warned about.
    "PE3", "PE5", "PE6", "PE7", "PEU", "PGF", "P2K", "P3G", "P5G", "MPO",
    # Free amino acids deposited as non-polymer components. Pin1's substrate is
    # a pSer/pThr-Pro motif, so a lone ALA or PRO is a peptide remnant or a
    # cryo additive, not a small-molecule binder.
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # N-/O-glycans. These are covalent post-translational modifications of the
    # protein, not small molecules that bind a pocket, so redocking them is not
    # the experiment. 2M9F/2M9J carry NAG.
    "NAG", "NDG", "BMA", "MAN", "BGC", "GAL", "GLA", "FUC", "FUL", "XYS",
    "SIA", "NGA", "A2G", "RAM",
    # Heavy-atom derivatives / phasing agents.
    "TAS", "HG", "AU", "PT", "PB", "IR", "OS", "SM", "EU", "GD", "YB", "XE",
    "KR", "LMR",
}

# Name patterns. A new PEG-type code we have never seen must still be caught,
# so the name is checked as well as the code.
#
# NO BARE HETEROCYCLE NAMES HERE. The first version carried `imidazole\s*$` to
# catch the buffer, and it rejected SGH, SFW, TKK and P5N -- four genuine
# phenyl-imidazole Pin1 ligands whose chemical names simply END in "imidazole".
# A pattern that describes a substructure will always match real medicinal
# chemistry; buffers must be matched as WHOLE names (ADDITIVE_EXACT_NAMES) or
# by code, never by a trailing fragment.
ADDITIVE_NAME_RE = re.compile(
    r"(polyethylene\s*glycol|\bpeg\b|ethylene\s*glycol|glycerol|propanediol|"
    r"butanediol|pentanediol|cryoprotectant|sulfate\s*ion|phosphate\s*ion|"
    r"acetate\s*ion|formate\s*ion|citrate|tartrate|malonate|"
    r"\btris\b|hepes|\bmes\b|bicine|tricine|"
    r"dimethyl\s*sulfoxide|mercaptoethanol|dithiothreitol|dithioerythritol|"
    r"\bion\b|detergent|monoolein|maltoside|spermidine|spermine|"
    r"trehalose|sucrose|"
    # Polyether chains: "...TRIDECAOXAHENTETRACONTANE-1,41-DIOL" and friends.
    r"(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca|trideca|"
    r"tetradeca|pentadeca)oxa|"
    r"(di|tri|tetra|penta|hexa|hepta|octa|nona)ethylene\s*glycol)",
    re.IGNORECASE,
)

# Buffers whose WHOLE name is the additive. Matched exactly (case-folded) so a
# substituted derivative bearing the same ring is not caught.
ADDITIVE_EXACT_NAMES = {
    "imidazole", "glycerol", "acetic acid", "formic acid", "urea",
    "1,2-ethanediol", "ethanol", "methanol", "2-propanol", "acetone",
    "dimethyl sulfoxide", "beta-mercaptoethanol", "acetonitrile",
}

# A polyether the name regex misses is still catchable in the structure:
# three or more consecutive -O-CH2-CH2-O- units is a PEG, whatever it is called.
POLYETHER_SMARTS = "OCCOCCOCCO"

_ENTRY_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined }
    struct { title }
    nonpolymer_entities {
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
      nonpolymer_comp {
        chem_comp { id name formula formula_weight type }
        pdbx_chem_comp_descriptor { type program descriptor }
      }
    }
  }
}
"""


def _post(url: str, payload: dict, retries: int = 4) -> dict:
    """POST JSON with retries -- RCSB rate-limits large GraphQL batches."""
    body = json.dumps(payload).encode()
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, data=body, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as fh:
                return json.loads(fh.read().decode())
        except Exception as exc:  # noqa: BLE001 - retry transient RCSB failures
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"POST {url} failed after {retries} tries: {last}")


def fetch_entry_ids() -> list[str]:
    """Every PDB entry mapped to the Pin1 UniProt accession."""
    payload = {
        "query": {
            "type": "terminal", "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers."
                             "reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": UNIPROT,
            },
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    d = _post(SEARCH_URL, payload)
    ids = [r["identifier"] for r in d["result_set"]]
    log.info("RCSB reports %d entries for %s", d["total_count"], UNIPROT)
    return ids


def _smiles_from_descriptors(descriptors: list[dict] | None) -> str | None:
    """Prefer a stereo-bearing canonical SMILES; fall back to any SMILES."""
    if not descriptors:
        return None
    ranked = sorted(
        (d for d in descriptors if d.get("type", "").upper().startswith("SMILES")),
        key=lambda d: (0 if "CANONICAL" in d.get("type", "").upper() else 1,
                       0 if d.get("program", "") == "OpenEye OEToolkits" else 1),
    )
    return ranked[0]["descriptor"] if ranked else None


_PEG_PATT = Chem.MolFromSmarts(POLYETHER_SMARTS)


def classify(comp_id: str, name: str, mw: float | None,
             smiles: str | None) -> tuple[str, int | None]:
    """('ligand' | rejection reason, heavy_atom_count).

    Checked in order code -> exact name -> name pattern -> structure -> mass,
    so the most specific evidence wins and the reason recorded is the real one.
    """
    heavy: int | None = None
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is not None:
        heavy = mol.GetNumHeavyAtoms()

    if comp_id in ADDITIVE_CODES:
        return "additive_code", heavy
    if name and name.strip().lower() in ADDITIVE_EXACT_NAMES:
        return "additive_name_exact", heavy
    if name and ADDITIVE_NAME_RE.search(name):
        return "additive_name", heavy
    if mol is not None and _PEG_PATT is not None and mol.HasSubstructMatch(_PEG_PATT):
        return "additive_polyether_structure", heavy
    if mw is None:
        return "no_formula_weight", heavy
    if mw < MIN_MW:
        return f"mw_below_{MIN_MW:g}", heavy
    if mw > MAX_MW:
        return f"mw_above_{MAX_MW:g}", heavy
    return "ligand", heavy


def survey() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the per-entry and per-(entry, ligand) tables."""
    ids = fetch_entry_ids()
    entries: list[dict] = []
    rows: list[dict] = []
    for i in range(0, len(ids), 40):
        batch = ids[i:i + 40]
        d = _post(GRAPHQL_URL, {"query": _ENTRY_QUERY, "variables": {"ids": batch}})
        if "errors" in d:
            raise RuntimeError(f"GraphQL errors: {d['errors'][:2]}")
        for e in d["data"]["entries"]:
            pdb_id = e["rcsb_id"]
            methods = [m["method"] for m in (e.get("exptl") or [])]
            res = (e.get("rcsb_entry_info") or {}).get("resolution_combined")
            resolution = float(res[0]) if res else None
            entries.append({
                "pdb_id": pdb_id,
                "method": "; ".join(methods),
                "resolution_a": resolution,
                "title": ((e.get("struct") or {}).get("title") or "")[:200],
                "n_nonpolymer_entities": len(e.get("nonpolymer_entities") or []),
            })
            for ent in (e.get("nonpolymer_entities") or []):
                comp = (ent.get("nonpolymer_comp") or {}).get("chem_comp") or {}
                desc = (ent.get("nonpolymer_comp") or {}).get("pdbx_chem_comp_descriptor")
                mw = comp.get("formula_weight")
                mw = float(mw) if mw is not None else None
                name = comp.get("name") or ""
                comp_id = comp.get("id") or ""
                smi = _smiles_from_descriptors(desc)
                verdict, heavy = classify(comp_id, name, mw, smi)
                rows.append({
                    "pdb_id": pdb_id,
                    "method": "; ".join(methods),
                    "resolution_a": resolution,
                    "comp_id": comp_id,
                    "comp_name": name[:160],
                    "formula": comp.get("formula"),
                    "formula_weight": mw,
                    "comp_type": comp.get("type"),
                    "smiles": smi,
                    "heavy_atoms": heavy,
                    "auth_asym_ids": ",".join(
                        (ent.get("rcsb_nonpolymer_entity_container_identifiers")
                         or {}).get("auth_asym_ids") or []),
                    "classification": verdict,
                    "tier": ("drug_like" if (mw or 0) >= FRAGMENT_MW else "fragment")
                            if verdict == "ligand" else "",
                })
        log.info("fetched %d/%d entries", min(i + 40, len(ids)), len(ids))
    return pd.DataFrame(entries), pd.DataFrame(rows)


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    entry_df, lig_df = survey()
    entry_df.to_csv(OUT_DIR / "pdb_entry_survey_1.csv", index=False)
    lig_df.to_csv(OUT_DIR / "pdb_ligand_survey_1.csv", index=False)

    log.info("entries: %d (%s)", len(entry_df),
             entry_df["method"].value_counts().to_dict())
    log.info("non-polymer components: %d rows", len(lig_df))
    log.info("classification: %s", lig_df["classification"].value_counts().to_dict())
    log.info("tier: %s", lig_df[lig_df.classification == "ligand"]["tier"].value_counts().to_dict())
    keep = lig_df[lig_df["classification"] == "ligand"]
    log.info("distinct candidate ligands: %d over %d entries",
             keep["comp_id"].nunique(), keep["pdb_id"].nunique())
    xray = keep[(keep["method"].str.contains("X-RAY")) & (keep["resolution_a"] <= 2.0)]
    log.info("X-ray <=2.0 A: %d distinct ligands over %d entries",
             xray["comp_id"].nunique(), xray["pdb_id"].nunique())


if __name__ == "__main__":
    main()
