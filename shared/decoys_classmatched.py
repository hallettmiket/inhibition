"""
Purpose: Covalent decoys drawn from the SAME warhead class as each active.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the frozen actives set; per-class ChEMBL substructure searches
Output: decoys_covalent_4.csv under immutable/inhibition/decoys/, plus provenance

WHY A SECOND GENERATOR RATHER THAN A PATCH (D0031).

`shared.decoys.build` property-matches against ONE generic 20k ChEMBL pool and
then filters for a warhead. That order is the defect. A generic sample of
drug-like molecules contains almost no sulfamate acetamides or nitro-chloro-
azines, so the rare classes are wiped out by the property filter before the
warhead filter ever runs, and the shortfall is topped up from whatever class
does survive. The result was 104 acrylamide / 5 naphthoquinone / 3
chloroacetamide decoys against actives that are chloroacetamides, sulfamates,
naphthoquinones and an sNAr azine.

D0020 says affinity is not comparable across warhead classes. A gate scoring one
class's actives against another class's decoys therefore measures the quantity
D0020 calls meaningless. Inverting the order — retrieve BY CLASS first, then
property-match within the class — is what makes the control test binding rather
than chemotype.

RETRIEVAL IS LOOSE, VERIFICATION IS STRICT, AND THEY ARE DIFFERENT PATTERNS.
The ChEMBL query casts wide (`CC(=O)CCl`). Membership is then decided by the
whole warhead group plus a successful adduct transform. Deciding class
membership from the retrieval pattern is exactly the error that let
cyclophosphamide and lomustine into the chloroacetamide bucket (D0028), and
naming a chemotype from a reactive-atom SMARTS is the recurring defect behind
D0025, D0028 and D0029.

A DECOY MUST BE DOCKABLE AS AN ADDUCT. Every candidate is put through
`covalent_adduct.to_adduct_form` with its class's SMARTS. One that fails cannot
be docked by the same protocol the approaches use, so it cannot sit in a control
for them (D0022, D0030).
"""

from __future__ import annotations

import json
import logging
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator as fpg

from . import covalent_adduct as cad
from . import descriptors as desc
from . import smiles as smi
from . import warhead_library as wl
from .decoys import (CHEMBL_API, MATCH_TOLERANCE, MAX_SIMILARITY_TO_ACTIVE,
                     OUT_DIR, pin1_actives_chembl_ids)

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)

# class_pools_4, and under append_only rather than immutable (read-only by
# project rule). A NEW directory on purpose: class_pools_3 caches the pools
# retrieved with the OLD chemotype queries, and a cache keyed only on class_id
# cannot tell that the query behind it has changed. Reusing the directory would
# have silently served the 3-molecule snar pool built from the narrow query
# after the query was relaxed -- which is exactly what happened once already.
CLASS_POOL_DIR = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"
                      "/decoys/class_pools_4")
MAX_PER_CLASS_FETCH = 3000

_REPO = Path(__file__).resolve().parent.parent
CHEMOTYPES = _REPO / "data" / "reference" / "decoy_chemotypes_3.csv"


def load_chemotypes() -> pd.DataFrame:
    """The decoy chemotype table: group SMARTS, retrieval query, and the
    library class used for the adduct check.

    SEPARATE FROM THE WARHEAD LIBRARY ON PURPOSE. The library's classes are
    T_4's ENUMERATION units — `naphthoquinone_c2` and `naphthoquinone_benzo`
    are two attachment points on one chemistry. A decoy has to match the
    active's CHEMISTRY, and Juglone (5-hydroxy-1,4-naphthoquinone) is a genuine
    naphthoquinone that sits at neither attachment point. Deciding decoy class
    membership from an enumeration unit rejected Juglone from its own class.
    """
    return pd.read_csv(CHEMOTYPES)


@dataclass
class ClassPool:
    """Everything retrieved for one warhead class, and what survived."""

    class_id: str
    query: str
    n_retrieved: int = 0
    n_group_match: int = 0
    n_adduct_ok: int = 0
    molecules: pd.DataFrame = field(default_factory=pd.DataFrame)
    unavailable_reason: str | None = None


def _get(url: str, retries: int = 3) -> dict:
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"ChEMBL request failed after {retries} tries: {last}")


def warhead_group_pattern(row) -> Chem.Mol | None:
    """The chemotype's whole reactive group, as SMARTS.

    MUST be SMARTS, not SMILES. An earlier version built this from the warhead
    library's `[*]C(=O)CCl` fragments with MolFromSmiles: in SMARTS `[*]` means
    "any atom", but in SMILES it is a dummy of atomic number 0 that matches only
    another dummy. That pattern matched 0 of 2,038 retrieved chloroacetamides
    and 0 of 3,963 naphthoquinones, and every chemotype was reported chemically
    unavailable while the pools held thousands of molecules.
    """
    return Chem.MolFromSmarts(str(row["group_smarts"]))


def verify_class(smiles: str, class_id: str, pattern: Chem.Mol, library) -> bool:
    """Membership: the chemotype's whole reactive group AND a valid adduct form.

    Two independent conditions. The group decides chemistry; the adduct decides
    dockability under the same protocol the approaches use (D0022, D0030). A
    molecule that passes only the first cannot sit in a control for them.
    """
    m = smi.to_mol(smiles)
    if m is None or pattern is None:
        return False
    if not m.HasSubstructMatch(pattern):
        return False
    try:
        cad.to_adduct_form(smiles, class_id, library=library)
    except Exception:  # noqa: BLE001 - AdductError and RDKit failures alike
        return False
    return True


def fetch_class_pool(class_id: str, query: str, *, force: bool = False,
                     max_n: int = MAX_PER_CLASS_FETCH) -> pd.DataFrame:
    """ChEMBL substructure search for one class, cached under immutable/."""
    CLASS_POOL_DIR.mkdir(parents=True, exist_ok=True)
    cache = CLASS_POOL_DIR / f"{class_id}.csv"
    if cache.is_file() and not force:
        try:
            df = pd.read_csv(cache)
        except pd.errors.EmptyDataError:
            # A chemotype that genuinely returns nothing caches an EMPTY file.
            # That is a real, reusable result — "ChEMBL has none of these" — so
            # it must read back as an empty pool, not crash the whole build.
            df = pd.DataFrame(columns=["chembl_id", "smiles"])
        log.info("[%s] cached class pool: %d molecules", class_id, len(df))
        return df

    excluded = pin1_actives_chembl_ids()
    rows, offset, limit = [], 0, 1000
    total_seen = 0
    while len(rows) < max_n:
        url = (f"{CHEMBL_API}/substructure/{urllib.parse.quote(query)}.json"
               f"?limit={limit}&offset={offset}")
        data = _get(url)
        total_seen = max(total_seen,
                         int(data.get("page_meta", {}).get("total_count", 0) or 0))
        mols = data.get("molecules", [])
        if not mols:
            break
        for m in mols:
            cid = m.get("molecule_chembl_id")
            s = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not s or cid in excluded:
                continue
            rows.append({"chembl_id": cid, "smiles": s})
        offset += limit
        total = data.get("page_meta", {}).get("total_count", 0)
        if offset >= total:
            break
        log.info("[%s] %d retrieved", class_id, len(rows))
    df = pd.DataFrame(rows).drop_duplicates("chembl_id")

    # NEVER CACHE "NOTHING EXISTS" WHEN CHEMBL SAYS OTHERWISE.
    #
    # A transient empty `molecules` list on a 200 response would otherwise be
    # written to the cache and read back forever as the finding "this chemotype
    # does not exist in the database". That happened to chloroacetamide — the
    # single most important class here, Sulfopin's own — and the run reported it
    # as chemically unavailable while ChEMBL held 4,430 of them. A retrieval
    # failure and a real absence must never look the same.
    reported = int(total_seen or 0)
    if df.empty and reported > 0:
        raise RuntimeError(
            f"[{class_id}] ChEMBL reports {reported:,} hits for {query!r} but "
            "none were retrieved. This is a fetch failure, not an absence; "
            "refusing to cache it as one.")
    df.to_csv(cache, index=False)
    log.info("[%s] cached %d molecules -> %s", class_id, len(df), cache)
    return df


def build_class_pool(chemotype: str, library, chemotypes=None, *,
                     force: bool = False) -> ClassPool:
    """Retrieve, verify and describe every usable decoy candidate for a chemotype."""
    ct = chemotypes if chemotypes is not None else load_chemotypes()
    rows = ct[ct["chemotype"] == chemotype]
    if rows.empty:
        return ClassPool(chemotype, "", unavailable_reason="unknown chemotype")
    row = rows.iloc[0]
    class_id = str(row["representative_class"])
    query = row.get("chembl_query")
    if not isinstance(query, str) or not query:
        return ClassPool(chemotype, "", unavailable_reason="no retrieval query")

    pool = ClassPool(chemotype, query)
    raw = fetch_class_pool(chemotype, query, force=force)
    pool.n_retrieved = len(raw)
    if raw.empty:
        pool.unavailable_reason = (
            f"ChEMBL substructure search for {query!r} returns nothing; this "
            "chemotype does not exist in the database")
        return pool

    patt = warhead_group_pattern(row)
    keep = []
    for _, r in raw.iterrows():
        s = smi.canonical(r["smiles"])
        if s is None:
            continue
        m = smi.to_mol(s)
        if m is None or patt is None or not m.HasSubstructMatch(patt):
            continue
        pool.n_group_match += 1
        try:
            cad.to_adduct_form(s, class_id, library=library)
        except Exception:  # noqa: BLE001
            continue
        pool.n_adduct_ok += 1
        keep.append({"chembl_id": r["chembl_id"], "canonical_smiles": s})

    pool.molecules = pd.DataFrame(keep)
    if pool.molecules.empty:
        pool.unavailable_reason = (
            f"{pool.n_retrieved} retrieved, {pool.n_group_match} carry the whole "
            f"warhead group, 0 produce a valid adduct")
    else:
        pool.molecules = desc.compute_frame(pool.molecules,
                                            smiles_col="canonical_smiles")
    return pool


def _fp(s: str):
    m = smi.to_mol(s)
    return _gen.GetFingerprint(m) if m is not None else None


def match_within_class(active: pd.Series, pool: ClassPool, active_fps: list,
                       *, n_per_active: int, used: set[str]) -> pd.DataFrame:
    """Property-match decoys to one active, drawing ONLY from its own class."""
    if pool.molecules.empty:
        return pd.DataFrame()
    a = desc.compute_frame(
        pd.DataFrame([{"canonical_smiles": active["canonical_smiles"]}]),
        smiles_col="canonical_smiles").iloc[0]

    m = pool.molecules[~pool.molecules["chembl_id"].isin(used)].copy()
    for prop, tol in MATCH_TOLERANCE.items():
        if prop not in m.columns:
            continue
        m = m[m[prop] == a[prop]] if tol == 0 else \
            m[(m[prop] - a[prop]).abs() <= tol]
    if m.empty:
        return m

    # Topological dissimilarity: a decoy too close to ANY active is more likely
    # a real binder than a control, and salting the set that way makes docking
    # look worse than it is.
    sims = []
    for s in m["canonical_smiles"]:
        fp = _fp(s)
        sims.append(max(DataStructs.BulkTanimotoSimilarity(fp, active_fps))
                    if fp is not None else 1.0)
    m = m.assign(max_sim_to_active=sims)
    m = m[m["max_sim_to_active"] <= MAX_SIMILARITY_TO_ACTIVE]

    # Closest on properties first, so the set is a control rather than a
    # convenience sample of whatever the search returned.
    m = m.assign(_d=((m["MW"] - a["MW"]).abs() / MATCH_TOLERANCE["MW"]
                     + (m["cLogP"] - a["cLogP"]).abs() / MATCH_TOLERANCE["cLogP"]
                     + (m["TPSA"] - a["TPSA"]).abs() / MATCH_TOLERANCE["TPSA"]))
    return m.nsmallest(n_per_active, "_d").drop(columns="_d")
