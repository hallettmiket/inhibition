"""
Purpose: Build property-matched decoys for the docking-enrichment gate.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: the frozen actives set; a cached ChEMBL molecule pool
Output: decoys_<stratum>.csv under immutable/inhibition/decoys/, plus a manifest

WHAT A DECOY HAS TO BE, AND WHY IT IS FIDDLY.

A decoy must be **physically similar** to an active (so the docking cannot
succeed by trivially preferring heavier or greasier molecules) and
**topologically dissimilar** (so it is unlikely to actually bind). Get the first
wrong and the enrichment is measuring molecular weight. Get the second wrong and
you have salted the decoy set with real binders, which makes docking look worse
than it is.

DUD-E style, built from ChEMBL because that source is already a dependency and
can be hash-pinned — a decoy set nobody can regenerate is not a control
(D0007's lesson applied here).

Anything with a recorded Pin1 (CHEMBL2288) activity is excluded from the pool
outright. Note the earlier trap: CHEMBL3391 is Threonine--tRNA ligase, NOT Pin1,
and two compounds were misattributed on exactly that confusion. The target id is
therefore a named constant here rather than a literal typed at each call site.
"""

from __future__ import annotations

import json
import logging
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator as fpg

from . import descriptors as desc
from . import smiles as smi
from .manifest import Manifest

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
PIN1_TARGET = "CHEMBL2288"          # Pin1. NOT CHEMBL3391 (Threonine--tRNA ligase).

_REPO_ROOT = Path(__file__).resolve().parent.parent
POOL_PATH = Path("/data/lab_vm/immutable/inhibition/decoys/chembl_pool.csv")
OUT_DIR = Path("/data/lab_vm/immutable/inhibition/decoys")

_gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)

# Properties a decoy must match, and how close is close enough. Windows are
# deliberately generous: a decoy set that matches too tightly becomes a set of
# near-analogs, which is the opposite of what is wanted.
MATCH_TOLERANCE = {
    "MW": 40.0,
    "cLogP": 1.0,
    "TPSA": 30.0,
    "rot_bonds": 3,
    "formal_charge": 0,      # exact
}
# Above this ECFP4 Tanimoto to ANY active, a candidate is refused as a decoy —
# it is too likely to be a real binder.
MAX_SIMILARITY_TO_ACTIVE = 0.35


class DecoyError(RuntimeError):
    """Decoy generation failed or produced an unusable set."""


@dataclass
class DecoySet:
    """Decoys plus the bookkeeping needed to defend them."""

    frame: pd.DataFrame
    per_active: dict[str, int]
    pool_size: int
    shortfalls: dict[str, int]


def _get(url: str, retries: int = 3) -> dict:
    """GET JSON with retries — the ChEMBL API rate-limits under bulk paging."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return json.loads(fh.read().decode("utf-8"))
        except Exception as exc:  # noqa: BLE001 - transient network is expected
            if attempt == retries - 1:
                raise DecoyError(f"ChEMBL request failed: {url}\n{exc}") from exc
            time.sleep(2 * (attempt + 1))
    raise DecoyError("unreachable")


def pin1_actives_chembl_ids() -> set[str]:
    """Every ChEMBL molecule id with a recorded Pin1 activity.

    These are excluded from the decoy pool. A decoy that is a known active is
    not a decoy.
    """
    ids: set[str] = set()
    offset, limit = 0, 1000
    while True:
        url = (f"{CHEMBL_API}/activity.json?target_chembl_id={PIN1_TARGET}"
               f"&limit={limit}&offset={offset}")
        data = _get(url)
        acts = data.get("activities", [])
        ids.update(a["molecule_chembl_id"] for a in acts if a.get("molecule_chembl_id"))
        if len(acts) < limit:
            break
        offset += limit
    log.info("excluded %d molecules with recorded Pin1 (%s) activity", len(ids), PIN1_TARGET)
    return ids


def fetch_pool(target_n: int = 20000, *, force: bool = False) -> pd.DataFrame:
    """Fetch (and cache) a drug-like ChEMBL molecule pool to draw decoys from.

    Cached under immutable/ so the decoy set is regenerable rather than a
    one-off artifact.
    """
    if POOL_PATH.is_file() and not force:
        df = pd.read_csv(POOL_PATH)
        log.info("using cached ChEMBL pool: %d molecules", len(df))
        return df

    excluded = pin1_actives_chembl_ids()
    rows, offset, limit = [], 0, 1000
    POOL_PATH.parent.mkdir(parents=True, exist_ok=True)
    while len(rows) < target_n:
        q = urllib.parse.urlencode({
            "molecule_properties__mw_freebase__gte": 150,
            "molecule_properties__mw_freebase__lte": 700,
            "limit": limit, "offset": offset,
        })
        data = _get(f"{CHEMBL_API}/molecule.json?{q}")
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
        if offset % 10000 == 0:
            log.info("pool: %d molecules fetched", len(rows))
    df = pd.DataFrame(rows).drop_duplicates("chembl_id")
    df.to_csv(POOL_PATH, index=False)
    log.info("cached ChEMBL pool: %d molecules -> %s", len(df), POOL_PATH)
    return df


def _fp(s: str):
    m = smi.to_mol(s)
    return _gen.GetFingerprint(m) if m is not None else None


def warhead_motifs() -> dict[str, Chem.Mol]:
    """Reactive-atom SMARTS per warhead class, from the warhead library."""
    from . import warhead_library as wl
    lib = wl.load()
    out = {}
    for _, r in lib.iterrows():
        p = Chem.MolFromSmarts(str(r["reactive_atom_smarts"]))
        if p is not None:
            out[str(r["class_id"])] = p
    return out


def classify_warheads(smiles: str, motifs: dict[str, Chem.Mol]) -> list[str]:
    """Which warhead classes' reactive-atom patterns a molecule matches."""
    m = smi.to_mol(smiles)
    if m is None:
        return []
    return [c for c, p in motifs.items() if m.HasSubstructMatch(p)]


def build(actives: pd.DataFrame, *, stratum: str, n_per_active: int = 50,
          pool: pd.DataFrame | None = None,
          require_warhead: bool = False,
          active_warhead_class: dict[str, str] | None = None) -> DecoySet:
    """Property-match decoys to each active.

    Parameters
    ----------
    actives : pandas.DataFrame
        Must carry ``name`` and ``canonical_smiles``.
    stratum : str
        ``covalent`` or ``non_covalent`` — recorded, and used in the filename.
    n_per_active : int
        Target decoys per active. A shortfall is REPORTED, never silently
        accepted: uneven decoy counts skew enrichment per active.
    require_warhead : bool
        Covalent stratum only. A decoy with no electrophile CANNOT be covalently
        docked — gnina needs a reactive atom to bond. Without this, ~90% of the
        decoy set is unrunnable and the surviving comparison is "electrophiles
        versus inert molecules", which docking wins trivially and which answers
        nothing.
    active_warhead_class : dict, optional
        ``{active_name: warhead_class}``. When given, decoys are preferentially
        drawn from the SAME warhead class as the active they match. That makes
        the control ask the interesting question — does docking prefer *our*
        chloroacetamide over *other* chloroacetamides — rather than the trivial
        one. Falls back to any electrophile when a class is too sparse, and the
        fallback is recorded.

    Returns
    -------
    DecoySet
    """
    pool = pool if pool is not None else fetch_pool()
    log.info("matching %d actives against a pool of %d", len(actives), len(pool))

    motifs = warhead_motifs() if require_warhead else {}
    if require_warhead:
        pool = pool.copy()
        pool["warhead_classes"] = [classify_warheads(s, motifs) for s in pool["smiles"]]
        before = len(pool)
        pool = pool[pool["warhead_classes"].map(len) > 0]
        log.info("warhead constraint: %d/%d pool molecules carry an electrophile",
                 len(pool), before)
        if pool.empty:
            raise DecoyError(
                "no electrophile-bearing molecules in the pool — a covalent "
                "decoy set cannot be built from it.")

    # Descriptors for the pool, once.
    pool = pool.copy()
    pool_desc = [desc.compute(s) for s in pool["smiles"]]
    for col in ("MW", "cLogP", "TPSA", "rot_bonds", "formal_charge"):
        pool[col] = [d.get(col) for d in pool_desc]
    pool = pool.dropna(subset=["MW", "cLogP", "TPSA", "rot_bonds", "formal_charge"])

    active_fps = [(r["name"], _fp(r["canonical_smiles"])) for _, r in actives.iterrows()]
    active_fps = [(n, f) for n, f in active_fps if f is not None]

    chosen: list[dict] = []
    per_active: dict[str, int] = {}
    shortfalls: dict[str, int] = {}
    used: set[str] = set()

    for _, act in actives.iterrows():
        a_desc = desc.compute(act["canonical_smiles"])
        if a_desc["MW"] is None:
            log.warning("active %r unparseable; skipped", act["name"])
            continue
        m = pool
        for prop, tol in MATCH_TOLERANCE.items():
            if tol == 0:
                m = m[m[prop] == a_desc[prop]]
            else:
                m = m[(m[prop] - a_desc[prop]).abs() <= tol]
        m = m[~m["chembl_id"].isin(used)]

        # Prefer decoys carrying the SAME warhead class as this active, so the
        # control asks whether docking prefers our electrophile over other
        # electrophiles of the same chemistry. Fall back to any electrophile
        # when the class is too sparse — recorded, not silent.
        want_class = (active_warhead_class or {}).get(act["name"])
        class_matched = 0
        if require_warhead and want_class:
            same = m[m["warhead_classes"].map(lambda cs: want_class in cs)]
            if len(same) >= n_per_active:
                m = same
                class_matched = len(same)
            elif not same.empty:
                m = pd.concat([same, m[~m.index.isin(same.index)]])
                class_matched = len(same)

        picked = 0
        for _, cand in m.iterrows():
            if picked >= n_per_active:
                break
            f = _fp(cand["smiles"])
            if f is None:
                continue
            # Topologically dissimilar to EVERY active, not just this one —
            # otherwise a decoy for one active can be a near-analog of another.
            if max((DataStructs.TanimotoSimilarity(f, af) for _, af in active_fps),
                   default=1.0) > MAX_SIMILARITY_TO_ACTIVE:
                continue
            rec = {"chembl_id": cand["chembl_id"], "smiles": cand["smiles"],
                   "matched_active": act["name"], "stratum": stratum, "label": 0}
            if require_warhead:
                rec["warhead_classes"] = "|".join(cand["warhead_classes"])
                rec["class_matched"] = bool(want_class and want_class in cand["warhead_classes"])
            chosen.append(rec)
            used.add(cand["chembl_id"])
            picked += 1
        if require_warhead and want_class and class_matched < n_per_active:
            log.info("active %r: only %d pool molecules share its warhead class %r; "
                     "topped up with other electrophiles",
                     act["name"], class_matched, want_class)

        per_active[act["name"]] = picked
        if picked < n_per_active:
            shortfalls[act["name"]] = n_per_active - picked
            log.warning("active %r: only %d/%d decoys matched", act["name"],
                        picked, n_per_active)

    frame = pd.DataFrame(chosen)
    return DecoySet(frame=frame, per_active=per_active, pool_size=len(pool),
                    shortfalls=shortfalls)


def write(ds: DecoySet, stratum: str) -> Path:
    """Write a decoy set and its manifest under immutable/."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"decoys_{stratum}_1.csv"
    if out.exists():
        raise DecoyError(f"{out} exists; decoys are immutable — bump the version")
    ds.frame.to_csv(out, index=False)

    mf = (Manifest(stage="build_decoys", approach="shared",
                   params={"stratum": stratum, "pool_size": ds.pool_size,
                           "tolerance": MATCH_TOLERANCE,
                           "max_similarity_to_active": MAX_SIMILARITY_TO_ACTIVE,
                           "per_active": ds.per_active})
          .add_input("chembl_pool", POOL_PATH)
          .add_output("decoys", out))
    if ds.shortfalls:
        mf.note(f"decoy shortfall for {len(ds.shortfalls)} active(s): {ds.shortfalls}. "
                "Uneven decoy counts skew per-active enrichment; interpret accordingly.")
    mf.write(Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"),
             filename=f"decoys_{stratum}_manifest.json")
    log.info("wrote %d decoys -> %s", len(ds.frame), out)
    return out
