"""
Purpose: Everything the integration GUI reads, in one place.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: the append-only frames, run manifests and decision records
Output: cached, GUI-shaped frames

THE GUI READS; IT DOES NOT OWN (D0008). Every function here is a rendering of
files in the repo or the append-only tree. Nothing is computed for display that
is not already a recorded fact, and the app must be safe to close and rebuild at
any moment.

WHY THE LOADING IS SEPARATE FROM THE RENDERING. Streamlit reruns the whole
script on every widget interaction. Parquet reads and decision parsing are
cached here so that interaction cost does not scale with the size of the
append-only tree.
"""

from __future__ import annotations

#: Mtime of THIS file at the moment it was imported. Frozen at import, so
#: comparing it with the file's current mtime is the only reliable way to tell
#: that a running process is executing stale code -- Streamlit re-runs the
#: script on every interaction but never re-imports helper modules.
LOADED_MTIME = __import__("os").stat(__file__).st_mtime

import json
from functools import lru_cache
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from shared import decisions as dec              # noqa: E402
from shared import io as dio                     # noqa: E402

DATA_ROOT = Path("/data/lab_vm/append_only/inhibition")

# Each approach, its experiment directory, its own rank metric, and the stratum
# whose enrichment gate governs it. The metrics are NOT comparable with each
# other — that is the whole reason integration presents rather than ranks.
APPROACHES = {
    "t1": {"experiment": "01_t1_de_novo", "name": "T_1 · de novo (DiffSBDD)",
           "metric": "vina_affinity", "stratum": "non_covalent",
           "mechanism": "reversible", "seed": "none (seed-free)"},
    "t2": {"experiment": "02_t2_atra_crem", "name": "T_2 · ATRA analogues (CReM)",
           "metric": "vina_affinity", "stratum": "non_covalent",
           "mechanism": "reversible", "seed": "ATRA"},
    "t3": {"experiment": "03_t3_reinvent", "name": "T_3 · decoration (LibInvent)",
           "metric": "affinity_kcal", "stratum": "covalent",
           "mechanism": "covalent Cys113", "seed": "sulfopin"},
    "t4": {"experiment": "04_t4_combinatorial", "name": "T_4 · combinatorial",
           "metric": "affinity_kcal", "stratum": "covalent",
           "mechanism": "covalent Cys113", "seed": "sulfopin"},
}

SHARED_AXES = ["MW", "cLogP", "TPSA", "QED", "SAscore", "novelty_external", "HAC"]


def load_frame(approach: str) -> tuple[pd.DataFrame | None, str]:
    """The approach's latest frame, or None with a reason."""
    cfg = APPROACHES[approach]
    try:
        df, path = dio.latest_frame(cfg["experiment"], approach)
        return df, path.name
    except Exception as exc:  # noqa: BLE001 - a missing stage is a normal state
        return None, f"no frame: {exc}"


def shortlist(approach: str) -> pd.DataFrame:
    """The approach's shortlisted candidates, with its own metric named."""
    df, _ = load_frame(approach)
    if df is None or "shortlist" not in df.columns:
        return pd.DataFrame()
    s = df[df["shortlist"].fillna(False)].copy()
    s["approach"] = approach
    s["metric_name"] = APPROACHES[approach]["metric"]
    s["metric_value"] = s.get(APPROACHES[approach]["metric"])
    return s


def all_shortlists() -> pd.DataFrame:
    """Every approach's shortlist, pooled for the score-free panels ONLY.

    Pooling rows is not pooling scores. `metric_value` carries a different
    quantity per approach and must never be sorted across them; what IS
    comparable is the shared physicochemical axes and structural identity, which
    is exactly what the convergence and axes panels use.
    """
    frames = [shortlist(a) for a in APPROACHES]
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def gate_verdicts() -> dict:
    """The enrichment gate token, or an empty dict if none has been written."""
    p = DATA_ROOT / "00_shared_substrate" / "enrichment_gate.token"
    if not p.is_file():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}


def manifests(approach: str, limit: int = 40) -> list[dict]:
    """Run manifests for one approach, newest first.

    `git.dirty` is surfaced by the caller as a warning: it means the recorded
    commit does not fully describe the code that ran, and the outputs are
    provisional. That must not be discoverable only by reading JSON.
    """
    d = DATA_ROOT / APPROACHES[approach]["experiment"]
    out = []
    for p in sorted(d.glob("*manifest*.json"), reverse=True)[:limit]:
        try:
            m = json.loads(p.read_text(encoding="utf-8"))
            m["_file"] = p.name
            out.append(m)
        except Exception:  # noqa: BLE001
            continue
    return out


def protocol_fingerprints() -> dict[str, set]:
    """Per-approach covalent protocol fingerprints, for the parity check.

    The optional within-covalent re-score is only defensible if T_3 and T_4
    docked through the identical protocol. "We both ran gnina" is not that
    claim, so the GUI compares the recorded fingerprints and DISABLES the
    comparison when they differ rather than showing numbers made under
    different rules.
    """
    out: dict[str, set] = {}
    for a in ("t3", "t4"):
        df, _ = load_frame(a)
        if df is None or "protocol_fingerprint" not in df.columns:
            out[a] = set()
            continue
        out[a] = set(df["protocol_fingerprint"].dropna().unique())
    return out


def decisions_all() -> list[dict]:
    return [d.to_dict() for d in dec.load()]


def decisions_for(approach: str) -> list[dict]:
    """An approach's own records PLUS the shared ones it inherits."""
    return [d.to_dict() for d in dec.by_approach(approach)]


def decisions_affecting(fragment: str) -> list[dict]:
    return [d.to_dict() for d in dec.affecting(fragment)]


def open_questions() -> list[dict]:
    """Records still `proposed`, plus anything marked unverified in the sources."""
    items = [d for d in decisions_all() if str(d.get("status")) == "proposed"]
    extra = []
    lock = REPO / "config" / "sources.lock.json"
    if lock.is_file():
        try:
            s = json.loads(lock.read_text(encoding="utf-8"))
            for k, v in (s.get("sources") or {}).items():
                if str(v.get("status", "")).lower() == "pending":
                    extra.append({"id": k, "title": "source pending",
                                  "status": "pending", "approach": "shared"})
        except Exception:  # noqa: BLE001
            pass
    for csv, col in ((REPO / "data" / "reference" / "warhead_classes_7.csv",
                      "structure_status"),):
        if csv.is_file():
            try:
                t = pd.read_csv(csv)
                if col in t.columns:
                    for _, r in t[t[col].isin(["NEEDS_DESIGN", "UNVERIFIED"])].iterrows():
                        extra.append({"id": r["class_id"],
                                      "title": f"warhead structure {r[col]}",
                                      "status": str(r[col]), "approach": "t4"})
            except Exception:  # noqa: BLE001
                pass
    return items + extra


# --------------------------------------------------------------------------
# cross-approach convergence
#
# WHY THIS PANEL EXISTS AND WHY IT CURRENTLY REPORTS NOTHING. The choreography
# offers "structural convergence" -- a molecule surfaced independently by more
# than one approach -- as a soft cross-validation that needs no shared metric.
# Measured on this build, it does not occur: exact overlap between every pair
# of approaches is ZERO (checked on InChIKey, so a canonicalisation difference
# is not hiding it), and the closest cross-approach shortlist pair is T_3~T_4
# at Tanimoto 0.455, below any usual similarity threshold and unsurprising
# since both descend from sulfopin.
#
# That absence is a result, not a missing feature, and it is worth showing
# rather than leaving as an empty panel a reader assumes is broken. So the
# lookup reports the nearest analogue and its similarity when there is no exact
# match, which distinguishes "another approach found this molecule" from
# "another approach found nothing like it" -- two very different statements
# that a blank panel would collapse into one.

_FP_RADIUS = 2
_FP_BITS = 2048


def _rdkit():
    """Imported lazily: the data layer stays importable without RDKit."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    return Chem, DataStructs, rdFingerprintGenerator


@lru_cache(maxsize=1)
def _convergence_index() -> dict:
    """Per approach: every RANKED molecule by InChIKey, plus shortlist prints.

    Exact matching runs over every ranked row, because "was this ranked
    elsewhere at all" is the question asked. Similarity runs only over the
    other approaches' SHORTLISTS -- comparing against ~11k ranked molecules
    would cost seconds per candidate view for a number nobody would trust
    beyond the shortlist anyway.
    """
    Chem, _, rdfp = _rdkit()
    gen = rdfp.GetMorganGenerator(radius=_FP_RADIUS, fpSize=_FP_BITS)
    idx: dict[str, dict] = {}
    for approach in APPROACHES:
        df, _ = load_frame(approach)
        if df is None or "canonical_smiles" not in df.columns:
            idx[approach] = {"exact": {}, "short": [], "n_ranked": 0}
            continue
        ranked = df[df["rank"].notna()] if "rank" in df.columns else df.iloc[0:0]
        exact: dict[str, tuple] = {}
        for _, r in ranked.drop_duplicates("canonical_smiles").iterrows():
            m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
            if m is None:
                continue
            exact.setdefault(Chem.MolToInchiKey(m),
                             (r.get("candidate_id"), r.get("rank")))
        short = []
        s = shortlist(approach)
        for _, r in s.drop_duplicates("canonical_smiles").iterrows():
            m = Chem.MolFromSmiles(str(r["canonical_smiles"]))
            if m is None:
                continue
            short.append((r.get("candidate_id"), r.get("rank"),
                          gen.GetFingerprint(m)))
        idx[approach] = {"exact": exact, "short": short,
                         "n_ranked": int(len(exact))}
    return idx


def cross_approach_ranks(smiles: str, home: str) -> list[dict]:
    """Where else this molecule appears, one row per OTHER approach.

    `exact=True` means the identical molecule (InChIKey) was ranked there and
    the rank is reported. `exact=False` reports the most similar molecule in
    that approach's shortlist instead, so the panel says something truthful
    either way. Returns [] if the SMILES cannot be parsed.
    """
    if not smiles:
        return []
    Chem, DataStructs, rdfp = _rdkit()
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return []
    key = Chem.MolToInchiKey(m)
    gen = rdfp.GetMorganGenerator(radius=_FP_RADIUS, fpSize=_FP_BITS)
    fp = gen.GetFingerprint(m)

    out = []
    idx = _convergence_index()
    for approach, cfg in APPROACHES.items():
        if approach == home:
            continue
        entry = idx.get(approach, {})
        hit = entry.get("exact", {}).get(key)
        if hit is not None:
            cid, rank = hit
            out.append({"approach": approach, "name": cfg["name"], "exact": True,
                        "candidate_id": cid, "rank": rank,
                        "n_ranked": entry.get("n_ranked", 0),
                        "similarity": 1.0})
            continue
        best_sim, best = 0.0, None
        for cid, rank, other in entry.get("short", []):
            sim = DataStructs.TanimotoSimilarity(fp, other)
            if sim > best_sim:
                best_sim, best = sim, (cid, rank)
        out.append({"approach": approach, "name": cfg["name"], "exact": False,
                    "candidate_id": best[0] if best else None,
                    "rank": best[1] if best else None,
                    "n_ranked": entry.get("n_ranked", 0),
                    "similarity": round(best_sim, 3)})
    return out


# --------------------------------------------------------------------------
# approach parameters
#
# VALUES ARE READ FROM THE CONFIG, NEVER TRANSCRIBED. Only the SELECTION of
# which keys are headline-worthy is curated here; every value comes from
# config/approaches/*.yaml at read time. Transcribing them into the interface
# would create a second source of truth that silently drifts from the one the
# runs actually used -- the failure mode D0033/D0034 came from twice over.
#
# A missing path is skipped rather than shown as blank, so a config that gains
# or loses a key degrades quietly instead of displaying an empty row.

CONFIG_DIR = REPO / "config" / "approaches"

_CONFIG_FILES = {
    "t1": "t1_de_novo.yaml",
    "t2": "t2_atra_crem.yaml",
    "t3": "t3_reinvent.yaml",
    "t4": "t4_combinatorial.yaml",
}

# (dotted path, human label). Ordered as a reader would want to read them:
# what it starts from, what it generates, how much, then what constrains it.
_HEADLINE = {
    "t1": [("approach.seed", "seed"),
           ("generation.engine", "generator"),
           ("generation.checkpoint", "checkpoint"),
           ("generation.pocket_mode", "pocket conditioning"),
           ("generation.n_samples", "molecules sampled"),
           ("generation.batch_size", "batch size"),
           ("filters.min_heavy_atoms", "min heavy atoms"),
           ("filters.max_heavy_atoms", "max heavy atoms")],
    "t2": [("approach.seed", "seed"),
           ("expansion.engine", "expansion engine"),
           ("expansion.fragment_db", "fragment database"),
           ("expansion.radius", "CReM context radius"),
           ("expansion.max_degree", "max edit degree"),
           ("expansion.frontier_cap", "frontier cap"),
           ("expansion.max_replacements", "max replacements per site"),
           ("labelling.dock_subset_max", "docked subset cap")],
    "t3": [("approach.seed", "seed"),
           ("generation.engine", "generator"),
           ("generation.mode", "mode"),
           ("generation.prior", "prior"),
           ("generation.scaffold_smiles", "fixed scaffold"),
           ("generation.n_smiles", "SMILES requested"),
           ("warhead.class_id", "fixed warhead class"),
           ("docking.ligand_form", "docked ligand form")],
    "t4": [("approach.seed", "seed"),
           ("enumeration.core_smiles", "fixed core"),
           ("warheads.library", "warhead library"),
           ("enumeration.rgroup_library", "R-group library"),
           ("enumeration.library_size", "enumerated library size"),
           ("ranking.per_class_quota", "shortlist quota per class"),
           ("ranking.min_docked_for_meaningful_rank", "min docked for a rank")],
}


def _dig(d: dict, path: str):
    cur = d
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


@lru_cache(maxsize=8)
def approach_config(approach: str) -> dict:
    """The approach's raw config, or {} if unreadable."""
    p = CONFIG_DIR / _CONFIG_FILES.get(approach, "")
    if not p.is_file():
        return {}
    try:
        import yaml
        return yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    except Exception:  # noqa: BLE001 - a broken config must not break the view
        return {}


def approach_parameters(approach: str) -> list[dict]:
    """Headline parameters for one approach, as [{parameter, value}].

    Long paths are shortened to their basename for display: the full value is
    in the config, and a 90-character absolute path in a four-column layout
    pushes everything else off screen.
    """
    cfg = approach_config(approach)
    if not cfg:
        return []
    rows = []
    for path, label in _HEADLINE.get(approach, []):
        val = _dig(cfg, path)
        if val is None:
            continue
        if isinstance(val, str) and "/" in val and len(val) > 40:
            val = ".../" + val.rsplit("/", 1)[-1]
        if isinstance(val, (list, tuple)):
            val = ", ".join(str(v) for v in val)
        rows.append({"parameter": label, "value": val})
    return rows
