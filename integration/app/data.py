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

# T₂ IS FIVE SEED NEIGHBOURHOODS, NOT ONE EXPERIMENT.
#
# The reseeding (issue #5) turned T₂ from "ATRA analogues" into "the CReM
# neighbourhood of a known binder", run once per seed into its own experiment
# directory. This dict was a single hardcoded `02_t2_atra_crem`, so du_xu and
# guo_pfizer finished, ranked, and remained invisible to the GUI — not because
# they were still computing, but because nothing here could reach them.
#
# RESOLVED FROM `shared.seeds`, NOT RESTATED HERE. `config/seeds.yaml` is the
# single source of the seed -> experiment mapping, and `shared/seeds.py` says
# in its own docstring that four copies of that lookup is four. A hardcoded
# dict here would be the fifth, and it would go stale the moment a seed is
# added or an experiment renamed -- the same pin-that-cannot-announce-itself
# failure as `warhead_classes_3.csv`.
#
# The degree-2 sample is appended separately: it is a DERIVED run of the ATRA
# seed into its own directory (`--experiment`), not a seed in seeds.yaml, and
# it answers a different question (does a second edit help?).
#: Short display labels. `seeds.yaml` names the seeds by their chemistry
#: ("Du-Xu-naphthalenecarboxamide"), which is right for a config and far too
#: long for a column header. The EXPERIMENT still comes from seeds.yaml — only
#: the presentation is owned here, so a renamed directory cannot go stale.
_T2_SHORT_LABEL = {
    "potter_astex": "Potter-Astex",
    "du_xu": "Du-Xu",
    "guo_pfizer": "Guo-Pfizer",
}


def _t2_seeds() -> dict:
    from shared import seeds as sd                  # noqa: PLC0415
    out = {}
    for key in sd.declared_for("t2"):
        try:
            rec = sd.resolve("t2", key, require_radius=False)
        except Exception:  # noqa: BLE001 - a malformed seed must not kill the GUI
            continue
        out[key] = {"experiment": rec["experiment"],
                    "label": _T2_SHORT_LABEL.get(key) or rec.get("name") or key}
    out["atra_degree2"] = {"experiment": "02_t2_atra_crem_degree2",
                           "label": "ATRA degree-2"}
    return out


T2_SEEDS = _t2_seeds()

# Each approach, its experiment directory, its own rank metric, and the stratum
# whose enrichment gate governs it. The metrics are NOT comparable with each
# other — that is the whole reason integration presents rather than ranks.
#
# `variants` names the seed dict for an approach that runs more than one
# neighbourhood; absent for approaches with a single experiment.
APPROACHES = {
    "t1": {"experiment": "01_t1_de_novo", "name": "T₁ · de novo (DiffSBDD)",
           "metric": "vina_affinity", "stratum": "non_covalent",
           "mechanism": "reversible", "seed": "none (seed-free)"},
    "t2": {"experiment": "02_t2_atra_crem", "name": "T₂ · CReM neighbourhood",
           "metric": "vina_affinity", "stratum": "non_covalent",
           # No static seed: T₂'s seed is whichever variant is active, and a
           # literal "ATRA" here is a wrong answer waiting to be displayed --
           # it was, under a Guo-Pfizer header. Fails to an obviously-not-a-seed
           # string rather than to a plausible one.
           "mechanism": "reversible", "seed": "(per variant)",
           "variants": T2_SEEDS, "default_variant": "atra"},
    "t3": {"experiment": "03_t3_reinvent", "name": "T₃ · decoration (LibInvent)",
           "metric": "affinity_kcal", "stratum": "covalent",
           "mechanism": "covalent Cys113", "seed": "sulfopin"},
    "t4": {"experiment": "04_t4_combinatorial", "name": "T₄ · combinatorial",
           "metric": "affinity_kcal", "stratum": "covalent",
           "mechanism": "covalent Cys113", "seed": "sulfopin"},
}

SHARED_AXES = ["MW", "cLogP", "TPSA", "QED", "SAscore", "novelty_external", "HAC"]


# WHICH VARIANT EACH APPROACH IS CURRENTLY SHOWING.
#
# Module state rather than a threaded argument, deliberately. Streamlit re-runs
# this script top to bottom on every interaction, so the sidebar sets this once
# per run and every panel reads the same value -- threading a `seed=` parameter
# through load_frame, shortlist, manifests, the convergence index and each
# panel would touch far more code for the same result. The cost is that the
# selection is global to a rerun, which is exactly the semantics wanted: the
# app shows ONE T₂ neighbourhood at a time and says which.
_ACTIVE_VARIANT: dict[str, str] = {}


def variants(approach: str) -> dict:
    """The approach's seed neighbourhoods, or {} if it has a single experiment."""
    return APPROACHES[approach].get("variants") or {}


def active_variant(approach: str) -> str | None:
    """The selected seed key, or None for a single-experiment approach."""
    v = variants(approach)
    if not v:
        return None
    key = _ACTIVE_VARIANT.get(approach) or APPROACHES[approach].get("default_variant")
    return key if key in v else next(iter(v))


def set_variant(approach: str, key: str) -> None:
    """Choose which seed neighbourhood the app shows for this approach."""
    if key in variants(approach):
        # The convergence index no longer needs clearing here: it is keyed on
        # `_frame_signature()`, which includes the active variant, so a seed
        # switch changes the key and a new index is built automatically. The
        # explicit clear was also incomplete -- it caught a seed switch but not
        # a new frame arriving from a finished docking run.
        _ACTIVE_VARIANT[approach] = key


def experiment_for(approach: str) -> str:
    """The append-only directory this approach is currently reading."""
    key = active_variant(approach)
    if key is None:
        return APPROACHES[approach]["experiment"]
    return variants(approach)[key]["experiment"]


def variant_label(approach: str) -> str | None:
    key = active_variant(approach)
    return variants(approach)[key]["label"] if key else None


def display_name(approach: str) -> str:
    """The approach's name with its active seed, for headers and selectors."""
    label = variant_label(approach)
    name = APPROACHES[approach]["name"]
    return f"{name} · seed {label}" if label else name


def variant_status(approach: str) -> list[dict]:
    """Every variant with whether it has a docked, ranked frame yet.

    A seed that is still docking must be visibly ABSENT rather than missing: a
    selector that silently omits it reads as "this seed does not exist", which
    is a different statement from "this seed has not finished".
    """
    metric = APPROACHES[approach]["metric"]
    out = []
    for key, v in variants(approach).items():
        try:
            df, path = dio.latest_frame(v["experiment"], approach)
        except Exception:  # noqa: BLE001
            out.append({"key": key, "label": v["label"], "ready": False,
                        "n": 0, "frame": None, "why": "no frame written yet"})
            continue
        docked = int(df[metric].notna().sum()) if metric in df.columns else 0
        ranked = "rank" in df.columns and df["rank"].notna().any()
        out.append({"key": key, "label": v["label"],
                    "ready": bool(docked and ranked), "n": docked,
                    "frame": path.name, "rows": len(df),
                    "why": "" if (docked and ranked)
                           else ("docked, not yet ranked" if docked
                                 else "still docking")})
    return out


def load_frame(approach: str) -> tuple[pd.DataFrame | None, str]:
    """The approach's latest frame for its ACTIVE variant, or None with a reason."""
    try:
        df, path = dio.latest_frame(experiment_for(approach), approach)
        return df, path.name
    except Exception as exc:  # noqa: BLE001 - a missing stage is a normal state
        return None, f"no frame: {exc}"


# THE SHORTLIST THE GUI SHOWS BY DEFAULT (PI decision, issue #1).
#
# `shortlist` is the raw metric top-N. `shortlist_synth` is the same quota
# rebuilt with structural-synthesizability failures REMOVED and the next-best
# passing candidates promoted in their place (scripts/reshortlist_synthesizable.py).
#
# Default to the synthesizable list. The instruction was explicit: a compound
# that fails these rules should not hold a top-25 slot. Showing the raw list
# with failures merely sorted to the bottom satisfies neither reading of that,
# and worse, it HIDES the promoted replacements entirely — they are not in the
# raw shortlist at all, so no amount of reordering can surface them.
SYNTH_COLUMN = "shortlist_synth"
BASE_COLUMN = "shortlist"


def shortlist_column(df: pd.DataFrame, prefer_synth: bool = True) -> str:
    """Which selection column to read, falling back when the rebuild is absent.

    An approach whose frame predates the rebuild has no `shortlist_synth`; it
    must still render, and it must not silently look as though it were filtered.
    Callers get the column name back so they can say which one they used.
    """
    if prefer_synth and SYNTH_COLUMN in df.columns:
        return SYNTH_COLUMN
    return BASE_COLUMN


def shortlist(approach: str, prefer_synth: bool = True) -> pd.DataFrame:
    """The approach's shortlisted candidates, with its own metric named.

    Carries `shortlist_column` and `display_rank` so the caller can label which
    list this is without re-deriving it.
    """
    df, _ = load_frame(approach)
    return _shape_shortlist(df, approach, prefer_synth)


def load_variant_frame(approach: str, key: str) -> tuple[pd.DataFrame | None, str]:
    """One named variant's frame, WITHOUT changing what the app is showing.

    The seed comparison needs several neighbourhoods at once, which the active
    variant deliberately cannot express. Reading them explicitly is correct;
    flipping the global selection in a loop would leave whichever seed happened
    to be last as the app's state.
    """
    v = variants(approach).get(key)
    if v is None:
        return None, f"no variant {key!r} for {approach}"
    try:
        df, path = dio.latest_frame(v["experiment"], approach)
        return df, path.name
    except Exception as exc:  # noqa: BLE001
        return None, f"no frame: {exc}"


def shortlist_variant(approach: str, key: str,
                      prefer_synth: bool = True) -> pd.DataFrame:
    """One named variant's shortlist, shaped exactly like `shortlist()`."""
    df, _ = load_variant_frame(approach, key)
    s = _shape_shortlist(df, approach, prefer_synth)
    if len(s):
        s["variant"] = key
        s["variant_label"] = variants(approach)[key]["label"]
    return s


def _shape_shortlist(df, approach: str, prefer_synth: bool) -> pd.DataFrame:
    """Shared shaping so the active-variant and explicit-variant paths agree.

    Two code paths producing "the shortlist" is how they drift, and a drift
    here would show one column set in the seed comparison and another in the
    shortlist panel for the same molecules.
    """
    if df is None:
        return pd.DataFrame()
    col = shortlist_column(df, prefer_synth)
    if col not in df.columns:
        return pd.DataFrame()
    s = df[df[col].fillna(False)].copy()
    s["approach"] = approach
    s["shortlist_column"] = col
    # The synthesizable list has its own contiguous 1..N ordering; the raw
    # `rank` is kept alongside as provenance, never overwritten.
    rank_col = "rank_synth" if (col == SYNTH_COLUMN
                                and "rank_synth" in s.columns) else "rank"
    # Cast here as well as at the source: a rank is an ordinal, and frames
    # written before `rank_shortlist` started emitting Int64 still carry
    # float64, which renders as `12.0`. Nullable so an unranked row stays <NA>
    # rather than becoming 0.
    s["display_rank"] = (pd.to_numeric(s[rank_col], errors="coerce").astype("Int64")
                         if rank_col in s.columns else pd.NA)
    s["metric_name"] = APPROACHES[approach]["metric"]
    s["metric_value"] = s.get(APPROACHES[approach]["metric"])
    return s


def shortlist_delta(approach: str) -> dict:
    """How the synthesizable list differs from the raw one, for labelling.

    Returns counts of promoted (in synth, not in raw) and dropped (in raw, not
    in synth) candidates. A panel that silently swapped one list for the other
    would be as misleading as the reordering this replaced.
    """
    df, _ = load_frame(approach)
    if df is None or SYNTH_COLUMN not in df.columns or BASE_COLUMN not in df.columns:
        return {"available": False, "promoted": 0, "dropped": 0, "kept": 0}
    raw = df[BASE_COLUMN].fillna(False)
    syn = df[SYNTH_COLUMN].fillna(False)
    return {"available": True,
            "promoted": int((syn & ~raw).sum()),
            "dropped": int((raw & ~syn).sum()),
            "kept": int((raw & syn).sum())}


def all_shortlists(prefer_synth: bool = True) -> pd.DataFrame:
    """Every approach's shortlist, pooled for the score-free panels ONLY.

    Pooling rows is not pooling scores. `metric_value` carries a different
    quantity per approach and must never be sorted across them; what IS
    comparable is the shared physicochemical axes and structural identity, which
    is exactly what the convergence and axes panels use.
    """
    frames = [shortlist(a, prefer_synth) for a in APPROACHES]
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
    d = DATA_ROOT / experiment_for(approach)
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
    # Resolved, not pinned: this named `warhead_classes_7.csv` while the
    # pipeline enumerated from `_10`, so the GUI listed outstanding warhead
    # structures from a library three versions behind the one in use.
    from shared import warhead_library as wl          # noqa: PLC0415
    for csv, col in ((wl.DEFAULT_LIBRARY, "structure_status"),):
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


def _frame_signature() -> tuple:
    """(approach, variant, newest frame name) for every approach.

    THE CACHE KEY THE INDEX ACTUALLY DEPENDS ON. Resolved with `dio.latest`
    rather than `load_frame` so it costs a directory glob, not a parquet read.
    """
    sig = []
    for a in APPROACHES:
        p = dio.latest(DATA_ROOT / experiment_for(a), f"D{a[-1]}", ".parquet")
        sig.append((a, active_variant(a), p.name if p is not None else None))
    return tuple(sig)


def _convergence_index() -> dict:
    """Per approach: every RANKED molecule by InChIKey, plus shortlist prints.

    Exact matching runs over every ranked row, because "was this ranked
    elsewhere at all" is the question asked. Similarity runs only over the
    other approaches' SHORTLISTS -- comparing against ~11k ranked molecules
    would cost seconds per candidate view for a number nobody would trust
    beyond the shortlist anyway.

    KEYED ON THE FRAMES IT WAS BUILT FROM. This was `@lru_cache(maxsize=1)` on
    a zero-argument function -- a cache keyed on NOTHING, in a process that
    outlives the data it caches. Rank an approach or finish a docking run and
    the index still holds the previous frames, with no indication: the
    Convergence panel and the dossier's cross-approach ranks would quietly
    describe a build that no longer exists. Streamlit keeps the process alive
    for days, so "restart it" is not a defence.

    That is catalogue entries #8 and #9 in `docs/how_this_project_breaks.md` --
    the class pool keyed on `class_id` without the query, the ligand prep keyed
    on `candidate_id` without the protonation. A cache keyed on less than its
    inputs is a cache that lies; this one was keyed on less than nothing.
    """
    return _convergence_index_for(_frame_signature())


@lru_cache(maxsize=4)
def _convergence_index_for(signature: tuple) -> dict:
    """The real build. `signature` is unused in the body -- it IS the key."""
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
            out.append({"approach": approach, "name": display_name(approach), "exact": True,
                        "candidate_id": cid, "rank": rank,
                        "n_ranked": entry.get("n_ranked", 0),
                        "similarity": 1.0})
            continue
        best_sim, best = 0.0, None
        for cid, rank, other in entry.get("short", []):
            sim = DataStructs.TanimotoSimilarity(fp, other)
            if sim > best_sim:
                best_sim, best = sim, (cid, rank)
        out.append({"approach": approach, "name": display_name(approach), "exact": False,
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

def pose_clusters() -> pd.DataFrame | None:
    """The newest pose-cluster table, or None if the stage has not run.

    Resolved by glob through `shared.outputs`, never by a pinned version --
    re-pinning a version literal is the defect this project has now written
    five times, and a test walks the AST for it.
    """
    try:
        from shared import outputs as sout
        p = sout.latest_path("blacksmith", "pose_clusters", "pose_clusters",
                             ".csv")
    except Exception:  # noqa: BLE001 - a stage that has not run is normal
        return None
    try:
        return pd.read_csv(p)
    except Exception as exc:  # noqa: BLE001
        log = __import__("logging").getLogger(__name__)
        log.warning("pose-cluster table %s could not be read: %s", p, exc)
        return None
