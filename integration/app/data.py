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

import json
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
