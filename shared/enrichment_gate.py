"""
Purpose: The docking-enrichment gate — does docking actually rank on THIS receptor?
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: actives + property-matched decoys, docked through the real protocols
Output: a graded verdict token consumed by every approach

THE CONTROL (adversary finding B2). A docking protocol can produce entirely
plausible scores while ordering known actives no better than chance. Nothing
crashes; the numbers look like numbers. So docking is not trusted to RANK
anything until it is shown to enrich on this receptor, with these protocols.

THE VERDICT IS GRADED, NOT BINARY (D0012). Validated Pin1 chemistry is scarce.
A confident PASS from a handful of actives manufactures precision; refusing to
proceed until the statistics are strong discards real signal that will never
become strong, because the scarcity is a property of the target. So this reports
evidence strength and lets a human weigh it:

    STRONG        enriches, and the actives set can support the claim
    WEAK          enriches, but within noise of not enriching
    UNDERPOWERED  too few independent chemotypes to tell — carry forward
                  WITH the uncertainty displayed. NOT a veto.
    FAIL          demonstrably anti-correlated with known actives

Only FAIL demotes dock_score to a displayed label.

LEAVE-ONE-CHEMOTYPE-OUT. Six actives that are three chemotypes must be reported
as three. Without this, a cluster of close analogs ranking well is counted as
several independent successes, which is the single easiest way to fool yourself
about enrichment.
"""

from __future__ import annotations

import json
import logging
import math
import random
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from rdkit import DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator as fpg

from . import smiles as smi
from .manifest import Manifest

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parent.parent
_gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)

TOKEN_DEFAULT = Path(
    "/data/lab_vm/append_only/inhibition/00_shared_substrate/enrichment_gate.token")


class EnrichmentGateError(RuntimeError):
    """The gate cannot be evaluated as configured."""


# ---------------------------------------------------------------------------
# Chemotype clustering — the honest denominator
# ---------------------------------------------------------------------------

def cluster_chemotypes(smiles: list[str], threshold: float = 0.4) -> list[int]:
    """Single-linkage ECFP4 clustering; returns a cluster id per molecule.

    The cluster COUNT, not the molecule count, is the effective sample size for
    an enrichment claim.
    """
    fps = [_gen.GetFingerprint(m) if (m := smi.to_mol(s)) is not None else None
           for s in smiles]
    n = len(fps)
    cid = [-1] * n
    nxt = 0
    for i in range(n):
        if cid[i] != -1 or fps[i] is None:
            continue
        cid[i] = nxt
        stack = [i]
        while stack:
            k = stack.pop()
            for j in range(n):
                if cid[j] != -1 or fps[j] is None:
                    continue
                if DataStructs.TanimotoSimilarity(fps[k], fps[j]) >= threshold:
                    cid[j] = nxt
                    stack.append(j)
        nxt += 1
    return cid


# ---------------------------------------------------------------------------
# Enrichment metrics
# ---------------------------------------------------------------------------

def roc_auc(scores: list[float], labels: list[int], *, higher_is_better: bool) -> float:
    """ROC-AUC via the Mann-Whitney statistic, ties counted as half."""
    s = scores if higher_is_better else [-x for x in scores]
    pos = [x for x, l in zip(s, labels) if l == 1]
    neg = [x for x, l in zip(s, labels) if l == 0]
    if not pos or not neg:
        raise EnrichmentGateError("ROC-AUC needs both actives and decoys")
    wins = sum((1.0 if p > q else 0.5 if p == q else 0.0) for p in pos for q in neg)
    return wins / (len(pos) * len(neg))


def enrichment_factor(scores: list[float], labels: list[int], *,
                      higher_is_better: bool, fraction: float = 0.01) -> float:
    """EF at a top fraction of the ranked list."""
    order = sorted(range(len(scores)), key=lambda i: scores[i],
                   reverse=higher_is_better)
    n_top = max(1, int(round(len(scores) * fraction)))
    hits_top = sum(labels[i] for i in order[:n_top])
    total_active_rate = sum(labels) / len(labels)
    if total_active_rate == 0:
        raise EnrichmentGateError("no actives present")
    return (hits_top / n_top) / total_active_rate


def bedroc(scores: list[float], labels: list[int], *, higher_is_better: bool,
           alpha: float = 20.0) -> float:
    """BEDROC (Truchon & Bayly 2007) — early recognition, exponentially weighted."""
    n = len(scores)
    order = sorted(range(n), key=lambda i: scores[i], reverse=higher_is_better)
    ranks = [i + 1 for i, idx in enumerate(order) if labels[idx] == 1]
    n_act = len(ranks)
    if n_act == 0 or n_act == n:
        raise EnrichmentGateError("BEDROC needs a mix of actives and decoys")
    ra = n_act / n
    s = sum(math.exp(-alpha * r / n) for r in ranks)
    rand = ra * (1 - math.exp(-alpha)) / (math.exp(alpha / n) - 1)
    factor = (ra * math.sinh(alpha / 2) /
              (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra)))
    return (s / rand) * factor + 1 / (1 - math.exp(alpha * (1 - ra)))


def bootstrap_ci(scores: list[float], labels: list[int], *, higher_is_better: bool,
                 n_boot: int = 2000, seed: int = 42) -> tuple[float, float]:
    """95% bootstrap CI on ROC-AUC, resampling actives and decoys separately.

    The interval is the point of the exercise. At small n the width is the
    finding — a 0.72 point estimate with a [0.41, 0.95] interval is not evidence
    of enrichment, and reporting only the 0.72 would imply it is.
    """
    rng = random.Random(seed)
    ai = [i for i, l in enumerate(labels) if l == 1]
    di = [i for i, l in enumerate(labels) if l == 0]
    vals: list[float] = []
    for _ in range(n_boot):
        sa = [rng.choice(ai) for _ in ai]
        sd = [rng.choice(di) for _ in di]
        try:
            vals.append(roc_auc([scores[i] for i in sa + sd],
                                [1] * len(sa) + [0] * len(sd),
                                higher_is_better=higher_is_better))
        except EnrichmentGateError:
            continue
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    return (vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1])


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

@dataclass
class GateResult:
    """One metric's enrichment result, with the power to interpret it."""

    metric: str
    stratum: str
    higher_is_better: bool
    n_actives: int
    n_decoys: int
    n_chemotypes: int
    roc_auc: float
    roc_auc_ci: tuple[float, float]
    ef_1pct: float
    bedroc: float
    per_chemotype_auc: dict[str, float] = field(default_factory=dict)
    # Set only when a per-candidate measurement error is supplied. Kept
    # SEPARATE from roc_auc_ci because the two answer different questions
    # and conflating them was a real defect (D0038).
    measurement_error: dict | None = None
    verdict: str = "UNDERPOWERED"
    reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = self.__dict__.copy()
        d["roc_auc_ci"] = list(self.roc_auc_ci)
        return d


def _verdict(res: GateResult, thresholds: dict) -> tuple[str, list[str]]:
    """Grade the evidence. Power is checked BEFORE the point estimates.

    THE POWER FLOOR GOVERNS **FAIL** TOO. It did not: the FAIL branch sat above
    the floor, so a damning verdict could be returned from evidence the same
    function would refuse to call STRONG. The MM-GBSA gate hit this with ONE
    active — ROC-AUC 0.140, CI [0.040, 0.240], graded FAIL, "demonstrably
    anti-correlated with known actives".

    With a single active that interval is not what it appears to be. ROC-AUC
    reduces to the fraction of decoys that one molecule beats, and bootstrapping
    resamples the same active every time, so the CI describes only decoy
    sampling. It looks tight precisely because the quantity that actually
    matters — variation between actives — has no way to enter it.

    A negative claim needs at least as much power as a positive one. Below the
    floor the honest verdict is UNDERPOWERED in BOTH directions.
    """
    reasons: list[str] = []
    min_ct = thresholds.get("min_independent_chemotypes_for_verdict", 6)
    min_act = thresholds.get("min_actives_for_verdict", 3)

    underpowered = []
    if res.n_actives < min_act:
        underpowered.append(
            f"{res.n_actives} active(s) < {min_act} required; with so few, "
            "ROC-AUC is the rank of individual molecules and its bootstrap CI "
            "reflects only decoy resampling")
    if res.n_chemotypes < min_ct:
        underpowered.append(
            f"{res.n_chemotypes} independent chemotypes < {min_ct} required; "
            "no verdict above UNDERPOWERED can be claimed regardless of the point estimates")
    if underpowered:
        if res.roc_auc < 0.5 and res.roc_auc_ci[1] < 0.5:
            underpowered.append(
                f"the point estimate (ROC-AUC {res.roc_auc:.3f}, CI upper bound "
                f"{res.roc_auc_ci[1]:.3f}) is BELOW chance and would grade FAIL "
                "with adequate power — reported, not claimed")
        return "UNDERPOWERED", underpowered

    if res.roc_auc < 0.5 and res.roc_auc_ci[1] < 0.5:
        return "FAIL", [
            f"ROC-AUC {res.roc_auc:.3f} with CI upper bound {res.roc_auc_ci[1]:.3f} "
            "below 0.5 — anti-correlated with known actives, not merely uninformative"]

    if res.roc_auc_ci[0] <= 0.5:
        reasons.append(
            f"ROC-AUC {res.roc_auc:.3f}, but the 95% CI [{res.roc_auc_ci[0]:.3f}, "
            f"{res.roc_auc_ci[1]:.3f}] includes 0.5 — within noise of not enriching")
        return "WEAK", reasons

    if (res.roc_auc >= thresholds.get("roc_auc_min", 0.70)
            and res.bedroc >= thresholds.get("bedroc_min", 0.30)):
        reasons.append(
            f"ROC-AUC {res.roc_auc:.3f} (CI excludes 0.5), BEDROC {res.bedroc:.3f}, "
            f"over {res.n_chemotypes} independent chemotypes")
        return "STRONG", reasons

    reasons.append(
        f"enriches (CI excludes 0.5) but below thresholds: ROC-AUC {res.roc_auc:.3f}, "
        f"BEDROC {res.bedroc:.3f}")
    return "WEAK", reasons


MEASUREMENT_DRAWS = 4000
MEASUREMENT_SEED = 20260729


def propagate_measurement_error(scored: pd.DataFrame, *, metric: str,
                                sem_col: str, higher_is_better: bool,
                                n_draws: int = MEASUREMENT_DRAWS,
                                seed: int = MEASUREMENT_SEED) -> dict:
    """How much of the ROC-AUC is explained by each candidate's OWN error bar.

    WHY THIS LIVES HERE AND IS SEEDED. This quantity was reported in a decision
    record and a manuscript draft while existing in no code at all -- computed
    ad hoc, unversioned, unseeded, untested, and never routed through the gate
    that grades every other metric. That is how a number becomes load-bearing
    without ever being checked (D0038).

    WHAT IT DOES AND DOES NOT MEASURE. It resamples each candidate's score from
    N(value, sem) and recomputes the metric, so it answers: given how precisely
    each candidate was measured, how much could the ordering move? It holds the
    SET of molecules fixed.

    It therefore CANNOT represent between-active variance, and with two actives
    there is none to estimate. It is not a substitute for the bootstrap CI,
    which resamples molecules and is the interval that speaks to
    generalisation. Report both or neither -- quoting only this one makes a
    two-molecule result look decided. `evaluate` returns them in separate
    fields for exactly that reason.
    """
    for col in ("label", metric, sem_col):
        if col not in scored.columns:
            raise EnrichmentGateError(f"scored frame missing {col!r}")
    sub = scored[scored[metric].notna() & scored[sem_col].notna()]
    labels = sub["label"].astype(int).to_numpy()
    if labels.sum() == 0 or labels.sum() == len(labels):
        raise EnrichmentGateError("need both actives and decoys")

    mu = sub[metric].astype(float).to_numpy()
    se = np.abs(sub[sem_col].astype(float).to_numpy())
    rng = np.random.default_rng(seed)
    draws = np.empty(n_draws, dtype=float)
    for i in range(n_draws):
        draws[i] = roc_auc(rng.normal(mu, se).tolist(), labels.tolist(),
                           higher_is_better=higher_is_better)
    point = roc_auc(mu.tolist(), labels.tolist(),
                    higher_is_better=higher_is_better)
    return {
        "point": round(float(point), 4),
        "mean": round(float(draws.mean()), 4),
        "ci": [round(float(np.percentile(draws, 2.5)), 4),
               round(float(np.percentile(draws, 97.5)), 4)],
        "p_above_chance": round(float((draws > 0.5).mean()), 4),
        "n_draws": int(n_draws),
        "seed": int(seed),
        "sem_column": sem_col,
        "holds_fixed": "the set of molecules; cannot represent between-active "
                       "variance, and with few actives there is none to "
                       "estimate. Not a substitute for the bootstrap CI.",
    }


def evaluate(scored: pd.DataFrame, *, metric: str, stratum: str,
             higher_is_better: bool, thresholds: dict | None = None,
             chemotype_threshold: float = 0.4, n_boot: int = 2000,
             sem_col: str | None = None) -> GateResult:
    """Score one metric on one stratum and grade it.

    Parameters
    ----------
    scored : pandas.DataFrame
        Needs ``canonical_smiles``, ``label`` (1 active / 0 decoy), and ``metric``.
    """
    for col in ("canonical_smiles", "label", metric):
        if col not in scored.columns:
            raise EnrichmentGateError(f"scored frame missing {col!r}")
    df = scored.dropna(subset=[metric]).reset_index(drop=True)
    scores = df[metric].astype(float).tolist()
    labels = df["label"].astype(int).tolist()

    actives = df[df.label == 1]
    cids = cluster_chemotypes(actives["canonical_smiles"].tolist(), chemotype_threshold)
    n_chemo = len(set(cids))

    res = GateResult(
        metric=metric, stratum=stratum, higher_is_better=higher_is_better,
        n_actives=int(sum(labels)), n_decoys=int(len(labels) - sum(labels)),
        n_chemotypes=n_chemo,
        roc_auc=roc_auc(scores, labels, higher_is_better=higher_is_better),
        roc_auc_ci=bootstrap_ci(scores, labels, higher_is_better=higher_is_better,
                                n_boot=n_boot),
        ef_1pct=enrichment_factor(scores, labels, higher_is_better=higher_is_better),
        bedroc=bedroc(scores, labels, higher_is_better=higher_is_better),
    )

    # Leave-one-chemotype-out: recompute AUC with each chemotype held out, so a
    # single dominant cluster cannot carry the result on its own.
    for c in sorted(set(cids)):
        keep_names = {n for n, k in zip(actives["canonical_smiles"], cids) if k != c}
        sub = df[(df.label == 0) | (df.canonical_smiles.isin(keep_names))]
        if sub.label.sum() == 0:
            continue
        try:
            res.per_chemotype_auc[f"without_cluster_{c}"] = roc_auc(
                sub[metric].astype(float).tolist(), sub["label"].astype(int).tolist(),
                higher_is_better=higher_is_better)
        except EnrichmentGateError:
            continue

    if sem_col and sem_col in scored.columns:
        try:
            res.measurement_error = propagate_measurement_error(
                scored, metric=metric, sem_col=sem_col,
                higher_is_better=higher_is_better)
        except EnrichmentGateError as exc:
            log.warning("[%s/%s] measurement-error propagation skipped: %s",
                        stratum, metric, exc)

    res.verdict, res.reasons = _verdict(res, thresholds or {})
    log.info("[%s/%s] %s — AUC %.3f CI[%.3f,%.3f] EF1%% %.1f BEDROC %.3f over %d chemotypes",
             stratum, metric, res.verdict, res.roc_auc, *res.roc_auc_ci,
             res.ef_1pct, res.bedroc, res.n_chemotypes)
    return res


MIN_DECOYS_PER_ACTIVE = 10


def filter_adequately_matched(actives: pd.DataFrame, decoys: pd.DataFrame, *,
                              min_decoys: int = MIN_DECOYS_PER_ACTIVE
                              ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Drop actives that could not be property-matched, and their decoys with them.

    WHY THIS IS NOT OPTIONAL. With a pooled decoy set, an active that found no
    property-matched decoys is still scored against decoys matched to *other,
    chemically different* actives. Docking then separates it on size and
    polarity alone — and that registers as ENRICHMENT. The control would be
    measuring its own imbalance and reporting it as success, which is worse
    than reporting failure.

    Dropping the decoys too matters for the same reason in reverse: decoys
    matched to a removed active are matched to nothing that remains, and would
    skew the retained pool toward properties no surviving active has.

    Parameters
    ----------
    actives, decoys : pandas.DataFrame
        ``decoys`` must carry ``matched_active``.
    min_decoys : int
        Below this, an active is excluded and reported.

    Returns
    -------
    (actives, decoys, excluded)
        ``excluded`` maps active name -> decoy count, for the token and the GUI.
    """
    counts = decoys["matched_active"].value_counts().to_dict()
    excluded = {str(n): int(counts.get(n, 0)) for n in actives["name"]
                if counts.get(n, 0) < min_decoys}
    if excluded:
        log.warning(
            "excluding %d active(s) with < %d property-matched decoys: %s. "
            "They are outside the decoy pool's chemical space (typically "
            "peptidic macrocycles), and scoring them against decoys matched to "
            "smaller actives would manufacture enrichment.",
            len(excluded), min_decoys, excluded)
    keep = ~actives["name"].isin(excluded)
    return (actives[keep].reset_index(drop=True),
            decoys[~decoys["matched_active"].isin(excluded)].reset_index(drop=True),
            excluded)


def load_thresholds() -> dict:
    """Gate thresholds from config/gates.yaml."""
    cfg = yaml.safe_load((_REPO_ROOT / "config" / "gates.yaml").read_text(encoding="utf-8"))
    return cfg["enrichment_gate"]["thresholds"]


def write_token(results: list[GateResult], token_path: Path | None = None) -> Path:
    """Write the gate token every approach reads before dock-based ranking.

    Records EVERY metric evaluated, not just the winner — the comparison is
    itself the evidence for D0011's covalent-metric choice.
    """
    p = token_path or TOKEN_DEFAULT
    p.parent.mkdir(parents=True, exist_ok=True)

    # MERGE, DO NOT REPLACE. This used to rebuild the payload from only the
    # current run's results, so `run_enrichment_gate.py covalent` deleted the
    # non_covalent verdict — the one T_1 and T_2 read before ranking. Nothing
    # errored; their ranking stage just reported UNGATED, which reads like "no
    # gate was ever run" rather than "another run erased it". A stratum absent
    # from THIS run keeps the verdict it already had.
    by_stratum: dict[str, dict] = {}
    if p.is_file():
        try:
            by_stratum = json.loads(p.read_text(encoding="utf-8")).get("strata", {})
            log.info("merging into existing token (strata present: %s)",
                     sorted(by_stratum))
        except Exception as exc:  # noqa: BLE001
            log.warning("existing token unreadable (%s); starting fresh", exc)
            by_stratum = {}
    # SUPERSEDE PER METRIC, NOT PER STRATUM (D0034). This used to pop the whole
    # stratum before re-adding, so a run that evaluated ONE metric erased every
    # other metric in that stratum. `run_mmgbsa_gate.py` writes a single
    # `mmgbsa_dG` result for `covalent`, and doing so deleted the
    # `affinity_kcal` verdict D0031 established — leaving the covalent stratum
    # claiming its recommended rank metric was MM-GBSA (ROC-AUC 0.140) when
    # every approach actually ranks on docking.
    #
    # This is the SAME defect as the stratum-level one described above, one
    # level down: that fix raised the merge granularity from token to stratum
    # and left the destructive replace in place underneath it. A run supersedes
    # exactly the metrics it computed and nothing else.
    for r in results:
        s = by_stratum.setdefault(r.stratum, {"metrics": {}})
        s["metrics"][r.metric] = r.to_dict()
    for stratum, s in by_stratum.items():
        ranked = [m for m, d in s["metrics"].items() if d["verdict"] != "FAIL"]
        best = max(ranked, key=lambda m: s["metrics"][m]["roc_auc"], default=None)
        s["recommended_rank_metric"] = best
        s["dock_score_ranks"] = bool(best)

    payload = {"generated_by": "shared.enrichment_gate", "strata": by_stratum}
    p.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    (Manifest(stage="enrichment_gate", approach="shared",
              params={"strata": list(by_stratum)})
     .add_config(_REPO_ROOT / "config" / "gates.yaml")
     .add_output("token", p)
     .note("Verdicts are graded (D0012). UNDERPOWERED is not a veto — the "
           "ranking carries forward with its uncertainty displayed.")
     .write(p.parent, filename="enrichment_gate_manifest.json"))
    log.info("wrote gate token -> %s", p)
    return p


def read_token(token_path: Path | None = None) -> dict:
    """Read the gate token. Approaches call this before ranking on dock_score."""
    p = token_path or TOKEN_DEFAULT
    if not p.is_file():
        raise EnrichmentGateError(
            f"no enrichment-gate token at {p}. Docking may not be used to RANK "
            "until the gate has run (control B2)."
        )
    return json.loads(p.read_text(encoding="utf-8"))
