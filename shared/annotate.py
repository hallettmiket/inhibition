"""
Purpose: The annotate-and-gate step every approach runs, in one place.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27
Input: an approach's generated frame (canonical_smiles + candidate_id)
Output: the same frame with descriptors, novelty and alert columns, and
        `rejected_at` stamped where a gate fired

WHY SHARED. T_1, T_2, T_3 and T_4 all need the same physicochemical axes, the
same novelty definition and the same structural-alert screen before anything is
docked or ranked. If each approach computed its own, the integration phase would
be pooling four different meanings of "SAscore" and "novelty" onto one plot —
which is precisely the failure the D^i contract exists to prevent.

STAMP, DO NOT DELETE. A candidate that fails a gate keeps its row and gets
`rejected_at` set to the gate's name. Downstream stages skip it, so gates still
throttle compute, but a later reweighting can see what was excluded and why.
Dropping rows makes a frame look cleaner than the run was, and makes "how many
did we lose to PAINS" unanswerable after the fact.

ALERT SCOPING IS OPTIONAL AND APPROACH-SPECIFIC. T_4 screens the R-group rather
than the whole molecule, because sulfopin's own covalent warhead trips BRENK on
every candidate and would reject the entire library for having the feature it
was designed around. Approaches with a fixed reactive core pass
`core_smarts`; approaches generating whole molecules (T_1, T_2) do not.
"""

from __future__ import annotations

import logging

import pandas as pd

from . import alerts, descriptors, novelty, smiles as smi

log = logging.getLogger(__name__)

GATE_ALERTS = "alert_gate"
GATE_UNPARSEABLE = "unparseable"


def annotate(df: pd.DataFrame, *, approach: str, core_smarts: str | None = None,
             smiles_col: str = "canonical_smiles",
             alert_limit: int | None = None) -> pd.DataFrame:
    """Add descriptors, novelty and alerts; stamp `rejected_at`.

    Parameters
    ----------
    df : pandas.DataFrame
        Generated candidates. Must carry ``smiles_col``.
    approach : str
        Approach tag, used for logging and candidate ids.
    core_smarts : str, optional
        When given, alerts are scored on the R-group outside this core as well
        as on the whole molecule (T_3/T_4). When None, whole-molecule only.
    alert_limit : int, optional
        Reject above this many total alerts. None leaves alerts annotated but
        ungated — the right default for an exploratory approach, where a hard
        alert cut discards chemotypes before anything has been measured.

    Returns
    -------
    pandas.DataFrame
        A copy with the annotation columns added. Row count is unchanged.
    """
    if smiles_col not in df.columns:
        raise KeyError(f"frame has no {smiles_col!r}")
    out = df.copy()
    n0 = len(out)

    if "rejected_at" not in out.columns:
        out["rejected_at"] = pd.NA
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [smi.candidate_id(s, prefix=approach)
                               for s in out[smiles_col]]

    bad = out["candidate_id"].isna()
    if bad.any():
        log.warning("[%s] %d candidate(s) have no InChIKey and cannot be keyed; "
                    "stamping %r", approach, int(bad.sum()), GATE_UNPARSEABLE)
        out.loc[bad, "rejected_at"] = GATE_UNPARSEABLE

    log.info("[%s] computing descriptors for %d candidates", approach, n0)
    out = descriptors.compute_frame(out, smiles_col=smiles_col)

    log.info("[%s] computing novelty against the external reference set", approach)
    out = novelty.novelty_frame(out, smiles_col=smiles_col)

    log.info("[%s] screening structural alerts (%s)", approach,
             "R-group scoped" if core_smarts else "whole molecule")
    out = alerts.screen_frame(out, smiles_col=smiles_col, core_smarts=core_smarts)

    if alert_limit is not None:
        col = "rgroup_alert_total" if core_smarts and "rgroup_alert_total" in out \
            else "whole_alert_total"
        fired = out[col].fillna(0) > alert_limit
        fired &= out["rejected_at"].isna()
        out.loc[fired, "rejected_at"] = GATE_ALERTS
        log.info("[%s] alert gate (>%d on %s): %d stamped",
                 approach, alert_limit, col, int(fired.sum()))
    else:
        log.info("[%s] alerts annotated but NOT gated — no alert_limit set",
                 approach)

    if len(out) != n0:
        raise RuntimeError(f"annotation changed row count {n0} -> {len(out)}")
    n_rej = int(out["rejected_at"].notna().sum())
    log.info("[%s] %d/%d stamped rejected and RETAINED", approach, n_rej, n0)
    return out
