"""
Purpose: Rank docked candidates and build a shortlist — shared by all four approaches.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: an approach's docked frame + its rank metric
Output: rank / percentile / ligand-efficiency columns, a shortlist flag, and the
        enrichment-gate verdict attached to every ranked row

WHY SHARED. Each approach ranks on its own metric — `vina_affinity` for T_1 and
T_2, `affinity_kcal` for T_3 and T_4 — and those numbers are NOT comparable with
each other. What must be identical is everything else: how ties are handled, how
a percentile is defined, what a quota means, and above all what gets attached to
a ranking before anyone reads it. Four approaches computing "top N" four ways is
how the integration phase ends up pooling four different meanings of "best".

NO RANKING HERE IS VALIDATED (D0031). Class-matched decoys put the covalent gate
at ROC-AUC 0.537 and the non-covalent gate at 0.535 — both indistinguishable
from chance. Every ranked row therefore carries `rank_validated = False` plus the
gate's own verdict and interval, and a shortlist is an ordering the pipeline
produced, NOT evidence that the molecules at the top bind. The GUI must display
the verdict beside the rank; a rank shown bare implies a confidence no gate here
supports.

GROUPING IS PER APPROACH, AND POST-REACTION FOR T_4 (D0029). T_4 ranks within
ADDUCT class, not warhead class: chloroacetamide, sulfamate_acetamide and
sulfonate_acetamide differ only in what leaves, so they give one identical
adduct, and quota'ing them separately spent 9 shortlist slots on 3 molecules.
T_3 has a single fixed warhead and T_1/T_2 have none, so they rank as one group.

STAMP, DO NOT DROP. Rows without a docking result keep null ranks rather than
sorting to the bottom — "did not dock" and "docked badly" are different facts.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

GATE_TOKEN = Path("/data/lab_vm/append_only/inhibition/00_shared_substrate"
                  "/enrichment_gate.token")

# Every rank metric in this project is lower-is-better (kcal/mol). Recorded
# explicitly so a future higher-is-better metric cannot be added silently.
LOWER_IS_BETTER = {"affinity_kcal", "vina_affinity"}


def load_gate_verdict(stratum: str, metric: str) -> dict:
    """The enrichment gate's verdict for one metric, or a loud placeholder.

    A missing or unreadable token is NOT treated as "fine". An unvalidated
    ranking presented without its verdict is the failure this function exists
    to prevent, so the absence is itself reported as the verdict.
    """
    if not GATE_TOKEN.is_file():
        log.warning("no enrichment gate token at %s — ranking will be marked "
                    "UNGATED", GATE_TOKEN)
        return {"verdict": "UNGATED", "reasons": ["no gate token found"]}
    try:
        tok = json.loads(GATE_TOKEN.read_text(encoding="utf-8"))
        m = tok["strata"][stratum]["metrics"][metric]
    except Exception as exc:  # noqa: BLE001
        log.warning("gate token has no %s/%s (%s) — marking UNGATED",
                    stratum, metric, exc)
        return {"verdict": "UNGATED",
                "reasons": [f"no gate entry for {stratum}/{metric}"]}
    return m


def attach_gate(df: pd.DataFrame, stratum: str, metric: str) -> pd.DataFrame:
    """Put the gate's verdict on every row, and mark the ranking unvalidated."""
    g = load_gate_verdict(stratum, metric)
    out = df.copy()
    out["gate_stratum"] = stratum
    out["gate_metric"] = metric
    out["gate_verdict"] = g.get("verdict", "UNGATED")
    out["gate_roc_auc"] = g.get("roc_auc")
    ci = g.get("roc_auc_ci") or [None, None]
    out["gate_roc_auc_ci_low"], out["gate_roc_auc_ci_high"] = ci[0], ci[1]
    out["gate_ef_1pct"] = g.get("ef_1pct")
    out["gate_n_chemotypes"] = g.get("n_chemotypes")
    # A ranking counts as validated only if the gate cleared it. Nothing here
    # currently does, and that must be visible in the data rather than in a
    # docstring somebody has to find.
    out["rank_validated"] = str(g.get("verdict", "UNGATED")) not in (
        "UNDERPOWERED", "UNGATED", "FAIL")
    log.info("gate %s/%s: %s (ROC-AUC %s) -> rank_validated=%s",
             stratum, metric, g.get("verdict"), g.get("roc_auc"),
             bool(out["rank_validated"].iloc[0]) if len(out) else "n/a")
    return out


def rank(df: pd.DataFrame, *, metric: str, group_col: str | None,
         min_docked: int, identity_col: str | None = None) -> pd.DataFrame:
    """Rank within each group; add percentile, group size and ligand efficiency.

    `group_col=None` ranks the whole frame as one group, which is correct for
    T_1, T_2 and T_3 — none of them varies the warhead.

    `identity_col` RANKS MOLECULES, NOT ROWS. T_4 carries one row per
    (R-group, warhead route) but several routes reach the same adduct, so the
    same molecule appears more than once inside a group. Ranking rows gave
    `acetamide_adduct` a top-3 of one molecule listed three times — D0029's
    defect surviving the fix for D0029, because merging the classes removed the
    duplication BETWEEN groups and left it WITHIN one. Rows sharing an identity
    receive the same rank, so a quota of 3 means three distinct molecules.
    """
    if metric not in LOWER_IS_BETTER:
        raise ValueError(
            f"{metric!r} is not a known rank metric; add it to LOWER_IS_BETTER "
            "with its direction rather than assuming one")
    out = df.copy()
    for col in ("rank", "group_n_docked", "group_percentile",
                "ligand_efficiency", "rank_is_selective"):
        out[col] = pd.NA
    out["rank_group"] = out[group_col] if group_col else "all"

    docked = out[metric].notna()
    log.info("%d of %d rows carry a docking result", int(docked.sum()), len(out))
    if not docked.any():
        raise SystemExit(f"no row carries a {metric} value; nothing to rank")

    for grp, idx in out[docked].groupby("rank_group").groups.items():
        sub = out.loc[idx]
        if identity_col:
            # One entry per MOLECULE. Several rows may reach the same product by
            # different synthetic routes; they are one candidate, not several.
            best = (sub.groupby(identity_col)[metric].min()
                    .sort_values(ascending=True))
            n = len(best)
            rank_of = {ident: r for r, ident in enumerate(best.index, start=1)}
            ranks_series = sub[identity_col].map(rank_of)
        else:
            ordered = sub.sort_values(metric, ascending=True)
            n = len(ordered)
            ranks_series = pd.Series(range(1, n + 1), index=ordered.index)
        out.loc[sub.index, "rank"] = ranks_series.reindex(sub.index)
        out.loc[sub.index, "group_n_docked"] = n
        out.loc[sub.index, "group_percentile"] = [
            round(100.0 * (n - r) / (n - 1), 2) if n > 1 else 100.0
            for r in ranks_series.reindex(sub.index)]
        out.loc[sub.index, "rank_is_selective"] = bool(n >= min_docked)
        if n < min_docked:
            log.warning("group %r has only %d successful docks (< %d): its rank "
                        "is not selective and is flagged", grp, n, min_docked)

    if "HAC" in out.columns:
        hac = pd.to_numeric(out["HAC"], errors="coerce")
        aff = pd.to_numeric(out[metric], errors="coerce")
        out["ligand_efficiency"] = (-aff / hac).where(hac > 0).round(4)
    else:
        log.warning("no HAC column; ligand efficiency not computed")
    return out


def shortlist(df: pd.DataFrame, *, quota: int,
              exclude_groups: set[str] | tuple[str, ...] = ()) -> pd.DataFrame:
    """Flag the top `quota` of each rank group, recording why each was taken.

    `exclude_groups` names groups that get no quota at all. It is generic on
    purpose: T_4 uses it for classes its reactivity triage flagged
    OUTSIDE_WINDOW, but deciding WHICH groups those are is the approach's
    business, not this module's. Default is to exclude nothing — D0019's
    position is that a flagged class carries forward with its flag attached
    rather than being silently dropped.
    """
    out = df.copy()
    out["shortlist"] = False
    out["shortlist_reason"] = pd.NA
    exclude = set(exclude_groups)

    ranked = out[out["rank"].notna()]
    for grp, sub in ranked.groupby("rank_group"):
        if grp in exclude:
            log.info("group %r excluded from the shortlist by the caller", grp)
            continue
        take = sub[sub["rank"] <= quota]
        out.loc[take.index, "shortlist"] = True
        notes = []
        if "reactivity_flag" in sub.columns:
            flag = str(sub["reactivity_flag"].iloc[0])
            if flag and flag not in ("IN_WINDOW", "nan", "<NA>"):
                notes.append(f"reactivity={flag}")
        if not bool(sub["rank_is_selective"].iloc[0]):
            notes.append(f"only {int(sub['group_n_docked'].iloc[0])} docked — "
                         "rank not selective")
        v = str(sub["gate_verdict"].iloc[0]) if "gate_verdict" in sub else "UNGATED"
        notes.append(f"gate={v}")
        out.loc[take.index, "shortlist_reason"] = (
            f"top-{quota} of {grp}" + (" (" + "; ".join(notes) + ")" if notes else ""))
    return out


def summarise(df: pd.DataFrame, metric: str) -> str:
    """A short human-readable block for the stage's stdout."""
    dk = df[df[metric].notna()]
    lines = [f"  {'group':24s} {'docked':>7s} {'best':>8s} {'median':>8s} {'best LE':>8s}"]
    lines.append("  " + "-" * (len(lines[0]) - 2))
    for grp, g in dk.groupby("rank_group"):
        le = pd.to_numeric(g["ligand_efficiency"], errors="coerce").max()
        lines.append(f"  {str(grp):24s} {len(g):7d} {g[metric].min():8.2f} "
                     f"{g[metric].median():8.2f} "
                     f"{(le if pd.notna(le) else float('nan')):8.3f}")
    v = str(df["gate_verdict"].iloc[0]) if "gate_verdict" in df else "UNGATED"
    auc = df["gate_roc_auc"].iloc[0] if "gate_roc_auc" in df else None
    lo = df.get("gate_roc_auc_ci_low", pd.Series([None])).iloc[0]
    hi = df.get("gate_roc_auc_ci_high", pd.Series([None])).iloc[0]
    lines.append("")
    lines.append(f"  GATE: {v}" + (f" — ROC-AUC {auc:.3f}" if auc is not None else "")
                 + (f" CI [{lo:.3f}, {hi:.3f}]" if lo is not None else ""))
    lines.append("  This ranking is an ordering the pipeline produced. It is NOT")
    lines.append("  evidence that the molecules at the top bind (D0031).")
    return "\n".join(lines)
