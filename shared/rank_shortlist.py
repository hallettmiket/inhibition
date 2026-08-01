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

RANKING IS SIZE-DECORRELATED (D0043's open item, decided in #9). The raw score
is partly a molecular-size sort, so `rank` is computed on the metric's residual
against heavy-atom count in RANK space, not on the metric itself. Ligand
efficiency is NOT the fix — re-measured on 2026-08-01 it is worse than the raw
score in five of six pools — and it stays a displayed column rather than a sort
key. `rank_raw_metric`, `rank_metric_used` and `rank_size_decorrelated` are
carried so that a shortlist change is attributable rather than mysterious.
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
# `size_decorrelated_score` is a residual of a lower-is-better metric, so it
# inherits the direction: below the fit is better than expected for the size.
LOWER_IS_BETTER = {"affinity_kcal", "vina_affinity", "size_decorrelated_score"}

SIZE_COL = "HAC"

# The enrichment gate's graded vocabulary (D0012), plus the absence case.
# Listed so an unrecognised verdict is noticed rather than silently bucketed.
GATE_VERDICTS = {"STRONG", "WEAK", "UNDERPOWERED", "FAIL", "UNGATED"}

# ONLY `STRONG` VALIDATES A RANKING. `enrichment_gate` defines the grades:
# STRONG "enriches, and the actives set can support the claim"; WEAK
# "enriches, but within noise of not enriching". A ranking that might not
# enrich is not a validated ranking, so WEAK does not clear this bar.
VALIDATING_VERDICTS = {"STRONG"}

# Below this many docked rows a decorrelation fit is not worth trusting, and
# the group falls back to the raw metric with `rank_size_decorrelated=False`
# stamped rather than silently ranking on a fit from a handful of points.
MIN_DECORRELATION_N = 30

# Rows per size stratum. A median is stable well below the threshold needed to
# attempt a fit at all, and conflating the two starved the binning: reusing
# MIN_DECORRELATION_N gave three strata over a 100-row frame and left most of
# the trend inside the bins.
MIN_PER_STRATUM = 10
MAX_STRATA = 40


def size_decorrelated_score(metric: pd.Series, size: pd.Series) -> pd.Series:
    """The metric with its monotone dependence on heavy-atom count removed.

    WHY THIS EXISTS. D0043 established that our rankings are partly a size
    sort, and that ligand efficiency is not the fix -- it over-corrects into a
    smallness sort. Re-measured on 2026-08-01 over post-D0047, pH-7.4 data, LE
    is WORSE than the raw score in five of six pools (T_1 -0.938 vs -0.617;
    T_2/du_xu -0.654 vs +0.119), which settles issue #9's item 4 with a number.

    WHY NOT A STRAIGHT LINE. Regressing the raw metric on heavy-atom count
    leaves T_4 at Spearman -0.247 against size: the dependence is monotone but
    not linear, and a line cannot remove it. Refitting in RANK space is not
    enough either -- it still leaves rho ~0.26 when the scatter around the
    trend varies along it, which is the normal situation here because score
    spread widens with molecule size. Both were tried and both were rejected by
    `test_decorrelation_actually_removes_the_size_dependence`.

    WHAT THIS DOES INSTEAD. Non-parametric local centring: bin candidates by
    heavy-atom count into equal-population strata and subtract the MEDIAN score
    of each stratum. That assumes nothing about the shape of the size-score
    relationship, and the median makes it robust to the outliers a docking
    score reliably produces.

    THE UNITS SURVIVE. The result is still kcal/mol -- "this much better than
    the typical molecule of its size" -- so it stays readable rather than
    becoming a plausible float in unstated units, which is the failure shape
    `docs/how_this_project_breaks.md` catalogues.
    """
    out = pd.Series(pd.NA, index=metric.index, dtype="Float64")
    m = pd.to_numeric(metric, errors="coerce")
    s = pd.to_numeric(size, errors="coerce")
    ok = m.notna() & s.notna() & (s > 0)
    n = int(ok.sum())
    if n < MIN_DECORRELATION_N:
        return out
    # Enough strata to track the trend, enough rows per stratum for the median
    # to mean anything. `duplicates="drop"` handles arms whose heavy-atom count
    # takes few distinct values -- T_4 enumerates from a fixed core.
    n_bins = max(2, min(MAX_STRATA, n // MIN_PER_STRATUM))
    try:
        strata = pd.qcut(s[ok].rank(method="first"), n_bins, duplicates="drop")
    except ValueError:                       # not enough distinct sizes to bin
        return out
    centred = m[ok] - m[ok].groupby(strata, observed=True).transform("median")
    out.loc[centred.index] = centred.astype("Float64")
    return out


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
    # AN ALLOWLIST, NOT A DENYLIST. This tested `verdict not in (UNDERPOWERED,
    # UNGATED, FAIL)`, so any verdict the list did not anticipate validated the
    # ranking by default. When the non-covalent gate moved to WEAK -- defined
    # in `enrichment_gate` as "enriches, but WITHIN NOISE of not enriching", at
    # ROC-AUC 0.599 with a CI spanning 0.311-0.874 and EF1% 0.0 -- T_1 and T_2
    # were silently stamped `rank_validated=True`, contradicting D0041, the
    # README, and this module's own docstring.
    #
    # Disguise #4 in `docs/how_this_project_breaks.md`: a guard that cannot
    # fail closed because the permissive branch is the default. Naming what
    # DOES validate means a new verdict string is refused until someone
    # decides it counts.
    verdict = str(g.get("verdict", "UNGATED"))
    if verdict not in GATE_VERDICTS:
        log.warning("unrecognised gate verdict %r; treating the ranking as "
                    "NOT validated", verdict)
    out["rank_validated"] = verdict in VALIDATING_VERDICTS
    log.info("gate %s/%s: %s (ROC-AUC %s) -> rank_validated=%s",
             stratum, metric, g.get("verdict"), g.get("roc_auc"),
             bool(out["rank_validated"].iloc[0]) if len(out) else "n/a")
    return out


def rank(df: pd.DataFrame, *, metric: str, group_col: str | None,
         min_docked: int, identity_col: str | None = None,
         decorrelate_size: bool = True) -> pd.DataFrame:
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

    # A RE-RANK INVALIDATES EVERYTHING DERIVED FROM THE PREVIOUS RANKING.
    # `shortlist_synth` / `rank_synth` are built by
    # `scripts/reshortlist_synthesizable.py` from the ranking in force when it
    # ran. Re-ranking used to leave them untouched, so after D0049's
    # decorrelation D1_27 carried a `shortlist_synth` whose members reached
    # rank 90 under the NEW ranking and overlapped the new top-25 by 12 of 25
    # -- and the GUI prefers that column, so it displayed a synthesizable list
    # built from a ranking that no longer existed.
    #
    # Dropping is right rather than recomputing: this module knows nothing
    # about synthesizability rules, and a stale filtered list is worse than an
    # absent one because only the absent one is visible. Re-run the reshortlist
    # script after ranking.
    stale = [c for c in ("shortlist_synth", "rank_synth") if c in out.columns]
    if stale:
        log.warning("dropping %s — derived from the previous ranking; re-run "
                    "scripts/reshortlist_synthesizable.py", ", ".join(stale))
        out = out.drop(columns=stale)

    for col in ("rank", "group_n_docked", "group_percentile",
                "ligand_efficiency", "rank_is_selective"):
        out[col] = pd.NA
    out["rank_group"] = out[group_col] if group_col else "all"

    docked = out[metric].notna()
    log.info("%d of %d rows carry a docking result", int(docked.sum()), len(out))
    if not docked.any():
        raise SystemExit(f"no row carries a {metric} value; nothing to rank")

    # SIZE-DECORRELATED BY DEFAULT (D0043's open item, decided in #9).
    # Fitted once over the whole frame rather than per group: the size bias is
    # a property of the scoring function, not of an adduct class, and a global
    # fit is both more stable and comparable across groups.
    rank_on = metric
    out["size_decorrelated_score"] = pd.NA
    if decorrelate_size and SIZE_COL in out.columns:
        out["size_decorrelated_score"] = size_decorrelated_score(
            out[metric], out[SIZE_COL])
        covered = int((out["size_decorrelated_score"].notna() & docked).sum())
        if covered == int(docked.sum()):
            rank_on = "size_decorrelated_score"
        else:
            # Ranking part of a group on a residual and the rest on kcal/mol
            # would silently mix two scales in one ordering. Fall back wholesale
            # and say so.
            log.warning("size decorrelation covers %d of %d docked rows "
                        "(need all); ranking on raw %s instead",
                        covered, int(docked.sum()), metric)
    elif decorrelate_size:
        log.warning("no %s column; ranking on raw %s without size "
                    "decorrelation", SIZE_COL, metric)

    out["rank_metric_used"] = rank_on
    out["rank_size_decorrelated"] = (rank_on == "size_decorrelated_score")
    log.info("ranking on %s (size-decorrelated=%s)",
             rank_on, rank_on == "size_decorrelated_score")

    for grp, idx in out[docked].groupby("rank_group").groups.items():
        sub = out.loc[idx]
        if identity_col:
            # One entry per MOLECULE. Several rows may reach the same product by
            # different synthetic routes; they are one candidate, not several.
            best = (sub.groupby(identity_col)[rank_on].min()
                    .sort_values(ascending=True))
            n = len(best)
            rank_of = {ident: r for r, ident in enumerate(best.index, start=1)}
            ranks_series = sub[identity_col].map(rank_of)
        else:
            ordered = sub.sort_values(rank_on, ascending=True)
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

    # THE PRE-DECORRELATION RANK IS KEPT, NOT DISCARDED. Changing the sort key
    # churns every shortlist, and a reader comparing this frame to an older one
    # must be able to see WHICH change moved a molecule rather than inferring
    # it. `rank_raw_metric` is what this frame would have ranked before #9's
    # decision; it is never the sort key, only the audit trail.
    out["rank_raw_metric"] = pd.NA
    if rank_on != metric:
        for grp, idx in out[docked].groupby("rank_group").groups.items():
            sub = out.loc[idx]
            if identity_col:
                best = (sub.groupby(identity_col)[metric].min()
                        .sort_values(ascending=True))
                rank_of = {ident: r for r, ident in enumerate(best.index, start=1)}
                raw = sub[identity_col].map(rank_of)
            else:
                ordered = sub.sort_values(metric, ascending=True)
                raw = pd.Series(range(1, len(ordered) + 1), index=ordered.index)
            out.loc[sub.index, "rank_raw_metric"] = raw.reindex(sub.index)
    else:
        out.loc[docked, "rank_raw_metric"] = out.loc[docked, "rank"]

    # ORDINALS AND COUNTS ARE INTEGERS. These are built by assigning into
    # columns seeded with pd.NA, which makes them object and then float64 the
    # moment a null is present, so a rank reads as `304.0` in the GUI and in
    # any CSV export. Nullable Int64 keeps "did not dock" distinct from a rank
    # of 0 while displaying as the ordinal it is.
    for col in ("rank", "group_n_docked", "rank_raw_metric"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")

    if SIZE_COL in out.columns:
        hac = pd.to_numeric(out[SIZE_COL], errors="coerce")
        aff = pd.to_numeric(out[metric], errors="coerce")
        # Still computed and still displayed -- D0043 rejects LE as a SORT KEY,
        # not as a reported quantity, and chemists read it.
        out["ligand_efficiency"] = (-aff / hac).where(hac > 0).round(4)
    else:
        log.warning("no %s column; ligand efficiency not computed", SIZE_COL)
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
