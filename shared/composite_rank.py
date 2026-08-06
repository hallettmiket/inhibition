"""
Purpose: combine the four ranking components into one queryable ranked list that never drops a candidate.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06
Input: a candidate frame carrying whichever of the four components have been measured for each row
Output: the same rows — all of them — plus a composite score, its uncertainty
        interval, and a per-row record of exactly which components contributed

STAGE 4 OF `docs/ranking_rationale.md`, which that document has listed as
"unbuilt" and "now the main gap" since 2026-08-05. This is the layer that turns
several partially-measured components into one ordering.

    THE CONTRACT IS RANK AND RETAIN. @tt8804: "does not matter if its top 5 or
    top 5000, we want a ranked list that we can query, we dont want to delete
    any candidates." Nothing in this module filters, drops, truncates or
    deduplicates. Every input row leaves with a row, and `rank_candidates`
    raises if that is ever not true. A floor, a flag or a verdict is expressed
    as a COLUMN, because in a rank-and-retain list a column is the only way to
    express one.

THE FOUR COMPONENTS, and the fact that governs the whole design:

    pose_quality   nac enrichment          all 5,769 candidates
    consensus      agreement of top poses  cheap; expected to cover all
    bpmd           metadynamics stability  ~50-100 molecules, ~3.5 h GPU each
    md_residence   pocket occupancy        similar

MOST CANDIDATES WILL CARRY ONE OR TWO OF FOUR. That is not an edge case to be
tidied away, it is the normal row, and it is the reason this module is not a
weighted sum.

--------------------------------------------------------------------------
WHY NOT IMPUTE, AND WHY NOT DROP
--------------------------------------------------------------------------

The three usual answers to a missing component are all wrong here, and each is
wrong in a way this project has already been bitten by:

  * **Impute zero** — the missing value becomes the worst measured value. This
    is precisely "not measured masquerading as measured badly", and @tt8804
    ruled it out by name: a candidate with no BPMD score must not rank below
    one with a bad BPMD score.
  * **Impute the cohort mean** — quieter and worse. It manufactures a
    measurement, writes it into a column indistinguishable from a real one, and
    is then read downstream as if someone had run the simulation. That is
    disguise #4 of `docs/how_this_project_breaks.md`: a value that is populated
    and plausible and computed from nothing.
  * **Drop the row** — forbidden by the contract above, and it would silently
    make the ranked list mean "the molecules we could afford to measure".

WHAT THIS MODULE DOES INSTEAD: a missing component contributes its full possible
RANGE rather than a point. Every component is mapped onto [0, 1] (higher =
better), so total ignorance about a component is the interval [0, 1] — which is
a true statement, needs no imputation, and has exactly the ordering behaviour
that was asked for:

    bad, measured   ->  contributes ~[0.0, 0.05]   drags the composite DOWN
    not measured    ->  contributes  [0.0, 1.0]    midpoint 0.5, width 1.0
    good, measured  ->  contributes ~[0.95, 1.0]   pushes the composite UP

so `bad < missing < good`, with no value invented anywhere. The composite is the
interval's MIDPOINT and the interval travels beside it in `composite_lo` /
`composite_hi` / `composite_width`. The 0.5 a missing component contributes is
not an imputed average; it is the midpoint of a bounded interval of total
ignorance, and it arrives with a width of 1.0 that announces itself. A reader
who takes the midpoint without the width has been given the means not to.

    `composite_width` IS THE CONFIDENCE READOUT. A candidate measured on all
    four has a narrow interval; a candidate with only pose quality has a wide
    one. That is the answer to @tt8804's "how confident are the poses" asked at
    the level of the whole ranking rather than one component.

AND THE RANK ITSELF CARRIES ITS OWN HONESTY. `n_indistinguishable` counts the
other candidates whose interval overlaps this one's. Where that number is large
the rank is an artefact of tie-breaking, and `rank_is_separated` says so per
row. This is the idiom `scripts/nac_rank.py::report` already uses when it prints
"1,239 molecules have an interval reaching this top-25 band" — a rank quoted
without it implies a precision no measurement here supports.

--------------------------------------------------------------------------
THE NORMALISATION, AND WHAT IT COSTS
--------------------------------------------------------------------------

The four components are an enrichment RATIO, a similarity, a free energy in
kJ/mol and an occupancy FRACTION. Nothing about those scales is commensurable,
and any mapping between them is a choice rather than a measurement, so the
choice is made here, in one place, and stated.

**Each component is mapped to its within-cohort empirical percentile**, computed
over the candidates that carry it.

Why that and not the alternatives:

  * **min-max** is destroyed by a single outlier — one 40x enrichment
    compresses everything else into the bottom decile. Docking scores and
    enrichment ratios produce such outliers reliably.
  * **z-score** is unbounded, so it cannot support the [0, 1] ignorance
    interval above, and it assumes a symmetry these distributions do not have.
  * **percentile** is robust by construction (an outlier occupies one rank slot,
    not the whole scale), is exactly bounded to [0, 1], and asserts nothing
    about the shape of any component's distribution — which is right, because
    no such assertion is supported for any of the four.

WHAT IT COSTS, stated rather than buried:

  * **Magnitude is discarded.** Two molecules a hair apart get adjacent
    percentiles as though the gap were real. That is why the interval is carried
    beside the score, and why pose quality is BANDED rather than scored (below).
  * **It is cohort-relative.** The same molecule scores differently in a
    different candidate set, and a component's percentile is taken over the
    candidates that CARRY it — so the worst of a pre-filtered 50-molecule BPMD
    cohort scores near 0 despite having been selected for measurement in the
    first place. `cohort_n_<component>` is written on every row so the number
    can be read as the within-cohort statement it is.

--------------------------------------------------------------------------
POSE QUALITY IS A BAND, NOT A SCORE (D0068)
--------------------------------------------------------------------------

`nac_criterion.enrichment` DOES NOT CONVERGE. D0068 measured the same 300
molecules falling from 2.91x to 0.96x at 10x the search effort, and the 15
crystallographic positives — never selected on score — falling identically
(2.27x -> 0.97x), so it is not winner's curse. The viable fraction is partly a
property of how hard you looked.

What survives the correction is that it is a **valid coarse filter**: top-300 vs
random warhead-matched inactives AUC 0.620 (p = 0.0017), and top-300 vs known
crystallographic binders AUC 0.438 (p = 0.21) — statistically indistinguishable
from known binders. What does not survive is the fine ordering: Spearman 0.364
between run counts, top-50 overlap 23/50, and at 200 runs 1,239 of 1,806
molecules' confidence intervals reach the top-25 band.

So this module treats it accordingly: pose quality enters as a **quantile band**,
not a percentile. `POSE_QUALITY_BANDS = 4` — the measurement supports "which
quarter of the distribution", and does not support more. Banding has two
consequences, both wanted:

  1. Molecules within a band TIE on this component. That is the honest reading
     of a metric whose intervals overlap that heavily, and it hands the ordering
     within a band to the components that can actually order.
  2. The band width (1 / POSE_QUALITY_BANDS) becomes this component's
     uncertainty interval, so a metric that cannot order finely visibly widens
     the composite rather than silently pretending to sharpen it.

**And nrun is part of the metric's definition, not a tuning knob** (D0068
consequence 2). `rank_candidates` REFUSES to score pose quality without it and
stamps `pose_quality_nrun` on every row. A 200-run and a 2,000-run enrichment
are different measurements and must never be ordered against each other.

--------------------------------------------------------------------------
DRUG-LIKENESS IS A FLOOR. PURCHASABILITY IS NOT AN INPUT.
--------------------------------------------------------------------------

@tt8804: "we want a baseline level of drug likeness but we dont want to be too
biased to it." A floor is not a ranking term, and in a list that may not drop
rows a floor can only be expressed as a FLAG — `druglike_floor_pass`. It has
weight zero by construction, it is not in `COMPONENTS`, and it cannot be given a
weight because the weight vector is validated against the component allowlist.
The threshold is deliberately permissive: covalent warheads depress QED on their
own, so a strict floor would rank on warhead class, which is D0020's failure
mode wearing a different hat.

**Purchasability must not affect the ranking at all.** It is refused as a
component by name (`_FORBIDDEN_COMPONENTS`) rather than merely omitted, because
omission is a default and this project's entire catalogue of defects is
defaults. Passthrough columns are untouched, so a purchasability column may
travel in the frame and be queried afterwards; it just cannot reach the score.

--------------------------------------------------------------------------
THE INTERFACE TO `shared/pose_consensus.py`
--------------------------------------------------------------------------

That module is being written separately and is NOT imported here — deliberately.
This layer takes consensus as DATA (a column of values, optionally a column of
spreads), so it has no build-time dependency on a module that may not exist yet,
and its absence is handled by the same mechanism that handles every other
missing measurement: the component is simply not declared, or is declared and
absent for a row, and the interval widens. See `CONSENSUS_CONTRACT` for the two
numbers this layer needs from it.

--------------------------------------------------------------------------
THE WEIGHTS ARE A CHOICE, AND ARE DECLARED AS ONE
--------------------------------------------------------------------------

No weighting of these four is validated on this target. BPMD and MD residence
have no discrimination measurement here at all; MD residence's only prior
measurement was that it is NOT reproducible (D0038/D0044); pose quality's is
AUC 0.672 at convergence (D0068) and is a filter rather than an order.

Fitting weights against our own 15 positives would be a knob turned until the
ranking looked reasonable, which is exactly what the rationale's
pre-registration clause forbids. So the default is EQUAL weights over the
declared components — a pre-registered default, not a fitted one — and `weights`
is an explicit parameter so any other choice is visible in the call site and in
`components_declared` rather than buried in a module constant.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class CompositeRankError(ValueError):
    """The ranking could not be built — never a silently degraded ranking."""


# ---------------------------------------------------------------------------
# The component registry. An ALLOWLIST with each component's direction declared
# explicitly, following `rank_shortlist.LOWER_IS_BETTER`, which is the guard
# that caught catalogue defect #4 and the only one that did. A component not
# named here is refused rather than assumed higher-is-better: "higher is better"
# is true of all four today and would be silently wrong for the first component
# quoted as a barrier height or a dissociation constant.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ComponentSpec:
    """One ranking component and everything needed to read its numbers.

    `banded` is not a display preference. It records that the underlying
    measurement supports a coarse ordering and not a fine one, and it changes
    both the score (a band midpoint) and the uncertainty (the band's width).
    """

    name: str
    higher_is_better: bool
    banded: bool
    n_bands: int
    requires: tuple[str, ...]
    note: str


# Quartiles. The measurement supports "which quarter", and D0068's numbers say
# it does not support more: Spearman 0.364 between run counts, top-50 overlap
# 23/50, and 1,239 of 1,806 intervals reaching one 25-molecule band.
POSE_QUALITY_BANDS = 4

COMPONENTS: dict[str, ComponentSpec] = {
    "pose_quality": ComponentSpec(
        name="pose_quality",
        higher_is_better=True,
        banded=True,
        n_bands=POSE_QUALITY_BANDS,
        requires=("nrun",),
        note="nac_criterion.enrichment. Does NOT converge (D0068) — a valid "
             "coarse filter, never a fine order, and meaningless without nrun.",
    ),
    "consensus": ComponentSpec(
        name="consensus",
        higher_is_better=True,
        banded=False,
        n_bands=0,
        requires=("n_poses", "tolerance_a"),
        note="pose_consensus.ConsensusResult.agreement. N and tolerance are "
             "part of its definition — pose_consensus.require_same_n says so "
             "and names THIS caller as the one with no other way to notice.",
    ),
    "bpmd": ComponentSpec(
        name="bpmd",
        higher_is_better=True,
        banded=False,
        n_bands=0,
        requires=(),
        note="bpmd.PoseStability.score — higher is harder to push the warhead "
             "out of the near-attack window. Its spread is spread_frac scaled "
             "by (1 + median_bias_to_escape_kj), because that is the factor "
             "the score multiplies the fraction by; spread_frac raw understates "
             "it for every pose that resisted a large bias.",
    ),
    "md_residence": ComponentSpec(
        name="md_residence",
        higher_is_better=True,
        banded=False,
        n_bands=0,
        requires=(),
        note="md_ensemble.residence_metrics frac_frames_engaged — occupancy, "
             "higher is longer in the pocket.",
    ),
}

# What this layer needs from `shared/pose_consensus.py`. Stated as prose rather
# than taken as an import: the two modules are built independently, a hard
# dependency would make this layer unusable until the other one lands, and
# absence is already handled correctly by the ignorance interval. Nothing here
# imports pose_consensus.
CONSENSUS_CONTRACT = """
Per candidate, from `pose_consensus.ConsensusResult`:

    value  -> `agreement`, the fraction of top-N pose pairs whose reactive
              regions agree within `tolerance_a`. HIGHER is better. Only the
              ORDER is read, never the units.
    bounds -> `agreement_jackknife`, already a (min, max) in the same units.
              Hand it over as a PAIR of columns rather than a symmetric spread:
              it is a leave-one-out range and is not symmetric about the value,
              and symmetrising it would invent precision on one side.

    provenance -> `n_poses` and `tolerance_a`, REQUIRED.

That last line is not bookkeeping. `pose_consensus.require_same_n` refuses to
compare agreements measured at different N or tolerance, and its docstring
names this module as the caller that "has no way to notice: both values are
populated, plausible, and in [0, 1]". So the requirement is enforced here, the
same way D0068's `nrun` is enforced for pose quality — declare the component
and you must declare what defines it.

Absence needs no signalling: an undeclared component is not part of the
ranking, and a declared component missing for a row widens that row's interval.
"""

# Refused BY NAME rather than merely left out of COMPONENTS. Omission is a
# default, and every defect in `docs/how_this_project_breaks.md` is a default.
# @tt8804: purchasability must not affect the ranking at all.
_FORBIDDEN_COMPONENTS = frozenset({
    "purchasable", "purchasability", "purchase", "in_stock", "availability",
    "available", "vendor", "supplier", "catalogue", "catalog", "price", "cost",
    "lead_time",
})

# TOTAL IGNORANCE ABOUT A COMPONENT, on the [0, 1] scale every component is
# mapped onto. Module-level so the claim "missing is not bad" is mutable and can
# therefore be mutation-tested: set this to (0.0, 0.0) and the ordering
# guarantee must visibly break.
ABSENT_INTERVAL: tuple[float, float] = (0.0, 1.0)

# QED floor. Deliberately permissive — a warhead depresses QED on its own, so a
# strict floor would rank warhead classes rather than molecules. A FLAG, never a
# term: it is not in COMPONENTS and so cannot be given a weight.
DRUGLIKE_FLOOR_QED = 0.35


# ---------------------------------------------------------------------------
# normalisation primitives
# ---------------------------------------------------------------------------

def _ecdf(sorted_values: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Empirical CDF of `sorted_values`, evaluated anywhere — not only at data.

    Midpoint convention: at a data point this returns exactly the average-rank
    percentile `(rank - 0.5) / n`, so a value and a value+0 spread agree, and
    ties receive one shared score rather than an order decided by row position.
    Evaluable off the data because a spread has to be pushed through the same
    transform as the value it belongs to.
    """
    n = sorted_values.size
    if n == 0:
        return np.full(np.shape(x), np.nan, dtype=float)
    lo = np.searchsorted(sorted_values, x, side="left")
    hi = np.searchsorted(sorted_values, x, side="right")
    return (lo + hi) / (2.0 * n)


@dataclass(frozen=True)
class ComponentInterval:
    """One component's contribution, on [0, 1], higher = better.

    `lo`/`hi` are the uncertainty interval and `present` says whether the row
    was measured at all. `pct` is the percentile of the VALUE ITSELF, ignoring
    uncertainty, and it exists because the two are not the same number:

        the composite uses the interval's MIDPOINT, so an asymmetric
        uncertainty deliberately shifts a component's contribution away from
        its value's percentile.

    That shift is wanted — an agreement whose jackknife says it could be much
    lower is worth less than one that could not, and averaging the asymmetry
    away would hide exactly the fragility `pose_consensus` computed it to
    expose. But it means `score` alone cannot be read as "where this molecule
    sits on this component", so `pct` is carried beside it and the two are
    separately auditable.

    `pct` is None for a BANDED component. Publishing a fine percentile next to
    a deliberately coarse band would hand a reader the precision D0068 measured
    as absent, in the same frame, one column over.
    """

    lo: pd.Series
    hi: pd.Series
    present: pd.Series
    pct: pd.Series | None


def component_interval(values: pd.Series, spec: ComponentSpec,
                       spread: pd.Series | None = None,
                       bounds: tuple[pd.Series, pd.Series] | None = None
                       ) -> ComponentInterval:
    """One component's per-row interval on [0, 1], higher = better.

    Rows where the component was not measured get `ABSENT_INTERVAL` and
    `present=False` — they are NOT given the cohort's worst value, and they are
    not dropped.

    Percentile rather than min-max or z-score: robust to the outliers every one
    of these components produces, exactly bounded (which the ignorance interval
    requires), and assumption-free about distribution shape. The cost is that
    magnitude is discarded and the score is cohort-relative; both are recorded
    in the module docstring and `cohort_n_<component>` on every row.

    UNCERTAINTY ARRIVES IN TWO SHAPES and both are supported, because forcing
    one into the other would invent precision:

      * `spread` — a symmetric ± in the component's own units (BPMD's
        replica-to-replica spread).
      * `bounds` — an explicit (lo, hi) in the component's own units, which is
        what `pose_consensus.agreement_jackknife` already is. A leave-one-out
        range is NOT symmetric about its value, and halving its width to fake a
        ± would overstate one side and understate the other.

    Either way the bounds are pushed through the SAME empirical CDF as the
    value, so the interval stays on the component's own percentile scale rather
    than being a unit-mixing addition after the fact.

    A BANDED component (pose quality, per D0068) is quantised to `n_bands`
    equal-population bands and its interval is the band it falls in. Its width
    is then the band width, so a metric that cannot order finely WIDENS the
    composite instead of pretending to sharpen it.
    """
    if spread is not None and bounds is not None:
        raise CompositeRankError(
            f"component {spec.name!r} was given both a spread and explicit "
            "bounds; they are two descriptions of one uncertainty and the "
            "wrong one would be used silently")
    v = pd.to_numeric(values, errors="coerce")
    present = v.notna()
    n = int(present.sum())

    lo = pd.Series(float(ABSENT_INTERVAL[0]), index=v.index, dtype=float)
    hi = pd.Series(float(ABSENT_INTERVAL[1]), index=v.index, dtype=float)
    if n == 0:
        log.warning("component %r is empty across the whole cohort — every row "
                    "is treated as unmeasured", spec.name)
        return ComponentInterval(lo, hi, present, None)

    # Direction is taken from the registry, never inferred from the data.
    sign = 1.0 if spec.higher_is_better else -1.0
    obs = (sign * v[present]).to_numpy(dtype=float)
    ordered = np.sort(obs)
    raw_pct = _ecdf(ordered, obs)

    if spec.banded:
        if spread is not None or bounds is not None:
            log.info("component %r is banded; its supplied uncertainty is not "
                     "used — the band width IS its uncertainty", spec.name)
        band = np.clip((raw_pct * spec.n_bands).astype(int), 0, spec.n_bands - 1)
        lo.loc[present] = band / spec.n_bands
        hi.loc[present] = (band + 1) / spec.n_bands
        # pct is deliberately NOT returned for a banded component: see
        # ComponentInterval. The band is the reading.
        return ComponentInterval(lo, hi, present, None)

    if bounds is not None:
        b_lo = pd.to_numeric(bounds[0], errors="coerce").reindex(v.index)[present]
        b_hi = pd.to_numeric(bounds[1], errors="coerce").reindex(v.index)[present]
        b_lo = b_lo.to_numpy(dtype=float)
        b_hi = b_hi.to_numpy(dtype=float)
        inverted = np.isfinite(b_lo) & np.isfinite(b_hi) & (b_lo > b_hi)
        if inverted.any():
            raise CompositeRankError(
                f"component {spec.name!r}: {int(inverted.sum())} rows have "
                "bounds with lo > hi. A pair the wrong way round produces a "
                "negative width and still ranks perfectly plausibly.")
        # A row missing a bound falls back to its own value — a point interval,
        # which is narrower than the truth but never wider, so it cannot
        # manufacture confidence it does not have on the optimistic side.
        b_lo = np.where(np.isfinite(b_lo), b_lo, v[present].to_numpy(dtype=float))
        b_hi = np.where(np.isfinite(b_hi), b_hi, v[present].to_numpy(dtype=float))
        # Negating for a lower-is-better component swaps which bound is which.
        obs_lo, obs_hi = (sign * b_lo, sign * b_hi) if sign > 0 else \
                         (sign * b_hi, sign * b_lo)
    else:
        if spread is None:
            log.warning("component %r supplied without an uncertainty — its "
                        "interval is a POINT, which claims a precision it "
                        "probably does not have", spec.name)
            s = np.zeros_like(obs)
        else:
            s = pd.to_numeric(spread, errors="coerce").reindex(v.index)[present]
            s = np.abs(s.to_numpy(dtype=float))
            bad = ~np.isfinite(s)
            if bad.any():
                log.warning("component %r: %d of %d rows have no usable spread; "
                            "those rows get a point interval", spec.name,
                            int(bad.sum()), n)
                s = np.where(bad, 0.0, s)
        obs_lo, obs_hi = obs - s, obs + s

    lo.loc[present] = _ecdf(ordered, obs_lo)
    hi.loc[present] = _ecdf(ordered, obs_hi)
    pct = pd.Series(np.nan, index=v.index, dtype=float)
    pct.loc[present] = raw_pct
    return ComponentInterval(lo, hi, present, pct)


# ---------------------------------------------------------------------------
# the ranked list
# ---------------------------------------------------------------------------

def _resolve_weights(declared: Sequence[str],
                     weights: Mapping[str, float] | None) -> dict[str, float]:
    """Equal by default; whatever is given is renormalised over the declared set.

    Renormalising over DECLARED rather than over all four is the distinction
    between "this component is not part of this ranking" and "this component
    could have been measured for you and was not". The first must not affect
    anyone; the second must widen that row's interval relative to rows that
    have it.
    """
    if weights is None:
        return {c: 1.0 / len(declared) for c in declared}
    unknown = set(weights) - set(declared)
    if unknown:
        raise CompositeRankError(
            f"weights name components that are not declared: {sorted(unknown)}. "
            "A weight on an undeclared component silently does nothing.")
    w = {c: float(weights.get(c, 0.0)) for c in declared}
    if any(x < 0 for x in w.values()):
        raise CompositeRankError(f"negative weight in {w}")
    total = sum(w.values())
    if total <= 0:
        raise CompositeRankError("weights sum to zero — nothing would be ranked")
    return {c: x / total for c, x in w.items()}


def _validate(df: pd.DataFrame, columns: Mapping[str, str],
              provenance: Mapping[str, Mapping[str, object]] | None) -> None:
    if not columns:
        raise CompositeRankError(
            "no components declared; a ranking with no components is not a "
            "ranking. Declare at least pose_quality.")
    for comp, col in columns.items():
        if comp in _FORBIDDEN_COMPONENTS:
            raise CompositeRankError(
                f"{comp!r} must not affect the ranking (@tt8804). Carry it as a "
                "passthrough column and query it afterwards.")
        if comp not in COMPONENTS:
            raise CompositeRankError(
                f"unknown component {comp!r}; known: {sorted(COMPONENTS)}. Add it "
                "to COMPONENTS with its direction rather than assuming one.")
        if col not in df.columns:
            raise CompositeRankError(f"component {comp!r} names column {col!r}, "
                                     f"which is not in the frame")
    for comp in columns:
        need = COMPONENTS[comp].requires
        got = dict((provenance or {}).get(comp, {}))
        missing = [k for k in need if got.get(k) is None]
        if missing:
            raise CompositeRankError(
                f"component {comp!r} requires provenance {list(need)} and is "
                f"missing {missing}. For pose_quality this is D0068: nrun is "
                "part of the metric's definition, not a tuning knob, and a "
                "200-run enrichment is not the same measurement as a 2,000-run "
                "one.")


def _indistinguishable(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """For each row, how many OTHER rows have an overlapping composite interval.

    Intervals are contiguous, so j overlaps i exactly when `lo_j <= hi_i` and
    `hi_j >= lo_i`, and the second set is a subset of the first — which turns
    an O(n^2) sweep into two binary searches per row. At 5,769 candidates the
    naive version is 33M comparisons for a number that gets recomputed on every
    query.

    This is the number that says whether a rank is real. Where it is large the
    ordering is tie-breaking, not measurement.
    """
    s_lo = np.sort(lo)
    s_hi = np.sort(hi)
    starts_before = np.searchsorted(s_lo, hi, side="right")
    ends_before = np.searchsorted(s_hi, lo, side="left")
    return np.maximum(starts_before - ends_before - 1, 0)


def rank_candidates(
    df: pd.DataFrame,
    *,
    columns: Mapping[str, str],
    spreads: Mapping[str, str | tuple[str, str]] | None = None,
    weights: Mapping[str, float] | None = None,
    provenance: Mapping[str, Mapping[str, object]] | None = None,
    druglike_col: str | None = "QED",
    druglike_floor: float = DRUGLIKE_FLOOR_QED,
) -> pd.DataFrame:
    """Combine the declared components into one ranked, complete candidate list.

    RANK AND RETAIN. The output has exactly the input's rows, in the input's
    order, with columns added. It never filters, drops, truncates or
    deduplicates, and it raises if the row count or index ever changes — the
    contract is asserted rather than trusted.

    `columns` maps component name -> column in `df`. Declaring a component makes
    it part of THIS ranking's definition; a row that is missing a declared
    component gets `ABSENT_INTERVAL` for it and a wider composite interval than
    a row that has it. An undeclared component affects nobody.

    `spreads` maps component name -> either a column of symmetric ± values, or
    a `(lo_col, hi_col)` pair of explicit bounds in the component's own units.
    Both shapes exist upstream — BPMD reports a replica spread,
    `pose_consensus` reports an asymmetric jackknife range — and forcing either
    into the other's shape would invent precision on one side.

    `provenance` supplies the parameters that are part of a component's
    DEFINITION rather than settings of it: `nrun` for pose quality (D0068),
    `n_poses` and `tolerance_a` for consensus (`pose_consensus.require_same_n`).
    They are required, stamped on every row, and the call is refused without
    them — two values measured at different settings are populated, plausible,
    and not comparable.

    Columns added, and why each is there:

        score_<c>, lo_<c>, hi_<c>, present_<c>
            per component, so the composite is auditable rather than a number
            that must be trusted. @tt8804 asked which components contributed;
            this answers it per row rather than per run.
        composite, composite_lo, composite_hi, composite_width
            the interval and its midpoint. `composite_width` is the confidence
            readout — wide means "mostly unmeasured", not "measured as middling".
        n_components, components_used, components_declared
            what was actually behind this row's number, and what could have
            been. A rank whose inputs are unknown is not auditable.
        composite_rank, n_indistinguishable, rank_is_separated
            the ordering, and whether it is real for this row.
        cohort_n_<c>
            how many candidates carried the component the percentile was taken
            over, because the score is a within-cohort statement.
        druglike_floor_pass
            a FLOOR expressed as a flag, weight zero. In a list that may not
            drop rows a flag is the only way to express one.

    Ties share a rank (`method="min"`), which they must: banding pose quality
    produces genuine ties by design, and breaking them on row position would
    invent an order out of the frame's storage layout.
    """
    if df.index.has_duplicates:
        raise CompositeRankError(
            "the frame's index has duplicates; per-row intervals could not be "
            "assigned unambiguously. Reset the index first.")
    _validate(df, columns, provenance)
    declared = sorted(columns)
    w = _resolve_weights(declared, weights)
    spreads = dict(spreads or {})
    unknown_spread = set(spreads) - set(declared)
    if unknown_spread:
        raise CompositeRankError(
            f"spreads name undeclared components: {sorted(unknown_spread)}")

    out = df.copy()
    out["components_declared"] = "+".join(declared)

    lo_tot = pd.Series(0.0, index=out.index, dtype=float)
    hi_tot = pd.Series(0.0, index=out.index, dtype=float)
    present_cols: list[pd.Series] = []

    for comp in declared:
        spec = COMPONENTS[comp]
        # `spreads[comp]` is either a column name (symmetric ±) or a pair of
        # column names (explicit lo, hi in the component's own units). Both
        # shapes occur in what the upstream modules actually produce, and
        # coercing one into the other would invent precision — see
        # `component_interval`.
        spec_unc = spreads.get(comp)
        sp_col, bnd_cols = None, None
        if isinstance(spec_unc, str):
            sp_col = spec_unc
        elif spec_unc is not None:
            bnd_cols = tuple(spec_unc)
            if len(bnd_cols) != 2:
                raise CompositeRankError(
                    f"uncertainty for {comp!r} must be a column name or a "
                    f"(lo, hi) pair of column names; got {spec_unc!r}")
        for c in ([sp_col] if sp_col else list(bnd_cols or ())):
            if c not in df.columns:
                raise CompositeRankError(
                    f"uncertainty column {c!r} for {comp!r} is not in the frame")
        ci = component_interval(
            df[columns[comp]], spec,
            spread=df[sp_col] if sp_col is not None else None,
            bounds=(df[bnd_cols[0]], df[bnd_cols[1]]) if bnd_cols else None)
        lo, hi, present = ci.lo, ci.hi, ci.present

        out[f"lo_{comp}"] = lo.where(present)
        out[f"hi_{comp}"] = hi.where(present)
        # The interval MIDPOINT, which is what feeds the composite — not the
        # value's own percentile. An asymmetric uncertainty moves the two
        # apart, deliberately; `pct_<comp>` carries the value's percentile
        # where publishing one is not itself misleading.
        out[f"score_{comp}"] = ((lo + hi) / 2.0).where(present)
        if ci.pct is not None:
            out[f"pct_{comp}"] = ci.pct.where(present)
        out[f"present_{comp}"] = present
        out[f"cohort_n_{comp}"] = int(present.sum())
        out[f"weight_{comp}"] = w[comp]
        for field, value in dict((provenance or {}).get(comp, {})).items():
            out[f"{comp}_{field}"] = value

        lo_tot += w[comp] * lo
        hi_tot += w[comp] * hi
        present_cols.append(present.rename(comp))

    pres = pd.concat(present_cols, axis=1)
    out["n_components"] = pres.sum(axis=1).astype("Int64")
    # The comparability class, spelled out per row. Two rows with the same
    # string were measured on the same things; two rows with different strings
    # were not, and comparing their composites means comparing an interval that
    # contains a measurement with one that contains ignorance.
    out["components_used"] = [
        "+".join(sorted(c for c in declared if row[c])) or "none"
        for _, row in pres.iterrows()]

    out["composite_lo"] = lo_tot
    out["composite_hi"] = hi_tot
    out["composite"] = (lo_tot + hi_tot) / 2.0
    out["composite_width"] = hi_tot - lo_tot

    out["composite_rank"] = (out["composite"]
                             .rank(ascending=False, method="min")
                             .astype("Int64"))
    n_ind = _indistinguishable(lo_tot.to_numpy(dtype=float),
                               hi_tot.to_numpy(dtype=float))
    out["n_indistinguishable"] = pd.Series(n_ind, index=out.index).astype("Int64")
    out["rank_is_separated"] = out["n_indistinguishable"] == 0

    # A FLOOR, NOT A TERM. Absent rather than True when it cannot be evaluated:
    # "not assessed" and "passed" are different facts, and defaulting to the
    # permissive one is catalogue defect #14 in miniature.
    if druglike_col and druglike_col in out.columns:
        q = pd.to_numeric(out[druglike_col], errors="coerce")
        out["druglike_floor_pass"] = (q >= druglike_floor).astype("boolean")
        out.loc[q.isna(), "druglike_floor_pass"] = pd.NA
        out["druglike_floor"] = druglike_floor
    else:
        if druglike_col:
            log.warning("no %r column; the drug-likeness floor is not assessed "
                        "and is reported as absent, not as passed", druglike_col)
        out["druglike_floor_pass"] = pd.Series(pd.NA, index=out.index,
                                               dtype="boolean")

    # THE CONTRACT, ASSERTED. @tt8804: "we dont want to delete any candidates."
    # Checked rather than trusted, because a silent drop would look exactly like
    # a shorter list.
    if len(out) != len(df) or not out.index.equals(df.index):
        raise CompositeRankError(
            f"rows were lost or reordered: {len(df)} in, {len(out)} out. This "
            "layer ranks and retains; it never filters.")
    log.info("ranked %d candidates on %s; median width %.3f, "
             "%d rows carry no component at all",
             len(out), "+".join(declared), float(out.composite_width.median()),
             int((out.n_components == 0).sum()))
    return out


def sorted_view(ranked: pd.DataFrame, identity_col: str | None = None
                ) -> pd.DataFrame:
    """The same rows in rank order — a VIEW, not a shortlist.

    Every row is still here. Sorting is separated from ranking so that
    `rank_candidates` can preserve the caller's row order for joins, and so that
    nothing in this module ever has an opportunity to return fewer rows than it
    was given.
    """
    keys = ["composite_rank"]
    if identity_col and identity_col in ranked.columns:
        keys.append(identity_col)      # deterministic within a tie
    out = ranked.sort_values(keys, kind="mergesort")
    if len(out) != len(ranked):
        raise CompositeRankError("sorting changed the row count")
    return out


def summarise(ranked: pd.DataFrame) -> str:
    """A short human-readable block, leading with what is NOT measured."""
    n = len(ranked)
    lines = [f"  {n} candidates ranked, {n} retained (this layer never drops)"]
    declared = str(ranked["components_declared"].iloc[0]) if n else ""
    lines.append(f"  components declared: {declared}")
    lines.append("")
    lines.append(f"  {'components measured':<32} {'n':>6} {'median width':>13}")
    lines.append("  " + "-" * 53)
    for used, g in ranked.groupby("components_used"):
        lines.append(f"  {used:<32} {len(g):>6} "
                     f"{float(g.composite_width.median()):>13.3f}")
    lines.append("")
    sep = int(ranked["rank_is_separated"].sum())
    lines.append(f"  {sep} of {n} rows have an interval disjoint from every "
                 f"other row's.")
    lines.append("  Where n_indistinguishable is large the rank is tie-breaking,")
    lines.append("  not measurement. A wide composite_width means MOSTLY")
    lines.append("  UNMEASURED — it does not mean measured as middling.")
    return "\n".join(lines)
