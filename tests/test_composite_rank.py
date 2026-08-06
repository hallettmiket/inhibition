"""
Purpose: the composite ranking is tested on synthetic components whose answer is known by arithmetic.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-06

WHY SYNTHETIC AND NOT REAL CANDIDATES. A real candidate's "correct" composite is
whatever this module computes, so testing against one would be circular. Every
cohort here is built so the expected ordering follows from the numbers put in —
which is what lets the test disagree with the implementation.

THE TWO CLAIMS THAT MATTER ARE MUTATION-TESTED, because a test that passes
against a broken implementation is worse than no test:

  * `test_MUTATION_absent_as_zero_breaks_the_missing_beats_bad_guarantee`
    flips `ABSENT_INTERVAL` to (0, 0) — the zero-imputation this module exists
    to refuse — and asserts the guarantee visibly inverts.
  * `test_MUTATION_purchasability_values_would_move_the_ranking_if_read`
    feeds the identical purchasability values in as a legitimate component and
    asserts the ranking DOES move, so the invariance test above it is not
    passing merely because the numbers were inert.
  * `test_MUTATION_removing_the_banding_breaks_the_tie_guarantee` widens
    pose quality to 1,000 bands and asserts the D0068 tie structure collapses.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd
import pytest

from shared import composite_rank as cr


# Every component's DEFINING parameters. Not settings — a 200-run enrichment
# and a 2,000-run one are different measurements (D0068), and so are agreements
# measured over 10 poses and over 40 (`pose_consensus.require_same_n`).
PROV = {"pose_quality": {"nrun": 2000},
        "consensus": {"n_poses": 20, "tolerance_a": 1.0}}


def frame(**cols: object) -> pd.DataFrame:
    """A candidate frame with an `ident` column, from equal-length columns."""
    n = len(next(iter(cols.values())))
    return pd.DataFrame({"ident": [f"c{i:04d}" for i in range(n)], **cols})


def rank(df: pd.DataFrame, **kw) -> pd.DataFrame:
    kw.setdefault("provenance", PROV)
    return cr.rank_candidates(df, **kw)


# --------------------------------------------------------------------------
# THE CONTRACT: rank and retain. Nothing is ever dropped.
# --------------------------------------------------------------------------

def test_every_input_row_survives_however_sparse():
    """@tt8804: "we dont want to delete any candidates."

    The cohort is deliberately as ragged as production will be: most rows carry
    only pose quality, a handful carry everything, and some carry nothing at
    all.
    """
    rng = np.random.default_rng(0)
    n = 500
    pq = rng.normal(2.0, 1.0, n)
    bp = np.where(rng.random(n) < 0.02, rng.random(n), np.nan)   # ~10 molecules
    md = np.where(rng.random(n) < 0.02, rng.random(n), np.nan)
    cons = np.where(rng.random(n) < 0.9, rng.random(n), np.nan)
    # 25 candidates that could not be measured on ANYTHING. They are the rows a
    # filter would take first and the ones the contract most needs to protect,
    # so they are put there deliberately rather than left to the random draw.
    naked = np.arange(0, n, 20)
    pq[naked] = np.nan
    bp[naked] = np.nan
    md[naked] = np.nan
    cons[naked] = np.nan
    df = frame(enrichment=pq, bpmd_score=bp, md_frac=md, consensus=cons)

    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "consensus",
                            "bpmd": "bpmd_score",
                            "md_residence": "md_frac"})
    assert len(out) == len(df)
    assert out.index.equals(df.index)
    assert list(out.ident) == list(df.ident)
    assert out.composite.notna().all(), "every row must carry a composite"
    assert out.composite_rank.notna().all()
    assert (out.n_components == 0).sum() == len(naked), \
        "the rows with nothing measured are exactly the ones a filter drops"
    assert set(out.loc[naked, "ident"]) == set(df.loc[naked, "ident"])


def test_a_row_with_no_component_at_all_is_retained_and_says_so():
    """Unmeasured on everything is a row, not an absence.

    It sits at the midpoint with the maximum possible width — the width is the
    statement, and it says "we know nothing", not "we measured middling".
    """
    df = frame(enrichment=[3.0, np.nan], bpmd_score=[0.9, np.nan])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    naked = out.iloc[1]
    assert naked.n_components == 0
    assert naked.components_used == "none"
    assert naked.composite == pytest.approx(0.5)
    assert naked.composite_width == pytest.approx(1.0)


def test_sorted_view_is_a_view_and_not_a_shortlist():
    rng = np.random.default_rng(1)
    df = frame(enrichment=rng.random(300), bpmd_score=rng.random(300))
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    view = cr.sorted_view(out, identity_col="ident")
    assert len(view) == len(out)
    assert set(view.ident) == set(out.ident)
    assert view.composite_rank.is_monotonic_increasing


# --------------------------------------------------------------------------
# all four components present
# --------------------------------------------------------------------------

def test_all_components_present_orders_by_the_components():
    """A candidate better on every component must rank first, and vice versa.

    Nothing subtle — this is the sanity floor the rest of the file stands on.
    """
    df = frame(enrichment=[4.0, 2.5, 1.0],
               consensus=[0.9, 0.5, 0.1],
               bpmd_score=[0.9, 0.5, 0.1],
               md_frac=[0.9, 0.5, 0.1])
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "consensus",
                            "bpmd": "bpmd_score",
                            "md_residence": "md_frac"})
    assert list(out.composite_rank) == [1, 2, 3]
    assert (out.n_components == 4).all()
    assert (out.components_used == "bpmd+consensus+md_residence+pose_quality").all()


def test_a_fully_measured_row_is_narrower_than_a_barely_measured_one():
    """`composite_width` is the confidence readout, so it must track evidence."""
    df = frame(enrichment=[2.0, 2.0, 2.0],
               consensus=[0.5, 0.5, np.nan],
               bpmd_score=[0.5, np.nan, np.nan],
               md_frac=[0.5, np.nan, np.nan])
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "consensus",
                            "bpmd": "bpmd_score",
                            "md_residence": "md_frac"})
    assert out.composite_width.iloc[0] < out.composite_width.iloc[1]
    assert out.composite_width.iloc[1] < out.composite_width.iloc[2]
    assert list(out.n_components) == [4, 2, 1]


# --------------------------------------------------------------------------
# only pose quality present — the normal production row
# --------------------------------------------------------------------------

def test_only_pose_quality_present_still_ranks_and_is_labelled():
    """The 5,700-odd candidates that will never see a GPU-week.

    They rank, they are retained, and every row states that its number came from
    one component of four. The width is the same for all of them because the
    same three things are unknown for all of them.
    """
    rng = np.random.default_rng(2)
    df = frame(enrichment=rng.random(400) * 5,
               consensus=[np.nan] * 400,
               bpmd_score=[np.nan] * 400,
               md_frac=[np.nan] * 400)
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "consensus",
                            "bpmd": "bpmd_score",
                            "md_residence": "md_frac"})
    assert (out.n_components == 1).all()
    assert (out.components_used == "pose_quality").all()
    # three absent components at weight 0.25 each, plus one quartile band
    expected = 3 * 0.25 * 1.0 + 0.25 * (1 / cr.POSE_QUALITY_BANDS)
    assert out.composite_width.round(6).nunique() == 1
    assert out.composite_width.iloc[0] == pytest.approx(expected)
    # and the ordering is the band ordering, nothing else
    assert out.composite.nunique() == cr.POSE_QUALITY_BANDS


def test_pose_quality_alone_still_separates_the_top_band_from_the_bottom():
    df = frame(enrichment=[0.1, 0.4, 1.2, 1.6, 2.4, 2.9, 4.0, 9.0])
    out = rank(df, columns={"pose_quality": "enrichment"})
    assert out.composite_rank.iloc[-1] == 1, "highest enrichment in the top band"
    assert out.composite_rank.iloc[0] > out.composite_rank.iloc[-1]


# --------------------------------------------------------------------------
# THE CENTRAL CLAIM: a missing component is not a bad one
# --------------------------------------------------------------------------

def bpmd_cohort() -> pd.DataFrame:
    """20 molecules with BPMD spanning the range, plus one never measured.

    Pose quality is identical for all 21, so it contributes the same band to
    everyone and the ordering is entirely BPMD's — which is what makes the
    expected answer arithmetic rather than opinion.
    """
    vals = list(np.linspace(0.0, 0.95, 20)) + [np.nan]
    return frame(enrichment=[2.0] * 21, bpmd_score=vals)


def test_a_missing_component_ranks_above_a_bad_one_and_below_a_good_one():
    """@tt8804's constraint, stated as an ordering and checked as one.

    With two components at equal weight, the unmeasured molecule contributes
    the midpoint of [0, 1] for BPMD, so it must sit exactly where the median of
    the measured cohort sits: above the ten worst, below the ten best. Not
    "roughly in the middle" — exactly ten, which is a claim the implementation
    can fail.
    """
    out = rank(bpmd_cohort(),
               columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    absent = out.iloc[20]
    measured = out.iloc[:20]

    assert absent.n_components == 1
    assert absent.components_used == "pose_quality"
    assert (measured.composite > absent.composite).sum() == 10
    assert (measured.composite < absent.composite).sum() == 10
    assert absent.composite_rank < measured.composite_rank.max(), \
        "not measured must not rank below measured badly"
    assert absent.composite_rank > measured.composite_rank.min()


def test_the_worst_measured_molecule_ranks_last_and_the_unmeasured_one_does_not():
    """The failure mode named in the brief, checked at the extremes."""
    out = rank(bpmd_cohort(),
               columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    worst_measured = out.iloc[0]        # bpmd 0.0
    absent = out.iloc[20]
    assert worst_measured.composite_rank == out.composite_rank.max()
    assert absent.composite_rank != out.composite_rank.max()


def test_a_missing_component_widens_rather_than_lowering():
    """The mechanism behind the claim above, checked directly.

    Absence must move the interval's WIDTH, not its floor. If absence lowered
    the composite the width would be unchanged and the midpoint would fall —
    which is imputation wearing an interval's clothes.
    """
    out = rank(bpmd_cohort(),
               columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    absent, measured = out.iloc[20], out.iloc[:20]
    assert absent.composite_width > measured.composite_width.max()
    assert absent.composite_lo <= measured.composite_lo.min()
    assert absent.composite_hi >= measured.composite_hi.max()


def test_MUTATION_absent_as_zero_breaks_the_missing_beats_bad_guarantee(monkeypatch):
    """Zero-imputation is the mistake; this proves the test above would catch it.

    Setting ABSENT_INTERVAL to (0, 0) is exactly "impute the worst value", the
    thing @tt8804 ruled out by name. Under it the unmeasured molecule falls to
    LAST — so the guarantee is genuinely load-bearing and the assertions above
    are not passing by accident.
    """
    monkeypatch.setattr(cr, "ABSENT_INTERVAL", (0.0, 0.0))
    out = rank(bpmd_cohort(),
               columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    absent = out.iloc[20]
    assert absent.composite_rank == out.composite_rank.max(), \
        "the mutation must reproduce the defect, or the guarantee is untested"
    assert absent.composite_width == pytest.approx(
        out.iloc[0].composite_width), "and the ignorance would vanish silently"


# --------------------------------------------------------------------------
# ties
# --------------------------------------------------------------------------

def test_identical_candidates_share_a_rank():
    """Breaking a genuine tie on row position invents an order out of storage."""
    df = frame(enrichment=[2.0, 2.0, 5.0], bpmd_score=[0.4, 0.4, 0.9])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    assert out.composite_rank.iloc[0] == out.composite_rank.iloc[1]
    assert out.composite.iloc[0] == pytest.approx(out.composite.iloc[1])


def test_banding_makes_pose_quality_tie_by_design_D0068():
    """D0068: enrichment is a coarse filter and cannot order what it selects.

    100 distinct enrichments must collapse to exactly `POSE_QUALITY_BANDS`
    distinct scores. Anything finer would be a rank order the measurement does
    not carry — Spearman 0.364 between run counts, top-50 overlap 23/50.
    """
    df = frame(enrichment=np.linspace(0.1, 8.0, 100))
    out = rank(df, columns={"pose_quality": "enrichment"})
    assert out["score_pose_quality"].nunique() == cr.POSE_QUALITY_BANDS
    assert out.composite_rank.nunique() == cr.POSE_QUALITY_BANDS
    assert out.groupby("composite_rank").size().tolist() == [25, 25, 25, 25]


def test_MUTATION_removing_the_banding_breaks_the_tie_guarantee(monkeypatch):
    """Widening pose quality to 1,000 bands is "treat it as a fine score".

    The tie structure the test above asserts then collapses to 100 distinct
    ranks, which is the precision D0068 measured as absent.
    """
    spec = dataclasses.replace(cr.COMPONENTS["pose_quality"], n_bands=1000)
    monkeypatch.setitem(cr.COMPONENTS, "pose_quality", spec)
    df = frame(enrichment=np.linspace(0.1, 8.0, 100))
    out = rank(df, columns={"pose_quality": "enrichment"})
    assert out.composite_rank.nunique() == 100, \
        "the mutation must reproduce the over-precision, or banding is untested"


def test_n_indistinguishable_flags_ranks_that_are_only_tie_breaking():
    """A rank whose interval overlaps 399 others is not a measurement.

    Nobody here has been measured on BPMD, so every interval is half the scale
    wide and the four pose-quality bands cannot be told apart. The rank column
    still populates — it has to, the list is queryable — and
    `rank_is_separated` is what stops it being read as a result.
    """
    df = frame(enrichment=np.linspace(0.1, 8.0, 400), bpmd_score=[np.nan] * 400)
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    assert (out.n_indistinguishable == 399).all()
    assert not out.rank_is_separated.any()
    assert out.composite_rank.notna().all(), "still ranked, just not separated"


def test_a_well_measured_candidate_can_be_genuinely_separated():
    """The flag must be able to say yes, or it is not carrying information."""
    df = frame(enrichment=[0.1, 0.2, 0.3, 9.0],
               bpmd_score=[0.1, 0.2, 0.3, 0.99],
               consensus=[0.1, 0.2, 0.3, 0.99],
               md_frac=[0.1, 0.2, 0.3, 0.99])
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "consensus",
                            "bpmd": "bpmd_score",
                            "md_residence": "md_frac"})
    assert bool(out.rank_is_separated.iloc[3]) is True
    assert out.composite_rank.iloc[3] == 1


def test_n_indistinguishable_matches_a_brute_force_count():
    """The searchsorted shortcut must agree with the O(n^2) definition."""
    rng = np.random.default_rng(3)
    n = 200
    df = frame(enrichment=rng.random(n) * 5,
               bpmd_score=np.where(rng.random(n) < 0.4, rng.random(n), np.nan))
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    lo = out.composite_lo.to_numpy()
    hi = out.composite_hi.to_numpy()
    brute = [(int(((lo <= hi[i]) & (hi >= lo[i])).sum()) - 1) for i in range(n)]
    assert out.n_indistinguishable.tolist() == brute


# --------------------------------------------------------------------------
# purchasability must not reach the score
# --------------------------------------------------------------------------

def purchasability_cohort():
    rng = np.random.default_rng(4)
    n = 120
    return frame(enrichment=rng.random(n) * 5,
                 bpmd_score=np.where(rng.random(n) < 0.5, rng.random(n), np.nan),
                 purchasable=rng.random(n))


def test_purchasability_is_refused_as_a_component_by_name():
    """Refused, not merely absent — omission is a default, and defaults are the
    shape every defect in this project has taken."""
    df = purchasability_cohort()
    with pytest.raises(cr.CompositeRankError, match="must not affect the ranking"):
        rank(df, columns={"pose_quality": "enrichment",
                          "purchasable": "purchasable"})


def test_permuting_purchasability_leaves_the_ranking_bit_identical():
    df = purchasability_cohort()
    cols = {"pose_quality": "enrichment", "bpmd": "bpmd_score"}
    base = rank(df, columns=cols)

    shuffled = df.copy()
    shuffled["purchasable"] = (
        df.purchasable.sample(frac=1.0, random_state=7).to_numpy())
    moved = rank(shuffled, columns=cols)

    assert not shuffled.purchasable.equals(df.purchasable), "sanity: it moved"
    pd.testing.assert_series_equal(base.composite, moved.composite)
    pd.testing.assert_series_equal(base.composite_rank, moved.composite_rank)


def test_MUTATION_purchasability_values_would_move_the_ranking_if_read():
    """The invariance test above is only meaningful if the values carry order.

    Feeding the identical purchasability numbers in as a legitimate component
    must change the ranking. If it did not, the test above would be passing
    because the column was inert rather than because it was ignored.
    """
    df = purchasability_cohort()
    cols = {"pose_quality": "enrichment", "bpmd": "bpmd_score"}
    base = rank(df, columns=cols)
    shuffled = df.copy()
    shuffled["purchasable"] = (
        df.purchasable.sample(frac=1.0, random_state=7).to_numpy())

    as_component = {**cols, "consensus": "purchasable"}
    a = rank(df, columns=as_component)
    b = rank(shuffled, columns=as_component)
    assert not a.composite_rank.equals(b.composite_rank), \
        "these values DO carry ordering information — so ignoring them is a result"


# --------------------------------------------------------------------------
# drug-likeness is a floor, not a term
# --------------------------------------------------------------------------

def test_druglike_floor_is_a_flag_and_does_not_move_the_ranking():
    """@tt8804: a baseline level, without being biased toward it.

    A candidate below the floor is FLAGGED and kept — which in a list that may
    not drop rows is the only way to express a floor at all.
    """
    rng = np.random.default_rng(5)
    n = 80
    df = frame(enrichment=rng.random(n) * 5,
               bpmd_score=rng.random(n),
               QED=rng.random(n))
    cols = {"pose_quality": "enrichment", "bpmd": "bpmd_score"}
    base = rank(df, columns=cols)

    permuted = df.copy()
    permuted["QED"] = df.QED.sample(frac=1.0, random_state=11).to_numpy()
    moved = rank(permuted, columns=cols)

    pd.testing.assert_series_equal(base.composite_rank, moved.composite_rank)
    assert base.druglike_floor_pass.notna().all()
    assert (~base.druglike_floor_pass.astype(bool)).any(), "some must fail"
    assert len(base) == n, "and every one of them is still here"


def test_druglike_floor_cannot_be_given_a_weight():
    df = frame(enrichment=[1.0, 2.0], QED=[0.2, 0.8])
    with pytest.raises(cr.CompositeRankError, match="unknown component"):
        rank(df, columns={"pose_quality": "enrichment", "QED": "QED"})


def test_an_unassessed_floor_is_absent_not_passed():
    """"Not assessed" and "passed" are different facts; defaulting to the
    permissive one is catalogue defect #14 in miniature."""
    df = frame(enrichment=[1.0, 2.0])
    out = rank(df, columns={"pose_quality": "enrichment"})
    assert out.druglike_floor_pass.isna().all()

    df2 = frame(enrichment=[1.0, 2.0], QED=[np.nan, 0.9])
    out2 = rank(df2, columns={"pose_quality": "enrichment"})
    assert pd.isna(out2.druglike_floor_pass.iloc[0])
    assert bool(out2.druglike_floor_pass.iloc[1]) is True


# --------------------------------------------------------------------------
# provenance: nrun is part of the metric's definition (D0068)
# --------------------------------------------------------------------------

def test_pose_quality_without_nrun_is_refused():
    """D0068 consequence 2. A 200-run and a 2,000-run enrichment are different
    measurements, and an enrichment without its run count cannot be ordered
    against anything."""
    df = frame(enrichment=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError, match="nrun"):
        cr.rank_candidates(df, columns={"pose_quality": "enrichment"},
                           provenance=None)
    with pytest.raises(cr.CompositeRankError, match="nrun"):
        cr.rank_candidates(df, columns={"pose_quality": "enrichment"},
                           provenance={"pose_quality": {"nrun": None}})


def test_nrun_is_stamped_on_every_row():
    df = frame(enrichment=[1.0, 2.0, 3.0])
    out = rank(df, columns={"pose_quality": "enrichment"})
    assert (out.pose_quality_nrun == 2000).all()


def test_consensus_without_its_n_and_tolerance_is_refused():
    """`pose_consensus.require_same_n` names THIS module as the caller that
    otherwise has no way to notice — "both values are populated, plausible, and
    in [0, 1]". So the requirement is enforced here rather than hoped for."""
    df = frame(enrichment=[1.0, 2.0], agreement=[0.4, 0.8])
    cols = {"pose_quality": "enrichment", "consensus": "agreement"}
    with pytest.raises(cr.CompositeRankError, match="n_poses"):
        cr.rank_candidates(df, columns=cols,
                           provenance={"pose_quality": {"nrun": 2000}})
    with pytest.raises(cr.CompositeRankError, match="tolerance_a"):
        cr.rank_candidates(
            df, columns=cols,
            provenance={"pose_quality": {"nrun": 2000},
                        "consensus": {"n_poses": 20}})
    out = rank(df, columns=cols)
    assert (out.consensus_n_poses == 20).all()
    assert (out.consensus_tolerance_a == 1.0).all()


# --------------------------------------------------------------------------
# uncertainty arrives in two shapes, and neither is forced into the other
# --------------------------------------------------------------------------

def test_an_asymmetric_jackknife_range_stays_asymmetric():
    """`pose_consensus.agreement_jackknife` is a leave-one-out (min, max).

    It is not symmetric about the value, and halving its width to fake a ±
    would overstate one side and understate the other. Here the middle
    candidate's range runs far below its value and barely above it, and the
    resulting interval must lean the same way.
    """
    df = frame(enrichment=[1.0, 2.0, 3.0, 4.0, 5.0],
               agreement=[0.10, 0.30, 0.50, 0.70, 0.90],
               agree_lo=[0.09, 0.29, 0.12, 0.69, 0.89],   # a long lower tail
               agree_hi=[0.11, 0.31, 0.51, 0.71, 0.91])
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "agreement"},
               spreads={"consensus": ("agree_lo", "agree_hi")})
    mid = out.iloc[2]
    # measured against the VALUE's own percentile, not the interval midpoint —
    # the midpoint is symmetric within its interval by construction, which is
    # exactly why `pct_<c>` is carried separately.
    below = mid["pct_consensus"] - mid["lo_consensus"]
    above = mid["hi_consensus"] - mid["pct_consensus"]
    assert below > 3 * above, "the asymmetry must survive the transform"
    assert 0.0 <= mid["lo_consensus"] <= mid["hi_consensus"] <= 1.0

    # and the long lower tail must PULL THE CONTRIBUTION DOWN. An estimate that
    # could be much lower is worth less than one that could not, and averaging
    # the asymmetry away would hide the fragility the jackknife exists to show.
    assert mid["score_consensus"] < mid["pct_consensus"]


def test_a_banded_component_publishes_no_finer_percentile_beside_its_band():
    """D0068's precision, withheld rather than merely unused.

    A fine `pct_pose_quality` sitting one column over from a deliberately
    coarse band is the knob D0068 measured as non-existent, handed to the
    reader in the same frame.
    """
    df = frame(enrichment=np.linspace(0.1, 8.0, 40), agreement=np.linspace(0, 1, 40))
    out = rank(df, columns={"pose_quality": "enrichment",
                            "consensus": "agreement"})
    assert "pct_pose_quality" not in out.columns
    assert "score_pose_quality" in out.columns
    assert "pct_consensus" in out.columns, "the unbanded one still publishes it"


def test_bounds_the_wrong_way_round_are_refused():
    """A pair reversed produces a negative width and still ranks plausibly."""
    df = frame(enrichment=[1.0, 2.0], agreement=[0.4, 0.8],
               agree_lo=[0.9, 0.7], agree_hi=[0.1, 0.9])
    with pytest.raises(cr.CompositeRankError, match="lo > hi"):
        rank(df, columns={"pose_quality": "enrichment",
                          "consensus": "agreement"},
             spreads={"consensus": ("agree_lo", "agree_hi")})


def test_a_spread_and_bounds_together_are_refused():
    """Two descriptions of one uncertainty; the wrong one would be used silently."""
    values = pd.Series([1.0, 2.0, 3.0])
    with pytest.raises(cr.CompositeRankError, match="both a spread and"):
        cr.component_interval(values, cr.COMPONENTS["bpmd"],
                              spread=pd.Series([0.1] * 3),
                              bounds=(pd.Series([0.9] * 3), pd.Series([1.1] * 3)))


def test_a_malformed_bounds_pair_is_refused():
    df = frame(enrichment=[1.0, 2.0], agreement=[0.4, 0.8], a=[0.1, 0.2])
    with pytest.raises(cr.CompositeRankError, match=r"\(lo, hi\) pair"):
        rank(df, columns={"pose_quality": "enrichment",
                          "consensus": "agreement"},
             spreads={"consensus": ("a",)})


def test_cohort_size_travels_with_the_score():
    """The percentile is a within-cohort statement, so the cohort is recorded."""
    df = frame(enrichment=[1.0, 2.0, 3.0, 4.0],
               bpmd_score=[0.5, np.nan, np.nan, 0.9])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    assert (out.cohort_n_pose_quality == 4).all()
    assert (out.cohort_n_bpmd == 2).all()


def test_declared_components_are_recorded_even_when_nobody_has_them():
    """Declared-but-unmeasured and undeclared are different, and the difference
    is visible: the first widens everyone's interval, the second affects no one."""
    df = frame(enrichment=[1.0, 5.0], bpmd_score=[np.nan, np.nan])
    with_bpmd = rank(df, columns={"pose_quality": "enrichment",
                                  "bpmd": "bpmd_score"})
    without = rank(df, columns={"pose_quality": "enrichment"})
    assert with_bpmd.composite_width.iloc[0] > without.composite_width.iloc[0]
    assert with_bpmd.components_declared.iloc[0] == "bpmd+pose_quality"
    assert without.components_declared.iloc[0] == "pose_quality"
    # and the ORDER is unchanged, because the absence is uniform
    assert list(with_bpmd.composite_rank) == list(without.composite_rank)


# --------------------------------------------------------------------------
# normalisation: robust, bounded, and honest about its direction
# --------------------------------------------------------------------------

def test_a_single_outlier_does_not_compress_the_rest_of_the_scale():
    """Why percentile and not min-max.

    Under min-max the 1e9 candidate would push every other score below 1e-8.
    Under a percentile it takes one rank slot and the rest keep their spread.
    """
    df = frame(enrichment=[1.0] * 6,
               bpmd_score=[1.0, 2.0, 3.0, 4.0, 5.0, 1e9])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    s = out["score_bpmd"]
    assert s.iloc[4] == pytest.approx(9 / 12)      # the value 5.0, second-best
    assert s.min() == pytest.approx(1 / 12)
    assert s.max() == pytest.approx(11 / 12)


def test_scores_stay_inside_the_unit_interval_whatever_the_spread():
    """The ignorance interval requires a bounded scale; a spread must not leave it."""
    df = frame(enrichment=[1.0, 2.0, 3.0],
               bpmd_score=[0.1, 0.5, 0.9], bpmd_spread=[99.0, 99.0, 99.0])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"},
               spreads={"bpmd": "bpmd_spread"})
    assert out["lo_bpmd"].min() >= 0.0
    assert out["hi_bpmd"].max() <= 1.0
    assert out.composite.between(0.0, 1.0).all()


def test_a_larger_spread_gives_a_wider_interval():
    """Uncertainty must survive into the output rather than being averaged away."""
    base = frame(enrichment=[1.0] * 5, bpmd_score=[0.1, 0.3, 0.5, 0.7, 0.9],
                 bpmd_spread=[0.0] * 5)
    wide = base.copy()
    wide["bpmd_spread"] = 0.3
    cols = {"pose_quality": "enrichment", "bpmd": "bpmd_score"}
    a = rank(base, columns=cols, spreads={"bpmd": "bpmd_spread"})
    b = rank(wide, columns=cols, spreads={"bpmd": "bpmd_spread"})
    assert (b.composite_width >= a.composite_width).all()
    assert b.composite_width.iloc[2] > a.composite_width.iloc[2]


def test_direction_comes_from_the_registry_not_from_the_data():
    """`rank_shortlist.LOWER_IS_BETTER` is the guard that caught catalogue #4.

    Declaring a component lower-is-better must invert its score. Checked on a
    fabricated spec so the assertion is about the mechanism, not about which
    way today's four components happen to point.
    """
    values = pd.Series([1.0, 2.0, 3.0])
    up = dataclasses.replace(cr.COMPONENTS["bpmd"], higher_is_better=True)
    down = dataclasses.replace(cr.COMPONENTS["bpmd"], higher_is_better=False)
    u = cr.component_interval(values, up)
    d = cr.component_interval(values, down)
    assert u.lo.iloc[2] > u.lo.iloc[0]
    assert d.lo.iloc[2] < d.lo.iloc[0]
    assert list(u.hi) == list(reversed(list(d.hi)))
    assert list(u.pct) == list(reversed(list(d.pct)))


def test_bounds_respect_the_declared_direction():
    """Negating for a lower-is-better component swaps which bound is which.

    Get this wrong and lo > hi silently, producing a negative width on a row
    that otherwise looks entirely normal.
    """
    values = pd.Series([1.0, 2.0, 3.0])
    b = (pd.Series([0.5, 1.5, 2.5]), pd.Series([1.5, 2.5, 3.5]))
    down = dataclasses.replace(cr.COMPONENTS["bpmd"], higher_is_better=False)
    ci = cr.component_interval(values, down, bounds=b)
    assert (ci.hi >= ci.lo).all(), "a reversed bound would go unnoticed"
    assert ci.lo.iloc[0] > ci.lo.iloc[2], "and the direction still inverts"


def test_all_registered_components_declare_a_direction():
    """The registry is an allowlist; every member states its direction."""
    for name, spec in cr.COMPONENTS.items():
        assert spec.name == name
        assert isinstance(spec.higher_is_better, bool)
        assert spec.note, f"{name} carries no explanation of what it is"


def test_a_component_that_nobody_carries_is_handled_not_crashed():
    df = frame(enrichment=[1.0, 2.0], md_frac=[np.nan, np.nan])
    out = rank(df, columns={"pose_quality": "enrichment",
                            "md_residence": "md_frac"})
    assert out["score_md_residence"].isna().all()
    assert (out.cohort_n_md_residence == 0).all()
    assert (~out["present_md_residence"]).all()


# --------------------------------------------------------------------------
# weights
# --------------------------------------------------------------------------

def test_default_weights_are_equal_over_the_declared_components():
    """Pre-registered, not fitted. No weighting of these four is validated on
    this target, and fitting one against 15 positives would be a knob."""
    df = frame(enrichment=[1.0, 2.0], bpmd_score=[0.2, 0.8])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    assert out.weight_pose_quality.iloc[0] == pytest.approx(0.5)
    assert out.weight_bpmd.iloc[0] == pytest.approx(0.5)


def test_weights_are_renormalised_and_recorded():
    df = frame(enrichment=[1.0, 2.0], bpmd_score=[0.2, 0.8])
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"},
               weights={"pose_quality": 1.0, "bpmd": 3.0})
    assert out.weight_pose_quality.iloc[0] == pytest.approx(0.25)
    assert out.weight_bpmd.iloc[0] == pytest.approx(0.75)


def test_a_weight_on_an_undeclared_component_is_refused():
    """A weight that silently does nothing is worse than a missing one."""
    df = frame(enrichment=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError, match="not declared"):
        rank(df, columns={"pose_quality": "enrichment"},
             weights={"pose_quality": 1.0, "bpmd": 1.0})


@pytest.mark.parametrize("bad", [{"pose_quality": -1.0}, {"pose_quality": 0.0}])
def test_degenerate_weights_are_refused(bad):
    df = frame(enrichment=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError):
        rank(df, columns={"pose_quality": "enrichment"}, weights=bad)


# --------------------------------------------------------------------------
# the guards that can actually fail
# --------------------------------------------------------------------------

def test_an_unknown_component_is_refused():
    df = frame(enrichment=[1.0, 2.0], mystery=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError, match="unknown component"):
        rank(df, columns={"pose_quality": "enrichment", "mystery": "mystery"})


def test_a_missing_column_is_refused_rather_than_silently_absent():
    df = frame(enrichment=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError, match="not in the frame"):
        rank(df, columns={"pose_quality": "enrichment", "bpmd": "no_such_col"})


def test_declaring_nothing_is_refused():
    df = frame(enrichment=[1.0, 2.0])
    with pytest.raises(cr.CompositeRankError, match="no components declared"):
        rank(df, columns={})


def test_a_duplicated_index_is_refused():
    df = frame(enrichment=[1.0, 2.0, 3.0])
    df.index = [0, 0, 1]
    with pytest.raises(cr.CompositeRankError, match="duplicates"):
        rank(df, columns={"pose_quality": "enrichment"})


def test_a_spread_for_an_undeclared_component_is_refused():
    df = frame(enrichment=[1.0, 2.0], s=[0.1, 0.1])
    with pytest.raises(cr.CompositeRankError, match="undeclared"):
        rank(df, columns={"pose_quality": "enrichment"}, spreads={"bpmd": "s"})


def test_summarise_leads_with_what_is_not_measured():
    rng = np.random.default_rng(6)
    n = 50
    df = frame(enrichment=rng.random(n),
               bpmd_score=np.where(rng.random(n) < 0.3, rng.random(n), np.nan))
    out = rank(df, columns={"pose_quality": "enrichment", "bpmd": "bpmd_score"})
    text = cr.summarise(out)
    assert "never drops" in text
    assert "pose_quality" in text
    assert "UNMEASURED" in text
