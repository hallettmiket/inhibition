"""`engagement_rank` -- ranking modes on warhead geometry alone.

The properties a caller relies on, not the numbers:

  * the score is a POSE property attributed to a group, never an average over a
    group -- 93% of shipped modes span more than half the anchor-quality scale,
    which is why every mode-level aggregate scored rho ~ 0.11 against the
    simulated pose's +0.652;
  * an unknown statistic or aggregation RAISES rather than falling back to a
    default, so a typo cannot silently produce a differently-ordered table;
  * the pose-count gate is not applied by default, because an engagement score
    is estimable in a group of one where a frequency is not;
  * `n_modes` travels with any ligand-level mean, because that denominator moves
    with docking depth (D0092).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from shared import engagement_rank as er


def _poses(spec):
    """spec: {(ident, mode): [(distance, angle), ...]} -> a per-pose frame."""
    rows = []
    for (ident, mode), pairs in spec.items():
        for i, (d, a) in enumerate(pairs):
            rows.append(dict(ident=ident, mode=mode, pose_idx=i, distance=d,
                             angle=a, mechanism="sn2_displacement", warhead_class="wc1"))
    return pd.DataFrame(rows)


def test_ideal_geometry_scores_higher_than_poor():
    ideal = er.pose_engagement([3.5], [180.0], ["sn2_displacement"])[0]
    poor = er.pose_engagement([4.15], [151.0], ["sn2_displacement"])[0]
    outside = er.pose_engagement([8.0], [90.0], ["sn2_displacement"])[0]
    assert ideal > poor > outside
    assert 0.0 <= outside and ideal <= 1.0


def test_unknown_statistic_raises_rather_than_defaulting():
    p = _poses({("m1", 0): [(3.5, 180.0)]})
    with pytest.raises(ValueError, match="unknown statistic"):
        er.mode_engagement(p, statistic="best_ever")


def test_unknown_ligand_aggregation_raises():
    m = pd.DataFrame([dict(ident="m1", mode=0, engagement=0.5)])
    with pytest.raises(ValueError, match="unknown aggregation"):
        er.rank_ligands(m, how="geometric_mean")


def test_spread_is_reported_so_a_mixture_is_visible():
    """A group whose poses disagree has no meaningful summary; the reader must
    be able to see that rather than infer it."""
    p = _poses({("m1", 0): [(3.5, 180.0), (8.0, 90.0)],      # a mixture
                ("m1", 1): [(3.5, 180.0), (3.5, 179.0)]})    # agreeing
    e = er.mode_engagement(p).set_index("mode")
    assert e.loc[0, "engagement_spread"] > 0.5
    assert e.loc[1, "engagement_spread"] < 0.05


def test_a_group_of_one_is_scored_not_dropped():
    """An engagement score is a pose property; a frequency is not. The pose-count
    gate must not be applied by default."""
    p = _poses({("m1", 0): [(3.5, 180.0)],
                ("m1", 1): [(4.1, 155.0)] * 40})
    r = er.rank_modes(er.mode_engagement(p), within=None)
    assert len(r) == 2
    assert r.iloc[0]["mode"] == 0, "the single well-anchored pose must rank first"


def test_min_poses_is_opt_in_and_needs_the_column():
    p = _poses({("m1", 0): [(3.5, 180.0)], ("m1", 1): [(3.6, 178.0)] * 20})
    assert len(er.rank_modes(er.mode_engagement(p), within=None, min_poses=12)) == 1
    with pytest.raises(ValueError, match="n_poses_mode is absent"):
        er.rank_modes(pd.DataFrame([dict(ident="m1", mode=0, engagement=0.4)]),
                      within=None, min_poses=12)


def test_ranking_is_ordered_best_first():
    p = _poses({("m1", 0): [(4.15, 151.0)] * 5,
                ("m1", 1): [(3.5, 180.0)] * 5,
                ("m1", 2): [(3.9, 165.0)] * 5})
    r = er.rank_modes(er.mode_engagement(p), within=None)
    assert list(r["mode"]) == [1, 2, 0]
    assert list(r.engagement_rank) == [1, 2, 3]


def test_ligand_mean_and_best_differ_and_both_carry_n_modes():
    m = pd.DataFrame([
        dict(ident="a", mode=0, engagement=0.9), dict(ident="a", mode=1, engagement=0.1),
        dict(ident="b", mode=0, engagement=0.5), dict(ident="b", mode=1, engagement=0.5)])
    mean = er.rank_ligands(m, how="mean")
    best = er.rank_ligands(m, how="best")
    assert set(mean.columns) >= {"n_modes", "ligand_aggregation"}
    assert mean.iloc[0].ligand_engagement == pytest.approx(0.5)
    assert best.iloc[0].ident == "a", "best must reward the one strong mode"
    assert mean.iloc[0].ident in ("a", "b")
    assert (mean.n_modes == 2).all()


def test_statistic_is_recorded_on_every_row():
    """Two statistics produce differently-ordered tables; a row that does not say
    which one made it is not comparable with one that used another."""
    p = _poses({("m1", 0): [(3.5, 180.0), (8.0, 90.0)] * 3})
    for stat in ("max", "mean", "q75_mean"):
        e = er.mode_engagement(p, statistic=stat)
        assert (e.engagement_statistic == stat).all()
    assert (er.mode_engagement(p, statistic="max").engagement.iloc[0] >
            er.mode_engagement(p, statistic="mean").engagement.iloc[0])


def test_unmapped_mechanism_raises_rather_than_scoring_zero():
    """A typo in the mechanism name used to return 0.0 for every pose, which
    ranks the whole molecule last with no error raised anywhere."""
    from shared import nac_criterion as nac
    with pytest.raises(ValueError, match="unknown mechanism"):
        nac.anchor_quality(3.5, 180.0, "sn2")
    with pytest.raises(ValueError, match="unknown mechanism"):
        er.pose_engagement([3.5], [180.0], ["michael"])


# --------------------------------------------------------------------------- #
#  the conservative support term
# --------------------------------------------------------------------------- #
def test_support_can_never_promote_a_pose_that_cannot_react():
    """THE SAFETY PROPERTY. `rank_score = anchor * support`, so a pose with no
    attack geometry scores zero however well supported it is. An ADDITIVE
    support term would let contacts rescue a pose that cannot react, which is
    exactly the influence the rule is meant not to have."""
    assert er.rank_score(0.0, 1.0 + er.MAX_SUPPORT) == 0.0
    assert er.rank_score(0.5, 1.0) == pytest.approx(0.5)
    assert er.rank_score(0.5, 1.15) == pytest.approx(0.575)


def test_support_is_bonus_only_and_bounded():
    """It may never demote a pose, and may never outweigh geometry."""
    s = er.receptor_support_atoms()
    far = er.support_factor(np.array([[999.0, 999.0, 999.0]]), ["O"], s)
    assert far == 1.0, "an unsupported pose must keep its geometry score exactly"
    near = er.support_factor(np.array([s[0] + np.array([3.0, 0, 0]),
                                       s[1] + np.array([3.0, 0, 0])]), ["O", "N"], s)
    assert near == pytest.approx(1.0 + er.MAX_SUPPORT)
    assert er.MAX_SUPPORT <= 0.25, (
        "a tie-breaker, not a weight: above ~0.25 a supported mediocre pose "
        "could outrank an unsupported excellent one")


def test_only_polar_ligand_atoms_earn_support():
    """Letting carbon satisfy an H-bond claim turns this into a proximity bonus,
    which is a shape preference by another name."""
    s = er.receptor_support_atoms()
    pos = np.array([s[0] + np.array([3.0, 0, 0])])
    assert er.support_factor(pos, ["C"], s) == 1.0
    assert er.support_factor(pos, ["O"], s) > 1.0


def test_support_atoms_are_the_two_serines_only():
    """The omissions are the design: His59 would reward pi-stacking, the leucines
    would reward greasiness, and the Arg loop would select for anions."""
    resi = {r for _, r, _, _ in er.SUPPORT_ATOMS}
    assert resi == {114, 115}
    assert 59 not in resi and 122 not in resi and 68 not in resi


def test_a_missing_support_atom_raises(tmp_path):
    """Silently losing one serine would lower the factor for every pose
    uniformly and read as a chemistry result."""
    p = tmp_path / "half.pdb"
    p.write_text("ATOM      1  OG  SER A 114      10.743   0.304  -2.297  1.00  0.00           O\n")
    with pytest.raises(ValueError, match="support atoms absent"):
        er.receptor_support_atoms(p)


def test_fraction_above_refuses_a_default_cutoff():
    m = pd.DataFrame([dict(ident="a", engagement=0.9)])
    with pytest.raises(ValueError, match="needs an explicit cutoff"):
        er.rank_ligands(m, how="fraction_above")


def test_fraction_above_and_best_order_differently():
    """The reason both exist. `best` sees one group; `fraction_above` sees how
    much of the molecule's repertoire is reaction-competent."""
    m = pd.DataFrame([
        dict(ident="a", engagement=0.9), dict(ident="a", engagement=0.8),
        dict(ident="a", engagement=0.1),
        dict(ident="b", engagement=0.95), dict(ident="b", engagement=0.02),
        dict(ident="b", engagement=0.01)])
    frac = er.rank_ligands(m, how="fraction_above", cutoff=0.5)
    best = er.rank_ligands(m, how="best")
    assert frac.iloc[0].ident == "a", "a has 2 of 3 competent, b has 1 of 3"
    assert best.iloc[0].ident == "b", "b holds the single strongest group"
    assert set(frac.columns) >= {"n_above", "n_modes"}
