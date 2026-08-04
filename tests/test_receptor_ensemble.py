"""
Purpose: the receptor ensemble ranks on the median, and its paths are keyed on the receptor.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-04
Input: shared.noncovalent_dock_run
Output: pass/fail

WHY THIS EXISTS. #6 item 6 decided the ensemble and D0052 pre-registered how
four receptors become one number, BEFORE any ensemble result was looked at.
That pre-registration is only worth something if the code cannot drift off it,
so the rule is pinned here rather than left in prose.

Two independent things are guarded:

1. THE COMBINATION RULE IS THE MEDIAN, NOT BEST-ACROSS. Best-across is a
   maximum over four correlated draws, so its upward bias grows with the width
   of a ligand's score distribution -- and width scales with conformational
   flexibility. Our pools differ ~2x on that axis (liu_2024_c3 averages 10.65
   rotatable bonds against du_xu's 4.81), so ranking on best-across would hand
   the flexible pool an advantage unrelated to binding. That is the artefact
   class D0049 just removed, reintroduced one level up and invisibly, because
   "we docked into an ensemble" reads as a refinement.

   `test_best_across_rewards_spread_and_the_median_does_not` demonstrates the
   bias on constructed data rather than asserting it, so the reasoning survives
   somebody deciding best-across "looks better" on a future pool.

2. POSE DIRECTORIES ARE KEYED ON THE RECEPTOR. Four receptors writing into one
   `poses_ph7.4/` would overwrite each other and `collect_modes` would parse
   whichever finished last, with every manifest recording success. Same defect
   the ligand-prep cache carries a tag to prevent.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from shared import noncovalent_dock_run as ncd


def _scores(tag: str, values: dict[str, float]) -> pd.DataFrame:
    return pd.DataFrame({"candidate_id": list(values),
                         "vina_affinity": list(values.values())})


# --- the paths -------------------------------------------------------------

def test_each_receptor_gets_its_own_pose_directory():
    work = Path("/w")
    dirs = {r.tag: ncd.pose_dir(work, r) for r in ncd.ENSEMBLE}
    assert len(set(dirs.values())) == len(ncd.ENSEMBLE), (
        f"two receptors share a pose directory: {dirs}")


def test_the_default_receptor_is_tagged_too():
    """No special case for 6VAJ.

    Exempting the default would leave the common path keyed on less than its
    inputs -- the anti-pattern itself, with an exception bolted on.
    """
    assert ncd.DEFAULT_RECEPTOR.tag in ncd.pose_dir(Path("/w")).name


def test_every_receptor_carries_its_own_box():
    """A box is coordinates in ONE structure's frame.

    Pairing 6VAJ's box with 3IKG's receptor docks into empty space beside the
    site and returns affinities that look completely ordinary.
    """
    boxes = {r.tag: r.box for r in ncd.ENSEMBLE}
    assert len(set(boxes.values())) == len(ncd.ENSEMBLE), (
        f"two receptors share a box file: {boxes}")
    for r in ncd.ENSEMBLE:
        assert r.reference_ligand, f"{r.tag} has no reference ligand recorded"


# --- the combination rule --------------------------------------------------

def test_the_rank_metric_is_the_median():
    per = {"6VAJ": _scores("6VAJ", {"a": -8.0}),
           "3IKG": _scores("3IKG", {"a": -6.0}),
           "3IKD": _scores("3IKD", {"a": -5.0}),
           "9INR": _scores("9INR", {"a": -5.0})}
    row = ncd.combine_ensemble(per).iloc[0]
    assert row[ncd.ENSEMBLE_MEDIAN] == pytest.approx(-5.5)
    assert row[ncd.ENSEMBLE_BEST] == pytest.approx(-8.0)
    assert row[ncd.ENSEMBLE_ARGBEST] == "6VAJ"
    assert row[ncd.ENSEMBLE_N] == 4


def test_argbest_names_the_receptor_not_a_position():
    """Derived from the column name, never an index into a list.

    `cols[idx]` is correct only while column order is guaranteed, and selection
    by position where an identity was meant is the second of the four disguises
    in `how_this_project_breaks.md`.
    """
    per = {"6VAJ": _scores("6VAJ", {"a": -4.0}),
           "3IKG": _scores("3IKG", {"a": -4.5}),
           "3IKD": _scores("3IKD", {"a": -9.9}),   # the winner, in the middle
           "9INR": _scores("9INR", {"a": -4.1})}
    assert ncd.combine_ensemble(per).iloc[0][ncd.ENSEMBLE_ARGBEST] == "3IKD"


def test_best_across_rewards_SPREAD_and_the_median_does_not():
    """The measured reason D0052 rejected best-across.

    Two ligands with the SAME central tendency: one scored consistently, one
    with a wide spread across receptors -- the signature of a flexible molecule
    like liu_2024_c3's pool. Best-across prefers the wide one purely because a
    maximum over correlated draws rewards width. The median does not separate
    them, which is the point.
    """
    tight = {"6VAJ": -7.0, "3IKG": -7.0, "3IKD": -7.0, "9INR": -7.0}
    wide = {"6VAJ": -9.5, "3IKG": -4.5, "3IKD": -7.0, "9INR": -7.0}
    per = {tag: pd.DataFrame({"candidate_id": ["tight", "wide"],
                              "vina_affinity": [tight[tag], wide[tag]]})
           for tag in tight}

    out = ncd.combine_ensemble(per).set_index("candidate_id")

    assert out.loc["tight", ncd.ENSEMBLE_MEDIAN] == pytest.approx(
        out.loc["wide", ncd.ENSEMBLE_MEDIAN]), (
        "the median must not separate two ligands with the same centre")
    assert out.loc["wide", ncd.ENSEMBLE_BEST] < out.loc["tight", ncd.ENSEMBLE_BEST], (
        "best-across should favour the wide one -- if it does not, this test "
        "no longer demonstrates the bias it exists to document")


def test_a_partial_ensemble_is_counted_not_hidden():
    """A median over two receptors is not the same quantity as over four.

    It would otherwise look exactly like everyone else's. Nothing is dropped;
    the count is carried so a reader can refuse it.
    """
    per = {"6VAJ": _scores("6VAJ", {"a": -8.0, "b": -7.0}),
           "3IKG": _scores("3IKG", {"a": -6.0}),          # b missing
           "3IKD": _scores("3IKD", {"a": -5.0}),          # b missing
           "9INR": _scores("9INR", {"a": -5.0, "b": -7.0})}
    out = ncd.combine_ensemble(per).set_index("candidate_id")
    assert out.loc["a", ncd.ENSEMBLE_N] == 4
    assert out.loc["b", ncd.ENSEMBLE_N] == 2
    assert pd.notna(out.loc["b", ncd.ENSEMBLE_MEDIAN]), "b must not be dropped"


def test_duplicate_candidates_from_one_receptor_are_refused():
    """Two receptors' poses in one directory is exactly the collision fixed."""
    per = {"6VAJ": pd.DataFrame({"candidate_id": ["a", "a"],
                                 "vina_affinity": [-8.0, -6.0]})}
    with pytest.raises(RuntimeError, match="duplicate"):
        ncd.combine_ensemble(per)


def test_collect_modes_stamps_the_receptor_when_given_one(tmp_path):
    (tmp_path / "t2_aaa_out.pdbqt").write_text(
        "MODEL 1\nREMARK VINA RESULT:      -8.3      0.000      0.000\nENDMDL\n"
        "MODEL 2\nREMARK VINA RESULT:      -8.1      1.200      2.400\nENDMDL\n")
    tagged = ncd.collect_modes(tmp_path, ncd.ENSEMBLE[1])
    assert tagged.loc[0, "receptor"] == "3IKG"
    # Absent, not defaulted: stamping 6VAJ onto a legacy directory would
    # assert provenance nobody established.
    assert "receptor" not in ncd.collect_modes(tmp_path).columns


def test_the_ensemble_metric_is_not_the_gated_one():
    """D0052 decision 2: a median over four receptors is a NEW metric.

    D0016/D0041's verdict is a 6VAJ number measured with `box_expanded.json`.
    If the ensemble metric were ever named `vina_affinity` it would inherit
    that verdict silently -- the same shape as the defect D0051 fixed, where an
    unanticipated verdict defaulted to validating the ranking.
    """
    assert ncd.ENSEMBLE_MEDIAN != "vina_affinity"
    assert ncd.ENSEMBLE_MEDIAN.startswith("vina_affinity_")
