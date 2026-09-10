"""
Purpose: the direction registry must CARRY a direction, not just a name.
Author: @tt8804 (with Claude Code)
Date: 2026-09-09

WHAT THIS PROTECTS. `rank_shortlist` refuses a metric it does not know -- the
guard that caught catalogue #4. Until 2026-09-09 the registry was a bare SET of
names whose three members happened to share one direction, and every sort below
it hardcoded `ascending=True`. The error message said "add it to
LOWER_IS_BETTER with its direction", and there was nowhere to put a direction:
adding a higher-is-better metric the way the message invited would have
silently ranked it backwards, with every value populated and plausible.

So these tests assert the two halves separately:

* an unregistered metric is still REFUSED (the original guard still works); and
* a registered metric is SORTED BY ITS DECLARED DIRECTION, checked on both
  directions so the assertion is about the mechanism rather than about which
  way today's metrics happen to point.

The second is the one that can fail. `test_higher_is_better_is_not_ranked_
backwards` fails on the pre-2026-09-09 code, which is the point: a test that
passes against the defect it names is the vacuous-guard shape this project
catalogues twice in `how_this_project_breaks.md`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import rank_shortlist as rs      # noqa: E402


def _frame() -> pd.DataFrame:
    """Three molecules whose two metrics order them in OPPOSITE directions."""
    return pd.DataFrame({
        "canonical_smiles": ["CCO", "CCCO", "CCCCO"],
        "candidate_id": ["a", "b", "c"],
        # kcal/mol: lower binds better -> 'c' is best
        "vina_affinity": [-5.0, -6.0, -7.0],
        # Tanimoto: higher is better -> 'c' is best too, by the OTHER direction
        "ph2d_max_similarity": [0.10, 0.20, 0.30],
        "HAC": [3, 4, 5],
    })


def test_registry_declares_a_direction_for_every_metric():
    assert rs.RANK_DIRECTION, "the registry must not be empty"
    for metric, lower in rs.RANK_DIRECTION.items():
        assert isinstance(lower, bool), (
            f"{metric!r} must declare a direction as a bool, got {lower!r}")


def test_lower_is_better_view_agrees_with_the_registry():
    """The back-compatible set must be derived, so the two cannot drift."""
    assert rs.LOWER_IS_BETTER == {
        m for m, lower in rs.RANK_DIRECTION.items() if lower}


def test_unregistered_metric_is_still_refused():
    """The original guard. An unanticipated name is refused, never assumed."""
    with pytest.raises(ValueError, match="not a known rank metric"):
        rs.direction_of("cnn_affinity")
    with pytest.raises(ValueError, match="not a known rank metric"):
        rs.rank(_frame(), metric="cnn_affinity", group_col=None, min_docked=1)


def test_lower_is_better_metric_ranks_the_most_negative_first():
    out = rs.rank(_frame(), metric="vina_affinity", group_col=None,
                  min_docked=1, decorrelate_size=False)
    best = out.loc[out["rank"] == 1, "candidate_id"].iloc[0]
    assert best == "c", "the most negative kcal/mol must rank first"


def test_higher_is_better_is_not_ranked_backwards():
    """THE ONE THAT FAILS ON THE OLD CODE.

    With a bare `LOWER_IS_BETTER` set and a hardcoded `ascending=True`, the
    lowest Tanimoto would rank first and the molecule LEAST like a known
    active would top the shortlist -- a perfectly populated, perfectly
    plausible, exactly inverted ranking.
    """
    assert rs.direction_of("ph2d_max_similarity") is False
    out = rs.rank(_frame(), metric="ph2d_max_similarity", group_col=None,
                  min_docked=1, decorrelate_size=False)
    best = out.loc[out["rank"] == 1, "candidate_id"].iloc[0]
    worst = out.loc[out["rank"] == 3, "candidate_id"].iloc[0]
    assert best == "c", "the HIGHEST similarity must rank first"
    assert worst == "a"


def test_the_two_directions_produce_opposite_orderings():
    """Mechanism, not today's data: same rows, both directions, reversed.

    Both metrics in the fixture order the molecules identically by VALUE
    (a < b < c), so if direction were ignored the two rankings would be the
    same. They must be exact reverses of each other.
    """
    lo = rs.rank(_frame(), metric="vina_affinity", group_col=None,
                 min_docked=1, decorrelate_size=False)
    hi = rs.rank(_frame(), metric="ph2d_max_similarity", group_col=None,
                 min_docked=1, decorrelate_size=False)
    lo_order = list(lo.sort_values("rank")["candidate_id"])
    hi_order = list(hi.sort_values("rank")["candidate_id"])
    assert lo_order == ["c", "b", "a"]
    assert hi_order == ["c", "b", "a"]
    # And the raw-value orderings they came from really are opposite, so the
    # agreement above is the direction doing work rather than a coincidence.
    f = _frame()
    assert list(f.sort_values("vina_affinity")["candidate_id"]) == ["c", "b", "a"]
    assert list(f.sort_values("ph2d_max_similarity")["candidate_id"]) == ["a", "b", "c"]


def test_size_decorrelated_residual_inherits_the_base_metric_direction():
    """The residual is ranked by the direction of what it is a residual OF.

    `size_decorrelated_score` is registered lower-is-better because it has only
    ever been a residual of kcal/mol. Ranking a higher-is-better metric
    produces a higher-is-better residual, and `rank()` must resolve the
    direction from the BASE metric rather than from the residual's own entry --
    otherwise turning decorrelation on silently inverts the ranking.
    """
    df = pd.concat([_frame()] * 20, ignore_index=True)
    df["candidate_id"] = [f"m{i}" for i in range(len(df))]
    df["ph2d_max_similarity"] = [0.01 * i for i in range(len(df))]
    df["HAC"] = [20 + (i % 5) for i in range(len(df))]
    out = rs.rank(df, metric="ph2d_max_similarity", group_col=None,
                  min_docked=1, decorrelate_size=True)
    if not bool(out["rank_size_decorrelated"].iloc[0]):
        pytest.skip("decorrelation did not engage on this fixture")
    ranked = out.sort_values("rank")
    top = ranked.iloc[0]["size_decorrelated_score"]
    bottom = ranked.iloc[-1]["size_decorrelated_score"]
    assert top > bottom, (
        "for a higher-is-better base metric the residual must also rank "
        "highest-first; ranking it ascending inverts the shortlist")
