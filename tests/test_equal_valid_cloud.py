"""The analysed pose cloud is a fixed number of VALID poses, on one denominator.

THE DEFECT (D0106). `consensus` = mode_size / n_poses with `n_poses` = the
number DOCKED, equal at 500 for every molecule. But `labels` is -1 for every
PoseBusters-invalid pose, so only valid poses can enter a mode -- and the valid
fraction ranges 0.812 to 0.982 across nac_v6's acrylamide/bdhi set. The
denominator was equal; the numerator's CEILING was not. A molecule at 81% valid
could not reach a consensus above 0.81 however tight its poses were.

WHY NO EXISTING TEST CAUGHT IT. Because the obvious invariant -- "every row's
n_poses is the same" -- was TRUE, and had been checked across all 34,059 rows.
The tests below assert the invariant that actually matters instead: the
denominator equals the size of the set a mode can be drawn FROM.

WHAT WOULD MAKE THESE PASS WHEN THEY SHOULD FAIL. If `consensus` were recomputed
from `n_poses` rather than read from the frame, the first test would be checking
arithmetic it had just done. So it compares `consensus` against `n_poses_mode`
and `n_poses_kept` -- three columns the screen writes independently -- and
requires the mode populations to SUM to the kept count, which no single-row
recomputation can satisfy by accident.
"""

from __future__ import annotations

import glob
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(REPO))

from shared import run_paths as rp                        # noqa: E402
from shared import target_config as tc                    # noqa: E402


def _topic() -> str:
    """The topic under test.

    Overridable so the guard can be pointed at a PRE-D0106 topic and shown to
    fail -- a test that has only ever been run against data satisfying it has
    not been demonstrated to fail at all (`how_this_project_breaks`, the two
    vacuous guards).
    """
    import os
    return os.environ.get("INHIBITION_TEST_TOPIC") or rp.topic()


def _aggs() -> pd.DataFrame:
    fs = sorted(glob.glob(str(rp.BLACKSMITH / _topic() / "agg_s*.csv")))
    if not fs:
        pytest.skip(f"topic {_topic()!r} has no aggregate tables yet")
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)


def test_the_config_target_is_on_the_valid_count_not_the_attempt_count():
    """`n_runs` must exceed the target, or the target is unreachable.

    Not a style check: `n_runs` sized BELOW `target_pb_valid` would make
    `target_met` False for every molecule in the run, and the warning that says
    so would fire 562 times and be read as noise.
    """
    target = int(tc.get("docking.target_pb_valid", default=0) or 0)
    if not target:
        pytest.skip("target_pb_valid is 0 -- the cap is disabled by design")
    n_runs = int(tc.get("docking.n_runs"))
    assert n_runs > target, (
        f"docking.n_runs ({n_runs}) must exceed docking.target_pb_valid "
        f"({target}); no molecule passes PoseBusters at 100%, so a run sized at "
        f"or below the target cannot reach it for anything.")
    # The worst pass rate MEASURED on this library is 0.812 (D0106).
    assert n_runs * 0.812 >= target, (
        f"n_runs {n_runs} x the worst measured pass rate 0.812 = "
        f"{n_runs * 0.812:.0f}, short of the {target} target. Molecules in the "
        f"tail will be flagged `target_met = False` and will not be comparable.")


def test_every_mode_shares_one_denominator_and_it_is_the_analysed_cloud():
    """`n_poses` must be the KEPT count, and mode populations must sum to it."""
    d = _aggs()
    for col in ("n_poses", "n_poses_mode", "n_poses_kept", "n_poses_pb_valid",
                "n_poses_docked"):
        assert col in d.columns, f"the screen wrote no {col!r} column"

    ok = d[d.status == "ok"] if "status" in d.columns else d
    if ok.empty:
        pytest.skip("no molecules completed yet")

    # n_poses IS the kept count, not the docked count. This is the assertion
    # the pre-D0106 frames would fail: there n_poses == n_poses_docked.
    bad = ok[ok.n_poses != ok.n_poses_kept]
    assert bad.empty, (
        f"{len(bad)} row(s) have n_poses != n_poses_kept, so the denominator is "
        f"not the analysed cloud. e.g. {bad.iloc[0].ident}: n_poses="
        f"{bad.iloc[0].n_poses}, kept={bad.iloc[0].n_poses_kept}, "
        f"docked={bad.iloc[0].n_poses_docked} (D0106)")

    # Mode populations partition the analysed cloud exactly -- three
    # independently written columns agreeing, which recomputation cannot fake.
    for parent, g in ok.groupby("parent_ident"):
        kept = int(g.n_poses_kept.iloc[0])
        assert int(g.n_poses_mode.sum()) == kept, (
            f"{parent}: mode populations sum to {int(g.n_poses_mode.sum())}, "
            f"not the {kept} analysed poses -- some pose is in two modes or none")
        assert abs(float(g.consensus.sum()) - 1.0) < 1e-6, (
            f"{parent}: consensus sums to {g.consensus.sum():.6f}, not 1.0; "
            f"the shares are not fractions of one cloud")


def test_kept_never_exceeds_valid_and_valid_never_exceeds_docked():
    """The three counts must nest. A kept pose that is not valid is the bug."""
    d = _aggs()
    ok = d[d.status == "ok"] if "status" in d.columns else d
    if ok.empty:
        pytest.skip("no molecules completed yet")
    assert (ok.n_poses_kept <= ok.n_poses_pb_valid).all(), (
        "some molecule kept more poses than PoseBusters passed -- the cap is "
        "selecting from the wrong array")
    assert (ok.n_poses_pb_valid <= ok.n_poses_docked).all()


def test_molecules_short_of_the_target_are_flagged_not_hidden():
    """A molecule that cannot reach the target must say so on its own rows.

    It is not an error -- some molecule will eventually fall below the tail we
    sized for. It is only an error for it to be INDISTINGUISHABLE from one that
    met the target, because their consensus values are then on different
    denominators with nothing saying which.
    """
    d = _aggs()
    ok = d[d.status == "ok"] if "status" in d.columns else d
    if ok.empty:
        pytest.skip("no molecules completed yet")
    assert "target_met" in ok.columns, "no `target_met` column to check"
    target = int(tc.get("docking.target_pb_valid", default=0) or 0)
    if not target:
        pytest.skip("cap disabled")
    short = ok[ok.n_poses_kept < target]
    assert (~short.target_met.astype(bool)).all() if len(short) else True, (
        f"{len(short[short.target_met.astype(bool)])} row(s) are short of the "
        f"{target} target but carry target_met = True")
    met = ok[ok.n_poses_kept >= target]
    assert met.target_met.astype(bool).all() if len(met) else True
