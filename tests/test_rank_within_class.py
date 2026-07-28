"""
Purpose: Tests for the shared ranking + quota shortlist used by all four approaches.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-27 (moved to shared.rank_shortlist 2026-07-28)
Input: synthetic candidate frames
Output: pytest pass/fail

The central claim is that ranking within class produces a different — and more
useful — shortlist than a global sort. `test_global_sort_would_starve_a_class`
states that as an executable assertion rather than a comment, using a frame
built so one class's raw affinities dominate.

The rest pin properties that are easy to invert by accident: which direction is
"better", that undocked rows stay null rather than sorting to the bottom, that
flags survive into the shortlist reason, and — since D0031 — that no ranking
claims to be validated when the gate says UNDERPOWERED.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import rank_shortlist as rs          # noqa: E402


def _frame() -> pd.DataFrame:
    """Three classes with deliberately unequal sizes and affinity ranges.

    Class C has the best raw affinities but only 4 successful docks; class B is
    small and mediocre; class A is large and middling. A global sort therefore
    excludes B entirely, which is the outcome the design exists to prevent.
    """
    rows = []
    for i in range(30):
        rows.append(dict(warhead_class="A", affinity_kcal=-9.0 + i * 0.1, HAC=30,
                         reactivity_flag="IN_WINDOW", candidate_id=f"a{i}"))
    for i in range(5):
        rows.append(dict(warhead_class="B", affinity_kcal=-7.0 + i * 0.1, HAC=20,
                         reactivity_flag="IN_WINDOW", candidate_id=f"b{i}"))
    for i in range(4):
        rows.append(dict(warhead_class="C", affinity_kcal=-11.0 + i * 0.1, HAC=44,
                         reactivity_flag="OUTSIDE_WINDOW", candidate_id=f"c{i}"))
    for i in range(3):          # enumerated, never docked
        rows.append(dict(warhead_class="A", affinity_kcal=None, HAC=30,
                         reactivity_flag="IN_WINDOW", candidate_id=f"x{i}"))
    return pd.DataFrame(rows)


@pytest.fixture
def ranked():
    return rs.rank(_frame(), metric="affinity_kcal", group_col="warhead_class",
                   min_docked=20)


def test_lower_affinity_ranks_first(ranked):
    """affinity_kcal is lower-is-better; rank 1 must be the minimum."""
    a = ranked[ranked.warhead_class == "A"].dropna(subset=["rank"])
    assert a.loc[a["rank"] == 1, "affinity_kcal"].iloc[0] == a.affinity_kcal.min()


def test_an_unknown_metric_is_refused_rather_than_assumed():
    """Direction must be declared, not guessed.

    Every metric here is lower-is-better. A higher-is-better one added silently
    would invert every ranking without erroring anywhere.
    """
    with pytest.raises(ValueError, match="direction"):
        rs.rank(_frame(), metric="cnn_affinity", group_col=None, min_docked=1)


def test_percentile_spans_the_class(ranked):
    assert ranked.loc[ranked.candidate_id == "a0", "group_percentile"].iloc[0] == 100.0
    assert ranked.loc[ranked.candidate_id == "a29", "group_percentile"].iloc[0] == 0.0


def test_undocked_rows_keep_null_rank(ranked):
    """"Did not dock" and "docked badly" are different facts."""
    x = ranked[ranked.candidate_id.str.startswith("x")]
    assert x["rank"].isna().all()
    assert x["group_percentile"].isna().all()


def test_small_classes_are_flagged_not_selective(ranked):
    assert bool(ranked.loc[ranked.candidate_id == "a0", "rank_is_selective"].iloc[0])
    assert not bool(ranked.loc[ranked.candidate_id == "b0", "rank_is_selective"].iloc[0])


def test_ligand_efficiency_is_affinity_per_heavy_atom(ranked):
    r = ranked.loc[ranked.candidate_id == "a0"].iloc[0]
    assert abs(r.ligand_efficiency - (-r.affinity_kcal / r.HAC)) < 1e-6


def test_every_class_gets_its_quota(ranked):
    short = rs.shortlist(ranked, quota=3)
    short = short[short.shortlist]
    assert len(short) == 9
    assert dict(short.warhead_class.value_counts()) == {"A": 3, "B": 3, "C": 3}


def test_global_sort_would_starve_a_class(ranked):
    """The design claim, as a test rather than a comment."""
    df = _frame().dropna(subset=["affinity_kcal"])
    global_top = set(df.nsmallest(9, "affinity_kcal").warhead_class)
    s = rs.shortlist(ranked, quota=3)
    within_top = set(s[s.shortlist].warhead_class)
    assert "B" not in global_top, "fixture no longer demonstrates the failure mode"
    assert within_top == {"A", "B", "C"}


def test_flags_travel_into_the_shortlist_reason(ranked):
    s = rs.shortlist(ranked, quota=3)
    short = s[s.shortlist]
    assert "OUTSIDE_WINDOW" in str(
        short[short.warhead_class == "C"].shortlist_reason.iloc[0])
    assert "not selective" in str(
        short[short.warhead_class == "B"].shortlist_reason.iloc[0])
    assert s.loc[~s.shortlist, "shortlist_reason"].isna().all()


def test_flagged_classes_can_be_excluded(ranked):
    """The caller decides which groups to drop; the module just honours it."""
    s = rs.shortlist(ranked, quota=3, exclude_groups={"C"})
    assert set(s[s.shortlist].warhead_class) == {"A", "B"}


def test_quota_larger_than_class_takes_the_whole_class(ranked):
    s = rs.shortlist(ranked, quota=50)
    short = s[s.shortlist]
    assert (short.warhead_class == "C").sum() == 4      # class C only has 4 docked
    assert short.affinity_kcal.notna().all()            # never an undocked row


def test_single_group_ranking():
    """T_1, T_2 and T_3 vary no warhead, so they rank as one group."""
    r = rs.rank(_frame(), metric="affinity_kcal", group_col=None, min_docked=20)
    assert set(r["rank_group"].unique()) == {"all"}
    docked = r.dropna(subset=["rank"])
    assert docked["rank"].min() == 1
    assert docked.loc[docked["rank"] == 1, "affinity_kcal"].iloc[0] == \
        docked.affinity_kcal.min()


def test_an_underpowered_gate_never_marks_a_ranking_validated(tmp_path):
    """D0031. No ranking in this project is validated, and that must be data.

    A rank shown without its verdict implies a confidence no gate here supports,
    so `rank_validated` is computed from the gate rather than asserted anywhere.
    """
    import json

    tok = tmp_path / "gate.token"
    tok.write_text(json.dumps({"strata": {"covalent": {"metrics": {
        "affinity_kcal": {"verdict": "UNDERPOWERED", "roc_auc": 0.537,
                          "roc_auc_ci": [0.346, 0.728], "ef_1pct": 0.0,
                          "n_chemotypes": 2}}}}}))
    orig, rs.GATE_TOKEN = rs.GATE_TOKEN, tok
    try:
        g = rs.attach_gate(_frame(), "covalent", "affinity_kcal")
    finally:
        rs.GATE_TOKEN = orig
    assert g["gate_verdict"].iloc[0] == "UNDERPOWERED"
    assert not g["rank_validated"].any()
    assert g["gate_roc_auc"].iloc[0] == 0.537


def test_a_missing_gate_is_reported_not_ignored(tmp_path):
    """An absent token must not read as 'fine'."""
    orig, rs.GATE_TOKEN = rs.GATE_TOKEN, tmp_path / "nope.token"
    try:
        g = rs.attach_gate(_frame(), "covalent", "affinity_kcal")
    finally:
        rs.GATE_TOKEN = orig
    assert g["gate_verdict"].iloc[0] == "UNGATED"
    assert not g["rank_validated"].any()


def test_one_stratum_run_does_not_erase_another(tmp_path):
    """Running the covalent gate must not delete the non-covalent verdict.

    It did: `write_token` rebuilt the payload from only the current run's
    results, so a covalent-only run silently removed the verdict T_1 and T_2
    read, and their ranking reported UNGATED as if no gate had ever run.
    """
    import json

    from shared import enrichment_gate as eg

    tok = tmp_path / "gate.token"
    tok.write_text(json.dumps({"strata": {"non_covalent": {"metrics": {
        "vina_affinity": {"verdict": "UNDERPOWERED", "roc_auc": 0.535}}}}}))

    r = eg.GateResult(metric="affinity_kcal", stratum="covalent",
                      higher_is_better=False, n_actives=2, n_decoys=81,
                      n_chemotypes=2, roc_auc=0.537, roc_auc_ci=(0.346, 0.728),
                      ef_1pct=0.0, bedroc=0.001, per_chemotype_auc={},
                      verdict="UNDERPOWERED", reasons=[])
    eg.write_token([r], token_path=tok)

    strata = json.loads(tok.read_text())["strata"]
    assert "covalent" in strata, "the run's own stratum must be written"
    assert "non_covalent" in strata, "the other stratum must survive the run"
    assert strata["non_covalent"]["metrics"]["vina_affinity"]["roc_auc"] == 0.535


def test_shared_identities_rank_as_one_molecule():
    """D0029, second act. A quota of 3 must mean three MOLECULES.

    T_4 carries one row per (R-group, warhead route), and the three SN2
    acetamides reach an identical adduct. Merging their classes removed the
    duplication BETWEEN groups; this pins that it is also gone WITHIN one.
    Ranking rows gave `acetamide_adduct` a top-3 of one molecule three times.
    """
    rows = []
    for i in range(5):                       # 5 molecules, each via 3 routes
        for route in ("chloroacetamide", "sulfamate_acetamide",
                      "sulfonate_acetamide"):
            rows.append(dict(adduct_class="acetamide_adduct", dock_id=f"d{i}",
                             warhead_class=route, affinity_kcal=-9.0 + i,
                             HAC=30, candidate_id=f"{route}_{i}"))
    df = pd.DataFrame(rows)

    r = rs.rank(df, metric="affinity_kcal", group_col="adduct_class",
                min_docked=1, identity_col="dock_id")
    # All three routes to one molecule share its rank.
    assert r[r.dock_id == "d0"]["rank"].nunique() == 1
    assert set(r["rank"].unique()) == {1, 2, 3, 4, 5}
    assert r["group_n_docked"].iloc[0] == 5, "group size counts molecules, not rows"

    s = rs.shortlist(r, quota=3)
    short = s[s.shortlist]
    assert short.dock_id.nunique() == 3, "a quota of 3 must be 3 distinct molecules"
    assert len(short) == 9, "every route to a shortlisted molecule is carried"


def test_a_negative_verdict_needs_the_same_power_as_a_positive_one():
    """The power floor must govern FAIL, not just STRONG.

    It did not. The FAIL branch sat above the floor, so a damning verdict could
    be returned from evidence the same function would refuse to call STRONG. The
    MM-GBSA gate hit this with ONE active: ROC-AUC 0.140, CI [0.040, 0.240],
    graded FAIL — "demonstrably anti-correlated with known actives".

    With one active, ROC-AUC is just that molecule's rank among the decoys and
    the bootstrap resamples the same active every time, so the interval looks
    tight precisely because between-active variation cannot enter it.
    """
    from shared import enrichment_gate as eg

    underpowered = eg.GateResult(
        metric="mmgbsa_dG", stratum="covalent", higher_is_better=False,
        n_actives=1, n_decoys=50, n_chemotypes=1,
        roc_auc=0.140, roc_auc_ci=(0.040, 0.240), ef_1pct=0.0, bedroc=0.0,
        per_chemotype_auc={})
    verdict, reasons = eg._verdict(underpowered, {})
    assert verdict == "UNDERPOWERED", (
        "a below-chance point estimate from ONE active must not be graded FAIL")
    assert any("1 active" in r for r in reasons)
    # The number is still reported — suppressing it would be the opposite error.
    assert any("BELOW chance" in r for r in reasons)


def test_fail_is_still_reachable_with_adequate_power():
    """The floor must not make FAIL unreachable — that would be its own defect."""
    from shared import enrichment_gate as eg

    powered = eg.GateResult(
        metric="affinity_kcal", stratum="covalent", higher_is_better=False,
        n_actives=12, n_decoys=500, n_chemotypes=7,
        roc_auc=0.21, roc_auc_ci=(0.11, 0.32), ef_1pct=0.0, bedroc=0.0,
        per_chemotype_auc={})
    verdict, _ = eg._verdict(powered, {})
    assert verdict == "FAIL"
