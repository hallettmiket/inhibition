"""
Purpose: Pin the measurement-error propagation the gate now owns (D0038).
Author: Mike Hallett (with Claude Code)
Date: 2026-07-29
Input: synthetic scored frames
Output: pass/fail

This quantity was reported in a decision record and a manuscript draft while
existing in no code anywhere in the repository -- computed ad hoc, unseeded,
untested, and never routed through the gate that grades every other metric. An
adversarial audit found it by grepping for an implementation and finding none.

These tests hold the three properties that make it safe to quote: it is
reproducible, it is kept separate from the bootstrap CI, and it widens as the
measurement gets worse rather than staying reassuringly narrow.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import enrichment_gate as eg  # noqa: E402


#: Canonical warhead class_ids, cycled over the synthetic actives.
#
# These fixtures declare `stratum="covalent"`, and since D0045 the covalent
# gate counts chemotypes by warhead class and REFUSES to fall back to
# structural clustering when the column is absent. That refusal is the point --
# a silent fallback would reintroduce the exact defect D0045 removed -- so the
# fixture supplies the column rather than the gate relaxing its guard.
#
# These tests are about measurement error, not chemotype counting, so the
# values only need to be real class_ids and distinct enough not to collapse the
# denominator by accident.
_WARHEADS = ["chloroacetamide", "sulfamate acetamide",
             "cinnamamide (aryl Michael acceptor; acrylamide-class)"]


def _frame(n_act=3, n_dec=30, sep=3.0, sem=0.2, seed=1) -> pd.DataFrame:
    """Actives separated from decoys by `sep`, each with error bar `sem`."""
    import numpy as np
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_act):
        rows.append({"canonical_smiles": f"C{'C'*i}O", "label": 1,
                     "warhead_class": _WARHEADS[i % len(_WARHEADS)],
                     "dG": -10.0 - sep + rng.normal(0, 0.5), "sem": sem})
    for i in range(n_dec):
        rows.append({"canonical_smiles": f"N{'C'*i}", "label": 0,
                     "warhead_class": None,
                     "dG": -10.0 + rng.normal(0, 0.5), "sem": sem})
    return pd.DataFrame(rows)


def test_is_reproducible_from_the_seed():
    df = _frame()
    a = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)
    b = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)
    assert a == b, "same seed must give the same interval"


def test_a_different_seed_moves_it_only_slightly():
    df = _frame()
    a = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False, seed=1)
    b = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False, seed=999)
    assert abs(a["mean"] - b["mean"]) < 0.05


def test_tiny_error_bars_barely_move_the_metric():
    """With precise measurements the interval should hug the point estimate."""
    df = _frame(sem=0.001)
    r = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)
    assert abs(r["mean"] - r["point"]) < 0.02
    assert r["ci"][1] - r["ci"][0] < 0.1


def test_large_error_bars_widen_it():
    """The whole reason to compute it: worse measurement, wider interval."""
    tight = eg.propagate_measurement_error(_frame(sem=0.1), metric="dG",
                                           sem_col="sem", higher_is_better=False)
    loose = eg.propagate_measurement_error(_frame(sem=8.0), metric="dG",
                                           sem_col="sem", higher_is_better=False)
    w_tight = tight["ci"][1] - tight["ci"][0]
    w_loose = loose["ci"][1] - loose["ci"][0]
    assert w_loose > w_tight * 3, f"{w_loose} should be much wider than {w_tight}"


def test_p_above_chance_is_low_for_a_below_chance_metric():
    """Actives placed WORSE than decoys must not read as enriching."""
    df = _frame(sep=-3.0)
    r = eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)
    assert r["point"] < 0.5
    assert r["p_above_chance"] < 0.05


def test_missing_columns_raise():
    df = _frame().drop(columns=["sem"])
    with pytest.raises(eg.EnrichmentGateError, match="sem"):
        eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)


def test_all_actives_or_all_decoys_raises():
    df = _frame()
    df["label"] = 1
    with pytest.raises(eg.EnrichmentGateError):
        eg.propagate_measurement_error(df, metric="dG", sem_col="sem",
                                       higher_is_better=False)


def test_evaluate_attaches_it_in_its_own_field():
    """Structural separation: two fields, because two different questions."""
    res = eg.evaluate(_frame(), metric="dG", stratum="covalent",
                      higher_is_better=False, sem_col="sem")
    assert res.measurement_error is not None
    assert "ci" in res.measurement_error
    assert res.measurement_error["holds_fixed"].startswith("the set of molecules")
    # roc_auc_ci must still be the BOOTSTRAP interval, untouched by propagation
    assert isinstance(res.roc_auc_ci, tuple)


def test_the_two_intervals_differ_when_the_result_is_marginal():
    """With clean separation both are [1,1]; the distinction only shows when the
    answer is close, which is exactly when someone would quote the narrower
    one. Two actives and heavy overlap is that case."""
    df = _frame(n_act=2, n_dec=40, sep=0.3, sem=2.0, seed=7)
    res = eg.evaluate(df, metric="dG", stratum="covalent",
                      higher_is_better=False, sem_col="sem")
    boot = tuple(res.roc_auc_ci)
    meas = tuple(res.measurement_error["ci"])
    assert boot != meas, f"bootstrap {boot} vs measurement {meas}"


def test_neither_interval_is_universally_the_wider_one():
    """A correction to how this was first described.

    On the real gate set the measurement interval was the NARROWER of the two
    ([0.356, 0.463] against a bootstrap [0.281, 0.500]), and it was tempting to
    generalise that into "the measurement interval always flatters the result".
    It does not: which is wider depends on the per-candidate error relative to
    the spread between candidates. Pinned so the claim stays specific to the
    data it came from.
    """
    tight = eg.evaluate(_frame(n_act=2, n_dec=40, sep=0.3, sem=0.05, seed=7),
                        metric="dG", stratum="covalent",
                        higher_is_better=False, sem_col="sem")
    loose = eg.evaluate(_frame(n_act=2, n_dec=40, sep=0.3, sem=4.0, seed=7),
                        metric="dG", stratum="covalent",
                        higher_is_better=False, sem_col="sem")

    def width(ci):
        return ci[1] - ci[0]

    tight_meas = width(tight.measurement_error["ci"])
    loose_meas = width(loose.measurement_error["ci"])
    # Precise measurement -> narrower than its bootstrap; poor measurement ->
    # wider. Both orderings occur, so neither may be assumed.
    assert tight_meas < width(tight.roc_auc_ci)
    assert loose_meas > width(loose.roc_auc_ci)


def test_evaluate_without_sem_col_leaves_it_none():
    """Absence must be visible, not silently defaulted to a narrow interval."""
    res = eg.evaluate(_frame(), metric="dG", stratum="covalent",
                      higher_is_better=False)
    assert res.measurement_error is None
