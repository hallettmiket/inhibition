"""
Purpose: Pin D0034 -- writing one metric must not erase the other metrics in
         its stratum, and writing one stratum must not erase the others.
Author: Mike Hallett (with Claude Code)
Date: 2026-07-28
Input: synthetic GateResult objects + a temp token path
Output: pass/fail

This defect has now appeared twice at two granularities. The first fix made a
run stop erasing OTHER strata; it kept erasing other metrics inside its own
stratum. Both directions are tested here so the next change has to break a
named test rather than a silent invariant.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from shared import enrichment_gate as eg  # noqa: E402


def _result(stratum: str, metric: str, auc: float) -> eg.GateResult:
    return eg.GateResult(
        metric=metric, stratum=stratum, higher_is_better=False,
        n_actives=5, n_decoys=50, n_chemotypes=6,
        roc_auc=auc, roc_auc_ci=(auc - 0.1, auc + 0.1),
        ef_1pct=0.0, bedroc=0.0, verdict="WEAK", reasons=[])


def _metrics(path: Path, stratum: str) -> dict:
    return json.loads(path.read_text())["strata"][stratum]["metrics"]


def test_writing_one_metric_keeps_the_others_in_that_stratum(tmp_path):
    """The D0034 failure: mmgbsa_dG erased affinity_kcal from `covalent`."""
    tok = tmp_path / "gate.token"
    eg.write_token([_result("covalent", "affinity_kcal", 0.537)], tok)
    eg.write_token([_result("covalent", "mmgbsa_dG", 0.260)], tok)

    m = _metrics(tok, "covalent")
    assert set(m) == {"affinity_kcal", "mmgbsa_dG"}
    assert m["affinity_kcal"]["roc_auc"] == 0.537


def test_writing_one_stratum_keeps_the_other(tmp_path):
    """The earlier fix, still held."""
    tok = tmp_path / "gate.token"
    eg.write_token([_result("non_covalent", "vina_affinity", 0.535)], tok)
    eg.write_token([_result("covalent", "affinity_kcal", 0.537)], tok)

    strata = json.loads(tok.read_text())["strata"]
    assert set(strata) == {"non_covalent", "covalent"}
    assert _metrics(tok, "non_covalent")["vina_affinity"]["roc_auc"] == 0.535


def test_rewriting_the_same_metric_replaces_it(tmp_path):
    """Supersession must still work -- a corrected value must win."""
    tok = tmp_path / "gate.token"
    eg.write_token([_result("covalent", "mmgbsa_dG", 0.140)], tok)
    eg.write_token([_result("covalent", "mmgbsa_dG", 0.260)], tok)

    m = _metrics(tok, "covalent")
    assert len(m) == 1
    assert m["mmgbsa_dG"]["roc_auc"] == 0.260


def test_recommended_metric_is_the_best_surviving_one(tmp_path):
    """With both metrics present, docking (0.537) must beat MM-GBSA (0.260)."""
    tok = tmp_path / "gate.token"
    eg.write_token([_result("covalent", "affinity_kcal", 0.537)], tok)
    eg.write_token([_result("covalent", "mmgbsa_dG", 0.260)], tok)

    cov = json.loads(tok.read_text())["strata"]["covalent"]
    assert cov["recommended_rank_metric"] == "affinity_kcal"
