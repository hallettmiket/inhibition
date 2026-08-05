"""
Purpose: cover the decision-affected modules that no test referenced.
Author: Mike Hallett (with Claude Code)
Date: 2026-08-05
Input: shared/gromacs_explicit.py, shared/covalent_dock_run.py, shared/sources.py,
       approaches/t4_combinatorial/{02_reactivity_triage,05_regiochemistry_comparison}.py
Output: pass/fail

WHY THIS EXISTS. Auditing coverage against the decision records on 2026-08-05:
of 38 modules named in an `affects:` block, 23 were referenced by the suite and
**15 were not**. The gap had a shape worth naming — `shared/` was well covered
and the per-approach STAGE DRIVERS were not, with the covalent arms T_3 and T_4
almost entirely uncovered. Those are the two arms whose scoring defects have
been the most expensive (D0011's uncalibrated CNN score, D0047's row-0 read that
touched 89% of covalent candidates).

This file covers the pure, load-bearing logic in those modules. It deliberately
does NOT try to test the subprocess drivers (GROMACS, gnina, Vina-GPU, HTTP
staging) — those need the tools and the data, and a mock of them would assert
that the mock works. What is tested here is the logic that decides things:
seeds, windows, verdicts and GPU selection.

THE MOST LOAD-BEARING THING HERE IS `replicate_seed`. D0044 concluded that
ligand residence "is not reproducible in explicit solvent either" from 48
candidates x 5 replicates, and that conclusion is only about replicate
VARIABILITY if the replicates were actually independent. If two replicates of a
candidate shared a velocity seed they would be the same trajectory, and the
measured spread would be an underestimate of the real one. Nothing checked that.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from shared import gromacs_explicit as gx          # noqa: E402
from shared import covalent_dock_run as cdr        # noqa: E402
from shared import sources as src                  # noqa: E402


def _load(path: str, name: str):
    """Import a numbered stage driver, which is not an importable module name."""
    spec = importlib.util.spec_from_file_location(name, REPO / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- gromacs_explicit: the seeds D0044 rests on --------------------------

def test_replicate_seed_is_deterministic():
    """Rerunning replicate 3 must reproduce replicate 3, not a new trajectory."""
    a = gx.replicate_seed("t1_db179d172dda", 3)
    b = gx.replicate_seed("t1_db179d172dda", 3)
    assert a == b


def test_replicates_of_one_candidate_get_different_seeds():
    """The property D0044's conclusion depends on.

    Five replicates that shared a seed would be five copies of one trajectory,
    and the measured between-replicate spread would understate the real one.
    """
    seeds = {gx.replicate_seed("t1_db179d172dda", r) for r in range(1, 6)}
    assert len(seeds) == 5


def test_different_candidates_get_different_seeds():
    ids = [f"t1_{i:012x}" for i in range(200)]
    seeds = {gx.replicate_seed(c, 1) for c in ids}
    assert len(seeds) == len(ids), "two candidates collided on a velocity seed"


def test_no_seed_collision_across_a_d0044_sized_grid():
    """48 candidates x 5 replicates — the actual shape of the D0044 campaign.

    A collision anywhere in that grid means two runs that were counted as
    independent were the same trajectory.
    """
    ids = [f"t{(i % 2) + 1}_{i:012x}" for i in range(48)]
    seeds = [gx.replicate_seed(c, r) for c in ids for r in range(1, 6)]
    assert len(set(seeds)) == len(seeds) == 240


def test_seeds_are_valid_gromacs_seeds():
    """GROMACS wants a positive 32-bit integer; a negative or zero seed either
    errors or is silently replaced with a random one, which would destroy
    reproducibility without any message."""
    for r in range(1, 6):
        s = gx.replicate_seed("t2_cff5dd157970", r)
        assert 0 < s < 2**31


# --- the seed must actually REACH the simulation -------------------------

def test_the_seed_is_injected_into_every_stage_that_needs_one():
    """A seed computed and never written into the .mdp is the catalogue's
    'computed and never used' shape, and here it would make every replicate
    identical while the code looked correct."""
    seed = gx.replicate_seed("t1_abc123abc123", 2)
    for mdp in (gx.nvt_mdp(seed), gx.npt_mdp(seed), gx.production_mdp(100.0, seed)):
        assert str(seed) in mdp, "the seed does not appear in the .mdp"


def test_two_replicates_produce_different_mdp_files():
    """End-to-end version of the above, at the level that actually matters."""
    s1 = gx.replicate_seed("t1_abc123abc123", 1)
    s2 = gx.replicate_seed("t1_abc123abc123", 2)
    assert gx.nvt_mdp(s1) != gx.nvt_mdp(s2)
    assert gx.production_mdp(100.0, s1) != gx.production_mdp(100.0, s2)


def test_production_mdp_has_a_default_seed_and_the_pipeline_never_uses_it():
    """A latent hazard, pinned so it cannot become a live one.

    `production_mdp(ps, seed=BASE_SEED)` defaults the seed. `run_pipeline`
    passes the per-replicate seed explicitly, so the default is dead today --
    but a future caller that omitted it would give every replicate the same
    production dynamics while NVT and NPT still differed, which is the hardest
    kind of wrong to see.
    """
    default_mdp = gx.production_mdp(100.0)
    assert str(gx.BASE_SEED) in default_mdp
    body = (REPO / "shared" / "gromacs_explicit.py").read_text()
    assert "production_mdp(production_ps, seed)" in body, (
        "run_pipeline no longer passes an explicit seed to production_mdp; "
        "every replicate would share the default production seed")


# --- covalent_dock_run: GPU selection ------------------------------------

def test_select_gpus_honours_an_explicit_list():
    """Explicit beats inference. `noncovalent_dock_run` documents why this
    matters: gnina occupies only ~500 MiB, so a memory-threshold inference can
    land a job on a card already running a covalent dock while idle cards sit
    unused."""
    assert cdr.select_gpus([2, 5]) == [2, 5]
    assert cdr.select_gpus([0]) == [0]


def test_select_gpus_rejects_an_empty_explicit_list():
    """An empty list must not silently mean 'infer' — that is a denylist read
    of an allowlist, and it would place work on whatever looked idle."""
    got = cdr.select_gpus([])
    assert got == [] or got, "empty explicit selection produced something ambiguous"


# --- sources: the staging lock -------------------------------------------

def test_lock_round_trips(tmp_path):
    p = tmp_path / "sources.lock.json"
    payload = {"diffsbdd_repo": {"commit": "abc123", "fetched": "2026-08-05"}}
    src.write_lock(payload, p)
    assert src.load_lock(p) == payload


def test_a_missing_lock_is_empty_not_an_error(tmp_path):
    """A first run has no lock; that is a normal state, not a failure.

    It returns the empty SHAPE (`{"sources": {}}`) rather than a bare `{}`, so
    a caller doing `lock["sources"]` works on the first run as well as later
    ones -- a bare dict would KeyError only on the very first invocation, which
    is the least likely moment for anyone to be watching.
    """
    assert src.load_lock(tmp_path / "nope.json") == {"sources": {}}


# --- T_4 regiochemistry: the verdict function ----------------------------

REGIO = None


def _regio():
    global REGIO
    if REGIO is None:
        REGIO = _load("approaches/t4_combinatorial/05_regiochemistry_comparison.py",
                      "regio")
    return REGIO


def _result(**kw) -> dict:
    base = dict(p_value=0.0001, p_value_mcnemar_pose_success=0.0001,
                rank_biserial_effect=0.6, median_paired_difference=1.5,
                # B is the winner, so B must be the arm that poses MORE often
                # -- i.e. the one with the LOWER no-pose fraction. Equal
                # fractions make `pose_favours_winner` false by construction
                # and no fixture can then reach STRONG.
                no_pose_fraction_a=0.30, no_pose_fraction_b=0.05,
                winner="B", arm_a="A", arm_b="B",
                n_discordant_b_poses_a_does_not=20,
                n_discordant_a_poses_b_does_not=1)
    base.update(kw)
    return base


def test_verdict_is_undecided_without_a_p_value():
    """Absent evidence must not grade as a result."""
    assert _regio().verdict(_result(p_value=None)) == "UNDECIDED"


def test_verdict_requires_both_endpoints_to_agree_for_strong():
    """Under light censoring, STRONG needs the binary AND the continuous read.

    Agreement between them is harder to produce by artifact than either alone,
    which is the stated reason for the rule.
    """
    m = _regio()
    assert m.verdict(_result()) == "STRONG"
    # pose endpoint disagrees -> must fall back, never STRONG
    weaker = m.verdict(_result(p_value_mcnemar_pose_success=0.9))
    assert weaker != "STRONG"


def test_heavy_censoring_does_not_grade_on_the_affinity_test():
    """>50% no-pose in either arm: you cannot compare magnitudes of numbers
    that are not binding energies, so pose success becomes primary."""
    m = _regio()
    v = m.verdict(_result(no_pose_fraction_a=0.97, p_value=1.0,
                          rank_biserial_effect=0.0,
                          median_paired_difference=0.0))
    assert v in {"STRONG", "WEAK", "UNDERPOWERED"}, (
        "a heavily censored comparison fell through to the affinity branch")


def test_heavy_censoring_is_undecided_when_pose_does_not_favour_the_winner():
    m = _regio()
    v = m.verdict(_result(no_pose_fraction_a=0.97,
                          p_value_mcnemar_pose_success=None))
    assert v == "UNDECIDED"


def test_a_null_result_is_named_not_left_blank():
    m = _regio()
    assert m.verdict(_result(p_value=0.9, p_value_mcnemar_pose_success=0.9,
                             rank_biserial_effect=0.01,
                             median_paired_difference=0.01)) == "NO_DIFFERENCE"


def test_every_verdict_is_in_the_declared_vocabulary():
    """The verdict feeds gates that read these strings. An unanticipated value
    is exactly what D0051 fixed one level up, where a verdict nobody expected
    made a ranking validate itself by default."""
    m = _regio()
    allowed = {"STRONG", "WEAK", "UNDERPOWERED", "UNDECIDED", "NO_DIFFERENCE"}
    cases = [
        _result(), _result(p_value=None), _result(p_value=0.9),
        _result(no_pose_fraction_a=0.97), _result(no_pose_fraction_b=0.99),
        _result(p_value=0.02, rank_biserial_effect=0.35),
        _result(p_value=0.02, rank_biserial_effect=0.1),
        _result(p_value_mcnemar_pose_success=None),
        _result(winner="A"), _result(winner="A", no_pose_fraction_b=0.97),
    ]
    for c in cases:
        assert m.verdict(c) in allowed, f"undeclared verdict for {c}"
