# Gates and controls

Four controls came out of an adversary audit of the spec. Each blocks a failure
that is **invisible as a crash** — the pipeline runs, output looks plausible,
and the result is wrong. They are declared in `config/gates.yaml`.

## 1. Docking-enrichment gate

Docking is not trusted to *rank* anything until it is shown to enrich known
actives over property-matched decoys **on this receptor**.

- Build property-matched decoys; spike in the verified covalent-Cys113 actives.
- Run the **exact downstream protocols** — a gate validating a different
  protocol than the one used downstream validates nothing.
- Require ROC-AUC ≥ 0.70, EF1% ≥ 5.0, BEDROC(α=20) ≥ 0.30, recovery ≥ 0.50.

!!! note "On failure the choreography does not stall"
    `dock_score` is **demoted to a displayed label** in every approach and the
    approaches lean on their other evidence plus the human. That is a config
    branch, not a code change.

**Status: not yet run.** No decoy set has been built.

## 2. Frozen external reference set

Novelty is computed against the published binder set, **never against the
seed**. Seed-relative novelty is circular: it measures edit distance from the
starting molecule, mechanically rewarding T_1 (no seed, so everything looks
novel) and penalizing T_2 for doing exactly what it was asked to do.

`novelty.py` has **no seed parameter**, by design.

## 3. Warhead-validity gate

In a prior real run, **6 of 16** warhead classes collapsed to an inert
formamide or sulfonamide once attached to the core. Invisible — docking
succeeds against an inert molecule and returns a plausible score.

This gate confirms each *attached* warhead is still a reactive electrophile of
its intended class, **before** covalent docking spends on it. If more than half
the classes die, it flags and stops for review.

## 4. Covalent-protocol parity

T_3 and T_4 import one pinned gnina covalent setup unchanged. If their recorded
protocol hashes differ, the within-covalent re-score in the GUI is **disabled**
rather than silently comparing incomparable numbers.

## Kill thresholds

Per-approach compute ceilings, in `config/gates.yaml`. The point is that each
approach stops rather than overrunning:

| Approach | Ceiling | Stop when |
|---|---|---|
| T_1 | 10,000 samples | valid+novel yield negligible, or gate FAIL |
| T_2 | CReM degree 1, frontier cap 200k | frontier exceeds cap |
| T_3 | 200 RL epochs, 500 final batch | — |
| T_4 | ≤36 MM-GBSA representatives | >half the warhead classes dead |
