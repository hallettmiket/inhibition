---
id: D0017
title: Vina-GPU agrees on scores but is not adopted yet
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/build_vina_gpu.sh
  - scripts/validate_vina_gpu.py
evidence:
  - 'Spearman rho 0.869 and Pearson 0.890 against CPU Vina over 248 ligands'
  - 'mean absolute difference 0.229 kcal/mol, mean bias +0.188'
  - 'ROC-AUC 0.535 (CPU) -> 0.433 (GPU); drift 0.101 against a 0.10 threshold'
  - 'both CIs contain 0.5 and overlap over 0.238-0.624; EF1% is 0.0 for both'
  - '346.9 s on GPU vs ~20 min on 20 CPU cores — roughly 3.5x on 248 ligands'
runbook: null
---

## Context

AutoDock Vina is CPU-only, making it the throughput wall for T_1 (5-10k
candidates) and T_2. Vina-GPU 2.1 was built and validated against the 248
ligands M3 had already scored with CPU Vina, on the identical receptor and box.

Adoption required four checks: rank agreement, score agreement, that the
enrichment ROC-AUC reproduces within 0.10, and that the graded verdict matches.

## Decision

**Not adopted yet — but the failing check is not informative, and the score
agreement is good.**

Three of four checks pass. Rank correlation is 0.869 and scores differ by 0.229
kcal/mol on average with a small positive bias; both engines return the same
graded verdict (UNDERPOWERED) and the same EF1% of 0.0.

The failure is ROC-AUC drift: 0.101 against a 0.10 threshold, missing by 0.001.
That check is close to meaningless on this data. **The baseline AUC is 0.535 —
chance.** When there is no signal to preserve, ranking is decided by noise, so
small score changes move the AUC freely. The two confidence intervals overlap
across 0.238-0.624 and *both contain 0.5*: the engines are not disagreeing about
which molecules bind, they are agreeing that neither can tell.

Adopting on a knife-edge threshold would be as wrong as rejecting on one. The
honest position is that **score agreement is demonstrated and enrichment
agreement is untested**, because this data set cannot test it.

## Consequences

Vina-GPU stays built, wrapped and documented but out of the pipeline. CPU Vina
remains the T_1/T_2 engine, so nothing downstream changes and D0016's baseline
stands.

To settle it properly, one of:

- **Re-validate on a set where docking actually enriches.** An AUC-reproduction
  test is only meaningful against a non-null baseline. If the expanded actives
  set lifts non-covalent enrichment above chance, re-run this script.
- **Raise `--search_depth`.** The run used 10; deeper Monte Carlo search should
  converge closer to CPU Vina and shrink the drift. Cheap to try and it tests
  whether the difference is convergence rather than implementation.
- **Adopt on score agreement alone**, accepting that enrichment equivalence is
  unproven, if T_1's throughput becomes the binding constraint. A 3.5x speedup
  on 248 ligands should widen on 10,000, where kernel setup amortises.

The speedup is real but was measured on a small batch; it is not yet a
projection for T_1 scale.
