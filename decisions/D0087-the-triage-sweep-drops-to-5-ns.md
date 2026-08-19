---
id: D0087
title: The triage sweep drops to 5 ns, to part-fund the finer pose split
date: 2026-08-19
status: proposed
approach: shared
decided_by: '@tt8804'
origin: user
supersedes: []
superseded_by: null
amends: D0085
affects:
  - config/target.yaml
  - decisions/D0085-the-triage-sweep-is-8-ns-and-the-survivor-bar-is-0.35-nm.md
evidence:
  - '@tt8804: "drop to 5 ns md sweeps"'
  - "D0085 measured the optimum as a PLATEAU, bootstrap 95% CI 4.3-9.5 ns -- 5 ns is inside the interval that produced 8 ns"
  - "truncation is ONE-SIDED: max@5 <= max@8, so a shorter sweep cannot drop a genuine survivor from the set, only admit extras"
  - "total variation in pass rate across 8-10 ns was 3.7%, a step function moving 0.60% per crossing -- the surface is flat, not peaked"
  - "cost, from D0085's own measurement of 4.10 min/ns over 251 completed sweeps: 32.8 min/mode at 8 ns, 20.5 min at 5 ns -- a 37% cut"
  - "the finer pose split (exp/4_election, 0.1 nm cut, no cap) raises the worklist from ~5 to ~25 modes per molecule, so 147 modes becomes ~735"
  - "net GPU: 147 x 32.8 min = 80 GPU-h before; 735 x 20.5 min = 251 GPU-h after -- the shorter sweep softens a 5x increase to 3.1x, it does not absorb it"
  - "the survivor bar is UNCHANGED at 0.35 nm; only the observation window moves"
runbook: null
---

# D0087 — the triage sweep drops to 5 ns

## The decision

`md.sweep_ps: 8000 -> 5000`. The survivor bar stays at 0.35 nm max ligand RMSD.
This amends [D0085](D0085-the-triage-sweep-is-8-ns-and-the-survivor-bar-is-0.35-nm.md)
on length only; every other part of that record stands.

## Why it is safe

D0085 did not measure 8 ns as a peak. It measured a **plateau**: bootstrap 95% CI
4.3-9.5 ns, pass rate a step function moving 0.60% per crossing, total variation
across 8-10 ns of 3.7%. 5 ns is inside the interval that produced 8, so this is
a move within the measurement rather than against it.

More importantly the error is one-sided. Max ligand RMSD is monotone in
observation time -- `max@5 <= max@8` for the same trajectory -- so a shorter
sweep can only ever let MORE modes through the 0.35 nm bar, never fewer. A mode
that would have survived 8 ns survives 5 ns by construction. What a shorter
window buys is extras, and an extra costs one 100 ns run that then rejects it,
while a miss costs a candidate. D0085 chose that asymmetry deliberately
(`capture_target: 0.95`) and this preserves it.

## Why now

The finer pose split measured in `exp/4_election` (0.1 nm cut, no cap on mode
count, ranked on attack angle) elects the validated pose 27/30 against 16/30 for
the shipped rule -- but it produces about 25 modes per molecule where the current
splitting produces 5.

The arithmetic, using D0085's own 4.10 min/ns over 251 completed sweeps:

| | modes | min/mode | GPU-h |
|---|---|---|---|
| now (8 ns, ~5 modes/mol) | 147 | 32.8 | 80 |
| finer split at 8 ns | ~735 | 32.8 | 402 |
| finer split at 5 ns | ~735 | 20.5 | 251 |

So this **softens a 5x increase to 3.1x. It does not pay for it.** The rest has
to come from the selection rule -- the enrichment floor (#71, still unmeasured)
or a top-k per molecule -- and that is a separate decision, not something this
one quietly assumes.

## What changes in the run

Nothing but the number, and every reader of it derives from config: the sweep
worker, the reports, the rail legend and the explainer deck all take the length
from `md.sweep_ps`, so the GUI says "5 ns" without a second edit. That wiring was
put in on 2026-08-18 after every one of those strings had been reading "10 ns"
against an 8 ns run.

**Not comparable across the change.** A sweep at 5 ns and a sweep at 8 ns produce
different max-RMSD values for the same mode, so nac_v5's sweep table cannot be
pooled with anything measured under this setting. A re-run gets a new topic.
