---
id: D0116
title: The gate's ROC-AUC did not scale from matched decoys to measured inactives
date: 2026-09-09
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/enrichment_gate.py
  - tests/test_enrichment_gate_auc.py
evidence:
  - 'roc_auc was O(actives x decoys): 34 x 361,354 = 12.3M Python comparisons per call'
  - 'bootstrap_ci calls it n_boot=2000 times -> 24.6 billion comparisons'
  - 'bootstrap_ci also drew each index with random.Random.choice: 722M calls before any AUC'
  - 'the gate was designed for property-matched decoys at 50 per active (~1,700), a 200x smaller population'
  - 'first full run ground for ~15 min without completing one confidence interval; killed'
  - 'vectorised form: one AUC at gate scale in under 1 s, verified equal to the definition to 1e-12 including ties'
runbook: null
---

## Context

Phase 1 ingested two datasets of **assayed** negatives so that enrichment could
be measured against molecules someone actually tested rather than against
property-matched decoys. PubChem AID 504891 contributes **34 actives and
361,354 measured inactives**.

Running the pharmacophore scorer through `enrichment_gate.evaluate()` on that
population did not finish.

## What was wrong

`roc_auc` computed the Mann-Whitney statistic as the definition literally
states it — every active against every decoy, in a Python generator
expression:

```python
wins = sum((1.0 if p > q else 0.5 if p == q else 0.0) for p in pos for q in neg)
```

That is 34 × 361,354 = **12.3 million** interpreted comparisons per call.
`bootstrap_ci` calls it 2,000 times: **24.6 billion**. And before any of that,
the bootstrap drew each resampled index individually with
`random.Random.choice` — 722 million calls per gate.

Nothing was incorrect. Every number it would eventually have produced was the
right number. It simply does not finish.

## Why it looked right

**Because it was right when it was written, for the population it was written
for.** The gate's decoys are property-matched at 50 per active
(`config/gates.yaml: n_per_active: 50`) — on the covalent stratum's lead-tier
actives that is on the order of 1,700 decoys, where the pairwise form is
~115 million comparisons across the whole bootstrap and completes in about a
minute. The code carried no constant to go stale and no threshold to
mis-trip; the cost lived entirely in the *shape* of the input.

Then the input changed by a factor of 200, in a different part of the project,
for a good reason, and nothing connected the two. This is catalogue **#19**'s
shape — `timeout=86400`, sized when a pool was 1,882 molecules and fatal when
the pool became 16,806 — and it fails the same way: **by consuming the run
rather than by raising.** A run that cannot finish reports nothing at all,
which reads as "still going" for as long as anyone is willing to wait.

It is worth being exact about what makes this family dangerous: there is no
wrong number to notice. The catalogue's usual tell — a populated, plausible
value computed from the wrong thing — is absent. What you get instead is
silence, and silence is indistinguishable from patience.

## Decision

`roc_auc` sorts the decoys once and binary-searches each active into them:

```python
below = np.searchsorted(neg_sorted, pos, side="left")    # strictly below
upto  = np.searchsorted(neg_sorted, pos, side="right")   # below or tied
wins  = below.sum() + 0.5 * (upto - below).sum()
```

This counts **precisely the same pairs**: `"left"` is the number of decoys
strictly below an active, and the gap to `"right"` is the number tied, which
is where the definition's half weight goes. O(n log n) instead of
O(actives × decoys).

The pairwise implementation is **kept**, renamed `_roc_auc_pairwise`, as the
reference the fast path is checked against. It is not called in anger.

`bootstrap_ci` resamples with a vectorised numpy RNG. **This deliberately
changes the RNG stream**, and that is stated rather than hidden: a bootstrap CI
is a Monte Carlo estimate whose endpoints were never a fixed quantity, so what
must be preserved is the *estimator*, not the digits. The test requires
agreement with the old sampler to within Monte Carlo error on a population
small enough for the old one to finish.

## What was checked

`tests/test_enrichment_gate_auc.py`, 15 tests, following this project's own
best-scoring habit: **run the new thing on something whose answer is already
written down.** Here the written-down answer is the previous implementation.

- Exact equality with `_roc_auc_pairwise` to 1e-12, across both directions and
  three population shapes including the real 34-active geometry.
- **A dedicated tie test**, because tie handling is the part a searchsorted
  form gets wrong if it uses one insertion side. Scores are rounded to one
  decimal to force many exact ties — which is also the realistic case, since a
  Tanimoto over a 15-molecule model set takes few distinct values. The fixture
  asserts it actually produced ties rather than assuming it did.
- All-tied scores must read exactly 0.5; perfect and inverted separation must
  read 1.0 and 0.0; the two directions must sum to 1.
- A speed regression at the true gate size (34 × 361,354) asserting one AUC in
  under a second, so this cannot silently regress to the shape it just left.

## Consequence

The gate is now runnable against measured inactives, which is what Phase 2.1
was blocked on for reasons unrelated to the GPU. Any scorer — not only this
one — can now be graded against 361,354 assayed negatives in minutes.
