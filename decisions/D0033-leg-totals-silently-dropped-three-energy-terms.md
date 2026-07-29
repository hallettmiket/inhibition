---
id: D0033
title: Every dG in the project was summed from a partial energy, and a plausible number hid it
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/mmgbsa.py
  - scripts/recompute_mmgbsa_totals.py
  - tests/test_mmgbsa_energy_terms.py
  - decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md
evidence:
  - "ENERGY_TERMS asked for '1-4VDW'/'1-4EEL'; sander prints '1-4 VDW'/'1-4 EEL' with a space"
  - "the token regex [A-Z0-9\\-]+ stopped at the space and stored 1-4 VDW's value under the key 'VDW', which nothing read"
  - "1-4 EEL collided with the already-set 'EEL' key under setdefault and was discarded"
  - 'CMAP was never in ENERGY_TERMS at all'
  - 'net effect: three terms contributed exactly 0.0 to every leg total in the project'
  - 'recomputing with the old logic reproduces every stored dG exactly, confirming this produced all of them'
  - 'gate set: shift +17.00 kcal/mol mean, range +7.47 to +38.00, sd 6.63 -- large and non-constant'
  - 'shift differs BY APPROACH: t3 -8.77 mean vs gate +17.00, a ~26 kcal/mol systematic gap'
  - 'T_4 ranking near-inverts: Spearman(original, corrected) = -0.735, 0 of 5 top candidates retained'
  - 'D0032 re-run: MM-GBSA ROC-AUC 0.140 -> 0.260, Sulfopin 44/51 -> 38/51, still below docking 0.440, verdict still UNDERPOWERED'
  - '134 candidates recomputed by parsing alone; no minimisation was rerun'
---

# Leg totals were summed from a partial energy

## The defect

`LegEnergies.total` summed a fixed tuple of term names out of a dict parsed
from sander's `FINAL RESULTS` block:

```python
sum(self.terms.get(k, 0.0) for k in ENERGY_TERMS)
```

Three terms never arrived in that dict.

sander prints `1-4 VDW =  448.2115` and `1-4 EEL = 4403.4840` **with a space**.
The parser's token pattern `[A-Z0-9\-]+` stopped at the space, so it captured
`VDW` and `EEL` rather than the full labels. `VDW` was a key nothing looked up.
`EEL` had already been set by the genuine `EEL` term earlier in the block, and
`setdefault` kept the first value, so the 1-4 electrostatic term was silently
discarded. `CMAP` was parsed correctly and then excluded, because it was never
listed in `ENERGY_TERMS`.

So every leg total in this project omitted 1-4 VDW, 1-4 EEL and CMAP.

## Why it survived

Because `.get(k, 0.0)` cannot tell a term that is genuinely zero from a term
that failed to parse. Both are worth nothing to the sum, and neither says so.

And because the answers still looked right. The omitted terms are enormous --
1-4 EEL alone is ~4400 kcal/mol per leg -- but they cancel almost entirely
between the complex and receptor legs. What survives the cancellation is a
residue of +7 to +38 kcal/mol on the gate set, which lands a dG in exactly the
range a covalent adduct is supposed to occupy. **A plausible -15 kcal/mol is
not evidence of a correct -15 kcal/mol**, and nothing in the pipeline was
checking the number against anything but expectation.

The check that would have caught it on the first run was available the whole
time: sander prints its own total on the `NSTEP` line. Code that re-derives a
quantity the input already states should compare the two.

## Impact

**Every dG the project has reported is wrong**, across D1_8, D2_8, D3_11,
D4_20 and the D0032 gate. The correction is not a constant offset:

| set | n | mean shift | range |
|---|---|---|---|
| t1 | 23 | +5.21 | -5.18 .. +19.77 |
| t2 | 24 | +0.70 | -4.21 .. +9.45 |
| t3 | 25 | **-8.77** | -21.77 .. -1.26 |
| t4 | 11 | +14.83 | -14.15 .. +28.00 |
| gate | 51 | **+17.00** | +7.47 .. +38.00 |

Two consequences follow, and the second is worse than the first.

**Within an approach, rankings moved.** Spearman between old and corrected
ordering: t2 0.922, t3 0.870, t1 0.739, gate 0.685, and **T_4 -0.735** -- a
near-inversion in which none of the previous top five remain. Any shortlist
selected on MM-GBSA has to be re-derived.

**Across approaches, the sets are not comparable at all.** t3 shifted -8.77
while the gate shifted +17.00, a systematic gap of ~26 kcal/mol between two
groups of numbers that were being read on the same scale. D0020 already warned
that dG is comparable only within a warhead class; this was a second,
undocumented incomparability sitting underneath that one.

## What it does to D0032

D0032's **numbers are superseded, its conclusion is not.** Re-running the gate
on the same 51 ligands with corrected totals:

| metric | ROC-AUC (old) | ROC-AUC (corrected) | Sulfopin rank |
|---|---|---|---|
| docking `affinity_kcal` | 0.440 | 0.440 (unaffected) | 29/51 |
| MM-GBSA `dG_kcal` | 0.140 | **0.260** | 44/51 -> **38/51** |

MM-GBSA improves but remains below chance and well below docking, and the
verdict is still UNDERPOWERED on one active. D0032's finding -- that MM-GBSA
does not rescue the ranking, and that a negative verdict needs as much power as
a positive one -- stands. Its reported figures do not, and are replaced by
these.

That the conclusion survived is luck, not vindication. A +17 kcal/mol
non-constant error on the scoring function under test could as easily have
manufactured an enrichment as destroyed one.

## The fix

1. `ENERGY_TERMS` now lists the labels sander actually prints, including
   `1-4 VDW`, `1-4 EEL` and `CMAP`. `RESTRAINT` stays excluded deliberately --
   it is 0.0 here and is not part of the physical energy.
2. `ENERGY_LINE_RX` matches the multi-word labels **before** falling back to a
   bare token, so `1-4 VDW` can never again be read as `VDW`.
3. `LegEnergies.total` **raises** on a missing term instead of treating it as
   zero. `CMAP` is exempt because a ligand-only leg has no protein backbone and
   legitimately lacks it.
4. `parse_energy_block()` is shared, so the ensemble rescorer and the
   single-structure scorer cannot drift apart over one output format.
5. `tests/test_mmgbsa_energy_terms.py` pins all of it, including a test that
   the summed terms equal sander's own stated total.

## Correcting the record cost nothing

The defect was in **reading** sander's output, not in producing it. Every term
needed was already in the `.min.out` files the original runs wrote, so all 134
candidates were corrected by one pass over text with no minimisation rerun.
Corrected values are written to `dG_corrected.jsonl` per approach and a
combined `dG_corrected_index.jsonl`; the original `result.json` files are left
untouched, because they are the evidence for what D0032 reported at the time.

One candidate (`t4_9b8c1d7fb439`) has no `FINAL RESULTS` block at all -- a
minimisation that failed and was previously being summed to a number anyway.
The stricter parser now refuses it.

## The general lesson

This is the fifth time this project has been bitten by the same shape of bug
(D0025, D0028, D0029, D0030, now this): **a value was derived from a name
rather than from the thing itself, and the mismatch was silent.** Previously it
was chemotypes identified by a reactive-atom SMARTS or a mechanism label; here
it is energy terms identified by a string that did not match what the tool
prints.

The rule that follows: `dict.get(key, default)` in a sum over a fixed key list
is a silent-failure generator. If the key list is a contract, violating it must
raise.
