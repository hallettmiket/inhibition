---
id: D0041
title: The gate issued its first real verdict, and it is that docking does not demonstrably enrich
date: 2026-07-30
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/reference/pin1_reference_binders_3.csv
  - scripts/build_noncovalent_decoys_2.py
  - scripts/run_enrichment_gate.py
  - shared/reference_set.py
  - decisions/D0040-the-residual-auc-was-never-significant.md
evidence:
  - '[non_covalent/vina_affinity] WEAK — AUC 0.599 CI[0.311,0.873] EF1% 0.0 BEDROC 0.114'
  - '6 actives / 353 decoys / 6 independent chemotypes — the gate floor is 6, met exactly'
  - 'first non-UNDERPOWERED verdict in the project'
  - 'the CI includes 0.5, so enrichment is not demonstrated'
  - 'measured null at 6 actives among 82: 95% range [0.263, 0.741]; 0.599 sits inside it'
  - 'Liu-2024-C3 (PMID 39229909, Kd 130 nM, co-crystal 9INR) is the sixth chemotype'
  - 'decoys_non_covalent_1 topped out at MW 478 against Liu-2024-C3 at MW 547 — zero property matches'
  - '7 of 366 ligands failed to dock'
---

# A verdict at last, and it is negative

## What changed

Every enrichment gate this project has run returned UNDERPOWERED: fewer than
6 independent chemotypes, so no verdict could be claimed regardless of the
point estimates. D0040 quantified how bad that was — at 2 actives the null
spans 79% of the AUC scale, and every score we had measured fell inside it.

A systematic ChEMBL2288 sweep found **Liu-2024-C3** (Kd 130 nM by SPR,
co-crystal **9INR**), a non-covalent scaffold structurally distinct from the
five we had. Checked against the gate's own `cluster_chemotypes` rather than
asserted: 5 clusters become 6. Six is the floor.

## The verdict

```
[non_covalent/vina_affinity] WEAK
  ROC-AUC 0.599  CI[0.311, 0.873]  EF1% 0.0  BEDROC 0.114
  6 actives / 353 decoys / 6 chemotypes
  - ROC-AUC 0.599, but the 95% CI [0.311, 0.873] includes 0.5
```

**WEAK**, not UNDERPOWERED. The gate has enough independent chemotypes to
speak, and what it says is that Vina affinity does not demonstrably separate
known Pin1 binders from property-matched decoys. The point estimate is above
chance; the interval contains chance. EF1% is 0.0 — not one active in the top
1%.

The independently measured null agrees. At 6 actives among 82 the 95% range of
AUC under random labelling is [0.263, 0.741], and 0.599 sits comfortably inside
it. Two different calculations, same conclusion.

## What this does and does not license

Under the gate's own vocabulary only FAIL demotes `dock_score` to a displayed
label, so **nothing is demoted**. Docking continues to be used as it has been:
to generate and filter poses, not to rank candidates on affinity. That was
already the project's position, and it now rests on a measurement rather than
on caution.

What it does license is a statement we could not previously make. "We do not
know whether docking enriches" has become "docking does not demonstrably
enrich, measured at the gate's own power floor." Those are different claims and
only the second is worth reporting.

## Two silent failures on the way here, both nearly fatal to the result

**The decoy pool could not reach the new active.** Liu-2024-C3 has MW 547;
`decoys_non_covalent_1.csv` topped out at MW 478. Zero property matches, and
`filter_adequately_matched` would have excluded it exactly as it excludes EGCG
and the peptidic macrocycles — the sixth chemotype present in the reference
file and absent from the gate, with the verdict still reading UNDERPOWERED for
a reason no log line would have given. The cached ChEMBL pool already held 887
molecules within MW ±50 and logP ±1.5 of it; the earlier build had simply never
reached that far up the mass range.

**The mechanism string matched nothing.** Liu-2024-C3 was first recorded as
`non_covalent_active_site` — accurate English. The gate selects strata by exact
equality against `non_covalent`, so the row was dropped with no error, no
warning, and a log line reading "5 actives" identical to the previous run. It
was caught only because the number was one lower than expected.

Both are the failure mode this project keeps meeting: a value derived from a
name or a default, failing silently. `reference_set.load` now validates the
mechanism column against `VALID_MECHANISMS` and names the offending rows.

## What is still out of reach

The covalent side cannot be brought to power from published data.

- `sulfamate_acetamide` has **6** decoys available and cannot be deepened:
  broadened to "α-carbonyl CH₂ with any O-leaving group", a 20,858-molecule
  ChEMBL pool yields **1**. Reddi-4d and 4g are unrecoverable.
- `snar_chloroazine` has **3**, recoverable to 14 or 59 only by broadening the
  chemotype definition, which changes what the gate measures — a non-activated
  chloroazine is a weaker electrophile than the active it would stand in for.
- The covalent literature is ~5 chemotypes deep and converged on
  chloroacetamide. Targeted searches for boronic acid, SuFEx, vinyl sulfone,
  cyanoacrylamide, fumarate, epoxide, aldehyde and nitrile warheads against
  Pin1 returned nothing with a resolvable structure.

So T_3 and T_4 stay UNDERPOWERED, and no amount of compute changes that. The
remaining routes are unpublished fragment-screening data, or validating the
covalent protocol on a target that has hundreds of actives and applying it to
Pin1 as an application rather than as its own validation.

## The lesson

The floor did its job for weeks by refusing to certify anything, and the moment
it could speak it said the thing we did not want to hear. That is what a power
floor is for. The failure worth naming is not the negative result — it is that
the two bugs above would each have produced a *quieter* outcome, another
UNDERPOWERED that looked exactly like the last one, and neither would have been
noticed by anything except a count being one lower than expected.
