---
id: D0065
title: Covalent candidates are ranked on warhead-presentation geometry, validated on two of three mechanisms
date: 2026-08-05
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/ranking_rationale.md
  - shared/nac_criterion.py
  - scripts/nac_screen.py
  - scripts/nac_rank.py
evidence:
  - '3IKD, 200 runs/molecule, 105 molecules, 0 unmeasurable. Positives = ligands crystallographically bonded to Cys113; negatives = warhead-matched compounds measured inactive in AID 504891 (qHTS Assay to Find Inhibitors of Pin1), shuffled.'
  - 'chloroacetamide (SN2): 9 pos vs 30 neg, enrichment 2.39x vs 0.82x, AUC 0.822, p=0.0020'
  - 'naphthoquinone/fumarate (Michael): 4 pos vs 30 neg, enrichment 2.39x vs 1.01x, AUC 0.800, p=0.0288'
  - 'snar_chloroazine (SNAr): 2 pos vs 30 neg, enrichment 2.66x vs 2.39x, AUC 0.558, p=0.41'
  - 'pooled on enrichment: 15 pos vs 90 neg, AUC 0.722, p=0.0031'
  - 'replicates across independent docking seeds: chloroacetamide AUC 0.872, 0.881, 0.852, 0.822 (mean 0.857)'
  - 'pooling RAW viable fractions instead of enrichment gives AUC 0.514 -- an artefact of differing window solid angles, not a result'
---

# Rank covalent candidates on whether they can present the warhead

## The decision

Covalent candidates are ranked by **enrichment**: the fraction of independent
docking runs in which the molecule reaches a mechanism-appropriate near-attack
conformation at Cys113, divided by the fraction an isotropically-approaching
nucleophile would reach by chance.

This replaces the docking score, MM-GBSA and MD residence for the covalent arms.
It is the first quantity measured on this target that separates actives from
inactives at all: five previous levels of theory did not (D0041, D0046, D0036,
D0038/D0044, D0057, D0061).

## What was measured

| class | mechanism | pos | neg | pos enrich | neg enrich | AUC | p |
|---|---|---|---|---|---|---|---|
| chloroacetamide | SN2 | 9 | 30 | 2.39× | 0.82× | **0.822** | **0.0020** |
| naphthoquinone/fumarate | Michael | 4 | 30 | 2.39× | 1.01× | **0.800** | **0.0288** |
| snar_chloroazine | SNAr | 2 | 30 | 2.66× | 2.39× | 0.558 | 0.41 |
| **pooled on enrichment** | — | 15 | 90 | 2.39× | 1.14× | **0.722** | **0.0031** |

Two mechanisms separate independently. It replicates across docking seeds
(chloroacetamide AUC 0.872 / 0.881 / 0.852 / 0.822).

## Why enrichment and not the raw viable fraction

**Because pooling the raw fractions gives AUC 0.514, and that number is an
artefact.** The SN2 window (≥150° from the leaving group) admits 6.7% of approach
directions; the perpendicular window admits 8.2%. Before the Bürgi–Dunitz
constraint was added they differed by 4.4×, and warhead-matched inactives scored
3–19% under SN2 against 71–81% under Michael — a statement about the windows, not
the molecules.

Dividing by each mechanism's own isotropic baseline removes it. This is D0020's
"rank within warhead class, not globally" reached from a new direction, and it is
load-bearing rather than cosmetic: T_3 is entirely acrylamide while T_4 spans
eight classes, so ranking them together on raw fractions would rank the
mechanisms.

## The controls, and why each is what it is

- **Positives are crystallographic, not annotated.** 17 ligands with an observed
  covalent bond to Cys113 at 1.6–1.9 Å, verified against `_struct_conn` (the set
  curated in issue #12 §A). "This molecule reacts with Cys113" is a structural
  fact here, not a label.
- **The HTS's own 34 actives were rejected as positives.** Read as chemistry
  rather than as labels, the 11 warhead-bearing ones are frequent hitters — two
  rhodanines, an azlactone, an arylidene barbiturate, a furfurylidene indandione,
  an embelin-like dihydroxyquinone, two naphthoquinone sulfonylimines — at 3–75 µM
  in a 387,000-compound screen. One is a cephalosporin whose warhead match is
  spurious. Validating a geometric criterion against compounds that hit
  everything would confirm nothing.
- **Negatives carry the same warhead as the positives.** A molecule with no
  warhead fails the criterion trivially for lack of a reactive atom, so an
  unmatched control would measure the presence of a warhead and nothing else
  (D0014 makes the same argument for covalent decoys).
- **Negatives are shuffled, not taken in file order.** A PubChem datatable is
  ordered by submission, which tracks depositor and therefore chemical series.

## What this does not establish

- **SNAr cannot be settled.** Only two SNAr ligands have ever been crystallised
  at Cys113. n = 2 is underpowered by construction — issue #12 §A's "3 verified
  chemotypes against a statistical floor of 6", arriving as a concrete inability
  to validate rather than as a projection.
- **15 positives supports "clearly works" or "clearly does not", not a marginal
  claim.** The rationale doc said this before the measurement and it still holds.
- **The negatives are HTS inactives**, weak evidence per compound. Single-
  concentration inactivity has many causes besides failing to bind.
- **Whether one draw of 30 negatives was lucky** is being measured separately by
  `scripts/nac_robustness.py`, against ten disjoint draws from a 300-per-class
  pool. Until that lands, the AUCs above carry only seed-replication support, not
  negative-draw support.
- **Nothing yet ranks among the molecules that pass.** Per-pose binding energies
  are now captured as stage 4's raw material; whether they add anything over
  frequency alone is untested.

## Consequences

- `scripts/nac_rank.py` scores all 5,769 warhead-bearing T_3/T_4 candidates on
  this quantity.
- **T_1 and T_2 are out of scope, not ranked last.** They carry no warhead, so
  "can it present its warhead" is undefined for them rather than false (D0043).
- Any verdict attached under a previous ranking is invalid, for the reason D0059
  gives: `attach_gate` keys on (stratum, metric) and the metric name does not
  change when its definition does.
