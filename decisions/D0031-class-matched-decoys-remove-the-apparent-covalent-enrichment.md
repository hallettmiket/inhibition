---
id: D0031
title: Class-matched decoys remove the apparent covalent enrichment
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/decoys_classmatched.py
  - scripts/build_covalent_decoys.py
  - scripts/run_enrichment_gate.py
  - data/reference/decoy_chemotypes_2.csv
  - decisions/D0014-covalent-decoys-must-carry-a-warhead.md
  - decisions/D0015-covalent-ranking-uses-affinity-not-cnnaffinity.md
  - decisions/D0028-enrichment-gate-remeasured-on-adduct-forms.md
evidence:
  - 'class-matched gate: affinity_kcal ROC-AUC 0.537, CI [0.346, 0.728], EF1% 0.0, BEDROC 0.001'
  - 'class-matched gate: cnn_affinity ROC-AUC 0.552, CI [0.358, 0.741], EF1% 0.0, BEDROC 0.001'
  - 'the same metric measured 0.815 on decoys_covalent_2 (D0015) and 0.718 on adduct forms (D0028)'
  - 'decoys_covalent_2 held 104 acrylamide decoys against ZERO acrylamide actives'
  - 'new set: 90 decoys, every one carrying its chemotype whole reactive group and producing a valid adduct'
  - 'Sulfopin 50/50 same-class decoys; Juglone 31; BJP-06-005-3 8; Tian-6a 0; Reddi-4d 0; Reddi-4g 1'
  - 'ChEMBL holds 4,430 chloroacetamides, 3,963 naphthoquinones, 41 sulfonate acetamides, 6 sulfamate acetamides, and 3 non-Pin1 nitro-chloropyrimidines'
  - 'max ECFP4 Tanimoto of any decoy to any active: 0.344 against a 0.35 cap'
  - 'both verdicts remain UNDERPOWERED: 2 chemotypes against a floor of 6'
---

# Class-matched decoys remove the apparent covalent enrichment

## The result

Scored against decoys of **its own chemotype**, covalent docking does not
separate known Pin1 actives from non-binders:

| decoy set | affinity_kcal ROC-AUC | EF1% |
|---|---|---|
| `decoys_covalent_2`, pre-reaction ligands (D0015) | 0.815, CI [0.667, 0.931] | 16.7 |
| `decoys_covalent_2`, adduct forms (D0028) | 0.718, CI [0.483, 0.944] | 19.0 |
| **`decoys_covalent_6`, class-matched adducts** | **0.537, CI [0.346, 0.728]** | **0.0** |

D0028 suspected that some unknown share of the separation was chemotype
rather than binding. It was most of it. The old set held **104 acrylamide
decoys against zero acrylamide actives**, so the gate was substantially
asking "can docking tell an acrylamide from a chloroacetamide" — a
question D0020 says is not even well-posed, since affinity is not
comparable across warhead classes.

`cnn_affinity` now scores 0.552 against affinity's 0.537. The two are
indistinguishable, so **D0015's preference for affinity has no remaining
empirical support** either. It is retained only on gnina's own warning
that CNN scoring is uncalibrated for covalent docking — a mechanistic
argument, not a measured one.

## What this does and does not license

**Does:** T_3's and T_4's docking rankings must be presented as
*unvalidated ordering*. There is no demonstrated ability to enrich for
binders within a chemotype on this receptor.

**Does not:** conclude that docking fails. Only two actives survive
class-matching, so the test has little power in either direction. The
honest statement is that the previously reported enrichment was inflated
by chemotype mismatch, and what remains is indistinguishable from chance
on a sample too small to resolve it. The 6-chemotype floor held both
verdicts at UNDERPOWERED throughout, which is now the third time it has
prevented a flattering point estimate from being promoted.

## How the set is built

Retrieval and verification are deliberately **different patterns**:

1. **Retrieve loosely** — one ChEMBL substructure search per chemotype
   (`NC(=O)CCl`), cached under `immutable/`.
2. **Verify strictly** — membership requires the chemotype's *whole*
   reactive group (`[NX3]C(=O)[CH2][Cl]`) **and** a successful adduct
   transform. A decoy that cannot be docked as an adduct cannot sit in a
   control for approaches that are (D0022, D0030).
3. **Property-match within the class**, then reject anything above 0.35
   ECFP4 Tanimoto to any active.
4. **Never top up across classes.** An active without enough same-class
   decoys is reported untestable and dropped.

Inverting the old order is the whole point. `shared.decoys.build`
property-matched against one generic 20k pool and filtered for a warhead
afterwards, so rare chemotypes were eliminated before the warhead filter
ran and the shortfall was topped up from whatever survived.

## Chemotype is not the same thing as a warhead class

Decoy membership is decided by a **chemotype** table
(`decoy_chemotypes_2.csv`), separate from the warhead library. The
library's classes are T_4's *enumeration units*: `naphthoquinone_c2` and
`naphthoquinone_benzo` are two attachment positions on one chemistry.
Juglone is a genuine 1,4-naphthoquinone whose hydroxyl sits at neither,
so classifying it by an enumeration unit rejected it from its own class.

## What the literature does not contain

Three actives cannot be tested, and not for want of trying:

- **Tian-6a** (sNAr): ChEMBL holds 44 nitro-chloropyrimidines, of which
  **41 are Pin1 compounds from the Tian 2025 paper itself**. Three
  remain, two of which produce a valid adduct.
- **Reddi-4d / 4g** (sulfamate acetamide): six exist in all of ChEMBL.
  The Reddi 2023 chemotype is close to unprecedented.
- **BJP-06-005-3**: a peptidomimetic far outside the property range of
  any available chloroacetamide; 8 matched.

These are limits of the published chemical record, not of the code.

## KPT-6566 is excluded

It is a self-immolative aryl-sulfonyl-acetate: the species that
alkylates Cys113 is a naphthoquinone **released** from it, not the
deposited molecule. Docking the parent as a naphthoquinone would score a
structure that never reaches the cysteine — the same "which species is
being scored" error as D0022 and D0030. Testing it needs the released
fragment as its own entry.

## Three defects found while building this, all silent

Each produced a plausible, wrong scientific claim rather than an error:

1. **SMILES parsed where SMARTS was meant.** The warhead group query was
   built with `MolFromSmiles`, where `[*]` is a dummy atom matching only
   another dummy, instead of `MolFromSmarts`, where it means any atom. It
   matched 0 of 2,038 chloroacetamides and 0 of 3,963 naphthoquinones,
   and every chemotype was reported chemically unavailable while the
   pools were full.
2. **A ketone query for an amide warhead.** `CC(=O)CCl` is
   C–C(=O)–C–Cl, a chloro-*ketone*. These warheads are acet*amides*,
   N–C(=O)–CH2–Cl. It returned 2,038 molecules of which 7 were amides,
   and made sulfamate acetamide look non-existent (0 hits) when the
   amide query finds 6.
3. **A transient empty response cached as an absence.** One ChEMBL call
   returned HTTP 200 with an empty `molecules` list. That was written to
   the cache and read back as "this chemotype does not exist in the
   database" — for chloroacetamide, Sulfopin's own class, of which
   ChEMBL holds 4,430.

The first two were caught by a self-check that **every active must
satisfy its own class test**: a control whose own positives fail it is
not a control. The third is now caught by refusing to cache an empty
result when ChEMBL reports a non-zero `total_count` — a retrieval
failure and a real absence must never look the same.

## Consequences

- `decoys_covalent_2` must not be used for a covalent gate again.
- D0015's metric choice stands on mechanism, not measurement.
- The integration GUI must show the covalent ranking with this verdict
  attached, not the D0015 figure.
- Strengthening this gate requires more actives, which means the
  literature, not more compute.
