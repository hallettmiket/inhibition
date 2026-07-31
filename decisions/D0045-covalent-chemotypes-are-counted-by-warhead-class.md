---
id: D0045
title: Covalent chemotypes are counted by warhead class, not by structural similarity — and that leaves the gate underpowered
date: 2026-07-31
status: accepted
approach: shared
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/enrichment_gate.py
  - shared/warhead_library.py
  - decisions/D0031-decoys-are-matched-by-warhead-class.md
  - decisions/D0041-the-first-verdict-docking-does-not-demonstrably-enrich.md
evidence:
  - 'DECIDED BEFORE THE COUNTS WERE COMPUTED ON NEW DATA — see "Why the order matters"'
  - 'covalent lead-tier actives: 8 molecules'
  - 'ECFP4 Butina @0.4 on those 8: 6 clusters -> would CLEAR the floor of 6'
  - 'canonical warhead class on those 8: 4 classes -> does NOT clear the floor'
  - 'free-text warhead_class nunique() on those 8: 6 -> would clear the floor, but is an artefact of prose'
  - 'the prose overcount: "chloroacetamide" and "chloroacetamide (N-methyl peptidomimetic)" are one warhead, two strings'
  - 'the 4 canonical classes: chloroacetamide, sulfamate_acetamide, cinnamamide, snar_chloroazine'
  - 'Sulfopin (chloroacetamide) and Reddi-4d/4g (sulfamate) fall in ONE ECFP4 cluster — shared sulfolane scaffold'
  - 'BJP-06-005-3 (chloroacetamide, same warhead as Sulfopin) falls in a DIFFERENT ECFP4 cluster'
  - 'adding ZL-Pin13 raised n_actives but added NO chemotype — it is a third chloroacetamide'
  - '164A10 has an unestablished warhead and contributes no class'
  - 'non_covalent keeps ECFP4: 7 lead actives -> 6 clusters (a non-covalent binder has no warhead)'
  - 'the Pin1 PDB survey holds >=4 further covalent chemistries absent from the reference set: aryl aldehyde (11 structures), maleate/fumarate ester (4), SuFEx (2), Mannich (1)'
---

# One definition to build the comparison, another to size it

## The incoherence

The gate has been using two different notions of "independent" in two places.

**Decoys are matched by warhead class** (D0031), because affinity is not
comparable across warheads — a chloroacetamide and a naphthoquinone are not
competing on the same scale, so the comparison is built *within* warhead class.

**Chemotypes were counted by whole-molecule ECFP4 similarity**, Butina at 0.4.

On this actives set the two are close to orthogonal, and they disagree in both
directions:

| | warhead says | structure says |
|---|---|---|
| Sulfopin vs Reddi-4d/4g | different (chloroacetamide vs sulfamate) | **same** — shared sulfolane scaffold |
| Sulfopin vs BJP-06-005-3 | same (both chloroacetamide) | **different** |

Building a comparison on one axis and sizing it on another is not defensible.
The floor exists to stop us claiming enrichment from what is effectively one
compound wearing several hats; it can only do that if "one compound" means the
same thing as the comparison it guards.

## The decision

**Covalent stratum: count by canonical warhead class.**
**Non-covalent stratum: keep ECFP4 clustering** — a reversible binder has no
warhead, and there is nothing else to use.

The gate token now records `chemotype_method` alongside `n_chemotypes`, because
two defensible definitions give 4 and 6 against a floor of 6. A bare count is
not interpretable without knowing which produced it.

## Why the order matters, and what it cost

**This was decided before the counts on the new data were computed.** That was
deliberate. The definition determines whether the covalent gate can return a
verdict at all, and choosing it after seeing which option clears the floor is
choosing the answer. D0031's own warning applies to whoever holds the pen: a
floor that has blocked every covalent verdict should not be relaxed by the
person who wants the verdict.

Having decided, the result is the unfavourable one:

* structural clustering would give **6** — exactly the floor, verdict unblocked
* canonical warhead class gives **4** — still UNDERPOWERED

So this decision **keeps the covalent gate blocked**. That is the honest
outcome and it is the reason the ordering mattered.

## A third count that is simply wrong

Counting `nunique()` on the reference set's free-text `warhead_class` column
also gives **6** — and it is an artefact. That column is prose written per row:

    "chloroacetamide"
    "chloroacetamide (N-methyl peptidomimetic)"

One warhead, two strings, two classes. Naive string uniqueness OVERCOUNTS, and
it overcounts in the direction that certifies a verdict. `warhead_library.
canonical_class()` maps prose to the controlled `class_id` vocabulary and
**raises** on anything unmapped, so a newly added active cannot become its own
chemotype by default — that would be the same failure in a new costume.

## What this says about the two actives just added

`pin1_reference_binders_4.csv` recovered ZL-Pin13 and 164A10, and it is worth
being precise about what that bought. `min_actives_for_verdict` (3) is
comfortably satisfied either way. But **ZL-Pin13 is a third chloroacetamide and
adds no chemotype**, and 164A10's warhead is unestablished so it adds none
either. The actives shortage eased; the chemotype shortage did not.

## What would clear the floor

Not more chloroacetamides. Two more independent warhead chemistries.

The Pin1 PDB survey (190 entries, 39 covalent) contains at least four
chemistries absent from our reference set — **aryl aldehyde** (11 structures,
a reversible-covalent mechanism we do not model at all), **maleate/fumarate
ester** (4), **SuFEx fluorosulfonyl** (2, unpublished deposits), and a
**Mannich base / aminoketone** (1). Curating those would take the count from 4
to 8 legitimately, without loosening anything.

That makes this decision and the PDB curation task a single piece of work
rather than two: the definition was tightened here, and the structures are what
pay for it.

## What did not change

No shortlist, rank or candidate is affected — this changes only how the gate
sizes its own evidence. The non-covalent stratum is untouched. D0041's
non-covalent verdict stands exactly as recorded.
