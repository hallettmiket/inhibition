---
id: D0014
title: Covalent decoys must carry a warhead
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/decoys.py
  - data/reference/warhead_classes_3.csv
evidence:
  - 'decoys_covalent_1: only 32 of 302 (10.6%) carried any electrophile'
  - 'gnina covalent docking requires --covalent_lig_atom_pattern to MATCH the ligand'
  - 'decoys_covalent_2: 294 decoys, 0 without an electrophile'
  - '119 of 294 share their active warhead class; the ChEMBL pool holds 13 chloroacetamides and 0 sulfamate acetamides'
  - 'peptidomimetic BJP-06-005-3 still matches 0 warhead-bearing decoys and is excluded by the >=10 filter'
runbook: null
---

> **REVISION NOTICE (D0031, 2026-07-28).** The decision below is right and is
> now enforced properly. Its IMPLEMENTATION property-matched first and
> filtered for a warhead afterwards, which eliminated rare chemotypes before
> the warhead filter ran and topped up the shortfall across classes: 104
> acrylamide decoys against zero acrylamide actives. `decoys_covalent_2` must
> not be used for a covalent gate again; see `decoys_covalent_6.csv`.

## Context

The first covalent decoy set was property-matched on size, greasiness,
polarity and charge — the standard DUD-E criteria — but never checked whether a
decoy could *react*. Covalent docking requires a reactive atom to bond to;
gnina matches a SMARTS pattern against the ligand and a molecule with no
electrophile cannot be scored at all.

Only 10.6% of that set carried an electrophile. Nine decoys in ten were
unrunnable, and the comparison that survived would have been "electrophiles
versus inert molecules" — which docking wins for reasons that have nothing to do
with Pin1.

## Decision

Covalent decoys must carry a warhead motif, enforced at generation. Where the
pool allows, they are drawn from the SAME warhead class as the active they are
matched to.

That second part is what makes the control interesting rather than trivial. The
question becomes *does docking prefer our chloroacetamide over other
chloroacetamides* — discrimination within a chemistry — instead of *does docking
prefer electrophiles over inert molecules*, which is not in doubt and not
informative.

## Consequences

294 decoys, none of them unrunnable. But class-matching is only partial: 119 of
294. The ChEMBL pool holds 13 chloroacetamides and no sulfamate acetamides in
the relevant property range, so most decoys are Michael acceptors topped up from
outside the active's class. That fallback is recorded per decoy
(`class_matched`) rather than hidden, and it caps how strong a covalent verdict
can honestly be.

A targeted ChEMBL substructure search per warhead class would deepen the pool
and raise the class-matched fraction. Worth doing if the covalent gate turns out
to be the binding constraint on T_3/T_4.

The peptidomimetic BJP-06-005-3 still matches nothing and remains excluded by
the minimum-decoy filter, as before.
