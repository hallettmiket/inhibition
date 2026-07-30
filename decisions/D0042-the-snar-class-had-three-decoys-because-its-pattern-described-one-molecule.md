---
id: D0042
title: The SNAr class had 3 decoys because its adduct pattern described one molecule, not a chemistry
date: 2026-07-30
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/reference/warhead_classes_9.csv
  - data/reference/decoy_chemotypes_3.csv
  - shared/decoys_classmatched.py
  - shared/warhead_library.py
  - scripts/build_covalent_decoys.py
  - scripts/run_enrichment_gate.py
evidence:
  - 'snar_chloroazine adduct-valid pool: 3 -> 1449 after relaxing the class test'
  - 'Tian 6a property-matched decoys: 0 -> 23; testable False -> True'
  - 'covalent gate: 2 actives / 2 chemotypes -> 3 actives / 3 chemotypes'
  - '[covalent/cnn_affinity] UNDERPOWERED AUC 0.656 CI[0.387,0.974] EF1% 0.0'
  - '[covalent/affinity_kcal] UNDERPOWERED AUC 0.518 CI[0.298,0.716] EF1% 0.0'
  - 'Tian adduct under the general class: attachment_idx 8, adduct SMILES identical to the narrow class'
  - 'generalising the SHARED class instead would make 30 of T_4 198 snar products degenerate'
  - 'sulfamate_acetamide fresh fetch: 6 retrieved, 6 adduct-valid, 0 and 1 property-matched — unrecoverable'
  - 'null 95% AUC range at 82 candidates: 2 actives [0.100,0.894], 3 actives [0.173,0.831]'
---

# A pattern that described one compound

## The symptom and the assumed cause

`snar_chloroazine` held **3** decoys against a gate minimum of 10, so Tian 6a —
one of only six covalent actives, and the sole representative of its chemotype —
could not enter the gate at all. The obvious reading was that ChEMBL simply does
not contain many nitro-activated chloropyrimidines, and that the chemotype was
too narrowly *defined* to be populated.

Both halves were wrong.

## The actual cause

The class's `adduct_attachment_smarts` was

```
[cH]1:[n]:[c]([NX3]):[c]([N+](=O)[O-]):[c]:[n]:1
```

a pyrimidine bearing **both** an amine and a nitro. That is not a description of
SNAr chemistry; it is a description of Tian 6a. Every other chloroazine failed
it. Of 175 genuine chloroazines in the cached ChEMBL pool, **0** passed —
including molecules whose chlorine sits between two ring nitrogens and which
undergo the identical displacement.

Chemically, chlorine on a ring carbon flanked by two ring nitrogens is
SNAr-activated by the ring nitrogens alone. The nitro raises the rate; it does
not create the mechanism. The pattern conflated a rate-enhancing substituent
with the reaction it enhances.

## The fix, and what it deliberately does not touch

A **separate** class, `snar_chloroazine_general`, with
`reactive_atom_smarts = [c]([Cl])([n])[n]` and
`adduct_attachment_smarts = [c;R]([n;R])[n;R]`.

Separate, not generalised in place, and the reason is measured: generalising the
shared class makes **30 of T_4's 198** enumerated SNAr products degenerate in
their attachment site. T_4's classes are ENUMERATION units; a decoy has to match
the active's CHEMISTRY. The project already documented exactly this split for
naphthoquinone. The new class carries `warhead_fragment_smiles = UNRESOLVED`, so
`enumerable()` can never select it — verified.

**The active is scored identically.** Under the general class Tian's adduct
attaches at atom index 8 with a byte-identical adduct SMILES. Nothing about the
molecule the gate scores has changed.

## Result

| | before | after |
|---|---|---|
| snar class pool, adduct-valid | 3 | **1449** |
| Tian property-matched decoys | 0 | **23** |
| covalent gate actives | 2 | **3** |
| covalent chemotypes | 2 | **3** |

```
[covalent/cnn_affinity]   UNDERPOWERED  AUC 0.656  CI[0.387, 0.974]  EF1% 0.0
[covalent/affinity_kcal]  UNDERPOWERED  AUC 0.518  CI[0.298, 0.716]  EF1% 0.0
  3 actives / 114 decoys / 3 chemotypes
```

Still UNDERPOWERED: 3 chemotypes against a floor of 6. What improved is the
null, from [0.100, 0.894] at 2 actives to [0.173, 0.831] at 3. Real, modest,
not a fix. EF1% is 0.0 on both metrics.

## What is now known to be unrecoverable

`sulfamate_acetamide` was refetched from ChEMBL rather than assumed: **6
molecules retrieved, all 6 adduct-valid**, and property matching to Reddi-4d and
4g yields **0** and **1**. The chemotype barely exists outside the Reddi paper.
Those two actives cannot be brought into the gate by any pattern change.

`Ieda-2019-(S)-2`, added earlier today, is excluded for a different and fixable
reason: its cinnamamide chemotype has no row in `decoy_chemotypes_3.csv` at all,
so no pool exists to match it against. That is a gap created by adding an active
without adding its chemotype definition, and it is open work.

## Three cache and write failures on the way, all the same shape

**The class-pool cache is keyed on `class_id` alone.** After the query was
relaxed, `build_class_pool("snar_chloroazine")` served the 3-molecule pool
retrieved with the OLD query and reported `class_pool_adduct_valid: 3` beside a
freshly measured 1449. A cache that cannot tell its inputs changed is D0033 and
the MD-length cache again. Now cached in a NEW directory, `class_pools_4/` under
`append_only`, so an old query's results cannot be served under a new query's
name.

**`warhead_library` still pointed at `warhead_classes_7.csv`.** The acrylamide
precedent correction written earlier the same day as `_8` was never loaded by
anything — a silent no-op, found only because the new class also failed to
appear. Writing a new version of a reference file is not the same as using it.

**An `unlink()` added so the builder could rewrite its output** would have
violated the append-only rule; the hook would have refused it. Output is
versioned to `_8` instead. The manifest write still collides with an existing
path and is not recorded — open work.

## The lesson

Three separate times today a correction was written and not connected: the
interaction split that no frame carried (D0039), the warhead file no loader
read, and the relaxed query no cache honoured. The pattern is not carelessness
about the fix; it is that a fix feels finished when the value is right, and the
value being right is the smaller half. The question that catches all three is
the same one: **which file does the consumer open, and does the change appear in
it?**
