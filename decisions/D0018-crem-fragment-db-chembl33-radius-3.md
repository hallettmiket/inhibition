---
id: D0018
title: CReM fragment DB — ChEMBL33 SA<=2 primary, Enamine secondary, radius 2
date: 2026-07-27
status: accepted
approach: t2
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/sources.yaml
  - config/approaches/t2_atra_crem.yaml
evidence:
  - 'published CReM DBs differ by SOURCE, SA filter and min fragment frequency — NOT by radius'
  - 'radius is a runtime argument to grow_mol/mutate_mol, not a property of the file'
  - 'chembl33_sa2_f5 is 281 MB; enamine2025_sa2_f5 is 777 MB'
  - "f5 matches the spec's stated --min-freq 5 default"
  - 'the DB contains radius1..radius5 tables — all radii ship in one file'
  - 'ATRA at radius 3, 4 and 5 yields ZERO mutations and ZERO grows'
  - 'ATRA at radius 2 yields 43 mutations / 38 grows; radius 1 yields 45 / 50'
  - 'control: benzoic acid at radius 3 yields 47 mutations — so radius 3 works in general, just not for this seed'
  - "CReM docstring: 'radius: radius of context which will be considered for replacement. Default: 3.'"
  - 'benzoic acid replacements by radius: 1449 / 593 / 584 / 530 / 526 for r=1..5 — MONOTONICALLY DECREASING'
  - 'benzoic acid by max_size at fixed radius 2: 40 / 204 / 593 / 593 for max_size=2/4/8/12 — INCREASING'
  - 'max_size defaults to 10 and was never set; ATRA has 22 heavy atoms'
runbook: docs/runbooks/adding_a_source.md
---

## Context

T_2 was blocked on "choosing the fragment-DB radius". That framing was wrong and
worth correcting: the published CReM databases differ by **source**, **synthetic
-accessibility filter** and **minimum fragment frequency**. The context
**radius is a runtime argument** to `grow_mol`/`mutate_mol`, not baked into the
file. Two separate decisions were hiding inside one.

## Decision

**Database — both staged, ChEMBL33 first.**

- **Primary: `chembl33_sa2_f5`** (281 MB). `f5` matches the spec's stated
  `--min-freq 5` default exactly. `sa2` restricts to synthetically accessible
  fragments, which serves condition (iii) *at enumeration time* rather than
  discovering unmakeable derivatives after labelling 10^4-10^5 of them.
- **Secondary: `enamine2025_sa2_f5`** (777 MB), staged now and unused for the
  first run. Enamine stock biases the neighbourhood toward fragments that can
  actually be bought — a different scientific stance from ChEMBL's "what has
  been published". Having it hash-pinned means that comparison is later a config
  change, not a fresh acquisition.

**Radius: 2** — revised from 3 on evidence, before any run.

Radius 3 was chosen first, on the reasoning that a larger radius demands more
surrounding context and so yields fewer but more chemically sensible
replacements. A smoke test against the actual seed refuted it:

| radius | ATRA mutate | ATRA grow |
|---|---|---|
| 1 | 45 | 50 |
| **2** | **43** | **38** |
| 3 | **0** | **0** |
| 4 | 0 | 0 |
| 5 | 0 | 0 |

**Radius 3 produces nothing at all for ATRA.** A control rules out a broken
setup: benzoic acid at radius 3 gives 47 mutations from the same database. The
cause is the seed itself — ATRA is a conjugated polyene with methyl branches, and
requiring three bonds of matching context finds no precedent for it anywhere in
ChEMBL33.

So radius 3 is not "conservative" here, it is inoperable: T_2 would have
enumerated an empty frontier and reported success. Radius 2 retains more context
than radius 1 while still being productive, and is the operative choice.

## What `radius` actually means — a naming trap worth recording

The natural reading is that a larger radius permits larger changes and so yields
MORE variants. It is the opposite. CReM's own docstring:

> `radius: radius of context which will be considered for replacement.`

The radius is the amount of surrounding structure that must **match** between
the query molecule and a database fragment before a swap is permitted. A larger
radius is a STRICTER constraint, so it yields FEWER replacements. Measured on
benzoic acid, where every radius works:

| radius | replacements |
|---|---|
| 1 | 1,449 |
| 2 | 593 |
| 3 | 584 |
| 4 | 530 |
| 5 | 526 |

The parameter that behaves the way the name suggests is **`max_size`** — heavy
atoms in the fragment being swapped in. At fixed radius 2: 40 / 204 / 593 / 593
for `max_size` 2 / 4 / 8 / 12. Two distinct knobs; conflating them is easy.

This also sharpens what happened with ATRA. Benzoic acid only drifts 593 → 526
from radius 2 to 5, a gentle decline. ATRA goes 43 → **0**. That is not a filter
tightening gradually — it is ATRA's conjugated polyene having NO 3-bond
precedent anywhere in ChEMBL33 while an ordinary aromatic has plenty. Radius 2
is not a mild relaxation for this seed; it is the last radius at which the seed
has any precedent at all.

A related consequence: at radius 1 ATRA gives 45 against 43 at radius 2 — nearly
identical, unlike benzoic acid's 2.4x jump. So the radius-1-vs-2 choice matters
far less for ATRA than it would for a typical scaffold.

## `max_size`: pinned, not inherited

`max_size` defaults to **10** heavy atoms and was never set explicitly. ATRA has
22 heavy atoms, so the default permits swapping nearly half the molecule — hard
to call a *derivative* neighbourhood. It is now pinned in the T_2 config rather
than inherited silently, because an unstated default is exactly the kind of
parameter that turns out to have been doing the work.

## Measured frontier with the pinned parameters

ATRA at radius 2, after pinning:

| operator | parameters | products |
|---|---|---|
| `mutate_mol` | `max_size 8, min_size 1, max_inc 4` | 203 |
| `grow_mol` | `max_atoms 8, min_atoms 1` | 1,679 |
| union | | **1,882** |

`mutate` at the unpinned defaults gives 239, so the pins cost about 15% of the
mutation products — a modest price for bounding what "derivative" means.

**`grow` dominates the frontier by roughly 8:1.** Worth knowing before tuning
either operator: the neighbourhood is mostly ATRA-plus-a-substituent, not
ATRA-with-a-piece-swapped. Degree-1 over ~1,882 products is a tractable starting
frontier, well inside the 200,000 cap.

A signature trap alongside the naming one: `grow_mol` does NOT accept
`max_size`/`min_size`. It sets them internally (it appends to hydrogens rather
than replacing a fragment) and takes `min_atoms`/`max_atoms`. Passing `max_size`
to it raises `TypeError`, so the two operators are configured separately.

## Consequences

The radius is the parameter that actually defines the neighbourhood, so changing
it changes what "degree-1 derivative of ATRA" means and invalidates comparison
with any earlier run.

**A general lesson for retargeting:** the usable radius is a property of the
SEED, not of the method. Any new seed needs this smoke test before its radius is
pinned — an unusual scaffold can silently produce an empty neighbourhood at a
radius that works fine for ordinary drug-like molecules. It is recorded in the T_2 config and pinned into each
run's manifest rather than passed ad hoc.

Running the Enamine database is a documented extension: swap the source in the
T_2 config, re-run, and compare. Both are hash-pinned so the comparison is
between two known inputs.

**Verification gap, stated plainly:** the file list and the SA/frequency
semantics come from the CReM download page. The radius guidance is a reading of
what the parameter does, not a quoted recommendation — the docs page gave none.
Worth checking against crem.readthedocs.io before any published result rests
on it.
