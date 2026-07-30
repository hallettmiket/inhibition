---
id: D0035
title: The sp2 junction gap that cost 32 ligands and an active was three missing angle lines
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - data/params/cys_gaff2_junction_3.frcmod
  - shared/mmgbsa.py
  - decisions/D0032-mmgbsa-gate-and-the-power-floor-on-negative-verdicts.md
evidence:
  - "tleap error was always the same single line: 'Could not find angle parameter for atom types: 2C - S - cc'"
  - 'junction v2 declared BONDS S-c2, S-c3, S-ca, S-cc, S-cd but only ANGLES 2C-S-c2 and 2C-S-c3'
  - 'so every ligand attaching through aromatic (ca) or conjugated (cc/cd) carbon failed at the last build step'
  - 'cost: 32 of 83 gate ligands, 5 of 7 T_4 adduct classes, and Juglone -- the second active'
  - 'fix is 3 lines, taken from gaff2.dat by the analogue rule the file already documented'
  - 'the 2 pre-existing entries match gaff2.dat exactly (c2-ss-c3, c3-ss-c3), confirming the rule before applying it'
  - 'Juglone now builds: tleap Errors = 0, all three prmtops written'
---

# Three missing angle lines

> **REVISION NOTICE (D0037, 2026-07-29).** This record verified the derivation
> rule against the ANGLE block only. The DIHE block violated the same rule --
> it applied GAFF2's sp3 analogue (1.00 / 3-fold / 0 deg) to every sp2
> attachment, where GAFF2 gives X-c2-ss-X 2.200/2-fold/180 and X-ca-ss-X
> 0.800/2-fold/180. Corrected in junction v5. The claim below that "a
> derivation rule that reproduces the existing values can be trusted to produce
> the missing ones" is sound, but it was applied to one section of the file and
> not the other.

## What was blocking

Every failure carried the same tleap message:

```
Could not find angle parameter for atom types: 2C - S - cc
```

`2C` is ff19SB's Cys CB, `S` its SG, `cc` a GAFF2 conjugated carbon. Junction
version 2 declared the **bonds** for all five carbon types a ligand might
attach through — `S-c2`, `S-c3`, `S-ca`, `S-cc`, `S-cd` — but only two of the
matching **angles**, `2C-S-c2` and `2C-S-c3`.

So a ligand attaching through an aromatic or conjugated carbon got a valid bond
parameter, passed every earlier stage, and died at the last step of topology
construction.

## What it cost

- **32 of 83** gate ligands
- **5 of 7** T_4 adduct classes (bdhi_c4, bdhi_c5, naphthoquinone_benzo,
  naphthoquinone_c2, snar_chloroazine)
- **Juglone** — the second active

The last item is the expensive one. D0032 graded the MM-GBSA gate UNDERPOWERED
on **one** active, and recorded that "restoring Juglone alone doubles the
actives and the chemotypes". The single most limiting fact about this build's
statistical power was three absent lines in a parameter file.

## The fix

The file already documented its own derivation rule: take the GAFF2 parameter
for the same geometry, substituting GAFF2's thioether sulfur `ss` for the
protein `S`, and `c3` for the protein CB (which is sp3). Applying it:

```
2C-S -ca    92.72     102.53     from gaff2 c3-ss-ca
2C-S -cc    94.13     101.06     from gaff2 c3-ss-cc
2C-S -cd    93.40     102.22     from gaff2 c3-ss-cd
```

The rule was **verified before being used**: the two entries already present
match `gaff2.dat` to the digit (`c2-ss-c3` = 93.32/101.26, `c3-ss-c3` =
92.67/99.64). A derivation rule that reproduces the existing values is one that
can be trusted to produce the missing ones.

Juglone now builds — `Errors = 0`, all three topologies written.

## The invariant that was missing

**The angle section must cover every carbon type the bond section covers.**

A bond parameter without its matching angle is not partial support. It is a
build that consumes a docking pose, a charge derivation and an antechamber run,
and then fails. Worse, it fails *per ligand*, so it reads as a property of that
molecule rather than as a systematic hole — which is why it was carried for
several decisions as "the sp2 junction gap" rather than fixed.

## A format trap, recorded because it cost three attempts

`tleap` parses the frcmod header loosely and will read prose there as parameter
specifications. Two drafts of this file failed with

```
Unknown keyword ("Could not find angle parameter" for every ligand attaching...
```

— the parser consuming the *explanation of the bug* as a parameter. Version 3
is therefore a **minimal diff** from version 2: one changed title line and
three added angle lines. All rationale lives here instead.

The general form: a file consumed by a strict external parser is not the place
for narrative. Put the reasoning in the decision record and keep the artefact
machine-shaped.

## What this does not settle

Restoring Juglone takes the covalent gate from one active to two. Two is still
below the `min_actives_for_verdict` floor of three that D0032 established, so
the verdict will remain UNDERPOWERED and *should*. This buys statistical power,
not a conclusion. The gate's own floor decides what may be claimed, and two
actives does not clear it.
