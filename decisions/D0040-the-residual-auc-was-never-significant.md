---
id: D0040
title: The residual's enrichment AUC was never significant, and size is not the explanation
date: 2026-07-30
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - decisions/D0037-the-junction-dihedral-was-the-sp3-analogue-and-dG-was-not-an-interaction-energy.md
  - shared/enrichment_gate.py
evidence:
  - 'exact permutation over all C(82,2) = 3321 relabellings of 2 actives among 82 candidates'
  - 'null 95% range of the AUC: [0.106, 0.894] -- identical for every score tested'
  - 'internal residual AUC 0.700, p = 0.188'
  - 'interaction energy AUC 0.237, p = 0.886'
  - 'full potential (the column previously displayed) AUC 0.350, p = 0.755'
  - 'heavy-atom count used alone AUC 0.469, p = 0.559'
  - 'corr(residual, heavy atoms) Pearson +0.223, Spearman +0.278'
  - 'the two actives carry 37 and 19 heavy atoms; decoys median 28, range 18-41'
  - 'D0037 quoted 0.831 and 0.181 from a partial 37-of-82 subset; the complete set gives 0.700 and 0.237'
---

# Two actives cannot support any of it

## What D0037 said

That the internal residual separated actives from decoys with an AUC of **0.831**
while the standard interaction energy managed **0.181** — the modelling artefact
outperforming the physics, and the physics running backwards. It was offered as
one of the more alarming results in the project, and the leading explanation was
that the residual encodes chemotype: with roughly one active per chemotype,
chemotype and label are perfectly confounded.

The open question left behind was whether the residual is a chemotype or
molecular-size proxy. That question has an answer, and a prior question has a
different one.

## The AUC is not distinguishable from random labelling

There are 2 actives among 82 candidates. Every possible way of choosing which
2 of those 82 molecules are the "actives" gives 3321 relabellings, and computing
the AUC under each gives the exact null distribution — no sampling, no
approximation.

**The null 95% range is [0.106, 0.894].** Every score tested falls inside it:

| score | AUC | exact permutation p |
|---|---|---|
| internal residual | 0.700 | 0.188 |
| full potential (what the dossier displayed) | 0.350 | 0.755 |
| interaction energy | 0.237 | 0.886 |
| heavy-atom count, used alone | 0.469 | 0.559 |

Nothing is significant. An AUC of 0.831 or 0.181 sits comfortably within what
relabelling two arbitrary molecules produces, so neither number was evidence
about anything. The comparison between them — artefact beats physics — was a
comparison of two draws from the same null.

Two secondary corrections fall out. The 0.831 and 0.181 were computed on a
**partial 37-of-82 subset**; on the complete set they are 0.700 and 0.237. And
the direction that looked most alarming, interaction energy running backwards at
0.181, is 0.237 on complete data with p = 0.886 — the strongest-looking result
in D0037 was the least supported.

## Size is not the explanation

The obvious mechanism was extensivity: internal energy scales with atom count,
so a larger ligand has more bonded terms and a larger residual, with no
reference to the pocket at all. If so, heavy-atom count alone should reproduce
the AUC.

It does not. Heavy-atom count scores **0.469** — chance — and correlates with
the residual at only Pearson +0.223 / Spearman +0.278. The two actives carry 37
and 19 heavy atoms against a decoy median of 28, i.e. one above and one below,
which is why size carries no signal here.

So the chemotype/size hypothesis is not supported. But that finding is close to
moot: at p = 0.188 there is no effect that needs explaining.

## What stands and what is withdrawn

**Stands — the interaction/residual split.** It rests on physics, not on this
AUC: single-trajectory MM-GBSA is defined by bonded-term cancellation, and the
link-atom caps break it. D0039 confirmed the split empirically by a route that
does not involve actives at all — T_1 and T_2, having no link atom, return a
residual of exactly 0.00, and T_3 returns 2.88. That control is unaffected by
anything here.

**Withdrawn — every enrichment statement in D0037.** Specifically: that the
residual separates actives from decoys better than the interaction energy; that
the interaction energy is anti-correlated with activity; and the chemotype
confounding argument, which explained an effect that was not established.

No shortlist, gate or rank changes, because none consumes these values.

## The gate already knew

`min_actives_for_verdict = 3`. The gate has returned UNDERPOWERED throughout and
refused to certify anything on 2 actives. It was right, and the prose written
around it was not: the power floor was enforced on the *verdict* and ignored in
the *commentary*, where the same numbers were narrated as findings.

That is the correction worth keeping. A quantity the gate declines to act on
should not appear in a decision record as a result, in the manuscript as an
observation, or in the dossier as a column a reader is invited to interpret. The
guard existed and did its job; what leaked past it was me describing the numbers
it had already rejected.

## What would change this

Three or more actives clears the gate's floor, and the practical target is
enough that a permutation null is narrower than the effect being claimed. At 2
actives the null spans 0.79 of the AUC range; nothing measured against it can
be informative. More actives is the single highest-value addition to this
project, and it is a data problem rather than a compute one.
