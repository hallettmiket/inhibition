---
id: D0037
title: The junction dihedral was the sp3 analogue, and the reported dG was never an interaction energy
date: 2026-07-29
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - data/params/cys_gaff2_junction_5.frcmod
  - shared/mmgbsa.py
  - shared/mmgbsa_ensemble.py
  - shared/enrichment_gate.py
  - decisions/D0035-the-sp2-junction-gap-was-three-missing-angles.md
  - decisions/D0036-ensemble-mmgbsa-is-precise-and-still-below-chance.md
evidence:
  - "the frcmod header states every term is GAFF2's ss analogue; the whole DIHE block said 'from parm19' and used 1.00/3-fold/0deg for all five carbon types"
  - 'GAFF2 gives X-c2-ss-X 2.200/2-fold/180, X-ca-ss-X 0.800/2-fold/180; only X-c3-ss-X matched what was in use'
  - "Juglone's built topology carried 2C-S-cc at PK 0.333, per 3, phase 0.0 -- the sp3 form on an sp2 attachment"
  - 'affected 31 of 82 gate ligands: every cc/cd/ca/c2 attachment, i.e. exactly the ligands D0035 restored'
  - 'GAFF2 has NO generic X-cc-ss-X or X-cd-ss-X, only specific 4-atom terms; cc/cd use cd-cc-ss-ca (2.430/2-fold/0) and are the least certain of the five'
  - 'link-atom cap prevents bonded-term cancellation: residual 9.77 +/- 11.38 kcal/mol against a decoy spread of 6.35'
  - 'gate ROC-AUC: full potential 0.425, standard interaction energy 0.181, the residual ALONE 0.831'
  - 'Juglone interaction energy +9.08 kcal/mol, i.e. unfavourable, for a compound with published covalent activity'
  - 'the measurement-error propagation quoted in D0036 existed in no code in the repository'
---

# Two defects the adversary found

Both were found by an adversarial audit that was asked to refute, not confirm.
Neither would have surfaced from the runs succeeding, because both produced
plausible numbers.

## 1. The junction dihedral was the wrong analogue

`cys_gaff2_junction_4.frcmod` states in its own header that every term is the
GAFF2 parameter for the same geometry with GAFF2's thioether sulfur `ss` in
place of the protein's `S`. The entire `DIHE` block violated that, and said so
in its own comments -- `generic, from parm19 X-CT-S-X` -- applying
1.00 kcal/mol, 3-fold, 0 degrees to all five attachment carbon types.

GAFF2's actual values:

| junction | in use | GAFF2 |
|---|---|---|
| `c3` (sp3) | 1.000 / 3 / 0 | 1.000 / 3 / 0 — correct |
| `c2` | 1.000 / 3 / 0 | **2.200 / 2 / 180** |
| `ca` | 1.000 / 3 / 0 | **0.800 / 2 / 180** |
| `cc`, `cd` | 1.000 / 3 / 0 | no generic exists |

Wrong barrier, wrong periodicity, wrong phase for four of five. The 2-fold
180-degree form is what enforces the conjugation plane at an sp2 carbon-sulfur
bond; a 3-fold 0-degree term does not.

It reached the topologies: Juglone carried `2C-S-cc` at PK 0.333, per 3, phase
0.0. That is 31 of the 82 gate ligands -- every `cc`/`cd`/`ca`/`c2` attachment,
which is precisely the set D0035 had just restored and D0036 built its
two-active result on.

**How it survived D0035.** That record argued the derivation rule could be
trusted because it reproduced the two entries already present. The argument is
sound. It was applied to the ANGLE block and never to the DIHE block, which
sits sixty lines below and contradicts the header in plain text.

**The fix and its weakest point.** Junction v5 uses GAFF2's values. For
`cc`/`cd` GAFF2 provides no generic form at all, only specific four-atom terms,
and our terminal atom is the protein `2C`, which none of them match. We use the
2-fold member of that set (`cd-cc-ss-ca`, 2.430/2/0). That is a judgement call,
it is the least defensible of the five, and it wants a sensitivity check before
any conclusion rests on the naphthoquinone chemotypes.

The 45 affected candidates rescored automatically, because the cache
fingerprint includes the frcmod name. Gate single-structure dG moved 0.350 to
0.425.

## 2. The reported dG was never an interaction energy

Single-trajectory MM/GBSA is DEFINED by the bonded terms cancelling between the
three legs, leaving dVDW + dEEL + dEGB + dESURF. The link-atom cap breaks that:
the complex has an S-C bond where the two legs have S-H and C-H, so BOND,
ANGLE, DIHED, 1-4 and CMAP do not cancel.

The remainder is not a rounding error. Across the gate set it is
**9.77 +/- 11.38 kcal/mol**, larger than the 6.35 kcal/mol spread of the decoys
the ranking has to discriminate within.

| quantity | gate ROC-AUC |
|---|---|
| full potential difference (what was reported as "dG") | 0.425 |
| standard MM/GBSA interaction energy | **0.181** |
| the non-cancelling internal residual, **alone** | **0.831** |

The artefact is the best classifier in the dataset, and the reported score was
propped up by it. On the correct quantity the two known actives rank 56th and
78th of 82, and Juglone's interaction energy is **+9.08 kcal/mol** --
unfavourable, for a compound with published covalent activity.

**Decision: report both, and name the remainder.** `dG_kcal` keeps its existing
meaning so nothing downstream shifts silently; `dG_interaction_kcal` is the
standard quantity; `dG_internal_residual_kcal` is the remainder, no longer
folded invisibly into a number people read as a binding energy.
`INTERACTION_TERMS` and `INTERNAL_TERMS` are asserted at import to partition
`ENERGY_TERMS` exactly.

**What 0.831 probably is.** A term that separates two molecules from eighty is
more plausibly encoding chemotype than binding, and with one active per
chemotype the two are perfectly confounded in this set. That is the confound
the class-matched decoys were built to remove from docking (D0031), reappearing
in a different term. It should not be read as signal until it is tested against
heavy-atom count and class.

## 3. Two smaller findings from the same audit

**A number in the manuscript existed in no code.** The measurement-error
propagation behind `P(AUC > 0.5) = 0.002` was computed ad hoc, unseeded,
untested, and never routed through the gate that grades every other metric. It
is now `enrichment_gate.propagate_measurement_error` -- seeded, versioned,
tested, and attached to `GateResult` in its own field, kept separate from the
bootstrap CI because the two answer different questions.

A related claim of ours did not survive its own test: the measurement interval
is not inherently the narrower, flattering one. Which interval is wider depends
on the per-candidate error relative to the spread between candidates. True of
our data, not true in general, and now pinned as a test.

**The Spearman 0.283 did not reproduce.** D0036 argued that single structures
were noise-dominated on the gate set, citing rho = 0.283 between one frame and
the trajectory mean. Recomputed across all 90 frames on the complete set, rho
runs 0.655 to 0.887, median 0.822. The original figure came from a partial set
of 37 of 82 ligands assembled in identifier order while the remaining runs were
still in flight -- the exact selection bias this project had refused to accept
in the gate one message earlier, and then committed in an analysis. The
argument that rested on it is withdrawn.

## The pattern, restated

D0033 was a value derived from a name that did not match. D0034 was a
destructive default. These are the same family: a rule verified on one section
of a file and assumed for the rest, and a quantity whose definition was never
stated so nobody checked which one was being computed. In every case the run
succeeded and the number looked reasonable.
