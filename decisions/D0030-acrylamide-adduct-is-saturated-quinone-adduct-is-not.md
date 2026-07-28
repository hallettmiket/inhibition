---
id: D0030
title: Acrylamide's adduct is saturated; the quinones' is not — one mechanism, two chemistries
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/covalent_adduct.py
  - shared/covalent_dock_run.py
  - data/reference/warhead_classes_5.csv
  - approaches/t3_reinvent/03_covalent_dock.py
  - approaches/t4_combinatorial/03_covalent_dock.py
  - decisions/D0022-dock-the-adduct-not-the-pre-reaction-ligand.md
evidence:
  - 'D0022 transformed only classes with a leaving group; all three michael_addition classes were passed through untouched'
  - '561 of 1,782 T_4 rows carried the resulting approximation note (acrylamide, naphthoquinone_c2, naphthoquinone_benzo)'
  - 'acrylamide is T_3''s ONLY warhead, so 100% of T_3 was affected'
  - 'docking the alkene gives Cys-S-CH=CH-C(=O)NR2, a planar vinyl thioether; the true adduct is Cys-S-CH2-CH2-C(=O)NR2 with two rotatable bonds'
  - 'thiol addition to 1,4-naphthoquinone gives a hydroquinone that re-oxidizes to the 2-thio-quinone, so the quinone sulfur genuinely sits on an sp2 carbon'
  - 'protocol fingerprint 67366274f425a371 -> a2854a6e6f7edc43'
  - 'T_3 smoke test: 12/12 docked on the saturated adduct'
---

# Acrylamide's adduct is saturated; the quinones' is not

## What was wrong

D0022 established that gnina must be handed the **post-reaction** ligand,
because gnina replaces an implicit hydrogen rather than performing
reaction chemistry. The transform it introduced strips the leaving
group — and Michael acceptors have no leaving group, so all three
`michael_addition` classes were passed through unchanged, with a note
recorded in `adduct_approximation`.

That note called the error "one hydrogen on a flexible position."
It is not, and the classes do not all share the same error.

## Two chemistries under one mechanism label

| class | product | attachment carbon | transform |
|---|---|---|---|
| `acrylamide` | `Cys-S-CH2-CH2-C(=O)NR2` | sp3 | **saturate the C=C** |
| `naphthoquinone_c2` | 2-thio-1,4-naphthoquinone | sp2 | none |
| `naphthoquinone_benzo` | 2-thio-1,4-naphthoquinone | sp2 | none |

**Quinones re-aromatize.** Thiol addition gives a hydroquinone that
re-oxidizes to the substituted quinone — the basis of quinone protein
arylation, and the isolated product in practice. The sulfur ends up on
an sp2 carbon with the ring alkene intact, so the untransformed molecule
already *is* the product. Saturating it would have modelled a transient
species. The old treatment was right here, for a reason nobody had
written down.

**Acrylamide cannot.** There is no aromatic ring to return to, so the
product is the saturated 3-thio-propanamide. Leaving the alkene in place
handed gnina an sp2 carbon and produced a **vinyl thioether**: planar,
conjugated, and rigid, where the real linker has two freely rotating
bonds. That is not a missing hydrogen — it is the wrong linker
flexibility, and it systematically penalises decorations that must bend
to reach a subsite.

Since acrylamide is **T_3's only warhead**, that bias would have applied
to all of T_3 — and T_3's entire question is which decoration best
complements a fixed covalent scaffold. The approach would have been
answering its own question with an instrument tilted against half the
answers.

## Decision

1. Acrylamide's adduct form saturates the acceptor C=C. Its attachment
   SMARTS becomes `[CH3][CH2][CX3](=O)[NX3]` — the terminal carbon of a
   propanamide, exactly parallel to how the acetamide classes bond the
   CH3 of `[*]C(=O)C`.
2. The quinones are unchanged; their note now states the
   re-aromatization assumption instead of claiming a missing hydrogen.
3. Whether a Michael acceptor saturates is declared in the library
   (`adduct_saturates_alkene`, `warhead_classes_5.csv`), not inferred.
   Both chemistries carry `mechanism: michael_addition` and need
   opposite treatment, so the mechanism label cannot decide it and
   neither can a SMARTS. It is chemistry, and it is written down.
4. The docking run itself moves to `shared/covalent_dock_run.py`. T_3
   and T_4 now execute one function rather than two scripts that began
   identical — the drift risk T_3's own config warns about, applied to
   the loop instead of just the transform.

## Consequences

The protocol fingerprint changes, so **T_4's existing docks are no
longer at parity with T_3's** and the integration GUI will correctly
refuse the within-covalent comparison until T_4 re-docks. T_4's
non-acrylamide classes will reproduce exactly — the docking is seeded
and deterministic and their SMARTS are untouched — so the re-run's real
content is acrylamide's 187 rows.

`05_regiochemistry_comparison`'s naphthoquinone result is **not**
affected. Both of its arms are quinones, neither is transformed, and its
STRONG verdict (benzo over c2 on pose success, p = 1.2e-14) stands. The
likely cause of c2's 96.3% pose failure is that c2 places the R-group
and the Cys113 sulfur on adjacent carbons of one rigid ring while benzo
places the R-group on the distal ring — a real steric difference, not a
modelling artefact.

## Why it was missed

D0022 was written while chasing a leaving group that was still present,
and it partitioned the classes on exactly that question:
`has_leaving_group`. Everything on the false side of that split was
treated as one case. Two of those three classes were fine, which made
the third look fine too.

The recurring lesson from D0025, D0028 and D0029 is that a chemotype
must not be identified by a reactive-atom feature. This is its
complement: **a chemotype must not be identified by the mechanism label
either.** `michael_addition` describes how the bond forms. It says
nothing about what the product looks like once it has.
