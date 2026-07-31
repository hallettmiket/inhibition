---
id: D0046
title: Docking cannot reproduce known Pin1 poses — 5% in production, and the failure is scoring, not sampling
date: 2026-07-31
status: accepted
approach: shared
decided_by: '@mhallet'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - decisions/D0041-the-first-verdict-docking-does-not-demonstrably-enrich.md
  - decisions/D0016-non-covalent-enrichment-is-indistinguishable-from-chance.md
  - shared/noncovalent_dock_run.py
  - approaches/t1_de_novo/03_dock.py
  - approaches/t2_atra_crem/03_dock.py
evidence:
  - 'redocking benchmark over 80 scored X-ray Pin1 ligand complexes, symmetry-corrected RMSD'
  - 'SELF-DOCK into each ligand own receptor, production 26 A box: 13/80 = 16.2% within 2.0 A, median 8.43 A'
  - 'CROSS-DOCK into 6VAJ (the actual T_1/T_2 protocol): 4/80 = 5.0% within 2.0 A, median 7.34 A'
  - 'tight-box control, top-1: 18/80 = 22.5%, median 4.04 A — the box is not the explanation'
  - 'tight-box control, BEST-OF-9: 44/80 = 55.0%, median 1.80 A — the search DOES find the pose'
  - 'self-dock best-of-9: 31/80 = 38.8%; cross-dock best-of-9: 7/80 = 8.8%'
  - 'literature norm for a well-behaved target on self-docking is 60-80% within 2 A'
  - 'CA superposition quality: median 0.38 A, max 0.54 A over PPIase residues 50-163'
  - 'protocol parity by import: Vina-GPU 2.1, SEARCH_DEPTH 20 (D0017), pH 7.4 prep, production run_vina_gpu unmodified'
  - 'RMSD by CalcRMS, NOT GetBestRMS: on a reference translated exactly 3.0 A, CalcRMS returns 3.0000 and GetBestRMS returns 0.0000'
  - '43 covalent complexes excluded and reported separately; 3 entries are 14-3-3 not Pin1; 2 RMSD uncomputable'
  - 'all numbers independently recomputed from redock_per_case_1.csv before acceptance'
---

# The scoring function declines to pick the answer it already found

## The control D0041 never had

D0041 measured non-covalent enrichment on this pocket at ROC-AUC 0.599, CI
[0.311, 0.874], EF1% 0.0 — indistinguishable from chance. It left open *why*,
and the leading hypothesis was decoy construction: property-matched molecules
that were never actually tested, with the AVE/analogue bias that implies.

That hypothesis was never tested against the simpler question. Before asking
"does docking rank binders above non-binders?", one should ask **"can docking
reproduce a binding pose we already know?"** If it cannot, the enrichment null
needs no exotic explanation.

It cannot.

## The numbers

80 scored X-ray Pin1 ligand complexes, symmetry-corrected RMSD, the production
protocol imported rather than reimplemented.

| arm | n | ≤2.0 Å | median RMSD |
|---|---|---|---|
| Self-dock, production 26 Å box | 80 | **16.2%** (13) | 8.43 Å |
| **Cross-dock into 6VAJ — what T_1/T_2 actually run** | 80 | **5.0%** (4) | 7.34 Å |
| Tight-box control, top-1 | 80 | 22.5% (18) | 4.04 Å |
| **Tight-box control, best-of-9** | 80 | **55.0%** (44) | **1.80 Å** |

A well-behaved target gives **60–80%** on self-docking. We get 16%. In
production, into 6VAJ, **5%**.

## The mechanism, which is the actionable part

**Best-of-9 is 55% while top-1 is 22.5%.**

Vina returns nine ranked modes and the pipeline consumes mode 1. The search
*generates* the crystallographic pose in more than half of cases, and the
scoring function then ranks it below a wrong pose roughly three times in four.

The engine finds the right answer and declines to pick it.

Two things follow immediately:

* **The box is not the explanation.** Shrinking from the 26 Å production box to
  a ligand-sized one moves top-1 only 16.2% → 22.5%. Worth having, not a fix.
* **Sampling is not the explanation either.** More search would not help; the
  pose is already in the returned set.

The failure is the ranking function, on this pocket.

## Why cross-docking is worse than self-docking

5.0% vs 16.2%, on the same ligands. Rigid-receptor cross-docking additionally
pays for 6VAJ's side-chain conformations being wrong for other ligands. Pin1's
PPIase site is shallow and solvent-exposed (D0016), which is exactly the
geometry where a single receptor snapshot fails hardest.

**6VAJ is what T_1 and T_2 dock into, so 5% is the number that describes
production**, and 16% is the generous upper bound the engine could reach if
every ligand were given its own cognate structure.

## What this does to D0041

D0041's null stands, and it now has a mechanism. It is much better explained as
**a scoring-function failure on this pocket** than as a decoy-construction
artefact, because the same function fails a strictly *easier* question — not
"rank this binder above other molecules" but "rank this binder's own true pose
above its own wrong poses."

**This makes a prediction, and that is the point of recording it.** Phase 2.1 of
issue #4 re-runs Vina through Gates B and C, which have *experimentally measured*
inactives rather than constructed decoys. If decoy construction were the cause,
that should rescue the enrichment. On this evidence it should not. Phase 2.1 is
now a test of this explanation rather than an open question.

It also raises the priority of anything orthogonal to Vina scoring — the
pharmacophore arm (#4 Phase 2.2) and the receptor-ensemble work (Phase 0.3c),
both of which attack the two failures measured here.

## Method notes that matter

**RMSD is `CalcRMS`, not `GetBestRMS` — and my own brief specified the wrong
one.** Both are symmetry-corrected, but `GetBestRMS` *superposes before
measuring*. On a reference translated by exactly 3.0 Å, `CalcRMS` returns
3.0000 and `GetBestRMS` returns 0.0000. Using `GetBestRMS` would have reported
near-zero for every pose including complete failures, and manufactured a
flawless benchmark. The instruction was overridden and the correct choice
validated against RDKit's own `CalcRMS` (max |diff| 0.0000 Å on 79/82 cases).

This belongs in `docs/how_this_project_breaks.md`: a value chosen because its
*name* sounded right — "best RMS" — rather than because of what it computes.

**Protocol parity by import, not by copying.** Every setting comes from
`shared/noncovalent_dock_run.py`; arm B calls the production `run_vina_gpu`
unmodified. A benchmark run under different settings measures nothing about
this pipeline.

## Honest accounting of what was excluded

From 190 entries and 382 non-polymer components: 253 additives filtered (PEGs,
sulfate, glycerol), 129 kept as ligands, **82 benchmarked, 80 scored**.

* **43 covalent complexes excluded and reported separately.** Redocking a
  covalent adduct non-covalently is a different experiment: the crystal holds
  the *adduct* while we would dock the pre-reaction compound, so the molecular
  graphs differ and a symmetry-corrected RMSD is not defined. Includes 6VAJ's
  own QT7.
* **3 entries are not Pin1** — 7OQ9, 7OQA and 8C3C are 14-3-3 proteins carrying
  a 5-residue Pin1 phosphopeptide; the ligand binds 14-3-3. This corrects the
  "190 Pin1 structures" figure used in issue #4.
* 1 partial density (3I6C), 2 RMSD uncomputable (3TCZ, 3TDB).
* Fragments (<150 Da) reported as their own tier, not pooled: **0%** top-1 at
  the production box.

## What would change this conclusion

* A scoring function that recovers poses here. That is what the pharmacophore
  arm and any rescoring proposal must be measured against **first** — pose
  recovery is a cheaper and stricter test than enrichment, and it is now
  available as a harness.
* Ensemble docking. If cross-docking into an ensemble of the 163 X-ray
  receptors substantially beats 5%, rigid-receptor error is the dominant term
  and the ensemble is the fix.
* An error in the covalent exclusion rule or the RMSD validation. Both should
  be audited by the adversary before this record is cited externally.

## What did not change

No shortlist, rank or candidate. This measures the protocol, not the
candidates. But it does mean that **every docking-derived ranking in this
project is produced by a function that picks the right pose 1 time in 20 on the
receptor it actually uses** — and any statement resting on those rankings
should carry that number.
