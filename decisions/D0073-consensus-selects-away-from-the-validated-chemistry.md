---
id: D0073
title: Consensus does not enrich for validated mechanisms — it depletes them, and the best-validated class has the worst pass rate
date: 2026-08-06
status: proposed
approach: shared
decided_by: null
review_requested_from: '@tt8804'
origin: adversary
supersedes: []
superseded_by: null
affects:
  - shared/pose_consensus.py
  - scripts/nac_consensus_all.py
  - docs/recap_2.0.0.md
  - decisions/D0070-consensus-preserves-rank-where-frequency-does-not.md
  - decisions/D0065-warhead-presentation-geometry-ranks-covalent-candidates.md
evidence:
  - '5,765 screened molecules with usable consensus; 397 clear the >= 0.90 single-mode bar (6.9% pass rate)'
  - 'validated-mechanism share FALLS across the filter: library 90.3% -> pool 77.8%. Fisher odds ratio 0.34, p = 2.29e-14'
  - 'pass rate by mechanism: sn2_ring_opening (BDHI) 16.6%, snar_displacement 13.9%, michael_addition 6.3%, sn2_displacement (chloroacetamide) 2.9%'
  - 'share lift (pool share / library share): BDHI 2.41x, SNAr 2.02x, Michael 0.92x, SN2-displacement 0.41x'
  - 'the two mechanisms consensus favours are exactly the two with no usable validation: BDHI has 0 crystallographic positives, SNAr has 2 and failed at AUC 0.558 (p = 0.41)'
  - 'the mechanism it most disfavours is the best-validated one: chloroacetamide, 9 positives, AUC 0.908 [0.857, 0.954]'
  - 'largely rigidity: pass rate is monotone in rotatable bonds — rotb 0-2 19.5%, 3-4 10.1%, 5-6 4.9%, 7+ 2.9%; rho(rotb, consensus) = -0.259, p = 3.8e-89'
  - 'median rotatable bonds by mechanism: BDHI 4.0, Michael 5.0, SNAr 5.0, chloroacetamide 7.0 — the best-validated class is the most flexible in this library'
  - 'rigidity is NOT the whole story: SNAr and Michael share a median rotb of 5.0 but pass at 13.9% vs 6.3%'
  - 'only 16 of the 397 are SN2-displacement (11 sulfamate acetamide, 5 chloroacetamide)'
runbook: null
---

# Consensus selects away from the chemistry we can validate

**Flagged for @tt8804 to review before it shapes the 2.1.0 ranking.**

## What prompted it

The observation on the table was that filtering 5,769 → 397 on consensus leaves
78% of the pool in a mechanism that cleared validation, which looks like the
filter doing useful work. It is worth checking whether that 78% is *selection* or
just *composition* — and it is composition, with the sign running the wrong way.

## The finding

**The library is 90.3% validated-mechanism. The consensus pool is 77.8%.**
Consensus does not enrich for validated chemistry; it **depletes** it, at odds
ratio 0.34 (p = 2.3×10⁻¹⁴).

| mechanism | library | pool | pass rate | lift | validation |
|---|---:|---:|---:|---:|---|
| SN2 ring-opening (BDHI) | 6.5% | 15.6% | **16.6%** | **2.41×** | **0 positives, never tested** |
| SNAr | 3.2% | 6.5% | 13.9% | 2.02× | 2 positives, **AUC 0.558, failed** |
| Michael addition | 80.5% | 73.8% | 6.3% | 0.92× | 4 positives, AUC 0.734 |
| SN2 displacement (chloroacetamide) | 9.7% | 4.0% | **2.9%** | **0.41×** | 9 positives, **AUC 0.908** |

The ordering is close to exactly inverted against the evidence: **the two
mechanisms consensus favours are the two we cannot validate, and the one it
penalises hardest is the one we validate best.**

## Why

Mostly the rigidity confound, now visible at library scale. Pass rate is monotone
in rotatable-bond count:

| rotatable bonds | n | pass rate |
|---|---:|---:|
| 0–2 | 41 | 19.5% |
| 3–4 | 2,435 | 10.1% |
| 5–6 | 2,346 | 4.9% |
| 7+ | 943 | 2.9% |

ρ(rotb, consensus) = −0.259, p = 3.8×10⁻⁸⁹. A rigid molecule has fewer ways to
sit, so its top poses agree — which is a statement about the molecule's
conformational freedom, not about the pocket.

The warhead classes differ systematically in flexibility: BDHI is a fused ring
system (median rotb 4.0) while chloroacetamide is an acyclic chain (median 7.0).
So a filter that rewards rigidity rewards BDHI and punishes chloroacetamide, and
the mechanism skew follows.

**Rigidity is not the whole explanation.** SNAr and Michael share a median rotb
of 5.0 yet pass at 13.9% vs 6.3%. Something mechanism-specific is left over and
is not yet accounted for.

## Why this matters for 2.1.0

The plan on the table was to build the shortlist on consensus plus the validated
class. **These two criteria pull against each other**, and the size of the
conflict was not visible before: consensus leaves only **16 SN2-displacement
molecules** in the whole pool, 5 of them chloroacetamide.

This does not make consensus wrong. D0070 stands — it is rank-stable where
enrichment is not, and that is a real property. But it does mean:

1. **A consensus cut applied library-wide will be dominated by chemistry we
   cannot defend.** Two-thirds of the classes it promotes have no working
   criterion behind them.
2. **Consensus should be applied within mechanism, not across it.** Comparing a
   BDHI's consensus to a chloroacetamide's compares their ring counts as much as
   their poses. Per-mechanism ranking, or explicit rigidity correction, is the
   minimum fix.
3. **Any "top N overall" list is a rigidity ranking wearing a geometry label.**
   That is the same defect class as the rest of 2.0.0 — a value taken by a
   correlate rather than by the thing itself.

## What this does not say

- It does **not** say consensus is measuring nothing. Rigid molecules genuinely
  do have better-determined poses; the problem is that the property is confounded
  with warhead chemistry in *this* library.
- It does **not** say BDHI is bad chemistry. It says BDHI is unvalidated, and
  that a filter preferring it moves the shortlist further from evidence.
- It does **not** establish the right correction. Rotatable-bond matching worked
  for the elevation cohort at n = 8 per group; whether it is the right instrument
  at library scale is untested.

## Open

- What is the residual mechanism effect after rigidity is accounted for?
  (SNAr 13.9% vs Michael 6.3% at equal median rotb.)
- Should the consensus bar be set **per mechanism** at a fixed pass rate rather
  than at a fixed 0.90 threshold?
- Does the chloroacetamide arm survive at all if only 5 molecules clear the bar,
  and does that argue for generating a more rigid chloroacetamide series rather
  than lowering the bar?
