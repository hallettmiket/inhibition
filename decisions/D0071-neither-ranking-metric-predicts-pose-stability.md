---
id: D0071
title: Neither enrichment nor consensus predicts pose stability, and the crystallographic anchor is what makes that readable
date: 2026-08-06
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/elevation_prereg.md
  - docs/elevation_results.md
  - scripts/elevation_run.py
  - scripts/elevation_analysis.py
  - shared/composite_rank.py
  - decisions/D0068-enrichment-depends-on-search-effort-and-energy-beats-it-at-convergence.md
  - decisions/D0070-consensus-preserves-rank-where-frequency-does-not.md
evidence:
  - 'PRE-REGISTERED before any simulation ran: docs/elevation_prereg.md, groups + readouts + readings all fixed in advance'
  - '37 molecules x 3 replicas, 111/111 tier-1 replicas succeeded, 0 failures'
  - 'tier 1 (|delta d| across 300 ps unrestrained equilibration, nm): A 0.277, B 0.198, D 0.204, V 0.203, REF 0.102'
  - 'tier-1 pre-registered contrasts all null: A-B p=0.130 (Holm 0.391, delta -0.469); B-D p=0.878 (delta -0.062); A-D p=0.442 (delta -0.250)'
  - 'tier-1 anchor contrasts all significant: A-REF p=0.0070 delta -0.781; B-REF p=0.0104 delta -0.750; D-REF p=0.0499 delta -0.594'
  - 'still in the 0.28-0.42 nm near-attack window at the start of production: A 0.08, B 0.21, D 0.08, V 0.53, REF 0.54'
  - 'POST-HOC: drift under 300 ps of dynamics is REF 0.049 nm vs A 0.226, B 0.159, D 0.174, V 0.173; minimisation contributes 0.018-0.055 nm in every group'
  - 'signed drift positive in 110 of 111 replicas; all 37 molecules have a positive mean'
  - 'group V is n = 5 and carries no significance claim, per the pre-registration'
runbook: null
---

# Neither metric predicts stability; the anchor is what makes that readable

## Context

Two metrics rank the same 5,769 candidates and disagree almost completely.
**Enrichment** — the fraction of docking runs reaching a near-attack
conformation — does not converge: D0068 measured the same molecules falling from
2.91× to 0.96× at 10× the search effort, with rank correlation ρ = −0.117 across
efforts. **Consensus** — whether a molecule's top-10 poses by energy agree — is
rank-stable at ρ = +0.568 (D0070). 397 molecules have a single binding mode at
≥ 0.90 agreement and only 4 of them clear enrichment > 5.70, so the enrichment
cut captures 1% of the well-aligned molecules and misses five chloroacetamides —
the one warhead class with a validated criterion.

The question is which of the two, if either, selects for something physical.
Because four groups and several possible readouts can be made to say almost
anything after the fact — and D0045 exists in this project because that has
happened — **the design, the readouts and the readings were written and
committed before any simulation ran** (`docs/elevation_prereg.md`).

## Decision

**Neither metric predicts pose stability.** The pre-registered reading is
**A ≈ B ≈ D**: all three pre-registered contrasts are null, with Holm-corrected
p of 0.391, 0.884 and 0.884. Its fixed conclusion stands — the ranking has no
physical support from this experiment, and no shortlist may be described as
selecting for stability on the strength of either metric.

**The anchor is what makes that a result rather than a dead assay.** All three
BDHI groups are significantly *less* stable than the 8 crystallographic Cys113
positives (p = 0.007, 0.010, 0.050; Cliff's δ = −0.78, −0.75, −0.59). The
measurement separates molecules known to react from generated candidates at
every group. So the null is a statement about the metrics, not about the
measurement — and the prereg's own wording, "BPMD is measuring something
orthogonal", is **too generous to the metrics and is amended here**: tier 1
measures pose survival, which is what it was built to measure.

The three point estimates all run in the direction of the *higher*-enrichment
group being the *less* stable one. **No claim is made from that**, since none is
significant; it is recorded so that a later run with more power knows which
direction to look.

## Consequences

- **Neither metric earns a place in the ranking on stability grounds.**
  Consensus keeps whatever independent support D0070 gave it for rank
  *stability across search effort*, which is a different property and is not
  affected by this result. Enrichment gains nothing here.
- **The screen still does its job as a filter.** Nothing in this experiment
  bears on concentration of active-like molecules, only on ordering.
- **A cohort where some group survived would test the metrics harder.** Every
  BDHI group's poses left the near-attack window; the contrasts are being drawn
  between degrees of failure. That is a limit on power, not a defect in the
  design, and it was visible only after the run.
- **Group V and the anchor agree on one post-hoc readout and not on the
  pre-registered one.** On "still in the near-attack window" V is 0.53 against
  REF's 0.54; on |Δd| V is 0.203 against REF's 0.102. Both are reported in
  `docs/elevation_results.md` and neither is allowed to stand in for the other.
  The prereg's *V ≈ REF* reading is therefore **not** claimed.
- **What would change this:** more power on the A-vs-B contrast (δ = −0.469 is a
  moderate effect that n = 8 cannot resolve), or a cohort drawn so that at least
  one group's poses survive.
