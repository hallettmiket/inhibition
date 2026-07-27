---
id: D0016
title: Non-covalent docking barely enriches on Pin1
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - docs/approaches/t1.md
  - docs/approaches/t2.md
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - 'non-covalent Vina: ROC-AUC 0.535, CI [0.215, 0.855], EF1% 0.0, BEDROC 0.083'
  - 'EF1% 0.0 — no known Pin1 binder reached the top 1% of the ranking'
  - '5 actives, 243 decoys, 5 independent chemotypes'
  - 'the CI spans from clearly-worse-than-chance to strongly-enriching'
  - 'Rev 3 section 3 already flags the PPIase pocket as shallow and solvent-exposed'
runbook: null
---

## Context

The enrichment gate ran the non-covalent stratum through the same Vina protocol
T_1 and T_2 will use, against property-matched decoys on the prepared 6VAJ
receptor.

## Decision

Record that **non-covalent docking shows essentially no enrichment on this
target**, and treat T_1 and T_2 dock-based rankings as weakly supported until
more actives are available.

ROC-AUC 0.535 is a coin flip. The interval [0.215, 0.855] is so wide it cannot
distinguish "docking works" from "docking is actively misleading". EF1% 0.0 and
BEDROC 0.083 agree: no known binder reached the top of the ranking.

This does NOT trip FAIL, which requires the interval's upper bound to fall below
0.5, so under D0012 the ranking carries forward with its uncertainty displayed
rather than being vetoed.

## Consequences

**T_1 and T_2 both rank on Vina.** On this evidence their dock-based shortlists
should not be presented as evidence-backed, and the GUI must show the gate
verdict beside them rather than a bare score. The covalent stratum is in
markedly better shape (AUC 0.815, CI excluding 0.5), so the two families are not
equally supported and should not be displayed as if they were.

This is a plausible property of the target rather than a defect in the setup.
Rev 3 section 3 already notes Pin1's PPIase pocket is shallow and
solvent-exposed — the regime where structure-based methods are weakest, and the
reason T_1's sanitise/filter stages were called load-bearing.

Before concluding docking cannot rank here, the honest next step is more
actives: 5 non-covalent actives over 5 chemotypes cannot separate a real null
from a small effect. If the result survives an expanded set, T_1 and T_2 should
lean on their other evidence and the human, exactly as the FAIL branch intends.
