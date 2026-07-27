---
id: D0012
title: Gates report evidence strength; they do not adjudicate
date: 2026-07-27
status: accepted
approach: shared
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/gates.yaml
  - shared/enrichment_gate.py
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - 'covalent stratum has 6 verified actives spanning ~3-4 independent chemotypes'
  - 'ROC-AUC standard error at n=6 is roughly +/-0.2 — a 0.70 PASS sits within noise of FAIL'
  - 'EF1% over 306 molecules resolves to the top 3 compounds, so it is quantised, not continuous'
  - 'the choreography already refuses an authoritative cross-approach numeric join (Rev 3 section 7)'
runbook: null
---

## Context

The enrichment gate was specified as a binary PASS/FAIL on ROC-AUC, EF1% and
BEDROC thresholds. Assessing the available actives showed those thresholds
cannot carry that weight: six actives, three or four independent chemotypes, and
an EF1% that is quantised to a handful of values.

The tempting responses were both wrong. Emitting a confident PASS from six
actives manufactures precision that is not there. Emitting FAIL, or refusing to
run until the statistics are strong, discards real signal — and the statistics
will *never* be strong here, because validated Pin1 chemistry is genuinely
scarce. That scarcity is a finding about the target, not a defect to engineer
around.

The PI's framing settles it: **this is not purely a statistical exercise.** The
choreography exists so PIs and researchers can bring their own priors and
expertise to bear on the ranking. False discovery is expected and accepted; the
job is to compute where the most promising leads are, not to prove significance.

## Decision

Gates **report evidence strength; they do not adjudicate.** Concretely:

- The enrichment gate emits a graded verdict — `STRONG`, `WEAK`,
  `UNDERPOWERED`, `FAIL` — not a binary, and always alongside confidence
  intervals, the actives count, and the **independent chemotype count**.
- Evaluation is **per chemotype**, leave-one-chemotype-out, so analog bias
  cannot inflate a result. Six actives that are three chemotypes get reported as
  three.
- `UNDERPOWERED` does **not** mean "discard". It means the ranking is carried
  forward with its uncertainty displayed, for a human to weigh. It is a label on
  the evidence, not a veto.
- Only `FAIL` — docking demonstrably anti-correlated with known actives —
  demotes `dock_score` to a displayed label.

This makes the gate consistent with the choreography's existing stance in Rev 3
section 7: present the evidence and its limits, let the human adjudicate. A gate
that silently vetoed on thin statistics would be the one component doing the
opposite of what everything else does.

## Consequences

The GUI must display gate verdicts with their power characteristics, not just a
green tick — a `WEAK` verdict shown as a pass would be worse than no gate. The
Open Questions panel already surfaces this class of limitation.

D0011 still stands: M3 decides the covalent rank metric empirically. But if the
comparison between CNNaffinity and Vina-style affinity comes back
`UNDERPOWERED`, that is a real answer — it means the choice should be made on
mechanistic grounds and the tool's own calibration warning, not on six data
points.
