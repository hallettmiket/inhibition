---
id: D0011
title: gnina CNN scoring is uncalibrated for covalent docking
date: 2026-07-27
status: proposed
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/covalent_protocol.py
  - config/choreography.yaml
  - docs/approaches/t3.md
  - docs/approaches/t4.md
evidence:
  - 'gnina v1.3.3 prints: "CNN scoring not yet calibrated for covalent docking. Recommend running with --cnn_scoring none"'
  - 'the Rev 3 spec makes gnina CNNaffinity T_3 rank metric and T_4 secondary metric'
  - 'Sulfopin covalent dock: CNNaffinity 4.64, CNNscore 0.6286, Vina-style affinity -2.35 kcal/mol'
  - 'the warning is emitted on every covalent run, not just edge cases'
runbook: null
---

## Context

M2 pinned the shared gnina covalent protocol and ran it end-to-end on Sulfopin,
the anchor whose true covalent pose 6VAJ resolves. It works. But gnina prints,
on every covalent run:

> CNN scoring not yet calibrated for covalent docking. Recommend running with
> `--cnn_scoring none`

The Rev 3 spec makes `CNNaffinity` **T_3's rank metric** and T_4's secondary
metric, and §7 offers an optional within-covalent-stratum re-score built on it.
The tool's own authors are saying that number is not calibrated for the mode we
are using it in.

This is the same shape as the LUMO finding in D0005: a metric the spec treats as
a ranking signal, which the evidence says is not one.

## Decision

**PROPOSED, not yet accepted — this changes a spec-level metric choice and
wants the PI's call.** Options, in the order I would rank them:

1. **Rank covalent candidates by gnina's Vina-style `affinity` (kcal/mol,
   lower better) and carry `CNNaffinity` as an advisory annotation.** Keeps the
   shared-protocol parity (S3) intact, uses a metric that is at least
   calibrated for what it measures, and costs nothing to adopt.
2. **Keep CNNaffinity but run the enrichment gate (M3) on the covalent stratum
   specifically**, and let the measured result decide. This is what the gate
   exists for, and it turns the question into an empirical one.
3. Keep CNNaffinity as the rank metric and note the caveat. Weakest option: it
   ranks on a number the tool says is uncalibrated.

Option 2 subsumes option 1 and is the honest route, at the cost of needing M3
before T_3 can rank anything.

Interim, already implemented: `dock()` detects the warning, logs it, and
returns `cnn_uncalibrated_for_covalent` so it reaches the manifest and the GUI
rather than scrolling past in a log.

## Consequences

If CNNaffinity is demoted, T_3's output contract changes (rank metric, and its
stated direction flips from higher-better to lower-better), and §7's
within-covalent re-score changes with it. The docs state the direction
explicitly in several places and would need updating together.

Nothing is blocked meanwhile: the protocol is pinned, parity holds, and both
numbers are recorded per dock regardless of which one ends up ranking.
