---
id: D0025
title: isolate_rgroup ignores its cap argument, producing false aldehyde alerts
date: 2026-07-27
status: proposed
approach: shared
decided_by: '@mhallet'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/alerts.py
  - approaches/t3_reinvent/02_annotate.py
  - approaches/t4_combinatorial/01_enumerate.py
evidence:
  - 'isolate_rgroup(cap="C"), cap="[H]" and cap="CC" all return the identical R-group SMILES — the argument has no effect'
  - 'T_3 candidate C=CC(=O)N(C(=O)NCc1ccccc1OC)C1CCS(=O)(=O)C1 has a UREA decoration; its isolated R-group reads COc1ccccc1CNC=O, a formamide'
  - 'BRENK correctly flags H-C(=O)-N as an aldehyde; the intact molecule has no aldehyde'
  - '2657 of 5396 T_3 candidates were stamped with "aldehyde, aldehyde", the single largest rejection reason'
  - 'two_tier defaults to max_rgroup_alerts=0, so a single artefactual alert rejects'
  - '5270/5396 (97.7%) of T_3 stamped rejected at the alert gate'
runbook: null
---

## Context

T_3 decorates a fixed scaffold at the sulfopin nitrogen, and LibInvent
frequently attaches its decoration through a carbonyl — ureas, carbamates,
amides. That is ordinary medicinal chemistry.

`shared.alerts.isolate_rgroup` cuts the R-group away from the core so alerts are
scored on the decoration rather than on the warhead every candidate is required
to carry. It takes a `cap` argument, documented as the group used to satisfy the
open valence left by the cut, defaulting to `"C"`.

**The argument is ignored.** `cap="C"`, `cap="[H]"` and `cap="CC"` all return
byte-identical output. The cut valence is filled with hydrogen regardless.

For a decoration attached through carbon this is harmless. For one attached
through a carbonyl it is not: an amide `>N–C(=O)–R` becomes, once severed from
the nitrogen and H-capped, a formamide `H–C(=O)–R`. BRENK then flags an aldehyde
that exists only in the fragment, never in the molecule.

## Decision (proposed, NOT yet applied)

Fix `isolate_rgroup` to honour `cap`, and re-run the affected gates.

Not applied tonight, deliberately. `shared/alerts.py` is also T_4's alert gate,
and T_4's enumeration, triage, docking, ranking and MM-GBSA have all been run
against its current behaviour. Changing it invalidates that gate and requires
re-running the chain. That is a considered change with a verification cost, not
a one-line patch to make at the end of a long session — and an unvalidated fix
to a control is worse than a recorded defect.

Two things need deciding together with it:

1. **What the correct cap is.** A methyl keeps a severed amide reading as a
   ketone rather than an aldehyde, which is closer to the truth, but no cap is
   right for every attachment chemistry. The alternative is to screen alerts on
   the INTACT molecule while subtracting those the core contributes — more
   faithful, more work.
2. **Whether `max_rgroup_alerts=0` is the right default.** Rejecting on a single
   alert is strict for a generative approach whose value is proposing chemistry
   nobody would have picked.

## Consequences

**T_3's current annotation (D3_3) is provisional and its 97.7% rejection rate
should not be believed.** The 126 candidates that passed are a biased sample —
they are the ones whose decoration happened not to attach through a carbonyl.

T_4 is affected in principle but far less in practice: its R-groups come from a
frequency-derived ChEMBL library and attach through carbon, not carbonyl. The
count of T_4 candidates rejected specifically for `aldehyde` should be checked
before assuming its gate was unaffected.

The general shape is the same one D0022 had: a control that ran, reported a
plausible number, and was measuring an artefact of how its input was prepared.
Both were found by looking at a specific molecule rather than at a summary
statistic.
