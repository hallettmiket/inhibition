---
id: D0025
title: isolate_rgroup ignores its cap argument, producing false aldehyde alerts
date: 2026-07-27
status: accepted
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

## Decision (APPLIED 2026-07-28 — option 3, PI's choice)

**Do not cut the molecule at all.** Screen it intact and attribute each alert by
the atoms RDKit reports it matched: wholly inside the excused region it is the
mechanism, wholly outside it belongs to the decoration, and a match spanning
both is charged to the decoration but reported separately as a boundary hit.

`isolate_rgroup` is left in place, its result retained for human inspection and
for nothing that decides anything. Its `cap` argument is documented as ignored
rather than fixed, because the function is no longer load-bearing.

TWO REGRESSIONS FOUND WHILE APPLYING THIS, both reintroducing the exact false
positive the two-tier design exists to prevent:

1. **Attributing by the core alone rejected all 1,683 T_4 survivors.** T_4 scopes
   against `N[CH]1CCS(=O)(=O)C1` — the sulfolane and its nitrogen — while the
   warhead hangs outside that pattern, so every warhead alert was charged to the
   decoration. The warhead must be excused alongside the core.
2. **Excusing only the reactive ATOMS still rejected 1,146.** chloroacetamide's
   reactive-atom SMARTS is `[CH2][Cl]`, so the carbonyl stayed exposed and
   `alpha_halo_carbonyl` straddled the boundary. The excused region is defined
   by `warhead_fragment_smiles`, the whole group.

Both were caught by diffing the new gate against the old on all 1,782 T_4
candidates rather than trusting it. Final effect on T_4: 1,683 -> 1,624 passing,
59 changed, all newly rejected on genuine decoration alerts.

A THIRD BUG, found in the same pass: re-running annotation over its own output
inherited the previous run's rejections, because `screen_frame` only stamps rows
whose `rejected_at` is still null. 1,114 T_3 rows were stamped `alerts` while
carrying `alert_gate_pass = True` — a frame disagreeing with itself. The stage
now clears the stamps it owns before re-applying them.

The original note, retained: `shared/alerts.py` is also T_4's alert gate,
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

**T_3's rejection rate barely moved — 97.7% to 77.2% — and that is the
interesting part.** The aldehyde and thiol artefacts are gone entirely, but the
dominant alert is now `acyclic_imide` (3,787 candidates), and it is REAL. The
T_3 scaffold nitrogen already carries the acrylamide carbonyl; LibInvent
frequently decorates it with a second acyl group, and an N with two acyl groups
is an acyclic imide — genuinely more electrophilic and hydrolytically labile
than the acrylamide alone.

So the finding is not "the gate was broken", it is **"LibInvent's preferred
decoration chemistry for this scaffold creates imides"**. That is a fact about
T_3 worth carrying to the panel, and it would have stayed hidden underneath the
aldehyde artefact. 1,233 of 5,396 now pass.

`max_rgroup_alerts = 0` remains the default and remains worth revisiting
separately: rejecting on a single alert is strict for an approach whose value is
proposing chemistry a person would not have picked.

T_4 is affected in principle but far less in practice: its R-groups come from a
frequency-derived ChEMBL library and attach through carbon, not carbonyl. The
count of T_4 candidates rejected specifically for `aldehyde` should be checked
before assuming its gate was unaffected.

The general shape is the same one D0022 had: a control that ran, reported a
plausible number, and was measuring an artefact of how its input was prepared.
Both were found by looking at a specific molecule rather than at a summary
statistic.
