---
id: D0019
title: Reactivity window flags classes, it does not reject them
date: 2026-07-27
status: accepted
approach: t4
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - approaches/t4_combinatorial/02_reactivity_triage.py
  - config/approaches/t4_combinatorial.yaml
  - integration/app/DECISIONS_TAB_SPEC.md
evidence:
  - 'window with ALL anchors spans 3.96 eV and excludes ZERO of 9 classes — a filter that does not filter'
  - 'dropping the 2 promiscuous quinone anchors narrows it to 1.58 eV and excludes 4 of 9'
  - 'adding the 8 measured-kinetics compounds barely moves it (1.61 eV, same 4 excluded)'
  - 'excluded: acrylamide, naphthoquinone_c2, naphthoquinone_benzo, snar_chloroazine'
  - 'acrylamide misses by 0.106 eV against a hand-chosen 0.5 eV tolerance'
  - 'every clean anchor and every kinetics compound is chloroacetamide/sulfamate/sulfonate — ONE chemical family'
  - 'D0005: computed reactivity vs measured Pin1 labelling correlate at only r = 0.396'
runbook: null
---

## Context

The reactivity window as first built excluded nothing. Anchored on all six
verified covalent actives, it spanned ~4 eV — bounded at one end by mild
chloroacetamides and at the other by highly reactive quinones — and every one of
the nine warhead classes fell inside. A filter that does not filter, reported as
though everything had passed a safety check.

Two corrections were approved: drop the promiscuous quinone anchors (juglone and
KPT-6566 both carry `promiscuity_flag = y`, so bounding a SAFETY window with them
admits reactivity nobody would accept in a lead), and calibrate against the
measured rate constants rather than computed LUMO alone.

Both were applied. The window narrows from 3.96 eV to 1.61 eV and starts
discriminating — and the measured kinetics barely move the bounds set by the
clean anchors, which is reassuring: computed and measured chemistry agree on
where the range sits.

But it now excludes four of nine classes, **including acrylamide — the warhead
chosen by the PI for T_3**.

## Decision

**The window FLAGS classes; it does not reject them.** Candidates outside the
window carry `reactivity_flag = OUTSIDE_WINDOW` and proceed to docking with
their evidence collected. `rejected_at` is not set.

Three reasons the exclusion is not strong enough to act on as a veto:

1. **The window is chemotype-narrow by construction.** Every clean anchor and
   every kinetics compound is chloroacetamide, sulfamate or sulfonate — one
   chemical family. A window built from one family will exclude other families
   almost by definition, whether or not they are genuinely unsafe. This is the
   "chloroacetamide-centric" caveat from the reference provenance, now biting in
   the opposite direction from before.
2. **LUMO is a weak proxy and we measured how weak.** D0005 found computed
   reactivity and Pin1 labelling correlate at r = 0.396. The window answers
   condition (ii) — is this electrophile in a precedented safety range — not
   whether the warhead will work.
3. **The margin is 0.106 eV**, against a 0.5 eV tolerance I chose myself. That
   is not a robust exclusion, and widening the tolerance after seeing which
   class it would readmit would be choosing the answer.

This is D0012 applied to a second gate: report evidence strength, let the human
adjudicate. The docking data for the flagged classes will exist regardless, so
the flag costs nothing and preserves the option.

## Consequences

T_4 keeps all nine classes through docking. The GUI must display
`reactivity_flag` beside those candidates — a flagged candidate presented
without it would imply a safety assessment it did not pass.

Note this does not bear on T_3 directly: T_3's condition (ii) rests on the
up-front expert warhead choice, not on this window. It is informative about that
choice, not a contradiction of it.

If a flagged class subsequently ranks well on docking and MM-GBSA, that is
exactly the case worth a chemist's attention — a warhead outside the
precedented range but performing — and the honest resolution is measured
kinetics on the specific compound, not a recomputed LUMO.
