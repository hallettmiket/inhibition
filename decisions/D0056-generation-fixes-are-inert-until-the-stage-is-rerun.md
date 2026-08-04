---
id: D0056
title: A fix in a generation stage is inert until the stage is re-run, and nothing announced it
date: 2026-08-04
status: accepted
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - scripts/check_frame_code_currency.py
  - tests/test_frame_code_currency.py
  - approaches/t3_reinvent/01_generate.py
evidence:
  - 'b776843 (2026-07-30) added the pocket size ceiling for T_2 and T_3; T_3 frame D3_2 was written 2026-07-27 23:20 and T_3 generation has never been re-run'
  - 'no T_3 frame (26 of them) has ever carried the heavy_atoms column that verify() computes'
  - 'measured: 2 of 5,396 T_3 molecules exceed 55 heavy atoms (max 59), both already rejected by alerts, neither on shortlist or shortlist_synth'
  - 'T_2 ATRA generated 2026-07-27 (before the ceiling); the four newer seeds 2026-07-31 (after)'
  - 'measured: no T_2 pool exceeds 49 heavy atoms (atra max 30) — the ceiling would never have fired'
  - 'all four generation frames record git.dirty = true, so their commit does not describe the code that ran'
  - 'measured: 0 candidate_id collisions across nine pools, 72,104 molecules, injective id->SMILES'
---

# The fix is in the code and not in the data

## What was found

Auditing the generation stages, two stages turned out to have been last run
BEFORE a fix that changes what they produce:

* **T_3's pocket size ceiling.** `verify()` computes `heavy_atoms` and stamps
  `exceeds_pocket_ceiling`. It was added on 2026-07-30. T_3's production frame
  `D3_2` was written on 2026-07-27 and T_3 has never been re-generated, so
  **no T_3 frame has ever carried that column** — across all 26 of them.
* **T_2 ATRA**, from the same commit. ATRA predates the ceiling; the four newer
  seeds postdate it.

Neither is visible from the code, which reads as though the ceiling is in
force, or from the frame, which is populated and plausible.

## The measured impact is small. The class is not.

Being honest about magnitude matters here, because the instinct on finding this
is to overstate it:

| | measured |
|---|---|
| T_3 molecules over the 55-atom ceiling | **2 of 5,396** (max 59) — both already rejected by `alerts`, on neither shortlist |
| T_2 molecules over the ceiling, any pool | **0** — the largest pool tops out at 49 heavy atoms |

So the ceiling would have changed two stamps on two already-excluded molecules,
and nothing else. **The defect is real and its consequences here are nearly
nil.** The reason to record it is that nothing detected it, and the next
instance may not be harmless.

**A related observation worth separating out:** the ceiling as configured (55,
the pocket limit) was never going to address the problem D0043 identified. D0043
is about T_3 shortlisting molecules at ~2x the generated median — a size *drift*
well inside the pocket limit, not molecules that cannot fit. The ceiling removes
what is physically impossible; it is not a lever on the drift. (Separately
measured 2026-08-04: T_3's shortlist median is now **31** heavy atoms against a
generated median of 25, where D0043 reported 39 vs 25 — the gap has narrowed
from 14 to 6 since D0047/D0049 landed.)

## The decision

`scripts/check_frame_code_currency.py` compares each generation stage's newest
frame — via the commit its manifest records — against later commits touching
that stage's source. `tests/test_frame_code_currency.py` fails the suite on
staleness nobody has looked at.

**It reports rather than simply failing.** Staleness is frequently the correct
state: re-running DiffSBDD or a 16,806-molecule CReM expansion to pick up a
comment change would be absurd. So known-stale stages are ACKNOWLEDGED with a
reason *and a measured impact*, and a second test asserts every acknowledgement
contains the word MEASURED and a number — because "we checked, it's fine" is
worth nothing without the figure attached. Allowlist, not denylist (D0051): a
stage nobody has looked at fails rather than passing unnoticed.

## What it found immediately

Two more stages nobody had examined — and **both dissolved on measurement**:

* **T_1** — `1c40ba4` (floor 10 + `size_class`) postdates `D1_3`'s commit, but
  the change is demonstrably in the data: `size_class` is on the frame, 1,376
  rows stamped `degenerate_too_small`, 82 `too_large`.
* **T_4** — `72cf331` (dangling attachment points, colliding candidate ids)
  postdates `D4_6`, but 0 of 1,782 SMILES carry an unfilled attachment point, 0
  fail to parse, and there are 0 candidate_id collisions across all nine pools.

## The caveat that makes this honest

**All four generation frames record `git.dirty = true`**, so the commit each one
names does not describe the code that actually ran, and the comparison
systematically OVER-reports. That is the right direction to be wrong in — it
asks a question rather than granting a pass — but it means a STALE verdict is a
prompt to measure, never a finding. Two of four flags dissolved on measurement.
Read the impact line, not the tag.
