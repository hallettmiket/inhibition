---
id: D0026
title: Excuse named alerts rather than raising the tolerated count
date: 2026-07-28
status: accepted
approach: shared
decided_by: '@mhallet'
origin: user
supersedes: []
superseded_by: null
affects:
  - shared/alerts.py
  - shared/annotate.py
  - config/approaches/t3_reinvent.yaml
  - approaches/t3_reinvent/02_annotate.py
evidence:
  - 'T_3: 3054 candidates carry exactly one attributable alert, and 2853 of those (93%) are acyclic_imide alone'
  - 'raising the count to 1 would also admit thioester, phenol_ester, Thiocarbonyl_group and isolated_alkene one-offs'
  - 'T_3 passing: 1233 (count 0) -> 4086 (count 0, imide excused); rejected 4163 -> 1310'
  - 'T_4 is barely sensitive: 1683 pass at count 0, 1728 at count 1 of 1782'
  - 'the gate only applies to approaches passing core_smarts — T_1 and T_2 have zero alert-gate rejections'
runbook: null
---

## Context

The decoration gate rejected on a single alert, failing 77% of T_3. The natural
response is to raise the tolerated count, and the PI asked whether we could.

The distribution says a count is the wrong instrument. Of the 3,054 T_3
candidates carrying exactly one alert, **2,853 carry `acyclic_imide` and nothing
else**. So "tolerate one alert" would have been a decision to accept imides,
taken by accident, while simultaneously admitting every other one-off alert —
including thioesters and phenol esters, which hydrolyse, and thiocarbonyls,
which are reactive. Those are worse liabilities than the thing the change was
actually meant to permit.

## Decision

**Keep the count at 0. Excuse alerts by NAME, and carry the excusal.**

For T_3, `acyclic_imide` is excused. The scaffold nitrogen already bears the
acrylamide carbonyl, and LibInvent frequently adds a second acyl group, which
makes a genuine imide. This is not a false positive: an N-acyl acrylamide is
more electrophilic and more hydrolytically labile than the acrylamide alone.

It is excused rather than rejected because it is a structural consequence of
*where T_3 decorates*, so gating on it discards 53% of the approach before
anything has been measured — the filter-before-you-measure trap that D0012 and
D0019 both refuse. The alert travels with the candidate as
`excused_alert_names`, and the GUI must display it: an imide shown without that
caveat implies a cleanliness it did not earn.

## Consequences

T_3 goes from 1,233 to 4,086 passing. Every one of the 2,853 imide-bearing
candidates carries a visible flag rather than a silent pass.

This changes nothing for T_1 and T_2, and the record should be explicit because
it was briefly misremembered: the gate lives on the two-tier path, which runs
only for approaches that pass a `core_smarts`. T_1 and T_2 are non-covalent,
have no fixed core, and take the whole-molecule path where nothing is
disqualifying by default. **Neither has ever had an alert-gate rejection.** T_1's
own lever is `--alert-limit`, still unset, so its alerts are annotated only.

T_4 is barely sensitive to the count (1,683 of 1,782 at 0, 1,728 at 1) and gets
no excusals; its decoration alerts are genuine and few.

If docking and MM-GBSA end up favouring the imides, the honest resolution is a
hydrolytic-stability measurement, not a recomputed alert.
