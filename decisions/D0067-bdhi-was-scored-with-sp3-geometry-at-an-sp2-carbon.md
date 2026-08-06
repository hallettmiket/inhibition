---
id: D0067
title: The BDHI classes were scored with sp3 backside geometry at an sp2 carbon, and read as unreactive
date: 2026-08-06
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/nac_criterion.py
  - scripts/nac_rank.py
evidence:
  - 'bdhi_c4 and bdhi_c5 are the only members of mechanism sn2_ring_opening; reactive_atom_smarts [CX3]([Br])=[NX2]'
  - 'RDKit on the class fragment [*]C1C(Br)=NOC1: attacked atom is C, hybridisation SP2, degree 3, non-aromatic'
  - 'production ranking, 374 BDHI candidates: median enrichment 0.00x, 67% (bdhi_c4) and 96% (bdhi_c5) scoring exactly zero'
  - 'their MEASURED median S-C-Br angles were 91.8 and 110.5 degrees -- near-perpendicular, and the SN2 window demands >=150'
  - 'median distances 3.67 and 3.64 A, inside the 2.8-4.2 A near-attack window -- the failure was purely angular'
  - 'the two VALIDATED classes are unaffected: chloroacetamide is sp3 (backside, unchanged), michael_addition was already perpendicular'
---

# BDHI was asked the wrong geometric question

## What happened

`shared/nac_criterion.py` maps a warhead's **mechanism** to a geometry.
`sn2_ring_opening` was mapped to backside attack anti to the leaving group — on
the strength of the name.

Its only members are `bdhi_c4` and `bdhi_c5`, the 3-bromo-4,5-dihydroisoxazoles.
Their attacked atom is the carbon of a **C=N**, which RDKit confirms is **sp2**,
degree 3. A thiolate does not displace bromide there by attacking from behind:
it adds perpendicular to the C=N plane and bromide leaves. That is
addition–elimination, and the ring-opening name describes the *outcome* rather
than the *trajectory*.

## How it presented

Exactly as a real chemical result would.

| | bdhi_c4 | bdhi_c5 |
|---|---|---|
| median enrichment | **0.00×** | **0.00×** |
| scoring exactly zero | 67% | 96% |
| median S–C–Br angle | **91.8°** | **110.5°** |
| median S···C distance | 3.67 Å | 3.64 Å |

The distances sit comfortably inside the 2.8–4.2 Å near-attack window, so nothing
was failing to reach the sulfur. **The search was finding near-perpendicular
approaches — the correct geometry for an sp2 centre — and the criterion was
scoring them dead because it wanted ≥150°.**

Left alone, this would have ranked 374 candidates last, across two of the nine
warhead classes, and the reason would have read as "BDHI cannot present its
warhead" rather than "we measured the wrong angle."

## The pattern, again

This is D0064's observation arriving a third time: *90° is dead for SN2 and ideal
for SNAr*. BDHI is a third sp2 case, and it was mislabelled because a mechanism
NAME was trusted over the attacked atom's hybridisation.

The general form is the project's signature defect in a new dress — a value taken
by **label** rather than by identity, producing a populated, plausible, wrong
answer that raises no exception.

## The fix

`sn2_ring_opening` now maps to `perpendicular_to_plane`. The class's SMARTS
`[CX3]([Br])=[NX2]` supplies exactly the three atoms a plane needs (C, Br, N).

**This is a defect fix and not a tuning.** Applying sp3 geometry to an sp2 centre
is wrong independent of what any data shows; the measured 91.8°/110.5° confirmed
the diagnosis but the chemistry decided it. Crucially it **does not touch either
validated class** — chloroacetamide is genuinely sp3 and keeps backside geometry,
`michael_addition` was already perpendicular — so D0065's AUCs stand unchanged.

A genuinely sp3 ring-opening (an epoxide, an aziridine) would need its own
mechanism label, because it would want the backside geometry this entry no longer
provides. Noted in the source so the next person adding one does not inherit this
mistake inverted.

## Consequence for the data on disk

The 374 BDHI rows already written are invalid. `append_only` forbids rewriting
them, so they stay, and `nac_rank.load_scored()` now resolves duplicates with
`keep="last"` in **version order** — the project's usual convention that the
highest integer is current. `--redo-classes` re-scores named classes ignoring the
done-set, and the new rows supersede.

Version-order sorting matters here and is not incidental: `nac_rank_s0_10.csv`
sorts before `nac_rank_s0_2.csv` lexically, which would have made "newest" wrong
in precisely the case this mechanism exists to handle.
