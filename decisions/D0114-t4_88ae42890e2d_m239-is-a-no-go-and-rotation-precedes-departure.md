---
id: D0114
title: t4_88ae42890e2d_m239 is a NO GO — the pocket holds it while the warhead rotates away, and the rotation leads the departure by ~20 ns
date: 2026-09-06
status: accepted
approach: shared
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - scripts/attack_sweep.py
  - docs/state_of_the_project.md
evidence:
  - '100 ns explicit-solvent run: mean ligand RMSD 0.352 nm, final 0.424, max 0.744, PBC guard NOT triggered -- the ligand never left the pocket'
  - 'warhead-to-Cys113 SG over the same run: engaged (2.8-3.5 A) in only 8.0% of frames, median 4.46 A, max 11.10 A'
  - 'by 10 ns block the distance goes 3.69 / 3.67 / 3.76 / 4.01 / 3.88 / 5.65 / 8.83 / 7.70 / 6.55 / 7.15 A -- departure at 50-60 ns and no return'
  - 'the OFF-NORMAL ANGLE degrades first: median 17.6 / 18.2 / 18.0 / 28.9 / 27.8 / 53.9 / 64.7 deg. The fraction within 30 deg falls 95% -> 55% at 30-40 ns while the distance is still ~4 A; the distance only fails at 50-60 ns'
  - 'at 70-100 ns the angle RECOVERS (26-36 deg, 45-62% within criterion) while the distance stays 6.5-7.7 A -- aligned but far, so the angle alone is not a criterion either'
  - 'over the whole run distance-only and distance+angle give the SAME 8.0%, because whenever it is close it is also aligned; the angle adds nothing to the SCORE and everything to the MECHANISM'
  - 'the 10 ns adaptive sweep had it holding to the cap at 38.3% engagement (common window) -- the strongest combination in the campaign at the time'
  - 'angle drift over the 10 ns holders: m239 was FLAT (14.9 -> 15.3 deg, p = 0.38), so its 10 ns sweep gave no warning of a rotation that happened at 30-40 ns; but t4_27141ff6418f_m113 drifts 28.6 -> 63.6 deg (p < 1e-4) and t4_75642786dc33_m143 47.9 -> 66.0 deg WITHIN their 10 ns, so the drift is detectable for some modes at sweep length'
  - 'the top-ranked mode t4_4b5557fe1591_m140 sits at 59.7 -> 64.3 deg throughout -- near the isotropic expectation -- while scoring 52% on the distance-only criterion'
  - 'across 1,597 rows the angle removes 50% of engaged frames; modes above 30% fall from 17 to 3; rank correlation between the two orderings is rho +0.827, so they agree broadly and disagree enormously in places'
  - 'the sharpest disagreement: t4_2bd5ba0aa666_m187 is #5 of 1,597 on distance-only (33.9%) and #1257 with the angle (0.0%) -- close to Cys113 and never once in attack geometry. t4_2ae730ebea05_m70 goes #6 -> #398, t4_bcfa4048fd32_m161 #6 -> #213'
runbook: null
---

# NO GO, and the reason is orientation rather than residence

@twu383, on reading the 100 ns result: *"no good, it rotates"*.

## The two readouts disagree, and both are right

| | |
|---|---|
| ligand RMSD | mean **0.352 nm**, max 0.744, PBC guard clean |
| warhead engagement | **8.0%** of the run, median distance 4.46 A |

The pocket holds the molecule for 100 ns -- this is the best residence figure the
campaign has produced, better than sulfopin's 0.803 nm reactant form on the same
receptor (D0107) -- while the warhead turns away from Cys113 and does not come
back. A low ligand RMSD with the reactive carbon 7-8 A from the sulfur is a
molecule that has settled into a comfortable, **non-reactive** pose.

`mdprio_report` renders it as *"Held, residence 1.000"*. That is true of the
ligand and must not be read as a hit. It is the case the two readouts exist to
separate, and the first time this campaign has caught it.

## The rotation LEADS the departure by about 20 ns

This is the part worth keeping.

| block | median distance | median off-normal | within 30 deg |
|---|---:|---:|---:|
| 0-30 ns | 3.7 A | **18 deg** | 88-95% |
| 30-40 ns | 4.0 A | **29 deg** | **55%** |
| 40-50 ns | 3.9 A | 28 deg | 55% |
| 50-60 ns | 5.7 A | 54 deg | 28% |
| 60-70 ns | 8.8 A | 65 deg | 0% |

The angle halves at **30-40 ns while the distance is still ~4 A**. The distance
does not fail until 50-60 ns. Whatever lets go, lets go of the orientation first.

**And the angle alone is not a criterion.** At 70-100 ns it recovers to 26-36 deg
while the distance sits at 6.5-7.7 A: aligned, and far. Both terms are needed,
and neither is sufficient.

## What this does and does not change about the gate

**It does not reverse D0111.** Over the whole run, distance-only and
distance-plus-angle give the *same* 8.0%, because whenever this warhead is
inside 3.5 A it is also inside 30 deg. The angle adds nothing to the SCORE. D0111
dropped it from the gate on 1.2 ns evidence and that conclusion survives.

**What it changes is what the angle is FOR.** It is not a second filter on the
same frames; it is the earlier-moving quantity, and its TREND over a long run is
a leading indicator of departure that the distance does not provide. A 1.2 ns
sweep cannot see it -- at 1.2 ns this molecule was 95% aligned and looked
excellent, which is exactly what the sweep reported.

That suggests a readout nobody has built: not occupancy at any instant, but
whether alignment is DEGRADING across a run. It costs nothing extra -- the angle
is already measured on every frame of every sweep -- and it is the natural thing
to test on the adaptive runs, where 167 modes now reach 10 ns.

## Status of the elevation queue

Two 100 ns runs, two distinct failure modes, neither supporting synthesis:

* `t4_215b12bd9b34_m184` -- **left the site** (mean RMSD 1.854 nm, 27% of frames
  with zero protein contact).
* `t4_88ae42890e2d_m239` -- **held and disengaged**, this record.

`t4_4b5557fe1591_m140` (51.7% engagement over the common window, held the full
10 ns cap) is the strongest remaining candidate. Before it is elevated it is
worth asking the question this record raises: what did its ANGLE do over its
10 ns, and was it already drifting?


---

## Addendum: the angle reorders the top, and that is a live question for D0111

Checked immediately after this record was written, because the obvious next
move was to elevate `t4_4b5557fe1591_m140` and this record says to look at its
angle first. It was worth looking.

**`m140` is not well aligned at all.** Its off-normal angle runs 59.7 -> 64.3 deg
across its 10 ns -- near the isotropic expectation for a randomly oriented
warhead -- while it scores **52%** on the distance-only criterion and sits top of
the sweep rail. Close, and pointing the wrong way.

Across all 1,597 scored rows:

| | distance only | with the angle |
|---|---:|---:|
| mean engagement | 5.2% | 2.6% |
| modes above 30% | **17** | **3** |

Rank correlation is **rho +0.827** -- the two orderings agree broadly and
disagree enormously in individual cases:

| mode | distance-only | with angle |
|---|---|---|
| `t4_2bd5ba0aa666_m187` | 33.9% (#5) | **0.0% (#1257)** |
| `t4_2ae730ebea05_m70` | 33.1% (#6) | 3.3% (#398) |
| `t4_bcfa4048fd32_m161` | 33.1% (#6) | 6.6% (#213) |
| `t4_7b02d1dc4fd2_m50` | 46.9% (#1) | 29.9% (#4) |

`m187` is the one to look at: **fifth of 1,597 on distance, and never once in
attack geometry.**

### What this means for D0111

D0111 dropped the angle from the gate on the grounds that at a 3.0 A cutoff it
moved the count of discriminating modes from 1 to 0 -- a statement about the
COUNT, which was true and is still true. It was not a statement about the
ORDERING, and the ordering is what the rail is for.

The D0110 objection stands and has not gone away: on a distance-selected set the
angle is not class-neutral, because BDHI's off-normal collapses to ~9 deg inside
3 A for steric reasons while acrylamide's does not. Putting the angle back into
the ranking would re-weight the campaign towards BDHI for a reason that is not
reaction competence.

So this is genuinely a decision and not a defect: **rank on distance and lose
molecules like `m187` into the top ten, or rank on both and take a known
class bias.** Both figures are already on every row
(`frac_attack_ready`, `frac_attack_ready_angle`), so the choice can be made
without recomputing anything, and either ordering can be shown.

Not decided here. Flagged because the next elevation was about to be chosen off
the distance-only ordering.
