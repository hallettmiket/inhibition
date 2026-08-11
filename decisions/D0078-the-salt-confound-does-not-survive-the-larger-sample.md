---
id: D0078
title: The salt confound does not survive the larger sample, and the near-independence of engagement and residence is reversed
date: 2026-08-11
status: proposed
approach: shared
decided_by: '@tt8804'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/gromacs_explicit.py
  - docs/gui_spec.md
  - scripts/mdprio_combine.py
evidence:
  - 'ion count is not the variable: `addions cpx Na+ 0` / `Cl- 0` neutralises only, so NO system has salt; every box is at ~0 M, not 0.15 M'
  - 'n=12 (2026-08-07): cations left 3/3, neutrals 0/9, Fisher p = 0.0045 — perfect separation'
  - 'n=58 (2026-08-11): cations left 21/32 (66%), neutrals 11/26 (42%), OR 2.60, Fisher p = 0.112 — NOT significant'
  - 'ligand charge and ion count are collinear by construction: 28/28 neutrals got 1 ion, 29/36 cations got 0'
  - 'the 7 cations that DID get an ion left 7/7 (100%) vs 14/25 (56%) for cations with none — direction OPPOSITE to the hypothesis, Fisher p = 0.066'
  - 'engagement vs residence_frac recomputed from 62 runs on rmsd.xvg at the projects own BOUND_NM = 1.2: rho = +0.684, p = 9.15e-10'
  - 'engagement vs max ligand RMSD, same 62 runs: rho = -0.732, p = 1.43e-11'
  - 'within cations rho(engagement, max RMSD) = -0.837; within neutrals -0.609 — the relation is not a charge artefact'
  - 'the cited rho = -0.007 appears in docs/gui_spec.md attributed to #46; #46 does not contain that number anywhere'
---

# The blocker was measured at n=12 and does not reproduce

#33 was raised by the `adversary` audit on 2026-08-07 and marked a blocker on the
weekend sweep. Its table is exact and its arithmetic is right:

| molecule | ions in box | outcome |
|---|---|---|
| `t4_9265b4bff789` | 0 | came unbound |
| `t4_9a973be6b946` | 0 | came unbound |
| `t4_da2e98512d02` | 0 | came unbound |
| all nine others | 1 | stayed bound |

Perfect separation, Fisher p = 0.0045. On twelve molecules.

We now have 58 completed 100 ns runs that can be matched to a charge class. The
separation is gone:

| charge at pH 7.4 | n | left | rate |
|---|---:|---:|---:|
| cation | 32 | 21 | 66% |
| neutral | 26 | 11 | 42% |

Odds ratio 2.60, **Fisher p = 0.112**. A real-looking direction, not a
significant one, and nothing like a blocker.

## Ion count was never the variable

`shared/gromacs_explicit.py` builds every system with

```
addions cpx Na+ 0
addions cpx Cl- 0
```

The `0` means *add exactly enough to neutralise*. **No system in this project has
ever contained added salt.** A box with one Na+ is not "salted" and a box with
zero is not unusually bare — both are at ~0 M against a physiological 0.15 M. The
difference the issue's table turns on is one ion in roughly eight thousand
waters, which screens nothing.

The issue says this itself, in a paragraph beginning "Separately:". That
paragraph is the real finding and the table is not.

## The two variables cannot be separated in this dataset, but the sign is wrong

Ligand charge and ion count are collinear by construction: a +1 ligand cancels
the receptor's −1, so tleap adds nothing.

| | 0 ions | 1 ion |
|---|---:|---:|
| neutral | 0 | 28 |
| cation | 29 | 7 |

Which is why **re-running the three cations with salt cannot settle it** — the
experiment #33 proposes changes salt and charge together, exactly as the original
comparison did.

The 7 cations that happened to receive an ion are the only place the two come
apart, and they go the wrong way:

| cations | n | left |
|---|---:|---:|
| 0 ions | 25 | 14 (56%) |
| 1 ion | 7 | **7 (100%)** |

Fisher p = 0.066, n = 7, so this is not proof of anything. But it is not evidence
*for* the ion hypothesis, and it is the only within-charge contrast available.

## What this does to the claim #33 said it blocked

#33's stated stake:

> The entire reason we built the short pre-screen is this claim: *how long a
> molecule stays put and whether it is aimed correctly are nearly unrelated.*

Recomputed on the project's own quantities — `explicit_frac_frames_engaged`
against `residence_frac` taken from each run's `rmsd.xvg` at `BOUND_NM = 1.2`,
62 runs:

| pair | rho | p | n |
|---|---:|---:|---:|
| engagement vs residence_frac | **+0.684** | 9.15e-10 | 62 |
| engagement vs max ligand RMSD | **−0.732** | 1.43e-11 | 62 |
| within cations, engagement vs max RMSD | −0.837 | 1e-8 | 32 |
| within neutrals, engagement vs max RMSD | −0.609 | 0.001 | 26 |

`docs/gui_spec.md` justifies the **combined / split held-left** toggle with
"engagement and residence are near-independent (ρ = −0.007, #46)". At n = 62
they are **strongly related in the opposite direction**, and the relation holds
inside each charge class, so it is not a charge artefact.

The cited −0.007 could not be traced. It appears in `gui_spec.md` attributed to
#46; **#46 does not contain that number**. A number in the spec with no reachable
derivation is the same defect this project keeps finding in its data paths, moved
into prose.

## Decision

1. **#33 is not a blocker.** The finding it rests on is a twelve-sample
   coincidence that does not reproduce at 58.
2. **The genuine defect it found is real and larger than it claimed**: no
   simulation in this project has physiological ionic strength. That affects all
   63 runs, not three, and it is a setup deficiency to fix going forward rather
   than a confound that invalidates a specific comparison.
3. **The salt experiment, if run, must vary salt alone** — a matched set of
   cations *and* neutrals at 0.15 M NaCl. Re-running only the three cations
   changes both variables again.
4. **The toggle's justification in `gui_spec.md` is withdrawn.** The toggle may
   still earn its keep — a molecule that engages 80% of the run and still leaves
   is worth seeing separately — but not on the grounds stated, and the number is
   corrected in the doc.

## What is NOT concluded

That cations behave the same as neutrals. 66% versus 42% at p = 0.112 is a
direction worth watching, and with `charge_ph74` already stamped on D4 it costs
nothing to keep watching. It is simply not established, and it was never the
ion count.
