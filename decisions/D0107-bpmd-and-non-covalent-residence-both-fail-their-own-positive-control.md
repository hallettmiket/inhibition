---
id: D0107
title: BPMD and 100 ns non-covalent residence both fail their own positive control, and neither may be used to rank covalent candidates
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: implementation
supersedes: []
superseded_by: null
affects:
  - shared/bpmd.py
  - scripts/md_residence_3ikd.py
  - docs/how_this_project_breaks.md
evidence:
  - 'BPMD at the 10 ns protocol: 7 of 7 runs ever completed have escaped = True, INCLUDING xtal:6VAJ:QT7, sulfopin in its own crystal pose'
  - 'frac_in_window: sulfopin 0.103 (n=1); t4_80fbed3bdf1e 0.080 +/- 0.045 (n=3); t4_a4908f1fc7b0 0.069 +/- 0.039 (n=3) -- candidates within ~0.5 SD of the positive'
  - 'every run that did NOT escape was 100-300 ps, too short for anything to leave; `escaped` is length-dependent and no 10 ns negative control exists'
  - 'escape occurs at 0.26-1.04 ns with bias_at_exit_kj = 0.0 in 3 of 6 replicates -- the crossing is thermal, so no barrier was measured'
  - '`escaped = any(d >= UNBOUND_NM)` is a TOUCH-ONCE test over 10 ns with no requirement to stay gone; t4_80fbed3bdf1e rep1 spends 33% of its post-crossing time back inside 0.6 nm'
  - 'max_cv_nm clusters at 1.599-1.659 nm across all 6 replicates because UPPER_WALLS sits at WALL_NM = 1.5 -- that number is the restraint, not the molecule'
  - 'the CV is warhead-to-SG distance, not ligand displacement: the same poses hold mean ligand RMSD 0.314 nm with 80% residence over 100 ns while the warhead crosses 1.0 nm'
  - '100 ns non-covalent residence, rx_6VAJ (sulfopin reactant form, 3IKD, same protocol): mean ligand RMSD 0.803 nm, max 5.681, final 3.797, residence 0.743'
  - 'the same measure on t4_80fbed3bdf1e_m10: mean 0.314, max 0.790, final 0.212, residence 0.803 -- the CANDIDATE is better than the positive on every column'
  - 'all three classify HELD, UNSTABLE under residence_tier (optimal < 0.45 nm, residence floor 0.95)'
  - 'BPMD start distances disagree with the recorded pose: start_distance_A = 3.338 while the first COLVAR frame is 4.63 / 7.00 / 6.83 A across replicates'
runbook: null
---

# Both physics readouts fail their positive control

Asked whether `t4_80fbed3bdf1e` survives BPMD, the answer looked like "no": all
three replicates escaped. **The crystal pose of sulfopin also escapes**, and so
does every 10 ns BPMD run this project has ever completed — 7 of 7. The
candidates sit within about half a standard deviation of the positive on
`frac_in_window`.

The same thing happened one tier up. The 100 ns non-covalent control had never
been run; D0072 called it *"the cheapest missing control"* and it stayed missing
for four weeks. Run now on 3IKD with sulfopin's reactant form, **the positive is
the worst of the three**: mean ligand RMSD 0.803 nm against the candidate's
0.314, ending 3.8 nm away and not returning.

## Why the numbers looked usable

**`escaped` is a legal, plausible boolean that answers a different question.**
It is `any(d >= 1.0 nm)` across 10 ns — a touch-once test with no requirement to
stay gone. The ligands cross within the first nanosecond with **zero accumulated
bias in half the replicates**, so metadynamics had deposited nothing and no
barrier was measured; and they come back, the parent spending 33% of its
post-crossing time back inside 0.6 nm. `max_cv_nm` looked like a measurement of
how far each ligand got and is 1.60-1.66 nm in all six runs because
`WALL_NM = 1.5` is where the `UPPER_WALLS` restraint sits.

**The CV is warhead-to-SG distance, not ligand displacement.** On a flexible arm
the warhead can swing past 1.0 nm while the ligand stays seated — which is
exactly what the unbiased 100 ns run shows for the same pose. The two readouts
were not contradicting each other; they measure different things and were read as
if they measured one.

**And for the 100 ns leg the failure is chemical, not mechanical.** Sulfopin's
potency comes from the covalent bond. Its *non-covalent* residence is genuinely
poor, so a non-covalent trajectory drifting away is the correct result for the
right reason. D0072 said as much in prose — *"a covalent inhibitor does not need
its warhead parked on the sulfur"* — and the pipeline went on ranking covalent
candidates on non-covalent residence anyway.

## Decision

1. **Neither `escaped` nor 100 ns non-covalent ligand RMSD may be used to rank
   or reject a covalent candidate.** They have no demonstrated discriminating
   power: the only ground truth available fails both.
2. **Tier 1 remains the one validated readout.** D0071 measured that warhead
   displacement over 300 ps of unrestrained equilibration separates
   crystallographic binders from generated candidates at p = 0.007, and the REF
   group's median of 0.102 nm is reproduced by molecules outside it (D0108).
3. **`frac_in_window` is the salvageable part of BPMD** — sulfopin ranks highest
   on it, and unlike `escaped` it is not a threshold crossing. It is not yet
   validated, and calling it a readout requires a negative control at 10 ns,
   which does not exist.
4. **`escaped` must not be reported as a verdict** while it is a touch-once test
   whose threshold sits below the wall that caps the CV. If it is kept it should
   require sustained departure.

## What this does NOT say

- **Not that BPMD is the wrong method.** The protocol was never calibrated on
  this target: 4 of 5 convergence replicates of the positive crashed
  (`PLMD::GridBase::index_t`), leaving n = 1, and the surviving replicate is on
  **6VAJ**, the receptor D0059 retired. A properly powered control on 3IKD could
  still establish it.
- **Not that the candidates are good.** They are indistinguishable from the
  positive on readouts that cannot distinguish anything. That is an absence of
  evidence in both directions.

## The habit this cost

`how_this_project_breaks` already carries the rule: *"before trusting a new
pipeline, run it on something whose answer is already written down in
`decisions/`."* BPMD was run and its verdict reported before that check. The
check took one query against a table already on disk and inverted the
conclusion. This is the fourth entry where the defence existed, was written
down, and was not applied.

## Loose end worth its own look

BPMD's recorded `start_distance_A` is the pose as handed in (3.338 Å) while the
first COLVAR frame — after the same 300 ps equilibration — is 4.63 to 7.00 Å
depending on replicate. **The pose BPMD biases is not the pose the column names**,
and a reader comparing `start_distance_A` against `max_cv_nm` is comparing two
different starting points. That is a labelling defect independent of everything
above.
