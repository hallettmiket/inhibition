---
id: D0111
title: Attack ready is occupancy within 3.5 A of Cys113, distance only, and the 100 ns gate is that occupancy plus a pose-stability test
date: 2026-09-02
status: accepted
approach: shared
decided_by: '@twu383'
origin: user
supersedes: []
superseded_by: null
affects:
  - config/target.yaml
  - scripts/attack_sweep.py
  - scripts/recompute_attack_ready.py
  - scripts/sweep_combine.py
  - scripts/sweep_gui_refresh.py
  - shared/sweep_state.py
  - tests/test_elevation_gate.py
evidence:
  - 'the readout used `nac_criterion.NAC_DIST_MAX` (4.2 A), the near-attack WINDOW, while the worklist selected modes at < 3.0 A -- so a mode at a trajectory median of 3.6 A scored 93% attack-ready'
  - 'measured over 93 completed 1.2 ns sweeps, modes scoring >1% ready: distance-only 1/15/39/64/78 at 3.0/3.2/3.4/3.6/4.2 A; distance+angle 0/7/31/46/58'
  - 'the angle is not what limits the answer: at 3.0 A it moves the count from 1 to 0, and median closest approach across the 93 is 3.37 A'
  - 'at 3.5 A distance-only over 98 sweeps: best occupancy 33.9%, median 2.1%; 0 of 98 reach the 60% bar; 71 of 98 pass the RMSD half'
  - '`sweep_combine` built its pose-rank map from the worklist `ident` (the MOLECULE) and looked it up with a mode id: 0 of 98 finished modes resolved, 18 of 98 were reading another pose trajectory, median disagreement 1.12 A, max 21.65 A, 7 flipped the pose-held verdict in both directions'
  - '`sweep_state.state()` carried a hand-listed allowlist of nine result columns, so `rmsd_max_a`, `rmsd_mean_a`, `elevate` and `pose_held` never reached the page'
  - '`sweep_gui_refresh.n_finished` summed `ok` rows per file without deduping on ident: a re-scored table took the count from 98 to 196 with nothing new run'
runbook: null
---

# What "attack ready" means, and what earns a 100 ns run

## 1. The readout disagreed with the selection

The nac_v8 worklist selects modes whose warhead sits under **3.0 A** from
Cys113 SG. The readout that graded them, `frac_attack_ready`, used
`nac_criterion.NAC_DIST_MAX` = **4.2 A** -- the near-attack *window*, a
different quantity that happens to be a distance in the same units.

Both halves of one campaign therefore meant different things by "close", and
the readout was the looser of the two. A mode whose warhead sat at a trajectory
median of 3.6 A scored **93% attack ready**: true under that definition, and
not what anyone reading the number took it to mean.

@twu383, 2026-09-02: *"we need to clearly establish that attack ready means
within angle and under 3 A at any time"*, then, after seeing what 3 A costs,
*"okay 3.5 then"* and *"just distance"*.

## 2. Why the angle came out of the gate

Measured across 93 completed trajectories before choosing, counting modes that
score above 1% ready:

| cutoff | distance only | distance + angle |
|---:|---:|---:|
| 3.0 A | 1 / 93 | 0 / 93 |
| 3.2 | 15 | 7 |
| 3.4 | 39 | 31 |
| 3.6 | 64 | 46 |
| 4.2 | 78 | 58 |

**The distance is the binding constraint, not the angle.** At 3.0 A the angular
term moves the count from 1 to 0 -- it removes modes without changing the
conclusion. Median closest approach across the 93 is 3.37 A, so these molecules
simply do not sit inside 3 A once dynamics relax them.

It is also the term **D0110** showed is *not class-neutral* on a
distance-selected set: BDHI's off-normal angle collapses from 68 to 9 degrees
inside 3 A for steric reasons (the reactive carbon is inside a rigid
dihydroisoxazole ring), while acrylamide's barely moves. Gating on it re-weights
the campaign towards BDHI for a reason that is not reaction competence.

**The angle is still measured.** Every row carries `frac_attack_ready_angle`
beside `frac_attack_ready`, and `attack_ready_uses_angle` says which one the
gate used. Dropped from the gate, not from the record.

## 3. The 100 ns gate

@twu383: *"<3.5 max RSMD or mean RMSD 3.0 (need to account for quick spikes but
overall low rmsd) + 60%+ warhead within 3.5 A goes to 100 ns"*.

Two conditions, **both** required:

| | test | why it is separate |
|---|---|---|
| pose held | `rmsd_max_a < 3.5` **OR** `rmsd_mean_a < 3.0` | the whole LIGAND stayed where it was docked |
| warhead engaged | `frac_attack_ready >= 0.60` within 3.5 A | the WARHEAD stayed in reach of Cys113 |

They are independent on purpose, because they disagree: a molecule can pivot its
warhead into place while its scaffold walks off, and it can sit rock-still in
the wrong orientation. Either condition alone elevates one of those.

The **OR** is the spike allowance that was explicitly asked for: a run touching
4 A for 20 ps and sitting at 2 A otherwise passes on the mean; one that never
exceeds 3.5 passes without needing the mean. 3.5 A = 0.35 nm, the same bar the
triage sweep has used since D0087.

**"At any time" was rejected in favour of occupancy.** Of the three modes that
ever reached under 3.0 A, two did so for exactly one frame out of 121. A single
qualifying frame is noise, and a gate that a single frame can pass is not a
measurement of residence.

### What it currently elevates: nothing

Over the 98 sweeps finished when this was written:

* best occupancy at 3.5 A **33.9%**, median **2.1%**
* **71 of 98** pass the RMSD half -- the poses mostly hold
* **0 of 98** reach 60% engagement
* therefore **0 of 98 elevate**

That is a real reading of ~2% of a 4,295-mode campaign, not a calibration
failure. The bar was set from the chemistry and whether anything clears it is
the experiment. It is recorded here so that if the bar is later moved, it is
moved knowingly rather than because nothing passed.

## 4. Three silent paths that fed the gate, all fixed

Found while wiring the gate, and the first is the one that mattered.

### `sweep_combine` was reading the wrong trajectory for 18% of modes

The block that computes the page's RMSD carries a comment saying **"POSE RANK IS
PART OF THE KEY"** -- written for catalogue #23, with the numbers that motivated
it. It then built the map from the worklist's `ident` column, which is the
**molecule** (`t4_215b12bd9b34`), and looked it up with a **mode**
(`t4_215b12bd9b34_m184`).

**0 of 98 finished modes resolved a pose rank.** Every lookup returned None and
`rep_dir` fell back to whichever sibling directory sorted first:

* 18 of 98 modes were reading another pose's trajectory
* median disagreement **1.12 A**, max **21.65 A**
* **7 flipped the pose-held verdict**, in both directions -- `t4_5f904942f66c_m28`
  is 13.09 A and was shown as 4.39; `t4_2bd5ba0aa666_m187`, the best mode in the
  campaign, is 1.89 A and was shown as 3.78

**Why it looked right.** The worklist carries both `ident` and `task_id`, both
populated, both plausible string ids. Only 18% differed because most molecules
have a single swept mode so far -- the error grows as more modes per molecule
finish. And the fix for this exact defect had already been written once; it was
keyed on the wrong column and so did nothing, under a comment asserting it did.

`sweep_assets` and `recompute_attack_ready` already keyed on `task_id`. This was
the last reader that did not.

### `sweep_state` dropped every new measurement on the way to the page

`state()` selected result columns through a hand-maintained allowlist of nine
names, so `rmsd_max_a`, `rmsd_mean_a`, `elevate` and `pose_held` never arrived.
That is catalogue #5, and it is *why* the defect above could bite: the page had
to recompute RMSD by a second route because the first route silently withheld
it. Two implementations of one quantity is how a page and a gate come to
disagree while both look right.

Now derived from `res.columns`, minus what `base` already provides.

### `n_finished` counted rows, not modes

`sweep_gui_refresh.n_finished` summed `ok` rows per file with no dedupe on
ident. Re-scoring rewrites every finished row into a new versioned file, which
took the count from 98 to **196** with nothing new run. Only a rebuild trigger,
so it cost nothing -- but it is a completion count that reads as progress and it
was wrong by a factor of two.

## 5. Re-scoring costs no GPU

`recompute_attack_ready.py` re-reads each mode's persisted `sweep_dense.pdb` and
`rmsd.xvg` and writes a **new versioned table**. Old rows are never edited: the
outputs root is append-only, and a rewritten row would make two definitions
indistinguishable in one file. Every row carries `attack_ready_max_a`,
`attack_ready_uses_angle` and the three `elevate_*` thresholds, so which
definition produced a number is a property of the row.

## 6. Guards

`tests/test_elevation_gate.py`, 10 tests, each verified to fail under a
deliberate mutation:

* both conditions required; neither elevates alone
* the spike allowance is an **OR** -- turning it into an AND fails the suite
* a **missing** reading never elevates (catalogue #30, #31: the guard that
  passes because what it inspects is absent)
* RMSD is reported in **Angstrom** -- GROMACS writes nm, and a missing factor of
  10 makes every run pass by a mile and read as a spectacular result
* an absent `rmsd.xvg` is absent, not zero
* every row records the definition that produced it
* the pose-rank map resolves **mode** ids, and the molecule-keyed version is
  asserted to resolve nothing -- so the test cannot pass vacuously
* `rescols` is derived from `res.columns`, not listed by hand


---

## 7. The green zone had to move with the gate

*Added 2026-09-02, @twu383: "update the rmsd plots to show the correct green
zones".*

Three drawers shaded `NAC_DIST_MIN..NAC_DIST_MAX` (2.8-4.2 A) -- the SCREEN's
near-attack window -- while the sweep is judged at 2.8-3.5:

| | was | now |
|---|---|---|
| `sweep_assets` figure | 2.8-4.2 "attack range" | 2.8-3.5 "attack ready" |
| `mdprio_report` figure | 2.8-4.2 "attack window" | 2.8-3.5 |
| `md_movie` viewer readout | `nac_lo=2.8, nac_hi=4.2` defaults | read from the criterion |

**A trace could sit inside the green zone for most of a run and score 0% engaged
directly beneath it.** The viewer was worse: `nac_lo`/`nac_hi` were keyword
defaults no caller ever overrode -- catalogue #32/#35 -- so it printed "(in
attack window)" for a pose at 4.0 A that the same page ranked as not engaged.

`shared/nac_criterion.attack_ready_window()` is now the single source. It is
deliberately SEPARATE from `NAC_DIST_MIN`/`NAC_DIST_MAX`, which are unchanged:
the screen's criterion and the sweep's bar are different quantities, and
collapsing them would have silently re-scored the screen. The wider window is
still drawn, as a dotted hairline, so the two pages do not describe different
physics.

Guarded by `test_the_plotted_band_is_the_gate_band` (the band must equal
`attack_ready_max_a()`, and `NAC_DIST_MAX` must still be 4.2) and
`test_no_plot_hardcodes_the_wider_window_as_the_ready_band`.

**And the refresher could not see the change.** Its staleness check compared each
page against `<ident>.pdb` only, but the figure is a base64-embedded
`<ident>.png` -- so regenerating 114 figures changed no `.pdb` and 113 pages
would have gone on serving the old band with nothing reporting anything wrong.
`_asset_mtime` now takes the newest of both. A staleness check scoped to one of
two inputs is a guard that cannot fail for the other.
