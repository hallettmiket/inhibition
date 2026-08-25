# Build plan — the next release

*A LIVING DOCUMENT. Written with @tt8804 as each pipeline step is discussed, one
at a time. Sections appear as we settle them; a section marked **DRAFT** has not
been agreed, and a section marked **OPEN** is a question we have not answered.*

*Started 2026-08-17. Companion to
[`pipeline_in_plain_language.md`](pipeline_in_plain_language.md), which
describes the pipeline as it is TODAY. This one describes what we are changing
it to, and why.*

---

## 0. Status, and the version question — **OPEN**

| step | change under discussion | section | state |
|---|---|---|---|
| 1. Generate | — | — | not discussed |
| **2. Dock** | **add PoseBusters as a validity gate; quota on valid poses** | **§1** | **DRAFT** |
| 3. Group into modes | — | — | not discussed |
| 4. Score (near-attack) | — | — | not discussed |
| 5. Rank | — | — | not discussed |
| 6. Shortlist | — | — | not discussed |
| 7. Triage sweep | — | — | not discussed |
| 8. Production MD | — | — | not discussed |
| 9. BPMD | — | — | not discussed |

**This is probably a MAJOR release, and we should say so before we start.** By
the CHANGELOG's own rule, MAJOR means *previously reported numbers are invalid
and must be re-measured*. Changing which poses exist changes the mode split, the
per-mode counts, `viable_fraction`, `enrichment`, the ranking and the shortlist —
so nothing from 3.1.0 stays quotable beside it. 3.1.0 is itself still in
progress, so the numbering (3.2.0 vs 4.0.0) is a decision for @tt8804, not a
consequence of this document.

**Consequence worth stating early:** every change we agree here costs a full
re-screen, and they should therefore be batched into one release rather than
landed one at a time. That is an argument for discussing all nine steps before
running anything.

---

## 1. Docking — PoseBusters as a validity gate — **DRAFT**

### 1.1 What we do now

Each in-scope molecule is docked with AutoDock-GPU at a **fixed `--nrun 500`**.
Whatever poses come back are what the molecule gets. The pose count is therefore
*variable and unremarked*: measured on `nac_v5`, molecules return anywhere from
~418 to ~486 poses from the same 500 runs, and nothing downstream records that
one molecule was scored on 418 observations and another on 486.

No physical validity check is applied to a pose at any point. A pose that
clashes into the protein is scored on its near-attack geometry exactly like one
that does not.

### 1.2 What @tt8804 proposes

> *"Add PoseBusters to this step since it is quite lightweight. Change 500 rounds
> of docking → 500 poses, to running docking checked by PoseBusters to get up to
> 500 PoseBusters-stable poses as a quota."*

So: **the quota moves from runs to valid poses.** Dock until 500 poses pass
PoseBusters, rather than docking 500 times and keeping whatever appears.

Two things this buys:

1. **Physical validity.** Poses that clash into the protein, or that carry
   absurd internal energy, stop entering the scoring population.
2. **An equal denominator.** Every molecule is scored on 500 poses rather than
   on however many its docking happened to return. `viable_fraction` is a
   proportion, and it is currently measured on populations that differ by ~15%
   between molecules with nothing recording it.

The second is arguably the bigger win and was not the stated motivation.

### 1.3 Measured: what PoseBusters actually rejects

`posebusters 0.6.5` is already installed in `dwi_cheminf`. No new environment.

Sample: **1,500 poses**, 25 molecules drawn at random from each of the three
warhead families, run in `dock` mode (ligand + protein, no crystal reference)
against the prepared 3IKD.

| | pass rate |
|---|---:|
| **overall** | **90.60%** |
| acrylamide | 91.6% |
| bdhi_c4 | 88.2% |
| bdhi_c5 | 92.0% |

Per molecule: **min 70%, median 90%, max 100%**; 65 of 75 molecules have at
least one failing pose.

**Only 2 of the 22 checks ever fail:**

| check | failure rate |
|---|---:|
| `minimum_distance_to_protein` | 7.53% |
| `internal_energy` | 2.00% |
| the other 20 | 0.00% |

### 1.4 Why 20 of the 22 checks can never fail here

This is structural and worth writing down, because it decides how much of
PoseBusters is actually doing work.

AutoDock-GPU varies **rigid-body placement and torsion angles only**. Bond
lengths, bond angles, ring conformations and stereochemistry are inherited
unchanged from the input conformer, which meeko/RDKit generated and
energy-minimised before docking began.

So `bond_lengths`, `bond_angles`, `aromatic_ring_flatness`,
`double_bond_flatness`, `all_atoms_connected`, `sanitization` and the rest are
**tests of our conformer generator, not of our docking** — and they are the same
answer for all 500 poses of a molecule. They will pass at 100% until we change
how ligands are prepared.

The two that vary with the pose are the two that fail:
`minimum_distance_to_protein` (where the pose sits relative to the protein) and
`internal_energy` (which torsions can change).

**Implication for the design:** we are not adding "22 physical checks". We are
adding **two**, plus twenty that certify the conformer generator once per
molecule. That is still worth having — but the honest description of this change
is *a protein-ligand clash gate and an internal-energy gate*, and the build
should say so rather than claim broader coverage than it has.

### 1.5 Compute cost

**The gate itself is cheap.** 232 ms/pose, measured sustained over 1,500 poses
(not a warm-up artefact — the rate held from the first batch to the last).

A full screen at the new quota is ~561 × 552 ≈ **310,000 poses ≈ 20 CPU-hours**,
or about **25 minutes wall** on the project's 50-worker cap
(`shared/compute.py`). Negligible against everything else in the pipeline.

**The docking increase is the real cost.**

| | |
|---|---:|
| runs needed for 500 valid poses at 90.6% | **~552** |
| GPU increase, on the mean | **+10.4%** |
| runs needed for the worst molecule sampled (70%) | **714** |
| GPU increase, worst case | **+43%** |

The spread matters more than the mean. **A fixed over-request sized on the
average will undershoot the quota for roughly half the library** — which would
reintroduce the unequal-denominator problem this change is partly meant to fix,
while charging us the extra GPU time anyway.

### 1.6 The risk we have not resolved — **OPEN**

PoseBusters flags a clash when any ligand heavy atom is closer to a protein
heavy atom than **0.75 × the sum of their van der Waals radii**
(`clash_cutoff: 0.75`, `radius_type: vdw`, confirmed in the installed
`dock.yml`).

For the reactive carbon against Cys113's sulfur:

```
0.75 × (1.70 + 1.80)  =  2.625 Å
```

**Our near-attack window starts at 2.8 Å.**

So a pose sitting at the ideal SN2 approach distance is **0.175 Å above being
declared a clash** — and `minimum_distance_to_protein` is precisely the check
doing three quarters of the rejecting. Compounding it: we dock into a
*reactive* receptor whose van der Waals parameters were deliberately softened
(`R_EQ_12 = 3.2`) to pull the warhead toward the sulfur. We built a receptor
that encourages close approach and are now proposing a filter that penalises it.

**If PoseBusters rejects attack-ready poses at a higher rate than ordinary
ones, this gate spends GPU time buying reaction-incompetent poses.** That is
the same shape as D0082, where a filter rejects the one molecule known to work.

**Measurement in flight.** All poses of the 18 molecules with the most viable
poses, joined to the screen's own `viable` flag by the established
(mode, pose_idx) correspondence, cross-tabulated against the PoseBusters
verdict. Result goes in §1.7.

### 1.6a A data defect this measurement surfaced — already fixed in code

The join above refused two molecules for a pose-count mismatch (461 SDF vs 457
rows; 474 vs 456), so I measured how widespread it is:

> **105 of 561 molecules (18.7%) have a persisted pose cloud that does not match
> their own scored rows.** Differences run from −42 to +39 poses; the SDF holds
> more in 43 cases and fewer in 62.

For those molecules the cloud cannot be joined to its own measurements at all,
which is #44's rule — *the cloud must come from the same run that produced the
scores* — broken in practice.

**The cause is found and already fixed on this branch** (`a7cca45`,
`nac_screen_v2.py`): the all-poses SDF was written under `if not
adest.exists()`, so a re-screened molecule kept the **previous** run's cloud
beside the current run's table. It is now always rewritten.

**But the `nac_v5` data on disk predates the fix**, so the 18.7% stands for
every number currently derived from those clouds. It resolves only on a
re-screen — which this build requires anyway. Worth stating explicitly because
it means *the current pose clouds are not a safe baseline to compare the
PoseBusters gate against*: §1.10's criterion 4 must be measured on a re-screen
with the fix in, not against what is on disk today.

### 1.7 Result of the attack-ready cross-tabulation — **PENDING**

*To be filled in.*

Decision rule, **set before the number is known**:

* **If the odds ratio is ≈ 1** (PoseBusters rejects attack-ready and ordinary
  poses at indistinguishable rates) — adopt the gate as proposed, all checks on.
* **If attack-ready poses fail materially more often** — do **not** adopt
  `minimum_distance_to_protein` as-is. Two options, in preference order:
  1. exempt the anchor atom pair (reactive atom ↔ Cys113 SG) from the clash
     check and keep the check for every other atom, since the close approach we
     want is *specifically* at that pair;
  2. drop `minimum_distance_to_protein` from the gate entirely and keep the
     other 21 — which, per §1.4, means keeping `internal_energy` and twenty
     conformer certifications.
* **If attack-ready poses fail materially LESS often** — that is evidence the
  criterion and the clash check agree, and worth recording as a small
  independent validation of both.

### 1.8 Implementation options — **DRAFT**

The quota loop interacts with two known AutoDock-GPU properties.

**(a) The `--nrun` ceiling.** Measured on this build: `--nrun 5000` corrupts its
own stack *and still writes a `.dlg`*; `--nrun 10000` exits 0 with no output.
`dock()` already refuses all three failure shapes. Any quota loop must therefore
cap the per-call `--nrun` well below 2,000 and top up with additional calls
rather than one large request.

**(b) Seeding.** Docking is seeded for reproducibility (#77) — the same seed
returns the same poses, so a top-up call must use a *different* seed or it
returns the identical cloud. The seed sequence has to be deterministic and
recorded, or we lose the reproducibility that #77 bought.

Three shapes:

| option | how | cost | risk |
|---|---|---|---|
| **A. Flat over-request** | one call at `--nrun ceil(500 / 0.906) = 552`, keep the first 500 valid | simplest, one call | undershoots for ~half the library — the molecules with the *worst* pose quality get the *smallest* populations, which is backwards |
| **B. Iterate to quota** | call at 500, gate, top up in batches of ~100 with successive seeds until 500 valid or a hard cap | exact quota | more calls, each with AutoDock start-up cost; needs the seed ledger |
| **C. Adaptive** | call at 500, measure this molecule's pass rate, size one top-up call from it | one extra call, near-exact | slight overshoot/undershoot; still needs the seed ledger |

**Recommendation: C**, with B's hard cap as a backstop. It costs one extra
AutoDock call per molecule, hits the quota for nearly everything, and its
failure mode (a molecule that still falls short after the top-up) is *recorded*
rather than silent — see §1.9.

### 1.9 What must be recorded — **DRAFT**

The current screen records `nrun: 500` as a constant. Under a quota it becomes a
per-molecule outcome, and a constant in the config would then be a lie.

Per molecule, the run must persist:

- `n_runs_requested` — total AutoDock runs actually spent
- `n_poses_returned` — before the gate
- `n_poses_valid` — after the gate
- `quota_met` — bool; **false is a result, not an error**
- `pb_fail_reasons` — counts per failing check
- `seeds` — the seed sequence used, for reproducibility

**A molecule that cannot reach 500 valid poses must be stamped, not dropped and
not silently short.** That is the same rule as `docked_species_ok` (#58): the
substitution we are guarding against is a molecule quietly entering the ranking
on a smaller population than everything it is ranked against.

### 1.10 Pre-registered failure criteria

This change is **rejected** if any of the following is true after
implementation:

1. **It selects against reaction competence.** Attack-ready poses are rejected
   materially more often than ordinary poses and the fix in §1.7 does not
   resolve it.
2. **The quota is not actually met.** More than 10% of molecules end the screen
   below the 500-valid-pose quota, leaving the unequal denominators the change
   was meant to remove.
3. **It costs more than +25% GPU on the library mean.** At that point the same
   GPU time is better spent on the sweep, where §"What we can and cannot claim"
   in the plain-language doc says the real uncertainty lies.
4. **It changes nothing measurable.** If the ranking produced with the gate is
   indistinguishable from the ranking without it (rank correlation ≳ 0.99 on the
   same seeds), then we have spent GPU time and added a dependency for no
   change in what the pipeline recommends — and the honest thing is to record
   that and not adopt it.

Criterion 4 deserves emphasis. **It is the most likely outcome**: 90.6% of poses
already pass, the rejects are concentrated in a check that may be firing on the
poses we want, and mode assignment is a clustering over hundreds of poses that
is unlikely to move much when 9% are removed. We should measure it explicitly
rather than assume the gate helped because it sounds like it should.

---

## 2-9. Remaining steps — **not yet discussed**

Placeholders, so the order of work is visible. Each gets its own section, with
the same structure: what we do now, what changes, what it costs, what would make
us reject it.

---

## Appendix — measurements behind this document

| what | where |
|---|---|
| PoseBusters pass rate, 1,500 poses, 3 families | §1.3; raw CSV in the session scratchpad |
| per-check failure breakdown | §1.3 |
| 232 ms/pose timing | §1.5 |
| clash threshold 2.625 Å vs the 2.8 Å window | §1.6, from the installed `dock.yml` |
| attack-ready cross-tabulation | §1.7, pending |

**Reproduce:** the sampling scripts are in the session scratchpad and should be
promoted to `exp/` before this document is used to justify a decision. A
measurement that cannot be re-run is a claim, not a result.
